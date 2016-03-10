/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <utility>
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <cstdlib>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>

#include "graph.hh"
#include "dimacs.hh"

namespace po = boost::program_options;

using std::make_tuple;
using std::max;
using std::set_union;
using std::move;
using std::vector;
using std::atomic;
using std::mutex;
using std::unique_lock;
using std::thread;
using std::condition_variable;
using std::cv_status;
using std::string;

using std::cout;
using std::endl;
using std::boolalpha;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::seconds;
using std::chrono::steady_clock;
using std::chrono::time_point;

using boost::dynamic_bitset;

struct Stats
{
    atomic<unsigned long long> nodes_created{ 0 };
    atomic<unsigned long long> nodes_reused{ 0 };
    atomic<unsigned long long> bdds_created{ 0 };
    atomic<unsigned long long> bdds_pruned{ 0 };
};

struct Node
{
    unsigned long long score; // number of accepted vertices so far (i.e. path length)
    dynamic_bitset<> clique;  // vertices accepted so far (i.e. one path of best length)
    dynamic_bitset<> state;   // vertices which remain to be selectable
};

struct Level
{
    vector<Node> nodes;       // children of a particular BDD level

    /**
     * Add node to level, updating stats, merging it if it is a duplicate.
     */
    void add_node(const Graph & graph, Node && node, Stats & stats, const int domination)
    {
        for (auto & d : nodes) {
            if (d.state == node.state) {
                if (node.score > d.score) {
                    d.score = node.score;
                    d.clique = node.clique;
                }
                ++stats.nodes_reused;
                return;
            }
            else if ((domination > 0) && d.score <= node.score && d.state.is_subset_of(node.state)) {
                d = node;
                ++stats.nodes_reused;
                return;
            }
            else if ((domination > 0) && node.score <= d.score && node.state.is_subset_of(d.state)) {
                ++stats.nodes_reused;
                return;
            }

            if (domination > 1)
            {
                auto remainder = node.state - d.state;
                unsigned long long missing_score = 0;
                for (auto v = remainder.find_first() ; v != dynamic_bitset<>::npos ; v = remainder.find_next(v)) {
                    remainder.reset(v);
                    missing_score += graph.weight(v);
                    if (node.score + missing_score > d.score)
                        break;
                }

                if (node.score + missing_score <= d.score) {
                    ++stats.nodes_reused;
                    return;
                }
            }

            if (domination > 1)
            {
                auto remainder = d.state - node.state;
                unsigned long long missing_score = 0;
                for (auto v = remainder.find_first() ; v != dynamic_bitset<>::npos ; v = remainder.find_next(v)) {
                    remainder.reset(v);
                    missing_score += graph.weight(v);
                    if (d.score + missing_score > node.score)
                        break;
                }

                if (d.score + missing_score <= node.score) {
                    d = node;
                    ++stats.nodes_reused;
                    return;
                }
            }
        }

        nodes.push_back(move(node));
        ++stats.nodes_created;
    }
};

struct BDD
{
    unsigned height;          // how many levels in this BDD have been assigned?
    unsigned n_unassigned;    // how many variables are left to assign?

    vector<Level> levels;     // nodes at each level in the BDD

    explicit BDD(unsigned h, unsigned n) :
        height(h),
        n_unassigned(n),
        levels(h + 1)
    {
    }
};

struct Incumbent
{
    atomic<unsigned long long> value{ 0 };

    mutex clique_mutex;
    dynamic_bitset<> clique;

    void update(const Node & node)
    {
        while (true) {
            unsigned long long current_value = value;
            if (node.score > current_value) {
                if (value.compare_exchange_strong(current_value, node.score)) {
                    unique_lock<mutex> lock(clique_mutex);
                    clique = node.clique;
                    cout << "incumbent = " << node.score << endl;
                    cout << "candidate =";
                    for (auto i = clique.find_first() ; i != dynamic_bitset<>::npos ; i = clique.find_next(i))
                        cout << " " << i + 1; // dimacs files are 1-indexed
                    cout << endl;
                    break;
                }
            }
            else
                break;
        }
    }
};

void solve(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats, atomic<bool> & abort,
        const int domination);

/**
 * Select a vertex to branch on, from level. May return -1, if there are no
 * vertices to branch on in the restricted case.
 */
int select_branch_vertex_in_most_states(const Graph & graph, const Level & level)
{
    vector<int> counts(graph.size());
    for (auto & n : level.nodes)
         for (auto v = n.state.find_first() ; v != dynamic_bitset<>::npos ; v = n.state.find_next(v))
            ++counts[v];

    int best = -1;

    for (unsigned i = 0 ; i < counts.size() ; ++i)
        if (counts[i] != 0 && (-1 == best || counts[i] < counts[best]))
            best = i;

    return best;
}

/**
 * Build level of bdd by branching on branch_vertex.
 */
void build_level(const Graph & graph, BDD & bdd, int level, int branch_vertex, Stats & stats, const int domination)
{
    bdd.levels[level].nodes.reserve(bdd.levels[level - 1].nodes.size() * 2);

    // For each node on the previous level...
    for (auto & n : bdd.levels[level - 1].nodes) {
        // Does this node's state allow us to accept the branch vertex?
        if (n.state[branch_vertex]) {
            // Yes, the new level of the BDD gets an accept node.
            Node accept = n;
            accept.score += graph.weight(branch_vertex);
            accept.clique.set(branch_vertex);
            accept.state &= graph.neighbourhood(branch_vertex);
            bdd.levels[level].add_node(graph, move(accept), stats, domination);
        }

        // The new level of the BDD always gets a reject node.
        Node reject = n;
        reject.state.set(branch_vertex, false);
        bdd.levels[level].add_node(graph, move(reject), stats, domination);
    }
}

/**
 * Solve a relaxed form of the BDD, to get an upper bound on its actual best
 * value.
 */
unsigned solve_relaxed(const Graph & graph, BDD & bdd, Stats & stats, atomic<bool> & abort, const int domination)
{
    if (abort.load())
        return 0;

    // Build up each level...
    unsigned level = 1;
    for ( ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex_in_most_states(graph, bdd.levels[level - 1]);
        if (-1 == branch_vertex)
            break;

        build_level(graph, bdd, level, branch_vertex, stats, domination);

        // Is our newly created level too wide?
        auto & level_nodes = bdd.levels[level].nodes;
        if (level_nodes.size() > bdd.n_unassigned + 1) {
            // Yup, start merging together states.
            sort(level_nodes.begin(), level_nodes.end(), [] (const Node & a, const Node & b) {
                    return make_tuple(a.score, a.state.count()) > make_tuple(b.score, b.state.count());
                    });

            Node merged = { 0, dynamic_bitset<>(graph.size(), 0), dynamic_bitset<>(graph.size(), 0) };
            while (level_nodes.size() >= bdd.n_unassigned + 1) {
                Node & n = level_nodes[level_nodes.size() - 1];
                merged.score = max(merged.score, n.score);
                merged.state |= n.state;
                merged.clique |= n.clique;
                level_nodes.pop_back();
            }
            level_nodes.push_back(move(merged));
        }
    }

    // The bound is found by looking at the score of every leaf node and
    // picking the best.
    unsigned long long bound = 0;
    for (auto & n : bdd.levels[level - 1].nodes)
        bound = max(bound, n.score);
    return bound;
}

/**
 * Solve a restricted form of the BDD, in the hopes of getting a good incumbent quickly.
 */
void solve_restricted(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats, atomic<bool> & abort, const int domination)
{
    if (abort.load())
        return;

    int last_level = 0;

    // Build up each level...
    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex_in_most_states(graph, bdd.levels[level - 1]);
        if (-1 == branch_vertex)
            break;

        build_level(graph, bdd, level, branch_vertex, stats, domination);
        last_level = level;

        // Is our newly created level too wide?
        if (bdd.levels[level].nodes.size() > bdd.n_unassigned + 1) {
            // Yes, just drop some bad nodes.
            sort(bdd.levels[level].nodes.begin(), bdd.levels[level].nodes.end(), [] (const Node & a, const Node & b) {
                    return make_tuple(a.score, a.state.count()) > make_tuple(b.score, b.state.count());
                    });

            bdd.levels[level].nodes.resize(bdd.n_unassigned + 1);
        }
    }

    // We may have found a better solution.
    for (auto & n : bdd.levels[last_level].nodes)
        incumbent.update(n);
}

/**
 * Given a level of nodes in a BDD, solve them exactly by creating new BDDs.
 */
void solve_level_by_branching(const Graph & graph, const BDD & bdd, const vector<Node> & nodes,
        Incumbent & incumbent, Stats & stats, atomic<bool> & abort,
        const int domination)
{
    vector<unsigned long long> bounds(nodes.size(), 0);

    for (unsigned n_index = 0 ; n_index < nodes.size() ; ++n_index) {
        const auto & n = nodes[n_index];

        if (abort.load())
            return;

        BDD relaxed_bdd(n.state.count(), bdd.n_unassigned);

        // Start by getting a relaxed bound, to see if it's worth continuing.
        relaxed_bdd.levels[0] = { Level{ { n } } };
        bounds[n_index] = solve_relaxed(graph, relaxed_bdd, stats, abort, domination);

        // Are we any good?
        if (bounds[n_index] > incumbent.value) {
            // Yes, now get a restricted solution.
            BDD restricted_bdd(n.state.count(), bdd.n_unassigned);
            restricted_bdd.levels[0] = { Level{ { n } } };
            solve_restricted(graph, restricted_bdd, incumbent, stats, abort, domination);
        }
    }

    for (unsigned n_index = 0 ; n_index < nodes.size() ; ++n_index) {
        const auto & n = nodes[n_index];

        if (bounds[n_index] > incumbent.value) {
            ++stats.bdds_created;

            BDD child_bdd(n.state.count(), bdd.n_unassigned);
            child_bdd.levels[0] = { Level{ { n } } };

            solve(graph, child_bdd, incumbent, stats, abort, domination);
        }
        else
            ++stats.bdds_pruned;
    }
}

/**
 * Solve a BDD.
 */
void solve(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats, atomic<bool> & abort, const int domination)
{
    if (abort.load())
        return;

    // Build up each level...
    unsigned level = 1;
    for ( ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex_in_most_states(graph, bdd.levels[level - 1]);
        if (-1 == branch_vertex)
            break;

        build_level(graph, bdd, level, branch_vertex, stats, domination);

        // Is our newly created level too wide?
        if (bdd.levels[level].nodes.size() > bdd.n_unassigned + 1) {
            // Yes, solve each of its nodes independently using new BDDs.
            solve_level_by_branching(graph, bdd, bdd.levels[level].nodes, incumbent, stats, abort, domination);
            return;
        }
    }

    // We built a complete BDD. Pick off its leaf nodes for new candidate
    // solutions.
    for (auto & n : bdd.levels[level - 1].nodes)
        incumbent.update(n);
}

BDD create_root_bdd(const Graph & graph, Stats & stats)
{
    BDD bdd(graph.size(), graph.size());

    // Initially, we can pick all vertices.
    dynamic_bitset<> all_vertices(graph.size());
    all_vertices.set();

    // Create BDD with root node.
    bdd.levels[0] = { Level{ { Node{ 0, dynamic_bitset<>(graph.size()), all_vertices } } } };

    ++stats.bdds_created;

    return bdd;
}

int main(int argc, char * argv[])
{
    try {
        po::options_description display_options{ "Program options" };
        display_options.add_options()
            ("help",                                  "Display help information")
            ("timeout",            po::value<int>(),  "Abort after this many seconds")
            ("domination",         po::value<int>(),  "Dominance detection power (0 for none)")
            ;

        po::options_description all_options{ "All options" };
        all_options.add_options()
            ("clique-file",      po::value<std::string>(), "Clique file")
            ;

        all_options.add(display_options);

        po::positional_options_description positional_options;
        positional_options
            .add("clique-file", 1)
            ;

        po::variables_map options_vars;
        po::store(po::command_line_parser(argc, argv)
                .options(all_options)
                .positional(positional_options)
                .run(), options_vars);
        po::notify(options_vars);

        /* --help? Show a message, and exit. */
        if (options_vars.count("help")) {
            std::cout << "Usage: " << argv[0] << " [options] clique-file" << std::endl;
            std::cout << std::endl;
            std::cout << display_options << std::endl;
            return EXIT_SUCCESS;
        }

        /* No input file specified? Show a message and exit. */
        if (! options_vars.count("clique-file")) {
            std::cout << "Usage: " << argv[0] << " [options] clique-file" << std::endl;
            return EXIT_FAILURE;
        }

        int timeout = 0;
        if (options_vars.count("timeout"))
            timeout = options_vars["timeout"].as<int>();

        int domination = 0;
        if (options_vars.count("domination"))
            domination = options_vars["domination"].as<int>();

        bool aborted = false;

        Graph graph = read_dimacs(options_vars["clique-file"].as<std::string>());

        thread timeout_thread;
        mutex timeout_mutex;
        condition_variable timeout_cv;
        atomic<bool> abort{ false };
        if (0 != timeout) {
            timeout_thread = thread([&] {
                    auto abort_time = steady_clock::now() + seconds(timeout);
                    {
                        /* Sleep until either we've reached the time limit,
                         * or we've finished all the work. */
                        unique_lock<mutex> guard(timeout_mutex);
                        while (! abort.load()) {
                            if (cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                                /* We've woken up, and it's due to a timeout. */
                                aborted = true;
                                break;
                            }
                        }
                    }
                    abort.store(true);
                    });
        }

        /* Start the clock */
        auto start_time = steady_clock::now();

        Stats stats;
        Incumbent incumbent;

        {
            BDD bdd = create_root_bdd(graph, stats);

            {
                // Start by building a restricted BDD, just to get a good incumbent quickly.
                BDD restricted_bdd = bdd;
                solve_restricted(graph, restricted_bdd, incumbent, stats, abort, domination);
            }

            solve(graph, bdd, incumbent, stats, abort, domination);
        }

        /* Clean up the timeout thread */
        if (timeout_thread.joinable()) {
            {
                unique_lock<mutex> guard(timeout_mutex);
                abort.store(true);
                timeout_cv.notify_all();
            }
            timeout_thread.join();
        }

        cout << "timeout = " << boolalpha << aborted << endl;
        cout << "solution =";
        for (auto i = incumbent.clique.find_first() ; i != dynamic_bitset<>::npos ; i = incumbent.clique.find_next(i))
            cout << " " << i + 1; // dimacs files are 1-indexed
        cout << endl;

        cout << "size = " << incumbent.value << endl;
        cout << "runtime = " << duration_cast<milliseconds>(steady_clock::now() - start_time).count() << "ms" << endl;
        cout << "nodes = " << stats.nodes_created << " + " << stats.nodes_reused << endl;
        cout << "branches = " << stats.bdds_created << " + " << stats.bdds_pruned << endl;

        // sanity check...
        for (auto i = incumbent.clique.find_first() ; i != dynamic_bitset<>::npos ; i = incumbent.clique.find_next(i))
            for (auto j = incumbent.clique.find_first() ; j != dynamic_bitset<>::npos ; j = incumbent.clique.find_next(j))
                if (i != j && ! graph.adjacent(i, j)) {
                    cout << "Oops! Horrific bug detected" << endl;
                    return EXIT_FAILURE;
                }

        return EXIT_SUCCESS;
    }
    catch (const po::error & e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Try " << argv[0] << " --help" << std::endl;
        return EXIT_FAILURE;
    }
    catch (const std::exception & e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}

