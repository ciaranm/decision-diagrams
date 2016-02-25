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

#include "graph.hh"
#include "dimacs.hh"

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
    void add_node(Node && node, Stats & stats, const bool dominate)
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
            else if (dominate && d.score <= node.score && d.state.is_subset_of(node.state)) {
                d = node;
                ++stats.nodes_reused;
                return;
            }
            else if (dominate && node.score <= d.score && node.state.is_subset_of(d.state)) {
                ++stats.nodes_reused;
                return;
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
                    break;
                }
            }
            else
                break;
        }
    }
};

void solve(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats, atomic<bool> & abort, const bool dominate);

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
void build_level(const Graph & graph, BDD & bdd, int level, int branch_vertex, Stats & stats, const bool dominate)
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
            bdd.levels[level].add_node(move(accept), stats, dominate);
        }

        // The new level of the BDD always gets a reject node.
        Node reject = n;
        reject.state.set(branch_vertex, false);
        bdd.levels[level].add_node(move(reject), stats, dominate);
    }
}

/**
 * Solve a relaxed form of the BDD, to get an upper bound on its actual best
 * value.
 */
unsigned solve_relaxed(const Graph & graph, BDD & bdd, Stats & stats, atomic<bool> & abort, const bool dominate)
{
    if (abort.load())
        return 0;

    // Build up each level...
    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex_in_most_states(graph, bdd.levels[level - 1]);
        build_level(graph, bdd, level, branch_vertex, stats, dominate);

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
    for (auto & n : bdd.levels[bdd.height].nodes)
        bound = max(bound, n.score);
    return bound;
}

/**
 * Solve a restricted form of the BDD, in the hopes of getting a good incumbent quickly.
 */
void solve_restricted(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats, atomic<bool> & abort, const bool dominate)
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

        build_level(graph, bdd, level, branch_vertex, stats, dominate);
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
 * Given a node in a BDD, solve it exactly by creating a new BDD.
 */
void solve_by_branching(const Graph & graph, const BDD & bdd, const Node & n, Incumbent & incumbent, Stats & stats, atomic<bool> & abort,
        const bool dominate)
{
    if (abort.load())
        return;

    BDD relaxed_bdd(n.state.count(), bdd.n_unassigned);

    // Start by getting a relaxed bound, to see if it's worth continuing.
    relaxed_bdd.levels[0] = { Level{ { n } } };
    unsigned bound = solve_relaxed(graph, relaxed_bdd, stats, abort, dominate);

    // Are we any good?
    if (bound > incumbent.value) {
        // Yes, now get an exact solution.
        ++stats.bdds_created;
        BDD child_bdd(n.state.count(), bdd.n_unassigned);
        child_bdd.levels[0] = { Level{ { n } } };

        solve(graph, child_bdd, incumbent, stats, abort, dominate);
    }
    else {
        ++stats.bdds_pruned;
    }
}

/**
 * Solve a BDD.
 */
void solve(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats, atomic<bool> & abort,
        const bool dominate)
{
    if (abort.load())
        return;

    // Start by building a restricted BDD, just to get a good incumbent quickly.
    BDD restricted_bdd = bdd;
    solve_restricted(graph, restricted_bdd, incumbent, stats, abort, dominate);

    // Build up each level...
    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex_in_most_states(graph, bdd.levels[level - 1]);
        build_level(graph, bdd, level, branch_vertex, stats, dominate);

        // Is our newly created level too wide?
        if (bdd.levels[level].nodes.size() > bdd.n_unassigned + 1) {
            // Yes, solve each of its nodes independently and in parallel using
            // new BDDs.
            for (auto & n : bdd.levels[level].nodes) {
                solve_by_branching(graph, bdd, n, incumbent, stats, abort, dominate);
            }
            return;
        }
    }

    // We built a complete BDD. Pick off its leaf nodes for new candidate
    // solutions.
    for (auto & n : bdd.levels[bdd.height].nodes)
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
    if (argc != 2 && (argc != 3 || string(argv[1]) != "-d")) {
        cout << "Usage: " << argv[0] << " [ -d ] file.clq" << endl;
        return EXIT_FAILURE;
    }

    int timeout = 14400;
    bool aborted = false;

    Graph graph = read_dimacs(argv[argc - 1]);

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

    BDD bdd = create_root_bdd(graph, stats);
    solve(graph, bdd, incumbent, stats, abort, argc == 3);

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

