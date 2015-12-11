/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <utility>
#include <atomic>

#include <boost/dynamic_bitset.hpp>

#include <cilk/cilk.h>

#include "graph.hh"
#include "dimacs.hh"

using std::make_tuple;
using std::max;
using std::set_union;
using std::move;
using std::vector;
using std::atomic;

using std::cout;
using std::endl;

using boost::dynamic_bitset;

struct Node
{
    unsigned score;           // number of accepted vertices so far
    dynamic_bitset<> state;   // vertices which remain to be selectable
};

struct Level
{
    vector<Node> nodes;       // children of a particular BDD level
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

struct Stats
{
    atomic<unsigned long long> nodes_created{ 0 };
    atomic<unsigned long long> nodes_reused{ 0 };
    atomic<unsigned long long> bdds_created{ 0 };
    atomic<unsigned long long> bdds_pruned{ 0 };
};

struct Incumbent
{
    atomic<unsigned> value{ 0 };

    void update(unsigned new_value)
    {
        while (true) {
            unsigned current_value = value;
            if (new_value > current_value) {
                if (value.compare_exchange_strong(current_value, new_value)) {
                    cout << "found " << new_value << endl;
                    break;
                }
            }
            else
                break;
        }
    }
};

void solve(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats);

/**
 * Add node to level, updating stats, merging it if it is a duplicate.
 */
void add_node(Level & level, Node && node, Stats & stats)
{
    for (auto & d : level.nodes)
        if (d.state == node.state) {
            d.score = max(d.score, node.score);
            ++stats.nodes_reused;
            return;
        }

    level.nodes.push_back(move(node));
    ++stats.nodes_created;
}

/**
 * Select a vertex to branch on, from level. May return -1, if there are no
 * vertices to branch on in the restricted case.
 */
int select_branch_vertex(const Graph & graph, const Level & level)
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
void build_level(const Graph & graph, BDD & bdd, int level, int branch_vertex, Stats & stats)
{
    bdd.levels[level].nodes.reserve(bdd.levels[level - 1].nodes.size() * 2);

    // For each node on the previous level...
    for (auto & n : bdd.levels[level - 1].nodes) {
        // Does this node's state allow us to accept the branch vertex?
        if (n.state[branch_vertex]) {
            // Yes, the new level of the BDD gets an accept node.
            Node accept = { n.score + 1, n.state };
            accept.score = n.score + 1;
            accept.state &= graph.neighbourhood(branch_vertex);
            add_node(bdd.levels[level], move(accept), stats);
        }

        // The new level of the BDD always gets a reject node.
        Node reject;
        reject.score = n.score;
        reject.state = n.state;
        reject.state.set(branch_vertex, false);
        add_node(bdd.levels[level], move(reject), stats);
    }
}

/**
 * Solve a relaxed form of the BDD, to get an upper bound on its actual best
 * value.
 */
unsigned solve_relaxed(const Graph & graph, BDD & bdd, Stats & stats)
{
    // Build up each level...
    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex(graph, bdd.levels[level - 1]);
        build_level(graph, bdd, level, branch_vertex, stats);

        // Is our newly created level too wide?
        auto & level_nodes = bdd.levels[level].nodes;
        if (level_nodes.size() > bdd.n_unassigned + 1) {
            // Yup, start merging together states.
            sort(level_nodes.begin(), level_nodes.end(), [] (const Node & a, const Node & b) {
                    return make_tuple(a.score, a.state.count()) > make_tuple(b.score, b.state.count());
                    });

            Node merged = { 0, dynamic_bitset<>(graph.size(), 0) };
            while (level_nodes.size() >= bdd.n_unassigned + 1) {
                Node & n = level_nodes[level_nodes.size() - 1];
                merged.score = max(merged.score, n.score);
                merged.state |= n.state;
                level_nodes.pop_back();
            }
            level_nodes.push_back(move(merged));
        }
    }

    // The bound is found by looking at the score of every leaf node and
    // picking the best.
    unsigned bound = 0;
    for (auto & n : bdd.levels[bdd.height].nodes)
        bound = max(bound, n.score);
    return bound;
}

/**
 * Solve a restricted form of the BDD, in the hopes of getting a good incumbent quickly.
 */
void solve_restricted(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats)
{
    int last_level = 0;

    // Build up each level...
    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex(graph, bdd.levels[level - 1]);
        if (-1 == branch_vertex)
            break;

        build_level(graph, bdd, level, branch_vertex, stats);
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
        incumbent.update(n.score);
}

/**
 * Given a node in a BDD, solve it exactly by creating a new BDD.
 */
void solve_by_branching(const Graph & graph, const BDD & bdd, const Node & n, Incumbent & incumbent, Stats & stats)
{
    BDD relaxed_bdd(n.state.count(), bdd.n_unassigned);

    // Start by getting a relaxed bound, to see if it's worth continuing.
    relaxed_bdd.levels[0] = { Level{ { Node{ n.score, n.state } } } };
    unsigned bound = solve_relaxed(graph, relaxed_bdd, stats);

    // Are we any good?
    if (bound > incumbent.value) {
        // Yes, now get an exact solution.
        ++stats.bdds_created;
        BDD child_bdd(n.state.count(), bdd.n_unassigned);
        child_bdd.levels[0] = { Level{ { Node{ n.score, n.state } } } };

        solve(graph, child_bdd, incumbent, stats);
    }
    else {
        ++stats.bdds_pruned;
    }
}

/**
 * Solve a BDD.
 */
void solve(const Graph & graph, BDD & bdd, Incumbent & incumbent, Stats & stats)
{
    // Start by building a restricted BDD, just to get a good incumbent quickly.
    BDD restricted_bdd = bdd;
    solve_restricted(graph, restricted_bdd, incumbent, stats);

    // Build up each level...
    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        --bdd.n_unassigned;

        // Pick a branch vertex, and build the next level.
        int branch_vertex = select_branch_vertex(graph, bdd.levels[level - 1]);
        build_level(graph, bdd, level, branch_vertex, stats);

        // Is our newly created level too wide?
        if (bdd.levels[level].nodes.size() > bdd.n_unassigned + 1) {
            // Yes, solve each of its nodes independently and in parallel using
            // new BDDs.
            for (auto & n : bdd.levels[level].nodes) {
                cilk_spawn solve_by_branching(graph, bdd, n, incumbent, stats);
            }
            cilk_sync;
            return;
        }
    }

    // We built a complete BDD. Pick off its leaf nodes for new candidate
    // solutions.
    for (auto & n : bdd.levels[bdd.height].nodes)
        incumbent.update(n.score);
}

BDD create_root_bdd(const Graph & graph, Stats & stats)
{
    BDD bdd(graph.size(), graph.size());

    // Initially, we can pick all vertices.
    dynamic_bitset<> all_vertices(graph.size());
    all_vertices.set();

    // Create BDD with root node.
    bdd.levels[0] = { Level{ { Node{ 0, all_vertices } } } };

    ++stats.bdds_created;

    return bdd;
}

int main(int, char * argv[])
{
    Graph graph = read_dimacs(argv[1]);

    Stats stats;
    Incumbent incumbent;

    BDD bdd = create_root_bdd(graph, stats);
    solve(graph, bdd, incumbent, stats);

    cout << "size=" << incumbent.value << " nodes=" << stats.nodes_created << " nodes_reused=" << stats.nodes_reused
        << " branches=" << stats.bdds_created << " pruned_branches=" << stats.bdds_pruned << endl;
}

