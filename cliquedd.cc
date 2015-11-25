/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <utility>

#include <boost/dynamic_bitset.hpp>

#include "graph.hh"
#include "dimacs.hh"

using std::make_tuple;
using std::max;
using std::set_union;
using std::move;
using std::vector;

using std::cout;
using std::endl;

using boost::dynamic_bitset;

struct Node
{
    int score;
    dynamic_bitset<> state;
};

struct Level
{
    vector<Node> nodes;
};

struct BDD
{
    unsigned full_max_width;
    unsigned relaxed_max_width;
    unsigned restricted_max_width;
    unsigned height;

    vector<Level> levels;
};

void add_node(Level & level, Node && node, unsigned long long & nodes_created, unsigned long long & nodes_reused)
{
    for (auto & d : level.nodes)
        if (d.state == node.state) {
            d.score = max(d.score, node.score);
            ++nodes_reused;
            return;
        }

    level.nodes.push_back(move(node));
    ++nodes_created;
}

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

void build_level(const Graph & graph, BDD & bdd, int level, int branch_vertex,
        unsigned long long & a_nodes_created, unsigned long long & a_nodes_reused)
{
    bdd.levels[level].nodes.reserve(bdd.levels[level - 1].nodes.size() * 2);

    for (auto & n : bdd.levels[level - 1].nodes) {
        if (n.state[branch_vertex]) {
            Node accept = { n.score + 1, n.state };
            accept.score = n.score + 1;
            accept.state &= graph.neighbourhood(branch_vertex);
            add_node(bdd.levels[level], move(accept), a_nodes_created, a_nodes_reused);
        }

        Node reject;
        reject.score = n.score;
        reject.state = n.state;
        reject.state.set(branch_vertex, false);
        add_node(bdd.levels[level], move(reject), a_nodes_created, a_nodes_reused);
    }
}

int solve_relaxed(const Graph & graph, BDD & bdd,
        unsigned long long & relaxed_nodes_created, unsigned long long & relaxed_nodes_reused)
{
    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        int branch_vertex = select_branch_vertex(graph, bdd.levels[level - 1]);
        build_level(graph, bdd, level, branch_vertex, relaxed_nodes_created, relaxed_nodes_reused);

        auto & level_nodes = bdd.levels[level].nodes;
        if (level_nodes.size() > bdd.relaxed_max_width) {
            sort(level_nodes.begin(), level_nodes.end(), [] (const Node & a, const Node & b) {
                    return make_tuple(a.score, a.state.count()) > make_tuple(b.score, b.state.count());
                    });

            Node merged = { 0, dynamic_bitset<>(graph.size(), 0) };
            while (level_nodes.size() >= bdd.relaxed_max_width) {
                Node & n = level_nodes[level_nodes.size() - 1];
                merged.score = max(merged.score, n.score);
                merged.state |= n.state;
                level_nodes.pop_back();
            }
            level_nodes.push_back(move(merged));
        }
    }

    int bound = 0;
    for (auto & n : bdd.levels[bdd.height].nodes)
        bound = max(bound, n.score);
    return bound;
}

void solve_restricted(const Graph & graph, BDD & bdd, int & incumbent, unsigned long long & restricted_nodes_created,
        unsigned long long & restricted_nodes_reused)
{
    int last_level = 0;

    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        int branch_vertex = select_branch_vertex(graph, bdd.levels[level - 1]);
        if (-1 == branch_vertex)
            break;

        build_level(graph, bdd, level, branch_vertex, restricted_nodes_created, restricted_nodes_reused);
        last_level = level;

        if (bdd.levels[level].nodes.size() > bdd.restricted_max_width) {
            sort(bdd.levels[level].nodes.begin(), bdd.levels[level].nodes.end(), [] (const Node & a, const Node & b) {
                    return make_tuple(a.score, a.state.count()) > make_tuple(b.score, b.state.count());
                    });

            bdd.levels[level].nodes.resize(bdd.restricted_max_width);
        }
    }

    int old_incumbent = incumbent;
    for (auto & n : bdd.levels[last_level].nodes)
        incumbent = max(incumbent, n.score);
    if (incumbent > old_incumbent)
        cout << "~~ " << incumbent << endl;
}

void solve(const Graph & graph, BDD & bdd, int & incumbent, unsigned long long & nodes_created, unsigned long long & nodes_reused,
        unsigned long long & relaxed_nodes_created, unsigned long long & relaxed_nodes_reused,
        unsigned long long & restricted_nodes_created, unsigned long long & restricted_nodes_reused,
        unsigned long long & bdds_created, unsigned long long & bdds_pruned)
{
    BDD restricted_bdd = bdd;
    solve_restricted(graph, restricted_bdd, incumbent, restricted_nodes_created, restricted_nodes_reused);

    for (unsigned level = 1 ; level <= bdd.height ; ++level) {
        int branch_vertex = select_branch_vertex(graph, bdd.levels[level - 1]);
        build_level(graph, bdd, level, branch_vertex, nodes_created, nodes_reused);

        if (bdd.levels[level].nodes.size() > bdd.full_max_width) {
            for (auto & n : bdd.levels[level].nodes) {
                BDD relaxed_bdd{ n.state.count(), n.state.count(), n.state.count(), n.state.count(), vector<Level>(n.state.count() + 1) };
                relaxed_bdd.levels[0] = { Level{ { Node{ n.score, n.state } } } };
                int bound = solve_relaxed(graph, relaxed_bdd, relaxed_nodes_created, relaxed_nodes_reused);

                if (bound > incumbent) {
                    ++bdds_created;
                    BDD child_bdd{ n.state.count(), n.state.count(), n.state.count(), n.state.count(), vector<Level>(n.state.count() + 1) };
                    child_bdd.levels[0] = { Level{ { Node{ n.score, n.state } } } };

                    solve(graph, child_bdd, incumbent, nodes_created, nodes_reused, relaxed_nodes_created,
                            relaxed_nodes_reused, restricted_nodes_created, restricted_nodes_reused,
                            bdds_created, bdds_pruned);
                }
                else {
                    ++bdds_pruned;
                }
            }
            return;
        }
    }

    int old_incumbent = incumbent;
    for (auto & n : bdd.levels[bdd.height].nodes)
        incumbent = max(incumbent, n.score);
    if (incumbent > old_incumbent)
        cout << "-- " << incumbent << endl;
}

int main(int, char * argv[])
{
    Graph graph = read_dimacs(argv[1]);

    BDD bdd{ unsigned(graph.size()), unsigned(graph.size()), unsigned(graph.size()), unsigned(graph.size()), vector<Level>(graph.size() + 1) };

    dynamic_bitset<> all_vertices(graph.size(), 0);
    for (int i = 0 ; i < graph.size() ; ++i)
        all_vertices.set(i);
    bdd.levels[0] = { Level{ { Node{ 0, all_vertices } } } };

    int incumbent = 0;

    unsigned long long nodes_created = 1, nodes_reused = 0, relaxed_nodes_created = 0, relaxed_nodes_reused = 0,
                  restricted_nodes_created = 0, restricted_nodes_reused = 0, bdds_created = 1, bdds_pruned = 0;

    solve(graph, bdd, incumbent, nodes_created, nodes_reused, relaxed_nodes_created, relaxed_nodes_reused,
            restricted_nodes_created, restricted_nodes_reused, bdds_created, bdds_pruned);

    cout << incumbent << " " << nodes_created << " " << nodes_reused << " " << bdds_created << " " << bdds_pruned << endl;
}

