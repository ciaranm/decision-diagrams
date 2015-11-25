/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <utility>

#include "graph.hh"
#include "dimacs.hh"

using std::make_tuple;
using std::max;
using std::set_union;
using std::move;

using std::set;
using std::vector;

using std::cout;
using std::endl;

struct Node
{
    int score;
    set<int> p;
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
        if (d.p == node.p) {
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
         for (auto & v : n.p)
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
        if (n.p.count(branch_vertex)) {
            Node accept;
            accept.score = n.score + 1;
            for (auto & v : n.p)
                if (graph.adjacent(v, branch_vertex))
                    accept.p.insert(v);
            add_node(bdd.levels[level], move(accept), a_nodes_created, a_nodes_reused);
        }

        Node reject;
        reject.score = n.score;
        reject.p = n.p;
        reject.p.erase(branch_vertex);
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
                    return make_tuple(a.score, a.p.size()) > make_tuple(b.score, b.p.size());
                    });

            Node merged = { 0, {} };
            while (level_nodes.size() >= bdd.relaxed_max_width) {
                Node & n = level_nodes[level_nodes.size() - 1];
                merged.score = max(merged.score, n.score);
                set_union(merged.p.begin(), merged.p.end(), n.p.begin(), n.p.end(), inserter(merged.p, merged.p.end()));
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
                    return make_tuple(a.score, a.p.size()) > make_tuple(b.score, b.p.size());
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
                BDD relaxed_bdd{ n.p.size(), n.p.size(), n.p.size(), n.p.size(), vector<Level>(n.p.size() + 1) };
                relaxed_bdd.levels[0] = { Level{ { Node{ n.score, n.p } } } };
                int bound = solve_relaxed(graph, relaxed_bdd, relaxed_nodes_created, relaxed_nodes_reused);

                if (bound > incumbent) {
                    ++bdds_created;
                    BDD child_bdd{ n.p.size(), n.p.size(), n.p.size(), n.p.size(), vector<Level>(n.p.size() + 1) };
                    child_bdd.levels[0] = { Level{ { Node{ n.score, n.p } } } };

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

    set<int> all_vertices;
    for (int i = 0 ; i < graph.size() ; ++i)
        all_vertices.insert(i);
    bdd.levels[0] = { Level{ { Node{ 0, all_vertices } } } };

    int incumbent = 0;

    unsigned long long nodes_created = 1, nodes_reused = 0, relaxed_nodes_created = 0, relaxed_nodes_reused = 0,
                  restricted_nodes_created = 0, restricted_nodes_reused = 0, bdds_created = 1, bdds_pruned = 0;

    solve(graph, bdd, incumbent, nodes_created, nodes_reused, relaxed_nodes_created, relaxed_nodes_reused,
            restricted_nodes_created, restricted_nodes_reused, bdds_created, bdds_pruned);

    cout << incumbent << " " << nodes_created << " " << nodes_reused << " " << bdds_created << " " << bdds_pruned << endl;
}

