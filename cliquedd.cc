/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <utility>

#include "graph.hh"
#include "dimacs.hh"

struct Node
{
    int score;
    std::set<int> p;
};

struct Level
{
    std::vector<Node> nodes;
};

struct BDD
{
    int height;
    std::vector<Level> levels;
};

void add_node(Level & level, Node && node, unsigned long long & nodes_created, unsigned long long & nodes_reused)
{
    for (auto & d : level.nodes)
        if (d.p == node.p) {
            d.score = std::max(d.score, node.score);
            ++nodes_reused;
            return;
        }

    level.nodes.push_back(std::move(node));
    ++nodes_created;
}

const int max_width = 20;

int select_branch_vertex(const Graph & graph, const Level & level)
{
    std::vector<int> counts(graph.size());

    for (auto & n : level.nodes)
        for (auto & v : n.p)
            ++counts[v];

    int best = -1;

    for (unsigned i = 0 ; i < counts.size() ; ++i)
        if (counts[i] != 0 && (-1 == best || counts[i] < counts[best]))
            best = i;

    if (-1 == best)
        throw 0;

    return best;
}

void build_level(const Graph & graph, BDD & bdd, int level, unsigned long long & a_nodes_created, unsigned long long & a_nodes_reused)
{
    bdd.levels[level].nodes.reserve(bdd.levels[level - 1].nodes.size() * 2);

    int branch_vertex = select_branch_vertex(graph, bdd.levels[level - 1]);

    for (auto & n : bdd.levels[level - 1].nodes) {
        if (n.p.count(branch_vertex)) {
            Node accept;
            accept.score = n.score + 1;
            for (auto & v : n.p)
                if (graph.adjacent(v, branch_vertex))
                    accept.p.insert(v);
            add_node(bdd.levels[level], std::move(accept), a_nodes_created, a_nodes_reused);
        }

        Node reject;
        reject.score = n.score;
        reject.p = n.p;
        reject.p.erase(branch_vertex);
        add_node(bdd.levels[level], std::move(reject), a_nodes_created, a_nodes_reused);
    }
}

int solve_relaxed(const Graph & graph, BDD & bdd,
        unsigned long long & relaxed_nodes_created, unsigned long long & relaxed_nodes_reused)
{
    for (int level = 1 ; level <= bdd.height ; ++level) {
        build_level(graph, bdd, level, relaxed_nodes_created, relaxed_nodes_reused);

        auto & level_nodes = bdd.levels[level].nodes;
        if (level_nodes.size() > max_width) {
            std::sort(level_nodes.begin(), level_nodes.end(), [] (const Node & a, const Node & b) {
                    return std::make_tuple(a.score, a.p.size()) > std::make_tuple(b.score, b.p.size());
                    });

            Node merged = { 0, {} };
            while (level_nodes.size() >= max_width) {
                Node & n = level_nodes[level_nodes.size() - 1];
                merged.score = std::max(merged.score, n.score);
                std::set_union(merged.p.begin(), merged.p.end(), n.p.begin(), n.p.end(), std::inserter(merged.p, merged.p.end()));
                level_nodes.pop_back();
            }
            level_nodes.push_back(std::move(merged));
        }
    }

    int bound = 0;
    for (auto & n : bdd.levels[bdd.height].nodes)
        bound = std::max(bound, n.score);
    return bound;
}

void solve(const Graph & graph, BDD & bdd, int & incumbent, unsigned long long & nodes_created, unsigned long long & nodes_reused,
        unsigned long long & relaxed_nodes_created, unsigned long long & relaxed_nodes_reused,
        unsigned long long & bdds_created, unsigned long long & bdds_pruned)
{
    for (int level = 1 ; level <= bdd.height ; ++level) {
        build_level(graph, bdd, level, nodes_created, nodes_reused);

        if (bdd.levels[level].nodes.size() > max_width) {
            for (auto & n : bdd.levels[level].nodes) {
                BDD relaxed_bdd{ int(n.p.size()), std::vector<Level>(n.p.size() + 1) };
                relaxed_bdd.levels[0] = { Level{ { Node{ n.score, n.p } } } };
                int bound = solve_relaxed(graph, relaxed_bdd, relaxed_nodes_created, relaxed_nodes_reused);

                if (bound > incumbent) {
                    ++bdds_created;
                    BDD child_bdd{ int(n.p.size()), std::vector<Level>(n.p.size() + 1) };
                    child_bdd.levels[0] = { Level{ { Node{ n.score, n.p } } } };

                    solve(graph, child_bdd, incumbent, nodes_created, nodes_reused, relaxed_nodes_created,
                            relaxed_nodes_reused, bdds_created, bdds_pruned);
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
        incumbent = std::max(incumbent, n.score);
    if (incumbent > old_incumbent)
        std::cout << "-- " << incumbent << std::endl;
}

int main(int, char * argv[])
{
    Graph graph = read_dimacs(argv[1]);

    BDD bdd{ graph.size(), std::vector<Level>(graph.size() + 1) };

    std::set<int> all_vertices;
    for (int i = 0 ; i < graph.size() ; ++i)
        all_vertices.insert(i);
    bdd.levels[0] = { Level{ { Node{ 0, all_vertices } } } };

    int incumbent = 0;
    unsigned long long nodes_created = 1, nodes_reused = 0, relaxed_nodes_created = 0, relaxed_nodes_reused = 0,
                  bdds_created = 1, bdds_pruned = 0;
    solve(graph, bdd, incumbent, nodes_created, nodes_reused, relaxed_nodes_created, relaxed_nodes_reused, bdds_created, bdds_pruned);
    std::cout << incumbent << " " << nodes_created << " " << nodes_reused << " " << bdds_created << " " << bdds_pruned << std::endl;
}

