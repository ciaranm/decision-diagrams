/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <set>
#include <list>
#include <iostream>

#include "graph.hh"
#include "dimacs.hh"

struct Node
{
    int score;
    std::set<int> p;
};

struct Level
{
    std::list<Node> nodes;
};

struct BDD
{
    int height;
    std::vector<Level> levels;
};

void add_node(Level & level, Node & node)
{
    for (auto & d : level.nodes)
        if (d.p == node.p) {
            d.score = std::max(d.score, node.score);
            return;
        }

    level.nodes.push_back(node);
}

const int max_width = 10;

void solve(const Graph & graph, BDD & bdd, int & incumbent)
{
    for (int level = 1 ; level <= bdd.height ; ++level) {
        int branch_vertex = *std::next(bdd.levels[0].nodes.begin()->p.begin(), level - 1);

        for (auto & n : bdd.levels[level - 1].nodes) {
            if (n.p.count(branch_vertex)) {
                Node accept;
                accept.score = n.score + 1;
                for (auto & v : n.p)
                    if (graph.adjacent(v, branch_vertex))
                        accept.p.insert(v);
                add_node(bdd.levels[level], accept);
            }

            Node reject;
            reject.score = n.score;
            reject.p = n.p;
            reject.p.erase(branch_vertex);
            add_node(bdd.levels[level], reject);
        }

        if (bdd.levels[level].nodes.size() > max_width) {
            for (auto & n : bdd.levels[level].nodes) {
                BDD child_bdd{ int(n.p.size()), std::vector<Level>(n.p.size() + 1) };
                child_bdd.levels[0] = { Level{ { Node{ n.score, n.p } } } };

                solve(graph, child_bdd, incumbent);
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
    solve(graph, bdd, incumbent);
    std::cout << incumbent << std::endl;
}

