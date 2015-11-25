/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "graph.hh"
#include <algorithm>

Graph::Graph(int size)
{
    if (0 != size)
        resize(size);
}

auto Graph::resize(int size) -> void
{
    _size = size;
    _adjacency.resize(_size);
    for (auto & r : _adjacency)
        r.resize(_size);
}

auto Graph::add_edge(int a, int b) -> void
{
    _adjacency[a][b] = 1;
    _adjacency[b][a] = 1;
}

auto Graph::adjacent(int a, int b) const -> bool
{
    return _adjacency[a][b];
}

auto Graph::size() const -> int
{
    return _size;
}

auto Graph::degree(int a) const -> int
{
    return _adjacency[a].count();
}

auto Graph::neighbourhood(int a) const -> const boost::dynamic_bitset<> &
{
    return _adjacency[a];
}

