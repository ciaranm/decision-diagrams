/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GUARD_GRAPH_HH
#define GUARD_GRAPH_HH 1

#include <vector>
#include <string>
#include <type_traits>
#include <cstdint>

#include <boost/dynamic_bitset.hpp>

/**
 * A graph, with an adjaceny matrix representation. We only provide the
 * operations we actually need.
 *
 * Indices start at 0.
 */
class Graph
{
    private:
        using AdjacencyMatrix = std::vector<boost::dynamic_bitset<> >;
        using Weights = std::vector<unsigned long long>;

        int _size = 0;
        AdjacencyMatrix _adjacency;
        Weights _weights;

    public:
        /**
         * \param initial_size can be 0, if resize() is called afterwards.
         */
        Graph(int initial_size);

        Graph(const Graph &) = default;

        explicit Graph() = default;

        /**
         * Number of vertices.
         */
        auto size() const -> int;

        /**
         * Change our size. Must be called before adding an edge, and must
         * not be called afterwards.
         */
        auto resize(int size) -> void;

        /**
         * Add an edge from a to b (and from b to a).
         */
        auto add_edge(int a, int b) -> void;

        /**
         * Are vertices a and b adjacent?
         */
        auto adjacent(int a, int b) const -> bool;

        /**
         * What is the degree of a given vertex?
         */
        auto degree(int a) const -> int;

        /**
         * What is the neighbourhood of a given vertex?
         */
        auto neighbourhood(int v) const -> const boost::dynamic_bitset<> &;

        auto set_weight(int v, unsigned long long w) -> void;

        auto weight(int v) const -> unsigned long long;
};

#endif
