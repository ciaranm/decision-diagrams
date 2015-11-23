/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <graph_file_error.hh>

GraphFileError::GraphFileError(const std::string & filename, const std::string & message) throw () :
    _what("Error reading graph file '" + filename + "': " + message)
{
}

auto GraphFileError::what() const throw () -> const char *
{
    return _what.c_str();
}

