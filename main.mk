BUILD_DIR := intermediate
TARGET_DIR := ./

boost_ldlibs := -lboost_regex -lboost_thread -lboost_system -lboost_program_options

override CXXFLAGS += -O3 -march=native -std=c++14 -I./ -W -Wall -g -ggdb3 -pthread -fcilkplus
override LDFLAGS += -pthread $(boost_ldlibs) -fcilkplus

TARGET := cliquedd

SOURCES := \
    graph.cc \
    dimacs.cc \
    graph_file_error.cc \
    cliquedd.cc

TGT_LDLIBS := $(boost_ldlibs)


