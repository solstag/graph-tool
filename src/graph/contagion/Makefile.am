## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_contagiondir = $(pythondir)/graph_tool/contagion

libgraph_tool_contagion_LTLIBRARIES = libgraph_tool_contagion.la

libgraph_tool_contagion_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_contagion_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_contagion_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_contagion_la_SOURCES = \
    graph_generalized_contagion.cc\
    graph_strategic_seeding.cc

libgraph_tool_contagion_la_include_HEADERS = \
    graph_generalized_contagion.hh\
    graph_strategic_seeding.hh
