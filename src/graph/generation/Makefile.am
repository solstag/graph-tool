## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CXXFLAGS = $(CXXFLAGS)

libgraph_tool_generationdir = $(pythondir)/graph_tool/generation

libgraph_tool_generation_LTLIBRARIES = libgraph_tool_generation.la

libgraph_tool_generation_la_includedir = $(pythondir)/graph_tool/include

libgraph_tool_generation_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_generation_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_generation_la_SOURCES = \
    graph_generation.cc

libgraph_tool_generation_la_include_HEADERS = \
    graph_generation.hh
