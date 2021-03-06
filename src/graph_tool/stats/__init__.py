#! /usr/bin/env python
# graph_tool.py -- a general graph manipulation python module
#
# Copyright (C) 2007 Tiago de Paula Peixoto <tiago@forked.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from .. dl_import import dl_import
dl_import("import libgraph_tool_stats")

from .. core import _degree, _prop
from numpy import *

__all__ = ["vertex_hist", "edge_hist", "vertex_average", "edge_average",
           "label_components", "label_parallel_edges", "remove_parallel_edges",
           "label_self_loops", "remove_self_loops", "remove_labeled_edges"]

def vertex_hist(g, deg, bins=[1], float_count=True):
    ret = libgraph_tool_stats.\
          get_vertex_histogram(g._Graph__graph, _degree(g, deg), bins)
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]

def edge_hist(g, eprop, bins=[1], float_count=True):
    ret = libgraph_tool_stats.\
          get_edge_histogram(g._Graph__graph, _prop("e", g, eprop), bins)
    return [array(ret[0], dtype="float64") if float_count else ret[0], ret[1]]

def vertex_average(g, deg):
    ret = libgraph_tool_stats.\
          get_vertex_average(g._Graph__graph, _degree(g, deg))
    return ret

def edge_average(g, eprop):
    ret = libgraph_tool_stats.\
          get_edge_average(g._Graph__graph, _prop("e", g, eprop))
    return ret

def label_components(g, vprop=None):
    if vprop == None:
        vprop = g.new_vertex_property("int32_t")
    libgraph_tool_stats.\
          label_components(g._Graph__graph, _prop("v", g, vprop))
    return vprop

def remove_labeled_edges(g, label):
    g.stash_filter(all=False, directed=True, reversed=True)
    libgraph_tool_stats.\
          remove_labeled_edges(g._Graph__graph, _prop("e", g, label))
    g.pop_filter(all=False, directed=True, reversed=True)

def label_parallel_edges(g, eprop=None):
    if eprop == None:
        eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_parallel_edges(g._Graph__graph, _prop("e", g, eprop))
    return eprop

def remove_parallel_edges(g):
    eprop = label_parallel_edges(g)
    remove_labeled_edges(g, eprop)

def label_self_loops(g, eprop=None):
    if eprop == None:
        eprop = g.new_edge_property("int32_t")
    libgraph_tool_stats.\
          label_self_loops(g._Graph__graph, _prop("e", g, eprop))
    return eprop

def remove_self_loops(g):
    eprop = label_self_loops(g)
    remove_labeled_edges(g, eprop)

