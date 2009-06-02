#! /usr/bin/env python
# graph-wars -- a general graph contagion dynamics thingy
#
# Copyright (C) 2006  Alexandre Hannud Abdo <abdo@member.fsf.org>
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
dl_import("import libgraph_tool_contagion")

from .. core import _degree, _prop, GraphError
from numpy import *
import random, sys

__all__ = ["strategic_seed", "generalized_contagion"]

def strategic_seed(g, strategy, number_of_seeds, property, rng_seed=None):
    """Seeds the state vertex propety of a graph in preparation to run a
    generalized contagion model"""
    if rng_seed is None:
        rng_seed = random.randint(0, sys.maxint)
    libgraph_tool_contagion.\
        strategic_seed(g._Graph__graph, strategy, number_of_seeds,
                       _prop("v", g, property), rng_seed)

def generalized_contagion(g, parameters, properties, iterations,
                          rng_seed=None, verbose=False):
    """Runs the generalized contagion model for a given number of iterations"""
    parameter_names = ("infection_rate", "recovery_rate",
                       "resusceptibility_rate",
                       "infection_dose_distribution", "timescale")
    property_names = ("sir_state","dose_threshold","doses", "availability")
    if any([name not in parameters for name in parameter_names]):
        raise GraphError("Missing parameters.")
    if any([name not in properties for name in property_names]):
        raise GraphError("Missing properties.")
    for name in property_names:
        properties[name] = _prop("v", g, properties[name])
    if rng_seed is None:
        rng_seed = random.randint(0, sys.maxint)
    return libgraph_tool_contagion.generalized_contagion\
        (g._Graph__graph, parameters, properties, iterations, rng_seed, verbose)
