// graph-wars -- a general graph contagion dynamics thingy
//
// Copyright (C) 2006  Alexandre Hannud Abdo <abdo@member.fsf.org>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <boost/python.hpp>

#include "graph.hh"
#include "graph_properties.hh"
#include "graph_filtering.hh"

#include "graph_strategic_seeding.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void
strategic_seed(GraphInterface& gi, string seeding_strategy, size_t
               number_of_seeds, boost::any property, size_t rng_seed)
{
    rng_t rng(static_cast<rng_t::result_type>(rng_seed));

    if( seeding_strategy == "random" )
        run_action<>()(gi, lambda::bind(random_seed(), lambda::_1, lambda::_2,
                                number_of_seeds, ref(rng)),
                     writable_vertex_scalar_properties())(property);
    else if( seeding_strategy == "top_bottom" )
    {
        bool reversed = gi.GetReversed();
        gi.SetReversed(!reversed);
        run_action<>()(gi, lambda::bind(top_bottom_seed(), lambda::_1,
                                        lambda::_2, number_of_seeds, ref(rng)),
                     writable_vertex_scalar_properties())(property);
        gi.SetReversed(reversed);
    }
    else if( seeding_strategy == "neighbor" )
        run_action<>()(gi, lambda::bind(neighbor_seed(), lambda::_1, lambda::_2,
                                        number_of_seeds, ref(rng)),
                     writable_vertex_scalar_properties())(property);
    else if( seeding_strategy == "walker" )
        run_action<>()(gi, lambda::bind(walker_seed(), lambda::_1, lambda::_2,
                                number_of_seeds, ref(rng)),
                     writable_vertex_scalar_properties())(property);
    else if( seeding_strategy == "starwalker" )
        run_action<>()(gi, lambda::bind(starwalker_seed(), lambda::_1,
                                        lambda::_2, number_of_seeds, ref(rng)),
                     writable_vertex_scalar_properties())(property);
    else if( seeding_strategy == "bfs" )
        run_action<>()(gi, lambda::bind(bfs_seed(), lambda::_1, lambda::_2,
                                        number_of_seeds, ref(rng),
                                        gi.GetVertexIndex()),
                     writable_vertex_scalar_properties())(property);
    else throw GraphException("error running dynamics: "
                              "requested invalid seeding strategy '"+
                              seeding_strategy + "'");
}
