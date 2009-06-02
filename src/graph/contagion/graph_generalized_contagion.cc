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
#include <boost/vector_property_map.hpp>

#include "graph.hh"
#include "graph_properties.hh"
#include "graph_filtering.hh"

#include "graph_generalized_contagion.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

struct doses_from_string
{
    template <class Graph, class DS, class D, class AD>
    void operator()(const Graph& g, DS doses_as_string, D doses,
                    AD accumulated_dose) const
    {
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            string s = doses_as_string[*vi];
            istringstream ist(s);
            double dose;
            while(ist >> dose)
                doses[*vi].push_back(dose);
            accumulated_dose[*vi] =
                accumulate(doses[*vi].begin(), doses[*vi].end(), 0.0);
        }
    }
};

struct doses_to_string
{
    template <class Graph, class DS, class D>
    void operator()(const Graph& g, DS doses_as_string, D doses) const
    {
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            {
                ostringstream ost;
                deque<double>::const_iterator dqi = doses[*vi].begin();
                ost << *dqi;
                for(++dqi; dqi != doses[*vi].end(); ++dqi) ost << " " << *dqi;
                doses_as_string[*vi] = ost.str();
            }
    }
};

python::list
generalized_contagion(GraphInterface& gi, python::object parameters,
                            python::object properties, size_t iterations,
                            size_t rng_seed, bool verbose)
{
    int timescale =
        python::extract<int>(parameters["timescale"])();
    double infection_rate =
        python::extract<double>(parameters["infection_rate"])();
    double recovery_rate =
        python::extract<double>(parameters["recovery_rate"])();
    double resusceptibility_rate =
        python::extract<double>(parameters["resusceptibility_rate"])();
    string infection_dose_distribution =
        python::extract<string>(parameters["infection_dose_distribution"])();

    typedef GraphInterface::vertex_index_map_t g_vim_t;
    checked_vector_property_map<string, g_vim_t> doses_as_string =
        any_cast<checked_vector_property_map<string, g_vim_t> >
        (python::extract<any>(properties["doses"]));
    checked_vector_property_map<int, g_vim_t> sir_state =
        any_cast<checked_vector_property_map<int, g_vim_t> >
        (python::extract<any>(properties["sir_state"]));
    checked_vector_property_map<double, g_vim_t> dose_threshold =
        any_cast<checked_vector_property_map<double, g_vim_t> >
        (python::extract<any>(properties["dose_threshold"]));
    checked_vector_property_map<int, g_vim_t> availability =
        any_cast<checked_vector_property_map<int, g_vim_t> >
        (python::extract<any>(properties["availability"]));

    checked_vector_property_map<deque<double>, g_vim_t> doses;
    checked_vector_property_map<double, g_vim_t> accumulated_dose;
    run_action<>()(gi, lambda::bind<void>(doses_from_string(), lambda::_1,
                                          doses_as_string, doses,
                                          accumulated_dose))();

    vector<vector<size_t> > state_members(3);
    rng_t rng(static_cast<rng_t::result_type>(rng_seed));
	run_action<>()
    (gi, lambda::bind<void>(
        advance_contagion
        <checked_vector_property_map<double, g_vim_t>,
         checked_vector_property_map<deque<double>, g_vim_t>,
         checked_vector_property_map<double, g_vim_t>,
         checked_vector_property_map<int, g_vim_t>,
         checked_vector_property_map<int, g_vim_t> >
        (&state_members, &dose_threshold, &doses, &accumulated_dose,
         &sir_state, &availability, rng),
     lambda::_1, infection_rate, recovery_rate, resusceptibility_rate,
     infection_dose_distribution, iterations, timescale))();

    run_action<>()(gi, lambda::bind<void>(doses_to_string(), lambda::_1,
                                          doses_as_string, doses))();
    
    python::list l;
    for(vector<vector<size_t> >::iterator i1 = state_members.begin();
        i1 != state_members.end(); ++i1)
    {
        python::list li;
        for(vector<size_t>::iterator i2 = i1->begin(); i2 != i1->end(); ++i2)
            li.append(*i2);
        l.append(li);
    }
    return l;
}

using namespace boost::python;

void
strategic_seed(GraphInterface& gi, string seeding_strategy, size_t
               number_of_seeds, boost::any property, size_t rng_seed);

BOOST_PYTHON_MODULE(libgraph_tool_contagion)
{
    def("generalized_contagion", &generalized_contagion);
    def("strategic_seed", &strategic_seed);
}
