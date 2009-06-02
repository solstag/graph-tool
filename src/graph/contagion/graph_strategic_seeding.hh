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

#ifndef GRAPH_STRATEGIC_SEEDING_HH
#define GRAPH_STRATEGIC_SEEDING_HH

#include <algorithm>
#include <ext/algorithm>
#include <numeric>
#include <vector>
#include <string>
#include <sstream>
#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <tr1/random>
typedef std::tr1::mt19937 rng_t;

namespace graph_tool
{
using namespace std;
using namespace boost;
using namespace boost::lambda;

struct random_seed
{
    typedef void result_type;
    template<class Graph, class SIRState, class NumberOfSeeds>
    void operator()(const Graph& g, SIRState& sir_state,
                    NumberOfSeeds number_of_seeds, rng_t& rng) const
    {
        if(number_of_seeds > HardNumVertices()(g))
            throw GraphException("error running dynamics: "
                                 "number of seeds too high:'" +
                                 lexical_cast<string>(number_of_seeds) + "'");
        enum states { S, I, R, NUMBER_OF_STATES };

        vector<typename graph_traits<Graph>::vertex_descriptor> descriptors;
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            descriptors.push_back(*vi);

        typedef tr1::uniform_int<size_t> idist_t;
        tr1::variate_generator<rng_t, idist_t> random(rng, idist_t());

        vector<typename graph_traits<Graph>::vertex_descriptor>
            sample(number_of_seeds);
        __gnu_cxx::random_sample(descriptors.begin(), descriptors.end(),
                                 sample.begin(), sample.end(), random);
        for(size_t i = 0; i < sample.size(); ++i) sir_state[sample[i]] = I ;
    }
};

struct top_bottom_seed
{
    typedef void result_type;
    template <class Graph, class SIRState, class NumberOfSeeds>
    void operator()(const Graph& g, SIRState sir_state,
                    NumberOfSeeds number_of_seeds, rng_t& rng) const
    {
        if(number_of_seeds > HardNumVertices()(g))
            throw GraphException("error running dynamics: "
                                 "number of seeds too high:'"+
                                 lexical_cast<string>(number_of_seeds) + "'");
        enum states { S, I, R, NUMBER_OF_STATES };
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;

        set<size_t> out_degrees;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            out_degrees.insert( out_degree(*vi, g) );
        }

        set<size_t>::reverse_iterator max_degree = out_degrees.rbegin();
        vector<typename graph_traits<Graph>::vertex_descriptor> descriptors;
        while(number_of_seeds)
        {
            for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
                if ( out_degree(*vi, g) == *max_degree )
                    descriptors.push_back(*vi);

            if(number_of_seeds >= descriptors.size())
            {
                for(size_t i=0; i < descriptors.size(); ++i)
                    sir_state[descriptors[i]] = I;
                number_of_seeds -= descriptors.size();
            }
            else
            {
                typedef tr1::uniform_int<size_t> idist_t;
                tr1::variate_generator<rng_t, idist_t> random(rng, idist_t());
                vector<typename graph_traits<Graph>::vertex_descriptor>
                    sample(number_of_seeds);
                __gnu_cxx::random_sample(descriptors.begin(), descriptors.end(),
                                         sample.begin(), sample.end(), random);
                for(size_t i = 0; i < sample.size(); ++i)
                    sir_state[sample[i]] = I;
                number_of_seeds -= sample.size();
            }

            descriptors.clear();
            ++max_degree;
        }
    }
};

struct neighbor_seed
{
    typedef void result_type;
    template<class Graph, class SIRState, class NumberOfSeeds>
    void operator()(const Graph& g, SIRState& sir_state,
                    NumberOfSeeds number_of_seeds, rng_t& rng) const
    {
        if(number_of_seeds > HardNumVertices()(g))
            throw GraphException("error running dynamics: "
                                 "number of seeds too high:'"+
                                 lexical_cast<string>(number_of_seeds) + "'");
        enum states { S, I, R, NUMBER_OF_STATES };

        vector<typename graph_traits<Graph>::vertex_descriptor> descriptors;
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            descriptors.push_back(*vi);

        typedef tr1::uniform_int<size_t> idist_t;
        tr1::variate_generator<rng_t, idist_t> random(rng, idist_t());
        while(number_of_seeds)
        {
            vector<typename graph_traits<Graph>::vertex_descriptor>
                sample(1), neighbor(1);
            __gnu_cxx::random_sample(descriptors.begin(), descriptors.end(),
                                     sample.begin(), sample.end(), random);
            if( out_degree(sample[0], g) == 0 ) continue;
            typename graph_traits<Graph>::adjacency_iterator ni, ni_end;
            tie(ni,ni_end) = adjacent_vertices(sample[0], g);
            __gnu_cxx::random_sample(ni, ni_end,
                                     neighbor.begin(), neighbor.end(), random);
            if(sir_state[neighbor[0]] != I )
            {
                sir_state[neighbor[0]] = I;
                --number_of_seeds;
            }
        }
    }
};

struct walker_seed
{
    typedef void result_type;
    template<class Graph, class SIRState, class NumberOfSeeds>
    void operator()(const Graph& g, SIRState& sir_state,
                    NumberOfSeeds number_of_seeds, rng_t& rng) const
    {
        if(number_of_seeds > HardNumVertices()(g))
            throw GraphException("error running dynamics: "
                                 "number of seeds too high:'"+
                                 lexical_cast<string>(number_of_seeds) + "'");
        enum states { S, I, R, NUMBER_OF_STATES };

        vector<typename graph_traits<Graph>::vertex_descriptor> descriptors;
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            descriptors.push_back(*vi);

        typedef tr1::uniform_int<size_t> idist_t;
        tr1::variate_generator<rng_t, idist_t> random(rng, idist_t());
        while(number_of_seeds)
        {
            vector<typename graph_traits<Graph>::vertex_descriptor>
                sample(1), neighbor(1);
            __gnu_cxx::random_sample(descriptors.begin(), descriptors.end(),
                                     sample.begin(), sample.end(), random);
            for (int count = 0; count < 5; ++count){
                while( out_degree( sample[0], g) == 0 )
                {
                    __gnu_cxx::random_sample(descriptors.begin(),
                                             descriptors.end(), sample.begin(),
                                             sample.end(), random);
                    count = 0;
                }
                typename graph_traits<Graph>::adjacency_iterator ni, ni_end;
                tie(ni,ni_end) = adjacent_vertices(sample[0], g);
                __gnu_cxx::random_sample(ni, ni_end, neighbor.begin(),
                                         neighbor.end(), random);
                sample[0] = neighbor[0];
            }
            if(sir_state[sample[0]] != I )
            {
                sir_state[sample[0]] = I;
                --number_of_seeds;
            }
        }
    }
};

struct starwalker_seed
{
    typedef void result_type;
    template<class Graph, class SIRState, class NumberOfSeeds>
    void operator()(const Graph& g, SIRState& sir_state,
                    NumberOfSeeds number_of_seeds, rng_t& rng) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        if(number_of_seeds > HardNumVertices()(g))
            throw GraphException("error running dynamics: "
                                 "number of seeds too high:'"+
                                 lexical_cast<string>(number_of_seeds) + "'");
        enum states { S, I, R, NUMBER_OF_STATES };

        vector<vertex_t> descriptors;
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            descriptors.push_back(*vi);

        typedef tr1::uniform_int<size_t> idist_t;
        tr1::variate_generator<rng_t, idist_t> random(rng, idist_t());
        while(number_of_seeds)
        {
            vector<vertex_t>
                sample(1), neighbor(1);
            __gnu_cxx::random_sample(descriptors.begin(), descriptors.end(),
                                     sample.begin(), sample.end(), random);
            for (int count = 0; count < 5; ++count){
                while( out_degree( sample[0], g) == 0 )
                {
                    __gnu_cxx::random_sample(descriptors.begin(),
                                             descriptors.end(), sample.begin(),
                                             sample.end(), random);
                    count = 0;
                }
                typename graph_traits<Graph>::adjacency_iterator ni, ni_end;
                tie(ni,ni_end) = adjacent_vertices(sample[0], g);
                __gnu_cxx::random_sample(ni, ni_end, neighbor.begin(),
                                         neighbor.end(), random);
                sample[0] = neighbor[0];
            }
            if(sir_state[sample[0]] != I )
            {
                sir_state[sample[0]] = I;
                --number_of_seeds;
            }

            vector<vertex_t> star;
            typename graph_traits<Graph>::adjacency_iterator ni, ni_end;
            for(tie(ni,ni_end) = adjacent_vertices(sample[0], g); ni != ni_end;
                ++ni)
                star.push_back(*ni);
            random_shuffle(star.begin(), star.end(), random);
            typename vector<vertex_t>::iterator starwalker = star.begin();
            for(int count = 0; (starwalker != star.end()) && (count < 5);
                ++starwalker)
            {
                if(number_of_seeds)
                {
                    if(sir_state[*starwalker] != I )
                    {
                        sir_state[*starwalker] = I;
                        --number_of_seeds;
                        ++count;
                    }
                }
            }
        }
    }
};

struct bfs_seed
{
    typedef void result_type;
    // this will abort the BFS search when no longer useful
    class bfs_stop_exception {};
    template <class SIRState, class NumberOfSeeds>
    struct bfs_seeder
    {
        typedef on_tree_edge event_filter;

        bfs_seeder(SIRState* sir_state, NumberOfSeeds* number_of_seeds)
            : _sir_state(sir_state), _number_of_seeds(*number_of_seeds) {}
        
        template <class Graph>
        void operator()(typename graph_traits<Graph>::edge_descriptor e,
                        const Graph& g) 
        {
            enum states { S, I, R, NUMBER_OF_STATES };
            typename graph_traits<Graph>::vertex_descriptor v = target(e,g);
            if ( ((*_sir_state)[v] != I) && _number_of_seeds)
            {
                (*_sir_state)[v] = I;
                --_number_of_seeds;
            }
            if (!_number_of_seeds)
                throw bfs_stop_exception();
        }
        
        SIRState* _sir_state;
        NumberOfSeeds& _number_of_seeds;
    };
    template<class Graph, class SIRState, class NumberOfSeeds, class IndexMap>
    void operator()(const Graph& g, SIRState& sir_state,
                    NumberOfSeeds number_of_seeds, rng_t& rng,
                    IndexMap vertex_index) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        if(number_of_seeds > HardNumVertices()(g))
            throw GraphException("error running dynamics: "
                                 "number of seeds too high:'"+
                                 lexical_cast<string>(number_of_seeds) + "'");
        enum states { S, I, R, NUMBER_OF_STATES };

        vector<vertex_t> descriptors;
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            descriptors.push_back(*vi);

        typedef tr1::uniform_int<size_t> idist_t;
        tr1::variate_generator<rng_t, idist_t> random(rng, idist_t());
        while(number_of_seeds)
        {
            vector<vertex_t>
                sample(1), neighbor(1);
            __gnu_cxx::random_sample(descriptors.begin(), descriptors.end(),
                                     sample.begin(), sample.end(), random);
            for (int count = 0; count < 5; ++count){
                while( out_degree( sample[0], g) == 0 )
                {
                    __gnu_cxx::random_sample(descriptors.begin(),
                                             descriptors.end(), sample.begin(),
                                             sample.end(), random);
                    count = 0;
                }
                typename graph_traits<Graph>::adjacency_iterator ni, ni_end;
                tie(ni,ni_end) = adjacent_vertices(sample[0], g);
                __gnu_cxx::random_sample(ni, ni_end, neighbor.begin(),
                                         neighbor.end(), random);
                sample[0] = neighbor[0];
            }
            if( (sir_state[sample[0]] != I) && number_of_seeds )
            {
                sir_state[sample[0]] = I;
                --number_of_seeds;
            }

            // BFS
            typedef tr1::unordered_map<vertex_t,default_color_type,
                                       DescriptorHash<IndexMap> > cmap_t;
            cmap_t cmap(0, DescriptorHash<IndexMap>(vertex_index));
            InitializedPropertyMap<cmap_t>
                color_map(cmap, color_traits<default_color_type>::white());

            try
            {
                bfs_seeder<SIRState, NumberOfSeeds> seeder(&sir_state,
                                                           &number_of_seeds);
                breadth_first_visit(g, sample[0],
                                    visitor(make_bfs_visitor(seeder)).
                                    color_map(color_map));
            }
            catch(bfs_stop_exception) {}
        }
    }
};

} //graph-tool namespace

#endif // GRAPH_STRATEGIC_SEEDING_HH
