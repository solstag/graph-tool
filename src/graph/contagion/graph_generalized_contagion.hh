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

#ifndef GRAPH_GENERALIZED_CONTAGION_HH
#define GRAPH_GENERALIZED_CONTAGION_HH

#include <tr1/unordered_map>
#include <algorithm>
#include <ext/algorithm>
#include <numeric>
#include <set>
#include <deque>
#include <vector>
#include <string>
#include <sstream>
#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>

#include <tr1/random>
typedef std::tr1::mt19937 rng_t;

namespace graph_tool
{
using namespace std;
using namespace boost;

// This will iterate over a random permutation of a random access sequence, by
// swapping the values of the sequence as it iterates
template <class RandomAccessIterator, class RNG,
          class RandomDist = tr1::uniform_int<size_t> >
class random_permutation_iterator : public
    std::iterator<input_iterator_tag, typename RandomAccessIterator::value_type>
{
public:
    random_permutation_iterator(RandomAccessIterator begin,
                                RandomAccessIterator end, RNG& rng)
        : _i(begin), _end(end), _rng(&rng)
    {
        if(_i != _end)
        {
            RandomDist random(0,  _end - _i - 1);
            std::iter_swap(_i, _i + random(*_rng));
        }
    }

    typename RandomAccessIterator::value_type operator*()
    {
        return *_i;
    }

    random_permutation_iterator& operator++()
    {
        ++_i;
        if(_i != _end)
        {
            RandomDist random(0,  _end - _i - 1);
            std::iter_swap(_i, _i + random(*_rng));
        }
        return *this;
    }

    bool operator==(const random_permutation_iterator& ri)
    {
        return _i == ri._i;
    }

    bool operator!=(const random_permutation_iterator& ri)
    {
        return _i != ri._i;
    }
private:
    RandomAccessIterator _i, _end;
    RNG* _rng;
};

// Gives us some choices (currently not) of infection doeses, implemented in C++
struct sample_infection_dose
{
    sample_infection_dose(string inf_dose_dist, rng_t& rng) : _rng(rng)
    {
        if( inf_dose_dist == "constant_one" ){ _sampler = constant_one(); }
          else throw GraphException("error running dynamics: requested invalid"
                                    "infection dose distribution '"+
                                    inf_dose_dist + "'");
    }
    double operator()()
    {
        return _sampler();
    }
    struct constant_one
    {
        double operator()() { return 1.0; }
    };
    function<double ()> _sampler;
    rng_t& _rng;
};

// to be used with 'for_each'
void push_bback(vector<size_t>& vec)
{
    vec.push_back(vec.back());
}

template <class DoseThreshold, class Doses, class AccumulatedDose,
          class SIRState, class Availability >
struct advance_contagion
{
    advance_contagion(vector<vector<size_t> >* sm, DoseThreshold* dt, Doses* ds,
                      AccumulatedDose* ad, SIRState* ss, Availability* av,
                      rng_t& rng)
        :state_members(*sm), dose_threshold(*dt), doses(*ds),
         accumulated_dose(*ad), sir_state(*ss), availability(*av), _rng(rng) {}

    template < class Graph, class InfectionRate, class RecoveryRate,
               class ResusceptibilityRate, class InfectionDoseDistribution >
    void operator()(const Graph& g, InfectionRate infection_rate, RecoveryRate
                    recovery_rate, ResusceptibilityRate resusceptibility_rate,
                    InfectionDoseDistribution infection_dose_distribution,
                    int max_iterations, int timescale ) const
    {
        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
        tr1::unordered_map<typename property_traits<SIRState>::key_type,
                           typename property_traits<SIRState>::value_type>
            next_state;
        enum states { S, I, R, NUMBER_OF_STATES };
        typename graph_traits<Graph>::vertex_iterator vi, vi_end;
        int stability_index = 0, max_stability = 100;
        int max_availability, rounds_in_turn;
        typedef tr1::uniform_real<double> rdist_t;
        tr1::variate_generator<rng_t, rdist_t> sampler_01(_rng, rdist_t());
        sample_infection_dose
            infection_dose_sampler(infection_dose_distribution, _rng);
        
        if(timescale > 0)
            { max_availability = timescale; rounds_in_turn = 1; }
        else if(timescale < 0)
            { max_availability = 1; rounds_in_turn = -timescale; }
        else if(timescale == 0)
            return;
        
        for_each(state_members.begin(), state_members.end(),
                 lambda::bind(&vector<size_t>::push_back, lambda::_1, 0));
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            ++state_members[sir_state[*vi]].back();
        
        typedef random_permutation_iterator<typename vector<vertex_t>::iterator,
                                            rng_t> random_vertex_iter;
        vector<vertex_t> all_vertices;
        for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
            all_vertices.push_back(*vi);

        for(int iterations = 0; iterations < max_iterations; )
        {
            random_vertex_iter
                svi(all_vertices.begin(), all_vertices.end(), _rng),
                svi_end(all_vertices.end(), all_vertices.end(), _rng);
            for (; svi != svi_end; ++svi)
            {
                accumulated_dose[*svi] -= doses[*svi].front();
                doses[*svi].pop_front();
                vector<vertex_t> neighbors;
                typename graph_traits<Graph>::adjacency_iterator ni, ni_end;
                for(tie(ni,ni_end)=adjacent_vertices(*svi, g); ni!=ni_end; ++ni)
                    neighbors.push_back(*ni);
                random_vertex_iter
                    rni(neighbors.begin(), neighbors.end(), _rng),
                    rni_end(neighbors.end(), neighbors.end(), _rng);
                rni = find_if(rni, rni_end, is_available(availability));
                switch(sir_state[*svi])
                {
                    case S:
                        if(rni != rni_end)
                        {
                            if((sir_state[*rni] == I) &&
                               (sampler_01() < infection_rate))
                            {
                                doses[*svi].push_back(infection_dose_sampler());
                                accumulated_dose[*svi]+=doses[*svi].back();
                                if(accumulated_dose[*svi] >=
                                   dose_threshold[*svi])
                                    next_state[*svi]=I;
                            }
                            else doses[*svi].push_back(0.0);

                            --availability[*rni];
                        }
                        else doses[*svi].push_back(0.0);
                        break;

                    case I:
                        if(rni != rni_end)
                        {
                            if((sir_state[*rni] == I) &&
                               (sampler_01() < infection_rate))
                            {
                                doses[*svi].push_back(infection_dose_sampler());
                                accumulated_dose[*svi]+=doses[*svi].back();
                            }
                            else doses[*svi].push_back(0.0);

                            --availability[*rni];
                        }
                        else doses[*svi].push_back(0.0);

                        if(accumulated_dose[*svi] < dose_threshold[*svi] &&
                           sampler_01() < recovery_rate )
                            next_state[*svi] =
                                (resusceptibility_rate > 1.0) ? S : R ;
                                // SIS or SIR model             ^^^^^
                        break;

                    case R:
                        if(rni != rni_end)
                            --availability[*rni];
                        doses[*svi].push_back(0.0);
                        if( sampler_01() < resusceptibility_rate )
                            next_state[*svi]=S;
                        break;

                    default:
                        throw GraphException
                                  ("error running dynamics: unknown state '" +
                                   lexical_cast<string>(sir_state[*svi]) + "'");
                        
                } // vertex state switch
            } // graph loop
            
            for_each(state_members.begin(), state_members.end(), push_bback);
            for(typeof(next_state.begin()) node_state = next_state.begin();
                node_state != next_state.end(); ++node_state)
            {
                --(state_members[ sir_state[node_state->first] ].back());
                sir_state[node_state->first] = node_state->second;
                ++(state_members[ sir_state[node_state->first] ].back());
            }
            next_state.clear();

            if(state_members[I].back() ==
               state_members[I][state_members[I].size()-2])
                ++stability_index;
            else stability_index = 0;
            
            ++iterations;
            
            if( !(iterations % rounds_in_turn) )
                for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
                    availability[*vi]=max_availability;
            
            if((state_members[I].back() == 0 && state_members[R].back() == 0) ||
               stability_index == max_stability)
                break;
            
        } // iteration loop
    } // operator()
    
    struct is_available
    {
        is_available(Availability& a) : _a(a){}
        template<class Vertex>
        bool operator()(const Vertex& v)
        {
            return _a[v] > 0;
        }
        Availability _a;
    };
    
    vector<vector<size_t> >& state_members;
    DoseThreshold& dose_threshold;
    Doses& doses;
    AccumulatedDose& accumulated_dose;
    SIRState& sir_state;
    Availability& availability;
    rng_t& _rng;
}; // advance_contagion

} //graph-tool namespace

#endif // GRAPH_GENERALIZED_CONTAGION_HH
