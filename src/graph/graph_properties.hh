// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2007  Tiago de Paula Peixoto <tiago@forked.de>
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

#ifndef GRAPH_PROPERTIES_HH
#define GRAPH_PROPERTIES_HH

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <tr1/memory>
#include <boost/python/object.hpp>

#include <boost/property_map.hpp>
#include <boost/dynamic_property_map.hpp>
#include "fast_vector_property_map.hh"
#include <boost/functional/hash.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/find.hpp>
#include <boost/lambda/bind.hpp>

#include "graph.hh"

// this file provides general functions for manipulating graph properties

namespace graph_tool
{
using namespace std;
using namespace boost;

// Metaprogramming
// ===============
//
// Metafunctions and data structures to deal with property maps

// Global Property Types
// ---------------------
// global property types. only these types are allowed in property maps
// Note: we must avoid a vector<bool> (and bools in general) since it is quite
//       broken, and use a vector<uint8_t> instead!
//       see: http://www.gotw.ca/publications/N1211.pdf

typedef mpl::vector<uint8_t, int32_t, int64_t, double, long double, string,
                    vector<uint8_t>, vector<int32_t>, vector<int64_t>,
                    vector<double>, vector<long double>, vector<string>,
                    python::object>
    value_types;

extern const char* type_names[]; // respective type names (defined in
                                 // graph_properties.cc)

// scalar types: types contained in value_types which are scalar
typedef mpl::vector<uint8_t, int32_t, int64_t, double, long double>
    scalar_types;

// integer_types: scalar types which are integer
typedef mpl::vector<uint8_t, int32_t, int64_t> integer_types;

// floating_types: scalar types which are floating point
typedef mpl::vector<double, long double> floating_types;

struct make_vector
{
    template <class ValueType> struct apply { typedef vector<ValueType> type; };
};

// scalar_vector_types: vector types with floating point values
typedef mpl::transform<scalar_types,make_vector>::type scalar_vector_types;

// integer_vector_types: vector types with floating point values
typedef mpl::transform<integer_types,make_vector>::type integer_vector_types;

// floating_vector_types: vector types with floating point values
typedef mpl::transform<floating_types,make_vector>::type floating_vector_types;

//
// Property Map Types
// ------------------

// metafunction to generate the correct property map type given a value type and
// an index map
struct property_map_type
{
    template <class ValueType, class IndexMap>
    struct apply
    {
        typedef checked_vector_property_map<ValueType,IndexMap> type;
    };
};

// metafunction to get the sequence of property map types of ValueTypes and
// IndexMap
struct property_map_types
{
   // this wraps an unary metafunction class Bind into a unary metafunction,
   // i.e., it is an identity operation. I'm not sure why it's necessary, but
   // using pure unary bind expressions below didn't work for me, and this fixed
   // it.
    template <class Bind>
    struct bind_wrap1
    {
        template <class T1> struct apply
        { typedef typename Bind::template apply<T1>::type type; };
    };

    template <class ValueTypes, class IndexMap,
              class IncludeIndexMap = mpl::bool_<true> >
    struct apply
    {
        typedef typename mpl::transform<
            ValueTypes,
            bind_wrap1<mpl::bind2<property_map_type,
                                  mpl::_1,
                                  IndexMap> >
            >::type scalar_properties;

        // put index map itself
        typedef typename mpl::if_<
            IncludeIndexMap,
            typename mpl::push_back<scalar_properties,IndexMap>::type,
            scalar_properties
            >::type type;
    };
};

// Property map manipulation
// =========================
//
// Functions which deal with several aspects of property map manipulation


// this functor tests whether or not a given boost::any object holds a type
// contained in a given type Sequence
template <class Sequence>
struct belongs
{
    struct get_type
    {
        get_type(const boost::any& val, bool& found)
            : _val(val), _found(found) {}

        template <class Type>
        void operator()(Type) const
        {
            const Type* ptr = any_cast<Type>(&_val);
            if (ptr != 0)
                _found = true;
        }

        const boost::any& _val;
        bool& _found;
    };

    bool operator()(const boost::any& prop)
    {
        bool found = false;
        mpl::for_each<Sequence>(get_type(prop, found));
        return found;
    }
};

// this will return the name of a given type
template <class TypeSequence = value_types,
          class NamedSequence = value_types>
class get_type_name
{
public:
    get_type_name(const char* names[] = type_names)
        : _type_names(type_names)
    {
        if (_all_names.empty())
        {
            using namespace lambda;
            mpl::for_each<TypeSequence>
                (lambda::bind<void>(get_all_names(), lambda::_1,
                                    var(_type_names), var(_all_names)));
        }
    }

    const string& operator()(const type_info& type) const
    {
        using namespace lambda;
        string* name;
        mpl::for_each<TypeSequence>
            (lambda::bind<void>(find_name(), lambda::_1, var(type),
                                var(_all_names), var(name)));
        return *name;
    }

    const vector<string>& all_names() const
    {
        return _all_names;
    }

private:
    struct find_name
    {
        template <class Type>
        void operator()(Type, const type_info& type,
                        vector<string>& all_names,
                        string*& name) const
        {
            size_t index = mpl::find<TypeSequence,Type>::type::pos::value;
            if (type == typeid(Type))
                name = &all_names[index];
        }
    };

    struct get_all_names
    {
        template <class Type>
        void operator()(Type, const char** type_names,
                        vector<string>& names) const
        {
            size_t index = mpl::find<NamedSequence,Type>::type::pos::value;
            names.push_back(type_names[index]);
        }
    };

    const char** _type_names;
    static vector<string> _all_names;
};

template <class TypeSequence, class NamedSequence>
vector<string> get_type_name<TypeSequence,NamedSequence>::_all_names;

//
// Extra Property Map Types
// ========================

// the following class wraps a generic property map, so it can be used as a
// property with a given Key and Value type. The keys and values are converted
// to the desired Key and Value type, which may cause a performance impact,
// since virtual functions are used. Should be used only when property map
// access time is not crucial
template <class Value, class Key>
class DynamicPropertyMapWrap
{
public:
    typedef Value value_type;
    typedef Value reference;
    typedef Key key_type;
    typedef read_write_property_map_tag category;

    template <class PropertyTypes>
    DynamicPropertyMapWrap(boost::any pmap, PropertyTypes)
    {
        ValueConverter* converter = 0;
        mpl::for_each<PropertyTypes>
            (lambda::bind<void>(choose_converter(), lambda::_1,
                                lambda::var(pmap), lambda::var(converter)));
        if (converter == 0)
            throw bad_lexical_cast();
        else
            _converter = tr1::shared_ptr<ValueConverter>(converter);
    }

    DynamicPropertyMapWrap() {}

    Value get(const Key& k) const
    {
        return (*_converter).get(k);
    }

    void put(const Key& k, const Value& val)
    {
        (*_converter).put(k, val);
    }

private:
    class ValueConverter
    {
    public:
        virtual Value get(const Key& k) = 0;
        virtual void put(const Key& k, const Value& val) = 0;
    };

    template <class PropertyMap>
    class ValueConverterImp: public ValueConverter
    {
    public:
        ValueConverterImp(PropertyMap pmap): _pmap(pmap) {}

        virtual Value get(const Key& k)
        {
            typedef typename property_traits<PropertyMap>::key_type key_t;
            return boost::get(_pmap, key_t(k));
        }

        virtual void put(const Key& k, const Value& val)
        {
            typedef typename property_traits<PropertyMap>::key_type key_t;
            typedef typename property_traits<PropertyMap>::value_type val_t;
            boost::put(_pmap, key_t(k), val_t(val));
        }

    private:
        PropertyMap _pmap;
    };

    struct choose_converter
    {
        template <class PropertyMap>
        void operator()(PropertyMap, boost::any& dmap,
                        ValueConverter*& converter) const
        {
            if (typeid(PropertyMap) == dmap.type())
                converter = new ValueConverterImp<PropertyMap>
                    (any_cast<PropertyMap>(dmap));
        }
    };

    tr1::shared_ptr<ValueConverter> _converter;
};

template <class Value, class Key, class ConvKey>
Value get(const graph_tool::DynamicPropertyMapWrap<Value,Key>& pmap,
          ConvKey k)
{
    Key key = k;
    return pmap.get(key);
}

template <class Value, class Key>
void put(graph_tool::DynamicPropertyMapWrap<Value,Key>& pmap,
         Key k, const Value& val)
{
    pmap.put(k, val);
}

// the following is hash functor which, will hash a vertex or edge descriptor
// based on its index
template <class IndexMap>
class DescriptorHash
    : public unary_function<typename IndexMap::key_type, size_t>
{
public:
    DescriptorHash() {}
    DescriptorHash(IndexMap index_map): _index_map(index_map) {}
    size_t operator()(typename IndexMap::key_type const& d) const
    {
        return hash_value(_index_map[d]);
    }
private:
    IndexMap _index_map;
};

// the following is a property map based on a hashed container, which uses the
// above hash function for vertex or edge descriptors
template <class IndexMap, class Value>
class HashedDescriptorMap
    : public put_get_helper<Value&, HashedDescriptorMap<IndexMap,Value> >
{
public:
    typedef DescriptorHash<IndexMap> hashfc_t;
    typedef tr1::unordered_map<typename IndexMap::key_type,Value,hashfc_t>
        map_t;
    typedef associative_property_map<map_t> prop_map_t;

    typedef typename property_traits<prop_map_t>::value_type value_type;
    typedef typename property_traits<prop_map_t>::reference reference;
    typedef typename property_traits<prop_map_t>::key_type key_type;
    typedef typename property_traits<prop_map_t>::category category;

    HashedDescriptorMap(IndexMap index_map)
        : _base_map(new map_t(0, hashfc_t(index_map))), _prop_map(*_base_map) {}
    HashedDescriptorMap(){}

    reference operator[](const key_type& k) { return _prop_map[k]; }
    const reference operator[](const key_type& k) const { return _prop_map[k]; }

private:
    shared_ptr<map_t> _base_map;
    prop_map_t _prop_map;
};


// this wraps a container as a property map which is automatically initialized
// with a given default value
template <class Container>
class InitializedPropertyMap
    : public put_get_helper<typename Container::value_type::second_type&,
                            InitializedPropertyMap<Container> >
{
public:
    typedef typename Container::value_type::second_type value_type;
    typedef value_type& reference;
    typedef typename Container::key_type key_type;
    typedef read_write_property_map_tag category;

    InitializedPropertyMap(Container& base_map, value_type def)
        : _base_map(&base_map), _default(def) {}
    InitializedPropertyMap(){}

    reference operator[](const key_type& k)
    {
        return get(k);
    }

    const reference operator[](const key_type& k) const
    {
        return get(k);
    }

    const reference get(const key_type& k) const
    {
        typename Container::iterator val;
        val = _base_map->find(k);
        if (val == _base_map->end())
            val = _base_map->insert(make_pair(k, _default)).first;
        return val->second;
    }

private:
    Container* _base_map;
    value_type _default;
};

// the following is a property map which always returns a constant value
template <class Value, class Key>
class ConstantPropertyMap
    : public put_get_helper<Value, ConstantPropertyMap<Value,Key> >
{
public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef Key key_type;
    typedef readable_property_map_tag category;

    ConstantPropertyMap(const value_type& c): _c(c) {}
    ConstantPropertyMap(){}

    const value_type& operator[](const key_type& k) const { return _c; }

private:
    value_type _c;
};

// this wraps an existing property map, but always converts its values to a
// given type
template <class PropertyMap, class Type>
class ConvertedPropertyMap
{
public:
    typedef Type value_type;
    typedef typename property_traits<PropertyMap>::value_type orig_type;
    typedef value_type reference;
    typedef typename property_traits<PropertyMap>::key_type key_type;
    typedef read_write_property_map_tag category;

    ConvertedPropertyMap(PropertyMap base_map)
        : _base_map(base_map) {}
    ConvertedPropertyMap(){}

    value_type get(const key_type& k) const
    {
        return boost::get(_base_map, k);
    }

    void put(const key_type& k, const Type& v)
    {
        put(_base_map, k, orig_type(v));
    }
private:
    PropertyMap _base_map;
};

template <class PropertyMap, class Type>
Type get(const graph_tool::ConvertedPropertyMap<PropertyMap,Type>& pmap,
         const typename property_traits<PropertyMap>::key_type& k)
{
    return pmap.get(k);
}

template <class PropertyMap, class Type>
void put(graph_tool::ConvertedPropertyMap<PropertyMap,Type> pmap,
         const typename property_traits<PropertyMap>::key_type& k,
         const typename property_traits<PropertyMap>::value_type& val)
{
    pmap.put(k, val);
}

} // graph_tool namespace

#endif
