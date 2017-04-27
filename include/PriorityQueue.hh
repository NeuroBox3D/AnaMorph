/*
 * This file is part of
 *
 * AnaMorph: a framework for geometric modelling, consistency analysis and surface
 * mesh generation of anatomically reconstructed neuron morphologies.
 * 
 * Copyright (c) 2013-2017: G-CSC, Goethe University Frankfurt - Queisser group
 * Author: Konstantin Mörschel
 * 
 * AnaMorph is free software: Redistribution and use in source and binary forms,
 * with or without modification, are permitted under the terms of the
 * GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 *
 * (3) Neither the name "AnaMorph" nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * (4) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Mörschel K, Breit M, Queisser G. Generating neuron geometries for detailed
 *   three-dimensional simulations using AnaMorph. Neuroinformatics (2017)"
 * "Grein S, Stepniewski M, Reiter S, Knodel MM, Queisser G.
 *   1D-3D hybrid modelling – from multi-compartment models to full resolution
 *   models in space and time. Frontiers in Neuroinformatics 8, 68 (2014)"
 * "Breit M, Stepniewski M, Grein S, Gottmann P, Reinhardt L, Queisser G.
 *   Anatomically detailed and large-scale simulations studying synapse loss
 *   and synchrony using NeuroBox. Frontiers in Neuroanatomy 10 (2016)"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PRIORITY_QUEUE
#define PRIORITY_QUEUE

#include "common.hh"
#include "IdQueue.hh"

template <typename Tkey, typename Tvalue>
class PriorityQueue {
    private:
        IdQueue                                                                             idq;
        std::map<std::pair<Tkey, uint32_t>, Tvalue>                                         Q;
        std::map<Tvalue, typename std::map< std::pair<Tkey, uint32_t>, Tvalue>::iterator>   value_key_map;

    public:
        void                        insert(Tkey key, Tvalue value);
        void                        insert(std::pair<Tkey, Tvalue> key_value_pair);
        std::pair<Tkey, Tvalue>     top();
        void                        deleteMin();
        bool                        changeKey(Tvalue value, Tkey new_key);
        void                        clear();
        bool                        empty();
        void                        checkHeap();
};

template<typename Tkey, typename Tvalue>
void
PriorityQueue<Tkey, Tvalue>::clear()
{
    this->Q.clear();
    this->value_key_map.clear();
}

template<typename Tkey, typename Tvalue>
bool
PriorityQueue<Tkey, Tvalue>::empty()
{
    return (this->Q.empty());
}

template<typename Tkey, typename Tvalue>
void
PriorityQueue<Tkey, Tvalue>::insert(
    Tkey    key,
    Tvalue  value)
{
    typename std::map<std::pair<Tkey, uint32_t>, Tvalue>::iterator qit;

    /* insert (key, value) pair and store iterator */
    qit = Q.insert( std::pair< std::pair<Tkey, uint32_t>, Tvalue>( std::pair<Tkey, uint32_t>(key, idq.getId()), value )).first;

    /* associate key with iterator in value_key_map */
    value_key_map.insert( std::pair<Tvalue, typename std::map< std::pair<Tkey, uint32_t>, Tvalue>::iterator> (value, qit) );
}

template<typename Tkey, typename Tvalue>
void
PriorityQueue<Tkey, Tvalue>::insert(
    std::pair<Tkey, Tvalue> key_value_pair)
{
    this->insert(key_value_pair.first, key_value_pair.second);
}

template<typename Tkey, typename Tvalue>
std::pair<Tkey, Tvalue>
PriorityQueue<Tkey, Tvalue>::top()
{
    std::pair< std::pair<Tkey, uint32_t>, Tvalue>  p(Q.begin()->first, Q.begin()->second);
    return (std::pair<Tkey, Tvalue>(p.first.first, p.second));
}

template<typename Tkey, typename Tvalue>
void
PriorityQueue<Tkey, Tvalue>::deleteMin()
{
    std::pair< std::pair<Tkey, uint32_t>, Tvalue>  p(Q.begin()->first, Q.begin()->second);

    /* free id */
    idq.freeId(p.first.second);

    Q.erase( Q.begin() );
    value_key_map.erase(p.second);
}

template<typename Tkey, typename Tvalue>
bool
PriorityQueue<Tkey, Tvalue>::changeKey(
        Tvalue  value,
        Tkey    new_key)
{
    uint32_t id;

    auto it = value_key_map.find(value);
    if ( it != value_key_map.end()) {
        debugl(1, "PriorityQueue()::changeKey(): value found => changing key..\n");
        debugTabInc();

        id          = it->second->first.second;
        debugl(1, "id: %d\n", id);

        /* erase old entry from q with obtained iterator. insert new entry into Q and save iterator */
        debugl(1, "erasing iterator from key_value_map Q.\n");
        Q.erase(it->second);

        debugl(1, "re-inserting into Q with correct key.\n");
        auto qpair  = Q.insert( { {new_key, id}, value } );

        debugl(1, "updating key in value_key_map.\n");
        /* only key has changed, value is unaffected => use iterator it to update entry in value_key_map: overwrite the
         * currently stored old iterator with qpair.first, which has been returned by the Q.insert(..) call above */
        it->second  = qpair.first;

        debugl(1, "done. returning..\n");
        debugTabDec();
        return true;
    }
    else {
        debugl(1, "PriorityQueue()::changeKey(): value not found..\n");
        return false;
    }
}


template<typename Tkey, typename Tvalue>
void
PriorityQueue<Tkey, Tvalue>::checkHeap()
{
    std::list<Tvalue>                                                                   Q_values, value_key_map_keys;
    typename std::map<std::pair<Tkey, uint32_t>, Tvalue>::iterator                      qit;
    typename std::map<uint32_t, typename std::map<Tkey, Tvalue>::iterator>::iterator    vit;

    for (qit = Q.begin(); qit != Q.end(); ++qit) {
        Q_values.push_back(qit->second.first);
    }

    for (vit = value_key_map.begin(); vit != value_key_map.end(); ++vit) {
        value_key_map_keys.push_back(vit->first);
    }

    Q_values.sort();
    value_key_map_keys.sort();
    if (Q_values != value_key_map_keys) {
        throw("PriorityQueue::checkHeap(): heap invariant violated. internal logic error.");
    }
}

#endif
