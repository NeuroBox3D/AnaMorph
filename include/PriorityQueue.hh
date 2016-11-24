/* --------------------------------------------------------------------------------
 * 
 *                              THIS FILE IS PART OF                               
 * 
 * AnaMorph: A Framework for Geometric Modelling, Consistency Analysis and Surface
 * Mesh Generation of Anatomically Reconstructed Neuron Morphologies
 * 
 * Web site: http://www.anamorph.net
 * 
 * Copyright (c) 2013-2014, Konstantin Mörschel.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 * 
 *    This product includes software developed by Konstantin Mörschel for
 *    AnaMorph (http://www.anamorph.net).
 * 
 * 4. Neither the name "AnaMorph" nor the names of its contributors may be
 *    used to endorse or promote products derived from this software without
 *    specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS OF ANAMORPH ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS OF ANAMORPH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * -------------------------------------------------------------------------------- */

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
        debugl(0, "PriorityQueue()::changeKey(): value found => changing key..\n");
        debugTabInc();

        id          = it->second->first.second;
        debugl(0, "id: %d\n", id);

        /* erase old entry from q with obtained iterator. insert new entry into Q and save iterator */
        debugl(0, "erasing iterator from key_value_map Q.\n");
        Q.erase(it->second);

        debugl(0, "re-inserting into Q with correct key.\n");
        auto qpair  = Q.insert( { {new_key, id}, value } );

        debugl(0, "updating key in value_key_map.\n");
        /* only key has changed, value is unaffected => use iterator it to update entry in value_key_map: overwrite the
         * currently stored old iterator with qpair.first, which has been returned by the Q.insert(..) call above */
        it->second  = qpair.first;

        debugl(0, "done. returning..\n");
        debugTabDec();
        return true;
    }
    else {
        debugl(0, "PriorityQueue()::changeKey(): value not found..\n");
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
