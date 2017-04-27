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

#include "common.hh"
#include "IdQueue.hh"

IdQueue::IdQueue(uint32_t _smallest_id)
: refill_count(1024), smallest_id(_smallest_id), next_id(_smallest_id), last_id(_smallest_id)
{}

IdQueue::IdQueue()
: refill_count(1024), smallest_id(0), next_id(0), last_id(0)
{}

void
IdQueue::clear(uint32_t smallest_id) {
    /* reset queue, there's no clear, so fresh copy with trivial constructor */
    this->q             =   std::priority_queue<
                                uint32_t,
                                std::deque<uint32_t>,
                                std::greater<uint32_t> > ();

    this->smallest_id   = smallest_id;
    this->next_id       = smallest_id;
    this->last_id       = smallest_id;
    this->refill_count  = 1024;
}

uint32_t
IdQueue::getId()
{
    uint32_t    id;

    debugl(5, "IdQueue::getID()\n");

    debugTabInc();
    /* if q is empty, refill it with larger ids. note that next_id == last_id here.. */
    if (this->q.empty()) {
        debugl(5, "q.size(): %10d, refilling id queue.. last_id: %10d, refill count: %10d \n", this->q.size(), this->last_id, this->refill_count);

        debugTabInc();
        uint32_t i;
        for (i = this->last_id + 1; i < (this->last_id + this->refill_count + 1) && i < UINT32_MAX; i++) {
            debugl(6, "pushing id %5d\n", i);
            this->q.push(i);
        }
        debugTabDec();

        /* value of i after the loop is the last pushed id 1 => update last_id = i - 1. note that it is possible
         * that next_id == last_id was (UINT32_MAX - 1) and no id has been pushed at all, where in this case last_id
         * retains its value of (UINT32_MAX - 1). */
        this->last_id = i - 1;

        /* if q is still empty here, no new id could be pushed and there was no "refilling". this happens if and only
         * if next_id == last_id == UINT32_MAX - 1 at the beginning of the call. throw exception to indicate overflow to
         * the caller. */
        if (q.empty()) {
            throw("Q could not be refilled, since UINT32_MAX - 1 values have already been used => overflow.");
        }

        debugl(5, "q.size(): %10d, last_id: %10d.\n", this->q.size(), this->last_id);
    }
    debugTabDec();

    /* save next_id in id, update next_id with min (top()) element of q and deleteMin (pop()) of q */
    id              = this->next_id;
    this->next_id   = this->q.top();
    this->q.pop();

    return id;
}

void
IdQueue::freeId(uint32_t id)
{
    debugl(5, "IdQueue::freeId()\n");
    if (id >= this->smallest_id) {
        /* if id < next_id, push next_id and set next_id = id */
        if (id < next_id) {
            this->q.push(this->next_id);
            this->next_id = id;
        }
        /* otherwise push id */
        else this->q.push(id);
    }
    else {
        debugl(1, "IdQueue::freeId(): WARNING: attempting to free id (%5d), which is smaller than this->smallest_id (%5d). ignoring..\n", id, this->smallest_id);
    }
}
