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

    debugl(4, "IdQueue::getID()\n");

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
    debugl(4, "IdQueue::freeId()\n");
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
        debugl(0, "IdQueue::freeId(): WARNING: attempting to free id (%5d), which is smaller than this->smallest_id (%5d). ignoring..\n", id, this->smallest_id);
    }
}
