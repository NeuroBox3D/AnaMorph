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

/* initialize static member predicates */
#if 0
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
neuron_vertex_true = std::function<bool(NeuronVertex const &v)>(
    [] (NeuronVertex const &v) -> bool
    {
        return true;
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
neuron_vertex_false = std::function<bool(NeuronVertex const &v)>(
    [] (NeuronVertex const &v) -> bool
    {
        return false;
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
is_cell_vertex = std::function<bool(NeuronVertex const &v)>(
    [] (NeuronVertex const &v) -> bool
    {
        return (v.getType() == SOMA_VERTEX ||
                v.getType() == AXON_VERTEX ||
                v.getType() == DENDRITE_VERTEX);
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
is_neurite_vertex = std::function<bool(NeuronVertex const &v)>(
    [] (NeuronVertex const &v) -> bool
    {
        return (v.getType() == SOMA_VERTEX ||
                v.getType() == AXON_VERTEX ||
                v.getType() == DENDRITE_VERTEX);
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
is_axon_vertex = std::function<bool(NeuronVertex const &v)>(
    [] (NeuronVertex const &v) -> bool
    {
        return (v.getType() == AXON_VERTEX);
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
is_dendrite_vertex = std::function<bool(NeuronVertex const &v)>(
    [] (NeuronVertex const &v) -> bool
    {
        return (v.getType() == DENDRITE_VERTEX);
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
neuron_edge_true = std::function<bool(NeuronEdge const &e)>(
    [] (NeuronEdge const &e) -> bool
    {
        return true;
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
neuron_edge_false = std::function<bool(NeuronEdge const &e)>(
    [] (NeuronEdge const &e) -> bool
    {
        return false;
    }
);

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::function<bool(typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
is_neurite_segment = std::function<bool(NeuronEdge const &e)>(
    [] (NeuronEdge const &e) -> bool
    {
        return (e.getType() == AXON_VERTEX || e.getType() == DENDRITE_SEGMENT);
    }
);
#endif

#include "aux.hh"

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of (neuron) vertex class hierarchy for CellNetworks                 
 *                        starting abstract class NeuronVertex
 *
 * ----------------------------------------------------------------------------------------------------------------- */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::NeuronVertex(
    CellNetwork * const                &network,
    std::list<CellSection> const       &sections,
    Tv const                           &data) : Graph<Tn, Tv, Te>::Vertex(network, data)
{
    this->network = network;
    if (!sections.empty()) {
        this->sections      = sections;
    }
    else {
        throw("CellNetwork::NeuronVertex::NeuronVertex(): given section list empty. neuron vertices must contain information about at least one section.");
    }
}

/* protected copy constructor: use with care: this invokes Graph::Vertex copy constructor, which in turn copies POINTERS
 * to the internal data. must be used only in the implementation and care must be taken so as not to create data races
 * with pointers. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::NeuronVertex(NeuronVertex const &v)
: Graph<Tn, Tv, Te>::Vertex(v), network(NULL)

{
    this->sections      = v.sections;
}

/* see protected copy constructor above */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::operator=(NeuronVertex const &v)
{
    Graph<Tn, Tv, Te>::Vertex::operator=(v);

    this->sections      = v.sections;

    return (*this);
}

/* virtual dtor */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::~NeuronVertex()
{
}


/* NOTE: this default-implementation returns the network id of the first section, which is guaranteed to exist since
 * there are no NeuronVertices with empty section list. however, this might have to be adapted to the semantics of the
 * vertex type in question (e.g. somas, synapses, ..) */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getNetworkCompartmentId() const
{
    return (this->sections.begin()->compartment_id());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::numSections() const
{
    return (this->sections.size());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getSectionMinRadius() const
{
    R max_r = -Aux::Numbers::inf<R>();

    for (auto &s : this->sections) {
        max_r = std::max(max_r, s.getRadius());
    }
    return max_r;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getSectionMaxRadius() const
{
    R min_r = Aux::Numbers::inf<R>();

    for (auto &s : this->sections) {
        min_r = std::min(min_r, s.getRadius());
    }
    return min_r;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
Vec3<R>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getSectionPositionCentroid() const
{
    Vec3<R> sum(0, 0, 0);

    for (auto &s : this->sections) {
        sum += s.getPosition();
    }
    return (sum / (R)(this->sections.size()));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getBoundingBox(
    Vec3<R> &bb_min,
    Vec3<R> &bb_max) const
{
    using namespace Aux::VecMat;

    R       inf = Aux::Numbers::inf<R>();
    Vec3<R> sp;
    R       max_r = -inf;

    bb_min = Vec3<R>( inf,  inf,  inf);
    bb_max = Vec3<R>(-inf, -inf, -inf);

    /* get vectors of min / max components of all section positions. also get maximum radius in the loop */
    for (auto &s : this->sections) {
        max_r   = std::max(max_r, s.radius());
        sp      = s.position();
        minVec3<R>(bb_min, bb_min, sp);
        maxVec3<R>(bb_max, bb_max, sp);
    }

    /* extend bounding box by onesVec3() * max_r */
    bb_min     -= onesVec3<R>() * max_r;
    bb_max     += onesVec3<R>() * max_r;

    /* extend bounding box by 2.5%, but no less than 1E-3, on each side. */
    Vec3<R> offset;
    maxVec3<R>(offset, Vec3<R>(1E-3, 1E-3, 1E-3), fabsVec3<R>(bb_max - bb_min)*0.025);
    bb_max         += offset;
    bb_min         -= offset;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::list<typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::CellSection>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getSections() const
{
    return (this->sections);
}

/*
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::setSections(std::list<CellSection> const &sections)
{
    if (!sections.empty()) {
        this->sections = sections;
    }
    else {
        throw("CellNetwork::NeuronVertex::setSections(): given section list empty. neuron vertices must contain information about at least one section.");
    }
}
*/

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::translate(Vec3<R> const &d)
{
    for (auto &c : this->sections) {
        c = CellSection(c.getCompartmentId(), c.getPosition() += d, c.getRadius());
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::scale(R const &x)
{
    for (auto &c : this->sections) {
        c = CellSection(c.getCompartmentId(), c.getPosition(), c.getRadius() * x);
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::isBranchingVertex() const
{
    return (this->outdeg() > 1);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::isTerminalVertex() const
{
    return (this->outdeg() == 0);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::isSimpleVertex() const
{
    return (this->outdeg() == 1 && this->indeg() == 1);
}

/* override Graph::Vertex methods for retrieving neighbours: in a CellNetwork, all (relevant) neighours are
 * Neuron{Vertex,Edge} objects */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::list<typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getOutNeighbours() const
{
    std::list<NeuronVertex *>       out_nbs;
    std::list<NeuronVertex const *> out_nbs_const;

    /* use version for NeuronVertex * to retrieve list of pointers to non-const NeuronVertex objects */
    this->getOutNeighbours(out_nbs);

    /* convert to list of pointers to const NeuronVertex and return */
    for (auto &nb : out_nbs) {
        out_nbs_const.push_back(nb);
    }
    return out_nbs_const;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getOutNeighbours(std::list<NeuronVertex *> &out_nbs) const
{
    /* clear return list */
    out_nbs.clear();

    std::list<typename Graph<Tn, Tv, Te>::Vertex *>     out_nbs_upcast;
    NeuronVertex                                       *nb_downcast;

    /* get neighbours "upcast" to Graph::Vetex * */
    Graph<Tn, Tv, Te>::Vertex::getOutNeighbours(out_nbs_upcast);

    for (auto &nb : out_nbs_upcast) {
        /* down-cast to neuron vertex and append to result list. */
        nb_downcast = dynamic_cast<NeuronVertex *>(nb);
        if (nb_downcast) {
            out_nbs.push_back(nb_downcast);
        }
        else {
            throw("CellNetwork::NeuronVertex::getOutNeighbours() (reimplemented from Graph::Vertex): failed to "\
                "down-cast neighbour to type NeuronVertex. all vertex types must derive from NeuronVertex. do "
                "not circumvent intended typing system by upcasting to Graph and inserting Vertex objects that don't "\
                "derive from NeuronVertex.");
        }
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getOutNeighbours(std::list<NeuronVertex const *> &out_nbs) const
{
    out_nbs = this->getOutNeighbours();
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::list<typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronVertex const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getInNeighbours() const
{
    std::list<NeuronVertex *>       in_nbs;
    std::list<NeuronVertex const *> in_nbs_const;

    /* use version for NeuronVertex * to retrieve list of pointers to non-const NeuronVertex objects */
    this->getInNeighbours(in_nbs);

    /* convert to list of pointers to const NeuronVertex and return */
    for (auto &nb : in_nbs) {
        in_nbs_const.push_back(nb);
    }
    return in_nbs_const;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getInNeighbours(std::list<NeuronVertex *> &in_nbs) const
{
    /* clear return list */
    in_nbs.clear();

    std::list<typename Graph<Tn, Tv, Te>::Vertex *>     in_nbs_upcast;
    NeuronVertex                                       *nb_downcast;

    /* get neighbours "upcast" to Graph::Vetex * */
    Graph<Tn, Tv, Te>::Vertex::getInNeighbours(in_nbs_upcast);

    for (auto &nb : in_nbs_upcast) {
        /* down-cast to neuron vertex and append to result list. */
        nb_downcast = dynamic_cast<NeuronVertex *>(nb);
        if (nb_downcast) {
            in_nbs.push_back(nb_downcast);
        }
        else {
            throw("CellNetwork::NeuronVertex::getInNeighbours() (reimplemented from Graph::Vertex): failed to "\
                "down-cast neighbour to type NeuronVertex. all vertex types must derive from NeuronVertex. do "
                "not circumvent intended typing system by upcasting to Graph and inserting Vertex objects that don't "\
                "derive from NeuronVertex.");
        }
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getInNeighbours(std::list<NeuronVertex const *> &in_nbs) const
{
    in_nbs = this->getInNeighbours();
}

/* override Graph::Vertex methods for retrieving neighbours: in a CellNetwork, all (relevant) neighours are
 * Neuron{Vertex,Edge} objects */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::list<typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getOutEdges() const
{
    std::list<NeuronEdge *>         out_edges;
    std::list<NeuronEdge const *>   out_edges_const;

    /* use version for NeuronEdge * to retrieve list of pointers to non-const NeuronEdge objects */
    this->getOutEdges(out_edges);

    /* convert to list of pointers to const NeuronVertex and return */
    for (auto &e : out_edges) {
        out_edges_const.push_back(e);
    }
    return out_edges_const;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getOutEdges(std::list<NeuronEdge *> &out_edges) const
{
    /* clear return list */
    out_edges.clear();

    std::list<typename Graph<Tn, Tv, Te>::Edge *>   out_edges_upcast;
    NeuronEdge                                     *e_downcast;

    /* get neighbours "upcast" to Graph::Edge * */
    Graph<Tn, Tv, Te>::Vertex::getOutEdges(out_edges_upcast);

    for (auto &e : out_edges_upcast) {
        /* down-cast to neuron vertex and append to result list. */
        e_downcast = dynamic_cast<NeuronEdge *>(e);
        if (e_downcast) {
            out_edges.push_back(e_downcast);
        }
        else {
            throw("CellNetwork::NeuronVertex::getOutEdges() (reimplemented from Graph::Vertex): failed to "\
                "down-cast edge type to NeuronEdge. all edge types must derive from NeuronEdge. do "
                "not circumvent intended typing system by upcasting to Graph and inserting Edge objects that don't "\
                "derive from NeuronEdge.");
        }
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getOutEdges(std::list<NeuronEdge const *> &out_edges) const
{
    out_edges = this->getOutEdges();
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::list<typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getInEdges() const
{
    std::list<NeuronEdge *>         in_edges;
    std::list<NeuronEdge const *>   in_edges_const;

    /* use version for NeuronVertex * to retrieve list of pointers to non-const NeuronVertex objects */
    this->getInEdges(in_edges);

    /* convert to list of pointers to const NeuronVertex and return */
    for (auto &e : in_edges) {
        in_edges_const.push_back(e);
    }
    return in_edges_const;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getInEdges(std::list<NeuronEdge *> &in_edges) const
{
    /* clear return list */
    in_edges.clear();

    std::list<typename Graph<Tn, Tv, Te>::Edge *>       in_edges_upcast;
    NeuronEdge                                         *e_downcast;

    /* get neighbours "upcast" to Graph::Edge * */
    Graph<Tn, Tv, Te>::Vertex::getInEdges(in_edges_upcast);

    for (auto &e : in_edges_upcast) {
        /* down-cast to neuron vertex and append to result list. */
        e_downcast = dynamic_cast<NeuronEdge *>(e);
        if (e_downcast) {
            in_edges.push_back(e_downcast);
        }
        else {
            throw("CellNetwork::NeuronVertex::getOutEdges() (reimplemented from Graph::Vertex): failed to "\
                "down-cast edge type to NeuronEdge. all edge types must derive from NeuronEdge. do "
                "not circumvent intended typing system by upcasting to Graph and inserting Edge objects that don't "\
                "derive from NeuronEdge.");
        }
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getInEdges(std::list<NeuronEdge const *> &in_edges) const
{
    in_edges = this->getInEdges();
}

/* methods to access neighbours of a certain type only. safe-guards down-casting */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename NbType>
std::list<NbType const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredOutNeighbours() const
{
    std::list<NbType const *>   ret;

    /* iterate through neighbour list containing const NeuronVertex * pointers, try to down-cast to Nbtype, and append
     * to result list if dynamic cast is successful */
    NbType const *dp;
    for (auto &nb : this->getOutNeighbours() ) {
        dp = dynamic_cast<NbType const *>(nb);
        if (dp) {
            ret.push_back(dp);
        }
    }
    return ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename NbType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredOutNeighbours(std::list<NbType const *> &out_nbs) const
{
    std::list<NbType const *> ret = this->getFilteredOutNeighbours<NbType>();
    out_nbs = ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename NbType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredOutNeighbours(std::list<NbType *> &out_nbs) const
{
    out_nbs.clear();

    /* iterate through neighbour list containing "upcast" NeuronVertex * pointers, try to down-cast to NbType, and
     * append to result list if dynamic cast is successful */
    std::list<NeuronVertex *> out_nbs_upcast;
    this->getOutNeighbours(out_nbs_upcast);

    NbType *dp;
    for (auto &nb : out_nbs_upcast) {
        dp = dynamic_cast<NbType *>(nb);
        if (dp) {
            out_nbs.push_back(dp);
        }
    }
}

/* analogous for in-neighbours: set out-neighbour versions */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename NbType>
std::list<NbType const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredInNeighbours() const
{
    std::list<NbType const *> ret;

    NbType const *dp;
    std::list<NeuronVertex const *> in_nbs_upcast = this->getInNeighbours();
    for (auto &nb : in_nbs_upcast) {
        dp = dynamic_cast<NbType const *>(nb);
        if (dp) {
            ret.push_back(dp);
        }
    }
    return ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename NbType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredInNeighbours(std::list<NbType const *> &in_nbs) const
{
    std::list<NbType const *> ret = this->getFilteredInNeighbours<NbType>();
    in_nbs = ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename NbType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredInNeighbours(std::list<NbType *> &in_nbs) const
{
    in_nbs.clear();

    std::list<NeuronVertex *> in_nbs_upcast;
    this->getInNeighbours(in_nbs_upcast);

    NbType *dp;
    for (auto &nb : in_nbs_upcast) {
        dp = dynamic_cast<NbType *>(nb);
        if (dp) {
            in_nbs.push_back(dp);
        }
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename NbType>
std::list<NbType const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredNeighbours() const
{
    std::list<NbType const *> nbs = this->getFilteredOutNeighbours<NbType>();
    std::list<NbType const *> tmp = this->getFilteredInNeighbours<NbType>();

    nbs.insert(nbs.end(), tmp.begin(), tmp.end());
    return nbs;
}

/* methods to access edges of a certain type only. safe-guards down-casting */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename EdgeType>
std::list<EdgeType const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredOutEdges() const
{
    std::list<EdgeType const *>   ret;

    /* iterate through neighbour list containing const NeuronVertex * pointers, try to down-cast to Nbtype, and append
     * to result list if dynamic cast is successful */
    EdgeType const *dp;
    for (auto &nb : this->getOutEdges() ) {
        dp = dynamic_cast<EdgeType const *>(nb);
        if (dp) {
            ret.push_back(dp);
        }
    }
    return ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename EdgeType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredOutEdges(std::list<EdgeType const *> &out_edges) const
{
    std::list<EdgeType const *> ret = this->getFilteredOutEdges<EdgeType>();
    out_edges = ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename EdgeType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredOutEdges(std::list<EdgeType *> &out_edges) const
{
    out_edges.clear();

    /* iterate through neighbour list containing "upcast" NeuronVertex * pointers, try to down-cast to EdgeType, and
     * append to result list if dynamic cast is successful */
    std::list<NeuronEdge *>   out_edges_upcast;
    this->getOutEdges(out_edges_upcast);

    EdgeType *dp;
    for (auto &e : out_edges_upcast) {
        dp = dynamic_cast<EdgeType *>(e);
        if (dp) {
            out_edges.push_back(dp);
        }
    }
}

/* analogous for in-neighbours: set out-neighbour versions */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename EdgeType>
std::list<EdgeType const *>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredInEdges() const
{
    std::list<EdgeType const *>   ret;

    EdgeType const *dp;
    for (auto &e : this->getInEdges() ) {
        dp = dynamic_cast<EdgeType const *>(e);
        if (dp) {
            ret.push_back(dp);
        }
    }
    return ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename EdgeType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredInEdges(std::list<EdgeType const *> &in_edges) const
{
    std::list<EdgeType const *> ret = this->getFilteredInEdges<EdgeType>();
    in_edges = ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename EdgeType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::getFilteredInEdges(std::list<EdgeType *> &in_edges) const
{
    in_edges.clear();

    std::list<NeuronEdge *> in_edges_upcast;
    this->getInEdges(in_edges_upcast);

    EdgeType *dp;
    for (auto &e : in_edges_upcast) {
        dp = dynamic_cast<EdgeType *>(e);
        if (dp) {
            in_edges.push_back(dp);
        }
    }
}

/* iterator and const_iterator getters */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::iterator()
{
    return neuron_iterator(this->network, this->g_vit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex::iterator() const
{
    return neuron_const_iterator(this->network, this->g_vit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of SomaVertex vertex class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::SomaVertex(
    CellNetwork * const                                    &network,
    Graph<uint32_t, CellSection, Common::UnitType> const   &section_graph,
    Tv const                                               &vertex_data,
    Tso const                                              &soma_data)
        /* pass dummy section to NeuronVertex constructor, which required at least one section. ugly, but alternative is
         * to reimplement the empty() check in all leafs of the hierarchy. */
        : NeuronVertex(network, { CellSection( 0, Vec3<R>(0, 0, 0), 0) }, vertex_data)
{
    this->section_graph = section_graph;
    this->soma_data     = soma_data;

    /* clear sections, copy sections from section_graph */
    this->sections.clear();
    if (section_graph.vertices.size() > 0) {
        for (auto &v : this->section_graph.vertices) {
            this->sections.push_back(*v);
        }
    }
    else {
        throw("SomaVertex::SomaVertex(): given section graph must not be empty, i.e. must contain at least one section for the soma.");
    }
}

/* private copy ctor, use with care! see NeuronVertex and Graph::Vertex copy ctor for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::SomaVertex(SomaVertex const &s) : NeuronVertex(s)
{
}

/* private assignemt operator, use with care! see NeuronVertex and Graph::Vertex assignment operators for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::operator=(SomaVertex const &s)
{
    NeuronVertex::operator=(s);
    return (*this);
}

/* private new and delete operators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::operator new(size_t size)
{
    if (size != sizeof(CellNetwork::SomaVertex)) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "CellNetwork::SomaVertex::operator new(): size does not match sizeof(SomaVertex). internal logic error.");
    }
    return ::operator new(size);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::operator delete(void *p)
{
    ::operator delete(p);
}

/* polymorph getType() method */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::getType() const
{
    return SOMA_VERTEX;
}

/* clone method invokes copy ctor */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::clone() const
{
    return (new SomaVertex(*this));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
Vec3<R>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::getSinglePointPosition() const
{
    /* for now: return centroid of all section positions */ 
    Vec3<R> c = Aux::VecMat::nullvec<R>();

    for (auto &section : this->sections) {
        c += section.position();
    }

    c /= (R)(this->sections.size());

    return c;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::getSinglePointRadius() const
{
    /* return max radius over all sections */
    R rmax = -Aux::Numbers::inf<R>();
    for (auto &section : this->sections) {
        rmax = std::max(rmax, section.radius());
    }
    return rmax;
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::iterator()
{
    return soma_iterator(this->network, this->g_vit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaVertex::iterator() const
{
    return soma_const_iterator(this->network, this->g_vit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of abstract NeuriteVertex vertex class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::NeuriteVertex(
    CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R> * const &network,
    CellSection const                  &section,
    Tnv const                          &neurite_vertex_data,
    Tv const                           &data)
        : NeuronVertex(network, { section }, data),
        neurite_data(neurite_vertex_data)
{
    this->soma.explicitlyInvalidate();
    this->neurite.explicitlyInvalidate();
}

/* private copy ctor, use with care! see NeuronVertex and Graph::Vertex copy ctor for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::NeuriteVertex(NeuriteVertex const &n) : NeuronVertex(n)
{
    this->soma      = n.soma;
    this->neurite   = n.neurite;
}

/* private assignemt operator, use with care! see NeuronVertex and Graph::Vertex assignment operators for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::operator=(NeuriteVertex const &n)
{
    NeuronVertex::operator=(n);
    this->soma      = n.soma;
    this->neurite   = n.neurite;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
Vec3<R>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getSinglePointPosition() const
{
    return (this->getPosition());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getSinglePointRadius() const
{
    return (this->getRadius());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
Vec3<R>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getPosition() const
{
    return (this->sections.begin()->position());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::setPosition(Vec3<R> const &p)
{
    this->sections.begin()->position() = p;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getRadius() const
{
    return (this->sections.begin()->radius());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::setRadius(R const &r)
{
    this->sections.begin()->radius() = r;
}

/* get soma this neurite vertex belongs to */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getSoma() const
{
    return this->soma;
}

/* get neurite_rootedge_iterator, representing the neurite, that this neurite vertex belongs to */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getNeurite() const
{
    return this->neurite;
}

/* customized predicates restricted to cell-tree edges */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::isNeuriteSimpleVertex() const
{
    return (
        (this->template getFilteredInEdges<NeuriteSegment>().size() + this->template getFilteredInEdges<NeuriteRootEdge>().size())  == 1 &&
        this->template getFilteredOutEdges<NeuriteSegment>().size() == 1);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::isNeuriteBranchingVertex() const
{
    return (this->template getFilteredOutEdges<NeuriteSegment>().size() > 1);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::isNeuriteTerminalVertex() const
{
    return (this->template getFilteredOutEdges<NeuriteSegment>().size() == 0);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::isNeuriteRootVertex() const
{
    return (this->template getFilteredInEdges<NeuriteRootEdge>().size() == 1);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getNeuriteParent() const
{
    std::list<NeuriteVertex *>  tmp;
    this->template getFilteredInNeighbours<NeuriteVertex>(tmp);
    if (tmp.size() != 1) {
        throw("CellNetwork::NeuriteVertex::getNeuriteParent(): vertex does not have exactly one in-neurite-segment.");
    }
    else {
        return (tmp.front()->iterator());
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getNeuriteFirstChild() const
{
    std::list<NeuriteVertex *>  tmp;
    this->template getFilteredOutNeighbours<NeuriteVertex>(tmp);
    if (tmp.empty()) {
        throw("CellNetwork::NeuriteVertex::getNeuriteFirstChild(): vertex does not have any \"children\", i.e. no NeuriteVertex as out-neighbour.");
    }
    else {
        return (tmp.front()->iterator());
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_segment_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::getNeuriteSegmentFirstChild() const
{
    std::list<NeuriteSegment *>  tmp;
    this->template getFilteredOutEdges<NeuriteSegment>(tmp);
    if (tmp.empty()) {
        throw("CellNetwork::NeuriteVertex::getNeuriteFirstChild(): vertex does not have any \"children\", i.e. no NeuriteVertex as out-neighbour.");
    }
    else {
        return (tmp.front()->iterator());
    }
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::iterator()
{
    return neurite_iterator(this->network, this->g_vit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertex::iterator() const
{
    return neurite_const_iterator(this->network, this->g_vit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of AxonVertex vertex class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::AxonVertex(
    CellNetwork * const                &network,
    CellSection const                  &section,
    Tax const                          &axon_data,
    Tnv const                          &neurite_vertex_data,
    Tv const                           &vertex_data)
        : NeuriteVertex(network, section, neurite_vertex_data, vertex_data),
        axon_data(axon_data)
{
}

/* private copy ctor, use with care! see NeuronVertex and Graph::Vertex copy ctor for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::AxonVertex(AxonVertex const &a) : NeuriteVertex(a)
{
    this->axon_data = a.axon_data;
}

/* private assignemt operator, use with care! see NeuronVertex and Graph::Vertex assignment operators for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::operator=(AxonVertex const &a)
{
    NeuriteVertex::operator=(a);
    this->axon_data = a.axon_data;
    
    return (*this);
}

/* private new and delete operators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::operator new(size_t size)
{
    if (size != sizeof(CellNetwork::AxonVertex)) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "CellNetwork::AxonVertex::operator new(): size does not match sizeof(AxonVertex). internal logic error.");
    }
    return ::operator new(size);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::operator delete(void *p)
{
    ::operator delete(p);
}

/* polymorph getType() method */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::getType() const
{
    return AXON_VERTEX;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::clone() const
{
    return (new AxonVertex(*this));
}

/* get axon_rootedge_iterator, representing the axon, that this axon vertex belongs to */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::getAxon() const
{
    AxonRootEdge *ar = dynamic_cast<AxonRootEdge *>(&(*this->neurite));
    if (ar) {
        return ar->iterator();
    }
    else {
        throw("AxonVertex::getAxon(): failed to down-cast NeuriteRootEdge to AxonRootEdge. internal logic error.");
    }
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::iterator()
{
    return axon_iterator(this->network, this->g_vit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonVertex::iterator() const
{
    return axon_const_iterator(this->network, this->g_vit);
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of DendriteVertex vertex class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::DendriteVertex(
    CellNetwork * const                &network,
    CellSection const                  &section,
    bool                                apical_dendrite,
    Tde const                          &dendrite_data,
    Tnv const                          &neurite_vertex_data,
    Tv const                           &vertex_data)
        : NeuriteVertex(network, section, neurite_vertex_data, vertex_data),
        apical_dendrite(apical_dendrite),
        dendrite_data(dendrite_data)
{
}

/* private copy ctor, use with care! see NeuronVertex and Graph::Vertex copy ctor for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::DendriteVertex(DendriteVertex const &d) : NeuriteVertex(d)
{
    this->apical_dendrite   = d.apical_dendrite;
    this->dendrite_data     = d.dendrite_data;
}

/* private assignemt operator, use with care! see NeuronVertex and Graph::Vertex assignment operators for details. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::operator=(DendriteVertex const &d)
{
    NeuriteVertex::operator=(d);
    this->apical_dendrite   = d.apical_dendrite;
    this->dendrite_data     = d.dendrite_data;

    return (*this);
}

/* private new and delete operators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::operator new(size_t size)
{
    if (size != sizeof(CellNetwork::DendriteVertex)) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "CellNetwork::DendriteVertex::operator new(): size does not match sizeof(DendriteVertex). internal logic error.");
    }
    return ::operator new(size);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::operator delete(void *p)
{
    ::operator delete(p);
}

/* polymorph getType() method */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::getType() const
{
    return DENDRITE_VERTEX;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertex *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::clone() const
{
    return (new DendriteVertex(*this));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::isApicalDendrite() const
{
    return (this->apical_dendrite);
}

/* get dendrite_rootedge_iterator, representing the dendrite, that this dendrite vertex belongs to */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::getDendrite() const
{
    DendriteRootEdge *dr = dynamic_cast<DendriteRootEdge *>(&(*this->neurite));
    if (dr) {
        return dr->iterator();
    }
    else {
        throw("DendriteVertex::getDendrite(): failed to down-cast NeuriteRootEdge to DendriteRootEdge. internal logic error.");
    }
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::iterator()
{
    return dendrite_iterator(this->network, this->g_vit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteVertex::iterator() const
{
    return dendrite_const_iterator(this->network, this->g_vit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of (neuron) edge class hierarchy for CellNetworks                 
 *                        starting with abstract class NeuronEdge
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::NeuronEdge(
    CellNetwork * const    &network,
    NeuronVertex           *v_src,
    NeuronVertex           *v_dst,
    Te const               &edge_data)
        : Graph<Tn, Tv, Te>::Edge(network, v_src, v_dst, edge_data, false)
{
    this->network       = network;
    this->v_neuron_src  = v_src;
    this->v_neuron_dst  = v_dst;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::NeuronEdge(NeuronEdge const &e) : Graph<Tn, Tv, Te>::Edge(e)
{
    this->network       = e.network;
    this->v_neuron_src  = e.v_neuron_src;
    this->v_neuron_dst  = e.v_neuron_dst;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::operator=(NeuronEdge const &e)
{
    Graph<Tn, Tv, Te>::Edge::operator=(e);
    this->network       = e.network;
    this->v_neuron_src  = e.v_neuron_src;
    this->v_neuron_dst  = e.v_neuron_dst;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::~NeuronEdge()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::updateVertexPointers()
{
    this->v_neuron_src = dynamic_cast<NeuronVertex *>(this->v_src);
    this->v_neuron_dst = dynamic_cast<NeuronVertex *>(this->v_dst);
}

/* virtual method getLength(). default implementation: return difference of getSinglePointPosition() for source and
 * destination vertex */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::getLength() const
{
    Vec3<R> psrc = this->v_neuron_src->getSinglePointPosition();
    Vec3<R> pdst = this->v_neuron_dst->getSinglePointPosition();
    return ( (psrc - pdst).len2() );
}

/* get source vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::getSourceVertex()
{
    typename Graph<Tn, Tv, Te>::vertex_iterator srcit = this->v_src->iterator();
    return neuron_iterator(this->network, Graph<Tn, Tv, Te>::getInternalIterator(srcit));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::getSourceVertex() const
{
    typename Graph<Tn, Tv, Te>::vertex_const_iterator srcit = this->v_src->iterator();
    return neuron_const_iterator(this->network, Graph<Tn, Tv, Te>::getInternalIterator(srcit));
}

/* get destination vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::getDestinationVertex()
{
    typename Graph<Tn, Tv, Te>::vertex_iterator dstit = this->v_dst->iterator();
    return neuron_iterator(this->network, Graph<Tn, Tv, Te>::getInternalIterator(dstit));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::getDestinationVertex() const
{
    typename Graph<Tn, Tv, Te>::vertex_const_iterator dstit = this->v_dst->iterator();
    return neuron_const_iterator(this->network, Graph<Tn, Tv, Te>::getInternalIterator(dstit));
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_edge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::iterator()
{
    return neuron_edge_iterator(this->network, this->g_eit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_edge_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdge::iterator() const
{
    return neuron_edge_const_iterator(this->network, this->g_eit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of (abstract) NeuriteRootEdge edge class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::NeuriteRootEdge(
    CellNetwork * const    &network,
    SomaVertex             *v_src,
    NeuriteVertex          *v_dst,
    Tnr const              &neurite_root_edge_data,
    Te const               &edge_data)
        : NeuronEdge(network, v_src, v_dst, edge_data), v_soma_src(v_src), v_neurite_dst(v_dst),
        neurite_root_edge_data(neurite_root_edge_data)
{
    this->neurite_root_edge_data = neurite_root_edge_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::NeuriteRootEdge(NeuriteRootEdge const &e) : NeuronEdge(e)
{
    this->v_soma_src    = e.v_soma_src;
    this->v_neurite_dst = e.v_neurite_dst;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuriteRootEdge &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::operator=(NeuriteRootEdge const &e)
{
    NeuronEdge::operator=(e);
    this->v_soma_src    = e.v_soma_src;
    this->v_neurite_dst = e.v_neurite_dst;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::~NeuriteRootEdge()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::updateVertexPointers()
{
    NeuronEdge::updateVertexPointers();
    this->v_soma_src    = dynamic_cast<SomaVertex *>(this->v_neuron_src);
    this->v_neurite_dst = dynamic_cast<NeuriteVertex *>(this->v_neuron_dst);
}

/* get source vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::getSourceVertex()
{
    return (this->v_soma_src->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::getSourceVertex() const
{
    return (this->v_soma_src->iterator());
}

/* get destination vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::getDestinationVertex()
{
    return (this->v_neurite_dst->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::getDestinationVertex() const
{
    return (this->v_neurite_dst->iterator());
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::iterator()
{
    return neurite_rootedge_iterator(this->network, this->g_eit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_rootedge_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdge::iterator() const
{
    return neurite_rootedge_const_iterator(this->network, this->g_eit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of AxonRootEdge edge class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::AxonRootEdge(
    CellNetwork * const    &network,
    SomaVertex             *v_src,
    AxonVertex             *v_dst,
    Tar const              &axon_root_edge_data,
    Tnr const              &neurite_root_edge_data,
    Te const               &edge_data)
        : NeuriteRootEdge(network, v_src, v_dst, neurite_root_edge_data, edge_data),
        v_axon_dst(v_dst)
{
    this->axon_root_edge_data = axon_root_edge_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::AxonRootEdge(AxonRootEdge const &e) : NeuriteRootEdge(e)
{
    this->v_axon_dst            = e.v_axon_dst;
    this->axon_root_edge_data   = e.axon_root_edge_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::AxonRootEdge &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::operator=(AxonRootEdge const &e)
{
    NeuriteRootEdge::operator=(e);
    this->v_axon_dst            = e.v_axon_dst;
    this->axon_root_edge_data   = e.axon_root_edge_data;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::~AxonRootEdge()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::updateVertexPointers()
{
    NeuriteRootEdge::updateVertexPointers();
    this->v_axon_dst    = dynamic_cast<AxonVertex *>(this->v_neurite_dst);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::getType() const
{
    return AXON_ROOT_EDGE;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::clone() const
{
    return (new AxonRootEdge(*this));
}

/* get destination vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::getDestinationVertex()
{
    return (this->v_axon_dst->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::getDestinationVertex() const
{
    return (this->v_axon_dst->iterator());
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::iterator()
{
    return axon_rootedge_iterator(this->network, this->g_eit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_rootedge_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdge::iterator() const
{
    return axon_rootedge_const_iterator(this->network, this->g_eit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of DendriteRootEdge edge class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::DendriteRootEdge(
    CellNetwork * const    &network,
    SomaVertex             *v_src,
    DendriteVertex         *v_dst,
    Tdr const              &dendrite_root_edge_data,
    Tnr const              &neurite_root_edge_data,
    Te const               &edge_data)
        : NeuriteRootEdge(network, v_src, v_dst, neurite_root_edge_data, edge_data),
        v_dendrite_dst(v_dst)
{
    this->dendrite_root_edge_data = dendrite_root_edge_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::DendriteRootEdge(DendriteRootEdge const &e) : NeuriteRootEdge(e)
{
    this->v_dendrite_dst            = e.v_dendrite_dst;
    this->dendrite_root_edge_data   = e.dendrite_root_edge_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::DendriteRootEdge &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::operator=(DendriteRootEdge const &e)
{
    NeuriteRootEdge::operator=(e);
    this->v_dendrite_dst            = e.v_dendrite_dst;
    this->dendrite_root_edge_data   = e.dendrite_root_edge_data;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::~DendriteRootEdge()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::updateVertexPointers()
{
    NeuriteRootEdge::updateVertexPointers();
    this->v_dendrite_dst = dynamic_cast<DendriteVertex *>(this->v_neurite_dst);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::getType() const
{
    return DENDRITE_ROOT_EDGE;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::clone() const
{
    return (new DendriteRootEdge(*this));
}

/* get destination vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::getDestinationVertex()
{
    return (this->v_dendrite_dst->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::getDestinationVertex() const
{
    return (this->v_dendrite_dst->iterator());
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::iterator()
{
    return dendrite_rootedge_iterator(this->network, this->g_eit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_rootedge_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdge::iterator() const
{
    return dendrite_rootedge_const_iterator(this->network, this->g_eit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of (abstract) NeuriteSegment edge class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::NeuriteSegment(
    CellNetwork * const    &network,
    NeuriteVertex          *v_src,
    NeuriteVertex          *v_dst,
    Tns const              &neurite_segment_data,
    Te const               &edge_data)
        : NeuronEdge(network, v_src, v_dst, edge_data),
        v_neurite_src(v_src), v_neurite_dst(v_dst),
        neurite_segment_data(neurite_segment_data)
{
    this->soma.explicitlyInvalidate();
    this->neurite.explicitlyInvalidate();
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::NeuriteSegment(NeuriteSegment const &e) : NeuronEdge(e)
{
    this->soma                  = e.soma;
    this->neurite               = e.neurite;
    this->v_neurite_src         = e.v_neurite_src;
    this->v_neurite_dst         = e.v_neurite_dst;
    this->neurite_segment_data  = e.neurite_segment_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuriteSegment &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::operator=(NeuriteSegment const &e)
{
    NeuronEdge::operator=(e);
    this->soma                  = e.soma;
    this->neurite               = e.neurite;
    this->v_neurite_src         = e.v_neurite_src;
    this->v_neurite_dst         = e.v_neurite_dst;
    this->neurite_segment_data  = e.neurite_segment_data;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::~NeuriteSegment()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::updateVertexPointers()
{
    NeuronEdge::updateVertexPointers();
    this->v_neurite_src = dynamic_cast<NeuriteVertex *>(this->v_neuron_src);
    this->v_neurite_dst = dynamic_cast<NeuriteVertex *>(this->v_neuron_dst);
}

/* getters for radius difference / ratio and length of neurite segment */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getAbsoluteRadiusDifference() const
{
    R rsrc = this->v_neuron_src->getSinglePointRadius();
    R rdst = this->v_neuron_dst->getSinglePointRadius();
    return ( std::abs(rsrc - rdst) );
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getRadiusRatio() const
{
    R rsrc = this->v_neuron_src->getSinglePointRadius();
    R rdst = this->v_neuron_dst->getSinglePointRadius();
    return ( std::max( rsrc / rdst, rdst / rsrc ) );
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getMinRadius() const
{
    return (std::min(this->v_neurite_src->getRadius(), this->v_neurite_dst->getRadius()));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getMaxRadius() const
{
    return (std::max(this->v_neurite_src->getRadius(), this->v_neurite_dst->getRadius()));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getAngle() const
{
    /* get "parent" of source vertex, might be a soma, so use neuron_iterator */
    auto plist = this->getSourceVertex()->getInNeighbours();
    if (plist.size() == 1) {
        NeuronVertex const &p   = *(plist.front());
        NeuriteVertex const &u  = *(this->getSourceVertex());
        NeuriteVertex const &v  = *(this->getDestinationVertex());

        Vec3<R> d_pu            = u.getPosition() - p.getSinglePointPosition();
        Vec3<R> d_uv            = v.getPosition() - u.getPosition();

        d_pu.normalize();
        d_uv.normalize();

        R  scalprod = d_pu * d_uv;
        R  angle;
        if (scalprod > 1.0 - 1E-13) {
            angle = 0;
        }
        else {
            if (std::is_same<R, double>::value) {
                angle = acos(scalprod);
            }
            else if (std::is_same<R, float>::value) {
                angle = acosf(scalprod);
            }
            else if (std::is_same<R, long double>::value) {
                angle = acosl(scalprod);
            }
            else {
                throw("NeuriteSegment::getAngle(): angle computation not supported for template type R.");
            }
        }
        return angle;
    }
    else {
        throw("NeuriteSegment::getAxon(): can't get angle if source vertex of neurite segment does not have exactly "\
            "one in-neighbour, i.e. a unique parent.");
    }
}

/* used in the check for secondary minimum distance violation: for a neurite segment (u, v), let rmax_u_nb be
*
*  1. the maximum radius of all vertices adjacent to u (including u itself) EXCEPT for v if u is NOT a neurite root
*  vertex.
*
*  2. 0, if u is a neurite root vertex
*
*  and let rmax_v_nb be 
*
*  1. the maximum radius of all vertices adjacent to v (including v itself) EXCEPT for u if v is NOT a neurite terminal
*  vertex
*
*  2. 0, if v is a neurite terminal vertex.
*
* a secondary minimum distance violation occurs if len(u,v) is less than rmax_u_nb + rmax_v_nb, because then two
* neighbouring canal segments will "trivially" intersect by the definition of volume used for intersection analysis. */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<R, R>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getSMDVRadii() const
{
    neurite_const_iterator  u = this->getSourceVertex();
    neurite_const_iterator  v = this->getDestinationVertex();
    neurite_const_iterator  p;
    R                       rmax_u_nb, rmax_v_nb;

    if (!u->isNeuriteRootVertex()) {
        p           = u->getNeuriteParent();
        rmax_u_nb   = std::max(u->getRadius(), p->getRadius());

        for (auto &nb : u->template getFilteredOutNeighbours<NeuriteVertex>()) {
            /* skip v */
            if (nb->iterator() != v) {
                rmax_u_nb = std::max( rmax_u_nb, nb->getRadius() );
            }
        }
    }
    else {
        rmax_u_nb = 0.0;
    }

    /* get rmax_v_nb: v is the child of u, skip parent u */
    if (!v->isNeuriteTerminalVertex()) {
        rmax_v_nb = v->getRadius();
        for (auto &nb : u->template getFilteredOutNeighbours<NeuriteVertex>()) {
            rmax_v_nb = std::max( rmax_v_nb, nb->getRadius() );
        }
    }
    else {
        rmax_v_nb = 0.0;
    }

    return std::pair<R, R>(rmax_u_nb, rmax_v_nb);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getLength() const
{
    Vec3<R> psrc = this->v_neuron_src->getSinglePointPosition();
    Vec3<R> pdst = this->v_neuron_dst->getSinglePointPosition();
    return ( (psrc - pdst).len2() );
}

/* get source vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getSourceVertex()
{
    return (this->v_neurite_src->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getSourceVertex() const
{
    return (this->v_neurite_src->iterator());
}

/* get destination vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getDestinationVertex()
{
    return (this->v_neurite_dst->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getDestinationVertex() const
{
    return (this->v_neurite_dst->iterator());
}

/* get soma / neurite */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getSoma() const
{
    return (this->soma);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::getNeurite() const
{
    return (this->neurite);
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_segment_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::iterator()
{
    return neurite_segment_iterator(this->network, this->g_eit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_segment_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegment::iterator() const
{
    return neurite_segment_const_iterator(this->network, this->g_eit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of AxonSegment edge class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::AxonSegment(
    CellNetwork * const    &network,
    AxonVertex             *v_src,
    AxonVertex             *v_dst,
    Tas const              &axon_segment_data,
    Tns const              &neurite_segment_data,
    Te const               &edge_data)
        : NeuriteSegment(network, v_src, v_dst, neurite_segment_data, edge_data), 
        v_axon_src(v_src), v_axon_dst(v_dst),
        axon_segment_data(axon_segment_data)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::AxonSegment(AxonSegment const &e) : NeuriteSegment(e)
{
    this->v_axon_src        = e.v_axon_src;
    this->v_axon_dst        = e.v_axon_dst;
    this->axon_segment_data = e.axon_segment_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::AxonSegment &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::operator=(AxonSegment const &e)
{
    NeuriteSegment::operator=(e);
    this->v_axon_src        = e.v_axon_src;
    this->v_axon_dst        = e.v_axon_dst;
    this->axon_segment_data = e.axon_segment_data;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::~AxonSegment()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::updateVertexPointers()
{
    NeuriteSegment::updateVertexPointers();
    this->v_axon_src = dynamic_cast<AxonVertex *>(this->v_neurite_src);
    this->v_axon_dst = dynamic_cast<AxonVertex *>(this->v_neurite_dst);
}

/* get type */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::getType() const
{
    return AXON_SEGMENT;
}

/* clone() method */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::clone() const
{
    return (new AxonSegment(*this));
}

/* get source vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::getSourceVertex()
{
    return (this->v_axon_src->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::getSourceVertex() const
{
    return (this->v_axon_src->iterator());
}

/* get destination vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::getDestinationVertex()
{
    return (this->v_axon_dst->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::getDestinationVertex() const
{
    return (this->v_axon_dst->iterator());
}

/* iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_segment_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::iterator()
{
    return axon_segment_iterator(this->network, this->g_eit);
}

/* get axon_rootedge_iterator, representing the axon, that this axon segment belongs to */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::getAxon() const
{
    AxonRootEdge *ar = dynamic_cast<AxonRootEdge *>(&(*this->neurite));
    if (ar) {
        return ar->iterator();
    }
    else {
        throw("AxonVertex::getAxon(): failed to down-cast NeuriteRootEdge to AxonRootEdge. internal logic error.");
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_segment_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegment::iterator() const
{
    return axon_segment_const_iterator(this->network, this->g_eit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of DendriteSegment edge class
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::DendriteSegment(
    CellNetwork * const    &network,
    DendriteVertex         *v_src,
    DendriteVertex         *v_dst,
    Tds const              &dendrite_segment_data,
    Tns const              &neurite_segment_data,
    Te const               &edge_data)
        : NeuriteSegment(network, v_src, v_dst, neurite_segment_data, edge_data),
        v_dendrite_src(v_src), v_dendrite_dst(v_dst),
        dendrite_segment_data(dendrite_segment_data)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::DendriteSegment(DendriteSegment const &e) : NeuriteSegment(e)
{
    this->v_dendrite_src        = e.v_dendrite_src;
    this->v_dendrite_dst        = e.v_dendrite_dst;
    this->dendrite_segment_data = e.dendrite_segment_data;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::DendriteSegment &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::operator=(DendriteSegment const &e)
{
    NeuriteSegment::operator=(e);
    this->v_dendrite_src        = e.v_dendrite_src;
    this->v_dendrite_dst        = e.v_dendrite_dst;
    this->dendrite_segment_data = e.dendrite_segment_data;

    return (*this);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::~DendriteSegment()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::updateVertexPointers()
{
    NeuriteSegment::updateVertexPointers();
    this->v_dendrite_src = dynamic_cast<DendriteVertex *>(this->v_neurite_src);
    this->v_dendrite_dst = dynamic_cast<DendriteVertex *>(this->v_neurite_dst);
}

/* get type */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::getType() const
{
    return DENDRITE_SEGMENT;
}

/* clone() method */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::NeuronEdge *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::clone() const
{
    return (new DendriteSegment(*this));
}

/* get source vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::getSourceVertex()
{
    return (this->v_dendrite_src->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::getSourceVertex() const
{
    return (this->v_dendrite_src->iterator());
}

/* get destination vertex, both const and non-const version */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::getDestinationVertex()
{
    return (this->v_dendrite_dst->iterator());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::getDestinationVertex() const
{
    return (this->v_dendrite_dst->iterator());
}

/* get dendrite_rootedge_iterator, representing the dendrite, that this dendrite segment belongs to */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_rootedge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::getDendrite() const
{
    DendriteRootEdge *dr = dynamic_cast<DendriteRootEdge *>(&(*this->neurite));
    if (dr) {
        return dr->iterator();
    }
    else {
        throw("DendriteVertex::getDendrite(): failed to down-cast NeuriteRootEdge to DendriteRootEdge. internal logic error.");
    }
}

/* get const / non-const dendrite_segment_iterators */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_segment_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::iterator()
{
    return dendrite_segment_iterator(this->network, this->g_eit);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_segment_const_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegment::iterator() const
{
    return dendrite_segment_const_iterator(this->network, this->g_eit);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of abstract generic CellNetworkAccessor class                                
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
CellNetworkAccessor(
    CellNetwork                    &C_arg,
    InternalType                   &int_ds_arg,
    ValuePred const                &p) : C(C_arg), int_ds(int_ds_arg), pred(p)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
IteratorType
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
begin()
{
    /* start iterating over internal data structure, cast to delegator type, check dynamic cast to Valuetype until
     * either an object of matching type has been found or end() has been reached. */
    typename InternalType::iterator     it;
    DelegatorType                      *dp;
    for (it = this->int_ds.begin(); it != this->int_ds.end(); ++it) {
        dp = dynamic_cast<DelegatorType *>(BaseType::getPtr(it));

        if (dp) {
            /* if type matches, return wrapped iterator, otherwise return end(), since the search has been negative,
             * just as if the id had not been present. */
            if (this->pred(*dp)) {
                return IteratorType( &(this->C), it);
            }
        }
        else {
            throw("CellNetwork::CellNetworkAccessor::begin(): failed to down-cast base type to delegator type. internal logic error.");
        }
    }

    /* no objects of matching type found. return iterator to end() */
    return IteratorType( &(this->C), this->int_ds.end());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
ConstIteratorType
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
begin() const
{
    /* start iterating over internal data structure, cast to delegator type, check dynamic cast to Valuetype until
     * either an object of matching type has been found or end() has been reached. */
    typename InternalType::const_iterator   it;
    DelegatorType const                    *dp;
    for (it = this->int_ds.begin(); it != this->int_ds.end(); ++it) {
        dp = dynamic_cast<DelegatorType const *>(BaseType::getPtr(it));

        if (dp) {
            /* if type matches, return wrapped iterator, otherwise return end(), since the search has been negative,
             * just as if the id had not been present. */
            if (this->pred(*dp)) {
                return ConstIteratorType( &(this->C), it);
            }
        }
        else {
            throw("CellNetwork::CellNetworkAccessor::begin(): failed to down-cast base type to delegator type. internal logic error.");
        }
    }

    /* no objects of matching type found. return iterator to end() */
    return ConstIteratorType( &(this->C), this->int_ds.end());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
IteratorType
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
end()
{
    /* for any specialized type (vertex or edge), the end() iterator is simply a wrapped version of the vertex / edge
     * end() operator, which in turn wraps X::end(). hence this is analogous to Graph::{Vertex, Edge}Accessor::end(),
     * but the CellNetworkIterator constructor takes the internal data structure int_ds as an additional argument, since 
     * it needs to "skip" over objects of non-matching type */
    return IteratorType( &(this->C), this->int_ds.end());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
ConstIteratorType
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
end() const
{
    /* see non-const end() version */
    return ConstIteratorType( &(this->C), this->int_ds.end());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
IteratorType
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
find(uint32_t id)
{
    /* locate id in internal data structure using InternalType::find(). note that this and InternalType::at() (including
     * their type signature of course) are requirements on the template type InternalType. std::map is used internally,
     * but other choices are possible. */
    typename InternalType::iterator searchit_int = this->int_ds.find(id);

    /* id internally found, check if object has matching type using ValuePred, just as in CellNetworkIterator */
    if (searchit_int != this->int_ds.end()) {
        /* down-cast searchit_int to delegator type, throw if downcast fails. */
        DelegatorType *dp = dynamic_cast<DelegatorType *>(BaseType::getPtr(searchit_int));
        if (dp) {
            /* if type matches, return wrapped iterator, otherwise return end(), since the search has been negative,
             * just as if the id had not been present. */
            if (this->pred(*dp)) {
                return IteratorType(&(this->C), searchit_int);
            }
            else {
                return (this->end());
            }
        }
        else {
            throw("CellNetwork::CellNetworkAccessor::find(): failed to down-cast base type to delegator type. internal logic error.");
        }
    }
    /* not found in internal data structure, return end() to indicate to the caller */
    else {
        return (this->end());
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
ConstIteratorType
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
find(uint32_t id) const
{
    /* locate id in internal data structure using InternalType::find(). note that this and InternalType::at() (including
     * their type signature of course) are requirements on the template type InternalType. std::map is used internally,
     * but other choices are possible. */
    typename InternalType::const_iterator searchit_int = this->int_ds.find(id);

    /* id internally found, check if object has matching type using ValuePred, just as in CellNetworkIterator */
    if (searchit_int != this->int_ds.end()) {
        /* down-cast searchit_int to delegator type, throw if downcast fails. */
        DelegatorType const *dp = dynamic_cast<DelegatorType const *>(BaseType::getPtr(searchit_int));
        if (dp) {
            /* if type matches, return wrapped iterator, otherwise return end(), since the search has been negative,
             * just as if the id had not been present. */
            if (this->pred(*dp)) {
                return ConstIteratorType(&(this->C), searchit_int);
            }
            else {
                return (this->end());
            }
        }
        else {
            throw("CellNetwork::CellNetworkAccessor::find(): failed to down-cast base type to delegator type. internal logic error.");
        }
    }
    /* not found in internal data structure, return end() to indicate to the caller */
    else {
        return (this->end());
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
ValueType &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
at(uint32_t id)
{
    /* see find() vs. find() const above, analogous procedure here. */

    /* create const reference to (this) in order to access const versions of find() and end() */
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType> const     &const_this = (*this);

    /* use const version of find() to retrieve result */
    ConstIteratorType constit = const_this.at(id);

    /* if object of matching id AND matching type has been found, return non-const reference */

    /* compare returned const iterator to const end() */
    if (constit != const_this.end()) {
        /* return non-const reference by using internals of ConstIteratorType. alternative: retrieve
         * value and cast away constness, but this is even uglier and more dangerous than this one here.
         * however, the method below assumes that IteratorType / ConstIteratorType are of type CellNetworkIterator
         * (or derived from the latter). */
        return *(ValueType::getPtr(constit.int_it));
    } 
    /* of nothing has been found, throw out-of-range exception */
    else {
        throw("CellNetworkAccessor<..>::at(uint32_t id): object of matching id and type not found. out of range.");
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
ValueType const &
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
at(uint32_t id) const
{
    /* use const-version of find() to retrieve object. if it has not been found, throw an exception */
    ConstIteratorType it = this->find(id);

    /* if object has been found, return const reference to object by dereferencing iterator */
    if (it != this->end()) {
        return (*it);
    }
    /* of nothing has been found, throw out-of-range exception */
    else {
        throw("CellNetworkAccessor<..>::at(uint32_t id) const: object of matching id and type not found. out of range.");
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
exists(uint32_t id) const
{
    /* use const-version of find() to check for existence */
    return (this->find(id) != this->end());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <
    typename ValueType, typename BaseType, typename DelegatorType, typename InternalType, typename ValuePred,
    typename IteratorType, typename  ConstIteratorType,typename DataType
>
size_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
    CellNetworkAccessor<
        ValueType,
        BaseType,
        DelegatorType,
        InternalType,
        ValuePred,
        IteratorType,
        ConstIteratorType,
        DataType
    >::
size() const
{
    /* use internal iterator architecture to compute size */
    ConstIteratorType it = this->begin();
    size_t size = 0;
    while (it != this->end()) {
        size++;
        ++it;
    }
    return size;
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of NeuronVertexAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertexAccessor::NeuronVertexAccessor(CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R> &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        NeuronVertex,
        typename Graph<Tn, Tv, Te>::Vertex,
        NeuronVertex,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
        NeuronVertexPred,
        neuron_iterator,
        neuron_const_iterator,
        Tv
    >(C, C.V, C.isNeuronVertex)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertexAccessor::erase(neuron_iterator it)
{
    /* "up"-convert neuron_iterator to Graph::vertex_iterator by passing pointer to CellNetwork, which is upcast to
     * Graph, and internal iterator */
    typename Graph<Tn, Tv, Te>::vertex_iterator vit(&(this->C), it.int_it);

    try {
        vit = this->C.vertices.erase(vit);
    }
    catch (GraphEx& ex) {
        /* FIXME: exception handling, wrap into CellNetwork exception class */
        throw;
    }
    /* convert back "down" by wrapping internal iterator in a  neuron_iterator. */
    return neuron_iterator(&(this->C), Graph<Tn, Tv, Te>::getInternalIterator(vit));
}

/* FIXME: implement me */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronVertexAccessor::erase(uint32_t id)
{
    return true;
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of NeuriteVertexAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertexAccessor::NeuriteVertexAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        NeuriteVertex,
        typename Graph<Tn, Tv, Te>::Vertex,
        NeuronVertex,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
        NeuriteVertexPred,
        neurite_iterator,
        neurite_const_iterator,
        Tv
    >(C, C.V, C.isNeuriteVertex)
{
}

/*
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertexAccessor::erase(neurite_iterator it)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteVertexAccessor::erase(uint32_t id)
{
}
*/

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of SomaAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>
::SomaAccessor::SomaAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        SomaVertex,
        typename Graph<Tn, Tv, Te>::Vertex,
        NeuronVertex,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
        SomaVertexPred,
        soma_iterator,
        soma_const_iterator,
        Tv
    >(C, C.V, C.isSomaVertex)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaAccessor::insert(
    Graph<uint32_t, CellSection, Common::UnitType> const   &section_graph,
    Tv const                                               &vertex_data,
    Tso const                                              &soma_data)
{
    SomaVertex *s   = new SomaVertex(&(this->C), section_graph, vertex_data, soma_data);
    auto rpair      = this->C.protectedVertexInsert(s);

    /* if vertex has been inserted successfully, convert returned Graph::vertex_iterator to soma_iterator using the
     * returned internal map iterator referring to Graph::V (this is not a Graph::vertex_iterator, which does not allow
     * access to internals). */
    if (rpair.second) {
        return soma_iterator( &(this->C), rpair.first);
    }
    else {
        throw("CellNetwork::SomaAccessor::insert(): failed to insert soma vertex.");
    }
}

/*
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::soma_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaAccessor::erase(soma_iterator it)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
SomaAccessor::erase(uint32_t id)
{
}
*/

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of AxonAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonAccessor::AxonAccessor(CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R> &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        AxonVertex,
        typename Graph<Tn, Tv, Te>::Vertex,
        NeuronVertex,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
        AxonVertexPred,
        axon_iterator,
        axon_const_iterator,
        Tv
    >(C, C.V, C.isAxonVertex)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonAccessor::insert(
    CellSection const      &section,
    Tax const              &axon_data,
    Tnv const              &neurite_vertex_data,
    Tv const               &vertex_data)
{
    AxonVertex *a   = new AxonVertex(&(this->C), section, axon_data, neurite_vertex_data, vertex_data);
    auto rpair      = this->C.protectedVertexInsert(a);

    /* if vertex has been inserted successfully, convert returned Graph::vertex_iterator to axon_iterator using the
     * returned internal map iterator referring to Graph::V (this is not a Graph::vertex_iterator, which does not allow
     * access to internals). */
    if (rpair.second) {
        return axon_iterator( &(this->C), rpair.first );
    }
    else {
        throw("CellNetwork::AxonAccessor::insert(): failed to insert axon vertex.");
    }
}

/*
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonAccessor::erase(axon_iterator it)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonAccessor::erase(uint32_t id)
{
}
*/

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of DendriteAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteAccessor::DendriteAccessor(CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R> &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        DendriteVertex,
        typename Graph<Tn, Tv, Te>::Vertex,
        NeuronVertex,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
        DendriteVertexPred,
        dendrite_iterator,
        dendrite_const_iterator,
        Tv
    >(C, C.V, C.isDendriteVertex)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteAccessor::insert(
    CellSection const      &section,
    bool                    apical_dendrite,
    Tde const              &dendrite_data,
    Tnv const              &neurite_vertex_data,
    Tv const               &vertex_data)
{
    DendriteVertex *d   = new DendriteVertex(&(this->C), section, apical_dendrite, dendrite_data, neurite_vertex_data, vertex_data);
    auto rpair          = this->C.protectedVertexInsert(d);

    /* if vertex has been inserted successfully, convert returned Graph::vertex_iterator to dendrite_iterator using the
     * returned internal map iterator referring to Graph::V (this is not a Graph::vertex_iterator, which does not allow
     * access to internals). */
    if (rpair.second) {
        return dendrite_iterator( &(this->C), rpair.first );
    }
    else {
        throw("CellNetwork::DendriteAccessor::insert(): failed to insert dendrite vertex.");
    }
}

/*
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteAccessor::erase(dendrite_iterator it)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteAccessor::erase(uint32_t id)
{
}
*/

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of NeuronEdgeAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdgeAccessor::NeuronEdgeAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        NeuronEdge,
        typename Graph<Tn, Tv, Te>::Edge,
        NeuronEdge,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
        NeuronEdgePred,
        neuron_edge_iterator,
        neuron_edge_const_iterator,
        Te
    >(C, C.E, C.isNeuronEdge)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neuron_edge_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdgeAccessor::erase(neuron_edge_iterator it)
{

    /* "up"-convert neuron_edge_iterator to Graph::edge_iterator by passing pointer to CellNetwork, which is upcast to
     * Graph, and internal iterator */
    typename Graph<Tn, Tv, Te>::edge_iterator eit(&(this->C), it.int_it);

    try {
        eit = this->C.edges.erase(eit);
    }
    catch (GraphEx& ex) {
        /* FIXME: exception handling, wrap into CellNetwork exception class */
        throw;
    }
    /* convert back "down" by wrapping internal iterator in a neuron_edge_iterator. */
    return neuron_edge_iterator(&(this->C), Graph<Tn, Tv, Te>::getInternalIterator(eit));
}

/* FIXME: implement me */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
bool
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuronEdgeAccessor::erase(uint32_t id)
{
    return true;
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of NeuriteRootEdgeAccessor
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdgeAccessor::NeuriteRootEdgeAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        NeuriteRootEdge,
        typename Graph<Tn, Tv, Te>::Edge,
        NeuronEdge,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
        NeuriteRootEdgePred,
        neurite_rootedge_iterator,
        neurite_rootedge_const_iterator,
        Te
    >(C, C.E, C.isNeuriteRootEdge)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<
    typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_rootedge_iterator,
    std::list<typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_segment_iterator>
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteRootEdgeAccessor::split(
    neurite_rootedge_iterator           nr_it,
    std::vector<
            std::pair<
                Vec3<R>,
                R
            >
        > const                        &intermediate_vertex_data)
{
    debugl(1, "CellNetwork::NeuriteRootEdgeAccessor::split(). neurite root edge: e = (u, v): %d = (%d, %d).\n",
            nr_it->id(),
            nr_it->getSourceVertex()->id(),
            nr_it->getDestinationVertex()->id());
    debugTabInc();

    if (intermediate_vertex_data.empty()) {
        throw("CellNetwork::NeuriteRootEdgeAccessor::split(): no information about intermediate vertices given. won't silently discard call due to likely erroneous input.");
    }

    /* insert intermediate points p1, .., pn in-between neurite root edge nr = (s, v), forming the edges (s, p1) - (p1,
     * p2) - .. - (pn, v). the data attached to the old neurite root edge nr is copied to the new neurite root edge
     * (s, p1) in this process. */
    uint32_t const                  n       = intermediate_vertex_data.size();
    soma_iterator                   s_it    = nr_it->getSourceVertex();
    neurite_iterator                v_it    = nr_it->getDestinationVertex();
    std::vector<neurite_iterator>   intermediate_vertices(n);
    neurite_rootedge_iterator       nr_new_it;

    /* depending on the type of neurite root edge (and equivalently, the type of v), insert Axon- or DendriteVertices
     * with default constructed instances for Tv and (Ta (x)or Td). */
    debugl(2, "adding %d intermediate vertices..\n", n);
    debugTabInc();
    if (v_it->getType() == this->C.AXON_VERTEX) {
        debugl(2, "(s, v) was an axon root edge => inserting intermediate axon vertices.\n");
        AxonVertex *v_downcast = dynamic_cast<AxonVertex *>(&(*v_it));
        if (v_downcast) {
            axon_iterator           p0a_it;
            axon_rootedge_iterator  s_p0a_it;
            Tar                     nra_data = (v_downcast->template getFilteredInEdges<AxonRootEdge>()).front()->axon_root_edge_data;

            debugl(2, "erasing old neurite root edge e = (u, v) = %d\n", nr_it->id());
            this->C.neuron_edges.erase(nr_it);

            for (uint32_t i = 0; i < n; i++) {
                if (i == 0) {
                    p0a_it = 
                        this->C.axon_vertices.insert(
                            /* FIXME: which value for compartment id? keep IdQueue for compartment ids.. */
                            CellSection(
                                0,
                                intermediate_vertex_data[i].first,
                                intermediate_vertex_data[i].second
                            )
                        );

                    /* also stored upcast neurite_iterator in intermediate vertices */
                    intermediate_vertices[0] = p0a_it;

                    /* insert new neurite root edge (s, p0d) using copied DendriteRootEdgeInfo */
                    auto rpair = 
                        this->C.axon_root_edges.insert(
                                s_it,
                                p0a_it,
                                nra_data);

                    if (!rpair.second) {
                        throw("CellNetwork::NeuriteRootEdgeAccessor::split(): couldn't insert fresh axon root segment"\
                            " in the process of splitting given one.");
                    }

                    /* save axon root edge iterator as neurite root edge iterator, used below. */
                    s_p0a_it    = rpair.first;
                    nr_new_it   = s_p0a_it;
                }
                else {
                    intermediate_vertices[i] =
                        this->C.axon_vertices.insert(
                            /* FIXME: which value for compartment id? keep IdQueue for compartment ids.. */
                            CellSection(
                                0,
                                intermediate_vertex_data[i].first,
                                intermediate_vertex_data[i].second
                            )
                        );
                }
            }
        }
        else {
            throw("NLM_CellNetwork::splitNeuriteRootEdge(): failed to down-cast source vertex to indicated type "\
                "DendriteVertex. internal logic error.");
        }
    }
    else if (v_it->getType() == this->C.DENDRITE_VERTEX) {
        debugl(2, "(s, v) was a dendrite root edge => inserting intermediate dendrite vertices.\n");
        DendriteVertex *v_downcast = dynamic_cast<DendriteVertex *>(&(*v_it));
        if (v_downcast) {
            dendrite_iterator           p0d_it;
            dendrite_rootedge_iterator  s_p0d_it;
            Tdr                         nrd_data = (v_downcast->template getFilteredInEdges<DendriteRootEdge>()).front()->dendrite_root_edge_data;

            debugl(2, "erasing old neurite root edge e = (u, v) = %d\n", nr_it->id());
            this->C.neuron_edges.erase(neuron_edge_iterator(nr_it));

            for (uint32_t i = 0; i < n; i++) {
                if (i == 0) {
                    p0d_it = 
                        this->C.dendrite_vertices.insert(
                            /* FIXME: which value for compartment id? keep IdQueue for compartment ids.. */
                            CellSection(
                                0,
                                intermediate_vertex_data[i].first,
                                intermediate_vertex_data[i].second
                            ),
                            /* apical dendrite flag from u / v */
                            v_downcast->isApicalDendrite()
                        );

                    /* also stored upcast neurite_iterator in intermediate vertices */
                    intermediate_vertices[0] = p0d_it;

                    /* insert new neurite root edge (s, p0d) using copied DendriteRootEdgeInfo */
                    auto rpair = 
                        this->C.dendrite_root_edges.insert(
                                s_it,
                                p0d_it,
                                nrd_data);

                    if (!rpair.second) {
                        throw("CellNetwork::NeuriteRootEdgeAccessor::split(): couldn't insert fresh dendrite root segment"\
                            " in the process of splitting given one.");
                    }

                    /* save dendrite root edge iterator as neurite root edge iterator, used below. */
                    s_p0d_it    = rpair.first;
                    nr_new_it   = s_p0d_it;
                }
                else {
                    intermediate_vertices[i] =
                        this->C.dendrite_vertices.insert(
                            /* FIXME: which value for compartment id? keep IdQueue for compartment ids.. */
                            CellSection(
                                0,
                                intermediate_vertex_data[i].first,
                                intermediate_vertex_data[i].second
                            ),
                            /* apical dendrite flag from u / v */
                            v_downcast->isApicalDendrite()
                        );
                }
            }
        }
        else {
            throw("NLM_CellNetwork::splitNeuriteRootEdge(): failed to down-cast source vertex to indicated type "\
                "DendriteVertex. internal logic error.");
        }
    }
    else {
        throw("NLM_CellNetwork::splitNeuriteRootEdge(): source and destination vertex of given neurite segment neither "\
            "both of type AxonVertex nor both of type DendriteVertex. internal logic error.");
    }
    debugTabDec();

    /* insert new neurite segments, on the neurite-segment abstraction level.. */
    debugl(2, "inserting new neurite segments (p1, p_2), (p1, p2), .., (p_n, v)\n");
    debugTabInc();

    std::list<neurite_segment_iterator> ns_list;
    for (uint32_t i = 0; i < n - 1; i++) {
        debugl(2, "inserting ns (%d, %d)\n",
                intermediate_vertices[i]->id(),
                intermediate_vertices[i+1]->id()
            );
        ns_list.push_back( (this->C.neurite_segments.insert(intermediate_vertices[i], intermediate_vertices[i+1])).first );
    }

    debugl(2, "inserting ns (%d, %d)\n", intermediate_vertices[n-1]->id(), v_it->id());
    ns_list.push_back( (this->C.neurite_segments.insert(intermediate_vertices[n-1], v_it)).first );

    debugTabDec();
    debugl(2, "done inserting new neurite segments..\n");

    /* network info needs update, since neurite information has changed. */
    this->C.network_info_initialized = false;

    debugTabDec();
    debugl(1, "CellNetwork::NeuriteRootEdgeAccessor::split(): done.\n");

    return std::pair<
            neurite_rootedge_iterator,
            std::list<neurite_segment_iterator>
        >(nr_new_it, ns_list);
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of AxonRootEdgeAccessor
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdgeAccessor::AxonRootEdgeAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        AxonRootEdge,
        typename Graph<Tn, Tv, Te>::Edge,
        NeuronEdge,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
        AxonRootEdgePred,
        axon_rootedge_iterator,
        axon_rootedge_const_iterator,
        Te
    >(C, C.E, C.isAxonRootEdge)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<
    typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_rootedge_iterator,
    bool
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonRootEdgeAccessor::insert(
    soma_iterator      &v_src_it,
    axon_iterator      &v_dst_it,
    Tar const          &axon_root_edge_data,
    Tnr const          &neurite_root_edge_data,
    Te const           &edge_data)
{
    /* check if both refer to (this) container */
    if (!v_src_it.checkContainer(this->C) || !v_src_it.sameContainer(v_dst_it)) {
        throw("CellNetwork::AxonRootEdgeAccessor::insert(): at least one input vertex does not refer to this (this) CellNetwork.");
    }

    std::pair<CellNetwork::axon_rootedge_iterator, bool> ret_pair;

    AxonRootEdge *a     = new AxonRootEdge(&(this->C), &(*v_src_it), &(*v_dst_it), axon_root_edge_data, neurite_root_edge_data, edge_data);
    auto rpair          = this->C.protectedEdgeInsert(a, v_src_it->id(), v_dst_it->id() );

    /* if vertex has been inserted successfully, convert returned Graph::vertex_iterator to dendrite_iterator using the
     * returned internal map iterator referring to Graph::V (this is not a Graph::vertex_iterator, which does not allow
     * access to internals). */
    if (rpair.second) {
        ret_pair.first  = axon_rootedge_iterator( &(this->C), rpair.first );
        ret_pair.second = true;
        return ret_pair;
    }
    else {
        ret_pair.first.explicitlyInvalidate();
        ret_pair.second = false;
        return ret_pair;
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of DendriteRootEdgeAccessor
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdgeAccessor::DendriteRootEdgeAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        DendriteRootEdge,
        typename Graph<Tn, Tv, Te>::Edge,
        NeuronEdge,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
        DendriteRootEdgePred,
        dendrite_rootedge_iterator,
        dendrite_rootedge_const_iterator,
        Te
    >(C, C.E, C.isDendriteRootEdge)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<
    typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_rootedge_iterator,
    bool
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteRootEdgeAccessor::insert(
    soma_iterator       v_src_it,
    dendrite_iterator   v_dst_it,
    Tdr const          &dendrite_root_edge_data,
    Tnr const          &neurite_root_edge_data,
    Te const           &edge_data)
{
    /* check if both refer to (this) container */
    if (!v_src_it.checkContainer(this->C) || !v_src_it.sameContainer(v_dst_it)) {
        throw("CellNetwork::DendriteRootEdgeAccessor::insert(): at least one input vertex does not refer to this (this) CellNetwork.");
    }

    std::pair<CellNetwork::dendrite_rootedge_iterator, bool> ret_pair;

    DendriteRootEdge *d = new DendriteRootEdge(&(this->C), &(*v_src_it), &(*v_dst_it), dendrite_root_edge_data, neurite_root_edge_data, edge_data);
    auto rpair          = this->C.protectedEdgeInsert(d, v_src_it->id(), v_dst_it->id() );

    /* if vertex has been inserted successfully, convert returned Graph::vertex_iterator to dendrite_iterator using the
     * returned internal map iterator referring to Graph::V (this is not a Graph::vertex_iterator, which does not allow
     * access to internals). */
    if (rpair.second) {
        ret_pair.first  = dendrite_rootedge_iterator( &(this->C), rpair.first );
        ret_pair.second = true;
        return ret_pair;
    }
    else {
        ret_pair.first.explicitlyInvalidate();
        ret_pair.second = false;
        return ret_pair;
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of NeuriteSegmentAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegmentAccessor::NeuriteSegmentAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        NeuriteSegment,
        typename Graph<Tn, Tv, Te>::Edge,
        NeuronEdge,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
        NeuriteSegmentPred,
        neurite_segment_iterator,
        neurite_segment_const_iterator,
        Te
    >(C, C.E, C.isNeuriteSegment)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<
    typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_segment_iterator,
    bool
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegmentAccessor::insert(
    neurite_iterator    v_src_it,
    neurite_iterator    v_dst_it,
    Tns const          &neurite_segment_data,
    Te const           &edge_data)
{
    /* check if v_src and v_dst are both of type axon / dendrite and forward call accordingly */
    if (v_src_it->getType() == AXON_VERTEX && v_dst_it->getType() == AXON_VERTEX) {
        /* cast down to axon vertices */
        AxonVertex *v_src_axon = dynamic_cast<AxonVertex *>(&(*v_src_it));
        AxonVertex *v_dst_axon = dynamic_cast<AxonVertex *>(&(*v_dst_it));

        if (v_src_axon && v_dst_axon) {
            auto rpair = this->C.axon_segments.insert(
                v_src_axon->iterator(),
                v_dst_axon->iterator(),
                Tas(),
                neurite_segment_data,
                edge_data);

            return std::pair<neurite_segment_iterator, bool>(rpair.first, rpair.second);
        }
        else {
            throw("CellNetwork::NeuriteSegmentAccessor::insert(): two given neurite vertex iterators both "\
                "indicate type AxonVertex, yet downcasting of at least one vertex to AxonVertex failed. internal "\
                "logic error.");
        }
    }
    else if (v_src_it->getType() == DENDRITE_VERTEX && v_dst_it->getType() == DENDRITE_VERTEX) {
        DendriteVertex *v_src_dendrite = dynamic_cast<DendriteVertex *>(&(*v_src_it));
        DendriteVertex *v_dst_dendrite = dynamic_cast<DendriteVertex *>(&(*v_dst_it));

        if (v_src_dendrite && v_dst_dendrite) {
            auto rpair = this->C.dendrite_segments.insert(
                v_src_dendrite->iterator(),
                v_dst_dendrite->iterator(),
                Tds(),
                neurite_segment_data,
                edge_data);

            return std::pair<neurite_segment_iterator, bool>(rpair.first, rpair.second);
        }
        else {
            throw("CellNetwork::NeuriteSegmentAccessor::insert(): two given neurite vertex iterators both "\
                "indicate type DenriteVertex, yet downcasting of at least one vertex to DendriteVertex failed. internal "\
                "logic error.");
        }
    }
    else {
        throw("CellNetwork::NeuriteSegmentAccessor::insert(): two given neurite vertex iterators neither both "\
            "refer to axon nor both to dendrite vertices. semantically invalid input data.");
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::list<typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_segment_iterator>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegmentAccessor::split(
    neurite_segment_iterator            ns_it,
    std::vector<
            std::pair<
                Vec3<R>,
                R
            >
        > const                        &intermediate_vertex_data)
{
    debugl(1, "CellNetwork::NeuriteSegmentAccessor::split(). neurite segment: e = (u, v): %d = (%d, %d).\n",
            ns_it->id(),
            ns_it->getSourceVertex()->id(),
            ns_it->getDestinationVertex()->id());
    debugTabInc();

    if (intermediate_vertex_data.empty()) {
        throw("CellNetwork::NeuriteSegmentAccessor::split(): no information about intermediate vertices given. won't silently discard call due to likely erroneous input.");
    }

    /* insert intermediate points p1, .., pn in-between neurite segments ns = (u, v), forming
     * the neurite segments (u, p1) - (p1, p2) - .. - (pn, v). the data attached to the neurite segment ns is lost in
     * this process. */
    uint32_t const                  n = intermediate_vertex_data.size();
    neurite_iterator                u = ns_it->getSourceVertex();
    neurite_iterator                v = ns_it->getDestinationVertex();
    std::vector<neurite_iterator>   intermediate_vertices(n);

    debugl(2, "erasing old edge e = (u, v) = %d\n", ns_it->id());
    this->C.neuron_edges.erase(neuron_edge_iterator(ns_it));

    /* depending on the type of ns (and equivalently, the types of u and v), insert Axon- or DendriteVertices with
     * default constructed instances for Tv and (Ta (x)or Td). */
    debugl(2, "adding %d intermediate vertices..\n", n);
    debugTabInc();
    if (u->getType() == this->C.AXON_VERTEX && v->getType() == this->C.AXON_VERTEX) {
        debugl(2, "(u, v) was an axon segment => inserting intermediate axon vertices.\n");
        for (uint32_t i = 0; i < n; i++) {
            intermediate_vertices[i] = this->C.axon_vertices.insert(
                    /* FIXME: which value for compartment id? keep IdQueue for compartment ids.. */
                    CellSection(
                        0,
                        intermediate_vertex_data[i].first,
                        intermediate_vertex_data[i].second
                    )
                );
        }
    }
    else if (u->getType() == this->C.DENDRITE_VERTEX && v->getType() == this->C.DENDRITE_VERTEX) {
        debugl(2, "(u, v) was a dendrite segment => inserting intermediate dendrite vertices.\n");
        DendriteVertex *u_downcast = dynamic_cast<DendriteVertex *>(&(*u));
        if (u_downcast) {
            for (uint32_t i = 0; i < n; i++) {
                intermediate_vertices[i] = this->C.dendrite_vertices.insert(
                        /* FIXME: which value for compartment id? keep IdQueue for compartment ids.. */
                        CellSection(
                            0,
                            intermediate_vertex_data[i].first,
                            intermediate_vertex_data[i].second
                        ),
                        /* apical dendrite flag from u / v */
                        u_downcast->isApicalDendrite()
                    );
            }
        }
        else {
            throw("NLM_CellNetwork::splitNeuriteSegment(): failed to down-cast source vertex to indicated type "\
                "DendriteVertex. internal logic error.");
        }
    }
    else {
        throw("NLM_CellNetwork::splitNeuriteSegment(): source and destination vertex of given neurite segment neither "\
            "both of type AxonVertex nor both of type DendriteVertex. internal logic error.");
    }
    debugTabDec();

    /* set neurite and soma information inside newly added neurite vertices */
    debugl(2, "setting soma / neurite information in intermediate vertices.\n");
    for (uint32_t i = 0; i < n; i++) {
        intermediate_vertices[i]->soma      = u->getSoma();
        intermediate_vertices[i]->neurite   = u->getNeurite();
    }

    /* insert new neurite segments, on the neurite-segment abstraction level.. */
    debugl(2, "inserting new neurite segments (u, p_1), (p1, p2), .., (p_n, v)\n");
    debugTabInc();

    std::list<neurite_segment_iterator> ret;
    debugl(2, "inserting ns (%d, %d)\n", u->id(), intermediate_vertices[0]->id());
    ret.push_back( (this->C.neurite_segments.insert(u, intermediate_vertices[0])).first );

    for (uint32_t i = 0; i < n - 1; i++) {
        debugl(2, "inserting ns (%d, %d)\n",
                intermediate_vertices[i]->id(),
                intermediate_vertices[i+1]->id()
            );
        ret.push_back( (this->C.neurite_segments.insert(intermediate_vertices[i], intermediate_vertices[i+1])).first );
    }

    debugl(2, "inserting ns (%d, %d)\n", intermediate_vertices[n-1]->id(), v->id());
    ret.push_back( (this->C.neurite_segments.insert(intermediate_vertices[n-1], v)).first );

    debugTabDec();
    debugl(2, "done inserting new neurite segments..\n");

    /* set neurite and soma information for all newly added neurite segments and return list of iterators */
    for (auto &ns : ret) {
        ns->soma    = u->getSoma();
        ns->neurite = u->getNeurite();
    }

    debugTabDec();
    debugl(1, "CellNetwork::NeuriteSegmentAccessor::split(): done.\n");
    return ret;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::neurite_iterator
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
NeuriteSegmentAccessor::collapse(
    neurite_segment_iterator            ns_it,
    std::pair<Vec3<R>, R> const        &collapsed_vertex_data)
{
    /* perform collapse of edge ns_it on Graph abstraction level, store return vertex_iterator */
    typename Graph<Tn, Tv, Te>::vertex_iterator g_it =
        this->C.Graph<Tn, Tv, Te>::edges.collapse(ns_it);

    /* convert returned Graph::vertex_iterator to neurite vertex iterator. */
    neurite_iterator c_it(&(this->C), Graph<Tn, Tv, Te>::getInternalIterator(g_it) );

    /* call polymorphic updateVertexPointers() method on all neighbours of collapsed vertex c */
    std::list<NeuriteSegment *> nbs, tmp;
    c_it->template getFilteredInEdges<NeuriteSegment>(nbs);
    c_it->template getFilteredOutEdges<NeuriteSegment>(tmp);
    nbs.insert(nbs.end(), tmp.begin(), tmp.end());
    
    for (auto &nb : nbs) {
        nb->updateVertexPointers();
    }

    /* set position / radius supplied by the caller */
    c_it->setPosition(collapsed_vertex_data.first);
    c_it->setRadius(collapsed_vertex_data.second);

    return c_it;
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of AxonSegmentAccessor 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegmentAccessor::AxonSegmentAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        AxonSegment,
        typename Graph<Tn, Tv, Te>::Edge,
        NeuronEdge,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
        AxonSegmentPred,
        axon_segment_iterator,
        axon_segment_const_iterator,
        Te
    >(C, C.E, C.isAxonSegment)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<
    typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::axon_segment_iterator,
    bool
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
AxonSegmentAccessor::insert(
    axon_iterator       v_src_it,
    axon_iterator       v_dst_it,
    Tas const          &axon_segment_data,
    Tns const          &neurite_segment_data,
    Te const           &edge_data)
{
    /* check if both refer to (this) container */
    if (!v_src_it.checkContainer(this->C) || !v_src_it.sameContainer(v_dst_it)) {
        throw("CellNetwork::AxonSegmentAccessor::insert(): at least one axon input vertex does not refer to this (this) CellNetwork.");
    }

    std::pair<CellNetwork::axon_segment_iterator, bool> ret_pair;

    AxonSegment *a      = new AxonSegment(&(this->C), &(*v_src_it), &(*v_dst_it), axon_segment_data, neurite_segment_data, edge_data);
    auto rpair          = this->C.protectedEdgeInsert(a, v_src_it->id(), v_dst_it->id() );

    /* if vertex has been inserted successfully, convert returned Graph::vertex_iterator to dendrite_iterator using the
     * returned internal map iterator referring to Graph::V (this is not a Graph::vertex_iterator, which does not allow
     * access to internals). */
    if (rpair.second) {
        ret_pair.first  = axon_segment_iterator( &(this->C), rpair.first );
        ret_pair.second = true;
        return ret_pair;
    }
    else {
        ret_pair.first.explicitlyInvalidate();
        ret_pair.second = false;
        return ret_pair;
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of DendriteSegmentAccessor
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegmentAccessor::DendriteSegmentAccessor(CellNetwork &C)
    /* call base CellNetworkAccessor constructor, which initializes the references */
    : CellNetworkAccessor<
        DendriteSegment,
        typename Graph<Tn, Tv, Te>::Edge,
        NeuronEdge,
        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
        DendriteSegmentPred,
        dendrite_segment_iterator,
        dendrite_segment_const_iterator,
        Te
    >(C, C.E, C.isDendriteSegment)
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<
    typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::dendrite_segment_iterator,
    bool
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
DendriteSegmentAccessor::insert(
    dendrite_iterator   v_src_it,
    dendrite_iterator   v_dst_it,
    Tds const          &dendrite_segment_data,
    Tns const          &neurite_segment_data,
    Te const           &edge_data)
{
    /* check if both refer to (this) container */
    if (!v_src_it.checkContainer(this->C) || !v_src_it.sameContainer(v_dst_it)) {
        throw("CellNetwork::DendriteSegmentAccessor::insert(): at least one input dendrite vertex does not refer to this (this) CellNetwork.");
    }

    std::pair<CellNetwork::dendrite_segment_iterator, bool> ret_pair;

    /* check if both are apical or basal dendrites */
    if (Aux::Logic::lequiv(v_src_it->isApicalDendrite(), v_dst_it->isApicalDendrite())) {
        DendriteSegment *d  = new DendriteSegment(&(this->C), &(*v_src_it), &(*v_dst_it), dendrite_segment_data, neurite_segment_data, edge_data);
        auto rpair          = this->C.protectedEdgeInsert(d, v_src_it->id(), v_dst_it->id() );

        /* if vertex has been inserted successfully, convert returned Graph::vertex_iterator to dendrite_iterator using the
         * returned internal map iterator referring to Graph::V (this is not a Graph::vertex_iterator, which does not allow
         * access to internals). */
        if (rpair.second) {
            ret_pair.first  = dendrite_segment_iterator( &(this->C), rpair.first );
            ret_pair.second = true;
            return ret_pair;
        }
        else {
            ret_pair.first.explicitlyInvalidate();
            ret_pair.second = false;
            return ret_pair;
        }
    }
    else {
        throw("CellNetwork::DendriteSegmentAccessor::insert(): input vertices invalid. one vertex belongs to a basal, the other to an apical dendrite.");
    }
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        implementation of CellNetwork interface 
 *
 * ----------------------------------------------------------------------------------------------------------------- */

// obviously unused as w.iterator in the return statement cannot be resolved
#if 0
/* two selection predicate generation template methods from namespace CellNetwork */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename VertexType>
std::function<bool(VertexType const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
genVertexSelectionPred(VertexType const &v)
{
    return std::function<bool(NeuronVertex const &v)>(
            [v] (NeuronVertex const &w) -> bool
            {
                return (v.iterator == w.iterator);
            }
        );
}
#endif

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template<typename EdgeType>
std::function<bool(EdgeType const &)>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
genEdgeSelectionPred(EdgeType const &e)
{
    return std::function<bool(EdgeType const &v)>(
            [e] (EdgeType const &f) -> bool
            {
                return (e.iterator == f.iterator);
            }
        );
}

/* protected methods */

/* check network topology */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
checkTopology()
{
    /* verify that every soma hast at most one axon */
    for (auto &s : this->soma_vertices) {
        if (s.template getFilteredOutNeighbours<AxonVertex>().size() > 1) {
            throw("CellNetwork::checkTopology(): discovered soma with more than one axon stem.");
        }
    }

    /* verify that the subgraph consisting only of 
     *
     *  1. soma and neurite vertices and
     *
     *  2. neurite root edges and neurite segments
     *
     *  is a forest where all trees are rooted in soma vertices. only synapse edges (or later maybe other specialized
     *  morphology information edges) inter-connect different cell trees. a simple test uses the fact that the
     *  soma/neurite subgraph is a forest with the desired properties if and only if every soma is a source and every
     *  neurite vertex has in-degree exactly one. in particular, there can't be two soma-trees that can be
     *  interconnected by a neurite segment and no neurite vertex acting as the "root" of a tree in the forest (this is
     *  sometimes found in the original files from NeuroMorpho.org, where unconnected neurite fragment from some other
     *  cell is encoded in a file). */

    for (auto &nv : neurite_vertices) {
        /* get "in-degree" of neurite vertex when restricting graph to soma and neurite vertices. */
        if (nv.template getFilteredInEdges<NeuriteRootEdge>().size() + 
            nv.template getFilteredInEdges<NeuriteSegment>().size()     != 1) 
        {
            throw("CellNetwork::checkTopology(): discovered neurite vertex with in-degree != 1 when restricted to "\
                "soma / neurite vertices=> cell network does not exhibit required structure: the sub-graph "\
                "consisting of soma and neurite vertices must be a forest of neurite trees rooted in the somas. "\
                "semantic error.\n");
        }
    }
}

/* initialize network information */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
initializeNetworkInfo()
{
    if (!this->network_info_initialized) {
        debugl(0, "CellNetwork::initializeNetworkInfo().\n");
        debugTabInc();
        /* for all cells (i.e. somas), traverse each neurite individually and store iterators to soma and respective
         * neurite root edge in neurite vertices / neurite segment.  this makes later access much more efficient and
         * comfortable. */
        debugl(1, "iterating over all cells (somas)..\n");
        debugTabInc();
        for (auto &s : this->soma_vertices) {
            debugl(1, "processing soma %d\n", s.id());
            std::list<NeuriteRootEdge *>    s_neurite_root_edges;

            s.template getFilteredOutEdges<NeuriteRootEdge>(s_neurite_root_edges);

            debugl(1, "processing all neurites connected to soma %d\n", s.id());
            debugTabInc();
            for (auto &nre : s_neurite_root_edges) {
                debugl(2, "processing neurite root edge %d with neurite root vertex r %d. retrieving connected component..\n", nre->id(), nre->getDestinationVertex()->id());

                /* extract neurite root vertex r of neurite starting with edge nre */
                neurite_iterator r      = nre->getDestinationVertex();

                debugl(2, "retrieving all neurite vertices / edges reachable from neurite roor vertex %d.\n", r->id());
                /* get all neurite vertices / segments reachable from r */
                std::list<NeuriteVertex *>  r_cc_nv;
                std::list<NeuriteSegment *> r_cc_ns;

                this->getNeuriteConnectedComponent(r, r_cc_nv, r_cc_ns);

                /* store info inside all reachable neurite vertices */
                debugl(2, "assigning soma / neurite information to all %d neurite vertices reachable from r(%d).\n", r_cc_nv.size(), r->id());
                debugTabInc();
                for (auto &nv : r_cc_nv) {
                    debugl(3, "processing neurite vertex %d.\n", nv->id());
                    nv->soma        = s.iterator();
                    nv->neurite     = nre->iterator();
                }
                debugTabDec();
                debugl(2, "done assigning neurite vertex info.\n");

                /* store info inside all reachable neurite edges */
                debugl(2, "assigning soma / neurite information to all %d neurite segments reachable from r(%d).\n", r_cc_ns.size(), r->id());
                debugTabInc();
                for (auto &ns : r_cc_ns) {
                    debugl(3, "processing neurite segment %d.\n", ns->id());
                    ns->soma    = s.iterator();
                    ns->neurite = nre->iterator();
                }
                debugTabDec();
            }
            debugTabDec();
            debugl(1, "done with soma %d.\n", s.id());
        }
        debugTabDec();
        debugl(0, "done processing somas.\n");

        this->network_info_initialized = true;

        debugTabDec();
        debugl(0, "CellNetwork::initializeNetworkInfo(): done.\n");
    }
}

/* protected casting methods */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <typename VertexType>
VertexType *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
downcastVertex(NeuronVertex *v)
{
    return (dynamic_cast<VertexType *>(v));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <typename EdgeType>
EdgeType *
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
downcastEdge(NeuronEdge *e)
{
    return (dynamic_cast<EdgeType *>(e));
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getVertexType(neuron_const_iterator const &v_it) const
{
    return (v_it->getType());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
uint32_t
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getEdgeType(neuron_edge_const_iterator const &e_it) const
{
    return (e_it->getType());
}

/* public method interface */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
CellNetwork(const std::string& _network_name)
: Graph<Tn, Tv, Te>(),
network_name(_network_name), network_info_initialized(false),
global_coordinate_displacement(Aux::VecMat::nullvec<R>()), global_coordinate_scaling_factor(1.0),
neuron_vertices(*this), soma_vertices(*this),
neurite_vertices(*this), axon_vertices(*this), dendrite_vertices(*this),
neuron_edges(*this), neurite_root_edges(*this), axon_root_edges(*this), dendrite_root_edges(*this),
neurite_segments(*this), axon_segments(*this), dendrite_segments(*this)
{}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
~CellNetwork()
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
template <typename VertexType, typename EdgeType>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
traverseBreadthFirst(
    neuron_iterator                                     vstart_it,
    bool const                                         &directed,
    bool const                                         &return_paths,
    std::list<
            std::tuple<
                VertexType *,
                std::list<NeuronVertex *>,
                R
            >
        >                                              &reachable_vertices_info,
    std::list<EdgeType *>                              &reachable_edges,
    std::function<bool(NeuronVertex const &v)> const   &vertex_pred,
    std::function<bool(NeuronEdge const &e)> const     &edge_pred,
    std::function<bool(VertexType const &v)> const     &vertex_selection_pred,
    std::function<bool(EdgeType const &e)> const       &edge_selection_pred,
    std::function<bool(VertexType const &v)> const     &vertex_termination_pred,
    std::function<bool(EdgeType const &e)> const       &edge_termination_pred,
    int32_t                                             tid_arg)
{
    debugl(0, "CellNetwork::getConnectedComponent(): starting %s traversal from node %2d.\n",
        (directed) ? "directed" : "undirected", vstart_it->id());
    debugTabInc();

    /* queue consisting of tuples (iterator, path, distance), where iterator identifies the currently visited vertex v,
     * path is path from vstart to v encoded as a list of NeuronVertex pointers and dist is the d(vstart, v) */
    std::list<
            std::tuple<
                neuron_iterator,
                std::list<NeuronVertex *>,
                R
            >
        >                                       Q;

    uint32_t                                    tid;
    neuron_iterator                             v_it, nb_it;
    VertexType                                 *v_downcast;
    std::list<NeuronVertex *>                   v_path;
    R                                           v_dist, nb_dist;
    std::list<NeuronEdge *>                     tmp_edges, v_filtered_out_edges, v_filtered_in_edges;

    /* if tid_arg == -1, get fresh traversal id, otherwise use tid (or throw). */
    if (tid_arg == -1) {
        tid = this->getFreshTraversalId();  
    }
    else if (tid_arg >= 0) {
        tid = tid_arg;
    }
    else {
        throw("CellNetwork::traverseBreadthFirst(): supplied traversal id is neither -1 (to indicate that fresh "\
            "traversal id should be fetched internaly) nor >= 0.\n");
    }

    /* initialize Q. note that tuple constructor is explicit (for whatever reason) and one is forced to use the ugly
     * explicit construction instead of the much more elegant initialization list(which works for e.g. pairs).
     * more precisely, this would be logical, but for some rather shallow technical reason this is not allowed for
     * tuples.. but for pairs
     *
     * Q = { { vstart_it, { &(*vstart_it) }, 0) } };
     * */
    Q = { std::tuple<neuron_iterator, std::list<NeuronVertex *>, R>(vstart_it, { &(*vstart_it) }, 0) };

    /* set vstart to enqueued. */
    vstart_it->setTraversalState(tid, this->TRAV_ENQUEUED);

    while (!Q.empty()) {
        /* extract info and dequeue front() tuple */
        v_it    = std::get<0>(Q.front());
        v_path  = std::get<1>(Q.front());
        v_dist  = std::get<2>(Q.front());
        Q.pop_front();

        /* if v is both of matching type VertexType and satisfied the vertex selection predicate, insert it into
         * reachable_vertices_info with empty path and distance 0. */
        v_downcast = dynamic_cast<VertexType *>(&(*v_it));

        /* NOTE: moderately "subtle" lazy evaluation issue here.. */
        if (v_downcast && vertex_selection_pred(*v_downcast)) {
            /* return path to v if desired by the caller */
            if (return_paths) {
                reachable_vertices_info.push_back(
                        std::tuple<VertexType *, std::list<NeuronVertex *>, R>
                            (v_downcast, v_path, v_dist)
                    );
            }
            /* don't return path to v, init second element of tuple to empty list. */
            else {
                reachable_vertices_info.push_back(
                        std::tuple<VertexType *, std::list<NeuronVertex *>, R>
                            (v_downcast, {}, v_dist)
                    );
            }
        }

        /* if v is of matching type VertexType and satisfies vertex_termination_pred, abort the traversal */
        v_downcast = dynamic_cast<VertexType *>(&(*v_it));
        if (v_downcast && vertex_termination_pred(*v_downcast)) {
            break;
        }

        /* termination criterion not satisfied. proceed with the traversal: clear data from previous vertex */
        v_filtered_out_edges.clear();
        v_filtered_in_edges.clear();

        /* for both directed and undirected version, get all of v's OUT edges and extract those satisfying edge_pred */
        v_it->template getOutEdges(tmp_edges);
        for (auto &e : tmp_edges) {
            if (edge_pred(*e)) {
                v_filtered_out_edges.push_back(e);
            }
        }

        /* get the destination vertex of all filtered out edges, i.e. the neighbour connected to v by this edge */
        for (auto &e : v_filtered_out_edges) {
            nb_it = e->getDestinationVertex();

            /* consider neighbour only if vertex predicate returns true */
            if (vertex_pred(*nb_it)) {
                /* if type of neighbour matches VertexType, the edge e connects two vertices satisfying
                 * vertex_pred. if additionally, e itself satisfies edge_pred AND edge_selection_pred, append it to
                 * reachable_edges return list, yet only if its traversal state is TRAV_UNSEEN. */
                VertexType *nb_downcast = dynamic_cast<VertexType *>(&(*nb_it));
                if (nb_downcast) {
                    EdgeType *e_downcast = dynamic_cast<EdgeType *>(e);
                    /* NOTE: moderately "subtle" lazy evaluation issue here.. */
                    if (e_downcast &&
                        edge_selection_pred(*e_downcast) &&
                        e_downcast->getTraversalState(tid) == this->TRAV_UNSEEN)
                    {
                        reachable_edges.push_back(e_downcast);

                        /* set traversal state of edge to DONE to prevent double inserts. other than that, edge
                         * traversal states are not used for anything. in particular, the edge is not blocked (or
                         * similar) afterwards. */
                        e_downcast->setTraversalState(tid, this->TRAV_DONE);
                    }
                }

                /* extent path to v to get path to neighbour */
                std::list<NeuronVertex *> nb_path = v_path;
                nb_path.push_back(&(*nb_it));

                /* calculate distance to neighbour */
                nb_dist = v_dist + e->getLength();

                /* if neighbour has not been traversed yet, enqueue it */
                if (nb_it->getTraversalState(tid) == this->TRAV_UNSEEN) {
                    Q.push_back(
                        std::tuple<neuron_iterator, std::list<NeuronVertex *>, R>
                            (nb_it, nb_path, nb_dist)
                        );
                }
            }
        }

        /* for undirected version, also process v's IN edges in a completely symmetric way. */
        if (!directed) {
            v_it->template getInEdges(tmp_edges);
            for (auto &e : tmp_edges) {
                if (edge_pred(*e)) {
                    v_filtered_in_edges.push_back(e);
                }
            }

            for (auto &e : v_filtered_in_edges) {
                nb_it = e->getSourceVertex();
                if (vertex_pred(*nb_it)) {
                    VertexType *nb_downcast = dynamic_cast<VertexType *>(&(*nb_it));
                    if (nb_downcast) {
                        EdgeType *e_downcast = dynamic_cast<EdgeType *>(e);
                        /* NOTE: moderately "subtle" lazy evaluation issue here.. */
                        if (e_downcast &&
                            edge_selection_pred(*e_downcast) &&
                            e_downcast->getTraversalState(tid) == this->TRAV_UNSEEN)
                        {
                            reachable_edges.push_back(e_downcast);
                            e_downcast->setTraversalState(tid, this->TRAV_DONE);
                        }
                    }

                    std::list<NeuronVertex *> nb_path = v_path;
                    nb_path.push_back(&(*nb_it));
                    nb_dist = v_dist + e->getLength();
                    if (nb_it->getTraversalState(tid) == this->TRAV_UNSEEN) {
                        Q.push_back(
                            std::tuple<neuron_iterator, std::list<NeuronVertex *>, R>
                                (nb_it, nb_path, nb_dist)
                            );
                    }
                }
            }
        }

        /* vertex v is done. */
        v_it->setTraversalState(tid, this->TRAV_DONE);
    }

    debugTabDec();
    debugl(0, "CellNetwork::getConnectedComponent(): %s traversal from node %2d finished.\n",
        (directed) ? "directed" : "undirected", vstart_it->id());
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
std::pair<uint32_t, R>
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuriteSubTreeDepths(neurite_iterator n_it)
{
    /* compute depth from directed traversal of sub-tree rooted in n_it. first, retrieve all reachable neurite leafs
     * along with paths to them. then the depth is given as the number of edges in the longest returned path. */
    auto neurite_vertex_pred    = std::function<bool(NeuronVertex const &v)>(
            [] (NeuronVertex const &v) -> bool
            {
                return (v.getType() == SOMA_VERTEX ||
                        v.getType() == AXON_VERTEX ||
                        v.getType() == DENDRITE_VERTEX);
            }
        );

    auto neurite_segment_pred   = std::function<bool(NeuronEdge const &e)>(
            [] (NeuronEdge const &e) -> bool
            {
                return (e.getType() == AXON_SEGMENT ||
                        e.getType() == DENDRITE_SEGMENT);
            }
        );

    auto neurite_leaf_pred      = std::function<bool(NeuriteVertex const &v)>(
            [] (NeuriteVertex const &v) -> bool
            {
                return (v.template getFilteredOutNeighbours<NeuriteVertex>().size() == 0);
            }
        );

    auto vertex_false_pred      = std::function<bool(NeuriteVertex const &v)>(
            [] (NeuriteVertex const &v) -> bool
            {
                return false;
            }
        );

    auto edge_false_pred        = std::function<bool(NeuriteSegment const &e)>(
            [] (NeuriteSegment const &e) -> bool
            {
                return false;
            }
        );

    std::list<std::tuple<NeuriteVertex *, std::list<NeuronVertex *>, R>>    reachable_vertices_info;
    std::list<NeuriteSegment *>                                             reachable_edges;

    this->traverseBreadthFirst<NeuriteVertex, NeuriteSegment>(
        /* start traversal from soma s */
        n_it,
        /* directed traversal, return paths */
        true, true,
        /* return lists */
        reachable_vertices_info, reachable_edges,
        /* predicates for neurite veritces / segments */
        neurite_vertex_pred, neurite_segment_pred,
        /* selection predicates for leaf neurite vertices, no edges at all */
        neurite_leaf_pred, edge_false_pred,
        /* trivial termination predicates always returning false */
        vertex_false_pred, edge_false_pred);

    /* get size of longest path to any neurite leaf */
    if (!reachable_vertices_info.empty()) {
        uint32_t    integral_depth  = 0, itmp;
        R           chord_depth     = 0, rtmp;
        for (auto &vi_tuple : reachable_vertices_info) {
            itmp = (std::get<1>(vi_tuple)).size() - 1;
            if (itmp > integral_depth) {
                integral_depth = itmp;
            }

            rtmp = (std::get<2>(vi_tuple));
            if (rtmp > chord_depth) {
                chord_depth = rtmp;
            }
        }
        return std::pair<uint32_t, R>(integral_depth, chord_depth);
    }
    /* no neurite leaf (other than n_it) reachable. depths are zero */
    else {
        return std::pair<uint32_t, R>(0, 0);
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getDirectedPath(
    const neuron_iterator&             u_it,
    const neuron_iterator&             v_it,
    std::list<NeuronVertex *>  &path,
    R                          &pathlen)
{
    auto neuron_vertex_true_pred    = std::function<bool(NeuronVertex const &v)>(
            [] (NeuronVertex const &v) -> bool
            {
                return true;
            }
        );

    auto neuron_vertex_false_pred   = std::function<bool(NeuronVertex const &v)>(
            [] (NeuronVertex const &v) -> bool
            {
                return false;
            }
        );

    auto neuron_edge_true_pred  = std::function<bool(NeuronEdge const &e)>(
            [] (NeuronEdge const &e) -> bool
            {
                return true;
            }
        );
    
    auto neuron_edge_false_pred = std::function<bool(NeuronEdge const &e)>(
            [] (NeuronEdge const &e) -> bool
            {
                return false;
            }
        );


    auto v_pred                     = std::function<bool(NeuronVertex const &v)> (
            [u_it] (NeuronVertex const &v) -> bool
            {
                return (neuron_const_iterator(u_it) == v.iterator());
            }
        );

    std::list<std::tuple<NeuronVertex *, std::list<NeuronVertex *>, R>>     reachable_vertices_info;
    std::list<NeuronEdge *>                                                 reachable_edges;

    this->traverseBreadthFirst<NeuronVertex, NeuronEdge>(
        /* start traversal from u */
        u_it,
        /* directed traversal, return paths */
        true, true,
        /* return lists */
        reachable_vertices_info, reachable_edges,
        /* use all neuron edges / vertices */
        neuron_vertex_true_pred, neuron_edge_true_pred,
        /* selection predicates: get v and no edges */
        v_pred, neuron_edge_false_pred,
        /* trivial termination predicates always returning false */
        neuron_vertex_false_pred, neuron_edge_false_pred);

    /* if v has been reached, copy info */
    if (reachable_vertices_info.size() == 1) {
    }
    /* if v is not reachable, indicate with empty path and distance infinity */
    else if (reachable_vertices_info.empty()) {
        path    = {};
        pathlen = Aux::Numbers::inf<R>(); 
    }
    /* throw: ambiguous, at least two paths from u to v found */
    else {
        throw("CellNetwork::getDirectedPath(): ambiguous result. more than two paths from u to v found.");
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuriteConnectedComponent(
    neurite_iterator                v_it,
    std::list<NeuriteVertex *>     &reachable_neurite_vertices,
    std::list<NeuriteSegment *>    &reachable_neurite_segments)
{
    auto neurite_vertex_pred        = std::function<bool(NeuronVertex const &v)>(
            [] (NeuronVertex const &v) -> bool
            {
                return (v.getType() == AXON_VERTEX || v.getType() == DENDRITE_VERTEX);
            }
        );

    auto neurite_segment_pred       = std::function<bool(NeuronEdge const &e)>(
            [] (NeuronEdge const &e) -> bool
            {
                return (e.getType() == AXON_SEGMENT || e.getType() == DENDRITE_SEGMENT);
            }
        );

    auto neurite_segment_true_pred  = std::function<bool(NeuriteSegment const &e)>(
            [] (NeuriteSegment const &e) -> bool
            {
                return true;
            }
        );

    auto neurite_segment_false_pred = std::function<bool(NeuriteSegment const &e)>(
            [] (NeuriteSegment const &e) -> bool
            {
                return false;
            }
        );

    auto neurite_vertex_true_pred   = std::function<bool(NeuriteVertex const &v)>(
            [] (NeuriteVertex const &v) -> bool
            {
                return true;
            }
        );

    auto neurite_vertex_false_pred  = std::function<bool(NeuriteVertex const &v)>(
            [] (NeuriteVertex const &v) -> bool
            {
                return false;
            }
        );
    

    std::list<std::tuple<NeuriteVertex *, std::list<NeuronVertex *>, R>>    reachable_vertices_info;
    std::list<NeuriteSegment *>                                             reachable_edges;

    this->traverseBreadthFirst<NeuriteVertex, NeuriteSegment>(
        /* start traversal from u */
        v_it,
        /* undirected traversal, don't return paths */
        false, false,
        /* return lists */
        reachable_vertices_info, reachable_edges,
        /* use all neurite vertices / segments */
        neurite_vertex_pred, neurite_segment_pred,
        /* selection predicates: get all neurite vertices and edge */
        neurite_vertex_true_pred, neurite_segment_true_pred,
        /* trivial termination predicates always returning false */
        neurite_vertex_false_pred, neurite_segment_false_pred);

    /* assemble list of all reachable neurite vertices */
    reachable_neurite_vertices.clear();
    for (auto &v_info : reachable_vertices_info) {
        reachable_neurite_vertices.push_back(std::get<0>(v_info));
    }

    /* copy reachable edges to return list */
    reachable_neurite_segments = reachable_edges;
}

/* statistical methods */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
typename CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::CellNetworkStatistics
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
computeStatistics() const
{
    CellNetworkStatistics stat;

    stat.neurite_tree_total_length      = this->getTotalNeuriteTreeLength();
    //stat.morphological_diameter         = this->getMorphologicalDiameter();

    this->getNeuronVerticesDistanceStat(
        stat.neuron_vertices_dist_min,
        stat.neuron_vertices_dist_max,
        stat.neuron_vertices_dist_avg,
        stat.neuron_vertices_dist_sigma);

    this->getNeuronVerticesCoordinateStat(
        stat.neuron_vertices_xmin,
        stat.neuron_vertices_xmax,
        stat.neuron_vertices_xavg,
        stat.neuron_vertices_xsigma,
        //
        stat.neuron_vertices_ymin,
        stat.neuron_vertices_ymax,
        stat.neuron_vertices_yavg,
        stat.neuron_vertices_ysigma,
        //
        stat.neuron_vertices_zmin,
        stat.neuron_vertices_zmax,
        stat.neuron_vertices_zavg,
        stat.neuron_vertices_zsigma);

    this->getNeuriteRadiusStat(
        stat.neurite_radius_min,
        stat.neurite_radius_max,
        stat.neurite_radius_avg,
        stat.neurite_radius_sigma);

    this->getNeuriteSegmentLengthStat(
        stat.neurite_segment_length_min,
        stat.neurite_segment_length_max,
        stat.neurite_segment_length_avg,
        stat.neurite_segment_length_sigma,
        //
        stat.neurite_segment_reduced_length_min,
        stat.neurite_segment_reduced_length_max,
        stat.neurite_segment_reduced_length_avg,
        stat.neurite_segment_reduced_length_sigma);

    return stat;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getTotalNeuriteTreeLength() const
{
    R len = 0;
    for (auto &ns : this->neurite_segments) {
        len += ns.getLength();
    }
    return len;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
R
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getMorphologicalDiameter(soma_iterator s_it)
{
    debugl(0, "CellNetwork::getMorphologicalDiameter(): called on soma s with id: %d\n", s_it->id());
    debugTabInc();

    /* NOTE: this assumes that the CellNetwork::checkTopology() has been called to ensure that (this) CellNetwork forms
     * a forest of properly directed trees when restricted to soma and neurite vertices (cell vertices) as well as
     * neurite root and neurite segment edges (cell edges).
     *
     * with this restriction in place, the morphological diameter of the cell with the given soma is simply the diameter
     * of the soma / neurite vertex sub-tree, i.e. the longest path between any two vertices in said tree. */

    /* find diameter of cell rooted in soma referred to by soma_it */
    auto cell_vertex_pred   = std::function<bool(NeuronVertex const &v)>(
            [] (NeuronVertex const &v) -> bool
            {
                return (v.getType() == SOMA_VERTEX ||
                        v.getType() == AXON_VERTEX ||
                        v.getType() == DENDRITE_VERTEX);
            }
        );

    auto cell_edge_pred     = std::function<bool(NeuronEdge const &e)>(
            [] (NeuronEdge const &e) -> bool
            {
                return (e.getType() == AXON_ROOT_EDGE ||
                        e.getType() == DENDRITE_ROOT_EDGE ||
                        e.getType() == AXON_SEGMENT ||
                        e.getType() == DENDRITE_SEGMENT);
            }
        );

    auto neuron_vertex_true_pred    = std::function<bool(NeuronVertex const &v)>(
            [] (NeuronVertex const &v) -> bool
            {
                return true;
            }
        );

    auto neuron_vertex_false_pred   = std::function<bool(NeuronVertex const &v)>(
            [] (NeuronVertex const &v) -> bool
            {
                return false;
            }
        );

    auto neuron_edge_true_pred  = std::function<bool(NeuronEdge const &e)>(
            [] (NeuronEdge const &e) -> bool
            {
                return true;
            }
        );
    
    auto neuron_edge_false_pred = std::function<bool(NeuronEdge const &e)>(
            [] (NeuronEdge const &e) -> bool
            {
                return false;
            }
        );

    std::list<std::tuple<NeuronVertex *, std::list<NeuronVertex *>, R>> reachable_vertices_info;
    std::list<NeuronEdge *>                                             reachable_edges;

    this->traverseBreadthFirst<NeuronVertex, NeuronEdge>(
        /* start traversal from soma s */
        s_it,
        /* undirected traversal, don't return paths */
        false, false,
        /* addresses of return lists */
        reachable_vertices_info, reachable_edges,
        /* predicates for cell veritces / edges */
        cell_vertex_pred, cell_edge_pred,
        /* trivial selection predicates always returning true */
        neuron_vertex_true_pred, neuron_edge_true_pred,
        /* trivial termination predicates always returning false */
        neuron_vertex_false_pred, neuron_edge_false_pred);

    debugl(0, "reachable_vertices_info.size(): %zu.\n", reachable_vertices_info.size());
    debugl(0, "reachable_edges.size(): %zu.\n", reachable_edges.size());

    if (!reachable_vertices_info.empty()) {
        /* find vertex u of maximum distance to s */
        R               su_dist = -Aux::Numbers::inf<R>();
        neuron_iterator u_it;

        for (auto &vi_tuple : reachable_vertices_info) {
            if (std::get<2>(vi_tuple) > su_dist) {
                su_dist = std::get<2>(vi_tuple);
                u_it    = std::get<0>(vi_tuple)->iterator();
            }
        }
        debugl(0, "neuron vertex of maximum distance to soma s (%d): vertex u (%d) with d(s,u) %10.5f\n", s_it->id(), u_it->id(), su_dist);

        /* start traversal from u, find vertex v of maximum distance to u. d(u, v) is the diameter of the cell tree, i.e.
         * the morphological diameter of the cell rooted in s */
        reachable_vertices_info.clear();
        reachable_edges.clear();
        this->traverseBreadthFirst<NeuronVertex, NeuronEdge>(
            /* start traversal from soma iterator sit */
            u_it,
            /* undirected traversal, don't return paths */
            false, false,
            /* addresses of return lists */
            reachable_vertices_info, reachable_edges,
            /* predicates for cell veritces / edges */
            cell_vertex_pred, cell_edge_pred,
            /* trivial selection predicates always returning true */
            neuron_vertex_true_pred, neuron_edge_true_pred,
            /* trivial termination predicates always returning false */
            neuron_vertex_false_pred, neuron_edge_false_pred);

        /* find vertex u of maximum distance to s */
        R               uv_dist = -Aux::Numbers::inf<R>();
        neuron_iterator v_it;

        for (auto &vi_tuple : reachable_vertices_info) {
            if (std::get<2>(vi_tuple) > uv_dist) {
                uv_dist = std::get<2>(vi_tuple);
                v_it    = std::get<0>(vi_tuple)->iterator();
            }
        }
        debugl(0, "neuron vertex of maximum distance to u (%d): vertex v(%d) with d(u, v) = %10.5f\n", u_it->id(), v_it->id(), uv_dist);
        debugl(0, "<=> diameter of cell tree, i.e. the morphological diameter, is %10.5f\n", uv_dist);

        debugTabDec();
        debugl(0, "CellNetwork::getMorphologicalDiameter(): done.\n");

        return uv_dist;
    }
    /* no neuron vertex reachable from s => return 0 */
    else {
        return 0;
    }
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuronVerticesDistanceStat(
    R  &v_dist_min,
    R  &v_dist_max,
    R  &v_dist_avg,
    R  &v_dist_sigma) const
{
    std::vector<R>  distances;
    R               d;

    for (auto vit = this->neuron_vertices.begin(); vit != this->neuron_vertices.end(); ++vit) {
        auto vit_next = vit;
        ++vit_next;
        for (auto wit = vit_next; wit != this->neuron_vertices.end(); ++wit) {
            d = (wit->getSinglePointPosition() - vit->getSinglePointPosition()).len2();
            distances.push_back(d);
        }
    }

    Aux::Stat::computeMinMaxAvgSigma(
        distances,
        v_dist_min,
        v_dist_max,
        v_dist_avg,
        v_dist_sigma);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuronVerticesCoordinateStat(
    R  &v_xmin,
    R  &v_xmax,
    R  &v_xavg,
    R  &v_xsigma,
    //
    R  &v_ymin,
    R  &v_ymax,
    R  &v_yavg,
    R  &v_ysigma,
    //
    R  &v_zmin,
    R  &v_zmax,
    R  &v_zavg,
    R  &v_zsigma) const
{
    std::vector<R>  x_components, y_components, z_components;
    Vec3<R>         p;

    for (auto &v : this->neuron_vertices) {
        p = v.getSinglePointPosition();
        x_components.push_back(p[0]);
        y_components.push_back(p[1]);
        z_components.push_back(p[2]);
    }

    Aux::Stat::computeMinMaxAvgSigma(
        x_components,
        v_xmin,
        v_xmax,
        v_xavg,
        v_xsigma);

    Aux::Stat::computeMinMaxAvgSigma(
        y_components,
        v_ymin,
        v_ymax,
        v_yavg,
        v_ysigma);

    Aux::Stat::computeMinMaxAvgSigma(
        z_components,
        v_zmin,
        v_zmax,
        v_zavg,
        v_zsigma);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuriteRadiusStat(
    R  &n_rmin,
    R  &n_rmax,
    R  &n_ravg,
    R  &n_rsigma) const
{
    std::vector<R>  values;
    for (auto &nv : this->neurite_vertices) {
        values.push_back(nv.getRadius());
    }
    Aux::Stat::computeMinMaxAvgSigma(
        values,
        n_rmin,
        n_rmax,
        n_ravg,
        n_rsigma);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuriteSegmentRadiusStat(
    R  &ns_abs_radius_diff_min,
    R  &ns_abs_radius_diff_max,
    R  &ns_abs_radius_diff_avg,
    R  &ns_abs_radius_diff_sigma,
    //
    R  &ns_radius_ratio_min,
    R  &ns_radius_ratio_max,
    R  &ns_radius_ratio_avg,
    R  &ns_radius_ratio_sigma) const
{
    std::vector<R>  abs_radius_diff_values;
    std::vector<R>  radius_ratio_values;
    for (auto &ns : this->neurite_segments) {
        abs_radius_diff_values.push_back(ns.getAbsoluteRadiusDifference());
        radius_ratio_values.push_back(ns.getRadiusRatio());

    }

    Aux::Stat::computeMinMaxAvgSigma(
        abs_radius_diff_values,
        ns_abs_radius_diff_min,
        ns_abs_radius_diff_max,
        ns_abs_radius_diff_avg,
        ns_abs_radius_diff_sigma);

    Aux::Stat::computeMinMaxAvgSigma(
        radius_ratio_values,
        ns_radius_ratio_min,
        ns_radius_ratio_max,
        ns_radius_ratio_avg,
        ns_radius_ratio_sigma);
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuriteSegmentLengthStat(
    R  &ns_len_min,
    R  &ns_len_max,
    R  &ns_len_avg,
    R  &ns_len_sigma,
    //
    R  &ns_rlen_min,
    R  &ns_rlen_max,
    R  &ns_rlen_avg,
    R  &ns_rlen_sigma) const
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
getNeuriteSegmentAngleStat(
    R  &ns_angle_min,
    R  &ns_angle_max,
    R  &ns_angle_avg,
    R  &ns_angle_sigma) const
{
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
transformToCentroidSystem()
{
    /* compute centroid of all vertex positions */
    Vec3<R> centroid = Aux::VecMat::nullvec<R>();
    for (auto &v : this->neuron_vertices) {
        centroid += v.getSinglePointPosition();
    }
    centroid /= (R)(this->vertices.size());

    printf("centroid = (%f, %f, %f)\n", centroid[0], centroid[1], centroid[2]);

    /* apply to all vertices */
    for (auto &v : this->neuron_vertices) {
        for (auto &cs : v.sections) {
            cs.position() -= centroid;
        }
    }

    /* update somas individually, since these store the original "graph" used to encode the soma in the input files */
    for (auto &s : this->soma_vertices) {
        for (auto &v : s.section_graph.vertices) {
            v.vertex_data.position() -= centroid;
        }
    }

    /* update global displacement vector */
    this->global_coordinate_displacement += centroid;
}

template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
transformToSomaSystem(soma_const_iterator const &s_it)
{
    /* compute centroid of all vertex positions */
    Vec3<R> d = s_it->getSinglePointPosition();

    debugl(0, "soma position (displacement vector) : (%5.4f, %5.4f, %5.4f)\n", d[0], d[1], d[2]);

    /* apply displacement to all vertices */
    for (auto &v : this->neuron_vertices) {
        for (auto &cs : v.sections) {
            cs.position() -= d;
        }
    }

    /* update somas individually, since these store the original "graph" used to encode the soma in the input files */
    for (auto &s : this->soma_vertices) {
        for (auto &v : s.section_graph.vertices) {
            v.vertex_data.position() -= d;
        }
    }

    /* update global displacement vector */
    this->global_coordinate_displacement += d;
}

/* initialize network from standardized SWC file as used by NeuroMorpho.org */
template <
    typename Tn, typename Tv, typename Te, typename Tso, typename Tnv, typename Tax, typename Tde,
    typename Tns, typename Tas, typename Tds, typename Tnr, typename Tar, typename Tdr, typename R
>
void
CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>::
readFromNeuroMorphoSWCFile(
        std::string     filename,
        bool const     &check_coincident_positions)
{
    using Common::UnitType;
    using Aux::Alg::listContains;

    debugl(0, "CellNetwork::readFromNeuroMorphoSWCFile():\n");
    debugTabInc();

    /* define sorting function as lambda */
    auto sort_id = [] (SWCNode const &a, SWCNode const &b) -> bool {
        return (a.compartment_id < b.compartment_id);
    };

    /* try to open input file */
    std::ifstream           f;
    std::string             line_string;
    char                   *line;
    char                   *linestart;
    std::vector<SWCNode>    swc_nodes;
    std::vector<bool>       node_traversed;
    Vec3<R>                 p;
    R                       r;
    uint32_t                compartment_id, compartment_type;
    int32_t                 parent_id;
    uint32_t                nmatch;

    debugl(1, "trying to open input file \"%s\". \n", filename.c_str() );

    f.open(filename, std::ifstream::in);
    if ( !f.is_open() ) {
        //printf("CellNetwork::readFromNeuroMorphoSWCFile(): cannot open SWC file \"%s\" for reading.\n", filename.c_str() );
        throw("CellNetwork::readFromNeuroMorphoSWCFile(): unable to open SWC file for reading.");
    }
    debugl(1, "successfully opened. parsing input file one line at a time..\n");

    debugTabInc();
    while ( f.good() ) {
        std::getline(f, line_string);
        if (line_string.size() == 0) {
            debugl(3, "empty line..\n");
            continue;
        }

        line        = new char [ line_string.size() + 1];
        linestart   = line;
        std::strcpy(line, line_string.c_str() );

        while (isspace(*linestart)) {
            ++linestart;
        }

        /* skip comments */
        if (*linestart == '#') {
            debugl(3, "comment found.. ignoring.\n");
            delete[] line;
            continue;
        }
        else {
            debugl(3, "line: \"%s\".\n", line);
        }

        if (std::is_same<double, R>::value) {
            nmatch = sscanf(linestart, "%u %u %lf %lf %lf %lf %d\n", &compartment_id, &compartment_type, &p[0], &p[1], &p[2], &r, &parent_id);
        }
        if (std::is_same<float, R>::value) {
            nmatch = sscanf(linestart, "%u %u %f %f %f %f %d\n", &compartment_id, &compartment_type, &p[0], &p[1], &p[2], &r, &parent_id);
        }
        else {
            nmatch = sscanf(linestart, "%u %u %lf %lf %lf %lf %d\n", &compartment_id, &compartment_type, &p[0], &p[1], &p[2], &r, &parent_id);
        }

        if (nmatch != 7) {
            throw("CellNetwork::readFromNeuroMorphoSWCFile(): line is not a comment, yet necessary information could not be matched. syntax error..");
        }

        /* decrement indices, we count fromt 0 */
        compartment_id--;
        /* if parent_id is not -1 (soma), decrement as well */
        if (parent_id != -1) {
            parent_id--;
        }
        //printf("(%u,%u,%f,%f,%f,%f,%d)\n", section_idx, section_type, x, y, z, r, parent_idx);

        swc_nodes.push_back(
                SWCNode(
                    compartment_id,
                    compartment_type,
                    parent_id,
                    { { compartment_id, p, r } }
                )
            );

        delete [] line;
    }
    debugTabDec();
    debugl(1, "input file parsed into SWCNode array.\n");
    
    /* sort nodes by compartment id. */
    std::sort(swc_nodes.begin(), swc_nodes.end(), sort_id);

    /* resize and init traversal flag array */
    node_traversed.assign(swc_nodes.size(), false);

    debugl(1, "checking for consecutive ids.\n");
    /* check whether all ids are consecutive and hence also pairwise distinct */
    for (uint32_t i = 0; i < swc_nodes.size(); i++) {
        if (swc_nodes[i].compartment_id != i) {
            throw("CellNetwork::readFromNeuroMorphoSWCFile(): compartment ids not consecutive. semantic error.\n");
        }
    }

    debugl(1, "computing preliminary tree structure.\n");
    /* compute tree structure among nodes and get all soma nodes, that is: all nodes with parent_id -1. other nodes may
     * refer encode information for the same soma and will be merged below */
    // FIXME: The assertion "soma nodes are nodes with parent_id -1" is not true in general.
    //        Root nodes (parent_id == -1) are not required to be soma nodes (and vice-versa).
    //        I get a hard-to-track RTE for a geometry where this is the root node is a dendritic node.
    std::list<uint32_t> soma_root_node_ids;
    for (auto &node : swc_nodes) {
        /* get parent of node, add node as child of parent */
        if (node.parent_id == -1) {
            soma_root_node_ids.push_back( node.compartment_id );
        }
        else {
            swc_nodes[node.parent_id].child_ids.push_back(node.compartment_id);
        }
    }

    /* sort all child_id lists. */
    for (auto &node : swc_nodes) {
        node.child_ids.sort();
    }

    debugl(1, "gathering soma information from sub-graphs of soma info nodes.\n");

    /* for all soma nodes, which should be the roots of all subtrees (topology will be checked in detail below), get all
     * descendants whose compartment type indicates soma as well. these nodes (including the root soma node) shall be
     * referred to as soma info nodes. in the order of appearance of a depth-first traversal, encode their sections in a
     * graph (soma_sectiongraph), which is then stored in the SomaVertex in (this) CellNetwork. this graph exactly
     * mirrors the data contained in the SWC file, ambiguous and ill-formed as it may be in many cases.
     *
     * remember ids of all soma info nodes excluding the root soma node, mark all of them as traversed and disconnect
     * them from the soma root node by erasing the respective ids from the soma root node's children list. */
    std::list<uint32_t>                             soma_info_node_ids;

    /* map associating soma root node ids to soma graphs and a matching iterator */
    std::map<
        uint32_t,
        Graph<uint32_t, CellSection, UnitType>
    >                                               soma_sectiongraph_map;
    typename std::map<
        uint32_t,
        Graph<uint32_t, CellSection, UnitType>
    >::iterator                                     sgm_it;

    /* stack containing pairs (id, Graph::vertex_iterator) where the id refers to swc_nodes array and graph_iterator is
     * the iterator to the respective vertex in soma_sectiongraph */
    std::list<
            std::pair<
                uint32_t,
                typename Graph<uint32_t, CellSection, UnitType>::vertex_iterator
            >
        >                                           S;

    /* traverse connected components of all soma root nodes */
    debugTabInc();
    for (auto &sid : soma_root_node_ids) {
        debugl(2, "processing soma root node %2d.\n", sid);
        SWCNode &soma_root_node = swc_nodes[sid];

        /* reset data */
        soma_info_node_ids.clear();
        
        /* create new soma graph for soma_root_node with id sid */
        auto sgmpair = soma_sectiongraph_map.insert({ sid, Graph<uint32_t, CellSection, UnitType>() });

        if (!sgmpair.second) {
            throw("CellNetwork::readFromNeuroMorphoSWCFile(): failed to insert fresh graph into soma graph map for "\
                "current root soma node.");
        }

        /* get reference to fresh soma graph from map iterator */
        Graph<uint32_t, CellSection, UnitType> &soma_sectiongraph = (sgmpair.first)->second;
        
        /* insert soma vertex into fresh soma section graph */
        auto somagraph_sit = soma_sectiongraph.vertices.insert(soma_root_node.sections.front());

        /* store id of inserted vertex representing soma root node in the graph's data. this is used to quickly identify
         * the root (or source) in soma_sectiongraph */
        soma_sectiongraph.data() = somagraph_sit->id();

        debugl(2, "initializing traversal stack with all direct soma info node children of current soma root node %d.\n", soma_root_node.compartment_id);
        /* init stack to contain pairs (id, soma_sectiongraph iterator) for all non-root soma info nodes directly
         * connected to the root soma node. */
        S.clear();
        debugTabInc();
        for (auto &child_id : soma_root_node.child_ids) {
            SWCNode &c = swc_nodes[child_id];
            if (c.compartment_type == SOMA_COMPARTMENT) {
                /* insert vertex for c and edge (soma_root_node, c) into soma_sectiongraph */
                auto somagraph_c_it = soma_sectiongraph.vertices.insert(c.sections.front());
                auto rpair          = soma_sectiongraph.edges.insert(somagraph_sit, somagraph_c_it);
                if (rpair.second) {
                    debugl(2, "adding soma info node %d to _bottom_ of stack.\n", c.compartment_id);
                    S.push_front( { c.compartment_id, somagraph_c_it } );
                    soma_info_node_ids.push_back(c.compartment_id);
                }
                else {
                    throw("CellNetwork::readFromNeuroMorphoSWCFile(): failed to add edge (r, c) into soma info "\
                        "graph for soma root node r and its child c. duplicate edge in SWC file?");
                }
            }
        }
        debugTabDec();
        debugl(2, "traversal stack initialized. starting depth-first traversal..\n");

        /* stack has been initialized => depth-first traversal of all soma info node sub-trees for current
         * soma_root_node. */
        debugTabInc();
        while ( !S.empty()) {
            /* get top() soma info node id / soma section graph iterator */
            SWCNode &v          = swc_nodes[S.back().first];
            auto somagraph_v_it = S.back().second;

            S.pop_back();

            debugl(2, "stack top() current soma info node v: %2d. processing children..\n", v.compartment_id);

            /* scan through children of v */
            debugTabInc();
            for (auto child_it = v.child_ids.begin(); child_it != v.child_ids.end(); ) {
                SWCNode c = swc_nodes[*child_it];

                debugl(2, "current child c = %2d of v = %2d.\n", c.compartment_id, v.compartment_id);

                /* if child is a soma info node, handle it */
                debugTabInc();
                if (c.compartment_type == SOMA_COMPARTMENT) {
                    debugl(2, "child is SOMA compartment => adding child c and edge (v,c) to soma graph.\n");

                    /* add vertex c and edge (v, c) to soma section graph */
                    auto somagraph_c_it     = soma_sectiongraph.vertices.insert(c.sections.front());
                    auto somagraph_vc_eit   = soma_sectiongraph.edges.insert(somagraph_v_it, somagraph_c_it);

                    /* although it might be reasonable to assume that the soma info data is a tree, there seems to be
                     * little regularity in the data from NeuroMorpho.org. therefore check if child has already been
                     * traversed before pushing onto stack. */
                    if (!node_traversed[c.compartment_id]) {
                        debugl(2, "child c not yet traversed. pushing onto stack..\n");
                        S.push_back( { c.compartment_id, somagraph_c_it } );

                        /* append child node to soma info node list */
                        soma_info_node_ids.push_back(c.compartment_id);
                    }

                    ++child_it;
                }
                /* redirect all children with compartment type other than soma to the soma root node. */
                else {
                    debugl(2, "child is non-SOMA compartment => redirecting edge (v = %d, c = %d) to soma root node %d as edge (%d, c = %d).\n",
                        v.compartment_id, c.compartment_id, soma_root_node.compartment_id, soma_root_node.compartment_id, c.compartment_id);

                    soma_root_node.child_ids.push_back(*child_it); 
                    child_it = v.child_ids.erase(child_it);
                }
                debugTabDec();
            }
            debugTabDec();

            /* all v's children have been processed (redirected to root soma node, skipped or pushed onto the stack in
             * case of soma info nodes). mark v as traversed. */
            node_traversed[v.compartment_id] = true;
        }
        debugTabDec();

        debugl(2, "marking all non-root soma info nodes as traversed.\n");
        /* all sub-trees of soma info vertices attached to root soma_node have been traversed. mark all info nodes as
         * traversed. */
        for (auto &vid : soma_info_node_ids) {
            node_traversed[vid] = true;
        }

        /* disconnect all non-root soma info nodes from root soma node */
        debugl(2, "removing all (non-root) soma info node children of soma root node.\n");
        debugTabInc();
        for (auto child_it = soma_root_node.child_ids.begin(); child_it != soma_root_node.child_ids.end(); ) {
            if (swc_nodes[*child_it].compartment_type == SOMA_COMPARTMENT) {
                debugl(2,"removing child node %2d of soma root node %2d\n", *child_it, soma_root_node.compartment_id);
                child_it = soma_root_node.child_ids.erase(child_it);
            }
            else ++child_it;
        }
        debugTabDec();
    }
    debugTabDec();

    debugl(1, "all soma root nodes processed. number of somas: %3zu\n", soma_root_node_ids.size());

    /* traverse the preliminary swc "tree" breadth first and insert all vertices / edges into (this) CellNetwork */

    /* queue containing pairs (id, CellNetwork::neuron_iterator) where the id refers to swc_nodes array and
     * neuron_iterator is the iterator to the respective NeuronVertex in (this) CellNetwork */
    std::list<
            std::pair<
                uint32_t,
                neuron_iterator
            >
        >                       Q;

    /* iterators for all possible situations.. */
    neuron_iterator             nit;

    soma_iterator               n_sit;
    axon_iterator               n_ait;
    dendrite_iterator           n_dit;

    soma_iterator               c_sit;
    axon_iterator               c_ait;
    dendrite_iterator           c_dit;

    axon_rootedge_iterator      arit;
    dendrite_rootedge_iterator  drit;
    axon_segment_iterator       asit;
    dendrite_segment_iterator   dsit;

    debugl(1, "traversing preliminary forest breadth-first and constructing CellNetwork..\n");
    debugTabInc();
    for (auto &sid : soma_root_node_ids) {
        debugl(2, "traversing connected component of soma %2d.\n", sid);

        /* clear Q */
        Q.clear();

        /* insert soma root node with section graph retrieved from soma_sectiongraph_map, implicitly convert returned
         * soma_iterator to neuron_iterator when enqueueing pair (sid, soma_iterator) into Q. */
        sgm_it = soma_sectiongraph_map.find(sid);
        if (sgm_it != soma_sectiongraph_map.end()) {
            n_sit = this->soma_vertices.insert(sgm_it->second);
            Q.push_back( {sid, n_sit} );
        }
        else {
            throw("CellNetwork::readFromNeuroMorphoSWCFile(): no soma graph ground for current soma root "\
                "node. internal logic error.");
        }

        /* breadth-first traversal. while queue non-empty, keep going */
        debugTabInc();
        while (!Q.empty()) {
            /* get current front() node */
            SWCNode &n  = swc_nodes[Q.front().first];
            nit         = Q.front().second;
            Q.pop_front();

            debugl(2, "current node n: %d.\n", n.compartment_id);

            /* mark node as traversed */
            node_traversed[n.compartment_id] = true;

            debugl(2, "inspecting all of n's neighbours and adding corresponding edges to CellNetwork..\n");
            /* examine all children of n, that is to say all out-going edge (n, c) */
            debugTabInc();
            for (auto &child_id : n.child_ids) {
                SWCNode &c = swc_nodes[child_id];
                debugl(2, "current neighbour c = %2d of n = %2d. adding c as vertex in CellNetwork..\n",
                    c.compartment_id, n.compartment_id);

                debugTabInc();

                /* explicitly invalidate child iterators */
                c_sit.explicitlyInvalidate();
                c_ait.explicitlyInvalidate();
                c_dit.explicitlyInvalidate();

                debugTabInc();
                /* add child node node c as vertex */
                switch (c.compartment_type) {
                    case SOMA_COMPARTMENT:
                        debugl(2, "c is soma vertex => inserting..\n");
                        /* insert soma vertex, use default values for Tv and Ts, since SWC files do not provide any such
                         * additional information. information must be annotated to vertices later on. same for other vertex
                         * types below.
                         *
                         * for the special case of the soma, we need the section graph computed from all soma info nodes
                         * above. it has been inserted into soma_sectiongraph_map with the root soma node's compartment id
                         * as key => find graph in map and insert soma vertex. */
                        sgm_it = soma_sectiongraph_map.find(c.compartment_id);
                        if (sgm_it != soma_sectiongraph_map.end()) {
                            /* insert soma with graph extracted from soma graph map iterator */
                            c_sit = this->soma_vertices.insert(sgm_it->second);
                        }
                        else {
                            throw("CellNetwork::readFromNeuroMorphoSWCFile(): no soma graph ground for current soma root "\
                                "node. internal logic error.");
                        }
                        break;

                    /* axon node */
                    case AXON_COMPARTMENT:
                        debugl(2, "c is axon => inserting..\n");
                        if (c.sections.size() == 1) {
                            c_ait = this->axon_vertices.insert(c.sections.front());
                        }
                        else {
                            throw("CellNetwork::readFromNeuroMorphoSWCFile(): discovered axon SWC node containing not "\
                                "exactly one CellSection. internal logic error.");
                        }
                        break;

                    /* apical dendrite node */
                    case APICAL_DENDRITE_COMPARTMENT:
                        debugl(2, "c is apical dendrite => inserting..\n");
                        if (c.sections.size() == 1) {
                            c_dit = this->dendrite_vertices.insert(c.sections.front(), true);
                        }
                        else {
                            throw("CellNetwork::readFromNeuroMorphoSWCFile(): discovered basal dendrite SWC node "\
                                "containing not exactly one CellSection. internal logic error.");
                        }
                        break;

                    /* basal dendrite node */
                    case BASAL_DENDRITE_COMPARTMENT:
                        debugl(2, "c is basal dendrite => inserting..\n");
                        if (c.sections.size() == 1) {
                            c_dit = this->dendrite_vertices.insert(c.sections.front(), false);
                        }
                        else {
                            throw("CellNetwork::readFromNeuroMorphoSWCFile(): discovered apical dendrite SWC node "\
                                "containing not exactly one CellSection. internal logic error.");
                        }
                        break;

                    default:
                        debugl(2, "c has unsupported section type. skipping n and its sub-tree.\n");
                        /* unsupported compartment type. skip node with continue, i.e. don't examine out-going edges */
                        printf("CellNetwork::readFromNeuroMorphoSWCFile(): WARNING: SWC node with unsupported compartment"\
                               " type %2u: skipping node's sub-tree.\n", c.compartment_type);
                        continue;
                }
                debugTabDec();

                /* insert edge (n, c) into CellNetwork, the type of which depends on the compartment types of n and c.
                 * "down" convert neuron_iterator for n as needed. */
                debugl(2, "inserting edge (n, c) into CellNetwork..\n");

                debugTabInc();
                /* axon root edge */
                if      (n.compartment_type == SOMA_COMPARTMENT &&
                        c.compartment_type == AXON_COMPARTMENT)
                {
                    debugl(2, "edge (n, c) is axon root edge => adding..\n");

                    /* "down"-convert neuron iterator for n to soma_iterator and insert edge.. analogous for other types
                     * below. */
                    n_sit = soma_iterator(nit.network, nit.int_it);
                    this->axon_root_edges.insert(n_sit, c_ait);
                }
                /* apical dendrite root edge */
                else if (n.compartment_type == SOMA_COMPARTMENT &&
                        c.compartment_type == APICAL_DENDRITE_COMPARTMENT)
                {
                    debugl(2, "edge (n, c) is apical dendrite root edge => adding..\n");

                    n_sit = soma_iterator(nit.network, nit.int_it);
                    this->dendrite_root_edges.insert(n_sit, c_dit);
                }
                /* basal dendrite root edge */
                else if (n.compartment_type == SOMA_COMPARTMENT &&
                        c.compartment_type == BASAL_DENDRITE_COMPARTMENT)
                {
                    debugl(2, "edge (n, c) is basal dendrite root edge => adding..\n");

                    n_sit = soma_iterator(nit.network, nit.int_it);
                    this->dendrite_root_edges.insert(n_sit, c_dit);
                }
                /* axon segment */
                else if (n.compartment_type == AXON_COMPARTMENT &&
                        c.compartment_type == AXON_COMPARTMENT)
                {
                    debugl(2, "edge (n, c) is axon segment => adding..\n");

                    n_ait = axon_iterator(nit.network, nit.int_it);
                    this->axon_segments.insert(n_ait, c_ait);
                }
                /* apical dendrite segment */
                else if (n.compartment_type == APICAL_DENDRITE_COMPARTMENT &&
                        c.compartment_type == APICAL_DENDRITE_COMPARTMENT)
                {
                    debugl(2, "edge (n, c) is apical dendrite segment => adding..\n");

                    n_dit = dendrite_iterator(nit.network, nit.int_it);
                    this->dendrite_segments.insert(n_dit, c_dit);
                }
                /* basal dendrite segment */
                else if (n.compartment_type == BASAL_DENDRITE_COMPARTMENT &&
                        c.compartment_type == BASAL_DENDRITE_COMPARTMENT)
                {
                    debugl(2, "edge (n, c) is basal dendrite segment => adding..\n");

                    n_dit = dendrite_iterator(nit.network, nit.int_it);
                    this->dendrite_segments.insert(n_dit, c_dit);
                }
                /* (yet) unsupported edge type (not contained in NeuroMorpho SWC files as of 01/2014). throw exception.. */
                else {
                    debugl(2, "edge (n, c) is unsupported. throwing..\n");
                    throw("CellNetwork::readFromNeuroMorphoSWCFile(): discovered unsupported edge type during "\
                        "traversal of preliminary swc tree and construction of CellNetwork.");
                }

                /* in SWC files from NeuroMorpho, the encoded graph is really a tree and thus an already traversed node
                 * should never be found here. for file types containing synapses, this is not necessarily true.
                 * following the general scheme: enqueue child only if it has not been traversed yet. issue warning: */
                if (!node_traversed[child_id]) {
                    debugl(2, "current neighbour c = %2d of n = %2d. enqueueing c.. \n", c.compartment_id, n.compartment_id);

                    /* enqueue pair (c.compartment_id, "up"-converted neuron_iterator) in Q. the iterator to be used
                     * depends on the compartment type of c, as above when n was "down"-converted */
                    switch (c.compartment_type) {
                        case SOMA_COMPARTMENT:
                            Q.push_back( { child_id, c_sit } );
                            break;

                        case AXON_COMPARTMENT:
                            Q.push_back( { child_id, c_ait } );
                            break;

                        case BASAL_DENDRITE_COMPARTMENT:
                        case APICAL_DENDRITE_COMPARTMENT:
                            Q.push_back( { child_id, c_dit } );
                            break;

                        default:
                            throw("CellNetwork::readFromNeuroMorphoSWCFile(): unsupported compartment type discovered "\
                                "while about to enqueue child. this must never happen, since the node should have "\
                                "already been skipped above. internal logic error.");
                    };
                }
                else {
                    printf("CellNetwork::readFromNeuroMorphoSWCFile(): WARNING: discovered already traversed node in "\
                        "neighbours of current node in bread-first traversal. should not happen for trees.\n");
                }
                debugTabDec();

                debugTabDec();
            }
            debugTabDec();
        }
        debugTabDec();
    }
    debugTabDec();
    debugl(1, "tree (or a subset of it in case of unsupported compartment types) converted to CellNetwork. "\
        "NOTE: topology remains to be checked.\n");

    debugl(1, "checking whether all nodes of supported compartment type have been traversed.\n");
    std::list<uint32_t> supported_compartment_types = {
            SOMA_COMPARTMENT,
            AXON_COMPARTMENT,
            APICAL_DENDRITE_COMPARTMENT,
            BASAL_DENDRITE_COMPARTMENT
        };

    /* check if all nodes of supported compartment type have been traversed. */
    for (uint32_t i = 0; i < swc_nodes.size(); i++) {
        if (listContains<uint32_t>(supported_compartment_types, swc_nodes[i].compartment_type)) {
            if (node_traversed[i] == false) {
                throw("CellNetwork::readFromNeuroMorphoSWCFile(): discovered un-traversed ndoe with supported "\
                    " compartment type after CellNetwork construction.\n");
            }
        }
    }

    /* check for coincident positions in all cell sections. */
    bool coincident_positions = false;
    if (check_coincident_positions) {
        debugl(1, "checking for coincident vertices.\n");
        for (auto vit = this->neuron_vertices.begin(); vit != this->neuron_vertices.end(); ++vit) {
            auto vit_next = vit;
            ++vit_next;
            for (auto wit = vit_next; wit != this->neuron_vertices.end(); ++wit) {
                for (auto &vs : vit->sections) {
                    for (auto &ws : wit->sections) {
                        Vec3<R> pvs = vs.position();
                        Vec3<R> pws = ws.position();

                        if ( (pvs - pws).len2() < 1E-8) {
                            coincident_positions = true;
                            printf("CellNetwork::readFromNeuroMorphoSWCFile(): compartments %5d <-> %5d: distance < 1E-8 = 0.01pm => coincident positions.\n",
                                vs.compartment_id() + 1, ws.compartment_id() + 1);
                            throw("CellNetwork::readFromNeuroMorphoSWCFile(): discovered two distinct compartments with coincident positions. semantic error.\n");
                        }
                    }
                }
            }
        }
    }

    debugl(1, "checking for zero radii.\n");
    /* check for zero radii */
    R r_vs;
    bool zero_radii = false;
    for (auto &v : this->neuron_vertices) {
        for (auto &vs : v.sections) {
            r_vs = vs.radius();    
            if (r_vs < 1E-8) {
                zero_radii = true;
                printf("CellNetwork::readFromNeuroMorphoSWCFile(): compartment: %5d. radius %12.5e < 1E-8 = 0.01pm. biologically infeasible.\n",
                    vs.compartment_id() + 1, r_vs);
                throw("CellNetwork::readFromNeuroMorphoSWCFile(): discovered compartment with (close-to) zero radius. biologically infeasible. semantic error.");
            }
        }
    }

    debugl(1, "checking topology of created CellNetwork.\n");
    this->checkTopology();
     
    if (coincident_positions) {
        throw("SWCParser(): cell contains compartments with unrealistically close-to coincident positions. semantic error.\n");
    }
    if (zero_radii) {
        throw("SWCParser(): cell contains compartments with unrealistically close-to zero radii. semantic error.\n");
    }

    /* finally, initialize basic network information (store iterators to soma and neurite in neurite vertices, etc.) */
    this->initializeNetworkInfo();

    debugTabDec();
    debugl(0, "CellNetwork::readFromNeuroMorphoSWCFile(): done.\n");
}
