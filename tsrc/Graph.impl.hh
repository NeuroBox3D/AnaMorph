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

#include "aux.hh"

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        graph vertex class implementation ....                                                        
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::resetTraversalStates()
{
    /* clear id queue for traversal ids, get first id, reset all states with that id */
    this->traversal_idq.clear();
    uint32_t first_id = this->traversal_idq.getId();

    for (auto &v : this->vertices) {
        v.setTraversalState(first_id, TRAV_UNSEEN);
    }

    for (auto &f : this->edges) {
        f.setTraversalState(first_id, TRAV_UNSEEN);
    }
}

template <typename Tg, typename Tv, typename Te>
std::pair<typename std::map<uint32_t, typename Graph<Tg, Tv, Te>::VertexPointerType>::iterator, bool>
Graph<Tg, Tv, Te>::protectedVertexInsert(Vertex *v)
{
    debugl(2, "Graph::protectedVertexInsert()\n");
    debugTabInc();

    /* erase all topological information from v */
    v->in_edges.clear();
    v->out_edges.clear();

    /* init data as in vertex constructor */
    v->graph                = this;
    v->current_traversal_id = 0;
    v->traversal_state      = TRAV_UNSEEN;
    v->blocked              = false;

    /* insert vertex just as in VertexAccessor::insert().. */
    std::pair<typename std::map<uint32_t, VertexPointerType >::iterator, bool> pair, ret;

    /* get fresh id for new vertex, allocate new vertex, insert pair (id, vertex) into map */
    uint32_t v_id   = this->V_idq.getId();
    pair            = this->V.insert( { v_id, VertexPointerType(v) } );
    if (!pair.second) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::protectedVertexInsert(): vertex with fresh id from idq already present in vertex map. this must never happen..");
    }
    else {
        /* the (private) constructor Vertex(Graph * const &, const Tv & const &data) does not set the internal iterator
         * Vertex::g_vit, since it has no way of knowing it in advance. set it now to bring the new
         * vertex into a consistent state. */
        v->g_vit    = pair.first;

        /* prepare and return pair */
        ret.first   = v->g_vit;
        ret.second  = true;

        debugTabDec();
        debugl(2, "Graph::protectedVertexInsert(): done.\n");

        return ret;
    }
}

template <typename Tg, typename Tv, typename Te>
std::pair<typename std::map<uint32_t, typename Graph<Tg, Tv, Te>::EdgePointerType>::iterator, bool>
Graph<Tg, Tv, Te>::protectedEdgeInsert(
    Edge       *e,
    uint32_t    v_src_id,
    uint32_t    v_dst_id)  
{
    debugl(2, "Graph::protectedEdgeInsert()\n");
    debugTabInc();

    /* locate two vertices by id and check if they match with the ones stored in e */
    vertex_iterator v_src_it, v_dst_it;

    v_src_it = this->vertices.find(v_src_id);
    v_dst_it = this->vertices.find(v_dst_id);

    if (v_src_it != this->vertices.end() && v_dst_it != this->vertices.end()) {
        /* compare pointers */
        if ( v_src_it == e->getSourceVertex() && v_dst_it == e->getDestinationVertex() ) {
            /* pointers match, check if edge already exists via iterators */
            std::pair<typename std::map<uint32_t, EdgePointerType >::iterator, bool> ret;

            if (this->edges.exists(v_src_it, v_dst_it)) {
                ret.first = this->E.end();
                ret.second = false;
            }

            /* insert edge as in EdgeAccessor::insert(). first, init edge data as in Edge::Edge constructor */
            e->graph                = this; 
            e->blocked              = false;
            e->current_traversal_id = 0;
            e->traversal_state      = TRAV_UNSEEN;

            uint32_t    edge_id;
            edge_id     = this->E_idq.getId();
            auto rpair  = this->E.insert( {edge_id, EdgePointerType(e) } );
            if (rpair.second) {
                /* set Edge::m_fit iterator, which is required for Edge to be in a consistent state and has not
                 * been set by the (private) Edge ctor, just as for Graph::Vertex */
                e->g_eit = rpair.first;

                /* topology information update */
                bool all_inserted = 
                    Aux::Alg::listSortedInsert(v_src_it->out_edges, e, false) &&
                    Aux::Alg::listSortedInsert(v_dst_it->in_edges, e, false);

                if (all_inserted) {
                    ret.first   = e->g_eit;
                    ret.second  = true;

                    debugTabDec();
                    debugl(2, "Graph::protectedEdgeInsert(): done.\n");

                    return ret;
                }
                else {
                    throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::protectedEdgeInsert(): at least one (sorted) insertion of edge pointer into {out/in} lists of {v_src/v_dst} failed. internal logic error.");
                }
            }
            else {
                throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::protectedEdgeInsert(): new Edge with fresh id from idq already present in Edge map. this must never happen..");
            }
        }
        else {
            throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::protectedEdgeInsert(): source and destination vertices found by id, yet pointers do not match those stored in given (upcast) Edge object. internal logic error.");
        }
    }
    else {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::protectedEdgeInsert(): at least one of the two given ids was not found. internal logic error.");
    }
}

/* protected (static) methods to extract internal iterators from Graph::{vertex,edge}_iterators */
template <typename Tg, typename Tv, typename Te>
typename std::map<uint32_t, typename Graph<Tg, Tv, Te>::VertexPointerType>::const_iterator
Graph<Tg, Tv, Te>::getInternalIterator(vertex_const_iterator const &it)
{
    return it.int_it;
}

template <typename Tg, typename Tv, typename Te>
typename std::map<uint32_t, typename Graph<Tg, Tv, Te>::VertexPointerType>::iterator
Graph<Tg, Tv, Te>::getInternalIterator(vertex_iterator const &it)
{
    return it.int_it;
}

template <typename Tg, typename Tv, typename Te>
typename std::map<uint32_t, typename Graph<Tg, Tv, Te>::EdgePointerType>::const_iterator
Graph<Tg, Tv, Te>::getInternalIterator(edge_const_iterator const &it)
{
    return it.int_it;
}

template <typename Tg, typename Tv, typename Te>
typename std::map<uint32_t, typename Graph<Tg, Tv, Te>::EdgePointerType>::iterator
Graph<Tg, Tv, Te>::getInternalIterator(edge_iterator const &it)
{
    return it.int_it;
}

/* graph vertex ctors */
/* NOTE: Vertex::g_vit is not set after these constructors have been called and thus the vertex is
 * not in a consistent state. the caller (which has to be a friend of Graph::Vertex, since the
 * following ctor are private) needs to take care to set the iterator correctly */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Vertex::Vertex() : vertex_data()
{
    this->graph                 = NULL;
    this->current_traversal_id  = 0;
    this->traversal_state       = TRAV_UNSEEN;
    this->blocked               = false;
}

template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Vertex::Vertex(
    Graph<Tg, Tv, Te> * const  &graph,
    Tv const                   &data) : vertex_data(data)
{
    this->graph                 = graph;
    this->current_traversal_id  = 0;
    this->blocked               = false;
    this->traversal_state       = TRAV_UNSEEN;
}

/* private copy ctor */
/* NOTE: this constructor is protected (i.e. not publicly accessible) and is only used internally for
 * the implementation of Graph methods.  anway: careful when changing the internals, pointers are
 * copied!  if this is constructor is invoked to change the id of a edge within a graph, that's ok.
 * if however this is used to copy vertices from one graph to another, the caller must take care to
 * adjust the adjacency / incidence information accordingly to avoid UB. */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Vertex::Vertex(const Vertex &x) : vertex_data(x.vertex_data)
{    
    this->graph                 = x.graph;
    this->g_vit                 = x.g_vit;
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->blocked               = x.blocked;
    this->in_edges              = x.in_edges;
    this->out_edges             = x.out_edges;
}

/* private assignment operator. same as for private copy ctor: to be used only internally and only
 * with care. */
template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::Vertex &
Graph<Tg, Tv, Te>::Vertex::operator=(const Vertex &x)
{
    this->graph                 = x.graph;
    this->g_vit                 = x.g_vit;
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->blocked               = x.blocked;
    this->in_edges              = x.in_edges;
    this->out_edges             = x.out_edges;
    this->vertex_data           = x.vertex_data;

    return (*this);
}

/* virtual destructor to allow polymorphic inheritance */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Vertex::~Vertex()
{
}

template <typename Tg, typename Tv, typename Te>
void *
Graph<Tg, Tv, Te>::Vertex::operator new(size_t size)
{
    if (size != sizeof(Graph::Vertex)) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::Vertex::operator new(): size does not match sizeof(Vertex). internal logic error.");
    }
    return ::operator new(size);
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::operator delete(void *p)
{
    ::operator delete(p);
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::replaceAdjacentVertices(const std::map<Vertex *, Vertex*> &replace_map)
{
    typename std::map<Vertex *, Vertex *>::const_iterator    mit;

    /* replace pointers in all out-going and in-coming edges */
    for (auto &e_out : this->out_edges) {
        if ( (mit = replace_map.find(e_out->v_dst)) != replace_map.end() ) {
            e_out->v_dst = mit->second;
        }
    }
    for (auto &e_in : this->in_edges) {
        if ( (mit = replace_map.find(e_in->v_src)) != replace_map.end() ) {
            e_in->v_src = mit->second;
        }
    }
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::Vertex *
Graph<Tg, Tv, Te>::Vertex::getPtr(typename std::map<uint32_t, VertexPointerType >::const_iterator const &it)
{
    return (it->second);
}


/* retrieve out and in neighbour vertices */
template <typename Tg, typename Tv, typename Te>
std::list<const typename Graph<Tg, Tv, Te>::Vertex *>
Graph<Tg, Tv, Te>::Vertex::getOutNeighbours() const
{
    std::list<const Vertex *> ret;
    this->getOutNeighbours(ret);
    return ret;
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getOutNeighbours(std::list<Vertex *> &out_nbs) const
{
    out_nbs.clear();

    for (auto &e : this->out_edges) {
        out_nbs.push_back( &(*e->getDestinationVertex()) );
    }
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getOutNeighbours(std::list<const Vertex *> &out_nbs) const
{
    out_nbs.clear();

    for (auto &e : this->out_edges) {
        out_nbs.push_back( &(*e->getDestinationVertex()) );
    }
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getOutNeighboursIndices(std::list<uint32_t> &out_nbs) const
{
    out_nbs.clear();

    for (auto &e : this->out_edges) {
        out_nbs.push_back( e->getDestinationVertex()->id() );
    }
}


template <typename Tg, typename Tv, typename Te>
std::list<const typename Graph<Tg, Tv, Te>::Vertex *>
Graph<Tg, Tv, Te>::Vertex::getInNeighbours() const
{
    std::list<const Vertex *> ret;
    this->getInNeighbours(ret);
    return ret;
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getInNeighbours(std::list<Vertex *> &in_nbs) const
{
    in_nbs.clear();

    for (auto &e : this->in_edges) {
        in_nbs.push_back( &(*e->getSourceVertex()) );
    }
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getInNeighbours(std::list<const Vertex *> &in_nbs) const
{
    in_nbs.clear();

    for (auto &e : this->in_edges) {
        in_nbs.push_back( &(*e->getSourceVertex()) );
    }
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getInNeighboursIndices(std::list<vertex_iterator> &in_nbs) const
{
    in_nbs.clear();

    for (auto &e : this->in_edges) {
        in_nbs.push_back( e->getDestinationVertex()->iterator() );
    }
}


/* retrieve in and out edges */
template <typename Tg, typename Tv, typename Te>
std::list<const typename Graph<Tg, Tv, Te>::Edge *>
Graph<Tg, Tv, Te>::Vertex::getOutEdges() const
{
    std::list<const Edge *> ret;
    this->getOutEdges(ret);
    return ret;
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getOutEdges(std::list<Edge *> &out_edges) const
{
    out_edges.clear();
    for (auto &e : this->out_edges) {
        out_edges.push_back(e);
    }
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getOutEdges(std::list<const Edge *> &out_edges) const
{
    out_edges.clear();
    for (auto &e : this->out_edges) {
        out_edges.push_back(e);
    }
}


template <typename Tg, typename Tv, typename Te>
std::list<const typename Graph<Tg, Tv, Te>::Edge *>
Graph<Tg, Tv, Te>::Vertex::getInEdges() const
{
    std::list<const Edge *> ret;
    this->getInEdges(ret);
    return ret;
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getInEdges(std::list<Edge *> &in_edges) const
{
    in_edges.clear();
    for (auto &e : this->in_edges) {
        in_edges.push_back(e);
    }
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::getInEdges(std::list<const Edge *> &in_edges) const
{
    in_edges.clear();
    for (auto &e : this->in_edges) {
        in_edges.push_back(e);
    }
}


/* get and set traversal states for given traversal id */
template <typename Tg, typename Tv, typename Te>
uint8_t
Graph<Tg, Tv, Te>::Vertex::getTraversalState(const uint32_t &traversal_id)
{
    if (this->current_traversal_id != traversal_id) {
        this->current_traversal_id  = traversal_id;
        this->traversal_state       = TRAV_UNSEEN;
    }
    return this->traversal_state;
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Vertex::setTraversalState(
    const uint32_t &traversal_id,
    const uint8_t  &state)
{
    this->current_traversal_id  = traversal_id;
    this->traversal_state       = state;
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        graph edge class implementation ....                                                        
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* ctors */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Edge::Edge()
: graph(NULL), v_src(NULL), v_dst(NULL), blocked(false), current_traversal_id(0),
  traversal_state(TRAV_UNSEEN), edge_data()
{}

/* NOTE: private ctor only to be used by the internal implementation, which default constructs the
 * internal iterator Graph::g_eit. the (internal) caller must set the iterator (once it is known) to
 * bring the Edge into a consistent state. */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Edge::Edge
(
    Graph<Tg, Tv, Te>* _graph,
    Vertex*            _v_src,
    Vertex*            _v_dst,
    const Te&          data,
    const bool         _blocked
)
: graph(_graph), v_src(_v_src), v_dst(_v_dst), blocked(_blocked), current_traversal_id(0),
  traversal_state(TRAV_UNSEEN), edge_data()
{}

/* private copy ctor */
/* NOTE: this constructor is private (i.e. not publicly accessible) and is only used internally for
 * the implementation of Graph methods.  anway: careful when changing the internals, pointers are
 * copied!  if this is constructor is invoked to change the id of a edge within a graph, that's ok.
 * if however this is used to copy vertices from one graph to another, the caller must take care to
 * adjust the adjacency / incidence information accordingly to avoid UB. */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Edge::Edge(const Edge &x)
{
    this->graph                 = x.graph;
    this->g_eit                 = x.g_eit;
    this->blocked               = x.blocked;
    this->v_src                 = x.v_src;
    this->v_dst                 = x.v_dst;
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->edge_data             = x.edge_data;
}

/* private assignment operator. same as for private copy ctor: to be used only internally and only
 * with care. */
template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::Edge &
Graph<Tg, Tv, Te>::Edge::operator=(const Graph::Edge &x)
{
    this->graph                 = x.graph;
    this->g_eit                 = x.g_eit;
    this->blocked               = x.blocked;
    this->v_src                 = x.v_src;
    this->v_dst                 = x.v_dst;
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->edge_data             = x.edge_data;

    return (*this);
}

/* virtual destructor to allow polymorphic inheritance */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Edge::~Edge()
{
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::Edge *
Graph<Tg, Tv, Te>::Edge::getPtr(typename std::map<uint32_t, EdgePointerType >::const_iterator const &it)
{
    return (it->second);
}

/* iterate through adjacent vertices and replace ids using the given map */
template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Edge::replaceVertices(const std::map<Vertex *, Vertex*> &replace_map)
{
    typename std::map<Vertex *, Vertex *>::const_iterator mit;

    /* replace source and destination vertex if found in replace map */
    if ( (mit = replace_map.find(this->v_src)) != replace_map.end() ) {
        this->v_src = mit->second;
    }
    if ( (mit = replace_map.find(this->v_dst)) != replace_map.end() ) {
        this->v_dst = mit->second;
    }
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::operator<(const Edge &b) const
{
    std::pair<
            Graph<Tg, Tv, Te>::Vertex *,
            Graph<Tg, Tv, Te>::Vertex *
        >                                   e_ptr(this->v_src, this->v_dst),
                                            b_ptr(b.v_src, b.v_dst);

    return Vertex::ptr_less(e_ptr.first, b_ptr.first) ||
            (!Vertex::ptr_less(b_ptr.first, e_ptr.first) && Vertex::ptr_less(e_ptr.second, b_ptr.second));
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::operator>(const Edge &b) const
{
    return (b < (*this));
}

template <typename Tg, typename Tv, typename Te>
uint32_t
Graph<Tg, Tv, Te>::Edge::id() const
{
    return this->g_eit->first;
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_iterator
Graph<Tg, Tv, Te>::Edge::iterator() const
{
    return edge_iterator(this->graph, this->g_eit );
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_const_iterator
Graph<Tg, Tv, Te>::Edge::const_iterator() const
{
    return edge_iterator(this->graph, this->g_eit );
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::isBlocked() const
{
    return (this->blocked == 1);
}

template <typename Tg, typename Tv, typename Te>
Te &
Graph<Tg, Tv, Te>::Edge::operator*()
{
    return (this->edge_data);
}

template <typename Tg, typename Tv, typename Te>
Te const &
Graph<Tg, Tv, Te>::Edge::operator*() const
{
    return (this->edge_data);
}

template <typename Tg, typename Tv, typename Te>
Te *
Graph<Tg, Tv, Te>::Edge::operator->()
{
    return &(this->edge_data);
}

template <typename Tg, typename Tv, typename Te>
Te const *
Graph<Tg, Tv, Te>::Edge::operator->() const
{
    return &(this->edge_data);
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::contains(const vertex_const_iterator &v) const
{
    return (this->contains(&(*v)));
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::contains(const Vertex * const &v) const
{
    return (v == this->v_src || v == this->v_dst);
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::contains(const uint32_t &v_id) const
{
    return (v_id == this->v_src->id() || v_id == this->v_dst->id());
}

/* check if edge has a specific edge-neighbour */
template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::gotNeighbour(const edge_const_iterator &f) const
{
    /* edge f is a neighbour of (this) edge if it shares (at least) one vertex with it */
    return (this->v_src == f->v_src || this->v_src == f->v_dst || this->v_dst == f->v_src || this->v_dst == f->v_dst);
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::gotNeighbour(const Edge * const &f) const
{
    /* call iterator version */
    return (this->gotNeighbour(f->iterator()));
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::Edge::gotNeighbour(const uint32_t &f_id) const
{
    /* search for edge, use edge_const_iterator version */
    Graph<Tg, Tv, Te> const &const_graph    = (*this->graph);
    edge_const_iterator f                   = const_graph.edges.find(f_id);
    
    if (f != const_graph.edges.end()) {
        /* call edge_const_iterator version */
        return (this->gotNeighbour(f));
    }
    /* edge with id f_id not found. throw */
    else {
        throw("Graph::Edge::gotNeighbour() (id version): no edge with specified id found. out of range.");
    }
}

/* get and set traversal state for a given traversal id */
template <typename Tg, typename Tv, typename Te>
uint8_t
Graph<Tg, Tv, Te>::Edge::getTraversalState(const uint32_t &traversal_id)
{
    if (this->current_traversal_id != traversal_id) {
        this->current_traversal_id  = traversal_id;
        this->traversal_state       = TRAV_UNSEEN;
    }
    return this->traversal_state;
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::Edge::setTraversalState(
    const uint32_t &traversal_id,
    const uint8_t  &state)
{
    this->current_traversal_id  = traversal_id;
    this->traversal_state       = state;
}

/* getters for source and destination vertex */
template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::Edge::getSourceVertex()
{
    return (this->v_src->iterator());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_const_iterator
Graph<Tg, Tv, Te>::Edge::getSourceVertex() const
{
    return (this->v_src->iterator());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::Edge::getDestinationVertex()
{
    return (this->v_dst->iterator());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_const_iterator
Graph<Tg, Tv, Te>::Edge::getDestinationVertex() const
{
    return (this->v_dst->iterator());
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        graph class implementation ....                                                        
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* graph ctors */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Graph(Tg const &graph_data) : vertices(*this) , edges(*this)
{
    this->graph_data = graph_data;
}

/* copy ctor */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::Graph(const Graph<Tg, Tv, Te> &X) : vertices(*this), edges(*this) {
    /* use assignment operator. although this initializes all members with the default ctor and
     * immediately overwrites them again, this was deemed preferrable to copying the code of
     * Graph::operator=() and have virtually the same piece of code twice.  the alternative: put the
     * code in the copy ctor and use a swap trick in operator=(), which was deemed too inefficient,
     * since (partially quite large) graphes are often assigned, but not so much copy-constructed.
     * furthermore, the default construction is very little work. */
    this->operator=(X);
}

/* assignment operator: special care has to be taken, since the topology information inside a graph
 * uses iterators and pointers, which must of course not be copied during assignment to prevent
 * pointer aliasing and races. */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>&
Graph<Tg, Tv, Te>::operator=(const Graph<Tg, Tv, Te>&X)
{
    /* clear all data */
    this->clear();

    /* assign idqs, ids used in X will be reused in (this) graph, although this must not be relied
     * upon by the caller: its just more efficient and easiert to realize, since ids can be reused
     * during the construction of topological information for (this) graph */
    this->V_idq             = X.V_idq;
    this->E_idq             = X.E_idq;
    this->traversal_idq     = X.traversal_idq;

    /* copy graph data */
    this->graph_data        = X.graph_data;

    /* copy vertex map by value: this will temporarily copy pointers referring to data from X into
     * (this) graph, which is not desirable in the final result. however, all ids are copied
     * correctly. in the next step, "deep copy" all elements in this->V, i.e. allocate new vertices
     * and copy-assign the vertices from X. note that the link between (this) graph and X still isn't
     * broken, since the Vertex objects contains pointers, which still refer to components of X. */
    this->V                 = X.V;

    /* deep copy */
    Vertex *v_new;
    typename std::map<uint32_t, VertexPointerType >::iterator vit;
    for (vit = this->V.begin(); vit != this->V.end(); ++vit) {
        /* make a copy of the Vertex object currently pointed to by vit, which is a Vertex object
         * allocated by X, with the private copy ctor of Vertex */
        v_new           = new Vertex(*(vit->second)); 

        /* store pointer to new copy in iterator */
        vit->second     = VertexPointerType(v_new);

        /* although the vertex has been deep-copied, it still contains pointers referring to X.
         * remove entirely. all topological information for (this) graph will be constructed will be
         * constructed below */

        /* firstly, clear adjacency and incidence information */
        v_new->out_edges.clear();
        v_new->in_edges.clear();

        /* set graph pointer to (this) */
        v_new->graph    = this;

        /* set iterator to vit */
        v_new->g_vit    = vit;
    }

    /* now the edges of X are added using the integer _id_ based version (neither copy ctor nor
     * iterator based version). */
    edge_iterator       eit;
    bool                e_inserted;
    for (auto &xe : X.edges) {
        auto rpair  = this->edges.insert(xe.getSourceVertex()->id(), xe.getDestinationVertex()->id(), xe.edge_data);
        eit         = rpair.first;
        e_inserted  = rpair.second;

        if (!e_inserted) {
            throw("Graph::operator=(): could not insert edge into (this) mesh. internal logic error.");
        }

        /* copy the rest of the information from xf */
        eit->blocked                = xe.blocked;
        eit->current_traversal_id   = xe.current_traversal_id;
        eit->traversal_state        = xe.traversal_state;
        //eit->data                   = xe.data;
    }

    return (*this);
}

/* dtor */
template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::~Graph()
{
    /* delete all allocated vertices */
    for (auto &v : this->vertices) {
        delete (&v);
    }

    /* delete all allocated edges */
    for (auto &e : this->edges) {
        delete (&e);
    }
}

/* access to data object of template type Tg */
template <typename Tg, typename Tv, typename Te>
Tg &
Graph<Tg, Tv, Te>::data()
{
    return (this->graph_data);
}

template <typename Tg, typename Tv, typename Te>
Tg const &
Graph<Tg, Tv, Te>::data() const
{
    return (this->graph_data);
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::clear()
{
    /* delete all allocated vertices */
    for (auto &v : this->vertices) {
        delete (&v);
    }

    /* delete all allocated edges */
    for (auto &e : this->edges) {
        delete (&e);
    }

    /* clear vertex / edge maps */
    this->V.clear();
    this->E.clear();

    /* clear id queues */
    this->V_idq.clear();
    this->E_idq.clear();
    this->traversal_idq.clear();
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::clearEdges()
{
    /* delete all allocated edges */
    for (auto &e : this->edges) {
        delete (&e);
    }

    /* clear edges map and edge id queue.*/
    this->E.clear();
    this->E_idq.clear();

    /* all vertices are now isolated, simply clear all {in,out}_edges lists for all vertices */
    for (auto &vit : this->V) {
        vit.second->in_edges.clear();
        vit.second->out_edges.clear();
    }
}

/* renumber vertices and edges consecutively, starting from vertex_start_id for vertices and
 * vedge_start_id for edges. */
template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::renumberConsecutively(
    uint32_t vertex_start_id,
    uint32_t edge_start_id)
{
    debugl(1, "Graph::renumberConsecutively(): vertex_start_id: %5d, edge_start_id: %5d.\n", vertex_start_id, edge_start_id);
    debugTabInc();

    /* renumber indices of vertices and edges in a consecutive fashion. since the internal maps
     * Graph::V and Graph::E store pointers to the data, the data itself is not touched by
     * manipulation of these maps and hence all topological information inside Vertex and Edge
     * objects remains valid EXCEPT the internal interators, which have to be reset to bring the
     * entire graph back to a consistent state.
     *
     * in general, all iterators are invalidated by this method */
    std::map<uint32_t, VertexPointerType >    vertices_swap; 
    std::map<uint32_t, EdgePointerType >      edges_swap;
    
    /* swap vertices and edges with vertices_swap / edges_swap in-place */
    this->V.swap(vertices_swap);
    this->E.swap(edges_swap);

    /* reset id queues for vertices / edges to start at vertex_start_id / edge_start_id */
    this->V_idq = IdQueue(vertex_start_id);
    this->E_idq = IdQueue(edge_start_id);

    /* iterate through swap arrays and insert Vertex and Edge shared pointers into now empty 
     * maps this->V and this->E with correct ids */
    typename std::map<uint32_t, VertexPointerType >::iterator vit;
    Vertex *v;
    for (uint32_t current_vertex_id = vertex_start_id; !vertices_swap.empty(); current_vertex_id++) {
        vit         = vertices_swap.begin();
        v           = Vertex::getPtr(vit);

        /* insert pointer to vertex (vit->second) into this->V with correct id, store new
         * iterator, which is the first element of the result (iterator, bool) pair. */
        auto rpair  = ( this->V.insert( {current_vertex_id, VertexPointerType(vit->second) } ) );

        /* vertex could not be inserted */
        if (!rpair.second) {
            throw("Graph::renumberConsecutively(): failed to insert renumbered vertex into (this) graph. internal logic error.");
        }
        /* update iterator inside vertex */
        v->g_vit    = rpair.first;

        /* erase vertex from vertices_swap */
        vertices_swap.erase(vit);

        /* pop() id from id queue. this is guaranteed to be consecutive since IdQueue object has
         * been initialized above and returned consecutive ids as long as no IdQuee::freeId() is
         * ever issued */
        this->V_idq.getId();
    }

    /* same for all edges */
    typename std::map<uint32_t, EdgePointerType >::iterator eit;
    Edge *e;

    for (uint32_t current_edge_id = edge_start_id; !edges_swap.empty(); current_edge_id++) {
        eit         = edges_swap.begin();
        e           = Edge::getPtr(eit);
        auto rpair  = ( this->E.insert( {current_edge_id, EdgePointerType(eit->second) } ) );

        /* edge could not be inserted */
        if (!rpair.second) {
            throw("Graph::renumberConsecutively(): failed to insert renumbered edge in (this) graph. internal logic error.");
        }
        /* update iterator inside edge */
        e->g_eit    = rpair.first;

        edges_swap.erase(eit);
        this->E_idq.getId();
    }

    debugTabDec();
    debugl(1, "Graph::renumberConsecutively(). done.\n");
}

template <typename Tg, typename Tv, typename Te>
uint32_t
Graph<Tg, Tv, Te>::numVertices() const
{
    return this->vertices.size();
}

template <typename Tg, typename Tv, typename Te>
uint32_t
Graph<Tg, Tv, Te>::numEdges() const
{
    return this->edges.size();
}

template <typename Tg, typename Tv, typename Te>
uint32_t
Graph<Tg, Tv, Te>::getFreshTraversalId()
{
    uint32_t tid = this->traversal_idq.getId();

    /* reset after 2^23 - 1 = 8388607 traversals, since traversal_id is dimensioned to be a 23 bit
     * unsigned integer.  this way, at least noone else will have to mess around with the global
     * resetTraversalStates() method.. */
    if (tid == 8388607u) { 
        /* reset traversal states, which clear()s the id queue for traversal ids */
        this->resetTraversalStates();
        tid = this->traversal_idq.getId();
    }
    return tid;
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::copyAppend(const Graph<Tg, Tv, Te>&B)
{
    std::map<uint32_t, vertex_iterator>                     new_vertex_its;
    typename std::map<uint32_t, vertex_iterator>::iterator  idit;
    vertex_iterator                                         new_it;

    debugl(1, "Graph::appendCopy()\n");
    debugTabInc();

    /* add all vertices of B to (this) graph, store iterators to new vertices */
    for (auto &B_v : B.vertices) {
        new_it = this->vertices.insert(*B_v);

        if ( (idit = new_vertex_its.find( B_v.id() )) != new_vertex_its.end()) {
            throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::appendCopy(): duplicate entry in new_vertex_ids map, should be impossible.");
        }
        new_vertex_its[B_v.id()] = new_it;  

        debugl(1, "added vertex %5d from b under new id %5d.\n", B_v.id(), new_it->id() );
    }

    /* add all edges of B to (this) graph while taking care to replace old vertex ids from B with the
     * proper new ids. */
    for (auto &B_e : B.edges) {
        this->edges.insert(
            new_vertex_its[B_e.getSourceVertex()->id()],
            new_vertex_its[B_e.getDestinationVertex()->id()],
            *B_e);
    }
    debugTabDec();
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::moveAppend(
    Graph<Tg, Tv, Te>              &B,
    std::list<vertex_iterator>     *update_vits)
{
    debugl(2, "Graph::moveAppend()\n");
    debugTabInc();

    uint32_t                            new_id;
    Vertex                             *v;
    typename std::map<
            uint32_t,
            VertexPointerType
        >::iterator                     v_newit;

    Edge                               *e;
    typename std::map<
            uint32_t,
            EdgePointerType
        >::iterator                     e_newit;
    bool                                inserted;
        
    /* if update_vits != NULL, generate the list of vertex pointers from the given iterators.
     * although iterators are invalidated, pointers remain valid during this method, since no
     * vertices or edges are allocated or deleted during moving. */
    std::list<Vertex *> update_list;
    if (update_vits) {
        update_list.resize(update_vits->size());

        auto vit = update_vits->begin();
        auto lit = update_list.begin();
        while (vit != update_vits->end()) {
            *lit = &(**vit);
            ++vit, ++lit;
        }
    }

    /* add all vertices of B to (this) graph, store iterators to new vertices */
    auto B_vit = B.V.begin();
    while (B_vit != B.V.end()) {
        /* get fresh vertex id */
        new_id      = this->V_idq.getId();

        /* insert pair (new_id, pointer B_vit->second) into this->V */
        auto rpair  = this->V.insert( { new_id, VertexPointerType(B_vit->second) } );

        /* if vertex has been successfully inserted */
        if (rpair.second) {
            /* the vertex pointer has been moved to (this) graph: update graph pointer Vertex::graph
             * and iterator Vertex::g_vit */
            /* get pointer to vertex */
            v_newit     = rpair.first;
            v           = v_newit->second;
            v->graph    = this;
            v->g_vit    = v_newit;

            /* erase B_vit from B.V, call returns iterator to next vertex to be processed (or end()) */
            B_vit       = B.V.erase(B_vit); 
        }
        /* vertex insertion failed */
        else {
            throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::moveAppend(): insertion of vertex pointer into this->V with fresh id failed. internal logic error.");
        }
    }

    /* move all edges of B to (this) graph in very much the same way */
    auto B_eit = B.E.begin();
    while (B_eit != B.E.end()) {
        new_id      = this->E_idq.getId();

        auto rpair  = this->E.insert( { new_id, EdgePointerType(B_eit->second) } );
        inserted    = rpair.second;
        if (inserted) {
            e_newit     = rpair.first;
            e           = e_newit->second;
            e->graph    = this;
            e->g_eit    = e_newit;
            B_eit       = B.E.erase(B_eit); 
        }
        else {
            throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::moveAppend(): insertion of edge pointer into this->E with fresh id failed. internal logic error.");
        }
    }

    /* clear all information from B (B.V and B.F are empty, yet id queues etc are still set */
    if (!B.E.empty() || !B.V.empty()) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::moveAppend(): internal vertex and edge maps of appended graph not empty at the end of process. internal logic error.");
    }
    B.clear();

    /* if update_vits != NULL, the list update_list of vertex pointers corresponding to the list of
     * given old vertex iterators (referring to B before the append) has been generated above: use
     * it to update the iterators in-place, i.e. an old iterator referring to a vertex of B is
     * replaced with an iterator referring to this moved vertex as a component of (this) graph. */
    if (update_vits) {
        auto vit = update_vits->begin();
        auto lit = update_list.begin();
        while (vit != update_vits->end()) {
            *vit= (*lit)->iterator();
            ++vit, ++lit;
        }
    }
    debugTabDec();
    debugl(2, "Graph::moveAppend(): done.\n");

}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::deleteConnectedComponent(vertex_iterator vstart_it)
{
    debugl(1, "Graph::deleteConnectedComponent(). start vertex: %5d\n", vstart_it->id());
    debugTabInc();

    Vertex                 *v;
    std::queue<Vertex *>    Q;
    std::list<Vertex *>     v_nbs;
    std::list<Vertex *>     cc_vertices;

    const uint32_t          traversal_id = this->getFreshTraversalId();

    /* traverse cc of vstart_it, burning it down as we go.. */
    Q.push( &(*vstart_it) );
    vstart_it->setTraversalState(traversal_id, TRAV_ENQUEUED);

    debugl(2, "finding all vertices of connected component.\n");
    debugTabInc();
    while (!Q.empty()) {
        v = Q.front();
        Q.pop();

        debugl(3, "current vertex %5d\n", v->id() );

        /* get out neighbours of v */
        v->getOutNeighbours(v_nbs);

        debugTabInc();
        /* iterator over all neighbours u and enqueue them if they haven't been seen yet */
        for (auto &u : v_nbs) {
            if (u->getTraversalState(traversal_id) == TRAV_UNSEEN) {
                debugl(3, "yet unseen neighbour %5d => enqueueing..\n", u->id() );
                Q.push(u);
                u->setTraversalState(traversal_id, TRAV_ENQUEUED);
            }
        }
        debugTabDec();

        /* v is done */
        debugl(3, "vertex %5d done.\n", v->id() );
        v->setTraversalState(traversal_id, TRAV_DONE);
        cc_vertices.push_back(v);
    }
    debugTabDec();

    debugl(2, "traversal of connected component completed. deleting..\n");

    debugTabInc();
    /* delete all vertices of cc, which in turn deletes all edges containing any such vertex. */
    for (auto &v : cc_vertices) {
        debugl(3, "deleting vertex %5d\n", v->id() );
        this->vertices.erase(v->iterator()) ;
        debugl(3, "vertex deleted.\n");
    }
    debugTabDec();

    debugTabDec();
    debugl(1, "Graph::deleteConnectedComponent(). done.\n");
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::getConnectedComponentBreadthFirst(
    Graph::vertex_iterator      vstart_it,    
    const uint32_t             &traversal_id,
    std::list<Vertex *>        *cc_vertices,
    std::list<Edge *>          *cc_edges)
{
    if (!vstart_it.checkContainer(*this) || vstart_it == this->vertices.end()) {
        throw("Graph::getConnectedComponent(): invalid start vertex iterator.");
    }

    debugl(1, "Graph::getConnectedComponent(). start vertex: %5d\n", vstart_it->id());
    debugTabInc();

    /* only process if one of the given list pointers isn't NULL and the start vertex's traversal
     * state is set to UNSEEN. */
    if ( (cc_vertices || cc_edges) && vstart_it->getTraversalState(traversal_id) == TRAV_UNSEEN) {
        std::queue<Vertex *>    Q;

        Vertex                 *v;
        std::list<Vertex *>     v_out_neighbours;
        std::list<Edge *>       v_out_edges;

        /* traverse cc of vstart_it, burning it down as we go.. */
        Q.push( &(*vstart_it) );
        vstart_it->setTraversalState(traversal_id, TRAV_ENQUEUED);

        debugl(2, "traversing connected component of start vertex %5d.\n", Q.front()->id());
        debugTabInc();
        while (!Q.empty()) {
            v = Q.front();
            Q.pop();

            debugl(3, "current vertex %5d\n", v->id() );

            /* add v and all its incident edges to the result lists if desired by the caller */
            if (cc_vertices) {
                cc_vertices->push_back(v);
            }
            if (cc_edges) {
                /* add all out-edges of v to list */
                v->getOutEdges(v_out_edges);
                cc_edges->insert(cc_edges->end(), v_out_edges.begin(), v_out_edges.end());
            }

            /* get all out-neighbours of v */
            v->getOutNeighbours(v_out_neighbours);

            debugTabInc();
            /* iterate over all vertex neighbours u and enqueue them if they haven't been seen yet. */
            for (auto &u : v_out_neighbours) {
                if (u->getTraversalState(traversal_id) == TRAV_UNSEEN) {
                    debugl(3, "yet unseen neighbour %5d => enqueueing..\n", u->id() );
                    Q.push(u);
                    u->setTraversalState(traversal_id, TRAV_ENQUEUED);
                }
            }
            debugTabDec();

            /* v is done. */
            debugl(3, "vertex %5d done.\n", v->id() );
            v->setTraversalState(traversal_id, TRAV_DONE);
        }
        debugTabDec();
        debugl(2, "traversal of connected component completed.\n");
    }

    debugTabDec();
    debugl(1, "Graph::getConnectedComponent(). done.\n");
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::checkEdge(
    const std::string&              fn,
    const vertex_const_iterator    &u_it,
    const vertex_const_iterator    &v_it) const
{
    if (!u_it.checkContainer(*this) || !u_it.sameContainer(v_it) || u_it == v_it) {
        throw GraphEx(GRAPH_LOGIC_ERROR, fn + "given edge (u, v) invalid: u and v not from (this) graph or identical.");
    }
    else if (!u_it->gotOutNeighbour(v_it)) {
        throw GraphEx(GRAPH_LOGIC_ERROR, fn + "given edge (u, v) invalid: both vertices from (this) graph, yet edge does not exist.");
    }
}

template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::checkEdge(
    const char                     *fn,
    const vertex_const_iterator    &u_it,
    const vertex_const_iterator    &v_it) const
{
    this->checkEdge(std::string(fn), u_it, v_it);
}


template <typename Tg, typename Tv, typename Te>
void
Graph<Tg, Tv, Te>::checkInternalConsistency() const
{
    debugl(0, "Graph::checkInternalConsistency()\n");
    debugTabInc();

    /* FIXME: adapt this to graph */

    debugTabDec();
    debugl(0, "Graph::checkInternalConsistency(): graph internally consistent.\n");
}

template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::VertexAccessor::VertexAccessor(Graph<Tg, Tv, Te> &g) : graph(g) 
{
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::VertexAccessor::begin()
{
    return Graph::vertex_iterator( &(this->graph), this->graph.V.begin());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_const_iterator
Graph<Tg, Tv, Te>::VertexAccessor::begin() const
{
    /* the seemingly returned vertex_iterator will be implicitly converted to a
     * vertex_const_iterator to match the return type */
    return Graph::vertex_iterator( &(this->graph), this->graph.V.begin());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::VertexAccessor::end()
{
    return Graph::vertex_iterator( &(this->graph), this->graph.V.end());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_const_iterator
Graph<Tg, Tv, Te>::VertexAccessor::end() const
{
    return Graph::vertex_const_iterator( &(this->graph), this->graph.V.end());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::VertexAccessor::find(uint32_t id)
{
    return Graph::vertex_iterator( &(this->graph), this->graph.V.find(id));
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_const_iterator
Graph<Tg, Tv, Te>::VertexAccessor::find(uint32_t id) const
{
    return Graph::vertex_const_iterator( &(this->graph), this->graph.V.find(id));
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::Vertex &
Graph<Tg, Tv, Te>::VertexAccessor::at(uint32_t id)
{
    try {
        return *(this->graph.V.at(id));
    }
    catch (std::out_of_range &e) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::VertexAccessor::at(): vertex index out of range.");
    }
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::VertexAccessor::exists(uint32_t id) const
{
    return (this->graph.V.find(id) != this->graph.V.end());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::VertexAccessor::insert(Tv const &data)
{
    debugl(2, "Graph::VertexAccessor::insert().\n");
    debugTabInc();

    std::pair<typename std::map<uint32_t, VertexPointerType >::iterator, bool> pair;

    /* get fresh id for new vertex, allocate new vertex, insert pair (id, vertex) into map */
    uint32_t v_id   = this->graph.V_idq.getId();
    Vertex *v       = new Vertex(&(this->graph), data);
    pair            = this->graph.V.insert( { v_id, VertexPointerType(v) } );
    if (!pair.second) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "vertex with fresh id from idq already present in vertex map. this must never happen..");
    }
    else {
        /* the (private) constructor Vertex(Graph * const &, const Tv & const &data) does not set the internal iterator
         * Vertex::g_vit, since it has no way of knowing it in advance. set it now to bring the new
         * vertex into a consistent state. */
        v->g_vit    = pair.first;

        debugTabDec();
        debugl(2, "Graph::VertexAccessor::insert(). done.\n");

        return (v->iterator());
    }
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::VertexAccessor::erase(vertex_iterator it)
{
    debugl(2, "Graph::VertexAccessor::erase(vertex_iterator it): vertex id: %6d. deleting %2zu and %2zu incident out and in edges, respectively.\n", it->id(), it->in_edges.size(), it->out_edges.size() );
    debugTabInc();

    if (!it.checkContainer( this->graph)) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::VertexAccessor::erase(): given input iterator refers to different container.");
    }

    /* delete all edges that the vertex is incident to. NOTE: the erase call will alter
     * vit->in_edges and vit->out_edges: always take the begin() iterator, which is guaranteed always to be valid
     * as long as the list is non-empty. */
    while (!it->in_edges.empty()) {
        this->graph.edges.erase( (it->in_edges.front())->iterator() );
    }
    while (!it->out_edges.empty()) {
        this->graph.edges.erase( (it->out_edges.front())->iterator() );
    }

    debugl(3, "freeing id and erasing vertex from internal vertex map..\n");

    /* free id */
    this->graph.V_idq.freeId(it->id());

    debugl(3, "deleting (deallocating) vertex object..\n");
    /* delete allocated vertex object */
    delete &(*it);

    debugTabDec();
    debugl(2, "Graph::VertexAccessor::erase(). erase()ing and returning vertex_iterator to next vertex.\n");

    /* erase vertex from internal map this->graph.V and return vertex_iterator to the next element
     * (by wrapping the result of std::map::erae in the vertex_iterator constructor). */
    return Graph::vertex_iterator( &(this->graph), this->graph.V.erase(it.int_it) );
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::VertexAccessor::erase(uint32_t id)
{
    /* try to get iterator and call iterator version of erase */
    Graph::vertex_iterator vit;
    if ( (vit = this->find(id)) != this->end() ) {
        this->erase(vit);
        return true;

    }
    else {
        return false;
    }
}

template <typename Tg, typename Tv, typename Te>
size_t
Graph<Tg, Tv, Te>::VertexAccessor::size() const
{
    return (this->graph.V.size());
}

/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                        graph edge accessor class implementation..                                                    
 *
 * ----------------------------------------------------------------------------------------------------------------- */

template <typename Tg, typename Tv, typename Te>
Graph<Tg, Tv, Te>::EdgeAccessor::EdgeAccessor(Graph<Tg, Tv, Te> &g) : graph(g) 
{
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::begin()
{
    return Graph::edge_iterator( &(this->graph), this->graph.E.begin() );
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_const_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::begin() const
{
    /* edge_iterator will be implicitly converted to edge_const_iterator */
    return Graph::edge_iterator( &(this->graph), this->graph.E.begin() );
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::end()
{
    return Graph::edge_iterator( &(this->graph), this->graph.E.end());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_const_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::end() const
{
    /* edge_iterator will be implicitly converted to edge_const_iterator */
    return Graph::edge_iterator( &(this->graph), this->graph.E.end());
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::find(uint32_t id)
{
    return Graph::edge_iterator( &(this->graph), this->graph.E.find(id));
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_const_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::find(uint32_t id) const
{
    return Graph::edge_const_iterator( &(this->graph), this->graph.E.find(id));
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::find(
    vertex_const_iterator const    &v_src_it,
    vertex_const_iterator const    &v_dst_it)
{
    /* check whether the two vertex iterators refer to (this) graph! */
    if (!( v_src_it.checkContainer(this->graph) && v_src_it.sameContainer(v_dst_it) ) ) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::exists(): invalid input iterators: not both referring to (this) graph.");
    }

    /* is edge (v_src, v_dst) is an in-edge of v_dst */
    for (auto &e_out : v_src_it->getOutEdges()) {
        if (e_out->getDestinationVertex()->const_iterator() == v_dst_it) {
            return (e_out->iterator());
        }
    }
    return this->end();
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_const_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::find(
    vertex_const_iterator const    &v_src_it,
    vertex_const_iterator const    &v_dst_it) const
{
    /* check whether the two vertex iterators refer to (this) graph! */
    if (!( v_src_it.checkContainer(this->graph) && v_src_it.sameContainer(v_dst_it) ) ) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::exists(): invalid input iterators: not both referring to (this) graph.");
    }

    /* is edge (v_src, v_dst) is an in-edge of v_dst */
    for (auto &e_out : v_src_it->getOutEdges()) {
        if (e_out->getDestinationVertex()->const_iterator() == v_dst_it) {
            return (e_out->const_iterator());
        }
    }
    return this->end();
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::Edge &
Graph<Tg, Tv, Te>::EdgeAccessor::at(uint32_t id)
{
    try {
        return *(this->graph.E.at(id));
    }
    catch (std::out_of_range &e) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::at(): vertex index out of range.");
    }
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::Edge const &
Graph<Tg, Tv, Te>::EdgeAccessor::at(uint32_t id) const
{
    try {
        return *(this->graph.E.at(id));
    }
    catch (std::out_of_range &e) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::at(): vertex index out of range.");
    }
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::EdgeAccessor::exists(uint32_t id) const
{
    return (this->graph.E.find(id) != this->graph.E.end());
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::EdgeAccessor::exists(
    vertex_const_iterator const    &v_src_it,
    vertex_const_iterator const    &v_dst_it) const
{
    /* check whether the two vertex iterators refer to (this) graph! */
    if (!( v_src_it.checkContainer(this->graph) && v_src_it.sameContainer(v_dst_it) ) ) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::exists(): invalid input iterators: not both referring to (this) graph.");
    }

    /* check if &(*v_dst_it) is contained in the neighbours v_dst_it's in-neighbours (and the
     * symmetric case as an additional check for internal consistency) */
    const Vertex *v_src = &(*v_src_it), *v_dst = &(*v_dst_it);

    /* edge (v_src, v_dst) is an in-edge of v_dst */
    if (Aux::Alg::listContains<const Vertex *>(v_src->getOutNeighbours(), v_dst)) {
        /* for consistency, check whether v_dst is also an out-neighbour of v_src */
        if (!Aux::Alg::listContains<const Vertex *>(v_dst->getInNeighbours(), v_src)) {
            throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::exists(): v_dst is out-neighbour of v_src, yet v_src not in-neighbour of v_dst. internal logic error.");
        }
        else return true;
    }
    else return false;
}

template <typename Tg, typename Tv, typename Te>
std::pair<typename Graph<Tg, Tv, Te>::edge_iterator, bool>
Graph<Tg, Tv, Te>::EdgeAccessor::insert(
    vertex_iterator     v_src_it,
    vertex_iterator     v_dst_it,
    Te const           &data)
{
    debugl(2, "Graph::EdgeAccessor::insert().\n");
    debugTabInc();

    std::pair<
            typename Graph<Tg, Tv, Te>::edge_iterator,
            bool
        >               ret;

    if (this->exists(v_src_it, v_dst_it)) {
        ret.first.explicitlyInvalidate();
        ret.second = false;
    }

    uint32_t    edge_id;
    Edge       *e;

    edge_id     = this->graph.E_idq.getId();
    e           = new Edge (&(this->graph), &(*v_src_it), &(*v_dst_it), data);
    auto rpair  = this->graph.E.insert( {edge_id, EdgePointerType(e) } );

    if (rpair.second) {
        /* set Edge::m_fit iterator, which is required for Edge to be in a consistent state and has not
         * been set by the (private) Edge ctor, just as for Graph::Vertex */
        e->g_eit = rpair.first;

        /* topology information update */
        bool all_inserted = 
            Aux::Alg::listSortedInsert(v_src_it->out_edges, e, false) &&
            Aux::Alg::listSortedInsert(v_dst_it->in_edges, e, false);

        if (all_inserted) {
            ret.first   = e->iterator();
            ret.second  = true;

            debugTabDec();
            debugl(2, "Graph::EdgeAccessor::insert().\n");

            return ret;
        }
        else {
            throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::insert(): at least one (sorted) insertion of edge pointer into {out/in} lists of {v_src/v_dst} failed. internal logic error.");
        }
    }
    else {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::insert(): new Edge with fresh id from idq already present in Edge map. this must never happen..");
    }
}

template <typename Tg, typename Tv, typename Te>
std::pair<typename Graph<Tg, Tv, Te>::edge_iterator, bool>
Graph<Tg, Tv, Te>::EdgeAccessor::insert(
    uint32_t    v_src_id,
    uint32_t    v_dst_id,
    Te const   &data)
{
    std::pair<
            typename Graph<Tg, Tv, Te>::edge_iterator,
            bool
        >               ret;
    vertex_iterator     v_src_it, v_dst_it, vend;

    v_src_it    = this->graph.vertices.find(v_src_id);
    v_dst_it    = this->graph.vertices.find(v_dst_id);
    vend        = this->graph.vertices.end();

    /* note: self-edges (v, v) ARE allowed here.. */
    if (v_src_it == vend || v_dst_it == vend) {
        ret.first.explicitlyInvalidate();
        ret.second = false;
        return ret;
    }
    else return (this->insert(v_src_it, v_dst_it, data));
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::edge_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::erase(edge_iterator it)
{
    debugl(2, "Graph::EdgeAccessor::erase(): erasing edge (%5d, %5d) with id %5d\n", it->v_src->id(), it->v_dst->id(), it->id());
    debugTabInc();

    /* check if edge belongs to (this) graph */
    if (!it.checkContainer(this->graph)) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::erase(): given iterator does not refer to (this) graph.");
    }

    /* get pointer to edge */
    Edge *e = &(*it);

    /* remove edge from {out, in}-edge lists of {src, dst} vertex */
    bool all_erased = 
        Aux::Alg::removeFirstOccurrenceFromList(e->v_src->out_edges, e) &&
        Aux::Alg::removeFirstOccurrenceFromList(e->v_dst->in_edges, e);

    if (!all_erased) {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::erase(): source or destination vertex does not contain given edge in its incidence lists. invalid input iterator?");
    }

    if (
        Aux::Alg::listContains(e->v_src->out_edges, e) ||
        Aux::Alg::listContains(e->v_dst->in_edges, e)
       )
    {
        throw GraphEx(GRAPH_LOGIC_ERROR, "Graph::EdgeAccessor::erase(): {src,dst} vertex {out,in} edge lists still contains edge to be deleted after removing first occurrence => edge has been in list multiple times. internal logic error.");
    }

    /* free edge_id */
    this->graph.E_idq.freeId( it->id() );

    /* delete allocated edge object */
    delete &(*it);

    debugTabDec();
    debugl(2, "Graph::EdgeAccessor::erase(): done.\n");

    /* return edge_iterator to next element by wrapping return iterator of map::erase inside a edge_iterator */
    return Graph::edge_iterator( &(this->graph), this->graph.E.erase(it.int_it));
}

template <typename Tg, typename Tv, typename Te>
bool
Graph<Tg, Tv, Te>::EdgeAccessor::erase(uint32_t id)
{
    /* locate edge with given id. if found, remove it with iterator version of erase */
    edge_iterator it = this->find(id);
    if ( it != this->end() ) {
        this->erase(it);
        return true;
    }
    /* edge not found, return false */
    else {
        return false;
    }
}

template <typename Tg, typename Tv, typename Te>
typename Graph<Tg, Tv, Te>::vertex_iterator
Graph<Tg, Tv, Te>::EdgeAccessor::collapse(edge_iterator it)
{
    debugl(0, "Graph::EdgeAccessor::collapse(): edge e = (u, v) | %d = (%d, %d).\n",
        it->id(),
        it->getSourceVertex()->id(),
        it->getDestinationVertex()->id());
    debugTabInc();

    vertex_iterator u_it = it->getSourceVertex();
    vertex_iterator v_it = it->getDestinationVertex();

    debugl(1, "erasing e..\n");
    /* delete edge it */
    this->graph.edges.erase(it);

    debugl(1, "remaining edge-neighbours of u: \n");
    debugTabInc();
    std::list<Edge *> nbs, tmp;

    u_it->getInEdges(nbs);
    u_it->getOutEdges(tmp);
    nbs.insert(nbs.end(), tmp.begin(), tmp.end());

#ifdef __DEBUG__
    for (auto &nb : nbs) {
        debugl(0, "%d = (%d, %d)\n", nb->id(), nb->getSourceVertex()->id(), nb->getDestinationVertex()->id());
    }
#endif

    debugTabDec();

    debugl(0, "remaining edge-neighbours of v: \n");
    debugTabInc();

    nbs.clear(), tmp.clear();
    v_it->getInEdges(nbs);
    v_it->getOutEdges(tmp);
    nbs.insert(nbs.end(), tmp.begin(), tmp.end());

#ifdef __DEBUG__
    for (auto &nb : nbs) {
        debugl(0, "%d = (%d, %d)\n", nb->id(), nb->getSourceVertex()->id(), nb->getDestinationVertex()->id());
    }
#endif

    debugTabDec();

    /* wlog, choose u to be the vertex that "survives" the collapse, i.e. v is deleted => redirect all in- and out-edges
     * of v to u by replacing destination / source vertices. */
    for (auto &out_e : v_it->out_edges) {
        out_e->v_src = &(*u_it);
        u_it->out_edges.push_back(out_e);
    }
    v_it->out_edges.clear();

    for (auto &in_e : v_it->in_edges) {
        in_e->v_dst = &(*u_it);
        u_it->in_edges.push_back(in_e);
    }
    v_it->in_edges.clear();

    /* delete v, which is isolated now. */
    if (v_it->outdeg() + v_it->indeg() != 0) {
        throw("Graph::EdgeAccessor::collapse(): vertex to be deleted not isolated after redirecting all edges."\
            "internal logic error.");
    }
    this->graph.vertices.erase(v_it);

    /* sort edges of u */
    u_it->in_edges.sort(Edge::ptr_less);
    u_it->out_edges.sort(Edge::ptr_less);

    debugl(0, "neighbours of collapsed vertex: \n");
    debugTabInc();
    u_it->getInEdges(nbs);
    u_it->getOutEdges(tmp);
    nbs.insert(nbs.end(), tmp.begin(), tmp.end());

#ifdef __DEBUG__
    for (auto &nb : nbs) {
        debugl(0, "%d = (%d, %d)\n", nb->id(), nb->getSourceVertex()->id(), nb->getDestinationVertex()->id());
    }
#endif

    debugTabDec();

    debugTabDec();
    debugl(0, "Graph::EdgeAccessor::collapse(). done\n");

    /* return iterator */
    return u_it;
}

template <typename Tg, typename Tv, typename Te>
size_t
Graph<Tg, Tv, Te>::EdgeAccessor::size() const
{
    return (this->graph.E.size());
}
