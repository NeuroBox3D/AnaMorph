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

#ifndef GRAPH_H
#define GRAPH_H

#include "common.hh"
#include "IdQueue.hh"
#include "Octree.hh"

enum graph_error_types {
    GRAPH_NOERROR,
    GRAPH_NOT_FOUND,
    GRAPH_LOGIC_ERROR,
    GRAPH_IO_ERROR,
    GRAPH_NONMANIFOLD,
    GRAPH_INVALID_ID,
    GRAPH_CSG_EDGE_CASE,
    GRAPH_CSG_COMPLEX_EDGE,
    GRAPH_CSG_NUM_ISEC_POLY,
    GRAPH_CSG_TRIANGULATION,
    GRAPH_CSG_FATAL
};

struct GraphEx : public std::runtime_error {
    const uint32_t      error_type;
    const std::string   error_msg;

    GraphEx(const uint32_t &type, const std::string &msg) :
        std::runtime_error(msg),
        error_type(type),
        error_msg(msg)
    {
    }
};

struct GraphEx_OutOfRange : public GraphEx {
};

template<
    typename    Tg,
    typename    Tv,
    typename    Te
>
class Graph {
    public:
        /* ----------------- iterator related declarations ---------------------- */
        template <
            typename ValueType,
            typename InternalIteratorType
        >
        class GraphIterator : public std::iterator<std::bidirectional_iterator_tag, ValueType>
        {
            friend class Graph<Tg, Tv, Te>;

            protected:
                Graph<Tg, Tv, Te>          *graph;
                InternalIteratorType        int_it;

                /* private ctor, can only be called by Graph */
                GraphIterator(
                    Graph<Tg, Tv, Te>      *g,
                    InternalIteratorType    it)
                {
                    this->graph     = g;
                    this->int_it    = it;
                }


            public:
                GraphIterator()
            	: graph(NULL)
                {
                }

                GraphIterator(const GraphIterator &x)
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                }

                ~GraphIterator()
                {
                }

                inline GraphIterator &
                operator=(const GraphIterator &x)
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                    return (*this);
                }

                inline GraphIterator &
                operator++()
                {
                    ++this->int_it;
                    return (*this);
                }

                inline GraphIterator
                operator++(int)
                {
                    GraphIterator tmp(*this);
                    this->operator++();
                    return tmp;
                }

                inline GraphIterator &
                operator--()
                {
                    --this->int_it;
                    return (*this);
                }

                inline GraphIterator
                operator--(int)
                {
                    GraphIterator tmp(*this);
                    this->operator--();
                    return tmp;
                }

                inline bool
                operator==(const GraphIterator &x) const
                {
                    return (this->int_it == x.int_it);
                }

                inline bool
                operator!=(const GraphIterator &x) const
                {
                    return !(this->operator==(x));
                }

                /*! \brief checks if (this) iterator refers to a given graph container.
                 *
                 * \param   g   graph container in question.
                 * \return      true iff iterator refers to m. always returns false for explicitly
                 *              invalidated iterators.
                 * */
                inline bool
                checkContainer(Graph<Tg, Tv, Te> const &g) const
                {
                    return (this->graph == &g);
                }

                /*! \brief checks if (this) iterator refers to the same graph container as given
                 * iterator.
                 *
                 * \param   x   input iterator.
                 * \return      true iff both iterators (this) and x refer to the same graph
                 *              container. returns false if any iterator has been explicitly
                 *              invalidated.
                 * */
                inline bool
                sameContainer(const GraphIterator &x) const
                {
                    return (!this->explicitlyInvalid() && !x.explicitlyInvalid() && this->graph == x.graph);
                }

                /*! \brief explicitly invalidates (this) iterator.
                 *
                 * \sideeffect (this) iterator is permanently invalidated.
                 *
                 * \par foobar\n\n lbla
                 *
                 * \note all other iterators referring to the same graph component remain unaffected.
                 * it is therefore generally possible that other iterators referring to the same
                 * object are not marked as explicitly invalid.
                 *
                 * \note it is the client's responsibility to maintain information about iterator
                 * validity.  this method does not absolve the client from that responsibility.
                 *
                 * */
                inline void
                explicitlyInvalidate()
                {
                    this->graph = NULL;
                }

                /*! \brief check if (this) iterator has been explicitly invalidated.
                 *
                 * \returns true iff (this) iterator has been explicitly invalidated.
                 *
                 * \note this method must not be used to check for iterator validity in the
                 * general case: an iterator can in general be invalid although this method 
                 * return false. using the method in this way will result in undefined behaviour.
                 * the main use case: an algorithm or method receives an iterator (or
                 * a container with iterators) per reference, manipulates the graph(es) the passed
                 * iterators refer to, and updates the iterators accordingly.
                 * for example, a vertex might be deleted, which can be reflected by explicitly
                 * invalidating the corresponding iterator.
                 *
                 * \note it is the client's responsibility to maintain information about iterator validity.
                 * this method does not absolve the client from that responsibility.
                 *
                 * */
                inline bool
                explicitlyInvalid() const
                {
                    return (this->graph == NULL);
                }

                inline ValueType &
                operator*() const
                {
                    return *(ValueType::getPtr(this->int_it));
                }

                inline ValueType *
                operator->() const
                {
                    return (ValueType::getPtr(this->int_it));
                }
        };
        /* forward declaration of classes Vertex and Edge */
        class   Vertex;
        class   Edge;

    protected:
        /* typedefs for pointer typed internally used. could also be std::shared_ptr<..> with minor modifications */
        typedef Graph<Tg, Tv, Te>::Vertex *      VertexPointerType;    
        typedef Graph<Tg, Tv, Te>::Edge *        EdgePointerType;    

    public:
        /* NOTE: typdefs are not treated as full types by either the standard or compilers. for
         * more type safety, declare vertex_iterator as template specialization of GraphIterator
         * instead of typedef */
        /*
        typedef
            GraphIterator<
                Vertex,
                std::map<uint32_t, VertexPointerType >
            >   vertex_iterator;
        */

        class   vertex_iterator :
            public GraphIterator<
                Graph<Tg, Tv, Te>::Vertex,
                typename std::map<uint32_t, VertexPointerType>::iterator
            >
        {
            public:
                vertex_iterator()
                {
                }

                vertex_iterator(
                    Graph<Tg, Tv, Te>                                          *g,
                    typename std::map<uint32_t, VertexPointerType >::iterator   it)
                {
                    this->graph     = g;
                    this->int_it    = it;
                    //this->val       = Vertex::getPtr(it);
                }

                vertex_iterator(const vertex_iterator &x)
                    : GraphIterator<
                        Graph<Tg, Tv, Te>::Vertex,
                        typename std::map<uint32_t, VertexPointerType>::iterator
                      >()
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                }

               ~vertex_iterator()
               {
               }
        };

        class   vertex_const_iterator :
            public GraphIterator<
                const Graph<Tg, Tv, Te>::Vertex,
                typename std::map<uint32_t, VertexPointerType >::const_iterator
            >
        {
            public:
                vertex_const_iterator()
                {
                }

                vertex_const_iterator(
                    Graph<Tg, Tv, Te>                                                  *g,
                    typename std::map<uint32_t, VertexPointerType >::const_iterator     it)
                {
                    this->graph     = g;
                    this->int_it    = it;
                    //this->val       = Vertex::getPtr(it);
                }
                    
                vertex_const_iterator(const vertex_const_iterator &x)
                    :  GraphIterator<
                        const Graph<Tg, Tv, Te>::Vertex,
                        typename std::map<uint32_t, VertexPointerType>::const_iterator
                    > ()
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

                /* provide copy constructor with non-const vertex_iterator argument to allow implicit
                 * conversion of non-const vertex_iterator to const vertex_const_iterator */
                vertex_const_iterator(const vertex_iterator &x)
                    :  GraphIterator<
                        const Graph<Tg, Tv, Te>::Vertex,
                        typename std::map<uint32_t, VertexPointerType>::const_iterator
                    > ()
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

               ~vertex_const_iterator()
               {
               }
        };

        /* NOTE: same as for vertex_iterator above. typedefs are not treated as full new types, but
         * rather as a sort macro, which is undesirable for type checking. fully specialize
         * GraphIterator template for const and non-const face iterators */
        /*
        typedef GraphIterator<
                Edge,
                std::map<uint32_t, EdgePointerType >
            > edge_iterator;
        */

        class   edge_iterator :
            public GraphIterator<
                Graph<Tg, Tv, Te>::Edge,
                typename std::map<uint32_t, EdgePointerType>::iterator
            >
        {
            public:
                edge_iterator()
                {
                }

                edge_iterator(
                    Graph<Tg, Tv, Te>                                          *g,
                    typename std::map<uint32_t, EdgePointerType >::iterator     it)
                {
                    this->graph     = g;
                    this->int_it    = it;
                    //this->val       = Vertex::getPtr(it);
                }

                edge_iterator(const edge_iterator &x)
                    :  GraphIterator<
                        Graph<Tg, Tv, Te>::Edge,
                        typename std::map<uint32_t, EdgePointerType>::iterator
                    > ()
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                }

               ~edge_iterator()
               {
               }
        };

        class edge_const_iterator :
            public GraphIterator<
                const Graph<Tg, Tv, Te>::Edge,
                typename std::map<uint32_t, EdgePointerType>::const_iterator
            >
        {
            public:
                edge_const_iterator()
                {
                }

                edge_const_iterator(
                    Graph<Tg, Tv, Te>                                              *g,
                    typename std::map<uint32_t, EdgePointerType >::const_iterator   it)
                {
                    this->graph     = g;
                    this->int_it    = it;
                    //this->val       = Edge::getPtr(it);
                }

                edge_const_iterator(const edge_const_iterator &x)
                    : GraphIterator<
                        const Graph<Tg, Tv, Te>::Edge,
                        typename std::map<uint32_t, EdgePointerType>::const_iterator
                    > ()
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

                /* provide copy constructor with non-const vertex_iterator argument to allow implicit
                 * conversion of non-const vertex_iterator to const vertex_const_iterator */
                edge_const_iterator(const edge_iterator &x)
                    : GraphIterator<
                        const Graph<Tg, Tv, Te>::Edge,
                        typename std::map<uint32_t, EdgePointerType>::const_iterator
                    > ()
                {
                    this->graph     = x.graph;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

               ~edge_const_iterator()
                {
                }
        };

        /* forward declaration of Graph::VertexAccessor */
        class VertexAccessor;

        /* vertex class */
        class Vertex {
            friend class Graph<Tg, Tv, Te>;
            /* NOTE: in C++11, all nested classes are friends of the containing class. furthermore,
             * they act as an implementation detail of the containing class and have access to all
             * data that the containing class has access to.
             * So: Graph::Vertex has access to the private data of Graph, but by default not vice
             * versa. Since Graph has been declared a friend of Graph::Vertex above, all nested
             * classes of Graph have access to Graph::Vertex as well, rendering the following 
             * delcarations obsolete */

            /* friend class Graph::VertexAccessor; */
            /* can't declare friend with typedef vertex_iterator, so write full template
             * argument list. logically, this is equivalent to:
             *
             * friend class vertex_iterator; */
            /*
            friend class Graph::GraphIterator<
                            Vertex,
                            std::map<uint32_t, Vertex>
                         >;
            */
            /* vertex_const_iterator is inherited from GraphIterator: this works */
            /* friend class vertex_const_iterator; */

            protected:
                Graph<Tg, Tv, Te>                  *graph;
                typename std::map<
                    uint32_t,
                    VertexPointerType >::iterator   g_vit;

                uint32_t                            blocked : 1,  current_traversal_id : 23, traversal_state : 8;

                std::list<Edge *>                   in_edges, out_edges;

                /* protected ctors */
                                                    Vertex();
                                                    Vertex(
                                                        Graph<Tg, Tv, Te> * const  &graph,
                                                        Tv const                   &data = Tv());

                /* protected copy ctor and assignment operator. to be used only internally and only with care. */
                                                    Vertex(const Vertex &x);
                Vertex                             &operator=(const Vertex &x);

                virtual                            ~Vertex();

                /* Vertex objects must not be publically allocated with new() or delete() => protected
                 * new and delete operators to prevent this at compile time. */
                static void                        *operator new(size_t size);
                static void                         operator delete(void *p);

                /* other private methods */
                void                                replaceAdjacentVertices(const std::map<Vertex *, Vertex*> &replace_map);

            public:
                Tv                                  vertex_data;

                /* static getPtr() method required by iterator */
                static Vertex *                     getPtr(typename std::map<uint32_t, VertexPointerType >::const_iterator const &it);

                uint32_t
                id() const
                {
                    return this->g_vit->first;
                }


                vertex_iterator
                iterator() const
                {
                    return vertex_iterator(this->graph, this->g_vit);
                }

                vertex_const_iterator
                const_iterator() const
                {
                    return vertex_iterator(this->graph, this->g_vit);
                }

                uint32_t
                indeg() const
                {
                    return (this->getInNeighbours().size());
                }

                uint32_t
                outdeg() const
                {
                    return (this->getOutNeighbours().size());
                }

                bool
                gotOutNeighbour(const vertex_const_iterator &v) const
                {
                    return this->gotOutNeighbour(&(*v));
                }

                bool
                gotOutNeighbour(const Vertex * const &v) const
                {
                    for (auto &e : this->out_edges) {
                        if (v->id() == e->getDestinationVertex()->id()) {
                            return true;
                        }
                    }
                    return false;
                }

                bool
                gotOutNeighbour(const uint32_t &v_id) const
                {
                    Graph<Tg, Tv, Te> const &const_graph    = (*this->graph);
                    vertex_const_iterator vit               = this->graph->vertices.find(v_id);

                    if (vit != const_graph.vertices.end()) {
                        return (this->gotOutNeighbour(vit));
                    }
                    else {
                        throw("Graph::Vertex::gotOutNeighbour(const uint32_t &): specified vertex id not found.");
                    }
                }

                bool
                gotInNeighbour(const vertex_const_iterator &v) const
                {
                    return this->gotInNeighbour(&(*v));
                }

                bool
                gotInNeighbour(const Vertex * const &v) const
                {
                    for (auto &e : this->in_edges) {
                        if (v->id() == e->getSourceVertex()->id()) {
                            return true;
                        }
                    }
                    return false;
                }

                bool
                gotInNeighbour(const uint32_t &v_id) const
                {
                    Graph<Tg, Tv, Te> const &const_graph    = (*this->graph);
                    vertex_const_iterator vit               = this->graph->vertices.find(v_id);

                    if (vit != const_graph.vertices.end()) {
                        return (this->gotInNeighbour(vit));
                    }
                    else {
                        throw("Graph::Vertex::gotInNeighbour(const uint32_t &): specified vertex id not found.");
                    }
                }

                Tv &
                operator*()
                {
                    return (this->vertex_data);
                }

                Tv const &
                operator*() const
                {
                    return (this->vertex_data);
                }

                Tv *
                operator->()
                {
                    return &(this->vertex_data);
                }
                
                Tv const *
                operator->() const
                {
                    return &(this->vertex_data);
                }

                std::list<const Vertex *>           getOutNeighbours() const;
                void                                getOutNeighbours(std::list<Vertex *> &out_nbs) const;
                void                                getOutNeighbours(std::list<const Vertex *> &out_nbs) const;
                void                                getOutNeighboursIndices(std::list<uint32_t> &out_nbs) const;

                std::list<const Vertex *>           getInNeighbours() const;
                void                                getInNeighbours(std::list<Vertex *> &in_nbs) const;
                void                                getInNeighbours(std::list<const Vertex *> &in_nbs) const;
                void                                getInNeighboursIndices(std::list<vertex_iterator> &out_nbs) const;

                std::list<const Edge *>             getOutEdges() const;
                void                                getOutEdges(std::list<Edge *> &out_edges) const;
                void                                getOutEdges(std::list<const Edge *> &out_edges) const;

                std::list<const Edge *>             getInEdges() const;
                void                                getInEdges(std::list<Edge *> &in_edges) const;
                void                                getInEdges(std::list<const Edge *> &in_edges) const;
                
                uint8_t                             getTraversalState(const uint32_t &traversal_id);
                void                                setTraversalState(const uint32_t &traversal_id, const uint8_t &state);
        };

        /* forward declaration of EdgeAccessor */
        class EdgeAccessor;

    public: 
        /* face class */
        class Edge {
            friend class Graph<Tg, Tv, Te>;
            /* NOTE: in C++11, all nested classes are friends of the containing class. furthermore,
             * they act as an implementation detail of the containing class and have access to all
             * data that the containing class has access to.
             * So: Graph::Vertex has access to the private data of Graph, but by default not vice
             * versa. Since Graph has been declared a friend of Graph::Vertex above, all nested
             * classes of Graph have access to Graph::Vertex as well, rendering the following 
             * delcarations obsolete */

            /* friend class Graph::EdgeAccessor; */
            /* C++ can't declare friend with typedef edge_iterator, so write full template
             * argument list. logically, this is equivalent to
             *
             * friend class edge_iterator; */
            /*
            friend class Graph::GraphIterator<
                            Edge,
                            std::map<uint32_t, Edge>
                         >; */
            /* edge_const_iterator is inherited from GraphIterator: this works */
            /* friend class edge_const_iterator; */

            protected:
                Graph<Tg, Tv, Te>                  *graph;
                typename std::map<
                    uint32_t,
                    EdgePointerType
                >::iterator                         g_eit;

                Vertex                             *v_src, *v_dst;
                uint32_t                            blocked : 1,  current_traversal_id : 23, traversal_state : 8;

                /* private ctors, can only be called by Graph and nested classes */
                                                    Edge(); 
                                                    Edge(
                                                        Graph<Tg, Tv, Te>  *graph,
                                                        Vertex             *v_src,
                                                        Vertex             *v_dst,
                                                        Te const           &data = Te(),
                                                        bool const          blocked = false);

                                                    Edge(const Edge &x);
                Edge                               &operator=(const Edge &b);
                virtual                            ~Edge();

                void                                replaceVertices(const std::map<Vertex *, Vertex*> &replace_map);
                bool                                operator<(const Edge &b) const;
                bool                                operator>(const Edge &b) const;

            public:
                Te                                  edge_data;

                /* static getPtr() method required by iterator */
                static Edge *                       getPtr(typename std::map<uint32_t, EdgePointerType >::const_iterator const &it);

                uint32_t                            id() const;
                edge_iterator                       iterator() const;
                edge_const_iterator                 const_iterator() const;

                bool                                isBlocked() const;

                Te                                 &operator*();
                Te const                           &operator*() const;
                Te                                 *operator->();
                Te const                           *operator->() const;

                bool                                contains(vertex_const_iterator const &v) const;
                bool                                contains(const Vertex * const &v) const;
                bool                                contains(uint32_t const &v_id) const;

                bool                                gotNeighbour(edge_const_iterator const &e) const;
                bool                                gotNeighbour(const Edge * const &e) const;
                bool                                gotNeighbour(uint32_t const &e_id) const;

                uint8_t                             getTraversalState(uint32_t const &traversal_id);
                void                                setTraversalState(uint32_t const &traversal_id, const uint8_t &state);

                vertex_iterator                     getSourceVertex();
                vertex_const_iterator               getSourceVertex() const;

                vertex_iterator                     getDestinationVertex();
                vertex_const_iterator               getDestinationVertex() const;
        };

        /* ------------------ Accessor classes for vertices / faces ------------- */
        class VertexAccessor {
            friend class Graph<Tg, Tv, Te>;

            private:
                Graph<Tg, Tv, Te>                   &graph;

                /* private constructor, can only be called by Graph */
                                                    VertexAccessor(Graph<Tg, Tv, Te> &g);

                /* not assignable, not copy-constructible */
                                                    VertexAccessor(const VertexAccessor &x) = delete;
                VertexAccessor                     &operator=(const VertexAccessor &x)      = delete;

            public:
                vertex_iterator                     begin();
                vertex_const_iterator               begin() const;

                vertex_iterator                     end();
                vertex_const_iterator               end() const;

                vertex_iterator                     find(uint32_t id);
                vertex_const_iterator               find(uint32_t id) const;

                Graph::Vertex                      &at(uint32_t id);
                const Graph::Vertex                &at(uint32_t id) const;

                bool                                exists(uint32_t id) const;

                vertex_iterator                     insert(Tv const &data);

                vertex_iterator                     erase(vertex_iterator it);
                bool                                erase(uint32_t id);

                size_t                              size() const;
        };

        class EdgeAccessor {
            friend class Graph<Tg, Tv, Te>;

            private:
                Graph<Tg, Tv, Te>                  &graph;

                /* private constructor, can only be called by Graph */
                                                    EdgeAccessor(Graph<Tg, Tv, Te> &g);
                /* not assignable, not copy-constructible */
                EdgeAccessor                       &operator=(EdgeAccessor const &x)    = delete;
                                                    EdgeAccessor(EdgeAccessor const &x) = delete;

                /* private insertion method that can be used to insert an allocated Edge object already containing its
                 * data. this method cannot be used by clients, but it is useful to insert previously allocated derived
                 * Edge class objects, which have been upcast to a pointer to Edge for this purpose (as in
                 * CellNetwork).
                 *
                 * all topological information inside the given vertex v is cleared (in_edges / out_edges). */
                std::pair<edge_iterator, bool>      privateInsert(Edge *e);


            public:
                edge_iterator                       begin();
                edge_const_iterator                 begin() const;

                edge_iterator                       end();
                edge_const_iterator                 end() const;

                edge_iterator                       find(uint32_t id);
                edge_const_iterator                 find(uint32_t id) const;

                edge_iterator                       find(
                                                        vertex_const_iterator const    &v_src_it,
                                                        vertex_const_iterator const    &v_dst_it);
                edge_const_iterator                 find(
                                                        vertex_const_iterator const    &v_src_it,
                                                        vertex_const_iterator const    &v_dst_it) const;

                Graph::Edge                        &at(uint32_t id);
                Graph::Edge const                  &at(uint32_t id) const;

                bool                                exists(uint32_t id) const;
                bool                                exists(
                                                        vertex_const_iterator const    &v_src_it,
                                                        vertex_const_iterator const    &v_dst_it) const;


                std::pair<edge_iterator, bool>      insert(
                                                        vertex_iterator     v_src_it,
                                                        vertex_iterator     v_dst_it,
                                                        Te const           &data = Te());

                std::pair<edge_iterator, bool>      insert(
                                                        uint32_t    v_src_id,
                                                        uint32_t    v_dst_id,
                                                        Te const    &data = Te());
                                                        
                edge_iterator                       erase(edge_iterator it);
                bool                                erase(uint32_t id);

                vertex_iterator                     collapse(edge_iterator it);

                size_t                              size() const;
        };


    protected:
        /* id queues */
        IdQueue                             V_idq;
        IdQueue                             E_idq;
        IdQueue                             traversal_idq;

        /* data object of template type Tg */
        Tg                                  graph_data;

        /* vertex and face maps */
        std::map<
                uint32_t,
                VertexPointerType
            >                               V;
        std::map<
                uint32_t,
                EdgePointerType
            >                               E;

        /* globally reset the traversal states of all vertices and faces to TRAV_UNSEEN and reset
         * the traversal id queue. */
        void                                resetTraversalStates();

        /* protected insertion method that can be used to insert an allocated {Vertex, Edge} object already containing its
         * data. this method cannot be used by clients, but it is useful to insert previously allocated derived
         * {Vertex, Edge} class objects, which have been upcast to a pointer to {Vertex, Edge} for this purpose (as in
         * CellNetwork).
         *
         * all topological information inside the given vertex are erased (incident edges), the edge insertion method
         * takes two ids along with the pointer to verify that the allocated edge is indeed topologically intact. */
        std::pair<
                typename std::map<uint32_t, VertexPointerType>::iterator,
                bool
            >                               protectedVertexInsert(Vertex *v);

        std::pair<
                typename std::map<uint32_t, EdgePointerType>::iterator,
                bool
            >                               protectedEdgeInsert(
                                                Edge       *e,
                                                uint32_t    v_src_id,
                                                uint32_t    v_dst_id);

        /* protected static methods to extract internal std::map<..,..> iterators from Graph::{vertex,edge}_iterators */
        static typename std::map<
                uint32_t,
                VertexPointerType
            >::const_iterator               getInternalIterator(vertex_const_iterator const &it);

        static typename std::map<
                uint32_t,
                VertexPointerType
            >::iterator                     getInternalIterator(vertex_iterator const &it);

        static typename std::map<
                uint32_t,
                EdgePointerType
            >::const_iterator               getInternalIterator(edge_const_iterator const &it);

        static typename std::map<
                uint32_t,
                EdgePointerType
            >::iterator                     getInternalIterator(edge_iterator const &it);

    public:
        /* Graph public interface */

        /* the only two publically accessible members are immutable instances of the accessor classes for vertices and
         * faces. they are both non-copy-constructible and non-assignable. */
        VertexAccessor                      vertices;
        EdgeAccessor                        edges;

        /* ctor, copy ctor */
                                            Graph(Tg const &graph_data = Tg());
                                            Graph(Graph<Tg, Tv, Te> const &X);
        /* assignment operator */                                    
        Graph<Tg, Tv, Te>                  &operator=(Graph<Tg, Tv, Te> const &X);

        /* dtor */
                                           ~Graph();

        /* data object of template type Tg */
        Tg                                 &data();
        Tg const                           &data() const;


        /* clear all data, clear faces only */
        void                                clear();
        void                                clearEdges();
        /* renumber vertices and edges consecutively from the given start ids onwards */
        void                                renumberConsecutively(uint32_t vertex_start_id = 0, uint32_t face_start_id = 0);

        /* general info / statistics */
        uint32_t                            numVertices() const;
        uint32_t                            numEdges() const;

        /* traveral state is publically accessible for traversal methods. */
        enum GraphTraversalStates {
            TRAV_UNSEEN      = 0,
            TRAV_ENQUEUED    = 1,
            TRAV_DONE        = 2,
            TRAV_BLOCKED     = 3
        };

        uint32_t                            getFreshTraversalId();

        /* append another graph: add vertices / faces and offset indices. not that no topological connection
         * between (this) graph and the appended graph is performed, it simply computes the union of two
         * distinct graphes inside one object. note also that the old ids (and of course all
         * iterators) to B are NO LONGER valid in the appended graph, i.e. ids from B cannot be used
         * to refer to the copy of B that has been appended to graph in question. */
        void                                copyAppend(Graph<Tg, Tv, Te> const &B);

        /* move version: append graph B to (this) graph by moving B's contents. B is empty afterwards.
         * since pointers are used internaly, this does not require significant additional memory
         * and can be done "in-place".  the parameter update_vits is a pointer to a list containing
         * vertex iterators referring to B, which will be updated in-place with iterators referring
         * to the respective vertices of (this) graph after the append. */
        void                                moveAppend(
                                                Graph<Tg, Tv, Te>              &B,
                                                std::list<vertex_iterator>     *update_vits = NULL);

        /* delete the connected component (i.e. all reachable vertices / faces) of the vertex identified by iterator
         * vstart_it */
        void                                deleteConnectedComponent(vertex_iterator vstart_it);

        /* get all vertices / faces / edges of the connected component of a given start vertex for a
         * given traversal id. note that "boundary conditions" can be implemeneted by manually
         * setting the traversal state of selected vertices for the supplied traversal id before
         * calling this method. example: given a "ring" of vertices whose deletion separates the
         * graph in two connected components, one can set the traversal state of all these ring
         * vertices to DONE, which causes them to act as a boundary that the traversal cannot pass
         * */
        void                                getConnectedComponentBreadthFirst(
                                                Graph::vertex_iterator      vstart_it,    
                                                const uint32_t             &traversal_id,
                                                std::list<Vertex *>        *cc_vertices = NULL,
                                                std::list<Edge *>          *cc_faces    = NULL);

        /* ----------- EDGE related methods, which don't have their own Iterator / Accessor scheme
         * (yet) ------- */
        void                                checkEdge(
                                                const std::string&              fn,
                                                const vertex_const_iterator    &u_it,
                                                const vertex_const_iterator    &v_it) const;

        void                                checkEdge(
                                                const char                     *fn,
                                                const vertex_const_iterator    &u_it,
                                                const vertex_const_iterator    &v_it) const;

        /* ----------------- consistency check method  ----------------- */
        /* debug method: check internal data consistency */
        void                                checkInternalConsistency() const;

        /* -----------------  I/O  ----------------- */
        /* graph format would be nice.. serializer for usually data?
        void                                readFromObjFile(const char *filename);
        void                                writeObjFile(const char *jobname, bool include_materials = false);
        */

        /* NOTE: In the C++11 standard, nested classes are automatically "friends" of the containing
         * class, but not vice versa. the declarations below are therefore obsolete */
        /* "friend class vertex_iterator;" cannot be used since C++ does not allow friend
         * declaration through typdefs. use full template argument list */
        /*
        friend class Graph::VertexAccessor;
        friend class Graph::EdgeAccessor;
        */
        /*
        friend class Graph::GraphIterator<
                        Vertex,
                        std::map<uint32_t, Vertex>
                     >;
        friend class vertex_const_iterator;
        */

        /* "friend class edge_iterator;" cannot be used since C++ does not allow friend
         * declaration through typdefs. use full template argument list */
        /*
        friend class Graph::GraphIterator<
                        Edge,
                        std::map<uint32_t, Edge>
                     >;
        friend class edge_const_iterator;
        */
};

/* include header for template implementation */
#include  "../tsrc/Graph_impl.hh"

#endif
