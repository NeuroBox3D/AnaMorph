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

#ifndef CELL_NETWORK_H
#define CELL_NETWORK_H

#include "common.hh"
#include "Graph.hh"

template <
    typename Tn,
    typename Tv,
    typename Te,
    typename Tso   = Common::UnitType,
    typename Tnv    = Common::UnitType,
    typename Tax    = Common::UnitType,
    typename Tde    = Common::UnitType,
    typename Tns    = Common::UnitType,
    typename Tas    = Common::UnitType,
    typename Tds    = Common::UnitType,
    typename Tnr    = Common::UnitType,
    typename Tar    = Common::UnitType,
    typename Tdr    = Common::UnitType,
    // typename Tsyn = Common::UnitType
    typename R      = double
>
class CellNetwork : public Graph<Tn, Tv, Te> {
    /* public data types */
    public:
        /* derived vertex type hierarchy */
        enum NeuronVertexTypes {
            SOMA_VERTEX,
            AXON_VERTEX,
            DENDRITE_VERTEX
        };

        class CellSection {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            private:
                uint32_t            _compartment_id;
                Vec3<R>             _pos;
                R                   _radius;

            public:
                                    CellSection() : _compartment_id(0), _pos(), _radius(0)
                                    {
                                    }

                                    CellSection(
                                        uint32_t        compartment_id,
                                        Vec3<R> const   &pos,
                                        R const        &radius)
                                    : _compartment_id(compartment_id), _pos(pos), _radius(radius)
                                    {}
            /* NOTE: implicit copy constructor and assignment operator suffice here, as there are no non-trivial
             * members.*/
                uint32_t           &compartment_id()
                                    {
                                        return this->_compartment_id;
                                    }

                uint32_t const     &compartment_id() const
                                    {
                                        return this->_compartment_id;
                                    }

                Vec3<R>            &position()
                                    {
                                        return (this->_pos);
                                    }

                Vec3<R> const      &position() const
                                    {
                                        return (this->_pos);
                                    }

                R                  &radius()
                                    {
                                        return (this->_radius);
                                    }

                R const            &radius() const
                                    {
                                        return (this->_radius);
                                    }
        };

    private:
        /* forward declarations of predicates */
        struct NeuronVertexPred;
        struct NeuriteVertexPred;
        struct SomaVertexPred;
        struct AxonVertexPred;
        struct DendriteVertexPred;

    public:

        /* forward declarations */
        class NeuronVertexAccessor;
        class neuron_iterator;
        class neuron_const_iterator;

        class NeuronEdge;

        /* just a for Graph::Vertex in Graph, this class cannot be allocated or freed, assigned to or copied by a client. access to
         * objects of this type can only be provided by CellNetwork. */
        class NeuronVertex : public Graph<Tn, Tv, Te>::Vertex {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            friend struct NeuronVertexPred;
            friend struct NeuriteVertexPred;
            friend struct SomaVertexPred;
            friend struct AxonVertexPred;
            friend struct DendriteVertexPred;

            protected:
                /* pointer to cell network this vertex belongs to. upcast version of this is contained in
                 * Graph::Vertex::graph. */
                CellNetwork                *network;

                /* protected members */
                std::list<CellSection>      sections;

            protected:
                /* protected constructor */
                                                NeuronVertex(
                                                    CellNetwork * const                &network,
                                                    std::list<CellSection> const       &sections,
                                                    Tv const                           &data = Tv());

                /* copy constructor and assignment operator, to be used only internally and only with care, since it is
                 * at least questionable semantically to "assign" to vertices in a topological object such as a graph:
                 * what happens to the connectivity if a vertex is assigned another? the data inside the vertex is
                 * assignable, the vertex itself is handled only internally. */
                                                NeuronVertex(NeuronVertex const &v);
                NeuronVertex                   &operator=(NeuronVertex const &v);

                /* virtual destructor for polymorphism */
                virtual                        ~NeuronVertex();

                virtual uint32_t                getType() const = 0;
                /* virtual clone method that creates derived vertex object and returns upcast pointer to abstract class
                 * NeuronVertex*/
                virtual NeuronVertex           *clone() const = 0;

            public:
                virtual uint32_t                getNetworkCompartmentId() const;

                virtual Vec3<R>                 getSinglePointPosition() const = 0;
                virtual R                       getSinglePointRadius() const = 0;
                virtual void                    scaleSinglePointRadius(R scale) = 0;

                uint32_t                        numSections() const;

                R                               getSectionMinRadius() const;
                R                               getSectionMaxRadius() const;
                Vec3<R>                         getSectionPositionCentroid() const;
                virtual void                    getBoundingBox(
                                                    Vec3<R> &bb_min,
                                                    Vec3<R> &bb_max) const;

                std::list<CellSection>          getSections() const;
                //void                            setSections(std::list<CellSection> const &sections);

                void                            translate(Vec3<R> const &d);
                void                            scale(R const &x);

                bool                            isBranchingVertex() const;
                bool                            isTerminalVertex() const;
                bool                            isSimpleVertex() const;

                /* override Graph::Vertex methods for retrieving neighbours: in a CellNetwork, all (relevant) neighours
                 * are Neuron{Vertex,Edge} objects */
                std::list<NeuronVertex const *> getOutNeighbours() const;
                void                            getOutNeighbours(std::list<NeuronVertex *> &out_nbs) const;
                void                            getOutNeighbours(std::list<const NeuronVertex *> &out_nbs) const;

                std::list<NeuronVertex const *> getInNeighbours() const;
                void                            getInNeighbours(std::list<NeuronVertex *> &in_nbs) const;
                void                            getInNeighbours(std::list<const NeuronVertex *> &in_nbs) const;

                std::list<NeuronEdge const *>   getOutEdges() const;
                void                            getOutEdges(std::list<NeuronEdge *> &out_edges) const;
                void                            getOutEdges(std::list<const NeuronEdge *> &out_edges) const;

                std::list<NeuronEdge const *>   getInEdges() const;
                void                            getInEdges(std::list<NeuronEdge *> &in_edges) const;
                void                            getInEdges(std::list<const NeuronEdge *> &in_edges) const;

                /* templated neighbour access method versions that allow comfortable and safe (albeit not statically
                 * type-safe) access to the {in,out}-neighbours of a specific type, e.g. all SomaVertex or AxonVertex
                 * neighbours */
                template <typename NbType>  
                std::list<NbType const *>       getFilteredOutNeighbours() const;
                template <typename NbType>      
                void                            getFilteredOutNeighbours(std::list<NbType const *> &out_nbs) const;
                template <typename NbType>      
                void                            getFilteredOutNeighbours(std::list<NbType *> &out_nbs) const;

                template <typename NbType>      
                std::list<NbType const *>       getFilteredInNeighbours() const;
                template <typename NbType>      
                void                            getFilteredInNeighbours(std::list<NbType const *> &in_nbs) const;
                template <typename NbType>      
                void                            getFilteredInNeighbours(std::list<NbType *> &in_nbs) const;

                template <typename NbType>      
                std::list<NbType const *>       getFilteredNeighbours() const;

                /* same for edge types */
                template <typename EdgeType>  
                std::list<EdgeType const *>     getFilteredOutEdges() const;
                template <typename EdgeType>      
                void                            getFilteredOutEdges(std::list<EdgeType const *> &out_nbs) const;
                template <typename EdgeType>      
                void                            getFilteredOutEdges(std::list<EdgeType *> &out_nbs) const;

                template <typename EdgeType>      
                std::list<EdgeType const *>     getFilteredInEdges() const;
                template <typename EdgeType>      
                void                            getFilteredInEdges(std::list<EdgeType const *> &in_nbs) const;
                template <typename EdgeType>      
                void                            getFilteredInEdges(std::list<EdgeType *> &in_nbs) const;
            
                /* iterators */
                neuron_iterator                 iterator();
                neuron_const_iterator           iterator() const;
        };

        class SomaAccessor;
        class soma_iterator;
        class soma_const_iterator;

        class SomaVertex : public NeuronVertex {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            friend class SomaAccessor;

            private:
                /* private constructors */
                                            SomaVertex(
                                                CellNetwork * const                                    &network,
                                                Graph<uint32_t, CellSection, Common::UnitType> const   &section_graph,
                                                Tv const                                               &vertex_data = Tv(),
                                                Tso const                                              &soma_data = Tso());

                                            SomaVertex(const SomaVertex &s);
                SomaVertex                 &operator=(const SomaVertex &s);

                /* vertex objects must not be publically allocated with new() or delete() => private new and delete
                 * operators to prevent this at compile time. */
                static void                *operator new(size_t size);
                static void                 operator delete(void *p);

                uint32_t                    getType() const;
                virtual NeuronVertex       *clone() const final;

                /* soma encoding has special semantics: input data is parsed as-is and stored as a graph whose vertices
                 * are CellSections. the graph data is the id of the CellSection stored in the root node, i.e. the node
                 * with parent_id -1 for SWC files.
                 *
                 * NOTE: SomaVertex also contains the CellSection list NeuronVertex::sections, which must always be
                 * synchronized with the graph representation used here. this could be achieved through a
                 * reimplementation of all NeuronVertex methods. instead, the graph vertices (i.e. CellSections) are
                 * dumped into NeuronVertex::sections. no write access is possible at the NeuronVertex level, so this is
                 * fine for the moment. 
                 *
                 * NOTE: soma encoding is semantically ill-defined and ambiguous for SWC files. contact NeuroMorpho
                 * admins for more information about encoding / discussion of special extension for soma encodings. */
                Graph<uint32_t, CellSection, Common::UnitType>  section_graph;

            public:
                /* public member storing special data annotated to somas only */
                Tso                         soma_data;

                virtual Vec3<R>             getSinglePointPosition() const final;
                virtual R                   getSinglePointRadius() const final;
                virtual void                scaleSinglePointRadius(R scale) final;

                /* iterators */
                soma_iterator               iterator();
                soma_const_iterator         iterator() const;
        };

        /* forward declarations */
        class NeuriteVertexAccessor;
        class neurite_iterator;
        class neurite_const_iterator;
        class neurite_rootedge_iterator;
        class neurite_segment_iterator;

        /* abstract class for neurite vertices, i.e. axon and dendrite vertices */
        class NeuriteVertex : public NeuronVertex {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            /* protected members */
            protected:
                soma_iterator               soma;
                neurite_rootedge_iterator   neurite;

            /* methods */
            protected:
                                            NeuriteVertex(
                                                CellNetwork * const    &network,
                                                CellSection const      &section,
                                                Tnv const              &neurite_vertex_data = Tnv(),
                                                Tv const               &data = Tv());

                                            NeuriteVertex(NeuriteVertex const &n);
                NeuriteVertex              &operator=(NeuriteVertex const &s);

                virtual uint32_t            getType() const = 0;
                virtual NeuronVertex       *clone() const = 0;

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tnv                         neurite_data;

                virtual Vec3<R>             getSinglePointPosition() const;
                virtual R                   getSinglePointRadius() const;
                virtual void                scaleSinglePointRadius(R scale);

                virtual Vec3<R>             getPosition() const;
                virtual void                setPosition(Vec3<R> const &p);

                virtual R                   getRadius() const;
                virtual void                setRadius(R const &r);

                soma_iterator               getSoma() const;
                neurite_rootedge_iterator   getNeurite() const;

                /* define edge type predicates restricted to cell-tree edges, i.e. neurite root edges and neurite
                 * segment. new predicate for neurite root vertices. */
                bool                        isNeuriteSimpleVertex() const;
                bool                        isNeuriteBranchingVertex() const;
                bool                        isNeuriteTerminalVertex() const;
                bool                        isNeuriteRootVertex() const;

                neurite_iterator            getNeuriteParent() const;
                neurite_iterator            getNeuriteFirstChild() const;
                neurite_segment_iterator    getNeuriteSegmentFirstChild() const;

                /* iterators */
                neurite_iterator            iterator();
                neurite_const_iterator      iterator() const;
        };

        /* forward declarations */
        class AxonAccessor;
        class axon_iterator;
        class axon_const_iterator;
        class axon_rootedge_iterator;

        class AxonVertex : public NeuriteVertex {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            friend class AxonAccessor;

            private:
                                                    AxonVertex(
                                                        CellNetwork * const    &network,
                                                        CellSection const      &section,
                                                        Tax const              &axon_data = Tax(),
                                                        Tnv const              &neurite_vertex_data = Tnv(),
                                                        Tv const               &vertex_data = Tv());

                                                    AxonVertex(AxonVertex const &n);
                AxonVertex                         &operator=(AxonVertex const &s);

                /* vertex objects must not be publicly allocated with new() or delete() => private new and delete
                 * operators to prevent this at compile time. */
                static void                        *operator new(size_t size);
                static void                         operator delete(void *p);

                uint32_t                            getType() const final;
                NeuronVertex                       *clone() const final;

                /* member storing special data annotated to axons only */
                //axon_rootedge_iterator              axon;

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tax                                 axon_data;

                axon_rootedge_iterator              getAxon() const;

                /* iterators */
                axon_iterator                       iterator();
                axon_const_iterator                 iterator() const;

        };

        /* forward declarations */
        class DendriteAccessor;
        class dendrite_iterator;
        class dendrite_const_iterator;
        class dendrite_rootedge_iterator;

        class DendriteVertex : public NeuriteVertex {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            friend class DendriteAccessor;

            private:
                bool                                apical_dendrite;

                                                    DendriteVertex(
                                                        CellNetwork * const    &network,
                                                        CellSection const      &section,
                                                        bool                    apical_dendrite,
                                                        Tde const              &dendrite_data = Tv(),
                                                        Tnv const              &neurite_vertex_data = Tnv(), 
                                                        Tv const               &vertex_data = Tv());

                                                    DendriteVertex(DendriteVertex const &n);
                DendriteVertex                     &operator=(DendriteVertex const &s);

                /* vertex objects must not be publicly allocated with new() or delete() => private new and delete
                 * operators to prevent this at compile time. */
                static void                        *operator new(size_t size);
                static void                         operator delete(void *p);

                uint32_t                            getType() const final;
                NeuronVertex                       *clone() const final;

                /* member storing special data annotated to dendrites only */
                //dendrite_rootedge_iterator          dendrite;

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tde                                 dendrite_data;

                bool                                isApicalDendrite() const;

                dendrite_rootedge_iterator          getDendrite() const;

                /* iterators */
                dendrite_iterator                   iterator();
                dendrite_const_iterator             iterator() const;
        };

        /* derived edge type hierarchy */
        enum NeuronEdgeTypes {
            AXON_ROOT_EDGE,
            DENDRITE_ROOT_EDGE,
            AXON_SEGMENT,
            DENDRITE_SEGMENT
        };

    private:
        /* forward delcarations of edge predicates */
        struct NeuronEdgePred;
        struct AxonRootEdgePred;
        struct DendriteRootEdgePred;
        struct NeuriteSegmentPred;
        struct AxonSegmentPred;
        struct DendriteSegmentPred;

    public:

        class NeuronEdgeAccessor;
        class neuron_edge_iterator;
        class neuron_edge_const_iterator;

        class NeuronEdge : public Graph<Tn, Tv, Te>::Edge {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            friend struct NeuronEdgePred;
            friend struct AxonRootEdgePred;
            friend struct DendriteRootEdgePred;
            friend struct NeuriteSegmentPred;
            friend struct AxonSegmentPred;
            friend struct DendriteSegmentPred;

            friend class NeuronEdgeAccessor;

            protected:
                /* pointer to cell network this vertex belongs to. upcast version of this is contained in
                 * Graph::Edge::graph. */
                CellNetwork                *network;

                /* NeuronVertex pointers to source and destination vertices. upcast versions of these are contained in
                 * Graph::Edge::{v_src, v_dst} */
                NeuronVertex               *v_neuron_src, *v_neuron_dst;

            protected:
                /* protected constructor */
                                            NeuronEdge(
                                                CellNetwork * const    &network,
                                                NeuronVertex           *v_src,
                                                NeuronVertex           *v_dst,
                                                Te const               &edge_data = Te());

                /* copy constructor and assignment operator, to be used only internally and only with care, see
                 * NeuronVertex copy ctor */
                                            NeuronEdge(NeuronEdge const &e);
                NeuronEdge                 &operator=(NeuronEdge const &e);

                /* virtual destructor for polymorphism */
                virtual                    ~NeuronEdge();

                virtual uint32_t            getType() const = 0;

                /* virtual clone method */
                virtual NeuronEdge         *clone() const = 0;

                /* virtual update pointer method, which must be called to update down-cast pointers if the edge is
                 * changed on ANY level of abstraction */
                virtual void                updateVertexPointers();

            public:
                const CellNetwork           *getNetwork() const;
                /* get edge of length for traversals. this might involve "simplifying" complex edge types and provide an
                 * approximate length (e.g. for axon,dendrite root edges and synapses). */
                virtual R                   getLength() const;

                /* override get{Source,Destination}Vertex() to return specialized iterators */
                neuron_iterator             getSourceVertex();
                neuron_const_iterator       getSourceVertex() const;

                neuron_iterator             getDestinationVertex();
                neuron_const_iterator       getDestinationVertex() const;

                neuron_edge_iterator        iterator();
                neuron_edge_const_iterator  iterator() const;
        };

        class NeuriteRootEdgeAccessor;
        //class neurite_rootedge_iterator;
        class neurite_rootedge_const_iterator;

        class NeuriteRootEdge : public NeuronEdge {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            friend class NeuriteRootEdgeAccessor;

            protected:
                SomaVertex                         *v_soma_src;
                NeuriteVertex                      *v_neurite_dst;

                /* protected constructor */
                                                    NeuriteRootEdge(
                                                        CellNetwork * const    &network,
                                                        SomaVertex             *v_src,
                                                        NeuriteVertex          *v_dst,
                                                        Tnr const              &neurite_root_edge_data = Tnr(),
                                                        Te const               &edge_data = Te());

                /* copy constructor and assignment operator, to be used only internally and only with care, see
                 * NeuronVertex copy ctor */
                                                    NeuriteRootEdge(NeuriteRootEdge const &v);
                NeuriteRootEdge                    &operator=(NeuriteRootEdge const &v);

                /* virtual destructor for polymorphism */
                                                   ~NeuriteRootEdge();

                virtual uint32_t                    getType() const = 0;

                /* virtual clone method */
                virtual NeuronEdge                 *clone() const = 0;

                /* virtual update pointer method, which must be called to update down-cast pointers if the edge is
                 * changed on ANY level of abstraction */
                virtual void                        updateVertexPointers();

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tnr                                 neurite_root_edge_data;

                //virtual R                           getLength() const final;

                /* override get{Source,Destination}Vertex() to return specialized iterators */
                soma_iterator                       getSourceVertex();
                soma_const_iterator                 getSourceVertex() const;

                neurite_iterator                    getDestinationVertex();
                neurite_const_iterator              getDestinationVertex() const;

                /* iterators */
                neurite_rootedge_iterator           iterator();
                neurite_rootedge_const_iterator     iterator() const;
        };
            
        class AxonRootEdgeAccessor;
        //class axon_rootedge_iterator;
        class axon_rootedge_const_iterator;

        class AxonRootEdge : public NeuriteRootEdge {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;
            friend class AxonRootEdgeAccessor;

            protected:
                AxonVertex                         *v_axon_dst;

                /* protected constructor */
                                                    AxonRootEdge(
                                                        CellNetwork * const    &network,
                                                        SomaVertex             *v_src,
                                                        AxonVertex             *v_dst,
                                                        Tar const              &axon_root_edge_data = Tar(),
                                                        Tnr const              &neurite_root_edge_data = Tnr(),
                                                        Te const               &edge_data = Te());

                /* copy constructor and assignment operator, to be used only internally and only with care, see
                 * NeuronVertex copy ctor */
                                                    AxonRootEdge(AxonRootEdge const &v);
                AxonRootEdge                       &operator=(AxonRootEdge const &v);

                /* virtual destructor for polymorphism */
                                                   ~AxonRootEdge();

                virtual uint32_t                    getType() const final;

                /* virtual clone method */
                virtual NeuronEdge                 *clone() const final;

                /* virtual update pointer method, which must be called to update down-cast pointers if the edge is
                 * changed on ANY level of abstraction */
                virtual void                        updateVertexPointers();

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tar                                 axon_root_edge_data;

                //virtual R                           getLength() const final;

                /* override getDestinationVertex() to return specialized iterators */
                axon_iterator                       getDestinationVertex();
                axon_const_iterator                 getDestinationVertex() const;

                /* iterators */
                axon_rootedge_iterator              iterator();
                axon_rootedge_const_iterator        iterator() const;
                
        };

        class DendriteRootEdgeAccessor;
        //class dendrite_rootedge_iterator;
        class dendrite_rootedge_const_iterator;

        class DendriteRootEdge : public NeuriteRootEdge {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;
            friend class DendriteRootEdgeAccessor;

            protected:
                DendriteVertex                     *v_dendrite_dst;

                /* protected constructor */
                                                    DendriteRootEdge(
                                                        CellNetwork * const    &network,
                                                        SomaVertex             *v_src,
                                                        DendriteVertex         *v_dst,
                                                        Tdr const              &dendrite_root_edge_data,
                                                        Tnr const              &neurite_root_edge_data,
                                                        Te const               &edge_data);

                /* copy constructor and assignment operator, to be used only internally and only with care, see
                 * NeuronVertex copy ctor */
                                                    DendriteRootEdge(DendriteRootEdge const &v);
                DendriteRootEdge                   &operator=(DendriteRootEdge const &v);

                /* virtual destructor for polymorphism */
                                                   ~DendriteRootEdge();

                virtual uint32_t                    getType() const final;

                /* virtual clone method */
                virtual NeuronEdge                 *clone() const final;

                /* virtual update pointer method, which must be called to update down-cast pointers if the edge is
                 * changed on ANY level of abstraction */
                virtual void                        updateVertexPointers();

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tdr                                 dendrite_root_edge_data;

                //virtual R                           getLength() const final;

                /* override getDestinationVertex() to return specialized iterators */
                dendrite_iterator                   getDestinationVertex();
                dendrite_const_iterator             getDestinationVertex() const;

                /* iterators */
                dendrite_rootedge_iterator          iterator();
                dendrite_rootedge_const_iterator    iterator() const;
                
        };

        class NeuriteSegmentAccessor;
        //class neurite_segment_iterator;
        class neurite_segment_const_iterator;

        class NeuriteSegment : public NeuronEdge {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;

            friend class NeuriteSegmentAccessor;

            protected:
                NeuriteVertex                  *v_neurite_src;
                NeuriteVertex                  *v_neurite_dst;

                /* protected constructor */
                                                NeuriteSegment(
                                                    CellNetwork * const    &network,
                                                    NeuriteVertex          *v_src,
                                                    NeuriteVertex          *v_dst,
                                                    Tns const              &neurite_segment_data,
                                                    Te const               &edge_data);

                /* copy constructor and assignment operator, to be used only internally and only with care, see
                 * NeuronVertex copy ctor */
                                                NeuriteSegment(NeuriteSegment const &v);
                NeuriteSegment                 &operator=(NeuriteSegment const &v);

                virtual uint32_t                getType() const = 0;

                /* virtual destructor for polymorphism */
                virtual                        ~NeuriteSegment();

                /* virtual clone method */
                virtual NeuronEdge             *clone() const = 0;

                /* virtual update pointer method, which must be called to update down-cast pointers if the edge is
                 * changed on ANY level of abstraction */
                virtual void                    updateVertexPointers();

                soma_iterator                   soma;
                neurite_rootedge_iterator       neurite;

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tns                             neurite_segment_data;

                /* get length of neurite segment */
                virtual R                       getLength() const final;
                R                               getAbsoluteRadiusDifference() const;
                R                               getRadiusRatio() const;
                R                               getMinRadius() const;
                R                               getMaxRadius() const;
                R                               getAngle() const;
                std::pair<R, R>                 getSMDVRadii() const;

                /* override get{Source,Destination}Vertex() to return neurite_segment_iterator */
                neurite_iterator                getSourceVertex();
                neurite_const_iterator          getSourceVertex() const;

                neurite_iterator                getDestinationVertex();
                neurite_const_iterator          getDestinationVertex() const;

                /* get soma / neurite root edge (representing the neurite) this neurite segment belongs to */
                soma_iterator                   getSoma() const;
                neurite_rootedge_iterator       getNeurite() const;

                /* iterators */
                neurite_segment_iterator        iterator();
                neurite_segment_const_iterator  iterator() const;
        };

        class AxonSegmentAccessor;
        class axon_segment_iterator;
        class axon_segment_const_iterator;

        class AxonSegment : public NeuriteSegment {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;
            friend class AxonSegmentAccessor;

            protected:
                AxonVertex                         *v_axon_src;
                AxonVertex                         *v_axon_dst;

                /* protected constructor */
                                                    AxonSegment(
                                                        CellNetwork * const    &network,
                                                        AxonVertex             *v_src,
                                                        AxonVertex             *v_dst,
                                                        Tas const              &axon_segment_data,
                                                        Tns const              &neurite_segment_data,
                                                        Te const               &edge_data);

                /* copy constructor and assignment operator, to be used only internally and only with care, see
                 * NeuronVertex copy ctor */
                                                    AxonSegment(AxonSegment const &v);
                AxonSegment                        &operator=(AxonSegment const &v);

                /* virtual destructor for polymorphism */
                                                   ~AxonSegment();

                virtual uint32_t                    getType() const final;

                /* virtual clone method */
                virtual NeuronEdge                 *clone() const final;

                /* virtual update pointer method, which must be called to update down-cast pointers if the edge is
                 * changed on ANY level of abstraction */
                virtual void                        updateVertexPointers();

                //axon_rootedge_iterator              axon;

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tas                                 axon_segment_data;

                axon_iterator                       getSourceVertex();
                axon_const_iterator                 getSourceVertex() const;

                axon_iterator                       getDestinationVertex();
                axon_const_iterator                 getDestinationVertex() const;

                axon_rootedge_iterator              getAxon() const;

                /* iterators */
                axon_segment_iterator               iterator();
                axon_segment_const_iterator         iterator() const;
        };

        class DendriteSegmentAccessor;
        class dendrite_segment_iterator;
        class dendrite_segment_const_iterator;

        class DendriteSegment : public NeuriteSegment {
            friend class CellNetwork<Tn, Tv, Te, Tso, Tnv, Tax, Tde, Tns, Tas, Tds, Tnr, Tar, Tdr, R>;
            friend class DendriteSegmentAccessor;

            protected:
                DendriteVertex                     *v_dendrite_src;
                DendriteVertex                     *v_dendrite_dst;

                /* protected constructor */
                                                    DendriteSegment(
                                                        CellNetwork * const    &network,
                                                        DendriteVertex         *v_src,
                                                        DendriteVertex         *v_dst,
                                                        Tds const              &dendrite_segment_data,
                                                        Tns const              &neurite_segment_data,
                                                        Te const               &edge_data);

                /* copy constructor and assignment operator, to be used only internally and only with care, see
                 * NeuronVertex copy ctor */
                                                    DendriteSegment(DendriteSegment const &v);
                DendriteSegment                    &operator=(DendriteSegment const &v);

                /* virtual destructor for polymorphism */
                virtual                            ~DendriteSegment();

                virtual uint32_t                    getType() const final;

                /* virtual clone method */
                virtual NeuronEdge                 *clone() const final;

                /* virtual update pointer method, which must be called to update down-cast pointers if the edge is
                 * changed on ANY level of abstraction */
                virtual void                        updateVertexPointers();

                //dendrite_rootedge_iterator          dendrite;

            public:
                /* public member storing annotated data specific to level of abstraction */
                Tds                                 dendrite_segment_data;

                dendrite_iterator                   getSourceVertex();
                dendrite_const_iterator             getSourceVertex() const;

                dendrite_iterator                   getDestinationVertex();
                dendrite_const_iterator             getDestinationVertex() const;

                dendrite_rootedge_iterator          getDendrite() const;

                /* iterators */
                dendrite_segment_iterator           iterator();
                dendrite_segment_const_iterator     iterator() const; 
        };

        /* once synapses are supported, they can likely be modelled as edges */

    /* predicate functors used to internally identify derived types from base type pointers / references in the above
     * vertex / edge hierarchy */
    private:
        /* predicates for vertex types */
        struct NeuronVertexPred {
            bool
            operator()(NeuronVertex const &v) const
            {
                return true;
            }
        } isNeuronVertex;

        struct NeuriteVertexPred {
            bool
            operator()(NeuronVertex const &v) const
            {
                /* call polymorph getType() method to extract derived type id */
                uint32_t vtype = v.getType();

                return (vtype == AXON_VERTEX || vtype == DENDRITE_VERTEX);
            }
        } isNeuriteVertex;

        struct SomaVertexPred {
            bool
            operator()(NeuronVertex const &v) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (v.getType() == SOMA_VERTEX);
            }
        } isSomaVertex;

        struct AxonVertexPred {
            bool
            operator()(NeuronVertex const &v) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (v.getType() == AXON_VERTEX);
            }
        } isAxonVertex;

        struct DendriteVertexPred {
            bool
            operator()(NeuronVertex const &v) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (v.getType() == DENDRITE_VERTEX);
            }
        } isDendriteVertex;

        /* predicates for edge types */
        struct NeuronEdgePred {
            bool
            operator()(NeuronEdge const &e) const
            {
                return true;
            }
        } isNeuronEdge;

        struct NeuriteRootEdgePred {
            bool
            operator()(NeuronEdge const &e) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (e.getType() == AXON_ROOT_EDGE || e.getType() == DENDRITE_ROOT_EDGE);
            }
        } isNeuriteRootEdge;

        struct AxonRootEdgePred {
            bool
            operator()(NeuronEdge const &e) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (e.getType() == AXON_ROOT_EDGE);
            }
        } isAxonRootEdge;

        struct DendriteRootEdgePred {
            bool
            operator()(NeuronEdge const &e) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (e.getType() == DENDRITE_ROOT_EDGE);
            }
        } isDendriteRootEdge;

        struct NeuriteSegmentPred {
            bool
            operator()(NeuronEdge const &e) const
            {
                /* call polymorph getType() method to extract derived type id */
                uint32_t etype = e.getType();

                return (etype == AXON_SEGMENT || etype == DENDRITE_SEGMENT);
            }
        } isNeuriteSegment;

        struct AxonSegmentPred {
            bool
            operator()(NeuronEdge const &e) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (e.getType() == AXON_SEGMENT);
            }
        } isAxonSegment;

        struct DendriteSegmentPred {
            bool
            operator()(NeuronEdge const &e) const
            {
                /* call polymorph getType() method to extract derived type id */
                return (e.getType() == DENDRITE_SEGMENT);
            }
        } isDendriteSegment;

    /* public declarations of iterators / accessors */
    public:
        /* forward declaration of CellNetworkAccessor */
        template<
            typename  ValueType,
            typename  BaseType,
            typename  DelegatorType,
            typename  InternalType,
            typename  ValuePred,
            typename  IteratorType,
            typename  ConstIteratorType,
            typename  DataType
        >
        class CellNetworkAccessor;

        /* generic iterator that is specialized and derived to generate specific iterator types and provide a safe-guard
         * for dynamic down-casting of Graph::Vertex or CellNetwork::NeuronVertex to more specific vertex types, e.g. to
         * SomaVertex. the down-casting is necessary since the graph stored Graph::Vertex pointers internally and at
         * this level, more specific static type information is lost. */
        template <
            typename    ValueType,
            typename    BaseType,
            typename    DelegatorType,
            typename    InternalType,
            typename    InternalIteratorType,
            typename    ValuePred
        >
        class CellNetworkIterator : public Graph<Tn, Tv, Te>::template GraphIterator<BaseType, InternalIteratorType> {
            friend class CellNetwork;
            template<
                typename    AValueType,
                typename    ABaseType,
                typename    ADelegatorType,
                typename    AInternalType,
                typename    AValuePred,
                typename    AIteratorType,
                typename    AConstIteratorType,
                typename    ADataType
            > friend class CellNetworkAccessor;

            protected:
                /* NOTE: CellNetworkIterator has inherited Graph<Tn, Tv, Te> *graph; and IteratorType::iterator int_it
                 * from Graph::GraphIterator */
                InternalType   *int_ds;   
                ValuePred      *pred;
                CellNetwork    *network;

            private:
                /* private ctors, can only be called by CellNetwork. */
                CellNetworkIterator()
            	: int_ds(NULL), pred(NULL), network(NULL)
                {
                    this->explicitlyInvalidate();
                }

                CellNetworkIterator(
                    CellNetwork            *network,
                    InternalType           *ds,
                    InternalIteratorType    it,
                    ValuePred * const      &p)
                        /* call Graph::GraphIterator constructor in initializer list */
                        : Graph<Tn, Tv, Te>::template GraphIterator<BaseType, InternalIteratorType>(network, it)
                {
                    this->int_ds    = ds;
                    this->pred      = p;
                    this->network   = network;
                }

            public:
                CellNetworkIterator(CellNetworkIterator const &x)
            	/* call Graph::GraphIterator constructor in initializer list */
                : Graph<Tn, Tv, Te>::template GraphIterator<BaseType, InternalIteratorType>(x),
				  int_ds(x.int_ds), pred(NULL), network(x.network)
                {}

                CellNetworkIterator &
                operator=(CellNetworkIterator const &x)
                {
                    /* call operator= of base class Graph::GraphIterator */
                    Graph<Tn, Tv, Te>::template GraphIterator<BaseType, InternalIteratorType>::operator=(x);
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;

                    return (*this);
                }

               ~CellNetworkIterator()
                {
                }


                /* increment internal iterator until type id matches */
                CellNetworkIterator &
                operator++()
                {
                    /* check if this->int_it is already this->int_ds->end(). throw in case */
                    if (this->int_it != this->int_ds->end()) {
                        DelegatorType *dp;

                        /* increment internal iterator at least once */
                        ++this->int_it;

                        /* while iterator hasn't reached end() */
                        while (this->int_it != this->int_ds->end()) {
                            /* down-cast to delegator type, throw if downcast fails. */
                            dp = dynamic_cast<DelegatorType *>(BaseType::getPtr(this->int_it));
                            if (dp) {
                                /* if type matches, return, otherwise increment the internal iterator
                                 * and restart the loop if end() hasn't been reached yet */
                                if ( this->pred->operator()(*dp) ) {
                                    break;
                                }
                                else {
                                    ++this->int_it;
                                }
                            }
                            else {
                                throw("CellNetwork::CellNetworkIterator::operator++(): failed to down-cast base type to delegator type. internal logic error.");
                            }
                        }
                        return (*this);
                    }
                    else {
                        throw("CellNetwork::CellNetworkIterator::operator++(): out of range. won't increment beyond end().");
                    }
                }

                CellNetworkIterator &
                operator++(int)
                {
                    CellNetworkIterator tmp(*this);
                    this->operator++();
                    return tmp;
                }

                /* NOTE: we need to check whether there's a matching type to the "left" of the current element that is
                 * reachable with operator--(). this code will return wrapped version of internal begin(), which is in
                 * general of false type and will produce segfaults on dereferencing..
                 *
                 * however, that is the case for all containers? the caller needs to check against begin() when using
                 * operator--, maybe like this:
                 *
                 *  it = end(); while (it != begin) { --it; foo(it); bar(it); }
                 *
                 * care must be taken not to abort before begin() has been seen and not to use anything beyond begin(),
                 * i.e. --begin(), which is undefined..
                 *
                 * so this might be a general problem and not specific to this special iterator here..*/
                CellNetworkIterator &
                operator--()
                {
                    DelegatorType *dp;

                    /* check un-decremented version of iterator against begin(), decrement _before_ the loop body to be
                     * able to deal with begin() as well (necessary since there are no relational operators such as <=
                     * in this case) */
                    while (this->int_it != this->int_ds->begin()) {
                        /* decrement iterator before loop body */
                        --this->int_it;

                        /* down-cast to delegator type, throw() if downcast fails. */
                        dp = dynamic_cast<DelegatorType *>(BaseType::getPtr(this->int_it));
                        if (dp) {
                            /* if type matches, return, otherwise decrement the internal iterator
                             * and restart the loop if end() hasn't been reached yet */
                            if ( this->pred->operator()(*dp) ) {
                                return (*this);
                            }
                        }
                        else {
                            throw("CellNetwork::CellNetworkIterator::operator-- failed to down-cast base type to delegator type. internal logic error.");
                        }
                    }
                    /* int_it is not this->int_ds->begin() and the type does not match. throw out of range exception */
                    throw("CellNetwork::CellNetworkIterator::operator--(): out of range, won't decrement further than begin().");
                }

                CellNetworkIterator
                operator--(int)
                {
                    CellNetworkIterator tmp(*this);
                    this->operator--();
                    return tmp;
                }

                bool
                operator==(const CellNetworkIterator &x) const
                {
                    if (this->explicitlyInvalid() || x.explicitlyInvalid()) {
                        return false;
                    }
                    else {
                        return (this->int_it == x.int_it);
                    }
                }

                bool
                operator!=(const CellNetworkIterator &x) const
                {
                    return !(this->operator==(x));
                }

                ValueType &
                operator*() const
                {
                    return *(this->operator->());
                }

                ValueType *
                operator->() const
                {
                    /* directly down-cast base type to value type */
                    ValueType   *vp = dynamic_cast<ValueType *>(BaseType::getPtr(this->int_it)); 

                    /* if cast has been successful, return ValueType pointer, otherwise throw(). */
                    if (vp) {
                        return vp;
                    }
                    else {
                        throw("CellNetworkIterator::operator->(): down-casting to value type failed. internal logic error.");
                    }
                }

                CellNetwork *
                getContainer()
                {
                    return (this->network);
                }
        };

        /* vertex iterator types for different levels of abstraction, derived from CellNetworkIterator */

        /* NOTE: forward declarations of iterators have been made above. */

        /* neuron vertex iterator, const and non-const */
        class neuron_iterator :
            public CellNetworkIterator<
                NeuronVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                NeuronVertexPred
            > 
        {
            public:
                neuron_iterator()
                    : CellNetworkIterator<
                          NeuronVertex,
                          typename Graph<Tn, Tv, Te>::Vertex,
                          NeuronVertex,
                          std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                          typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                          NeuronVertexPred
                      >()
                {
                }

                neuron_iterator(
                    CellNetwork                                                                            *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator    it)
                        : CellNetworkIterator<
                              NeuronVertex,
                              typename Graph<Tn, Tv, Te>::Vertex,
                              NeuronVertex,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                              NeuronVertexPred
                          >(C, &(C->V), it, &(C->isNeuronVertex) )
                {
                }

                neuron_iterator(neuron_iterator const &x)
                        : CellNetworkIterator<
                              NeuronVertex,
                              typename Graph<Tn, Tv, Te>::Vertex,
                              NeuronVertex,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                              NeuronVertexPred
                          >(x)
                {
                }

               ~neuron_iterator()
               {
               }

                /* copy constructors for all iterators referring to derived Vertex types inheriting from NeuronVertex, 
                 * since these iterators should be up-convertible just as the underlying types */
                neuron_iterator(soma_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_iterator(neurite_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_iterator(axon_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_iterator(dendrite_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        class neuron_const_iterator :
            public CellNetworkIterator<
                const NeuronVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                NeuronVertexPred
            > 
        {
            public:
                neuron_const_iterator()
                {
                }

                neuron_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator  it)
                        : CellNetworkIterator<
                              const NeuronVertex,
                              typename Graph<Tn, Tv, Te>::Vertex,
                              NeuronVertex,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                              NeuronVertexPred
                          >(C, &(C->V), it, &(C->isNeuronVertex) )
                {
                }

                neuron_const_iterator(neuron_const_iterator const &x)
                        : CellNetworkIterator<
                              const NeuronVertex,
                              typename Graph<Tn, Tv, Te>::Vertex,
                              NeuronVertex,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                              NeuronVertexPred
                          >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                neuron_const_iterator(neuron_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~neuron_const_iterator()
                {
                }

                /* copy constructors for all iterators referring to derived Vertex types inheriting from NeuronVertex, 
                 * since these iterators should be up-convertible just as the underlying types */
                neuron_const_iterator(soma_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_const_iterator(neurite_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_const_iterator(axon_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_const_iterator(dendrite_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        /* soma vertex iterator, const and non-const */
        class soma_iterator : 
            public CellNetworkIterator<
                SomaVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                SomaVertexPred
            > 
        {
            public:
                soma_iterator()
                    : CellNetworkIterator<
                            SomaVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            SomaVertexPred
                    >()
                {
                }

                soma_iterator(
                    CellNetwork                                                                            *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator    it)
                        : CellNetworkIterator<
                            SomaVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            SomaVertexPred
                        >(C, &(C->V), it, &(C->isSomaVertex) )
                {
                }

                soma_iterator(soma_iterator const &x)
                        : CellNetworkIterator<
                            SomaVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            SomaVertexPred
                        >(x)
                {
                }

               ~soma_iterator()
               {
               }
        };

        class soma_const_iterator :
            public CellNetworkIterator<
                const SomaVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                SomaVertexPred
            > 
        {
            public:
                soma_const_iterator()
                {
                }

                soma_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator  it)
                        : CellNetworkIterator<
                                const SomaVertex,
                                typename Graph<Tn, Tv, Te>::Vertex,
                                NeuronVertex,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                                SomaVertexPred
                        >(C, &(C->V), it, &(C->isSomaVertex) )
                {
                }

                soma_const_iterator(soma_const_iterator const &x)
                        : CellNetworkIterator<
                                const SomaVertex,
                                typename Graph<Tn, Tv, Te>::Vertex,
                                NeuronVertex,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                                SomaVertexPred
                        >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                soma_const_iterator(soma_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~soma_const_iterator()
                {
                }
        };

        /* neurite vertex iterator, const and non-const */
        class neurite_iterator : 
            public CellNetworkIterator<
                NeuriteVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                NeuriteVertexPred
            > 
        {
            public:
                neurite_iterator()
                    : CellNetworkIterator<
                            NeuriteVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            NeuriteVertexPred
                      >()
                {
                }

                neurite_iterator(
                    CellNetwork                                                                            *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator    it)
                        : CellNetworkIterator<
                            NeuriteVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            NeuriteVertexPred
                        >(C, &(C->V), it, &(C->isNeuriteVertex) )
                {
                }

                neurite_iterator(neurite_iterator const &x)
                        : CellNetworkIterator<
                            NeuriteVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            NeuriteVertexPred
                        >(x)
                {
                }

               ~neurite_iterator()
               {
               }

                /* copy constructors for axon_iterator and dendrite_iterator, which should be up-convertible as the
                 * underlying types are */
                neurite_iterator(axon_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neurite_iterator(dendrite_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        class neurite_const_iterator :
            public CellNetworkIterator<
                const NeuriteVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                NeuriteVertexPred
            > 
        {
            public:
                neurite_const_iterator()
                {
                }

                neurite_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator  it)
                        : CellNetworkIterator<
                                const NeuriteVertex,
                                typename Graph<Tn, Tv, Te>::Vertex,
                                NeuronVertex,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                                NeuriteVertexPred
                          >(C, &(C->V), it, &(C->isNeuriteVertex) )
                {
                }

                neurite_const_iterator(neurite_const_iterator const &x)
                        : CellNetworkIterator<
                                const NeuriteVertex,
                                typename Graph<Tn, Tv, Te>::Vertex,
                                NeuronVertex,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                                NeuriteVertexPred
                          >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                neurite_const_iterator(neurite_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~neurite_const_iterator()
                {
                }

                /* copy constructors for axon_const_iterator and dendrite_const_iterator, which should be up-convertible as the
                 * underlying types are */
                neurite_const_iterator(axon_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neurite_const_iterator(dendrite_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        /* axon vertex iterator, const and non-const */
        class axon_iterator : 
            public CellNetworkIterator<
                AxonVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                AxonVertexPred
            > 
        {
            public:
                axon_iterator()
                    : CellNetworkIterator<
                            AxonVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            AxonVertexPred
                      >()
                {
                }

                axon_iterator(
                    CellNetwork                                                                            *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator    it)
                        : CellNetworkIterator<
                            AxonVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            AxonVertexPred
                        >(C, &(C->V), it, &(C->isAxonVertex) )
                {
                }

                axon_iterator(axon_iterator const &x)
                        : CellNetworkIterator<
                            AxonVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            AxonVertexPred
                        >(x)
                {
                }

               ~axon_iterator()
               {
               }
        };

        class axon_const_iterator :
            public CellNetworkIterator<
                const AxonVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                AxonVertexPred
            > 
        {
            public:
                axon_const_iterator()
                {
                }

                axon_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator  it)
                        : CellNetworkIterator<
                                const AxonVertex,
                                typename Graph<Tn, Tv, Te>::Vertex,
                                NeuronVertex,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                                AxonVertexPred
                          >(C, &(C->V), it, &(C->isAxonVertex) )
                {
                }

                axon_const_iterator(axon_const_iterator const &x)
                        : CellNetworkIterator<
                                const AxonVertex,
                                typename Graph<Tn, Tv, Te>::Vertex,
                                NeuronVertex,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                                AxonVertexPred
                          >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                axon_const_iterator(axon_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~axon_const_iterator()
                {
                }
        };

        /* dendrite vertex iterator, const and non-const */
        class dendrite_iterator : 
            public CellNetworkIterator<
                DendriteVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                DendriteVertexPred
            > 
        {
            public:
                dendrite_iterator()
                    : CellNetworkIterator<
                            DendriteVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            DendriteVertexPred
                      >()
                {
                }

                dendrite_iterator(
                    CellNetwork                                                                            *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator    it)
                        : CellNetworkIterator<
                            DendriteVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            DendriteVertexPred
                        >(C, &(C->V), it, &(C->isDendriteVertex) )
                {
                }

                dendrite_iterator(dendrite_iterator const &x)
                        : CellNetworkIterator<
                            DendriteVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::iterator,
                            DendriteVertexPred
                        >(x)
                {
                }

               ~dendrite_iterator()
               {
               }
        };

        class dendrite_const_iterator :
            public CellNetworkIterator<
                const DendriteVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                DendriteVertexPred
            > 
        {
            public:
                dendrite_const_iterator()
                {
                }

                dendrite_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator  it)
                        : CellNetworkIterator<
                            const DendriteVertex,
                            typename Graph<Tn, Tv, Te>::Vertex,
                            NeuronVertex,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                            DendriteVertexPred
                        >(C, &(C->V), it, &(C->isDendriteVertex) )
                {
                }

                dendrite_const_iterator(dendrite_const_iterator const &x)
                    : CellNetworkIterator<
                        const DendriteVertex,
                        typename Graph<Tn, Tv, Te>::Vertex,
                        NeuronVertex,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>::const_iterator,
                        DendriteVertexPred
                    >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                dendrite_const_iterator(dendrite_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~dendrite_const_iterator()
                {
                }
        };

        /* edge iterator types for different levels of abstraction, derived from CellNetworkIterator */

        /* NOTE: forward declarations of iterators have been made above*/

        /* neuron edge iterator, const and non-const */
        class neuron_edge_iterator :
            public CellNetworkIterator<
                NeuronEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                NeuronEdgePred
            > 
        {
            public:
                neuron_edge_iterator()
                    : CellNetworkIterator<
                        NeuronEdge,
                        typename Graph<Tn, Tv, Te>::Edge,
                        NeuronEdge,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                        NeuronEdgePred
                      >()
                {
                }

                neuron_edge_iterator(
                    CellNetwork                                                                        *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator  it)
                        :   CellNetworkIterator<
                                NeuronEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                NeuronEdgePred
                            >(C, &(C->E), it, &(C->isNeuronEdge) )
                {
                }

                neuron_edge_iterator(neuron_edge_iterator const &x)
                        :   CellNetworkIterator<
                                NeuronEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                NeuronEdgePred
                            >(x)
                {
                }

               ~neuron_edge_iterator()
               {
               }

                /* copy constructors for all iterators referring to derived Vertex types inheriting from NeuronVertex, 
                 * since these iterators should be up-convertible just as the underlying types */
                neuron_edge_iterator(neurite_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_edge_iterator(axon_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_edge_iterator(dendrite_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_edge_iterator(axon_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_edge_iterator(dendrite_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                /* extra copy constructor for neurite segment iterators */
                neuron_edge_iterator(neurite_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        class neuron_edge_const_iterator :
            public CellNetworkIterator<
                const NeuronEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                NeuronEdgePred
            > 
        {
            public:
                neuron_edge_const_iterator()
                {
                }

                neuron_edge_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator    it)
                        : CellNetworkIterator<
                            const NeuronEdge,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            NeuronEdgePred
                        >(C, &(C->E), it, &(C->isNeuronEdge) )
                {
                }

                neuron_edge_const_iterator(neuron_edge_const_iterator const &x)
                        : CellNetworkIterator<
                            const NeuronEdge,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            NeuronEdgePred
                        >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                neuron_edge_const_iterator(neuron_edge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~neuron_edge_const_iterator()
                {
                }

                /* copy constructors for all iterators referring to derived Vertex types inheriting from NeuronVertex, 
                 * since these iterators should be up-convertible just as the underlying types */
                neuron_edge_const_iterator(axon_rootedge_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_edge_const_iterator(dendrite_rootedge_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_edge_const_iterator(axon_segment_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neuron_edge_const_iterator(dendrite_segment_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                /* extra copy constructor for neurite segment iterators */
                neuron_edge_const_iterator(neurite_segment_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        /* neurite root edge iterator, const and non-const */
        class neurite_rootedge_iterator :
            public CellNetworkIterator<
                NeuriteRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                NeuriteRootEdgePred
            > 
        {
            public:
                neurite_rootedge_iterator()
                    : CellNetworkIterator<
                        NeuriteRootEdge,
                        typename Graph<Tn, Tv, Te>::Edge,
                        NeuronEdge,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                        NeuriteRootEdgePred
                      >()
                {
                }

                neurite_rootedge_iterator(
                    CellNetwork                                                                        *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator  it)
                        :   CellNetworkIterator<
                                NeuriteRootEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                NeuriteRootEdgePred
                            >(C, &(C->E), it, &(C->isNeuriteRootEdge) )
                {
                }

                neurite_rootedge_iterator(neurite_rootedge_iterator const &x)
                        :   CellNetworkIterator<
                                NeuriteRootEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                NeuriteRootEdgePred
                            >(x)
                {
                }

               ~neurite_rootedge_iterator()
                {
                }

                /* copy constructors for all iterators referring to derived edge types inheriting from NeuriteSegment, 
                 * since these iterators should be up-convertible just as the underlying types */
                neurite_rootedge_iterator(axon_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neurite_rootedge_iterator(dendrite_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        class neurite_rootedge_const_iterator :
            public CellNetworkIterator<
                const NeuriteRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                NeuriteRootEdgePred
            > 
        {
            public:
                neurite_rootedge_const_iterator()
                {
                }

                neurite_rootedge_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator    it)
                        : CellNetworkIterator<
                            const NeuriteRootEdge,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            NeuriteRootEdgePred
                        >(C, &(C->E), it, &(C->isNeuriteRootEdge) )
                {
                }

                neurite_rootedge_const_iterator(neurite_rootedge_const_iterator const &x)
                        : CellNetworkIterator<
                            const NeuriteRootEdge,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            NeuriteRootEdgePred
                        >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                neurite_rootedge_const_iterator(neurite_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~neurite_rootedge_const_iterator()
                {
                }

                /* copy constructors for all iterators referring to derived edge types inheriting from NeuriteSegment, 
                 * since these iterators should be up-convertible just as the underlying types */
                neurite_rootedge_const_iterator(axon_rootedge_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neurite_rootedge_const_iterator(dendrite_rootedge_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }
        };

        /* axon root edge iterator, const and non-const */
        class axon_rootedge_iterator :
            public CellNetworkIterator<
                AxonRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                AxonRootEdgePred
            > 
        {
            public:
                axon_rootedge_iterator()
                    : CellNetworkIterator<
                        AxonRootEdge,
                        typename Graph<Tn, Tv, Te>::Edge,
                        NeuronEdge,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                        AxonRootEdgePred
                      >()
                {
                }

                axon_rootedge_iterator(
                    CellNetwork                                                                        *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator  it)
                        :   CellNetworkIterator<
                                AxonRootEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                AxonRootEdgePred
                            >(C, &(C->E), it, &(C->isAxonRootEdge) )
                {
                }

                axon_rootedge_iterator(axon_rootedge_iterator const &x)
                        :   CellNetworkIterator<
                                AxonRootEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                AxonRootEdgePred
                            >(x)
                {
                }

               ~axon_rootedge_iterator()
               {
               }
        };

        class axon_rootedge_const_iterator :
            public CellNetworkIterator<
                const AxonRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                AxonRootEdgePred
            > 
        {
            public:
                axon_rootedge_const_iterator()
                {
                }

                axon_rootedge_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator    it)
                        : CellNetworkIterator<
                              const AxonRootEdge,
                              typename Graph<Tn, Tv, Te>::Edge,
                              NeuronEdge,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                              AxonRootEdgePred
                          >(C, &(C->E), it, &(C->isAxonRootEdge) )
                {
                }

                axon_rootedge_const_iterator(axon_rootedge_const_iterator const &x)
                        : CellNetworkIterator<
                              const AxonRootEdge,
                              typename Graph<Tn, Tv, Te>::Edge,
                              NeuronEdge,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                              AxonRootEdgePred
                          >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                axon_rootedge_const_iterator(axon_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~axon_rootedge_const_iterator()
                {
                }
        };

        /* axon root edge iterator, const and non-const */
        class dendrite_rootedge_iterator :
            public CellNetworkIterator<
                DendriteRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                DendriteRootEdgePred
            > 
        {
            public:
                dendrite_rootedge_iterator()
                    : CellNetworkIterator<
                        DendriteRootEdge,
                        typename Graph<Tn, Tv, Te>::Edge,
                        NeuronEdge,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                        DendriteRootEdgePred
                      >()
                {
                }

                dendrite_rootedge_iterator(
                    CellNetwork                                                                        *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator  it)
                        :   CellNetworkIterator<
                                DendriteRootEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                DendriteRootEdgePred
                            >(C, &(C->E), it, &(C->isDendriteRootEdge) )
                {
                }

                dendrite_rootedge_iterator(dendrite_rootedge_iterator const &x)
                        :   CellNetworkIterator<
                                DendriteRootEdge,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                DendriteRootEdgePred
                            >(x)
                {
                }

               ~dendrite_rootedge_iterator()
               {
               }
        };

        class dendrite_rootedge_const_iterator :
            public CellNetworkIterator<
                const DendriteRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                DendriteRootEdgePred
            > 
        {
            public:
                dendrite_rootedge_const_iterator()
                {
                }

                dendrite_rootedge_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator    it)
                        : CellNetworkIterator<
                              const DendriteRootEdge,
                              typename Graph<Tn, Tv, Te>::Edge,
                              NeuronEdge,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                              DendriteRootEdgePred
                          >(C, &(C->E), it, &(C->isDendriteRootEdge) )
                {
                }

                dendrite_rootedge_const_iterator(dendrite_rootedge_const_iterator const &x)
                        : CellNetworkIterator<
                              const DendriteRootEdge,
                              typename Graph<Tn, Tv, Te>::Edge,
                              NeuronEdge,
                              std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                              typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                              DendriteRootEdgePred
                          >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                dendrite_rootedge_const_iterator(dendrite_rootedge_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~dendrite_rootedge_const_iterator()
                {
                }
        };

        /* neurite segment iterator, const and non-const */
        class neurite_segment_iterator :
            public CellNetworkIterator<
                NeuriteSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                NeuriteSegmentPred
            > 
        {
            public:
                neurite_segment_iterator()
                    : CellNetworkIterator<
                        NeuriteSegment,
                        typename Graph<Tn, Tv, Te>::Edge,
                        NeuronEdge,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                        NeuriteSegmentPred
                      >()
                {
                }

                neurite_segment_iterator(
                    CellNetwork                                                                        *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator  it)
                        :   CellNetworkIterator<
                                NeuriteSegment,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                NeuriteSegmentPred
                            >(C, &(C->E), it, &(C->isNeuriteSegment) )
                {
                }

                neurite_segment_iterator(neurite_segment_iterator const &x)
                        :   CellNetworkIterator<
                                NeuriteSegment,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                NeuriteSegmentPred
                            >(x)
                {
                }

               ~neurite_segment_iterator()
               {
               }

                /* copy constructors for all iterators referring to derived edge types inheriting from NeuriteSegment, 
                 * since these iterators should be up-convertible just as the underlying types */
                neurite_segment_iterator(axon_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neurite_segment_iterator(dendrite_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                /* cast operator that allows implicit conversion to Graph<Tn, Tv, Te>::edge_iterator */
                operator typename Graph<Tn, Tv, Te>::edge_iterator()
                {
                    return typename Graph<Tn, Tv, Te>::edge_iterator(
                        this->network,
                        this->int_it
                    );
                }

                operator typename Graph<Tn, Tv, Te>::edge_const_iterator()
                {
                    return typename Graph<Tn, Tv, Te>::edge_const_iterator(
                        this->network,
                        this->int_it
                    );
                }
        };

        class neurite_segment_const_iterator :
            public CellNetworkIterator<
                const NeuriteSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                NeuriteSegmentPred
            > 
        {
            public:
                neurite_segment_const_iterator()
                {
                }

                neurite_segment_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator    it)
                        : CellNetworkIterator<
                            const NeuriteSegment,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            NeuriteSegmentPred
                        >(C, &(C->E), it, &(C->isNeuriteSegment) )
                {
                }

                neurite_segment_const_iterator(neurite_segment_const_iterator const &x)
                        : CellNetworkIterator<
                            const NeuriteSegment,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            NeuriteSegmentPred
                        >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                neurite_segment_const_iterator(neurite_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~neurite_segment_const_iterator()
                {
                }

                /* copy constructors for all iterators referring to derived edge types inheriting from NeuriteSegment, 
                 * since these iterators should be up-convertible just as the underlying types */
                neurite_segment_const_iterator(axon_segment_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                neurite_segment_const_iterator(dendrite_segment_const_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_it    = x.int_it;
                    this->int_ds    = x.int_ds;
                }

                /* cast operator that allows implicit conversion to Graph<Tn, Tv, Te>::edge_const_iterator */
                operator typename Graph<Tn, Tv, Te>::edge_const_iterator()
                {
                    return typename Graph<Tn, Tv, Te>::edge_const_iterator(
                        this->network,
                        this->int_it
                    );
                }
        };

        /* axon segment iterator, const and non-const */
        class axon_segment_iterator :
            public CellNetworkIterator<
                AxonSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                AxonSegmentPred
            > 
        {
            public:
                axon_segment_iterator()
                    : CellNetworkIterator<
                        AxonSegment,
                        typename Graph<Tn, Tv, Te>::Edge,
                        NeuronEdge,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                        AxonSegmentPred
                      >()
                {
                }

                axon_segment_iterator(
                    CellNetwork                                                                        *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator  it)
                        :   CellNetworkIterator<
                                AxonSegment,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                AxonSegmentPred
                            >(C, &(C->E), it, &(C->isAxonSegment) )
                {
                }

                axon_segment_iterator(axon_segment_iterator const &x)
                        :   CellNetworkIterator<
                                AxonSegment,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                AxonSegmentPred
                            >(x)
                {
                }

               ~axon_segment_iterator()
               {
               }
        };

        class axon_segment_const_iterator :
            public CellNetworkIterator<
                const AxonSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                AxonSegmentPred
            > 
        {
            public:
                axon_segment_const_iterator()
                {
                }

                axon_segment_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator    it)
                        : CellNetworkIterator<
                            const AxonSegment,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            AxonSegmentPred
                        >(C, &(C->E), it, &(C->isAxonSegment) )
                {
                }

                axon_segment_const_iterator(axon_segment_const_iterator const &x)
                        : CellNetworkIterator<
                            const AxonSegment,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            AxonSegmentPred
                        >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                axon_segment_const_iterator(axon_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~axon_segment_const_iterator()
                {
                }
        };

        /* dendrite segment iterator, const and non-const */
        class dendrite_segment_iterator :
            public CellNetworkIterator<
                DendriteSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                DendriteSegmentPred
            > 
        {
            public:
                dendrite_segment_iterator()
                    : CellNetworkIterator<
                        DendriteSegment,
                        typename Graph<Tn, Tv, Te>::Edge,
                        NeuronEdge,
                        std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                        typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                        DendriteSegmentPred
                      >()
                {
                }

                dendrite_segment_iterator(
                    CellNetwork                                                                        *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator  it)
                        :   CellNetworkIterator<
                                DendriteSegment,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                DendriteSegmentPred
                            >(C, &(C->E), it, &(C->isDendriteSegment) )
                {
                }

                dendrite_segment_iterator(dendrite_segment_iterator const &x)
                        :   CellNetworkIterator<
                                DendriteSegment,
                                typename Graph<Tn, Tv, Te>::Edge,
                                NeuronEdge,
                                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::iterator,
                                DendriteSegmentPred
                            >(x)
                {
                }

               ~dendrite_segment_iterator()
               {
               }
        };

        class dendrite_segment_const_iterator :
            public CellNetworkIterator<
                const DendriteSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                DendriteSegmentPred
            > 
        {
            public:
                dendrite_segment_const_iterator()
                {
                }

                dendrite_segment_const_iterator(
                    CellNetwork                                                                                *C,
                    typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator    it)
                        : CellNetworkIterator<
                            const DendriteSegment,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            DendriteSegmentPred
                        >(C, &(C->E), it, &(C->isDendriteSegment) )
                {
                }

                dendrite_segment_const_iterator(dendrite_segment_const_iterator const &x)
                        : CellNetworkIterator<
                            const DendriteSegment,
                            typename Graph<Tn, Tv, Te>::Edge,
                            NeuronEdge,
                            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                            typename std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>::const_iterator,
                            DendriteSegmentPred
                        >(x)
                {
                }

                /* provide copy constructor with non-const argument to allow implicit conversion of non-const to const
                 * iterator */
                dendrite_segment_const_iterator(dendrite_segment_iterator const &x)
                {
                    this->graph     = x.graph;
                    this->network   = x.network;
                    this->int_ds    = x.int_ds;
                    this->int_it    = x.int_it;
                }

               ~dendrite_segment_const_iterator()
                {
                }
        };
        /*
        typedef CellNetworkIterator<
            NeuronEdge,
            typename Graph<Tn, Tv, Te>::Edge,
            NeuronEdge,
            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
            NeuronEdgePred
        > neuron_edge_iteartor;

        typedef CellNetworkIterator<
            NeuriteSegment,
            typename Graph<Tn, Tv, Te>::Edge,
            NeuronEdge,
            std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
            NeuriteSegmentPred
        > neurite_segment_iterator;

        */
          
        /* generic cell network accessor template. derived specializations for vertices (neuron /
         * neurite / soma / axon / ..) and edeges (neurite segment) of increasing depth in the
         * inheritance trees are declared below. uses CellNetworkIterator from above.
         *
         * it relies and depends on the internally used data structure for vertices / edges
         *
         *      std::map<uint32_t, {Vertex, Edge, Soma, NeuriteSegment, ...}PointerType>
         *
         * it contains the code shared between all accessors for specific vertex / edge types. since
         * explicit down-casting should be avoided as far as possible, an individual accessor will
         * be defined for every specialied vertex / edge type, e.g. for somata / axon vertices /
         * dendrite vertices / neurite segment edges / soma-neurite edges / ..
         *
         * the class CellNetworkAccessor provides a generic interface and is instantiated with
         * specializations of the above "skipping" iterator CellNetworkIterator as iterator types.
         * these iterator types use the generic this->vertices map containing pointers to vertices
         * internally, but enables access only to objects of matching specialized type (e.g.
         * somata), which are returned in a properly down-cast fashion.
         *
         * the overall structure is very similar to {Vertex,Edge}Accessor in Graph (or
         * {Vertex,Face}Accessor in MEsh)): begin(), end(), exists(), at(), find(), etc.. */
        template<
            typename ValueType,
            typename BaseType,
            typename DelegatorType,
            typename InternalType,
            typename ValuePred,
            typename IteratorType,
            typename ConstIteratorType,
            typename DataType
        >
        class CellNetworkAccessor {
            friend class CellNetwork;

            protected:
                CellNetwork                            &C;
                InternalType                           &int_ds;
                ValuePred const                        &pred;

                /* protected constructor */
                                                        CellNetworkAccessor(
                                                            CellNetwork     &C,
                                                            InternalType                   &X,
                                                            ValuePred const                &p);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                CellNetworkAccessor                    &operator=(CellNetworkAccessor const &x)             = delete;
                                                        CellNetworkAccessor(CellNetworkAccessor const &x)   = delete;

            public:
                IteratorType                            begin();
                ConstIteratorType                       begin() const;

                IteratorType                            end();
                ConstIteratorType                       end() const;

                IteratorType                            find(uint32_t id);
                ConstIteratorType                       find(uint32_t id) const;

                ValueType                              &at(uint32_t id);
                ValueType const                        &at(uint32_t id) const;

                bool                                    exists(uint32_t id) const;

                /* NOTE: since "skipping" iterators are used internally to scan through a map of base class types, it is
                 * not obvious how many objects of matching types there are => size() becomes an O(n) operation, just as
                 * in some std::list implementations. shouldn't be too problematic. the only alternative would be to
                 * init to zero and update on new insertions, but the accessor class cannot guarantee that these updates
                 * are performed faithfully. so for now: the O(n) version. */
                size_t                                  size() const;

                /*
                virtual std::pair<IteratorType, bool>   insert(DataType const &data) = 0;
                virtual IteratorType                    erase(IteratorType it) = 0;
                virtual bool                            erase(uint32_t id);
                */
        };


        /* derived CellNetworkAccessor specialization for abstract type NeuronVertex, using neuron_(const_)iterator.
         * does NOT provide insert functionality, since abstract and semantically incomplete type NeuronVertex cannot be
         * created without knowledge which concrete type should be created.
         *
         * however, erase functionality is only provided at this most abstract level, since semantics of edge erasing
         * does not depend on the particular type of edge being erased. */
        class NeuronVertexAccessor :
            public CellNetworkAccessor<
                NeuronVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                NeuronVertexPred,
                neuron_iterator,
                neuron_const_iterator,
                Tv>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            NeuronVertexAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            NeuronVertexAccessor(NeuronVertexAccessor const &x) = delete;
                NeuronVertexAccessor                       &operator=(NeuronVertexAccessor const &x) = delete;

            public:
                neuron_iterator                             erase(neuron_iterator it);
                bool                                        erase(uint32_t id);
        };

        /* derived CellNetworkAccessor specialization for abstract type NeuriteVertex, using neurite_(const_)iterator
         * this does NOT provide insert functionality, since abstract and semantically incomplete type NeuronVertex
         * cannot be created without knowledge which concrete type should be created. */
        class NeuriteVertexAccessor :
            public CellNetworkAccessor<
                NeuriteVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                NeuriteVertexPred,
                neurite_iterator,
                neurite_const_iterator,
                Tv>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            NeuriteVertexAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            NeuriteVertexAccessor(NeuriteVertexAccessor const &x) = delete;
                NeuriteVertexAccessor                       &operator=(NeuriteVertexAccessor const &x) = delete;

            public:
                /*
                neurite_iterator                            erase(neurite_iterator it);
                bool                                        erase(uint32_t id);
                */
        };

        /* derived CellNetworkAccessor specialization for somas, using soma_(const_)iterator.
         * provides insert functionality. */
        class SomaAccessor :
            public CellNetworkAccessor<
                SomaVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                SomaVertexPred,
                soma_iterator,
                soma_const_iterator,
                Tv>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                        SomaAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                        SomaAccessor(SomaAccessor const &x) = delete;
                SomaAccessor                           &operator=(SomaAccessor const &x) = delete;

            public:
                soma_iterator                           insert(
                                                            Graph<uint32_t, CellSection, Common::UnitType> const   &section_graph,
                                                            Tv const                                               &vertex_data = Tv(),
                                                            Tso const                                              &soma_data = Tso());

                /*
                soma_iterator                           erase(soma_iterator it);
                bool                                    erase(uint32_t id);
                */
        };

        /* derived CellNetworkAccessor specialization for axon vertices, using axon_(const_)iterator
         * provides insert functionality. */
        class AxonAccessor :
            public CellNetworkAccessor<
                AxonVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                AxonVertexPred,
                axon_iterator,
                axon_const_iterator,
                Tv>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                        AxonAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                        AxonAccessor(AxonAccessor const &x) = delete;
                AxonAccessor                           &operator=(AxonAccessor const &x) = delete;

            public:
                axon_iterator                           insert(
                                                            CellSection const  &sections,
                                                            Tax const          &axon_data = Tax(),
                                                            Tnv const          &neurite_vertex_data = Tnv(),
                                                            Tv const           &vertex_data = Tv());

                /*
                axon_iterator                           erase(axon_iterator it);
                bool                                    erase(uint32_t id);
                */
        };

        /* derived CellNetworkAccessor specialization for dendrite vertices, using dendrite_(const_)iterator
         * provides insert functionality. */
        class DendriteAccessor :
            public CellNetworkAccessor<
                DendriteVertex,
                typename Graph<Tn, Tv, Te>::Vertex,
                NeuronVertex,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::VertexPointerType>,
                DendriteVertexPred,
                dendrite_iterator,
                dendrite_const_iterator,
                Tv>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            DendriteAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            DendriteAccessor(DendriteAccessor const &x) = delete;
                DendriteAccessor                           &operator=(DendriteAccessor const &x) = delete;

            public:
                dendrite_iterator                           insert(
                                                                CellSection const  &sections,
                                                                bool                apical_dendrite,
                                                                Tde const          &dendrite_data = Tde(),
                                                                Tnv const          &neurite_vertex_data = Tnv(),
                                                                Tv const           &vertex_data = Tv());

                /*
                dendrite_iterator                           erase(dendrite_iterator it);
                bool                                        erase(uint32_t id);
                */
        };

        /* derived CellNetworkAccessor specialization for abstract type NeuronEdge, using neuron_edge_(const_)iterator.
         * does NOT provde insert functionality, since abstract and semantically incomplete type NeuronVertex cannot be
         * created without knowledge which concrete type should be created.
         *
         * however, erase functionality is only provided at this most abstract level, since semantics of edge erasing
         * does not depend on the particular type of edge being erased. */
        class NeuronEdgeAccessor :
            public CellNetworkAccessor<
                NeuronEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                NeuronEdgePred,
                neuron_edge_iterator,
                neuron_edge_const_iterator,
                Te>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            NeuronEdgeAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            NeuronEdgeAccessor(NeuronEdgeAccessor const &x) = delete;
                NeuronEdgeAccessor                         &operator=(NeuronEdgeAccessor const &x) = delete;

            public:
                neuron_edge_iterator                        erase(neuron_edge_iterator it);
                bool                                        erase(uint32_t id);
        };

        /* derived CellNetworkAccessor specialization for type NeuriteRootEdge, using neurite_rootedge_(const_)iterator.
         * does NOT provde insert functionality, since abstract and semantically incomplete type NeuronVertex cannot be
         * created without knowledge which concrete type should be created.  erase() functionality is provided in
         * NeuronEdgeAccessor only. */
        class NeuriteRootEdgeAccessor : 
            public CellNetworkAccessor<
                NeuriteRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                NeuriteRootEdgePred,
                neurite_rootedge_iterator,
                neurite_rootedge_const_iterator,
                Te>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            NeuriteRootEdgeAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            NeuriteRootEdgeAccessor(NeuriteRootEdgeAccessor const &x) = delete;
                NeuriteRootEdgeAccessor                    &operator=(NeuriteRootEdgeAccessor const &x) = delete;

            public:

                std::pair<
                        neurite_rootedge_iterator,
                        std::list<neurite_segment_iterator>
                    >                                       split(
                                                                neurite_rootedge_iterator   nr_it,
                                                                std::vector<
                                                                        std::pair<
                                                                            Vec3<R>,
                                                                            R
                                                                        >
                                                                    > const                &intermediate_vertex_data);

        };


        /* derived CellNetworkAccessor specialization for type AxonRootEdge, using axon_rootedge_(const_)iterator.
         * provides neither insert nor erase functionality */
        class AxonRootEdgeAccessor :
            public CellNetworkAccessor<
                AxonRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                AxonRootEdgePred,
                axon_rootedge_iterator,
                axon_rootedge_const_iterator,
                Te>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            AxonRootEdgeAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            AxonRootEdgeAccessor(AxonRootEdgeAccessor const &x) = delete;
                AxonRootEdgeAccessor                       &operator=(AxonRootEdgeAccessor const &x) = delete;

            public:
                std::pair<axon_rootedge_iterator, bool>     insert(
                                                                soma_iterator  &v_src_it,
                                                                axon_iterator  &v_dst_it,
                                                                Tar const      &axon_root_edge_data = Tar(),
                                                                Tnr const      &neurite_root_edge_data = Tnr(), 
                                                                Te const       &edge_data = Te());

        };

        /* derived CellNetworkAccessor specialization for type DendriteRootEdge, using
         * dendrite_rootedge_(const_)iterator.  provides insert functionality. erase() functionality is provided in
         * NeuronEdgeAccessor only. */
        class DendriteRootEdgeAccessor :
            public CellNetworkAccessor<
                DendriteRootEdge,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                DendriteRootEdgePred,
                dendrite_rootedge_iterator,
                dendrite_rootedge_const_iterator,
                Te>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            DendriteRootEdgeAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            DendriteRootEdgeAccessor(DendriteRootEdgeAccessor const &x) = delete;
                DendriteRootEdgeAccessor                   &operator=(DendriteRootEdgeAccessor const &x) = delete;

            public:
                std::pair<dendrite_rootedge_iterator, bool> insert(
                                                                soma_iterator       v_src_it,
                                                                dendrite_iterator   v_dst_it,
                                                                Tdr const          &dendrite_root_edge_data = Tdr(),
                                                                Tnr const          &neurite_root_edge_data = Tnr(), 
                                                                Te const           &edge_data = Te());
        };

        /* derived CellNetworkAccessor specialization for abstract type NeuriteSegment, using
         * neurite_segment_(const_)iterator.  does provde insert functionality by checking to given neurite_iterators
         * for types and forwarding to DendriteSegmentAccessor/AxonSegmentAccessor::insert. since no type information is
         * available, default-constructed instances of Tas and Tds are used as dendrite/axon specific annotated data.
         * erase() functionality is provided in NeuronEdgeAccessor only. */
        class NeuriteSegmentAccessor :
            public CellNetworkAccessor<
                NeuriteSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                NeuriteSegmentPred,
                neurite_segment_iterator,
                neurite_segment_const_iterator,
                Te>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            NeuriteSegmentAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            NeuriteSegmentAccessor(NeuronEdgeAccessor const &x) = delete;
                NeuriteSegmentAccessor                     &operator=(NeuriteSegmentAccessor const &x) = delete;

            public:
                /* insert neurite segment (u, v) for two neurite vertices (u, v), which are checked in the
                 * implementation to be of matching type */
                std::pair<neurite_segment_iterator, bool>   insert(
                                                                neurite_iterator    v_src_it,
                                                                neurite_iterator    v_dst_it,
                                                                Tns const          &neurite_segment_data = Tns(),
                                                                Te const           &edge_data = Te());

                std::list<neurite_segment_iterator>         split(
                                                                neurite_segment_iterator            ns_it,
                                                                std::vector<
                                                                        std::pair<
                                                                            Vec3<R>,
                                                                            R
                                                                        >
                                                                    > const                        &intermediate_vertex_data);

                neurite_iterator                            collapse(
                                                                neurite_segment_iterator            ns_it,
                                                                std::pair<Vec3<R>, R> const        &collapsed_vertex_data);

        };

        /* derived CellNetworkAccessor specialization for type AxonSegment, using axon_segment_(const_)iterator.
         * provides insert functionality. erase() functionality is provided in NeuronEdgeAccessor only. */
        class AxonSegmentAccessor :
            public CellNetworkAccessor<
                AxonSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                AxonSegmentPred,
                axon_segment_iterator,
                axon_segment_const_iterator,
                Te>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            AxonSegmentAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            AxonSegmentAccessor(AxonSegmentAccessor const &x) = delete;
                AxonSegmentAccessor                        &operator=(AxonSegmentAccessor const &x) = delete;

            public:
                std::pair<axon_segment_iterator, bool>      insert(
                                                                axon_iterator       v_src_it,
                                                                axon_iterator       v_dst_it,
                                                                Tas const          &axon_segment_data = Tas(),
                                                                Tns const          &neurite_segment_data = Tns(),
                                                                Te const           &edge_data = Te());
        };

        /* derived CellNetworkAccessor specialization for type DendriteSegment, using dendrite_segment_(const_)iterator.
         * provides insert functionality. erase() functionality is provided in NeuronEdgeAccessor only. */
        class DendriteSegmentAccessor :
            public CellNetworkAccessor<
                DendriteSegment,
                typename Graph<Tn, Tv, Te>::Edge,
                NeuronEdge,
                std::map<uint32_t, typename Graph<Tn, Tv, Te>::EdgePointerType>,
                DendriteSegmentPred,
                dendrite_segment_iterator,
                dendrite_segment_const_iterator,
                Te>
        {
            friend class CellNetwork;

            private:
                /* private ctor */
                                                            DendriteSegmentAccessor(CellNetwork &C);

                /* accessors are neither copyable nor assignable and exist but once for every CellNetwork */
                                                            DendriteSegmentAccessor(DendriteSegmentAccessor const &x) = delete;
                DendriteSegmentAccessor                    &operator=(DendriteSegmentAccessor const &x) = delete;

            public:
                std::pair<dendrite_segment_iterator, bool>  insert(
                                                                dendrite_iterator   v_src_it,
                                                                dendrite_iterator   v_dst_it,
                                                                Tds const          &dendrite_segment_data = Tds(),
                                                                Tns const          &neurite_segment_data = Tns(), 
                                                                Te const           &edge_data = Te());
        };

    protected:
        /* struct used for SWC parsing */
        struct SWCNode {
            uint32_t                compartment_id;
            uint32_t                compartment_type;
            int32_t                 parent_id;
            std::list<CellSection>  sections;
            std::list<uint32_t>     child_ids;

            SWCNode(
                    uint32_t                        _compartment_id,
                    uint32_t                        _compartment_type,
                    int32_t                         _parent_id,
                    std::list<CellSection> const&   _sections)
            : compartment_id(_compartment_id), compartment_type(_compartment_type),
              parent_id(_parent_id), sections(_sections)
            {}
        };

        /* enum for compartment types */
        enum SWCCompartmentTypes {
            SOMA_COMPARTMENT            = 1,
            AXON_COMPARTMENT            = 2,
            BASAL_DENDRITE_COMPARTMENT  = 3,
            APICAL_DENDRITE_COMPARTMENT = 4
        };

    /* protected members */
    protected:
        std::string                                 network_name;
        bool                                        network_info_initialized;

        Vec3<R>                                     global_coordinate_displacement;
        R                                           global_coordinate_scaling_factor;

    /* "access declaration", achieved with using delcarations in C++11, to make vertex / edge accessor from Graph<Tn,
     * Tv, Te> private to the public: it must not be possible to insert generic Graph::Vertex into a cell network */
    /*
    private:
        using Graph<Tn, Tv, Te>::vertices;
        using Graph<Tn, Tv, Te>::edges;
    */

    /* public members */
    public:
        /* the various accessors */
        NeuronVertexAccessor                        neuron_vertices;
        SomaAccessor                                soma_vertices;
        NeuriteVertexAccessor                       neurite_vertices;
        AxonAccessor                                axon_vertices; 
        DendriteAccessor                            dendrite_vertices; 

        NeuronEdgeAccessor                          neuron_edges;
        NeuriteRootEdgeAccessor                     neurite_root_edges;
        AxonRootEdgeAccessor                        axon_root_edges;
        DendriteRootEdgeAccessor                    dendrite_root_edges;
        NeuriteSegmentAccessor                      neurite_segments;
        AxonSegmentAccessor                         axon_segments;
        DendriteSegmentAccessor                     dendrite_segments;

        /* static predicates (stored in std::function objects) useful for traversal calls */
        static std::function<bool(NeuronVertex const &)>    neuron_vertex_true;
        static std::function<bool(NeuronVertex const &)>    neuron_vertex_false;
        static std::function<bool(NeuronVertex const &)>    is_cell_vertex;
        static std::function<bool(NeuronVertex const &)>    is_neurite_vertex;
        static std::function<bool(NeuronVertex const &)>    is_axon_vertex;
        static std::function<bool(NeuronVertex const &)>    is_dendrite_vertex;

        static std::function<bool(NeuronEdge const &)>      neuron_edge_true;
        static std::function<bool(NeuronEdge const &)>      neuron_edge_false;
        static std::function<bool(NeuronEdge const &)>      is_neurite_segment;

        /* template methods returning predicates that match only for the input vertex / edge of the given template
         * type VertexType / EdgeType */
        template<typename VertexType>
        static std::function<bool(VertexType const &)>      genVertexSelectionPred(VertexType const &v);

        template<typename EdgeType>
        static std::function<bool(EdgeType const &)>        genEdgeSelectionPred(EdgeType const &e);

    /* protected methods */
    protected:
        void                                        checkTopology();
        void                                        initializeNetworkInfo();

        /* down-cast of Neuron{Veretx,Edge} to a specialized type */
        template<typename VertexType>
        VertexType                                 *downcastVertex(NeuronVertex *v);

        template<typename EdgeType>
        EdgeType                                   *downcastEdge(NeuronEdge *e);

        uint32_t                                    getVertexType(neuron_const_iterator const &v_it) const;
        uint32_t                                    getEdgeType(neuron_edge_const_iterator const &e_it) const;

    protected:
        /* struct storing statistical information about CellNetwork */
        struct CellNetworkStatistics {
            R       neurite_tree_total_length;
            R       morphological_diameter;

            R       neuron_vertices_dist_min;
            R       neuron_vertices_dist_max;
            R       neuron_vertices_dist_avg;
            R       neuron_vertices_dist_sigma;

            R       neuron_vertices_xmin;
            R       neuron_vertices_xmax;
            R       neuron_vertices_xavg;
            R       neuron_vertices_xsigma;

            R       neuron_vertices_ymin;
            R       neuron_vertices_ymax;
            R       neuron_vertices_yavg;
            R       neuron_vertices_ysigma;

            R       neuron_vertices_zmin;
            R       neuron_vertices_zmax;
            R       neuron_vertices_zavg;
            R       neuron_vertices_zsigma;

            R       neurite_radius_min;
            R       neurite_radius_max;
            R       neurite_radius_avg;
            R       neurite_radius_sigma;

            R       neurite_segment_length_min;
            R       neurite_segment_length_max;
            R       neurite_segment_length_avg;
            R       neurite_segment_length_sigma;

            R       neurite_segment_reduced_length_min;
            R       neurite_segment_reduced_length_max;
            R       neurite_segment_reduced_length_avg;
            R       neurite_segment_reduced_length_sigma;

            R       neurite_segment_angle_min;
            R       neurite_segment_angle_max;
            R       neurite_segment_angle_avg;
            R       neurite_segment_angle_sigma;

            R       neurite_segment_abs_radius_diff_max;
            R       neurite_segment_abs_radius_diff_avg;
            R       neurite_segment_abs_radius_diff_sigma;

            R       neurite_segment_radius_ratio_max;
            R       neurite_segment_radius_ratio_avg;
            R       neurite_segment_radius_ratio_sigma;
        };

    /* public method interface */
    public:
                                                    CellNetwork(const std::string& network_name = std::string(""));
        /* FIXME: implement copy ctor, assignment operator using "deep-copy" clone()-ing and topology update, similar to
         * Graph and Mesh. */
                                                    CellNetwork(CellNetwork const &X) = delete;
        CellNetwork                                &operator=(CellNetwork const &X) = delete;

        /* virtual destructor for polymorphism */
        virtual                                    ~CellNetwork();

        /* topology related methods */
        /* traverse CellNetwork breadth-first starting at vstart_it, but
         *
         *  1.  consider only edges for which edge_pred returns true AND for which both endpoints satisfy vetex
         *  predicate vertex_pred. fill return list reachable_edges only with reachable edges whose type matches the
         *  template type EdgeType.
         *
         *  2.  consider only vertices for which vertex_pred returns true and fill return list reachable_vertices only
         *  with reachable vertices whose type matches the template type VertexType.
         *
         *  3.  terminate traversal whenever the predicate vertex_termination_pred returns true for the currently
         *  visited vertex, which must be of matching type VertexType.
         *
         *  4.  terminate traversal whenever the predicate edge_termination_pred returns true for the currently
         *  inspected edge, which must be of matching type EdgeType.
         *
         *  5. the directed flag is used to control whether the direction of all (internally always directed) edges is
         *  taken into account or not.  => if directed == false, the neighbours of a vertex are all in- and
         *  out-neighbours and hence an undirected version of (this) CellNetwork (with possibly doubled undirected
         *  edges) will be traversed. NOTE: the return list reachable_edges still contains pointers to directed edges.
         *  it is the callers responsibility to interpret the information correctly if undirected == true: the direction
         *  is meaningless and the edge should be interpreted as a 2-subset stored in arbitrary order rather than an
         *  ordered pair.  for example, an undirected path connecting two neurite vertices from different neurites will
         *  contain "reversely" directed edges from the start neurite vertex to the soma and "correctly" directed edges
         *  from the soma to the destinatino neurite vertex. 
         *
         *  NOTE: the filter functionality can be "disabled" by choosing the defaults: base classes Neuron{Vertex,Edge}
         *  as types and the trivial predicates (which just return true). these defaults result in a normal traversal of
         *  the entire CellNetwork, where all vertex neighbours / edges being are considered. choosing vertex / edge
         *  types further down the hierarchy and choosing specialized predicates will result in increasing filtering and
         *  different behaviour. */
        template <typename VertexType = NeuronVertex, typename EdgeType = NeuronEdge>
        void                                        traverseBreadthFirst(
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
                                                        std::function<bool(NeuronVertex const &v)> const   &vertex_pred,// =           [] (NeuronVertex const &v) -> bool { return true; },
                                                        std::function<bool(NeuronEdge const &e)> const     &edge_pred,// =           [] (NeuronVertex const &v) -> bool { return true; },
                                                        std::function<bool(VertexType const &v)> const     &vertex_selection_pred,
                                                        std::function<bool(EdgeType const &e)> const       &edge_selection_pred,
                                                        std::function<bool(VertexType const &v)> const     &vertex_termination_pred,
                                                        std::function<bool(EdgeType const &e)> const       &edge_termination_pred,
                                                        int32_t                                             tid_arg = -1);// =    [] (VertexType const &v)   -> bool { return true; });
#if 0
                                                        old function pointer type signature: statless => cannot use capturing lambdas or functors of all sorts..

                                                        bool                                          (&vertex_pred)(NeuronVertex const &v),// =           [] (NeuronVertex const &v) -> bool { return true; },
                                                        bool                                          (&edge_pred)(NeuronEdge const &e),// =               [] (NeuronEdge const &e)   -> bool { return true; },
                                                        bool                                          (&vertex_selection_pred)(VertexType const &v),// =   [] (NeuronVertex const &v) -> bool { return true; },
                                                        bool                                          (&edge_selection_pred)(EdgeType const &e),// =       [] (NeuronEdge const &e)   -> bool { return true; },
                                                        bool                                          (&vertex_termination_pred)(VertexType const &v),// = [] (VertexType const &v)   -> bool { return true; },
                                                        bool                                          (&edge_termination_pred)(EdgeType const &e),
#endif

        /* get both integral depth (graph-theoretic depth) and chord-depth (chord-length of longest path to a neurite
         * leaf) of neurite sub-tree rooted in neurite vertex n_it */
        std::pair<uint32_t, R>                      getNeuriteSubTreeDepths(neurite_iterator n_it);

        /* get directed path and path lengths from vertex u to v */
        void                                        getDirectedPath(
                                                        const neuron_iterator&      u_it,
                                                        const neuron_iterator&      v_it,
                                                        std::list<NeuronVertex *>  &path,
                                                        R                          &pathlen);

        /* get all neurite vertices / segments reachable from v in an undirected traversal when only traversing neurite
         * segments and neurite vertices (in particular, no somas). */
        void                                        getNeuriteConnectedComponent(
                                                        neurite_iterator                v_it,
                                                        std::list<NeuriteVertex *>     &reachable_neurite_vertices,
                                                        std::list<NeuriteSegment *>    &reachable_neurite_segments);

        /* statistical methods */
        CellNetworkStatistics                       computeStatistics() const;

        R                                           getTotalNeuriteTreeLength() const;

        /* get morphological diameter of cell rooted with soma vertex referred to by iterator sit */
        R                                           getMorphologicalDiameter(soma_iterator s_it);

        void                                        getNeuronVerticesDistanceStat(
                                                        R  &v_dist_min,
                                                        R  &v_dist_max,
                                                        R  &v_dist_avg,
                                                        R  &v_dist_sigma) const;

        void                                        getNeuronVerticesCoordinateStat(
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
                                                        R  &v_zsigma) const;


        void                                        getNeuriteRadiusStat(
                                                        R  &n_rmin,
                                                        R  &n_rmax,
                                                        R  &n_ravg,
                                                        R  &n_rsigma) const;

        void                                        getNeuriteSegmentRadiusStat(
                                                        R  &ns_abs_radius_diff_min,
                                                        R  &ns_abs_radius_diff_max,
                                                        R  &ns_abs_radius_diff_avg,
                                                        R  &ns_abs_radius_diff_sigma,
                                                        //
                                                        R  &ns_radius_ratio_min,
                                                        R  &ns_radius_ratio_max,
                                                        R  &ns_radius_ratio_avg,
                                                        R  &ns_radius_ratio_sigma) const;

        void                                        getNeuriteSegmentLengthStat(
                                                        R  &ns_len_min,
                                                        R  &ns_len_max,
                                                        R  &ns_len_avg,
                                                        R  &ns_len_sigma,
                                                        //
                                                        R  &ns_rlen_min,
                                                        R  &ns_rlen_max,
                                                        R  &ns_rlen_avg,
                                                        R  &ns_rlen_sigma) const;

        void                                        getNeuriteSegmentAngleStat(
                                                        R  &ns_angle_min,
                                                        R  &ns_angle_max,
                                                        R  &ns_angle_avg,
                                                        R  &ns_angle_sigma) const;


        /* transform to centroid coordinate system */
        void                                        transformToCentroidSystem();

        /* transform to coordinate system whose is the position of the soma given by s_it */
        void                                        transformToSomaSystem(soma_const_iterator const &s_it);

        Vec3<R>                                     getGlobalCoordinateDisplacement() const;
        /* -----------------------------------  I / O  ------------------------------------------------------------- */
        /* initialize network from standardized SWC file as used by NeuroMorpho.org */
        void                                        readFromNeuroMorphoSWCFile(
                                                        std::string     filename,
                                                        bool const     &check_coincident_positions = false);
};

#include "../tsrc/CellNetwork.impl.hh"

#endif
