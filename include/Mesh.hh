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

#ifndef MESH_H
#define MESH_H

#include "common.hh"
#include "Vec3.hh"
#include "BoundingBox.hh"
#include "IdQueue.hh"
#include "Octree.hh"

enum mesh_error_types {
    MESH_NOERROR,
    MESH_NOT_FOUND,
    MESH_LOGIC_ERROR,
    MESH_IO_ERROR,
    MESH_NONMANIFOLD,
    MESH_INVALID_ID
};

struct MeshEx : public std::runtime_error {
    const uint32_t      error_type;
    const std::string   error_msg;

    MeshEx(const uint32_t &type, const std::string &msg) :
        std::runtime_error(msg),
        error_type(type),
        error_msg(msg)
    {
    }
};

struct MeshEx_OutOfRange : public MeshEx {
};

template<
    typename    Tm,
    typename    Tv,
    typename    Tf,
    typename    R
>
class Mesh {
    public:
        /* ----------------- iterator related declarations ---------------------- */
        template <typename ValueType, typename InternalType>
        class MeshIterator : public std::iterator<std::bidirectional_iterator_tag, ValueType>
        {
            friend class Mesh;

            protected:
                typename InternalType::iterator     int_it;
                //ValueType                          *val;
                Mesh<Tm, Tv, Tf, R>                *mesh;

                /* private ctor, can only be called by Mesh */
                MeshIterator(
                    Mesh<Tm, Tv, Tf, R>                *m,
                    typename InternalType::iterator     it)
                : mesh(m), int_it(it)
                {}


            public:
                MeshIterator()
                : mesh(NULL)
                {}

                MeshIterator(const MeshIterator &x)
                : mesh(x.mesh), int_it(x.int_it)
                {}

               ~MeshIterator()
                {
                }

                inline MeshIterator &
                operator=(const MeshIterator &x)
                {
                    this->mesh      = x.mesh;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                    return (*this);
                }

                inline MeshIterator &
                operator++()
                {
                    ++this->int_it;
                    //this->val = ValueType::getPtr(this->int_it);
                    return (*this);
                }

                inline MeshIterator
                operator++(int)
                {
                    MeshIterator tmp(*this);
                    this->operator++();
                    return tmp;
                }

                inline MeshIterator &operator--()
                {
                    --this->int_it;
                    //this->val = ValueType::getPtr(this->int_it);
                    return (*this);
                }

                inline MeshIterator
                operator--(int)
                {
                    MeshIterator tmp(*this);
                    this->operator--();
                    return tmp;
                }

                inline bool
                operator==(const MeshIterator &x) const
                {
                    return (this->int_it == x.int_it);
                }

                inline bool
                operator!=(const MeshIterator &x) const
                {
                    return !(this->operator==(x));
                }

                /*
                inline bool
                operator<(const MeshIterator &x) const
                {
                    return (this->int_it < x.int_it);
                }

                inline bool
                operator>(const MeshIterator &x) const
                {
                    return !((*this) == x || (*this) < x);
                }

                inline bool
                operator>=(const MeshIterator &x) const
                {
                    return !( (*this) < x);
                }

                inline bool
                operator<=(const MeshIterator &x) const
                {
                    return !( (*this) > x);
                }
                */

                /*! \brief checks if (this) iterator refers to a given mesh container.
                 *
                 * \param   m   mesh container in question.
                 * \return      true iff iterator refers to m. always returns false for explicitly
                 *              invalidated iterators.
                 * */
                inline bool
                checkContainer(const Mesh<Tm, Tv, Tf, R> &m) const
                {
                    return ( this->mesh == &m);
                }

                Mesh<Tm, Tv, Tf, R> *
                getContainer()
                {
                    return (this->mesh);
                }

                /*! \brief checks if (this) iterator refers to the same mesh container as given
                 * iterator.
                 *
                 * \param   x   input iterator.
                 * \return      true iff both iterators (this) and x refer to the same mesh
                 *              container. returns false if any iterator has been explicitly
                 *              invalidated.
                 * */
                inline bool
                sameContainer(const MeshIterator &x) const
                {
                    return (!this->explicitlyInvalid() && !x.explicitlyInvalid() && this->mesh == x.mesh);
                }

                /*! \brief explicitly invalidates (this) iterator.
                 *
                 * \sideeffect (this) iterator is permanently invalidated.
                 *
                 * \par foobar\n\n lbla
                 *
                 * \note all other iterators referring to the same mesh component remain unaffected.
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
                    this->mesh = NULL;
                }

                /*! \brief check if (this) iterator has been explicitly invalidated.
                 *
                 * \returns true iff (this) iterator has been explicitly invalidated.
                 *
                 * \note this method must not be used to check for iterator validity in the
                 * general case: an iterator can in general be invalid although this method 
                 * return false. using the method in this way will result in undefined behaviour.
                 * the main use case: an algorithm or method receives an iterator (or
                 * a container with iterators) per reference, manipulates the mesh(es) the passed
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
                    return (this->mesh == NULL);
                }

                inline ValueType &
                operator*() const
                {
                    return *(ValueType::getPtr(this->int_it));
                    //return *(this->int_it->second.get());
                }

                inline ValueType *
                operator->() const
                {
                    return (ValueType::getPtr(this->int_it));
                    //return (this->int_it->second.get());
                }
        };
        /* forward declaration of classes Vertex and Face */
        class   Vertex;
        class   Face;

    private:
        /* typedefs for pointer typed internally used. could also be std::shared_ptr<..> with minor modifications */
        typedef Mesh<Tm, Tv, Tf, R>::Vertex *   VertexPointerType;    
        typedef Mesh<Tm, Tv, Tf, R>::Face *     FacePointerType;    

    public:
        /* NOTE: typdefs are not treated as full types by either the standard or compilers. for
         * more type safety, declare vertex_iterator as template specialization of MeshIterator
         * instead of typedef */
        /*
        typedef
            MeshIterator<
                Vertex,
                std::map<uint32_t, VertexPointerType >
            >   vertex_iterator;
        */

        class   vertex_iterator :
            public MeshIterator<
                Mesh<Tm, Tv, Tf, R>::Vertex,
                std::map<uint32_t, VertexPointerType>
            >
        {
            public:
                vertex_iterator()
                {
                }

                vertex_iterator(
                    Mesh<Tm, Tv, Tf, R>                                        *m,
                    typename std::map<uint32_t, VertexPointerType >::iterator   it)
                {
                    this->mesh      = m;
                    this->int_it    = it;
                    //this->val       = Vertex::getPtr(it);
                }

                vertex_iterator(const vertex_iterator &x)
                    : MeshIterator<
                        Mesh<Tm, Tv, Tf, R>::Vertex,
                        std::map<uint32_t, VertexPointerType>
                      >()
                {
                    this->mesh      = x.mesh;
                    this->int_it    = x.int_it;
                }

               ~vertex_iterator()
               {
               }
        };

        class   vertex_const_iterator :
            public MeshIterator<
                const Mesh<Tm, Tv, Tf, R>::Vertex,
                std::map<uint32_t, VertexPointerType >
            >
        {
            public:
                vertex_const_iterator()
                {
                }

                vertex_const_iterator(
                    Mesh<Tm, Tv, Tf, R>                                        *m,
                    typename std::map<uint32_t, VertexPointerType >::iterator   it)
                {
                    this->mesh      = m;
                    this->int_it    = it;
                    //this->val       = Vertex::getPtr(it);
                }
                    
                vertex_const_iterator(const vertex_const_iterator &x)
                    :  MeshIterator<
                        const Mesh<Tm, Tv, Tf, R>::Vertex,
                        std::map<uint32_t, VertexPointerType >
                    > ()
                {
                    this->mesh      = x.mesh;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

                /* provide copy constructor with non-const vertex_iterator argument to allow implicit
                 * conversion of non-const vertex_iterator to const vertex_const_iterator */
                vertex_const_iterator(const vertex_iterator &x)
                    :  MeshIterator<
                        const Mesh<Tm, Tv, Tf, R>::Vertex,
                        std::map<uint32_t, VertexPointerType >
                    > ()
                {
                    this->mesh      = x.mesh;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

               ~vertex_const_iterator()
               {
               }
        };

        /* NOTE: same as for vertex_iterator above. typedefs are not treated as full new types, but
         * rather as a sort macro, which is undesirable for type checking. fully specialize
         * MeshIterator template for const and non-const face iterators */
        /*
        typedef MeshIterator<
                Face,
                std::map<uint32_t, FacePointerType >
            > face_iterator;
        */

        class   face_iterator :
            public MeshIterator<
                Mesh<Tm, Tv, Tf, R>::Face,
                std::map<uint32_t, FacePointerType >
            >
        {
            public:
                face_iterator()
                {
                }

                face_iterator(
                    Mesh<Tm, Tv, Tf, R>                                        *m,
                    typename std::map<uint32_t, FacePointerType >::iterator     it)
                {
                    this->mesh      = m;
                    this->int_it    = it;
                    //this->val       = Vertex::getPtr(it);
                }

                face_iterator(const face_iterator &x)
                    :  MeshIterator<
                        Mesh<Tm, Tv, Tf, R>::Face,
                        std::map<uint32_t, FacePointerType >
                    > ()
                {
                    this->mesh      = x.mesh;
                    this->int_it    = x.int_it;
                }

               ~face_iterator()
               {
               }
        };

        class face_const_iterator :
            public MeshIterator<
                const Mesh<Tm, Tv, Tf, R>::Face,
                std::map<uint32_t, FacePointerType >
            >
        {
            public:
                face_const_iterator()
                {
                }

                face_const_iterator(
                    Mesh<Tm, Tv, Tf, R>                                        *m,
                    typename std::map<uint32_t, FacePointerType >::iterator     it)
                {
                    this->mesh      = m;
                    this->int_it    = it;
                    //this->val       = Face::getPtr(it);
                }

                face_const_iterator(const face_const_iterator &x)
                    : MeshIterator<
                        const Mesh<Tm, Tv, Tf, R>::Face,
                        std::map<uint32_t, FacePointerType >
                    > ()
                {
                    this->mesh      = x.mesh;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

                /* provide copy constructor with non-const vertex_iterator argument to allow implicit
                 * conversion of non-const vertex_iterator to const vertex_const_iterator */
                face_const_iterator(const face_iterator &x)
                    : MeshIterator<
                        const Mesh<Tm, Tv, Tf, R>::Face,
                        std::map<uint32_t, FacePointerType >
                    > ()
                {
                    this->mesh      = x.mesh;
                    this->int_it    = x.int_it;
                    //this->val       = x.val;
                }

               ~face_const_iterator()
                {
                }
        };

        /* forward declaration of Mesh::VertexAccessor */
        class VertexAccessor;

        /* vertex class */
        class Vertex {
            friend class Mesh<Tm, Tv, Tf, R>;
            /* NOTE: in C++11, all nested classes are friends of the containing class. furthermore,
             * they act as an implementation detail of the containing class and have access to all
             * data that the containing class has access to.
             * So: Mesh::Vertex has access to the private data of Mesh, but by default not vice
             * versa. Since Mesh has been declared a friend of Mesh::Vertex above, all nested
             * classes of Mesh have access to Mesh::Vertex as well, rendering the following 
             * delcarations obsolete */

            /* friend class Mesh::VertexAccessor; */
            /* can't declare friend with typedef vertex_iterator, so write full template
             * argument list. logically, this is equivalent to:
             *
             * friend class vertex_iterator; */
            /*
            friend class Mesh::MeshIterator<
                            Vertex,
                            std::map<uint32_t, Vertex>
                         >;
            */
            /* vertex_const_iterator is inherited from MeshIterator: this works */
            /* friend class vertex_const_iterator; */

            private:
                Mesh<Tm, Tv, Tf, R>                *mesh;
                typename std::map<
                    uint32_t,
                    VertexPointerType >::iterator   m_vit;

                Vec3<R>                             position;
                uint32_t                            current_traversal_id : 23, traversal_state : 8;
                Tv                                  data;

                std::list<Vertex *>                 adjacent_vertices;
                std::list<Face *>                   incident_faces;

                /* private ctors */
                                                    Vertex();
                                                    Vertex(
                                                        Mesh<Tm, Tv, Tf, R>    *mesh,
                                                        const                   Vec3<R> &pos,
                                                        const Tv               *data = NULL);
                                                    Vertex(const Vertex &x);

                /* private assignment operator */
                Vertex                             &operator=(const Vertex &x);

                /* Vertex objects must not be publically allocated with new() or delete() => private
                 * new and delete operators to prevent this at compile time. */
                static void                        *operator new(size_t size);
                static void                         operator delete(void *p);

                /* other private methods */
                void                                replaceAdjacentVertices(const std::map<Vertex *, Vertex*> &replace_map);

                /* static getPtr() method required by iterator */
                static Vertex *                     getPtr(typename std::map<uint32_t, VertexPointerType >::iterator it);

                /* private methods to insert / delete adjacent vertices / incident faces, which
                 * abstract from the internally used lists (used to be set for incident faces). this
                 * avoids sort() / unique() calls all over the place. isolated calls to these
                 * methods may leave the Vertex in an inconistent state: they are to be used only by
                 * the implementation and with care to retain data structure invariants. */
                void                                insertAdjacentVertex(Vertex *v);
                /* note that this deletes only ONE occurrence of the adjacent vertex. if this vertex
                 * is contained in more than one incident faces, this->isNeighbour(v) will still
                 * return true. this returns true iff v has been removed. */
                bool                                eraseAdjacentVertex(Vertex *v);
                void                                insertIncidentFace(Face *f);
                /* this removes the incident face f if present, whereupon true is returned,
                 * otherwise false. */
                bool                                eraseIncidentFace(Face *f);

            public:
                uint32_t
                id() const
                {
                    return this->m_vit->first;
                }

                Vec3<R>
                pos() const
                {
                    return this->position;
                }

                Vec3<R> &
                pos()
                {
                    return this->position;
                }

                vertex_iterator
                iterator() const
                {
                    return vertex_iterator(this->mesh, this->m_vit);
                }

                vertex_const_iterator
                const_iterator() const
                {
                    return vertex_iterator(this->mesh, this->m_vit);
                }

                uint32_t
                deg() const
                {
                    std::list<Vertex *> vstar;
                    this->getVertexStar(vstar);
                    return (vstar.size());
                }

                bool
                gotNeighbour(const vertex_const_iterator &v) const
                {
                    return this->gotNeighbour(&(*v));
                }

                bool
                gotNeighbour(const Vertex * const &v) const
                {
                    for (auto &u : this->adjacent_vertices) {
                        if (u && v->id() == u->id()) {
                            return true;
                        }
                    }
                    return false;
                }

                bool
                gotNeighbour(const uint32_t &v_id) const
                {
                    for (auto &u : this->adjacent_vertices) {
                        if (u && u->id() == v_id) {
                            return true;
                        }
                    }
                    return false;
                }

                Tv &
                operator*()
                {
                    return (this->data);
                }

                Tv *
                operator->()
                {
                    return &(this->data);
                }

                /*
                Vec3 &
                operator=(const Vec3 &p)
                {
                    this->position = p;
                    return (this->position);
                }

                Vec3 &
                operator+=(const Vec3 &p)
                {
                    this->position += p;
                    return (this->position);
                }

                Vec3 &
                operator-=(const Vec3 &p)
                {
                    this->position -= p;
                    return (this->position);
                }

                Vec3 &
                operator*=(const double &x)
                {
                    this->position *= x;
                    return (this->position);
                }

                Vec3 &
                operator/=(const double &x)
                {
                    this->position /= x;
                    return (this->position);
                }
                */

                void                                getVertexStar(std::list<Vertex *> &vstar) const;
                void                                getVertexStarIndices(std::list<uint32_t> &vstar) const;
                void                                getVertexStarIterators(std::list<vertex_iterator> &vstar) const;
                
                // way more efficient:
                const std::list<Face*>&             getFaceStar() const;
                void                                getFaceStar(std::list<Face *> &fstar) const;
                void                                getFaceStar(std::list<Face const *> &fstar) const;
                void                                getFaceStarIndices(std::list<uint32_t> &fstar) const;
                void                                getFaceStarIterators(std::list<face_iterator> &fstar) const;
                void                                getFaceStarIterators(std::list<face_const_iterator> &fstar) const;

                BoundingBox<R>                      getBoundingBox() const;
                uint8_t                             getTraversalState(const uint32_t &traversal_id);
                void                                setTraversalState(const uint32_t &traversal_id, const uint8_t &state);

                bool                                isIsolated() const;
                bool                                isManifold() const;
        };

        /* forward declaration of FaceAccessor */
        class FaceAccessor;

    public: 
        /* face class */
        class Face {
            friend class Mesh<Tm, Tv, Tf, R>;
            /* NOTE: in C++11, all nested classes are friends of the containing class. furthermore,
             * they act as an implementation detail of the containing class and have access to all
             * data that the containing class has access to.
             * So: Mesh::Vertex has access to the private data of Mesh, but by default not vice
             * versa. Since Mesh has been declared a friend of Mesh::Vertex above, all nested
             * classes of Mesh have access to Mesh::Vertex as well, rendering the following 
             * delcarations obsolete */

            /* friend class Mesh::FaceAccessor; */
            /* C++ can't declare friend with typedef face_iterator, so write full template
             * argument list. logically, this is equivalent to
             *
             * friend class face_iterator; */
            /*
            friend class Mesh::MeshIterator<
                            Face,
                            std::map<uint32_t, Face>
                         >; */
            /* face_const_iterator is inherited from MeshIterator: this works */
            /* friend class face_const_iterator; */

            private:
                Mesh<Tm, Tv, Tf, R>                *mesh;
                typename std::map<
                    uint32_t,
                    FacePointerType
                >::iterator                         m_fit;

                std::array<Vertex *,4>              vertices;
                uint32_t                            quad : 1,  current_traversal_id : 23, traversal_state : 8;
                Tf                                  data;

                /* private ctors, can only be called by Mesh and nested classes */
                                                    Face(); 
                                                    Face(
                                                        Mesh<Tm, Tv, Tf, R>    *mesh,
                                                        bool                    quad,
                                                        Vertex                 *vi,
                                                        Vertex                 *vj,
                                                        Vertex                 *vk,
                                                        Vertex                 *vl,
                                                        const Tf               *data = NULL);

                                                    Face(const Face &x);
                Face                               &operator=(const Face &b);

                /* static getPtr() method required by iterator */
                static Face *                       getPtr(typename std::map<uint32_t, FacePointerType >::iterator it);

                void                                replaceVertices(const std::map<Vertex *, Vertex*> &replace_map);
                bool                                operator<(const Face &b) const;
                bool                                operator>(const Face &b) const;

            public:
                uint32_t                            id() const;
                face_iterator                       iterator() const;
                face_const_iterator                 const_iterator() const;
                bool                                isTri() const;
                bool                                isQuad() const;

                Tf                                 &operator*();
                Tf                                 *operator->();

                uint8_t                             getTraversalState(const uint32_t &traversal_id);
                void                                setTraversalState(const uint32_t &traversal_id, const uint8_t &state);

                bool                                contains(const vertex_const_iterator &v) const;
                bool                                contains(const Vertex * const &v) const;
                bool                                contains(const uint32_t &v_id) const;

                bool                                gotNeighbour(const face_const_iterator &f) const;
                bool                                gotNeighbour(const Face * const &v) const;
                bool                                gotNeighbour(const uint32_t &v_id) const;

                void                                getIterators(std::list<vertex_iterator> &l) const;
                void                                getVertices(std::list<Vertex *> &l) const;
                void                                getIndices(uint32_t& i1, uint32_t& i2, uint32_t& i3, uint32_t& i4) const;
                std::vector<uint32_t>               getIndices() const;
                void                                getPositions(std::list<Vec3<R>> &l) const;

                Vec3<R>                             getNormal() const;
                BoundingBox<R>                      getBoundingBox() const;
                void                                getEdges(std::list<std::pair<vertex_iterator, vertex_iterator> > &edge_list) const;
                void                                getFaceNeighbours(std::vector<Face*>& nb_faces) const;
                void                                getFaceNeighbourhood(uint32_t max_depth, std::vector<Face*> &face_neighbourhood);
                void                                invertOrientation();

                /* triangle-specific methods */
                inline void                         checkTri(const char *fn) const
                                                    {
                                                        // we do not want to init a string every time we call this
                                                        // only if we really need it
                                                        if (!this->isTri())
                                                            throw MeshEx(MESH_LOGIC_ERROR, std::string(fn) + "triangle-specific method called on non-triangle face.");
                                                    }

                inline void                         checkTri(const std::string& fn) const
                                                    {
                                                        if ( !this->isTri() ) {
                                                            throw MeshEx(MESH_LOGIC_ERROR, fn + "triangle-specific method called on non-triangle face.");
                                                        }
                                                    }
                inline void                         checkQuad(const char *fn) const
                                                    {
                                                        if ( !this->isQuad() ) {
                                                            throw MeshEx(MESH_LOGIC_ERROR, std::string(fn) + "quad-specific method called on non-quad face.");
                                                        }
                                                    }

                inline void                         checkQuad(const std::string& fn) const
                                                    {
                                                        if ( !this->isQuad() ) {
                                                            throw MeshEx(MESH_LOGIC_ERROR, fn + "quad-specific method called on non-quad face.");
                                                        }
                                                    }

                inline void                         checkTriQuad(const char *fn) const
                                                    {
                                                        if ( !this->isQuad() && !this->isTri()) {
                                                            throw MeshEx(MESH_LOGIC_ERROR, std::string(fn) + "supplied face is neither quad nor triangle. general case intentionally unsupported right now => internal logic error.");
                                                        }
                                                    }

                inline void                         checkTriQuad(const std::string& fn) const
                                                    {
                                                        if ( !this->isQuad() && !this->isTri()) {
                                                            throw MeshEx(MESH_LOGIC_ERROR, fn + "supplied face is neither quad nor triangle. general case intentionally unsupported right now => internal logic error.");
                                                        }
                                                    }

                void                                getTriIterators(vertex_iterator &v0, vertex_iterator &v1, vertex_iterator &v2) const;
                void                                getTriVertices(Vertex *&v0, Vertex *&v1, Vertex *&v2) const;
                void                                getTriIndices(uint32_t &v0_id, uint32_t &v1_id, uint32_t &v2_id) const;
                void                                getTriPositions(Vec3<R> &v0, Vec3<R> &v1, Vec3<R> &v2) const;
                void                                getTriInfo(uint32_t &v0_id, Vec3<R> &v0, uint32_t &v1_id, Vec3<R> &v1, uint32_t &v2_id, Vec3<R> &v2) const;

                Vec3<R>                             getTriNormal() const;
                /* default ar def: radius of circumcircle / twice radius incircle */
                R                                   getTriAspectRatio() const;
                /* old ar def: length of longest / length of shortes edge. */
                R                                   getTriAspectRatioMaxMin() const;
                R                                   getTriArea() const;
                R                                   getTriSignedVolume() const;

                void                                getTriShortestEdge(vertex_iterator &u_it, vertex_iterator &v_it) const;
                void                                getTriShortestEdge(uint32_t &u_id, uint32_t &v_id) const;
                uint8_t                             getTriEdgeNum(uint32_t u_id, uint32_t v_id) const;
                bool                                getTriEdgeOrientation(uint32_t u, uint32_t v) const;

                void                                getTriOrthonormalBase2d(Vec3<R> &e_x, Vec3<R> &e_y) const;
                void                                getTriProjectedPointBarycentricCoordinates(Vec3<R> x, R &x_s, R &x_t) const;
                bool                                triContainsPoint(Vec3<R> x, R &x_s, R &x_t, R const &eps = 1E-15) const;

                /* useful static methods */
                static bool                         getTriEdgeOrientationIndices(
                                                        uint32_t v0_id,
                                                        uint32_t v1_id,
                                                        uint32_t v2_id,
                                                        uint32_t u,
                                                        uint32_t v);

                static bool                         getTwoTrianglesSharedEdgeAndRemainingVertices(
                                                        vertex_iterator     A_v0,
                                                        vertex_iterator     A_v1,
                                                        vertex_iterator     A_v2,
                                                        vertex_iterator     B_v0,
                                                        vertex_iterator     B_v1,
                                                        vertex_iterator     B_v2,
                                                        vertex_iterator    &fst_shared_vertex,
                                                        vertex_iterator    &snd_shared_vertex,
                                                        vertex_iterator    &A_remaining_vertex,
                                                        vertex_iterator    &B_remaining_vertex);

                static bool                         getTwoTrianglesCommonVertexAndRemainingVertices(
                                                        vertex_iterator     A_v0,
                                                        vertex_iterator     A_v1,
                                                        vertex_iterator     A_v2,
                                                        vertex_iterator     B_v0,
                                                        vertex_iterator     B_v1,
                                                        vertex_iterator     B_v2,
                                                        vertex_iterator    &shared_vertex,
                                                        vertex_iterator    &A_fst_remaining_vertex,
                                                        vertex_iterator    &A_snd_remaining_vertex,
                                                        vertex_iterator    &B_fst_remaining_vertex,
                                                        vertex_iterator    &B_snd_remaining_vertex);

                static vertex_iterator              getTriRemainingVertex(
                                                        vertex_iterator     v0_it,
                                                        vertex_iterator     v1_it,
                                                        vertex_iterator     v2_it,
                                                        vertex_iterator     u_it,
                                                        vertex_iterator     v_it);

                static void                         getTriSharedAndRemainingVertex(
                                                        vertex_iterator     u1_it,
                                                        vertex_iterator     v1_it,
                                                        vertex_iterator     u2_it,
                                                        vertex_iterator     v2_it,
                                                        vertex_iterator    &shared_vit,
                                                        vertex_iterator    &remaining_vit);

                /* quad-specific methods */
                void                                getQuadIterators(vertex_iterator &v0, vertex_iterator &v1, vertex_iterator &v2, vertex_iterator &v3) const;
                void                                getQuadVertices(Vertex *&v0, Vertex *&v1, Vertex *&v2, Vertex *&v3) const;
                void                                getQuadIndices(uint32_t &v0_id, uint32_t &v1_id, uint32_t &v2_id, uint32_t &v3_id) const;
                void                                getQuadPositions(Vec3<R> &v0, Vec3<R> &v1, Vec3<R> &v2, Vec3<R> &v3) const;
                void                                getQuadInfo(uint32_t &v0_id, Vec3<R> &v0, uint32_t &v1_id, Vec3<R> &v1, uint32_t &v2_id, Vec3<R> &v2, uint32_t &v3_id, Vec3<R> &v3) const;

        };

        /* ------------------ Accessor classes for vertices / faces ------------- */
        class VertexAccessor {
            friend class Mesh<Tm, Tv, Tf, R>;

            private:
                Mesh<Tm, Tv, Tf, R>        &mesh;

                /* private constructor, can only be called by Mesh */
                                            VertexAccessor(Mesh<Tm, Tv, Tf, R> &m);

                /* not assignable, not copy-constructible */
                VertexAccessor             &operator=(const VertexAccessor &x)      = delete;
                                            VertexAccessor(const VertexAccessor &x) = delete;


            public:
                vertex_iterator             begin();
                vertex_const_iterator       begin() const;

                vertex_iterator             end();
                vertex_const_iterator       end() const;

                vertex_iterator             find(const uint32_t &id);
                vertex_const_iterator       find(const uint32_t &id) const;

                Mesh::Vertex               &at(const uint32_t &id);
                bool                        exists(const uint32_t &id) const;

                vertex_iterator             insert(const Vec3<R> &v);
                        
                vertex_iterator             erase(vertex_iterator it);
                bool                        erase(const uint32_t &id);

                size_t                      size() const;
                bool                        empty() const;
        };

        class FaceAccessor {
            friend class Mesh<Tm, Tv, Tf, R>;

            private:
                Mesh<Tm, Tv, Tf, R>        &mesh;

                /* private constructor, can only be called by Mesh */
                                            FaceAccessor(Mesh<Tm, Tv, Tf, R> &m);
                /* not assignable, not copy-constructible */
                FaceAccessor               &operator=(const FaceAccessor &x) = delete;
                                            FaceAccessor(const FaceAccessor &x) = delete;


            public:
                face_iterator               begin();
                face_const_iterator         begin() const;

                face_iterator               end();
                face_const_iterator         end() const;

                face_iterator               find(const uint32_t &id);
                face_const_iterator         find(const uint32_t &id) const;

                Mesh::Face                 &at(const uint32_t &id);
                bool                        exists(const uint32_t &id) const;

                face_iterator               insert(const std::list<uint32_t> &id_list);

                face_iterator               insert(vertex_iterator v0_it, vertex_iterator v1_it, vertex_iterator v2_it);
                face_iterator               insert(const uint32_t &v0_id, const uint32_t &v1_id, const uint32_t &v2_id);

                face_iterator               insert(vertex_iterator v0_it, vertex_iterator v1_it, vertex_iterator v2_it, vertex_iterator v3_it);
                face_iterator               insert(const uint32_t &v0_id, const uint32_t &v1_id, const uint32_t &v2_id, const uint32_t &v3_id);

                face_iterator               erase(face_iterator it);
                bool                        erase(const uint32_t &id);

                size_t                      size() const;
                bool                        empty() const;
        };


    private:
        /* structs used for the template types of the employed Octree */
        struct
        Mesh_OctreeInfo {
        };

        struct
        Mesh_OctreeNodeInfo {
            Vec3<R>                 cube_min, cube_max;
            std::list<Face *>       face_list;
            std::list<Vertex *>     vertex_list;
        };

        struct
        Mesh_OctreeIsecNodeInfo {
            Vec3<R>                 cube_min, cube_max;
            std::list<Face *>       A_facelist;
            std::list<Face *>       B_facelist;
        };

        /* id queues */
        IdQueue                             V_idq;
        IdQueue                             F_idq;
        IdQueue                             traversal_idq;

        /* vertex and face maps */
        std::map<
                uint32_t,
                VertexPointerType
            >                               V;
        std::map<
                uint32_t,
                FacePointerType
            >                               F;

        /* data object of template type Tm */
        Tm                                  data;

        /* bounding box and octree */
        BoundingBox<R>                      bb;
        Octree<
                Mesh_OctreeInfo,
                Mesh_OctreeNodeInfo,
                R
            >                              *O;
        bool                                octree_updated;

        static void
        partitionOctree(
            Octree<Mesh_OctreeInfo, Mesh_OctreeNodeInfo, R>    &O,
            const Mesh<Tm, Tv, Tf, R>                          &m,
            OctreeNode<Mesh_OctreeNodeInfo>                    *n,
            std::list<Face *>                                  &rec_facelist,
            std::list<Vertex *>                                &rec_vertexlist,
            uint32_t                                            rec_depth,
            uint32_t                                            max_elements,
            uint32_t                                            max_rec_depth);

        void
        findOctreeRecursive(
            Octree<Mesh_OctreeInfo, Mesh_OctreeNodeInfo, R>    *O,
            OctreeNode<Mesh_OctreeNodeInfo>                    *n,
            Vec3<R>                                             aabb_min,
            Vec3<R>                                             aabb_max,
            std::list<Vertex *>                                *vertex_list,
            std::list<Face *>                                  *face_list);

        void
        getOctreeLeafCubes(
            Octree<Mesh_OctreeInfo, Mesh_OctreeNodeInfo, R>    *O,
            OctreeNode<Mesh_OctreeNodeInfo>                    *n,
            std::list<std::pair<Vec3<R>, Vec3<R>>>             *cube_list);

        /* globally reset the traversal states of all vertices and faces to TRAV_UNSEEN and reset
         * the traversal id queue. */
        void                                resetTraversalStates();

    public:
        /* Mesh public interface */

        /* the only two publically accessible members are immutable instances of the accessor classes for vertices and
         * faces. they are both non-copy-constructible and non-assignable. */
        VertexAccessor                      vertices;
        FaceAccessor                        faces;

        /* ctor, copy ctor */
                                            Mesh();
                                            Mesh(const Mesh<Tm, Tv, Tf, R> &X);
        /* assignment operator */                                    
        Mesh<Tm, Tv, Tf, R>                &operator=(const Mesh<Tm, Tv, Tf, R> &X);

        /* dtor */
                                           ~Mesh();


        /* clear all data, clear faces only */
        void                                clear();
        void                                clearFaces();
        /* renumber vertices and faces consecutively from the given start ids onwards */
        void                                renumberConsecutively(uint32_t vertex_start_id = 0, uint32_t face_start_id = 0);

        /* general info / statistics */
        uint32_t                            numVertices() const;
        uint32_t                            numEdges() const;
        uint32_t                            numFaces() const;
        uint32_t                            numObtuseTriangles() const;
        R                                   getTotalArea() const;
        R                                   getTotalVolume() const;
        void                                getAvgAspectRatio(R &ar_avg, R &ar_sigma, R *ar_max = NULL, R *ar_min = NULL) const;
        BoundingBox<R>                      getBoundingBox() const;

        /* invert orientation of all faces */
        void                                invertOrientation();

        /* traveral state is publically accessible for traversal methods. */
        enum MeshTraversalStates {
            TRAV_UNSEEN      = 0,
            TRAV_ENQUEUED    = 1,
            TRAV_DONE        = 2,
            TRAV_BLOCKED     = 3
        };
        uint32_t                            getFreshTraversalId();

        /* append another mesh: add vertices / faces and offset indices. not that no topological connection
         * between (this) mesh and the appended mesh is performed, it simply computes the union of two
         * distinct meshes inside one object. note also that the old ids (and of course all
         * iterators) to B are NO LONGER valid in the appended mesh, i.e. ids from B cannot be used
         * to refer to the copy of B that has been appended to mesh in question. */
        void                                copyAppend(const Mesh<Tm, Tv, Tf, R> &B);

        /* move version: append mesh B to (this) mesh by moving B's contents. B is empty afterwards.
         * since pointers are used internally, this does not require significant additional memory
         * and can be done "in-place".  the parameter update_vits is a pointer to a list containing
         * vertex iterators referring to B, which will be updated in-place with iterators referring
         * to the respective vertices of (this) mesh after the append. */
        void                                moveAppend(
                                                Mesh<Tm, Tv, Tf, R>            &B,
                                                std::list<vertex_iterator>     *update_vits = NULL);

        /* delete the connected component (i.e. all reachable vertices / faces) of the vertex identified by iterator
         * vstart_it */
        void                                deleteConnectedComponent(vertex_iterator vstart_it);

        /* get all vertices / faces / edges of the connected component of a given start vertex for a
         * given traversal id. note that "boundary conditions" can be implemeneted by manually
         * setting the traversal state of selected vertices for the supplied traversal id before
         * calling this method. example: given a "ring" of vertices whose deletion separates the
         * mesh in two connected components, one can set the traversal state of all these ring
         * vertices to DONE, which causes them to act as a boundary that the traversal cannot pass
         * */
        void                                getConnectedComponent(
                                                Mesh::vertex_iterator   vstart_it,    
                                                const uint32_t         &traversal_id,
                                                std::list<Vertex *>    *cc_vertices = NULL,
                                                std::list<Face *>      *cc_faces    = NULL);

        /* delete all connected components containing isolated vertices or border edges */
        void                                deleteBorderCCsAndIsolatedVertices();

        /* ----------------- location routines using the octree ----------------- */
        void                                updateOctree();
        void                                findVertices(
                                                BoundingBox<R> const   &search_box,
                                                std::list<Vertex *>    &vertex_list);

        void                                findFaces(
                                                BoundingBox<R> const  &search_box,
                                                std::list<Face *>      &face_list);

        void                                find(
                                                BoundingBox<R> const   &search_box,
                                                std::list<Vertex *>    *vertex_list,
                                                std::list<Face *>      *face_list);

        /* ----------------- selection-related methods -------------------------- */
        void                                selectNonManifoldVertices(std::list<Vertex *> &vlist) const;
        void                                selectIsolatedVertices(std::list<Vertex *> &vlist) const;
        void                                invertVertexSelection(std::list<Vertex *> &vlist) const;
        void                                invertFaceSelection(std::list<Face *> &flist) const;

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

        void                                getFacesIncidentToEdge(
                                                const vertex_const_iterator    &u_it,
                                                const vertex_const_iterator    &v_it,
                                                Face**                         incident_faces,
                                                size_t&                        sizeInOut) const;

        void                                getFacesIncidentToEdge(
                                                uint32_t                u,
                                                uint32_t                v,
                                                Face**                  incident_faces,
                                                size_t&                 sizeInOut) const;
        /* version for manifold meshes, where each edge is incident to exactly two faces.
         * throws exception if this is violated */
        void                                getFacesIncidentToManifoldEdge(
                                                const vertex_const_iterator    &u_it,
                                                const vertex_const_iterator    &v_it,
                                                face_iterator                  &fst_face,
                                                face_iterator                  &snd_face) const;

        face_iterator                       getOtherFaceIncidentToManifoldEdge(
                                                const vertex_const_iterator    &u_it,
                                                const vertex_const_iterator    &v_it,
                                                const face_const_iterator      &known_face_it) const;

        bool                                getTriCommonEdge(
                                                const face_const_iterator      &A_it,
                                                const face_const_iterator      &B_it,
                                                vertex_iterator                &u_it,
                                                vertex_iterator                &v_it) const;

        bool                                getTriCommonEdge(
                                                uint32_t    A_id,
                                                uint32_t    B_id,
                                                uint32_t   &u,
                                                uint32_t   &v) const;

        void                                getEdgeBoundingBox(
                                                uint32_t    u,
                                                uint32_t    v,
                                                Vec3<R>    &aabb_min,
                                                Vec3<R>    &aabb_max) const;

        face_iterator                       getSharedTri(
                                                vertex_const_iterator const    &u1_it,
                                                vertex_const_iterator const    &v1_it,
                                                vertex_const_iterator const    &u2_it,
                                                vertex_const_iterator const    &v2_it) const;

        /* ---------------------- topological / geometric modifications ------------------------- */
        void                                scale(R const &r);
        void                                translate(Vec3<R> const &d);
        /* given an edge and a relative coordinate lambda, insert a new vertex on that edge, delete both
         * incident triangles and the four resultling new triangles */
        void                                splitEdge(
                                                vertex_iterator     u_it,
                                                vertex_iterator     v_it,
                                                R const            &lambda);

        /* given an edge and a relative coordinate lambda, insert a new vertex on that edge, delete both
         * incident triangles and the four resultling new triangles */
        void                                splitEdge(
                                                uint32_t            u_it,
                                                uint32_t            v_it,
                                                R const            &lambda);

        void                                splitEdge(
                                                vertex_iterator         u_it,
                                                vertex_iterator         v_it,
                                                std::vector<R> const   &lambda_values);

        void                                split_face_with_center(
                                                uint32_t face_id,
                                                const Vec3<R>& splitPos);

        void                                splitQuadIntoTwoTris(face_iterator quad);
        void                                splitQuadIntoTwoTrisShorterDiagonal(face_iterator quad);
        void                                triangulateQuads();

        /* given an edge {u,v} and a new vertex position, collapse it, i.e.
         *      1. delete the two incident triangles, create a new vertex 
         *      2. create a new vertex at the given position
         *      3. replace the indices of u and v in all other faces formlery incident to u or v
         *      with the index of the new vertex
         *  
         * if pos == NULL, default to the midpoint between u and v. the return value indicates whether
         * the collapse has beeen performed. it is performed only if it is topologically safe, i.e. the
         * topology of the mesh is not affected. if the mesh was manifold before, it is manifold after. 
         * it is safe to collapse the edge {u,v} if |nbs(u) \cap nbs{v}| == 2. if not, false is
         * returned. */
        bool                                collapseTriEdge(
                                                vertex_iterator     u_it,
                                                vertex_iterator     v_it,
                                                vertex_iterator    *vnew_it = NULL,
                                                Vec3<R>            *new_pos = NULL);

        /* given two vertices u and v which are topologically unconnected in the following sense:
         *
         * 1. neither one is a neighbour of the other, i.e. the edge {u, v} does not exist.
         * 2. neither one is contained in a face that is incident to the other.
         *
         * this methods "merges" these two vertices together to produce a vertex w, whose
         *
         * 1. vertex star is the union of the vertex stars of u and v
         * 2. face star is the union of the face stars of u and v, where all ocurances of u and v in
         * incident faces are replaced by w. */
        vertex_iterator                     mergeUnrelatedVertices(
                                                vertex_iterator     u_it,
                                                vertex_iterator     v_it,
                                                const Vec3<R>      *new_pos = NULL);


        /* ----------------- topology / geometry check methods  ----------------- */
        void                                checkTopology();
        void                                checkGeometry(
                                                R const &eps_v = 1E-3,
                                                R const &eps_T = 1E-3);

        /* debug method: check internal data consistency */
        void                                checkInternalConsistency() const;

        /* -----------------  I/O  ----------------- */
        void                                readFromObjFile(const char *filename);
        void                                writeObjFile(const char *jobname);


        /* NOTE: In the C++11 standard, nested classes are automatically "friends" of the containing
         * class, but not vice versa. the declarations below are therefore obsolete */
        /* "friend class vertex_iterator;" cannot be used since C++ does not allow friend
         * declaration through typdefs. use full template argument list */
        /*
        friend class Mesh::VertexAccessor;
        friend class Mesh::FaceAccessor;
        */
        /*
        friend class Mesh::MeshIterator<
                        Vertex,
                        std::map<uint32_t, Vertex>
                     >;
        friend class vertex_const_iterator;
        */

        /* "friend class face_iterator;" cannot be used since C++ does not allow friend
         * declaration through typdefs. use full template argument list */
        /*
        friend class Mesh::MeshIterator<
                        Face,
                        std::map<uint32_t, Face>
                     >;
        friend class face_const_iterator;
        */
};

/* include header for template implementation */
#include "../tsrc/Mesh.impl.hh"

#endif
