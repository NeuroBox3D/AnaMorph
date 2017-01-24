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
 *                       mesh vertex class implementation..                                                        
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* mesh vertex ctors */
/* NOTE: Vertex::m_vit is not set after these constructors have been called and thus the vertex is
 * not in a consistent state. the caller (which has to be a friend of Mesh::Vertex, since the
 * following ctor are private) needs to take care to set the iterator correctly */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Vertex::Vertex()
{
    this->mesh                  = NULL;
    this->position              = Aux::VecMat::nullvec<R>;
    this->current_traversal_id  = 0;
    this->traversal_state       = TRAV_UNSEEN;
}

template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Vertex::Vertex(
    Mesh               *mesh,
    Vec3<R> const      &pos,
    Tv const           *data)
{
    this->mesh                  = mesh;
    this->position              = pos;
    this->current_traversal_id  = 0;
    this->traversal_state       = TRAV_UNSEEN;
    if (data) {
        this->data              = *data;
    }
}

/* private copy ctor */
/* NOTE: this constructor is private (i.e. not publicly accessible) and is only used internally for
 * the implementation of Mesh methods.  anway: careful when changing the internals, pointers are
 * copied!  if this is constructor is invoked to change the id of a face within a mesh, that's ok.
 * if however this is used to copy vertices from one mesh to another, the caller must take care to
 * adjust the adjacency / incidence information accordingly to avoid UB. */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Vertex::Vertex(const Vertex &x)
{    
    this->mesh                  = x.mesh;
    this->m_vit                 = x.m_vit;
    this->position              = x.position;
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->adjacent_vertices     = x.adjacent_vertices;
    this->incident_faces        = x.incident_faces;
    this->data                  = x.data;
}

/* private assignment operator. same as for private copy ctor: to be used only internally and only
 * with care. */
template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::Vertex &
Mesh<Tm, Tv, Tf, R>::Vertex::operator=(const Vertex &x)
{
    this->mesh                  = x.mesh;
    this->m_vit                 = x.m_vit;
    this->position              = x.position;
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->adjacent_vertices     = x.adjacent_vertices;
    this->incident_faces        = x.incident_faces;
    this->data                  = x.data;

    return (*this);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void *
Mesh<Tm, Tv, Tf, R>::Vertex::operator new(size_t size)
{
    if (size != sizeof(Mesh::Vertex)) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::Vertex::operator new(): size does not match sizeof(Vertex). internal logic error.");
    }
    return ::operator new(size);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::operator delete(void *p)
{
    ::operator delete(p);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::replaceAdjacentVertices(const std::map<Vertex *, Vertex*> &replace_map)
{
    typename std::list<Vertex *>::iterator                   nbit;
    typename std::map<Vertex *, Vertex *>::const_iterator    mit;

    for (nbit = this->adjacent_vertices.begin(); nbit != this->adjacent_vertices.end(); ++nbit) {
        /* if nb is is found in index change map, replace it */
        if ( (mit = replace_map.find(*nbit)) != replace_map.end() ) {
            *nbit = mit->second;
        }
    }

    this->adjacent_vertices.sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::Vertex *
Mesh<Tm, Tv, Tf, R>::Vertex::getPtr(
    typename std::map<uint32_t, VertexPointerType >::iterator it)
{
    return (it->second);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::insertAdjacentVertex(Vertex *v)
{
    Aux::Alg::listSortedInsert(this->adjacent_vertices, v, true);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Vertex::eraseAdjacentVertex(Vertex *v)
{
    return Aux::Alg::removeFirstOccurrenceFromList(this->adjacent_vertices, v);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::insertIncidentFace(Face *f)
{
    Aux::Alg::listSortedInsert(this->incident_faces, f, false);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Vertex::eraseIncidentFace(Face *f)
{
    return Aux::Alg::removeFirstOccurrenceFromList(this->incident_faces, f);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getVertexStar(
    std::list<Vertex *> &vstar) const
{
    /* NOTE: assigning will assign to or destroy all values stored in vstar when it is given to this
     * method. since pointers have trivial operator=(), this is semantically equivalent to
     * clear()ing vstar and copying the list with std::copy. */
    vstar = this->adjacent_vertices;
    vstar.sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});
    vstar.unique([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() == y->id());});
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getVertexStarIndices(
    std::list<uint32_t> &vstar) const
{
    vstar.clear();
    for (auto &nb : this->adjacent_vertices) {
        vstar.push_back(nb->id());
    }
    vstar.sort();
    vstar.unique();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getVertexStarIterators(
    std::list<vertex_iterator> &vstar) const
{
    auto av_copy = this->adjacent_vertices;
    av_copy.sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});
    av_copy.unique([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() == y->id());});

    vstar.clear();
    for (auto &nb : av_copy) {
        vstar.push_back(nb->iterator());
    }
}


template <typename Tm, typename Tv, typename Tf, typename R>
const std::list<typename Mesh<Tm, Tv, Tf, R>::Face*>&
Mesh<Tm, Tv, Tf, R>::Vertex::getFaceStar() const
{
    return incident_faces;
}


template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getFaceStar(std::list<Face *> &fstar) const
{
    fstar.clear();
    std::copy(this->incident_faces.begin(), this->incident_faces.end(), std::back_inserter(fstar) );
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getFaceStar(std::list<Face const *> &fstar) const
{
    fstar.clear();
    std::copy(this->incident_faces.begin(), this->incident_faces.end(), std::back_inserter(fstar) );
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getFaceStarIndices(std::list<uint32_t> &fstar) const
{
    fstar.clear();
    for (auto &f : this->incident_faces) {
        fstar.push_back(f->id());
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getFaceStarIterators(std::list<face_iterator> &fstar) const
{
    fstar.clear();
    for (auto &f : this->incident_faces) {
        fstar.push_back(f->iterator());
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::getFaceStarIterators(std::list<face_const_iterator> &fstar) const
{
    fstar.clear();
    for (auto &f : this->incident_faces) {
        fstar.push_back(f->iterator());
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
BoundingBox<R>
Mesh<Tm, Tv, Tf, R>::Vertex::getBoundingBox() const
{
    /* offset vector in every dimension */
    Vec3<R> offset(1E-3, 1E-3, 1E-3);

    /* initialize bounding box */
    Vec3<R> bb_min = this->pos(), bb_max = this->pos();
    bb_min -= offset;
    bb_max += offset;

    return BoundingBox<R>( { bb_min, bb_max });
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint8_t
Mesh<Tm, Tv, Tf, R>::Vertex::getTraversalState(const uint32_t &traversal_id)
{
    if (this->current_traversal_id != traversal_id) {
        this->current_traversal_id  = traversal_id;
        this->traversal_state       = TRAV_UNSEEN;
    }
    return this->traversal_state;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Vertex::setTraversalState(
    const uint32_t &traversal_id,
    const uint8_t  &state)
{
    this->current_traversal_id  = traversal_id;
    this->traversal_state       = state;
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Vertex::isIsolated() const
{
    return (this->deg() == 0);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Vertex::isManifold() const
{
    std::list<Mesh<Tm, Tv, Tf, R>::Face const *>    fstar;
    std::list<Mesh<Tm, Tv, Tf, R>::Face const *>    dual_circle;
    Mesh<Tm, Tv, Tf, R>::Face const *               last_on_circle;
    bool                                            got_neighbour;
    Mesh<Tm, Tv, Tf, R>::vertex_iterator            e_x, e_y;

    /* get face star */
    this->getFaceStar(fstar);

    /* if face star consists of < 3 face, the vertex is non-manifold */
    if (fstar.size() < 3) {
        return false;
    }

    /* get off first face and push it as first dual "vertex", i.e. primal face, in  dual_circle */
    dual_circle.push_back(fstar.front());
    fstar.pop_front();

    /* always take last element from dual circle, search for neighbour  */
    while ( !fstar.empty() ) {
        last_on_circle  = dual_circle.back();
        got_neighbour   = false;

        /* scan remaining incident faces for neighbour of last_on_circle */
        for (auto ifit = fstar.begin(); ifit != fstar.end(); ++ifit) {
            /* if *ifit is a neighbour, remove it from list and append it to dual cycle */
            if ( (this->mesh->getTriCommonEdge(last_on_circle->iterator(), (*ifit)->iterator(), e_x, e_y)) == true) {
                /* one of e_x and e_y must be v itself */
                if (e_x->id() != this->id() && e_y->id() != this->id() ) {
                    debugl(1, "Mesh::Vertex::isManifold(): face in face star of vertex %d shares common edge with other face from face star, but vertex %d is not contained in in.\n", this->id(), this->id() );
                    throw("Mesh::Vertex::isManifold(): face in face star of vertex shares common edge with other face from face star, but vertex v is not contained in in.\n");
                }

                /* got the neighbour, delete *ifit and push_back on dual circle */
                got_neighbour = true;

                dual_circle.push_back(*ifit);
                fstar.erase(ifit);
                break;
            }
        }

        /* incident_faces non-empty and no closed cycle => non-manifold vertex. */
        if (!got_neighbour) {
            return false;
        }
    }

    /* check if last element of v_dual_circle shares edge with first element of v_dual_circle */
    if (!this->mesh->getTriCommonEdge(dual_circle.front()->iterator(), dual_circle.back()->iterator(), e_x, e_y)) {
        return false;
    }
    else {
        /* (this) vertex is regular or "manifold" */
        return true;
    }
}


/* ----------------------------------------------------------------------------------------------------------------- *
 *
 *                       mesh face class implementation..                                                        
 *
 * ----------------------------------------------------------------------------------------------------------------- */

/* iterate through adjacent vertices and replace ids using the given map */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::replaceVertices(const std::map<Vertex *, Vertex*> &replace_map)
{
    typename std::map<Vertex *, Vertex *>::const_iterator mit;

    /* for all four vertex pointers: if not NULL, search in replace_map and replace if found */
    for (int i = 0; i < 4; i++) {
        if (this->vertices[i]) {
            if ( (mit = replace_map.find(this->vertices[i])) != replace_map.end() ) {
                this->vertices[i] = mit->second;
            }
        }
    }
}


/* meshface ctors */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Face::Face()
{
    this->mesh                  = NULL;
    this->current_traversal_id  = 0;
    this->vertices.fill(NULL);
    this->quad                  = false;
    this->current_traversal_id  = 0;
    this->traversal_state       = TRAV_UNSEEN;
}

/* NOTE: private ctor only to be used by the internal implementation, which default constructs the
 * internal iterator Mesh::m_fit. the (internal) caller must set the iterator (once it is known) to
 * bring the Face into a consistent state. */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Face::Face(
    Mesh       *mesh,
    bool       quad,
    Vertex    *vi,
    Vertex    *vj,
    Vertex    *vk,
    Vertex    *vl,
    const Tf  *data)
{
    this->mesh                  = mesh;
    this->quad                  = quad;
    this->vertices[0]           = vi;
    this->vertices[1]           = vj;
    this->vertices[2]           = vk;
    this->vertices[3]           = vl;
    this->current_traversal_id  = 0;
    this->traversal_state       = TRAV_UNSEEN;
    if (data) {
        this->data              = *data;
    }
}

/* private copy ctor */
/* NOTE: this constructor is private (i.e. not publicly accessible) and is only used internally for
 * the implementation of Mesh methods.  anway: careful when changing the internals, pointers are
 * copied!  if this is constructor is invoked to change the id of a face within a mesh, that's ok.
 * if however this is used to copy vertices from one mesh to another, the caller must take care to
 * adjust the adjacency / incidence information accordingly to avoid UB. */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Face::Face(const Face &x)
{
    this->mesh                  = x.mesh;
    this->m_fit                 = x.m_fit;
    this->quad                  = x.quad;
    for (int i = 0; i < 4; i++) {
        this->vertices[i]       = x.vertices[i];
    }
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->data                  = x.data;
}

/* private assignment operator. same as for private copy ctor: to be used only internally and only
 * with care. */
template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::Face &
Mesh<Tm, Tv, Tf, R>::Face::operator=(const Mesh::Face &x)
{
    this->mesh                  = x.mesh;
    this->m_fit                 = x.m_fit;
    this->quad                  = x.quad;
    for (int i = 0; i < 4; i++) {
        this->vertices[i]       = x.vertices[i];
    }
    this->current_traversal_id  = x.current_traversal_id;
    this->traversal_state       = x.traversal_state;
    this->data                  = x.data;

    return (*this);
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::Face *
Mesh<Tm, Tv, Tf, R>::Face::getPtr(typename std::map<uint32_t, FacePointerType >::iterator it)
{
    return (it->second);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::operator<(const Face &b) const
{
    this->checkTriQuad("Mesh::Face::operator<()");

    /* triangles are always smaller than quads */
    if (this->isTri() && b.isQuad()) {
        return true;
    }
    else if (this->isQuad() && b.isTri()) {
        return false;
    }
    /* either both triangles or both quads */
    else {
        /* get two lists of indices, get minimum, cycle until minimum values are in front() of both
         * lists */

        if (this->isTri() && b.isTri())
        {
        	uint32_t i[4];
            getIndices(i[0], i[1], i[2], i[3]);
            if (i[1] <= i[0] && i[1] <= i[2])
            {
            	i[3] = i[0];	// use i3 as tmp
            	i[0] = i[1];
            	i[1] = i[2];
            	i[2] = i[3];
            }
            else if (i[2] <= i[0] && i[2] <= i[1])
            {
            	i[3] = i[2];
				i[2] = i[1];
				i[1] = i[0];
				i[0] = i[3];
            }

        	uint32_t j[4];
            b.getIndices(j[0], j[1], j[2], j[3]);
            if (j[1] <= j[0] && j[1] <= j[2])
            {
            	j[3] = j[0];	// use j3 as tmp
            	j[0] = j[1];
            	j[1] = j[2];
            	j[2] = j[3];
            }
            else if (j[2] <= j[0] && j[2] <= j[1])
            {
            	j[3] = j[2];
				j[2] = j[1];
				j[1] = j[0];
				j[0] = j[3];
            }

            /* component-wise compare */
            return (i[0] < j[0] && i[1] < j[1] && i[2] < j[2]);
        }
        else if (this->isQuad() && b.isQuad()) {
            throw MeshEx(MESH_LOGIC_ERROR, "Face::operator<(): don't use me on quads yet .__^");
        }
        else throw MeshEx(MESH_LOGIC_ERROR, "Face::operator<(): one of the two faces is neither a quad nor a triangle. general case intentionally unsupported right now => internal logic error.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::operator>(const Face &b) const
{
    this->checkTriQuad("Mesh::Face::operator>()");

    /* triangles are always smaller than quads */
    if (this->isTri() && b.isQuad()) {
        return false;
    }
    else if (this->isQuad() && b.isTri()) {
        return true;
    }
    /* either both triangles or both quads */
    else {
        /* get two lists of indices, get minimum, cycle until minimum values are in front() of both
         * lists */
        if (this->isTri() && b.isTri())
        {
        	uint32_t i[4];
			getIndices(i[0], i[1], i[2], i[3]);
			if (i[1] <= i[0] && i[1] <= i[2])
			{
				i[3] = i[0];	// use i3 as tmp
				i[0] = i[1];
				i[1] = i[2];
				i[2] = i[3];
			}
			else if (i[2] <= i[0] && i[2] <= i[1])
			{
				i[3] = i[2];
				i[2] = i[1];
				i[1] = i[0];
				i[0] = i[3];
			}

			uint32_t j[4];
			b.getIndices(j[0], j[1], j[2], j[3]);
			if (j[1] <= j[0] && j[1] <= j[2])
			{
				j[3] = j[0];	// use j3 as tmp
				j[0] = j[1];
				j[1] = j[2];
				j[2] = j[3];
			}
			else if (j[2] <= j[0] && j[2] <= j[1])
			{
				j[3] = j[2];
				j[2] = j[1];
				j[1] = j[0];
				j[0] = j[3];
			}

			/* component-wise compare */
			return (i[0] > j[0] && i[1] > j[1] && i[2] > j[2]);
        }
        else if (this->isQuad() && b.isQuad()) {
            throw MeshEx(MESH_LOGIC_ERROR, "Face::operator>(): don't use me on quads yet .__^");
        }
        else throw MeshEx(MESH_LOGIC_ERROR, "Face::operator>(): one of the two faces is neither a quad nor a triangle. general case intentionally unsupported right now => internal logic error.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint8_t
Mesh<Tm, Tv, Tf, R>::Face::getTraversalState(const uint32_t &traversal_id)
{
    if (this->current_traversal_id != traversal_id) {
        this->current_traversal_id  = traversal_id;
        this->traversal_state       = TRAV_UNSEEN;
    }
    return this->traversal_state;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::setTraversalState(
    const uint32_t &traversal_id,
    const uint8_t  &state)
{
    this->current_traversal_id  = traversal_id;
    this->traversal_state       = state;
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint32_t
Mesh<Tm, Tv, Tf, R>::Face::id() const
{
    return this->m_fit->first;
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::Face::iterator() const
{
    return face_iterator(this->mesh, this->m_fit );
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_const_iterator
Mesh<Tm, Tv, Tf, R>::Face::const_iterator() const
{
    return face_iterator(this->mesh, this->m_fit );
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::isTri() const
{
    return (this->quad == 0);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::isQuad() const
{
    return (this->quad == 1);
}

template <typename Tm, typename Tv, typename Tf, typename R>
Tf &
Mesh<Tm, Tv, Tf, R>::Face::operator*()
{
    return (this->data);
}

template <typename Tm, typename Tv, typename Tf, typename R>
Tf *
Mesh<Tm, Tv, Tf, R>::Face::operator->()
{
    return &(this->data);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::contains(const vertex_const_iterator &v) const
{
    this->checkTriQuad("Mesh::Face::contains()");
    return (this->contains(&(*v)));
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::contains(const Vertex * const &v) const
{
    this->checkTriQuad("Mesh::Face::contains()");

    for (auto &u : this->vertices) {
        if (u && v == u) {
            return true;
        }
    }
    return false;
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::contains(const uint32_t &v_id) const
{
    this->checkTriQuad("Mesh::Face::contains()");

    for (auto &u : this->vertices) {
        if (u && v_id == u->id()) {
            return true;
        }
    }
    return false;
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::gotNeighbour(const face_const_iterator &f) const
{
    return this->gotNeighbour(&(*f));
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::gotNeighbour(const Face * const &f) const
{
    std::vector<Face*> f_nbs;
    f_nbs.reserve(4);
    this->getFaceNeighbours(f_nbs);

    for (auto &g : f_nbs) {
        if (g && g == f) {
            return true;
        }
    }
    return false;
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::gotNeighbour(const uint32_t &f_id) const
{
    std::vector<Face*> f_nbs;
    f_nbs.reserve(4);
    this->getFaceNeighbours(f_nbs);

    for (auto &g : f_nbs) {
        if (g->id() == f_id) {
            return true;
        }
    }
    return false;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getVertices(
    std::list<Mesh::Vertex *> &l) const
{
    this->checkTriQuad("Mesh::Face::getVertices()");

    l.clear();
    l.push_back(this->vertices[0]);
    l.push_back(this->vertices[1]);
    l.push_back(this->vertices[2]);
    if (this->isQuad()) {
        l.push_back(this->vertices[3]);
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
std::vector<uint32_t>
Mesh<Tm, Tv, Tf, R>::Face::getIndices() const
{
    this->checkTriQuad("Mesh::Face::getIndices()");

    std::vector<uint32_t> l;
    l.push_back(this->vertices[0]->id());
    l.push_back(this->vertices[1]->id());
    l.push_back(this->vertices[2]->id());
    if (this->isQuad()) {
        l.push_back(this->vertices[3]->id());
    }

    return l;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getIndices(uint32_t& i1, uint32_t& i2, uint32_t& i3, uint32_t& i4) const
{
    this->checkTriQuad("Mesh::Face::getIndices()");

    i1 = vertices[0]->id();
    i2 = vertices[1]->id();
    i3 = vertices[2]->id();
    if (this->isQuad())
        i4 = vertices[3]->id();
    else
    	i4 = i3;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getPositions(std::list<Vec3<R>> &l) const
{
    this->checkTriQuad("Mesh::Face::getIndices()");

    l.clear();
    l.push_back(this->vertices[0]->pos());
    l.push_back(this->vertices[1]->pos());
    l.push_back(this->vertices[2]->pos());
    if (this->isQuad()) {
        l.push_back(this->vertices[3]->pos());
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
Vec3<R>
Mesh<Tm, Tv, Tf, R>::Face::getNormal() const
{
    this->checkTriQuad("Mesh::Face::getNormal()");

    if (this->isQuad()) {
        Vec3<R> v0, v1, v2, v3, n1, n2, n;
        this->getQuadPositions(v0, v1, v2, v3);

        n1  = (v1 - v0).cross(v2 - v0);
        n1.normalize();
        n2  = (v2 - v0).cross(v3 - v0);
        n2.normalize();
        n   = n1 + n2;
        n.normalize();

        return n;
    }
    /* triangle */
    else {
        Vec3<R> v0, v1, v2, n;
        this->getTriPositions(v0, v1, v2);

        n = (v1 - v0).cross(v2 - v0);
        n.normalize();
        
        return n;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
BoundingBox<R>
Mesh<Tm, Tv, Tf, R>::Face::getBoundingBox() const
{ 
    using Aux::VecMat::minVec3;
    using Aux::VecMat::maxVec3;
    using Aux::VecMat::fabsVec3;

    this->checkTriQuad("Mesh::Face::getBoundingBox()");

    debugl(5, "Mesh::Face::getBoundingBox()\n");
    //debugTabInc();

    Vec3<R> v_i, v_j, v_k, v_l;

    if (this->isTri()) {
        this->getTriPositions(v_i, v_j, v_k);
    }
    else {
        this->getQuadPositions(v_i, v_j, v_k, v_l);
    }

    Vec3<R> bb_min  = v_i;
    minVec3(bb_min, v_i, v_j);
    minVec3(bb_min, bb_min, v_k);
    if (this->isQuad()) {
        minVec3(bb_min, bb_min, v_l);
    }

    Vec3<R> bb_max  = v_i;
    maxVec3(bb_max, bb_max, v_j);
    maxVec3(bb_max, bb_max, v_k);
    if (this->isQuad()) {
    	maxVec3(bb_max, bb_max, v_l);
    }

    /* extend bounding box by 2.5%, but no less than 1E-3, on each side. */
    return (BoundingBox<R>(bb_min, bb_max).extend(0.025, Vec3<R>(1E-3, 1E-3, 1E-3)));

    // code never gets here!?s
    debugl(5, "aabb (widened): (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f)\n",
            bb_min[0], bb_min[1], bb_min[2],
            bb_max[0], bb_max[1], bb_max[2]);

    debugTabDec();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getEdges(
    std::list<std::pair<vertex_iterator, vertex_iterator> > &edge_list) const
{
    this->checkTriQuad("Mesh::Face::getEdges");

    /* clear return list, handle two cases */
    edge_list.clear();
    if (this->isTri()) {
        edge_list.push_back( { this->vertices[0]->iterator(), this->vertices[1]->iterator() });
        edge_list.push_back( { this->vertices[1]->iterator(), this->vertices[2]->iterator() });
        edge_list.push_back( { this->vertices[2]->iterator(), this->vertices[0]->iterator() });
    }
    else {
        edge_list.push_back( { this->vertices[0]->iterator(), this->vertices[1]->iterator() });
        edge_list.push_back( { this->vertices[1]->iterator(), this->vertices[2]->iterator() });
        edge_list.push_back( { this->vertices[2]->iterator(), this->vertices[3]->iterator() });
        edge_list.push_back( { this->vertices[3]->iterator(), this->vertices[0]->iterator() });
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getFaceNeighbours(
    std::vector<Face*>& nb_faces) const
{
    this->checkTriQuad("Mesh::Face::getFaceNeighbours()");

    Face* thisFace = Face::getPtr(this->m_fit);
    Face* tmp[2];
    if (this->isTri()) {
        vertex_const_iterator   v_it[3];
        for (int i = 0; i < 3; i++) {
            v_it[i] = this->vertices[i]->iterator();
        }

        /* get all faces incident to edges (0, 1), (1, 2) and (2,0) */
        for (size_t i = 0; i < 3; ++i)
        {
            size_t nIF = 2;
            this->mesh->getFacesIncidentToEdge(v_it[i], v_it[(i+1)%3], tmp, nIF);
            for (size_t j = 0; j < nIF; ++j)
                if (tmp[j] != thisFace)
                    nb_faces.push_back(tmp[j]);
        }
    }
    /* same for quads */
    else {
        vertex_const_iterator   v_it[4];
        for (int i = 0; i < 4; i++) {
            v_it[i] = this->vertices[i]->iterator();
        }

        for (size_t i = 0; i < 4; ++i)
        {
            size_t nIF = 2;
            this->mesh->getFacesIncidentToEdge(v_it[i], v_it[(i+1)%4], tmp, nIF);
            for (size_t j = 0; j < nIF; ++j)
                if (tmp[j] != thisFace)
                    nb_faces.push_back(tmp[j]);
        }
    }

    // not required, is it?
    /*
    std::sort(nb_faces.begin(), nb_faces.end(), Mesh::Face::ptr_less);
    auto newEnd = std::unique(nb_faces.begin(), nb_faces.end());
    nb_faces.erase(newEnd, nb_faces.end());
    */
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getFaceNeighbourhood(
    uint32_t            max_depth,
    std::vector<Face*>  &face_neighbourhood)
{
    std::queue< std::pair<Face*, uint32_t> >    Q;
    Face                                       *f;
    uint32_t                                    f_depth;
    std::vector<Face*>                          f_neighbours;
    f_neighbours.reserve(3);

    /* get new traversal id */
    const uint32_t                              traversal_id = this->mesh->getFreshTraversalId();

    // each new element can have two new neighbors at most (triangles)
    // this leads to a maximum of 3*2^n-2 triangles in the neighborhood
    // be careful if depth parameter is chosen ridiculously high
    face_neighbourhood.clear();
    face_neighbourhood.reserve(3*pow(2,std::min(max_depth,7u))-2);

    /* initialize Q: we start at (this) face in depth 0 */
    Q.push({this, 0});

    //printf("go..\n");
    while (!Q.empty())
    {
        //printf("Q.size(): %d\n", Q.size());
        /* get front element and dequeue */
        f       = Q.front().first;
        f_depth = Q.front().second;
        Q.pop();

        /* append to face_list, set f's traversal state to done */
        face_neighbourhood.push_back(f);
        f->setTraversalState(traversal_id, TRAV_DONE);

        /* if max_depth has not yet been reached, inspect all neighbours of f, if they haven't been
         * enqueue yet, enqueue them. */
        if (f_depth <= max_depth) {
            f->getFaceNeighbours(f_neighbours);
            for (Face *nb : f_neighbours) {
                if (nb->getTraversalState(traversal_id) == TRAV_UNSEEN) {
                    Q.push({nb, f_depth + 1});
                    nb->setTraversalState(traversal_id, TRAV_ENQUEUED);
                }
            }
            f_neighbours.clear();
        }
    }

    // why would this be necessary!?
    /*
    std::sort(face_neighbourhood.begin(), face_neighbourhood.end(), Mesh::Face::ptr_less);
    auto newEnd = std::unique(face_neighbourhood.begin(), face_neighbourhood.end());
    face_neighbourhood.erase(newEnd, face_neighbourhood.end());
    */
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::invertOrientation()
{
    this->checkTriQuad("Mesh::Face::invertOrientation()");

    /* for triangles, vertices[3] contains NULL, reverse only first three elements of array */
    if (this->isTri()) {
        std::reverse(this->vertices.begin(), this->vertices.begin() + 3);
    }
    /* for quads, rever entire vertices array */
    else {
        std::reverse(this->vertices.begin(), this->vertices.end()); 
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
R
Mesh<Tm, Tv, Tf, R>::Face::getTriAspectRatio() const
{
    this->checkTri("Mesh::Face::getTriAspectRatio():");

    Vec3<R> v0, v1, v2;
    R       a, b, c, s;

    this->getTriPositions(v0, v1, v2);
    a       = (v1 - v0).len2();
    b       = (v2 - v1).len2();
    c       = (v0 - v2).len2();
    s       = 0.5*(a + b + c);

    return ( (a*b*c) / ( 8.0*(s - a)*(s - b)*(s - c) ));
}

template <typename Tm, typename Tv, typename Tf, typename R>
R
Mesh<Tm, Tv, Tf, R>::Face::getTriAspectRatioMaxMin() const
{
    this->checkTri("Mesh::Face::getTriAspectRatioMaxMin():");

    Vec3<R> v0, v1, v2;
    R       l01, l12, l20, len_min, len_max;

    this->getTriPositions(v0, v1, v2);
    l01     = (v1 - v0).len2();
    l12     = (v2 - v1).len2();
    l20     = (v0 - v2).len2();

    len_min = l01;
    if (l12 < len_min) {
        len_min = l12;
    }
    if (l20 < len_min) {
        len_min = l20;
    }

    len_max = l01;
    if (l12 > len_max) {
        len_max = l12;
    }
    if (l20 > len_max) {
        len_max = l20;
    }

    return (len_max / len_min);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriIterators(
    vertex_iterator    &v0,
    vertex_iterator    &v1,
    vertex_iterator    &v2) const
{
    this->checkTri("Mesh::Face::getTriIterators()");
    v0  = this->vertices[0]->iterator();
    v1  = this->vertices[1]->iterator();
    v2  = this->vertices[2]->iterator();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriVertices(
    Mesh::Vertex *&v0,
    Mesh::Vertex *&v1,
    Mesh::Vertex *&v2) const
{
    this->checkTri("Mesh::Face::getTriVertices()");
    v0  = this->vertices[0];
    v1  = this->vertices[1];
    v2  = this->vertices[2];
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriIndices(
    uint32_t &v0_id,
    uint32_t &v1_id,
    uint32_t &v2_id) const
{
    this->checkTri("Mesh::Face::getTriIndices()");
    v0_id   = this->vertices[0]->id();
    v1_id   = this->vertices[1]->id();
    v2_id   = this->vertices[2]->id();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriPositions(
    Vec3<R> &v0,
    Vec3<R> &v1,
    Vec3<R> &v2) const
{
    this->checkTri("Mesh::Face::getTriPositions()");
    v0  = this->vertices[0]->pos();
    v1  = this->vertices[1]->pos();
    v2  = this->vertices[2]->pos();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriInfo(
    uint32_t   &v0_id,
    Vec3<R>    &v0, 
    uint32_t   &v1_id,
    Vec3<R>    &v1, 
    uint32_t   &v2_id,
    Vec3<R>    &v2) const
{
    this->getTriIndices(v0_id, v1_id, v2_id);
    this->getTriPositions(v0, v1, v2);
}

template <typename Tm, typename Tv, typename Tf, typename R>
Vec3<R>
Mesh<Tm, Tv, Tf, R>::Face::getTriNormal() const
{
    this->checkTri("Mesh::Face::getTriNormal():");

    Vec3<R> v0, v1, v2, n;
    this->getTriPositions(v0, v1, v2);

    n = (v1 - v0).cross(v2 - v0);
    n.normalize();
    
    return n;
}

template <typename Tm, typename Tv, typename Tf, typename R>
R
Mesh<Tm, Tv, Tf, R>::Face::getTriArea() const
{
    this->checkTri("Mesh::Face::getTriVertices():");

    Vec3<R> v0, v1, v2;
    this->getTriPositions(v0, v1, v2);
    return ( ((v1 - v0).cross(v2 - v0)).len2() / 2.0 );
}

template <typename Tm, typename Tv, typename Tf, typename R>
R
Mesh<Tm, Tv, Tf, R>::Face::getTriSignedVolume() const
{
    this->checkTri("Mesh::Face::getTriVertices():");

    Vec3<R> v0, v1, v2;
    this->getTriPositions(v0, v1, v2);
    return ( (v0 * v1.cross(v2)) / 6.0 );
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriShortestEdge(
    uint32_t &u_id,
    uint32_t &v_id) const
{
    this->checkTri("Mesh::Face::getTriShortestEdge():");

    Vec3<R>     v0, v1, v2;
    uint32_t    v0_id, v1_id, v2_id;

    this->getTriInfo(v0_id, v0, v1_id, v1, v2_id, v2);

    R       l01, l12, l20, len_min;

    l01     = (v1 - v0).len2();
    l12     = (v2 - v1).len2();
    l20     = (v0 - v2).len2();

    len_min = l01;
    len_min = std::min(len_min, l12);
    len_min = std::min(len_min, l20);

    if (len_min == l01) {
        u_id = v0_id; 
        v_id = v1_id;
    }
    else if (len_min == l12) {
        u_id = v1_id;
        v_id = v2_id;
    }
    else {
        u_id = v2_id;
        v_id = v0_id;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriShortestEdge(
    vertex_iterator &u_it,
    vertex_iterator &v_it) const
{
    this->checkTri("Mesh::Face::getTriShortestEdge():");

    vertex_iterator v0, v1, v2;
    R               l01, l12, l20, len_min;

    this->getTriIterators(v0, v1, v2);

    l01     = (v1->pos() - v0->pos()).len2();
    l12     = (v2->pos() - v1->pos()).len2();
    l20     = (v0->pos() - v2->pos()).len2();

    len_min = l01;
    len_min = std::min(len_min, l12);
    len_min = std::min(len_min, l20);

    if (len_min == l01) {
        u_it = v0;
        v_it = v1;
    }
    else if (len_min == l12) {
        u_it = v1;
        v_it = v2;
    }
    else {
        u_it = v2;
        v_it = v0;
    }
}

/* get edge number in {0,1,2} of triangle, given an edge (u_id, v_id) */
template <typename Tm, typename Tv, typename Tf, typename R>
uint8_t
Mesh<Tm, Tv, Tf, R>::Face::getTriEdgeNum(
    uint32_t    u_id,
    uint32_t    v_id) const
{
    this->checkTri("Mesh::Face::getTriEdgeNum():");

    uint32_t v0_id, v1_id, v2_id;
    this->getTriIndices(v0_id, v1_id, v2_id);

    if ( (u_id == v0_id && v_id == v1_id) || (u_id == v1_id && v_id == v0_id) ) {
        return 0;
    }
    else if ( (u_id == v1_id && v_id == v2_id) || (u_id == v2_id && v_id == v1_id) ) {
        return 1;
    }
    else if ( (u_id == v2_id && v_id == v0_id) || (u_id == v0_id && v_id == v2_id) ) {
        return 2;
    }
    else {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::Face::getTriEdgeNum(): given triangle does not contain oriented edge (u_id, v_id).");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::getTriEdgeOrientation(
    uint32_t u,
    uint32_t v) const
{
    this->checkTri("Mesh::Face::getTriEdgeOrientation():");

    uint32_t v0, v1, v2;
    this->getTriIndices(v0, v1, v2); 

    /* check various cases */
    if (u == v0 && v == v1) {
        return true;
    }
    else if (u == v1 && v == v0) {
        return false;
    }
    else if (u == v1 && v == v2) {
        return true;
    }
    else if (u == v2 && v == v1) {
        return false;
    }
    else if (u == v2 && v == v0) {
        return true;
    }
    else if (u == v0 && v == v2) {
        return false;
    }
    else {
        throw MeshEx(MESH_NOT_FOUND, "Mesh:getEdgeOrientationInTri(): given edge is not contained in supplied triangle.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriOrthonormalBase2d(
    Vec3<R> &e_x,
    Vec3<R> &e_y) const
{
    /* compute u, v and normal vector as above */
    this->checkTri("Mesh::getTriOrthonormalBase2d():");

    Vec3<R> v0, v1, v2, u, v, n;
    this->getTriPositions(v0, v1, v2);
    u       = v1 - v0;
    v       = v2 - v0;
    n       = u.cross(v);

    /* normalize both n and u */
    n.normalize();
    u.normalize();

    /* now e_x is normalized u, e_y is n x e_x */
    e_x     = u;
    e_y     = n.cross(e_x);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriProjectedPointBarycentricCoordinates(
    Vec3<R> x,
    R      &x_s,
    R      &x_t) const
{
    this->checkTri("Mesh::Face::getTriProjectedPointBarycentricCoordinates()");

    Vec3<R> v0, v1, v2;
    this->getTriPositions(v0, v1, v2);
    Aux::Geometry::computeBaryCoordsOfProjectedPoint(x, v0, v1, v2, x_s, x_t);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::triContainsPoint(
    Vec3<R>     x,
    R          &x_s,
    R          &x_t,
    R const    &eps) const
{
    this->checkTri("Mesh::Face::triContainsPoint()");

    R       s, t;
    Vec3<R> v0, v1, v2;

    /* get vertices and compute bary coords of x */
    this->getTriPositions(v0, v1, v2);
    Aux::Geometry::computeBaryCoordsOfProjectedPoint(x, v0, v1, v2, s, t);

    /* check barycentric coordiantes */
    if ( !std::isfinite(s) || !std::isfinite(t) || s < -eps || t < -eps || s + t > 1.0 + 2.0*eps) {
        x_s = s;
        x_t = t;
        return false;
    }
    /* point is in, return true and sert barycentric coordinates for further checking */
    else {
        x_s = s;
        x_t = t;
        return true;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getQuadIterators(
    vertex_iterator    &v0,
    vertex_iterator    &v1,
    vertex_iterator    &v2,
    vertex_iterator    &v3) const
{
    this->checkQuad("Mesh::Face::getQuadIterators()");

    v0  = this->vertices[0]->iterator();
    v1  = this->vertices[1]->iterator();
    v2  = this->vertices[2]->iterator();
    v3  = this->vertices[3]->iterator();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getQuadVertices(
    Vertex    *&v0,
    Vertex    *&v1,
    Vertex    *&v2,
    Vertex    *&v3) const
{
    this->checkQuad("Mesh::Face::getQuadVertices()");

    v0  = this->vertices[0];
    v1  = this->vertices[1];
    v2  = this->vertices[2];
    v3  = this->vertices[3];
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getQuadIndices(
    uint32_t &v0_id,
    uint32_t &v1_id,
    uint32_t &v2_id,
    uint32_t &v3_id) const
{
    this->checkQuad("Mesh::Face::getQuadIndices()");

    v0_id   = this->vertices[0]->id();
    v1_id   = this->vertices[1]->id();
    v2_id   = this->vertices[2]->id();
    v3_id   = this->vertices[3]->id();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getQuadPositions(
    Vec3<R> &v0,
    Vec3<R> &v1,
    Vec3<R> &v2,
    Vec3<R> &v3) const
{
    this->checkQuad("Mesh::Face::getQuadPositions()");

    v0  = this->vertices[0]->pos();
    v1  = this->vertices[1]->pos();
    v2  = this->vertices[2]->pos();
    v3  = this->vertices[3]->pos();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getQuadInfo(
    uint32_t   &v0_id,
    Vec3<R>    &v0,
    uint32_t   &v1_id,
    Vec3<R>    &v1,
    uint32_t   &v2_id,
    Vec3<R>    &v2,
    uint32_t   &v3_id,
    Vec3<R>    &v3) const

{
    this->checkQuad("Mesh::Face::getQuadInfo()");

    this->getQuadIndices(v0_id, v1_id, v2_id, v3_id);
    this->getQuadPositions(v0, v1, v2, v3);
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::getTriEdgeOrientationIndices(
    uint32_t v0,
    uint32_t v1,
    uint32_t v2,
    uint32_t u,
    uint32_t v)
{
    /* check various cases */
    if (u == v0 && v == v1) {
        return true;
    }
    else if (u == v1 && v == v0) {
        return false;
    }
    else if (u == v1 && v == v2) {
        return true;
    }
    else if (u == v2 && v == v1) {
        return false;
    }
    else if (u == v2 && v == v0) {
        return true;
    }
    else if (u == v0 && v == v2) {
        return false;
    }
    else {
        throw MeshEx(MESH_NOT_FOUND, "Mesh:getEdgeOrientationInTri(): given edge is not contained in supplied triangle.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::getTwoTrianglesSharedEdgeAndRemainingVertices(
    vertex_iterator     A_v0,
    vertex_iterator     A_v1,
    vertex_iterator     A_v2,
    vertex_iterator     B_v0,
    vertex_iterator     B_v1,
    vertex_iterator     B_v2,
    vertex_iterator    &fst_shared_vertex,
    vertex_iterator    &snd_shared_vertex,
    vertex_iterator    &A_remaining_vertex,
    vertex_iterator    &B_remaining_vertex)
{
    /* init pointers lists by taking address of references returned by
     * vertex_iterator::operator*(). this does not create copies */
    std::vector<Vertex *>   A_vertices = { &(*A_v0), &(*A_v1), &(*A_v2) };
    std::vector<Vertex *>   B_vertices = { &(*B_v0), &(*B_v1), &(*B_v2) };
    std::vector<Vertex *>   AB_shared_vertices;
    std::vector<Vertex *>   A_remaining_vertices, B_remaining_vertices;

    /* first, check if triangles share exactly one common edge */
    auto vrtCmp = [] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());};
    std::sort(A_vertices.begin(), A_vertices.end(), vrtCmp);
    std::sort(B_vertices.begin(), B_vertices.end(), vrtCmp);
    std::set_intersection(
            A_vertices.begin(),
            A_vertices.end(),
            B_vertices.begin(),
            B_vertices.end(),
            std::inserter(AB_shared_vertices, AB_shared_vertices.begin()),
            vrtCmp
        );

    /* if both triangles share exactly two vertices, i.e. a common edge */
    if (AB_shared_vertices.size() == 2) {
        /* get remaining vertices by set difference. the created lists should have size == 1
         * always */
        std::sort(AB_shared_vertices.begin(), AB_shared_vertices.end(), vrtCmp);
        std::set_difference(
                A_vertices.begin(),
                A_vertices.end(),
                AB_shared_vertices.begin(),
                AB_shared_vertices.end(),
                std::inserter(A_remaining_vertices, A_remaining_vertices.begin()),
                vrtCmp
            );

        if (A_remaining_vertices.size() != 1) {
            throw("(static) Mesh::getTwoTrianglesSharedEdgeAndRemainingVertices(): two triangles share one edge, but A_remaining_vertices vector has size != 1. should be impossible..");
        }

        std::set_difference(
                B_vertices.begin(),
                B_vertices.end(),
                AB_shared_vertices.begin(),
                AB_shared_vertices.end(),
                std::inserter(B_remaining_vertices, B_remaining_vertices.begin()),
                vrtCmp
            );

        if (B_remaining_vertices.size() != 1) {
            throw("(static) Mesh::getTwoTrianglesSharedEdgeAndRemainingVertices(): two triangles share one edge, but B_remaining_vertices vector has size != 1. should be impossible..");
        }

        /* write return variables */
        fst_shared_vertex   = AB_shared_vertices.front()->iterator();
        snd_shared_vertex   = AB_shared_vertices.back()->iterator();
        A_remaining_vertex  = A_remaining_vertices.front()->iterator();
        B_remaining_vertex  = B_remaining_vertices.front()->iterator();

        return true;
    }
    /* not edge-neighbours, return false */
    else return false;
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::Face::getTwoTrianglesCommonVertexAndRemainingVertices(
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
    vertex_iterator    &B_snd_remaining_vertex)
{
    /* init pointers lists by taking address of references returned by vertex_iterator::operator*().
     * this does not create copies */
    std::vector<Vertex *>   A_vertices = { &(*A_v0), &(*A_v1), &(*A_v2) };
    std::vector<Vertex *>   B_vertices = { &(*B_v0), &(*B_v1), &(*B_v2) };
    std::vector<Vertex *>   AB_shared_vertices;
    std::vector<Vertex *>   A_remaining_vertices, B_remaining_vertices;

    auto vrtCmp = [] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());};
    std::sort(A_vertices.begin(), A_vertices.end(), vrtCmp);
    std::sort(B_vertices.begin(), B_vertices.end()), vrtCmp;
    std::set_intersection(
            A_vertices.begin(),
            A_vertices.end(),
            B_vertices.begin(),
            B_vertices.end(),
            std::inserter(AB_shared_vertices, AB_shared_vertices.begin()),
            vrtCmp
        );

    if (AB_shared_vertices.size() == 1) {
        /* get remaining vertices for A and B */
        std::sort(AB_shared_vertices.begin(), AB_shared_vertices.end(), vrtCmp);
        std::set_difference(
                A_vertices.begin(),
                A_vertices.end(),
                AB_shared_vertices.begin(),
                AB_shared_vertices.end(),
                std::inserter(A_remaining_vertices, A_remaining_vertices.begin()),
                vrtCmp
            );

        if (A_remaining_vertices.size() != 2) {
            throw("(static) Mesh::getTwoTrianglesCommonVertexAndRemainingVertices(): two triangles share only one vertex, but A_remaining_vertices vector has size != 2. should be impossible..");
        }

        std::set_difference(
                B_vertices.begin(),
                B_vertices.end(),
                AB_shared_vertices.begin(),
                AB_shared_vertices.end(),
                std::inserter(B_remaining_vertices, B_remaining_vertices.begin()),
                vrtCmp
            );

        if (B_remaining_vertices.size() != 2) {
            throw("(static) Mesh::getTwoTrianglesCommonVertexAndRemainingVertices(): two triangles share only one vertex, but B_remaining_vertices vector has size != 2. should be impossible..");
        }

        /* assign return values */
        shared_vertex           = AB_shared_vertices.front()->iterator();
        A_fst_remaining_vertex  = A_remaining_vertices.front()->iterator();
        A_snd_remaining_vertex  = A_remaining_vertices.back()->iterator();

        B_fst_remaining_vertex  = B_remaining_vertices.front()->iterator();
        B_snd_remaining_vertex  = B_remaining_vertices.back()->iterator();

        return true;
    }
    else return false;
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
Mesh<Tm, Tv, Tf, R>::Face::getTriRemainingVertex(
    vertex_iterator     v0_it,
    vertex_iterator     v1_it,
    vertex_iterator     v2_it,
    vertex_iterator     u_it,
    vertex_iterator     v_it)
{
    if ((u_it->id() == v1_it->id() && v_it->id() == v2_it->id())
        || (v_it->id() == v1_it->id() && u_it->id() == v2_it->id()))
        return v0_it;

    if ((u_it->id() == v0_it->id() && v_it->id() == v2_it->id())
        || (v_it->id() == v0_it->id() && u_it->id() == v2_it->id()))
        return v1_it;

    if ((u_it->id() == v0_it->id() && v_it->id() == v1_it->id())
        || (v_it->id() == v0_it->id() && u_it->id() == v1_it->id()))
        return v2_it;

    throw MeshEx(MESH_LOGIC_ERROR, "(static) getTriRemainingVertex(): seems like given edge is not part of given triangle.");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::Face::getTriSharedAndRemainingVertex(
    vertex_iterator     u1_it,
    vertex_iterator     v1_it,
    vertex_iterator     u2_it,
    vertex_iterator     v2_it,
    vertex_iterator    &shared_vit,
    vertex_iterator    &remaining_vit)
{
   if (u1_it == u2_it) {
       shared_vit       = u1_it;
       remaining_vit    = v2_it;
   }
   else if (u1_it == v2_it) {
       shared_vit       = u1_it;
       remaining_vit    = u2_it;
   }
   else if (v1_it == u2_it) {
       shared_vit       = v1_it;
       remaining_vit    = v2_it;
   }
   else if (v1_it == v2_it) {
       shared_vit       = v1_it;
       remaining_vit    = u2_it;
   }
   else {
       throw MeshEx(MESH_LOGIC_ERROR, "(static void) Mesh::getTriSharedAndRemainingVertex(): given to edges do not share a vertex.");
   }
}
/* mesh ctors */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Mesh() : vertices(*this) , faces(*this)
{
    this->O                 = NULL;
    this->octree_updated    = false;
}

/* copy ctor */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::Mesh(const Mesh &X) : vertices(*this), faces(*this) {
    /* default init */
    this->O                 = NULL;
    this->octree_updated    = false;

    /* use assignment operator. although this initializes all members with the default ctor and
     * immediately overwrites them again, this was deemed preferable to copying the code of
     * Mesh::operator=() and have virtually the same piece of code twice.  the alternative: put the
     * code in the copy ctor and use a swap trick in operator=(), which was deemed too inefficient,
     * since (partially quite large) meshes are often assigned, but not so much copy-constructed.
     * furthermore, the default construction is very little work. */
    this->operator=(X);
}

/* assignment operator: special care has to be taken, since the topology information inside a mesh
 * uses iterators and pointers, which must of course not be copied during assignment to prevent
 * pointer aliasing and races. */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R> &
Mesh<Tm, Tv, Tf, R>::operator=(const Mesh &X)
{
	//static int calls = 0;
    //std::cout << "Mesh::operator= call " << ++calls << std::endl;

    /* clear all data */
    this->clear();

    /* assign idqs, ids used in X will be reused in (this) mesh. this can be relied upon by the caller. */
    this->V_idq             = X.V_idq;
    this->F_idq             = X.F_idq;
    this->traversal_idq     = X.traversal_idq;

    /* copy vertex map by value: this will temporarily copy pointers referring to data from X into
     * (this) mesh, which is not desirable in the final result. however, all ids are copied
     * correctly. in the next step, "deep copy" all elements in this->V, i.e. allocate new vertices
     * and copy-assign the vertices from X. note that the link between (this) mesh and X still isn't
     * broken, since the Vertex objects contains pointers, which still refer to components of X. */
    this->V                 = X.V;

    /* deep copy */
    Vertex *v_new;
    typename std::map<uint32_t, VertexPointerType >::iterator vit;
    for (vit = this->V.begin(); vit != this->V.end(); ++vit) {

        /* make a copy of the Vertex object currently pointed to by vit, which is a Vertex object
         * allocated by X, with the private copy ctor of Vertex */
        v_new       = new Vertex(*(vit->second)); 

        /* store pointer to new copy in iterator */
        vit->second = VertexPointerType(v_new);

        /* although the vertex has been deep-copied, it still contains pointers referring to X.
         * remove entirely. all topological information for (this) mesh will be constructed will be
         * constructed below */

        /* firstly, clear adjacency and incidence information */
        v_new->adjacent_vertices.clear();
        v_new->incident_faces.clear();

        /* set mesh pointer to (this) */
        v_new->mesh     = this;

        /* set iterator to vit */
        v_new->m_vit    = vit;
    }

    /* now the faces of X are added using the integer _id_ based version, not copy ctor or iterator
     * based version. through this process, all topological information (pointers, etc) is
     * constructed correctly for (this) mesh. */

    /* ids of vertices are copied, ids of faces are invalidated during the process. why not simply iterator of X.F
     * instead of the accessor, extract the id of the face as well, locate vertices via id and generate face map entry
     * in (this) mesh. this way, both vertex and face ids will be intact, which is very convenient in some applications:
     * the caller can copy a mesh and re-use ids, although pointers and iterators are invalidated. in case of moving,
     * ids are invalidated, by pointers stay intact. easy association of mesh components simplified complex processed
     * significantly.. */
    face_iterator       fit;
	uint32_t i1, i2, i3, i4;
    for (auto &xf : X.faces)
    {
        xf.getIndices(i1, i2, i3, i4);
        if (i3 == i4) // triangle
        	fit = faces.insert(i1, i2, i3);
        else // quadrilateral
        	fit = faces.insert(i1, i2, i3, i4);

        /* copy the rest of the information from xf */
        fit->traversal_state    = xf.traversal_state;
    }

    /* copy mesh boudning box */
    this->bb                = X.bb;

    /* we could copy the entire octee by value, but well... */
    this->O                 = NULL;
    this->octree_updated    = false;

    return (*this);
}

/* dtor */
template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::~Mesh()
{
    if (this->O != NULL) {
        delete this->O;
    }

    /* delete all allocated vertices */
    for (auto &v : this->vertices) {
        delete (&v);
    }

    /* delete all allocated faces */
    for (auto &f : this->faces) {
        delete (&f);
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::clear()
{
    /* delete all allocated vertices */
    for (auto &v : this->vertices) {
        delete (&v);
    }

    /* delete all allocated faces */
    for (auto &f : this->faces) {
        delete (&f);
    }

    /* clear vertex / face maps */
    this->V.clear();
    this->F.clear();

    /* clear id queues */
    this->V_idq.clear();
    this->F_idq.clear();
    this->traversal_idq.clear();

    /* delete octree */
    if (this->O != NULL) {
        delete this->O;
    }
    this->O                 = NULL;
    this->octree_updated    = false;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::clearFaces()
{
    /* delete all allocated faces */
    for (auto &f : this->faces) {
        delete (&f);
    }

    /* clear faces map and face id queue.*/
    this->F.clear();
    this->F_idq.clear();

    /* since there are no isolated edges, simply clear all adjacency and incidence information in
     * all vertices. */
    for (auto &vit : this->V) {
        vit.second->adjacent_vertices.clear();
        vit.second->incident_faces.clear();
    }
}

/* renumber vertices and faces consecutively, starting from vertex_start_id for vertices and
 * vface_start_id for faces. */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::renumberConsecutively(
    uint32_t vertex_start_id,
    uint32_t face_start_id)
{
    debugl(2, "Mesh::renumberConsecutively(): vertex_start_id: %5d, face_start_id: %5d.\n", vertex_start_id, face_start_id);
    debugTabInc();

    /* renumber indices of vertices and faces in a consecutive fashion. since the internal maps
     * Mesh::V and Mesh::E store pointers to the data, the data itself is not touched by
     * manipulation of these maps and hence all topological information inside Vertex and Face
     * objects remains valid EXCEPT the internal interators, which have to be reset to bring the
     * entire mesh back to a consistent state.
     *
     * in general, all iterators are invalidated by this method */
    std::map<uint32_t, VertexPointerType >    vertices_swap; 
    std::map<uint32_t, FacePointerType >      faces_swap;
    
    /* swap vertices and faces with vertices_swap / faces_swap in-place */
    this->V.swap(vertices_swap);
    this->F.swap(faces_swap);

    /* reset id queues for vertices / faces to start at vertex_start_id / face_start_id */
    this->V_idq = IdQueue(vertex_start_id);
    this->F_idq = IdQueue(face_start_id);

    /* iterate through swap arrays and insert Vertex and Face shared pointers into now empty 
     * maps this->V and this->F with correct ids */
    typename std::map<uint32_t, VertexPointerType >::iterator   vit, vnew_it;
    Vertex *v;
    for (uint32_t current_vertex_id = vertex_start_id; !vertices_swap.empty(); current_vertex_id++) {
        vit         = vertices_swap.begin();
        v           = Vertex::getPtr(vit);

        /* insert pointer to vertex (vit->second) into this->V with correct id, store new
         * iterator, which is the first element of the result (iterator, bool) pair. */
        vnew_it     = ( this->V.insert( {current_vertex_id, VertexPointerType(vit->second) } ) ).first;

        /* update iterator in Vertex object */
        v->m_vit    = vnew_it;

        /* erase vertex from vertices_swap */
        vertices_swap.erase(vit);

        /* pop() id from id queue. this is guaranteed to be consecutive since IdQueue object has
         * been initialized above and returned consecutive ids as long as no IdQuee::freeId() is
         * ever issued */
        this->V_idq.getId();
    }

    /* same for all faces */
    typename std::map<uint32_t, FacePointerType >::iterator fit, fnew_it;
    Face *f;

    for (uint32_t current_face_id = face_start_id; !faces_swap.empty(); current_face_id++) {
        fit         = faces_swap.begin();
        f           = Face::getPtr(fit);
        fnew_it     = ( this->F.insert( {current_face_id, FacePointerType(fit->second) } ) ).first;
        f->m_fit    = fnew_it;
        faces_swap.erase(fit);
        this->F_idq.getId();
    }

    debugTabDec();
    debugl(2, "Mesh::renumberConsecutively(). done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint32_t
Mesh<Tm, Tv, Tf, R>::numVertices() const
{
    return this->vertices.size();
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint32_t
Mesh<Tm, Tv, Tf, R>::numFaces() const
{
    return (this->faces.size());
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint32_t
Mesh<Tm, Tv, Tf, R>::numEdges() const
{
    uint32_t two_E = 0;
    for (auto &v : this->vertices) {
        two_E += v.deg();
    }
    return (two_E / 2);
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint32_t
Mesh<Tm, Tv, Tf, R>::numObtuseTriangles() const
{
    Vec3<R>     v0, v1, v2;
    R           a, b, c;
    R           aa, bb, cc;

    uint32_t    nobtuse_tris = 0;
    for (auto &f : this->faces) {
        f.checkTri("Mesh::numObtuseTriangles():");
        f.getTriPositions(v0, v1, v2);

        a       = (v1 - v0).len2();
        b       = (v2 - v1).len2();
        c       = (v0 - v2).len2();
        aa      = a*a;
        bb      = b*b;
        cc      = c*c;

        if (aa + bb < cc || bb + cc < aa || cc + aa < bb) {
            nobtuse_tris++;
        }
    }
    return nobtuse_tris;
}

/* area, volume, ar statistics */
template <typename Tm, typename Tv, typename Tf, typename R>
R     
Mesh<Tm, Tv, Tf, R>::getTotalArea() const
{
    R      area = 0.0;
    for (auto &f : this->faces) {
        f.checkTri("Mesh::getTotalArea():");
        area += f.getTriArea();
    }
    return area;
}

template <typename Tm, typename Tv, typename Tf, typename R>
R     
Mesh<Tm, Tv, Tf, R>::getTotalVolume() const
{
    R      volume = 0.0;
    for (auto &f : this->faces) {
        f.checkTri("Mesh::getTotalVolume():");
        volume += f.getTriSignedVolume();
    }
    return fabs(volume);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::getAvgAspectRatio(
    R       &ar_avg,
    R       &ar_sigma,
    R       *ar_max,
    R       *ar_min) const
{
    using Aux::Numbers::inf;

    std::vector<R     > ratios;

    for (auto &f : this->faces) {
        f.checkTri("Mesh::getAvgAspectRatio():");
        ratios.push_back( f.getTriAspectRatio() );
    }

    R      ar_min_tmp, ar_max_tmp;
    Aux::Stat::computeMinMaxAvgSigma<R     >(
            ratios, 
            ar_min_tmp,
            ar_max_tmp,
            ar_avg,
            ar_sigma);

    if (ar_max) *ar_max = ar_max_tmp;
    if (ar_min) *ar_min = ar_min_tmp;
}

template <typename Tm, typename Tv, typename Tf, typename R>
BoundingBox<R>
Mesh<Tm, Tv, Tf, R>::getBoundingBox() const
{
    using Aux::Numbers::inf;
    using Aux::VecMat::minVec3;
    using Aux::VecMat::maxVec3;
    
    if (this->octree_updated) {
        return (this->bb);
    }
    else {
        BoundingBox<R> bb;

        for (auto &v : this->vertices) {
            bb.update(v.getBoundingBox());
        }

        for (auto &f : this->faces) {
            bb.update(f.getBoundingBox());
        }

        return (bb.extend(0.025, Vec3<R>(1E-3, 1E-3, 1E-3)));
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::invertOrientation()
{
    for (auto &f : this->faces) {
        f.invertOrientation();
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::copyAppend(const Mesh &B)
{
    std::map<uint32_t, vertex_iterator>                     new_vertex_its;
    typename std::map<uint32_t, vertex_iterator>::iterator  idit;
    vertex_iterator                                         new_it;

    debugl(4, "Mesh::appendCopy()\n");
    debugTabInc();
    /* add all vertices of B to (this) mesh, store iterators to new vertices */
    for (auto &B_v : B.vertices) {
        new_it  = this->vertices.insert( B_v.pos() );

        if ( (idit = new_vertex_its.find( B_v.id() )) != new_vertex_its.end()) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::appendCopy(): duplicate entry in new_vertex_ids map, should be impossible.");
        }
        new_vertex_its[B_v.id()] = new_it;  

        debugl(5, "added vertex %5d from b under new id %5d. position: \n", B_v.id(), new_it->id() );
        B_v.pos().print_debugl(5);
    }

    /* add all faces of B to (this) mesh while taking care to replace old vertex ids from B with the
     * proper new ids */
    uint32_t v0_id, v1_id, v2_id, v3_id;
    for (auto &B_f : B.faces) {
        if (B_f.isQuad()) {
            B_f.getQuadIndices(v0_id, v1_id, v2_id, v3_id);
            this->faces.insert(
                new_vertex_its[v0_id],
                new_vertex_its[v1_id],
                new_vertex_its[v2_id],
                new_vertex_its[v3_id]);
        }
        else if (B_f.isTri()) {
            B_f.getTriIndices(v0_id, v1_id, v2_id);
            this->faces.insert(
                new_vertex_its[v0_id],
                new_vertex_its[v1_id],
                new_vertex_its[v2_id]);
        }
        else throw MeshEx(MESH_LOGIC_ERROR, "Mesh::copyAppend(): found face that is neither quad nor triangle. general case intentionally unsupported right now => internal logic error.");
    }

    this->octree_updated    = false;

    debugTabDec();
    debugl(4, "Mesh::appendCopy(): done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::moveAppend(
    Mesh                           &B,
    std::list<vertex_iterator>     *update_vits)
{
    uint32_t                            new_id;
    Vertex                             *v;
    typename std::map<
            uint32_t,
            VertexPointerType
        >::iterator                     v_newit;
    Face                               *f;
    typename std::map<
            uint32_t,
            FacePointerType
        >::iterator                     f_newit;
    bool                                inserted;
        

    debugl(2, "Mesh::moveAppend()\n");
    debugTabInc();

    /* if update_vits != NULL, generate the list of vertex pointers from the given iterators.
     * although iterators are invalidated, pointers remain valid during this method, since no
     * vertices or faces are allocated or deleted during moving. */
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

    /* add all vertices of B to (this) mesh, store iterators to new vertices */
    std::pair<
            typename std::map<uint32_t, VertexPointerType>::iterator,
            bool
        > v_rpair;
    auto B_vit = B.V.begin();
    while (B_vit != B.V.end()) {
        /* get fresh vertex id */
        new_id      = this->V_idq.getId();

        /* insert pair (new_id, pointer B_vit->second) into this->V */
        v_rpair     = this->V.insert( { new_id, VertexPointerType(B_vit->second) } );
        inserted    = v_rpair.second;
        if (inserted) {
            /* the vertex pointer has been moved to (this) mesh: update mesh pointer Vertex::mesh
             * and iterator Vertex::m_vit */
            /* get pointer to vertex */
            v_newit     = v_rpair.first;
            v           = v_newit->second;
            v->mesh     = this;
            v->m_vit    = v_newit;

            /* erase B_vit from B.V */
            B_vit       = B.V.erase(B_vit); 
        }
        else {
            debugTabDec();
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::moveAppend(): insertion of vertex pointer into this->V with fresh id failed. internal logic error.");
        }
    }

    /* move all faces of B to (this) mesh in very much the same way */
    std::pair<
            typename std::map<uint32_t, FacePointerType >::iterator,
            bool
        > f_rpair;
    auto B_fit = B.F.begin();
    while (B_fit != B.F.end()) {
        new_id      = this->F_idq.getId();
        auto it = this->F.find(new_id);
        if ( it != this->F.end()) {
            debugTabDec();
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::moveAppend(): fresh id already taken in this->F. internal logic error.\n");
        }
        f_rpair     = this->F.insert( { new_id, FacePointerType(B_fit->second) } );
        inserted    = f_rpair.second;
        if (!inserted) {
            debugTabDec();
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::moveAppend(): insertion of face pointer into this->F with fresh id failed. internal logic error.");
        }
        else {
            f_newit     = f_rpair.first;
            f           = f_newit->second;
            f->mesh     = this;
            f->m_fit    = f_newit;
            B_fit       = B.F.erase(B_fit); 
        }
    }

    /* mesh octree needs update */
    this->octree_updated    = false;

    /* clear all information from B (B.V and B.F are empty, yet id queues etc are still set */
    if (!B.F.empty() || !B.V.empty()) {
        debugTabDec();
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::moveAppend(): internal vertex and face maps of appended mesh not empty at the end of process. internal logic error.");
    }
    B.clear();

    /* if update_vits != NULL, the list update_list of vertex pointers corresponding to the list of
     * given old vertex iterators (referring to B before the append) has been generated above: use
     * it to update the iterators in-place, i.e. an old iterator referring to a vertex of B is
     * replaced with an iterator referring to this moved vertex as a component of (this) mesh. */
    if (update_vits) {
        auto vit = update_vits->begin();
        auto lit = update_list.begin();
        while (vit != update_vits->end()) {
            *vit= (*lit)->iterator();
            ++vit, ++lit;
        }
    }
    debugTabDec();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::deleteConnectedComponent(vertex_iterator vstart_it)
{
    debugl(2, "Mesh::deleteConnectedComponent(). start vertex: %5d\n", vstart_it->id());
    debugTabInc();

    Vertex                 *v;
    std::queue<Vertex *>    Q;
    std::list<Vertex *>     v_nbs;
    std::list<Vertex *>     cc_vertices;

    const uint32_t          traversal_id = this->getFreshTraversalId();

    /* traverse cc of vstart_it, burning it down as we go.. */
    Q.push( &(*vstart_it) );

    debugl(3, "finding all vertices of connected component.\n");
    debugTabInc();
    while (!Q.empty()) {
        v = Q.front();
        Q.pop();

        debugl(4, "current vertex %5d\n", v->id() );

        /* get vertex star of v */
        v->getVertexStar(v_nbs);

        debugTabInc();
        /* iterator over all neighbours u and enqueue them if they haven't been seen yet */
        for (auto &u : v_nbs) {
            if (u->getTraversalState(traversal_id) == TRAV_UNSEEN) {
                debugl(4, "yet unseen neighbour %5d => enqueueing..\n", u->id() );
                Q.push(u);
                u->setTraversalState(traversal_id, TRAV_ENQUEUED);
            }
        }
        debugTabDec();

        /* v is done */
        debugl(4, "vertex %5d done.\n", v->id() );
        v->setTraversalState(traversal_id, TRAV_DONE);
        cc_vertices.push_back(v);
    }
    debugTabDec();

    debugl(3, "traversal of connected component completed. deleting..\n");

    debugTabInc();
    /* delete all vertices of cc, which in turn deletes all faces containing any such vertex. */
    for (auto &v : cc_vertices) {
        debugl(4, "deleting vertex %5d\n", v->id() );
        this->vertices.erase(v->iterator()) ;
        debugl(4, "vertex deleted.\n");
    }
    debugTabDec();

    this->octree_updated = false;

    debugTabDec();
    debugl(2, "Mesh::deleteConnectedComponent(). done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::getConnectedComponent(
    Mesh::vertex_iterator   vstart_it,    
    const uint32_t         &traversal_id,
    std::list<Vertex *>    *cc_vertices,
    std::list<Face *>      *cc_faces)
{
    debugl(2, "Mesh::getConnectedComponent(). start vertex: %5d\n", vstart_it->id());
    debugTabInc();

    if (!vstart_it.checkContainer(*this) || vstart_it == this->vertices.end()) {
        throw("Mesh::getConnectedComponent(): invalid start vertex iterator.");
    }

    /* only process if one of the given list pointers isn't NULL and the start vertex's traversal
     * state is set to UNSEEN. */
    if ( (cc_vertices || cc_faces) && vstart_it->getTraversalState(traversal_id) == TRAV_UNSEEN) {
        std::queue<Vertex *>    Q;

        Vertex                 *v;
        std::list<Vertex *>     v_vstar;
        std::list<Face *>       v_fstar;

        /* traverse cc of vstart_it, burning it down as we go.. */
        Q.push( &(*vstart_it) );

        debugl(3, "traversing connected component of start vertex %5d.\n", Q.front()->id());
        debugTabInc();
        while (!Q.empty()) {
            v = Q.front();
            Q.pop();

            debugl(4, "current vertex %5d\n", v->id() );

            /* add v and all its incident faces to the result lists if desired by the caller */
            if (cc_vertices) {
                cc_vertices->push_back(v);
            }
            if (cc_faces) {
                /* get face star of v */
                v->getFaceStar(v_fstar);

                /* a face incident incident to v is appended to cc_faces list iff this face hasn't
                 * been seen yet. in this case, set traversal state to DONE to prevent multiple
                 * insertions of the same face. */
                for (auto v_if : v_fstar) {
                    if (v_if->getTraversalState(traversal_id) == TRAV_UNSEEN) {
                        cc_faces->push_back(v_if);
                        v_if->setTraversalState(traversal_id, TRAV_DONE);
                    }
                }
            }

            /* get vertex star of v */
            v->getVertexStar(v_vstar);

            debugTabInc();
            /* iterate over all vertex neighbours u and enqueue them if they haven't been seen yet. */
            for (auto &u : v_vstar) {
                if (u->getTraversalState(traversal_id) == TRAV_UNSEEN) {
                    debugl(4, "yet unseen neighbour %5d => enqueueing..\n", u->id() );
                    Q.push(u);
                    u->setTraversalState(traversal_id, TRAV_ENQUEUED);
                }
            }
            debugTabDec();

            /* v is done. */
            debugl(4, "vertex %5d done.\n", v->id() );
            v->setTraversalState(traversal_id, TRAV_DONE);
        }
        debugTabDec();
        debugl(3, "traversal of connected component completed.\n");
    }

    debugTabDec();
    debugl(2, "Mesh::getConnectedComponent(). done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::deleteBorderCCsAndIsolatedVertices()
{
    debugl(2, "Mesh::deleteBorderCCsAndIsolatedVertices()\n");
    debugTabInc();

    bool done = false;

restart_loop:
    while (!done) {
        debugl(3, "(re)entering loop..\n");
        done = true;

        Vertex                 *s, *v;
        std::list<Vertex *>     v_vstar;
        std::list<Face *>       v_fstar;

        std::queue<Vertex *>    Q;

        uint8_t                 u_state;
        uint32_t                ccs, cc_size;

        const uint32_t          traversal_id = this->getFreshTraversalId();

        ccs         = 0;
        debugl(3, "deleting all isolated vertices..\n");

        debugTabInc();

        /* assemble list of isolated vertices */
        std::list<Vertex *> isolated_vertices;
        for (auto &v : this->vertices) {
            v.getVertexStar(v_vstar);

            if (v_vstar.empty()) {
                v.getFaceStar(v_fstar);
#ifdef __DEBUG__
                debugl(3, "found isolated vertex: %d. vertex_star.size(): %d, face_star.size(): %d: incident faces..:\n", v.id(), v_vstar.size(), v_fstar.size());
                debugTabInc();
                for (auto &f : v_fstar) {
                    debugl(3, "%5d : (%5d, %5d, %5d)\n", f->id(), f->vertices[0]->id(), f->vertices[1]->id(), f->vertices[2]->id());
                }
                debugTabDec();
#endif
                isolated_vertices.push_back(&v);
            }
        }

        /* delete all isolated vertices */
        for (auto &v : isolated_vertices) {
            debugl(3, "deleting isolated vertex: %d\n", v->id());
            this->vertices.erase(v->iterator());
        }

        debugTabDec();

        /* traverse all remaining connected components */
        debugl(3, "analysing all connected components.\n");

        debugTabInc();
        for (Vertex &vref : this->vertices) {
            if (vref.getTraversalState(traversal_id) == TRAV_UNSEEN) {
                /* new connected component with root s */
                s       = &vref;
                cc_size = 0;
                ccs++;

                debugl(4, "new connected component with root vertex %5d\n", s->id() );

                /* traverse cc */
                Q.push(s);
                while (!Q.empty()) {
                    v = Q.front();
                    Q.pop();
                    cc_size++;

                    /* get neighbours of v */
                    v_vstar.clear();
                    v->getVertexStar(v_vstar);

                    /* inspect all neighbours. note: since all isolated vertices have been deleted
                     * above, this loop is always non-trivial */
                    for (auto &u : v_vstar) {
                        /* get traversal state of current neighbour u */
                        u_state = u->getTraversalState(traversal_id);

                        /* we handle the edge in in the direction {v -> u}. if u is not done (but
                         * can be either enqueued or unseen), otherwise u has already handled the
                         * edge in direction {u -> v} */
                        if (u_state != TRAV_DONE) {
                            Face* vu_incident_faces[2];

                            /* get faces containing edge {v, u} */
                            size_t nIF = 2;
                            this->getFacesIncidentToEdge(v->iterator(), u->iterator(), vu_incident_faces, nIF);

                            /* in case of border edge {v, u}, delete connected component of v and
                             * restart the loop */
                            if (nIF == 1) {
                                this->deleteConnectedComponent( v->iterator() );
                                done = false;
                                goto restart_loop;
                            }
                            /* isolated edge, this should be impossible */
                            else if (nIF == 0) {
                                debugl(1, "Mesh::deleteBorderCCsAndIsolatedVertices(): isolated edge (%d, %d) found.\n", v->id(), u->id());
                                throw MeshEx(MESH_LOGIC_ERROR, "Mesh::deleteBorderCCsAndIsolatedVertices(): isolated edge found. this should never happen due to internal invariants => internal logic error.");
                            }
                            /* two or more faces, case irrelevant here */
                        }
                        
                        /* enqueue u if it hasn't been seen yet */
                        if (u_state == TRAV_UNSEEN) {
                            Q.push(u);
                            u->setTraversalState(traversal_id, TRAV_ENQUEUED);
                        }
                    }

                    /* v is done */
                    v->setTraversalState(traversal_id, TRAV_DONE);
                }

                debugl(4, "finished traversal of connected component with root vertex %5d. size: %5d\n", s, cc_size);
            }
        }
        debugTabDec();
    }

    /* octree no longer up to date */
    this->octree_updated = false;

    debugTabDec();
    debugl(2, "Mesh::deleteBorderCCsAndIsolatedVertices(): done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::updateOctree()
{
    using namespace Aux::Timing;
    using Aux::Numbers::inf;
    using Aux::VecMat::minVec3;
    using Aux::VecMat::maxVec3;
    using Aux::VecMat::fabsVec3;

    debugl(1, "Mesh::updateOctree():..\n");
    debugTabInc();

    if (!this->octree_updated) {
        R                       rootcube_len;
        Mesh_OctreeInfo         tree_info;
        Mesh_OctreeNodeInfo     root_info;
        std::list<Face *>       root_face_list;
        std::list<Vertex *>     root_vertex_list;

        tick(15);

        /* compute AABB fore entire mesh */
        this->bb    = BoundingBox<R>();

        for (auto &f : this->faces) {
            this->bb.update(f.getBoundingBox());
        }

        this->bb.extend(0.025, Vec3<R>(1E-3, 1E-3, 1E-3));

        /* compute root cube length with mesh bounding box */
        Vec3<R> aabb_min = this->bb.min();
        Vec3<R> aabb_max = this->bb.max();

        debugl(1, "bounding box for entire mesh: (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f)\n", 
                aabb_min[0], aabb_min[1], aabb_min[2],
                aabb_max[0], aabb_max[1], aabb_max[2]);

        /* fresh Octree .. */
        if (O) {
            delete O;
        }

        rootcube_len = std::max(aabb_max[0] - aabb_min[0], aabb_max[1] - aabb_min[1]);
        rootcube_len = std::max(rootcube_len, aabb_max[2] - aabb_min[2]);

        root_info.cube_min      = aabb_min;
        root_info.cube_max[0]   = aabb_min[0] + rootcube_len;
        root_info.cube_max[1]   = aabb_min[1] + rootcube_len;
        root_info.cube_max[2]   = aabb_min[2] + rootcube_len;

        debugl(1, "root cube:: (%5.4f, %5.4f, %5.4f) len: (%5.4f, %5.4f, %5.4f)\n", 
                root_info.cube_min[0], root_info.cube_min[1], root_info.cube_min[2],
                root_info.cube_max[0], root_info.cube_max[1], root_info.cube_max[2]);


        /* init face list containing pointers to all mesh faces, which is passed to the partitioning
         * function */
        for (auto &f : this->faces) {
            root_face_list.push_back( &f );
        }

        for (auto &v : this->vertices) {
            root_vertex_list.push_back( &v );
        }

        /* topoloy information is updated -> bounding box is up to date as well */
        O = new Octree<Mesh_OctreeInfo, Mesh_OctreeNodeInfo, R>(
                tree_info,
                root_info,
                aabb_min,
                rootcube_len);

        /* call recursive partitioning algorithm */
        debugl(2, "calling recursive partitioning / construction algorithm.\n", time);
        Mesh::partitionOctree(*O, (*this), &(O->root), root_face_list, root_vertex_list, 0, 64, 9);

        debugl(2, "recursive octree construction done. time: %10.5f\n", tack(15));

        /* octree has been updated */
        this->octree_updated = true;
    }

    debugTabDec();
    debugl(1, "Mesh::updateOctree(). done.\n");
}

/* locate vertices in or very nearly in bounding box (octree is not fully constructed due to massive
 * overhead, to the vicinity might be returned as well, but not much..) */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::findVertices(
    BoundingBox<R> const   &search_box,
    std::list<Vertex *>    &vertex_list)
{
    using namespace Aux::Geometry::IntersectionTestResults;

    Vec3<R> aabb_min = search_box.min(), aabb_max = search_box.max();

    debugl(2, "Mesh::findVertices(): input bb (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f)(\n",
            aabb_min[0], aabb_min[1], aabb_min[2],
            aabb_max[0], aabb_max[1], aabb_max[2]);
    debugTabInc();

    if (!octree_updated) {
        this->updateOctree();
    }

    /* recursive octree location traversal, set face list argument to NULL since only faces need to
     * be found */
    Mesh::findOctreeRecursive(O, &(O->root), aabb_min, aabb_max, &vertex_list, NULL);

    /* we want unique lists */
    vertex_list.sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});
    vertex_list.unique([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() == y->id());});

    /* delete vertices whose bounding box does not intersect the given search box */
    for (auto it = vertex_list.begin(); it != vertex_list.end(); ) {
        if (search_box && (*it)->getBoundingBox()) {
            ++it;
        }
        else {
            it = vertex_list.erase(it);
        }
    }

    debugTabDec();
    debugl(2, "Mesh::findVertices(): done. %d vertices found.\n", vertex_list.size() );
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::findFaces(
    BoundingBox<R> const   &search_box,
    std::list<Face *>      &face_list)
{
    using namespace Aux::Geometry::IntersectionTestResults;

    Vec3<R> aabb_min = search_box.min(), aabb_max = search_box.max();

    debugl(2, "Mesh::findFaces(): input bb (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f)(\n",
            aabb_min[0], aabb_min[1], aabb_min[2],
            aabb_max[0], aabb_max[1], aabb_max[2]);
    debugTabInc();

    if (aabb_min >= aabb_max) {
        debugl(1, "Mesh::findFaces(): invalid bounding box: (%5.4f, %5.4f, %20.10e) - (%5.4f, %5.4f, %20.10e)\n",
                aabb_min[0], aabb_min[1], aabb_min[2],
                aabb_max[0], aabb_max[1], aabb_max[2]);

        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::findFaces(): invalid input bounding boxes: minimum value >= maximum value for some component.");
    }

    /* clear face_list */
    face_list.clear();

    if (!octree_updated) {
        this->updateOctree();
    }

    /* recursive octree location traversal, set vertex list argument to NULL since only faces need
     * to be found */
    Mesh::findOctreeRecursive(O, &(O->root), aabb_min, aabb_max, NULL, &face_list);

    /* we want unique lists */
    face_list.sort([] (const Face* x, const Face* y) -> bool {return (x->id() < y->id());});
    face_list.unique([] (const Face* x, const Face* y) -> bool {return (x->id() == y->id());});

    /* delete faces whose bounding box does not intersect the given search box */
    for (auto it = face_list.begin(); it != face_list.end(); ) {
        if (search_box && (*it)->getBoundingBox()) {
            ++it;
        }
        else {
            it = face_list.erase(it);
        }
    }

    debugTabDec();
    debugl(2, "Mesh::findFaces(): done. %d faces found.\n", face_list.size() );
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::find(
    BoundingBox<R> const   &search_box,
    std::list<Vertex *>    *vertex_list,
    std::list<Face *>      *face_list)
{
    using namespace Aux::Geometry::IntersectionTestResults;

    Vec3<R> aabb_min = search_box.min(), aabb_max = search_box.max();

    debugl(2, "Mesh::findFaces(): input bb (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f)(\n",
            aabb_min[0], aabb_min[1], aabb_min[2],
            aabb_max[0], aabb_max[1], aabb_max[2]);
    debugTabInc();

    if (aabb_min >= aabb_max) {
        debugl(1, "Mesh::findFaces(): invalid bounding box: (%5.4f, %5.4f, %20.10e) - (%5.4f, %5.4f, %20.10e)\n",
                aabb_min[0], aabb_min[1], aabb_min[2],
                aabb_max[0], aabb_max[1], aabb_max[2]);

        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::findFaces(): invalid input bounding boxes: minimum value >= maximum value for some component.");
    }

    /* clear input lists */
    if (vertex_list)    vertex_list->clear();
    if (face_list)      face_list->clear();

    if (!octree_updated) {
        this->updateOctree();
    }

    /* if not both vertex_list and face_list are NULL, call recursive octree location traversal */
    if (vertex_list || face_list) {
        Mesh::findOctreeRecursive(O, &(O->root), aabb_min, aabb_max, vertex_list, face_list);

        /* we want unique lists */
        vertex_list->sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});
        vertex_list->unique([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() == y->id());});

        /* delete vertices whose bounding box does not intersect the given search box */
        Vec3<R> v_bb_min, v_bb_max;
        for (auto it = vertex_list->begin(); it != vertex_list->end(); ) {
            if (search_box && (*it)->getBoundingBox()) {
                ++it;
            }
            else {
                it = vertex_list->erase(it);
            }
        }

        face_list->sort([] (const Face* x, const Face* y) -> bool {return (x->id() < y->id());});
        face_list->unique([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() == y->id());});

        /* delete faces whose bounding box does not intersect the given search box */
        for (auto it = face_list->begin(); it != face_list->end(); ) {
            if (search_box && (*it)->getBoundingBox()) {
                ++it;
            }
            else {
                it = face_list->erase(it);
            }
        }

    }

    debugTabDec();
    debugl(2, "Mesh::findFaces(): done. %d faces found.\n", face_list->size() );
}

/* static recursive partitioning function for Octree construction */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::partitionOctree(
    Octree<Mesh_OctreeInfo, Mesh_OctreeNodeInfo, R>    &O,
    Mesh const                                         &m,
    OctreeNode<Mesh_OctreeNodeInfo>                    *n,
    std::list<Face *>                                  &rec_facelist,
    std::list<Vertex *>                                &rec_vertexlist,
    uint32_t                                            rec_depth,
    uint32_t                                            max_elements,
    uint32_t                                            max_rec_depth)
{
    using namespace Aux::Geometry::IntersectionTestResults;

    debugl(5, "Mesh::partitionOctree() (recursive)\n");
    debugTabInc();

    debugl(5, "depth %4d. facelist.size(): %6ld\n", rec_depth, rec_facelist.size() );
    /* leaf reached. save face_list */
    if (rec_depth >= max_rec_depth) {
        n->data.face_list   = rec_facelist;
        n->data.vertex_list = rec_vertexlist;
        debugl(5, "max depth reached => creating leaf.\n");
    }
    /* partition face list with current cube, recursive call */
    else {
        /* if rec_facelist contains more than max_faces, process it */
        if (rec_facelist.size() + rec_vertexlist.size() > max_elements) {
            debugl(5, "max depth not yet reached => partitioning %8ld faces and %8ld vertices among children..\n", rec_facelist.size(), rec_vertexlist.size() );
            uint32_t                            i;
            Vec3<R>                             nc_min, nc_max, nc_m;
            Vec3<R>                             xdisp, ydisp, zdisp;

            std::list<Face *>                  *children_facelists[8];
            std::list<Vertex *>                *children_vertexlists[8];
            OctreeNode<Mesh_OctreeNodeInfo>    *children[8];

            /* split node n, which is currently still a leaf */
            O.splitLeaf(n);

            /* get child pointers and allocate face lists for children */
            for (i = 0; i < 8; i++) {
                children[i]             = n->getChild(i);
                children_facelists[i]   = new std::list<Face *>;
                children_vertexlists[i] = new std::list<Vertex *>;
            }

            /* compute midpoint of n's cube and init children's cube min / max vectors.
             * incidentally, this does not need to be a cube here.. but who cares. the OtId thing
             * has been designed for DMC and is not used here. space.. well.. o_O */
            nc_min                      = n->data.cube_min;
            nc_max                      = n->data.cube_max;
            nc_m                        = (nc_min + nc_max) * 0.5;

            xdisp                       = Vec3<R>( (nc_max[0] - nc_min[0]) / 2.0, 0.0, 0.0);
            ydisp                       = Vec3<R>( 0.0, (nc_max[1] - nc_min[1]) / 2.0, 0.0);
            zdisp                       = Vec3<R>( 0.0, 0.0, (nc_max[2] - nc_min[2]) / 2.0);

            /* child 0: n_min and midpoint. offsetting from there with displacement vectors */
            children[0]->data.cube_min  = nc_min;
            children[0]->data.cube_max  = nc_m;

            children[1]->data.cube_min  = nc_min + xdisp;
            children[1]->data.cube_max  = nc_m + xdisp;

            children[2]->data.cube_min  = nc_min + xdisp + ydisp;
            children[2]->data.cube_max  = nc_m + xdisp + ydisp;

            children[3]->data.cube_min  = nc_min + ydisp;
            children[3]->data.cube_max  = nc_m + ydisp;

            children[4]->data.cube_min  = nc_min + zdisp;
            children[4]->data.cube_max  = nc_m + zdisp;

            children[5]->data.cube_min  = nc_min + xdisp + zdisp;
            children[5]->data.cube_max  = nc_m + xdisp + zdisp;

            /* special case: child 6 has min corner m and max corner nc_max */
            children[6]->data.cube_min  = nc_m;
            children[6]->data.cube_max  = nc_max;

            children[7]->data.cube_min  = nc_min + ydisp + zdisp;
            children[7]->data.cube_max  = nc_m + ydisp + zdisp;

            /* iterate over all faces in rec_facelist and partition them into the child lists.
             * recursive call */
            Vec3<R> f_aabb_min, f_aabb_max;
            for (auto &f : rec_facelist) {
                auto f_bb   = f->getBoundingBox();
                f_aabb_min  = f_bb.min();
                f_aabb_max  = f_bb.max();

                /* check for intersections with all eight child cubes. */
                for (i = 0; i < 8; i++) {
                    /* mesh face bb intersecting child cube -> push intro children list */
                    if (Aux::Geometry::simpleIntersect2AABB(
                            f_aabb_min, 
                            f_aabb_max,
                            children[i]->data.cube_min,
                            children[i]->data.cube_max) == INTERSECTION)
                    {
                        children_facelists[i]->push_back(f);
                    }
                }
            }

            /* iterate over all vertices in the list, get the vector, construct a tiny AABB and
             * distribute among children. AABB length: 2.0*1E-8 */
            Vec3<R> vpos;
            Vec3<R> v_aabb_min, v_aabb_max;
            for (auto &v : rec_vertexlist) {
                vpos        = v->pos();
                auto v_bb   = v->getBoundingBox();
                v_aabb_min  = v_bb.min();
                v_aabb_max  = v_bb.max();

                /* check for intersections with all eight child cubes. */
                for (i = 0; i < 8; i++) {
                    /* mesh face bb intersecting child cube -> push intro children list */
                    if (Aux::Geometry::simpleIntersect2AABB(
                            v_aabb_min,
                            v_aabb_max,
                            children[i]->data.cube_min,
                            children[i]->data.cube_max) == INTERSECTION)
                    {
                        children_vertexlists[i]->push_back(v);
                    }
                }
            }

            /* clear recursion lists, since they are not needed anymore */
            rec_facelist.clear();
            rec_vertexlist.clear();

            /* issue recursive calls */
            for (i = 0; i < 8; i++) {
                debugl(5, "recursive call for child %02d..\n", i);

                debugTabInc();
                /* here we go.. */
                partitionOctree(O, m, children[i], *(children_facelists[i]), *(children_vertexlists[i]), rec_depth + 1, max_elements, max_rec_depth);
                debugTabDec();

                debugl(5, "finished recursive call for child %02d..\n", i);

                /* delete child lists allocated above */
                delete children_facelists[i];
                delete children_vertexlists[i];
            }
        }
        /* otherwise, facelist contains <= max_faces faces and will be turned into a leaf */
        else {
            debugl(5, "numer of vertices and faces small enough => creating leaf. \n");
            n->data.face_list   = rec_facelist;
            n->data.vertex_list = rec_vertexlist;
        }
    }

    debugTabDec();
    debugl(5, "Mesh::partitionOctree() (recursive): done.\n");
}


/* static recursive octree location traversal: given a bounding box, find all leafs in the octree
 * that intersect the bounding box and return the union of their face lists */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::findOctreeRecursive(
    Octree<Mesh_OctreeInfo, Mesh_OctreeNodeInfo, R>    *O,
    OctreeNode<Mesh_OctreeNodeInfo>                    *n,
    Vec3<R>                                             aabb_min,
    Vec3<R>                                             aabb_max,
    std::list<Vertex *>                                *vertex_list,
    std::list<Face *>                                  *face_list)
{
    using namespace Aux::Geometry::IntersectionTestResults;

    debugl(5, "Mesh::findOctreeRecursive() (recursive)\n");
    debugTabInc();

    /* n is a leaf, append n's facelist to argument face_list */
    if (n->isLeaf()) {
        debugl(5, "given bb (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f), n->bb: (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f): leaf => appending %5d faces and / or %d vertices.\n",
                aabb_min[0], aabb_min[1], aabb_min[2], 
                aabb_max[0], aabb_max[1], aabb_max[2], 
                n->data.cube_min[0], n->data.cube_min[1], n->data.cube_min[2],
                n->data.cube_max[0], n->data.cube_max[1], n->data.cube_max[2],
                n->data.vertex_list.size(), n->data.face_list.size());

        if (vertex_list)    vertex_list->insert( vertex_list->end(), n->data.vertex_list.begin(), n->data.vertex_list.end() );
        if (face_list)      face_list->insert( face_list->end(), n->data.face_list.begin(), n->data.face_list.end() );
    }
    /* otherwise, get all children of n, and call recursively if given bounding box and cube
     * intersect */
    else {
        debugl(5, "given bb (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f), n->bb: (%5.4f, %5.4f, %5.4f) - (%5.4f, %5.4f, %5.4f): no leaf => recursive calls..\n",
                aabb_min[0], aabb_min[1], aabb_min[2], 
                aabb_max[0], aabb_max[1], aabb_max[2], 
                n->data.cube_min[0], n->data.cube_min[1], n->data.cube_min[2],
                n->data.cube_max[0], n->data.cube_max[1], n->data.cube_max[2]);

        OctreeNode<Mesh_OctreeNodeInfo>    *child;
        uint32_t                            i;
        for (i = 0; i < 8; i++) {
            child = n->getChild(i);
            if (child && Aux::Geometry::simpleIntersect2AABB(
                        child->data.cube_min,
                        child->data.cube_max,
                        aabb_min,
                        aabb_max) == INTERSECTION)
            {
                debugl(5, "child %2d's bb intersects => calling.\n", i);

                debugTabInc();
                Mesh::findOctreeRecursive(O, child, aabb_min, aabb_max, vertex_list, face_list);
                debugTabDec();

                debugl(5, "child %2d's call done.\n", i);
            }
        }
    }

    debugTabDec();
    debugl(5, "Mesh::findOctreeRecursive() (recursive): done.\n");
}

/* DEBUG: for debugging, output cubes of all leafs */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::getOctreeLeafCubes(
    Octree<Mesh_OctreeInfo, Mesh_OctreeNodeInfo, R>    *O,
    OctreeNode<Mesh_OctreeNodeInfo>                    *n,
    std::list<std::pair<Vec3<R>, Vec3<R>>>             *cube_list)
{
    debugl(5, "Mesh::getOctreeLeafCubes() (recursive)\n");
    debugTabInc();

    /* n is a leaf, append n's facelist to argument face_list */
    if (n->isLeaf()) {
        debugl(5, "get_leaf_cubes: leaf => appending cube to list.\n");
        cube_list->push_back( {n->data.cube_min, n->data.cube_max} );
    }
    /* otherwise, get all children of n, and call recursively if given bounding box and cube
     * intersect */
    else {
        OctreeNode<Mesh_OctreeNodeInfo>    *child;
        uint32_t                            i;
        debugl(5, "get_leaf_cubes: no leaf => recursive calls for children\n");
        for (i = 0; i < 8; i++) {
            child = n->getChild(i);
            debugl(5, "child %2d's bb intersects => calling.\n", i);

            debugTabInc();
            Mesh::getOctreeLeafCubes(O, child, cube_list);
            debugTabDec();

            debugl(5, "child %2d's call done.\n", i);
        }
    }

    debugTabDec();
    debugl(5, "Mesh::getOctreeLeafCubes() (recursive): done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::selectNonManifoldVertices(std::list<Vertex *> &vlist) const
{
    for (auto &v : this->vertices) {
        if (!v.isManifold()) {
            vlist.push_back(&v);
        }
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::selectIsolatedVertices(std::list<Vertex *> &vlist) const
{
    for (auto &v : this->vertices) {
        if (!v.isIsolated()) {
            vlist.push_back(&v);
        }
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::invertVertexSelection(std::list<Vertex *> &vlist) const
{
    class id_sort_cmp {
        public:
            bool
            operator()(
                Vertex * const &x,
                Vertex * const &y)
            {
                return (x->id() < y->id());
            }
    };

    class cmp_functor {
        public:
            bool
            operator()(
                Vertex * const                                 &x,
                std::pair<uint32_t, VertexPointerType> const   &ypair)
            {
                return (x->id() < ypair.first);
            }

            bool operator()(
                std::pair<uint32_t, VertexPointerType> const   &xpair,
                Vertex * const                                 &y)
            {
                return (xpair.first < y->id());
            }
    };

    std::list<std::pair<uint32_t, Vertex *>> result;

    vlist.sort(id_sort_cmp());
    std::set_difference(
        this->V.begin(), this->V.end(),
        vlist.begin(), vlist.end(),
        std::back_inserter(result), cmp_functor() );

    vlist.clear();
    for (auto &vp : result) {
        vlist.push_back(vp.second);
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::invertFaceSelection(std::list<Face *> &flist) const
{
    class id_sort_cmp {
        public:
            bool
            operator()(
                Face * const &x,
                Face * const &y)
            {
                return (x->id() < y->id());
            }
    };

    class cmp_functor {
        public:
            bool operator()(
                Face * const                               &x,
                std::pair<uint32_t, FacePointerType> const &ypair)
            {
                return (x->id() < ypair.first);
            }

            bool operator()(
                std::pair<uint32_t, FacePointerType> const &xpair,
                Face * const                               &y)
            {
                return (xpair.first < y->id());
            }
    };

    std::list<std::pair<uint32_t, Face *> > result;

    flist.sort(id_sort_cmp());
    std::set_difference(
        this->F.begin(), this->F.end(),
        flist.begin(), flist.end(),
        std::back_inserter(result), cmp_functor() );

    flist.clear();
    for (auto &fp : result) {
        flist.push_back(fp.second);
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::checkEdge(
    const std::string&              fn,
    const vertex_const_iterator    &u_it,
    const vertex_const_iterator    &v_it) const
{
    if (!u_it.checkContainer( *this ) || !u_it.sameContainer(v_it) || u_it == v_it) {
        throw MeshEx(MESH_LOGIC_ERROR, fn + "given edge (u, v) invalid: u and v not from (this) mesh or identical.");
    }
    else if (!u_it->gotNeighbour(v_it)) {
        throw MeshEx(MESH_LOGIC_ERROR, fn + "given edge (u, v) invalid: both vertices from (this) mesh, yet edge does not exist.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::checkEdge(
    const char                     *fn,
    const vertex_const_iterator    &u_it,
    const vertex_const_iterator    &v_it) const
{
    if (!u_it.checkContainer( *this ) || !u_it.sameContainer(v_it) || u_it == v_it) {
        throw MeshEx(MESH_LOGIC_ERROR, std::string(fn) + "given edge (u, v) invalid: u and v not from (this) mesh or identical.");
    }
    else if (!u_it->gotNeighbour(v_it)) {
        throw MeshEx(MESH_LOGIC_ERROR, std::string(fn) + "given edge (u, v) invalid: both vertices from (this) mesh, yet edge does not exist.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::getFacesIncidentToEdge(
    const vertex_const_iterator    &u_it,
    const vertex_const_iterator    &v_it,
    Face**                          incident_faces,
    size_t&                         sizeInOut) const
{
    this->checkEdge("Mesh::getFacesIncidentToEdge()", u_it, v_it);

    // this is extremely slow due to recurring call of new operator for creation of lists
#if 0
    std::list<Face *> ufaces, vfaces;

    u_it->getFaceStar(ufaces);
    v_it->getFaceStar(vfaces);

    /* clear() result, sort and set intersect into result */
    incident_faces.clear();
    ufaces.sort(Mesh::Face::ptr_less);
    vfaces.sort(Mesh::Face::ptr_less);

    /* use back_inserter to compile result into incident_faces */
    std::set_intersection(ufaces.begin(), ufaces.end(), vfaces.begin(), vfaces.end(), std::back_inserter(incident_faces) );
#endif

    // we do it brute force without sorting and are still faster
    size_t sz = 0;
    const std::list<Face*>& uFaces = u_it->getFaceStar();
    const std::list<Face*>& vFaces = v_it->getFaceStar();
    typename std::list<Face*>::const_iterator itU = uFaces.begin();
    typename std::list<Face*>::const_iterator itV;
    typename std::list<Face*>::const_iterator itUend = uFaces.end();
    typename std::list<Face*>::const_iterator itVend = vFaces.end();
    for (; itU != itUend; ++itU)
    {
        for (itV = vFaces.begin(); itV != itVend; ++itV)
        {
            if ((*itU)->id() == (*itV)->id())
            {
                if (++sz > sizeInOut)
                    throw MeshEx(MESH_LOGIC_ERROR, "Found more incident faces to edge than expected.");
                incident_faces[sz-1] = *itU;
                break;
            }
        }
    }

    if (!sz)
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getFacesIncidentToEdge(): input edge invalid or non-existent.");
    sizeInOut = sz;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::getFacesIncidentToEdge(
    uint32_t                u,
    uint32_t                v,
    Face**                  incident_faces,
    size_t&                 sizeInOut) const
{
    /* get face with ids u and v and both lists of incident faces */
    vertex_const_iterator u_it, v_it;
    u_it = this->vertices.find(u);
    v_it = this->vertices.find(v);

    if (u_it == this->vertices.end() || v_it == this->vertices.end()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getFacesIncidentToEdge(): at least one of the given vertex ids is out of range.");
    }
    else this->checkEdge("Mesh::getFacesIncidentToEdge()", u_it, v_it);

    /* call version with iterators as parameters */
    this->getFacesIncidentToEdge(u_it, v_it, incident_faces, sizeInOut);
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::getFacesIncidentToManifoldEdge(
    const vertex_const_iterator    &u_it,
    const vertex_const_iterator    &v_it,
    face_iterator                  &fst_face,
    face_iterator                  &snd_face) const
{
    this->checkEdge("Mesh::getFacesIncidentToManifoldEdge()", u_it, v_it);

    Face* uv_faces[2];
    size_t nIF = 2;
    this->getFacesIncidentToEdge(u_it, v_it, uv_faces, nIF);

    if (nIF != 2) {
        throw MeshEx(MESH_LOGIC_ERROR, "getFacesIncidentToManifoldEdge(): specified edge not a manifold edge.");
    }
    else {
        fst_face = uv_faces[0]->iterator();
        snd_face = uv_faces[1]->iterator();
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::getOtherFaceIncidentToManifoldEdge(
    const vertex_const_iterator    &u_it,
    const vertex_const_iterator    &v_it,
    const face_const_iterator      &known_face_it) const
{
    this->checkEdge("Mesh::getOtherFaceIncidentToManifoldEdge()", u_it, v_it);

    Face* nbfaces[2];
    size_t nIF = 2;
    this->getFacesIncidentToEdge(u_it, v_it, nbfaces, nIF);

    /* if there are not exactly two faces, throw exception */
    if (nIF != 2) {
        if (nIF == 0) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getFacesIncidentToManifoldEdge(): edge has no incident faces.");
        }
        else {
            throw MeshEx(MESH_NONMANIFOLD, "Mesh::getFacesIncidentToManifoldEdge(): edge hat > 2 incident faces: non-manifold edge.");
        }
    }
    else {
        if ( known_face_it->id() == nbfaces[0]->id() ) {
            return nbfaces[1]->iterator();
        }
        else if (known_face_it->id() == nbfaces[1]->id() ) {
            return nbfaces[0]->iterator();
        }
        else {
            throw MeshEx(MESH_NOT_FOUND, "Mesh::getOtherFaceIncidentToManifoldEdge(): given known face id is not incident to specified edge.");
        }
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::getTriCommonEdge(
    const face_const_iterator      &A_it,
    const face_const_iterator      &B_it,
    vertex_iterator                &u_it,
    vertex_iterator                &v_it) const
{
    /* get vertex ids */
    Vertex *A_v0, *A_v1, *A_v2;
    Vertex *B_v0, *B_v1, *B_v2;

    /* check iterators */
    if (!A_it.checkContainer(*this) || !A_it.sameContainer(B_it)) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getTriCommonEdge(): input face iterators refer to different containers.");
    }
    else {
        A_it->checkTri("Mesh::getTriCommonEdge()");
        B_it->checkTri("Mesh::getTriCommonEdge()");
    }

    A_it->getTriVertices(A_v0, A_v1, A_v2);
    B_it->getTriVertices(B_v0, B_v1, B_v2);

    std::vector<Vertex *> A_vertices = { A_v0, A_v1, A_v2 };
    std::vector<Vertex *> B_vertices = { B_v0, B_v1, B_v2 };
    std::vector<Vertex *> AB_shared_vertices;

    auto vrtCmp = [] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());};
    std::sort(A_vertices.begin(), A_vertices.end(), vrtCmp);
    std::sort(B_vertices.begin(), B_vertices.end(), vrtCmp);
    std::set_intersection(
            A_vertices.begin(),
            A_vertices.end(),
            B_vertices.begin(),
            B_vertices.end(),
            std::inserter(AB_shared_vertices, AB_shared_vertices.begin()),
            vrtCmp
        );

    if (AB_shared_vertices.size() == 2) {
        /* write return variables */
        u_it = AB_shared_vertices.front()->iterator();
        v_it = AB_shared_vertices.back()->iterator();
        return true;
    }
    else if (AB_shared_vertices.size() == 0 || AB_shared_vertices.size() == 1) {
        return false;
    }
    else {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getTriCommonEdge(): two triangles sharing more than two vertices found. this must never happen.");
    }
}


template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::getTriCommonEdge(
    uint32_t    A_id,
    uint32_t    B_id,
    uint32_t   &u,
    uint32_t   &v) const
{
    face_const_iterator     A_it, B_it;
    vertex_iterator         u_it, v_it;
    bool                    ret;

    A_it = this->faces.find(A_id);
    B_it = this->faces.find(B_id);

    if (A_it == this->faces.end() || B_it == this->faces.end()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getTriCommonEdge(): at least one input id out of range.");
    }
    ret = this->getTriCommonEdge(A_it, B_it, u_it, v_it);

    if (ret) {
        u = u_it->id();
        v = v_it->id();
    }
    return ret;
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::getEdgeBoundingBox(
    uint32_t    u_id,
    uint32_t    v_id,
    Vec3<R>    &aabb_min,
    Vec3<R>    &aabb_max) const
{
    using namespace Aux::VecMat;

    Mesh::vertex_const_iterator uit, vit;
    Vec3<R>                     u, v, offset;
    std::list<uint32_t>         u_nbs;

    debugl(3, "Mesh::getEdgeBoundingBox(): ..\n");

    uit = this->vertices.find(u_id);
    vit = this->vertices.find(v_id);

    uit->getFaceStarIndices(u_nbs);

    if (!Aux::Alg::listContains<uint32_t>(u_nbs, v_id)) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getEdgeBoundingBox(): neither (u,v) nor (v,u) for the given vertices u, v is not an edge of the mesh.");
    }

    minVec3(aabb_min, u, v);
    maxVec3(aabb_max, u, v);

    /* enlarge bounding box slightly to be save */
    /* if the difference in one coord is exactly zero, then a relative offset won't help.  use a
     * minimum absolute value of 1E-3 to prevent this. the offset vector below always has positive
     * components.. */
    maxVec3<R>(offset, Vec3<R>(1E-3, 1E-3, 1E-3), fabsVec3<R>(aabb_max - aabb_min)*0.025);
    aabb_min   -= offset;
    aabb_max   += offset;

    debugl(3, "Mesh::getEdgeBoundingBox(): done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::getSharedTri(
    const vertex_const_iterator    &u1_it,
    const vertex_const_iterator    &v1_it,
    const vertex_const_iterator    &u2_it,
    const vertex_const_iterator    &v2_it) const
{
    /* find all faces incident to (u1, v1) and all faces incident to (u2, v2), compute intersection,
     * which must contain exactly one triangle. if not, throw exception */

    Face* u1_v1_faces[2];
    Face* u2_v2_faces[2];
    size_t nIF1 = 2;
    size_t nIF2 = 2;
    this->getFacesIncidentToEdge(u1_it, v1_it, u1_v1_faces, nIF1);
    this->getFacesIncidentToEdge(u2_it, v2_it, u2_v2_faces, nIF2);

    Face* f = NULL;
    size_t commonFaces = 0;
    for (size_t i = 0; i < nIF1; ++i)
    {
        for (size_t j = 0; j < nIF2; ++j)
        {
            if (u1_v1_faces[i] == u2_v2_faces[j])
            {
               ++commonFaces;
               f = u1_v1_faces[i];
               break;
            }
        }
    }

    /* must contain exactly one triangle */
    if (commonFaces != 1) {
        debugl(1, "Mesh::getSharedTri(): edges (%5d, %5d) / (%5d, %5d) share %du faces, not exactly one.\n",
           u1_it->id(), v1_it->id(), u2_it->id(), v2_it->id(), commonFaces);
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getSharedTri(): given two edges are not contained in exactly one face.");
    }

    if (f->isQuad()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::getSharedTri(): face containing two given edges is a quad, not a triangle");
    }
    else return f->iterator();
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::scale(R const &r)
{
    for (auto &v : this->vertices) {
        v.pos() *= r;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::translate(Vec3<R> const &d)
{
    for (auto &v : this->vertices) {
        v.pos() += d;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::splitEdge(
    vertex_iterator     u_it,
    vertex_iterator     v_it,
    R const            &lambda)
{
    vertex_iterator w_it;
    Vec3<R>         u_pos, v_pos, w_pos;

    face_iterator   uv_fst_tri_it, uv_snd_tri_it;
    bool            uv_fst_tri_orientation;

    vertex_iterator fst_rem_vertex, snd_rem_vertex;
    vertex_iterator fst_v0_it, fst_v1_it, fst_v2_it, snd_v0_it, snd_v1_it, snd_v2_it;

    /* check edge */
    this->checkEdge("Mesh::splitEdge", u_it, v_it);

    /* compute and add new vertex, lambda must be in [0,1] */
    if (lambda <= 0.0 || lambda >= 1.0) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::splitEdge(): given lambda must be in (0, 1).");
    }

    u_pos   = u_it->pos();
    v_pos   = v_it->pos();
    w_pos   = u_pos + (v_pos - u_pos)*lambda;
    w_it    = this->vertices.insert(w_pos);
    
    /* get two common faces of edge {u, v}, create four new faces and delete the old two */
    this->getFacesIncidentToManifoldEdge(u_it, v_it, uv_fst_tri_it, uv_snd_tri_it);

    /* we need the "remaining" vertices in uv_fst_tri and uv_snd_tri, that is to say the vertex of
     * both triangles that is neither u nor v. additionally, we need the orientation of the edge
     * {u, v} in, say, the first triangle uv_fst_tri. the rest is then uniquely determined */
    uv_fst_tri_it->getTriIterators(fst_v0_it, fst_v1_it, fst_v2_it);
    uv_snd_tri_it->getTriIterators(snd_v0_it, snd_v1_it, snd_v2_it);

    fst_rem_vertex          = Mesh::Face::getTriRemainingVertex(fst_v0_it, fst_v1_it, fst_v2_it, u_it, v_it);
    snd_rem_vertex          = Mesh::Face::getTriRemainingVertex(snd_v0_it, snd_v1_it, snd_v2_it, u_it, v_it);
    uv_fst_tri_orientation  = uv_fst_tri_it->getTriEdgeOrientation(u_it->id(), v_it->id());

    /* delete old two triangles */
    this->faces.erase(uv_fst_tri_it);
    this->faces.erase(uv_snd_tri_it);

    /* if orientation is positive, new triangles are:
     *
     * (u, w, fst_rem_vertex), (w, v, fst_rem_vertex)
     *
     * and 
     *
     * (v, w, snd_rem_vertex), (w, u, snd_rem_vertex).
     *
     * and if orientation is negative, swap u and v in the above triangles */
    if (uv_fst_tri_orientation) {
        this->faces.insert(u_it, w_it, fst_rem_vertex);
        this->faces.insert(w_it, v_it, fst_rem_vertex);

        this->faces.insert(v_it, w_it, snd_rem_vertex);
        this->faces.insert(w_it, u_it, snd_rem_vertex);
    }
    else {
        this->faces.insert(v_it, w_it, fst_rem_vertex);
        this->faces.insert(w_it, u_it, fst_rem_vertex);

        this->faces.insert(u_it, w_it, snd_rem_vertex);
        this->faces.insert(w_it, v_it, snd_rem_vertex);
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::splitEdge(
    vertex_iterator         u_it,
    vertex_iterator         v_it,
    std::vector<R> const   &lambda_values)
{
    debugl(1, "Mesh::splitEdge(): e = (%d, %d)\n", u_it->id(), v_it->id() );
    debugTabInc();

    uint32_t const n = lambda_values.size();
    if (n == 0) {
        debugTabDec();
        throw("Mesh::splitEdge(): vector of split points (lambda_values) is empty. won't silently discard ineffective call.");
    }

    Vec3<R>         u_pos, v_pos, n_pos;

    face_iterator   uv_fst_tri_it, uv_snd_tri_it;
    bool            uv_fst_tri_orientation;

    vertex_iterator fst_rem_vertex, snd_rem_vertex;
    vertex_iterator fst_v0_it, fst_v1_it, fst_v2_it, snd_v0_it, snd_v1_it, snd_v2_it;

    /* check edge */
    this->checkEdge("Mesh::splitEdge", u_it, v_it);

    /* all lambda values must be in ]0,1[ */
    for (auto &lambda : lambda_values) {
        if (lambda <= 0.0 || lambda >= 1.0) {
            debugTabDec();
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::splitEdge(): givne lambda must be in (0, 1).");
        }
    }

    u_pos   = u_it->pos();
    v_pos   = v_it->pos();

    /* compute new vertex positions, add them and store iterators */
    std::vector<vertex_iterator> new_vertex_its(n);
    for (uint32_t i = 0; i < n; i++) {
        n_pos               = u_pos + (v_pos - u_pos)*lambda_values[i];
        new_vertex_its[i]   = this->vertices.insert(n_pos);
    }
    
    /* get two common faces of edge {u, v}, create four new faces and delete the old two */
    this->getFacesIncidentToManifoldEdge(u_it, v_it, uv_fst_tri_it, uv_snd_tri_it);

    /* we need the "remaining" vertices in uv_fst_tri and uv_snd_tri, that is to say the vertex of
     * both triangles that is neither u nor v. additionally, we need the orientation of the edge
     * {u, v} in, say, the first triangle uv_fst_tri. the rest is then uniquely determined */
    uv_fst_tri_it->getTriIterators(fst_v0_it, fst_v1_it, fst_v2_it);
    uv_snd_tri_it->getTriIterators(snd_v0_it, snd_v1_it, snd_v2_it);

    fst_rem_vertex          = Mesh::Face::getTriRemainingVertex(fst_v0_it, fst_v1_it, fst_v2_it, u_it, v_it);
    snd_rem_vertex          = Mesh::Face::getTriRemainingVertex(snd_v0_it, snd_v1_it, snd_v2_it, u_it, v_it);
    uv_fst_tri_orientation  = uv_fst_tri_it->getTriEdgeOrientation(u_it->id(), v_it->id());

    /* delete old two triangles */
    this->faces.erase(uv_fst_tri_it);
    this->faces.erase(uv_snd_tri_it);

    /* if orientation is positive, new triangles are:
     *
     * 1.   (u, new_vertex_its[0], fst_rem_vertex)
     *
     *      (new_vertex_its[0], u, snd_rem_vertex)
     *
     * 2.   (new_vertex_its[i], new_vertex_its[i+1], fst_rem_vertex)        for i = 0..n-2
     *
     *      (new_vertex_its[i + 1], new_vertex_its[i], snd_rem_vertex)
     *
     * 3.   (new_vertex_its[n-1], v, fst_rem_vertex) 
     *
     *      (v, new_vertex_its[n-1], snd_rem_vertex)
     *
     * and if orientation is negative, swap fst_rem_vertex and snd_rem_vertex in the above triangles. */
    if (!uv_fst_tri_orientation) {
        std::swap(fst_rem_vertex, snd_rem_vertex);
    }

    this->faces.insert( u_it, new_vertex_its[0], fst_rem_vertex);
    this->faces.insert( new_vertex_its[0], u_it, snd_rem_vertex);

    for (uint32_t i = 0; i < n - 1; i++) {
        this->faces.insert( new_vertex_its[i], new_vertex_its[i+1], fst_rem_vertex );
        this->faces.insert( new_vertex_its[i+1], new_vertex_its[i], snd_rem_vertex );
    }

    this->faces.insert(new_vertex_its[n-1], v_it, fst_rem_vertex);
    this->faces.insert(v_it, new_vertex_its[n-1], snd_rem_vertex);

    debugTabDec();
    debugl(1, "Mesh::splitEdge(): done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::splitEdge(
    uint32_t    u_id,
    uint32_t    v_id,
    R const    &lambda)
{
    Mesh::vertex_iterator u_it, v_it;

    u_it = this->vertices.find(u_id);
    v_it = this->vertices.find(v_id);

    if (u_it == this->vertices.end() || v_it == this->vertices.end()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::splitEdge() (uint32_t id version): at least one given input id for edge (u, v) is invalid.");
    }
    this->splitEdge(u_it, v_it, lambda);
}


template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::split_face_with_center
(
    uint32_t face_id,
    const Vec3<R>& splitPos
)
{
    face_iterator fit = faces.find(face_id);
    if (fit == faces.end())
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::split_face_with_center(): Given input id for face is invalid.");

    // get corner vertices
    std::vector<vertex_iterator> corners(3);
    if (fit->isTri())
        fit->getTriIterators(corners[0], corners[1], corners[2]);
    else if (fit->isQuad())
    {
        corners.resize(4);
        fit->getQuadIterators(corners[0], corners[1], corners[2], corners[3]);
    }
    else
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::split_face_with_center(): "
            "Face associated to input id is neither triangle not quadrilateral.");

    // remove face
    faces.erase(fit);

    // create new vertex
    vertex_iterator center = vertices.insert(splitPos);

    // create new faces
    size_t nCorner = corners.size();
    for (size_t i = 0; i < nCorner; ++i)
        faces.insert(corners[i], corners[(i+1)%nCorner], center);
}


/* split quad such that the sum of aspect ratios of the two resulting triangles is minimized */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::splitQuadIntoTwoTris(face_iterator quad)
{
    using Aux::Numbers::fmin3; 
    using Aux::Numbers::fmax3; 

    debugl(3, "Mesh::splitQuadIntoTwoTris().\n");

    if (!quad.checkContainer( *this )) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::splitQuadIntoTwoTris(): given face_iterator does not refer to (this) mesh.");
    }
    else if (!quad->isQuad()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::splitQuadIntoTwoTris(): face identified by given input face_iterator is not quad.");
    }
    
    Vec3<R>         v0, v1, v2, v3;
    vertex_iterator v0_it, v1_it, v2_it, v3_it;
    R               l01, l12, l23, l30, l02, l13;
    R               split_02_ar1, split_02_ar2, split_13_ar1, split_13_ar2;

    quad->getQuadIterators(v0_it, v1_it, v2_it, v3_it);
    quad->getQuadPositions(v0   , v1   , v2   , v3   );

    l01     = (v1 - v0).len2();
    l12     = (v2 - v1).len2();
    l23     = (v3 - v2).len2();
    l30     = (v0 - v3).len2();

    l02     = (v2 - v0).len2();
    l13     = (v3 - v1).len2();

    /* aspect ratio of triangle (0, 1, 2) */
    split_02_ar1    =   fmax3(l01, l12, l02) / fmin3(l01, l12, l02);

    /* aspect ratio of triangle (0, 2, 3) */
    split_02_ar2    =   fmax3(l02, l23, l30) / fmin3(l02, l23, l30);

    /* aspect ratio of triangle (0, 1, 3) */
    split_13_ar1    =   fmax3(l01, l13, l30) / fmin3(l01, l13, l30);

    /* aspect ratio of triangle (1, 2, 3) */
    split_13_ar2    =   fmax3(l12, l23, l13) / fmin3(l12, l23, l13);

    /* delete old quad, insert two new triangles  */
    this->faces.erase(quad);

    /* split along diagonal 0 - 2 */
    if ( (split_02_ar1 + split_02_ar2) < (split_13_ar1 + split_13_ar2) ) {
        this->faces.insert(v0_it, v1_it, v2_it);
        this->faces.insert(v0_it, v2_it, v3_it);
    }
    /* split along diagonal 1 - 3 */
    else {
        this->faces.insert(v0_it, v1_it, v3_it);
        this->faces.insert(v1_it, v2_it, v3_it);
    }
}

/* split quad along the shorter diagonal */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::splitQuadIntoTwoTrisShorterDiagonal(face_iterator quad)
{
    using Aux::Numbers::fmin3; 
    using Aux::Numbers::fmax3; 

    debugl(3, "Mesh::splitQuadIntoTwoTrisShorterDiagonal().\n");

    if (!quad.checkContainer( *this )) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::splitQuadIntoTwoTrisShorterDiagonal(): given face_iterator does not refer to (this) mesh.");
    }
    else if (!quad->isQuad()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::splitQuadIntoTwoTrisShorterDiagonal(): face identified by given input face_iterator is not quad.");
    }
    
    Vec3<R>         v0, v1, v2, v3;
    vertex_iterator v0_it, v1_it, v2_it, v3_it;
    R               diag_02_len, diag_13_len;

    quad->getQuadIterators(v0_it, v1_it, v2_it, v3_it);
    quad->getQuadPositions(v0   , v1   , v2   , v3   );

    diag_02_len = (v2 - v0).len2();
    diag_13_len = (v3 - v1).len2();

    /* delete old quad, insert two new triangles  */
    this->faces.erase(quad);

    /* split along diagonal 0 - 2 */
    if ( diag_02_len < diag_13_len ) {
        this->faces.insert(v0_it, v1_it, v2_it);
        this->faces.insert(v0_it, v2_it, v3_it);
    }
    /* split along diagonal 1 - 3 */
    else {
        this->faces.insert(v0_it, v1_it, v3_it);
        this->faces.insert(v1_it, v2_it, v3_it);
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::triangulateQuads()
{
    /* NOTE: face iterators are required NOT to be invalidated by inserting or deleting elements
     * other than the one the iterator points to.. this does the trick here, since splitQuad()
     * deletes the quad face and inserte two new triangles, possibly under totally different id's
     * waiting on the id queue, which returns the smallest currently free id.. */
    face_iterator fit = this->faces.begin(), fit_erase_copy;

    while (fit != this->faces.end()) {
        /* handle only quads */
        if (fit->isQuad()) {
            /* copy iterator, increment it, then destroy the quad being pointed to */
            fit_erase_copy = fit;
            ++fit;
            this->splitQuadIntoTwoTris(fit_erase_copy);
        }
        else ++fit;
    }
}
    
template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::collapseTriEdge(
    vertex_iterator     u_it,
    vertex_iterator     v_it,
    vertex_iterator    *vnew_it,
    Vec3<R>            *new_pos)
{
    using namespace Aux;

    debugl(1, "Mesh::collapseTriEdge(): u_id: %5d, v_id: %5d\n", u_it->id(), v_it->id());
    debugTabInc();

    Vec3<R>             w_pos;
    vertex_iterator     w_it;
    face_iterator       uv_fst_tri_it, uv_snd_tri_it;

    /* get two incident triangles to edge {u, v}. this might throw an exception if the edge does not
     * exist or is a non-manifold edge (which is also required) */
    this->getFacesIncidentToManifoldEdge(u_it, v_it, uv_fst_tri_it, uv_snd_tri_it);

    /* assert that the two faces are indeed triangles */
    uv_fst_tri_it->checkTri("Mesh::collapseTriEdge()");
    uv_snd_tri_it->checkTri("Mesh::collapseTriEdge()");

    /* check topological saftey of collapse. criterion: the set of common neighbours of both u and v
     * (which excludes u and v themselves) has size two. */
    std::list<uint32_t> u_nbs, v_nbs;
    std::list<uint32_t> shared_nbs;

    u_it->getVertexStarIndices(u_nbs);
    v_it->getVertexStarIndices(v_nbs);

    /* collapse topologically safe => proceed */
    if (Aux::Alg::listIntersection<uint32_t>(u_nbs, v_nbs, shared_nbs) == 2) {
        debugl(2, "collapse topologically safe. proceeding.\n");
        uint32_t fst_v0, fst_v1, fst_v2, snd_v0, snd_v1, snd_v2;

        uv_fst_tri_it->getTriIndices(fst_v0, fst_v1, fst_v2);
        uv_snd_tri_it->getTriIndices(snd_v0, snd_v1, snd_v2);

        debugl(3, "deleting common faces: %5d = (%5d, %5d, %5d), %5d = (%5d, %5d, %5d)\n", 
                uv_fst_tri_it->id(), fst_v0, fst_v1, fst_v2,
                uv_snd_tri_it->id(), snd_v0, snd_v1, snd_v2);

        /* delete two faces incident to manifold edge {u, v} */
        this->faces.erase(uv_fst_tri_it);
        this->faces.erase(uv_snd_tri_it);

        /* add new vertex w, representing {u, v} after the collapse. if new position is specified,
         * take it. otherwise default to average position of u and v */
        if (new_pos) {
            w_pos = *new_pos;
        }
        else {
            w_pos = (u_it->pos() + v_it->pos()) * 0.5;
        }
        w_it = this->vertices.insert(w_pos);
        
        /* write iterator w_it to new vertex w if desired by the caller */
        if (vnew_it) {
            *vnew_it = w_it;
        }

        debugl(3, "created new vertex w with id %5d and pos (%5.4f, %5.4f, %5.4f)\n",
                w_it->id(), w_pos[0], w_pos[1], w_pos[2]);

        /* the new vertex w is incident to all faces that u and v were incident to, excluding the
         * two deleted ones, which have already been deleted. similarly, as adjacent vertices, w has
         * all neighbours of u and v combined, again excluding the neighbour relations from the old
         * two triangles.  gather that information. sort(), but not unique() the adjacency lists,
         * because R      entries indicate "multiple" shared edges, i.e. non-bordre edges contained
         * in multiple faces */
        Vertex *u_vertex = &(*u_it);
        Vertex *v_vertex = &(*v_it);
        Vertex *w_vertex = &(*w_it);

        /* create map for vertex replacement */
        std::map<Vertex *, Vertex *> replace_map = { {u_vertex, w_vertex}, {v_vertex, w_vertex} };

        /* w's incident faces is the union of u's and v's. */
        debugl(3, "getting incident faces of u and v and computing union as new incident faces of w.\n");
        w_vertex->incident_faces = u_vertex->incident_faces;
        w_vertex->incident_faces.insert(
                w_vertex->incident_faces.end(),
                v_vertex->incident_faces.begin(),
                v_vertex->incident_faces.end());

        w_vertex->incident_faces.sort([] (const Face* x, const Face* y) -> bool {return (x->id() < y->id());});
        w_vertex->incident_faces.unique([] (const Face* x, const Face* y) -> bool {return (x->id() == y->id());});

        /* replace vertex pointers to u and v with pointers to v in all incident faces */
        debugl(3, "new vertex w's incident faces: replacing vertex pointers to u/v with pointers to w.\n");
        debugl(4, "union of faces with old ids:\n");
        debugTabInc();
        for (Face *f : w_vertex->incident_faces) {
            debugl(1, "%5d = (%5d, %5d, %5d)\n", f->id(), f->vertices[0]->id(), f->vertices[1]->id(), f->vertices[2]->id() );
            f->replaceVertices(replace_map);
        } 
        debugTabDec();

#ifdef __DEBUG__
        debugl(4, "union of faces with new ids:\n");
        debugTabInc();
        for (Face *f : w_vertex->incident_faces) {
            debugl(1, "%5d = (%5d, %5d, %5d)\n", f->id(), f->vertices[0]->id(), f->vertices[1]->id(), f->vertices[2]->id() );
        } 
        debugTabDec();
#endif

        /* w's adjacent vertices are the union of all vertices that were incident to u and v */
        debugl(3, "merging lists of incident vertices of u and v.. updating pointers.\n");
        w_vertex->adjacent_vertices = u_vertex->adjacent_vertices;
        w_vertex->adjacent_vertices.insert(w_vertex->adjacent_vertices.end(), v_vertex->adjacent_vertices.begin(), v_vertex->adjacent_vertices.end());
        w_vertex->adjacent_vertices.sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});

        /* replace pointers to u/v with pointers to w in all adjacent vertices */
        debugl(3, "new vertex w's adjacent vertices: replacing vertex pointers to u/v with pointers to w.\n");
        std::list<Vertex *> w_adj_vertices;
        w_vertex->getVertexStar(w_adj_vertices);
        for (Vertex *w_nb : w_adj_vertices) {
            w_nb->replaceAdjacentVertices(replace_map);
        }

        debugl(3, "manually clear()ing topological information in u/v and erase()ing u and v from mesh.\n");
        /* finally, delete the vertices u and v, which have now been entirely replaced by w.
         * manually reset adjacency / incidence data to prevent VertexAccessor::erase() from erasing
         * the faces u and v are (were) contained in */
        u_it->adjacent_vertices.clear();
        u_it->incident_faces.clear();
        v_it->adjacent_vertices.clear();
        v_it->incident_faces.clear();

        this->vertices.erase(u_it);
        this->vertices.erase(v_it);

        /* done */
        debugTabDec();
        debugl(1, "Mesh::collapseTriEdge(): done. edge collapsed.\n");

        return true;
    }
    /* collapse topologically unsafe. don't proceed, return false. */
    else {
        debugl(2, "collapse topologically unsafe => not going to proceeed.\n");
        debugTabDec();
        return false;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
Mesh<Tm, Tv, Tf, R>::mergeUnrelatedVertices(
    vertex_iterator     u_it,
    vertex_iterator     v_it,
    const Vec3<R>      *new_pos)
{
    using namespace Aux::Alg::Functional;

    debugl(3, "Mesh::mergeUnrelatedVertices(): u = %5d, v = %5d.\n", u_it->id(), v_it->id() );
    debugTabInc();

    /* firstly, verify that u and v satisfy the condition of being toplogically "unrelated" */
    if ( !u_it.checkContainer( *this ) || !u_it.sameContainer(v_it) || u_it == v_it)
    {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::mergeUnrelatedVertices(): vertices u and v not both referring to (this) or identical.");
    }

    std::list<Face *> u_fstar, v_fstar;
    u_it->getFaceStar(u_fstar);
    v_it->getFaceStar(v_fstar);

    bool uv_unrelated = true;
    for (auto &uf : u_fstar) {
        uv_unrelated &= !uf->contains(v_it);
    }
    for (auto &vf : v_fstar) {
        uv_unrelated &= !vf->contains(u_it);
    }

    if (!uv_unrelated) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::mergeUnrelatedVertices(): given input vertices not unrelated.");
    }

    /* perform the merge, add new vertex with average of the two positions or *new_pos if new_pos !=
     * NULL */
    Vec3<R> w_pos;
    if (new_pos) {
        w_pos = *new_pos;
    }
    else {
        w_pos = (u_it->pos() + v_it->pos()) * 0.5;
    }

    Mesh::vertex_iterator w_it = this->vertices.insert(w_pos);

    /* manually copy adjacent_vertices and incident from both u and v to w and clear info inside u
     * and v. */
    w_it->adjacent_vertices = u_it->adjacent_vertices;
    w_it->adjacent_vertices.insert(
            w_it->adjacent_vertices.end(),
            v_it->adjacent_vertices.begin(),
            v_it->adjacent_vertices.end());

    /* sort with custom static comparison function for pointers */
    w_it->adjacent_vertices.sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});

    /* clear adjacent_vertices of u and v */
    u_it->adjacent_vertices.clear();
    v_it->adjacent_vertices.clear();

    /* analogous procedure for faces */
    w_it->incident_faces = u_it->incident_faces;
    w_it->incident_faces.insert(
            w_it->incident_faces.end(),
            v_it->incident_faces.begin(),
            v_it->incident_faces.end());
    w_it->incident_faces.sort([] (const Face* x, const Face* y) -> bool {return (x->id() < y->id());});
    w_it->incident_faces.unique([] (const Face* x, const Face* y) -> bool {return (x->id() == y->id());});

    u_it->incident_faces.clear();
    v_it->incident_faces.clear();

    /* replace all occurrences of u and v in all incident faces of the new vertex w and in the
     * adjacency lists of all neighbours of w (which is the union of the old neighbours of u and v
     * and still think they're connected to u or v). */
    std::map<Vertex *, Vertex *> replace_map = { { &(*u_it), &(*w_it)}, { &(*v_it), &(*w_it)} };
    for (auto &f : w_it->incident_faces) {
        f->replaceVertices(replace_map);
    } 

    std::list<Vertex *> w_vstar;
    w_it->getVertexStar(w_vstar);
    for (auto &w_nb : w_vstar) {
        w_nb->replaceAdjacentVertices(replace_map);
    }

    /* erase u and v. since these are now isolated, this won't remove any faces. */
    this->vertices.erase(u_it);
    this->vertices.erase(v_it);

    debugTabDec();
    debugl(3, "Mesh::mergeUnrelatedVertices(): done.\n");

    /* return iterator to new vertex w. */
    return w_it;
}

/* the mesh is 2-manifold iff every edge is contained in exactly two vertices and all faces incident
 * to a vertex (the face star) forms a single closed circle in the dual graph */
template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::checkTopology()
{
    debugl(1, "Mesh::checkTopology()\n");
    debugTabInc();

    Vertex                 *s, *v;
    std::list<Vertex *>     v_vstar;
    std::list<Face *>       v_fstar;

    std::queue<Vertex *>    Q;

    uint8_t                 u_state;
    uint32_t                ccs, cc_size;

    const uint32_t          traversal_id = this->getFreshTraversalId();

    ccs         = 0;
    debugl(2, "checking for orphaned vertices..\n");

    debugTabInc();
    for (auto &v : this->vertices) {
        v.getVertexStar(v_vstar);
        v.getFaceStar(v_fstar);

        if (v_vstar.empty() || v_fstar.empty()) {
            debugl(1, "vertex %d: incident_faces: %2d, adjacent_vertices: %2d\n",
                v.id(), v_fstar.size(), v_vstar.size());
        }
    }
    debugTabDec();

    debugl(2, "checking if mesh is triangular..\n");
    debugTabInc();
    /* checking if all faces are triangles */
    for (auto &f : this->faces) {
        if (!f.isTri()) {
            debugl(1, "Mesh::checkTopology(): quad with id %d found.\n", f.id());
            throw("Mesh::checkTopology(): quad found.");
        }
    }
    debugTabDec();

    /* traverse all connected components, check manifold property for edges */
    debugl(2, "checking for non-manifold edges by traversing all connected components.\n");
    debugTabInc();
    for (Vertex &vref : this->vertices) {
        if (vref.getTraversalState(traversal_id) == TRAV_UNSEEN) {
            /* new connected component with root s */
            s       = &vref;
            cc_size = 0;
            ccs++;

            /* traverse cc */
            Q.push(s);

            debugl(1, "traversing new connected component with root vertex %5d\n", s->id() );
            debugTabInc();
            while (!Q.empty()) {
                v = Q.front();
                Q.pop();
                cc_size++;

                /* get neighbours of v */
                v_vstar.clear();
                v->getVertexStar(v_vstar);

                for (auto u : v_vstar) {
                    /* current neighbour u */
                    u_state = u->getTraversalState(traversal_id);

                    /* we handle the edge in in the direction {v -> u} if u is not done (but can be
                     * either enqueued or unseen), otherwise u has already handled the edge in
                     * direction {u -> v} */
                    if (u_state != TRAV_DONE) {
                        Face* vu_incident_faces[2];
                        size_t nIF = 2;
                        /* get faces containing edge {v, u} */
                        this->getFacesIncidentToEdge(v->iterator(), u->iterator(), vu_incident_faces, nIF);
                        if (nIF > 2) {
                            debugl(1, "Mesh::checkTopology(): FATAL: non-manifold edge (%d, %d) found. faces: \n", v->id() , u->id() );
                            throw("Mesh::checkTopology(): FATAL: non-manifold edge found.");
                        }
                        else if (nIF == 1) {
                            debugl(1, "Mesh::checkTopology(): FATAL: border edge (%d, %d) found!\n", v->id(), u->id() );
                            throw("Mesh::checkTopology(): FATAL: border edge found.");
                        }
                        else if (nIF == 0) {
                            debugl(1, "Mesh::checkTopology(): FATAL: isolated edge (%d, %d) found.\n", v->id(), u->id());
                            throw("Mesh::checkTopology(): FATAL: isolated edge found.");
                        }

                        /* get two faces, check if their are triangles, get orientation of shared
                         * edge, check for orientability */
                        else if (nIF == 2) {
                            Face   *t1 = vu_incident_faces[0];
                            Face   *t2 = vu_incident_faces[1];
                            bool    t1_uv_orient, t2_uv_orient;
                            
                            /* get orientations in both triangles */
                            t1_uv_orient = t1->getTriEdgeOrientation(u->id(), v->id() );
                            t2_uv_orient = t2->getTriEdgeOrientation(u->id(), v->id() );

                            if (t1_uv_orient == t2_uv_orient) {
                                debugl(1, "Mesh::checkTopology(): triangels %d / %d sharing edge (%d, %d) not compatible in orientation!\n", t1->id(), t2->id(), u->id(), v->id());
                                throw("Mesh::checkTopology(): triangles with incompatible orientation found.");
                            }

                        }
                        else {
                            throw("Mesh::checkTopology(): logical error: must-be-impossible case for edge_faces.size().");
                        }
                    }
                    
                    /* enqueue u if it hasn't been seen yet */
                    if (u_state == TRAV_UNSEEN) {
                        Q.push(u);
                        u->setTraversalState(traversal_id, TRAV_ENQUEUED);
                    }
                }

                /* v is done */
                v->setTraversalState(traversal_id, TRAV_DONE);
            }
            debugTabDec();

            debugl(1, "finished traversal of connected component with root vertex %5d. size: %5d\n", s->id(), cc_size);
        }
    }
    debugTabDec();

    if (ccs > 1) {
        debugl(2, "WARNING: mesh has %d > 1 connected components.\n", ccs);
    }

    std::list<Face *>                       v_dual_circle;
    Face                                   *last_on_circle;
    vertex_iterator                         e_x, e_y;
    bool                                    got_neighbour;

    /* check manifold property for vertices */
    debugl(2, "checking for non-manifold vertices.\n");
    debugTabInc();
    for (auto &v : this->vertices) {
        /* clear and get all incident faces */
        v_fstar.clear();
        v.getFaceStar(v_fstar);

        /* must be at least three elements */
        if (v_fstar.size() < 3) {
            debugl(1, "Mesh::checkTopology(): vertex %d got < 3 incident triangles.\n", v.id() );
            throw("Mesh::checkTopology(); vertex with < 3 incident triangles faces found.");
        }

        /* clear dual circle list */
        v_dual_circle.clear();

        /* get off first face and push it as first "vertex", i.e. face, in  v_dual_circle */
        v_dual_circle.push_front(v_fstar.front());
        v_fstar.pop_front();

        /* always take last element from dual circle, search for neighbour  */
        while ( !v_fstar.empty() ) {
            last_on_circle  = v_dual_circle.back();
            got_neighbour   = false;

            /* scan remaining incident faces for neighbour of last_on_circle */
            for (auto ifit = v_fstar.begin(); ifit != v_fstar.end(); ++ifit) {
                /* if *ifit is a neighbour, remove it from list and append it to dual cycle */
                if ( (this->getTriCommonEdge(last_on_circle->iterator(), (*ifit)->iterator(), e_x, e_y)) == true) {
                    /* one of e_x and e_y must be v itself */
                    if (e_x->id() != v.id() && e_y->id() != v.id() ) {
                        debugl(1, "Mesh::checkTopology(): face in face star of vertex %d shares common edge with other face from face star, but vertex %d is not contained in in.\n", v.id(), v.id() );
                        throw("Mesh::checkTopology(): face in face star of vertex shares common edge with other face from face star, but vertex v is not contained in in.\n");
                    }

                    /* got the neighbour, delete *ifit and push_back on dual circle */
                    got_neighbour = true;

                    v_dual_circle.push_back(*ifit);
                    v_fstar.erase(ifit);
                    break;
                }
            }

            /* incident_faces non-empty and no closed cycle => non-manifold vertex. */
            if (!got_neighbour) {
                debugl(1, "Mesh::checkTopology(): non-manifold vertex %d found.\n", v.id() );
                throw("Mesh::checkTopology(): non-manifold vertex found.");
            }
        }

        /* check if last element of v_dual_circle shares edge with first element of v_dual_circle */
        if (!this->getTriCommonEdge(v_dual_circle.front()->iterator(), v_dual_circle.back()->iterator(), e_x, e_y)) {
            debugl(1, "Mesh::checkTopology(): non-manifold vertex %d found. first and last face in v_dual_circle not neighbours => open path.\n", v.id() );
            throw("Mesh::checkTopology(): non-manifold vertex found. first and last face in v_dual_circle not neighbours => open path.\n");
        }

        /* we're getting paranoid here.. */
        if (v_dual_circle.front() == v_dual_circle.back()) {
            debugl(1, "Mesh::checkTopology(): front() and back() of dual circle identical. must not be.\n");
            throw("Mesh::checkTopology(): front() and back() of dual circle identical. must not be.");
        }
        /* the currently analysed vertex v is regular or "manifold" */
    }
    debugTabDec();

    debugTabDec();
    debugl(1, "Mesh::checkTopology(): done. mesh topology ok.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::checkGeometry(
    R const &eps_v,
    R const &eps_T)
{
    debugl(1, "Mesh::checkGeometry(): updating Octree..\n");
    debugTabInc();

    /* update octree if necessary */
    this->updateOctree();

    Vec3<R>                 v_pos, n_pos;
    R                       v_n_dist;
    std::list<Vertex *>     v_close_vertices;

    vertex_iterator         T_v0_it, T_v1_it, T_v2_it, n_v0_it, n_v1_it, n_v2_it;
    Vec3<R>                 T_v0, T_v1, T_v2, n_v0, n_v1, n_v2;
    std::list<Face *>       T_close_tris;

    /* vertices of shared edge and remainin vertices */
    vertex_iterator         x_it, u_it, v_it, y_it; //v_id declared above
    Vec3<R>                 x, u, v, y;

    /* for two triangles sharing exactly one vertex:
     * s:       shared vertex
     * a0,a1:   remaining vertices of tri 1
     * b0,b1:   remaining vertices of tri 2
     * */
    vertex_iterator         s_it, a0_it, a1_it, b0_it, b1_it;
    Vec3<R>                 s, a0, a1, b0, b1; 

    printf("Mesh::checkGeometry(): checking for coincident vertices with threshold %3.2e and self-intersecting faces..\n", eps_v);
    uint32_t    nchecks = 0;
    uint32_t    F       = this->numFaces();
    uint32_t    tenperc = F / 10;


    /* check for ''coincident'' or rather too close vertices */
    for (auto &v : this->vertices) {
        /* get bounding box of v */
        auto v_bb = v.getBoundingBox();

        /* search for vertices in the vicinity */
        v_close_vertices.clear();
        this->findVertices(v_bb, v_close_vertices);

        /* check all of them, excepting v itself */
        for (Vertex *n : v_close_vertices) {
            if (n->id() == v.id() ) continue;
            else {
                v_n_dist = ( n->pos() - v.pos() ).len2();
                if ( v_n_dist < eps_v) {
                    printf("Mesh::checkGeometry(): dist between vertices %d, %d = %4.2e < %4.2e\n", v.id(), n->id(), v_n_dist, eps_v);
                    throw("Mesh::checkGeometry(): two vertices with distance below threshold.");
                }
            }
        }
    }
    printf("Mesh::checkGeometry(): no coincident vertices (within threshold) found.\n");


    /* check all triangles for degenerated edges by computing the three altitudes (heights) and
     * comparing them to eps_T */
    R       a, b, c, sp, Rc;
    R       h_a, h_b, h_c;

    printf("Mesh::checkGeometry(): checking for faces with any altitude < %3.2e\n", eps_T);
    for (auto &T: this->faces) {
        T.getTriPositions(T_v0, T_v1, T_v2);

        a       = (T_v1 - T_v0).len2();
        b       = (T_v2 - T_v1).len2();
        c       = (T_v0 - T_v2).len2();
        sp      = 0.5*(a + b + c);
        Rc      = (a*b*c) / ( 4.0 * std::sqrt( sp*(a + b - sp)*(a + c - sp)*(b + c - sp) ) );

        h_a     = (0.5*b*c) / Rc;
        h_b     = (0.5*a*c) / Rc;
        h_c     = (0.5*a*b) / Rc;

        if (h_a < eps_T || h_b < eps_T || h_c < eps_T) {
            printf("Mesh::checkGeometry(): triangle %d, side length: %3.2e, %3.2e, %3.2e, AR: %3.5e  has altitudes (%3.2e, %3.2e, %3.2e), one of which < eps_T = %3.2e\n",
                T.id(), a, b, c, T.getTriAspectRatio(), h_a, h_b, h_c, eps_T);
            throw("Mesh::checkGeometry(): triangle has altitude < eps_T.");
        }
    }
    printf("Mesh::checkGeometry(): no such faces found.\n");

    printf("Mesh::checkGeometry(): checking for self-intersections..\n");
    /* check for self-intersections on all faces. */
    for (auto &T : this->faces) {
        if (nchecks % tenperc == 0) {
            printf("\t %3u%% done.\n", (uint32_t)(round( 100.0 * (R)nchecks/(R)F )) );
        }


        /* get iterators, vertex positions and bounding box of T */
        T.getTriIterators(T_v0_it, T_v1_it, T_v2_it);
        T.getTriPositions(T_v0, T_v1, T_v2);
        auto T_bb = T.getBoundingBox();

        /* locate all close triangles */
        T_close_tris.clear();
        this->findFaces(T_bb, T_close_tris);

        /* distinguish three cases:
         * 1. triangles sharing an edge                     => special test
         * 2. triangles not sharing an edge but a vertex.   => special test
         * 3. triangles not sharing anything                => generic tritri test
         * */
        for (auto n : T_close_tris) {
            /* skip T itself */
            if (n->id() == T.id() ) continue;
            else {
                /* get vertex ids of n */
                n->getTriIterators(n_v0_it, n_v1_it, n_v2_it);

                /* check if they share an edge */
                if (Mesh::Face::getTwoTrianglesSharedEdgeAndRemainingVertices(
                            T_v0_it,
                            T_v1_it,
                            T_v2_it,
                            n_v0_it,
                            n_v1_it,
                            n_v2_it,
                            /* shared vertices */
                            u_it, v_it,
                            /* remaining vertex in T */
                            x_it, 
                            /* remaining vertex in n */
                            y_it) == true)
                {
                    //printf("tri: %d, %d, sharing an edge %d, %d\n", T_id, n->id, u_id, v_id);
                    /* get vertices and check for self-intersection */
                    u = u_it->pos();
                    v = v_it->pos();
                    x = x_it->pos();
                    y = y_it->pos();
                    bool si = Aux::Geometry::triTriSharingEdge<R>(x, u, v, y);
                    if (si) {
                        printf("Mesh::checkGeometry(): tri: %d, %d sharing edge %d, %d are self-intersecting.\n",
                                T.id(), n->id(), u_it->id(), v_it->id());
                        throw("Mesh::checkGeometry(): two triangles sharing an edge are self-intersecting.");
                    }
                }
                /* check if they share exactly one vertex */
                else if (Mesh::Face::getTwoTrianglesCommonVertexAndRemainingVertices(
                            T_v0_it,
                            T_v1_it,
                            T_v2_it,
                            n_v0_it,
                            n_v1_it,
                            n_v2_it,
                            /* shared vertices */
                            s_it,
                            /* two vertices remaining for T */
                            a0_it, a1_it,
                            /* two vertices remaining for n */
                            b0_it, b1_it) == true)
                {
                    /* get vertices */
                    a0  = a0_it->pos();
                    a1  = a1_it->pos();  
                    s   = s_it->pos();   
                    b0  = b0_it->pos();  
                    b1  = b1_it->pos();  

                    //printf("tris %d %d sharing vertex %d => check\n", T_id, n->id, s_id);
                    bool si = Aux::Geometry::triTriSharingVertex<R>(a0, a1, s, b0, b1);
                    if (si) {
                        printf("Mesh::checkGeometry(): two triangles %d, %d sharing only vertex %d are self-intersecting.\n",
                                T.id(), n->id(), s_it->id() );
                        throw("Mesh::checkGeometry(): two triangles sharing only one vertex self-intersecting.\n");
                    }
                }
                /* otherwise generic Moeller tri tri test, which should be robust enough */
                else {
                    n->getTriPositions(n_v0, n_v1, n_v2);
                    bool si = Aux::Geometry::triTri3d(T_v0, T_v1, T_v2, n_v0, n_v1, n_v2);
                    if (si) {
                        printf("Mesh::checkGeometry(): two non-incident triangels %d, %d intersecting => self-intersection.\n", T.id(), n->id() );
                        throw("Mesh::checkGeometry(): two non-incident triangels intersecting => self-intersection.");
                    }
                }
            }
        }
        nchecks++;
    }
    printf("Mesh::checkGeometry(): no self-intersections found.\n");

    /* geometry ok */
    debugTabDec();
    debugl(1, "Mesh::checkGeometry(): geometry ok => no coincident vertices within threshold, no self-intersections.\n");
}

/* I/O */
struct obj_face {
    bool        quad;
    uint32_t    v_idx[4];
    uint32_t    n_idx[4];

    obj_face(
        bool        quad,
        bool        use_normals,
        uint32_t    v0,
        uint32_t    v1,
        uint32_t    v2,
        uint32_t    v3,
        uint32_t    n0 = 0,
        uint32_t    n1 = 0,
        uint32_t    n2 = 0,
        uint32_t    n3 = 0)
    {
        this->quad      = quad;
        this->v_idx[0]  = v0 - 1;
        this->v_idx[1]  = v1 - 1;  
        this->v_idx[2]  = v2 - 1;  
        if (quad) {
            this->v_idx[3]  = v3 - 1;
        }

        if (use_normals) {
            n_idx[0]    = n0 - 1;
            n_idx[1]    = n1 - 1;
            n_idx[2]    = n2 - 1;
            if (quad) {
                n_idx[3]    = n3 - 1;
            }
        }
        else {
            n_idx[0] = 0;
            n_idx[1] = 0;
            n_idx[2] = 0;
            n_idx[3] = 0;
        }
    }
};

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::readFromObjFile(const char *filename)
{
    debugl(1, "Mesh::readFromObjFile()");
    debugTabInc();

    /* clear mesh */
    this->clear();

    std::list<Vec3<R>>      vertices;
    std::list<Vec3<R>>      normals;
    std::list<obj_face>     faces;

    std::string             line_string;
    char                   *line;
    std::ifstream           f;

    f.open(filename, std::ifstream::in);
    if ( !f.is_open() ) {
        throw MeshEx(MESH_IO_ERROR, "Mesh::readFromObjFile(): can't open input file\n");
    }

    R           x, y, z;
    uint32_t    v0, v1, v2, v3, n0, n1, n2, n3;

    /* keep reading until eof */
    while ( f.good() ) {
        /* get line, convert to c string */
        std::getline(f, line_string);
        if (line_string.size() == 0) {
            debugl(3, "empty line..\n");
            continue;
        }

        line = new char [ line_string.size() + 1];
        std::strcpy(line, line_string.c_str() );

        /* parse line with sscanf, if not a comment */
        if (line[0] == '#') {
            //printf("comment: ");
        }
        /* object declaration */
        else if (line[0] == 'o') {
            //printf("comment: ");
        }
        /* dummy texture coordinate */
        else if (line[0] == 'v' && line[1] == 't') {
        }
        else if ( sscanf(line, "v %lf %lf %lf", &x, &y, &z) == 3 ) {
            //printf("vertex.");
            vertices.push_back( Vec3<R>(x, y, z) );
        }
        else if ( sscanf(line, "vn %lf %lf %lf", &x, &y, &z) == 3) {
            //printf("normal.");
            normals.push_back( Vec3<R>(x, y, z) );
        }
        else if ( sscanf(line, "f %u %u %u", &v0, &v1, &v2) == 3) {
            //printf("simple triangle.");
            faces.push_back( obj_face(false, false, v0, v1, v2, 0) );
        }
        else if ( sscanf(line, "f %u %u %u %u", &v0, &v1, &v2, &v3) == 4) {
            //printf("simple quad.");
            throw("Mesh::readFromObjFile(): QUAD!\n");
            //faces.push_back( obj_face(true, false, v0, v1, v2, v3) );
        }
        else if ( sscanf(line, "f %u//%u %u//%u %u//%u", &v0, &n0, &v1, &n1, &v2, &n2) == 6) {
            //printf("triangle with normals.");
            faces.push_back( obj_face(false, true, v0, v1, v2, 0, n0, n1, n2, 0) );
        }
        else if ( sscanf(line, "f %u//%u %u//%u %u//%u %u//%u", &v0, &n0, &v1, &n1, &v2, &n2, &v3, &n3) == 8) {
            //printf("quad with normals.");
            throw("Mesh::readFromObjFile(): QUAD!\n");
            //faces.push_back( obj_face(true, true, v0, v1, v2, v3, n0, n1, n2, n3) );
        }
        else {
            printf("Mesh::readFromObjFile(): unrecognized line: \"%s\".\n", line);
            throw("Mesh::readFromObjFile(): unrecognized line.\n");
        }

        delete[] line;
    }

    /* try to add all vertices and faces. */
    try {
        /* NOTE: this relies upon the fact that adding n vertices to an empty mesh will number them 0...(n-1) */
        debugl(2, "Adding %5d vertices\n", (uint32_t)vertices.size());
        while (!vertices.empty()) {
            this->vertices.insert(vertices.front());
            vertices.pop_front();
        }
        debugl(2, "done adding vertices.\n");

        debugl(2, "adding %5d faces..\n", (uint32_t)faces.size());
        obj_face *f;
        while (!faces.empty()) {
            f = &(faces.front());
            if (f->quad) {
                this->faces.insert(
                        f->v_idx[0],
                        f->v_idx[1],
                        f->v_idx[2],
                        f->v_idx[3]
                    );
            }
            /* obj reader only supports quads and triangles, so !f->quad <=> f is triangle */
            else {
                this->faces.insert(
                        f->v_idx[0],
                        f->v_idx[1],
                        f->v_idx[2]
                    );
            }

            /* pop_front added face */
            faces.pop_front();
        }
        debugl(2, "done adding faces.\n");
    }
    catch (MeshEx& err) {
        printf("caught exception: \"%s\".\n", err.error_msg.c_str() );
    }

    debugTabDec();
    debugl(1, "Mesh::readFromObjFile(): done reading mesh from obj: numVertices(): %d, numFaces(): %d, numEdges(): %d\n", this->numVertices(), this->numFaces(), this->numEdges() );
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::writeObjFile(const char *jobname)
{
    debugl(4, "Mesh::writeObjFile(): writing mesh as obj to outfile \"%s\".\n", jobname);
    debugTabInc();

    /* renumber vertices */
    this->renumberConsecutively();

    char obj_filename[512];
    snprintf(obj_filename, 512, "%s.obj", jobname);
    FILE *outfile = fopen(obj_filename, "w");

    if (!outfile) {
        debugl(1, "Mesh::writeObjFile(): can't open file \'%s\' for writing.\n", obj_filename);
        throw("Mesh::writeObjFile(): can't open output file for writing.");
    }

    fprintf(outfile, "# obj file automatically generated by AnaMorph for jobname: \"%s\".\n", jobname);
    fprintf(outfile, "o %s\n", jobname);

    /* write in all the vertices, preceeded by a comment */
    fprintf(outfile, "\n# %15ld vertices.\n", this->V.size());

    debugl(4, "writing %15ld vertices..\n", this->vertices.size());
    debugTabInc();

    Vec3<R> vpos;
    for (auto &v : this->vertices) {
        vpos = v.pos();
        debugl(5, "writing vertex %5d..\n", v.id());
        vpos.print_debugl(5);
        fprintf(outfile, "v %+.10e %+.10e %+.10e\n", vpos[0], vpos[1], vpos[2]);
    }

    debugTabDec();

    /* write a dummy texture coordinate to increase compatability with the somewhat ill-defined
     * wavefront format. some readers don't accept empty texture coordinates */
    fprintf(outfile, "\n# dummy texture coordinate to increase compatability with several programs importing .obj files.\n");
    fprintf(outfile, "vt 0.0 0.0\n");

    /* write in all the face normals, preceeded by a comment */
    fprintf(outfile, "\n# %15ld face normals.\n", this->faces.size());

    debugl(4, "writing %15ld face normals..\n", this->faces.size());
    debugTabInc();

    Vec3<R> n;
    for (auto &f : this->faces) {
        debugl(5, "writing face normal %5d..\n", f.id());
        n = f.getNormal();
        fprintf(outfile, "vn %+.10e %+.10e %+.10e\n", n[0], n[1], n[2]);
    }

    debugTabDec();

    /* newline, comment and then all faces */
    fprintf(outfile, "\n# %15ld faces.\n", this->F.size() );

    debugl(4, "writing %15ld faces..\n", this->faces.size());
    debugTabInc();

    uint32_t v0_id, v1_id, v2_id, v3_id;
    for (auto &f: this->faces) {
        if (f.isQuad()) {
            f.getQuadIndices(v0_id, v1_id, v2_id, v3_id);
            fprintf(outfile, "f %u//%u %u//%u %u//%u %u//%u\n",
                    v0_id + 1, 
                    v0_id + 1, 
                    v1_id + 1,
                    v1_id + 1,
                    v2_id + 1,
                    v2_id + 1,
                    v3_id + 1,
                    v3_id + 1);

            debugl(5, "quad  %u/%u/%u/%u\n",
                    v0_id, 
                    v1_id,
                    v2_id,
                    v3_id);
        }
        else if (f.isTri()) {
            f.getTriIndices(v0_id, v1_id, v2_id);
            fprintf(outfile, "f %u//%u %u//%u %u//%u\n",
                    v0_id + 1, 
                    v0_id + 1, 
                    v1_id + 1,
                    v1_id + 1,
                    v2_id + 1,
                    v2_id + 1);

            debugl(5, "tri  %u/%u/%u/\n",
                    v0_id, 
                    v1_id,
                    v2_id);
        }
        else {
            fclose(outfile);
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::writeObjFile(): discovered face that is neither triangle nor quad. internal logic error.");
        }
    }
    debugTabDec();

    fclose(outfile);

    debugTabDec();
    debugl(4, "Mesh::writeObjFile(): done.\n");
}

template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::VertexAccessor::VertexAccessor(Mesh<Tm, Tv, Tf, R> &m) : mesh(m) 
{
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::begin()
{
    return Mesh::vertex_iterator( &(this->mesh), this->mesh.V.begin());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_const_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::begin() const
{
    /* the seemingly returned vertex_iterator will be implicitly converted to a
     * vertex_const_iterator to match the return type */
    return Mesh::vertex_iterator( &(this->mesh), this->mesh.V.begin());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::end()
{
    return Mesh::vertex_iterator( &(this->mesh), this->mesh.V.end());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_const_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::end() const
{
    return Mesh::vertex_const_iterator( &(this->mesh), this->mesh.V.end());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::find(const uint32_t &id)
{
    return Mesh::vertex_iterator( &(this->mesh), this->mesh.V.find(id));
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_const_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::find(const uint32_t &id) const
{
    return Mesh::vertex_const_iterator( &(this->mesh), this->mesh.V.find(id));
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::Vertex &
Mesh<Tm, Tv, Tf, R>::VertexAccessor::at(const uint32_t &id)
{
    try {
        return *(this->mesh.V.at(id));
    }
    catch (std::out_of_range &e) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::VertexAccessor::at(): vertex index out of range.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::VertexAccessor::exists(const uint32_t &id) const
{
    return (this->mesh.V.find(id) != this->mesh.V.end());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::insert(const Vec3<R> &vpos)
{
    std::pair<
            typename std::map<uint32_t, VertexPointerType >::iterator,
            bool
        >                                                               pair;

    typename std::map<uint32_t, VertexPointerType >::iterator           vit;

    /* get fresh id for new vertex, allocate new vertex, insert pair (id, vertex) into map */
    uint32_t v_id   = this->mesh.V_idq.getId();
    Vertex *v       = new Vertex(&(this->mesh), vpos);
    pair            = this->mesh.V.insert( { v_id, VertexPointerType(v) } );
    if (!pair.second) {
        throw MeshEx(MESH_LOGIC_ERROR, "vertex with fresh id from idq already present in vertex map. this must never happen..");
    }
    else {
        /* the (private) constructor Vertex(const Vec3 &pos) does not set the internal iterator
         * Vertex::m_vit, since it has no way of knowing it in advance. set it now to bring the new
         * vertex into a consistent state. */
        vit                         = pair.first;
        v->m_vit                    = vit;

        /* mesh octree needs update */
        this->mesh.octree_updated   = false;

        /* return iterator */
        return (vit->second->iterator());
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
Mesh<Tm, Tv, Tf, R>::VertexAccessor::erase(vertex_iterator it)
{
    debugl(3, "Mesh::VertexAccessor::erase(vertex_iterator it): vertex id: %6d. deleting %2d incident faces..\n", it->id(), it->incident_faces.size() );
    debugTabInc();

    if (!it.checkContainer( this->mesh)) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::VertexAccessor::erase(): given input iterator refers to different container.");
    }

    /* delete all faces that the vertex pointed to by "it" is contained in.  NOTE: the erase call
     * will alter vit->incident_faces: always take the begin() iterator, which is guaranteed always
     * to be valid as long as set is non-empty */
    Face *incident_face;
    while (!it->incident_faces.empty()) {
        incident_face = *(it->incident_faces.begin());
        this->mesh.faces.erase( incident_face->iterator() );
    }

    /* now it->incident_faces should be empty */
    if (it->adjacent_vertices.size() > 0) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::deleteVertex(): vertex still got incident faces or adjacent vertices after deleteing all incident faces.. o_O");
    }

    debugl(4, "freeing id and erasing vertex from internal vertex map..\n");

    /* free id */
    this->mesh.V_idq.freeId(it->id());

    debugl(4, "deleting (deallocating) vertex object..\n");
    /* delete allocated vertex object */
    delete &(*it);

    /* mesh octree needs update */
    this->mesh.octree_updated = false;

    debugTabDec();
    debugl(3, "Mesh::VertexAccessor::erase(). erase()ing and returning vertex_iterator to next vertex.\n");

    /* erase vertex from internal map this->mesh.V and return vertex_iterator to the next element
     * (by wrapping the result of std::map::erase in the vertex_iterator constructor). */
    return Mesh::vertex_iterator( &(this->mesh), this->mesh.V.erase(it.int_it) );
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::VertexAccessor::erase(const uint32_t &id)
{
    /* try to get iterator and call iterator version of erase */
    Mesh::vertex_iterator vit;
    if ( (vit = this->find(id)) != this->end() ) {
        this->erase(vit);
        return true;

    }
    else {
        debugTabDec();
        return false;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
size_t
Mesh<Tm, Tv, Tf, R>::VertexAccessor::size() const
{
    return (this->mesh.V.size());
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::VertexAccessor::empty() const
{
    return (this->mesh.V.empty());
}

template <typename Tm, typename Tv, typename Tf, typename R>
Mesh<Tm, Tv, Tf, R>::FaceAccessor::FaceAccessor(Mesh<Tm, Tv, Tf, R> &m) : mesh (m) 
{
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::begin()
{
    return Mesh::face_iterator( &(this->mesh), this->mesh.F.begin() );
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_const_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::begin() const
{
    /* face_iterator will be implicitly converted to face_const_iterator */
    return Mesh::face_iterator( &(this->mesh), this->mesh.F.begin() );
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::end()
{
    return Mesh::face_iterator( &(this->mesh), this->mesh.F.end());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_const_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::end() const
{
    /* face_iterator will be implicitly converted to face_const_iterator */
    return Mesh::face_iterator( &(this->mesh), this->mesh.F.end());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::find(const uint32_t &id)
{
    return Mesh::face_iterator( &(this->mesh), this->mesh.F.find(id));
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_const_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::find(const uint32_t &id) const
{
    return Mesh::face_const_iterator( &(this->mesh), this->mesh.F.find(id));
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::Face &
Mesh<Tm, Tv, Tf, R>::FaceAccessor::at(const uint32_t &id)
{
    try {
        return *(this->mesh.F.at(id));
    }
    catch (std::out_of_range &e) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::at(): vertex index out of range.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::FaceAccessor::exists(const uint32_t &id) const
{
    return (this->mesh.F.find(id) != this->mesh.F.end());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::insert(
    const std::list<uint32_t> &id_list)
{
    if (id_list.size() == 3) {
        uint32_t    v0_id, v1_id, v2_id;
        auto        lit = id_list.begin();

        v0_id = *lit;
        v1_id = *++lit;
        v2_id = *++lit;

        return ( this->insert(v0_id, v1_id, v2_id) );
    }
    else if (id_list.size() == 4) {
        uint32_t    v0_id, v1_id, v2_id, v3_id;
        auto        lit = id_list.begin();

        v0_id = *lit;
        v1_id = *++lit;
        v2_id = *++lit;
        v3_id = *++lit;

        return ( this->insert(v0_id, v1_id, v2_id, v3_id) );
    }
    else {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert(): given index list contains neither 3 nor 4 indices. only triangles and quads supported.");
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::insert(
    vertex_iterator v0_it,
    vertex_iterator v1_it,
    vertex_iterator v2_it)
{
    uint32_t                        tri_id;
    Face                           *tri;
    Vertex                         *v0, *v1, *v2;
    std::pair<
        typename std::map<
            uint32_t,
            FacePointerType
        >::iterator,
        bool>                       rpair;

    /* at least check whether all iterators refer to (this) mesh! */
    if (!( v0_it.checkContainer(this->mesh) && v0_it.sameContainer(v1_it) && v0_it.sameContainer(v2_it) )) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert() (tri version): invalid input iterators: at least two iterators referring to different meshes found.");
    }

    /* check whether all three are distinct */
    if (v0_it == v1_it || v0_it == v2_it || v1_it == v2_it) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert() (tri version): at least two of the given three iterators are identical.");
    }

    /* get face pointers */
    v0      = &(*v0_it);
    v1      = &(*v1_it);
    v2      = &(*v2_it);

    /* get fresh id for new triangle */
    tri_id  = this->mesh.F_idq.getId();
    tri     = new Face(&(this->mesh), false, v0, v1, v2, NULL);

    /* insert into map, directly set iterator inside newly created Face */
    rpair   = this->mesh.F.insert( {tri_id, FacePointerType(tri) } );
    if (!rpair.second) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert(): new Face with fresh id from idq already present in Face map. this must never happen..");
    }

    /* set Face::m_fit iterator, which is required for Face to be in a consistent state and has not
     * been set by the (private) Face ctor, just as for Mesh::Vertex */
    tri->m_fit  = rpair.first;

    /* vertex ids can be in adjacent_vertices multiple times, for two vertices can be an edge of two
     * incident faces. when getAdjacentIndices/Vertices() is called, the unique() list is computed.
     * when deleting a face however, this is convenient: one simply deletes only ONE instance of
     * adjacent vertices. if the two vertices are an edge in two faces, they will still be adjacent
     * afterwards. if unique all the time, we'd have to check the surrounding faces every time a
     * face is deleted..  incident_faces is a set, however. same for quads below */
    v0->insertAdjacentVertex(v1);
    v0->insertAdjacentVertex(v2);
    v0->insertIncidentFace(tri);

    v1->insertAdjacentVertex(v0);
    v1->insertAdjacentVertex(v2);
    v1->insertIncidentFace(tri);

    v2->insertAdjacentVertex(v0);
    v2->insertAdjacentVertex(v1);
    v2->insertIncidentFace(tri);

    this->mesh.octree_updated = false;

    /* return iterator to newly inserted tri */
    return (tri->iterator());
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::insert(
    const uint32_t &v0_id, 
    const uint32_t &v1_id, 
    const uint32_t &v2_id)
{
    vertex_iterator v0_it, v1_it, v2_it, vend;

    v0_it   = this->mesh.vertices.find(v0_id);
    v1_it   = this->mesh.vertices.find(v1_id);
    v2_it   = this->mesh.vertices.find(v2_id);
    vend    = this->mesh.vertices.end();

    if (v0_it == vend || v1_it == vend || v2_it == vend) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert(): at least one of the given input ids out of range.");
    }
    else return (this->insert(v0_it, v1_it, v2_it));
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::insert(
    vertex_iterator v0_it,
    vertex_iterator v1_it,
    vertex_iterator v2_it,
    vertex_iterator v3_it)
{
    uint32_t                        quad_id;
    Face                           *quad;
    Vertex                         *v0, *v1, *v2, *v3;
    std::pair<
        typename std::map<
            uint32_t,
            FacePointerType >
        ::iterator,
        bool>                       rpair;

    /* check whether all three iterators refer to (this) mesh! */
    if (!( v0_it.checkContainer(this->mesh) && v0_it.sameContainer(v1_it) && v0_it.sameContainer(v2_it) && v0_it.sameContainer(v3_it) )) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert() (quad version): invalid input iterators: at least two iterators referring to different meshes found.");
    }

    /* check whether all four are distinct */
    if (v0_it == v1_it || v0_it == v2_it || v0_it == v3_it || v1_it == v2_it || v1_it == v3_it || v2_it == v3_it) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert() (quad version): at least two of the given frou iterators are identical.");
    }

    /* get face pointers */
    v0      = &(*v0_it);
    v1      = &(*v1_it);
    v2      = &(*v2_it);
    v3      = &(*v3_it);

    /* get fresh id for new triangle */
    quad_id = this->mesh.F_idq.getId();
    quad    = new Face( &(this->mesh), true, v0, v1, v2, v3);

    /* insert into map, directly set iterator inside newly created Face */
    rpair   = this->mesh.F.insert( { quad_id, FacePointerType(quad) } );
    if (!rpair.second) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert(): new Face with fresh id from idq already present in Face map. this must never happen..");
    }

    /* set Face::m_fit iterator, which is required for Face to be in a consistent state and has not
     * been set by the (private) Face ctor, just as for Mesh::Vertex */
    quad->m_fit = rpair.first;

    /* topology information update */
    v0->insertAdjacentVertex(v3);
    v0->insertAdjacentVertex(v1);
    v0->insertIncidentFace(quad);

    v1->insertAdjacentVertex(v0);
    v1->insertAdjacentVertex(v2);
    v1->insertIncidentFace(quad);

    v2->insertAdjacentVertex(v1);
    v2->insertAdjacentVertex(v3);
    v2->insertIncidentFace(quad);

    v3->insertAdjacentVertex(v2);
    v3->insertAdjacentVertex(v0);
    v3->insertIncidentFace(quad);

    /* mesh octree needs update */
    this->mesh.octree_updated = false;

    /* return iterator to newly inserted quad */
    return ( quad->iterator() );
}

template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::insert(
    const uint32_t &v0_id,
    const uint32_t &v1_id,
    const uint32_t &v2_id,
    const uint32_t &v3_id)
{
    vertex_iterator v0_it, v1_it, v2_it, v3_it, vend;

    v0_it   = this->mesh.vertices.find(v0_id);
    v1_it   = this->mesh.vertices.find(v1_id);
    v2_it   = this->mesh.vertices.find(v2_id);
    v3_it   = this->mesh.vertices.find(v3_id);
    vend    = this->mesh.vertices.end();

    if (v0_it == vend || v1_it == vend || v1_it == vend || v2_it == vend || v3_it == vend) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::insert(): at least one of the given input ids out of range.");
    }
    else return (this->insert(v0_it, v1_it, v2_it, v3_it));
}


template <typename Tm, typename Tv, typename Tf, typename R>
typename Mesh<Tm, Tv, Tf, R>::face_iterator
Mesh<Tm, Tv, Tf, R>::FaceAccessor::erase(face_iterator it)
{
    using Aux::Alg::removeFirstOccurrenceFromList;

    debugl(3, "Mesh::FaceAccessor::erase(): erasing face with it: %6d\n", it->id());
    debugTabInc();
    bool all_erased;

    it->checkTriQuad("Mesh::FaceAccessor::erase()");

    /* get pointer to face */
    Face *f = &(*it);

    /* quads */
    if ( f->isQuad() ) {
        Vertex *v_i, *v_j, *v_k, *v_l;
        f->getQuadVertices(v_i, v_j, v_k, v_l);

        /* erase face_id from std::SET of incident faces */
        all_erased = 
            (v_i->eraseIncidentFace(f)) &&
            (v_j->eraseIncidentFace(f)) &&
            (v_k->eraseIncidentFace(f)) && 
            (v_l->eraseIncidentFace(f));

        if (!all_erased) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::erase(): could not erase face pointer from the incident_faces set of at least one contained vertex. internal logic error.");
        }

        /* now remove vertex adjacencies created by this edge. the adjacency lists might contain
         * multiple duplicate entries, because an edge might be incident to several faces. remove
         * only ONE copy from the adjacency list. with a set, this would be very problematic */
        auto sortFct = [] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());};
        removeFirstOccurrenceFromList(v_i->adjacent_vertices, v_l);
        removeFirstOccurrenceFromList(v_i->adjacent_vertices, v_j);
        v_i->adjacent_vertices.sort(sortFct);

        removeFirstOccurrenceFromList(v_j->adjacent_vertices, v_i);
        removeFirstOccurrenceFromList(v_j->adjacent_vertices, v_k);
        v_j->adjacent_vertices.sort(sortFct);

        removeFirstOccurrenceFromList(v_k->adjacent_vertices, v_j);
        removeFirstOccurrenceFromList(v_k->adjacent_vertices, v_l);
        v_k->adjacent_vertices.sort(sortFct);

        removeFirstOccurrenceFromList(v_l->adjacent_vertices, v_k);
        removeFirstOccurrenceFromList(v_l->adjacent_vertices, v_i);
        v_l->adjacent_vertices.sort(sortFct);
    }
    /* same for triangles */
    else if (f->isTri()) {
        Vertex *v_i, *v_j, *v_k;
        f->getTriVertices(v_i, v_j, v_k);

        /* erase face_id from incidence set */
        all_erased = 
            (v_i->eraseIncidentFace(f)) &&
            (v_j->eraseIncidentFace(f)) &&
            (v_k->eraseIncidentFace(f));

        if (!all_erased) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::erase(): could not erase face pointer from the incident_faces set of at least one contained vertex. internal logic error.");
        }

        /* remove one occurrence of adjacencies from the face */
        auto sortFct = [] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());};
        removeFirstOccurrenceFromList(v_i->adjacent_vertices, v_k);
        removeFirstOccurrenceFromList(v_i->adjacent_vertices, v_j);
        v_i->adjacent_vertices.sort(sortFct);

        removeFirstOccurrenceFromList(v_j->adjacent_vertices, v_i);
        removeFirstOccurrenceFromList(v_j->adjacent_vertices, v_k);
        v_j->adjacent_vertices.sort(sortFct);

        removeFirstOccurrenceFromList(v_k->adjacent_vertices, v_j);
        removeFirstOccurrenceFromList(v_k->adjacent_vertices, v_i);
        v_k->adjacent_vertices.sort(sortFct);
    }
    else { 
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::FaceAccessor::erase(): supplied face is neither quad nor triangle. general case intentionally unsupported right now => internal logic error.");
    }

    /* free face_id */
    this->mesh.F_idq.freeId( it->id() );

    /* delete allocated face object */
    delete &(*it);

    /* mesh octree needs update */
    this->mesh.octree_updated = false;

    debugTabDec();

    /* return face_iterator to next element by wrapping return iterator of map::erase inside a face_iterator */
    return Mesh::face_iterator( &(this->mesh), this->mesh.F.erase(it.int_it));
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::FaceAccessor::erase(const uint32_t &id)
{
    /* locate face with given id. if found, remove it with iterator version of erase */
    face_iterator fit = this->find(id);
    if ( fit != this->end() ) {
        this->erase(fit);
        return true;
    }
    /* face not found, return false. */
    else {
        return false;
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
size_t
Mesh<Tm, Tv, Tf, R>::FaceAccessor::size() const
{
    return (this->mesh.F.size());
}

template <typename Tm, typename Tv, typename Tf, typename R>
bool
Mesh<Tm, Tv, Tf, R>::FaceAccessor::empty() const
{
    return (this->mesh.F.empty());
}

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::resetTraversalStates()
{
    /* clear id queue for traversal ids, get first id, reset all states with that id */
    this->traversal_idq.clear();
    uint32_t first_id = this->traversal_idq.getId();

    for (auto &v : this->vertices) {
        v.setTraversalState(first_id, TRAV_UNSEEN);
    }

    for (auto &f : this->faces) {
        f.setTraversalState(first_id, TRAV_UNSEEN);
    }
}

template <typename Tm, typename Tv, typename Tf, typename R>
uint32_t
Mesh<Tm, Tv, Tf, R>::getFreshTraversalId()
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

template <typename Tm, typename Tv, typename Tf, typename R>
void
Mesh<Tm, Tv, Tf, R>::checkInternalConsistency() const
{
    debugl(1, "Mesh::checkInternalConsistency()\n");
    /* iterate over internal vertex / face maps, check if mesh pointers and iterators inside the
     * objects are consistent. if shared pointers are used, check if the reference counter for each
     * shared_ptr is exactly one */
    debugTabInc();
    Vertex  *v;
    std::list<Vertex *> vptr_list;
    for (auto vit = this->V.begin(); vit != this->V.end(); ++vit) {
        v = vit->second;
        if (v->mesh != this || v->m_vit != vit || v->id() != vit->first) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::checkInternalConsistency(): found vertex with inconsistent mesh pointer or internal iterator. internal logic error.");
        }
        /* only if shared_ptr<Vertex> is VertexPointerType */
        /*
        if (vit->second.use_count() != 1) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::checkInternalConsistency(): found vertex whose shared_ptr in internal vertex map has use_count != 1. internal logic error.");
        }
        debugl(2, "vertex id %d smart ptr use count: %d\n", v->id(), vit->second.use_count());
        */
        vptr_list.push_back(v);
    }

    vptr_list.sort([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() < y->id());});
    vptr_list.unique([] (const Vertex* x, const Vertex* y) -> bool {return (x->id() == y->id());});
    if (vptr_list.size() != this->V.size()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::checkInternalConsistency(): at least one vertex (pointer) stored under two differend ids..\n");
    }

    Face  *f;
    std::list<Face *> fptr_list;
    for (auto fit = this->F.begin(); fit != this->F.end(); ++fit) {
        f = fit->second;
        if (f->mesh != this || f->m_fit != fit || f->id() != fit->first) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::checkInternalConsistency(): found face with inconsistent mesh pointer or internal iterator. internal logic error.");
        }
        /* only if shared_ptr<Face> is FacePointerType */
        /*
        if (fit->second.use_count() != 1) {
            throw MeshEx(MESH_LOGIC_ERROR, "Mesh::checkInternalConsistency(): found face whose shared_ptr in internal face map has use_count != 1. internal logic error.");
        }
        debugl(2, "face id %d smart ptr use count: %d\n", f->id(), fit->second.use_count());
        */
        fptr_list.push_back(f);
    }

    fptr_list.sort([] (const Face* x, const Face* y) -> bool {return (x->id() < y->id());});
    fptr_list.unique([] (const Face* x, const Face* y) -> bool {return (x->id() == y->id());});

    if (fptr_list.size() != this->F.size()) {
        throw MeshEx(MESH_LOGIC_ERROR, "Mesh::checkInternalConsistency(): at least one face (pointer) stored under two differend ids..\n");
    }
    debugTabDec();
    debugl(1, "Mesh::checkInternalConsistency(): mesh internally consistent.\n");
}
