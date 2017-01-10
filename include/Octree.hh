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

#ifndef OCTREE_H
#define OCTREE_H

#include <vector>
#include <list>

#include "CubeCommon.hh"

/* flags bit:
 *  0       : leaf bool
 *  1..3    : child type
 *  4       : last sibling
 *  */

/* 28.02.2013: old thingy ate up 50% up running time. reimplement with std::vector<uint8_t> for
 * testing */
class OtId {
    public:
        //uint8_t                 prefix_len; 
        uint8_t                 len;
        std::vector<uint8_t>    vec;

                                OtId();

        /* copy constructor necessary, since implicitly generated copy constructor copies pointers
         * instead of contents */
                                OtId(const OtId &x);
                               ~OtId();
        uint8_t                 getPrefixLength() const;
        uint8_t                 size() const;
        void                    clear();
        void                    push_front(uint8_t x);
        void                    push(uint8_t x);
        void                    pop();
        void                    set(int16_t pos, uint8_t v);
        uint8_t                 front();
        uint8_t                 back();
        void                    print() const;
        void                    print_nodebug(bool prefix = true) const;

        OtId                   &operator=(const OtId &x);
        uint8_t                 operator[](int16_t pos) const;
        bool                    operator<(const OtId &x) const;
        bool                    operator>(const OtId &x) const;
        bool                    operator==(const OtId &x) const;
};

/* specialization of std::hash<T> for T = OtId copied from boost hash_combine() */
namespace std {
    template <>
    struct hash<OtId> {
        size_t
        operator()(const OtId &x) const 
        {
            size_t                                  hash;
            std::hash<uint8_t>                      bytehasher;
            std::vector<uint8_t>::const_iterator    it;

            hash = 0;
            for (it = x.vec.begin(); it != x.vec.end(); ++it) {
                hash ^= bytehasher(*it) + 0x9e3779b9 + (hash << 6) + (hash >> 2); 
            }
            return hash; 
        }
    };
}

template <typename T>
struct OctreeNode {
    T                   data;
    uint32_t            flags;
    OctreeNode<T>      *first_child;
    union {
        OctreeNode<T>   *next_sibling;
        OctreeNode<T>   *parent;
    } uplink;


    OctreeNode()
    : flags(0), first_child(NULL)
    {}

    OctreeNode(T data, uint8_t child_type, bool leaf = false)
    {
        this->data              = data;
        this->flags             = 0;
        this->first_child       = NULL;
        this->uplink.parent     = NULL;

        this->setLastSibling(false);
        this->setChildType(child_type);
        this->setLeaf(leaf);
    }

    ~OctreeNode() {
        //printf("~OctreeNode()\n");
    }


    inline bool
    isLeaf() const
    {
        return (this->flags & 1);
    }

    inline void
    setLeaf(bool is_leaf = true)
    {
        if (is_leaf) {
            this->flags |= 1;
        }
        else {
            this->flags &= ~1;
        }
    }

    inline bool
    isLastSibling() const
    {
        return (this->flags & 16);
    }

    inline void
    setLastSibling(bool last_sibling = true)
    {
        if (last_sibling) {
            this->flags |= 16;
        }
        else {
            this->flags &= ~16;
        }
    }

    inline uint8_t
    getChildType() const
    {
        return ( (this->flags >> 1) & 7);
    }

    inline void
    setChildType(uint8_t child_type)
    {
        if (child_type < 8) {
            /* first, zero child type bits 1..3 */
            this->flags &= ~14;

            /* then OR the new type in */
            this->flags |= (child_type << 1);
        }
        else throw("OctreeNode::setChildType(): invalid child type given.\n");
    }

    OctreeNode<T> *
    getChild(uint8_t child_type)
    {
        OctreeNode<T> *tmp;
        tmp = this->first_child;
        if (tmp) {
            while (tmp != this && tmp->getChildType() != child_type) {
                tmp = tmp->uplink.next_sibling;
            }

            if (tmp != this) {
                return tmp;
            }
            else return NULL;
        }
        else return NULL;
    }

    OctreeNode<T> *
    getParent()
    {
        OctreeNode<T> *tmp;
        tmp = this;
        while (!tmp->isLastSibling()) {
            //debug("still not last sibling..\n");
            tmp = tmp->uplink.next_sibling;
        }
        //debug("last sibling reached. parent = %d\n", tmp->uplink.parent);
        return tmp->uplink.parent;
    }
};

template <typename tT, typename nT, typename R>
class Octree {
    private:
        void octree_generate_leaf_list_recursive(OctreeNode<nT> *n, std::list<OctreeNode<nT> *> &l);
        void octree_generate_leaf_vector_recursive(OctreeNode<nT> *n, std::vector<OctreeNode<nT> *> &v);
        void octree_free_recursive(OctreeNode<nT> *n);
            
    public:
        tT                  data;
        OctreeNode<nT>      root;
        Vec3<R>             root_vertex;
        R                   root_l;

        Octree(tT _data, nT root_data, Vec3<R> _root_vertex, R l)
        : data(_data), root_vertex(_root_vertex), root_l(l)
        {
            /* init octree. root_vertex is the lower left front corner, (x,y,z) an orthonormal base
             * to use, and roo_tl the cube edge length */

            /* init root node */
            this->root.data             = root_data;
            this->root.flags            = 0;
            this->root.first_child      = NULL;
            this->root.uplink.parent    = NULL;

            this->root.setLeaf(true);
            this->root.setChildType(0);
            this->root.setLastSibling(true);
        }

        ~Octree()
        {
            this->free();
        }

        void
        generateLeafList(std::list<OctreeNode<nT> *> &l) 
        {
            octree_generate_leaf_list_recursive(&(this->root), l);
        }

        void
        generateLeafVector(std::vector<OctreeNode<nT> *> &v) 
        {
            octree_generate_leaf_vector_recursive(&(this->root), v);
        }

        void
        free()
        {
            octree_free_recursive(&(this->root));
        }

        void
        freeNode(OctreeNode<nT> *n)
        {
            octree_free_recursive(n);
        }

        void
        freeNonLeafs()
        {
            octree_free_nonleafs_recursive(&(this->root));
        }

        uint32_t
        getDepth(OctreeNode<nT> *n) const
        {
            uint32_t depth = 0;
            OctreeNode<nT>   *node;
            node = n;
            if (node) {
                while ( node != &(this->root) ) {
                    node = node->getParent();
                    if (!node) {
                        throw("getParent() returned NULL before root\n");
                    }
                    depth++;
                }
                return depth;
            }
            else throw("Octree::getDepth(): NULL pointer received.\n");
        }

        void
        getCubeId(OctreeNode<nT> *n, OtId &cube_id) const
        {
            std::list<uint8_t>              list;
            std::list<uint8_t>::iterator    lit;
            OctreeNode<nT>                *tmp;

            /* reset cube_id to empty tuple */
            cube_id.clear();

            /* clear list */
            list.clear();

            /* walk back to the root and save reverse ordered child types in list, then push all list elements
             * into OtId structure, which only supports push() easily due to bitwise hacks, which
             * make it space efficient */
            for (tmp = n; tmp->getParent() != NULL; tmp = tmp->getParent()) {
                list.push_front(tmp->getChildType());
            }
        
            /* now push list contents into node id, in the right order */
            for (lit = list.begin(); lit != list.end(); ++lit) {
                cube_id.push(*lit);
            }
        }

        void
        getCubeCornerZeroAndLength(const OtId &cube_id, Vec3<R> &cube_vertex, R &length) const
        {
            using namespace CubeCommon;

            R  disp;
            uint8_t i;

            /* check if cube_id uses last prefix field. if so, initialize cube_vertex to the right
             * rootcube corner */

            if (cube_id.getPrefixLength() > 0) {
                disp        = this->root_l;
                cube_vertex = this->root_vertex;

                switch(cube_id[-1]) {
                    case V_LEFT_FRONT_DOWN:
                        break;

                    case V_RIGHT_FRONT_DOWN:
                        cube_vertex[0] += disp;
                        break;

                    case V_RIGHT_BACK_DOWN:
                        cube_vertex[0] += disp;
                        cube_vertex[1] += disp;
                        break;

                    case V_LEFT_BACK_DOWN:
                        cube_vertex[1] += disp;
                        break;

                    case V_LEFT_FRONT_UP:
                        cube_vertex[2] += disp;
                        break;

                    case V_RIGHT_FRONT_UP:
                        cube_vertex[0] += disp;
                        cube_vertex[2] += disp;
                        break;

                    case V_RIGHT_BACK_UP:
                        cube_vertex[0] += disp;
                        cube_vertex[1] += disp;
                        cube_vertex[2] += disp;
                        break;

                    case V_LEFT_BACK_UP:
                        cube_vertex[1] += disp;
                        cube_vertex[2] += disp;
                        break;
                }
            }
            else {
                cube_vertex = this->root_vertex;
            }

            disp        = this->root_l * 0.5;
            length      = this->root_l;

            for (i = 0; i < cube_id.size(); i++, disp *= 0.5, length *= 0.5) {
                switch(cube_id[i]) {
                    case V_LEFT_FRONT_DOWN:
                        break;

                    case V_RIGHT_FRONT_DOWN:
                        cube_vertex[0] += disp;
                        break;

                    case V_RIGHT_BACK_DOWN:
                        cube_vertex[0] += disp;
                        cube_vertex[1] += disp;
                        break;

                    case V_LEFT_BACK_DOWN:
                        cube_vertex[1] += disp;
                        break;

                    case V_LEFT_FRONT_UP:
                        cube_vertex[2] += disp;
                        break;

                    case V_RIGHT_FRONT_UP:
                        cube_vertex[0] += disp;
                        cube_vertex[2] += disp;
                        break;

                    case V_RIGHT_BACK_UP:
                        cube_vertex[0] += disp;
                        cube_vertex[1] += disp;
                        cube_vertex[2] += disp;
                        break;

                    case V_LEFT_BACK_UP:
                        cube_vertex[1] += disp;
                        cube_vertex[2] += disp;
                        break;
                }
            }
        }

        void
        getCubeCorners(const OtId &cube_id, Vec3<R> (&corners)[8]) const
        {
            uint32_t   i;
            Vec3<R>     cube_vertex;
            R           cube_length;
            OtId        uid;

            debug("-----------\n");
            for (i = 0; i < 8; i++) {
                uid         = this->getUniqueVertexId(cube_id, i);
                this->getCubeCornerZeroAndLength(uid, cube_vertex, cube_length);

                corners[i]  = cube_vertex;
            }
            debug("cube length: %f\n", cube_length);
            debug("-----------\n");
        }

        void
        getCubeCorners(OctreeNode<nT> *n, Vec3<R> (&corners)[8]) const
        {
            OtId cube_id;

            this->getCubeId(n, cube_id);
            this->getCubeCorners(cube_id, corners);
        }

        void
        getCubeCorner(const OtId &cube_id, uint8_t ctype, Vec3<R> &corner) const
        {
            R cube_length;

            OtId uid = this->getUniqueVertexId(cube_id, ctype);
            this->getCubeCornerZeroAndLength(uid, corner, cube_length);
        }

        bool
        getNeighbourId(uint8_t neighbour_type, const OtId &cube_id, OtId &neighbour_id, bool use_last_prefix_field = false) const
        {
            using namespace CubeCommon;

            uint8_t i, j;
            uint8_t nb1, nb2, nb3;

            /* we're working on a copy of cube_id, altering the components through the algorithm below */
            neighbour_id = cube_id;

            if (isFaceEnum(neighbour_type)) {
                //debug("Octree::getNeighbourId(): face neighbour requested: %d. id:\n", neighbour_type);
                //neighbour_id.print();

                /* backtrack the instances and mirror according to neighbour_type */
                for (j = neighbour_id.size(), i = j - 1; j >= 1 && isFaceVertex(neighbour_id[i], neighbour_type); j--, i--) {
                    //debug("(j, i) = (%d, %d). current component neighbour_id[i = j - 1 = %d] = %d is a face vertex for face %d\n", j, i, i, neighbour_id[i], neighbour_type);
                    neighbour_id.set(i, mirrorVertexFace(neighbour_id[i], neighbour_type));
                }

                /* if (j == 0) after the loop, we've reached the root cube and the child type there is still,
                 * say a F_DOWN for example, if we wanted a DOWN neighbour.
                 * in this case we got a little cube on the "down face fringe" of the root cube. 
                 * he's got no down neighbour neighbour inside this octree -> return false.  */
                if (j == 0) {
                    //debug("(j, i) = (%d, %d): reached root and all child types of face type. no neighbour in this octree.\n", j, i);
                    if (use_last_prefix_field) {
                        /* if face type if F_{RIGHT, BACK, UP}, then mirror the last prefix field
                         * (which is now 0) */
                        if (neighbour_type == F_RIGHT || neighbour_type == F_BACK || neighbour_type == F_UP) {
                            //debug("face type is F_{RIGHT, BACK, UP}: mirroring last prefix field.\n");
                            neighbour_id.set(-1, mirrorVertexFace(neighbour_id[-1], neighbour_type));
                        }
                    }
                    /* notice: return value is still false, since neighbour is not part of the
                     * octree */
                    return false;
                }
                /* otherwise j >= 1 and i >= 0 is the index of the first child type that does not
                 * satisfy isFaceVertex(..) (F_DOWN in the above example). this represents the first
                 * node on the way to the root who's not a DOWN child => 
                 * mirror this chile type using mirrorVertexFace(). */
                else {
                    //debug("first position of non-face-type child-type: (j, i) = (%d, %d). mirroring this one.. (type, mirrored_type) = (%d, %d)\n", j, i, neighbour_id[i], mirrorVertexFace(neighbour_id[i], neighbour_type) );
                    neighbour_id.set(i, mirrorVertexFace(neighbour_id[i], neighbour_type) );
                    //neighbour_id.print();
                    return true;
                }
            }
            else if (isEdgeEnum(neighbour_type)) {
                /* an "edge" neighbour can be reached by two face direction steps */
                //debug("Octree::getNeighbourId(): edge neighbour requested. splitting into multiple face neighbour calls.\n");
                getEdgeNeighbourDirections(neighbour_type, nb1, nb2);

                if (!use_last_prefix_field) {
                    if (!this->getNeighbourId(nb1, neighbour_id, neighbour_id, use_last_prefix_field)) return false;
                    if (!this->getNeighbourId(nb2, neighbour_id, neighbour_id, use_last_prefix_field)) return false;
                    return true;
                }
                else {
                    bool in = true;
                    if (!this->getNeighbourId(nb1, neighbour_id, neighbour_id, use_last_prefix_field)) in = false;
                    if (!this->getNeighbourId(nb2, neighbour_id, neighbour_id, use_last_prefix_field)) in = false;
                    return in;
                }

            }
            else if (isVertexEnum(neighbour_type)) {
                /* a "vertex" neighbour can be reached by three face direction steps */
                //debug("Octree::getNeighbourId(): vertex neighbour requested. splitting into multiple face neighbour calls.\n");
                getVertexNeighbourDirections(neighbour_type, nb1, nb2, nb3);

                if (!use_last_prefix_field) {
                    if (!this->getNeighbourId(nb1, neighbour_id, neighbour_id, use_last_prefix_field)) return false;
                    if (!this->getNeighbourId(nb2, neighbour_id, neighbour_id, use_last_prefix_field)) return false;
                    if (!this->getNeighbourId(nb3, neighbour_id, neighbour_id, use_last_prefix_field)) return false;
                    return true;
                }
                else {
                    bool in = true;
                    if (!this->getNeighbourId(nb1, neighbour_id, neighbour_id, use_last_prefix_field)) in = false;
                    if (!this->getNeighbourId(nb2, neighbour_id, neighbour_id, use_last_prefix_field)) in = false;
                    if (!this->getNeighbourId(nb3, neighbour_id, neighbour_id, use_last_prefix_field)) in = false;
                    return in;
                }
            }
            else throw("Octree::getNeighbourId(): invalid neighbour type given: not a face, edge or vertex enum.\n");
        }

        /* get unique id for a lattice vertex, i.e. a cube corner. we identify that cube corner by
         * the id of the LARGEST cube (cube of maximum edge length) in the octree that has the given
         * vertex as V_LEFT_FRONT_DOWN, i.e. vertex 0. 
         * to get to this cube, we first find a cube of equal size that has the vertex as
         * V_LEFT_FRONT_DOWN, and then go up the octree as far as possible, enlarging the cube.
         * we can go up if the current cube is the left fron down child of its parent, i.e. if its
         * child type is 0 */
        OtId
        getUniqueVertexId(const OtId& cube_id, uint8_t vtype) const {
            using namespace CubeCommon;

            OtId uid;

            /* start cube is the cube of euqal size that has the specified vertex as 
             * V_LEFT_FRONT_DOWN */
            switch (vtype) {
                /* already left front down vertex in given cube */
                case V_LEFT_FRONT_DOWN:
                    uid = cube_id;
                    break;

                /* right neighbour */
                case V_RIGHT_FRONT_DOWN:
                    this->getNeighbourId(F_RIGHT, cube_id, uid, true);
                    break;

                /* right and back */
                case V_RIGHT_BACK_DOWN:
                    this->getNeighbourId(E_RIGHT_BACK, cube_id, uid, true);
                    break;

                /* back neighbour */
                case V_LEFT_BACK_DOWN:
                    this->getNeighbourId(F_BACK, cube_id, uid, true);
                    break;

                /* up neighbour */
                case V_LEFT_FRONT_UP:
                    this->getNeighbourId(F_UP, cube_id, uid, true);
                    break;

                /* right and up */
                case V_RIGHT_FRONT_UP:
                    this->getNeighbourId(E_RIGHT_UP, cube_id, uid, true);
                    break;

                /* right, back and up */
                case V_RIGHT_BACK_UP:
                    this->getNeighbourId(V_RIGHT_BACK_UP, cube_id, uid, true);
                    break;

                /* back and up */
                case V_LEFT_BACK_UP:
                    this->getNeighbourId(E_BACK_UP, cube_id, uid, true);
                    break;
            }

            /* uid identifies the cube that has the vertex as V_LEFT_FRONT_DOWN. as long as
             * the cube is the left, front, down child of its parent, go up to the parent by
             * pop()ing from uid */
            while (uid.size() > 0 && uid.back() == getVertexNumber(V_LEFT_FRONT_DOWN) ) {
                uid.pop();
            }

            return uid;
        }

        bool
        isFringeCube(OctreeNode<nT> *n, const OtId &cube_id, bool edgecubes_only = false) const
        {
            using namespace CubeCommon;

            uint32_t i;
            bool     all_down, all_front, all_right, all_back, all_left, all_up;
            uint8_t  value;
            OtId    n_id;

            all_down    = true;
            all_front   = true;
            all_right   = true;
            all_back    = true;
            all_left    = true;
            all_up      = true;

            /* get n's node id */
            this->getCubeId(n, n_id);

            /* compare the first depth(n) components of n's and the given cube id */
            for (i = 0; i < n_id.size(); i++) {
                if (n_id[i] != cube_id[i]) {
                }
            }

            while (i < cube_id.size()) {
                value = cube_id[i];
                if (!isDownVertex(value)) {
                    all_down    = false;
                }
                if (!isFrontVertex(value)) {
                    all_front   = false;
                }
                if (!isRightVertex(value)) {
                    all_right   = false;
                }
                if (!isBackVertex(value)) {
                    all_back    = false;
                }
                if (!isLeftVertex(value)) {
                    all_left    = false;
                }
                if (!isUpVertex(value)) {
                    all_up      = false;
                }
                i++;
            }

            if (!edgecubes_only) {
                return (all_down || all_front || all_right || all_back || all_left || all_up);
            }
            else {
                return (
                        (all_down   && all_front)   ||
                        (all_down   && all_left)    ||
                        (all_down   && all_back)    ||
                        (all_down   && all_right)   ||

                        (all_front  && all_left)    ||
                        (all_front  && all_right)   ||
                        (all_back   && all_left)    ||
                        (all_back   && all_right)   ||

                        (all_up     && all_front)   ||
                        (all_up     && all_left)    ||
                        (all_up     && all_back)    ||
                        (all_up     && all_right)
                       );

            }
        }

        bool
        isFringeCube(const OtId &n_id, const OtId &cube_id, bool edgecubes_only = false) const
        {
            using namespace CubeCommon;

            uint32_t i;
            bool     all_down, all_front, all_right, all_back, all_left, all_up;
            uint8_t  value;

            all_down    = true;
            all_front   = true;
            all_right   = true;
            all_back    = true;
            all_left    = true;
            all_up      = true;

            /* compare the first depth(n) components of n's and the given cube id */
            for (i = 0; i < n_id.size(); i++) {
                if (n_id[i] != cube_id[i]) {
                }
            }

            while (i < cube_id.size()) {
                value = cube_id[i];
                if (!isDownVertex(value)) {
                    all_down    = false;
                }
                if (!isFrontVertex(value)) {
                    all_front   = false;
                }
                if (!isRightVertex(value)) {
                    all_right   = false;
                }
                if (!isBackVertex(value)) {
                    all_back    = false;
                }
                if (!isLeftVertex(value)) {
                    all_left    = false;
                }
                if (!isUpVertex(value)) {
                    all_up      = false;
                }
                i++;
            }

            if (!edgecubes_only) {
                return (all_down || all_front || all_right || all_back || all_left || all_up);
            }
            else {
                return (
                        (all_down   && all_front)   ||
                        (all_down   && all_left)    ||
                        (all_down   && all_back)    ||
                        (all_down   && all_right)   ||

                        (all_front  && all_left)    ||
                        (all_front  && all_right)   ||
                        (all_back   && all_left)    ||
                        (all_back   && all_right)   ||

                        (all_up     && all_front)   ||
                        (all_up     && all_left)    ||
                        (all_up     && all_back)    ||
                        (all_up     && all_right)
                       );

            }
        }

        inline void
        splitLeaf(OctreeNode<nT> *n)
        {
            OctreeNode<nT>    *children[8];
            uint32_t            i;

            if (n->isLeaf()) {
                /* create new children */
                for (i = 0; i < 8; i++) {
                    /* copy n's (which is the parent) data to children by default. this does not
                     * matter: if its a pointer, we can overwrite it, if not, space for data of type
                     * nT is set aside anyway.. */
                    children[i] = new OctreeNode<nT>(n->data, i);
                    children[i]->setLeaf(true);
                    children[i]->setChildType(i);
                    children[i]->setLastSibling(false);
                }

                /* child 7 is the last sibling */
                children[7]->setLastSibling(true);

                /* assign first child pointer for current node */
                n->first_child = children[0];

                /* do the sibling pointers */
                for (i = 0; i < 7; i++) {
                    children[i]->uplink.next_sibling = children[i + 1];
                }

                /* the last child has the pointer to its parent, which is n */
                children[7]->uplink.parent = n;
                
                /* now n aint no leaf no more o_O */
                n->setLeaf(false);
            }
            else throw("Octree::splitLeaf(). can't split non-leaf.\n"); 
        }

        void
        deleteLeaf(OctreeNode<nT> *n)
        {
            OctreeNode<nT>  *tmp;

            if (n->isLeaf()) {
                if (n != &(this->root)) {
                    tmp = n->getParent();
                    /* n is the first child */
                    if (tmp->first_child == n) {
                        //debug("n is the first child\n");
                        /* if n is the only child, the parent becomes a leaf */
                        if (n->isLastSibling()) {
                            //debug("and the last sibling, i.e. only child.\n");
                            tmp->first_child = NULL;
                            tmp->setLeaf();
                        }
                        /* otherwise, set the first child pointer of the parent to the successor
                         * sibling of n */
                        else {
                            //debug("and not the last sibling\n");
                            tmp->first_child = n->uplink.next_sibling;
                            //debug("first child of parent set to next sibling.\n");
                        }
                    }
                    /* n is not the first child */
                    else {
                        //debug("n is not the first child\n");
                        /* find the predecessor sibling of n */ 
                        for (tmp = tmp->first_child; tmp->uplink.next_sibling != n; tmp = tmp->uplink.next_sibling) {
                            ;
                        }

                        /* if n in the last sibling, set parent pointer of predecessor sibling */
                        if (n->isLastSibling()) {
                            //debug("but the last sibling\n");
                            tmp->setLastSibling();
                            tmp->uplink.parent = n->uplink.parent;
                            //debug("new parent pointer set\n");
                        }
                        /* otherwise set predecessor sibling pointer to skip n */
                        else {
                            //debug("and not the last sibling\n");
                            tmp->uplink.next_sibling = n->uplink.next_sibling;
                            //debug("predecessor pointer modified to skip n\n");
                        }
                    }

                    /* delete n */
                    delete n;
                }
                /* else: its the root, which is a leaf -> only one node tree. we don't delete the
                 * root */
            }
            else throw("Octree::deleteLeaf(): non-leaf node given.\n");
        }

        void
        disconnectChild(OctreeNode<nT> *n)
        {
            OctreeNode<nT>  *tmp;

            /* first, get n's parent */
            tmp = n->getParent();

            /* n is the first child */
            if (tmp->first_child == n) {
                //debug("n is the first child\n");
                /* if n is the only child, the parent becomes a leaf */
                if (n->isLastSibling()) {
                    //debug("and the last sibling, i.e. only child.\n");
                    tmp->first_child = NULL;
                    tmp->setLeaf();
                }
                /* otherwise, set the first child pointer of the parent to the successor
                 * sibling of n */
                else {
                    //debug("and not the last sibling\n");
                    tmp->first_child = n->uplink.next_sibling;
                    //debug("first child of parent set to next sibling.\n");
                }
            }
            /* n is not the first child */
            else {
                //debug("n is not the first child\n");
                /* find the predecessor sibling of n */ 
                for (tmp = tmp->first_child; tmp->uplink.next_sibling != n; tmp = tmp->uplink.next_sibling) {
                    ;
                }

                /* if n in the last sibling, set parent pointer of predecessor sibling */
                if (n->isLastSibling()) {
                    //debug("but the last sibling\n");
                    tmp->setLastSibling();
                    tmp->uplink.parent = n->uplink.parent;
                    //debug("new parent pointer set\n");
                }
                /* otherwise set predecessor sibling pointer to skip n */
                else {
                    //debug("and not the last sibling\n");
                    tmp->uplink.next_sibling = n->uplink.next_sibling;
                    //debug("predecessor pointer modified to skip n\n");
                }
            }
        }
};

template <typename tT, typename nT, typename R>
void
Octree<tT, nT, R>::octree_generate_leaf_list_recursive(OctreeNode<nT> *n, std::list<OctreeNode<nT> *> &l)
{
    if (n->isLeaf()) {
        l.push_back(n);
    }
    else {
        for (int i = 0; i < 8; i++) {
            if (n->children[i]) {
                octree_generate_leaf_list_recursive(n->children[i], l);
            }
        }
    }
}

template <typename tT, typename nT, typename R>
void
Octree<tT, nT, R>::octree_generate_leaf_vector_recursive(OctreeNode<nT> *n, std::vector<OctreeNode<nT> *> &v)
{
    if (n->isLeaf()) {
        v.push_back(n);
    }
    else {
        for (int i = 0; i < 8; i++) {
            if (n->children[i]) {
                octree_generate_leaf_vector_recursive(n->children[i], v);
            }
        }
    }
}

template <typename tT, typename nT, typename R>
void
Octree<tT, nT, R>::octree_free_recursive(OctreeNode<nT> *n)
{
    int i;
    OctreeNode<nT> *children[8];

    if (!n->isLeaf()) {
        /* first, get all children while pointers are still intact */
        for (i = 0; i < 8; i++) {
            children[i] = n->getChild(i);
        }
        
        /* now free them if they exist */
        for (i = 0; i < 8; i++) {
            if (children[i]) {
                octree_free_recursive(children[i]);
                delete children[i];
            }
        }
    }
}

#endif
