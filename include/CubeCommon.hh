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

#ifndef CUBE_COMMON_H
#define CUBE_COMMON_H

#include "Vec3.hh"

namespace CubeCommon {

    enum cube_parts {
        /* 8 vertices of the cube 0 - 4, 5 - 7 */
        V_LEFT_FRONT_DOWN   = 0,
        V_RIGHT_FRONT_DOWN  = 1,
        V_RIGHT_BACK_DOWN   = 2,
        V_LEFT_BACK_DOWN    = 3,
        V_LEFT_FRONT_UP     = 4,
        V_RIGHT_FRONT_UP    = 5,
        V_RIGHT_BACK_UP     = 6,
        V_LEFT_BACK_UP      = 7,

        /* 12 edges of the cube, each one has two representations */

        /* four edges on the "down" cube */
        E_DOWN_FRONT        = 8,
        E_FRONT_DOWN        = 8,

        E_DOWN_RIGHT        = 9,
        E_RIGHT_DOWN        = 9,

        E_DOWN_BACK         = 10,
        E_BACK_DOWN         = 10,

        E_DOWN_LEFT         = 11,
        E_LEFT_DOWN         = 11,

        /* four edges connecting "lower" square with "upper" square */
        E_FRONT_LEFT        = 12,
        E_LEFT_FRONT        = 12,

        E_FRONT_RIGHT       = 13,
        E_RIGHT_FRONT       = 13,

        E_BACK_RIGHT        = 14,
        E_RIGHT_BACK        = 14,

        E_BACK_LEFT         = 15,
        E_LEFT_BACK         = 15,

        /* four edges on the "up" cube */
        E_UP_FRONT          = 16,
        E_FRONT_UP          = 16,

        E_UP_RIGHT          = 17,
        E_RIGHT_UP          = 17,

        E_UP_BACK           = 18,
        E_BACK_UP           = 18,

        E_UP_LEFT           = 19,
        E_LEFT_UP           = 19,

        /* 6 faces of the cube */
        F_DOWN              = 20,
        F_FRONT             = 21,
        F_RIGHT             = 22,
        F_BACK              = 23,
        F_LEFT              = 24,
        F_UP                = 25
    };

    inline bool
    isDownVertex(uint32_t vtype)
    {
        return (vtype == V_LEFT_FRONT_DOWN || vtype == V_RIGHT_FRONT_DOWN || vtype == V_RIGHT_BACK_DOWN || vtype == V_LEFT_BACK_DOWN);
    }

    inline bool
    isUpVertex(uint32_t vtype)
    {
        return (vtype == V_LEFT_FRONT_UP || vtype == V_RIGHT_FRONT_UP || vtype == V_RIGHT_BACK_UP || vtype == V_LEFT_BACK_UP);
    }

    inline bool
    isRightVertex(uint32_t vtype)
    {
        return (vtype == V_RIGHT_FRONT_DOWN || vtype == V_RIGHT_BACK_DOWN || vtype == V_RIGHT_FRONT_UP || vtype == V_RIGHT_BACK_UP);
    }

    inline bool
    isLeftVertex(uint32_t vtype)
    {
        return (vtype == V_LEFT_FRONT_DOWN || vtype == V_LEFT_BACK_DOWN || vtype == V_LEFT_FRONT_UP || vtype == V_LEFT_BACK_UP);
    }

    inline bool
    isFrontVertex(uint32_t vtype)
    {
        return (vtype == V_LEFT_FRONT_DOWN || vtype == V_RIGHT_FRONT_DOWN || vtype == V_LEFT_FRONT_UP || vtype == V_RIGHT_FRONT_UP);
    }

    inline bool
    isBackVertex(uint32_t vtype)
    {
        return (vtype == V_LEFT_BACK_DOWN || vtype == V_RIGHT_BACK_DOWN || vtype == V_LEFT_BACK_UP || vtype == V_RIGHT_BACK_UP);
    }

    inline bool
    isFaceVertex(uint8_t vtype, uint8_t face_type) {
        switch (face_type) {
            case F_DOWN:
                return isDownVertex(vtype);

            case F_FRONT:
                return isFrontVertex(vtype);
                
            case F_RIGHT:
                return isRightVertex(vtype);

            case F_BACK:
                return isBackVertex(vtype);

            case F_LEFT:
                return isLeftVertex(vtype);

            case F_UP:
                return isUpVertex(vtype);

            default:
                throw("invalid face type given.\n");
        }
    }

    inline uint32_t
    mirrorVertexLeftRight(uint32_t vtype)
    {
        switch (vtype) {
            case V_LEFT_FRONT_DOWN:
                return V_RIGHT_FRONT_DOWN;

            case V_RIGHT_FRONT_DOWN:
                return V_LEFT_FRONT_DOWN;

            case V_RIGHT_BACK_DOWN:
                return V_LEFT_BACK_DOWN;

            case V_LEFT_BACK_DOWN:
                return V_RIGHT_BACK_DOWN;

            case V_LEFT_FRONT_UP:
                return V_RIGHT_FRONT_UP;

            case V_RIGHT_FRONT_UP:
                return V_LEFT_FRONT_UP;

            case V_RIGHT_BACK_UP:
                return V_LEFT_BACK_UP;

            case V_LEFT_BACK_UP:
                return V_RIGHT_BACK_UP;

            default:
                throw("mirrorVertexLeftRight(): invalid vertex type given.\n");

        }
    }

    inline uint32_t
    mirrorVertexDownUp(uint32_t vtype)
    {
        switch (vtype) {
            case V_LEFT_FRONT_DOWN:
                return V_LEFT_FRONT_UP;

            case V_RIGHT_FRONT_DOWN:
                return V_RIGHT_FRONT_UP;

            case V_RIGHT_BACK_DOWN:
                return V_RIGHT_BACK_UP;

            case V_LEFT_BACK_DOWN:
                return V_LEFT_BACK_UP;

            case V_LEFT_FRONT_UP:
                return V_LEFT_FRONT_DOWN;

            case V_RIGHT_FRONT_UP:
                return V_RIGHT_FRONT_DOWN;

            case V_RIGHT_BACK_UP:
                return V_RIGHT_BACK_DOWN;

            case V_LEFT_BACK_UP:
                return V_LEFT_BACK_DOWN;
            
            default:
                throw("mirrorVertexDownUp(): invalid vertex type given.\n");

        }
    }

    inline uint32_t
    mirrorVertexFrontBack(uint32_t vtype)
    {
        switch (vtype) {
            case V_LEFT_FRONT_DOWN:
                return V_LEFT_BACK_DOWN;

            case V_RIGHT_FRONT_DOWN:
                return V_RIGHT_BACK_DOWN;

            case V_RIGHT_BACK_DOWN:
                return V_RIGHT_FRONT_DOWN;

            case V_LEFT_BACK_DOWN:
                return V_LEFT_FRONT_DOWN;

            case V_LEFT_FRONT_UP:
                return V_LEFT_BACK_UP;

            case V_RIGHT_FRONT_UP:
                return V_RIGHT_BACK_UP;

            case V_RIGHT_BACK_UP:
                return V_RIGHT_FRONT_UP;

            case V_LEFT_BACK_UP:
                return V_LEFT_FRONT_UP;
            
            default:
                throw("mirrorVertexFrontBack(): invalid vertex type given.\n");

        }
    }

    inline uint8_t
    mirrorVertexFace(uint8_t vtype, uint8_t mirror_type) {
        switch (mirror_type) {
            case F_DOWN:
            case F_UP:
                return mirrorVertexDownUp(vtype);
                break;

            case F_RIGHT:
            case F_LEFT:
                return mirrorVertexLeftRight(vtype);
                break;

            case F_FRONT:
            case F_BACK:
                return mirrorVertexFrontBack(vtype);
                break;

            default:
                throw("mirrorVertexFace(): invalid mirror type (== face type) given.\n");
        }
    }

    inline uint8_t
    mirrorVertexCentre(uint8_t vtype)
    {
        switch (vtype) {
            case V_LEFT_FRONT_DOWN:
                return V_RIGHT_BACK_UP;

            case V_RIGHT_FRONT_DOWN:
                return V_LEFT_BACK_UP;

            case V_RIGHT_BACK_DOWN:
                return V_LEFT_FRONT_UP;

            case V_LEFT_BACK_DOWN:
                return V_RIGHT_FRONT_UP;

            case V_LEFT_FRONT_UP:
                return V_RIGHT_BACK_DOWN;

            case V_RIGHT_FRONT_UP:
                return V_LEFT_BACK_DOWN;

            case V_RIGHT_BACK_UP:
                return V_LEFT_FRONT_DOWN;

            case V_LEFT_BACK_UP:
                return V_RIGHT_FRONT_DOWN;

            default:
                throw("mirrorVertexCentre(): invalid vertex type given.\n");
        }
    }

    inline uint8_t
    mirrorEdgeCentre(uint8_t etype)
    {
        switch(etype) {
             /* four edges on the "down" cube */
            //case E_DOWN_FRONT:
            case E_FRONT_DOWN:
                return E_BACK_UP;

            //case E_DOWN_RIGHT:
            case E_RIGHT_DOWN:
                return E_LEFT_UP;

            //case E_DOWN_BACK:
            case E_BACK_DOWN:
                return E_FRONT_UP;

            //case E_DOWN_LEFT:
            case E_LEFT_DOWN:
                return E_RIGHT_UP;

             /* four edges connecting "lower" square with "upper" square */
            //case E_FRONT_LEFT:
            case E_LEFT_FRONT:
                return E_RIGHT_BACK;

            //case E_FRONT_RIGHT:
            case E_RIGHT_FRONT:
                return E_LEFT_BACK;

            //case E_BACK_RIGHT:
            case E_RIGHT_BACK:
                return E_LEFT_FRONT;

            //case E_BACK_LEFT:
            case E_LEFT_BACK:
                return E_RIGHT_FRONT;

            /* four edges on the "up" cube */
            //case E_UP_FRONT:
            case E_FRONT_UP:
                return E_BACK_DOWN;

            //case E_UP_RIGHT:
            case E_RIGHT_UP:
                return E_LEFT_DOWN;

            //case E_UP_BACK:
            case E_BACK_UP:
                return E_FRONT_DOWN;

            //case E_UP_LEFT:
            case E_LEFT_UP:
                return E_RIGHT_DOWN;

            default:
                throw("mirrorEdgeCentre(): invalid edge type giben.\n");
        }
    }

    /* get vertices, edges and faces from number [0..7], [0..12] and [0..6] */
#define getVertexEnumValue(vertex_num)  (vertex_num)
#define getVertexNumber(vertex_num)     (vertex_num)
#define getEdgeEnumValue(edge_num)      (8 + edge_num)
#define getEdgeNumber(edge_enum)        (edge_enum - 8)
#define getFaceEnumValue(face_num)      (20 + face_num)
    /* is value a vertex / edge / face enum value? */
#define isVertexEnum(v)                 (v <= 7)
#define isEdgeEnum(v)                   (v >= 8 && v <= 19)
#define isFaceEnum(v)                   (v >= 20 && v <= 25)

    /* NOTE: duplicate case values not allowed, so the second representation for earch edge is commented
     * out but still left there for documentation */
    /* NOTE: atype is always the smaller value, btype always the larger value, vertex types are numbered
     * from 0..7 */
    inline void
    getEdgeEndpointTypes(uint32_t etype, uint8_t &atype, uint8_t &btype)
    {
        switch (etype) {
            case E_DOWN_FRONT:
            //case E_FRONT_DOWN:
                atype = V_LEFT_FRONT_DOWN;
                btype = V_RIGHT_FRONT_DOWN;
                return;
            
            case E_DOWN_RIGHT:
            //case E_RIGHT_DOWN:
                atype = V_RIGHT_FRONT_DOWN;
                btype = V_RIGHT_BACK_DOWN;
                return;

            case E_DOWN_BACK:
            //case E_BACK_DOWN:
                atype = V_RIGHT_BACK_DOWN;
                btype = V_LEFT_BACK_DOWN;
                return;

            case E_DOWN_LEFT:
            //case E_LEFT_DOWN:
                atype = V_LEFT_FRONT_DOWN;
                btype = V_LEFT_BACK_DOWN;
                return;

                /* four edges connecting "lower" square with "upper" square */
            case E_FRONT_LEFT:
            //case E_LEFT_FRONT:
                atype = V_LEFT_FRONT_DOWN;
                btype = V_LEFT_FRONT_UP;
                return;

            case E_FRONT_RIGHT:
            //case E_RIGHT_FRONT:
                atype = V_RIGHT_FRONT_DOWN;
                btype = V_RIGHT_FRONT_UP;
                return;

            case E_BACK_RIGHT:
            //case E_RIGHT_BACK:
                atype = V_RIGHT_BACK_DOWN;
                btype = V_RIGHT_BACK_UP;
                return;

            case E_BACK_LEFT:
            //case E_LEFT_BACK:
                atype = V_LEFT_BACK_DOWN;
                btype = V_LEFT_BACK_UP;
                return;

                /* four edges on the "up" cube */
            case E_UP_FRONT:
            //case E_FRONT_UP:
                atype = V_LEFT_FRONT_UP;
                btype = V_RIGHT_FRONT_UP;
                return;

            case E_UP_RIGHT:
            //case E_RIGHT_UP:
                atype = V_RIGHT_FRONT_UP;
                btype = V_RIGHT_BACK_UP;
                return;

            case E_UP_BACK:
            //case E_BACK_UP:
                atype = V_RIGHT_BACK_UP;
                btype = V_LEFT_BACK_UP;
                return;

            case E_UP_LEFT:
            //case E_LEFT_UP:
                atype = V_LEFT_FRONT_UP;
                btype = V_LEFT_BACK_UP;
                return;

            default:
                throw("getEdgeEndpointTypes(): case default shall never be reached.\n");
        }
    }

    /* get the "directions" to the three neighbours that contain edge e of type etype in the cube.
     * example: E_DOWN_FRONT is contained in the
     *
     * down, front and down->front (==front->down)
     *
     * neighbour, where down->front means the front neighbour of the down neighbour ("diagonal"
     * neighbour). we need only two neighbour directions, since the third is implicitly given (in the
     * example: down and front, the third is down->front / front->down) */
    inline void
    getEdgeNeighbourDirections(uint8_t etype, uint8_t &nb1_type, uint8_t &nb2_type)
    {
        switch (etype) {
            case E_DOWN_FRONT:
            //case E_FRONT_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_FRONT;
                return;
            
            case E_DOWN_RIGHT:
            //case E_RIGHT_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_RIGHT;
                return;

            case E_DOWN_BACK:
            //case E_BACK_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_BACK;
                return;

            case E_DOWN_LEFT:
            //case E_LEFT_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_LEFT;
                return;

                /* four edges connecting "lower" square with "upper" square */
            case E_FRONT_LEFT:
            //case E_LEFT_FRONT:
                nb1_type = F_FRONT;
                nb2_type = F_LEFT;
                return;

            case E_FRONT_RIGHT:
            //case E_RIGHT_FRONT:
                nb1_type = F_FRONT;
                nb2_type = F_RIGHT;
                return;

            case E_RIGHT_BACK:
            //case E_BACK_RIGHT:
                nb1_type = F_RIGHT;
                nb2_type = F_BACK;
                return;

            case E_BACK_LEFT:
            //case E_LEFT_BACK:
                nb1_type = F_BACK;
                nb2_type = F_LEFT;
                return;

                /* four edges on the "up" cube */
            case E_FRONT_UP:
            //case E_UP_FRONT:
                nb1_type = F_FRONT;
                nb2_type = F_UP;
                return;

            case E_RIGHT_UP:
            //case E_UP_RIGHT:
                nb1_type = F_RIGHT;
                nb2_type = F_UP;
                return;

            case E_BACK_UP:
            //case E_UP_BACK:
                nb1_type = F_BACK;
                nb2_type = F_UP;
                return;

            case E_LEFT_UP:
            //case E_UP_LEFT:
                nb1_type = F_LEFT;
                nb2_type = F_UP;
                return;

            default:
                throw("getEdgeNeighbhourDirections(): case default shall never be reached.\n");
        }
    }
        
    /* get "directions" to neighbour identified by the given vertex. example: V_LEFT_FRONT_DOWN ->
     * direcitons: left, front, down */
    inline void
    getVertexNeighbourDirections(uint8_t vtype, uint8_t &nb1_type, uint8_t &nb2_type, uint8_t &nb3_type)
    {
        switch (vtype) {
            case V_LEFT_FRONT_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_FRONT;
                nb3_type = F_LEFT;
                return;

            case V_RIGHT_FRONT_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_FRONT;
                nb3_type = F_RIGHT;
                return;

            case V_RIGHT_BACK_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_RIGHT;
                nb3_type = F_BACK;
                return;

            case V_LEFT_BACK_DOWN:
                nb1_type = F_DOWN;
                nb2_type = F_BACK;
                nb3_type = F_LEFT;
                return;

            case V_LEFT_FRONT_UP:
                nb1_type = F_FRONT;
                nb2_type = F_LEFT;
                nb3_type = F_UP;
                return;

            case V_RIGHT_FRONT_UP:
                nb1_type = F_FRONT;
                nb2_type = F_RIGHT;
                nb3_type = F_UP;
                return;

            case V_RIGHT_BACK_UP:
                nb1_type = F_RIGHT;
                nb2_type = F_BACK;
                nb3_type = F_UP;
                return;

            case V_LEFT_BACK_UP:
                nb1_type = F_BACK;
                nb2_type = F_LEFT;
                nb3_type = F_UP;
                return;
            
            default:
                throw("mirrorVertexFrontBack(): invalid vertex type given.\n");
        }
    }

    /* get all neighbours (encoded as enum values) containing a given vertex, specified as enum value */
    inline void
    getNeighbourTypesContainingVertex(uint8_t vtype, uint8_t *nb_types)
    {
        switch (vtype) {
            case V_LEFT_FRONT_DOWN:
                nb_types[0] = F_LEFT;
                nb_types[1] = F_FRONT;
                nb_types[2] = F_DOWN;
                nb_types[3] = E_LEFT_FRONT;
                nb_types[4] = E_LEFT_DOWN;
                nb_types[5] = E_FRONT_DOWN;
                nb_types[6] = V_LEFT_FRONT_DOWN;
                return;

            case V_RIGHT_FRONT_DOWN:
                nb_types[0] = F_RIGHT;
                nb_types[1] = F_FRONT;
                nb_types[2] = F_DOWN;
                nb_types[3] = E_RIGHT_FRONT;
                nb_types[4] = E_RIGHT_DOWN;
                nb_types[5] = E_FRONT_DOWN;
                nb_types[6] = V_RIGHT_FRONT_DOWN;
                return;

            case V_RIGHT_BACK_DOWN:
                nb_types[0] = F_RIGHT;
                nb_types[1] = F_BACK;
                nb_types[2] = F_DOWN;
                nb_types[3] = E_RIGHT_BACK;
                nb_types[4] = E_RIGHT_DOWN;
                nb_types[5] = E_BACK_DOWN;
                nb_types[6] = V_RIGHT_BACK_DOWN;
                return;

            case V_LEFT_BACK_DOWN:
                nb_types[0] = F_LEFT;
                nb_types[1] = F_BACK;
                nb_types[2] = F_DOWN;
                nb_types[3] = E_LEFT_BACK;
                nb_types[4] = E_LEFT_DOWN;
                nb_types[5] = E_BACK_DOWN;
                nb_types[6] = V_LEFT_BACK_DOWN;
                return;

            case V_LEFT_FRONT_UP:
                nb_types[0] = F_LEFT;
                nb_types[1] = F_FRONT;
                nb_types[2] = F_UP;
                nb_types[3] = E_LEFT_FRONT;
                nb_types[4] = E_LEFT_UP;
                nb_types[5] = E_FRONT_UP;
                nb_types[6] = V_LEFT_FRONT_UP;
                return;

            case V_RIGHT_FRONT_UP:
                nb_types[0] = F_RIGHT;
                nb_types[1] = F_FRONT;
                nb_types[2] = F_UP;
                nb_types[3] = E_RIGHT_FRONT;
                nb_types[4] = E_RIGHT_UP;
                nb_types[5] = E_FRONT_UP;
                nb_types[6] = V_RIGHT_FRONT_UP;
                return;

            case V_RIGHT_BACK_UP:
                nb_types[0] = F_RIGHT;
                nb_types[1] = F_BACK;
                nb_types[2] = F_UP;
                nb_types[3] = E_RIGHT_BACK;
                nb_types[4] = E_RIGHT_UP;
                nb_types[5] = E_BACK_UP;
                nb_types[6] = V_RIGHT_BACK_UP;
                return;

            case V_LEFT_BACK_UP:
                nb_types[0] = F_LEFT;
                nb_types[1] = F_BACK;
                nb_types[2] = F_UP;
                nb_types[3] = E_LEFT_BACK;
                nb_types[4] = E_LEFT_UP;
                nb_types[5] = E_BACK_UP;
                nb_types[6] = V_LEFT_BACK_UP;
                return;
            
            default:
                throw("getVertexNeighbours(): invalid vertex type given.\n");
        }
    }

    /* get all neighbours (encoded as enum values) containing a given edge, specified as enum value */
    inline void
    getNeighbourTypesContainingEdge(uint8_t etype, uint8_t *nb_types)
    {
        switch (etype) {
            case E_DOWN_FRONT:
            //case E_FRONT_DOWN:
                nb_types[0] = F_DOWN;
                nb_types[1] = F_FRONT;
                nb_types[2] = E_DOWN_FRONT;
                return;
            
            case E_DOWN_RIGHT:
            //case E_RIGHT_DOWN:
                nb_types[0] = F_DOWN;
                nb_types[1] = F_RIGHT;
                nb_types[2] = E_DOWN_RIGHT;
                return;

            case E_DOWN_BACK:
            //case E_BACK_DOWN:
                nb_types[0] = F_DOWN;
                nb_types[1] = F_BACK;
                nb_types[2] = E_DOWN_BACK;
                return;

            case E_DOWN_LEFT:
            //case E_LEFT_DOWN:
                nb_types[0] = F_DOWN;
                nb_types[1] = F_LEFT;
                nb_types[2] = E_DOWN_LEFT;
                return;

                /* four edges connecting "lower" square with "upper" square */
            case E_FRONT_LEFT:
            //case E_LEFT_FRONT:
                nb_types[0] = F_FRONT;
                nb_types[1] = F_LEFT;
                nb_types[2] = E_FRONT_LEFT;
                return;

            case E_FRONT_RIGHT:
            //case E_RIGHT_FRONT:
                nb_types[0] = F_FRONT;
                nb_types[1] = F_RIGHT;
                nb_types[2] = E_FRONT_RIGHT;
                return;

            case E_RIGHT_BACK:
            //case E_BACK_RIGHT:
                nb_types[0] = F_RIGHT;
                nb_types[1] = F_BACK;
                nb_types[2] = E_RIGHT_BACK;
                return;

            case E_BACK_LEFT:
            //case E_LEFT_BACK:
                nb_types[0] = F_BACK;
                nb_types[1] = F_LEFT;
                nb_types[2] = E_BACK_LEFT;
                return;

                /* four edges on the "up" cube */
            case E_FRONT_UP:
            //case E_UP_FRONT:
                nb_types[0] = F_FRONT;
                nb_types[1] = F_UP;
                nb_types[2] = E_FRONT_UP;
                return;

            case E_RIGHT_UP:
            //case E_UP_RIGHT:
                nb_types[0] = F_RIGHT;
                nb_types[1] = F_UP;
                nb_types[2] = E_RIGHT_UP;
                return;

            case E_BACK_UP:
            //case E_UP_BACK:
                nb_types[0] = F_BACK;
                nb_types[1] = F_UP;
                nb_types[2] = E_BACK_UP; 
                return;

            case E_LEFT_UP:
            //case E_UP_LEFT:
                nb_types[0] = F_LEFT;
                nb_types[1] = F_UP;
                nb_types[2] = E_LEFT_UP; 
                return;

            default:
                throw("getEdgeNeighbhourDirections(): case default shall never be reached.\n");
        }   
    }

    template <typename R>
    inline Vec3<R>
    getAACubeCorner(
        Vec3<R> const  &pos,
        R const        &l,
        uint8_t         vtype)
    {
        Vec3<R> corner = pos;

        /* compute child i's cube vertex */
        switch (vtype) {
            case V_LEFT_FRONT_DOWN:
                break;

            case V_RIGHT_FRONT_DOWN:
                corner[0] += l;
                break;

            case V_RIGHT_BACK_DOWN:
                corner[0] += l;
                corner[1] += l;
                break;

            case V_LEFT_BACK_DOWN:
                corner[1] += l;
                break;

            case V_LEFT_FRONT_UP:
                corner[2] += l;
                break;

            case V_RIGHT_FRONT_UP:
                corner[0] += l;
                corner[2] += l;
                break;

            case V_RIGHT_BACK_UP:
                corner[0] += l;
                corner[1] += l;
                corner[2] += l;
                break;

            case V_LEFT_BACK_UP:
                corner[1] += l;
                corner[2] += l;
                break;
        }

        return corner;
    }

    template <typename R>
    inline void
    getAACubeCorners(
        Vec3<R> const  &pos,
        R const        &l,
        Vec3<R>       (&corners)[8])
    {
        corners[0]      = pos;

        corners[1]      = pos;
        corners[1][0]  += l;

        corners[2]      = corners[1];
        corners[2][1]  += l;

        corners[3]      = pos;
        corners[3][1]  += l;

        corners[4]      = pos;
        corners[4][2]  += l;

        corners[5]      = corners[4];
        corners[5][0]  += l;

        corners[6]      = corners[5];
        corners[6][1]  += l;

        corners[7]      = corners[4];
        corners[7][1]  += l;
    }

}

#endif
