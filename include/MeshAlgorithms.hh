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

#ifndef MESH_ALGORITHMS_H
#define MESH_ALGORITHMS_H


/* exception classes for Red-Blue Union algorithms */
namespace RedBlue_ExCodes {
    enum RedBlue_ExceptionCodes_enum {
        RB_INTERNAL_LOGIC_ERROR,
        RB_DISJOINT,
        RB_COMPLEX_EDGES,
        RB_NUMERICAL_EDGE_CASE,
        RB_TRIANGULATION,
        RB_ISECPOLY_NUM,
        RB_AFFECTED_CIRCLE_TRIVIAL
    };
};

struct RedBlue_Ex : public std::runtime_error {
    const uint32_t      type;
    const std::string   msg;
    bool                R_intact;
    bool                B_intact;

    RedBlue_Ex(
        const uint32_t     &type,
        const std::string  &msg,
        bool                R_intact,
        bool                B_intact)
            : std::runtime_error(msg), type(type), msg(msg), R_intact(R_intact), B_intact(B_intact)
    {
    }
};

/* struct to store information about intersecting edges during the RedBlue algorithm in all flavours. primarily needed
 * for exceptions, since this info is bundled and thrown if certain invariants do not hold. */
template <typename TR>
struct RedBlue_EdgeIsecInfo {
    /* edge red or blue */
    bool                    red;

    /* ids of the edge endpoints in the order that the intersection was performed, which is
     * important for the lambda value below */
    uint32_t                u_id, v_id;

    /* ids of the other mesh's intersected faces */
    std::vector<uint32_t>   isec_face_ids;

    /* corresponding lambda values, i.e. intersection points are u + (v-u)*lambda(k) for
     * in-range k */
    std::vector<TR>         edge_lambdas;

    RedBlue_EdgeIsecInfo(
        bool                            red,
        uint32_t                        u_id,
        uint32_t                        v_id,
        std::vector<uint32_t> const    &isec_face_ids,
        std::vector<TR> const          &edge_lambdas)
    {
        this->red           = red;
        this->u_id          = u_id;
        this->v_id          = v_id;
        this->isec_face_ids = isec_face_ids;
        if (!edge_lambdas.empty()) {
            this->edge_lambdas  = edge_lambdas;
        }
        else {
            throw("RedBlue_Ex_ComplexEdges(): no lambda values supplied for complex edge. internal logic error.");
        }
    }
};

struct RedBlue_Ex_InternalLogic : public RedBlue_Ex {
    RedBlue_Ex_InternalLogic(std::string msg)
        /* default to setting R and B to not-intact. internal logic errors are fatal anyway. */
        : RedBlue_Ex(RedBlue_ExCodes::RB_INTERNAL_LOGIC_ERROR, msg, false, false)
    {
    }
};

struct RedBlue_Ex_Disjoint : public RedBlue_Ex {
    RedBlue_Ex_Disjoint(std::string msg)
        : RedBlue_Ex(RedBlue_ExCodes::RB_DISJOINT, msg, true, true)
    {
    }
};

template <typename TR>
struct RedBlue_Ex_ComplexEdges : public RedBlue_Ex {
    std::list<RedBlue_EdgeIsecInfo<TR>> edge_isec_info;

    RedBlue_Ex_ComplexEdges(
            std::string                                 msg,
            std::list<RedBlue_EdgeIsecInfo<TR>> const  &einfo_vec)
                : RedBlue_Ex(RedBlue_ExCodes::RB_COMPLEX_EDGES, msg, true, true)
    {
        this->edge_isec_info = einfo_vec;
    }
};

struct RedBlue_Ex_NumericalEdgeCase : public RedBlue_Ex {
    RedBlue_Ex_NumericalEdgeCase(
        std::string const  &msg,
        bool                R_intact,
        bool                B_intact)
            : RedBlue_Ex(RedBlue_ExCodes::RB_NUMERICAL_EDGE_CASE, msg, R_intact, B_intact)
    {
    }
};

template <typename TR>
struct RedBlue_Ex_Triangulation : public RedBlue_Ex {
    RedBlue_Ex_Triangulation(
        std::string const  &msg,
        bool                R_intact,
        bool                B_intact)
            : RedBlue_Ex(RedBlue_ExCodes::RB_TRIANGULATION, msg, R_intact, B_intact)
    {
    }
};

struct RedBlue_Ex_NumIsecPoly : public RedBlue_Ex {
    RedBlue_Ex_NumIsecPoly(std::string msg)
        : RedBlue_Ex(RedBlue_ExCodes::RB_ISECPOLY_NUM, msg, true, true)
    {
    }
};

struct RedBlue_Ex_AffectedCircleTrivial : public RedBlue_Ex {
    bool        red;
    uint32_t    face_id;

    RedBlue_Ex_AffectedCircleTrivial(
        std::string msg,
        bool        red,
        uint32_t    face_id)
            : RedBlue_Ex(RedBlue_ExCodes::RB_AFFECTED_CIRCLE_TRIVIAL, msg, true, true), red(red), face_id(face_id)
    {
    }
};

namespace MeshAlg {
    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    generateOctSphere(
        Vec3<R>                 c,
        R const                &r,
        uint32_t                tessellation_depth,
        Mesh<Tm, Tv, Tf, R>    &S);

    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    generateIcoSphere(
        Vec3<R>                 c,
        R const                &r,
        uint32_t                tessellation_depth,
        Mesh<Tm, Tv, Tf, R>    &S);


    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    appendHalfSphereToCanalSurfaceMesh(
        Mesh<Tm, Tv, Tf, R>                            &M,
        Vec3<R>                                         render_vector,
        Vec3<R>                                         start,
        R const                                        &radius,
        Vec3<R>                                         direction,
        uint32_t                                        nphisegments,
        R const                                        &phi_offset,
        std::vector<
                typename Mesh<Tm, Tv, Tf, R>
                    ::vertex_iterator
            >                                           start_circle_its,
        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator   closing_vertex_it);

    /* for two given meshes X and Y, get pairs of potentially intersecting edges (from X, Y) / faces
     * (from Y, X) using a modified Octree-like construction and traversal.  an Octree is implicitly
     * constructed, but not stored, since it is not needed. instead, two lists given per reference
     * are filled with all (unique) result pairs in the process. */
    template <typename Tm, typename Tv, typename Tf, typename R>
    void 
    getPotentiallyIntersectingEdgeFacePairs(
        Mesh<Tm, Tv, Tf, R>                        &X,
        Mesh<Tm, Tv, Tf, R>                        &Y,
        std::list<
                std::pair<
                    std::pair<
                        typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                        typename Mesh<Tm, Tv, Tf, R>::Vertex *
                    >,
                    typename Mesh<Tm, Tv, Tf, R>::Face *
                >
            >                                      &X_edges_Y_faces_candidates,
        std::list<
                std::pair<
                    typename std::pair<
                        typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                        typename Mesh<Tm, Tv, Tf, R>::Vertex *
                    >,
                    typename Mesh<Tm, Tv, Tf, R>::Face *
                >
            >                                      &Y_edges_X_faces_candidates,
        uint32_t                                    max_components      = 128,
        uint32_t                                    max_recursion_depth = 7);

    /* red blue union algorithm and specializations for set operations */
    template <typename Tm, typename Tv, typename Tf, typename TR>
    void 
    RedBlueAlgorithm(
        Mesh<Tm, Tv, Tf, TR>                               &R,
        Mesh<Tm, Tv, Tf, TR>                               &B,
        const bool                                         &keep_red_outside_part,
        const bool                                         &keep_blue_outside_part,
        std::vector<
                typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
            >                                              *blue_update_its = NULL);

    template <typename Tm, typename Tv, typename Tf, typename TR>
    void
    RedBlueUnion(
        Mesh<Tm, Tv, Tf, TR>                               &R,
        Mesh<Tm, Tv, Tf, TR>                               &B,
        std::vector<
                typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
            >                                              *blue_update_its = NULL);
    
    template <typename Tm, typename Tv, typename Tf, typename TR>
    void
    RedBlueRedMinusBlue(
        Mesh<Tm, Tv, Tf, TR>                               &R,
        Mesh<Tm, Tv, Tf, TR>                               &B,
        std::vector<
                typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
            >                                              *blue_update_its = NULL);

    template <typename Tm, typename Tv, typename Tf, typename TR>
    void
    RedBlueIntersection(
        Mesh<Tm, Tv, Tf, TR>                               &R,
        Mesh<Tm, Tv, Tf, TR>                               &B,
        std::vector<
                typename Mesh<Tm, Tv, Tf, TR>::vertex_iterator
            >                                              *blue_update_its = NULL);

    
    /* greedy edge collapse post-processing */
    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    greedyEdgeCollapsePostProcessing(
        Mesh<Tm, Tv, Tf, R>    &M,
        R const                &alpha   = 1.75,
        R const                &lambda  = 0.125,
        R const                &mu      = 0.5,
        uint32_t                d       = 15);

    /* smoothing algorithms */
    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    simpleLaplacianSmoothing(
        Mesh<Tm, Tv, Tf, R>    &M,
        R const                &lambda,
        uint32_t                maxiter);

    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    HCLaplacianSmoothing(
        Mesh<Tm, Tv, Tf, R>    &M,
        R const                &alpha   = 0.4,
        R const                &beta    = 0.7,
        uint32_t                maxiter = 100);

    /* functions to allow partial flushing of a mesh to an obj file. NOTE: this does not provide paging functionality
     * that can be applied transparently by the user, i.e.: if a part of a mesh has been dumped, it is no longer part of
     * the mesh it came from and access to the respective mesh components is impossible as long as these have not been
     * explicitly reloaded from the file. 
     *
     * these functions are mainly used by the inductive cell meshing algorithm in NLM_CellNetwork: if the partial cell
     * mesh grows large, the merging times (typically) increase as ~O(nlog(n)), where n is the number of mesh
     * components. however, large parts of the mesh reside in memory and will never be reused after they have been
     * created, since no neurite path meshes remain to be merged in those parts. these can be identified by getting the
     * bounding boxes of all neurite canal segments that remain to be merged. all mesh components that lie properly
     * outside of the union of these bounding boxes can safely be flushed to disk. 
     *
     * in order not to meddle with the internal structure (e.g. vertex / face numbering), the following methods have
     * been designed to take care of the vertex indexing / numbering issue. obj files vertex lines semantics do not
     * specify vertex indices, but number the vertices (represented as single lines each) consecutively in order of
     * appearance. */
    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    partialFlushToObjFile(
        Mesh<Tm, Tv, Tf, R>                                        &M,
        std::pair<FILE **, std::string> const                      &obj_file_info,
        std::list<typename Mesh<Tm, Tv, Tf, R>::Face *>            &face_list,
        std::list<
                std::pair<
                    typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                    uint32_t
                >
            >                                                      &in_boundary_vertices,
        uint32_t const                                             &in_last_flush_vertex_id,
        std::list<
                std::pair<
                    typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                    uint32_t          
                >
            >                                                      &out_boundary_vertices,
        uint32_t                                                   &out_last_flush_vertex_id);

    /* class storing information about the flushing process to ease use of the above function. used in conjunction with
     * the wrapper overload of partialFlushToObjFile(..) below. */
    template<typename Tm, typename Tv, typename Tf, typename R>
    class MeshObjFlushInfo {
        public:
            std::string     filename;
            FILE *          obj_file; 
            std::list<
                    std::pair<
                        typename Mesh<Tm, Tv, Tf, R>::Vertex *,
                        uint32_t
                    >
                >                                                   last_boundary_vertices;
            uint32_t                                                last_flush_vertex_id;

            MeshObjFlushInfo() {
                this->obj_file  = NULL;
            }

            MeshObjFlushInfo(std::string const &filename) {
                this->filename                  = filename;
                this->obj_file  = fopen( (this->filename + ".obj").c_str(), "w");
                if (!this->obj_file) {
                    throw("MeshAlg::MeshObjFlushInfo::MeshObjFlushInfo(std:: string filemame): couldn't open given obj file for writing..\n");
                }
                this->last_flush_vertex_id  = 0;
                this->last_boundary_vertices.clear();
            }

            void    
            finalize()
            {
                if (this->obj_file) {
                    fclose(this->obj_file);
                }
                this->filename = std::string();
                this->last_boundary_vertices.clear();
                this->last_flush_vertex_id = 0;
            }
    };


    /* convenience wrapper function that takes care of everything, given only an initialized struct of type
     * MeshObjFlushInfo */
    template <typename Tm, typename Tv, typename Tf, typename R>
    void
    partialFlushToObjFile(
        Mesh<Tm, Tv, Tf, R>                                        &M,
        MeshObjFlushInfo<Tm, Tv, Tf, R>                            &M_flush_info,
        std::list<typename Mesh<Tm, Tv, Tf, R>::Face *>            &face_list);
}

#include "../tsrc/MeshAlgorithms.impl.hh"

#endif 
