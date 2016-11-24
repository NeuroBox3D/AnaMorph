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

#ifndef NLM_HH
#define NLM_HH

#include "Vec3.hh"
#include "BoundingBox.hh"
#include "MeshAlgorithms.hh"

/* forward declaration if NLM_CellNetwork, which is _not_ part of namespace NLM, but is build upon those components
 * declared and defined within namespace NLM. */
template <typename R> class NLM_CellNetwork;

namespace NLM {
    template <typename R> class CellNetworkInfo;
    template <typename R> class NeuronVertexInfo;
    template <typename R> class NeuronEdgeInfo;
    template <typename R> class SomaInfo;
    template <typename R> class NeuriteInfo;
    template <typename R> class AxonInfo;
    template <typename R> class DendriteInfo;
    template <typename R> class NeuriteSegmentInfo;
    template <typename R> class AxonSegmentInfo;
    template <typename R> class DendriteSegmentInfo;
    template <typename R> class NeuriteRootEdgeInfo;
    template <typename R> class AxonRootEdgeInfo;
    template <typename R> class DendriteRootEdgeInfo;


    enum NeuritePathParametrizationAlgos {
        CHORDLENGTH_PARAMETRIZATION,
        UNIFORM_PARAMETRIZATION,
        CENTRIPETAL_PARAMETRIZATION
    };

    /* Soma Sphere class */
    template <typename R>
    class SomaSphere {
        private:
            Vec3<R>             _centre;
            R                   _radius;

        public:
            SomaSphere() : _centre({0,0,0}), _radius(0)
            {
            }

            SomaSphere(
                Vec3<R> const  &centre,
                R const        &radius)
                    : _centre(centre), _radius(radius)
            {
            }

            Vec3<R> &
            centre()
            {
                return (this->_centre);
            }

            Vec3<R> const &
            centre() const
            {
                return (this->_centre);
            }

            R &
            radius()
            {
                return (this->_radius);
            }

            R const &
            radius() const
            {
                return (this->_radius);
            }

            BoundingBox<R>
            getBoundingBox() const
            {
                using Aux::VecMat::minVec3;
                using Aux::VecMat::maxVec3;
                using Aux::VecMat::onesVec3;
                using Aux::VecMat::fabsVec3;

                Vec3<R> rdisp   = Aux::VecMat::onesVec3<R>() * this->radius();
                Vec3<R> bb_min  = this->centre() - rdisp;
                Vec3<R> bb_max  = this->centre() + rdisp;
                Vec3<R> offset  = maxVec3<R>(
                                        Vec3<R>(1E-3, 1E-3, 1E-3),
                                        fabsVec3<R>(bb_max - bb_min)*0.025
                                    );
                bb_min         -= offset;
                bb_max         += offset;

                return BoundingBox<R>({ bb_min, bb_max });
            }

            bool
            contains(Vec3<R> const &p) const
            {
                return ( (p - this->centre()).len2() < this->radius() );
            }

            template <typename Tm, typename Tv, typename Tf>
            Mesh<Tm, Tv, Tf, R>
            generateMesh(uint32_t tessellation_depth = 4) const
            {
                Mesh<Tm, Tv, Tf, R> M_S;
                MeshAlg::generateIcoSphere(_centre, _radius, tessellation_depth, M_S);
                return M_S;
            }
    }; 

    /* NeuritePath class */
    template<typename R>
    class NeuritePath {
        friend class NLM_CellNetwork<R>;
        friend class
            CellNetwork<
                NLM::CellNetworkInfo<R>,
                NLM::NeuronVertexInfo<R>,
                NLM::NeuronEdgeInfo<R>,
                NLM::SomaInfo<R>,
                NLM::NeuriteInfo<R>,
                NLM::AxonInfo<R>,
                NLM::DendriteInfo<R>,
                NLM::NeuriteSegmentInfo<R>,
                NLM::AxonSegmentInfo<R>,
                NLM::DendriteSegmentInfo<R>,
                NLM::NeuriteRootEdgeInfo<R>,
                NLM::AxonRootEdgeInfo<R>,
                NLM::DendriteRootEdgeInfo<R>,
                R
            >;

        private:
            std::vector<
                    typename NLM_CellNetwork<R>::neurite_segment_iterator
                >                                   neurite_segments; 
            std::vector<R>                          neurite_vertex_parameters;
            std::vector<BLRCanalSurface<R> *>       canal_segments_magnified;

            NLM_CellNetwork<R>                     *C;

            bool                                    root_path;
            bool                                    geometry_updated;
            bool                                    bb_set;
            BoundingBox<R>                          bb;

        private:
            void                                    computeChordLengthParametrization();
            void                                    computeUniformParametrization();
            void                                    computeCentripetalParametrization();

        public:
                                                    //NeuritePath();
                                                    //NeuritePath(bool root_path = false);
                                                    NeuritePath(
                                                            typename NLM_CellNetwork<R>::neurite_segment_iterator   start_ns_it,
                                                            bool                                                    root_path);

                                                    NeuritePath(NeuritePath const &p);
            NeuritePath                            &operator=(NeuritePath const &p);
                                                   ~NeuritePath();

            uint32_t                                numVertices() const;
            uint32_t                                numEdges() const;

            uint32_t                                extend(typename NLM_CellNetwork<R>::neurite_segment_iterator ns_it);
            typename NLM_CellNetwork<R>::
                neurite_iterator                    operator[](uint32_t idx) const;

            R                                       getChordLength() const;

            void                                    updateGeometry(
                                                        std::function<
                                                                std::tuple<
                                                                    R,
                                                                    std::vector<R>,
                                                                    R
                                                                >(NLM::NeuritePath<R> const &P)
                                                            >  const                       &parametrization_algorithm);

            BoundingBox<R>                          getBoundingBox() const;

            Vec3<R>                                 findPermissibleRenderVector(
                                                        R lambda    = 0.2,
                                                        R mu        = 0.5) const;

        /* mesh generation */
            template<typename Tm, typename Tv, typename Tf>
            void                                    generateInitialSegmentMesh(
                                                        Mesh<Tm, Tv, Tf, R>                                    &M,
                                                        uint32_t                                                n_phi_segments,
                                                        Vec3<R>                                                 rvec,
                                                        R const                                                &phi_0,
                                                        R const                                                &arclen_dt,
                                                        bool                                                   &end_circle_offset,
                                                        std::vector<
                                                                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
                                                            >                                                  &end_circle_its,
                                                        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator          &closing_vertex_it,
                                                        R const                                                &radius_reduction_factor) const;


            template<typename Tm, typename Tv, typename Tf>
            void                                    appendTailMesh(
                                                        Mesh<Tm, Tv, Tf, R>                                    &M,
                                                        uint32_t                                                n_phi_segments,
                                                        Vec3<R>                                                 rvec,
                                                        R const                                                &phi_0,
                                                        R const                                                &arclen_dt,
                                                        bool const                                             &start_circle_offset,
                                                        std::vector<
                                                                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
                                                            > const                                            &start_circle_its,
                                                        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator const    &start_circle_closing_vertex_it,
                                                        bool                                                   *end_circle_offset   = NULL,
                                                        std::vector<
                                                                typename Mesh<Tm, Tv, Tf, R>::vertex_iterator
                                                            >                                                  *end_circle_its      = NULL,
                                                        typename Mesh<Tm, Tv, Tf, R>::vertex_iterator          *closing_vertex_it   = NULL) const;
                                                        

    };

    template<typename R>
    class CellNetworkInfo {
        friend class NLM_CellNetwork<R>;
    };

    template<typename R>
    class NeuronVertexInfo {
        friend class NLM_CellNetwork<R>;
    };

    template<typename R>
    class NeuronEdgeInfo {
        friend class NLM_CellNetwork<R>;
    };

    template<typename R>
    class SomaInfo {
        friend class NLM_CellNetwork<R>;

        private:
            /* information about the soma sphere */
            SomaSphere<R>                                       soma_sphere;

            /* neurite path graphs for every neurite: neurite paths are contained in the vertices and information about shared
             * vertices are contained in the edges. these graphs form trees for well-formed cells, so this list will be referred
             * to as the "path tree" associated to the cell tree rooted in the soma to which the instance of this object has
             * been annotated. */
            std::list<
                typename NLM_CellNetwork<R>::NeuritePathTree
            >                                                   neurite_path_trees;
    };

    template<typename R>
    class NeuriteInfo {
        friend class NLM_CellNetwork<R>;

        private:
            /* iterator to neurite path graph and list of pairs encoding the neurite paths and the indices within these
             * neurite paths that the neurite vertex is contained in. */
            typename NLM_CellNetwork<R>::NeuritePathTree   *npt;

            std::list<
                    std::pair<
                        typename NLM_CellNetwork<R>::
                            NeuritePathTree::
                            vertex_iterator,
                        uint32_t
                    >
                >                                           npt_info;
    };

    template<typename R>
    class AxonInfo : public NeuriteInfo<R> {
        friend class NLM_CellNetwork<R>;

        private:
    };

    template<typename R>
    class DendriteInfo : public NeuriteInfo<R> {
        friend class NLM_CellNetwork<R>;

        private:
    };

    template<typename R>
    class NeuriteSegmentInfo {
        friend class NLM_CellNetwork<R>;
        friend class NeuritePath<R>;

        private:
            /* pointer to neurite path tree, iterator of neurite path this neurite segment belongs to and the index of this
             * neurite segment in said neurite path. */
            typename NLM_CellNetwork<R>::NeuritePathTree   *npt;

            typename NLM_CellNetwork<R>::
                NeuritePathTree::vertex_iterator            npt_it;

            uint32_t                                        npt_ns_idx;

            /* magnified neurite segment canal surface */
            BLRCanalSurface<R>                                  canal_segment_magnified; 

            /* intersection-related members */
            bool                                                clean;

            /* radii of neighbours from smdv def */
            R                                                   v_src_rmax_nb;
            R                                                   v_dst_rmax_nb;

            bool                                                pmdv;
            bool                                                spmdv;
            bool                                                smdv;

            bool                                                reg;
            std::shared_ptr<
                    typename NLM_CellNetwork<R>::REG_IsecInfo
                >                                               reg_isec_info;

            bool                                                lsi;
            std::shared_ptr<
                    typename NLM_CellNetwork<R>::LSI_IsecInfo
                >                                               lsi_isec_info;


            bool                                                gsi;
            std::shared_ptr<
                    typename NLM_CellNetwork<R>::GSI_IsecInfo
                >                                               gsi_isec_info;

            bool                                                rc_sons;
            bool                                                ic_sons;

            bool                                                rc_nsns;
            bool                                                icrn_nsns;
            bool                                                icin_nsns;

            /* lists of intersection info structs: soma-neurite segment and neurite segment-neurite segmenta
             * */
            std::list<
                    typename NLM_CellNetwork<R>::
                        NSNS_IsecInfo *
                >                                           nsns_intersections;

            std::list<
                    typename NLM_CellNetwork<R>::
                        SONS_IsecInfo *
                >                                           sons_intersections;

        public:
            NeuriteSegmentInfo() : reg_isec_info(), lsi_isec_info(), gsi_isec_info()
            {
                this->npt           = NULL;
                this->npt_it.explicitlyInvalidate();
                this->npt_ns_idx    = 0;

                this->clean         = false;
                this->pmdv          = false;
                this->spmdv         = false;
                this->smdv          = false;
                this->reg           = false;
                //this->reg_isec_info = NULL;
                this->lsi           = false;
                //this->lsi_isec_info = NULL;
                this->gsi           = false;
                //this->gsi_isec_info = NULL;
                this->rc_nsns       = false;
                this->rc_sons       = false;

                this->icrn_nsns     = false;
                this->ic_sons       = false;
                this->icin_nsns     = false;
            }

            void
            resetIntersectionStatus()
            {
                /* empty shared pointers. this releases owernship and decrements reference count. usually, the allocated
                 * intersection info objects will be deleted if not explicitly referenced elsewhere. */
                this->reg_isec_info.reset();
                this->lsi_isec_info.reset();
                this->gsi_isec_info.reset();

                /* reset flags */
                this->clean         = false;
                this->pmdv          = false;
                this->spmdv         = false;
                this->smdv          = false;
                this->reg           = false;
                this->lsi           = false;
                this->gsi           = false;
                this->rc_nsns       = false;
                this->rc_sons       = false;
                this->icrn_nsns     = false;
                this->ic_sons       = false;
                this->icin_nsns     = false;
            }

            void
            updateREGStatus(typename NLM_CellNetwork<R>::REG_IsecInfo * const &reg_info = NULL)
            {
                if (reg_info) {
                    this->reg_isec_info.reset(reg_info);
                    this->reg = true;
                }
                else {
                    /* reset shared pointer (usually frees allocated REG_IsecInfo) and reset reg flag */
                    this->reg_isec_info.reset();
                    this->reg = false;
                }
            }
            void
            updateLSIStatus(typename NLM_CellNetwork<R>::LSI_IsecInfo * const &lsi_info)
            {
                if (lsi_info) {
                    this->lsi_isec_info.reset(lsi_info);
                    this->lsi = true;
                }
                else {
                    /* reset shared pointer (usually frees allocated LSI_IsecInfo) and reset lsi flag */
                    this->lsi_isec_info.reset();
                    this->lsi = false;
                }
            }

            void
            updateGSIStatus(typename NLM_CellNetwork<R>::GSI_IsecInfo * const &gsi_info)
            {
                if (gsi_info) {
                    this->gsi_isec_info.reset(gsi_info);
                    this->gsi = true;
                }
                else {
                    /* reset shared pointer (usually frees allocated GSI_IsecInfo) and reset gs flag */
                    this->gsi_isec_info.reset();
                    this->gsi = false;
                }
            }
    };

    template<typename R>
    class AxonSegmentInfo : public NeuriteSegmentInfo<R> {
        friend class NLM_CellNetwork<R>;

        private:
    };

    template<typename R>
    class DendriteSegmentInfo : public NeuriteSegmentInfo<R> {
        friend class NLM_CellNetwork<R>;

        private:
    };

    template<typename R>
    class NeuriteRootEdgeInfo {
        friend class NLM_CellNetwork<R>;
        friend class NeuritePath<R>;
        friend class
            CellNetwork<
                NLM::CellNetworkInfo<R>,
                NLM::NeuronVertexInfo<R>,
                NLM::NeuronEdgeInfo<R>,
                NLM::SomaInfo<R>,
                NLM::NeuriteInfo<R>,
                NLM::AxonInfo<R>,
                NLM::DendriteInfo<R>,
                NLM::NeuriteSegmentInfo<R>,
                NLM::AxonSegmentInfo<R>,
                NLM::DendriteSegmentInfo<R>,
                NLM::NeuriteRootEdgeInfo<R>,
                NLM::AxonRootEdgeInfo<R>,
                NLM::DendriteRootEdgeInfo<R>,
                R
            >;
        friend class
            CellNetwork<
                Common::UnitType,       
                Common::UnitType,       
                Common::UnitType,       
                NLM::CellNetworkInfo<R>,
                NLM::SomaInfo<R>,
                NLM::AxonInfo<R>,
                NLM::DendriteInfo<R>,
                NLM::AxonSegmentInfo<R>,
                NLM::DendriteSegmentInfo<R>,
                NLM::NeuriteRootEdgeInfo<R>,
                NLM::NeuriteRootEdgeInfo<R>,
                R
            >;

        private:
            /* initial cylinder for neurite root edge encoded as a bezier canal surface with linearly interpolated radii. */
            BLRCanalSurface<R> initial_cylinder;

        public:
            NeuriteRootEdgeInfo()
            {
            }
    };

    template <typename R>
    class AxonRootEdgeInfo : public NeuriteRootEdgeInfo<R> {
        friend class NLM_CellNetwork<R>;

        private:
    };

    template <typename R>
    class DendriteRootEdgeInfo : public NeuriteRootEdgeInfo<R> {
        friend class NLM_CellNetwork<R>;

        private:
    };

    /* structs to hold information about roots, focal points, stationary points, etc.. and
     * the corresponding values of univaraite and bivariate polynomials / polynomial systems */
    template <typename R>
    struct p2 {
        PolyAlg::RealInterval<R>    interval;
        R                           value;

        p2(
            PolyAlg::RealInterval<R> const &I,
            R const                        &v)
                : interval(I), value(v)
        {
        }
    };

    template <typename R>
    struct p3 {
        PolyAlg::RealRectangle<R>   rectangle;
        R                           value;

        p3(
            PolyAlg::RealRectangle<R> const    &I,
            R                         const    &v)
        {
            this->rectangle = I;
            this->value     = v;
        }
    };

}

#endif
