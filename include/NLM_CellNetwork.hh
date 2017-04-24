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

#ifndef NLM_CELL_NETWORK_HH
#define NLM_CELL_NETWORK_HH

#include "CellNetwork.hh"
#include "CanalSurface.hh"
#include "NLM.hh"

/* forward declarations */
template <typename R> class NLM_CellNetwork;

/* derived specialization of CellNetwork for the purpose of the non-linear goemetric model mainly used in the thesis */
template<typename R = double>
class NLM_CellNetwork
    : public CellNetwork<
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
    >
{    
    /* friends.. */
    friend class NLM::NeuritePath<R>;
    friend class NLM::NeuriteSegmentInfo<R>;

    /* abbreviate template base class type, the lengthy CellNetwork specialization, with NLM_CellNetworkBaseType */
    protected:
        typedef CellNetwork<
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
        > NLM_CellNetworkBaseType;


    /* settings used during analysis / mesh generation */
    protected:
        uint32_t        analysis_nthreads;
        R               analysis_univar_solver_eps;
        R               analysis_bivar_solver_eps;

        R               partition_filter_angle;
        R               partition_filter_max_ratio_ratio;
        std::function<
                typename NLM_CellNetwork<R>::neurite_segment_iterator(
                    NLM_CellNetwork<R> &,
                    typename NLM_CellNetwork<R>::neurite_iterator
                )
            >           partition_algo;

        std::function<
                std::tuple<
                    R,
                    std::vector<R>,
                    R
                >(NLM::NeuritePath<R> const &P)
            >           parametrization_algo;

        bool            meshing_flush;
        uint32_t        meshing_flush_face_limit;

        uint32_t        meshing_n_soma_refs;
        uint32_t        meshing_canal_segment_n_phi_segments;
        uint32_t        meshing_outer_loop_maxiter;
        uint32_t        meshing_inner_loop_maxiter;

        bool            meshing_preserve_crease_edges;
        R				meshing_cansurf_triangle_height_factor;

        R               meshing_radius_factor_decrement;
        R               meshing_complex_edge_max_growth_factor;

    /* pull iterators / vertices /edges from CellNetwork scope, since qualified access is cumbersome */
    public:
        using typename NLM_CellNetworkBaseType::neuron_iterator;
        using typename NLM_CellNetworkBaseType::neuron_const_iterator;

        using typename NLM_CellNetworkBaseType::soma_iterator;
        using typename NLM_CellNetworkBaseType::soma_const_iterator;

        using typename NLM_CellNetworkBaseType::neurite_iterator;
        using typename NLM_CellNetworkBaseType::neurite_const_iterator;
        
        using typename NLM_CellNetworkBaseType::axon_iterator;
        using typename NLM_CellNetworkBaseType::axon_const_iterator;

        using typename NLM_CellNetworkBaseType::dendrite_iterator;
        using typename NLM_CellNetworkBaseType::dendrite_const_iterator;

        using typename NLM_CellNetworkBaseType::neuron_edge_iterator;
        using typename NLM_CellNetworkBaseType::neuron_edge_const_iterator;

        using typename NLM_CellNetworkBaseType::neurite_rootedge_iterator;
        using typename NLM_CellNetworkBaseType::neurite_rootedge_const_iterator;

        using typename NLM_CellNetworkBaseType::axon_rootedge_iterator;
        using typename NLM_CellNetworkBaseType::axon_rootedge_const_iterator;

        using typename NLM_CellNetworkBaseType::dendrite_rootedge_iterator;
        using typename NLM_CellNetworkBaseType::dendrite_rootedge_const_iterator;

        using typename NLM_CellNetworkBaseType::neurite_segment_iterator;
        using typename NLM_CellNetworkBaseType::neurite_segment_const_iterator;

        using typename NLM_CellNetworkBaseType::axon_segment_iterator;
        using typename NLM_CellNetworkBaseType::axon_segment_const_iterator;

        using typename NLM_CellNetworkBaseType::dendrite_segment_iterator;
        using typename NLM_CellNetworkBaseType::dendrite_segment_const_iterator;

        using typename NLM_CellNetworkBaseType::NeuronVertex;
        using typename NLM_CellNetworkBaseType::SomaVertex;
        using typename NLM_CellNetworkBaseType::NeuriteVertex;
        using typename NLM_CellNetworkBaseType::AxonVertex;
        using typename NLM_CellNetworkBaseType::DendriteVertex;

        using typename NLM_CellNetworkBaseType::NeuronEdge;
        using typename NLM_CellNetworkBaseType::NeuriteRootEdge;
        using typename NLM_CellNetworkBaseType::AxonRootEdge;
        using typename NLM_CellNetworkBaseType::DendriteRootEdge;
        using typename NLM_CellNetworkBaseType::NeuriteSegment;
        using typename NLM_CellNetworkBaseType::AxonSegment;
        using typename NLM_CellNetworkBaseType::DendriteSegment;

        /* struct containing information about edges in neurite path graphs */
        struct NPTEdgeInfo
        {
            typename
                NLM_CellNetwork<R>::
                neurite_iterator        shared_vertex;
            uint32_t                    shared_vertex_src_path_idx;
        };

        /* typedef for neurite path graph */
        typedef Graph<
            /* iterator to neurite root edge is contained in neurite path tree */
            typename NLM_CellNetworkBaseType::neurite_rootedge_iterator,
            /* its vertices contain neurite paths */
            NLM::NeuritePath<R>,
            /* info class about edges in neurite path trees */
            NPTEdgeInfo
        > NeuritePathTree;
 
    /* settings class and methods to retrieve / update settings */
    public:
        struct Settings {
            uint32_t        analysis_nthreads;
            R               analysis_univar_solver_eps;
            R               analysis_bivar_solver_eps;

            /*
            R               partition_filter_angle;
            R               partition_filter_max_ratio_ratio;
            */
            std::function<
                    typename NLM_CellNetwork<R>::neurite_segment_iterator(
                        NLM_CellNetwork<R> &,
                        typename NLM_CellNetwork<R>::neurite_iterator
                    )
                >           partition_algo;

            std::function<
                    std::tuple<
                        R,
                        std::vector<R>,
                        R
                    >(NLM::NeuritePath<R> const &P)
                >           parametrization_algo;

            bool            meshing_flush;
            uint32_t        meshing_flush_face_limit;

            uint32_t        meshing_n_soma_refs;
            uint32_t        meshing_canal_segment_n_phi_segments;
            uint32_t        meshing_outer_loop_maxiter;
            uint32_t        meshing_inner_loop_maxiter;

            bool            meshing_preserve_crease_edges;
            R				meshing_cansurf_triangle_height_factor;

            R               meshing_radius_factor_decrement;
            R               meshing_complex_edge_max_growth_factor;
        };

        Settings            getSettings() const;
        void                updateSettings(Settings const &s);

    /* protected types / enums related to intersection analysis / classification in non-linear model. */
    protected:
        /* intersection types */
        enum IntersectionTypes {
            REG,
            LSI,
            GSI,
            SOSO,
            RC_SONS,
            IC_SONS,
            RC_NSNS,
            ICRN_NSNS,
            ICIN_NSNS,
            ICIN_NSNSA
        };

        /* intersection job types */
        enum NLM_ISEC_JOB_TYPES {
            JOB_REG             = 0,
            JOB_LSI             = 1,
            JOB_GSI             = 2,
            JOB_SONS            = 3,
            JOB_NS_NS_ADJ       = 4,
            JOB_NS_NS_NONADJ    = 5
        };

        enum NLM_ISEC_JOB_STATES {
            JOB_UNPROCESSED,
            JOB_IN_PROCESS,
            JOB_DONE
        };

        struct IsecJob {
            uint32_t            job_state;


            /* solver tolerances */
            R                   univar_solver_eps;
            R                   bivar_solver_eps;

            /* general result: true <=> intersection, false: ok */
            bool                result;

            IsecJob(
                R const    &univar_solver_eps,
                R const    &bivar_solver_eps) 
            {
                this->job_state         = JOB_UNPROCESSED;
                this->univar_solver_eps = univar_solver_eps;
                this->bivar_solver_eps  = bivar_solver_eps;
                this->result            = false;
            }

            virtual
           ~IsecJob()
            {
            }

            virtual uint32_t
            type() const = 0;

            void
            reset()
            {
                this->job_state         = JOB_UNPROCESSED;
                this->result            = false;
            }
        };

        struct REG_Job : public IsecJob {
            neurite_segment_const_iterator  ns_it;
            std::vector<NLM::p2<R>>         checkpoly_roots;

            REG_Job(
                neurite_segment_const_iterator const   &ns_it,
                R const                                &univar_solver_eps)
                    : IsecJob(univar_solver_eps, Aux::Numbers::inf<R>())
            {
                this->ns_it             = ns_it;
                this->checkpoly_roots   = checkpoly_roots;
            }    

            virtual uint32_t
            type() const final
            {
                return JOB_REG;
            }
        };

        struct LSI_Job : public IsecJob {
            neurite_segment_const_iterator  ns_it;
            std::vector<NLM::p2<R>>         lsi_neg_points;     

            LSI_Job(
                neurite_segment_const_iterator const   &ns_it,
                R const                                &univar_solver_eps)
                    : IsecJob(univar_solver_eps, Aux::Numbers::inf<R>())
            {
                this->ns_it = ns_it;
            }    

            virtual uint32_t
            type() const final
            {
                return JOB_LSI;
            }
        };

        struct GSI_Job : public IsecJob {
            neurite_segment_const_iterator  ns_it;
            std::vector<NLM::p3<R>>         gsi_stat_points;

            GSI_Job(
                neurite_segment_const_iterator const   &ns_it,
                R const                                &univar_solver_eps,
                R const                                &bivar_solver_eps)
                    : IsecJob(univar_solver_eps, bivar_solver_eps)
            {
                this->ns_it = ns_it;
            }    

            virtual uint32_t
            type() const final
            {
                return JOB_GSI;
            }
        };

        struct SONS_Job : public IsecJob {
            soma_const_iterator             s_it;
            neurite_segment_const_iterator  ns_it;
            std::vector<NLM::p2<R>>         isec_stat_points; 
            bool                            neurite_root_segment;

            SONS_Job(
                soma_const_iterator const              &s_it,
                neurite_segment_const_iterator const   &ns_it,
                R const                                &univar_solver_eps)
                    : IsecJob(univar_solver_eps, Aux::Numbers::inf<R>())
            {
                this->s_it                  = s_it;
                this->ns_it                 = ns_it;
                this->neurite_root_segment  = this->ns_it->getSourceVertex()->isNeuriteRootVertex();
            }

            virtual uint32_t
            type() const final
            {
                return JOB_SONS;
            }
        };

        struct NSNS_Job : public IsecJob {
            neurite_segment_const_iterator  ns_first_it;
            neurite_segment_const_iterator  ns_second_it;
            std::vector<NLM::p3<R>>         isec_stat_points;

            NSNS_Job(
                neurite_segment_const_iterator const   &ns_first_it,
                neurite_segment_const_iterator const   &ns_second_it,
                R const                                &univar_solver_eps,
                R const                                &bivar_solver_eps)
                    : IsecJob(univar_solver_eps, bivar_solver_eps)
            {
                this->ns_first_it   = ns_first_it;
                this->ns_second_it  = ns_second_it;
            }
        };

        struct NSNS_NonAdj_Job : public NSNS_Job {
            NSNS_NonAdj_Job(
                neurite_segment_const_iterator const   &ns_first_it,
                neurite_segment_const_iterator const   &ns_second_it,
                R const                                &univar_solver_eps,
                R const                                &bivar_solver_eps)
                    : NSNS_Job(ns_first_it, ns_second_it, univar_solver_eps, bivar_solver_eps)
            {
            }

            virtual uint32_t
            type() const final
            {
                return JOB_NS_NS_NONADJ;
            }
        };

        struct NSNS_Adj_Job : public NSNS_Job {
            bool fst_end_snd_start;

            NSNS_Adj_Job(
                neurite_segment_const_iterator const   &ns_first_it,
                neurite_segment_const_iterator const   &ns_second_it,
                bool                                    fst_end_snd_start,
                R const                                &univar_solver_eps,
                R const                                &bivar_solver_eps)
                    : NSNS_Job(ns_first_it, ns_second_it, univar_solver_eps, bivar_solver_eps)
            {
                this->fst_end_snd_start = fst_end_snd_start;
            }

            virtual uint32_t
            type() const final
            {
                return JOB_NS_NS_ADJ;
            }
        };

        /* intersection info for intersection types that are not defined per neurite segment *-NSNS and *-SONS */
        struct IsecInfo {
            protected:
                NLM_CellNetwork<R>         &C;

                /* protected ctor */
                IsecInfo(NLM_CellNetwork<R> &network) : C(network)
                {
                }
            
            public:
                virtual uint32_t            type() const = 0;

        };

        /* reg info. no further differentiation => can be constructed with a conversion constructor from REG_Job */
        struct REG_IsecInfo : public IsecInfo {
            neurite_segment_iterator    ns_it;
            std::vector<NLM::p2<R>>     checkpoly_roots;

            /* conversion constructor from REG_Job */
            REG_IsecInfo(
                NLM_CellNetwork<R> &network,
                REG_Job const      &reg_job) : IsecInfo(network)
            {
                /* convert to non-const iterator by locating neurite segment via id */
                this->ns_it             = this->C.neurite_segments.find(reg_job.ns_it->id());
                this->checkpoly_roots   = reg_job.checkpoly_roots;

                if (this->ns_it == this->C.neurite_segments.end()) {
                    throw("NLM_CellNetwork::REG_IsecInfo(REG_Job const &reg_job) (conversion ctor): failed to locate "\
                        "neurite segment contained in reg job via id. internal logic error.");
                }

            }

            virtual uint32_t
            type() const final
            {
                return REG;
            }
        };

        /* lsi isec info. no further differentiation => can be constructed with a conversion constructor from LSI_Job */
        struct LSI_IsecInfo : public IsecInfo {
            neurite_segment_iterator    ns_it;
            std::vector<NLM::p2<R>>     lsi_neg_points;

            /* conversion constructor from LSI_Job */
            LSI_IsecInfo(
                NLM_CellNetwork<R> &network,
                LSI_Job const      &lsi_job) : IsecInfo(network)
            {
                /* convert to non-const iterator by locating neurite segment via id */
                this->ns_it             = this->C.neurite_segments.find(lsi_job.ns_it->id());
                this->lsi_neg_points    = lsi_job.lsi_neg_points;

                if (this->ns_it == this->C.neurite_segments.end()) {
                    throw("NLM_CellNetwork::LSI_IsecInfo(LSI_Job const &lsi_job) (conversion ctor): failed to locate "\
                        "neurite segment contained in lsi job via id. internal logic error.");
                }
            }

            virtual uint32_t
            type() const final
            {
                return LSI;
            } 
        };

        /* gsi isec info. no further differentiation => can be constructed with a conversion constructor from GSI_Job */
        struct GSI_IsecInfo : public IsecInfo {
            neurite_segment_iterator    ns_it;
            std::vector<NLM::p3<R>>     gsi_stat_points;

            /* conversion constructor from GSI_Job */
            GSI_IsecInfo(
                NLM_CellNetwork<R> &network,
                GSI_Job const      &gsi_job) : IsecInfo(network)
            {
                /* convert to non-const iterator by locating neurite segment via id */
                this->ns_it             = this->C.neurite_segments.find(gsi_job.ns_it->id());
                this->gsi_stat_points   = gsi_job.gsi_stat_points;

                if (this->ns_it == this->C.neurite_segments.end()) {
                    throw("NLM_CellNetwork::GSI_IsecInfo(GSI_Job const &gsi_job) (conversion ctor): failed to locate "\
                        "neurite segment contained in gsi job via id. internal logic error.");
                }
            }

            virtual uint32_t
            type() const final
            {
                return GSI;
            } 
        };

        /* SOSO intersections, can be constructed with a conversion constructor from SOSO_Job */
        struct SOSO_IsecInfo : public IsecInfo {
            soma_iterator               s_it_first;
            soma_iterator               s_it_second;
        };

        /* sons-intersections. further distinguished into IC_SONS and RC_SONS with different analysis procedures =>
         * create this from a SONS_Job with a method that performs the necessary analysis. */
        struct SONS_IsecInfo : public IsecInfo {
            soma_iterator               s_it;
            neurite_segment_iterator    ns_it;
            std::vector<NLM::p2<R>>     isec_stat_points;

            SONS_IsecInfo(
                NLM_CellNetwork<R>             &network,
                SONS_Job const                 &sons_job)
                    : IsecInfo(network)
            {
                /* convert to non-const iterators */
                this->s_it              = this->C.soma_vertices.find(sons_job.s_it->id());
                this->ns_it             = this->C.neurite_segments.find(sons_job.ns_it->id());
                this->isec_stat_points  = sons_job.isec_stat_points;

                if (this->s_it == this->C.soma_vertices.end() || this->ns_it == this->C.neurite_segments.end()) {
                    throw("NLM_CellNetwork::SONS_IsecInfo(SONS_Job const &gsi_job) (conversion ctor): failed to locate "\
                        "soma iterator or neurite segment contained in sons job via id. internal logic error.");
                }
            }
        };

        struct RC_SONS_IsecInfo : public SONS_IsecInfo {
            RC_SONS_IsecInfo(
                NLM_CellNetwork<R>             &network,
                SONS_Job const                 &sons_job)
                    : SONS_IsecInfo(network, sons_job)
            {
            }

            virtual uint32_t
            type() const final
            {
                return RC_SONS;
            } 
        };

        struct IC_SONS_IsecInfo : public SONS_IsecInfo {
            std::list<neurite_iterator> path;
            R                           pathlen, pathlen_reduced;

            IC_SONS_IsecInfo(
                NLM_CellNetwork<R>             &network,
                SONS_Job const                 &sons_job)
                    : SONS_IsecInfo(network, sons_job)
            {
                /* compute path information for intra-cell soma-neurite segment info */
            }

            virtual uint32_t
            type() const final
            {
                return IC_SONS;
            } 
        };

        /* method to create a SONS_IsecInfo of matching specialized type (and gather respective specialized intersection
         * information) from a SONS_Job intersection job object */
        std::shared_ptr<SONS_IsecInfo> 
        generateSONSIsecInfo(SONS_Job const &sons_job)
        {
            /* get non-const versions of soma and neurite segment iterator stored in sons_job */
            auto s_it   = this->soma_vertices.find(sons_job.s_it->id());
            auto ns_it  = this->neurite_segments.find(sons_job.ns_it->id());

            if (s_it == this->soma_vertices.end() || ns_it == this->neurite_segments.end()) {
                throw("NLM_CellNetwork::generateSONSIsecInfo(): soma or neurite segment stored in sons job not found. "\
                    "internal logic error.");
            }

            /* analyse SONS intersection: if the neurite segment belongs to the cell rooted in the intersected soma =>
             * IC_SONS, otherwise RC_SONS. */
            if (sons_job.ns_it->getSoma() == s_it) {
                IC_SONS_IsecInfo *ic_sons_info
                    = new IC_SONS_IsecInfo(
                        (*this),
                        s_it,
                        ns_it,
                        sons_job.isec_stat_points);

                /* return shared pointer containing implicitly up-cast SONS_IsecInfo */
                return std::shared_ptr<SONS_IsecInfo>(ic_sons_info);
            }
            else {
                RC_SONS_IsecInfo *rc_sons_info
                    = new RC_SONS_IsecInfo(
                        (*this),
                        s_it,
                        ns_it,
                        sons_job.isec_stat_points);

                /* return shared pointer containing implicitly up-cast SONS_IsecInfo */
                return std::shared_ptr<SONS_IsecInfo>(rc_sons_info);
            } 
        }

        /* nsns-intersections. further distinguished into RC_NSNS, ICRN_NSNS, ICIN_NSNS and ICIN_NSNSa with different
         * analysis procedures => create this from a NSNS_Job with a method that performs the necessary analysis. */
        struct NSNS_IsecInfo : public IsecInfo {
            neurite_segment_iterator    ns_first_it, ns_second_it;
            std::vector<NLM::p3<R>>     isec_stat_points;

            NSNS_IsecInfo(
                NLM_CellNetwork<R>             &network,
                neurite_segment_iterator const &ns_first_it,
                neurite_segment_iterator const &ns_second_it,
                std::vector<NLM::p3<R>> const  &isec_stat_points) : IsecInfo(network)
            {
                this->ns_first_it       = ns_first_it;
                this->ns_second_it      = ns_second_it;
                this->isec_stat_points  = isec_stat_points;
            }
        };

        struct RC_NSNS_IsecInfo : public NSNS_IsecInfo {
            RC_NSNS_IsecInfo(
                NLM_CellNetwork<R>             &network,
                neurite_segment_iterator const &ns_first_it,
                neurite_segment_iterator const &ns_second_it,
                std::vector<NLM::p3<R>> const  &isec_stat_points)
                    : NSNS_IsecInfo(network, ns_first_it, ns_second_it, isec_stat_points)
            {
            }

            virtual uint32_t
            type() const final
            {
                return RC_NSNS;
            } 
        };

        struct ICRN_NSNS_IsecInfo : public NSNS_IsecInfo {
            ICRN_NSNS_IsecInfo(
                NLM_CellNetwork<R>             &network,
                neurite_segment_iterator const &ns_first_it,
                neurite_segment_iterator const &ns_second_it,
                std::vector<NLM::p3<R>> const  &isec_stat_points)
                    : NSNS_IsecInfo(network, ns_first_it, ns_second_it, isec_stat_points)
            {
                /* calculate intER-neurite path information */
            }


            virtual uint32_t
            type() const final
            {
                return ICRN_NSNS;
            }
        };

        struct ICIN_NSNS_IsecInfo : public NSNS_IsecInfo {
            bool                        adjacent;
            bool                        fst_end_snd_start;
            bool                        critical;
            std::list<neurite_iterator> path;
            R                           pathlen, pathlen_reduced;

            ICIN_NSNS_IsecInfo(
                NLM_CellNetwork<R>             &network,
                neurite_segment_iterator const &ns_first_it,
                neurite_segment_iterator const &ns_second_it,
                std::vector<NLM::p3<R>> const  &isec_stat_points,
                bool                            _adjacent            = false,
                bool                            _fst_end_snd_start   = false)
            : NSNS_IsecInfo(network, ns_first_it, ns_second_it, isec_stat_points),
              adjacent(_adjacent), fst_end_snd_start(_fst_end_snd_start), critical(false)
            {
                /* calculate inTRA-neurite path information */
            }

            virtual uint32_t
            type() const final
            {
                if (adjacent) {
                    return ICIN_NSNSA;
                }
                else {
                    return ICIN_NSNS;
                }
            }
        };

        /* method to create a NSNS_IsecInfo of matching specialized type (and gather respective specialized intersection
         * information) from an NSNS_Job intersection job object */
        std::shared_ptr<NSNS_IsecInfo> 
        generateNSNSIsecInfo(NSNS_Job const &nsns_job)
        {
            /* get non-const versions of neurite segment iterators stored in nsns job. */
            auto ns_first_it    = this->neurite_segments.find(nsns_job.ns_first_it->id());
            auto ns_second_it   = this->neurite_segments.find(nsns_job.ns_second_it->id());

            if (ns_first_it == this->neurite_segments.end() || ns_second_it == this->neurite_segments.end()) {
                throw("NLM_CellNetwork::generateNSNSIsecInfo(): at least one neurite segment stored in NSNS job not found. "\
                    "internal logic error.");
            }

            /* create specialized objects depending on the type of the NSNS intersection and wrap them in a shared ptr
             * containing an upcast NSNS_IsecInfo pointer to the created specialized object. */

            /* inter-cell */
            if (ns_first_it->getSoma() != ns_second_it->getSoma()) {
                RC_NSNS_IsecInfo *rc_nsns_info = 
                    new RC_NSNS_IsecInfo(
                        *this,
                        ns_first_it,
                        ns_second_it,
                        nsns_job.isec_stat_points
                    );
                return std::shared_ptr<NSNS_IsecInfo>(rc_nsns_info);
            }
            /* intra-cell */
            else {
                /* intra-cell inter-neurite */
                if (ns_first_it->getNeurite() != ns_second_it->getNeurite()) {
                    ICRN_NSNS_IsecInfo *icrn_nsns_isec_info = 
                        new ICRN_NSNS_IsecInfo(
                            *this,
                            ns_first_it,
                            ns_second_it,
                            nsns_job.isec_stat_points
                        );
                    return std::shared_ptr<NSNS_IsecInfo>(icrn_nsns_isec_info);
                }
                /* intra-cell-intra-neurite */
                else {
                    ICIN_NSNS_IsecInfo  *icin_nsns_isec_info;

                    /* check if neurite segments are adjacent */
                    if (ns_first_it->getSourceVertex() == ns_second_it->getSourceVertex()) {
                        icin_nsns_isec_info
                            = new ICIN_NSNS_IsecInfo(
                                *this,
                                ns_first_it,
                                ns_second_it,
                                nsns_job.isec_stat_points,
                                true,
                                false);
                    }
                    else if (ns_first_it->getDestinationVertex() == ns_second_it->getSourceVertex()) {
                        icin_nsns_isec_info
                            = new ICIN_NSNS_IsecInfo(
                                *this,
                                /* destination vertex of first ns is source vertex of second ns => order (first, second)
                                 * */
                                ns_first_it,
                                ns_second_it,
                                nsns_job.isec_stat_points,
                                true,
                                true);
                    }
                    else if (ns_first_it->getSourceVertex() == ns_second_it->getDestinationVertex()) {
                        icin_nsns_isec_info
                            = new ICIN_NSNS_IsecInfo(
                                *this,
                                /* order REVERSED compared to the above case: destination vertex of second ns is source
                                 * vertex of first ns => order (first, second) */
                                ns_second_it,
                                ns_first_it,
                                nsns_job.isec_stat_points,
                                true,
                                true);
                    }
                    else if (ns_first_it->getDestinationVertex() == ns_second_it->getDestinationVertex()) {
                        throw("NLM_CellNetwork::generateNSNSIsecInfo(): two neurite segments contained in NSNS_Job "\
                            "share common end vertex. internal logic error.");
                    }
                    /* non-adjacent. */
                    else {
                        icin_nsns_isec_info
                            = new ICIN_NSNS_IsecInfo(
                                *this,
                                ns_first_it,
                                ns_second_it,
                                nsns_job.isec_stat_points,
                                false,
                                false);

                    }

                    return std::shared_ptr<NSNS_IsecInfo>(icin_nsns_isec_info);
                }
            }
        }

        /* thread related types */
        enum ThreadState {
            THREAD_IDLE,
            THREAD_PROCESSING,
            THREAD_DONE
        };

        struct ThreadInfo {
            uint32_t                            thread_id;
            uint32_t                            thread_state;
            std::thread                         thread;
            std::mutex                          mutex;
            std::list<std::shared_ptr<IsecJob>> job_list;

            ThreadInfo(uint32_t id) : thread(), mutex()
            {
                thread_id       = id;
                thread_state    = THREAD_IDLE;
            }

            ThreadInfo(ThreadInfo const &ti) = delete;
            ThreadInfo &operator=(ThreadInfo const &ti) = delete;
            /* move constructor */

            ThreadInfo (ThreadInfo &&ti)
                /* mutex can't be moved, construct a fresh one.. */
                : mutex() {
                this->thread_id     = ti.thread_id; 
                this->thread_state  = ti.thread_state;
                /* move thread handler from ti to (this) with move constructor, this invalidates the thread in ti, yet
                 * leaves the (logical) thread process completely unaffected */
                this->thread        = std::move(ti.thread);
                this->job_list      = std::move(ti.job_list);
            }

            void
            reset()
            {
                thread_state    = THREAD_IDLE;
                job_list.clear();
            }
        };


    protected:
        /* partitioning of individual neurite sub-tree starting with neurite segment e = (u, v), where the (non-public)
         * caller must ensure that e's source vertex u is a neurite vertex of the cell-tree rooted in soma s.
         * selection_algorithm is used to end a neurite path or select a successor among the children of the currently
         * considered neurite branching vertex. */
        void                                        partitionNeuriteSubTree(
                                                        soma_iterator                                       s_it,
                                                        neurite_segment_iterator                            e_it,
                                                        std::function<
                                                                neurite_iterator(
                                                                    NLM_CellNetwork<R> &,
                                                                    neurite_iterator
                                                                )
                                                            >  const                                       &selection_algorithm);
        /* static intersection analysis methods */
        static bool                                 checkCanalSegmentRegularity(
                                                        BLRCanalSurface<3u, R> const   &Gamma,
                                                        R const                    &univar_solver_eps,
                                                        std::vector<NLM::p2<R>>    &checkpoly_roots);

        static bool                                 checkSomaNeuriteIntersection(
                                                        NLM::SomaSphere<R> const   &S,
                                                        BLRCanalSurface<3u, R> const   &Gamma,
                                                        bool                        neurite_root_segment,
                                                        R const                    &univar_solver_eps,
                                                        std::vector<NLM::p2<R>>    &isec_stat_points);

        /* local self intersection of neurite canal segment */
        static bool                                 checkNeuriteLocalSelfIntersection(
                                                        BLRCanalSurface<3u, R> const   &Gamma,
                                                        R const                    &univar_solver_eps,
                                                        std::vector<NLM::p2<R>>    &lsi_neg_points);

        /* global self intersection of neurite canal segment */
        static bool                                 checkNeuriteGlobalSelfIntersection(
                                                        BLRCanalSurface<3u, R> const   &Gamma,
                                                        R const                    &univar_solver_eps,
                                                        R const                    &bivar_solver_eps,
                                                        std::vector<NLM::p3<R>>    &gsi_stat_points);   

        /* neurite / neurite intersection for non-adjacent neurite segment canal surfaces */
        static bool                                 checkNeuriteNeuriteIntersection(
                                                        BLRCanalSurface<3u, R> const   &Gamma,
                                                        BLRCanalSurface<3u, R> const   &Delta,
                                                        R const                    &univar_solver_eps,
                                                        R const                    &bivar_solver_eps,
                                                        std::vector<NLM::p3<R>>    &isec_stat_points);

        /* same for adjacent neurite canal segments. if Gamma and Delta do not share their starting
         * point, then Gamma and Delta MUST be given in the order that  satisfies gamma(1.0) = delta(0.0),
         * i.e. the endpoint of Gamma's spine curve must be the starting point of Delta's spine curve. */
        static bool                                 checkAdjacentNeuriteNeuriteIntersection(
                                                        BLRCanalSurface<3u, R> const   &Gamma,
                                                        BLRCanalSurface<3u, R> const   &Delta,
                                                        R const                    &univar_solver_eps,
                                                        R const                    &bivar_solver_eps,
                                                        bool                        fst_end_snd_start,
                                                        std::vector<NLM::p3<R>>    &isec_stat_points);

        /* get list of pointers to all neurite paths of the cell network */
        void                                        getAllNeuritePaths(std::list<NLM::NeuritePath<R> *> &neurite_paths);
        void                                        getAllNeuritePaths(std::list<NLM::NeuritePath<R> const *> &neurite_paths) const;

        /* compute all intersection jobs for one full analysis cycle */
        void                                        computeFullAnalysisIntersectionJobs(std::list<std::shared_ptr<IsecJob>> &job_queue) const;

        /* thread-related methods */
        static void                                 startWorkerThread(ThreadInfo *tinfo);
       
        std::vector<ThreadInfo>                     thread_slots; 
        void                                        processIntersectionJobsMultiThreaded(
                                                        uint32_t const                             &nthreads,
                                                        std::list<std::shared_ptr<IsecJob>> const  &job_queue,
                                                        std::list<std::shared_ptr<IsecJob>>        &results);


        /* compute initial neurite root vertices as described in the thesis: all vertices inside the soma sphere are deleted.
         * if the neurite branches before leaving the soma sphere, an exception is thrown, but this case is extremely
         * rare in practice with the current soma sphere settings. */
        bool                                        nlm_neurite_root_vertices_adjusted;
    public: // we need this to be public for scaling radii
        void                                        computeInitialNeuriteRootVertices(); 
    protected:
        /* update basic information specific to NLM_CellNetwork */
        bool                                        nlm_network_info_updated;
        void                                        updateNLMNetworkInfo();

        /* update mdv information of neurite segment / all neurite segments */
        bool                                        updateNsMDVInformation(neurite_segment_iterator ns_it);                                        
    public:
        void                                        updateAllMDVInformation();


        /* cell-partitioning */

    public:
        /* methods returning pre-defined std::function wrapped selection algorithms for use during cell partitioning */
        static std::function<
                typename NLM_CellNetwork<R>::neurite_segment_iterator(
                    NLM_CellNetwork<R> &,
                    typename NLM_CellNetwork<R>::neurite_iterator
                )
            >                                       partition_select_max_chordal_depth(
                                                        R   filter_angle        = (2.0 * M_PI) / 3.0,
                                                        R   max_radius_ratio    = Aux::Numbers::inf<R>());

        static std::function<
                typename NLM_CellNetwork<R>::neurite_segment_iterator(
                    NLM_CellNetwork<R> &,
                    typename NLM_CellNetwork<R>::neurite_iterator
                )
            >                                       partition_select_min_angle(
                                                        R   filter_angle        = (2.0 * M_PI) / 3.0,
                                                        R   max_radius_ratio    = Aux::Numbers::inf<R>());

        static std::function<
                typename NLM_CellNetwork<R>::neurite_segment_iterator(
                    NLM_CellNetwork<R> &,
                    typename NLM_CellNetwork<R>::neurite_iterator
                )
            >                                       partition_select_simple_neurite_paths();

        /* methods returning pre-defined std::function wrapped parametrization algorithms for use during spine curve
         * computation by updateCellGeometry / NeuritePath::updateGeometry(). */
        static std::function<
                std::tuple<
                    R,
                    std::vector<R>,
                    R
                >(NLM::NeuritePath<R> const &P)
            >                                       parametrization_chord_length(),
                                                    parametrization_centripetal(),
                                                    parametrization_uniform();
    protected:
        /* partition an entire cell of the network, identified by the soma s */
        void                                        partitionCell(
                                                        soma_iterator                       s_it,
                                                        std::function<
                                                                neurite_segment_iterator(
                                                                    NLM_CellNetwork<R> &,
                                                                    neurite_iterator
                                                                )
                                                            >  const                       &selection_algorithm);

        /* update cell geometry for the entire cell C_s identified by soma s. this method assumes that C_s has been
         * properly partitioned, e.g. with partitionCell(). if the partitioning is incomplete, the geometry update will
         * be incomplete as well, which might be desirable in certain cases, yet it is not the intent of this method */
        void                                        updateCellGeometry(
                                                        soma_iterator   s_it,
                                                        std::function<
                                                                std::tuple<
                                                                    R,
                                                                    std::vector<R>,
                                                                    R
                                                                >(NLM::NeuritePath<R> const &P)
                                                            >  const                       &parametrization_algorithm);

    public:
                                                    NLM_CellNetwork(std::string network_name);
                                                   ~NLM_CellNetwork();

        /* -----------------------------------  I / O  ------------------------------------------------------------- */
        /* initialize network from standardized SWC file as used by NeuroMorpho.org. aside from a few small NLM specific
         * tasks, this just calls (and overrides) CellNetwork::readFromNeuroMorphoSWCFile(). */
        void                                        readFromNeuroMorphoSWCFile(
                                                        std::string     filename,
                                                        bool const     &check_coincident_positions = false);

        /* partition entire network */
        void                                        partitionNetwork();

        /* update geometry of entire network */
        void                                        updateNetworkGeometry();

        /* perform one full analysis iteration on the entire cell network. */
        bool                                        performFullAnalysis();

        /* mesh generation */
        template <typename Tm, typename Tv, typename Tf>
        void                                        renderCellNetwork(std::string filename);

        template <typename Tm, typename Tv, typename Tf>
        void                                        renderModellingMeshesIndividually(std::string filename) const;

        /* output visualization file for morphview */
        void                                        writeMorphViewFile(std::string filename) const;
};

#include "../tsrc/NLM_CellNetwork_impl.hh"

#endif
