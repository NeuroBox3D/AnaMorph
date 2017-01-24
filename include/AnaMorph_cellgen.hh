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

#ifndef ANAMORPH_CELLGEN_HH
#define ANAMORPH_CELLGEN_HH

#include "CLApplication.hh"
#include "NLM_CellNetwork.hh"

class AnaMorph_cellgen : private CLApplication {
    private:
        std::string         network_name;
        static const std::list<
                std::pair<
                    std::string,
                    uint32_t
                >
            >               cl_settings_info;

        static const std::list<
                std::pair<
                    std::string,
                    std::string
                >
            >               cl_mutex_switch_list;

        bool                ana;

        std::function<
                NLM_CellNetwork<double>::neurite_segment_iterator(
                    NLM_CellNetwork<double> &,
                    NLM_CellNetwork<double>::neurite_iterator
                )
            >               partition_algo;

        std::function<
                std::tuple<
                    double,
                    std::vector<double>,
                    double
                >(NLM::NeuritePath<double> const &P)
            >               parametrization_algo;

        uint32_t            ana_nthreads;
        double              ana_univar_solver_eps;
        double              ana_bivar_solver_eps;

        bool                meshing;
        bool                force_meshing;
        bool                meshing_individual_surfaces;

        bool                meshing_flush;
        uint32_t            meshing_flush_face_limit;

        uint32_t            meshing_canal_segment_n_phi_segments;
        uint32_t            meshing_outer_loop_maxiter;
        uint32_t            meshing_inner_loop_maxiter;

        bool                meshing_preserve_crease_edges;
        double 				meshing_cansurf_triangle_height_factor;

        double              meshing_radius_factor_decrement;
        double              meshing_complex_edge_max_growth_factor;


        bool                pc;
        double              pc_alpha;
        double              pc_beta;
        double              pc_gamma;

        bool                pp_gec;
        double              pp_gec_alpha;
        double              pp_gec_lambda;
        double              pp_gec_mu;
        uint32_t            pp_gec_d;

        bool                pp_hc;
        double              pp_hc_alpha;
        double              pp_hc_beta;
        uint32_t            pp_hc_maxiter;

                            AnaMorph_cellgen(AnaMorph_cellgen const &) = delete;
                            AnaMorph_cellgen(AnaMorph_cellgen const &&) = delete;
        AnaMorph_cellgen   &operator=(AnaMorph_cellgen const &) = delete;

    public:
                            AnaMorph_cellgen(int argc, char *argv[]);
        virtual				~AnaMorph_cellgen() {};
        virtual bool        processCommandLineArguments() final;
        virtual bool        run() final;
};

#endif
