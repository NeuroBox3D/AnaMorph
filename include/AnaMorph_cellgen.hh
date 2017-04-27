/*
 * This file is part of
 *
 * AnaMorph: a framework for geometric modelling, consistency analysis and surface
 * mesh generation of anatomically reconstructed neuron morphologies.
 * 
 * Copyright (c) 2013-2017: G-CSC, Goethe University Frankfurt - Queisser group
 * Author: Konstantin Mörschel
 * 
 * AnaMorph is free software: Redistribution and use in source and binary forms,
 * with or without modification, are permitted under the terms of the
 * GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works:
 * "Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph)."
 *
 * (3) Neither the name "AnaMorph" nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * (4) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Mörschel K, Breit M, Queisser G. Generating neuron geometries for detailed
 *   three-dimensional simulations using AnaMorph. Neuroinformatics (2017)"
 * "Grein S, Stepniewski M, Reiter S, Knodel MM, Queisser G.
 *   1D-3D hybrid modelling – from multi-compartment models to full resolution
 *   models in space and time. Frontiers in Neuroinformatics 8, 68 (2014)"
 * "Breit M, Stepniewski M, Grein S, Gottmann P, Reinhardt L, Queisser G.
 *   Anatomically detailed and large-scale simulations studying synapse loss
 *   and synchrony using NeuroBox. Frontiers in Neuroanatomy 10 (2016)"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

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

        uint32_t            meshing_n_soma_refs;
        double              scale_radius;
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
