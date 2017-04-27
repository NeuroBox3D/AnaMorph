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

#ifndef CL_APPLICATION_HH
#define CL_APPLICATION_HH

#include "common.hh"

enum CLAEx_Types {
    CLA_INVALID_CLARGS
};

struct CLAEx : public std::runtime_error {
    const uint32_t      error_type;
    const std::string   error_msg;

    CLAEx(
        const uint32_t &type,
        const std::string &msg) :
            std::runtime_error(msg),
            error_type(type),
            error_msg(msg)
    {
    }
};

class CLApplication {
    protected:
        std::list<std::pair<std::string, uint32_t>>                     cl_settings_info;
        std::list<std::pair<std::string, std::vector<std::string>>>     cl_settings;
        std::list<std::pair<std::string, std::string>>                  cl_mutex_switch_list;
        std::string                                                     usage_text;

        /* non-copy-constructive, non-assignable, non-movable */
                        CLApplication(CLApplication const &x)   = delete;
                        CLApplication(CLApplication const &&x)  = delete;
        CLApplication  &operator=(CLApplication const &x)       = delete;

    public:
                        CLApplication(
                            int                                                     argc,
                            char                                                   *argv[],
                            std::list<std::pair<std::string, uint32_t>> const      &cl_settings_info,
                            std::list<std::pair<std::string, std::string>> const   &cl_mutex_switch_list,
                            std::string const                                      &usage_text = "");

        virtual			~CLApplication() {};

        /* processing command line arguments and calling main loop is specific to application => pure virtual */
        virtual bool    processCommandLineArguments() = 0; 
        virtual bool    run() = 0;
};

#endif
