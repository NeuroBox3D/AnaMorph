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

#include "CLApplication.hh"
#include "aux.hh"

CLApplication::CLApplication(
    int                                                     argc,
    char                                                   *argv[],
    const std::list<std::pair<std::string, uint32_t> >&    _cl_settings_info,
    const std::list<std::pair<std::string, std::string> >& _cl_mutex_switch_list,
    const std::string&                                     _usage_text)
: cl_settings_info(_cl_settings_info), cl_mutex_switch_list(_cl_mutex_switch_list),
  usage_text(_usage_text)
{
    /* parse commandline arguments. arity of switches is given in (this->) cl_settings_info */
    /* parse string pairs (var, value) into list cl_settings */
    std::string                 s;
    std::vector<std::string>    s_args;
    int                         s_arity, i;

    i = 1;
    while (i < argc) {
        s       = std::string(argv[i]);

        //printf("i = %d, new switch: \"%s\".\n", i, s.c_str());

        /* s must start with a leading "-" */
        if (s[0] == '-') {
            bool s_valid = false;

            /* remove leading "-" */
            s.erase(s.begin());

            /* check if switch is supported by scanning cl_settings_info */
            for (auto clsi_pair : cl_settings_info) {
                if (s == clsi_pair.first) {
                    s_arity     = clsi_pair.second;
                    s_valid     = true;
                    break;
                }
            }

            if (s_valid) {
                /* get next arity arguments if possible and generate cl_settings element */
                if (argc - i > s_arity) {
                    s_args.clear();
                    for (int j = 1; j <= s_arity; j++) {
                        s_args.push_back(argv[i + j]);
                    }
                    this->cl_settings.push_back({ s, s_args });
                    i += s_arity + 1;
                }
                else {
                    throw CLAEx(CLA_INVALID_CLARGS, "ERROR: not enough arguments for command line switch \"" + s + "\".\n");
                }
            }
            else {
                throw CLAEx(CLA_INVALID_CLARGS, "ERROR: invalid command line switch \"" + s + "\".\n");
            }
        }
        else {
            throw CLAEx(CLA_INVALID_CLARGS, "ERROR: syntax error. expected new command line switch in place of \"" + s + "\".\n");
        }
    }

    /* sort command line info by switch name and check for double entries */
    this->cl_settings.sort(
        [] (
            std::pair<std::string, std::vector<std::string>> const &x,
            std::pair<std::string, std::vector<std::string>> const &y) -> bool
        {
            return (x.first < y.first);
        }
    );

    uint32_t cls_old_size = cl_settings.size();
    this->cl_settings.unique(
        [] (
            std::pair<std::string, std::vector<std::string>> const &x,
            std::pair<std::string, std::vector<std::string>> const &y) -> bool
        {
            return (x.first == y.first);
        }
    );
    
    if (cls_old_size != cl_settings.size()) {
        throw CLAEx(CLA_INVALID_CLARGS, "ERROR: no command line switch shall be supplied more than once.\n");
    }

    /* check if mutually exclusive switches have been set */
    std::list<std::string> switches;
    for (auto &cls_pair : this->cl_settings) {
        //printf("parsed switch: \"%s\"\n", cls_pair.first.c_str() );
        switches.push_back(cls_pair.first);
    }

    for (auto &mutexswitch_pair : this->cl_mutex_switch_list) {
        if (Aux::Alg::listContains(switches, mutexswitch_pair.first) &&
            Aux::Alg::listContains(switches, mutexswitch_pair.second))
        {
            throw CLAEx(CLA_INVALID_CLARGS, "ERROR: mutually exclusive command line switches \"" + mutexswitch_pair.first + "\" and \"" + mutexswitch_pair.second + "\" used simultaneously.\n");
        }
    }
}
