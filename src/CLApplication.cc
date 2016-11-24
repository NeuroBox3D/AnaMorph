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

#include "CLApplication.hh"
#include "aux.hh"

CLApplication::CLApplication(
    int                                                     argc,
    char                                                   *argv[],
    std::list<std::pair<std::string, uint32_t>> const      &cl_settings_info,
    std::list<std::pair<std::string, std::string>> const   &cl_mutex_switch_list,
    std::string const                                      &usage_text)
{
    this->cl_settings_info      = cl_settings_info;
    this->cl_mutex_switch_list  = cl_mutex_switch_list;
    this->usage_text            = usage_text;

    /* parse commandline arguments. arity of switches is given in (this->) cl_settings_info */
    /* parse string pairs (var, value) into list cl_settings */
    std::string                 s;
    std::vector<std::string>    s_args;
    int                         s_arity, i;
    bool                        s_valid;

    i = 1;
    while (i < argc) {
        s       = std::string(argv[i]);
        s_valid = false;

        //printf("i = %d, new switch: \"%s\".\n", i, s.c_str());

        /* s must start with a leading "-" */
        if (s[0] == '-') {
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
