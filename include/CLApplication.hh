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
