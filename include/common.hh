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

#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cstdint>

#include <limits>

#include <errno.h>
#include <math.h>
#include <complex>

#include <array>
#include <cstring>
#include <ctime>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <algorithm>
// author forgot <numeric> which s sometimes included implicitly
#include <numeric>
#include <stdexcept>

/* stream io */
#include <iostream>
#include <fstream>
#include <ios>

/* debugging and auxiliary functions */
#include "debug.hh"

#ifndef __WIN32__
    #include <unistd.h>
    #include <sys/types.h>
    #include <sys/socket.h>
    #include <netdb.h>
    #include <arpa/inet.h>
    #include <fcntl.h>
#endif

namespace Common {
    /* constants */
#ifdef __WIN32__
        #include <windows.h>
        #define M_PI           3.14159265358979323846  /* pi */
        #define M_PI_2         1.57079632679489661923  /* pi/2 */
        #define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif 

    const double    twopi  = 2.0*M_PI;
    const double    sqrt2  = sqrt(2.0);

    //const double    INF    = std::numeric_limits<double>::infinity();


/* NOTE: clang seems to use UINT32_MAX from stdint.h by default, g++ doesn't however. this macro will work.. */
#ifndef UINT32_MAX
    #define UINT32_MAX std::numeric_limits<uint32_t>::max()
#endif

    /* unit type (of non-empty size, although it is empty), which is used as default template parameter here and there
     * */
    class UnitType {};
}

#endif
