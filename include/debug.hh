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

#ifndef DEBUG_H
#define DEBUG_H

#include "common.hh"

/* --------- debug functionality ---------- */
enum DEBUG_COMPONENT {
    DBG_GLOBAL          = 0,
    DBG_POLYNOMIAL      = 1,
    DBG_POLYSOLVERS     = 2,
    DBG_DMC             = 3,
    DBG_CORE            = 4,
    DBG_NEURITE_PATH    = 5
};

void        initDebug();

uint32_t    getDebugComponent();
void        setDebugComponent(uint32_t comp);
void        setDebugLevel(uint32_t level);
void        setMaxDebugLevel(uint32_t max_level);
void        setDebugTab(uint32_t tab);
uint32_t    getDebugTab();

void        disableComponentDebug(uint32_t comp);
void        enableComponentDebug(uint32_t comp);

void        debugTabIncrement();
void        debugTabDecrement();

int         debugprintf(std::string filename, int line, std::string fmt, ...);
int         debugprintf_wlevel(uint32_t level, std::string filename, int line, std::string fmt, ...);

#ifdef __DEBUG__
    #define debug(...) debugprintf(__FILE__, __LINE__, __VA_ARGS__)
    #define debugl(level, ...) debugprintf_wlevel(level, __FILE__, __LINE__, __VA_ARGS__)
    #define debugTabInc() debugTabIncrement()
    #define debugTabDec() debugTabDecrement()
#else
    #define debug(...)
    #define debugl(level, ...)
    #define debugTabInc()
    #define debugTabDec()
#endif

#endif
