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

#include "common.hh"

#include "Polynomial.hh"

#include "Octree.hh"
#include "Mesh.hh"
#include "MeshAlgorithms.hh"

#include "CanalSurface.hh"

#include "CellNetwork.hh"
#include "NLM_CellNetwork.hh"

std::string const usage_text = 
"--------------------------------------------------------------------------------\n"
" AnaMorph: a framework for geometric modelling, consistency analysis and surface\n"
" mesh generation of anatomically reconstructed neuron morphologies.\n"
"\n"
" Copyright (c) 2013-2017: G-CSC, Goethe University Frankfurt - Queisser group\n"
" Created by Konstantin Mörschel.\n"
"\n"
" AnaMorph is free software: Redistribution and use in source and binary forms,\n"
" with or without modification, are permitted under the terms of the\n"
" GNU Lesser General Public License version 3 (as published by the\n"
" Free Software Foundation) with the following additional attribution\n"
" requirements (according to LGPL/GPL v3 §7):\n"
"\n"
" (1) The following notice must be displayed in the Appropriate Legal Notices\n"
" of covered and combined works:\n"
" \"Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph).\"\n"
"\n"
" (2) The following notice must be displayed at a prominent place in the\n"
" terminal output of covered works:\n"
" \"Based on AnaMorph (https://github.com/NeuroBox3D/AnaMorph).\"\n"
"\n"
" (3) Neither the name \"AnaMorph\" nor the names of its contributors may be\n"
" used to endorse or promote products derived from this software without\n"
" specific prior written permission.\n"
"\n"
" (4) The following bibliography is recommended for citation and must be\n"
" preserved in all covered files:\n"
" \"Mörschel K, Breit M, Queisser G. Generating neuron geometries for detailed\n"
"   three-dimensional simulations using AnaMorph. Neuroinformatics (2017)\"\n"
" \"Grein S, Stepniewski M, Reiter S, Knodel MM, Queisser G.\n"
"   1D-3D hybrid modelling – from multi-compartment models to full resolution\n"
"   models in space and time. Frontiers in Neuroinformatics 8, 68 (2014)\"\n"
" \"Breit M, Stepniewski M, Grein S, Gottmann P, Reinhardt L, Queisser G.\n"
"   Anatomically detailed and large-scale simulations studying synapse loss\n"
"   and synchrony using NeuroBox. Frontiers in Neuroanatomy 10 (2016)\"\n"
"\n"
" This program is distributed in the hope that it will be useful,\n"
" but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
" See the GNU Lesser General Public License for more details.\n"
" You should have received a copy of the GNU Lesser General Public License\n"
" along with this program. If not, see <http://www.gnu.org/licenses/>.\n"
"--------------------------------------------------------------------------------\n"\
"\n"\
"am_meshstat: generate mesh statistics.\n"\
"\n"\
"Usage: am_meshstat <OBJ_FILE>\n"\
"\n";

using namespace std;

int main(int argc, char *argv[])
{
    std::string meshname;

    if (argc == 2) {
        meshname = std::string(argv[1]);
    }
    else {
        printf("%s", usage_text.c_str());
        return EXIT_FAILURE;
    }
    try {
        double      area, volume, ar_avg, ar_sigma, ar_max;
        uint32_t    nobtuse_tris;
        int         nvertices, nfaces, nedges, chi;

        Mesh<bool, bool, bool, double> M;
        M.readFromObjFile(meshname.c_str());

        /* statistics */
        area            = M.getTotalArea();
        volume          = M.getTotalVolume();
        nvertices       = M.numVertices();
        nfaces          = M.numFaces();
        nedges          = M.numEdges();
        chi             = nvertices - nedges + nfaces;
        M.getAvgAspectRatio(ar_avg, ar_sigma, &ar_max);
        nobtuse_tris    = M.numObtuseTriangles();

        printf("Mesh: \"%s\"\n\n"\
               "nvertices:      %8d\n"\
               "nedges:         %8d\n"\
               "nfaces:         %8d\n"\
               "chi:            %8d\n"\
               "area:           %14.5f\n"\
               "volume:         %14.5f\n"\
               "ar_avg:         %14.5f\n"\
               "ar_sigma:       %14.5f\n"\
               "ar_max:         %14.5f\n"\
               "obtuse tris:    %8d\n",
                meshname.c_str(), nvertices, nedges, nfaces, chi, 
                area, volume, ar_avg, ar_sigma, ar_max,
                nobtuse_tris);

        fflush(stdout);
    }
    catch (const char *err) {
        printf("caught string err: \"%s\". shutting down..\n", err);
    }
    catch (std::string& err) {
        printf("caught string err: \"%s\". shutting down..\n", err.c_str());
    }
    catch (MeshEx& ex) {
        printf("caught MeshEx. error msg: \"%s\". shutting down..\n", ex.error_msg.c_str());
    }
    catch (...) {
        printf("caught unhandled exception. shutting down..\n");
    }
}
