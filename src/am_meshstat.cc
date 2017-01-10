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

#include "common.hh"

#include "Polynomial.hh"

#include "Octree.hh"
#include "Mesh.hh"
#include "MeshAlgorithms.hh"

#include "CanalSurface.hh"

#include "CellNetwork.hh"
#include "NLM_CellNetwork.hh"

std::string const usage_text = 
"--------------------------------------------------------------------------------\n"\
"\n"\
"AnaMorph: A Framework for Geometric Modelling, Consistency Analysis and Surface\n"\
"Mesh Generation of Anatomically Reconstructed Neuron Morphologies\n"\
"\n"\
"Web site: http://www.anamorph.net\n"\
"\n"\
"Copyright (c) 2013-2014, Konstantin Mörschel.\n"\
"All rights reserved.\n"\
"\n"\
"Redistribution and use in source and binary forms, with or without\n"\
"modification, are permitted provided that the following conditions are met:\n"\
"1. Redistributions of source code must retain the above copyright\n"\
"   notice, this list of conditions and the following disclaimer.\n"\
"2. Redistributions in binary form must reproduce the above copyright\n"\
"   notice, this list of conditions and the following disclaimer in the\n"\
"   documentation and/or other materials provided with the distribution.\n"\
"3. All advertising materials mentioning features or use of this software\n"\
"   must display the following acknowledgement:\n"\
"\n"\
"   This product includes software developed by Konstantin Mörschel for\n"\
"   AnaMorph (http://www.anamorph.net).\n"\
"\n"\
"4. Neither the name \"AnaMorph\" nor the names of its contributors may be\n"\
"   used to endorse or promote products derived from this software without\n"\
"   specific prior written permission.\n"\
"\n"\
"THIS SOFTWARE IS PROVIDED BY THE CONTRIBUTORS OF ANAMORPH ''AS IS'' AND ANY\n"\
"EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"\
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\n"\
"DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS OF ANAMORPH BE LIABLE FOR ANY\n"\
"DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n"\
"(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"\
"LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\n"\
"ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"\
"(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n"\
"SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n"\
"\n"\
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
