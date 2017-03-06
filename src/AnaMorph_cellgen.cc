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
 *    This product includes software developed by f for
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

#include "AnaMorph_cellgen.hh"
#include "NLM_CellNetwork.hh"
#include "CellNetworkAlg.hh"
#include "MeshAlgorithms.hh"

/* initialize static command line switch info */
const std::list<
    std::pair<
        std::string,
        uint32_t
    >
> AnaMorph_cellgen::cl_settings_info = 
    { 
        { "h",                                      0 },
        { "i",                                      1 },
        { "analysis",                               0 },
        { "no-analysis",                            0 },
        { "meshing",                                0 },
        { "no-meshing",                             0 },
        { "force-meshing",                          0 },
        { "meshing-individual-surfaces",            0 },
        { "cellnet-pc",                             3 },
        { "no-cellnet-pc",                          0 },
        { "scale-radius",                           1 },
        { "cellnet-partition-strategy",             3 },
        { "cellnet-parametrization-strategy",       1 },
        { "ana-nthreads",                           1 },
        { "ana-univar-eps",                         1 },
        { "ana-bivar-eps",                          1 },
        { "no-mesh-pp",                             0 },
        { "mesh-pp-gec",                            4 },
        { "no-mesh-pp-gec",                         0 },
        { "mesh-pp-hc",                             3 },
        { "no-mesh-pp-hc",                          0 },
        { "meshing-soma-refs",                      1 },
        { "meshing-cansurf-angularsegments",        1 },
        { "meshing-triangle-height",                1 },
        { "meshing-outerloop-maxiter",              1 },
        { "meshing-innerloop-maxiter",              1 },
        { "preserve-crease-edges",                  0 },
        { "meshing-flush",                          1 },
        { "no-meshing-flush",                       0 },
        { "meshing-merging-initial-radiusfactor",   1 },
        { "meshing-merging-radiusfactor-decrement", 1 },
        { "meshing-complexedge-max-growthfactor",   1 },
        { "debug-lvl",                              2 }
    };

const std::list<
    std::pair<
        std::string,
        std::string
    >
> AnaMorph_cellgen::cl_mutex_switch_list = 
    {
        { "analysis",       "no-analysis" },
        { "meshing",        "no-meshing" },
        { "cellnet-pc",     "no-cellnet-pc" },
        { "mesh-pp-gec",    "no-mesh-pp-gec" },
        { "mesh-pp-hc",     "no-mesh-pp-hc" },
        { "no-mesh-pp",     "mesh-pp-gec"},
        { "no-mesh-pp",     "mesh-pp-hc"},
        { "meshing-flush",  "no-meshing-flush" },
        { "no-analysis",    "meshing" },
        { "no-analysis",    "force-meshing" },
        { "no-analysis",    "meshing-individual-surfaces"           },
        { "no-analysis",    "meshing-cansurf-angularsegments",      },
        { "no-analysis",    "meshing-outerloop-maxiter",            },
        { "no-analysis",    "meshing-innerloop-maxiter",            },
        { "no-analysis",    "meshing-flush",                        },
        { "no-analysis",    "no-meshing-flush",                     },
        { "no-analysis",    "meshing-merging-initial-radiusfactor", },
        { "no-analysis",    "meshing-merging-radiusfactor-decrement"},
        { "no-analysis",    "meshing-complexedge-max-growthfactor"  },
    };

const std::string usage_string = 
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
"am_cellgen: command line interface for the non-linear geometric modelling and\n"\
"analysis approach.\n"\
"\n"\
"Usage: am_cellgen -i <NETWORKNAME> [OPTIONS]\n"\
"\n"\
"Command line switches:\n"\
" -i <NETWORKNAME>               MANDATORY: specify input cell network name. \n"\
"                                the network name is the name of the SWC\n"\
"                                morphology file without suffix (that is,\n"\
"                                without \".swc\").\n"\
"                                for convenience, the extensions \".swc\", \".obj\"\n"\
"                                and \".amv\" are stripped from the input\n"\
"                                string if present.\n"\
"                                EXAMPLE: the network name for \"ri05.CNG.swc\"\n"\
"                                if \"ri05.CNG\". using \"ri05.CNG.obj\" would yield\n"\
"                                the same network name.\n"\
"\n"\
" -analysis                  \n"\
" -no-analysis                   enable / disable geometric analysis. once a\n"\
"                                single full geometric analysis run completes,\n"\
"                                a result visualization file suitable for the\n"\
"                                AnaMorph morphology viewer (which allows the\n"\
"                                user to view the morphology and potential\n"\
"                                geometric errors three-dimensionally), is\n"\
"                                written to the file\n"\
"\n"\
"                                \"<CELLNETWORK>.amv\"\n"\
"\n"\
"                                NOTE: if cell network preconditioning is enabled\n"\
"                                (-cellnet-pc), the preconditioned cell network\n"\
"                                is analysed, not the originally encoded one.\n"\
"                                DEFAULT: enabled.\n"\
"\n"\
" -ana-nthreads <n>              number of worker threads used during geometric\n"\
"                                analysis. work is distributed approximately\n"\
"                                evenly among all threads. if the host CPU\n"\
"                                supports hyper-threading, a good choice is 2*n,\n"\
"                                otherwise n, where n is the number of physical\n"\
"                                cores of the CPU. n must be > 0.\n"\
"                                DEFAULT: 1. \n"\
"\n"\
" -ana-univar-eps <eps>          floating point tolerance value used for the\n"\
"                                numerical solver (Bezier Clipping) for\n"\
"                                univariate polynomial root finding problems. \n"\
"                                <eps> must be in [1E-11, 1E-3].\n"\
"                                DEFAULT: 1E-6.\n"\
"\n"\
" -ana-bivar-eps <eps>           floating point tolerance value used for the\n"\
"                                numerical solver (Bivariate Linear Clipping) for\n"\
"                                bivariate polynomial root finding problems. \n"\
"                                <eps> must be in [1E-11, 1E-3].\n"\
"                                DEFAULT: 1E-4.\n"\
"\n"\
" -cellnet-pc <alpha> <beta> <gamma>\n"\
" -no-cellnet-pc\n"\
"                                enable / disable cell network preconditioning.\n"\
"                                requires \"-analysis\" to be effective.\n"\
"                                if enabled: <alpha>, <beta> and <gamma> are the\n"\
"                                parameters for the (preliminary) preconditioning\n"\
"                                algorithm for cell networks, all of which must\n"\
"                                be > 0.0. the preconditioned cell network is \n"\
"                                subjected to geometric analysis.\n"\
"                                DEFAULT: enabled, <alpha> = 3.0, <beta> = 1.5\n"\
"                                <gamma> = 10.0.\n"\
"\n"\
" -cellnet-partition-strategy <strategy> <filter_angle> <filter_radius_ratio>\n"\
"                                specify cell network partitioning strategy.\n"\
"                                possible values for <strategy> (without quotes):\n"\
"\n"\
"                                    1. \"max-chordal-depth\": continue neurite\n"\
"                                    path with the neurite segment (u, v) whose\n"\
"                                    end neurite vertex v roots the sub-tree\n"\
"                                    of maximum chordal depth among all possible\n"\
"                                    choices. the chordal depth of a neurite sub-\n"\
"                                    tree with root v is the maximum chord length\n"\
"                                    of all paths from v to a neurite terminal\n"\
"                                    vertex (i.e. leaf).\n"\
"\n"\
"                                    2. \"min-angle\": continue neurite path with\n"\
"                                    the neurite segment of minimum angle.\n"\
"\n"\
"                                    3. \"simple-neurite-paths\": neurite path\n"\
"                                    always ends at the next neurite branching\n"\
"                                    vertex.\n"\
"\n"\
"                                the parameters <filter_angle>, given in degrees\n"\
"                                from [0, 180], and <filter_radius_ratio>, from \n"\
"                                [1.0, oo(, are used to narrow the set of\n"\
"                                possible neurite segments to choose from:\n"\
"                                among all possible neurite segments, only those\n"\
"                                having both an angle less than or equal to\n"\
"                                <filter_angle> and a radius radio of less than\n"\
"                                or equal to <filter_radius_ratio> are\n"\
"                                taken into consideration. if this set is empty,\n"\
"                                the neurite path ends at the terminal vertex.\n"\
"\n"\
"                                DEFAULT: <strategy> = \"max-chordal-depth\",\n"\
"                                         <filter_angle> = \"90.0\",\n"\
"                                         <filter_radius_ratio> = \"10.0\".\n"\
"\n"\
" -cellnet-parametrization-strategy <strategy>\n"\
"                                specify neurite path parametrization strategy.\n"\
"                                possible values for <strategy> (without quotes):\n"\
"\n"\
"                                    1. \"chord-length\": use chord-length\n"\
"                                    parametrization and matching scaling\n"\
"                                    factors.\n"\
"\n"\
"                                    2. \"uniform\": use uniform parametrization\n"\
"                                    and matching scaling factors.\n"\
"\n"\
"                                DEFAULT: <strategy> = \"chord-length\".\n"\
"\n"\
" -meshing                        \n"\
" -no-meshing                    enable / disable inductive cell meshing, which\n"\
"                                generates a consistent triangular output mesh\n"\
"                                representing the plasma membrane of the input\n"\
"                                cell network (the cell network union mesh).\n"\
"                                can be enabled only if geometric analysis is\n"\
"                                enabled as well (-analysis). inductive cell\n"\
"                                meshing is performed only if geometric analysis\n"\
"                                has confirmed the input cell network as clean.\n"\
"                                in order to attempt meshing of unclean networks,\n"\
"                                use -force-meshing. the mesh output file name is\n"\
"\n"\
"                                \"<CELLNETWORK>.obj\"\n"\
"                                \n"\
"                                DEFAULT: enabled.\n"\
"\n"\
" -force-meshing                 force inductive cell meshing, even for unclean\n"\
"                                networks. can be used only if geometric\n"\
"                                analysis is enabled as well (-analysis). \n"\
"                                DEFAULT: don't force meshing.\n"\
"\n"\
" -meshing-individual-surfaces   render all modelling surfaces individually\n"\
"                                into a generally inconsistent mesh.\n"\
"                                no merging is performed: every modelling\n"\
"                                surface produces a connected component of the\n"\
"                                resulting mesh. can be enabled only if geometric\n"\
"                                analysis is enabled as well (-analysis).\n"\
"                                the mesh output file name is\n"\
"\n"\
"                               \"<CELLNETWORK>_individual_modelling_surfaces.obj\"\n"\
"\n"\
"                                DEFAULT: no individual surface meshing.\n"\
"\n"\
" -meshing-flush <flush_face_limit>                    \n"\
" -no-meshing-flush                        \n"\
"                                enable / disable mesh flushing.\n"\
"                                requires -meshing and -analysis.\n"\
"                                during inductive cell network meshing, large\n"\
"                                parts of the partially completed cell network\n"\
"                                mesh are no longer needed and can already be\n"\
"                                flushed to disk to save significant amounts of\n"\
"                                RAM and computation time. the flushing process\n"\
"                                incrementally completes the mesh obj file\n"\
"                                (\"<CELLNETWORK>.obj\"). if enabled, the integral\n"\
"                                parameter <flush_face_limit>, which must be in\n"\
"                                [1024,oo[, defines the face count limit that\n"\
"                                triggers another partial flushing to disk. that\n"\
"                                is, if the partially completed mesh reaches\n"\
"                                the face count <flush_face_limit>, all \n"\
"                                definitely no longer needed mesh components are\n"\
"                                identified and written to disk.\n"\
"                                if disabled, the entire mesh is kept in RAM.\n"\
"                                NOTE: for large cell networks, this may result\n"\
"                                in intolerably high memory consumption, which\n"\
"                                can cause the system to kill the application\n"\
"                                or even freeze (occurred on win32) if the\n"\
"                                amount of available RAM is exceeded.\n"\
"                                DEFAULT: enabled, <flush_face_limit> = 100000.\n"\
"\n"\
" -meshing-soma-refs <n>         defines the number of refinements performed on an\n"
"                                icosahedron to represent the soma sphere,\n"
"                                default value: 3.\n"
"\n"
" -scale-radius <scale>          Defines a scaling factor used to scale the radius\n"
"                                after pre-conditioning. Can be used to create a\n"
"                                cell-in-cell representation of the ER.\n"
"                                default value: 1.0.\n"
"\n"
" -meshing-cansurf-angularsegments <nseg>      \n"\
"                                defines angular discretization precision used\n"\
"                                during neurite path / neurite canal segment \n"\
"                                meshing. during canal surface meshing, circles\n"\
"                                are approximated by ring polygons with <nseg>\n"\
"                                corners on the circle, where all <nseg> sides\n"\
"                                are chosen equilaterally. <nseg> must be\n"\
"                                in [3, 36].\n"\
"                                NOTE: vertex / face count of canal surfaces are\n"\
"                                O(<nseg>^2): choosing high values for <nseg>\n"\
"                                will produce meshes of correspondingly high\n"\
"                                vertex / face count.\n"\
"                                DEFAULT: <nseg> = 6.\n"\
"\n"\
" -meshing-triangle-height <h>   Defines the factor by which the optimal height of\n"
"                                canal surface triangles is multiplied. This argument\n"
"                                is useful to reduce the number of faces in the output\n"
"                                grid - at the cost of decreasing element quality.\n"
"                                DEFAULT: <h> = 1.0.\n"
"                                \n"
" -meshing-outerloop-maxiter <n> defines the maximum number of iterations\n"\
"                                performed in the outer meshing loop before\n"\
"                                the start radius factor for the currently \n"\
"                                processed neurite path mesh is decreased\n"\
"                                by the decrement defined by \n"\
"                                -meshing-merging-radiusfactor-decrement.\n"\
"                                <n> must be positive.\n"\
"                                DEFAULT: <n> = 16.\n"\
"                                 \n"\
" -meshing-innerloop-maxiter <n> defines the maximum number of iterations\n"\
"                                performed in the inner meshing loop before\n"\
"                                breaking the latter and re-randomizing in\n"\
"                                a fresh iteration of the outer meshing loop.\n"\
"                                <n> must be positive.\n"\
"                                DEFAULT: <n> = 8.\n"\
"                                \n"\
" -preserve-crease-edges         If this option is passed then the triangles along\n"
"                                the dendritic and axonal elements will be arranged\n"
"                                in such a way that there is no torsion in the neurites.\n"
"                                This prevents concavities along the neurites, but also\n"
"                                means the triangles in the neurites will no longer\n"
"                                be perfectly equilateral.\n"
"                                \n"
" -meshing-merging-radiusfactor-decrement <d>\n"\
"                                defines the decrement applied to the radius\n"\
"                                factor if \n"\
"                                    i) one of certain geometric invariants\n"\
"                                    used by the Red-Blue-Algorithm is violated,\n"\
"                                    which indicates that the radius factor must\n"\
"                                    be further reduced. \n"\
"\n"\
"                                    ii) the number of outer meshing loop\n"\
"                                    iterations for the currently processed\n"\
"                                    neurite path is a multiple of the value\n"\
"                                    set with -meshing-outerloop-maxiter.\n"\
"\n"\
"                                if f denotes the current factor, then the new\n"\
"                                factor f' = f - <d>. <d> must be a floating\n"\
"                                point value in ]0.0, 0.1] (max. 10% reduction in\n"\
"                                one step).\n"\
"                                DEFAULT: <d> = 0.01.\n"\
"\n"\
" -meshing-complexedge-max-growthfactor <c>\n"\
"                                defines the maximum growth factor for complexly\n"\
"                                intersecting mesh edges during merging.\n"\
"                                complexly intersecting edges must be split to\n"\
"                                ensure geometric and topological invariants used\n"\
"                                by the Red-Blue meshing algorithm.\n"\
"                                in certain situations, splitting complex edges\n"\
"                                creates even more complex edges, a process that\n"\
"                                can lead to exponential amplification, resulting\n"\
"                                in extreme edge counts and numerical problems.\n"\
"                                the maximum growth factor <c> is used to avoid\n"\
"                                such cases and trigger another iteration of the\n"\
"                                outer meshing loop. if m denotes the initial\n"\
"                                count of complexly intersecting edges of the\n"\
"                                original red / blue mesh pair (R, B),\n"\
"                                then _every new_ splitting step is allowed to\n"\
"                                produce at most <c>*m new complexly intersecting\n"\
"                                edges in the process of splitting and thus removing\n"\
"                                all old complexly intersecting edges.\n"\
"                                splitting continues until either the maximum\n"\
"                                number of inner meshing iterations is reached,\n"\
"                                another (possibly more severe) exception is\n"\
"                                thrown, the limit <c>*m is exceeded or until\n"\
"                                there are no more complexly intersecting edges\n"\
"                                left in both red and blue mesh.\n"\
"                                oftentimes, the splitting process dies out after\n"\
"                                a few iterations. in relatively rare cases\n"\
"                                however, it is necessary to stop a runaway\n"\
"                                process. <c> must be > 0.0, but a value close to\n"\
"                                1 is advisable in light of the above explanation.\n"\
"                                DEFAULT: <c> = 2.0.\n"\
"                                \n"\
" -no-mesh-pp                    disable cell network union mesh post-processing\n"\
"                                entirely. equivalent to\n"\
"                                -no-mesh-pp-gec -no-mesh-pp-hc\n"\
"\n"\
" -mesh-pp-gec <alpha> <lambda> <mu> <n>\n"\
" -no-mesh-pp-gec                enable / disable stage 1 of the post-processing\n"\
"                                chain for the cell network union mesh: improved\n"\
"                                greedy edge collapsing.\n"\
"                                <alpha>, <lambda> and <mu> are the parameters\n"\
"                                passed to the improved greedy edge collapse \n"\
"                                algorithm. <alpha> is a triangle aspect ratio\n"\
"                                and must hence be >= 1.0. <lambda> and <mu>\n"\
"                                are area scaling factors and must be > 0.0,\n"\
"                                n is the depth of the search for neighbors and\n"\
"                                must be >1.\n"\
"\n"\
"                                DEFAULT: enabled, <alpha> = 1.5,\n"\
"                                <lambda> = 0.125, <mu> = 0.5, <n> = 5.\n"\
"\n"\
"                                NOTE: if either -mesh-pp-gec (stage 1) or \n"\
"                                -mesh-pp-hc (stage2) is enabled, post-processing\n"\
"                                is performed as configured, where stage 1 always\n"\
"                                precedes stage 2. the mesh obtained after post-\n"\
"                                processing is written to the file\n"\
"\n"\
"                                \"<CELLNETWORK>_post_processed.obj\".\n"\
"\n"\
" -mesh-pp-hc <alpha> <beta> <maxiter>\n"\
" -no-mesh-pp-hc                 enable / disable stage 2 of the post-processing\n"\
"                                chain for the cell network union mesh:\n"\
"                                HC Laplacian smoothing.\n"\
"                                the algorithm is performed for <maxiter>\n"\
"                                (which must be > 0) iterations. <alpha> and\n"\
"                                <beta> are the key parameters for HC Laplacian\n"\
"                                smoothing, where <beta> must be greater than\n"\
"                                <alpha> to ensure a converging \"diffusion\",\n"\
"                                i.e. smoothing, process. \n"\
"                                DEFAULT: enabled, alpha=0.4, beta=0.7, maxiter=10\n"\
"\n"\
"                                NOTE: if either -mesh-pp-gec (stage 1) or \n"\
"                                -mesh-pp-hc (stage2) is enabled, post-processing\n"\
"                                is performed as configured, where stage 1 always\n"\
"                                precedes stage 2. the mesh obtained after post-\n"\
"                                processing is written to the file\n"\
"\n"\
"                                \"<CELLNETWORK>_post_processed.obj\".\n"\
"\n"\
" -debug-lvl <cmp> <lvl>         Enable debugging for component <cmp>\n"\
"                                and set debug level to <lvl>.\n"\
"                                Debug component 0 is global debugging.\n"\
"\n"\
" -h                             display this help message.\n";



AnaMorph_cellgen::AnaMorph_cellgen(
    int     argc,
    char   *argv[])
        : CLApplication(
            argc,
            argv,
            AnaMorph_cellgen::cl_settings_info,
            AnaMorph_cellgen::cl_mutex_switch_list,
            usage_string)
{
    /* default settings */
    this->ana                                       = true;
    this->ana_nthreads                              = 1;
    this->ana_univar_solver_eps                     = 1E-6;
    this->ana_bivar_solver_eps                      = 1E-4;

    this->partition_algo                            = NLM_CellNetwork<double>::partition_select_max_chordal_depth(
                                                          M_PI / 2.0,
                                                          10.0);
    this->parametrization_algo                      = NLM_CellNetwork<double>::parametrization_chord_length();

    this->meshing                                   = true;
    this->force_meshing                             = false;
    this->meshing_individual_surfaces               = false;

    this->meshing_flush                             = true;
    this->meshing_flush_face_limit                  = 100000;

    this->meshing_n_soma_refs                       = 3;
    this->scale_radius                              = 1.0;
    this->meshing_canal_segment_n_phi_segments      = 6;
    this->meshing_outer_loop_maxiter                = 16;
    this->meshing_inner_loop_maxiter                = 8;

    this->meshing_preserve_crease_edges             = false;
    this->meshing_cansurf_triangle_height_factor	= 1.0;

    this->meshing_radius_factor_decrement           = 0.01;
    this->meshing_complex_edge_max_growth_factor    = 2.0;

    this->pc                                        = true;
    this->pc_alpha                                  = 3.0;
    this->pc_beta                                   = 1.5;
    this->pc_gamma                                  = 10.0;

    this->pp_gec                                    = true;
    this->pp_gec_alpha                              = 1.5;
    this->pp_gec_lambda                             = 0.125;
    this->pp_gec_mu                                 = 0.5;
    this->pp_gec_d                                  = 5;

    this->pp_hc                                     = true;
    this->pp_hc_alpha                               = 0.4;
    this->pp_hc_beta                                = 0.7;
    this->pp_hc_maxiter                             = 10;
}

bool
AnaMorph_cellgen::processCommandLineArguments()
{
    using Aux::Alg::stou;

    if (this->cl_settings.empty()) {
        printf("%s", this->usage_text.c_str());
        return false;
    }

    else for (auto &cls_pair : this->cl_settings) {
        std::string &s                      = cls_pair.first;
        std::vector<std::string> &s_args    = cls_pair.second;

        if (s == "h") {
            printf("%s", this->usage_text.c_str());
            return false;
        }
        if (s == "i") {
            this->network_name = s_args.front();
        }
        else if (s == "analysis") {
            this->ana = true;
        }
        else if (s == "no-analysis") {
            this->ana = false;
        }
        else if (s == "meshing") {
            this->meshing = true;
        }
        else if (s == "no-meshing") {
            this->meshing = false;
        }
        else if (s == "force-meshing") {
            this->force_meshing = true;
        }
        else if (s == "meshing-individual-surfaces") {
            this->meshing_individual_surfaces = true;
        }
        else if (s == "cellnet-pc") {
            try {
                this->pc        = true;
                this->pc_alpha  = std::stod(s_args[0]);
                this->pc_beta   = std::stod(s_args[1]);
                this->pc_gamma  = std::stod(s_args[2]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: at least one argument to switch \"cellnet-pc\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: arguments to switch \"cellnet-pc\" could not be converted to double precision floating point values.\n");
                return false;
            }

            /* check values */
            if (this->pc_alpha <= 0.0 || this->pc_beta <= 0.0 || this->pc_gamma <= 0.0) {
                printf("ERROR: at least one argument to switch \"cellnet-pc\" invalid. alpha, beta and gamma must be in ]0, oo(.\n");
                return false;
            }
        }
        else if (s == "no-cellnet-pc") {
            this->pc = false;
        }

        else if (s == "scale-radius") {
            try {
                scale_radius = std::stod(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"scale-radius\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: argument to switch \"scale-radius\" could not be converted to a double.\n");
                return false;
            }

            // check value
            if (scale_radius < 0.1) {
                printf("ERROR: Radius scale must not be less than 0.1 as this would result in too big a geometry.\n");
                return false;
            }
        }

        else if (s == "cellnet-partition-strategy") {
            double filter_angle, max_radius_ratio;
            try {
                filter_angle        = (std::stod(s_args[1]) / 180.0) * M_PI;
                max_radius_ratio    = std::stod(s_args[2]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: filter angle or maximum radius ratio argument to switch \"cellnet-partition-strategy\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to \"cellnet-partition-strategy\".\n");
                return false;
            }

            /* check filter_angle and max_radius_ratio range */
            if (filter_angle < 0.0 || filter_angle > M_PI) {
                printf("ERROR: filter angle argument to switch \"cellnet-partition-strategy\" must be from interval [0, 180].\n");
                return false;
            }

            if (max_radius_ratio < 1.0) {
                printf("ERROR: maximum radius ratio argument to switch \"cellnet-partition-strategy\" must be from interval [1.0, oo(.\n");
                return false;
            }

            if (s_args[0] == "max-chordal-depth") {
                this->partition_algo = NLM_CellNetwork<double>::partition_select_max_chordal_depth(filter_angle, max_radius_ratio);
            }
            else if (s_args[0] == "min-angle") {
                this->partition_algo = NLM_CellNetwork<double>::partition_select_min_angle(filter_angle, max_radius_ratio);
            }
            else if (s_args[0] == "simple-neurite-paths") {
                this->partition_algo = NLM_CellNetwork<double>::partition_select_simple_neurite_paths();
            }
            else {
                printf("ERROR: argument to switch \"cellnet-partition-strategy\" invalid. possible choices: \"max-chordal-depth\" (default), \"min-angle\" and \"simple-neurite-paths\".\n");
                return false;
            }
        }
        else if (s == "cellnet-parametrization-strategy") {
            if (s_args[0] == "chord-length") {
                this->parametrization_algo = NLM_CellNetwork<double>::parametrization_chord_length();
            }
            else if (s_args[0] == "uniform") {
                this->parametrization_algo = NLM_CellNetwork<double>::parametrization_uniform();
            }
            else {
                printf("ERROR: argument to switch \"cellnet-parametrization-strategy\" invalid. possible choices: \"chord-length\" (default), \"uniform\".\n");
                return false;
            }
        }
        else if (s == "ana-nthreads") {
            try {
                this->ana_nthreads = stou(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"ana-nthreads\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: argument to switch \"ana-nthreads\" could not be converted to an unsigned integer.\n");
                return false;
            }

            /* check value */
            if (this->ana_nthreads == 0) {
                printf("ERROR: number of analysis threads must be >= 1\n");
                return false;
            }
        }
        else if (s == "ana-univar-eps") {
            try {
                this->ana_univar_solver_eps = std::stod(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"ana-univar-eps\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: argument to switch \"ana-univar-eps\" could not be converted to a double precision floating point value.\n");
                return false;
            }

            /* check value */
            if (this->ana_univar_solver_eps < 1E-11 || this->ana_univar_solver_eps > 1E-3) {
                printf("ERROR: invalid argument to switch \"ana-univar-eps\". univariate solver tolerance must be in [1E-3, 1E-11].\n");
                return false;
            }
        }
        else if (s == "ana-bivar-eps") {
            try {
                this->ana_bivar_solver_eps = std::stod(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"ana-bivar-eps\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: argument to switch \"ana-bivar-eps\" could not be converted to a double precision floating point value.\n");
                return false;
            }
            if (this->ana_bivar_solver_eps < 1E-11 || this->ana_bivar_solver_eps > 1E-3) {
                printf("ERROR: invalid argument to switch \"ana-bivar-eps\". bivariate solver tolerance must be in [1E-3, 1E-11].\n");
                return false;
            }
        }
        else if (s == "no-mesh-pp") {
            this->pp_gec    = false;
            this->pp_hc     = false;
        }
        else if (s == "mesh-pp-gec") {
            try {
                this->pp_gec        = true;
                this->pp_gec_alpha  = std::stod(s_args[0]); 
                this->pp_gec_lambda = std::stod(s_args[1]); 
                this->pp_gec_mu     = std::stod(s_args[2]); 
                this->pp_gec_d      = stou(s_args[3]); 
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: at least one argument to switch \"mesh-pp-gec\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to switch \"mesh-pp-gec\".\n");
                return false;
            }

            /* check values */
            if (this->pp_gec_alpha < 1.0) {
                printf("ERROR: \"alpha\" parameter to switch \"mesh-pp-gec\" is a triangle aspect ratio and must be in [1.0, oo(.\n");
                return false;
            }

            if (this->pp_gec_lambda <= 0.0 || this->pp_gec_mu <= 0.0) {
                printf("ERROR: \"lambda\" and \"mu\" parameters to switch \"mesh-pp-gec\" are relative area factors and must be in ]0.0, oo(.\n");
                return false;
            }

            if (this->pp_gec_d <= 1) {
                printf("ERROR: face neighbourhood size parameter \"d\" to switch \"mesh-pp-gec\" must be > 1.\n");
                return false;
            }
        }
        else if (s == "no-mesh-pp-gec") {
            this->pp_gec = false;
        }
        else if (s == "mesh-pp-hc") {
            try {
                this->pp_hc         = true;
                this->pp_hc_alpha   = std::stod(s_args[0]); 
                this->pp_hc_beta    = std::stod(s_args[1]); 
                this->pp_hc_maxiter = stou(s_args[2]); 
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: at least one argument to switch \"mesh-pp-hc\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to switch \"mesh-pp-hc\".\n");
                return false;
            }

            if (this->pp_hc_alpha < 0.0) {
                printf("ERROR: \"alpha\" parameter to switch \"mesh-pp-hc\" must be > 0.\n"); 
                return false;
            }

            if (this->pp_hc_beta <= this->pp_hc_alpha) {
                printf("ERROR: \"beta\" parameter (%5.4f) to switch \"mesh-pp-hc\" must be greater than alpha parameter (%5.4f).\n", this->pp_hc_beta, this->pp_hc_alpha); 
                return false;
            }

            if (this->pp_hc_maxiter == 0) {
                printf("ERROR: \"maxiter\" parameter to switch \"mesh-pp-hc\" must be > 0.\n");
                return false;
            }
        }
        else if (s == "no-mesh-pp-hc") {
            this->pp_hc = false;
        }
        else if (s == "meshing-soma-refs") {
            try {
                meshing_n_soma_refs = stou(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"meshing-soma-refs\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: argument to switch \"meshing-soma-refs\" could not be converted to an unsigned integer.\n");
                return false;
            }

            // check value
            if (this->meshing_canal_segment_n_phi_segments > 10) {
                printf("ERROR: number of soma refinements parameter specified in switch \"meshing-soma-refs\" must be in [0,10].\n");
                return false;
            }
        }
        else if (s == "meshing-cansurf-angularsegments") {
            try {
                this->meshing_canal_segment_n_phi_segments = stou(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"meshing-cansurf-angularsegments\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: argument to switch \"meshing-cansurf-angularsegments\" could not be converted to an unsigned integer.\n");
                return false;
            }

            /* check value */
            if (this->meshing_canal_segment_n_phi_segments < 3 || this->meshing_canal_segment_n_phi_segments > 64) {
                printf("ERROR: angular canal surface discretization parameter specified in switch \"meshing-cansurf-angularsegments\" must be in [12,36].\n");
                return false;
            }
        }
        else if (s == "meshing-outerloop-maxiter") {
            try {
                this->meshing_outer_loop_maxiter = stou(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"meshing-outerloop-maxiter\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to switch \"meshing-outerloop-maxiter\".\n");
                return false;
            }

            /* check value */
            if (this->meshing_outer_loop_maxiter == 0) {
                printf("ERROR: iteration parameter to switch \"meshing-outerloop-maxiter\" must be > 0.\n");
                return false;
            }
        }
        else if (s == "meshing-innerloop-maxiter") {
            try {
                this->meshing_inner_loop_maxiter = stou(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"meshing-innerloop-maxiter\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to switch \"meshing-innerloop-maxiter\".\n");
                return false;
            }

            /* check value */
            if (this->meshing_inner_loop_maxiter == 0) {
                printf("ERROR: iteration parameter to switch \"meshing-innerloop-maxiter\" must be > 0.\n");
                return false;
            }
        }
        else if (s == "preserve-crease-edges") {
            this->meshing_preserve_crease_edges = true;
        }
        else if (s == "meshing-triangle-height") {
            try {
                this->meshing_cansurf_triangle_height_factor = std::stod(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"meshing-triangle-height\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: argument to switch \"meshing-triangle-height\" could not be converted to a double precision floating point value.\n");
                return false;
            }
            if (this->meshing_cansurf_triangle_height_factor <= 0) {
                printf("ERROR: invalid argument to switch \"meshing-triangle-height\". factor must be greater than 0.\n");
                return false;
            }
        }

        else if (s == "meshing-flush") {
            try {
                this->meshing_flush             = true;
                this->meshing_flush_face_limit  = stou(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: at least one argument to switch \"meshing-flush\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to switch \"meshing-flush\".\n");
                return false;
            }

            /* check value */
            if (this->meshing_flush_face_limit < 1024) {
                printf("ERROR: flush trigger face count parameter to switch \"meshing-flush\" must be in [1024, oo].\n");
                return false;
            }
        }
        else if (s == "no-meshing-flush") {
            this->meshing_flush = false;
        }
        else if (s == "meshing-merging-radiusfactor-decrement") {
            try {
                this->meshing_radius_factor_decrement = std::stod(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"meshing-merging-radiusfactor-decrement\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to switch \"meshing-merging-radiusfactor-decrement\".\n");
                return false;
            }

            /* check value */
            if (this->meshing_radius_factor_decrement <= 0.0 || this->meshing_radius_factor_decrement > 0.1) {
                printf("ERROR: radius factor decrement parameter to switch \"meshing-merging-initial-radiusfactor\" must be in )0.0, 0.1].\n");
                return false;
            }
        }
        else if (s == "meshing-complexedge-max-growthfactor") {
            try {
                this->meshing_complex_edge_max_growth_factor = std::stod(s_args[0]);
            }
            catch (std::out_of_range& ex) {
                printf("ERROR: argument to switch \"meshing-complexedge-max-growthfactor\" out of range.\n");
                return false;
            }
            catch (...) {
                printf("ERROR: invalid arguments to switch \"meshing-complexedge-max-growthfactor\".\n");
                return false;
            }
            
            /* check value */
            if (this->meshing_complex_edge_max_growth_factor <= 0.0) {
                printf("ERROR: complex edge growth factor parameter to switch \"meshing-complexedge-max-growthfactor\" must be > 0.0.\n");
                return false;
            }
        }
        else if (s == "debug-lvl")
        {
            try
            {
                enableComponentDebug(stou(s_args[0]));
                setMaxDebugLevel(stou(s_args[1]));
            }
            catch (std::out_of_range& ex)
            {
                printf("ERROR: One of the arguments to switch \"debug-lvl\" is out of range.\n");
                return false;
            }
            catch (...)
            {
                printf("ERROR: Invalid argument to switch \"debug-lvl\".\n");
                return false;
            }
        }
        else {
            printf("ERROR: unknown command line switch \"%s\".\n%s",
                s.c_str(), this->usage_text.c_str());
            return false;
        }
    }

    /* further checks on successfully parsed command line arguments */
    if (this->network_name == "") {
        printf("ERROR: no input file name given.\n");
        return false;
    }
    /* remove .swc suffix from network_name */
    else {
        size_t  network_name_last_dot_index = this->network_name.find_last_of(".");
        if (network_name_last_dot_index != std::string::npos) {
            std::string network_name_extension = this->network_name.substr(network_name_last_dot_index, std::string::npos);
            if (    network_name_extension != ".swc" &&
                    network_name_extension != ".amv" &&
                    network_name_extension != ".obj" &&
                    network_name_extension != ".CNG")
            {
                printf("ERROR: input file name invalid.\n");
                return false;
            }
            else {
                if (network_name_extension != ".CNG") {
                    this->network_name = this->network_name.substr(0, network_name_last_dot_index);
                }
            }
        }
    }

    return true;
}

bool
AnaMorph_cellgen::run()
{
    try {
        /* process command line arguments and return false if an error has occurred */
        if (!this->processCommandLineArguments()) {
            return false;
        }

        /* try to open input file */
        printf("AnaMorph cell generator (non-linear geometric modelling). swc input file name: \"%s.swc\"\n", this->network_name.c_str());

        /* analysis and mesh generation */
        if (this->ana) {
            printf("reading network from input swc file \"%s.swc\"..", this->network_name.c_str());fflush(stdout);
            NLM_CellNetwork<double> C(this->network_name);
            C.readFromNeuroMorphoSWCFile( this->network_name + ".swc", false);

            printf("done.\n"\
                "\t neuron vertices: %6zu   somas:              %6zu  axon vertices:   %6zu  dendrite vertices:   %6zu\n"\
                "\t neuron edges:    %6zu   neurite root edges: %6zu  axon root edges: %6zu  dendrite root edges: %6zu\n"\
                "\t axon segments:   %6zu   dendrite segments:  %6zu\n\n",
                    C.neuron_vertices.size(), C.soma_vertices.size(), C.axon_vertices.size(), C.dendrite_vertices.size(),
                    C.neuron_edges.size(), C.neurite_root_edges.size(), C.axon_root_edges.size(), C.dendrite_root_edges.size(), C.axon_segments.size(), C.dendrite_segments.size());

            /* retrieve, update and store settings inside network */
            NLM_CellNetwork<double>::Settings C_settings = C.getSettings();

            C_settings.analysis_nthreads                        = this->ana_nthreads;
            C_settings.analysis_univar_solver_eps               = this->ana_univar_solver_eps;
            C_settings.analysis_bivar_solver_eps                = this->ana_bivar_solver_eps;

            C_settings.partition_algo                           = this->partition_algo;
            C_settings.parametrization_algo                     = this->parametrization_algo;

            C_settings.meshing_flush                            = this->meshing_flush;
            C_settings.meshing_flush_face_limit                 = this->meshing_flush_face_limit;

            C_settings.meshing_n_soma_refs                      = this->meshing_n_soma_refs;
            C_settings.meshing_canal_segment_n_phi_segments     = this->meshing_canal_segment_n_phi_segments;
            C_settings.meshing_outer_loop_maxiter               = this->meshing_outer_loop_maxiter;
            C_settings.meshing_inner_loop_maxiter               = this->meshing_inner_loop_maxiter;

            C_settings.meshing_preserve_crease_edges            = this->meshing_preserve_crease_edges;
            C_settings.meshing_cansurf_triangle_height_factor	= this->meshing_cansurf_triangle_height_factor;

            C_settings.meshing_radius_factor_decrement          = this->meshing_radius_factor_decrement;
            C_settings.meshing_complex_edge_max_growth_factor   = this->meshing_complex_edge_max_growth_factor;

            C.updateSettings(C_settings);
            
            // This does not seem to be necessary and is really annoying when trying to
            // match original 1d positions to 3d positions generated with AnaMorph.
            /*
            // transform cell network to centroid system
            printf("transforming network coordinate system to soma 0 as origin.. ");fflush(stdout);
            C.transformToSomaSystem(C.soma_vertices.begin());
            printf("done.\n");
            */

            /* apply preconditioning algorithm */
            if (this->pc) {
                printf("applying cell network preconditioning. parameters:\n"\
                    "\t alpha = %5.4f\n\t beta = %5.4f\n\t gamma = %5.4f\n",
                    this->pc_alpha, this->pc_beta, this->pc_gamma);
                fflush(stdout);

                CellNetworkAlg::preliminaryPreconditioning(C, this->pc_alpha, this->pc_beta, this->pc_gamma);
                printf("done.\n\n");
            }

            // possibly scale radius (useful to create a cell-in-cell ER)
            if (scale_radius != 1.0)
            {
                std::cout << "scaling radii using factor " << scale_radius << "." << std::endl;
                CellNetworkAlg::scale_radii(C, scale_radius);
            }

            /* partition cell network and update geometry */
            printf("partitioning cell network.. ");fflush(stdout);
            C.partitionNetwork();
            printf("done.\n");

            printf("updating cell network geometry.. ");fflush(stdout);
            C.updateNetworkGeometry();
            printf("done.\n");

            /* perform full analysis */
            printf("performing single full geometric analysis iteration.. ");fflush(stdout);
            bool clean = C.performFullAnalysis();
            printf("done.\n");

            /* output morphology viewer file */
            printf("writing AnaMorph visualization file \"%s.amv\".. ", this->network_name.c_str());fflush(stdout);
            C.writeMorphViewFile( std::string(network_name) + ".amv");
            printf("done.\n");

            /* render cell network mesh */
            if (clean || this->force_meshing) {
                printf("rendering cell network to consistent mesh \"%s.obj\".\n", network_name.c_str());
                if (this->force_meshing) {
                    printf("\t NOTE: meshing forced in spite of potentially unclean network.\n");fflush(stdout);
                }

                C.renderCellNetwork<bool, bool, bool>(network_name);

                printf("done.\n\n");
            }
         
            if (this->meshing_individual_surfaces) {
                std::string ims_filename = std::string(network_name + "_individual_modelling_surfaces");

                printf("rendering cell network modelling surfaces individually to output mesh \"%s.obj\".\n", ims_filename.c_str());fflush(stdout);
                /* render geometric modelling surfaces individually and output mesh */
                C.renderModellingMeshesIndividually<bool, bool, bool>(std::string(network_name + "_individual_modelling_surfaces"));

                printf("done.\n");
            }
        }

        /* mesh-post-processing */
        if (this->pp_gec || this->pp_hc) {
            printf("post-processing union mesh \"%s.obj\".\n", this->network_name.c_str() );
            /* reload mesh to ram */
            Mesh<bool, bool, bool, double> M_cell;
            try {
                M_cell.readFromObjFile( (network_name + ".obj").c_str());
                if (this->pp_gec) {
                    printf("\t stage 1: improved edge-collapse algorithm. parameters:\n"\
                        "\t\t alpha:  %5.4f\n"\
                        "\t\t lambda: %5.4f\n"\
                        "\t\t mu:     %5.4f\n"\
                        "\t\t d:      %5d\n",
                        this->pp_gec_alpha, this->pp_gec_lambda, this->pp_gec_mu, this->pp_gec_d);

                    MeshAlg::greedyEdgeCollapsePostProcessing(
                        M_cell,
                        this->pp_gec_alpha,
                        this->pp_gec_lambda,
                        this->pp_gec_mu,
                        this->pp_gec_d);
                }

                if (this->pp_hc) {
                    printf("\t stage 2: HC Laplacian smoothing. parameters:\n"\
                        "\t\t alpha:   %5.4f\n"\
                        "\t\t beta:    %5.4f\n"\
                        "\t\t maxiter: %5d\n",
                        this->pp_hc_alpha, this->pp_hc_beta, this->pp_hc_maxiter);

                    MeshAlg::HCLaplacianSmoothing(
                        M_cell,
                        this->pp_hc_alpha,
                        this->pp_hc_beta,
                        this->pp_hc_maxiter);
                }

                M_cell.writeObjFile( (this->network_name + "_post_processed").c_str() );
            }
            catch (MeshEx& e) {
                if (e.error_type == MESH_IO_ERROR) {
                    printf("\t ERROR: could not open mesh obj file for post-processing. skipping..\n");
                }
                else throw;
            }
            printf("done.\n\n");
        }
        printf("all tasks performed. shutting down..\n");
        return EXIT_SUCCESS;
    }
    catch (char const *x) {
        printf("ERROR: Caught string exception at top level. message: \"%s\"\n", x);
        return  EXIT_FAILURE;
    }
    catch (std::string& x) {
        printf("ERROR: Caught string exception at top level. message: \"%s\"\n", x.c_str());
        return  EXIT_FAILURE;
    }
    catch (GraphEx& e) {
        printf("ERROR: Caught GraphEx exception at top level. message: \"%s\".\n", e.error_msg.c_str());
        return  EXIT_FAILURE;
    }
    catch (MeshEx& e) {
        printf("ERROR: Caught MeshEx exception at top level. message: \"%s\".\n", e.error_msg.c_str());
        return  EXIT_FAILURE;
    }
    catch (std::runtime_error& e) {
        printf("ERROR: Caught std::runtime_error exception at top level. message: \"%s\".\n", e.what());
        return  EXIT_FAILURE;
    }
    catch (...) {
        printf("ERROR: Caught unhandled exception at top level.\n");
        return EXIT_FAILURE;
    }
}
