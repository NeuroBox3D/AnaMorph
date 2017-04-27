
# AnaMorph
Framework for creating 3d neuronal morphologies from point/diameter descriptions.


## Installing ##
1. Clone this repository using git:

		git clone https://github.com/NeuroBox3D/AnaMorph.git

2. Compile using cmake and make:
   In the following, we suppose the AnaMorph project has been cloned to $AnaMorphHome.

		cd $AnaMorphHome
		mkdir build && cd build
		cmake ..
		make
   
   AnaMorph makes heavy use of the C++11 standard. Up-to-date compilers are recommended.
   Once built successfully, the binary *am_cellgen* is present in \$AnaMorphHome/bin. You may want to add this to your \$PATH.

## Basic Usage ##
The most basic way to use AnaMorph is this: Suppose we have a valid SWC file *someCell.swc* in the current directory.
	
	am_cellgen -i someCell.swc

If all is well, then AnaMorph will output the final surface mesh as *someCell_post-processed.obj*.
Often, however, AnaMorph's consistency analysis will find some intersections in the morphology. To force creation of the (potentially self-intersecting) surface, one can use

	am_cellgen -i someCell.swc -force-meshing

## Advanced usage ##
AnaMorph has several parameters adjustable on the command line. Some of them are used fairly frequently,  e.g. *-meshing-cansurf-angularsegments* for the **surface resolution**:

	am_cellgen -i someCell.swc -force-meshing -meshing-cansurf-angularsegments 6

will produce neurites with a resolution of 6 nodes per cross-section (the default setting).

**Preconditioning** can be fine-tuned using the *-cellnet-pc* directive:

	am_cellgen -i someCell.swc -force-meshing -cellnet-pc 3.52 1.76 9999999

The first two parameters will cause short edges to be collapsed. Smaller parameters allow for shorter edges, bigger parameters increase the smoothing effect of preconditioning.

The **parametrization strategy** chosen for the spline interpolation is set using the *-cellnet-parametrization-strategy* argument:

	am_cellgen -i someCell.swc -force-meshing -cellnet-parametrization-strategy "centripetal"

The **post-processing** step of smoothing can be deactivated using:

	am_cellgen -i someCell.swc -force-meshing -no-mesh-pp-hc

A comprehensive list of adjustable parameters is contained in the file $AnaMorphHome/doc/ReadMe. It is also displayed if you simply call

	am_cellgen

without any arguments.
