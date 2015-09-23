Vessel Tools
------------

Introduction
------------
These tools were developed to carry out image processing and analysis on vascular network images created in a specific way - by staining of the luminal epithelium.  Since the stain does not label the interior of vessels, and does not always label the vessel wall uniformly or completely, special segmentation tools were needed.

The tools were almost all developed as C++ programs employing the ITK library, at least for reading and writing tiff files, and also for certain filters.  Each tool is built as an executable that can be invoked from the command line.  For ease of use in each case a GUI has been written using Qt.  The user specifies file names, input parameters etc. in the GUI, and then invokes the executable program that does the work by clicking a button.  The program's completion code indicates that the run was successful, or that an error was encountered.

Besides the basic set of programs needed to extract the network topology, there are many tools for a range of associated purposes.  In most cases the names are self-explanatory.  Our group makes extensive use of Amira and CMGUI, and for this reason there are several tools for operating on network files in these formats.

For assistance with building or using the tools please contact g.bogle@auckland.ac.nz

Building the tools
------------------
The "engine" program builds are set up using cmake.  For example, with the MS Visual Studio 2010 C++ compiler on Windows, the setup command is:

cmake -G "Visual Studio 10 Win64"

The CMakeLists.txt file associated with each program contains instructions to copy the resulting .exe file to vessel-tools\bin\exec.

The Qt GUI programs are built with Qt Creator, from the .pro files, and the .exe files are copied to vessel-tools\bin.

Processing an image
-------------------
Starting with an 8-bit greyscale image, raw.tif, the sequence of steps to extract the network topology is as follows:

- smoother
This tool creates a smoothed image that serves as the background for local thresholding.  The input parameter is the "radius" of the cube used for voxel averaging.
raw.tif -> sm.tif

- thresholder
Local thresholding is specified by two parameters.
raw.tif, sm.tif -> thresh.tif

- filler
The main filling tool operates on the binary thresholded image, with three parameters.
thresh.tif -> fill.tif

- holefiller
Remaining holes in the segmented image are removed.
fill.tif -> seg.tif

- connecter
Connected objects in the segmented image are identified and saved as separate images.  Only the largest will be of interest.
seg.tif -> conn.tif

- thinner
The segmented, connected image is skeletonised.
conn.tif -> skel.tif

- topology
Information in the segmented and skeleton images is combined to generate a description of network topology as a set of nodes and edges, together with vessel diameters.  The network description is conveyed by two CMGUI files.
conn.tif, skel.tif -> topo.exelem, topo.exnode
