Execute RUN_PARAMS to run the demo. 

Select an energy with the index iEnergy. The diretory 'inDir' should point to the
parameterized .obj file. If the .obj file contains a cut UV map, set the second argument
of readObjCut to true.

For users running an OS other than Windows, the mex files have to be compiled.
meshIsometricenergyC.cpp, AKVFC.cpp  require the header only library 'eigen'

http://eigen.tuxfamily.org/index.php?title=Main_Page

and OpenMP.

For improved performance replace the '\' linear solve with a solver that supports symbolic
factorization like Pardiso in the deformMesh.m file.