### The following Matlab scripts perform the mentioned tasks
1. `PIVlab_commandline.m` performs PIV analysis to get the velocity field 
2. `CalculatePlotFlow_Shear.m` extracts velocity data and performs various operations using [`pivmat`](http://www.fast.u-psud.fr/pivmat/versions/pivmat4.10.zip) library
3. `MovingWindow.m ` is used to extract velocity field within the moving placode domain
4. `AnalyzeNematic.m` is used to analyze the anisotropy data obtained from _OrientationJ_ analysis
5. `RemoveAmnio_UbiHis_Byn_RemoveEdge.m` is used to separate velocity field in the _byn_ region to get mean sliding velocity of the placode
6. `Index.m` is used to find winding number of the nematic field of the image anisotropy obtained from _OrientationJ_ analysis
7. `PlotDefects.m` is used to plot the nematic field and topological defects
