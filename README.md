urdme
=====

Using Comsol Multiphysics 3.5
=====
1. Start Comsol and open a model file, e.g. urdme-1.2/examples/mincde/coli.mph.2. From Comsol, start Matlab:     File > Client/Server/MATLAB > Connect to MATLAB3. Initialize the Matlab environment.¥ Change MatlabÕs working directory to the folder for the URDME model you wish to simulate. At the Matlab command prompt type	>> cd urdme-1.2/examples/mincde/4. Export the model geometry from Comsol to Matlab.¥ Update the Model data: Solve > Update Model3
¥ Export the data:File > Export > FEM Structure as ÕfemÕ5. Simulate the model. At the Matlab command prompt type:         >> umod = urdme(fem,ÕmincdeÕ)6. Visualize the results. At the Matlab command prompt type:	>> postplot(umod.comsol,ÕTetdataÕ,ÕMinD mÕ)
Using Comsol Multiphysics 4.x
====

1. Start the Comsol interface to Matlab (ÓLiveLinkÓ), ./comsol server matlab in Unix-based systems.2. Change MatlabÕs working directory to the folder for the URDME model you wish to simulate. At the Matlab command prompt type	>> cd urdme-1.2/examples/mincde/3. Load the Comsol geometry into Matlab:     	>> fem = mphload(Õcoli.mphÕ)4. Simulate the model. At the Matlab command prompt type:    	 >> umod = urdme(fem,ÕmincdeÕ);5. Visualize the results. At the Matlab command prompt type:	>> umod.comsol.result.create(Õres1Õ,ÕPlotGroup3DÕ);	>> umod.comsol.result(Õres1Õ).set(ÕtÕ,Õ900Õ);	>> umod.comsol.result(Õres1Õ).feature.create(Õsurf1Õ, ÕSurfaceÕ);	>> umod.comsol.result(Õres1Õ).feature(Õsurf1Õ).set(ÕexprÕ, ÕMinD mÕ); >> mphplot(umod.comsol,Õres1Õ);6. Optionally, save the output to a mph-file for further observations in the Comsol GUI: >> mphsave(umod.comsol,Õcoli output.mphÕ)