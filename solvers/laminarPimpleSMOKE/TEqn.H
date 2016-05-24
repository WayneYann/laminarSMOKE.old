{
    if(energyEquation == true)
    {
		radiation->correct();

		massDiffusionInEnergyEquation *= 0.;
		if (iMassDiffusionInEnergyEquation == true)
		{
			for (label i=0; i<Y.size(); i++)
			       massDiffusionInEnergyEquation -= CpSpecies[i]*( fvc::reconstruct(J[i]) & fvc::grad(T));
		}

		fvScalarMatrix TEqn
		(
			  cp*fvm::ddt(rho,T)
			+ cp*fvm::div(phi, T)		    
			  ==
	    		  fvm::laplacian(lambda,T)
	   	 	+ massDiffusionInEnergyEquation
			+ DpDt
			+ radiation->divq(T)
	    		+ cp*fvOptions(rho, T)
		);

		TEqn.relax();
		fvOptions.constrain(TEqn);
		TEqn.solve(mesh.solver("T"));
		fvOptions.correct(T);
    }

    Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

    // Spark
    spark.SetSpark(runTime.value(), T);
}
