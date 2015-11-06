{
	if(zMixEquation == true)
	{
		volScalarField& Dmixi = Dmix[zMixDmixIndex];

		fvScalarMatrix zMixEqn
		(
		    	fvm::div(phi, zMix) 
		      - fvm::laplacian(rho*Dmixi, zMix) 
		        == 
		     	fvOptions(rho, zMix)
		);

		// Solve
		zMixEqn.relax();
		fvOptions.constrain(zMixEqn);
		zMixEqn.solve(mesh.solver("zMix"));
	}
}
