{	
	// Mass diffusion
	forAll (Y,i)
	{
		volScalarField& Xi = X[i];
		volScalarField& Dmixi = Dmix[i];
		dimensionedScalar MWi("MWi", dimensionSet(1,0,0,0,-1,0,0),scalar(thermodynamicsMapXML->MW()[i+1]) ); 

		J[i] = ( fvc::interpolate ( -rho*Dmixi*MWi/MWmix ) ) * ( fvc::interpolate ( fvc::grad(Xi) ) & mesh.Sf() );
	}

	// Soret effect
	forAll (Y,i)
	{
		if (soretEffect == true)
		{
			if (soretEffectList[i] == true)
			{
				J[i] += ( fvc::interpolate ( -rho*Dsoret[indexSoret[i]]/T ) ) *
					( fvc::interpolate (  fvc::grad(T) ) & mesh.Sf () ) ;
			}
		}
	}

	// Thermophoretic effect
	forAll (Y,i)
	{
		if (thermophoreticEffect == true)
		{
			if (thermophoreticEffectList[i] == true)
			{
				volScalarField& Yi = Y[i];

				J[i] += ( fvc::interpolate ( -0.55*mu/T*Yi ) ) *
					( fvc::interpolate ( fvc::grad(T) ) & mesh.Sf () ) ;
			}
		}
	}

	// Correction flux
	Jc *= 0.;
	forAll (Y,i)
	{
		Jc -= J[i];
	}

	// Corrected fluxes
	forAll (Y,i)
	{
		J[i] += Jc*fvc::interpolate(Y[i]);
	}
}
