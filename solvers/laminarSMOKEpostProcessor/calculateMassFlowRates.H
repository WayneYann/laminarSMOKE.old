{
	std::vector<double> mass_flow_rates;
	std::vector<std::string> boundary_names;
 
	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		scalar mass_flow_rate = sum(phi.boundaryField()[patchi]); 
		mass_flow_rates.push_back(mass_flow_rate);
		boundary_names.push_back(patches[patchi].name());
	}

	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "                                 Mass flow rates                                          " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	for(int i=0;i<mass_flow_rates.size();i++)
	{
		std::cout << std::setw(25) << std::left << boundary_names[i];
		std::cout << std::setw(18) << std::right << mass_flow_rates[i];
		std::cout << std::endl;
	}
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(25) << std::left << "Net";
	std::cout << std::setw(18) << std::right << std::accumulate(mass_flow_rates.begin(), mass_flow_rates.end(), 0.);
	std::cout << std::endl;
}

volScalarField mu
(
	IOobject
	(
		"mu",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(1, -1, -1, 0, 0)
);

volVectorField soot_v_thermo
(
	IOobject
	(
		"soot::v_thermo",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(0, 1, -1, 0, 0)
);

{
	unsigned int ns = Y.size();
	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);
	Eigen::VectorXd y_eigen(ns);
	Eigen::VectorXd x_eigen(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	      scalarField& rhoCells = rho.internalField();
	      scalarField& muCells = mu.internalField();

	forAll(TCells, celli)
	{	
		double mw;

		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			y[i+1] = Y[i].internalField()[celli];
		const double sum = y.SumElements();
		for(unsigned int i=0;i<ns;i++)
			y[i+1] /= sum;
		

		// Molecular weight
		thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

		// Eigen vectors
		for(unsigned int i=0;i<ns;i++)
		{
			y_eigen(i) = y[i+1];
			x_eigen(i) = x[i+1];
		}

		// Density
		rhoCells[celli] = pCells[celli]*mw/PhysicalConstants::R_J_kmol/TCells[celli];

		// Viscosity
		transportMapXML->SetPressure(pCells[celli]);
		transportMapXML->SetTemperature(TCells[celli]);
		transportMapXML->DynamicViscosity(muCells[celli],x);

	}


	forAll(T.boundaryField(), patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
		const fvPatchScalarField& pp  = p.boundaryField()[patchi];
		      fvPatchScalarField& prho  = rho.boundaryField()[patchi];
		      fvPatchScalarField& pmu   = mu.boundaryField()[patchi];


		forAll(pT, facei)
		{
			double mw;

			// Extract the mass fractions
			for(unsigned int i=0;i<ns;i++)
				y[i+1] = Y[i].boundaryField()[patchi][facei];
			const double sum = y.SumElements();
			for(unsigned int i=0;i<ns;i++)
				y[i+1] /= sum;

			// Molecular weight
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

			// Eigen vectors
			for(unsigned int i=0;i<ns;i++)
			{
				y_eigen(i) = y[i+1];
				x_eigen(i) = x[i+1];
			}

			// Density
			prho[facei] = pp[facei]*mw/PhysicalConstants::R_J_kmol/pT[facei];

			// Viscosity
			transportMapXML->SetPressure(pp[facei]);
			transportMapXML->SetTemperature(pT[facei]);
			transportMapXML->DynamicViscosity(pmu[facei],x);
		}
	}

	Info << "Thermophoretic velocity" << endl;
	soot_v_thermo = -0.55*mu/rho/T*fvc::grad(T);
	soot_v_thermo.write();
}
