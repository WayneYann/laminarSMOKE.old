
volScalarField h
(
	IOobject
	(
		"h",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(0, 2, -2, 0, 0)
);

{
	unsigned int ns = Y.size();
	OpenSMOKE::OpenSMOKEVectorDouble massFractions(ns);
	OpenSMOKE::OpenSMOKEVectorDouble moleFractions(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	      scalarField& rhoCells = rho.internalField();

	forAll(TCells, celli)
	{
	}
 
	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		const fvPatchScalarField& pT  = T.boundaryField()[patchi];
		const fvPatchScalarField& pp  = p.boundaryField()[patchi];
		      fvPatchScalarField& ph  = h.boundaryField()[patchi];

		forAll(ph, facei)
		{
			thermodynamicsMapXML->SetPressure(pp[facei]);
			thermodynamicsMapXML->SetTemperature(pT[facei]);
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			double MWmix;
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions,MWmix,massFractions);

			thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(ph[facei], moleFractions);		// [J/kmol]
			ph[facei] /= MWmix;	
		}
	}
}
