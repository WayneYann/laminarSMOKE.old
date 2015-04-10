
volScalarField omega_O
(
	IOobject
	(
		"omega_O",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(0, 0, 0, 0, 0)
);

volScalarField omega_H
(
	IOobject
	(
		"omega_H",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(0, 0, 0, 0, 0)
);

volScalarField omega_C
(
	IOobject
	(
		"omega_C",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionSet(0, 0, 0, 0, 0)
);

{
	unsigned int ns = Y.size();
	OpenSMOKE::OpenSMOKEVectorDouble massFractions(ns);
	OpenSMOKE::OpenSMOKEVectorDouble moleFractions(ns);

	const scalarField& TCells = T.internalField();
	const scalarField& pCells = p.internalField();
	      scalarField& rhoCells = rho.internalField();

	const unsigned int jC = thermodynamicsMapXML->IndexOfElement("C");
	const unsigned int jO = thermodynamicsMapXML->IndexOfElement("O");
	const unsigned int jH = thermodynamicsMapXML->IndexOfElement("H");

	forAll(TCells, celli)
	{
	}
 
	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		fvPatchScalarField& pomega_O  = omega_O.boundaryField()[patchi];
		fvPatchScalarField& pomega_H  = omega_H.boundaryField()[patchi];
		fvPatchScalarField& pomega_C  = omega_C.boundaryField()[patchi];

		forAll(pomega_H, facei)
		{
			pomega_H[facei] = 0.;
			pomega_O[facei] = 0.;
			pomega_C[facei] = 0.;
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryField()[patchi][facei];

			// Molecular weight
			double MWmix;
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions, MWmix, massFractions);

			if (jC>0)
			{
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					pomega_C[facei] += thermodynamicsMapXML->atomic_composition()(i,jC-1)*moleFractions[i+1]*12./MWmix;
			}
			if (jO>0)
			{
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					pomega_O[facei] += thermodynamicsMapXML->atomic_composition()(i,jO-1)*moleFractions[i+1]*16./MWmix;
			}
			if (jH>0)
			{
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					pomega_H[facei] += thermodynamicsMapXML->atomic_composition()(i,jH-1)*moleFractions[i+1]*1./MWmix;
			}	
		}
	}
}
