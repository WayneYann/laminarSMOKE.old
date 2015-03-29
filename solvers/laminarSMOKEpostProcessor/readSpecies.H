PtrList<volScalarField> Y;
Y.resize(thermodynamicsMapXML->NumberOfSpecies());

PtrList<volScalarField> X;
if (iMoleFractions == true)
	X.resize(thermodynamicsMapXML->NumberOfSpecies());

PtrList<volScalarField> C;
if (iConcentrations == true)
	C.resize(thermodynamicsMapXML->NumberOfSpecies());

{
	Info<< "Reading species..." << endl;

	for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	{
		IOobject header
		(
			thermodynamicsMapXML->NamesOfSpecies()[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (header.headerOk())
		{
			std::cout << " * " << thermodynamicsMapXML->NamesOfSpecies()[i] << " -> from file..." << std::endl;
			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::NO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			volScalarField Ydefault
			(
				IOobject
				(
					"Ydefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					Ydefault
				)
			);
		}

		if (iConcentrations == true)
		{
			X.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"x::" + thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					Y[i]
				)
			);
		}

		if (iConcentrations == true)
		{
			C.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"c::" + thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, -3, 0, 0, 1)
				)
			);
		}
	}
}
