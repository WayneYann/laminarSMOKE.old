
Info << endl;
Info<< " * Reading chemical species mole fractions or mass fractions..." << endl;

PtrList<volScalarField> Y(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> X(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> Dmix(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> CpSpecies(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> FormationRates(outputFormationRatesIndices.size());

#if STEADYSTATE == 1
PtrList<volScalarField> sourceImplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
PtrList<volScalarField> sourceExplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
#endif

PtrList<volScalarField> Dsoret(transportMapXML->iThermalDiffusionRatios().size());
if (soretEffect == true)
{
	for(unsigned int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
	{ 
		unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
		indexSoret[index-1] = i;
		soretEffectList[index-1] = true;

		std::cout << "Soret effect: species " 	<< i+1 << "/" << transportMapXML->iThermalDiffusionRatios().size() 
							<< " " << thermodynamicsMapXML->NamesOfSpecies()[index-1] << " (" << index << ")" << std::endl;

		Dsoret.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"gas::Dsoret_"+ thermodynamicsMapXML->NamesOfSpecies()[index-1],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("Dsoret", dimensionSet(0, 2, -1, 0, 0), 1e-10)
			)
		);
	}
}

// Loop over all the species in the kinetic mechanism
for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{
	Info << "   Reading species: " << thermodynamicsMapXML->NamesOfSpecies()[i] << " > ";

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// Sets of diffusivities and specific heats of species
	{
		Dmix.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"gas::Dmix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("Dmix", dimensionSet(0, 2, -1, 0, 0), 1e-10)
			)
		);

		CpSpecies.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "thermo::Cp_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                                ),
                                mesh,
                                dimensionSet(0, 2, -2, -1, 0)
                        )
                );

		X.set
                (
                        i,
                        new volScalarField
                        (
                                IOobject
                                (
                                "thermo::X_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
                                mesh.time().timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                                ),
                                mesh,
                                dimensionSet(0, 0, 0, 0, 0)
                        )
                );
	}

	// check if field exists and can be read
	if (header.headerOk())
	{
		Info << thermodynamicsMapXML->NamesOfSpecies()[i] << endl;
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
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);
	}
	else
	{
		Info << "Ydefault" << endl;

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
					IOobject::AUTO_WRITE
				),
				Ydefault
			)
		);
	}
}

#if STEADYSTATE == 1

for (int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{
	sourceImplicit.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"SI::" + thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("SI", dimensionSet(1, -3, -1, 0, 0), 0.0)
		)
	);

	sourceExplicit.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"SE::" + thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionedScalar("SE", dimensionSet(1, -3, -1, 0, 0), 0.0)
		)
	);
}

sourceImplicit.set
(
	thermodynamicsMapXML->NumberOfSpecies(),
	new volScalarField
	(
		IOobject
		(
			"SI::T",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("SIT", dimensionSet(1, -1, -3, -1, 0), 0.0)
	)
);

sourceExplicit.set
(
	thermodynamicsMapXML->NumberOfSpecies(),
	new volScalarField
	(
		IOobject
		(
			"SE::T",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("SET", dimensionSet(1, -1, -3, 0, 0), 0.0)
	)
);

#endif

if (outputFormationRatesIndices.size() != 0)
{
	for (int i=0;i<outputFormationRatesIndices.size();i++)
	{
		FormationRates.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"R::" + thermodynamicsMapXML->NamesOfSpecies()[outputFormationRatesIndices(i)],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("R", dimensionSet(1, -3, -1, 0, 0), 0.0)
			)
		);
	}
}

PtrList<volScalarField> sootFields(18);
if (sootOnTheFlyFill == true)
{
	sootFields.set
	(
		0,
		new volScalarField
		(
			IOobject
			(
				"soot::fv_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		1,
		new volScalarField
		(
			IOobject
			(
				"soot::fv_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		2,
		new volScalarField
		(
			IOobject
			(
				"soot::rho_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		3,
		new volScalarField
		(
			IOobject
			(
				"soot::rho_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		4,
		new volScalarField
		(
			IOobject
			(
				"soot::N_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		5,
		new volScalarField
		(
			IOobject
			(
				"soot::N_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, -3, 0, 0, 0)
		)
	);

	sootFields.set
	(
		6,
		new volScalarField
		(
			IOobject
			(
				"soot::x_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		7,
		new volScalarField
		(
			IOobject
			(
				"soot::x_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		8,
		new volScalarField
		(
			IOobject
			(
				"soot::omega_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		9,
		new volScalarField
		(
			IOobject
			(
				"soot::omega_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		10,
		new volScalarField
		(
			IOobject
			(
				"pah::omega_1_2",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		11,
		new volScalarField
		(
			IOobject
			(
				"pah::omega_3_4",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

	sootFields.set
	(
		12,
		new volScalarField
		(
			IOobject
			(
				"pah::omega_more_4",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(0, 0, 0, 0, 0)
		)
	);

        sootFields.set
	(
		13,
		new volScalarField
		(
			IOobject
			(
				"soot::R_large",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);

        sootFields.set
	(
		14,
		new volScalarField
		(
			IOobject
			(
				"soot::R_small",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);

        sootFields.set
	(
		15,
		new volScalarField
		(
			IOobject
			(
				"pah::R_1_2",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);

        sootFields.set
	(
		16,
		new volScalarField
		(
			IOobject
			(
				"pah::R_3_4",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);

        sootFields.set
	(
		17,
		new volScalarField
		(
			IOobject
			(
				"pah::R_more_4",
				mesh.time().timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionSet(1, -3, -1, 0, 0)
		)
	);
}

#include "createFluxes.H"

// Create laminarSMOKEthermoClass object
std::vector<double> W(thermodynamicsMapXML->NumberOfSpecies());
std::vector<std::string> species_names(thermodynamicsMapXML->NumberOfSpecies());

for (unsigned int i=0;i<species_names.size();i++)
{
	species_names[i] = thermodynamicsMapXML->NamesOfSpecies()[i];
	W[i] = thermodynamicsMapXML->MW()[i+1];
}

autoPtr<Foam::laminarSMOKEthermoClass> myLaminarSMOKEthermo(Foam::laminarSMOKEthermoClass::New(T,p,Y, W, species_names));

