if (homogeneousReactions == true && (speciesEquations == true || energyEquation == true))
{
	
	if (implicitSourceTerm == true && jacobianUpdate == 1)
	{
		Info<< "Jacobian evaluation (continous)... " << endl;
		
	//	std::ofstream* fOut[thermodynamicsMapXML->NumberOfSpecies()];
	//	for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	//	{
	//		fOut[i] = new std::ofstream(thermodynamicsMapXML->NamesOfSpecies()[i].c_str(), std::ios::out);
	//		fOut[i]->setf(std::ios::scientific);
	//	}

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
	
			Eigen::VectorXd J(NC+1);
			OpenSMOKE::OpenSMOKEVectorDouble Source(thermodynamicsMapXML->NumberOfSpecies()+1);	
			OpenSMOKE::OpenSMOKEVectorDouble y(thermodynamicsMapXML->NumberOfSpecies()+1);

			forAll(TCells, celli)
			{
				for(int i=0;i<NC;i++)
					y[i+1] = Y[i].internalField()[celli];
				y[NC+1] = TCells[celli];

				linear_model.reactionSourceTerms(	*thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source);


				if (sparseJacobian == false)
					linear_model.reactionJacobian( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J );
				else
					linear_model.reactionJacobianSparse( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J, energyEquation, speciesEquations);

				for(int i=0;i<NC+1;i++)
					sourceImplicit[i].internalField()[celli] = J(i);
			
				for(int i=0;i<NC+1;i++)
					sourceExplicit[i].internalField()[celli] = Source[i+1] - J(i)*y[i+1];

		//		if (TCells[celli] > 2200. && TCells[celli] < 2220.)
		//		{
		//			for(int i=0;i<NC;i++)
		//				*fOut[i] << Source[i+1] << " " << J(i) << " " << y[i+1] << " " << J(i)*y[i+1] << " " << Source[i+1] - J(i)*y[i+1] << std::endl;
		//		}
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

	//	for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	//		fOut[i]->close();
	//
	//	delete[] fOut;
		
	}
	else if (implicitSourceTerm == true && jacobianUpdate != 1)
	{
		Info<< "Jacobian evaluation (discrete)... " << endl;

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
	
			Eigen::VectorXd J(NC+1);
			OpenSMOKE::OpenSMOKEVectorDouble Source(thermodynamicsMapXML->NumberOfSpecies()+1);	
			OpenSMOKE::OpenSMOKEVectorDouble y(thermodynamicsMapXML->NumberOfSpecies()+1);

			forAll(TCells, celli)
			{
				for(int i=0;i<NC;i++)
					y[i+1] = Y[i].internalField()[celli];
				y[NC+1] = TCells[celli];

				linear_model.reactionSourceTerms(	*thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source);

				if (jacobianCounter == jacobianUpdate)
				{
					if (sparseJacobian == false)
						linear_model.reactionJacobian( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J );
					else
						linear_model.reactionJacobianSparse( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J, energyEquation, speciesEquations);

					Jstore[celli] = J;
				}

				for(int i=0;i<NC+1;i++)
					sourceImplicit[i].internalField()[celli] = Jstore[celli](i);
			
				for(int i=0;i<NC+1;i++)
					sourceExplicit[i].internalField()[celli] = Source[i+1] - Jstore[celli](i)*y[i+1];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;

		if (jacobianCounter == jacobianUpdate)
			jacobianCounter = 0;

		jacobianCounter++;
	}
	else
	{
		Info<< "Source terms evaluation... " << endl;

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();	
		{
			const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();

			const scalarField& TCells = T.internalField();
			const scalarField& pCells = p.internalField(); 
	
			OpenSMOKE::OpenSMOKEVectorDouble Source(thermodynamicsMapXML->NumberOfSpecies()+1);	
			OpenSMOKE::OpenSMOKEVectorDouble y(thermodynamicsMapXML->NumberOfSpecies()+1);

			forAll(TCells, celli)
			{
				for(int i=0;i<NC;i++)
					y[i+1] = Y[i].internalField()[celli];
				y[NC+1] = TCells[celli];

				linear_model.reactionSourceTerms( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source );

				for(int i=0;i<NC+1;i++)
					sourceImplicit[i].internalField()[celli] = 0.;
			
				for(int i=0;i<NC+1;i++)
					sourceExplicit[i].internalField()[celli] = Source[i+1];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;	
	}
}
