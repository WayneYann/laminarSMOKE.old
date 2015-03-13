{
	if (implicitSourceTerm == true && jacobianUpdate == 1)
	{
		Info<< "Jacobian evaluation (continous)... " << endl;

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
					linear_model.reactionJacobianSparse( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J );

				for(int i=0;i<NC+1;i++)
					sourceImplicit[i].internalField()[celli] = J(i);
			
				for(int i=0;i<NC+1;i++)
					sourceExplicit[i].internalField()[celli] = Source[i+1] - J(i)*y[i+1];
			}
		}
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
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
						linear_model.reactionJacobianSparse( *thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], J );

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

				linear_model.reactionSourceTerms(	*thermodynamicsMapXML, *kineticsMapXML, y, pCells[celli], Source);

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
