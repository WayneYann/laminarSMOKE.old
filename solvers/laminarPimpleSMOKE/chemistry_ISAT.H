{
	//- Initial conditions
	scalarField& TCells = T.internalField();
	scalarField& QCells = Q.internalField();
	const scalarField& rhoCells = rho.internalField();
	const scalarField& vCells = mesh.V();

	if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+1;
		
		// Min and max values
		Eigen::VectorXd yMin(NEQ); for(unsigned int i=0;i<NEQ;i++) yMin(i) = 0.;  yMin(NC) = 200.;
		Eigen::VectorXd yMax(NEQ); for(unsigned int i=0;i<NEQ;i++) yMax(i) = 1.;  yMax(NC) = 6000.;
		Eigen::VectorXd y0(NEQ);
		Eigen::VectorXd yf(NEQ);

		Info <<" * Solving homogeneous chemistry (OpenSMOKE solver + ISAT)... "<<endl;
		{			
			unsigned int counter = 0;
			unsigned int nAddHOM = 0;
			unsigned int nGrowHOM = 0;
			unsigned int nRetHOM = 0;

			double cpuTimeRet      = 0.;
			double cpuTimeDI       = 0.;
			double cpuTimeGrowth   = 0.;
			double cpuTimeAddition = 0.;

			
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
				//- Solving for celli:	
				{
					for(unsigned int i=0;i<NC;i++)
						y0(i) = Y[i].internalField()[celli];
					y0(NC) = TCells[celli];
					
					if (constPressureBatchReactor == true)
					{
						// ISAT Algorithm
						{
							for(unsigned int i=0;i<NEQ;i++)
								phiISAT_HOM(i) = y0(i)*scalingFactors_ISAT(i);

							chemComp *phi0base = NULL;
							if(isat_HOM->retrieve(phiISAT_HOM, phi0base)) 
							{					
								double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
								
								// makes interpolation
								isat_HOM->interpol(phiISAT_HOM, RphiISAT_HOM, phi0base);
					
								//check negative value
								for(unsigned int i=0;i<NEQ;i++)
									yf(i) = std::max(RphiISAT_HOM(i), 0.)/scalingFactors_ISAT(i);

								nRetHOM++;	
								
								double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
								
								cpuTimeRet += (t2-t1);
							} 
							else 
							{			
								// Direct integration
								{
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
		
									// Set reactor
									batchReactorHomogeneousConstantPressure.SetReactor(thermodynamicPressure);
									batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
						
									// Set initial conditions
									odeSolverConstantPressure.SetInitialConditions(t0, y0);

									// Additional ODE solver options
									if (nGrowHOM == 0 && nAddHOM == 0)
									{
										// Set linear algebra options
										odeSolverConstantPressure.SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
										odeSolverConstantPressure.SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

										// Set relative and absolute tolerances
										odeSolverConstantPressure.SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
										odeSolverConstantPressure.SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

										// Set minimum and maximum values
										odeSolverConstantPressure.SetMinimumValues(yMin);
										odeSolverConstantPressure.SetMaximumValues(yMax);
									}

									OdeSMOKE::OdeStatus status = odeSolverConstantPressure.Solve(tf);
									odeSolverConstantPressure.Solution(yf);

									// Move the solution from DI to ISAT
									for(unsigned int i=0;i<NEQ;i++)
										RphiISAT_HOM(i) = std::max(yf(i), 0.)*scalingFactors_ISAT(i);

									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

									cpuTimeDI += (t2-t1);
								}

								// Growth
								if(isat_HOM->grow(phiISAT_HOM, RphiISAT_HOM, phi0base)) 
								{
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();

									phi0base->growEOA(phiISAT_HOM);
									nGrowHOM++;

									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

									cpuTimeGrowth += (t2-t1);
								} 
								// Addition
								else
								{
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
									
									// compute mapping gradient
									calcMappingGradient(	phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, scalingFactors_ISAT, 
												luSolver_ISAT, (tf-t0), numberSubSteps_ISAT, &odeSolverConstantPressure);
							
									// add a new leaf 
									bool flag = isat_HOM->add(phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, phi0base); 
						
									if(flag == false)	
										Info << "ISAT Error - Addition process failed..." << endl;
									
									nAddHOM++;
									
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
									
									cpuTimeAddition += (t2-t1);
								}
							}	
				
							// check if is necessary to balance
							isat_HOM->cleanAndBalance();
						}
					}
					else
					{
						// Set reactor pressure
						batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], thermodynamicPressure, rhoCells[celli]);
						batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);
						
						// Set initial conditions
						odeSolverConstantVolume.SetInitialConditions(t0, y0);

						// Additional ODE solver options
						if (celli == 0)
						{
							// Set linear algebra options
							odeSolverConstantVolume.SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
							odeSolverConstantVolume.SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

							// Set relative and absolute tolerances
							odeSolverConstantVolume.SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
							odeSolverConstantVolume.SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

							// Set minimum and maximum values
							odeSolverConstantVolume.SetMinimumValues(yMin);
							odeSolverConstantVolume.SetMaximumValues(yMax);
						}
						
						// Solve
						OdeSMOKE::OdeStatus status = odeSolverConstantVolume.Solve(tf);
						odeSolverConstantVolume.Solution(yf);
					}
				}
			
				// Check mass fractions
				normalizeMassFractions(yf, celli, massFractionsTol);

				// Assign mass fractions
				for(int i=0;i<NC;i++)
					Y[i].internalField()[celli] = yf(i);

				//- Allocating final values: temperature
				if (energyEquation == true)
					TCells[celli] = yf(NC);

				if (counter%(int(0.20*mesh.nCells())+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

				counter++;

				// Output
				if (runTime.outputTime())
				{
					if (constPressureBatchReactor == true)
					{
						QCells[celli] = batchReactorHomogeneousConstantPressure.QR();
						if (outputFormationRatesIndices.size() != 0)
						{
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
                                       	      		                                                   thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
						}
					}
					else
					{
						QCells[celli] = batchReactorHomogeneousConstantVolume.QR();
						if (outputFormationRatesIndices.size() != 0)
						{
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
						}
					}
				}

			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;

			if(isatCheck == true) 
			{
				Info << endl;
				Info << " ********* ISAT HOM stats **********" << endl;
				
				Info << "   Direct Integration : " << isat_HOM->nAdd()+isat_HOM->nGrow()  << " (" << nAddHOM+nGrowHOM << ")" << " (" << (nAddHOM+nGrowHOM)/double(mesh.nCells())*100. << "%)" << endl;
				Info << "      Add             : " << isat_HOM->nAdd()  << " (" << nAddHOM  << ")" << " (" << nAddHOM/double(mesh.nCells())*100. << "%)" << endl;
				Info << "      Grow            : " << isat_HOM->nGrow() << " (" << nGrowHOM << ")" << " (" << nGrowHOM/double(mesh.nCells())*100. << "%)" << endl;
				Info << "   Retrieve           : " << isat_HOM->nUse()  << " (" << nRetHOM  << ")" << " (" << nRetHOM/double(mesh.nCells())*100. << "%)" << endl;
				Info << endl;				
		
				const double cpuTimeIntegration = cpuTimeDI + cpuTimeGrowth + cpuTimeAddition;
				Info << "   CPU Integration  : " << cpuTimeIntegration  << " (" << cpuTimeIntegration/(tEnd-tStart)*100. << "%)" << endl;
				Info << "     CPU DI         : " << cpuTimeDI           << " (" << cpuTimeDI/(tEnd-tStart)*100.          << "%)" << endl;				
				Info << "     CPU Growth     : " << cpuTimeGrowth       << " (" << cpuTimeGrowth/(tEnd-tStart)*100.         << "%)" << endl;
				Info << "     CPU Addition   : " << cpuTimeAddition     << " (" << cpuTimeAddition/(tEnd-tStart)*100.    << "%)" << endl;
				Info << "   CPU Retrieve     : " << cpuTimeRet          << " (" << cpuTimeRet/(tEnd-tStart)*100.         << "%)" << endl;
				Info << endl;

				Info << "      BTS  : " << isat_HOM->nBTS()  << endl;
				Info << "      MRU  : " << isat_HOM->nMRU()  << endl;
				Info << "      MFU  : " << isat_HOM->nMFU()  << endl << endl;
				Info << endl;
			}
		}
	}
	else if (homogeneousReactions == true && odeParameterBatchReactorHomogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
		const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+1;
		
		// Min and max values
		OpenSMOKE::OpenSMOKEVectorDouble y0(NEQ);
		OpenSMOKE::OpenSMOKEVectorDouble yf(NEQ);

		Info <<" * Solving homogeneous chemistry (external ODE solver)... "<<endl;
		{			
			unsigned int counter = 0;
			
			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
				//- Solving for celli:	
				{
					for(unsigned int i=0;i<NC;i++)
						y0[i+1] = Y[i].internalField()[celli];
					y0[NC+1] = TCells[celli];
					
					if (constPressureBatchReactor == true)
					{
						// Set reactor pressure
						batchReactorHomogeneousConstantPressure.SetReactor(thermodynamicPressure);
						batchReactorHomogeneousConstantPressure.SetCheckMassFractions(true);
						batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
						
						if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
						{
							#if OPENSMOKE_USE_DVODE == 1
							typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DVODE> Solver;
							SolveOpenSourceSolvers<Solver>(dvodeSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);
							#endif
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
						{
							#if OPENSMOKE_USE_ODEPACK == 1
							typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODE> Solver;
							SolveOpenSourceSolvers<Solver>(dlsodeSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);					
							#endif
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
						{
							#if OPENSMOKE_USE_ODEPACK == 1
							typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODA> Solver;
							SolveOpenSourceSolvers<Solver>(dlsodaSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
							#endif							
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
						{
							#if OPENSMOKE_USE_DASPK == 1
							typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DASPK> Solver;
							SolveOpenSourceSolvers<Solver>(daspkSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);								
							#endif							
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
						{
							#if OPENSMOKE_USE_RADAU == 1
							typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_RADAU5> Solver;
							SolveOpenSourceSolvers<Solver>(radauSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);						
							#endif							
						}			
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
						{
							#if OPENSMOKE_USE_MEBDF == 1
							typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_MEBDF> Solver;
							SolveOpenSourceSolvers<Solver>(mebdfSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
							#endif							
						}	
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
						{
							#if OPENSMOKE_USE_SUNDIALS == 1
							typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE> Solver;
							SolveOpenSourceSolvers<Solver>(cvodeSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
							#endif							
						}																							
					}
					else
					{
						// Set reactor data
						batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], thermodynamicPressure, rhoCells[celli]);
						batchReactorHomogeneousConstantVolume.SetCheckMassFractions(true);
						batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);
						
						if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
						{
							#if OPENSMOKE_USE_DVODE == 1
							typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DVODE> Solver;
							SolveOpenSourceSolvers<Solver>(dvodeSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);
							#endif
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
						{
							#if OPENSMOKE_USE_ODEPACK == 1
							typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODE> Solver;
							SolveOpenSourceSolvers<Solver>(dlsodeSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);					
							#endif
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
						{
							#if OPENSMOKE_USE_ODEPACK == 1
							typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODA> Solver;
							SolveOpenSourceSolvers<Solver>(dlsodaSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
							#endif							
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
						{
							#if OPENSMOKE_USE_DASPK == 1
							typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DASPK> Solver;
							SolveOpenSourceSolvers<Solver>(daspkSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);								
							#endif							
						}
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
						{
							#if OPENSMOKE_USE_RADAU == 1
							typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_RADAU5> Solver;
							SolveOpenSourceSolvers<Solver>(radauSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);						
							#endif							
						}			
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
						{
							#if OPENSMOKE_USE_MEBDF == 1
							typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_MEBDF> Solver;
							SolveOpenSourceSolvers<Solver>(mebdfSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
							#endif							
						}	
						else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
						{
							#if OPENSMOKE_USE_SUNDIALS == 1
							typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CVODE> Solver;
							SolveOpenSourceSolvers<Solver>(cvodeSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
							#endif							
						}		
					}
				}
				
				// Check mass fractions
				normalizeMassFractions(yf, celli, massFractionsTol);

				// Assign mass fractions
				for(int i=0;i<NC;i++)
					Y[i].internalField()[celli] = yf[i+1];

				//- Allocating final values: temperature
				if (energyEquation == true)
					TCells[celli] = yf[NC+1];

				if (counter%(int(0.20*mesh.nCells())+1) == 0)
					Info <<"   Accomplished: " << counter << "/" << mesh.nCells() << endl;

				counter++;

				// Output
				if (runTime.outputTime())
				{
					if (constPressureBatchReactor == true)
					{
						QCells[celli] = batchReactorHomogeneousConstantPressure.QR();
						if (outputFormationRatesIndices.size() != 0)
						{
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantPressure.R()[outputFormationRatesIndices[i]+1] *
                                       	      		                                                   thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
						}
					}
					else
					{
						QCells[celli] = batchReactorHomogeneousConstantVolume.QR();
						if (outputFormationRatesIndices.size() != 0)
						{
							for (int i=0;i<outputFormationRatesIndices.size();i++)
								FormationRates[i].internalField()[celli] = batchReactorHomogeneousConstantVolume.R()[outputFormationRatesIndices[i]+1] *
                                           	                                                           thermodynamicsMapXML->MW()[outputFormationRatesIndices[i]+1];
						}
					}
				}

			}
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			
			Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per reactor)" << endl;
		}
	}
    
    Info<< " * T gas min/max (after chemistry) = " << min(T).value() << ", " << max(T).value() << endl;
}
