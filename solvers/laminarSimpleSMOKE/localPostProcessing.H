{
	if (minMaxCounter == minMaxUpdate)
	{
		for (label j=0; j<Y.size(); j++)
    		{
			Info << " * Max " << Y[j].name() << " = " << max(Y[j]).value() << endl;		
		}
		
		minMaxCounter = 0;
		
	}

	minMaxCounter++;
}

if (iPolimiSoot == true)
{
	if (sootCounter == sootUpdate)
	{
		const double large_bin_mass = largeBinMass(sootAnalyzer, rho, Y);
		const double small_bin_mass = smallBinMass(sootAnalyzer, rho, Y);
		const double large_bin_volume = largeBinVolume(sootAnalyzer, rho, Y);
		const double small_bin_volume = smallBinVolume(sootAnalyzer, rho, Y);
		const double large_bin_particles = largeBinNumberParticles(sootAnalyzer, rho, Y);
		const double small_bin_particles = smallBinNumberParticles(sootAnalyzer, rho, Y);
		const double soot_dimer_mass = sootDimerMass(sootAnalyzer, rho, Y);
		const double pah_1_2_rings_mass = pah_1_2_rings_Mass(sootAnalyzer, rho, Y);
		const double pah_3_4_rings_mass = pah_3_4_rings_Mass(sootAnalyzer, rho, Y);
		const double pah_more_than_4_rings_mass = pah_more_than_4_rings_Mass(sootAnalyzer, rho, Y);
		
		*fSootIntegrals << std::setw(20) << large_bin_mass;
		*fSootIntegrals << std::setw(20) << large_bin_volume;
		*fSootIntegrals << std::setw(20) << large_bin_particles;		
		*fSootIntegrals << std::setw(20) << small_bin_mass;
		*fSootIntegrals << std::setw(20) << small_bin_volume;
		*fSootIntegrals << std::setw(20) << small_bin_particles;
		*fSootIntegrals << std::setw(20) << soot_dimer_mass;		
		*fSootIntegrals << std::setw(20) << pah_1_2_rings_mass;
		*fSootIntegrals << std::setw(20) << pah_3_4_rings_mass;
		*fSootIntegrals << std::setw(20) << pah_more_than_4_rings_mass;
		*fSootIntegrals << std::endl;

		Info << "Soot mass:            " << large_bin_mass << endl;
		Info << "Soot volume:          " << large_bin_volume << endl;
		Info << "Soot particles:       " << large_bin_particles << endl;		
		Info << "Precursors mass:      " << small_bin_mass << endl;
		Info << "Precursors volume:    " << small_bin_volume << endl;
		Info << "Precursors particles: " << small_bin_particles << endl;
		Info << "Dimer mass:           " << soot_dimer_mass << endl;		
		Info << "PAH(340) mass:        " << pah_1_2_rings_mass << endl;
		Info << "PAH(400) mass:        " << pah_3_4_rings_mass << endl;
		Info << "PAH(500) mass:        " << pah_more_than_4_rings_mass << endl;

		sootCounter = 0;
	}

	sootCounter++;
}
