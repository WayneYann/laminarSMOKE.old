#include "soot/OpenSMOKE_PolimiSoot_Analyzer.h"

// Returns the integral mass of soot particles (large BINS) [kg] 
double largeBinMass(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_large_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_large_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral volume of soot particles (large BINs) [m3] 
double largeBinVolume(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_large_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_large_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/soot_analyzer->bin_density_large()[i]).value();
	}
	return sum;
}

// Returns the integral number of soot particles (large BINs) [#] 
double largeBinNumberParticles(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_large_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_large_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/soot_analyzer->bin_density_large()[i]/soot_analyzer->bin_V_large()[i]).value();
	}
	return sum;
}

// Returns the integral mass of soot particles (small BINS) [kg] 
double smallBinMass(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_small_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_small_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral volume of soot particles (small BINs) [m3] 
double smallBinVolume(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_small_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_small_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/soot_analyzer->bin_density_small()[i]).value();
	}
	return sum;
}

// Returns the integral number of soot particles (small BINs) [#] 
double smallBinNumberParticles(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_small_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_small_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]/soot_analyzer->bin_density_large()[i]/soot_analyzer->bin_V_small()[i]).value();
	}
	return sum;
}

// Returns the integral mass of dimer particles [kg] 
double sootDimerMass(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->soot_dimer_indices_global().size();i++)
	{
		const unsigned int j = soot_analyzer->soot_dimer_indices_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral mass of PAHs with 1 or 2 aromatic rings (340 nm) [kg] 
double pah_1_2_rings_Mass(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->pah_1_2_rings_indices_global().size();i++)
	{
		const unsigned int j = soot_analyzer->pah_1_2_rings_indices_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral mass of PAHs with 3 or 4 aromatic rings (400 nm) [kg] 
double pah_3_4_rings_Mass(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->pah_3_4_rings_indices_global().size();i++)
	{
		const unsigned int j = soot_analyzer->pah_3_4_rings_indices_global()[i];
		sum += fvc::domainIntegrate(rho*Y[j]).value();
	}
	return sum;
}

// Returns the integral mass of PAHs with more than 4 rings (500 nm) [kg] 
double pah_more_than_4_rings_Mass(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, volScalarField& rho, PtrList<volScalarField>& Y)
{
	return smallBinMass(soot_analyzer, rho, Y);
}

// Returns the total formation rate of soot particles (large BINs) [kg/m3/s]
double largeBinMassFormationRates(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, const OpenSMOKE::OpenSMOKEVectorDouble& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_large_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_large_global()[i]+1;
		sum += R_times_W[j];
	}
	return sum;
}

// Returns the total formation rate of soot particles (small BINs) [kg/m3/s]
double smallBinMassFormationRates(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, const OpenSMOKE::OpenSMOKEVectorDouble& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->bin_indices_small_global().size();i++)
	{
		const unsigned int j = soot_analyzer->bin_indices_small_global()[i]+1;
		sum += R_times_W[j];
	}
	return sum;
}

// Returns the total formation rate of PAHs with 1 or 2 aromatic rings (340 nm) [kg/m3/s] 
double pah12RingsFormationRates(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, const OpenSMOKE::OpenSMOKEVectorDouble& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->pah_1_2_rings_indices_global().size();i++)
	{
		const unsigned int j = soot_analyzer->pah_1_2_rings_indices_global()[i]+1;
		sum += R_times_W[j];
	}
	return sum;
}

// Returns the total formation rate of PAHs with 3 or 4 aromatic rings (400 nm) [kg/m3/s] 
double pah34RingsFormationRates(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, const OpenSMOKE::OpenSMOKEVectorDouble& R_times_W)
{
	double sum = 0.;
	for(int i=0;i<soot_analyzer->pah_3_4_rings_indices_global().size();i++)
	{
		const unsigned int j = soot_analyzer->pah_3_4_rings_indices_global()[i]+1;
		sum += R_times_W[j];
	}
	return sum;
}

// Returns the total formation rate of PAHs with more that aromatic rings (500 nm) [kg/m3/s] 
double pahMoreThan4RingsFormationRates(OpenSMOKE::PolimiSoot_Analyzer* soot_analyzer, const OpenSMOKE::OpenSMOKEVectorDouble& R_times_W)
{
	largeBinMassFormationRates(soot_analyzer, R_times_W);
}
