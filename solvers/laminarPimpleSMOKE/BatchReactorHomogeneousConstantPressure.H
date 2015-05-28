
#ifndef BatchReactorHomogeneousConstantPressure_H
#define	BatchReactorHomogeneousConstantPressure_H

class BatchReactorHomogeneousConstantPressure
{
public:

	BatchReactorHomogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
												OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap);

	void SetReactor( const double P0 );
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }
	double QR() const { return QR_; }

	unsigned int NumberOfEquations() const { return NE_; }

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	void SetDRG(OpenSMOKE::DRG* drg) { drg_ = drg; drgAnalysis_ = true; }
	void SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega );

private:

	double T_;
	double P0_;
	double rho_;

	double cTot_;
	double MW_;
	double CpMixMass_;
	double QR_;

	unsigned int NC_;
	unsigned int NE_;
	unsigned int NR_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN<double>&			kineticsMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	OpenSMOKE::OpenSMOKEVectorDouble r_;
	
	bool checkMassFractions_;
	bool energyEquation_;

	OpenSMOKE::DRG* drg_;
	bool drgAnalysis_;
};

BatchReactorHomogeneousConstantPressure::BatchReactorHomogeneousConstantPressure
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap ):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		NE_ = NC_+1;
		QR_ = 0.;
		
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		ChangeDimensions(NR_, &r_, true);
		
		checkMassFractions_ = false;
		energyEquation_ = true;
		drgAnalysis_ = false;
	}

void BatchReactorHomogeneousConstantPressure::SetReactor( const double P0 )
{
	P0_    = P0;
}

void BatchReactorHomogeneousConstantPressure::SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega )
{
	omega_ = omega;
}

int BatchReactorHomogeneousConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	if (drgAnalysis_ == false)
	{
		// Recover mass fractions
		if (checkMassFractions_ == true)
		{	for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = max(y[i], 0.);
		}
		else
		{
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = y[i];
		}

		// Recover temperature
		T_ = y[NC_+1];

		// Calculates the pressure and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
		cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
	    	rho_ = cTot_*MW_;
		Product(cTot_, x_, &c_);

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P0_);

		// Calculates kinetics
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P0_);
		kineticsMap_.KineticConstants();
		kineticsMap_.ReactionRates(c_);
		kineticsMap_.FormationRates(&R_);
		for (unsigned int i=1;i<=NC_;++i)	
			dy[i] = thermodynamicsMap_.MW()[i]*R_[i]/rho_;

	    	// Energy equation
	    	dy[NC_+1] = 0.;     
	    	if (energyEquation_ == true)
	    	{
			double CpMixMolar; 
			thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
			CpMixMass_ = CpMixMolar / MW_;
			QR_ = kineticsMap_.HeatRelease(R_);
			
			dy[NC_+1]  = QR_ / (rho_*CpMixMass_);
		}
	}

	else if (drgAnalysis_ == true)
	{	
		// Recover mass fractions
		if (checkMassFractions_ == true)
		{	
			for (unsigned int i=0;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omega_[j] = max(y[i+1], 0.);
			}	
		}
		else
		{
			for (unsigned int i=0;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omega_[j] = y[i+1];
			}
		}

		// Recover temperature
		unsigned int index_T = drg_->number_important_species()+1;
		T_ = y[index_T];

		// Calculates the pressure and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
		cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
	    	rho_ = cTot_*MW_;
		Product(cTot_, x_, &c_);

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P0_);

		// Calculates kinetics
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P0_);
		kineticsMap_.KineticConstants();
		kineticsMap_.ReactionRates(c_);

		// Reaction rates
		kineticsMap_.GetReactionRates(&r_);

		// Remove useless reactions
		for (unsigned int i=0;i<drg_->indices_unimportant_reactions().size();++i)
			r_[drg_->indices_unimportant_reactions()[i]+1] = 0.;

		// Formation rates
		kineticsMap_.stoichiometry().FormationRatesFromReactionRates(&R_, r_);

		// Equations
		for (unsigned int i=0;i<drg_->number_important_species();++i)	
		{
			const unsigned int j = drg_->indices_important_species()[i]+1;
			dy[i+1] = thermodynamicsMap_.MW()[j]*R_[j]/rho_;
		}	

	    	// Energy equation
	    	dy[index_T] = 0.;     
	    	if (energyEquation_ == true)
	    	{
			double CpMixMolar; 
			thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
			CpMixMass_ = CpMixMolar / MW_;
			QR_ = kineticsMap_.HeatRelease(R_);
			
			dy[index_T]  = QR_ / (rho_*CpMixMass_);
		}
	}
	
	return 0;
}

int BatchReactorHomogeneousConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
}

#endif // BatchReactorHomogeneousConstantPressure_H
