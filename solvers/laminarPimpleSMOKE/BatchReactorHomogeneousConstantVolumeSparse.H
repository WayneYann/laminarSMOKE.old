
#ifndef BatchReactorHomogeneousConstantVolumeSparse_H
#define	BatchReactorHomogeneousConstantVolumeSparse_H

class BatchReactorHomogeneousConstantVolumeSparse
{
public:

	BatchReactorHomogeneousConstantVolumeSparse(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap);

	void SetReactor( const double V0, const double P0, const double rho0);
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; }

	Eigen::SparseMatrix<double>& Jomega() { return Jomega_; }

	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }
	double QR() const { return QR_; }

	unsigned int NumberOfEquations() { return NE_; }

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

	virtual void SparseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::SparseMatrix<double> &J);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

private:

	double T_;
	double P_;
	double V0_;
	double rho0_;

	double cTot_;
	double MW_;
	double CvMixMass_;
	double QR_;

	unsigned int NC_;
	unsigned int NE_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN<double>&		kineticsMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	OpenSMOKE::OpenSMOKEVectorDouble U_;
	OpenSMOKE::OpenSMOKEVectorDouble u_;

	Eigen::SparseMatrix<double> Jomega_;
	Eigen::VectorXd JT_row_;
	Eigen::VectorXd JT_col_;
	OpenSMOKE::OpenSMOKEVectorDouble dy0_;

	
	bool checkMassFractions_;
	bool energyEquation_;
};

BatchReactorHomogeneousConstantVolumeSparse::BatchReactorHomogeneousConstantVolumeSparse(	
	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NE_ = NC_+1;
		QR_ = 0.;
		
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		ChangeDimensions(NC_, &U_, true);
		ChangeDimensions(NC_, &u_, true);

		ChangeDimensions(NC_+1, &dy0_, true);
		JT_row_.resize(NC_);
		JT_col_.resize(NC_+1);

		checkMassFractions_ = false;
		energyEquation_ = true;
	}

void BatchReactorHomogeneousConstantVolumeSparse::SetReactor( const double V0, const double P0, const double rho0)
{
	P_    = P0;
	V0_   = V0;
	rho0_ = rho0;
}

int BatchReactorHomogeneousConstantVolumeSparse::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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
	for (unsigned int i = 1; i <= NC_; ++i)
		c_[i] = omega_[i] * rho0_ / thermodynamicsMap_.MW()[i];

	// Calculates the pressure
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	cTot_ = rho0_ / MW_;
	P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P_);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P_);

	kineticsMap_.ReactionRates(c_);
	kineticsMap_.FormationRates(&R_);


	// Recovering residuals
	for (unsigned int i=1;i<=NC_;++i)	
		dy[i] = thermodynamicsMap_.MW()[i]*R_[i]/rho0_;

	// Energy equation
	dy[NC_+1] = 0.;   
	if (energyEquation_ == true)
	{
		double CpMixMolar; 
		thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
		CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
		thermodynamicsMap_.uMolar_Species(U_);

		QR_ = kineticsMap_.HeatRelease(R_);
		dy[NC_ + 1] = -OpenSMOKE::Dot(R_, U_) / (rho0_*CvMixMass_);
	}																						
     
	return 0;
}

void BatchReactorHomogeneousConstantVolumeSparse::SparseAnalyticalJacobian(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, Eigen::SparseMatrix<double> &J)
{
	if (energyEquation_ == false)
		std::cout << "Sparse without energy: not yet implemented!!!" << std::endl;

	// Recover mass fractions
	for (unsigned int i = 1; i <= NC_; ++i)
		omega_[i] = y[i];

	// Recover temperature
	T_ = y[NC_ + 1];

	// Recover matrix Jomega
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	cTot_ = rho0_ / MW_;
	P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;
	kineticsMap_.jacobian_sparsity_pattern_map()->Jacobian(omega_, T_, P_, Jomega_);

	// Row
	{
		Equations(t, y, dy0_);
		ElementByElementDivision(U_, thermodynamicsMap_.MW(), &u_);
		OpenSMOKE::OpenSMOKEVectorDouble cv_species(NC_);
		thermodynamicsMap_.cpMolar_Species(cv_species);
		for (unsigned int i = 1; i <= NC_; ++i)
		{
			cv_species[i] -= PhysicalConstants::R_J_kmol;
			cv_species[i] /= thermodynamicsMap_.MW()[i];
		}

		for (int k = 0; k < Jomega_.outerSize(); ++k)
		{
			double sum = 0.;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Jomega_, k); it; ++it)
			{
				const unsigned int index = it.row() + 1;
				it.valueRef() *= thermodynamicsMap_.MW()[index] / rho0_;
				sum += it.value()*u_[index];
			}
			JT_row_(k) = -(sum + dy0_[NC_ + 1] * cv_species[k + 1]) / CvMixMass_;
		}
	}

	// Column
	{
		OpenSMOKE::OpenSMOKEVectorDouble y1 = y;
		OpenSMOKE::OpenSMOKEVectorDouble dy1(NC_+1);

		const double eps = 0.1;
		y1[NC_+1] = y[NC_+1] + eps;
		
		Equations(t, y1, dy1);

		for (unsigned int i = 1; i <= NC_+1; ++i)
		{
			JT_col_(i-1) = (dy1[i] - dy0_[i]) / eps;
		}
	}

	// Semi-analytical Jacobian
	for (unsigned int k = 0; k < NC_; ++k)
	{
		Eigen::SparseMatrix<double>::InnerIterator it_local(Jomega_, k);
		for (Eigen::SparseMatrix<double>::InnerIterator it(J, k); it; ++it)
		{
			if (it.row() == NC_)
				it.valueRef() = JT_row_(it.col());
			else
			{
				it.valueRef() = it_local.value();
				++it_local;
			}
		}
	}

	for (Eigen::SparseMatrix<double>::InnerIterator it(J, NC_); it; ++it)
	{
		it.valueRef() = JT_col_(it.row());
	}
}

int BatchReactorHomogeneousConstantVolumeSparse::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
}

#endif // BatchReactorHomogeneousConstantVolumeSparse_H
