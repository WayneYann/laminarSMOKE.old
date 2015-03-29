
class PolimiSootAnalyzer
{
public:
	
	PolimiSootAnalyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>* thermodynamicsMapXML);
	void Analysis(const double T, const double P_Pa, const double rhoGas, const Eigen::VectorXd &omegaGas, const Eigen::VectorXd &xGas);
	void Distribution();
	void WriteBinData();
void WriteDistribution(std::ofstream& fSootDistribution, const double index);
void WriteDistributionLabel(std::ofstream &fSoot);

	double fv_small() const { return fv_small_;}
	double rho_small() const { return rho_small_;}
	double N_small() const { return N_small_;}
	double omega_small() const { return omega_small_;}
	double x_small() const { return x_small_;}
	double h_over_c_small() const { return h_over_c_small_;}
	double o_over_c_small() const { return o_over_c_small_;}
	double o_over_h_small() const { return o_over_h_small_;}

	double fv_large() const { return fv_large_;}
	double rho_large() const { return rho_large_;}
	double N_large() const { return N_large_;}
	double omega_large() const { return omega_large_;}
	double x_large() const { return x_large_;}
	double h_over_c_large() const { return h_over_c_large_;}
	double o_over_c_large() const { return o_over_c_large_;}
	double o_over_h_large() const { return o_over_h_large_;}

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermo_;

	std::string minimum_bin_;
	unsigned int bin_index_zero_;
	unsigned int bin_index_final_;
	double bin_density_A_;
	double bin_density_B_;
	double Df_;

	unsigned int iC_;
	unsigned int iO_;
	unsigned int iH_;
	unsigned nspecies_;

	std::vector<double> bin_density_;
	std::vector<unsigned int> bin_indices_;
	std::vector<std::string> bin_names_;
	std::vector<double> bin_mw_;
	std::vector<double> bin_m_;
	std::vector<double> bin_ds_;
	std::vector<double> bin_V_;
	std::vector<double> bin_c_;
	std::vector<double> bin_h_;
	std::vector<double> bin_o_;
	std::vector<double> bin_h_over_c_;
	std::vector<double> bin_o_over_c_;
	std::vector<double> bin_o_over_h_;
	std::vector<double> bin_indices_large_;
	std::vector<double> bin_indices_small_;
	std::vector<double> bin_np_;
	std::vector<double> bin_dc_;
	std::vector<double> bin_d_;

	std::vector<double> bin_omega_;
	std::vector<double> bin_x_;
	std::vector<double> bin_fv_;
	std::vector<double> bin_rho_;
	std::vector<double> bin_N_;

	std::vector<double> bin_baskets_;

	std::vector< std::vector<unsigned int> >bin_baskets_indices_;



	std::vector<double> bin_baskets_d_;
	std::vector<double> bin_baskets_mw_;
	std::vector<double> bin_baskets_log10d_;
	std::vector<double> bin_baskets_dlog10d_;
	std::vector<double> dN_over_dlog10d_;

	std::vector<double> bin_baskets_V_;
	std::vector<double> bin_baskets_log10V_;
	std::vector<double> bin_baskets_dlog10V_;
	std::vector<double> dN_over_dlog10V_;

	std::vector<double> bin_baskets_m_;
	std::vector<double> bin_baskets_log10m_;
	std::vector<double> bin_baskets_dlog10m_;
	std::vector<double> dN_over_dlog10m_;

	std::vector<double> bin_baskets_N_;
	std::vector<double> bin_baskets_fv_;
	std::vector<double> bin_baskets_rho_;
	std::vector<double> bin_baskets_x_;
	std::vector<double> bin_baskets_omega_;

	double fv_small_;
	double rho_small_;
	double N_small_;
	double omega_small_;
	double x_small_;
	double h_over_c_small_;
	double o_over_c_small_;
	double o_over_h_small_;

	double fv_large_;
	double rho_large_;
	double N_large_;
	double omega_large_;
	double x_large_;
	double h_over_c_large_;
	double o_over_c_large_;
	double o_over_h_large_;

	double dmean_N_large_;
	double d32_N_large_;
	double dvariance_N_;
	double dstd_N_;
};

PolimiSootAnalyzer::PolimiSootAnalyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>* thermodynamicsMapXML) : 
		thermo_(*thermodynamicsMapXML)
	{

		minimum_bin_     = "BIN5";
		bin_index_zero_  = 10;
		bin_index_final_ = 20;
		bin_density_A_   = 1500.;
		bin_density_B_   = 1700.;
		Df_ = 1.8;

		iC_ = thermo_.IndexOfElement("C")-1;
		iH_ = thermo_.IndexOfElement("H")-1;
		iO_ = thermo_.IndexOfElement("O")-1;

		nspecies_ = thermo_.NumberOfSpecies();

		double min_mw_soot = 1.e32;
		for(int i=0;i<nspecies_;i++)
			if (thermo_.NamesOfSpecies()[i].compare(0, minimum_bin_.size(), minimum_bin_) == 0)
			{
				if (thermo_.MW()[i+1] < min_mw_soot)	min_mw_soot = thermo_.MW()[i+1];
			}
		

		unsigned int BIN12_index = 0;
		for(int i=0;i<nspecies_;i++)
			if (thermo_.NamesOfSpecies()[i].compare(0, 3, "BIN") == 0)
			{
				const double nc = thermo_.atomic_composition()(i,iC_);
				const double nh = thermo_.atomic_composition()(i,iH_);
				const double no = thermo_.atomic_composition()(i,iO_);

				// BIN density
				{
					
					

					const double index = std::log(nc/24.)/std::log(2.) + 1.;

					std::cout << thermo_.NamesOfSpecies()[i] << " " << nc << " " << nh << " " << no << std::endl;
					
					if (index<=bin_index_zero_) 
					{
						bin_density_.push_back(bin_density_A_);
					}
					else
					{
						const double m = (bin_density_B_-bin_density_A_)/(bin_index_final_-bin_index_zero_);
						const double c = bin_density_A_-m*bin_index_zero_;
						bin_density_.push_back(c+m*index);
					}
				}

				if (thermo_.NamesOfSpecies()[i].compare(0, 5, "BIN12") == 0 && BIN12_index == 0)
					BIN12_index = bin_indices_.size();

				bin_indices_.push_back(i);						// Index of bin in the gas phase kinetic scheme [-]
				bin_names_.push_back(thermo_.NamesOfSpecies()[i]);			// Name of bin in the gas phase kinetic scheme [-]
				bin_mw_.push_back(thermo_.MW()[i+1]);					// Molecular weight [kg/kmol]
				bin_m_.push_back(thermo_.MW()[i+1]/PhysicalConstants::Nav_kmol);		// Mass of particle [kg]
				
				bin_ds_.push_back( std::pow(6./PhysicalConstants::pi*thermo_.MW()[i+1]/(bin_density_[bin_density_.size()-1]/1000.)
							/(PhysicalConstants::Nav_mol), 1./3.)*1.e-2 );	// Diameter of particle [m]
				
				bin_V_.push_back( PhysicalConstants::pi/6.*std::pow(bin_ds_[bin_ds_.size()-1],3.) );	// Volume of particle [m3]
				bin_c_.push_back(nc);	// C
				bin_h_.push_back(nh);// H
				bin_o_.push_back(no);// O
				bin_h_over_c_.push_back(nh/nc);	// Ratio H/C
				bin_o_over_c_.push_back(no/nc);// Ratio O/C
				
				if (nh>0)	bin_o_over_h_.push_back(no/nh);	// Ratio O/H
				else    	bin_o_over_h_.push_back(0.);	// Ratio O/H
				
				if ( thermo_.MW()[i+1] >=  min_mw_soot)	bin_indices_large_.push_back(bin_indices_.size()-1);
				else                            	bin_indices_small_.push_back(bin_indices_.size()-1);

				// Collisional diameter and diameter
				{
					int iCollisional = false;
					const double index = std::log(nc / 24.) / std::log(2.) + 1;

					if (index > 12)
						iCollisional = true;
					
					if (iCollisional == true)
					{
						bin_np_.push_back(bin_m_[bin_m_.size()-1] / bin_m_[BIN12_index]);
						bin_dc_.push_back(std::sqrt(5. / 3.)*bin_ds_[BIN12_index] * pow(bin_np_[bin_np_.size()-1] / std::pow(1. + 2. / Df_, Df_ / 2.), 1. / Df_));
						bin_d_.push_back(bin_dc_[bin_dc_.size()-1]);
					}
					else
					{
						bin_np_.push_back(0.);
						bin_dc_.push_back(0.);
						bin_d_.push_back(bin_ds_[bin_ds_.size()-1]);
					}
				}
			}

			WriteBinData();



			{
				// Memory allocation
				bin_omega_.resize(bin_indices_.size());
				bin_x_.resize(bin_indices_.size());
				bin_fv_.resize(bin_indices_.size());
				bin_rho_.resize(bin_indices_.size());
				bin_N_.resize(bin_indices_.size());

				//if (iBin_ == true)
				{
					for(int i=0;i<bin_indices_.size();i++)
					{
						bool iNew = true;
						for(int k=0;k<bin_baskets_.size();k++)
						if (bin_c_[i] == bin_baskets_[k])
						{
							iNew = false;
							break;
						}
						if (iNew == true)
							bin_baskets_.push_back(bin_c_[i]);
					}

					std::sort(bin_baskets_.begin(), bin_baskets_.end());
					bin_baskets_indices_.resize(bin_baskets_.size());


					for(int i=0;i<bin_indices_.size();i++)
						for(int k=0;k<bin_baskets_.size();k++)
						{
							if (bin_c_[i] == bin_baskets_[k])
								bin_baskets_indices_[k].push_back(i);
						}

					bin_baskets_d_.resize(bin_baskets_.size());
					bin_baskets_mw_.resize(bin_baskets_.size());
					bin_baskets_log10d_.resize(bin_baskets_.size());
					bin_baskets_dlog10d_.resize(bin_baskets_.size());
					dN_over_dlog10d_.resize(bin_baskets_.size());

					bin_baskets_V_.resize(bin_baskets_.size());
					bin_baskets_log10V_.resize(bin_baskets_.size());
					bin_baskets_dlog10V_.resize(bin_baskets_.size());
					dN_over_dlog10V_.resize(bin_baskets_.size());

					bin_baskets_m_.resize(bin_baskets_.size());
					bin_baskets_log10m_.resize(bin_baskets_.size());
					bin_baskets_dlog10m_.resize(bin_baskets_.size());
					dN_over_dlog10m_.resize(bin_baskets_.size());

					bin_baskets_N_.resize(bin_baskets_.size());
					bin_baskets_fv_.resize(bin_baskets_.size());
					bin_baskets_rho_.resize(bin_baskets_.size());
					bin_baskets_x_.resize(bin_baskets_.size());
					bin_baskets_omega_.resize(bin_baskets_.size());
				}

				for(int k=0;k<bin_baskets_.size();k++)
				{
					bin_baskets_mw_[k] = 0.;
					bin_baskets_d_[k]  = 0.;
					bin_baskets_m_[k]  = 0.;
					bin_baskets_V_[k]  = 0.;

					for(int i=0;i<bin_baskets_indices_[k].size();i++)
					{	
						int j = bin_baskets_indices_[k][i];
						bin_baskets_d_[k] += bin_d_[j];
						bin_baskets_mw_[k] += bin_mw_[j];
						bin_baskets_V_[k] += bin_V_[j];
						bin_baskets_m_[k] += bin_m_[j];
					}
				
					bin_baskets_d_[k]  /= double(bin_baskets_indices_[k].size());
					bin_baskets_mw_[k] /= double(bin_baskets_indices_[k].size());
					bin_baskets_V_[k]  /= double(bin_baskets_indices_[k].size());
					bin_baskets_m_[k]  /= double(bin_baskets_indices_[k].size());
				
					bin_baskets_log10d_[k] = log10(bin_baskets_d_[k]);
					bin_baskets_log10V_[k] = log10(bin_baskets_V_[k]);
					bin_baskets_log10m_[k] = log10(bin_baskets_m_[k]);
				}

				// Intervals
				bin_baskets_dlog10d_[0] = ((bin_baskets_log10d_[1]+bin_baskets_log10d_[0])/2. - bin_baskets_log10d_[0])*2.;
				bin_baskets_dlog10V_[0] = ((bin_baskets_log10V_[1]+bin_baskets_log10V_[0])/2. - bin_baskets_log10V_[0])*2.;
				bin_baskets_dlog10m_[0] = ((bin_baskets_log10m_[1]+bin_baskets_log10m_[0])/2. - bin_baskets_log10m_[0])*2.;

				for(int k=1;k<bin_baskets_.size()-1;k++)
				{
					bin_baskets_dlog10d_[k] = (bin_baskets_log10d_[k+1]+bin_baskets_log10d_[k])/2. - (bin_baskets_log10d_[k]+bin_baskets_log10d_[k-1])/2.;
					bin_baskets_dlog10V_[k] = (bin_baskets_log10V_[k+1]+bin_baskets_log10V_[k])/2. - (bin_baskets_log10V_[k]+bin_baskets_log10V_[k-1])/2.;
					bin_baskets_dlog10m_[k] = (bin_baskets_log10m_[k+1]+bin_baskets_log10m_[k])/2. - (bin_baskets_log10m_[k]+bin_baskets_log10m_[k-1])/2.;
				}

				int k = bin_baskets_.size()-1;
				bin_baskets_dlog10d_[k] = ((bin_baskets_log10d_[k]+bin_baskets_log10d_[k-1])/2. - bin_baskets_log10d_[k-1])*2.;
				bin_baskets_dlog10V_[k] = ((bin_baskets_log10V_[k]+bin_baskets_log10V_[k-1])/2. - bin_baskets_log10V_[k-1])*2.;
				bin_baskets_dlog10m_[k] = ((bin_baskets_log10m_[k]+bin_baskets_log10m_[k-1])/2. - bin_baskets_log10m_[k-1])*2.;
		}
	}

void PolimiSootAnalyzer::Analysis(const double T, const double P_Pa, const double rhoGas, const Eigen::VectorXd &omegaGas, const Eigen::VectorXd &xGas)
{
	//if (iBin_ == true)
	{
		double small_eps = 1e-20;

		for(int i=0;i<bin_indices_.size();i++)
		{
			int j = bin_indices_[i];
			bin_omega_[i]	= omegaGas(j);				// mass fraction
			bin_x_[i]	= xGas(j);				// mole fraction
			bin_rho_[i]	= rhoGas*omegaGas(j);			// density [kg_soot/m3]
			bin_fv_[i]	= bin_rho_[i]/bin_density_[i];		// volume fraction [m3_soot/m3]
			bin_N_[i]	= bin_fv_[i]/bin_V_[i];			// [1/m3]
		}

		fv_small_	= 0.;
		rho_small_	= 0.;
		N_small_	= 0.;
		omega_small_	= 0.;
		x_small_	= 0.;
		h_over_c_small_ = 0.;
		o_over_c_small_ = 0.;
		o_over_h_small_ = 0.;

		for(int i=0;i<bin_indices_small_.size();i++)
		{
			int j = bin_indices_small_[i];
			fv_small_	+= bin_fv_[j];
			rho_small_	+= bin_rho_[j];
			N_small_	+= bin_N_[j];
			omega_small_	+= bin_omega_[j];
			x_small_	+= bin_x_[j];
			h_over_c_small_ += bin_omega_[j]*bin_h_over_c_[j];
			o_over_c_small_ += bin_omega_[j]*bin_o_over_c_[j];
			o_over_h_small_ += bin_omega_[j]*bin_o_over_h_[j];
		}

		double denominator = omega_small_;
		if ( denominator < small_eps)	denominator = 1.e32;
		h_over_c_small_ /= denominator;
		o_over_c_small_ /= denominator;
		o_over_h_small_ /= denominator;

		//if (iSoot_ == true)
		{
			fv_large_	= 0.;
			rho_large_	= 0.;
			N_large_	= 0.;
			omega_large_	= 0.;
			x_large_	= 0.;
			h_over_c_large_ = 0.;
			o_over_c_large_ = 0.;
			o_over_h_large_ = 0.;

			for(int i=0;i<bin_indices_large_.size();i++)
			{
				int j = bin_indices_large_[i];
				fv_large_	+= bin_fv_[j];
				rho_large_	+= bin_rho_[j];
				N_large_	+= bin_N_[j];
				omega_large_	+= bin_omega_[j];
				x_large_	+= bin_x_[j];
				h_over_c_large_ += bin_omega_[j]*bin_h_over_c_[j];
				o_over_c_large_ += bin_omega_[j]*bin_o_over_c_[j];
				o_over_h_large_ += bin_omega_[j]*bin_o_over_h_[j];
			}
			
			double denominator = omega_large_;
			if ( denominator < small_eps)	denominator = 1.e32;
			h_over_c_large_ /= denominator;
			o_over_c_large_ /= denominator;
			o_over_h_large_ /= denominator;
		}

	}
}

void PolimiSootAnalyzer::Distribution()
{
	//if (iBin_ == true)
	{
		for(int k=0;k<bin_baskets_.size();k++)
		{
			bin_baskets_N_[k]   = 0.;
			bin_baskets_fv_[k]  = 0.;
			bin_baskets_rho_[k] = 0.;
			bin_baskets_omega_[k] = 0.;
			bin_baskets_x_[k] = 0.;
			
			for(int i=0;i<bin_baskets_indices_[k].size();i++)
			{	
				int j = bin_baskets_indices_[k][i];
				bin_baskets_N_[k]	+= bin_N_[j];
				bin_baskets_fv_[k]	+= bin_fv_[j];
				bin_baskets_rho_[k]	+= bin_rho_[j];
				bin_baskets_omega_[k]	+= bin_omega_[j];
				bin_baskets_x_[k]	+= bin_x_[j];
			}

			dN_over_dlog10d_[k] = bin_baskets_N_[k] / bin_baskets_dlog10d_[k];
			dN_over_dlog10V_[k] = bin_baskets_N_[k] / bin_baskets_dlog10V_[k];
			dN_over_dlog10m_[k] = bin_baskets_N_[k] / bin_baskets_dlog10m_[k];
		}
	}

	//if (iBin_ == true)
	{
		double small_eps = 1e-16;

		double m0 = 0.;
		double m1 = 0.;
		double m2 = 0.;
		double m3 = 0.;
		for(int k=0;k<bin_baskets_.size();k++)
		{
			m0 += bin_baskets_N_[k];
			m1 += bin_baskets_N_[k]*bin_baskets_d_[k];
			m2 += bin_baskets_N_[k]*bin_baskets_d_[k]*bin_baskets_d_[k];
			m3 += bin_baskets_N_[k]*bin_baskets_d_[k]*bin_baskets_d_[k]*bin_baskets_d_[k];
		}

		if (m0 < small_eps)	m0 = 1e32;
		if (m2 < small_eps)	m2 = 1e32;

		// Mean diameters [m]
		dmean_N_large_	= m1/m0;
		d32_N_large_	= m3/m2;

		// Variance [m2]
		dvariance_N_ = 0.;
		for(int k=0;k<bin_baskets_.size();k++)
			dvariance_N_ += bin_baskets_N_[k]*std::pow(bin_baskets_d_[k]-dmean_N_large_,2.);
		dvariance_N_	/= m0;

		// Standard deviation [m]
		dstd_N_ = std::sqrt(dvariance_N_);
	}
}

void PolimiSootAnalyzer::WriteDistribution(std::ofstream& fSootDistribution, const double coordinate)
{
	if (fv_large_ > 1e-14)
	{
		for (int k=0;k<bin_baskets_.size();k++)
		{
			fSootDistribution << std::setw(20) << std::left << 1;
			fSootDistribution << std::setw(20) << std::left << coordinate;

			fSootDistribution << std::setw(20) << std::left << bin_baskets_d_[k];
			fSootDistribution << std::setw(20) << std::left << bin_baskets_m_[k];
			fSootDistribution << std::setw(20) << std::left << bin_baskets_V_[k];

			fSootDistribution << std::setw(20) << std::left << std::log10(bin_baskets_d_[k]);
			fSootDistribution << std::setw(20) << std::left << std::log10(bin_baskets_m_[k]);
			fSootDistribution << std::setw(20) << std::left << std::log10(bin_baskets_V_[k]);

			fSootDistribution << std::setw(20) << std::left << std::log10(bin_baskets_dlog10d_[k]);
			fSootDistribution << std::setw(20) << std::left << std::log10(bin_baskets_dlog10m_[k]);
			fSootDistribution << std::setw(20) << std::left << std::log10(bin_baskets_dlog10V_[k]);

			fSootDistribution << std::setw(20) << std::left << bin_baskets_fv_[k];
			fSootDistribution << std::setw(20) << std::left << bin_baskets_x_[k];
			fSootDistribution << std::setw(20) << std::left << bin_baskets_omega_[k];
			fSootDistribution << std::setw(20) << std::left << bin_baskets_rho_[k];
			fSootDistribution << std::setw(20) << std::left << bin_baskets_N_[k] / 1.e6; // [#/cm3]
			fSootDistribution << std::setw(20) << std::left << bin_baskets_N_[k] / bin_baskets_dlog10d_[k] / 1.e6; // [#/cm3]

			fSootDistribution << std::setw(20) << std::left << bin_baskets_fv_[k]/(fv_large_+fv_small_);
			fSootDistribution << std::setw(20) << std::left << bin_baskets_x_[k]/(x_large_+x_small_);
			fSootDistribution << std::setw(20) << std::left << bin_baskets_omega_[k]/(omega_large_+omega_small_);
			fSootDistribution << std::setw(20) << std::left << bin_baskets_rho_[k]/(rho_large_+rho_small_);
			fSootDistribution << std::setw(20) << std::left << bin_baskets_N_[k]/(N_large_+N_small_);

			fSootDistribution << std::endl;
		}
	}
}

void PolimiSootAnalyzer::WriteDistributionLabel(std::ofstream &fSoot)
{
    fSoot	<< std::setw(20) << std::left << "x[cm](1)"
			<< std::setw(20) << std::left << "bin(2)";

    fSoot	<< std::setw(20) << std::left << "d[m](3)"
			<< std::setw(20) << std::left << "m[?](4)"
			<< std::setw(20) << std::left << "V[?](5)"
			<< std::setw(20) << std::left << "log10d(6)"
			<< std::setw(20) << std::left << "log10m(7)"
			<< std::setw(20) << std::left << "log10V(8)"
			<< std::setw(20) << std::left << "dlog10d(9)"
			<< std::setw(20) << std::left << "dlog10m(10)"
			<< std::setw(20) << std::left << "dlog10V(11)";

    fSoot	<< std::setw(20) << std::left << "fv(12)"
			<< std::setw(20) << std::left << "x(13)"
			<< std::setw(20) << std::left << "y(14)"
			<< std::setw(20) << std::left << "rho[kg/m3](15)"
			<< std::setw(20) << std::left << "N[#/cm3](16)"
			<< std::setw(20) << std::left << "dN/dlog10[#/cm3](17)";

    fSoot	<< std::setw(20) << std::left << "fvN(18)"
			<< std::setw(20) << std::left << "xN(19)"
			<< std::setw(20) << std::left << "yN(20)"
			<< std::setw(20) << std::left << "rhoN[](21)"
			<< std::setw(20) << std::left << "NN[](22)";

	fSoot << std::endl << std::endl;
}

void PolimiSootAnalyzer::WriteBinData()
{
	std::ofstream fOutput;
	fOutput.open("BinSummary.out", std::ios::out);
	fOutput.setf(std::ios::scientific);

	fOutput << std::setw(5)  << std::left << "#"; 
	fOutput << std::setw(5)  << std::left << "#";
	fOutput << std::setw(10) << std::left << "Name";	
	fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << "MW[kg/kmol]";	// [kg/kmol]
	fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "d[nm]";			// [nm]
	fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "dsph[nm]";		// [nm]
	fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "dcol[nm]";		// [nm]
	fOutput << std::setw(16) << std::scientific << std::left << "m[mug]";						// [mug] 
	fOutput << std::setw(16) << std::scientific << std::left << "V[cm3]";						// [cm3]
	fOutput << std::setw(12) << std::scientific << std::left << "np[-]";						// [-]
	fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "C";
	fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "H"; 
	fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "O";
	fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "H/C"; 
	fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "O/C"; 
	fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "O/H";
	fOutput << std::endl;

	for(int i=0;i<bin_indices_small_.size();i++)
	{
		int j = bin_indices_small_[i];
		fOutput << std::setw(5)  << std::left << j+1; 
		fOutput << std::setw(5)  << std::left << i+1; 
		fOutput << std::setw(10) << std::left << bin_names_[j].c_str();	
		fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << bin_mw_[j];			// [kg/kmol]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_d_[j] * 1e9;		// [nm]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
		fOutput << std::setw(16) << std::scientific << std::left << bin_m_[j] * 1e9;					// [mug] 
		fOutput << std::setw(16) << std::scientific << std::left << bin_V_[j] * 1e6;					// [cm3]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(2) << bin_np_[j];			// [-]
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_c_[j];
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_h_[j]; 
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_o_[j];
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_h_over_c_[j]; 
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_c_[j]; 
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_h_[j];
		fOutput << std::endl;
	}
	fOutput << std::endl;
	
	//if (iSoot_ == true)
	{
		for(int i=0;i<bin_indices_large_.size();i++)
		{
			int j = bin_indices_large_[i];
			fOutput << std::setw(5)  << std::left << j+1; 
			fOutput << std::setw(5)  << std::left << i+1; 
			fOutput << std::setw(10) << std::left << bin_names_[j].c_str();	
			fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << bin_mw_[j];			// [kg/kmol]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_d_[j] * 1e9;		// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
			fOutput << std::setw(16) << std::scientific << std::left << bin_m_[j]*1e9;						// [mug] 
			fOutput << std::setw(16) << std::scientific << std::left << bin_V_[j]*1e6;						// [cm3]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(2) << bin_np_[j];			// [-]
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_c_[j]; 
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_h_[j]; 
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_o_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_h_over_c_[j]; 
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_c_[j]; 
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_h_[j];
			fOutput << std::endl;
		}
	}

	fOutput.close();
}
