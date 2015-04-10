{
	std::vector<double> mass_flow_rates;
	std::vector<std::string> boundary_names;
 
	const fvPatchList& patches = mesh.boundary();
        forAll(patches, patchi)
	{
		scalar mass_flow_rate = sum(phi.boundaryField()[patchi]); 
		mass_flow_rates.push_back(mass_flow_rate);
		boundary_names.push_back(patches[patchi].name());
	}

	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "                                 Mass flow rates                                          " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	for(int i=0;i<mass_flow_rates.size();i++)
	{
		std::cout << std::setw(25) << std::left << boundary_names[i];
		std::cout << std::setw(18) << std::right << mass_flow_rates[i];
		std::cout << std::endl;
	}
	std::cout << "------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::setw(25) << std::left << "Net";
	std::cout << std::setw(18) << std::right << std::accumulate(mass_flow_rates.begin(), mass_flow_rates.end(), 0.);
	std::cout << std::endl;
}
