label jacobianUpdate = 1;
label propertiesUpdate = 1;
Switch implicitSourceTerm = true;
Switch sparseJacobian  = true;

const dictionary& steadyStateDictionary = solverOptions.subDict("SteadyState");
{
	sparseJacobian     = Switch(steadyStateDictionary.lookup(word("sparseJacobian")));
	jacobianUpdate     = readLabel(steadyStateDictionary.lookup("jacobianUpdate"));
	propertiesUpdate   = readLabel(steadyStateDictionary.lookup("propertiesUpdate"));
	implicitSourceTerm = Switch(steadyStateDictionary.lookup(word("implicitSourceTerm")));
}

label jacobianCounter = jacobianUpdate;
label propertiesCounter = propertiesUpdate;

std::vector<Eigen::VectorXd> Jstore;
if (jacobianUpdate != 1)
{
	Jstore.resize(mesh.nCells());
	for(int i=0;i<mesh.nCells();i++)
		Jstore[i].resize(thermodynamicsMapXML->NumberOfSpecies()+1);
}
