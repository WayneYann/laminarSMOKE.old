label jacobianUpdate = 1;
Switch implicitSourceTerm = true;

const dictionary& steadyStateDictionary = solverOptions.subDict("SteadyState");
{
	jacobianUpdate  = readLabel(steadyStateDictionary.lookup("jacobianUpdate"));
	implicitSourceTerm = Switch(steadyStateDictionary.lookup(word("implicitSourceTerm")));
}

label jacobianCounter = jacobianUpdate;
std::vector<Eigen::VectorXd> Jstore;

if (jacobianUpdate != 1)
{
	Jstore.resize(mesh.nCells());
	for(int i=0;i<mesh.nCells();i++)
		Jstore[i].resize(thermodynamicsMapXML->NumberOfSpecies()+1);
}
