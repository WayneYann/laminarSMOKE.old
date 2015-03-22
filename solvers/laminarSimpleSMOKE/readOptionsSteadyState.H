enum species_order_policy_enum {SPECIES_ORDER_POLICY_CONSTANT, SPECIES_ORDER_POLICY_SWEEP, SPECIES_ORDER_POLICY_RANDOM_SHUFFLE, SPECIES_ORDER_POLICY_ROTATE };

label jacobianUpdate = 1;
label propertiesUpdate = 1;
Switch implicitSourceTerm = true;
Switch sparseJacobian  = true;
species_order_policy_enum species_order_policy = SPECIES_ORDER_POLICY_CONSTANT;

const dictionary& steadyStateDictionary = solverOptions.subDict("SteadyState");
{
	sparseJacobian     = Switch(steadyStateDictionary.lookup(word("sparseJacobian")));
	jacobianUpdate     = readLabel(steadyStateDictionary.lookup("jacobianUpdate"));
	propertiesUpdate   = readLabel(steadyStateDictionary.lookup("propertiesUpdate"));
	implicitSourceTerm = Switch(steadyStateDictionary.lookup(word("implicitSourceTerm")));

	implicitSourceTerm = Switch(steadyStateDictionary.lookup(word("implicitSourceTerm")));

	word order_policy(steadyStateDictionary.lookup("orderSpecies"));
	if (order_policy == "constant")		species_order_policy = SPECIES_ORDER_POLICY_CONSTANT;
	else if (order_policy == "sweep")	species_order_policy = SPECIES_ORDER_POLICY_SWEEP;
	else if (order_policy == "random")	species_order_policy = SPECIES_ORDER_POLICY_RANDOM_SHUFFLE;
	else if (order_policy == "rotate")	species_order_policy = SPECIES_ORDER_POLICY_ROTATE;
	else
	{
		Info << "Wrong orderPolicy option: constant | sweep | random | rotate" << endl;
		abort();
	}
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

std::vector<int> species_order(thermodynamicsMapXML->NumberOfSpecies());
for(int i=0;i<species_order.size();i++)
{
	species_order[i] = i;
}

