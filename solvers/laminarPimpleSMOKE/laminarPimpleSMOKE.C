/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laminarSMOKE

Description
    Solver for reactive flows with detailed kinetic mechanisms.

\*---------------------------------------------------------------------------*/

// This is not a steady state simulation
#define STEADYSTATE 0

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
#include "soot/OpenSMOKE_PolimiSoot_Analyzer.h"

// Reactor utilities
#include "reactors/utilities/Utilities"

// ODE solvers
#include "math/multivalue-ode-solvers/MultiValueSolver"
#include "ode/ODE_Parameters.h"

// OpenFOAM
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "interpolation.H"
#include "radiationModel.H"

// Additional include files
#include "sparkModel.H"
#include "utilities.H"
#include "laminarSMOKEthermoClass.H"

// Homogeneous reactors
#include "BatchReactorHomogeneousConstantPressure.H"
#include "BatchReactorHomogeneousConstantPressure_ODE_Interface.H"
#include "BatchReactorHomogeneousConstantVolume.H"
#include "BatchReactorHomogeneousConstantVolume_ODE_Interface.H"

// ISAT
#if OPENSMOKE_USE_ISAT == 1
    #include "ISAT.h"
    #include "mappingGradient.h"
    #include "numericalJacobian4ISAT.H"
#endif

// Soot
#include "sootUtilities.H"

template<typename Solver, typename OdeBatch>
void SolveOpenSourceSolvers(OdeBatch& ode, const double t0, const double tf, const OpenSMOKE::OpenSMOKEVectorDouble& y0, OpenSMOKE::OpenSMOKEVectorDouble& yf, const OpenSMOKE::ODE_Parameters& parameters)
{
	Solver o(ode);
	o.SetDimensions(y0.Size());
	o.SetAbsoluteTolerance(parameters.absolute_tolerance());
	o.SetRelativeTolerance(parameters.relative_tolerance());
	o.SetAnalyticalJacobian(false);
	o.SetInitialValues(t0, y0.GetHandle());
	o.Solve(tf);
	o.Solution(yf.GetHandle());
}		

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "readGravitationalAcceleration.H"
	#include "createBasicFields.H"
	#include "readOptions.H"
	#include "createChemicalFields.H"
	#include "createFvOptions.H"
	#include "createRadiationModel.H"
	#include "memoryAllocation.H"
	#include "properties.H"
	#include "createAdditionalFields.H"
	#include "initContinuityErrs.H"
	#include "readTimeControls.H"
	#include "compressibleCourantNo.H"
	#include "setInitialDeltaT.H"
	pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
       		#include "readTimeControls.H"
	        #include "compressibleCourantNo.H"
		#include "setDeltaT.H"


		runTime++;
		Info<< "Time = " << runTime.timeName() << nl << endl;

		if (momentumEquations == false)
		{		
			if (strangAlgorithm == STRANG_MOMENTUM_TRANSPORT_REACTION)
			{
				#include "Policy_TransportReaction.H"
			}
			else if (strangAlgorithm == STRANG_MOMENTUM_REACTION_TRANSPORT)
			{
				#include "Policy_ReactionTransport.H"
			}
		}
		else
		{
			if (strangAlgorithm == STRANG_MOMENTUM_TRANSPORT_REACTION)
			{
				#include "Policy_MomentumTransportReaction.H"
			}
			else if (strangAlgorithm == STRANG_MOMENTUM_REACTION_TRANSPORT)
			{
				#include "Policy_MomentumReactionTransport.H"
			}		
		}

		// Passive scalars
            	#include "zMixEqn.H"
		#include "tauEqn.H"

		// Local post processing
		#include "localPostProcessing.H"
	
		runTime.write();

		
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
