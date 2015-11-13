#ifndef laminarSMOKEthermoClass_H
#define laminarSMOKEthermoClass_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "fvMatricesFwd.H"
#include "Switch.H"
#include <vector>
#include <string>

namespace Foam
{

class laminarSMOKEthermoClass : public IOdictionary
{
public:

	//- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

	// Constructors
        laminarSMOKEthermoClass(const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names);

        //- Return a reference to the selected radiation model
        static autoPtr<laminarSMOKEthermoClass> New(const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names);

	virtual ~laminarSMOKEthermoClass();

public:

	const volScalarField& T() const { return T_; }
	const volScalarField& p() const { return p_; }
	const PtrList<volScalarField>& Y() const { return Y_; }
	const volScalarField& Y(const label i) const { return Y_[i]; }
	double W(const label i) const { return W_[i]; }
	unsigned int species_index(const std::string name) const;
	
private:

	const volScalarField& T_;
	const volScalarField& p_;
	const PtrList<volScalarField>& Y_;
	std::vector<double> W_;
	std::vector<std::string> species_names_;
};

Foam::autoPtr<laminarSMOKEthermoClass> laminarSMOKEthermoClass::New ( const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names )
{
    IOobject radIO
    (
        "laminarSMOKEthermoClass",
        T.time().constant(),
        T.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    );

    return autoPtr<laminarSMOKEthermoClass>(new laminarSMOKEthermoClass(T,p,Y,W,species_names));
}

Foam::IOobject laminarSMOKEthermoClass::createIOobject ( const fvMesh& mesh ) const
{
    IOobject io
    (
        "laminarSMOKEthermoClass",
        mesh.time().constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


laminarSMOKEthermoClass::laminarSMOKEthermoClass(const volScalarField& T, const volScalarField& p, const PtrList<volScalarField>& Y, const std::vector<double>& W, const std::vector<std::string>& species_names) :
	IOdictionary
	    (
		IOobject
		(
		    "laminarSMOKEthermoClass",
		    T.time().constant(),
		    T.mesh(),
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		)
	    ),
	T_(T),
	p_(p),
	Y_(Y),
	W_(W),
	species_names_(species_names)
{
	Info << "Created object laminarSMOKEthermoClass" << endl;
}

unsigned int laminarSMOKEthermoClass::species_index(const std::string name) const
{
	for (unsigned int i=0;i<species_names_.size();i++)
	{
		if (name == species_names_[i])
			return i;
	}
	
	Info << "Species " << name << " not available!" << endl;
	abort(FatalError);

	return 0; 
}

laminarSMOKEthermoClass::~laminarSMOKEthermoClass()
{}

}

#endif
