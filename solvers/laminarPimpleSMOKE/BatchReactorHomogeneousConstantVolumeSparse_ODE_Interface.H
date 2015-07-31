#ifndef BatchReactorHomogeneousConstantVolumeSparse_ODE_Interface_H
#define	BatchReactorHomogeneousConstantVolumeSparse_ODE_Interface_H

#include "math/OpenSMOKEVector.h"
#include "ode/OpenSMOKE_OdeSystemObject.h"

#include "ode/OpenSMOKE_OdeSystemObject.h"
#include "ODESystemVirtualClassWithOpenSMOKEVectors.H"
#include "BatchReactorHomogeneousConstantVolumeSparse.H"

namespace OpenSMOKE
{
	class BatchReactorHomogeneousConstantVolumeSparse_ODE_OpenSMOKE
	{

	public:

		BatchReactorHomogeneousConstantVolumeSparse_ODE_OpenSMOKE() {};

		void SetReactor(BatchReactorHomogeneousConstantVolumeSparse* reactor) 
		{
			reactor_ = reactor; 
		}

	protected:

		unsigned int ne_;

		void MemoryAllocation()
		{
			OpenSMOKE::ChangeDimensions(ne_, &y_, true);
			OpenSMOKE::ChangeDimensions(ne_, &dy_, true);
		}

		virtual void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd& DY)
		{
			y_.CopyFrom(Y.data());
			reactor_->Equations(t, y_, dy_);
			dy_.CopyTo(DY.data());
		}

		void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::SparseMatrix<double> &J)
		{
			y_.CopyFrom(Y.data());
			reactor_->SparseAnalyticalJacobian(t, y_, J);
		}

		virtual void Print(const double t, const Eigen::VectorXd &Y)
		{
		}

	private:

		BatchReactorHomogeneousConstantVolumeSparse* reactor_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble dy_;
	};
	
}

#endif	// BatchReactorHomogeneousConstantVolumeSparse_ODE_Interface_H
