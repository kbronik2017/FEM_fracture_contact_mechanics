#pragma once
#include "FECore/FESolidDomain.h"
#include "FECore/FEModel.h"
#include "FEElasticDomain.h"

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FEElasticSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEElasticSolidDomain(FEMesh* pm, FEMaterial* pmat) : FESolidDomain(FE_SOLID_DOMAIN, pm, pmat) {}

	//! \todo Do I really use this?
	FEElasticSolidDomain& operator = (FEElasticSolidDomain& d) { m_Elem = d.m_Elem; m_pMesh = d.m_pMesh; return (*this); }

	//! initialize class
	bool Initialize(FEModel& fem);

	//! Init elements
	void InitElements();

	//! reset element data
	void Reset();

	//! Unpack solid element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public: // overrides from FEElasticDomain

	// update stresses
	void UpdateStresses(FEModel& fem);

	// update the element stress
	void UpdateElementStress(int iel, double dt);

	//! calculates the residual
//	void Residual(FESolver* psolver, vector<double>& R);

	//! intertial forces for dynamic problems
	void InertialForces(FEGlobalVector& R, vector<double>& F);

	//! internal stress forces
	void InternalForces(FEGlobalVector& R);

	//! body forces
	void BodyForce(FEGlobalVector& R, FEBodyForce& BF);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

	//! calculates inertial stiffness
	void MassMatrix(FESolver* psolver, double scale);
	///////////
	 void ConductionMatrix(FESolver* pnls) ;
	void  CapacitanceMatrix(FESolver* pnls, double dt) ;
	void ElementCapacitance(FESolidElement &el, matrix &ke, double dt);
	void ElementConduction(FESolidElement& el, matrix& ke);
	//void UnpackLM1(FEElement& el, vector<int>& lm);

	void HeatResidual(FEGlobalVector& R);

	void HeatElementResidual( FESolidElement& el, vector<double>& fe);


	void ConvectiveHeatStiffnessMatrix(FESolver* psolver){};

	



	///////////
	//! body force stiffness
	void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf);

public:
	// --- S T I F F N E S S ---

	//! calculates the solid element stiffness matrix
	virtual void ElementStiffness(FEModel& fem, int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(FESolidElement& el, matrix& ke);

	//! material stiffness component
	void ElementMaterialStiffness(FESolidElement& el, matrix& ke);

	//! calculates the solid element mass matrix
	void ElementMassMatrix(FESolidElement& el, matrix& ke, double a);

	//! calculates the stiffness matrix due to body forces 
	void ElementBodyForceStiffness(FEBodyForce& bf, FESolidElement& el, matrix& ke);

	void ConvectiveHeatElementStiffness(FESurfaceElement& el, matrix& ke, double hc){};



	// --- R E S I D U A L ---

	//! Calculates the internal stress vector for solid elements
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);

	//! Calculatess external body forces for solid elements
//	void ElementBodyForce(FEModel& fem, FESolidElement& elem, vector<double>& fe);

	//! Calculatess external body forces for solid elements
	void ElementBodyForce(FEBodyForce& BF, FESolidElement& elem, vector<double>& fe);

	// --- so far just default value (To do we need to change the PreView)
	double	m_k=1.0;	//!< heat conductivity
	double	m_c=1.0;	//!< heat capacitance
	//double	m_rho;	//!< density
	double m_Q_H = 0.0;
	double Capacitance() { return m_c; }
	void Conductivity(double D[3][3]);

private:
	//! Helper function for setting the material point's local coordinate system
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp, FEElasticMaterial* pme);
};
