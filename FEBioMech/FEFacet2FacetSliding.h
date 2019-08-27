#pragma once

#include "FEContactInterface.h"
#include "FEContactSurface.h"
#include "FECore/FEClosestPointProjection.h"
#include "FECore/vector.h"
//-----------------------------------------------------------------------------
//! Contact surface for facet-to-facet sliding interfaces
//class FEFacetSlidingSurface : public FEContactSurface


//-----------------------------------------------------------------------------
//! Sliding interface with facet to facet integration

//! This class is similar to the sliding interface except that it uses
//! a Gaussian quadrature rule in stead of a nodal integration rule
//! as its base class does.
//
//class FEFacet2FacetSliding : public FEContactInterface

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#include "FEContactSurface.h"
//#include "FEContactInterface.h"
//#include "FECore/FEClosestPointProjection.h"


//-----------------------------------------------------------------------------
class FEFacetSlidingSurface : public FEContactSurface
{
public:
	//! constructor
	//friend class FEFacet2FacetSliding;
	FEFacetSlidingSurface(FEMesh* pm = 0) : FEContactSurface(pm) {}

	//! Initializes data structures
	bool Init();

	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	//! evaluate net contact force
	vec3d GetContactForce();

	//! evaluate net contact area
	double GetContactArea();

	

	double GetFractureEnergy();
	double GetFrictionEnergy();
	double GetCriticalStrainEnergyReleaseRate();
	// Serialize data to 
	void Serialize(DumpFile& ar);

public:
	void GetNodalContactGap(int nface, double* pg);
	//void GetNodelTemperature(int nface, double* pg);

	void GetNodalContactOpeningGap(int nface, double* pg);


	void GetNodalContactPressure(int nface, double* pg);
	void GetNodalContactTraction(int nface, vec3d* pt);


public:

	//FEFacet2FacetSliding  sel ;
	vector<double>				m_gap;	//!< gap function at nodes
	vector<double>				D_temp;	//!< gap function at nodes
	vector<double>				D_R;	//!< gap function at nodes
	vector<double>				D_R_C;	//!< gap function at nodes
	vector<double>				 MAX;
	vector<double>				m_gap_p;	//!< gap function at nodes(previous time step)
	vector<bool>				DAM;
	vector<vec2d>				m_tan_gap_P;//!< tangential gap function at nodes(previous time step)
	vector<vec2d>				m_tan_gap;//!< tangential gap function at nodes
	//double DTCM = 0;
	//double  = 0;
	vector<double>				DTCM;	//!< gap function at nodes
	vector<double>				DT;	//!< gap function at nodes

	vector<double>				m_MSAT;	//!< temperature at nodes(Master,Slave)
	//vector<double>				m_SLT;	//!< temperature at nodes(Slave)
	//vector<double>				mt1_gap;	//!< gap function at nodes(tangential ONE)
	//vector<double>				mt2_gap;	//!< gap function at nodes(tangential TWO)
	vector<vec3d>				m_nu;	//!< master normal at slave node
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec2d>				m_rsp;	//!< natural coordinates at previous time step
	vector<vec2d>				m_rsp_f;	//!< natural coordinates at previous time step
	vector<vec2d>				m_Lt_P_C;
	vector<double>				m_Lm;	//!< Lagrange multipliers for contact pressure
	//vector<double>				m_Lmtemp;	//!< Lagrange multipliers for contact pressure
	vector<mat2d>				m_M;	//!< surface metric tensor
	vector<vec2d>				m_Lt;	//!< Lagrange multipliers for tangential

	vector<vec2d>				m_Lt_C;	//!< Lagrange multipliers for tangential

	vector<vec2d>				m_Lt_P;	//!< Lagrange multipliers for tangential(previous time step)
	vector<vec2d>				m_pt;	//!<net contact pressure for tangential
	vector<vec2d>				m_pt_P;	//!<net contact pressure for tangential(previous time step)
	vector<vec2d>				m_trans_S;  //  transformation needs to be done in curvilinear  
	vector<vec2d>				m_control;  //  transformation needs to be done in curvilinear
	vector<double>				m_off;	//!< gap offset (= shell thickness)
	vector<double>				m_eps;	//!< normal penalty factors
	vector<double>				m_Ln;	//!< net contact pressure normal    ///double max_normal = M_TN;
	vector<double>				m_Ln_max;
	vector<double>				m_Ln_P;	//!< net contact pressure normal(previous time step)
	vector<double>				m_Lnp;	//!< net contact pressure normal
	double TFractureEnergy;
	double TFrictionEnergy;
	double TCriticalStrainEnergyReleaseRate;
	
	//double  SetFractureEnergy;
	//double  SetFrictionEnergy;
	//FEMesh& meshS ;
	//friend class FEFacet2FacetSliding;
};

//-----------------------------------------------------------------------------
//! This class implements a sliding interface

//! The FESlidingInterface class defines an interface between two surfaces.
//! These surfaces define the slave and master surfaces of the contact
//! interface.
//double CalcFractureEnergy(){ return FractureEnergy; }

class FEFacet2FacetSliding : public FEContactInterface
{

	//friend class FEFacetSlidingSurface;
public:
	//! constructor
	FEFacet2FacetSliding(FEModel* pfem);
	//FEFacet2FacetSliding();
	//! destructor
	virtual ~FEFacet2FacetSliding(){}

	//! Initializes sliding interface
	bool Init();

	//! interface activation
	void Activate();

	//double CalcFractureEnergy();
	//double  CalcFrictionEnergy();
	//! update interface data
	void Update(int niter);

	//! projects slave nodes onto master nodes
	void ProjectSurface(FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, bool bupseg, bool bmove = false);

	//! shallow copy
	void ShallowCopy(DumpStream& dmp, bool bsave);

	//! calculate penalty value
	double Penalty() { return m_eps; }

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface() { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEStiffnessMatrix& K);

protected:
	//! calculate auto penalty factor
	void CalcAutoPenalty(FEFacetSlidingSurface& s);

	//! calculate the nodal force of a slave node
	//void ContactNodalForce(int m, FEFacetSlidingSurface& ss, FESurfaceElement& mel, vector<double>& fe);
	void FEFacet2FacetSliding::ContactNodalForce(int m, FEFacetSlidingSurface& ss, FESurfaceElement& mel, vector<double>& fe);

	//! calculate the stiffness contribution of a single slave node
	void ContactNodalStiffness(int m, FEFacetSlidingSurface& ss, FESurfaceElement& mel, matrix& ke);

	

	//! map the frictional data from the old element to the new element
	void MapFrictionData(int inode, FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, FESurfaceElement& sn, FESurfaceElement& so, vec3d& q);
	void MapTangentialComponent(int inode, FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, FESurfaceElement &en, FESurfaceElement &eo, vec3d &q);

private:
	//! copy constructor hidden
	FEFacet2FacetSliding(FEFacet2FacetSliding & si){}

public:
	FEFacetSlidingSurface	m_ss;	//!< slave surface
	FEFacetSlidingSurface	m_ms;	//!< master surface
	
	bool			m_btwo_pass;	//!< two pass algorithm flag


	//int				m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmax;
	int				m_naugmin;	//!< minimum nr of augmentations
	double			m_gtol;		//!< gap tolerance
	double			m_atol;		//!< augmentation tolerance
	double			m_ktmult;	//!< multiplier for tangential stiffness
	double			m_knmult;	//!< multiplier for normal stiffness


	double M_TN;
	double M_Tt;
	double δ_TN;
	double δ_Tt;
	double α;
	double β;
	double			m_stol;		//!< search tolerance

	bool			m_bautopen;	//!< auto penalty calculation
	double			m_eps;		//!< penalty scale factor 

	bool			m_breloc;	//!< initial node relocation

	double			m_mu;		//!< friction coefficient
	double			m_epsf;		//!< penalty scale factor for friction


	double			max_pent;		//!< maximum allowed penetration (for calculate the friction-Stick)

	double			fric_tol;   //!< friction tolerance
	int				m_nsegup;	//!< segment update parameter
	//heat capacity m^2/s^2K
	double C_heat;

	//default  heat transfer coefficients
	double γ1_transfer;
	double γ2_transfer;
	// default  thermal conductivity
	double κ_conductivity;
	// heat quantities
	double R[6];

	// contact reference temperatures
	double T0_Contact;
	// // contact surface  temperatures at contact  approach time (TC1 , TC2)
	double TCN_Contact;
	double κ_us;
	/// fluidity parameter
	double η;

private:
	bool	m_bfirst;	//!< flag to indicate the first time we enter Update
	double	m_normg0;	//!< initial gap norm

	

public:
	DECLARE_PARAMETER_LIST();
};

