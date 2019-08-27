#include "stdafx.h"
#include "FEFacet2FacetSliding.h"
#include "FEStiffnessMatrix.h"
#include "FECore/FEModel.h"
#include "FECore/FEClosestPointProjection.h"
#include "FECore/log.h"
////////////////


//#include "FESlidingInterface.h"
#include "FEElasticShellDomain.h"



//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FEFacet2FacetSliding, FEContactInterface)
ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL, "laugon");
ADD_PARAMETER(m_atol, FE_PARAM_DOUBLE, "tolerance");
ADD_PARAMETER(m_eps, FE_PARAM_DOUBLE, "penalty");
ADD_PARAMETER(m_bautopen, FE_PARAM_BOOL, "auto_penalty");
ADD_PARAMETER(m_btwo_pass, FE_PARAM_BOOL, "two_pass");
ADD_PARAMETER(m_gtol, FE_PARAM_DOUBLE, "gaptol");
ADD_PARAMETER(m_mu, FE_PARAM_DOUBLE, "fric_coeff");
ADD_PARAMETER(m_epsf, FE_PARAM_DOUBLE, "fric_penalty");
ADD_PARAMETER(m_naugmin, FE_PARAM_INT, "minaug");
ADD_PARAMETER(m_naugmax, FE_PARAM_INT, "maxaug"); 
ADD_PARAMETER(m_stol, FE_PARAM_DOUBLE, "search_tol");
ADD_PARAMETER(max_pent, FE_PARAM_DOUBLE, "max_penetration");

ADD_PARAMETER(fric_tol, FE_PARAM_DOUBLE, "fri_tolerance");

ADD_PARAMETER(M_TN, FE_PARAM_DOUBLE, "n_interf_strength");
ADD_PARAMETER(M_Tt, FE_PARAM_DOUBLE, "t_interf_strength");
ADD_PARAMETER(δ_TN, FE_PARAM_DOUBLE, "n_char_op");
ADD_PARAMETER(δ_Tt, FE_PARAM_DOUBLE, "t_char_op");

	ADD_PARAMETER(C_heat, FE_PARAM_DOUBLE, "heat_capacity");
	ADD_PARAMETER(γ1_transfer, FE_PARAM_DOUBLE, "s_heat__coeff");
	ADD_PARAMETER(γ2_transfer, FE_PARAM_DOUBLE, "m_heat__coeff");
	ADD_PARAMETER(T0_Contact, FE_PARAM_DOUBLE, "c_ini_temperature");
	ADD_PARAMETER(TCN_Contact, FE_PARAM_DOUBLE, "c_spat_temperature");
	ADD_PARAMETER(κ_conductivity, FE_PARAM_DOUBLE, "thermal_conductivity");
	ADD_PARAMETER(η, FE_PARAM_DOUBLE, "flu_parameter");
	ADD_PARAMETER(κ_us, FE_PARAM_DOUBLE, "user_def"); 
	ADD_PARAMETER(α, FE_PARAM_DOUBLE, "user_def_alpfa"); 
	ADD_PARAMETER(β, FE_PARAM_DOUBLE, "user_def_beta");

ADD_PARAMETER(m_ktmult, FE_PARAM_DOUBLE, "ktmult");
ADD_PARAMETER(m_knmult, FE_PARAM_DOUBLE, "knmult");
ADD_PARAMETER(m_breloc, FE_PARAM_BOOL, "node_reloc");
ADD_PARAMETER(m_nsegup, FE_PARAM_INT, "seg_up");
END_PARAMETER_LIST();




//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix


//! build the matrix profile for use in the stiffness matrix
void FEFacet2FacetSliding::BuildMatrixProfile(FEStiffnessMatrix& K)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	const int LMSIZE = 7 * (FEElement::MAX_NODES + 1);
	vector<int> lm(LMSIZE);
	int npass = (m_btwo_pass ? 2 : 1);
	for (int np = 0; np < npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0 ? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0 ? m_ms : m_ss);
		for (int j = 0; j < ss.Nodes(); ++j)
		{
			FEElement* pe = ss.m_pme[j];
			if (pe != 0)
			{
				FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pe);
				int* en = &me.m_node[0];

				int n = me.Nodes();
				lm.assign(LMSIZE, -1);

				lm[0] = ss.Node(j).m_ID[DOF_X];
				lm[1] = ss.Node(j).m_ID[DOF_Y];
				lm[2] = ss.Node(j).m_ID[DOF_Z];
				lm[3] = ss.Node(j).m_ID[DOF_RU];
				lm[4] = ss.Node(j).m_ID[DOF_RV];
				lm[5] = ss.Node(j).m_ID[DOF_RW];
				lm[6] = ss.Node(j).m_ID[DOF_T];

				for (int k = 0; k < n; ++k)
				{
					vector<int>& id = mesh.Node(en[k]).m_ID;
					lm[7 * (k + 1)] = id[DOF_X];
					lm[7 * (k + 1) + 1] = id[DOF_Y];
					lm[7 * (k + 1) + 2] = id[DOF_Z];
					lm[7 * (k + 1) + 3] = id[DOF_RU];
					lm[7 * (k + 1) + 4] = id[DOF_RV];
					lm[7 * (k + 1) + 5] = id[DOF_RW];
					lm[7 * (k + 1) + 6] = id[DOF_T];
				}

				K.build_add(lm);
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FEFacetSlidingSurface::Init()
{
	
	int i = 0, j = 0, n = 0;

	// always intialize base class first!
	if (FEContactSurface::Init() == false) return false;

	// make sure the sibling surface has been set
	assert(m_pSibling);

	// get the number of nodes
	int nn = Nodes();

	// int el = Elements();
	//SetFractureEnergy=0.0;
	//SetFrictionEnergy=0.0;
	// allocate other surface data
	m_gap.assign(nn, 0.0);	// gap funtion
	D_temp.assign(nn, 1.0);
	D_R.assign(nn, 1.0);
	D_R_C.assign(nn, 1.0);
	MAX.assign(nn, 0.0);
	m_gap_p.assign(nn, 0.0);
	DAM.assign(nn, false);
	
	//double DTCM = 0;
	//double  = 0;
	DTCM.assign(nn, 0.0);
	DT.assign(nn, 0.0);


	//	mt1_gap.assign(nn, 0.0);	// gap funtion(tangential)
	//	mt2_gap.assign(nn, 0.0);	// gap funtion(tangential)
	//m_nu.resize(nn);			// node normal 
	m_pme.assign(nn, static_cast<FESurfaceElement*>(0));		// penetrated master element
	//m_rs.resize(nn);			// natural coords of projected slave node on master element
	//m_rsp.resize(nn);
	//m_rsp_f.resize(nn);
	m_nu.assign(nn, vec3d(0, 0, 0));
	m_rs.assign(nn, vec2d(0, 0));
	m_rsp.assign(nn, vec2d(0, 0));
	m_rsp_f.assign(nn, vec2d(0, 0));

	//m_Lt_P_C.resize(nn);
	m_Lm.assign(nn, 0.0);
	TFractureEnergy=0.0;
	TFrictionEnergy=0.0;
	TCriticalStrainEnergyReleaseRate = 0.0;
	
	//m_Lmtemp.assign(nn, 0.0);
	//m_M.resize(nn);
	m_M.assign(nn, mat2d(0, 0, 0, 0));
	m_Lt.assign(nn, vec2d(0, 0));
	m_Lt_C.assign(nn, vec2d(0, 0));
	m_Lt_P_C.assign(nn, vec2d(0, 0));
	m_trans_S.assign(nn, vec2d(0, 0));
	m_control.assign(nn, vec2d(0, 0));
	m_Lt_P.assign(nn, vec2d(0, 0));
	m_pt.assign(nn, vec2d(0, 0));
	m_pt_P.assign(nn, vec2d(0, 0));
	m_tan_gap_P.assign(nn, vec2d(0, 0));
	m_tan_gap.assign(nn, vec2d(0, 0));
	m_off.assign(nn, 0.0);
	m_eps.assign(nn, 1.0);
	m_Ln.assign(nn, 0.0); 
	m_Ln_max.assign(nn, 0.0);
	m_Ln_P.assign(nn, 0.0);
	m_Lnp.assign(nn, 0.0);
	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	// note that we force rigid shells to have zero thickness
	FEMesh& m = *m_pMesh;
	vector<double> tag(m.Nodes());
	zero(tag);
	for (int nd = 0; nd<m.Domains(); ++nd)
	{
		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m.Domain(nd));
		if (psd)
		{
			for (i = 0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				n = el.Nodes();
				for (j = 0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
			}
		}
	}
	for (i = 0; i<nn; ++i) m_off[i] = tag[m_node[i]];

	return true;
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_Lm;
		//	dmp << m_Lmtemp;
		dmp << m_gap;
		dmp << m_Lt;
		dmp << m_Ln;
		dmp << m_Lnp;
		dmp << m_pt;

	}
	else
	{
		zero(m_pme);
		dmp >> m_Lm;
		dmp >> m_gap;
		dmp >> m_Lt;
		dmp >> m_Ln;
		dmp >> m_Lnp;
		dmp >> m_pt;
		//	dmp >> m_Lmtemp;
	}
}

//-----------------------------------------------------------------------------
//! 
vec3d FEFacetSlidingSurface::traction(int inode)
{
	vec3d t(0, 0, 0);
	if (m_pme[inode])
	{
		FESurfaceElement& el = *m_pme[inode];
		double Tn = m_Ln[inode];
		double Tnp = -m_Lnp[inode];
		double T1 = m_pt[el.m_lnode[inode]][0];
		double T2 = m_pt[el.m_lnode[inode]][1];
		//if (m_mu*m_epsf > 0){
		double TT1 = m_Lt[inode][0];
		double TT2 = m_Lt[inode][1];
		//}
		double r = m_rs[inode][0];
		double s = m_rs[inode][1];

		vec3d tn = m_nu[inode] * Tn, tt;
		vec3d tnp = m_nu[inode] * Tnp;
		vec3d e[2];
		ContraBaseVectors(el, r, s, e);
		tt = e[0] * T1 + e[1] * T2 + e[0] * TT1 + e[1] * TT2;
		t = tn + tt + tnp;
	}

	return t;
}



double FEFacetSlidingSurface::GetFrictionEnergy(){

	int n, i;
	const int MN = FEElement::MAX_NODES;
	double T1[MN], T2[MN], T1C[MN], T2C[MN], ΔtangentialgapX[MN], ΔtangentialgapY[MN];

	// initialize FrictionEnergy

	double E = 0;


	// loop over all elements of the surface
	for (n = 0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nseln = el.Nodes();

		// nodal contact pressures 
		for (i = 0; i<nseln; ++i) {

			//Δnormalgap[i] = m_gap[el.m_lnode[i]] - m_gap_p[el.m_lnode[i]];
			ΔtangentialgapX[i] = m_tan_gap[el.m_lnode[i]][0] - m_tan_gap_P[el.m_lnode[i]][0];
			ΔtangentialgapY[i] = m_tan_gap[el.m_lnode[i]][1] - m_tan_gap_P[el.m_lnode[i]][1];
			////  Trapezoidal rule

			//Tn[i] = fabs(((m_Ln[el.m_lnode[i]] + m_Ln_P[el.m_lnode[i]]) / 2) * Δnormalgap[i]);
			T1[i] = fabs(((m_Lt[el.m_lnode[i]][0] + m_Lt_P[el.m_lnode[i]][0]) / 2) * ΔtangentialgapX[i]);
			T2[i] = fabs(((m_Lt[el.m_lnode[i]][1] + m_Lt_P[el.m_lnode[i]][1]) / 2)*ΔtangentialgapY[i]);

			T1C[i] = fabs(((m_Lt_C[el.m_lnode[i]][0] + m_Lt_P_C[el.m_lnode[i]][0]) / 2) * ΔtangentialgapX[i]);
			T2C[i] = fabs(((m_Lt_C[el.m_lnode[i]][1] + m_Lt_P_C[el.m_lnode[i]][1]) / 2)*ΔtangentialgapY[i]);




		}
		int nint = el.GaussPoints();

		// evaluate the contact force for that element
		for (i = 0; i<nint; ++i)
		{
			// area in reference configuration
			vec3d g0[2], g[2];
			double r = el.gr(i);
			double s = el.gs(i);
			CoBaseVectors0(el, r, s, g0);
			double A = (g0[0] ^ g0[1]).unit();
			// traction components at integration point
			double t1 = el.eval(T1, i);
			double t2 = el.eval(T2, i);
			double t1C = el.eval(T1C, i);
			double t2C = el.eval(T2C, i);


			// unit normal vector
			vec3d normal = SurfaceNormal(el, i);
			// contravariant basis in spatial frame
			ContraBaseVectors(el, r, s, g);
			// Piola traction
			double t = (g[0] * t1)*g[0] + (g[1] * t2)*g[1] + (g[0] * t1C)*g[0] + (g[1] * t2C)*g[1];
			//double totall = t.norm();
			// gauss weight
			double w = el.GaussWeights()[i];
			// contact force
			E += fabs(t*(w*A));
		}
	}
	TFrictionEnergy += E;
	

	return TFrictionEnergy;

}
////
double FEFacetSlidingSurface::GetCriticalStrainEnergyReleaseRate(){

	int n, i;
	const int MN = FEElement::MAX_NODES;
	double Tn[MN], T1[MN], T2[MN], Δnormalgap[MN], ΔtangentialgapX[MN], ΔtangentialgapY[MN];

	// initialize FractureEnergy

	double E = 0;


	// loop over all elements of the surface
	for (n = 0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nseln = el.Nodes();

		// nodal contact pressures 
		for (i = 0; i<nseln; ++i) {

			Δnormalgap[i] = m_gap[el.m_lnode[i]] - m_gap_p[el.m_lnode[i]];
			ΔtangentialgapX[i] = m_tan_gap[el.m_lnode[i]][0] - m_tan_gap_P[el.m_lnode[i]][0];
			ΔtangentialgapY[i] = m_tan_gap[el.m_lnode[i]][1] - m_tan_gap_P[el.m_lnode[i]][1];
			////  Trapezoidal rule

			Tn[i] = fabs(((m_Ln[el.m_lnode[i]] + m_Ln_P[el.m_lnode[i]]) / 2) * Δnormalgap[i]);
			T1[i] = fabs(((m_pt[el.m_lnode[i]][0] + m_pt_P[el.m_lnode[i]][0]) / 2) * ΔtangentialgapX[i]);
			T2[i] = fabs(((m_pt[el.m_lnode[i]][1] + m_pt_P[el.m_lnode[i]][1]) / 2)*ΔtangentialgapY[i]);
		}
		int nint = el.GaussPoints();

		// evaluate the contact force for that element
		for (i = 0; i<nint; ++i)
		{
			// area in reference configuration
			vec3d g0[2], g[2];
			double r = el.gr(i);
			double s = el.gs(i);
			CoBaseVectors0(el, r, s, g0);
			double A = (g0[0] ^ g0[1]).unit();
			// traction components at integration point
			double t1 = el.eval(T1, i);
			double t2 = el.eval(T2, i);
			double t3 = el.eval(Tn, i);


			// unit normal vector
			vec3d normal = SurfaceNormal(el, i);
			// contravariant basis in spatial frame
			ContraBaseVectors(el, r, s, g);
			// Piola traction
			double t = (g[0] * t1)*g[0] + (g[1] * t2)*g[1] + (normal*t3)*normal;
			//double totall = t.norm();
			// gauss weight
			double w = el.GaussWeights()[i];
			// contact force
			E += fabs(t);
		}
	}
	TCriticalStrainEnergyReleaseRate += E;


	return TCriticalStrainEnergyReleaseRate;

}

//// StrainEnergyReleaseRate*area

double FEFacetSlidingSurface::GetFractureEnergy(){

	int n, i;
	const int MN = FEElement::MAX_NODES;
	double Tn[MN], T1[MN], T2[MN], Δnormalgap[MN], ΔtangentialgapX[MN], ΔtangentialgapY[MN];

	// initialize FractureEnergy
	
	double E=0;


	// loop over all elements of the surface
	for (n = 0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nseln = el.Nodes();

		// nodal contact pressures 
		for (i = 0; i<nseln; ++i) {

			Δnormalgap[i] = m_gap[el.m_lnode[i]] - m_gap_p[el.m_lnode[i]];
			ΔtangentialgapX[i] = m_tan_gap[el.m_lnode[i]][0] - m_tan_gap_P[el.m_lnode[i]][0];
			ΔtangentialgapY[i] = m_tan_gap[el.m_lnode[i]][1] - m_tan_gap_P[el.m_lnode[i]][1];
			////  Trapezoidal rule

			Tn[i] = fabs(((m_Ln[el.m_lnode[i]] + m_Ln_P[el.m_lnode[i]])/2) * Δnormalgap[i]);
			T1[i] = fabs(((m_pt[el.m_lnode[i]][0] + m_pt_P[el.m_lnode[i]][0])/2) * ΔtangentialgapX[i]);
			T2[i] = fabs(((m_pt[el.m_lnode[i]][1] + m_pt_P[el.m_lnode[i]][1])/2)*ΔtangentialgapY[i]);
		}
		int nint = el.GaussPoints();

		// evaluate the contact force for that element
		for (i = 0; i<nint; ++i)
		{
			// area in reference configuration
			vec3d g0[2], g[2];
			double r = el.gr(i);
			double s = el.gs(i);
			CoBaseVectors0(el, r, s, g0);
			double A = (g0[0] ^ g0[1]).unit();
			// traction components at integration point
			double t1 = el.eval(T1, i);
			double t2 = el.eval(T2, i);
            double t3 = el.eval(Tn, i);


			// unit normal vector
			vec3d normal = SurfaceNormal(el, i);
			// contravariant basis in spatial frame
			ContraBaseVectors(el, r, s, g);
			// Piola traction
			double t = (g[0] * t1)*g[0] + (g[1] * t2)*g[1] + (normal*t3)*normal;
			//double totall = t.norm();
			// gauss weight
			double w = el.GaussWeights()[i];
			// contact force
			E += fabs(t*(w*A));
		}
	}
	TFractureEnergy += E;
	

	return TFractureEnergy;

}




//-----------------------------------------------------------------------------
vec3d FEFacetSlidingSurface::GetContactForce()
{

	int n, i;
	const int MN = FEElement::MAX_NODES;
	double Tn[MN], T1[MN], T2[MN], Tnp[MN], TT1[MN], TT2[MN];

	// initialize contact force
	vec3d f(0, 0, 0);



	// loop over all elements of the surface
	for (n = 0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nseln = el.Nodes();

		// nodal contact pressures and frictional tractions
		for (i = 0; i<nseln; ++i) {
			Tn[i] = m_Ln[el.m_lnode[i]];
			Tnp[i] = -m_Lnp[el.m_lnode[i]];
			TT1[i] = m_Lt[el.m_lnode[i]][0];
			TT2[i] = m_Lt[el.m_lnode[i]][1];
			//T1[i] = m_Lt[el.m_lnode[i]][0];
			//T2[i] = m_Lt[el.m_lnode[i]][1];

			////////////
			T1[i] = m_pt[el.m_lnode[i]][0];
			T2[i] = m_pt[el.m_lnode[i]][1];
			/////////////
			/// change here tom.
		}
		int nint = el.GaussPoints();

		// evaluate the contact force for that element
		for (i = 0; i<nint; ++i)
		{
			// area in reference configuration
			vec3d g0[2], g[2];
			double r = el.gr(i);
			double s = el.gs(i);
			CoBaseVectors0(el, r, s, g0);
			double A = (g0[0] ^ g0[1]).unit();
			// traction components at integration point
			double t1 = el.eval(T1, i);
			double t2 = el.eval(T2, i);
			double t5 = el.eval(TT1, i);
			double t6 = el.eval(TT2, i);
			double t3 = el.eval(Tn, i);
			double t4 = el.eval(Tnp, i);
			// unit normal vector
			vec3d normal = SurfaceNormal(el, i);
			// contravariant basis in spatial frame
			ContraBaseVectors(el, r, s, g);
			// Piola traction
			vec3d t = g[0] * t1 + g[1] * t2 + normal*t3 + normal*t4 + g[0] * t5 + g[1] * t6;
			// gauss weight
			double w = el.GaussWeights()[i];
			// contact force
			f += t*(w*A);
		}
	}
	//felog.printf("ContactForce");
	//felog.printf("ContactForce  f.x ,f.y ,f.z are:  %.3f  %.3f  %.3f\n", f.x, f.y, f.z);
	
	return f;

}

//-----------------------------------------------------------------------------
double FEFacetSlidingSurface::GetContactArea()
{
	const int MN = FEElement::MAX_NODES;
	double Tn[MN];

	// initialize contact area
	double a = 0;

	// loop over all elements of the primary surface
	for (int n = 0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nint = el.GaussPoints();

		int nseln = el.Nodes();

		// nodal contact pressures
		for (int i = 0; i<nseln; ++i) {
			Tn[i] = m_Ln[el.m_lnode[i]];
			//Tn[i] = m_Ln[el.m_lnode[i]];
		}

		// evaluate the contact force for that element
		for (int i = 0; i<nint; ++i)
		{
			// get data for this integration point
			double Ln = el.eval(Tn, i);
			double s = (Ln < 0) ? 1 : 0;



			//

			// get the base vectors
			vec3d g[2];
			CoBaseVectors(el, i, g);

			// normal (magnitude = area)
			vec3d normal = g[0] ^ g[1];

			// gauss weight
			double w = el.GaussWeights()[i];

			// contact force
			a += normal.norm()*(w*s);
		}
	}

	return a;
}

//


//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::Serialize(DumpFile& ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_gap;
		ar << m_nu;
		ar << m_rs;
		ar << m_rsp_f;
		ar << m_rsp;
		ar << m_Lm;
		//ar << m_Lmtemp;
		ar << m_M;
		ar << m_Lt;
		ar << m_off;
		ar << m_eps;
		ar << m_Ln;
		ar << m_Lnp;
		ar << m_pt;
	}
	else
	{
		// read the contact data
		// Note that we do this after Init() (called in FESurface::Serialize) since this data gets 
		// initialized to zero there
		ar >> m_gap;
		ar >> m_nu;
		ar >> m_rs;
		ar >> m_rsp;
		ar >> m_Lm;
		ar >> m_rsp_f;
		//ar >> m_Lmtemp;
		ar >> m_M;
		ar >> m_Lt;
		ar >> m_off;
		ar >> m_eps;
		ar >> m_Ln;
		ar >> m_Lnp;
		ar >> m_pt;
	}
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::GetNodalContactGap(int nface, double* pg)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.m_lnode.size();
	//for (int j = 0; j< ne; ++j) pg[j] = m_nu[f.m_lnode[j]] * m_gap[f.m_lnode[j]] + mt1_gap[f.m_lnode[j]] + mt1_gap[f.m_lnode[j]];
	for (int j = 0; j< ne; ++j) pg[j] = m_gap[f.m_lnode[j]];
}

//			ΔtangentialgapX[i] = m_tan_gap[el.m_lnode[i]][0] - m_tan_gap_P[el.m_lnode[i]][0];
//ΔtangentialgapY[i] = m_tan_gap[el.m_lnode[i]][1] - m_tan_gap_P[el.m_lnode[i]][1];

void FEFacetSlidingSurface::GetNodalContactOpeningGap(int nface, double* pg)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.m_lnode.size();
	//for (int j = 0; j< ne; ++j) pg[j] = m_nu[f.m_lnode[j]] * m_gap[f.m_lnode[j]] + mt1_gap[f.m_lnode[j]] + mt1_gap[f.m_lnode[j]];
	for (int j = 0; j< ne; ++j) pg[j] = sqrt(m_gap[f.m_lnode[j]] * m_gap[f.m_lnode[j]] + m_tan_gap[f.m_lnode[j]][0] * m_tan_gap[f.m_lnode[j]][0] + m_tan_gap[f.m_lnode[j]][1] * m_tan_gap[f.m_lnode[j]][1]);
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::GetNodalContactPressure(int nface, double* pg)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.Nodes();
	//for (int j = 0; j<ne; ++j) pg[j] = m_Ln[f.m_lnode[j]];
	for (int j = 0; j<ne; ++j) pg[j] = m_Ln[f.m_lnode[j]] + m_pt[f.m_lnode[j]][0] + m_pt[f.m_lnode[j]][1] - m_Lnp[f.m_lnode[j]] + m_Lt[f.m_lnode[j]][0] + m_Lt[f.m_lnode[j]][1];
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::GetNodalContactTraction(int nface, vec3d* tn)
{
	
	FESurfaceElement& e = Element(nface);

	int nseln = e.Nodes();
	//	int n, i;
	const int MN = FEElement::MAX_NODES;
	double Tn[MN], T1[MN], T2[MN], Tnp[MN], TT1[MN], TT2[MN];
	// nodal contact pressures and frictional tractions
	for (int i = 0; i<nseln; ++i) {
		Tn[i] = m_Ln[e.m_lnode[i]];
		Tnp[i] = -m_Lnp[e.m_lnode[i]];
		TT1[i] = m_Lt[e.m_lnode[i]][0];
		TT2[i] = m_Lt[e.m_lnode[i]][1];
		//T1[i] = m_Lt[el.m_lnode[i]][0];
		//T2[i] = m_Lt[el.m_lnode[i]][1];
		/////////////////

		////////////
		T1[i] = m_pt[e.m_lnode[i]][0];
		T2[i] = m_pt[e.m_lnode[i]][1];
		/////////////
		/// change here tom.
	}
	int nint = e.GaussPoints();

	// evaluate the contact force for that element
	for (int i = 0; i<nint; ++i)
	{
		// area in reference configuration
		vec3d g0[2], g[2];
		double r = e.gr(i);
		double s = e.gs(i);
		CoBaseVectors0(e, r, s, g0);
		double A = (g0[0] ^ g0[1]).unit();
		// traction components at integration point
		double t1 = e.eval(T1, i);
		double t2 = e.eval(T2, i);
		double t5 = e.eval(TT1, i);
		double t6 = e.eval(TT2, i);
		double t3 = e.eval(Tn, i);
		double t4 = e.eval(Tnp, i);
		// unit normal vector
		vec3d n = SurfaceNormal(e, i);
		// contravariant basis in spatial frame
		ContraBaseVectors(e, r, s, g);
		// Piola traction
		tn[i] = g[0] * t1 + g[1] * t2 + n*t3 + n*t4 + g[0] * t5 + g[1] * t6;
		// gauss weight
		//double w = el.GaussWeights()[i];
		// contact force
		//f += t*(w*A);

	}


}



// FESlidingInterface
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
//! constructor
FEFacet2FacetSliding::FEFacet2FacetSliding(FEModel* pfem) : FEContactInterface(pfem), m_ss(&pfem->GetMesh()), m_ms(&pfem->GetMesh())
{
	static int count = 1;

	m_mu = 0;
	m_epsf = 0;
	max_pent = 0.01;
	fric_tol = 0;
	// default value for heat capacity m^2/s^2K
	C_heat = 1.0;

	//default value for  heat transfer coefficients
	γ1_transfer = 1.0;
	γ2_transfer = 1.0;
	// default value for thermal conductivity
	κ_conductivity = 1.0;
	// default value for heat quantities R0,R1,R2,R3,R4,R5
	//R[5] = { 0 };
	memset(R, 0, sizeof(R));
	// default value for contact reference temperatures
	T0_Contact = 1.0;
	// default value for contact surface  temperatures at contact  approach time (TC1 , TC2)
	TCN_Contact = 1.0;

	/// fluidity parameter
	η = 1.0;
	m_naugmin = 0;
	m_naugmax = 20;

	m_gtol = 1;

	// default values
	 M_TN= 1.0;
	 M_Tt = 1.0;
	 δ_TN= 0.001;
	 δ_Tt = 0.001;
	  α = 10.0;
	  β = 1.0;
	m_stol = 0.01;
	κ_us = 1.0;
	m_ktmult = 1;
	m_knmult = 1;
	//FractureEnergy = 0.0;
	//FrictionEnergy = 0.0;
	m_breloc = false;

	m_nsegup = 0;	// always do segment updates
	m_bautopen = false;	// don't use auto-penalty
	m_btwo_pass = false; // don't use two-pass
	 m_nID = count++;

	// set the siblings
	m_ms.SetSibling(&m_ss);
	m_ss.SetSibling(&m_ms);
};

//-----------------------------------------------------------------------------
//! Calculates the auto penalty factor

void FEFacet2FacetSliding::CalcAutoPenalty(FEFacetSlidingSurface& s)
{
	int i, k, m;

	// zero penalty values
	zero(s.m_eps);

	// get the mesh
	FEMesh& mesh = *s.GetMesh();

	// get the node element list for this surface
	FENodeElemList NEL;
	NEL.Create(s);

	// loop over all surface elements
	FEElement *pe;
	for (i = 0; i<s.Elements(); ++i)
	{
		// get the next face
		FESurfaceElement& face = s.Element(i);

		// grab the element this face belongs to
		pe = mesh.FindElementFromID(face.m_nelem);
		assert(pe);

		// we need a measure for the modulus
		double K = AutoPenalty(face, s);

		// calculate the facet area
		double area = s.FaceArea(face);

		// calculate the volume element
		double vol = mesh.ElementVolume(*pe);

		// set the auto calculation factor
		double eps = 0;
		if (vol != 0){ eps = K*area / vol; }
		//double eps = K*area / vol;

		// distribute values over nodes
		for (k = 0; k<face.Nodes(); ++k)
		{
			m = face.m_lnode[k];
			s.m_eps[m] += eps;
		}
	}

	// scale values according to valence (TODO: Why are we doing this?)
	for (i = 0; i<s.Nodes(); ++i) if (NEL.Valence(i) != 0) { s.m_eps[i] /= NEL.Valence(i); }
		
}

//-----------------------------------------------------------------------------
//! Initializes the sliding interface data

bool FEFacet2FacetSliding::Init()
{
	// set data
	m_bfirst = true;
	m_normg0 = 0.0;

	// create the surfaces
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Activate()
{
	// don't forget to call the base class
	FEContactInterface::Activate();
	//felog.printf("(Cohesive) force is being calculated... ");
	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms, true, m_breloc);
	if (m_bautopen) CalcAutoPenalty(m_ss);

	// for two-pass algorithms we repeat the previous
	// two steps with master and slave switched
	if (m_btwo_pass)
	{
		ProjectSurface(m_ms, m_ss, true, m_breloc);
		if (m_bautopen) CalcAutoPenalty(m_ms);
	}

}

//-----------------------------------------------------------------------------
//!  Projects the slave surface onto the master surface.
//!  That is for each slave node we determine the closest
//!  master element and the projection of the slave node onto
//!  this master element.

//! \todo this function needs to identify the different types of node-slave contact:
//!   1/ first contact
//!   2/ crossing of element boundary
//!	  3/ contact termination 
//!			either by failure to find master segment or when g < tolerance

void FEFacet2FacetSliding::ProjectSurface(FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, bool bupseg, bool bmove)
{
	// slave node projection
	double r = 0, s = 0;
	vec3d q(0, 0, 0);

	FEClosestPointProjection cpp(ms);
	cpp.SetTolerance(m_stol);
	cpp.Init();

	// loop over all slave nodes
	for (int i = 0; i<ss.Nodes(); ++i)
	{
		// get the node
		FENode& node = ss.Node(i);

		// get the nodal position
		vec3d x = node.m_rt;

		// get the global node number
		int m = ss.m_node[i];

		// get the previous master element (if any)
		FESurfaceElement* pme = ss.m_pme[i];

		// If the node is in contact, let's see if the node still is 
		// on the same master element
		if (pme != 0)
		{
			FESurfaceElement& mel = *pme;

			r = ss.m_rs[i][0];
			s = ss.m_rs[i][1];

			q = ms.ProjectToSurface(mel, x, r, s);
			ss.m_rs[i][0] = r;
			ss.m_rs[i][1] = s;

			// we only check when we can update the segments
			// otherwise, we just stick with this element, even
			// if the node is no longer inside it.
			if (bupseg)
			{
				if (!ms.IsInsideElement(mel, r, s, m_stol) && bupseg)
				{
					// see if the node might have moved to another master element
					FESurfaceElement* pold = pme;
					ss.m_rs[i] = vec2d(0, 0);
					pme = cpp.Project(x, q, ss.m_rs[i]);

					if (pme == 0)
					{
						// nope, if has genuinly left contact
						int* n = &pold->m_node[0];
						//						log.printf("node %d has left element (%d, %d, %d, %d)\n", m+1, n[0]+1, n[1]+1, n[2]+1, n[3]+1);
					}
					else if (pme != 0)
					{
						// the node has moved to another master segment.
						// If friction is active we need to translate the frictional
						// data to the new master segment.
						FESurfaceElement& eo = *pold;
						FESurfaceElement& en = *pme;
						MapTangentialComponent(i, ss, ms, en, eo, q);
						if (m_mu*m_epsf > 0) { MapFrictionData(i, ss, ms, en, eo, q); }
					}
				}
			}
		}
		else if (bupseg)
		{
			// get the master element
			// don't forget to initialize the search for the first node!
			ss.m_rs[i] = vec2d(0, 0);
			pme = cpp.Project(x, q, ss.m_rs[i]);
			if (pme )
			{
				// the node has come into contact so make sure to initialize
				// the previous natural coordinates for friction.
				ss.m_rsp[i] = ss.m_rs[i];
				ss.m_rsp_f[i] = ss.m_rs[i];
				
			}
		}

		// if we found a master element, update the gap and normal data
		ss.m_pme[i] = pme;
		if (pme != 0)
		{
			FESurfaceElement& mel = *ss.m_pme[i];

			r = ss.m_rs[i][0];
			s = ss.m_rs[i][1];

			// if this is a new contact, copy the current coordinates
			// to the previous ones
			ss.m_M[i] = ss.Metric0(mel, r, s);

			// the slave normal is set to the master element normal
			ss.m_nu[i] = ss.SurfaceNormal(mel, r, s);

			/////
			ss.m_gap_p[i] = ss.m_gap[i];

			// calculate gap
			ss.m_gap[i] = -(ss.m_nu[i] * (x - q)) + ss.m_off[i];
			if (bmove && (ss.m_gap[i]>0))
			{
				node.m_r0 = node.m_rt = q + ss.m_nu[i] * ss.m_off[i];
				ss.m_gap[i] = 0;
			}

			// TODO: what should we do if the gap function becomes
			// negative? setting the Lagrange multipliers to zero
			// might make the system unstable.
			/*			if (ss.gap[i] < 0)
			{
			ss.Lm[i] = 0;
			ss.Lt[i][0] = 0;
			ss.Lt[i][1] = 0;
			ss.pme[i] = 0;
			}
			*/
		}
		else
		{
			// TODO: Is this a good criteria for out-of-contact?
			//		 perhaps this is not even necessary.
			// since the node is not in contact, we set the gap function 
			// and Lagrangian multiplier to zero
			ss.m_gap[i] = 0;
			ss.m_Lm[i] = 0;
			ss.m_Lt[i][0] = ss.m_Lt[i][1] = 0;
		}
	}
}

//-----------------------------------------------------------------------------
//! updates sliding interface data
//! niter is the number of Newton iterations.
void FEFacet2FacetSliding::Update(int niter)
{
	// should we do a segment update or not?
	// TODO: check what happens when m_nsegup == -1 and m_npass = 2;
	// We have to make sure that in this case, both surfaces get at least
	// one pass!
	bool bupdate = (m_bfirst || (m_nsegup == 0) ? true : (niter <= m_nsegup));

	// project slave surface onto master surface
	// this also calculates the nodal gap functions
	ProjectSurface(m_ss, m_ms, bupdate, m_breloc);
	
	if (m_btwo_pass) ProjectSurface(m_ms, m_ss, bupdate, m_breloc);

	// Update the net contact pressures
	UpdateContactPressures();
	
	// set the first-entry-flag to false
	m_bfirst = false;
}

//-----------------------------------------------------------------------------

//void FEFacet2FacetSliding::ContactForces(FEGlobalVector& R)
// with temperature evaluation on contact interfaces

void FEFacet2FacetSliding::ContactForces(FEGlobalVector& R)
{
	

	int j, k, l, m, n, np;
	int nseln = 0, nmeln = 0, ndof = 0;

	// element contact force vector
	vector<double> fe;

	// the lm array for this force vector
	vector<int> lm;

	// the en array
	vector<int> en;

	// the elements LM vectors
	vector<int> sLM;
	vector<int> mLM;

	const int MN = FEElement::MAX_NODES;
	vec3d r0[MN];
	double w[MN];
	double* Gr, *Gs;
	double detJ[MN];
	vec3d dxr, dxs;

	// do two-pass
	int npass = (m_btwo_pass ? 2 : 1);
	for (np = 0; np<npass; ++np)
	{
		// pick the slave and master surfaces
		FEFacetSlidingSurface& ss = (np == 0 ? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0 ? m_ms : m_ss);

		// loop over all slave facets
		int ne = ss.Elements();
		for (j = 0; j<ne; ++j)
		{
			// get the slave element
			FESurfaceElement& sel = ss.Element(j);
			nseln = sel.Nodes();

			// get the element's LM array
			ss.UnpackLM(sel, sLM);

			// nodal coordinates
			for (int i = 0; i<nseln; ++i) r0[i] = ss.GetMesh()->Node(sel.m_node[i]).m_r0;

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (n = 0; n<nseln; ++n)
			{
				Gr = sel.Gr(n);
				Gs = sel.Gs(n);

				// calculate jacobian
				// note that we are integrating over the reference surface
				dxr = dxs = vec3d(0, 0, 0);
				for (k = 0; k<nseln; ++k)
				{
					dxr.x += Gr[k] * r0[k].x;
					dxr.y += Gr[k] * r0[k].y;
					dxr.z += Gr[k] * r0[k].z;

					dxs.x += Gs[k] * r0[k].x;
					dxs.y += Gs[k] * r0[k].y;
					dxs.z += Gs[k] * r0[k].z;
				}

				// jacobians
				detJ[n] = (dxr ^ dxs).norm();

				// integration weights
				w[n] = sel.GaussWeights()[n];
			}

			// loop over slave element nodes (which are the integration points as well)
			// and calculate the contact nodal force
			for (n = 0; n<nseln; ++n)
			{
				// get the local node number
				m = sel.m_lnode[n];

				// see if this node's constraint is active
				// that is, if it has a master element associated with it
				// TODO: is this a good way to test for an active constraint
				// The rigid wall criteria seems to work much better.
				if (ss.m_pme[m] != 0)
				{
					// This node is active and could lead to a non-zero
					// contact force.
					// get the master element
					FESurfaceElement& mel = *ss.m_pme[m];
					ms.UnpackLM(mel, mLM);

					// calculate the degrees of freedom
					nmeln = mel.Nodes();
					//ndof = 3*(nmeln+1);
					//  X,Y,Z and T (mortar) 
					ndof = 4 * (nmeln + 1);
					// test if sequence is empty
					if (!fe.empty() == true){ fe.clear(); }

					fe.resize(ndof);

					// calculate the nodal force
					ContactNodalForce(m, ss, mel, fe);

					// multiply force with weights
					for (l = 0; l<ndof; ++l) fe[l] *= detJ[n] * w[n];

					// fill the lm array
					//lm.resize(3*(nmeln+1));
					// add Temperature dof
					//(mortar) 
					// test if sequence is empty
					if (!lm.empty() == true){ lm.clear(); }
					lm.resize(4 * (nmeln + 1));
					lm[0] = sLM[n * 3];
					lm[1] = sLM[n * 3 + 1];
					lm[2] = sLM[n * 3 + 2];
					lm[3] = sLM[10 * nseln + n];
					for (l = 0; l < nmeln; ++l)
					{
						lm[4 * (l)+4] = mLM[l * 3];
						lm[4 * (l)+5] = mLM[l * 3 + 1];
						lm[4 * (l)+6] = mLM[l * 3 + 2];
						lm[4 * (l)+7] = mLM[10 * nmeln + l];

					}

					//for (l=0; l<nmeln; ++l)
					//{
					//lm[3*(l+1)  ] = mLM[l*3  ];
					//lm[3*(l+1)+1] = mLM[l*3+1];
					//lm[3*(l+1)+2] = mLM[l*3+2];
					//}

					// fill the en array
					//(mortar)
					// test if sequence is empty
					if (!en.empty() == true){ en.clear(); }
					en.resize(nmeln + 1);
					en[0] = sel.m_node[n];
					for (l = 0; l<nmeln; ++l) en[l + 1] = mel.m_node[l];

					// assemble into global force vector
					R.Assemble(en, lm, fe);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculates the contact force on a slave node.


void FEFacet2FacetSliding::ContactNodalForce(int m, FEFacetSlidingSurface& ss, FESurfaceElement& mel, vector<double>& fe)
{
	int k, l;

	vec3d dxr, dxs;

	// normal force
	double tn, Ln;

	// gap function
	double gap;

	// tangents
	vec3d tau1, tau2;
	const int MAXMN = FEElement::MAX_NODES;
	FEAnalysis& step = *(GetFEModel()->GetCurrentStep());
	//double dt = pstep->m_dt;
	double ∆time = step.m_dt;
	double F1 = 0;
	if (∆time > 0){ F1 = (C_heat*TCN_Contact) / (T0_Contact*∆time); }
	double F2 = F1*TCN_Contact;
	double R0 = (γ1_transfer*γ2_transfer) / (γ1_transfer + γ2_transfer + F1);
	R[0] = R0 / γ2_transfer;
	R[1] = R0 / γ1_transfer;
	R[2] = F1*R[0];
	R[3] = F2*R[0];
	R[4] = F1*R[1];
	R[5] = F2*R[1];
	// get the slave element nodal temperaturesss.GetMesh()->Node(se.m_node[i])m_MSAT
	double Stempo = ss.Node(m).m_T;
	double Mtmpo[MAXMN];

	// temperatures vectors
	//double a[4 * (MAXMN + 1)], b[4 * (MAXMN + 1)], c[4 * (MAXMN + 1)], d1[4 * (MAXMN + 1)], d2[4 * (MAXMN + 1)], e[4 * (MAXMN + 1)], f[4 * (MAXMN + 1)];
	//double NTE[4 * (MAXMN + 1)];

	std::vector<double>	a(4 * (MAXMN + 1));
	std::vector<double>	b(4 * (MAXMN + 1));
	std::vector<double>	c(4 * (MAXMN + 1));
	std::vector<double>	d1(4 * (MAXMN + 1));
	std::vector<double>	d2(4 * (MAXMN + 1));
	std::vector<double>	e(4 * (MAXMN + 1));
	std::vector<double>	f(4 * (MAXMN + 1));
	std::vector<double>	NTE(4 * (MAXMN + 1));





	// master element nodes
	vec3d rtm[MAXMN];

	// shape function values
	double H[MAXMN], Hr[MAXMN], Hs[MAXMN];

	// contact vectors
	//double N[3*(MAXMN+1)], N1[3*(MAXMN+1)], N2[3*(MAXMN+1)];
	//double T1[3*(MAXMN+1)], T2[3*(MAXMN+1)], D1[3*(MAXMN+1)], D2[3*(MAXMN+1)];

	// surface metrics
	double A[2][2], M[2][2], K[2][2];
	double detA;

	double eps, scale = Penalty();

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	double Tt[2];
	//add Temperature dof
	int nmeln, ndof;

	// gap function
	gap = ss.m_gap[m];

	// normal penalty
	eps = ss.m_eps[m] * scale;

	// get slave node normal force
	Ln = ss.m_Lm[m];
	tn = Ln + eps*gap;
	tn = MBRACKET(tn);

	// get the slave node normal
	vec3d& nu = ss.m_nu[m];

	nmeln = mel.Nodes();
	//ndof = 3*(1 + nmeln);


	// metric tensors
	mat2d Mk = ss.m_M[m];
	mat2d Mki = Mk.inverse();

	// get the master element node positions
	for (k = 0; k<nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;
	// get the master elements nodal temperatures
	for (k = 0; k<nmeln; ++k) Mtmpo[k] = mesh.Node(mel.m_node[k]).m_T;


	// isoparametric coordinates of the projected slave node
	// onto the master element
	double r = ss.m_rs[m][0];
	double s = ss.m_rs[m][1];

	

	// get the coordinates at the previous step
	double rp = ss.m_rsp[m][0];
	double sp = ss.m_rsp[m][1];
	double rpf = ss.m_rsp_f[m][0];
	double spf = ss.m_rsp_f[m][1];
	// get the master shape function values at this slave node
	mel.shape_fnc(H, r, s);
	mel.shape_deriv(Hr, Hs, r, s);
	/////////////////   modified part begin (mortar Thermodynamics  case!)
	// add Temperature dof
	ndof = 4 * (1 + nmeln);

	// add Temperature dof
	//double  N1[4 * (MAXMN + 1)], N2[4 * (MAXMN + 1)];
	//double T1[4 * (MAXMN + 1)], T2[4 * (MAXMN + 1)], D1[4 * (MAXMN + 1)], D2[4 * (MAXMN + 1)];


	std::vector<double>	N1(4 * (MAXMN + 1));
	std::vector<double>	N2(4 * (MAXMN + 1));
	std::vector<double>	T1(4 * (MAXMN + 1));
	std::vector<double>	T2(4 * (MAXMN + 1));
	std::vector<double>	D1(4 * (MAXMN + 1));
	std::vector<double>	D2(4 * (MAXMN + 1));





	// set up the a vectors
	a[0] = a[1] = a[2] = 0;
	a[3] = 1;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		a[4 * (l)+4] = 0;
		a[4 * (l)+5] = 0;
		a[4 * (l)+6] = 0;
		a[4 * (l)+7] = -H[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	b[0] = b[1] = b[2] = 0;
	b[3] = R[0];
	for (l = 0; l < nmeln; ++l)
	{

		b[4 * (l)+4] = 0;
		b[4 * (l)+5] = 0;
		b[4 * (l)+6] = 0;
		b[4 * (l)+7] = H[l] * R[1];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the d1 vectors
	d1[0] = d1[1] = d1[2] = 0;
	d1[3] = 0;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		d1[4 * (l)+4] = 0;
		d1[4 * (l)+5] = 0;
		d1[4 * (l)+6] = 0;
		d1[4 * (l)+7] = -Hr[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the d2 vectors
	d2[0] = d2[1] = d2[2] = 0;
	d2[3] = 0;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		d2[4 * (l)+4] = 0;
		d2[4 * (l)+5] = 0;
		d2[4 * (l)+6] = 0;
		d2[4 * (l)+7] = -Hs[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the f vectors
	f[0] = f[1] = f[2] = 0;
	f[3] = 0;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		f[4 * (l)+4] = 0;
		f[4 * (l)+5] = 0;
		f[4 * (l)+6] = 0;
		f[4 * (l)+7] = H[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the e vectors
	e[0] = e[1] = e[2] = 0;
	e[3] = 1;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		e[4 * (l)+4] = 0;
		e[4 * (l)+5] = 0;
		e[4 * (l)+6] = 0;
		e[4 * (l)+7] = 0;
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the c vectors
	c[0] = c[1] = c[2] = 0;
	c[3] = Stempo;

	double averTem = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		c[4 * (l)+4] = 0;
		c[4 * (l)+5] = 0;
		c[4 * (l)+6] = 0;
		c[4 * (l)+7] = Mtmpo[l];

	}

	// calculate contact vectors for normal traction(Thermodynamics  case!)
	NTE[0] = nu.x;
	NTE[1] = nu.y;
	NTE[2] = nu.z;
	NTE[3] = 0;


	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l<nmeln; ++l)
	{
		//N[3 * (2 * l + 1)] = 0;
		//N[3 * (2 * l + 1) + 1] = 0;
		//N[3 * (2 * l + 1) + 2] = 0;
		NTE[4 * (l)+4] = -H[l] * nu.x;
		NTE[4 * (l)+5] = -H[l] * nu.y;
		NTE[4 * (l)+6] = -H[l] * nu.z;
		NTE[4 * (l)+7] = 0;

	}
	for (k = 0; k < ndof; ++k) { fe[k] = 0; }


	if (gap > 0){
		for (l = 0; l < ndof; ++l) { fe[l] = tn*NTE[l]; }




		//if ((m_mu*m_epsf > 0) && (gap <= max_pent))
		if ((m_mu*m_epsf > 0))
		{

			m_friction = true;
			// Lagrangian traction
			double Lt[2];
			Lt[0] = ss.m_Lt[m][0];
			Lt[1] = ss.m_Lt[m][1];

			// calculate contact vector for tangential traction
			// only if both the friction coefficient and friction
			// penalty factor are non-zero

			// get the master shape function derivative values at this slave node

			//////////////////////////  modified part begin
			tau1 = tau2 = vec3d(0, 0, 0);
			for (k = 0; k < nmeln; ++k)
			{
				tau1.x += Hr[k] * rtm[k].x;
				tau1.y += Hr[k] * rtm[k].y;
				tau1.z += Hr[k] * rtm[k].z;

				tau2.x += Hs[k] * rtm[k].x;
				tau2.y += Hs[k] * rtm[k].y;
				tau2.z += Hs[k] * rtm[k].z;
			}

			// set up the Ti vectors
			T1[0] = tau1.x; T2[0] = tau2.x;
			T1[1] = tau1.y; T2[1] = tau2.y;
			T1[2] = tau1.z; T2[2] = tau2.z;
			T1[3] = 0; T2[3] = 0;


			for (k = 0; k < nmeln; ++k)
			{
				T1[4 * (k)+4] = -H[k] * tau1.x;
				T1[4 * (k)+5] = -H[k] * tau1.y;
				T1[4 * (k)+6] = -H[k] * tau1.z;
				T1[4 * (k)+7] = 0;

				T2[4 * (k)+4] = -H[k] * tau2.x;
				T2[4 * (k)+5] = -H[k] * tau2.y;
				T2[4 * (k)+6] = -H[k] * tau2.z;
				T2[4 * (k)+7] = 0;
			}

			// set up the Ni vectors
			N1[0] = N2[0] = 0;
			N1[1] = N2[1] = 0;
			N1[2] = N2[2] = 0;
			N1[3] = N2[4] = 0;
			for (k = 0; k < nmeln; ++k)
			{
				N1[4 * (k)+4] = -Hr[k] * nu.x;
				N1[4 * (k)+5] = -Hr[k] * nu.y;
				N1[4 * (k)+6] = -Hr[k] * nu.z;
				N1[4 * (k)+7] = 0;

				N2[4 * (k)+4] = -Hs[k] * nu.x;
				N2[4 * (k)+5] = -Hs[k] * nu.y;
				N2[4 * (k)+6] = -Hs[k] * nu.z;
				N1[4 * (k)+7] = 0;
			}

			// calculate metric tensor
			M[0][0] = tau1*tau1; M[0][1] = tau1*tau2;
			M[1][0] = tau2*tau1; M[1][1] = tau2*tau2;

			// calculate curvature tensor
			K[0][0] = 0; K[0][1] = 0;
			K[1][0] = 0; K[1][1] = 0;

			double Grr[FEElement::MAX_NODES];
			double Grs[FEElement::MAX_NODES];
			double Gss[FEElement::MAX_NODES];
			mel.shape_deriv2(Grr, Grs, Gss, r, s);
			for (k = 0; k < nmeln; ++k)
			{
				K[0][0] += (nu*rtm[k])*Grr[k];
				K[0][1] += (nu*rtm[k])*Grs[k];
				K[1][0] += (nu*rtm[k])*Grs[k];
				K[1][1] += (nu*rtm[k])*Gss[k];
			}

			// setup A matrix
			A[0][0] = M[0][0] + gap*K[0][0];
			A[0][1] = M[0][1] + gap*K[0][1];
			A[1][0] = M[1][0] + gap*K[1][0];
			A[1][1] = M[1][1] + gap*K[1][1];

			detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

			// setup Di vectors
			for (k = 0; k < ndof; ++k)
			{
				if (detA != 0){
					D1[k] = (1 / detA)*(A[1][1] * (T1[k] + gap*N1[k]) - A[0][1] * (T2[k] + gap*N2[k]));
					D2[k] = (1 / detA)*(A[0][0] * (T2[k] + gap*N2[k]) - A[0][1] * (T1[k] + gap*N1[k]));
				}
				//tntt[k] = D1[k];  // exp(-(((D1[k])*(D1[k])) + ((D2[k])*(D2[k]))) / (eps*eps));
				//tttn[k] = D2[k];
			}
			double x1 = (r - rpf);
			double x2 = (s - spf);
			Tt[0] = Lt[0] + m_epsf*(Mk[0][0] * x1 + Mk[0][1] * x2);
			Tt[1] = Lt[1] + m_epsf*(Mk[1][0] * x1 + Mk[1][1] * x2);
			//if (isinf(Tt[0]) != 0){ Tt[0] = 0; }
			//if (isinf(Tt[1]) != 0){ Tt[1] = 0; }

			double TMT = Tt[0] * (Mki[0][0] * Tt[0] + Mki[0][1] * Tt[1]) + Tt[1] * (Mki[1][0] * Tt[0] + Mki[1][1] * Tt[1]);
			//assert(TMT >= 0);

			double phi = 0;
			if (TMT >= 0){ 
				phi = sqrt(TMT) - m_mu*Ln;
			}

			// b. return map
			if ((phi > 0) && (TMT > 0))
			{
				Tt[0] = m_mu*Ln*Tt[0] / sqrt(TMT);
				Tt[1] = m_mu*Ln*Tt[1] / sqrt(TMT);
			}
			// note that in the case of non perfect sliding,
			// the calculation of  Mechanical dissipation could be slightly different 
			// Mechanical dissipation


			double DM = 0;
			//double DM_α = fabs(Tt[0] * (r - rp));
			//double DM_β = fabs(Tt[1] * (s - sp));
			//if (isinf(x1) != 0){ x1 = 0; }
			//if (isinf(x2) != 0){ x2 = 0; }

			double DM_α = fabs(Tt[0] * x1);
			double DM_β = fabs(Tt[1] * x2);
			DM = DM_α + DM_β;
			if (∆time > 0) { DM = DM / ∆time; }
			// Temperature force vector
			for (l = 0; l < ndof; ++l){
				fe[l] += -R0 * tn*(a[l] * c[l])*a[l] + DM*b[l] + tn*((R[2] * (e[l] * c[l]) - R[3])*e[l] - (R[4] * (f[l] * c[l]) - R[5])*f[l]);

			}

			///Temperature part  end   
			// tangential force vector
			for (l = 0; l < ndof; ++l) fe[l] -= (Tt[0] * D1[l] + Tt[1] * D2[l]);

		}

	}
////////////////////////////////////////////////////////////////////////////////////////////////////////
	if ((m_mu*m_epsf > 0) && (gap <= 0) && (fabs(gap) <= (fric_tol*δ_TN)))
		
	{

		m_friction = true;
		// Lagrangian traction

		/////////////////////////////////
		double max_normal = M_TN* ss.D_R[m];
		double max_tangential = M_Tt* ss.D_R[m];

		double mCL_normal = δ_TN *ss.D_R_C[m];
		double mCL_tangential = δ_Tt * ss.D_R_C[m];
		//eps = mCL_normal = mCL_tangential;
		double gap_n = fabs(gap);

		double Utn_α = -((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
		double Utn_β = -((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
		//double Utt_α = (exp(-(((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
		//double Utt_β = (exp(-(((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));


		double t_N_α = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_α);
		double t_N_β = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_β);
		double tracn = t_N_α + t_N_β;


		double LNX = -(MBRACKET(-tracn));

		/////////////////////////////////////

		double Lt[2];
		Lt[0] = ss.m_trans_S[m][0];
		Lt[1] = ss.m_trans_S[m][1];

		// calculate contact vector for tangential traction
		// only if both the friction coefficient and friction
		// penalty factor are non-zero

		// get the master shape function derivative values at this slave node
		//mel.shape_deriv(Hr, Hs, r, s);
		//////////////////////////  modified part begin
		tau1 = tau2 = vec3d(0, 0, 0);
		for (k = 0; k < nmeln; ++k)
		{
			tau1.x += Hr[k] * rtm[k].x;
			tau1.y += Hr[k] * rtm[k].y;
			tau1.z += Hr[k] * rtm[k].z;

			tau2.x += Hs[k] * rtm[k].x;
			tau2.y += Hs[k] * rtm[k].y;
			tau2.z += Hs[k] * rtm[k].z;
		}

		// set up the Ti vectors
		T1[0] = tau1.x; T2[0] = tau2.x;
		T1[1] = tau1.y; T2[1] = tau2.y;
		T1[2] = tau1.z; T2[2] = tau2.z;
		T1[3] = 0; T2[3] = 0;


		for (k = 0; k < nmeln; ++k)
		{
			T1[4 * (k)+4] = -H[k] * tau1.x;
			T1[4 * (k)+5] = -H[k] * tau1.y;
			T1[4 * (k)+6] = -H[k] * tau1.z;
			T1[4 * (k)+7] = 0;

			T2[4 * (k)+4] = -H[k] * tau2.x;
			T2[4 * (k)+5] = -H[k] * tau2.y;
			T2[4 * (k)+6] = -H[k] * tau2.z;
			T2[4 * (k)+7] = 0;
		}

		// set up the Ni vectors
		N1[0] = N2[0] = 0;
		N1[1] = N2[1] = 0;
		N1[2] = N2[2] = 0;
		N1[3] = N2[4] = 0;
		for (k = 0; k < nmeln; ++k)
		{
			N1[4 * (k)+4] = -Hr[k] * nu.x;
			N1[4 * (k)+5] = -Hr[k] * nu.y;
			N1[4 * (k)+6] = -Hr[k] * nu.z;
			N1[4 * (k)+7] = 0;

			N2[4 * (k)+4] = -Hs[k] * nu.x;
			N2[4 * (k)+5] = -Hs[k] * nu.y;
			N2[4 * (k)+6] = -Hs[k] * nu.z;
			N1[4 * (k)+7] = 0;
		}

		// calculate metric tensor
		M[0][0] = tau1*tau1; M[0][1] = tau1*tau2;
		M[1][0] = tau2*tau1; M[1][1] = tau2*tau2;

		// calculate curvature tensor
		K[0][0] = 0; K[0][1] = 0;
		K[1][0] = 0; K[1][1] = 0;

		double Grr[FEElement::MAX_NODES];
		double Grs[FEElement::MAX_NODES];
		double Gss[FEElement::MAX_NODES];
		mel.shape_deriv2(Grr, Grs, Gss, r, s);
		for (k = 0; k < nmeln; ++k)
		{
			K[0][0] += (nu*rtm[k])*Grr[k];
			K[0][1] += (nu*rtm[k])*Grs[k];
			K[1][0] += (nu*rtm[k])*Grs[k];
			K[1][1] += (nu*rtm[k])*Gss[k];
		}

		// setup A matrix
		A[0][0] = M[0][0] + gap*K[0][0];
		A[0][1] = M[0][1] + gap*K[0][1];
		A[1][0] = M[1][0] + gap*K[1][0];
		A[1][1] = M[1][1] + gap*K[1][1];

		detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

		// setup Di vectors
		for (k = 0; k < ndof; ++k)
		{
			if (detA != 0){
				D1[k] = (1 / detA)*(A[1][1] * (T1[k] + gap*N1[k]) - A[0][1] * (T2[k] + gap*N2[k]));
				D2[k] = (1 / detA)*(A[0][0] * (T2[k] + gap*N2[k]) - A[0][1] * (T1[k] + gap*N1[k]));
			}
			//tntt[k] = D1[k];  // exp(-(((D1[k])*(D1[k])) + ((D2[k])*(D2[k]))) / (eps*eps));
			//tttn[k] = D2[k];
		}


		// calculate friction tractions
		// a. calculate trial state
		double x1 = (r - rp);
		double x2 = (s - sp);

		Tt[0] = Lt[0] + m_epsf*(Mk[0][0] * x1 + Mk[0][1] * x2);
		Tt[1] = Lt[1] + m_epsf*(Mk[1][0] * x1 + Mk[1][1] * x2);

		//if (isinf(Tt[0]) != 0){ Tt[0] = 0; }
		//if (isinf(Tt[1]) != 0){ Tt[1] = 0; }

		double TMT = Tt[0] * (Mki[0][0] * Tt[0] + Mki[0][1] * Tt[1]) + Tt[1] * (Mki[1][0] * Tt[0] + Mki[1][1] * Tt[1]);
		//assert(TMT >= 0);


		// b. return map

		LNX = fabs(LNX);
		//double phi = sqrt(TMT) - m_mu*Ln;
		double phi = 0;
		
		if (TMT >= 0){
			phi = sqrt(TMT) - m_mu*LNX;
		}
		// b. return map
		if ((phi > 0) && (TMT > 0))
		{
			//Tt[0] = m_mu*Ln*Tt[0] / sqrt(TMT);
			//Tt[1] = m_mu*Ln*Tt[1] / sqrt(TMT);
			Tt[0] = m_mu*LNX*Tt[0] / sqrt(TMT);
			Tt[1] = m_mu*LNX*Tt[1] / sqrt(TMT);
		}

		// tangential force vector
		for (l = 0; l < ndof; ++l) fe[l] -= (Tt[0] * D1[l] + Tt[1] * D2[l]);


		///Temperature part begin
		// calculate Temperature force vector  and Mechanical dissipation  if 
		// only if both the friction coefficient and friction
		// penalty factor are non-zero
		// note that in the case of non perfect sliding!,(|Tt[0] * D1[l]| + |Tt[1] * D2[l]|)/∆time
		// the calculatation of  Mechanical dissipation could be slightly different 
		// Mechanical dissipation
		double DM = 0;
		//double DM_α = fabs(Tt[0] * (r - rp));
		//double DM_β = fabs(Tt[1] * (s - sp));
		//if (isinf(x1) != 0){ x1 = 0; }
		//if (isinf(x2) != 0){ x2 = 0; }

		double DM_α = fabs(Tt[0] * x1);
		double DM_β = fabs(Tt[1] * x2);
		DM = DM_α + DM_β;
		if (∆time > 0) { DM = DM / ∆time; }
		// Temperature force vector
		for (l = 0; l < ndof; ++l){ fe[l] += -R0*LNX*(a[l] * c[l])*a[l] + DM*b[l] + LNX*((R[2] * (e[l] * c[l]) - R[3])*e[l] - (R[4] * (f[l] * c[l]) - R[5])*f[l]); }

		// tangential force vector
		for (l = 0; l < ndof; ++l) fe[l] -= (Tt[0] * D1[l] + Tt[1] * D2[l]);


	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
if (gap <= 0){
	////////////////////
	
	
	
	double max_normal = M_TN* ss.D_R[m];
	double max_tangential = M_Tt* ss.D_R[m];

	double mCL_normal = δ_TN * ss.D_R_C[m];
	double mCL_tangential = δ_Tt * ss.D_R_C[m];
	//eps = mCL_normal = mCL_tangential;
	//mel.shape_deriv(Hr, Hs, r, s);
	//////////////////////////  modified part begin
	tau1 = tau2 = vec3d(0, 0, 0);
	for ( k = 0 ; k < nmeln; ++k)
	{
		tau1.x += Hr[k] * rtm[k].x;
		tau1.y += Hr[k] * rtm[k].y;
		tau1.z += Hr[k] * rtm[k].z;

		tau2.x += Hs[k] * rtm[k].x;
		tau2.y += Hs[k] * rtm[k].y;
		tau2.z += Hs[k] * rtm[k].z;
	}

	// set up the Ti vectors
	T1[0] = tau1.x; T2[0] = tau2.x;
	T1[1] = tau1.y; T2[1] = tau2.y;
	T1[2] = tau1.z; T2[2] = tau2.z;
	T1[3] = 0; T2[3] = 0;


	for ( k = 0 ; k < nmeln; ++k)
	{
		T1[4 * (k)+4] = -H[k] * tau1.x;
		T1[4 * (k)+5] = -H[k] * tau1.y;
		T1[4 * (k)+6] = -H[k] * tau1.z;
		T1[4 * (k)+7] = 0;

		T2[4 * (k)+4] = -H[k] * tau2.x;
		T2[4 * (k)+5] = -H[k] * tau2.y;
		T2[4 * (k)+6] = -H[k] * tau2.z;
		T2[4 * (k)+7] = 0;
	}

	// set up the Ni vectors
	N1[0] = N2[0] = 0;
	N1[1] = N2[1] = 0;
	N1[2] = N2[2] = 0;
	N1[3] = N2[4] = 0;
	for ( k = 0; k < nmeln; ++k)
	{
		N1[4 * (k)+4] = -Hr[k] * nu.x;
		N1[4 * (k)+5] = -Hr[k] * nu.y;
		N1[4 * (k)+6] = -Hr[k] * nu.z;
		N1[4 * (k)+7] = 0;

		N2[4 * (k)+4] = -Hs[k] * nu.x;
		N2[4 * (k)+5] = -Hs[k] * nu.y;
		N2[4 * (k)+6] = -Hs[k] * nu.z;
		N1[4 * (k)+7] = 0;
	}

	// calculate metric tensor
	M[0][0] = tau1*tau1; M[0][1] = tau1*tau2;
	M[1][0] = tau2*tau1; M[1][1] = tau2*tau2;

	// calculate curvature tensor
	K[0][0] = 0; K[0][1] = 0;
	K[1][0] = 0; K[1][1] = 0;

	double Grr[FEElement::MAX_NODES];
	double Grs[FEElement::MAX_NODES];
	double Gss[FEElement::MAX_NODES];
	mel.shape_deriv2(Grr, Grs, Gss, r, s);
	for ( k = 0; k < nmeln; ++k)
	{
		K[0][0] += (nu*rtm[k])*Grr[k];
		K[0][1] += (nu*rtm[k])*Grs[k];
		K[1][0] += (nu*rtm[k])*Grs[k];
		K[1][1] += (nu*rtm[k])*Gss[k];
	}

	// setup A matrix
	A[0][0] = M[0][0] + gap*K[0][0];
	A[0][1] = M[0][1] + gap*K[0][1];
	A[1][0] = M[1][0] + gap*K[1][0];
	A[1][1] = M[1][1] + gap*K[1][1];

	detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	// setup Di vectors
	for ( k = 0; k < ndof; ++k)
	{
		if (detA != 0){
			D1[k] = (1 / detA)*(A[1][1] * (T1[k] + gap*N1[k]) - A[0][1] * (T2[k] + gap*N2[k]));
			D2[k] = (1 / detA)*(A[0][0] * (T2[k] + gap*N2[k]) - A[0][1] * (T1[k] + gap*N1[k]));
		}

	}
	//eps = mCL_normal = mCL_tangential;
	double gap_n = fabs(gap);


	double Utn_α = -((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
	double Utn_β = -((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
	double Utt_α = (exp(-(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
	double Utt_β = (exp(-(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));

	
	double t_N_α = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_α);
	double t_N_β = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_β);
	double tracn = t_N_α + t_N_β;
	// contact traction( cohesive case! mode II and III)
	double t_T_α = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_α;
	double t_T_β = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_β;
	//double tract = t_T_α + t_T_β;
	double tn_x = -(MBRACKET(-tracn));
	double tt_α = -(MBRACKET(-t_T_α));
	double tt_β = -(MBRACKET(-t_T_β));


	for ( k = 0; k < ndof; ++k) {

		fe[k] += (tn_x)*NTE[k];
		fe[k] += ((tt_α*D1[k]) + (tt_β*D2[k]));

		}
	


	}
	
	

}


//-----------------------------------------------------------------------------

//void FEFacet2FacetSliding::ContactStiffness(FESolver* psolver)
void FEFacet2FacetSliding::ContactStiffness(FESolver* psolver)
{
	

	int j, k, l, n, m, np;
	int nseln = 0, nmeln = 0, ndof = 0;

	matrix ke;

	const int MAXMN = FEElement::MAX_NODES;
	//vector<int> lm(3*(MAXMN + 1));
	//vector<int> en(MAXMN+1);
	//const int MN = FEElement::MAX_NODES;
	//double *Gr, *Gs, w[6];
	//vec3d r0[6];
	//double detJ[6];
	//vec3d dxr, dxs;
	vec3d r0[MAXMN];
	double w[MAXMN];
	double* Gr, *Gs;
	double detJ[MAXMN];
	vec3d dxr, dxs;

	vector<int> lm;
	vector<int> en;
	vector<int> sLM;
	vector<int> mLM;

	// do two-pass
	int npass = (m_btwo_pass ? 2 : 1);
	for (np = 0; np<npass; ++np)
	{
		// get the master and slave surface
		FEFacetSlidingSurface& ss = (np == 0 ? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0 ? m_ms : m_ss);

		// loop over all slave elements
		int ne = ss.Elements();
		for (j = 0; j<ne; ++j)
		{
			// unpack the slave element
			FESurfaceElement& se = ss.Element(j);
			nseln = se.Nodes();

			// get the element's LM array
			ss.UnpackLM(se, sLM);

			// get the nodal coordinates
			for (int i = 0; i<nseln; ++i) r0[i] = ss.GetMesh()->Node(se.m_node[i]).m_r0;

			// get all the metrics we need 
			for (n = 0; n<nseln; ++n)
			{
				Gr = se.Gr(n);
				Gs = se.Gs(n);

				// calculate jacobian
				dxr = dxs = vec3d(0, 0, 0);
				for (k = 0; k<nseln; ++k)
				{
					dxr.x += Gr[k] * r0[k].x;
					dxr.y += Gr[k] * r0[k].y;
					dxr.z += Gr[k] * r0[k].z;

					dxs.x += Gs[k] * r0[k].x;
					dxs.y += Gs[k] * r0[k].y;
					dxs.z += Gs[k] * r0[k].z;
				}

				detJ[n] = (dxr ^ dxs).norm();
				w[n] = se.GaussWeights()[n];
			}

			// loop over all integration points (that is nodes)
			for (n = 0; n<nseln; ++n)
			{
				m = se.m_lnode[n];

				// see if this node's constraint is active
				// that is, if it has a master element associated with it
				if (ss.m_pme[m] != 0)
				{
					// get the master element
					FESurfaceElement& me = *ss.m_pme[m];

					// get the masters element's LM array
					ms.UnpackLM(me, mLM);

					nmeln = me.Nodes();
					// 
					//ndof = 3*(nmeln+1);
					// add Temperature dof
					ndof = 4 * (nmeln + 1);
					// calculate the stiffness matrix
					//ke.zero();
					ke.resize(ndof, ndof);
					ContactNodalStiffness(m, ss, me, ke);

					// muliply with weights
					for (k = 0; k<ndof; ++k)
						for (l = 0; l<ndof; ++l) ke[k][l] *= detJ[n] * w[n];


					// test if sequence is empty
					if (!lm.empty() == true){ lm.clear(); }

					// add Temperature dof
					//(mortar) 
					lm.resize(4 * (nmeln + 1));
					lm[0] = sLM[n * 3];
					lm[1] = sLM[n * 3 + 1];
					lm[2] = sLM[n * 3 + 2];
					lm[3] = sLM[10 * nseln + n];
					for (l = 0; l < nmeln; ++l)
					{
						lm[4 * (l)+4] = mLM[l * 3];
						lm[4 * (l)+5] = mLM[l * 3 + 1];
						lm[4 * (l)+6] = mLM[l * 3 + 2];
						lm[4 * (l)+7] = mLM[10 * nmeln + l];

					}

					// fill the lm array
					//lm[0] = sLM[n*3  ];
					//lm[1] = sLM[n*3+1];
					//lm[2] = sLM[n*3+2];

					//for (k=0; k<nmeln; ++k)
					//{
					//	lm[3*(k+1)  ] = mLM[k*3  ];
					//	lm[3*(k+1)+1] = mLM[k*3+1];
					//lm[3*(k+1)+2] = mLM[k*3+2];
					//}

					// create the en array
					//(mortar)
					// test if sequence is empty
					if (!en.empty() == true){ en.clear(); }
					en.resize(nmeln + 1);
					en[0] = se.m_node[n];
					for (k = 0; k<nmeln; ++k) en[k + 1] = me.m_node[k];

					// assemble stiffness matrix

					psolver->AssembleStiffness(en, lm, ke);

				}
			}
		}
	}
}

//-----------------------------------------------------------------------------

void FEFacet2FacetSliding::ContactNodalStiffness(int m, FEFacetSlidingSurface& ss, FESurfaceElement& mel, matrix& ke)
{

	int i, j, k, l;

	const int MAXMN = FEElement::MAX_NODES;
	FEAnalysis& step = *(GetFEModel()->GetCurrentStep());
	//double dt = pstep->m_dt;
	double ∆time = step.m_dt;
	double F1 = 0;
	if (∆time > 0 && T0_Contact){ F1 = (C_heat*TCN_Contact) / (T0_Contact*∆time); }
	double F2 = F1*TCN_Contact;
	double R0 = (γ1_transfer*γ2_transfer) / (γ1_transfer + γ2_transfer + F1);
	R[0] = R0 / γ2_transfer;
	R[1] = R0 / γ1_transfer;
	R[2] = F1*R[0];
	R[3] = F2*R[0];
	R[4] = F1*R[1];
	R[5] = F2*R[1];
	// get the slave element nodal temperaturesss.GetMesh()->Node(se.m_node[i])m_MSAT
	double Stempo = ss.Node(m).m_T;
	double Mtmpo[MAXMN];

	// temperatures vectors
	//double a[4 * (MAXMN + 1)], b[4 * (MAXMN + 1)], c[4 * (MAXMN + 1)], d1[4 * (MAXMN + 1)], d2[4 * (MAXMN + 1)], e[4 * (MAXMN + 1)], f[4 * (MAXMN + 1)];
	//double NTE[4 * (MAXMN + 1)];


	std::vector<double>	a(4 * (MAXMN + 1));
	std::vector<double>	b(4 * (MAXMN + 1));
	std::vector<double>	c(4 * (MAXMN + 1));
	std::vector<double>	d1(4 * (MAXMN + 1));
	std::vector<double>	d2(4 * (MAXMN + 1));
	std::vector<double>	e(4 * (MAXMN + 1));
	std::vector<double>	f(4 * (MAXMN + 1));
	std::vector<double>	NTE(4 * (MAXMN + 1));

	
	//vector<int> lm(3*(MAXMN+1));
	//vector<int> en(MAXMN + 1);


	//vector<int> lm;
	//vector<int> en;

	vec3d dxr, dxs;
	double H[MAXMN], Hr[MAXMN], Hs[MAXMN];

	//double N[3*(MAXMN+1)], T1[3*(MAXMN+1)], T2[3*(MAXMN+1)];
	//double N1[3*(MAXMN+1)], N2[3*(MAXMN+1)], D1[3*(MAXMN+1)], D2[3*(MAXMN+1)];
	//double Nb1[3*(MAXMN+1)], Nb2[3*(MAXMN+1)];

	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// nr of element nodes and degrees of freedom 
	int nmeln = mel.Nodes();
	//int ndof = 3*(1 + nmeln);

	// penalty factor
	double scale = Penalty();
	double eps = ss.m_eps[m] * scale;

	// nodal coordinates
	vec3d rt[MAXMN];
	for (j = 0; j<nmeln; ++j) rt[j] = mesh.Node(mel.m_node[j]).m_rt;
	// get the master elements nodal temperatures
	for (k = 0; k<nmeln; ++k) Mtmpo[k] = mesh.Node(mel.m_node[k]).m_T;

	// slave node natural coordinates in master element
	double r = ss.m_rs[m][0];
	double s = ss.m_rs[m][1];

	// slave gap
	double gap = ss.m_gap[m];

	// lagrange multiplier
	double Lm = ss.m_Lm[m];

	// get slave node normal force
	double tn = Lm + eps*gap;
	tn = MBRACKET(tn);

	// get the slave node normal
	vec3d& nu = ss.m_nu[m];

	// get the master shape function values and the derivatives at this slave node
	mel.shape_fnc(H, r, s);
	mel.shape_deriv(Hr, Hs, r, s);

	// get the tangent vectors
	vec3d tau[2];
	ss.CoBaseVectors(mel, r, s, tau);


	/////////////////   modified part begin (mortar Thermodynamics  case!)

	// temperatures vectors


	//double N1[4 * (MAXMN + 1)], N2[4 * (MAXMN + 1)];
	//double T1[4 * (MAXMN + 1)], T2[4 * (MAXMN + 1)], D1[4 * (MAXMN + 1)], D2[4 * (MAXMN + 1)];
	//double Nb1[4 * (MAXMN + 1)], Nb2[4 * (MAXMN + 1)];
	//vector<int> lm(4 * (MAXMN + 1));

	//double T11[4 * (MAXMN + 1)], T12[4 * (MAXMN + 1)], T21[4 * (MAXMN + 1)], T22[4 * (MAXMN + 1)];	// Tab matrix
	//double N11[4 * (MAXMN + 1)], N12[4 * (MAXMN + 1)], N21[4 * (MAXMN + 1)], N22[4 * (MAXMN + 1)];	// Nab matrix
	//double P1[4 * (MAXMN + 1)], P2[4 * (MAXMN + 1)];	// P arrays
	//double Tb11[4 * (MAXMN + 1)], Tb21[4 * (MAXMN + 1)], Tb12[4 * (MAXMN + 1)], Tb22[4 * (MAXMN + 1)]; // Tbar matrix
	//double Pb1[4 * (MAXMN + 1)], Pb2[4 * (MAXMN + 1)]; // Pbar array

	std::vector<double>	N1(4 * (MAXMN + 1));
	std::vector<double>	N2(4 * (MAXMN + 1));
	std::vector<double>	T1(4 * (MAXMN + 1));
	std::vector<double>	T2(4 * (MAXMN + 1));
	std::vector<double>	D1(4 * (MAXMN + 1));
	std::vector<double>	D2(4 * (MAXMN + 1));


	std::vector<double>	Nb1(4 * (MAXMN + 1));
	std::vector<double>	Nb2(4 * (MAXMN + 1));
	std::vector<double>	T11(4 * (MAXMN + 1));
	std::vector<double>	T12(4 * (MAXMN + 1));
	std::vector<double>	T21(4 * (MAXMN + 1));
	std::vector<double>	T22(4 * (MAXMN + 1));


	std::vector<double>	N11(4 * (MAXMN + 1));
	std::vector<double>	N12(4 * (MAXMN + 1));
	std::vector<double>	N21(4 * (MAXMN + 1));
	std::vector<double>	N22(4 * (MAXMN + 1));
	std::vector<double>	P1(4 * (MAXMN + 1));
	std::vector<double>	P2(4 * (MAXMN + 1));


	std::vector<double>	Tb11(4 * (MAXMN + 1));
	std::vector<double>	Tb21(4 * (MAXMN + 1));
	std::vector<double>	Tb12(4 * (MAXMN + 1));
	std::vector<double>	Tb22(4 * (MAXMN + 1));
	std::vector<double>	Pb1(4 * (MAXMN + 1));
	std::vector<double>	Pb2(4 * (MAXMN + 1));






	// add Temperature dof
	int ndof = 4 * (1 + nmeln);

	// add Temperature dof
	//	double  N1[4 * (MAXMN + 1)] = { 0 }, N2[4 * (MAXMN + 1)] = { 0 };
	//double T1[4 * (MAXMN + 1)] = { 0 }, T2[4 * (MAXMN + 1)] = { 0 }, D1[4 * (MAXMN + 1)] = { 0 }, D2[4 * (MAXMN + 1)] = { 0 };

	// set up the a vectors
	a[0] = a[1] = a[2] = 0;
	a[3] = 1;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		a[4 * (l)+4] = 0;
		a[4 * (l)+5] = 0;
		a[4 * (l)+6] = 0;
		a[4 * (l)+7] = -H[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	b[0] = b[1] = b[2] = 0;
	b[3] = R[0];
	for (l = 0; l < nmeln; ++l)
	{

		b[4 * (l)+4] = 0;
		b[4 * (l)+5] = 0;
		b[4 * (l)+6] = 0;
		b[4 * (l)+7] = H[l] * R[1];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the d1 vectors
	d1[0] = d1[1] = d1[2] = 0;
	d1[3] = 0;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		d1[4 * (l)+4] = 0;
		d1[4 * (l)+5] = 0;
		d1[4 * (l)+6] = 0;
		d1[4 * (l)+7] = -Hr[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the d2 vectors
	d2[0] = d2[1] = d2[2] = 0;
	d2[3] = 0;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		d2[4 * (l)+4] = 0;
		d2[4 * (l)+5] = 0;
		d2[4 * (l)+6] = 0;
		d2[4 * (l)+7] = -Hs[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the f vectors
	f[0] = f[1] = f[2] = 0;
	f[3] = 0;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		f[4 * (l)+4] = 0;
		f[4 * (l)+5] = 0;
		f[4 * (l)+6] = 0;
		f[4 * (l)+7] = H[l];
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the e vectors
	e[0] = e[1] = e[2] = 0;
	e[3] = 1;
	//a[2] = a[5] = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		e[4 * (l)+4] = 0;
		e[4 * (l)+5] = 0;
		e[4 * (l)+6] = 0;
		e[4 * (l)+7] = 0;
		//a[3 * (2 * l + 1) + 3] = -H[l] ;
		//a[3 * (2 * l + 1) + 4] = -H[l] ;
		//a[3 * (2 * l + 1) + 5] = -H[l] ;

	}
	// set up the c vectors
	c[0] = c[1] = c[2] = 0;
	c[3] = Stempo;

	double averTem = 0;
	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l < nmeln; ++l)
	{

		c[4 * (l)+4] = 0;
		c[4 * (l)+5] = 0;
		c[4 * (l)+6] = 0;
		c[4 * (l)+7] = Mtmpo[l];

	}

	// calculate contact vectors for normal traction(Thermodynamics  case!)
	NTE[0] = nu.x;
	NTE[1] = nu.y;
	NTE[2] = nu.z;
	NTE[3] = 0;


	//  temperature degrees of freedom + Mechanic degrees of freedom
	for (l = 0; l<nmeln; ++l)
	{
		//N[3 * (2 * l + 1)] = 0;
		//N[3 * (2 * l + 1) + 1] = 0;
		//N[3 * (2 * l + 1) + 2] = 0;
		NTE[4 * (l)+4] = -H[l] * nu.x;
		NTE[4 * (l)+5] = -H[l] * nu.y;
		NTE[4 * (l)+6] = -H[l] * nu.z;
		NTE[4 * (l)+7] = 0;

	}

	//	ndof = 4 * (1 + nmeln);
	// set up the Ti vectors
	T1[0] = tau[0].x; T2[0] = tau[1].x;
	T1[1] = tau[0].y; T2[1] = tau[1].y;
	T1[2] = tau[0].z; T2[2] = tau[1].z;
	T1[3] = 0; T2[3] = 0;


	for (k = 0; k < nmeln; ++k)
	{
		T1[4 * (k)+4] = -H[k] * tau[0].x;
		T1[4 * (k)+5] = -H[k] * tau[0].y;
		T1[4 * (k)+6] = -H[k] * tau[0].z;
		T1[4 * (k)+7] = 0;

		T2[4 * (k)+4] = -H[k] * tau[1].x;
		T2[4 * (k)+5] = -H[k] * tau[1].y;
		T2[4 * (k)+6] = -H[k] * tau[1].z;
		T2[4 * (k)+7] = 0;
	}

	// set up the Ni vectors
	N1[0] = N2[0] = 0;
	N1[1] = N2[1] = 0;
	N1[2] = N2[2] = 0;
	N1[3] = N2[4] = 0;
	for (k = 0; k < nmeln; ++k)
	{
		N1[4 * (k)+4] = -Hr[k] * nu.x;
		N1[4 * (k)+5] = -Hr[k] * nu.y;
		N1[4 * (k)+6] = -Hr[k] * nu.z;
		N1[4 * (k)+7] = 0;

		N2[4 * (k)+4] = -Hs[k] * nu.x;
		N2[4 * (k)+5] = -Hs[k] * nu.y;
		N2[4 * (k)+6] = -Hs[k] * nu.z;
		N1[4 * (k)+7] = 0;
	}

	// calculate metric tensor
	mat2d M;
	M[0][0] = tau[0] * tau[0]; M[0][1] = tau[0] * tau[1];
	M[1][0] = tau[1] * tau[0]; M[1][1] = tau[1] * tau[1];

	// calculate reciprocal metric tensor
	mat2d Mi = M.inverse();

	// calculate curvature tensor
	double K[2][2] = { 0 };
	double Grr[FEElement::MAX_NODES];
	double Grs[FEElement::MAX_NODES];
	double Gss[FEElement::MAX_NODES];
	mel.shape_deriv2(Grr, Grs, Gss, r, s);
	for (k = 0; k<nmeln; ++k)
	{
		K[0][0] += (nu*rt[k])*Grr[k];
		K[0][1] += (nu*rt[k])*Grs[k];
		K[1][0] += (nu*rt[k])*Grs[k];
		K[1][1] += (nu*rt[k])*Gss[k];
	}

	// setup A matrix A = M + gK
	double A[2][2];
	A[0][0] = M[0][0] + gap*K[0][0];
	A[0][1] = M[0][1] + gap*K[0][1];
	A[1][0] = M[1][0] + gap*K[1][0];
	A[1][1] = M[1][1] + gap*K[1][1];

	// calculate determinant of A
	double detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	// setup Di vectors
	for (k = 0; k<ndof; ++k)
	{
		if (detA != 0){
			D1[k] = (1 / detA)*(A[1][1] * (T1[k] + gap*N1[k]) - A[0][1] * (T2[k] + gap*N2[k]));
			D2[k] = (1 / detA)*(A[0][0] * (T2[k] + gap*N2[k]) - A[0][1] * (T1[k] + gap*N1[k]));
		}
		//tntt[k] = D1[k];  // exp(-(((D1[k])*(D1[k])) + ((D2[k])*(D2[k]))) / (eps*eps));
		//tttn[k] = D2[k];
	}

	// setup Nbi vectors
	for (k = 0; k<ndof; ++k)
	{
		Nb1[k] = N1[k] - K[0][1] * D2[k];
		Nb2[k] = N2[k] - K[0][1] * D1[k];
	}

	for (k = 0; k < ndof; ++k){
		for (l = 0; l < ndof; ++l){
			ke[k][l] =  0 ;

		}
	}

	// --- N O R M A L   S T I F F N E S S ---
	double sum=0;

	if (gap > 0){
		for (k = 0; k < ndof; ++k){
			for (l = 0; l < ndof; ++l)
			{

				sum = Mi[0][0] * Nb1[k] * Nb1[l] + Mi[0][1] * (Nb1[k] * Nb2[l] + Nb2[k] * Nb1[l]) + Mi[1][1] * Nb2[k] * Nb2[l];
				sum *= gap;
				sum -= D1[k] * N1[l] + D2[k] * N2[l] + N1[k] * D1[l] + N2[k] * D2[l];
				sum += K[0][1] * (D1[k] * D2[l] + D2[k] * D1[l]);
				sum *= tn*m_knmult;

				sum += eps*HEAVYSIDE(Lm + eps*gap)*NTE[k] * NTE[l];

				ke[k][l] = sum;
			}
		}
		// --- T A N G E N T I A L   S T I F F N E S S ---
		// We only calculate the tangential stiffness if friction is enabled. We also
		// make sure that the gap >= 0, i.e. the node is actually in contact, otherwise
		// I've noticed that the solution can diverge quickly.
			//if ((m_mu*m_epsf > 0) && (gap <= max_pent))
			if ((m_mu*m_epsf > 0))
		{
			// get the traction multipliers

			///////////////////   modified part begin

			double Lt[2];
			Lt[0] = ss.m_Lt[m][0];
			Lt[1] = ss.m_Lt[m][1];

			// get the metric tensor and its inverse
			mat2d& Mk = ss.m_M[m];
			mat2d Mki = Mk.inverse();

			// get the previous isoparameteric coordinates
			//double rp = ss.m_rsp[m][0];
			//double sp = ss.m_rsp[m][1];
			double rpf = ss.m_rsp_f[m][0];
			double spf = ss.m_rsp_f[m][1];
			/////////////////////


			// get the traction
			double Tt[2];
			// a. trial state
			double x1 = (r - rpf);
			double x2 = (s - spf);

			Tt[0] = Lt[0] + m_epsf*(Mk[0][0] * x1 + Mk[0][1] * x2);
			Tt[1] = Lt[1] + m_epsf*(Mk[1][0] * x1 + Mk[1][1] * x2);
			//Tt[0] = Lt[0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp);
			//Tt[1] = Lt[1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp);

			//if (isinf(Tt[0]) != 0){ Tt[0] = 0; }
			//if (isinf(Tt[1]) != 0){ Tt[1] = 0; }

			double TMT = Tt[0] * (Mki[0][0] * Tt[0] + Mki[0][1] * Tt[1]) + Tt[1] * (Mki[1][0] * Tt[0] + Mki[1][1] * Tt[1]);

			// calculate the normalized traction vector



			/////////////////////
			// calculate the covariant version
			double Pt[2] = { 0, 0 };

			if (TMT > 0){ 
				Pt[0] = Tt[0] / sqrt(TMT);
				Pt[1] = Tt[1] / sqrt(TMT);
			}
			double Ptc[2];
			Ptc[0] = Mki[0][0] * Pt[0] + Mki[0][1] * Pt[1];
			Ptc[1] = Mki[1][0] * Pt[0] + Mki[1][1] * Pt[1];

			//b. return map
			bool bstick = true;
			double phi = 0;
			
			if (TMT >= 0){ 
				phi = sqrt(TMT) - m_mu*tn;
			}
				


			if ((phi > 0) && (TMT > 0))
			{
				Tt[0] = m_mu*Tt[0] / sqrt(TMT);
				Tt[1] = m_mu*Tt[1] / sqrt(TMT);
				bstick = false;
			}

			// get the previous isoparameteric coordinates
			//	double rp = ss.m_rsp[m][0];
			//	double sp = ss.m_rsp[m][1];
			vec3d pt = tau[0] * Tt[0] + tau[1] * Tt[1];
			pt.unit();
			for (k = 0; k < nmeln; ++k)
			{
				T11[4 * (k)+4] = -Hr[k] * tau[0].x;
				T11[4 * (k)+5] = -Hr[k] * tau[0].y;
				T11[4 * (k)+6] = -Hr[k] * tau[0].z;
				T11[4 * (k)+7] = 0;

				T12[4 * (k)+4] = -Hs[k] * tau[0].x;
				T12[4 * (k)+5] = -Hs[k] * tau[0].y;
				T12[4 * (k)+6] = -Hs[k] * tau[0].z;
				T12[4 * (k)+7] = 0;

				T21[4 * (k)+4] = -Hr[k] * tau[1].x;
				T21[4 * (k)+5] = -Hr[k] * tau[1].y;
				T21[4 * (k)+6] = -Hr[k] * tau[1].z;
				T21[4 * (k)+7] = 0;

				T22[4 * (k)+4] = -Hs[k] * tau[1].x;
				T22[4 * (k)+5] = -Hs[k] * tau[1].y;
				T22[4 * (k)+6] = -Hs[k] * tau[1].z;
				T22[4 * (k)+7] = 0;



				P1[4 * (k)+4] = -Hr[k] * pt.x;
				P1[4 * (k)+5] = -Hr[k] * pt.y;
				P1[4 * (k)+6] = -Hr[k] * pt.z;
				P1[4 * (k)+7] = 0;

				P2[4 * (k)+4] = -Hs[k] * pt.x;
				P2[4 * (k)+5] = -Hs[k] * pt.y;
				P2[4 * (k)+6] = -Hs[k] * pt.z;
				P2[4 * (k)+7] = 0;
			}
			if (nmeln == 4)
			{
				N12[4] = N21[4] = -0.25*nu.x;
				N12[5] = N21[5] = -0.25*nu.y;
				N12[6] = N21[6] = -0.25*nu.z;
				N12[7] = N21[7] = 0;
				N12[8] = N21[8] = 0.25*nu.x;
				N12[9] = N21[9] = 0.25*nu.y;
				N12[10] = N21[10] = 0.25*nu.z;
				N12[11] = N21[11] = 0;
				N12[12] = N21[12] = -0.25*nu.x;
				N12[13] = N21[13] = -0.25*nu.y;
				N12[14] = N21[14] = -0.25*nu.z;
				N12[15] = N21[15] = 0;
				N12[16] = N21[16] = 0.25*nu.x;
				N12[17] = N21[17] = 0.25*nu.y;
				N12[18] = N21[18] = 0.25*nu.z;
				N12[19] = N21[19] = 0;
			}
			else if (nmeln == 6){
				N12[4] = N21[4] = -4*nu.x;
				N12[5] = N21[5] = -4*nu.y;
				N12[6] = N21[6] = -4*nu.z;
				N12[7] = N21[7] = 0;
				N12[8] = N21[8] = 0;
				N12[9] = N21[9] = 0;
				N12[10] = N21[10] = 0;
				N12[11] = N21[11] = 0;
				N12[12] = N21[12] = 0;
				N12[13] = N21[13] = 0;
				N12[14] = N21[14] = 0;
				N12[15] = N21[15] = 0;
				N12[16] = N21[16] = 4*nu.x;
				N12[17] = N21[17] = 4*nu.y;
				N12[18] = N21[18] = 4*nu.z;
				N12[19] = N21[19] = 0;
				N12[20] = N21[20] = -4*nu.x;
				N12[21] = N21[21] = -4*nu.y;
				N12[22] = N21[22] = -4*nu.z;
				N12[23] = N21[23] = 0;
				N12[24] = N21[24] = 4*nu.x;
				N12[25] = N21[25] = 4*nu.y;
				N12[26] = N21[26] = 4*nu.z;
				N12[27] = N21[27] = 0;


				N11[4] = -4 * nu.x;
				N11[5] = -4 * nu.y;
				N11[6] = -4 * nu.z;
				N11[7] = 0;
				N11[8] = -4 * nu.x;
				N11[9] = -4 * nu.y;
				N11[10] = -4 * nu.z;
				N11[11] = 0;
				N11[12] = 0;
				N11[13] = 0;
				N11[14] = 0;
				N11[15] = 0;
				N11[16] = 8 * nu.x;
				N11[17] = 8 * nu.y;
				N11[18] = 8 * nu.z;
				N11[19] = 0;
				N11[20] = 0;
				N11[21] = 0;
				N11[22] = 0;
				N11[23] = 0;
				N11[24] = 0;
				N11[25] = 0;
				N11[26] = 0;
				N11[27] = 0;


				N22[4] = -4 * nu.x;
				N22[5] = -4 * nu.y;
				N22[6] = -4 * nu.z;
				N22[7] = 0;
				N22[8] = 0;
				N22[9] = 0;
				N22[10] = 0;
				N22[11] = 0;
				N22[12] = -4 * nu.x;
				N22[13] = -4 * nu.y;
				N22[14] = -4 * nu.z;
				N22[15] = 0;
				N22[16] = 0;
				N22[17] = 0;
				N22[18] = 0;
				N22[19] = 0;
				N22[20] = 0;
				N22[21] = 0;
				N22[22] = 0;
				N22[23] = 0;
				N22[24] = 8 * nu.x;
				N22[25] = 8 * nu.y;
				N22[26] = 8 * nu.z;
				N22[27] = 0;


			}
			vec3d g12(0, 0, 0);
			if (nmeln == 4)
			{
				const double Grs[4] = { 0.25, -0.25, 0.25, -0.25 };
				g12 = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3];
			}
			else if (nmeln == 6)
			{
				const double Grs[6] = { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };
				g12 = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3] + rt[4] * Grs[4] + rt[5] * Grs[5];
			}

			double gt1 = g12*tau[0];
			double gt2 = g12*tau[1];
			double gp = g12*pt;

			for (k = 0; k < ndof; ++k)
			{
				Tb11[k] = T11[k] - gt1*D2[k];
				Tb12[k] = T12[k] - gt1*D1[k];

				Tb21[k] = T21[k] - gt2*D2[k];
				Tb22[k] = T22[k] - gt2*D1[k];

				Pb1[k] = P1[k] - gp*D2[k];
				Pb2[k] = P2[k] - gp*D1[k];
			}

			// raise the indices of A
			double Ac[2][2];
			for (k = 0; k < 2; ++k)
				for (l = 0; l < 2; ++l)
				{
				Ac[k][l] = 0;
				for (i = 0; i < 2; ++i)
					for (j = 0; j < 2; ++j) Ac[k][l] += Mki[k][i] * Mki[l][j] * A[i][j];
				}

			vec3d Hrs[2][2] = { { vec3d(0, 0, 0), vec3d(0, 0, 0) }, { vec3d(0, 0, 0), vec3d(0, 0, 0) } };
			if (nmeln == 4)
			{
				const double Grs[4] = { 0.25, -0.25, 0.25, -0.25 };
				Hrs[0][1] = Hrs[1][0] = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3];
			}
			else if (nmeln == 6)
			{
				const double Grs[6] = { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };
				Hrs[0][1] = Hrs[1][0] = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3] + rt[4] * Grs[4] + rt[5] * Grs[5];
			}
			///////////

			double kij;
			for (i = 0; i < ndof; ++i)
				for (j = 0; j < ndof; ++j)
				{
				// KT1
				kij = T11[i] * D1[j] + T12[i] * D2[j];
				kij += D1[i] * T11[j] + D2[i] * T12[j];
				kij -= (Hrs[0][1] * tau[0])*D1[i] * D2[j] + (Hrs[1][0] * tau[0])*D2[i] * D1[j];
				kij += Tb11[i] * D1[j] + Tb21[i] * D2[j];
				kij += D1[i] * Tb11[j] + D2[i] * Tb21[j];
				kij += gap*(N11[i] * D1[j] + N12[i] * D2[j] + D1[i] * N11[j] + D2[i] * N12[j]);
				kij -= NTE[i] * Nb1[j] - Nb1[i] * NTE[j];
				kij -= T1[i] * Mi[0][0] * Tb11[j] + T1[i] * Mi[0][1] * Tb21[j] + T2[i] * Mi[1][0] * Tb11[j] + T2[i] * Mi[1][1] * Tb21[j];
				kij -= Tb11[i] * Mi[0][0] * T1[j] + Tb21[i] * Mi[0][1] * T1[j] + Tb11[i] * Mi[1][0] * T2[j] + Tb21[i] * Mi[1][1] * T2[j];

				ke[i][j] += m_ktmult*(Tt[0] * Ac[0][0] + Tt[1] * Ac[1][0])*kij;

				// KT2
				kij = T21[i] * D1[j] + T22[i] * D2[j];
				kij += D1[i] * T21[j] + D2[i] * T22[j];
				kij -= (Hrs[0][1] * tau[1])*D1[i] * D2[j] + (Hrs[1][0] * tau[1])*D2[i] * D1[j];
				kij += Tb12[i] * D1[j] + Tb22[i] * D2[j];
				kij += D1[i] * Tb12[j] + D2[i] * Tb22[j];
				kij += gap*(N21[i] * D1[j] + N22[i] * D2[j] + D1[i] * N21[j] + D2[i] * N22[j]);
				kij -= NTE[i] * Nb2[j] - Nb2[i] * NTE[j];
				kij -= T1[i] * Mi[0][0] * Tb12[j] + T1[i] * Mi[0][1] * Tb22[j] + T2[i] * Mi[1][0] * Tb12[j] + T2[i] * Mi[1][1] * Tb22[j];
				kij -= Tb12[i] * Mi[0][0] * T1[j] + Tb22[i] * Mi[0][1] * T1[j] + Tb12[i] * Mi[1][0] * T2[j] + Tb22[i] * Mi[1][1] * T2[j];

				ke[i][j] += m_ktmult*(Tt[0] * Ac[0][1] + Tt[1] * Ac[1][1])*kij;

				double tempt = tn;
				// kdirect
				if ((tempt != 0) && (!bstick) && (η > 0) && (TMT > 0) && (κ_us > 0) && (∆time>0)){
					double W_1 = m_epsf*∆time*(Tt[0]) / tempt*η*sqrt(TMT);
					double W_2 = m_epsf*∆time*(Tt[1]) / tempt*η*sqrt(TMT);
					//ke[i][j]_α
					ke[i][j] += (D1[i])*e[j] * W_1*(1. / 1 + ((m_epsf*∆time) / tempt*η))*tempt*κ_us*m_mu;

					//ke[i][j]_β
					ke[i][j] += (D2[i])*e[j] * W_2*(1. / 1 + ((m_epsf*∆time) / tempt*η))*tempt*κ_us*m_mu;


				}
				if (bstick)
				{
					kij = Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j] + Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j];
					ke[i][j] += m_epsf*kij;
				}
				else
				{
					kij = (1.0 - Ptc[0] * Pt[0])*(Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j]);
					kij += (-Ptc[0] * Pt[1])*(Mk[1][0] * D1[i] * D1[j] + Mk[1][1] * D1[i] * D2[j]);
					kij += (-Ptc[1] * Pt[0])*(Mk[0][0] * D2[i] * D1[j] + Mk[0][1] * D2[i] * D2[j]);
					kij += (1.0 - Ptc[1] * Pt[1])*(Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j]);

					if (TMT > 0) {
						ke[i][j] += m_ktmult*m_epsf*m_mu*(tn / sqrt(TMT))*kij;
					}

				}
					///Temperature part begin
					// calculate Temperature force vector  and Mechanical dissipation  if 
					// only if both the friction coefficient and friction
					// penalty factor are non-zero

					//if (isinf(x1) != 0){ x1 = 0; }
					//if (isinf(x2) != 0){ x2 = 0; }
					double DM_α = fabs(Tt[0] * x1);
					double DM_β = fabs(Tt[1] * x2);

					if (∆time > 0){


						/// conduction stiffness

						ke[i][j] += R0* tempt*a[i] * a[j] - R0*(a[i] * c[i])*a[i] * NTE[j] * eps*HEAVYSIDE(gap);

						//ke[i][j]_α
						ke[i][j] += R0* tempt*(a[i] * c[i])*(d1[i] * D1[j]) - R0* tempt*(d1[i] * c[i])*a[i] * (D1[j]);

						//ke[i][j]_β
						ke[i][j] += R0* tempt*(a[i] * c[i])*(d2[i] * D2[j]) - R0* tempt*(d2[i] * c[i])*a[i] * (D2[j]);


						/// dissipation stiffness

						//ke[i][j]_α
						ke[i][j] += (R[1] / ∆time) *DM_α*(d1[i] * D1[j]) + sqrt(Tt[0] * Tt[0])*b[i] * (D1[j]) / ∆time;

						//ke[i][j]_β
						ke[i][j] += (R[1] / ∆time) *DM_β*(d2[i] * D2[j]) + sqrt(Tt[1] * Tt[1])*b[i] * (D2[j]) / ∆time;




						/// trapped debris(heat sinks) stiffness

						ke[i][j] += R[2] * tempt*e[i] * e[j] + R[4] * tempt*f[i] * f[j] + (-R[2] * (c[i] * e[i]) + R[3])*e[i] * NTE[j] * eps*HEAVYSIDE(gap) + (-R[4] * (c[i] * f[i]) + R[5])*f[i] * NTE[j] * eps*HEAVYSIDE(gap);
						//ke[i][j]_α
						ke[i][j] += (R[4] * (c[i] * f[i]) - R[5])*tempt*(d1[i] * D1[j]);

						//ke[i][j]_β
						ke[i][j] += (R[4] * (c[i] * f[i]) - R[5])*tempt*(d2[i] * D2[j]);

					}


				
					///Temperature part end
				

				}
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////////////

	if ((m_mu*m_epsf > 0) && (gap <= 0) && (fabs(gap) <= (fric_tol*δ_TN)))
		
	{
		///////////////////   modified part begin
		double max_normal = M_TN* ss.D_R[m];
		double max_tangential = M_Tt* ss.D_R[m];

		double mCL_normal = δ_TN * ss.D_R_C[m];
		double mCL_tangential = δ_Tt * ss.D_R_C[m];


		// get the metric tensor and its inverse
		mat2d& Mk = ss.m_M[m];
		mat2d Mki = Mk.inverse();

		// get the previous isoparameteric coordinates
		double rp = ss.m_rsp[m][0];
		double sp = ss.m_rsp[m][1];
		double gap_n = fabs(gap);



		double Utn_α = -((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
		double Utn_β = -((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
		double Utt_α = (exp(-(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
		double Utt_β = (exp(-(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));

		double t_N = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)));
		double t_N_α = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_α);
		double t_N_β = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_β);
		double tracn = t_N_α + t_N_β;
		// contact traction( cohesive case! mode II and III)
		double t_T = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)));
		double t_T_α = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_α;
		double t_T_β = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_β;
		double tn_x = -(MBRACKET(-tracn));
		double tt_α = -(MBRACKET(-t_T_α));
		double tt_β = -(MBRACKET(-t_T_β));
	     t_N = -(MBRACKET(-t_N));
		 t_N_α = -(MBRACKET(-t_N_α));
		 t_N_β = -(MBRACKET(-t_N_β));
		 t_T = -(MBRACKET(-t_T));
		//double tract = t_T_α + t_T_β;

		//double tractt1 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gapx / eps))*exp(-(gapx / eps)))*Utxt;
		//double tractt2 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gapx / eps))*exp(-(gapx / eps)))*Utyt;
		 //double tractt = tractt1 + tractt2;
		 //double abs_α = (Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)) / sqrt(((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))));
		 //double abs_β = (Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)) / sqrt(((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))));



		double abs_α = ((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)) + 1) / (sqrt(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) + 1);
		double abs_β = ((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)) + 1) / (sqrt(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) + 1);

		 if (abs_α < 1){ abs_α = -1; }

		 if (abs_β < 1){ abs_β = -1; }




		 double DUtn_α = ((-2 * (((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*exp(-((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential));
		 double DUtn_β = ((-2 * (((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*exp(-((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential));
		 double DUtt_α = (exp(-(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*(abs_α / mCL_tangential) + ((-2 * (((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*(exp(-(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
		 double DUtt_β = (exp(-(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*(abs_β / mCL_tangential) + ((-2 * (((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*(exp(-(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));



		double Lt[2];
		Lt[0] = ss.m_trans_S[m][0];
		Lt[1] = ss.m_trans_S[m][1];


		/////////////////////
		//double rpf = ss.m_rsp_f[m][0];
		//double spf = ss.m_rsp_f[m][1];

		// get the traction
		double Tt[2];
		// a. trial state
		double x1 = (r - rp);
		double x2 = (s - sp);

		Tt[0] = Lt[0] + m_epsf*(Mk[0][0] * x1 + Mk[0][1] * x2);
		Tt[1] = Lt[1] + m_epsf*(Mk[1][0] * x1 + Mk[1][1] * x2);
		//Tt[0] = Lt[0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp);
		//Tt[1] = Lt[1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp);

		//if (isinf(Tt[0]) != 0){ Tt[0] = 0; }
		//if (isinf(Tt[1]) != 0){ Tt[1] = 0; }

		double TMT = Tt[0] * (Mki[0][0] * Tt[0] + Mki[0][1] * Tt[1]) + Tt[1] * (Mki[1][0] * Tt[0] + Mki[1][1] * Tt[1]);

		// calculate the normalized traction vector
		double Pt[2] = { 0, 0 };

		if (TMT > 0){
			Pt[0] = Tt[0] / sqrt(TMT);
			Pt[1] = Tt[1] / sqrt(TMT);
		}


		/////////////////////
		// calculate the covariant version
		
		double Ptc[2];
		Ptc[0] = Mki[0][0] * Pt[0] + Mki[0][1] * Pt[1];
		Ptc[1] = Mki[1][0] * Pt[0] + Mki[1][1] * Pt[1];

		//b. return map
		bool bstick = true;
		double phi = 0;

		if (TMT >= 0){
			phi = sqrt(TMT) - m_mu*fabs(tn_x);
		}
		
		if ((phi > 0) && (TMT > 0))
		{
			Tt[0] = m_mu*Tt[0] / sqrt(TMT);
			Tt[1] = m_mu*Tt[1] / sqrt(TMT);
			bstick = false;
		}

		// get the previous isoparameteric coordinates
		//	double rp = ss.m_rsp[m][0];
		//	double sp = ss.m_rsp[m][1];
		vec3d pt = tau[0] * Tt[0] + tau[1] * Tt[1];
		pt.unit();
		for (k = 0; k < nmeln; ++k)
		{
			T11[4 * (k)+4] = -Hr[k] * tau[0].x;
			T11[4 * (k)+5] = -Hr[k] * tau[0].y;
			T11[4 * (k)+6] = -Hr[k] * tau[0].z;
			T11[4 * (k)+7] = 0;

			T12[4 * (k)+4] = -Hs[k] * tau[0].x;
			T12[4 * (k)+5] = -Hs[k] * tau[0].y;
			T12[4 * (k)+6] = -Hs[k] * tau[0].z;
			T12[4 * (k)+7] = 0;

			T21[4 * (k)+4] = -Hr[k] * tau[1].x;
			T21[4 * (k)+5] = -Hr[k] * tau[1].y;
			T21[4 * (k)+6] = -Hr[k] * tau[1].z;
			T21[4 * (k)+7] = 0;

			T22[4 * (k)+4] = -Hs[k] * tau[1].x;
			T22[4 * (k)+5] = -Hs[k] * tau[1].y;
			T22[4 * (k)+6] = -Hs[k] * tau[1].z;
			T22[4 * (k)+7] = 0;

			//if (nmeln == 4)
			//{
			//	N12[4 * (k)+4] = N21[4 * (k)+4] = -0.25*nu.x;
			//	N12[4 * (k)+5] = N21[4 * (k)+5] = -0.25*nu.y;
			//	N12[4 * (k)+6] = N21[4 * (k)+6] = -0.25*nu.z;
			//	N12[4 * (k)+7] = N21[4 * (k)+7] = 0;
			//}
			//else if (nmeln == 6) assert(false);

			P1[4 * (k)+4] = -Hr[k] * pt.x;
			P1[4 * (k)+5] = -Hr[k] * pt.y;
			P1[4 * (k)+6] = -Hr[k] * pt.z;
			P1[4 * (k)+7] = 0;

			P2[4 * (k)+4] = -Hs[k] * pt.x;
			P2[4 * (k)+5] = -Hs[k] * pt.y;
			P2[4 * (k)+6] = -Hs[k] * pt.z;
			P2[4 * (k)+7] = 0;
		}
		if (nmeln == 4)
		{
			N12[4] = N21[4] = -0.25*nu.x;
			N12[5] = N21[5] = -0.25*nu.y;
			N12[6] = N21[6] = -0.25*nu.z;
			N12[7] = N21[7] = 0;
			N12[8] = N21[8] = 0.25*nu.x;
			N12[9] = N21[9] = 0.25*nu.y;
			N12[10] = N21[10] = 0.25*nu.z;
			N12[11] = N21[11] = 0;
			N12[12] = N21[12] = -0.25*nu.x;
			N12[13] = N21[13] = -0.25*nu.y;
			N12[14] = N21[14] = -0.25*nu.z;
			N12[15] = N21[15] = 0;
			N12[16] = N21[16] = 0.25*nu.x;
			N12[17] = N21[17] = 0.25*nu.y;
			N12[18] = N21[18] = 0.25*nu.z;
			N12[19] = N21[19] = 0;
		}
		else if (nmeln == 6){
			N12[4] = N21[4] = -4 * nu.x;
			N12[5] = N21[5] = -4 * nu.y;
			N12[6] = N21[6] = -4 * nu.z;
			N12[7] = N21[7] = 0;
			N12[8] = N21[8] = 0;
			N12[9] = N21[9] = 0;
			N12[10] = N21[10] = 0;
			N12[11] = N21[11] = 0;
			N12[12] = N21[12] = 0;
			N12[13] = N21[13] = 0;
			N12[14] = N21[14] = 0;
			N12[15] = N21[15] = 0;
			N12[16] = N21[16] = 4 * nu.x;
			N12[17] = N21[17] = 4 * nu.y;
			N12[18] = N21[18] = 4 * nu.z;
			N12[19] = N21[19] = 0;
			N12[20] = N21[20] = -4 * nu.x;
			N12[21] = N21[21] = -4 * nu.y;
			N12[22] = N21[22] = -4 * nu.z;
			N12[23] = N21[23] = 0;
			N12[24] = N21[24] = 4 * nu.x;
			N12[25] = N21[25] = 4 * nu.y;
			N12[26] = N21[26] = 4 * nu.z;
			N12[27] = N21[27] = 0;


			N11[4] = -4 * nu.x;
			N11[5] = -4 * nu.y;
			N11[6] = -4 * nu.z;
			N11[7] = 0;
			N11[8] = -4 * nu.x;
			N11[9] = -4 * nu.y;
			N11[10] = -4 * nu.z;
			N11[11] = 0;
			N11[12] = 0;
			N11[13] = 0;
			N11[14] = 0;
			N11[15] = 0;
			N11[16] = 8 * nu.x;
			N11[17] = 8 * nu.y;
			N11[18] = 8 * nu.z;
			N11[19] = 0;
			N11[20] = 0;
			N11[21] = 0;
			N11[22] = 0;
			N11[23] = 0;
			N11[24] = 0;
			N11[25] = 0;
			N11[26] = 0;
			N11[27] = 0;


			N22[4] = -4 * nu.x;
			N22[5] = -4 * nu.y;
			N22[6] = -4 * nu.z;
			N22[7] = 0;
			N22[8] = 0;
			N22[9] = 0;
			N22[10] = 0;
			N22[11] = 0;
			N22[12] = -4 * nu.x;
			N22[13] = -4 * nu.y;
			N22[14] = -4 * nu.z;
			N22[15] = 0;
			N22[16] = 0;
			N22[17] = 0;
			N22[18] = 0;
			N22[19] = 0;
			N22[20] = 0;
			N22[21] = 0;
			N22[22] = 0;
			N22[23] = 0;
			N22[24] = 8 * nu.x;
			N22[25] = 8 * nu.y;
			N22[26] = 8 * nu.z;
			N22[27] = 0;


		}
		vec3d g12(0, 0, 0);
		if (nmeln == 4)
		{
			const double Grs[4] = { 0.25, -0.25, 0.25, -0.25 };
			g12 = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3];
		}
		else if (nmeln == 6)
		{
			const double Grs[6] = { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };
			g12 = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3] + rt[4] * Grs[4] + rt[5] * Grs[5];
		}

		double gt1 = g12*tau[0];
		double gt2 = g12*tau[1];
		double gp = g12*pt;

		for (k = 0; k < ndof; ++k)
		{
			Tb11[k] = T11[k] - gt1*D2[k];
			Tb12[k] = T12[k] - gt1*D1[k];

			Tb21[k] = T21[k] - gt2*D2[k];
			Tb22[k] = T22[k] - gt2*D1[k];

			Pb1[k] = P1[k] - gp*D2[k];
			Pb2[k] = P2[k] - gp*D1[k];
		}

		// raise the indices of A
		double Ac[2][2];
		for (k = 0; k < 2; ++k)
			for (l = 0; l < 2; ++l)
			{
			Ac[k][l] = 0;
			for (i = 0; i < 2; ++i)
				for (j = 0; j < 2; ++j) Ac[k][l] += Mki[k][i] * Mki[l][j] * A[i][j];
			}

		vec3d Hrs[2][2] = { { vec3d(0, 0, 0), vec3d(0, 0, 0) }, { vec3d(0, 0, 0), vec3d(0, 0, 0) } };
		if (nmeln == 4)
		{
			const double Grs[4] = { 0.25, -0.25, 0.25, -0.25 };
			Hrs[0][1] = Hrs[1][0] = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3];
		}
		else if (nmeln == 6)
		{
			const double Grs[6] = { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };
			Hrs[0][1] = Hrs[1][0] = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3] + rt[4] * Grs[4] + rt[5] * Grs[5];
		}
		/////////=
		//  m_ktmult is a user defined control function(default value m_ktmult = 1;)
		double kij;
		for (i = 0; i < ndof; ++i)
			for (j = 0; j < ndof; ++j)
			{
			// KT_α
			kij = T11[i] * D1[j] + T12[i] * D2[j];
			kij += D1[i] * T11[j] + D2[i] * T12[j];
			kij -= (Hrs[0][1] * tau[0])*D1[i] * D2[j] + (Hrs[1][0] * tau[0])*D2[i] * D1[j];
			kij += Tb11[i] * D1[j] + Tb21[i] * D2[j];
			kij += D1[i] * Tb11[j] + D2[i] * Tb21[j];
			kij += gap*(N11[i] * D1[j] + N12[i] * D2[j] + D1[i] * N11[j] + D2[i] * N12[j]);
			kij -= NTE[i] * Nb1[j] - Nb1[i] * NTE[j];
			kij -= T1[i] * Mi[0][0] * Tb11[j] + T1[i] * Mi[0][1] * Tb21[j] + T2[i] * Mi[1][0] * Tb11[j] + T2[i] * Mi[1][1] * Tb21[j];
			kij -= Tb11[i] * Mi[0][0] * T1[j] + Tb21[i] * Mi[0][1] * T1[j] + Tb11[i] * Mi[1][0] * T2[j] + Tb21[i] * Mi[1][1] * T2[j];

			ke[i][j] += m_ktmult*(Tt[0] * Ac[0][0] + Tt[1] * Ac[1][0])*kij;

			// KT_β
			kij = T21[i] * D1[j] + T22[i] * D2[j];
			kij += D1[i] * T21[j] + D2[i] * T22[j];
			kij -= (Hrs[0][1] * tau[1])*D1[i] * D2[j] + (Hrs[1][0] * tau[1])*D2[i] * D1[j];
			kij += Tb12[i] * D1[j] + Tb22[i] * D2[j];
			kij += D1[i] * Tb12[j] + D2[i] * Tb22[j];
			kij += gap*(N21[i] * D1[j] + N22[i] * D2[j] + D1[i] * N21[j] + D2[i] * N22[j]);
			kij -= NTE[i] * Nb2[j] - Nb2[i] * NTE[j];
			kij -= T1[i] * Mi[0][0] * Tb12[j] + T1[i] * Mi[0][1] * Tb22[j] + T2[i] * Mi[1][0] * Tb12[j] + T2[i] * Mi[1][1] * Tb22[j];
			kij -= Tb12[i] * Mi[0][0] * T1[j] + Tb22[i] * Mi[0][1] * T1[j] + Tb12[i] * Mi[1][0] * T2[j] + Tb22[i] * Mi[1][1] * T2[j];

			ke[i][j] += m_ktmult*(Tt[0] * Ac[0][1] + Tt[1] * Ac[1][1])*kij;

			double  tempt = fabs(tn_x);
			// kdirect
			if ((tempt != 0) && (!bstick) && (η > 0) && (TMT > 0) && (κ_us > 0) && (∆time>0)){
				double W_1 = m_epsf*∆time*(Tt[0]) / tempt*η*sqrt(TMT);
				double W_2 = m_epsf*∆time*(Tt[1]) / tempt*η*sqrt(TMT);
				//ke[i][j]_α
				ke[i][j] += (D1[i])*e[j] * W_1*(1. / 1 + ((m_epsf*∆time) / tempt*η))*tempt*κ_us*m_mu;

				//ke[i][j]_β
				ke[i][j] += (D2[i])*e[j] * W_2*(1. / 1 + ((m_epsf*∆time) / tempt*η))*tempt*κ_us*m_mu;


			}
			if (bstick)
			{
				kij = Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j] + Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j];
				ke[i][j] += m_epsf*kij;
			}
			else
			{
				kij = (1.0 - Ptc[0] * Pt[0])*(Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j]);
				kij += (-Ptc[0] * Pt[1])*(Mk[1][0] * D1[i] * D1[j] + Mk[1][1] * D1[i] * D2[j]);
				kij += (-Ptc[1] * Pt[0])*(Mk[0][0] * D2[i] * D1[j] + Mk[0][1] * D2[i] * D2[j]);
				kij += (1.0 - Ptc[1] * Pt[1])*(Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j]);

				if (TMT > 0) { ke[i][j] += m_ktmult*m_epsf*m_mu*(tn / sqrt(TMT))*kij; }


			}
				///Temperature part begin
				// calculate Temperature force vector  and Mechanical dissipation  if 
				// only if both the friction coefficient and friction


			//if (isinf(x1) != 0){ x1 = 0; }
			//if (isinf(x2) != 0){ x2 = 0; }
			double DM_α = fabs(Tt[0] * x1);
			double DM_β = fabs(Tt[1] * x2);

			if (∆time > 0){

				double dtntrail = -HEAVYSIDE(((gap_n / mCL_normal) * max_normal  * exp(1 - (gap_n / mCL_normal))))*(((max_normal / mCL_normal)*exp(1 - (gap_n / mCL_normal)) + (gap_n / mCL_normal) * max_normal *(-1 / mCL_normal)*exp(1 - (gap_n / mCL_normal))));
				if (dtntrail < 0){ dtntrail = -dtntrail; }
				double ∆tnn = dtntrail*(exp(Utt_α) + exp(Utn_β));
				double ∆tnt = t_N*DUtn_α*(Mk[0][0] * NTE[i] * D1[i] + Mk[0][1] * NTE[i] * D2[j]) + t_N*DUtn_β*(Mk[1][0] * NTE[i] * D1[j] + Mk[1][1] * NTE[i] * D2[j]);

				/// conduction stiffness

				ke[i][j] += R0* tempt*a[i] * a[j] - R0*(a[i] * c[i])*a[i] * (NTE[j] * ∆tnn + ∆tnt);

				//ke[i][j]_α
				ke[i][j] += R0* tempt*(a[i] * c[i])*(d1[i] * D1[j]) - R0* tempt*(d1[i] * c[i])*a[i] * (D1[j]);

				//ke[i][j]_β
				ke[i][j] += R0* tempt*(a[i] * c[i])*(d2[i] * D2[j]) - R0* tempt*(d2[i] * c[i])*a[i] * (D2[j]);

				/// dissipation stiffness

				//ke[i][j]_α
				ke[i][j] += (R[1] / ∆time) *DM_α*(d1[i] * D1[j]) + sqrt(Tt[0] * Tt[0])*b[i] * (D1[j]) / ∆time;

				//ke[i][j]_β
				ke[i][j] += (R[1] / ∆time) *DM_β*(d2[i] * D2[j]) + sqrt(Tt[1] * Tt[1])*b[i] * (D2[j]) / ∆time;


				/// trapped debris(heat sinks) stiffness

				ke[i][j] += R[2] * tempt*e[i] * e[j] + R[4] * tempt*f[i] * f[j] + (-R[2] * (c[i] * e[i]) + R[3])*e[i] * (NTE[j] * ∆tnn + ∆tnt) + (-R[4] * (c[i] * f[i]) + R[5])*f[i] * (NTE[j] * ∆tnn + ∆tnt);
				//ke[i][j]_α
				ke[i][j] += (R[4] * (c[i] * f[i]) - R[5])*tempt*(d1[i] * D1[j]);

				//ke[i][j]_β
				ke[i][j] += (R[4] * (c[i] * f[i]) - R[5])*tempt*(d2[i] * D2[j]);

			}
				

			}

			


	}



///////////////////////////////////////////////////////////////////////////////////////////////	
	if (gap <= 0){
		////////////////////
		
		//test
		
		double max_normal = M_TN* ss.D_R[m];
		double max_tangential = M_Tt* ss.D_R[m];

		double mCL_normal = δ_TN * ss.D_R_C[m];
		double mCL_tangential = δ_Tt * ss.D_R_C[m];
		
		//////////////////
		// get the traction multipliers
		//	double Lt[2];
		//	Lt[0] = ss.m_Lt[m][0];
		//	Lt[1] = ss.m_Lt[m][1];

		// get the metric tensor and its inverse
		mat2d& Mk = ss.m_M[m];
		mat2d Mki = Mk.inverse();

		// get the previous isoparameteric coordinates
		double rp = ss.m_rsp[m][0];
		double sp = ss.m_rsp[m][1];


		//pt.unit();
		for ( k = 0; k < nmeln; ++k)
		{
			T11[4 * (k)+4] = -Hr[k] * tau[0].x;
			T11[4 * (k)+5] = -Hr[k] * tau[0].y;
			T11[4 * (k)+6] = -Hr[k] * tau[0].z;
			T11[4 * (k)+7] = 0;

			T12[4 * (k)+4] = -Hs[k] * tau[0].x;
			T12[4 * (k)+5] = -Hs[k] * tau[0].y;
			T12[4 * (k)+6] = -Hs[k] * tau[0].z;
			T12[4 * (k)+7] = 0;

			T21[4 * (k)+4] = -Hr[k] * tau[1].x;
			T21[4 * (k)+5] = -Hr[k] * tau[1].y;
			T21[4 * (k)+6] = -Hr[k] * tau[1].z;
			T21[4 * (k)+7] = 0;

			T22[4 * (k)+4] = -Hs[k] * tau[1].x;
			T22[4 * (k)+5] = -Hs[k] * tau[1].y;
			T22[4 * (k)+6] = -Hs[k] * tau[1].z;
			T22[4 * (k)+7] = 0;

			//if (nmeln == 4)
			//{
			//	N12[4 * (k)+4] = N21[4 * (k)+4] = -0.25*nu.x;
			//	N12[4 * (k)+5] = N21[4 * (k)+5] = -0.25*nu.y;
			//	N12[4 * (k)+6] = N21[4 * (k)+6] = -0.25*nu.z;
			//	N12[4 * (k)+7] = N21[4 * (k)+7] = 0;
			//}
			//else if (nmeln == 6) assert(false);

			//P1[4 * (k)+4] = -Hr[k] * pt.x;
			//P1[4 * (k)+5] = -Hr[k] * pt.y;
			///P1[4 * (k)+6] = -Hr[k] * pt.z;
			//P1[4 * (k)+7] = 0;

			//P2[4 * (k)+4] = -Hs[k] * pt.x;
			//P2[4 * (k)+5] = -Hs[k] * pt.y;
			//P2[4 * (k)+6] = -Hs[k] * pt.z;
			//P2[4 * (k)+7] = 0;
		}
		if (nmeln == 4)
		{
			N12[4] = N21[4] = -0.25*nu.x;
			N12[5] = N21[5] = -0.25*nu.y;
			N12[6] = N21[6] = -0.25*nu.z;
			N12[7] = N21[7] = 0;
			N12[8] = N21[8] = 0.25*nu.x;
			N12[9] = N21[9] = 0.25*nu.y;
			N12[10] = N21[10] = 0.25*nu.z;
			N12[11] = N21[11] = 0;
			N12[12] = N21[12] = -0.25*nu.x;
			N12[13] = N21[13] = -0.25*nu.y;
			N12[14] = N21[14] = -0.25*nu.z;
			N12[15] = N21[15] = 0;
			N12[16] = N21[16] = 0.25*nu.x;
			N12[17] = N21[17] = 0.25*nu.y;
			N12[18] = N21[18] = 0.25*nu.z;
			N12[19] = N21[19] = 0;
		}
		else if (nmeln == 6){
			N12[4] = N21[4] = -4 * nu.x;
			N12[5] = N21[5] = -4 * nu.y;
			N12[6] = N21[6] = -4 * nu.z;
			N12[7] = N21[7] = 0;
			N12[8] = N21[8] = 0;
			N12[9] = N21[9] = 0;
			N12[10] = N21[10] = 0;
			N12[11] = N21[11] = 0;
			N12[12] = N21[12] = 0;
			N12[13] = N21[13] = 0;
			N12[14] = N21[14] = 0;
			N12[15] = N21[15] = 0;
			N12[16] = N21[16] = 4 * nu.x;
			N12[17] = N21[17] = 4 * nu.y;
			N12[18] = N21[18] = 4 * nu.z;
			N12[19] = N21[19] = 0;
			N12[20] = N21[20] = -4 * nu.x;
			N12[21] = N21[21] = -4 * nu.y;
			N12[22] = N21[22] = -4 * nu.z;
			N12[23] = N21[23] = 0;
			N12[24] = N21[24] = 4 * nu.x;
			N12[25] = N21[25] = 4 * nu.y;
			N12[26] = N21[26] = 4 * nu.z;
			N12[27] = N21[27] = 0;


			N11[4] = -4 * nu.x;
			N11[5] = -4 * nu.y;
			N11[6] = -4 * nu.z;
			N11[7] = 0;
			N11[8] = -4 * nu.x;
			N11[9] = -4 * nu.y;
			N11[10] = -4 * nu.z;
			N11[11] = 0;
			N11[12] = 0;
			N11[13] = 0;
			N11[14] = 0;
			N11[15] = 0;
			N11[16] = 8 * nu.x;
			N11[17] = 8 * nu.y;
			N11[18] = 8 * nu.z;
			N11[19] = 0;
			N11[20] = 0;
			N11[21] = 0;
			N11[22] = 0;
			N11[23] = 0;
			N11[24] = 0;
			N11[25] = 0;
			N11[26] = 0;
			N11[27] = 0;


			N22[4] = -4 * nu.x;
			N22[5] = -4 * nu.y;
			N22[6] = -4 * nu.z;
			N22[7] = 0;
			N22[8] = 0;
			N22[9] = 0;
			N22[10] = 0;
			N22[11] = 0;
			N22[12] = -4 * nu.x;
			N22[13] = -4 * nu.y;
			N22[14] = -4 * nu.z;
			N22[15] = 0;
			N22[16] = 0;
			N22[17] = 0;
			N22[18] = 0;
			N22[19] = 0;
			N22[20] = 0;
			N22[21] = 0;
			N22[22] = 0;
			N22[23] = 0;
			N22[24] = 8 * nu.x;
			N22[25] = 8 * nu.y;
			N22[26] = 8 * nu.z;
			N22[27] = 0;


		}
		vec3d g12(0, 0, 0);
		if (nmeln == 4)
		{
			const double Grs[4] = { 0.25, -0.25, 0.25, -0.25 };
			g12 = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3];
		}
		else if (nmeln == 6)
		{
			const double Grs[6] = { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };
			g12 = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3] + rt[4] * Grs[4] + rt[5] * Grs[5];
		}

		double gt1 = g12*tau[0];
		double gt2 = g12*tau[1];
		//double gp = g12*pt;///

		for ( k = 0; k < ndof; ++k)
		{
			Tb11[k] = T11[k] - gt1*D2[k];
			Tb12[k] = T12[k] - gt1*D1[k];

			Tb21[k] = T21[k] - gt2*D2[k];
			Tb22[k] = T22[k] - gt2*D1[k];

			//Pb1[k] = P1[k] - gp*D2[k];
			//Pb2[k] = P2[k] - gp*D1[k];
		}

		// raise the indices of A
		double Ac[2][2];
		for ( k = 0; k < 2; ++k)
			for ( l = 0; l < 2; ++l)
			{
			Ac[k][l] = 0;
			for ( i = 0; i < 2; ++i)
				for ( j = 0; j < 2; ++j) Ac[k][l] += Mki[k][i] * Mki[l][j] * A[i][j];
			}

		vec3d Hrs[2][2] = { { vec3d(0, 0, 0), vec3d(0, 0, 0) }, { vec3d(0, 0, 0), vec3d(0, 0, 0) } };
		if (nmeln == 4)
		{
			const double Grs[4] = { 0.25, -0.25, 0.25, -0.25 };
			Hrs[0][1] = Hrs[1][0] = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3];
		}
		else if (nmeln == 6)
		{
			const double Grs[6] = { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };
			Hrs[0][1] = Hrs[1][0] = rt[0] * Grs[0] + rt[1] * Grs[1] + rt[2] * Grs[2] + rt[3] * Grs[3] + rt[4] * Grs[4] + rt[5] * Grs[5];
		}
		/////////=

		///////////////////
		double gap_n = fabs(gap);

		// --- NORMAL and T A N G E N T I A L   S T I F F N E S S ---

		// the normal opening traction(coupled mode I cracks), the tangential opening traction (coupled mode II,III cracks)




		//  and the linearization (Δtn.δg+Δtt.δξ  etc.)
		double Utn_α = -((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
		double Utn_β = -((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
		double Utt_α = (exp(-(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
		double Utt_β = (exp(-(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));

		double t_N = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)));
		double t_N_α = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_α);
		double t_N_β = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_β);
		double tracn = t_N_α + t_N_β;
		
		double t_T = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)));
		double t_T_α = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_α;
		double t_T_β = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_β;
		double tn_x = -(MBRACKET(-tracn));
		double tt_α = -(MBRACKET(-t_T_α));
		double tt_β = -(MBRACKET(-t_T_β));
		 t_N = -(MBRACKET(-t_N));
		 t_N_α = -(MBRACKET(-t_N_α));
		 t_N_β = -(MBRACKET(-t_N_β));
		 t_T = -(MBRACKET(-t_T));


		 double abs_α = ((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)) + 1) / (sqrt(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) + 1);
		 double abs_β = ((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)) + 1) / (sqrt(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) + 1);
		 
		 if (abs_α < 1){ abs_α = -1; }
		 
		 if (abs_β < 1){ abs_β = -1; }

		 


	     double DUtn_α = ((-2 * (((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*exp(-((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential));
		 double DUtn_β = ((-2 * (((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*exp(-((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential));

		 //double DUtt_α = (exp(-(((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*(abs_α / mCL_tangential) + ((-2 * (((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*(exp(-(((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
		 //double DUtt_β = (exp(-(((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*(abs_β / mCL_tangential) + ((-2 * (((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*(exp(-(((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));


		 double DUtt_α = (exp(-(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*(abs_α / mCL_tangential) + ((-2 * (((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*(exp(-(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[m][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
		 double DUtt_β = (exp(-(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*(abs_β / mCL_tangential) + ((-2 * (((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))))) / (mCL_tangential*mCL_tangential))*(exp(-(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[m][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));

		// the reason for having extra term(HEAVISIDE function)  in front of the equation is: (in the case of complete separation)TN--->0 , TT--->0  ===>  stiffness--->0       
		double dtn_dn = -HEAVYSIDE(((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal))))*(((max_normal / mCL_normal)*exp(1 - (gap_n / mCL_normal)) + (gap_n / mCL_normal) *  max_normal*(-1 / mCL_normal)*exp(1 - (gap_n / mCL_normal))));
		double dtt_dn = -HEAVYSIDE((2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal))))*(2 * sqrt(0.5*exp(1))*max_tangential)*(((1 / mCL_normal)*exp(-(gap_n / mCL_normal))) + ((1 + (gap_n / mCL_normal))*(-1 / mCL_normal)*exp(-(gap_n / mCL_normal))));


		
		for ( k = 0; k < ndof; ++k){
			for (l = 0; l < ndof; ++l) { ke[k][l] +=  (dtn_dn*(exp(Utn_α) + exp(Utn_β)))*NTE[k] * NTE[l] +  (dtt_dn*((D1[k] * Utt_α) + (D2[k] * Utt_β))* NTE[k]); }
		}
		
		double tntotal = tn_x;
	///   problem zone	 ss.D_R[m]
		
		// tn.ΔUtxn.δg+tt.ΔUtxn.δξ
		double ktrail=0;
		for (i = 0; i < ndof; ++i){
			for ( j = 0; j < ndof; ++j)
			{
				ktrail = t_N*DUtn_α*(Mk[0][0] * NTE[i] * D1[i] + Mk[0][1] * NTE[i] * D2[j]);
				ktrail += t_N*DUtn_β*(Mk[1][0] * NTE[i] * D1[j] + Mk[1][1] * NTE[i] * D2[j]);
				ktrail += t_T*DUtt_α*(Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j]);
				ktrail += t_T*DUtt_β*(Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j]);
			//kij = Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j] + Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j];
			//ke[i][j] += m_epsf*kij;
				ke[i][j] +=  ktrail;
			}
		}
		
		
		// ---    S T I F F N E S S ---
		 
		for ( k = 0; k < ndof; ++k){
			for ( l = 0; l < ndof; ++l)
			{
				sum = 0;
				sum = Mi[0][0] * Nb1[k] * Nb1[l] + Mi[0][1] * (Nb1[k] * Nb2[l] + Nb2[k] * Nb1[l]) + Mi[1][1] * Nb2[k] * Nb2[l];
				sum *= gap_n;
				sum -= D1[k] * N1[l] + D2[k] * N2[l] + N1[k] * D1[l] + N2[k] * D2[l];
				sum += K[0][1] * (D1[k] * D2[l] + D2[k] * D1[l]);
				//sum *= tn*m_knmult;
				//sum *= tntotal*m_knmult;
				sum *= tntotal;
				//sum += eps*HEAVYSIDE(Lm + eps*gap)*N[k] * N[l];
				//sum += dtn*N[k] * N[l];
				ke[k][l] +=  sum;
			}
			
		}
		

		//  m_ktmult is a user defined control function(default value m_ktmult = 1;)
		double kij;

		for ( i = 0; i < ndof; ++i){
			for ( j = 0; j < ndof; ++j)
			{
				// KT_α
				kij = T11[i] * D1[j] + T12[i] * D2[j];
				kij += D1[i] * T11[j] + D2[i] * T12[j];
				kij -= (Hrs[0][1] * tau[0])*D1[i] * D2[j] + (Hrs[1][0] * tau[0])*D2[i] * D1[j];
				kij += Tb11[i] * D1[j] + Tb21[i] * D2[j];
				kij += D1[i] * Tb11[j] + D2[i] * Tb21[j];
				kij += gap_n*(N11[i] * D1[j] + N12[i] * D2[j] + D1[i] * N11[j] + D2[i] * N12[j]);
				kij -= NTE[i] * Nb1[j] - Nb1[i] * NTE[j];
				kij -= T1[i] * Mi[0][0] * Tb11[j] + T1[i] * Mi[0][1] * Tb21[j] + T2[i] * Mi[1][0] * Tb11[j] + T2[i] * Mi[1][1] * Tb21[j];
				kij -= Tb11[i] * Mi[0][0] * T1[j] + Tb21[i] * Mi[0][1] * T1[j] + Tb11[i] * Mi[1][0] * T2[j] + Tb21[i] * Mi[1][1] * T2[j];

				//ke[i][j] += ((tt_α* Ac[0][0] + tt_β* Ac[1][0])*kij)* m_ktmult;
				ke[i][j] +=  ((tt_α* Ac[0][0] + tt_β* Ac[1][0])*kij);

				// KT_β
				kij = T21[i] * D1[j] + T22[i] * D2[j];
				kij += D1[i] * T21[j] + D2[i] * T22[j];
				kij -= (Hrs[0][1] * tau[1])*D1[i] * D2[j] + (Hrs[1][0] * tau[1])*D2[i] * D1[j];
				kij += Tb12[i] * D1[j] + Tb22[i] * D2[j];
				kij += D1[i] * Tb12[j] + D2[i] * Tb22[j];
				kij += gap_n*(N21[i] * D1[j] + N22[i] * D2[j] + D1[i] * N21[j] + D2[i] * N22[j]);
				kij -= NTE[i] * Nb2[j] - Nb2[i] * NTE[j];
				kij -= T1[i] * Mi[0][0] * Tb12[j] + T1[i] * Mi[0][1] * Tb22[j] + T2[i] * Mi[1][0] * Tb12[j] + T2[i] * Mi[1][1] * Tb22[j];
				kij -= Tb12[i] * Mi[0][0] * T1[j] + Tb22[i] * Mi[0][1] * T1[j] + Tb12[i] * Mi[1][0] * T2[j] + Tb22[i] * Mi[1][1] * T2[j];

				///ke[i][j] += ((tt_α* Ac[0][1] + tt_β* Ac[1][1])*kij)* m_ktmult;

				ke[i][j] +=  ((tt_α* Ac[0][1] + tt_β* Ac[1][1])*kij);


				}
			}


		}

	
}


//-----------------------------------------------------------------------------

//bool FEFacet2FacetSliding::Augment(int naug)
bool FEFacet2FacetSliding::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;

	int i;
	double Ln;
	double Lt[2];
	bool bconv = true;
	mat2d Mi;

	// penalty factor
	double eps, scale = Penalty();

	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	double normL0 = 0;
	for (i = 0; i<m_ss.Nodes(); ++i)	normL0 += m_ss.m_Lm[i] * m_ss.m_Lm[i];
	for (i = 0; i<m_ms.Nodes(); ++i)	normL0 += m_ms.m_Lm[i] * m_ms.m_Lm[i];

	// b. tangential component
	if (m_mu*m_epsf > 0)
	{
		for (i = 0; i<m_ss.Nodes(); ++i)
		{
			if (m_ss.m_pme[i])
			{
				Lt[0] = m_ss.m_Lt[i][0];
				Lt[1] = m_ss.m_Lt[i][1];
				mat2d& M = m_ss.m_M[i];
				Mi = M.inverse();
				normL0 += Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
			}
		}

		for (i = 0; i<m_ms.Nodes(); ++i)
		{
			if (m_ms.m_pme[i])
			{
				Lt[0] = m_ms.m_Lt[i][0];
				Lt[1] = m_ms.m_Lt[i][1];
				mat2d& M = m_ms.m_M[i];
				Mi = M.inverse();
				normL0 += Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
			}
		}
	}
	normL0 = sqrt(normL0);

	// --- c a l c u l a t e   c u r r e n t   n o r m s ---
	// a. normal component
	double normL1 = 0;	// force norm
	double normg1 = 0;	// gap norm
	int N = 0;
	for (i = 0; i<m_ss.Nodes(); ++i)
	{
		eps = m_ss.m_eps[i] * scale;

		// update Lagrange multipliers
		Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;

		if (m_ss.m_gap[i] > 0)
		{
			normg1 += m_ss.m_gap[i] * m_ss.m_gap[i];
			++N;
		}
	}

	for (i = 0; i<m_ms.Nodes(); ++i)
	{
		eps = m_ms.m_eps[i] * scale;

		// update Lagrange multipliers
		Ln = m_ms.m_Lm[i] + eps*m_ms.m_gap[i];
		Ln = MBRACKET(Ln);

		normL1 += Ln*Ln;
		if (m_ms.m_gap[i] > 0)
		{
			normg1 += m_ms.m_gap[i] * m_ms.m_gap[i];
			++N;
		}
	}
	if (N == 0) N = 1;

	// b. tangential component
	if (m_mu*m_epsf > 0)
	{
		double r, s, rp, sp;
		for (i = 0; i<m_ss.Nodes(); ++i)
		{
			if (m_ss.m_pme[i])
			{
				r = m_ss.m_rs[i][0];
				s = m_ss.m_rs[i][1];
				rp = m_ss.m_rsp_f[i][0];
				sp = m_ss.m_rsp_f[i][1];

				Ln = m_ss.m_Lm[i];

				mat2d& Mk = m_ss.m_M[i];
				Mi = Mk.inverse();

				Lt[0] = m_ss.m_Lt[i][0] + m_epsf*(Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
				Lt[1] = m_ss.m_Lt[i][1] + m_epsf*(Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));

				double TMT = Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0 && TMT> 0)
				{
					Lt[0] = m_mu*Ln*Lt[0] / sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1] / sqrt(TMT);
				}

				normL1 += Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
			}
		}

		for (i = 0; i<m_ms.Nodes(); ++i)
		{
			if (m_ms.m_pme[i])
			{
				r = m_ms.m_rs[i][0];
				s = m_ms.m_rs[i][1];
				rp = m_ms.m_rsp_f[i][0];
				sp = m_ms.m_rsp_f[i][1];
				Ln = m_ms.m_Lm[i];

				mat2d& Mk = m_ms.m_M[i];
				Mi = Mk.inverse();

				Lt[0] = m_ms.m_Lt[i][0] + m_epsf*(Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
				Lt[1] = m_ms.m_Lt[i][1] + m_epsf*(Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));

				double TMT = Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0 && TMT> 0)
				{
					Lt[0] = m_mu*Ln*Lt[0] / sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1] / sqrt(TMT);
				}

				normL1 += Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
			}
		}
	}

	normL1 = sqrt(normL1);
	if (N != 0){ normg1 = sqrt(normg1 / N); }

	if (naug == 0) m_normg0 = 0;

	// calculate and print convergence norms
	double lnorm = 0, gnorm = 0;
	if (normL1 != 0) lnorm = fabs(normL1 - normL0) / normL1; else lnorm = fabs(normL1 - normL0);
	if (normg1 != 0) gnorm = fabs(normg1 - m_normg0) / normg1; else gnorm = fabs(normg1 - m_normg0);

	felog.printf(" sliding interface # %d\n", m_nID);
	felog.printf("                        CURRENT        REQUIRED\n");
	felog.printf("    normal force : %15le", lnorm);
	if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n");
	felog.printf("    gap function : %15le", gnorm);
	if (m_gtol > 0) felog.printf("%15le\n", m_gtol); else felog.printf("       ***\n");

	// check convergence
	bconv = true;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	// to do need to be changed  gnorm > m_gtol!!
	if ((m_gtol > 0) && (gnorm > m_gtol)) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if (m_naugmax <= naug) bconv = true;

	if (bconv == false)
	{
		// we did not converge so update multipliers
		for (i = 0; i<m_ss.Nodes(); ++i)
		{
			eps = m_ss.m_eps[i] * scale;

			// update Lagrange multipliers
			Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
			m_ss.m_Lm[i] = MBRACKET(Ln);

			if ((m_mu*m_epsf > 0) && (m_ss.m_pme[i]))
			{
				// update the metrics
				FESurfaceElement& mel = *m_ss.m_pme[i];

				double r = m_ss.m_rs[i][0], s = m_ss.m_rs[i][1];
				double rp = m_ss.m_rsp_f[i][0], sp = m_ss.m_rsp_f[i][1];

				 Ln = m_ss.m_Lm[i];

				mat2d Mk = m_ss.m_M[i];
				mat2d Mki = Mk.inverse();


				// update traction multipliers
				// a. trial state
				Lt[0] = m_ss.m_Lt[i][0] + m_epsf*(Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
				Lt[1] = m_ss.m_Lt[i][1] + m_epsf*(Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));

				double TMT = Lt[0] * (Mki[0][0] * Lt[0] + Mki[0][1] * Lt[1]) + Lt[1] * (Mki[1][0] * Lt[0] + Mki[1][1] * Lt[1]);
				//assert(TMT >= 0);

				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0 && TMT> 0)
				{
					Lt[0] = m_mu*Ln*Lt[0] / sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1] / sqrt(TMT);
				}

				m_ss.m_M[i] = m_ss.Metric0(mel, r, s);

				m_ss.m_Lt_P[i][0] = m_ss.m_Lt[i][0];
				m_ss.m_Lt_P[i][0] = m_ss.m_Lt[i][0];
				m_ss.m_Lt[i][0] = Lt[0];
				m_ss.m_Lt[i][1] = Lt[1];
			}
		}

		for (i = 0; i<m_ms.Nodes(); ++i)
		{
			eps = m_ms.m_eps[i] * scale;

			// update Lagrange multipliers
			Ln = m_ms.m_Lm[i] + eps*m_ms.m_gap[i];
			m_ms.m_Lm[i] = MBRACKET(Ln);

			if ((m_mu*m_epsf > 0) && (m_ms.m_pme[i]))
			{
				// update the metrics
				FESurfaceElement& mel = *m_ms.m_pme[i];

				double r = m_ms.m_rs[i][0], s = m_ms.m_rs[i][1];
				double rp = m_ms.m_rsp_f[i][0], sp = m_ms.m_rsp_f[i][1];

				 Ln = m_ms.m_Lm[i];

				mat2d Mk = m_ms.m_M[i];
				mat2d Mki = Mk.inverse();

				// update traction multipliers
				// a. trial state
				//double Lt[2];
				Lt[0] = m_ms.m_Lt[i][0] + m_epsf*(Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
				Lt[1] = m_ms.m_Lt[i][1] + m_epsf*(Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));

				double TMT = Lt[0] * (Mki[0][0] * Lt[0] + Mki[0][1] * Lt[1]) + Lt[1] * (Mki[1][0] * Lt[0] + Mki[1][1] * Lt[1]);
				//assert(TMT >= 0);

				double phi = sqrt(TMT) - m_mu*Ln;

				// b. return map
				if (phi > 0 && TMT> 0)
				{
					Lt[0] = m_mu*Ln*Lt[0] / sqrt(TMT);
					Lt[1] = m_mu*Ln*Lt[1] / sqrt(TMT);
				}

				m_ms.m_M[i] = m_ms.Metric0(mel, r, s);

				m_ms.m_Lt_P[i][0] = m_ms.m_Lt[i][0];
				m_ms.m_Lt_P[i][0] = m_ms.m_Lt[i][0];
				m_ms.m_Lt[i][0] = Lt[0];
				m_ms.m_Lt[i][1] = Lt[1];
			}
		}
	}

	if (bconv)
	{
		//m_ss.m_rsp = m_ss.m_rs;
		//m_ms.m_rsp = m_ms.m_rs;
		m_ss.m_rsp_f = m_ss.m_rs;
		m_ms.m_rsp_f = m_ms.m_rs;

	}

	// store the last gap norm
	m_normg0 = normg1;

	return bconv;
}



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//! This function transforms friction data between two master segments

void FEFacet2FacetSliding::MapFrictionData(int inode, FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, FESurfaceElement &en, FESurfaceElement &eo, vec3d &q)
{
	// first we find the projection of the old point on the new segment
	double r = ss.m_rs[inode][0];
	double s = ss.m_rs[inode][1];
	double rp = ss.m_rsp_f[inode][0], ro = rp;
	double sp = ss.m_rsp_f[inode][1], so = sp;
	vec3d xn = ms.Local2Global(eo, rp, sp);
	vec3d qn;
	qn = ms.ProjectToSurface(en, xn, rp, sp);
	ss.m_rsp_f[inode][0] = rp;
	ss.m_rsp_f[inode][1] = sp;

	// next, we transform the frictional traction
	// since these tractions are decomposed in the local 
	// element coordinate system, we have to do a coordinate transformation
	// note that this transformation needs to be done in curvilinear
	// coordinates since the base vectors may not be orthonormal. Also
	// note that we are doing this in the reference configuration
	vec3d to[2], tn[2];
	ms.ContraBaseVectors0(eo, ro, so, to);
	ms.CoBaseVectors0(en, r, s, tn);

	double Lt[2];
	Lt[0] = ss.m_Lt[inode][0];
	Lt[1] = ss.m_Lt[inode][1];

	vec3d t;
	t = to[0] * Lt[0] + to[1] * Lt[1];

	Lt[0] = t*tn[0];
	Lt[1] = t*tn[1];
	//!< Lagrange multipliers for tangential(previous time step)
	
	ss.m_Lt_P[inode][0] = ss.m_Lt[inode][0];
	ss.m_Lt_P[inode][0] = ss.m_Lt[inode][0];
	ss.m_Lt[inode][0] = Lt[0];
	ss.m_Lt[inode][1] = Lt[1];
}
////////////////////
/////////////////////
void FEFacet2FacetSliding::MapTangentialComponent(int inode, FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, FESurfaceElement &en, FESurfaceElement &eo, vec3d &q)
{
	// first we find the projection of the old point on the new segment
	double r = ss.m_rs[inode][0];
	double s = ss.m_rs[inode][1];
	double rp = ss.m_rsp[inode][0], ro = rp;
	double sp = ss.m_rsp[inode][1], so = sp;
	vec3d xn = ms.Local2Global(eo, rp, sp);
	vec3d qn;
	qn = ms.ProjectToSurface(en, xn, rp, sp);
	ss.m_rsp[inode][0] = rp;
	ss.m_rsp[inode][1] = sp;

	// next, we transform the frictional traction
	// since these tractions are decomposed in the local 
	// element coordinate system, we have to do a coordinate transformation
	// note that this transformation needs to be done in curvilinear
	// coordinates since the base vectors may not be orthonormal. Also
	// note that we are doing this in the reference configuration
	vec3d to[2], tn[2];
	ms.ContraBaseVectors0(eo, ro, so, to);
	ms.CoBaseVectors0(en, r, s, tn);

	double Lt[2];
	Lt[0] = ss.m_control[inode][0];
	Lt[1] = ss.m_control[inode][1];
	
	vec3d t;
	t = to[0] * Lt[0] + to[1] * Lt[1];

	Lt[0] = t*tn[0];
	Lt[1] = t*tn[1];
	//!< Lagrange multipliers for tangential(previous time step)
	//ss.m_Lt_P[inode][0] = ss.m_Lt[inode][0];
	//ss.m_Lt_P[inode][0] = ss.m_Lt[inode][0];
	//ss.m_Lt[inode][0] = Lt[0];
	//ss.m_Lt[inode][1] = Lt[1];

	//    ss.m_trans_S should be used later
	ss.m_trans_S[inode][0] = Lt[0];
	ss.m_trans_S[inode][1] = Lt[1];
}
//////////////////
//-----------------------------------------------------------------------------

void FEFacet2FacetSliding::UpdateContactPressures()
{
	int npass = (m_btwo_pass ? 2 : 1);
	
	bool temp = false;
	//double max_normal = M_TN* ss.D_R[m];
	//double max_tangential = M_Tt* ss.D_R[m];
	

	
	for (int np = 0; np<npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0 ? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0 ? m_ms : m_ss);

		// loop over all nodes of the primary surface
		for (int n = 0; n<ss.Nodes(); ++n)
		{
			// get the normal tractions at the integration points
			double gapP = ss.m_gap_p[n];
			double gap = ss.m_gap[n];
			double eps = m_eps*ss.m_eps[n];
			double max_normal = M_TN* ss.D_R[n];
			//felog.printf("max_normal-------------------        %15le\n", max_normal);

			double max_tangential = M_Tt* ss.D_R[n];
			double mCL_normal = δ_TN * ss.D_R_C[n];

			//felog.printf("mCL_normal-------------------        %15le\n", mCL_normal);

			double mCL_tangential = δ_Tt * ss.D_R_C[n];
			//double eps = m_eps;
			////////////////////

			//////
			//gap = pt.m_gap;
			//eps = m_epsn*pt.m_eps;
			//if (eps == 0){
				//felog.printf("Division by zero!....................penalty value==0");
				//break;
			//}
			if (gap > 0){
				//double gapeps = ss.m_eps[n] * m_gtol;
				//tn = Lm + eps*gap;
				
				//ss.m_Lnp[n] = MBRACKET(ss.m_Lmtemp[n] + eps*gap);
				ss.m_Lnp[n] = MBRACKET(ss.m_Lm[n] + eps*gap);
			}
			/////////////////////////////////
			if (gap <= 0){
				//test
				

				//////////////////
				mat2d& Mk = ss.m_M[n];
				mat2d Mki = Mk.inverse();

				// get the master element node positions
				//for (k = 0; k<nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;
				
				// isoparametric coordinates of the projected slave node
				// onto the master element
				double r = ss.m_rs[n][0];
				double s = ss.m_rs[n][1];

				// get the coordinates at the previous step
				double rp = ss.m_rsp[n][0];
				double sp = ss.m_rsp[n][1];

				ss.m_tan_gap_P[n][0] = ss.m_tan_gap[n][0];
				ss.m_tan_gap_P[n][1] = ss.m_tan_gap[n][1];
				ss.m_control[n][0] = ss.m_trans_S[n][0] + (Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
				ss.m_control[n][1] = ss.m_trans_S[n][1] + (Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));
				//eps = mCL_normal = mCL_tangential;

				
				double gap_n = fabs(gap);
				double Utn_α = -((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
				double Utn_β = -((ss.m_trans_S[n][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[n][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
				double Utt_α = (exp(-(((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
				double Utt_β = (exp(-(((ss.m_trans_S[n][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[n][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ss.m_trans_S[n][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ss.m_trans_S[n][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));


				double t_N_α = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_α);
				double t_N_β = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_β);
				double tracn = t_N_α + t_N_β;
				
				// contact traction( cohesive case! mode II and III)
				double t_T_α = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_α;
				double t_T_β = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_β;
				//double tract = t_T_α + t_T_β;
				double tn_x = -(MBRACKET(-tracn));
				//double tn_x = -(MBRACKET(-tnX));
				//double tt_α = -(MBRACKET(-t_T_α));
				//double tt_β = -(MBRACKET(-t_T_β));

				//double tn_x = -(MBRACKET(-tracn));

				
				
				ss.m_tan_gap[n][0] = ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp);
				ss.m_tan_gap[n][1] = ss.m_trans_S[n][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp);

				double tractt1 = t_T_α;
				double tractt2 = t_T_β;
				//double tractt = tractt1 + tractt2;

				//!< (previous time step)
				ss.m_pt_P[n][0] = ss.m_pt[n][0];
				ss.m_pt_P[n][1] = ss.m_pt[n][1];
				ss.m_Ln_P[n] = ss.m_Ln[n];

				//pt.m_Ln = -(MBRACKET(-tracs));
				ss.m_pt[n][0] =  (-(MBRACKET(-tractt1)));
				ss.m_pt[n][1] =  (-(MBRACKET(-tractt2)));
				ss.m_Ln[n] =  tn_x;

				////////////////////
				//////////////////^^^^^^^^^^^^^^^^^^^^^^^^
				if (m_mu*m_epsf > 0) {

					double Tt_0 = ss.m_trans_S[n][0] + m_epsf*(Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
					double Tt_1 = ss.m_trans_S[n][1] + m_epsf*(Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));

					//if (isinf(Tt[0]) != 0){ Tt[0] = 0; }
					//if (isinf(Tt[1]) != 0){ Tt[1] = 0; }

					double TMT = Tt_0 * (Mki[0][0] * Tt_0 + Mki[0][1] * Tt_1) + Tt_1 * (Mki[1][0] * Tt_0 + Mki[1][1] * Tt_1);
					//assert(TMT >= 0);


					// b. return map

					double LNX = fabs(tn_x);
					//double phi = sqrt(TMT) - m_mu*Ln;
					double phi = 0;

					if (TMT >= 0){
						phi = sqrt(TMT) - m_mu*LNX;
					}
					// b. return map
					if ((phi > 0) && (TMT > 0))
					{
						//Tt[0] = m_mu*Ln*Tt[0] / sqrt(TMT);
						//Tt[1] = m_mu*Ln*Tt[1] / sqrt(TMT);
						Tt_0 = m_mu*LNX*Tt_0 / sqrt(TMT);
						Tt_1 = m_mu*LNX*Tt_1 / sqrt(TMT);
					}
					ss.m_Lt_P_C[n][0] = ss.m_Lt_C[n][0];
					ss.m_Lt_P_C[n][1] = ss.m_Lt_C[n][1];
						ss.m_Lt_C[n][0] = Tt_0;
						ss.m_Lt_C[n][1] = Tt_1;

				}
				/////////////////^^^^^^^^^^^^^^^^^^^^^^^^^^^


				double gabTX = sqrt(((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[n][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))));
				double gabTy = sqrt(((ss.m_trans_S[n][1] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ss.m_trans_S[n][1] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))));
				double gabTxp = ss.m_tan_gap_P[n][0];
				double gabTyp = ss.m_tan_gap_P[n][1];



				if ((sqrt((gap_n*gap_n) + β*(gabTX*gabTX) + β*(gabTy*gabTy))>sqrt((gapP*gapP) + β*(gabTxp*gabTxp) + β*(gabTyp*gabTyp)))){
					if (sqrt((gap_n*gap_n) + (gabTX*gabTX) + (gabTy*gabTy)) > sqrt((mCL_normal*mCL_normal) + (mCL_tangential*mCL_tangential))){
						
						//ss.DTCM[n] = sqrt(mCL_normal*mCL_normal + β*mCL_tangential*mCL_tangential);
						double tnX = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_α)*exp(Utn_β);
						//ss.DTCM[n] = sqrt((gap_n*gap_n) + β*(gabTX*gabTX) + β*(gabTy*gabTy));//


						if ((mCL_normal != 0) && (mCL_tangential != 0)) {
							ss.DTCM[n] = (sqrt((gap_n*gap_n) / (mCL_normal*mCL_normal) + β*(gabTX*gabTX) / (mCL_tangential*mCL_tangential) + β*(gabTy*gabTy) / (mCL_tangential*mCL_tangential)));
						}
						//ss.DTCM[n] = (α - 1)*ss.DTCM[n];
						//felog.printf("Zero/cyclic stresses\n");
						//ss.DT[n] = sqrt((tnX*tnX) / (max_normal*max_normal) + β*(t_T_α*t_T_α) / (max_tangential*max_tangential) + β*(t_T_β*t_T_β) / (max_tangential*max_tangential));
						//ss.DT[n] = sqrt((tn_x*tn_x)  + β*(t_T_α*t_T_α)  + β*(t_T_β*t_T_β) );//
						if ((max_normal != 0) && (max_tangential != 0)) {
							ss.DT[n] = sqrt((tnX*tnX) / (max_normal*max_normal) + β*(t_T_α*t_T_α) / (max_tangential*max_tangential) + β*(t_T_β*t_T_β) / (max_tangential*max_tangential));
						}
						//felog.printf("ss.DT[n]-------------------        %15le\n", ss.DT[n]);
						ss.DAM[n] = true;
						//felog.printf("ss.DTCM[n]-------------------      %15le\n", ss.DTCM[n]);

					}
				}

				//if ((ss.DAM[n]) && ((sqrt((tn_x*tn_x) / (max_normal*max_normal) + (t_T_α*t_T_α) / (max_tangential*max_tangential) + (t_T_β*t_T_β) / (max_tangential*max_tangential)) == 0))){
				if ((ss.DAM[n]) && ((sqrt((gap_n*gap_n) / (mCL_normal*mCL_normal) + β*(gabTX*gabTX) / (mCL_tangential*mCL_tangential) + β*(gabTy*gabTy)/(mCL_tangential*mCL_tangential))) < 1)){

					//felog.printf("Zero/cyclic stresses\n");


					//DR = DT;
					
					//ss.D_R[n] = 1 - ss.DT[n];
					//ss.D_R[n] = ss.DT[n] / (max_normal + β*max_tangential);//
					ss.D_R[n] = ss.D_R[n]*ss.DT[n];
					//ss.D_R_C[n] = (1 - ss.DT[n])*ss.DTCM[n];
					//ss.D_R_C[n] = ss.DTCM[n] / sqrt(mCL_normal*mCL_normal + β*mCL_tangential*mCL_tangential);//
					ss.D_R_C[n] = ss.D_R_C[n]*ss.DTCM[n];
					ss.DAM[n] = false;

				}




				//////////////////////////////
				//ss.m_Ln_max[n] = δ_TN * ss.D_R_C[n];
				
			}
			//////
			//if (gap < 0){ ss.m_Ln[n] = MBRACKET(ss.m_Lm[n] + eps*gap); }
			FESurfaceElement* pme = ss.m_pme[n];
			if (m_btwo_pass && pme)
				
			{

				int me = pme->Nodes();
				//if (gap > 0){
					if (me < 6)
					{
						double ti[6] = { 0 }, tii[6] = { 0 };
						double  t_N_α[6] = { 0 }, t_N_β[6] = { 0 }, t_T_α[6] = { 0 }, t_T_β[6] = { 0 }, tTα[6] = { 0 }, tTβ[6] = { 0 };
						double Tt_0_m = 0, Tt_1_m = 0;
						for (int j = 0; j < me; ++j) {
							int k = pme->m_lnode[j];
							gap = ms.m_gap[k];

							eps = m_eps*ms.m_eps[k];

							if (gap > 0){
								ti[j] = MBRACKET(ms.m_Lm[k] + m_eps*ms.m_eps[k] * ms.m_gap[k]);

							}

							if (gap <= 0){

								mat2d& Mk = ms.m_M[k];
								mat2d Mki = Mk.inverse();

								double r = ms.m_rs[k][0];
								double s = ms.m_rs[k][1];

								// get the coordinates at the previous step
								double rp = ms.m_rsp[k][0];
								double sp = ms.m_rsp[k][1];

								ms.m_tan_gap_P[k][0] = ms.m_tan_gap[k][0];
								ms.m_tan_gap_P[k][1] = ms.m_tan_gap[k][1];
								ms.m_control[k][0] = ms.m_trans_S[k][0] + (Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
								ms.m_control[k][1] = ms.m_trans_S[k][1] + (Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));
								ms.m_tan_gap[k][0] = ms.m_trans_S[k][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp);
								ms.m_tan_gap[k][1] = ms.m_trans_S[k][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp);


								//eps = mCL_normal = mCL_tangential;
								double gap_n = fabs(gap);
								double Utn_α = -((ms.m_trans_S[k][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ms.m_trans_S[k][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
								double Utn_β = -((ms.m_trans_S[k][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ms.m_trans_S[k][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp))) / (mCL_tangential*mCL_tangential);
								double Utt_α = (exp(-(((ms.m_trans_S[k][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ms.m_trans_S[k][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ms.m_trans_S[k][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))*((ms.m_trans_S[k][0] + Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp)))) / mCL_tangential));
								double Utt_β = (exp(-(((ms.m_trans_S[k][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ms.m_trans_S[k][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / (mCL_tangential*mCL_tangential)))*((sqrt(((ms.m_trans_S[k][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))*((ms.m_trans_S[k][1] + Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp)))) / mCL_tangential));


								t_N_α[j] = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_α);
								t_N_β[j] = -((gap_n / mCL_normal) * max_normal * exp(1 - (gap_n / mCL_normal)))*exp(Utn_β);
								// contact traction( cohesive case! mode II and III)
								t_T_α[j] = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_α;
								t_T_β[j] = -(2 * sqrt(0.5*exp(1))*max_tangential *(1 + (gap_n / mCL_normal))*exp(-(gap_n / mCL_normal)))* Utt_β;
								//double tract = t_T_α + t_T_β;

								tTα[j] = -(MBRACKET(-t_T_α[j]));
								tTβ[j] = -(MBRACKET(-t_T_β[j]));

								tii[j] = -(MBRACKET(-(t_N_α[j] + t_N_β[j])));

								//////////////////^^^^^^^^^^^^^^^^^^^^^^^^
								if (m_mu*m_epsf > 0) {

									 Tt_0_m = ms.m_trans_S[k][0] + m_epsf*(Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
									 Tt_1_m = ms.m_trans_S[k][1] + m_epsf*(Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));

									//if (isinf(Tt[0]) != 0){ Tt[0] = 0; }
									//if (isinf(Tt[1]) != 0){ Tt[1] = 0; }

									 double TMT = Tt_0_m * (Mki[0][0] * Tt_0_m + Mki[0][1] * Tt_1_m) + Tt_1_m * (Mki[1][0] * Tt_0_m + Mki[1][1] * Tt_1_m);
									//assert(TMT >= 0);


									// b. return map

									double LNX = fabs(tii[j]);
									//double phi = sqrt(TMT) - m_mu*Ln;
									double phi = 0;

									if (TMT >= 0){
										phi = sqrt(TMT) - m_mu*LNX;
									}
									// b. return map
									if ((phi > 0) && (TMT > 0))
									{
										//Tt[0] = m_mu*Ln*Tt[0] / sqrt(TMT);
										//Tt[1] = m_mu*Ln*Tt[1] / sqrt(TMT);
										Tt_0_m = m_mu*LNX*Tt_0_m / sqrt(TMT);
										Tt_1_m = m_mu*LNX*Tt_1_m / sqrt(TMT);
									}
									ms.m_Lt_P_C[k][0] = ms.m_Lt_C[k][0];
									ms.m_Lt_P_C[k][1] = ms.m_Lt_C[k][1];
									ms.m_Lt_C[k][0] = Tt_0_m;
									ms.m_Lt_C[k][1] = Tt_1_m;

								}
								/////////////////^^^^^^^^^^^^^^^^^^^^^^^^^^^


							}
						}
						double tnn[6], ttα[6], ttβ[6];
						pme->project_to_nodes(tii, tnn);
						pme->project_to_nodes(tTα, ttα);
						pme->project_to_nodes(tTβ, ttβ);
						// now evaluate the traction at the intersection point
						double Lnn = pme->eval(tnn, ss.m_rs[n][0], ss.m_rs[n][1]);
						double Lnttα = pme->eval(ttα, ss.m_rs[n][0], ss.m_rs[n][1]);
						double Lnttβ = pme->eval(ttβ, ss.m_rs[n][0], ss.m_rs[n][1]);
						ss.m_Ln[n] += MBRACKET(Lnn);
						ss.m_pt[n][0] += MBRACKET(Lnttα);
						ss.m_pt[n][1] += MBRACKET(Lnttβ);


						// project the data to the nodes
						double tn[6];
						pme->project_to_nodes(ti, tn);
						// now evaluate the traction at the intersection point
						double Ln = pme->eval(tn, ss.m_rs[n][0], ss.m_rs[n][1]);
						ss.m_Lnp[n] += MBRACKET(Ln);
					}
				



			}

		}
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ShallowCopy(DumpStream& dmp, bool bsave)
{
	m_ss.ShallowCopy(dmp, bsave);
	m_ms.ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Serialize(DumpFile& ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);
}


