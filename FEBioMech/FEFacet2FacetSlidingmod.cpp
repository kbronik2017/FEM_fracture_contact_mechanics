#include "stdafx.h"
#include "FEFacet2FacetSliding.h"
#include "FEStiffnessMatrix.h"
#include "FECore/FEModel.h"
#include "FECore/FEClosestPointProjection.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_PARAMETER_LIST(FEFacet2FacetSliding, FEContactInterface)
	ADD_PARAMETER(m_epsn     , FE_PARAM_DOUBLE, "penalty"      );
	ADD_PARAMETER(m_bautopen , FE_PARAM_BOOL  , "auto_penalty" );
	ADD_PARAMETER(m_blaugon  , FE_PARAM_BOOL  , "laugon"       );
	ADD_PARAMETER(m_atol     , FE_PARAM_DOUBLE, "tolerance"    );
	ADD_PARAMETER(m_btwo_pass, FE_PARAM_BOOL  , "two_pass"     );
	ADD_PARAMETER(m_gtol     , FE_PARAM_DOUBLE, "gaptol"       );
	ADD_PARAMETER(m_naugmin  , FE_PARAM_INT   , "minaug"       );
	ADD_PARAMETER(m_naugmax  , FE_PARAM_INT   , "maxaug"       );
	ADD_PARAMETER(m_knmult   , FE_PARAM_DOUBLE, "knmult"       );
	ADD_PARAMETER(m_stol     , FE_PARAM_DOUBLE, "search_tol"   );
	ADD_PARAMETER(m_srad     , FE_PARAM_DOUBLE, "search_radius");
	ADD_PARAMETER(m_dxtol    , FE_PARAM_DOUBLE, "dxtol"        );
	ADD_PARAMETER(m_mu       , FE_PARAM_DOUBLE, "fric_coeff"   );
	ADD_PARAMETER(m_epsf     , FE_PARAM_DOUBLE, "fric_penalty" );
	ADD_PARAMETER(m_nsegup   , FE_PARAM_INT   , "seg_up"       );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//int FEFacet2FacetSliding::maxcohetraction = 20;

 double FEFacet2FacetSliding::rootcontrol=1.000000;
FEFacetSlidingSurface::Data::Data()
{
	m_gap = 0.0;
	m_Lm  = 0.0;
	m_Ln  = 0.0;
	m_eps = 1.0;
	//int nn = m_ms.;
	m_nu = vec3d(0,0,0);
	m_rs = vec2d(0,0);
	m_pme = (FESurfaceElement*)0;
	// get the number of nodes

	//m_M.resize(nn);
}

//-----------------------------------------------------------------------------
// FEFacetSlidingSurface
//-----------------------------------------------------------------------------

bool FEFacetSlidingSurface::Init()
{
	// initialize surface data first
	if (FEContactSurface::Init() == false) return false;
	int nn =Nodes();
	//int mk=
	//int em = Elements();
	for (int n = 0; n<Elements(); ++n)
	{
		//m_Lnt[n] = 1;
		FESurfaceElement& el = Element(n);
		//FESurfaceElement& se = ss.Element(i);
		int nint = el.GaussPoints();
		m_Lnt[el.m_nID].resize(nint);
		m_M_t[el.m_nID].resize(nint);
		m_rs_n[el.m_nID].resize(nint);
		m_rsp_n[el.m_nID].resize(nint);
		m_Lnt[el.m_nID].assign(nint, vec2d(0, 0));
		//m_Lnt[n].assign(nint, vec2d(0, 0));
		//m_Lnt[n].assign(nint, vec2d(0, 0));

	}
	//m_Lnt.resize(nn);
	//m_M_t.resize(nn);
	//m_rs_n.resize(nn);
	// natural coords of projected slave node on master element
	//m_rsp_n.resize(nn);
	
	//m_Lnt.assign
	// allocate data structures
	const int NE = Elements();
	m_Data.resize(NE);
	for (int i=0; i<NE; ++i)
	{
		FESurfaceElement& el = Element(i);
		int nint = el.GaussPoints();
		m_Data[i].resize(nint);
	}

	return true;
}

//-----------------------------------------------------------------------------
vec3d FEFacetSlidingSurface::GetContactForce()
{
	// initialize contact force
	vec3d f(0,0,0);
	const int MN = FEElement::MAX_NODES;
	double Tn[MN],T1[MN],T2[MN];
	vec3d e[2];
	// loop over all elements of the primary surface
	for (int n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nint = el.GaussPoints();
		
		// evaluate the contact force for that element
		for (int i=0; i<nint; ++i) 
		{
			Data& pt = m_Data[n][i];

			T1[i] = m_Lnt[el.m_nID][i][0];
			T2[i] = m_Lnt[el.m_nID][i][1];
			// unit vector
			vec3d n = SurfaceNormal(el, i);
			// gauss weight
			double w = el.GaussWeights()[i];
			// area in reference configuration
			vec3d g0[2], g[2];
			double r = el.gr(i);
			double s = el.gs(i);
			CoBaseVectors0(el, r, s, g0);
			double A = (g0[0] ^ g0[1]).unit();

			double rx = m_rs_n[el.m_nID][i][0];
			double sx = m_rs_n[el.m_nID][i][1];
			ContraBaseVectors(el, rx, sx, e);
			////////////////

			double t1 = el.eval(T1, i);
			double t2 = el.eval(T2, i);
			//double t3 = el.eval(Tn, i);
			// unit normal vector
			//vec3d n = SurfaceNormal(el, i);
			// contravariant basis in spatial frame
			//ContraBaseVectors(el, r, s, g);
			// Piola traction
			//double t = e[0] * t1 + e[1] * t2;
			vec3d t = e[0] * t1 + e[1] * t2;

			// contact traction( cohesive case!)
			//double tn = -((g / eps) * maxcohetraction * exp(1 - (g / eps)));
			//tn = MBRACKET(-tn);
			////////////////////////////////
			/*
			for (int k = 0; k<nseln; ++k)
			{
				fe[3 * k] = Hs[k] * nu.x;
				fe[3 * k + 1] = Hs[k] * nu.y;
				fe[3 * k + 2] = Hs[k] * nu.z;
			}

			for (int k = 0; k<nmeln; ++k)
			{
				fe[3 * (k + nseln)] = -Hm[k] * nu.x;
				fe[3 * (k + nseln) + 1] = -Hm[k] * nu.y;
				fe[3 * (k + nseln) + 2] = -Hm[k] * nu.z;
			}

			*/
			//////////////////////////////
			

			// contact force
			f += n*(w*pt.m_Ln*A) + t*(w*A);
		}
	}
	felog.printf("smart part just test  GetContactForce");
	felog.printf("f.x,f.y,f.z are:  %u \n", f.x,f.y,f.z);
	return f;
}




//-----------------------------------------------------------------------------
double FEFacetSlidingSurface::GetContactArea()
{
	// initialize contact area
	double a = 0;
	
	// loop over all elements of the primary surface
	for (int n=0; n<Elements(); ++n)
	{
		FESurfaceElement& el = Element(n);
		int nint = el.GaussPoints();
		
		// evaluate the contact force for that element
		for (int i=0; i<nint; ++i)
		{
			// get data for this integration point
			Data& data = m_Data[n][i];
            double s = (data.m_Ln > 0) ? 1 : 0;
            
			// get the base vectors
			vec3d g[2];
			CoBaseVectors(el, i, g);
            
			// normal (magnitude = area)
			vec3d n = g[0] ^ g[1];
            
			// gauss weight
			double w = el.GaussWeights()[i];
            
			// contact force
			a += n.norm()*(w*s);
		}
	}
	
	return a;
}

//-----------------------------------------------------------------------------
//! \todo Originally, we only copied Lmd, gap, Ln and reset pme to zero.
//!       Need to check if this achieves the same
void FEFacetSlidingSurface::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				dmp << d.m_gap;
				dmp << d.m_nu;
				dmp << d.m_rs;
				dmp << d.m_Lm;
				dmp << d.m_eps;
				dmp << d.m_Ln;
			}
		}
	}
	else
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				dmp >> d.m_gap;
				dmp >> d.m_nu;
				dmp >> d.m_rs;
				dmp >> d.m_Lm;
				dmp >> d.m_eps;
				dmp >> d.m_Ln;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::Serialize(DumpFile& ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsSaving())
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				ar << d.m_gap;
				ar << d.m_nu;
				ar << d.m_rs;
				ar << d.m_Lm;
				ar << d.m_eps;
				ar << d.m_Ln;
			}
		}
	}
	else
	{
		for (int i=0; i<(int) m_Data.size(); ++i)
		{
			vector<Data>& di = m_Data[i];
			int nint = (int) di.size();
			for (int j=0; j<nint; ++j)
			{
				Data& d = di[j];
				ar >> d.m_gap;
				ar >> d.m_nu;
				ar >> d.m_rs;
				ar >> d.m_Lm;
				ar >> d.m_eps;
				ar >> d.m_Ln;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::GetNodalContactGap(int nface, double* gn)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	double gi[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k) gi[k] = m_Data[nface][k].m_gap;

	//for (int k=0; k<ni; ++k) if (gi[k] < 0) gi[k] = 0;
	el.project_to_nodes(gi, gn);

	//for (int k=0; k<ne; ++k) if (gn[k] < 0) gn[k] = 0;
	for (int k = 0; k<ne; ++k) if (gn[k] < 0) gn[k] = -gn[k];
	felog.printf("smart part just test  GetNodalContactGap");
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::GetNodalContactPressure(int nface, double* pn)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	double ti[FEElement::MAX_INTPOINTS];
	for (int k=0; k<ni; ++k)
	{
		double L = m_Data[nface][k].m_Ln;
		double t1= m_Lnt[el.m_nID][k][0];
		double t2 =m_Lnt[el.m_nID][k][1];
		ti[k] = L + t1 + t2;// + pf->m_epsn*gi[k];
		//ti[k] = (ti[k]>=0?ti[k] : 0);		
	}

	el.project_to_nodes(ti, pn);
	//for (int k=0; k<ni; ++k)
		//pn[k] = (pn[k]>=0?pn[k] : 0);	
	felog.printf("smart part just test  GetNodalContactPressure");
}

//-----------------------------------------------------------------------------
void FEFacetSlidingSurface::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& el = Element(nface);
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	const int MN = FEElement::MAX_NODES;
	double  t1, t2;
	vec3d t,tt;
	vec3d e[2];
	const int MFI = FEElement::MAX_INTPOINTS;
	double tix[MFI], tiy[MFI], tiz[MFI];
	for (int k=0; k<ni; ++k)
	{
		FEFacetSlidingSurface::Data& pt = m_Data[nface][k];
		double gi = pt.m_gap;
		double Li = pt.m_Ln;
		vec3d  ti = pt.m_nu;
		t1= m_Lnt[el.m_nID][k][0];
		t2= m_Lnt[el.m_nID][k][1];
		
		double rx = m_rs_n[el.m_nID][k][0];
		double sx = m_rs_n[el.m_nID][k][1];
		ContraBaseVectors(el, rx, sx, e);
		// Piola traction
		tt = e[0] * t1 + e[1] * t2;
		t = ti*(Li)+tt;
		//if (gi > 0) t = ti*(Li); else t = vec3d(0,0,0);
		tix[k] = t.x; tiy[k] = t.y; tiz[k] = t.z;
	}

	// project traction to nodes
	const int MFN = FEElement::MAX_NODES;
	double tnx[MFN], tny[MFN], tnz[MFN];
	el.project_to_nodes(tix, tnx);
	el.project_to_nodes(tiy, tny);
	el.project_to_nodes(tiz, tnz);

	// store data
	for (int k=0; k<ne; ++k)
	{
		tn[k].x = tnx[k];
		tn[k].y = tny[k];
		tn[k].z = tnz[k];
	}
	felog.printf("smart part just test  GetNodalContactTraction");
}

//-----------------------------------------------------------------------------
// FEFacet2FacetSliding
//-----------------------------------------------------------------------------




FEFacet2FacetSliding::FEFacet2FacetSliding(FEModel* pfem) : FEContactInterface(pfem), m_ss(&pfem->GetMesh()), m_ms(&pfem->GetMesh())
{
	static int ncount = 1;
	m_nID = ncount++;
	//int maxcohetraction = m_naugmin;
	// default parameters
	m_epsn = 1.0;
	m_knmult = 1.0;
	m_ktmult = 1.0;
	m_stol = 0.01;
	m_btwo_pass = false;
	m_bautopen = false;
	m_nsegup = 0;	// always do segment updates
	GcModeI = 0.0;
	m_atol = 0.01;
	m_gtol = 0;
	m_naugmin = 0;
	m_naugmaxmod =10;
	m_naugmax = 20;
	m_srad = 1.0;




	m_dxtol = 0;

	// Note that friction has not been implemented yet
	m_mu = 0;
	m_epsf = 0;

	m_ss.SetSibling(&m_ms);
	m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------
//! build the matrix profile for use in the stiffness matrix
void FEFacet2FacetSliding::BuildMatrixProfile(FEStiffnessMatrix& K)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	vector<int> lm(6*FEElement::MAX_NODES*2);

	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0? m_ms : m_ss);

		for (int j=0; j<ss.Elements(); ++j)
		{
			FESurfaceElement& se = ss.Element(j);
			int nint = se.GaussPoints();
			int* sn = &se.m_node[0];
			for (int k=0; k<nint; ++k)
			{
				FEFacetSlidingSurface::Data& pt = ss.m_Data[j][k];
				FESurfaceElement* pe = pt.m_pme;
				if (pe != 0)
				{
					FESurfaceElement& me = dynamic_cast<FESurfaceElement&> (*pe);
					int* mn = &me.m_node[0];

					assign(lm, -1);

					int nseln = se.Nodes();
					int nmeln = me.Nodes();

					for (int l=0; l<nseln; ++l)
					{
						vector<int>& id = mesh.Node(sn[l]).m_ID;
						lm[6*l  ] = id[DOF_X];
						lm[6*l+1] = id[DOF_Y];
						lm[6*l+2] = id[DOF_Z];
						lm[6*l+3] = id[DOF_RU];
						lm[6*l+4] = id[DOF_RV];
						lm[6*l+5] = id[DOF_RW];
					}

					for (int l=0; l<nmeln; ++l)
					{
						vector<int>& id = mesh.Node(mn[l]).m_ID;
						lm[6*(l+nseln)  ] = id[DOF_X];
						lm[6*(l+nseln)+1] = id[DOF_Y];
						lm[6*(l+nseln)+2] = id[DOF_Z];
						lm[6*(l+nseln)+3] = id[DOF_RU];
						lm[6*(l+nseln)+4] = id[DOF_RV];
						lm[6*(l+nseln)+5] = id[DOF_RW];
					}

					K.build_add(lm);
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! Initialization routine
bool FEFacet2FacetSliding::Init()
{
	m_bfirst = true;
	m_normg0 = 0.0;

	// initialize surface data
	if (m_ss.Init() == false) return false;
	if (m_ms.Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Activate()
{
	// don't forget the base class
	FEContactInterface::Activate();

	// calculate penalty factors
	if (m_bautopen) CalcAutoPenalty(m_ss);

	// project slave surface onto master surface
	ProjectSurface(m_ss, m_ms, true);

	if (m_btwo_pass) 
	{
		ProjectSurface(m_ms, m_ss, true);
		if (m_bautopen) CalcAutoPenalty(m_ms);
	}

	// check friction parameters
	// since friction has not been implemented yet
	if ((m_mu != 0) || (m_epsf != 0))
	{
		felog.printbox("WARNING", "Friction has NOT been implemented yet for facet-to-facet contact\ninterfaces. Friction parameters are ignored.");
		m_mu = 0;
		m_epsf = 0;
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::CalcAutoPenalty(FEFacetSlidingSurface& s)
{
	// get the mesh
	FEMesh& m = GetFEModel()->GetMesh();

	// loop over all surface elements
	for (int i=0; i<s.Elements(); ++i)
	{
		// get the surface element
		FESurfaceElement& el = s.Element(i);

		// find the element this face belongs to
		FEElement* pe = m.FindElementFromID(el.m_nelem);
		assert(pe);

		// get the area of the surface element
		double A = s.FaceArea(el);

		// get the volume of the volume element
		double V = m.ElementVolume(*pe);

		// calculate a modulus
		double K = AutoPenalty(el, s);

		// calculate penalty
		double eps = K*A/V;

		// assign to integation points of surface element
		int nint = el.GaussPoints();
		for (int j=0; j<nint; ++j) 
		{
			FEFacetSlidingSurface::Data& pt = s.m_Data[i][j];
			pt.m_eps = eps;
		}
	}
}

//-----------------------------------------------------------------------------
//! In this function we project the integration points to the master surface,
//! calculate the projection's natural coordinates and normal vector
//
void FEFacet2FacetSliding::ProjectSurface(FEFacetSlidingSurface &ss, FEFacetSlidingSurface &ms, bool bsegup)
{
	FEClosestPointProjection cpp(ms);
	cpp.SetTolerance(m_stol);
	cpp.Init();

	// loop over all slave elements
	for (int i=0; i<ss.Elements(); ++i)
	{
		// get the slave element
		FESurfaceElement& se = ss.Element(i);
		int nn = se.Nodes();

		// get nodal coordinates
		vec3d re[FEElement::MAX_NODES];
		for (int l=0; l<nn; ++l) re[l] = ss.GetMesh()->Node(se.m_node[l]).m_rt;

		// loop over all its integration points
		int nint = se.GaussPoints();
		for (int j=0; j<nint; ++j)
		{
			// get the integration point data
			FEFacetSlidingSurface::Data& pt = ss.m_Data[i][j];

			// calculate the global coordinates of this integration point
			double* H = se.H(j);

			vec3d x(0,0,0), q, qt;
			for (int k=0; k<nn; ++k) x += re[k]*H[k];

			// see if the point still projects to the same element
			if (pt.m_pme)
			{
				// update projection to master element
				FESurfaceElement& mel = *pt.m_pme;
				q = ms.ProjectToSurface(mel, x, pt.m_rs[0], pt.m_rs[1]);

				//// tangential projection
				// onto the master element
				double rx = ss.m_rs_n[se.m_nID][j][0];
				double sx = ss.m_rs_n[se.m_nID][j][1];

				// get the coordinates at the previous step
				double rp = ss.m_rsp_n[se.m_nID][j][0];
				double sp = ss.m_rsp_n[se.m_nID][j][1];

				qt = ms.ProjectToSurface(mel, x, rx, sx);
				ss.m_rs_n[se.m_nID][j][0]=rx;
				ss.m_rs_n[se.m_nID][j][1]=sx;
				///////
				// see if the projection is still in the element
				if (bsegup && (!ms.IsInsideElement(mel, pt.m_rs[0], pt.m_rs[1], m_stol)))
				{
					// if not, do a new search
					FESurfaceElement* pold = pt.m_pme;
					pt.m_rs = vec2d(0,0);
					ss.m_rs_n[se.m_nID][j] = vec2d(0, 0);
					FESurfaceElement* pme = cpp.Project(x, q, pt.m_rs);
					FESurfaceElement* pmetxy;
					pmetxy = cpp.Project(x, qt, ss.m_rs_n[se.m_nID][j]);
					pt.m_pme = pme;
					if (pmetxy == 0)
					{
						// nope, if has genuinly left contact
						int* n = &pold->m_node[0];
						//						log.printf("node %d has left element (%d, %d, %d, %d)\n", m+1, n[0]+1, n[1]+1, n[2]+1, n[3]+1);
					}
					else 
					{
						// the node has moved to another master segment.
						// If friction is active we need to translate the frictional
						// data to the new master segment.
						FESurfaceElement& eo = *pold;
						FESurfaceElement& en = *pme;
						MapTangentialComponent(j, ss, ms, en, eo, qt);
					}

				}
			}
			if (bsegup)
			{
				// find the master segment this element belongs to
				pt.m_rs = vec2d(0,0);
				FESurfaceElement* pme = cpp.Project(x, q, pt.m_rs);
				pt.m_pme = pme;
				FESurfaceElement* pmet ;
				pmet = cpp.Project(x, qt, ss.m_rs_n[se.m_nID][j]);
				if (pmet)
				{
					// the node has come into contact so make sure to initialize
					// the previous natural coordinates for friction.
					ss.m_rsp_n[se.m_nID][j] = ss.m_rs_n[se.m_nID][j];
				}
			}

			// update normal and gap at integration point
			if (pt.m_pme)
			{
				double r = pt.m_rs[0];
				double s = pt.m_rs[1];

				//FESurfaceElement& mel = *ss.;
				//ss.m_M[i] = ss.Metric0(mel, r, s);
				double rx = ss.m_rs_n[se.m_nID][j][0];
				double sx = ss.m_rs_n[se.m_nID][j][1];
				ss.m_M_t[se.m_nID][j] = ss.Metric0(*pt.m_pme, rx, sx);

				// the slave normal is set to the master element normal
				pt.m_nu = ms.SurfaceNormal(*pt.m_pme, r, s);

				// calculate gap
				pt.m_gap = -pt.m_nu*(x - q);
			}
			else
			{
				// since the node is not in contact, we set the gap and Lagrange multiplier to zero
				//pt.m_gap = 0;
				pt.m_Lm = 0;
				ss.m_Lnt[se.m_nID][j][0] = 0;
				ss.m_Lnt[se.m_nID][j][1] = 0;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Update(int niter)
{
	FEModel& fem = *GetFEModel();

	// should we do a segment update or not?
	// TODO: check what happens when m_nsegup == -1 and m_npass = 2;
	// We have to make sure that in this case, both surfaces get at least
	// one pass!
	bool bupdate = (m_bfirst || (m_nsegup == 0)? true : (niter <= m_nsegup));

	// project slave surface to master surface
	ProjectSurface(m_ss, m_ms, bupdate);
	if (m_btwo_pass) ProjectSurface(m_ms, m_ss, bupdate);


	//bool test = true;
	//test begin
	//for (int i = 0; i < m_ss.Nodes(); ++i)
	//{
		//if (m_ss.m_Data[i]. < 0){
		//	test = false;
		//}

	//}

	//if (test){// Update the net contact pressures

		// Update the net contact pressures
		UpdateContactPressures();
	//}
	m_bfirst = false;
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ShallowCopy(DumpStream& dmp, bool bsave)
{
	m_ss.ShallowCopy(dmp, bsave);
	m_ms.ShallowCopy(dmp, bsave);
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::ContactForces(FEGlobalVector& R)
{
	vector<int> sLM, mLM, LM, en;
	vector<double> fe;

	////////////////////////////////////////////////
	//vec3d dxr, dxs;

	// normal force
	//double tn, Ln;

	// gap function
	///double gap;

	// tangents
	vec3d tau1, tau2;

	// max nr of master element nodes
	

	// master element nodes
	

	// shape function values
	



	// surface metrics
	double A[2][2], M[2][2], K[2][2];
	double detA;
	///////////////////////////////////////////////

	//FEClosestPointProjection ttt;
	//0.0148933197058=gap(infinitesimal)[δgap]/search tolerance
	double tolerance = 0.0148933197058*m_stol;
	//double tolerance = m_stol*5.5988986099716 * (10 ^ (-10));
	double infgap = (5.551115 / 151)*(exp(-18));
	felog.printf("    tolerance value : %1  \n", tolerance);
	felog.printf("    tolerance value : %15le\n", tolerance);
	const int MELN = FEElement::MAX_NODES;  //const int MAXMN = FEElement::MAX_NODES;
	double detJ[MELN], w[MELN], *Hs, H[MELN], Hm[MELN];  //double H[MAXMN], Hr[MAXMN], Hs[MAXMN];
	// contact vectors
	double N[3 * (MELN + 1)], N1[3 * (MELN + 1)], N2[3 * (MELN + 1)];
	double T1[3 * (MELN + 1)], T2[3 * (MELN + 1)], D1[3 * (MELN + 1)], D2[3 * (MELN + 1)];

	vec3d r0[MELN];  //vec3d rtm[MAXMN];

	int npass = (m_btwo_pass?2:1);
	for (int np=0; np<npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		for (int i=0; i<ss.Elements(); ++i)
		{
			FESurfaceElement& se = ss.Element(i);
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// get the element's LM vector
			ss.UnpackLM(se, sLM);

			// nodal coordinates
			for (int j=0; j<nseln; ++j) r0[j] = ss.GetMesh()->Node(se.m_node[j]).m_r0;

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (int j=0; j<nint; ++j)
			{
				double* Gr = se.Gr(j);
				double* Gs = se.Gs(j);

				// calculate jacobian
				// note that we are integrating over the reference surface
				vec3d dxr, dxs;
				for (int k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				// jacobians
				detJ[j] = (dxr ^ dxs).norm();

				// integration weights
				w[j] = se.GaussWeights()[j];
			}

			// loop over all integration points
			for (int j=0; j<nint; ++j)
			{
				// get integration point data
				FEFacetSlidingSurface::Data& pt = ss.m_Data[i][j];

				// get the master element
				FESurfaceElement* pme = pt.m_pme;
				//FESurfaceElement* psle = pt.;
				if (pme)
				{
					FESurfaceElement& me = *pme;

					int nmeln = me.Nodes();
					ms.UnpackLM(me, mLM);

					// calculate degrees of freedom
					int ndof = 3*(nseln + nmeln);

					// build the LM vector
					LM.resize(ndof);
					for (int k=0; k<nseln; ++k)
					{
						LM[3*k  ] = sLM[3*k  ];
						LM[3*k+1] = sLM[3*k+1];
						LM[3*k+2] = sLM[3*k+2];
					}

					for (int k=0; k<nmeln; ++k)
					{
						LM[3*(k+nseln)  ] = mLM[3*k  ];
						LM[3*(k+nseln)+1] = mLM[3*k+1];
						LM[3*(k+nseln)+2] = mLM[3*k+2];
					}

					// build the en vector
					en.resize(nseln+nmeln);
					for (int k=0; k<nseln; ++k) en[k] = se.m_node[k];
					for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

					// calculate shape functions
					Hs = se.H(j);

					double r = pt.m_rs[0];
					double s = pt.m_rs[1];
					//double rp = pt.[0];
					//double sp = pt.m_rs[1];
					//me.shape_fnc(Hm, r, s);

					/////////////////////////////////

					mat2d &Mk = ss.m_M_t[se.m_nID][j];
					mat2d Mki = Mk.inverse();

					// get the master element node positions
					//for (k = 0; k<nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;

					// isoparametric coordinates of the projected slave node
					// onto the master element
					double rx = ss.m_rs_n[se.m_nID][j][0];
					double sx = ss.m_rs_n[se.m_nID][j][1];

					// get the coordinates at the previous step
					double rp = ss.m_rsp_n[se.m_nID][j][0];
					double sp = ss.m_rsp_n[se.m_nID][j][1];

					// get the master shape function values at this slave node
					me.shape_fnc(H, rx, sx);
					/////////////////////////////////






					// get normal vector
					vec3d nu = pt.m_nu;

					// gap function
					double g = pt.m_gap;
					
					// lagrange multiplier
					double Lm = pt.m_Lm;

					// penalty value
					double eps = m_epsn*pt.m_eps;

					felog.printf("    penalty value : %15le\n", eps);

					felog.printf("    gap function : %15le\n", g );

					//felog.printf("    normal force : %15le", lnorm);

					if (eps == 0){
						felog.printf("Division by zero!....................penalty value==0");
						break;
						exit(1);
					}
					int maxcohetraction = m_naugmax;
					felog.printf("maxcohetraction:  %u \n", maxcohetraction);
					//if (eps = !0){
					felog.printf("maxcohetraction:  %u \n", maxcohetraction);
					felog.printf("rootcontrol:  %u \n", rootcontrol);
					
						double temp;
						
						double cohesivetolerance = (tolerance*((maxcohetraction / eps)*exp(1 - (tolerance / eps))) + ((tolerance / eps) * maxcohetraction*(-1 / eps)*exp(1 - (tolerance / eps))));
						double crtl = ((infgap / eps) * maxcohetraction * exp(1 - (infgap / eps)));
						if (g < 0){ g = -g; }
						// contact traction( cohesive case! mode I)
						double tn = -((g / eps) * maxcohetraction * exp(1 - (g / eps)));
						// contact traction( cohesive case! mode II and III)
						double tt = -(2*sqrt(0.5*exp(1))*maxcohetraction *(1 + (g / eps))*exp(-(g / eps)));

						if ((g> rootcontrol*eps) && ((-tn) < crtl)){
							tn = 0;
						}
						tn = -(MBRACKET(-tn));
						tt = -(MBRACKET(-tt));
						temp = tn;
					//}
						felog.printf("    contact traction value before : %15le\n", tn);
						// the projection tolerance( set with value of search tolerance)
						
						/*
						double	m_tol;	//!< projection tolerance
						*/
						felog.printf("    cohesivetolerance : %15le\n", cohesivetolerance);

						
						//Mix Mode begin

						// --- T A N G E N T I A L   T R A C T I O N ---


							// get the master shape function derivative values at this slave node
						me.shape_deriv(Hm, Hs, rx, sx);

							// get the tangent vectors
							tau1 = tau2 = vec3d(0, 0, 0);
							for (int k = 0; k<nmeln; ++k)
							{
								tau1.x += Hm[k] * r0[k].x;
								tau1.y += Hm[k] * r0[k].y;
								tau1.z += Hm[k] * r0[k].z;

								tau2.x += Hs[k] * r0[k].x;
								tau2.y += Hs[k] * r0[k].y;
								tau2.z += Hs[k] * r0[k].z;
							}

							// set up the Ti vectors
							T1[0] = tau1.x; T2[0] = tau2.x;
							T1[1] = tau1.y; T2[1] = tau2.y;
							T1[2] = tau1.z; T2[2] = tau2.z;

							for (int k = 0; k<nmeln; ++k)
							{
								T1[(k + 1) * 3] = -H[k] * tau1.x;
								T1[(k + 1) * 3 + 1] = -H[k] * tau1.y;
								T1[(k + 1) * 3 + 2] = -H[k] * tau1.z;

								T2[(k + 1) * 3] = -H[k] * tau2.x;
								T2[(k + 1) * 3 + 1] = -H[k] * tau2.y;
								T2[(k + 1) * 3 + 2] = -H[k] * tau2.z;
							}

							// set up the Ni vectors
							N1[0] = N2[0] = 0;
							N1[1] = N2[1] = 0;
							N1[2] = N2[2] = 0;

							for (int k = 0; k<nmeln; ++k)
							{
								N1[(k + 1) * 3] = -Hm[k] * nu.x;
								N1[(k + 1) * 3 + 1] = -Hm[k] * nu.y;
								N1[(k + 1) * 3 + 2] = -Hm[k] * nu.z;

								N2[(k + 1) * 3] = -Hs[k] * nu.x;
								N2[(k + 1) * 3 + 1] = -Hs[k] * nu.y;
								N2[(k + 1) * 3 + 2] = -Hs[k] * nu.z;
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
							me.shape_deriv2(Grr, Grs, Gss, rx, sx);
							for (int k = 0; k<nmeln; ++k)
							{
								K[0][0] += (nu*r0[k])*Grr[k];
								K[0][1] += (nu*r0[k])*Grs[k];
								K[1][0] += (nu*r0[k])*Grs[k];
								K[1][1] += (nu*r0[k])*Gss[k];
							}

							// setup A matrix
							A[0][0] = M[0][0] + g *K[0][0];
							A[0][1] = M[0][1] + g *K[0][1];
							A[1][0] = M[1][0] + g *K[1][0];
							A[1][1] = M[1][1] + g *K[1][1];

							detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];
							double *tntt, *tttn;
							// setup Di vectors
							for (int k = 0; k<ndof; ++k)
							{
								D1[k] = (1 / detA)*(A[1][1] * (T1[k] + g *N1[k]) - A[0][1] * (T2[k] + g *N2[k]));
								D2[k] = (1 / detA)*(A[0][0] * (T2[k] + g *N2[k]) - A[0][1] * (T1[k] + g *N1[k]));
								// T A N G E N T I A L   T R A C T I O N  trail
								tntt[k] = D1[k];  // exp(-(((D1[k])*(D1[k])) + ((D2[k])*(D2[k]))) / (eps*eps));
								tttn[k] = D2[k];  //(exp(-(((D1[k])*(D1[k])) + ((D2[k])*(D2[k]))) / (eps*eps)))*((sqrt((((D1[k])*(D1[k])) + ((D2[k])*(D2[k])))) / eps));
							}
							double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
							double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
							double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
							double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
							//cut here !
								/*
							// calculate friction tractions
							// a. calculate trial state
							Tt[0] = Lt[0] + m_epsf*(Mk[0][0] * (r - rp) + Mk[0][1] * (s - sp));
							Tt[1] = Lt[1] + m_epsf*(Mk[1][0] * (r - rp) + Mk[1][1] * (s - sp));

							double TMT = Tt[0] * (Mki[0][0] * Tt[0] + Mki[0][1] * Tt[1]) + Tt[1] * (Mki[1][0] * Tt[0] + Mki[1][1] * Tt[1]);
							assert(TMT >= 0);

							double phi = sqrt(TMT) - m_mu*Ln;

							// b. return map
							if (phi > 0)
							{
								Tt[0] = m_mu*Ln*Tt[0] / sqrt(TMT);
								Tt[1] = m_mu*Ln*Tt[1] / sqrt(TMT);
							}

							// tangential force vector
							for (l = 0; l<ndof; ++l) fe[l] -= test*(Tt[0] * D1[l] + Tt[1] * D2[l]);
						}





						//Mix Mode end

						*/


					
					// contact traction smart contact part
						/*
						if ((-g)> rootcontrol*eps && ((-tn)<cohesivetolerance)){
						felog.printf("smart part  back to original....................part ContactForces");
						//FEFacetSlidingSurface::Data::Data();
						
					
						
						double g = pt.m_gap;
						 tn = Lm + eps*g;
						tn = MBRACKET(tn);
					}
					*/
						felog.printf("   contact traction value after : %15le\n", tn);
						felog.printf("   contact T A N G E N T I A L  traction value after : %15le\n", tt);
					// calculate the force vector
					fe.resize(ndof);

					for (int k=0; k<nseln; ++k)
					{
						fe[3*k  ] = Hs[k]*nu.x;
						fe[3*k+1] = Hs[k]*nu.y;
						fe[3*k+2] = Hs[k]*nu.z;
					}

					for (int k=0; k<nmeln; ++k)
					{
						fe[3*(k+nseln)  ] = -Hm[k]*nu.x;
						fe[3*(k+nseln)+1] = -Hm[k]*nu.y;
						fe[3*(k+nseln)+2] = -Hm[k]*nu.z;
					}
					double tempGC = 0;
					//for (int k=0; k<ndof; ++k) fe[k] *= tn*detJ[j]*w[j];
					// force mixed mode
					// tn*N
					for (int k = 0; k<ndof; ++k) fe[k] *= ((tn* Utxn) + (tn*Utyn))*detJ[j] * w[j];
					// tt*(D1+D2)
					for (int k = 0; k<ndof; ++k) fe[k] += ((tt*tntt[k] * Utxt) + (tt*tttn[k] * Utyt))*detJ[j] * w[j];
					for (int m = 0; m<ndof; ++m) tempGC+=fe[m];
					
					felog.printf("Gc(I) local value  %15le \n", tempGC);
					GcModeI = GcModeI + tempGC;
					felog.printf("GcMode(I) total value   %15le \n", GcModeI);
					 temp = 0;
					//test begin
					if (tn == 0){
						//ss.GetMesh()->Node(se.m_node[j]).m_r0
						//ss.GetMesh()->Node(se.m_node[j]).m_pt = 0;
						//ss.GetMesh()->Node(me.m_node[j]).m_at.x = -ss.GetMesh()->Node(me.m_node[j]).m_ap.x;
						//ss.GetMesh()->Node(me.m_node[j]).m_at.y = -ss.GetMesh()->Node(me.m_node[j]).m_ap.y;
						//ss.GetMesh()->Node(me.m_node[j]).m_at.z = -ss.GetMesh()->Node(me.m_node[j]).m_ap.z;

							}
					else{

						// assemble the global residual
						R.Assemble(en, LM, fe);
					}
					//test end
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// cohesive Part of force

vec3d FEFacet2FacetSliding::GetCohesiveContactForce(){
	vec3d f(0, 0, 0);
	return f;

}
//




void FEFacet2FacetSliding::ContactStiffness(FESolver* psolver)
{
	vector<int> sLM, mLM, LM, en;
	//0.0148933197058=gap(infinitesimal)[δgap]/search tolerance
	//double tolerance = 0.0148933197058*m_stol;
	double tolerance = m_stol*5.5988986099716 * (10 ^ (-10));
	double infgap = (5.551115 / 151)*(exp(-18));
	const int MN = FEElement::MAX_NODES;
	const int ME = 3*MN*2;
	double N[ME], T1[ME], T2[ME], N1[ME] = {0}, N2[ME] = {0}, D1[ME], D2[ME], Nb1[ME], Nb2[ME];
	matrix ke;
	vec3d tau1, tau2;
	// surface metrics
	double A[2][2], M[2][2], K[2][2];
	double detA;
	// get the mesh
	FEMesh* pm = m_ss.GetMesh();

	// see how many reformations we've had to do so far
	int nref = psolver->m_nref;

	// get the "size" of the model
	// We need this to scale the insertion distance
	double R = GetFEModel()->GetMesh().GetBoundingBox().radius();
	double dxtol = R*m_dxtol;

	// set higher order stiffness mutliplier
	double knmult = m_knmult;
	if (m_knmult < 0)
	{
		int ni = int(-m_knmult);
		if (nref >= ni)
		{
			knmult = 1; 
			felog.printf("Higher order stiffness terms included.\n");
		}
		else knmult = 0;
	}

	double detJ[MN], w[MN], *Hs, Hm[MN], Hmr[MN], Hms[MN],H[MN] ;
	vec3d r0[MN];

	int npass = (m_btwo_pass?2:1);
	for (int np=0; np < npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0? m_ms : m_ss);

		// loop over all slave elements
		for (int i=0; i<ss.Elements(); ++i)
		{
			FESurfaceElement& se = ss.Element(i);
			int nseln = se.Nodes();
			int nint = se.GaussPoints();

			// get the element's LM vector
			ss.UnpackLM(se, sLM);

			// nodal coordinates
			for (int j=0; j<nseln; ++j) r0[j] = ss.GetMesh()->Node(se.m_node[j]).m_r0;

			// we calculate all the metrics we need before we
			// calculate the nodal forces
			for (int j=0; j<nint; ++j)
			{
				double* Gr = se.Gr(j);
				double* Gs = se.Gs(j);

				// calculate jacobian
				// note that we are integrating over the reference surface
				vec3d dxr, dxs;
				for (int k=0; k<nseln; ++k)
				{
					dxr.x += Gr[k]*r0[k].x;
					dxr.y += Gr[k]*r0[k].y;
					dxr.z += Gr[k]*r0[k].z;

					dxs.x += Gs[k]*r0[k].x;
					dxs.y += Gs[k]*r0[k].y;
					dxs.z += Gs[k]*r0[k].z;
				}

				// jacobians
				detJ[j] = (dxr ^ dxs).norm();

				// integration weights
				w[j] = se.GaussWeights()[j];
			}

			// loop over all integration points
			for (int j=0; j<nint; ++j)
			{
				// get integration point data
				FEFacetSlidingSurface::Data& pt = ss.m_Data[i][j];

				// get the master element
				FESurfaceElement* pme = pt.m_pme;
				if (pme)
				{
					FESurfaceElement& me = *pme;

					int nmeln = me.Nodes();
					ms.UnpackLM(me, mLM);

					// calculate degrees of freedom
					int ndof = 3*(nseln + nmeln);

					// build the LM vector
					LM.resize(ndof);
					for (int k=0; k<nseln; ++k)
					{
						LM[3*k  ] = sLM[3*k  ];
						LM[3*k+1] = sLM[3*k+1];
						LM[3*k+2] = sLM[3*k+2];
					}

					for (int k=0; k<nmeln; ++k)
					{
						LM[3*(k+nseln)  ] = mLM[3*k  ];
						LM[3*(k+nseln)+1] = mLM[3*k+1];
						LM[3*(k+nseln)+2] = mLM[3*k+2];
					}

					// build the en vector
					en.resize(nseln+nmeln);
					for (int k=0; k<nseln; ++k) en[k] = se.m_node[k];
					for (int k=0; k<nmeln; ++k) en[k+nseln] = me.m_node[k];

					// calculate shape functions
					Hs = se.H(j);
					double r = pt.m_rs[0];
					double s = pt.m_rs[1];
					//me.shape_fnc(Hm, r, s);

					/////////////////////////////////

					mat2d &Mk = ss.m_M_t[se.m_nID][j];
					mat2d Mki = Mk.inverse();

					// get the master element node positions
					//for (k = 0; k<nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;

					// isoparametric coordinates of the projected slave node
					// onto the master element
					double rx = ss.m_rs_n[se.m_nID][j][0];
					double sx = ss.m_rs_n[se.m_nID][j][1];

					// get the coordinates at the previous step
					double rp = ss.m_rsp_n[se.m_nID][j][0];
					double sp = ss.m_rsp_n[se.m_nID][j][1];

					// get the master shape function values at this slave node
					me.shape_fnc(H, rx, sx);
					/////////////////////////////////




					// get normal vector
					vec3d nu = pt.m_nu;

					// gap function
					double g = pt.m_gap;
					
					// lagrange multiplier
					double Lm = pt.m_Lm;

					// penalty value
					double eps = m_epsn*pt.m_eps;

					double temp;
					// contact traction( cohesive case!)
					//felog.printf("smart part just test");
					if (eps == 0){
						felog.printf("Division by zero!....................penalty value==0");
						break;
					}
					int maxcohetraction = m_naugmax;
					double cohesivetolerance = (tolerance*((maxcohetraction / eps)*exp(1 - (tolerance / eps))) + ((tolerance / eps) * maxcohetraction*(-1 / eps)*exp(1 - (tolerance / eps))));
					if (g < 0){ g = -g; }
					double crtl = ((infgap / eps) * maxcohetraction * exp(1 - (infgap / eps)));
					double tn = -((g / eps) * maxcohetraction * exp(1 - (g / eps)));
					// contact traction( cohesive case! mode II and III)
					double tt = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (g / eps))*exp(-(g / eps)));
					if ((g> rootcontrol*eps) && ((-tn) < crtl)){
						tn = 0;
					}
					tn =-(MBRACKET(-tn));
					tt = -(MBRACKET(-tt));
					temp = tn;
					// the directional derivative of traction-displacement function in the direction of search tolerance with respect to initial tolerance
					
					double dtn = -HEAVYSIDE(((g / eps) * maxcohetraction * exp(1 - (g / eps))))*(((maxcohetraction / eps)*exp(1 - (g / eps)) + (g / eps) * maxcohetraction*(-1 / eps)*exp(1 - (g / eps))));
					double dtt = -HEAVYSIDE((2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (g / eps))*exp(-(g / eps))))*(2 * sqrt(0.5*exp(1))*maxcohetraction)*(((1 / eps)*exp(-(g / eps))) + ((1 + (g / eps))*(-1 / eps)*exp(-(g / eps))));




					//Mix Mode begin

					// --- T A N G E N T I A L   T R A C T I O N ---


					// get the master shape function derivative values at this slave node
					me.shape_deriv(Hm, Hs, rx, sx);

					// get the tangent vectors
					tau1 = tau2 = vec3d(0, 0, 0);
					for (int k = 0; k<nmeln; ++k)
					{
						tau1.x += Hm[k] * r0[k].x;
						tau1.y += Hm[k] * r0[k].y;
						tau1.z += Hm[k] * r0[k].z;

						tau2.x += Hs[k] * r0[k].x;
						tau2.y += Hs[k] * r0[k].y;
						tau2.z += Hs[k] * r0[k].z;
					}

					// set up the Ti vectors
					T1[0] = tau1.x; T2[0] = tau2.x;
					T1[1] = tau1.y; T2[1] = tau2.y;
					T1[2] = tau1.z; T2[2] = tau2.z;

					for (int k = 0; k<nmeln; ++k)
					{
						T1[(k + 1) * 3] = -H[k] * tau1.x;
						T1[(k + 1) * 3 + 1] = -H[k] * tau1.y;
						T1[(k + 1) * 3 + 2] = -H[k] * tau1.z;

						T2[(k + 1) * 3] = -H[k] * tau2.x;
						T2[(k + 1) * 3 + 1] = -H[k] * tau2.y;
						T2[(k + 1) * 3 + 2] = -H[k] * tau2.z;
					}

					// set up the Ni vectors
					N1[0] = N2[0] = 0;
					N1[1] = N2[1] = 0;
					N1[2] = N2[2] = 0;

					for (int k = 0; k<nmeln; ++k)
					{
						N1[(k + 1) * 3] = -Hm[k] * nu.x;
						N1[(k + 1) * 3 + 1] = -Hm[k] * nu.y;
						N1[(k + 1) * 3 + 2] = -Hm[k] * nu.z;

						N2[(k + 1) * 3] = -Hs[k] * nu.x;
						N2[(k + 1) * 3 + 1] = -Hs[k] * nu.y;
						N2[(k + 1) * 3 + 2] = -Hs[k] * nu.z;
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
					me.shape_deriv2(Grr, Grs, Gss, rx, sx);
					for (int k = 0; k<nmeln; ++k)
					{
						K[0][0] += (nu*r0[k])*Grr[k];
						K[0][1] += (nu*r0[k])*Grs[k];
						K[1][0] += (nu*r0[k])*Grs[k];
						K[1][1] += (nu*r0[k])*Gss[k];
					}

					// setup A matrix
					A[0][0] = M[0][0] + g *K[0][0];
					A[0][1] = M[0][1] + g *K[0][1];
					A[1][0] = M[1][0] + g *K[1][0];
					A[1][1] = M[1][1] + g *K[1][1];

					detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];
					double *tntt, *tttn;
					// setup Di vectors
					for (int k = 0; k<ndof; ++k)
					{
						D1[k] = (1 / detA)*(A[1][1] * (T1[k] + g *N1[k]) - A[0][1] * (T2[k] + g *N2[k]));
						D2[k] = (1 / detA)*(A[0][0] * (T2[k] + g *N2[k]) - A[0][1] * (T1[k] + g *N1[k]));
						// T A N G E N T I A L   T R A C T I O N  trail
						tntt[k] = D1[k];  // exp(-(((D1[k])*(D1[k])) + ((D2[k])*(D2[k]))) / (eps*eps));
						tttn[k] = D2[k];  //(exp(-(((D1[k])*(D1[k])) + ((D2[k])*(D2[k]))) / (eps*eps)))*((sqrt((((D1[k])*(D1[k])) + ((D2[k])*(D2[k])))) / eps));
					}
					double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
					double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
					double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
					double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
					//cut here !
					
					///the directional derivative of the trails
					double DUtxn = ((-2 * (((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))))) / (eps*eps))*exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
					double DUtyn = ((-2 * (((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))))) / (eps*eps))*exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
					double DUtxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*(1 / eps) + ((-2 * (((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))))) / (eps*eps))*(exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
					double DUtyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*(1 / eps) + ((-2 * (((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))))) / (eps*eps))*(exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));

					
					//double tn
					
					//((tn*tntt[k] * Utxn) + (tn*tttn[k] * Utyn) + (tt*tntt[k] * Utxt) + (tt*tttn[k] * Utyt))

					/* smart contact part later
					// contact traction
					*/

					/*
					if ((-g)> rootcontrol*eps && ((-tn) < cohesivetolerance)){
						felog.printf("smart part  back to original....................part ContactStiffness");
						tn = Lm + eps*g;
						tn = MBRACKET(tn);

						 dtn = eps*HEAVYSIDE(Lm + eps*g);



						// smart contact part later
						// define buffer layer for penalty insertion
						if ((dtn < 1e-7) && (g < 0) && (dxtol != 0))
						{
							if (dxtol < 0) dtn = eps*exp(-g / dxtol);
							else if (-g <= dxtol) dtn = eps*(1 + g / dxtol);
						}
					} 
					*/

					// calculate the N-vector
					for (int k=0; k<nseln; ++k)
					{
						N[3*k  ] = Hs[k]*nu.x;
						N[3*k+1] = Hs[k]*nu.y;
						N[3*k+2] = Hs[k]*nu.z;
					}

					for (int k=0; k<nmeln; ++k)
					{
						N[3*(k+nseln)  ] = -Hm[k]*nu.x;
						N[3*(k+nseln)+1] = -Hm[k]*nu.y;
						N[3*(k+nseln)+2] = -Hm[k]*nu.z;
					}

					// --- N O R M A L + T A N G E N T I A L    S T I F F N E S S ---

					// create the stiffness matrix mix mode
					ke.resize(ndof, ndof);

					// add the first order term (= D(tn)*dg )
					// Δtn.δg+Δtt.δξ
					for (int k=0; k<ndof; ++k)
						for (int l = 0; l<ndof; ++l) ke[k][l] = (dtn*((Utxn) + (Utyn))*N[k] * N[l] * detJ[j] * w[j]) + (dtt*((tntt[l] * Utxt) + (tttn[l] * Utyt))* N[k] * detJ[j] * w[j]);

					double tntotal = tn*((Utxn)+(Utyn));
					double tttotal = tt*((Utxt)+(Utyt));
					// tn.ΔUtxn.δg+tt.ΔUtxn.δξ
					double ktrail;
					for (i = 0; i < ndof; ++i)
						for (j = 0; j < ndof; ++j)
						{
						ktrail = tn*DUtxn*(Mk[0][0] * N[i] * tntt[j] + Mk[0][1] * N[i] * tttn[j]) + tn*DUtyn*(Mk[1][0] * N[i] * tntt[j] + Mk[1][1] * N[i] * tttn[j]) + tt*DUtxt*(Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j]) + tt*DUtyt*(Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j]);
						//kij = Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j] + Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j];
						//ke[i][j] += m_epsf*kij;
						ke[i][j] += ktrail*detJ[j] * w[j] * m_ktmult;
						}
					// add the higher order terms (= tn*D(dg) )
					if (knmult > 0)
					{
						// calculate the master shape fncs derivatives
						me.shape_deriv(Hmr, Hms, r, s);

						// get the master nodes
						vec3d rt[MN];
						for (int k = 0; k < nmeln; ++k) rt[k] = ms.GetMesh()->Node(me.m_node[k]).m_rt;

						// get the tangent vectors
						vec3d tau1(0, 0, 0), tau2(0, 0, 0);
						for (int k = 0; k < nmeln; ++k)
						{
							tau1.x += Hmr[k] * rt[k].x;
							tau1.y += Hmr[k] * rt[k].y;
							tau1.z += Hmr[k] * rt[k].z;

							tau2.x += Hms[k] * rt[k].x;
							tau2.y += Hms[k] * rt[k].y;
							tau2.z += Hms[k] * rt[k].z;
						}

						// set up the Ti vectors
						for (int k = 0; k < nseln; ++k)
						{
							T1[k * 3] = Hs[k] * tau1.x; T2[k * 3] = Hs[k] * tau2.x;
							T1[k * 3 + 1] = Hs[k] * tau1.y; T2[k * 3 + 1] = Hs[k] * tau2.y;
							T1[k * 3 + 2] = Hs[k] * tau1.z; T2[k * 3 + 2] = Hs[k] * tau2.z;
						}

						for (int k = 0; k < nmeln; ++k)
						{
							T1[(k + nseln) * 3] = -Hm[k] * tau1.x;
							T1[(k + nseln) * 3 + 1] = -Hm[k] * tau1.y;
							T1[(k + nseln) * 3 + 2] = -Hm[k] * tau1.z;

							T2[(k + nseln) * 3] = -Hm[k] * tau2.x;
							T2[(k + nseln) * 3 + 1] = -Hm[k] * tau2.y;
							T2[(k + nseln) * 3 + 2] = -Hm[k] * tau2.z;
						}

						// set up the Ni vectors
						for (int k = 0; k < nmeln; ++k)
						{
							N1[(k + nseln) * 3] = -Hmr[k] * nu.x;
							N1[(k + nseln) * 3 + 1] = -Hmr[k] * nu.y;
							N1[(k + nseln) * 3 + 2] = -Hmr[k] * nu.z;

							N2[(k + nseln) * 3] = -Hms[k] * nu.x;
							N2[(k + nseln) * 3 + 1] = -Hms[k] * nu.y;
							N2[(k + nseln) * 3 + 2] = -Hms[k] * nu.z;
						}

						// calculate metric tensor
						mat2d M;
						M[0][0] = tau1*tau1; M[0][1] = tau1*tau2;
						M[1][0] = tau2*tau1; M[1][1] = tau2*tau2;

						// calculate reciprocal metric tensor
						mat2d Mi = M.inverse();

						// calculate curvature tensor
						double K[2][2] = { 0 };
						double Grr[MN];
						double Gss[MN];
						double Grs[MN];
						me.shape_deriv2(Grr, Grs, Gss, r, s);
						for (int k = 0; k < nmeln; ++k)
						{
							K[0][0] += (nu*rt[k])*Grr[k];
							K[0][1] += (nu*rt[k])*Grs[k];
							K[1][0] += (nu*rt[k])*Grs[k];
							K[1][1] += (nu*rt[k])*Gss[k];
						}

						// setup A matrix A = M + gK
						double A[2][2];
						A[0][0] = M[0][0] + g*K[0][0];
						A[0][1] = M[0][1] + g*K[0][1];
						A[1][0] = M[1][0] + g*K[1][0];
						A[1][1] = M[1][1] + g*K[1][1];

						// calculate determinant of A
						double detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

						// setup Di vectors
						for (int k = 0; k < ndof; ++k)
						{
							D1[k] = (1 / detA)*(A[1][1] * (T1[k] + g*N1[k]) - A[0][1] * (T2[k] + g*N2[k]));
							D2[k] = (1 / detA)*(A[0][0] * (T2[k] + g*N2[k]) - A[0][1] * (T1[k] + g*N1[k]));
						}

						// setup Nbi vectors
						for (int k = 0; k < ndof; ++k)
						{
							Nb1[k] = N1[k] - K[0][1] * D2[k];
							Nb2[k] = N2[k] - K[0][1] * D1[k];
						}

						// add it to the stiffness
						double sum;
						// calculate the normalized traction vector
						vec3d pt = tau1 * (Utxn)+tau2 * (Utyn)+tau1 * (Utxt)+tau2 * (Utyt);
						pt.unit();

						/////////////////////////////

						// we need to define additional arrays for the tangential
						// contribution of the contact stiffness
						double T11[3 * (MN + 1)] = { 0 }, T12[3 * (MN + 1)] = { 0 }, T21[3 * (MN + 1)] = { 0 }, T22[3 * (MN + 1)] = { 0 };	// Tab matrix
						double N11[3 * (MN + 1)] = { 0 }, N12[3 * (MN + 1)] = { 0 }, N21[3 * (MN + 1)] = { 0 }, N22[3 * (MN + 1)] = { 0 };	// Nab matrix
						double P1[3 * (MN + 1)] = { 0 }, P2[3 * (MN + 1)] = { 0 };	// P arrays
						double Tb11[3 * (MN + 1)], Tb21[3 * (MN + 1)], Tb12[3 * (MN + 1)], Tb22[3 * (MN + 1)]; // Tbar matrix
						double Pb1[3 * (MN + 1)], Pb2[3 * (MN + 1)]; // Pbar array

						for (int k = 0; k < nmeln; ++k)
						{
							T11[3 * (k + 1)] = -Hmr[k] * tau1.x;
							T11[3 * (k + 1) + 1] = -Hmr[k] * tau1.y;
							T11[3 * (k + 1) + 2] = -Hmr[k] * tau1.z;

							T12[3 * (k + 1)] = -Hms[k] * tau1.x;
							T12[3 * (k + 1) + 1] = -Hms[k] * tau1.y;
							T12[3 * (k + 1) + 2] = -Hms[k] * tau1.z;

							T21[3 * (k + 1)] = -Hmr[k] * tau2.x;
							T21[3 * (k + 1) + 1] = -Hmr[k] * tau2.y;
							T21[3 * (k + 1) + 2] = -Hmr[k] * tau2.z;

							T22[3 * (k + 1)] = -Hms[k] * tau2.x;
							T22[3 * (k + 1) + 1] = -Hms[k] * tau2.y;
							T22[3 * (k + 1) + 2] = -Hms[k] * tau2.z;

							if (nmeln == 4)
							{
								N12[3 * (k + 1)] = N21[3 * (k + 1)] = -0.25*nu.x;
								N12[3 * (k + 1) + 1] = N21[3 * (k + 1) + 1] = -0.25*nu.y;
								N12[3 * (k + 1) + 2] = N21[3 * (k + 1) + 2] = -0.25*nu.z;
							}
							else if (nmeln == 6) assert(false);

							P1[3 * (k + 1)] = -Hmr[k] * pt.x;
							P1[3 * (k + 1) + 1] = -Hmr[k] * pt.y;
							P1[3 * (k + 1) + 2] = -Hmr[k] * pt.z;

							P2[3 * (k + 1)] = -Hms[k] * pt.x;
							P2[3 * (k + 1) + 1] = -Hms[k] * pt.y;
							P2[3 * (k + 1) + 2] = -Hms[k] * pt.z;
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
						else
						{
							assert(false);
						}
						double gt1 = g12*tau1;
						double gt2 = g12*tau2;
						double gp = g12*pt;

						for (int k = 0; k < ndof; ++k)
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
						for (int k = 0; k < 2; ++k)
							for (int l = 0; l < 2; ++l)
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

						// calculate stiffness matrix
						// NOTE: I think I need Mi (iso Mki) for KT1 and KT2. I only need Mk (iso M) for the direct stiffnesses

						// tn.Utxn.Δδg+tt.Utxn.Δδξ

						double kij;

						for (i = 0; i < ndof; ++i)
							for (j = 0; j < ndof; ++j)
							{
							// KT1
							kij = T11[i] * D1[j] + T12[i] * D2[j];
							kij += D1[i] * T11[j] + D2[i] * T12[j];
							kij -= (Hrs[0][1] * tau1)*D1[i] * D2[j] + (Hrs[1][0] * tau1)*D2[i] * D1[j];
							kij += Tb11[i] * D1[j] + Tb21[i] * D2[j];
							kij += D1[i] * Tb11[j] + D2[i] * Tb21[j];
							kij += g*(N11[i] * D1[j] + N12[i] * D2[j] + D1[i] * N11[j] + D2[i] * N12[j]);
							kij -= N[i] * Nb1[j] - Nb1[i] * N[j];
							kij -= T1[i] * Mi[0][0] * Tb11[j] + T1[i] * Mi[0][1] * Tb21[j] + T2[i] * Mi[1][0] * Tb11[j] + T2[i] * Mi[1][1] * Tb21[j];
							kij -= Tb11[i] * Mi[0][0] * T1[j] + Tb21[i] * Mi[0][1] * T1[j] + Tb11[i] * Mi[1][0] * T2[j] + Tb21[i] * Mi[1][1] * T2[j];
							
							ke[i][j] += ((tt*(Utxt)* Ac[0][0] + tt*(Utyt)* Ac[1][0])*kij)*detJ[j] * w[j] * m_ktmult;

							// KT2
							kij = T21[i] * D1[j] + T22[i] * D2[j];
							kij += D1[i] * T21[j] + D2[i] * T22[j];
							kij -= (Hrs[0][1] * tau2)*D1[i] * D2[j] + (Hrs[1][0] * tau2)*D2[i] * D1[j];
							kij += Tb12[i] * D1[j] + Tb22[i] * D2[j];
							kij += D1[i] * Tb12[j] + D2[i] * Tb22[j];
							kij += g*(N21[i] * D1[j] + N22[i] * D2[j] + D1[i] * N21[j] + D2[i] * N22[j]);
							kij -= N[i] * Nb2[j] - Nb2[i] * N[j];
							kij -= T1[i] * Mi[0][0] * Tb12[j] + T1[i] * Mi[0][1] * Tb22[j] + T2[i] * Mi[1][0] * Tb12[j] + T2[i] * Mi[1][1] * Tb22[j];
							kij -= Tb12[i] * Mi[0][0] * T1[j] + Tb22[i] * Mi[0][1] * T1[j] + Tb12[i] * Mi[1][0] * T2[j] + Tb22[i] * Mi[1][1] * T2[j];

							ke[i][j] += ((tt*(Utxt)* Ac[0][1] + tt*(Utyt)* Ac[1][1])*kij)*detJ[j] * w[j] * m_ktmult;

							// kdirect
							//if (bstick)
							//{
							//kij = Mk[0][0] * D1[i] * D1[j] + Mk[0][1] * D1[i] * D2[j] + Mk[1][0] * D2[i] * D1[j] + Mk[1][1] * D2[i] * D2[j];
							//ke[i][j] += m_epsf*kij;
							//}

							}
					


						for (int k=0; k<ndof; ++k)
							for (int l=0; l<ndof; ++l)
							{
							// 
								sum = Mi[0][0]*Nb1[k]*Nb1[l]+Mi[0][1]*(Nb1[k]*Nb2[l]+Nb2[k]*Nb1[l])+Mi[1][1]*Nb2[k]*Nb2[l];
								sum *= g;
								sum -= D1[k]*N1[l]+D2[k]*N2[l]+N1[k]*D1[l]+N2[k]*D2[l];
								sum += K[0][1]*(D1[k]*D2[l]+D2[k]*D1[l]);
								sum *= tntotal*knmult;

								ke[k][l] += sum*detJ[j]*w[j];
							}
					}
					//test begin
					if (tn == 0){
						//ss.GetMesh()->Node(se.m_node[j]).m_pt = 0;

						//ss.GetMesh()->Node(me.m_node[j]).m_at.x = -ss.GetMesh()->Node(me.m_node[j]).m_ap.x;
						//ss.GetMesh()->Node(me.m_node[j]).m_at.y = -ss.GetMesh()->Node(me.m_node[j]).m_ap.y;
						//ss.GetMesh()->Node(me.m_node[j]).m_at.z = -ss.GetMesh()->Node(me.m_node[j]).m_ap.z;



					}
					else{

						// assemble the global residual
						psolver->AssembleStiffness(en, LM, ke);
					}
					//test end
				}
			}
		}
	}
}


void FEFacet2FacetSliding::MapTangentialComponent(int inode, FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, FESurfaceElement &en, FESurfaceElement &eo, vec3d &q)
{
	
	FESurfaceElement& se = ss.Element(inode);
	//// tangential projection
	// onto the master element
	double rx = ss.m_rs_n[se.m_nID][inode][0];
	double sx = ss.m_rs_n[se.m_nID][inode][1];

	// get the coordinates at the previous step
	double rp = ss.m_rsp_n[se.m_nID][inode][0], ro = rp;
	double sp = ss.m_rsp_n[se.m_nID][inode][1], so = sp;
	
	
	
	// first we find the projection of the old point on the new segment
	//double r = ss.m_rs[inode][0];
	//double s = ss.m_rs[inode][1];
	//double rp = ss.m_rsp[inode][0], ro = rp;
	//double sp = ss.m_rsp[inode][1], so = sp;
	vec3d xn = ms.Local2Global(eo, rp, sp);
	vec3d qn;
	qn = ms.ProjectToSurface(en, xn, rp, sp);
	ss.m_rsp_n[se.m_nID][inode][0] = rp;
	ss.m_rsp_n[se.m_nID][inode][1]= sp;
	//ss.m_rsp[inode][0] = rp;
	///ss.m_rsp[inode][1] = sp;

	// next, we transform the frictional traction
	// since these tractions are decomposed in the local 
	// element coordinate system, we have to do a coordinate transformation
	// note that this transformation needs to be done in curvilinear
	// coordinates since the base vectors may not be orthonormal. Also
	// note that we are doing this in the reference configuration
	vec3d to[2], tn[2];
	ms.ContraBaseVectors0(eo, ro, so, to);
	ms.CoBaseVectors0(en, rx, sx, tn);

	double Lt[2];

	Lt[0] = ss.m_Lnt[se.m_nID][inode][0];
	Lt[1] = ss.m_Lnt[se.m_nID][inode][1];

	//Lt[0] = ss.m_Lt[inode][0];
	//Lt[1] = ss.m_Lt[inode][1];

	vec3d t;
	t = to[0] * Lt[0] + to[1] * Lt[1];

	Lt[0] = t*tn[0];
	Lt[1] = t*tn[1];

	ss.m_Lnt[se.m_nID][inode][0] = Lt[0];
	ss.m_Lnt[se.m_nID][inode][1] = Lt[1];
	//ss.m_Lt[inode][0] = Lt[0];
	//ss.m_Lt[inode][1] = Lt[1];
}



//-----------------------------------------------------------------------------
// --- N O R M A L + T A N G E N T I A L   Pressures
void FEFacet2FacetSliding::UpdateContactPressures()
{
	int npass = (m_btwo_pass?2:1);
	//0.0148933197058=gap(infinitesimal)[δgap]/search tolerance
	//double tolerance = 0.0148933197058*m_stol;
	double tolerance = m_stol*5.5988986099716 * (10 ^ (-10));
	double infgap = (5.551115 / 151)*(exp(-18));

	const int MELN = FEElement::MAX_NODES;  //const int MAXMN = FEElement::MAX_NODES;
	//double detJ[MELN], w[MELN], *Hs, H[MELN], Hm[MELN];  //double H[MAXMN], Hr[MAXMN], Hs[MAXMN];
	// contact vectors
	//double N[3 * (MELN + 1)], N1[3 * (MELN + 1)], N2[3 * (MELN + 1)];
	//double T1[3 * (MELN + 1)], T2[3 * (MELN + 1)], D1[3 * (MELN + 1)], D2[3 * (MELN + 1)];

	//vec3d r0[MELN];  //vec3d rtm[MAXMN];



	for (int np=0; np<npass; ++np)
	{
		FEFacetSlidingSurface& ss = (np == 0? m_ss : m_ms);
		FEFacetSlidingSurface& ms = (np == 0? m_ms : m_ss);
		
		// loop over all elements of the primary surface
		for (int n=0; n<ss.Elements(); ++n)
		{
			FESurfaceElement& el = ss.Element(n);
			int elemn = el.m_nID;
			int nint = el.GaussPoints();
			
			// get the normal tractions at the integration points
			double gap, eps;
			//mat2d Mk, Mki;
			for (int i=0; i<nint; ++i) 
			{
				FEFacetSlidingSurface::Data& pt = ss.m_Data[n][i];

				gap = pt.m_gap;
				eps = m_epsn*pt.m_eps;
				if (eps == 0){
					felog.printf("Division by zero!....................penalty value==0");
					break;
				}

				/////////////////////////////////

				mat2d& Mk=ss.m_M_t[el.m_nID][i];
				mat2d Mki= Mk.inverse();

				// get the master element node positions
				//for (k = 0; k<nmeln; ++k) rtm[k] = mesh.Node(mel.m_node[k]).m_rt;

				// isoparametric coordinates of the projected slave node
				// onto the master element
				 double rx = ss.m_rs_n[el.m_nID][i][0];
				 double sx = ss.m_rs_n[el.m_nID][i][1];

				// get the coordinates at the previous step
				 double rp = ss.m_rsp_n[el.m_nID][i][0];
				 double sp = ss.m_rsp_n[el.m_nID][i][1];

				double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
				double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
				double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
				double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
				//cut here !

				int maxcohetraction = m_naugmax;
				// the directional derivative of traction-displacement function in the direction of search tolerance with respect to initial tolerance
				double cohesivetolerance = (tolerance*((maxcohetraction / eps)*exp(1 - (tolerance / eps))) + ((tolerance / eps) * maxcohetraction*(-1 / eps)*exp(1 - (tolerance / eps))));
				double crtl = ((infgap / eps) * maxcohetraction * exp(1 - (infgap / eps)));
				if (gap < 0){ gap = -gap; }
				double tracs1 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utxn;
				double tracs2 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utyn;
				double tracs = tracs1 + tracs2;
				double tractt1 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utxt;
				double tractt2 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utyt;
				double tractt = tractt1 + tractt2;
				if ((gap> rootcontrol*eps) && ((-tracs) < crtl)){
					tracs = 0;
				}
				ss.m_Lnt[el.m_nID][i][0] = -(MBRACKET(-tractt1));
				ss.m_Lnt[el.m_nID][i][1] = -(MBRACKET(-tractt2));
				pt.m_Ln = -(MBRACKET(-tracs));
				felog.printf("pt.m_Ln:  %u \n", pt.m_Ln);
				//pt.m_Ln = tracs;
				/*smart contact part
				double tn = -((g / eps) * maxcohetraction * exp(1 - (g / eps)));


				pt.m_Ln = MBRACKET(pt.m_Lm + eps*gap);
				ti[j] = MBRACKET(md.m_Lm + m_epsn*md.m_eps*md.m_gap);

				*/
				double tracm = 0;
				
				/*
				if ((-gap)> rootcontrol*eps && (-tracs) < cohesivetolerance){
					felog.printf("smart part  back to original....................part UpdateContactPressures");
					pt.m_Ln = MBRACKET(pt.m_Lm + eps*gap);
				}

				*/
				//felog.printf("test ....................part UpdateContactPressures");
				FESurfaceElement* pme = pt.m_pme;
				if (m_btwo_pass && pme)
				{
					// get master element data
					vector<FEFacetSlidingSurface::Data>& mdv = ms.m_Data[pme->m_lid];
					int mint = pme->GaussPoints();
					double ti[FEElement::MAX_NODES], tt1[FEElement::MAX_NODES], tt2[FEElement::MAX_NODES];
					for (int j=0; j<mint; ++j) 
					{
						FEFacetSlidingSurface::Data& md = mdv[j];
						gap = md.m_gap;
						eps = m_epsn*md.m_eps;
						mat2d& Mk = ms.m_M_t[pme->m_nID][i];
						mat2d Mki = Mk.inverse();
						double rx = ms.m_rs_n[pme->m_nID][i][0];
						double sx = ms.m_rs_n[pme->m_nID][i][1];

						// get the coordinates at the previous step
						double rp = ms.m_rsp_n[pme->m_nID][i][0];
						double sp = ms.m_rsp_n[pme->m_nID][i][1];

						double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
						double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
						double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
						double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
						//cut here !

						if (eps == 0){
							felog.printf("Division by zero!....................penalty value==0");
							break;
							exit(1);
						}
						if (md.m_gap < 0){ md.m_gap = -md.m_gap; }
						 //tracm = -((md.m_gap / m_epsn*md.m_eps) * maxcohetraction * exp(1 - (md.m_gap / m_epsn*md.m_eps)));
						double tracs1 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utxn;
						double tracs2 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utyn;
						double tracm = tracs1 + tracs2;
						double tracttm1 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utxt;
						double tracttm2 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utyt;
						double tractt = tractt1 + tractt2;

						 if ((gap> rootcontrol*eps) && ((-tracm) < crtl)){
							 tracm = 0;
						 }
							
						 tt1[i] = -(MBRACKET(-tracttm1));
						 tt2[i] = -(MBRACKET(-tracttm2));
						 ti[j] = -(MBRACKET(-tracm));
						 //ti[j] = tracm;
							/*
							if ((-gap)> rootcontrol*eps && (-tracm) < cohesivetolerance){
							felog.printf("smart part  back to original....................part two_pass");
							ti[j] = MBRACKET(md.m_Lm + m_epsn*md.m_eps*md.m_gap);
						}
							*/
					}
					// project the data to the nodes
					double tn[FEElement::MAX_NODES], ttx[FEElement::MAX_NODES], tty[FEElement::MAX_NODES];
					pme->project_to_nodes(ti, tn);
					pme->project_to_nodes(tt1, ttx);
					pme->project_to_nodes(tt2, tty);
					//double rx = ss.m_rs_n[el.m_nID][i][0];
					//double sx = ss.m_rs_n[el.m_nID][i][1];
					// now evaluate the traction at the intersection point
					double Ln = pme->eval(tn, pt.m_rs[0], pt.m_rs[1]);

					double L1 = pme->eval(ttx, ss.m_rs_n[el.m_nID][i][0], ss.m_rs_n[el.m_nID][i][1]);
					double L2 = pme->eval(tty, ss.m_rs_n[el.m_nID][i][0], ss.m_rs_n[el.m_nID][i][1]);
					(ss.m_Lnt[el.m_nID][i][0]) += -(MBRACKET(-L1));
					(ss.m_Lnt[el.m_nID][i][1]) += -(MBRACKET(-L2));
					pt.m_Ln += -(MBRACKET(-Ln));
					
					//pt.m_Ln += Ln;
					/*
					if ((-gap)> rootcontrol*eps && (-tracm) < cohesivetolerance){
						felog.printf("smart part  back to original....................part two_pass intersection point");
						pt.m_Ln += MBRACKET(Ln);
					}  */
				}
				felog.printf("pt.m_Ln if two pass:  %u \n", pt.m_Ln);
				// test UpdateContactPressures
				//if (pt.m_Ln < 0){ pt.m_Ln = 0;}
			}
		}
	}
	
}

//-----------------------------------------------------------------------------
bool FEFacet2FacetSliding::Augment(int naug)
{
	// make sure we need to augment
	if (!m_blaugon) return true;
	int maxcohetraction = m_naugmax;
	bool bconv = true;
	double normL0 = 0;
	// --- c a l c u l a t e   i n i t i a l   n o r m s ---
	// a. normal component
	int NS = (int) m_ss.m_Data.size();
	int NM = (int) m_ms.m_Data.size();
	double Ln;
	double Lt[2];
	//double scale = Penalty();
	////////////////
	// b. tangential component

	// loop over all elements of the master and slave  surfaces
	for (int n = 0; n<m_ss.Elements(); ++n)
	{
		FESurfaceElement& el = m_ss.Element(n);
		int elemn = el.m_nID;
		int nint = el.GaussPoints();
		for (int i = 0; i < nint; ++i)
		{
			Lt[0] = m_ss.m_Lnt[el.m_nID][i][0];
			Lt[1] = m_ss.m_Lnt[el.m_nID][i][1];
			mat2d& Mk = m_ss.m_M_t[el.m_nID][i];
			mat2d Mi = Mk.inverse();
			normL0 += Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
		}
		}

	for (int n = 0; n<m_ms.Elements(); ++n)
	{
		FESurfaceElement& el = m_ms.Element(n);
		int elemn = el.m_nID;
		int nint = el.GaussPoints();
		for (int i = 0; i < nint; ++i)
		{
			Lt[0] = m_ms.m_Lnt[el.m_nID][i][0];
			Lt[1] = m_ms.m_Lnt[el.m_nID][i][1];
			mat2d& Mk = m_ms.m_M_t[el.m_nID][i];
			mat2d Mi = Mk.inverse();
			normL0 += Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]);
		}
	}

	////////////////

	
	for (int i=0; i<NS; ++i)
	{
		vector<FEFacetSlidingSurface::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FEFacetSlidingSurface::Data& ds = sd[j];
			normL0 += ds.m_Lm*ds.m_Lm;
		}
	}
	for (int i=0; i<NM; ++i)
	{
		vector<FEFacetSlidingSurface::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FEFacetSlidingSurface::Data& dm = md[j];
			normL0 += dm.m_Lm*dm.m_Lm;
		}
	}
	normL0 = sqrt(normL0);

	// --- c a l c u l a t e   c u r r e n t   n o r m s ---
	// a. normal component
	double normL1 = 0;	// force norm
	double normg1 = 0;	// gap norm
	int N = 0;


	////////////////
	// b. tangential component

	// loop over all elements of the master and slave  surfaces
	for (int n = 0; n<m_ss.Elements(); ++n)
	{
		FESurfaceElement& el = m_ss.Element(n);
		int elemn = el.m_nID;
		int nint = el.GaussPoints();
		for (int i = 0; i < nint; ++i)
		{
			FEFacetSlidingSurface::Data& pt = m_ss.m_Data[n][i];

			double gap = pt.m_gap;
			double eps = m_epsn*pt.m_eps;
			double rx = m_ss.m_rs_n[el.m_nID][i][0];
			double sx = m_ss.m_rs_n[el.m_nID][i][1];
			mat2d& Mk = m_ss.m_M_t[el.m_nID][i];
			mat2d Mi = Mk.inverse();
			// get the coordinates at the previous step
			double rp = m_ss.m_rsp_n[el.m_nID][i][0];
			double sp = m_ss.m_rsp_n[el.m_nID][i][1];

			double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
			double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
			double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
			double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
			//cut here !
			double tracttm1 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utxt;
			double tracttm2 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utyt;
			Lt[0] = -(MBRACKET(-tracttm1));
			Lt[1] = -(MBRACKET(-tracttm2));

			normL1 += -(Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]));

			//eps = m_epsn*pt.m_eps*scale;
			double tracs1 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utxn;
			double tracs2 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utyn;
			double tracs = tracs1 + tracs2;
			// update Lagrange multipliers
			//Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
			Ln = -(MBRACKET(-tracs));

			normL1 += Ln*Ln;

			if (pt.m_gap > 0)
			{
				normg1 += pt.m_gap * pt.m_gap;
				++N;
			}
		}

	}

	for (int n = 0; n<m_ms.Elements(); ++n)
	{
		FESurfaceElement& el = m_ms.Element(n);
		int elemn = el.m_nID;
		int nint = el.GaussPoints();
		for (int i = 0; i < nint; ++i)
		{
			FEFacetSlidingSurface::Data& pt = m_ms.m_Data[n][i];

			double gap = pt.m_gap;
			double eps = m_epsn*pt.m_eps;
			double rx = m_ms.m_rs_n[el.m_nID][i][0];
			double sx = m_ms.m_rs_n[el.m_nID][i][1];
			mat2d& Mk = m_ms.m_M_t[el.m_nID][i];
			mat2d Mi = Mk.inverse();
			// get the coordinates at the previous step
			double rp = m_ms.m_rsp_n[el.m_nID][i][0];
			double sp = m_ms.m_rsp_n[el.m_nID][i][1];

			double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
			double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
			double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
			double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
			//cut here !
			double tracttm1 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utxt;
			double tracttm2 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utyt;
			Lt[0] = -(MBRACKET(-tracttm1));
			Lt[1] = -(MBRACKET(-tracttm2));

			normL1 += -(Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]));

			//eps = m_epsn*pt.m_eps*scale;
			double tracs1 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utxn;
			double tracs2 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utyn;
			double tracs = tracs1 + tracs2;
			// update Lagrange multipliers
			//Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
			Ln = -(MBRACKET(-tracs));

			normL1 += Ln*Ln;

			if (pt.m_gap > 0)
			{
				normg1 += pt.m_gap * pt.m_gap;
				++N;
			}


		}
	}
	////////////////


	/*
	for (int i=0; i<NS; ++i)
	{
		vector<FEFacetSlidingSurface::Data>& sd = m_ss.m_Data[i];
		for (int j=0; j<(int)sd.size(); ++j)
		{
			FEFacetSlidingSurface::Data& ds = sd[j];

			// penalty value
			double eps = m_epsn*ds.m_eps;

			// update Lagrange multipliers
			double Ln = ds.m_Lm + eps*ds.m_gap;
			Ln = MBRACKET(Ln);

			normL1 += Ln*Ln;

			if (ds.m_gap > 0)
			{
				normg1 += ds.m_gap*ds.m_gap;
				++N;
			}
		}
	}	

	for (int i=0; i<NM; ++i)
	{
		vector<FEFacetSlidingSurface::Data>& md = m_ms.m_Data[i];
		for (int j=0; j<(int)md.size(); ++j)
		{
			FEFacetSlidingSurface::Data& dm = md[j];

			// penalty value
			double eps = m_epsn*dm.m_eps;

			// update Lagrange multipliers
			double Ln = dm.m_Lm + eps*dm.m_gap;
			Ln = MBRACKET(Ln);

			normL1 += Ln*Ln;
			if (dm.m_gap > 0)
			{
				normg1 += dm.m_gap*dm.m_gap;
				++N;
			}
		}
	}


	*/

	if (N == 0) N=1;

	normL1 = sqrt(normL1);
	normg1 = sqrt(normg1 / N);

	if (naug == 0) m_normg0 = 0;

	// calculate and print convergence norms
	double lnorm = 0, gnorm = 0;
	if (normL1 != 0) lnorm = fabs(normL1 - normL0)/normL1; else lnorm = fabs(normL1 - normL0);
	if (normg1 != 0) gnorm = fabs(normg1 - m_normg0)/normg1; else gnorm = fabs(normg1 - m_normg0);

	felog.printf(" sliding interface # %d\n", m_nID);
	felog.printf("                        CURRENT        REQUIRED\n");
	felog.printf("    normal force : %15le", lnorm);
	if (m_atol > 0) felog.printf("%15le\n", m_atol); else felog.printf("       ***\n");
	felog.printf("    gap function : %15le", gnorm);
	if (m_gtol > 0) felog.printf("%15le\n", m_gtol); else felog.printf("       ***\n");

	// check convergence
	bconv = true;
	if ((m_atol > 0) && (lnorm > m_atol)) bconv = false;
	if ((m_gtol > 0) && (gnorm > m_gtol)) bconv = false;
	
	//if (m_naugminmod > naug) bconv = false;

	if (m_naugmin > naug) bconv = false;

	if (m_naugmaxmod <= naug) bconv = true;
	//if (m_naugmax <= naug) bconv = true;
		
	if (bconv == false)
	{
		// we did not converge so update multipliers


		for (int n = 0; n < m_ss.Elements(); ++n)
		{
			FESurfaceElement& el = m_ss.Element(n);
			int elemn = el.m_nID;
			int nint = el.GaussPoints();
			for (int i = 0; i < nint; ++i)
			{
				FEFacetSlidingSurface::Data& pt = m_ss.m_Data[n][i];

				double gap = pt.m_gap;
				double eps = m_epsn*pt.m_eps;
				double rx = m_ss.m_rs_n[el.m_nID][i][0];
				double sx = m_ss.m_rs_n[el.m_nID][i][1];
				mat2d& Mk = m_ss.m_M_t[el.m_nID][i];
				mat2d Mi = Mk.inverse();
				// get the coordinates at the previous step
				double rp = m_ss.m_rsp_n[el.m_nID][i][0];
				double sp = m_ss.m_rsp_n[el.m_nID][i][1];

				double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
				double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
				double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
				double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
				//cut here !
				double tracttm1 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utxt;
				double tracttm2 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utyt;
				Lt[0] = -(MBRACKET(-tracttm1));
				Lt[1] = -(MBRACKET(-tracttm2));

				//normL1 += -(Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]));

				//eps = m_epsn*pt.m_eps*scale;
				double tracs1 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utxn;
				double tracs2 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utyn;
				double tracs = tracs1 + tracs2;
				// update Lagrange multipliers
				//Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
				//Ln = -(MBRACKET(-tracs));

				m_ss.m_Lnt[el.m_nID][i][0] = Lt[0];
				m_ss.m_Lnt[el.m_nID][i][1] = Lt[1];
				pt.m_Ln = -(MBRACKET(-tracs));
			}

		}

		for (int n = 0; n < m_ms.Elements(); ++n)
		{
			FESurfaceElement& el = m_ms.Element(n);
			int elemn = el.m_nID;
			int nint = el.GaussPoints();
			for (int i = 0; i < nint; ++i)
			{
				FEFacetSlidingSurface::Data& pt = m_ms.m_Data[n][i];

				double gap = pt.m_gap;
				double eps = m_epsn*pt.m_eps;
				double rx = m_ms.m_rs_n[el.m_nID][i][0];
				double sx = m_ms.m_rs_n[el.m_nID][i][1];
				mat2d& Mk = m_ms.m_M_t[el.m_nID][i];
				mat2d Mi = Mk.inverse();
				// get the coordinates at the previous step
				double rp = m_ms.m_rsp_n[el.m_nID][i][0];
				double sp = m_ms.m_rsp_n[el.m_nID][i][1];

				double Utxn = exp(-((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp))) / (eps*eps));
				double Utyn = exp(-((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp))) / (eps*eps));
				double Utxt = (exp(-(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))*((Mk[0][0] * (rx - rp) + Mk[0][1] * (sx - sp)))) / eps));
				double Utyt = (exp(-(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / (eps*eps)))*((sqrt(((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))*((Mk[1][0] * (rx - rp) + Mk[1][1] * (sx - sp)))) / eps));
				//cut here !
				double tracttm1 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utxt;
				double tracttm2 = -(2 * sqrt(0.5*exp(1))*maxcohetraction *(1 + (gap / eps))*exp(-(gap / eps)))*Utyt;
				Lt[0] = -(MBRACKET(-tracttm1));
				Lt[1] = -(MBRACKET(-tracttm2));

				//normL1 += -(Lt[0] * (Mi[0][0] * Lt[0] + Mi[0][1] * Lt[1]) + Lt[1] * (Mi[1][0] * Lt[0] + Mi[1][1] * Lt[1]));

				//eps = m_epsn*pt.m_eps*scale;
				double tracs1 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utxn;
				double tracs2 = -((gap / eps) * maxcohetraction * exp(1 - (gap / eps)))*Utyn;
				double tracs = tracs1 + tracs2;
				// update Lagrange multipliers
				//Ln = m_ss.m_Lm[i] + eps*m_ss.m_gap[i];
				//Ln = -(MBRACKET(-tracs));

				m_ms.m_Lnt[el.m_nID][i][0] = Lt[0];
				m_ms.m_Lnt[el.m_nID][i][1] = Lt[1];
				pt.m_Ln = -(MBRACKET(-tracs));
			}

		}
		/*
		for (int i=0; i<NS; ++i)
		{
			vector<FEFacetSlidingSurface::Data>& sd = m_ss.m_Data[i];
			for (int j=0; j<(int)sd.size(); ++j)
			{
				FEFacetSlidingSurface::Data& ds = sd[j];

				// penalty value
				double eps = m_epsn*ds.m_eps;

				// update Lagrange multipliers
				double Ln = ds.m_Lm + eps*ds.m_gap;
				ds.m_Lm = MBRACKET(Ln);
			}
		}	

		for (int i=0; i<NM; ++i)
		{
			vector<FEFacetSlidingSurface::Data>& md = m_ms.m_Data[i];
			for (int j=0; j<(int)md.size(); ++j)
			{
				FEFacetSlidingSurface::Data& dm = md[j];

				// penalty value
				double eps = m_epsn*dm.m_eps;

				// update Lagrange multipliers
				double Ln = dm.m_Lm + eps*dm.m_gap;
				dm.m_Lm = MBRACKET(Ln);
			}
		}
		*/
	}

	// store the last gap norm
	m_normg0 = normg1;

	return bconv;
}

//-----------------------------------------------------------------------------
void FEFacet2FacetSliding::Serialize(DumpFile &ar)
{
	// store contact data
	FEContactInterface::Serialize(ar);

	// store contact surface data
	m_ms.Serialize(ar);
	m_ss.Serialize(ar);
}
