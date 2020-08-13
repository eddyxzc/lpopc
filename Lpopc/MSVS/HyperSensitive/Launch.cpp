// Copyright (C) 2014-2015 Xue Zhichen, Wang Yujie,Wang Na 
// All Rights Reserved.
// This file is a part of LPOPC , published under the Eclipse Public License.
// Author:Xue Zhichen 7/14 2015   13:16
// Email:eddy_lpopc@163.com
// 
#include "Launch.hpp"
#include "LpLpopcApplication.hpp"

const double PI= datum::pi;
using namespace Lpopc;
 void Launchoe2rv( vec& oe,double mu,vec& ri,vec& vi);
 void Launchrv2oe(vec& oe, vec rv, vec vv, double mu);
 void DLaunchrv2oe(const vec& rv,const vec& vv,double mu,double Re,mat& doe);
 double earthRadius         = 6378145;
	double gravParam           = 3.986012e14;
	double initialMass         = 301454;
	double earthRotRate        = 7.29211585e-5;
	double seaLevelDensity     = 1.225;
	double densityScaleHeight  = 7200;
	double g0                  = 9.80665;

	struct struct_scales 
	{
		double length;
		double speed;
		double time;
		double acceleration;
		double mass;
		double force;
		double area;
		double volume;
		double density;
		double gravparam;
	}scales={
		 earthRadius,
		sqrt(gravParam/scales.length),
		scales.length/scales.speed,
		 scales.speed/scales.time,
		 initialMass,
		scales.mass*scales.acceleration,
		 scales.length*scales.length,
		scales.area*scales.length,
		scales.mass/scales.volume,
		scales.acceleration*scales.length*scales.length
	};

	double omega=earthRotRate*scales.time;

	struct CONSTANTS_struct
	{
		double omega_matrix[9];
		double mu;
		double cd;
		double sa;
		double rho0;
		double H;
		double Re;
		double g0;
		double thrust_srb;
		double thrust_first;
		double thrust_second;
		double ISP_srb;
		double ISP_first;
		double ISP_second;
	}CONSTANTS={
		{0,1*omega,0,-1*omega,0,0,0,0},
		gravParam/scales.gravparam,
		0.5,
		4*PI/scales.area,
		seaLevelDensity/scales.density,
		densityScaleHeight/scales.length,
		earthRadius/scales.length,g0/scales.acceleration
	};
	
int main(void)
{
	double lat0 = 28.5*PI/180;               // Geocentric Latitude of Cape Canaveral
	double x0 = CONSTANTS.Re* cos (lat0);     // x component of initial position
	double z0 = CONSTANTS.Re* sin (lat0);      // z component of initial position
	double 	y0 = 0;

	vec r0(3, 1);
	r0(0) = x0; r0(1) = y0;r0(2) =z0;
	
	mat omega_matrix(CONSTANTS.omega_matrix,3,3);

	mat v0 = omega_matrix*r0;

	double bt_srb = 75.2/scales.time;
	double bt_first = 261.0/scales.time;
	double bt_second = 700.0/scales.time;

	double t0 = 0.0/scales.time;
	double t1 = 75.2/scales.time;
	double t2 = 150.4/scales.time;
	double t3 = 261/scales.time;
	double t4 = 961/scales.time;

	double m_tot_srb     = 19290/scales.mass;
	double m_prop_srb    = 17010/scales.mass;
	double m_dry_srb     = m_tot_srb-m_prop_srb;
	double m_tot_first   = 104380/scales.mass;
	double m_prop_first  = 95550/scales.mass;
	double m_dry_first   = m_tot_first-m_prop_first;
	double m_tot_second  = 19300/scales.mass;
	double m_prop_second = 16820/scales.mass;
	double m_dry_second  = m_tot_second-m_prop_second;
	double m_payload     = 4164/scales.mass;
	double thrust_srb    = 628500/scales.force;
	double thrust_first  = 1083100/scales.force;
	double thrust_second = 110094/scales.force;
	double mdot_srb      = m_prop_srb/bt_srb;
	double ISP_srb       = thrust_srb/(CONSTANTS.g0*mdot_srb);
	double mdot_first    = m_prop_first/bt_first;
	double ISP_first     = thrust_first/(CONSTANTS.g0*mdot_first);
	double mdot_second   = m_prop_second/bt_second;
	double ISP_second     = thrust_second/(CONSTANTS.g0*mdot_second);

	double af = 24361140/scales.length;
	double ef = 0.7308;
	double incf = 28.5*PI/180;
	double Omf = 269.8*PI/180;
	double omf = 130.5*PI/180;
	double nuguess = 0;
	double cosincf = cos(incf);
	double cosOmf = cos(Omf);
	double cosomf = cos(omf);

	vec oe_matrix(6);
	oe_matrix[0] = af;
	oe_matrix[1] = ef;
	oe_matrix[2] = incf;
	oe_matrix[3] = Omf;
	oe_matrix[4] = omf;
	oe_matrix[5] = nuguess;
	vec vout,rout;
	Launchoe2rv(oe_matrix,CONSTANTS.mu,rout,vout);
	double 	m10 = m_payload+m_tot_second+m_tot_first+9*m_tot_srb;
	double m1f = m10-(6*mdot_srb+mdot_first)*t1;
	double m20 = m1f-6*m_dry_srb;
	double m2f = m20-(3*mdot_srb+mdot_first)*(t2-t1);
	double m30 = m2f-3*m_dry_srb;
	double m3f = m30-mdot_first*(t3-t2);
	double m40 = m3f-m_dry_first;
	double m4f = m_payload;

	CONSTANTS.thrust_srb    = thrust_srb;
	CONSTANTS.thrust_first  = thrust_first;
	CONSTANTS.thrust_second = thrust_second;
	CONSTANTS.ISP_srb       = ISP_srb;
	CONSTANTS.ISP_first     = ISP_first;
	CONSTANTS.ISP_second    = ISP_second;

	double rmin = -2*CONSTANTS.Re;
	double rmax = -rmin;
	double vmin = -10000/scales.speed;
	double vmax = -vmin;
	//Step 1
	shared_ptr<LpopcApplication> app(new LpopcApplication(console_print));

	//Step2
	shared_ptr<Phase> Phase1 (new Phase(1,7,3,0,1,0));
	Phase1->SetTimeMin(t0,t1);
	Phase1->SetTimeMax(t0,t1);
	Phase1->SetStateMin(r0(0),rmin,rmin);
	Phase1->SetStateMax(r0(0),rmax,rmax);

	Phase1->SetStateMin(r0(1),rmin,rmin);
	Phase1->SetStateMax(r0(1),rmax,rmax);

	Phase1->SetStateMin(r0(2),rmin,rmin);
	Phase1->SetStateMax(r0(2),rmax,rmax);

	Phase1->SetStateMin(v0(0),vmin,vmin);
	Phase1->SetStateMax(v0(0),vmax,vmax);

	Phase1->SetStateMin(v0(1),vmin,vmin);
	Phase1->SetStateMax(v0(1),vmax,vmax);

	Phase1->SetStateMin(v0(2),vmin,vmin);
	Phase1->SetStateMax(v0(2),vmax,vmax);

	Phase1->SetStateMin(m10,m1f,m1f);
	Phase1->SetStateMax(m10,m10,m10);

	Phase1->SetcontrolMin(-1);
	Phase1->SetcontrolMax(1);

	Phase1->SetcontrolMin(-1);
	Phase1->SetcontrolMax(1);

	Phase1->SetcontrolMin(-1);
	Phase1->SetcontrolMax(1);
	//No parameter

	Phase1->SetpathMin(1);
	Phase1->SetpathMax(1);

	Phase1->SetTimeGuess(t0);
	Phase1->SetTimeGuess(t1);

	Phase1->SetStateGuess(1,r0(0));
	Phase1->SetStateGuess(1,r0(0));

	Phase1->SetStateGuess(2,r0(1));
	Phase1->SetStateGuess(2,r0(1));

	Phase1->SetStateGuess(3,r0(2));
	Phase1->SetStateGuess(3,r0(2));

	Phase1->SetStateGuess(4,v0(0));
	Phase1->SetStateGuess(4,v0(0));

	Phase1->SetStateGuess(5,v0(1));
	Phase1->SetStateGuess(5,v0(1));

	Phase1->SetStateGuess(6,v0(2));
	Phase1->SetStateGuess(6,v0(2));

	Phase1->SetStateGuess(7,m10);
	Phase1->SetStateGuess(7,m1f);

	Phase1->SetControlGuess(1,0);
	Phase1->SetControlGuess(1,0);

	Phase1->SetControlGuess(2,1);
	Phase1->SetControlGuess(2,1);

	Phase1->SetControlGuess(3,0);
	Phase1->SetControlGuess(3,0);
	//Phase 2
	shared_ptr<Phase> Phase2 (new Phase(2,7,3,0,1,0));
	Phase2->SetTimeMin(t1,t2);
	Phase2->SetTimeMax(t1,t2);

	Phase2->SetStateMin(rmin,rmin,rmin);
	Phase2->SetStateMax(rmax,rmax,rmax);

	Phase2->SetStateMin(rmin,rmin,rmin);
	Phase2->SetStateMax(rmax,rmax,rmax);

	Phase2->SetStateMin(rmin,rmin,rmin);
	Phase2->SetStateMax(rmax,rmax,rmax);

	Phase2->SetStateMin(vmin,vmin,vmin);
	Phase2->SetStateMax(vmax,vmax,vmax);

	Phase2->SetStateMin(vmin,vmin,vmin);
	Phase2->SetStateMax(vmax,vmax,vmax);

	Phase2->SetStateMin(vmin,vmin,vmin);
	Phase2->SetStateMax(vmax,vmax,vmax);

	Phase2->SetStateMin(m2f,m2f,m2f);
	Phase2->SetStateMax(m20,m20,m20);

	Phase2->SetcontrolMin(-1);
	Phase2->SetcontrolMax(1);

	Phase2->SetcontrolMin(-1);
	Phase2->SetcontrolMax(1);

	Phase2->SetcontrolMin(-1);
	Phase2->SetcontrolMax(1);
	//No parameter

	Phase2->SetpathMin(1);
	Phase2->SetpathMax(1);

	Phase2->SetTimeGuess(t1);
	Phase2->SetTimeGuess(t2);

	Phase2->SetStateGuess(1,r0(0));
	Phase2->SetStateGuess(1,r0(0));

	Phase2->SetStateGuess(2,r0(1));
	Phase2->SetStateGuess(2,r0(1));

	Phase2->SetStateGuess(3,r0(2));
	Phase2->SetStateGuess(3,r0(2));

	Phase2->SetStateGuess(4,v0(0));
	Phase2->SetStateGuess(4,v0(0));

	Phase2->SetStateGuess(5,v0(1));
	Phase2->SetStateGuess(5,v0(1));

	Phase2->SetStateGuess(6,v0(2));
	Phase2->SetStateGuess(6,v0(2));

	Phase2->SetStateGuess(7,m20);
	Phase2->SetStateGuess(7,m2f);

	Phase2->SetControlGuess(1,0);
	Phase2->SetControlGuess(1,0);

	Phase2->SetControlGuess(2,1);
	Phase2->SetControlGuess(2,1);

	Phase2->SetControlGuess(3,0);
	Phase2->SetControlGuess(3,0);


	shared_ptr<Phase> Phase3 (new Phase(3,7,3,0,1,0));
	Phase3->SetTimeMin(t2,t3);
	Phase3->SetTimeMax(t2,t3);

	Phase3->SetStateMin(rmin,rmin,rmin);
	Phase3->SetStateMax(rmax,rmax,rmax);

	Phase3->SetStateMin(rmin,rmin,rmin);
	Phase3->SetStateMax(rmax,rmax,rmax);

	Phase3->SetStateMin(rmin,rmin,rmin);
	Phase3->SetStateMax(rmax,rmax,rmax);

	Phase3->SetStateMin(vmin,vmin,vmin);
	Phase3->SetStateMax(vmax,vmax,vmax);

	Phase3->SetStateMin(vmin,vmin,vmin);
	Phase3->SetStateMax(vmax,vmax,vmax);

	Phase3->SetStateMin(vmin,vmin,vmin);
	Phase3->SetStateMax(vmax,vmax,vmax);

	Phase3->SetStateMin(m3f,m3f,m3f);
	Phase3->SetStateMax(m30,m30,m30);

	Phase3->SetcontrolMin(-1);
	Phase3->SetcontrolMax(1);

	Phase3->SetcontrolMin(-1);
	Phase3->SetcontrolMax(1);

	Phase3->SetcontrolMin(-1);
	Phase3->SetcontrolMax(1);
	//No parameter

	Phase3->SetpathMin(1);
	Phase3->SetpathMax(1);

	Phase3->SetTimeGuess(t2);
	Phase3->SetTimeGuess(t3);

	Phase3->SetStateGuess(1,rout(0));
	Phase3->SetStateGuess(1,rout(0));

	Phase3->SetStateGuess(2,rout(1));
	Phase3->SetStateGuess(2,rout(1));

	Phase3->SetStateGuess(3,rout(2));
	Phase3->SetStateGuess(3,rout(2));

	Phase3->SetStateGuess(4,vout(0));
	Phase3->SetStateGuess(4,vout(0));

	Phase3->SetStateGuess(5,vout(1));
	Phase3->SetStateGuess(5,vout(1));

	Phase3->SetStateGuess(6,vout(2));
	Phase3->SetStateGuess(6,vout(2));

	Phase3->SetStateGuess(7,m30);
	Phase3->SetStateGuess(7,m3f);

	Phase3->SetControlGuess(1,0);
	Phase3->SetControlGuess(1,0);

	Phase3->SetControlGuess(2,1);
	Phase3->SetControlGuess(2,1);

	Phase3->SetControlGuess(3,0);
	Phase3->SetControlGuess(3,0);

	shared_ptr<Phase> Phase4 (new Phase(4,7,3,0,1,5));
	Phase4->SetTimeMin(t3,t3);
	Phase4->SetTimeMax(t3,t4);

	Phase4->SetStateMin(rmin,rmin,rmin);
	Phase4->SetStateMax(rmax,rmax,rmax);

	Phase4->SetStateMin(rmin,rmin,rmin);
	Phase4->SetStateMax(rmax,rmax,rmax);

	Phase4->SetStateMin(rmin,rmin,rmin);
	Phase4->SetStateMax(rmax,rmax,rmax);

	Phase4->SetStateMin(vmin,vmin,vmin);
	Phase4->SetStateMax(vmax,vmax,vmax);

	Phase4->SetStateMin(vmin,vmin,vmin);
	Phase4->SetStateMax(vmax,vmax,vmax);

	Phase4->SetStateMin(vmin,vmin,vmin);
	Phase4->SetStateMax(vmax,vmax,vmax);

	Phase4->SetStateMin(m4f,m4f,m4f);
	Phase4->SetStateMax(m40,m40,m40);

	Phase4->SetcontrolMin(-1);
	Phase4->SetcontrolMax(1);

	Phase4->SetcontrolMin(-1);
	Phase4->SetcontrolMax(1);

	Phase4->SetcontrolMin(-1);
	Phase4->SetcontrolMax(1);
	//No parameter

	Phase4->SetpathMin(1);
	Phase4->SetpathMax(1);

	Phase4->SeteventMin(af);
	Phase4->SeteventMin(ef);
	Phase4->SeteventMin(incf);
	Phase4->SeteventMin(Omf);
	Phase4->SeteventMin(omf);

	Phase4->SeteventMax(af);
	Phase4->SeteventMax(ef);
	Phase4->SeteventMax(incf);
	Phase4->SeteventMax(Omf);
	Phase4->SeteventMax(omf);

	Phase4->SetTimeGuess(t3);
	Phase4->SetTimeGuess(t4);

	Phase4->SetStateGuess(1,rout(0));
	Phase4->SetStateGuess(1,rout(0));

	Phase4->SetStateGuess(2,rout(1));
	Phase4->SetStateGuess(2,rout(1));

	Phase4->SetStateGuess(3,rout(2));
	Phase4->SetStateGuess(3,rout(2));

	Phase4->SetStateGuess(4,vout(0));
	Phase4->SetStateGuess(4,vout(0));

	Phase4->SetStateGuess(5,vout(1));
	Phase4->SetStateGuess(5,vout(1));

	Phase4->SetStateGuess(6,vout(2));
	Phase4->SetStateGuess(6,vout(2));

	Phase4->SetStateGuess(7,m40);
	Phase4->SetStateGuess(7,m4f);

	Phase4->SetControlGuess(1,0);
	Phase4->SetControlGuess(1,0);

	Phase4->SetControlGuess(2,1);
	Phase4->SetControlGuess(2,1);

	Phase4->SetControlGuess(3,0);
	Phase4->SetControlGuess(3,0);

	shared_ptr<Linkage> Link1 (new Linkage(1,1,2));
	Link1->SetLinkMin(0);
	Link1->SetLinkMin(0);
	Link1->SetLinkMin(0);

	Link1->SetLinkMin(0);
	Link1->SetLinkMin(0);
	Link1->SetLinkMin(0);
	Link1->SetLinkMin(-6*m_dry_srb);

	Link1->SetLinkMax(0);
	Link1->SetLinkMax(0);
	Link1->SetLinkMax(0);

	Link1->SetLinkMax(0);
	Link1->SetLinkMax(0);
	Link1->SetLinkMax(0);
	Link1->SetLinkMax(-6*m_dry_srb);

	shared_ptr<Linkage> Link2 (new Linkage(2,2,3));

	Link2->SetLinkMin(0);
	Link2->SetLinkMin(0);
	Link2->SetLinkMin(0);

	Link2->SetLinkMin(0);
	Link2->SetLinkMin(0);
	Link2->SetLinkMin(0);
	Link2->SetLinkMin(-3*m_dry_srb);

	Link2->SetLinkMax(0);
	Link2->SetLinkMax(0);
	Link2->SetLinkMax(0);

	Link2->SetLinkMax(0);
	Link2->SetLinkMax(0);
	Link2->SetLinkMax(0);
	Link2->SetLinkMax(-3*m_dry_srb);

	shared_ptr<Linkage> Link3(new Linkage(3,3,4));
	Link3->SetLinkMin(0);
	Link3->SetLinkMin(0);
	Link3->SetLinkMin(0);

	Link3->SetLinkMin(0);
	Link3->SetLinkMin(0);
	Link3->SetLinkMin(0);
	Link3->SetLinkMin(-m_dry_first);

	Link3->SetLinkMax(0);
	Link3->SetLinkMax(0);
	Link3->SetLinkMax(0);

	Link3->SetLinkMax(0);
	Link3->SetLinkMax(0);
	Link3->SetLinkMax(0);
	Link3->SetLinkMax(-m_dry_first);

	//Step3
	shared_ptr<FunctionWrapper> userfun(new LaunchFunction());
	shared_ptr<OptimalProblem> optpro(new OptimalProblem(4,3,userfun));
	optpro->AddPhase(Phase1);
	optpro->AddPhase(Phase2);
	optpro->AddPhase(Phase3);
	optpro->AddPhase(Phase4);

	optpro->AddLinkage(Link1);
	optpro->AddLinkage(Link2);
	optpro->AddLinkage(Link3);
	
		
		try
		{
			//Step 4
			app->SetOptimalControlProblem(optpro);
			app->SolveOptimalProblem();
		}
		catch (LpopcException & e)
		{
			std::cout << "got an error" << std::endl;
		}
		catch (...)
		{
			std::cout << "got an error" << std::endl;
		}

	
	return 0;
}

void Launchoe2rv( vec& oe,double mu, vec& ri, vec& vi )
{
	assert(oe.n_elem==6);
	double a=oe(0); 
	double e=oe(1);
	double i=oe(2); 
	double Om=oe(3); 
	double om=oe(4); 
	double nu=oe(5);
	double p = a*(1-e*e);
	double r = p/(1+e*cos(nu));

	vec rv(3);
	rv[0] = r*cos(nu);
	rv[1] = r*sin(nu);
	rv[2] = 0;
	vec vv(3);
	vv[0] = -sin(nu);
	vv[1] = e + cos(nu);
	vv[2] = 0;
	 vv*= sqrt(mu/p);

	double cO = cos(Om); 
	double sO = sin(Om);
	double co = cos(om); 
	double so = sin(om);
	double ci = cos(i);  
	double si = sin(i);

	mat m_R;
	m_R << cO*co - sO*so*ci << -cO*so - sO*co*ci << sO*si << endr
		<< sO*co + cO*so*ci << -sO*so + cO*co*ci << -cO*si << endr
		<< so*si << co*si << ci << endr;
	/*m_R(3, 3, cO*co - sO*so*ci, sO*co + cO*so*ci, so*si,
				-cO*so-sO*co*ci,-sO*so+cO*co*ci,co*si,
				sO*si,                  -cO*si,                  ci);*/
	ri = m_R*rv;
	vi = m_R*vv;
}

void Launchrv2oe( vec& oe,vec rv,vec vv,double mu )
{
	vec K(3);
	K[0] = 0.0; K[1] = 0.0; K[2] = 1.0;
	vec hv=cross(rv,vv);
	vec nv=cross(K,hv);
	double n=sqrt(dot(nv,nv));
	double h2=dot(hv,hv);
	double v2=dot(vv,vv);
	double r=sqrt(dot(rv,rv));

	vec ev=rv*(v2-mu/r)-vv*dot(rv,vv);
	ev*=(1.0/mu);
	double p=h2/mu;
	 double e=sqrt(dot(ev,ev));//eccentricity
	 double a=p/(1-e*e);//semimajor axis
	 double i  = acos(hv(2)/sqrt(h2));//inclination
	 double Om1=acos(nv(0)/n);//RAAN
	 double eps=datum::eps;
	 if (nv(1)<0-eps)
	 {
		 Om1=2*PI-Om1;
	 }
	 double Om2=acos(dot(nv,ev)/n/e);

	 if (ev(2)<0)
	 {
		 Om2=2*PI-Om2;
	 }

	 double nu=acos(dot(ev,rv)/e/r);
	 if (dot(rv,vv)<0)
	 {
		 nu=2*PI-nu;
	 }
	 oe(0)=a;
	 oe(1)=e;
	 oe(2)=i;
	 oe(3)=Om1;
	 oe(4)=Om2;
	 oe(5)=nu;
}



void LaunchFunction::MayerCost( SolCost&mySolcost,double& mayer )
{
	mat xf=mySolcost.terminal_state_;
	if (mySolcost.phase_num_==4)
	{
		mayer=-xf(6);
	}else
	{
		mayer=0.0;
	}
}

void LaunchFunction::DerivMayer( SolCost&mySolcost,rowvec& deriv_mayer )
{
}

void LaunchFunction::LagrangeCost( SolCost& mySolcost,vec& langrange )
{
	vec t=mySolcost.time_;
	langrange=zeros(t.n_elem,1);
}

void LaunchFunction::DerivLagrange( SolCost& mySolcost,mat& deriv_langrange )
{
}

void LaunchFunction::DaeFunction( SolDae& mySolDae,mat& stateout,mat& pathout )
{
	vec t=mySolDae.time_;
	mat x=mySolDae.state_;
	mat u=mySolDae.contol_;

	int iphase=mySolDae.phase_num_;
	mat r=x.cols(0,2);
	mat v=x.cols(3,5);
	vec m=x.col(6);
	
	mat rad = sqrt(sum(r % r, 1));
	mat omega_matrix(CONSTANTS.omega_matrix, 3, 3);
	mat omegacrossr=r*trans(omega_matrix);
	mat vrel=v-omegacrossr;
	mat speedrel=sqrt(sum(vrel%(vrel),1));
	mat altitude=rad-CONSTANTS.Re;
	mat ret=-altitude/CONSTANTS.H;
	mat rho=exp (ret)*CONSTANTS.rho0;
	mat bc = rho / (m * 2)*(CONSTANTS.sa*CONSTANTS.cd);
	mat bcspeed=bc%(speedrel);
	mat bcspeedmat=repmat(bcspeed,1,3);
	mat Drag=bcspeedmat*(-1.0);
	Drag=Drag%vrel;//% means elem *
	mat muoverradcubed=ones(rad.n_rows,rad.n_cols)*CONSTANTS.mu;
	muoverradcubed = muoverradcubed / (pow(rad, 3));
	mat muoverradcubedMat=repmat(muoverradcubed,1,3);
	mat grav =-muoverradcubedMat %(r);
	mat T_tot,mdot;
	if (iphase==1)
	{ 
		mat T_srb,T_first,m1dot,m2dot;
		T_srb = ones(t.n_elem,1)*(6*CONSTANTS.thrust_srb);
		T_first =ones(t.n_elem,1)* (CONSTANTS.thrust_first);
		T_tot = T_srb+T_first;
		m1dot=zeros(T_tot.n_rows,T_tot.n_cols);
		m2dot=m1dot;
		m1dot -=T_srb/(CONSTANTS.g0*CONSTANTS.ISP_srb);
		m2dot -=T_first/(CONSTANTS.g0*CONSTANTS.ISP_first);
		mdot = m1dot+m2dot;

	}else if (iphase==2)
	{
		mat T_srb,T_first,m1dot,m2dot;
		T_srb = ones(t.n_elem,1)*(3*CONSTANTS.thrust_srb);
		T_first =ones(t.n_elem,1)* (CONSTANTS.thrust_first);
		T_tot = T_srb+T_first;
		m1dot=zeros(T_tot.n_rows,T_tot.n_cols);
		m2dot=m1dot;
		m1dot -=T_srb/(CONSTANTS.g0*CONSTANTS.ISP_srb);
		m2dot -=T_first/(CONSTANTS.g0*CONSTANTS.ISP_first);
		mdot = m1dot+m2dot;
	}else if( iphase==3){
		mat T_first;
		T_first =ones(t.n_elem,1)* CONSTANTS.thrust_first;
		T_tot = T_first;
		mdot=zeros(T_tot.n_rows,T_tot.n_cols);
		mdot -=T_first/(CONSTANTS.g0*CONSTANTS.ISP_first);
	}else if (iphase==4)
	{	vec T_second;
		T_second = ones(t.n_elem,1)*CONSTANTS.thrust_second;
		T_tot=T_second;
		mdot=zeros(T_tot.n_rows,T_tot.n_cols);
		mdot -=T_second/(CONSTANTS.g0*CONSTANTS.ISP_second);
	}

	pathout = sum(u%(u),1);
	mat Toverm = T_tot/(m);
	mat Tovermmat = repmat(Toverm,1,3);

	mat thrust = Tovermmat%(u);

	mat rdot = v;
	mat vdot = thrust+Drag+grav;
	stateout = zeros(t.n_elem, x.n_cols);// allocate memroy!!!!!
	stateout.cols(0,rdot.n_cols-1)=rdot;
	stateout.cols(rdot.n_cols,rdot.n_cols+vdot.n_cols-1)=vdot;
	stateout.col(rdot.n_cols+vdot.n_cols)=mdot;
	}

void LaunchFunction::DerivDae( SolDae& mySolDae,mat& deriv_state,mat& deriv_path )
{
}

void LaunchFunction::EventFunction( SolEvent& mySolEvent,vec& eventout )
{
	vec xf=mySolEvent.terminal_state_;
	int iphase=mySolEvent.phase_num_;
	if (iphase==4)
	{
		vec oe=zeros(6,1);
		Launchrv2oe(oe, xf.subvec(0, 2), xf.subvec(3, 5), CONSTANTS.mu);
		eventout=oe(span(0,4));
	}
}

void LaunchFunction::DerivEvent( SolEvent& mySolEvent,mat& deriv_event )
{
}

void LaunchFunction::LinkFunction( SolLink& mySolLink,vec& linkageout )
{
	vec x0_right=mySolLink.right_state_;
	vec xf_left=mySolLink.left_state_;
	linkageout=x0_right-xf_left;
}

void LaunchFunction::DerivLink( SolLink& mySolLink,mat& derive_link )
{
}
