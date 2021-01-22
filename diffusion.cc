/**
 * @file diffusion.cc
 * @authors Luca Maccione, Daniele Gaggero, Pedro de la Torre
 * @email luca.maccione@desy.de
 * @email daniele.gaggero@sissa.it
 * @email pedro.delatorreluque@fysik.su.se
 */


#include "diffusion.h"

#include "geometry.h"
#include "grid.h"
#include "input.h"
#include "sources.h"
#include "errorcode.h"

#include "bfield.h"

using namespace std;

TDiffusionCoefficient::TDiffusionCoefficient(TGrid* Coord, Input* in, TSource* SourceTerm, TBField* bf, std::vector<int>& nuclei, bool El) {

	_fCoordinates = Coord; 

	zt = in->zt;
	nrn_sn = SourceTerm->GetSource(8.5,0.0,0.0);
	index_radial = in->index_radial;
	set_profile = in->set_profile;

	delta = in->delta;
	delta_H = in->delta_H;
	delta_L = in->delta_L;
	rho_H = in->rho_H;
	rho_L = in->rho_L;
	
	if (in->Diff_table == true){
	  if (in->VariableD == true or in->OneZone_interp){
	    DiffFile_nameHalo = in->Halo_table; 
	    DiffFile_nameEDisc= in->Disk_table;}
 
	  else  DiffFile_name = in->Disk_table;
	}
	
}

double TDiffusionCoefficient::GetProfile(double x, double y, double zeta, TSource* SourceTerm) {

	double radial = 1;
	double result = 1.;

	switch(set_profile) {
	case Constant :
		return 1.0;
		break;

	case Exp :
		return exp(fabs(zeta)/zt);
		break;
	case Qtau :
		radial = SourceTerm->GetSource(x,y,zeta);
		radial /= nrn_sn;
		return pow(radial, index_radial);
		break;
	case Test :
		return 1.;
		break;		
	default :
		return -1;
	}

}


double TDiffusionCoefficient::GetXDerivative(double x, double y, double z, TSource* SourceTerm) {

	double delta_x = 0.00002; //Coord->GetDeltaX();
	double xup = x + 0.5*delta_x;
	double xdo = x - 0.5*delta_x;
	return ( GetProfile(xup,y,z, SourceTerm) - GetProfile(xdo,y,z, SourceTerm) )/delta_x ;
}

double TDiffusionCoefficient::GetYDerivative(double x, double y, double z, TSource* SourceTerm) {

	double delta_y = 0.00002; //Coord->GetDeltaY();
	double yup = y + 0.5*delta_y;
	double ydo = y - 0.5*delta_y;
	return ( GetProfile(x,yup,z, SourceTerm) - GetProfile(x,ydo,z, SourceTerm) )/delta_y;
}

double TDiffusionCoefficient::GetZDerivative(double x, double y, double z, TSource* SourceTerm) {

	double delta_z = 0.00002; //Coord->GetDeltaZ();
	double zup = z + 0.5*delta_z;
	double zdo = z - 0.5*delta_z;
	return ( GetProfile(x,y,zup, SourceTerm) - GetProfile(x,y,zdo, SourceTerm) )/delta_z;
}

// 2D

TDiffusionCoefficient2D::TDiffusionCoefficient2D(TGrid* Coord, Input* in, TSource* SourceTerm, TBField* bf, std::vector<int>& nuclei, bool El) : TDiffusionCoefficient(Coord, in, SourceTerm, bf, nuclei, El){
  
							   
  delta = in->delta;
  
  vector<double> pp;
  vector<double> beta;
  vector<double> Ek;
  if (El) {
    pp = Coord->GetMomentumEl();
    
    beta = Coord->GetBetaEl();
    Ek = Coord->GetEk();
  }
  else {
    pp = Coord->GetMomentum();
    beta = Coord->GetBeta();
    Ek = Coord->GetEk();
  }
  
  DRAGONmomentum = pp;
  
  if (in->VariableD == false) { 
    
    if (in->Diff_table == false) {
      
      delta = in->delta;
      
      if (in->BreakDiff == false) cout << "building diff. coeff. with constant delta" << endl;
      else { cout << "building diff. coeff. with break in delta at ";
	if (in->LowE_b == true) cout << in->rho_L << endl;
	else if (in->Double_b == true) cout << in->rho_L << " and at " << in->rho_H << endl;
	else cout << in->rho_H << endl;      }
      
      for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){ 
	
	int Z = -1000;  // nucleus charge
	int A = -1000;  // nucleus mass
	Utility::id_nuc(nuclei[nuloop], A, Z);   
	
	if (A==0){ 
	  A=1;
	  Z=1;
	}
	
	std::vector<double> thisD (pp.size(), 0.);
	
	for (int i = 0; i < pp.size(); ++i) {
	  
	  if ( delta > 0.0 ){
	    	      
	    if (in->BreakDiff == true){

	      if (in->LowE_b == true){  // Low energy break
		s_L = in->s_L; 
		Deltadelta_L = in->delta_L - delta;
		R_L = in->rho_L; 
		if(in->D_ref_rig>R_L){
		  //smoothly broken PL!
		  long double K = 1.e22;
		  thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(R_L/in->D_ref_rig, delta)*pow(pp[i]*A/(fabs(Z)*R_L), in->delta_L) *(1-((1 + tanh(K*((pp[i]*A/fabs(Z))-R_L)))/2.)) + ((1 + tanh(K*((pp[i]*A/fabs(Z))-R_L)))/2.) * in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta);	     
		  //Sharply_broken:
		  //if(pp[i]<=rho_H*fabs(Z)/A) thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(R_L/in->D_ref_rig, delta)*pow(pp[i]*A/(fabs(Z)*R_L), in->delta_L);
		  //if(pp[i]>rho_H*fabs(Z)/A)  thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta);  
		}
		else 
		  thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta) * pow(1 + pow(pp[i]*A/(fabs(Z)*R_L),Deltadelta_L/s_L),s_L);
	      }

	      else if (in->Double_b == true){ //Doubly broken diff coeff.
		s_H = in->s_H; //0.04;
		Deltadelta_H = in->delta_H - delta; // -0.14
		R_H = in->rho_H; // 312.
		s_L = in->s_L; //0.04;
		Deltadelta_L = in->delta_L - delta;
		R_L = in->rho_L;

		thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta) * pow(1 + pow(pp[i]*A/(fabs(Z)*R_L),Deltadelta_L/s_L),s_L) * pow(1 + pow(pp[i]*A/(fabs(Z)*R_L), (-1*Deltadelta_H)/s_H), -s_H);
		//following Weinrich et al, A&A 639, A131 (2020); https://doi.org/10.1051/0004-6361/202037875
		
	      }
	      
	      else {  // High energy break *** Smooth diffusion coefficient in rigidity, from Genolini et. al, PRL 119, 241101 (2017) *** //
		s_H = in->s_H; //0.04;
		Deltadelta_H = delta - in->delta_H; // 0.14
		R_H = in->rho_H; // 312.
		if(in->D_ref_rig>R_H){
		  //smoothly broken!
		  double K = 1.e22;
		  thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(R_H/in->D_ref_rig, delta)*pow(pp[i]*A/(fabs(Z)*R_H), delta)*(1-((1 + tanh(K*((pp[i]*A/fabs(Z))-R_H)))/2.)) + ((1 + tanh(K*((pp[i]*A/fabs(Z))-R_H)))/2.) * in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), in->delta_H);
		  //sharply broken PL
		  //if(pp[i]<=rho_H*fabs(Z)/A) thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(R_H/in->D_ref_rig, delta)*pow(pp[i]*A/(fabs(Z)*R_H), delta);
		  //if(pp[i]>rho_H*fabs(Z)/A)  thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta_H);
		}
		else 
		  thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta)/pow(1 + pow(pp[i]*A/(fabs(Z)*R_H),Deltadelta_H/s_H),s_H);  
	      } //End of High energy break (default)
	      
	    } //end of if BreakDiff

	    else // NO BREAK!	
	      thisD[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta);  
	    
	  }//End of delta>0
	  
	  else {
	    const double lgEk = log10(Ek[i]);
	    thisD[i] = in->D0*pow(10.0,-0.8990+1.1781*lgEk-0.1800*lgEk*lgEk+0.0141*lgEk*lgEk*lgEk);
	  }
	  
	} //End Energy (momentum) loop    
	sp.push_back(thisD);
	
	
      } // End of the nuclei loop
      
    } 
    
    else{ //if (Diff_table = true)
      

      if (in->OneZone_interp){ //Trick to parametrize a two-zone problem into a one-zone problem, following: Tomassetti, 2012  	arXiv:1204.4492
	
	std::cout << "Building one-zone TRICKY for two-zone diffusion... " << std::endl;

	datafileDiff_H.open(DiffFile_nameHalo.c_str());
	datafileDiff_ED.open(DiffFile_nameEDisc.c_str());
	
	if (!datafileDiff_ED){
	  std::cout << "Extended disc diffusion file not found!! --> Check the folder data in your setup for " << DiffFile_nameEDisc << "\n" << std::endl;
	  exit(-1);}
	
	if (!datafileDiff_H){
	  std::cout << "Halo diffusion file not found!! --> Check the folder data in your setup for " << DiffFile_nameHalo << "\n" << std::endl;
	  exit(-1);}
	
	std::cout << std::endl << "Reading diffusion coefficients: " << "\n    Extended Disc -> " << DiffFile_nameEDisc << "\n    Halo -> " << DiffFile_nameHalo << std::endl;
	
	double moment = 0.;
	double Dval, rlL; 	     	
	
	while (datafileDiff_H >> moment >> rlL >> Dval) {
	  
	  std::cout << "halo diffusion: " <<  moment << " " << Dval<< std::endl;
	  
	  DiffmomVector.push_back(moment); //   CHECK UNITS!!! The DragonEnergy is in GeV
	  Diff_vect_H.push_back(Dval/kpc/kpc*Myr); //*pow(10.,-28)); //Protons!
	  
	} // end of the while loop
	datafileDiff_H.close(); 
	
	moment = 0.; Dval = 0;	     
	
	while (datafileDiff_ED >> moment >> rlL >> Dval) {
	  
	  std::cout << moment << " " << Dval<< std::endl;
	  Diff_vect_ED.push_back(Dval/kpc/kpc*Myr); //*pow(10.,-28)); //Protons! 
	  
	} // end of the while loop
	datafileDiff_ED.close(); 
	
	std::cout << "...diffusion coefficients successfully read from files" << std::endl << std::endl;
	


	double ipP, dist0 = 100000.2;
	for (int ipd = 0; ipd < DiffmomVector.size(); ipd++) {
	  double dist = fabs(Diff_vect_ED[ipd] - Diff_vect_H[ipd]); //Find the intersection point... roughly
	  if (dist < dist0){
	    dist0 = dist;
	    ipP = ipd;}
	}
	
	double intersec_D = Diff_vect_H[ipP+1];
	double intersec_p = DiffmomVector[ipP+1];
	std::cout << "Intersection momentum: " << intersec_p << ", pseudo_D0: " << intersec_D*kpc*kpc/Myr << std::endl;
	double delta_H, delta_D, local_slope, epsi = 0.15, Xi = 0.36;

	Diff_vect.push_back(intersec_D);
	for (int ipd = 1; ipd < DiffmomVector.size(); ipd++) { //calculation of the effective slope
	    delta_H = log10(Diff_vect_H[ipd]/Diff_vect_H[ipd-1])/log10(DiffmomVector[ipd]/DiffmomVector[ipd-1]);
	    delta_D = log10(Diff_vect_ED[ipd]/Diff_vect_ED[ipd-1])/log10(DiffmomVector[ipd]/DiffmomVector[ipd-1]);
	    std::cout << "\nMomentum:" << DiffmomVector[ipd] << ", slope halo: " << delta_H << ", slope disk: " << delta_D << std::endl;
	    local_slope = delta_D + ( (delta_H - delta_D)/( 1 + (epsi/(Xi*(1-epsi)))* pow(DiffmomVector[ipd]/intersec_p, delta_H-delta_D) ));	  
	    Diff_vect.push_back(Diff_vect.back() * pow(DiffmomVector[ipd]/DiffmomVector[ipd-1], local_slope) );
	    
	}

	
	std::vector<double> thisDsp (pp.size(), 0.);
	double Dffinter;
	for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
	  int Z = -1000;  // nucleus charge
	  int A = -1000;  // nucleus mass
	  Utility::id_nuc(nuclei[nuloop], A, Z); 
	  
	  double AZ = (double)A/(double)fabs(Z);
	  
	  for (int ie = 0; ie < pp.size(); ie++) {
	   
	    Dffinter = GetDInterp(pp[ie]*AZ, AZ, "One_zone");     //D(p,A,Z) = D(p',A=1,Z=1) --> p' = p *A/fabs(Z) p of the nuclei is p'
	  
	    if (Dffinter > 100.)  
	      std::cout << "\n\nSomething strange happens";

	    std::cout << "Diff for nuclei " << nuclei[nuloop] <<  " is " << Dffinter*kpc*kpc/Myr << " for momentum: " << pp[ie] << std::endl;
	    thisDsp[ie] = Dffinter;
	  }	  
	  sp.push_back(thisDsp);
	}    

       	
      } //end OneZone_interp
      

      else{//Normal one zone model with Diff. read from table
	datafileDiff.open(DiffFile_name.c_str());
	if (!datafileDiff){
	  std::cout << "Diffusion file not found!! --> Check the folder data in your setup for the file " << DiffFile_name << "\n" << std::endl;
	  exit(-1);}
	
	std::cout << std::endl << "Reading diffusion coefficients from " << DiffFile_name << " ... " << std::endl;
	
	double moment = 0.;
	double Dval, rlL; 	     
	
	DiffmomVector.clear();	
	
	while (datafileDiff >> moment >> rlL >> Dval) {
	  
	  std::cout << moment << " " << Dval << std::endl;
	  
	  DiffmomVector.push_back(moment); //   CHECK UNITS!!! The DragonEnergy is in GeV
	  Diff_vect.push_back(Dval/kpc/kpc*Myr); //Protons!
	  
	} // end of the while loop
	datafileDiff.close(); 
	
	std::cout << "...diffusion coefficients successfully read from file" << std::endl << std::endl;
	
	for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
	  
	  int Z = -1000;  // nucleus charge
	  int A = -1000;  // nucleus mass
	  Utility::id_nuc(nuclei[nuloop], A, Z); 
	  
	  std::vector<double> thissp (pp.size(), 0.);
	  
	  for (int ie = 0; ie < pp.size(); ie++) {
	    double AZ = (double)A/(double)fabs(Z);
	    //std::cout << A << "  " << Z << "  " << AZ << std::endl;
	    
	    double Dffinter = GetDInterp(pp[ie]*AZ, AZ, "One_zone");    
	    
	    
	    if (Dffinter > 100.)  
	      std::cout << "\n\nSomething strange happens";
	    
	    thissp[ie] = Dffinter;
	  }
	  
	  sp.push_back(thissp);
	}
      }
    } //end of Diff_table
  }
  
  else { //Beginning of variable Diffusion coefficient
    
    if (in->TwoZoneModel == false) { 
      
      double delta_A = in->delta_A;
      double delta_B = in->delta_B;
      
      cout << "building diff. coeff. with variable delta(R) = Ar + B " << endl;
      cout << delta_A << "r + " << delta_B << endl;
      
      vector<double> r_vec = Coord->GetR();
      for (int i = 0; i < r_vec.size(); i++)
	cout << r_vec[i] << " " ;
      cout << endl;
      vector<double> z_vec = Coord->GetZ();
      for (int i = 0; i < z_vec.size(); i++)
	cout << z_vec[i] << " " ;
      cout << endl;
      dimr = r_vec.size();
      cout << dimr << endl;
      dimz = z_vec.size();
      cout << dimz << endl;
      unsigned int irsun = (unsigned int) ((in->robs-r_vec.front())/(r_vec.back()-r_vec.front())*(double)(dimr-1));
      unsigned int izsun = (unsigned int) ((in->zobs-z_vec.front())/(z_vec.back()-z_vec.front())*(double)(dimz-1));
      
      
      for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
	int Z = -1000;  // nucleus charge
	int A = -1000;  // nucleus mass
	Utility::id_nuc(nuclei[nuloop], A, Z);    
	std::vector<double> mysp (pp.size(), 0.);
	std::vector<double> myspe (pp.size()*z_vec.size()*r_vec.size(), 0.);
	int ipzr = 0;
	for (int ip = 0; ip < pp.size(); ++ip) {
	  //cout << ip << endl;
	  for (int iz = 0; iz < z_vec.size(); iz++) {
	    for (int ir = 0; ir < r_vec.size(); ir++) {
	      
	      double current_delta = delta_A*r_vec[ir]+delta_B;
	      
	      if (r_vec[ir] > in->diff_threshold)
		current_delta = delta_A*in->diff_threshold + delta_B;
	      
	      double current_delta_hi = current_delta;
	      
	      if(in->D_ref_rig<=rho_H){
		if(pp[ip]<=rho_H*fabs(Z)/A) myspe[ipzr] = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]*A/(fabs(Z)*in->D_ref_rig), current_delta);
		if(pp[ip]>rho_H*fabs(Z)/A)  myspe[ipzr] = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]*A/(fabs(Z)*rho_H), current_delta_hi)*pow(rho_H/in->D_ref_rig, current_delta);
	      }
	      
	      if(in->D_ref_rig>rho_H){
		if(pp[ip]<=rho_H*fabs(Z)/A) myspe[ipzr] = in->D0*pow(beta[ip], in->etaT)*pow(rho_H/in->D_ref_rig, current_delta_hi)*pow(pp[ip]*A/(fabs(Z)*rho_H), current_delta);
		if(pp[ip]>rho_H*fabs(Z)/A)  myspe[ipzr] = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]*A/(fabs(Z)*in->D_ref_rig), current_delta_hi);
	      }	
	      
	      
	      if (ir == irsun && iz == izsun) {
		if(in->D_ref_rig<=rho_H){
		  if(pp[ip]<=rho_H*fabs(Z)/A) mysp[ip] = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]*A/(fabs(Z)*in->D_ref_rig), current_delta);
		  if(pp[ip]>rho_H*fabs(Z)/A)  mysp[ip] = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]*A/(fabs(Z)*rho_H), current_delta_hi)*pow(rho_H/in->D_ref_rig, current_delta);
		}
		if(in->D_ref_rig>rho_H){
		  if(pp[ip]<=rho_H*fabs(Z)/A) mysp[ip] = in->D0*pow(beta[ip], in->etaT)*pow(rho_H/in->D_ref_rig, current_delta_hi)*pow(pp[ip]*A/(fabs(Z)*rho_H), current_delta);
		  if(pp[ip]>rho_H*fabs(Z)/A)  mysp[ip] = in->D0*pow(beta[ip], in->etaT)*pow(pp[ip]*A/(fabs(Z)*in->D_ref_rig), current_delta_hi);
		}
	      }
	      ipzr += 1;
	   }
	  }
	  
	} //End of energy loop
	
	sp.push_back(mysp);
	spectrum_extended.push_back(myspe);
      } //End of nuclei loop
    }
    
    else if (in->TwoZoneModel == true)  { 
      
      if (in->Diff_table == false){
	
	double slope_disk = 0.15; //in->delta;
	double slope_halo = 0.67; //slope_disk + in->Delta_delta;
	double L_connecting = 0.05; //in->L_connecting;
	double disk_size = 1.701; //in->disk_size;
	
	cout << "Building two-zone model..." << endl;
	cout << "delta_disk = " << slope_disk << "; delta_halo = " << slope_halo << endl;
	
	vector<double> r_vec = Coord->GetR();
	vector<double> z_vec = Coord->GetZ();
	//cout << dimr << endl;
	dimz = z_vec.size();
	//cout << dimz << endl;
	unsigned int irsun = (unsigned int) ((in->robs-r_vec.front())/(r_vec.back()-r_vec.front())*(double)(dimr-1));
	unsigned int izsun = (unsigned int) ((in->zobs-z_vec.front())/(z_vec.back()-z_vec.front())*(double)(dimz-1));
	for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
	  int Z = -1000;  // nucleus charge
	  int A = -1000;  // nucleus mass
	  Utility::id_nuc(nuclei[nuloop], A, Z); 
	  
	  std::vector<double> thissp (pp.size(), 0.);
	  std::vector<double> thisspe (pp.size()*z_vec.size()*r_vec.size(), 0.);
	  
	  int ipzr = 0;
	  
	  for (int ip = 0; ip < pp.size(); ++ip) {
	    for (int iz = 0; iz < z_vec.size(); iz++) {
	      for (int ir = 0; ir < r_vec.size(); ir++) {
		
		/*
		  double slope = 0.5 * slope_disk * (1.0 - tanh( (fabs(z_vec[iz]) - 0.5*disk_size)/L_connecting  ) ) + 0.5 * slope_halo * (1.0 + tanh( (fabs(z_vec[iz]) - 0.5*disk_size)/L_connecting  ) );
		  
		  //if (fabs(z_vec[iz]) <= disk_size) slope = slope_disk;
		  //else slope = slope_halo;
		  
		  
		  thisspe[ipzr] = in->D0 * pow(beta[ip], in->etaT) * pow(pp[ip]*A/(fabs(Z)*in->D_ref_rig), slope) ;	
		  
		  if (ir == irsun && iz == izsun) 
		  thissp[ip] = in->D0 * pow(beta[ip], in->etaT) * pow(pp[ip]*A/(fabs(Z)*in->D_ref_rig), slope) ;
		*/
				
		
		// Parametrisation from Feng, Tomassetti, Oliva,  Phys.Rev.D 94 (2016) 12, 123007  DOI:"10.1103/PhysRevD.94.123007"
		double _disk = (1.92/kpc/kpc*Myr*1.e28) * pow(beta[ip], -0.4) * pow(pp[ip]*A/(fabs(Z)*0.25), slope_disk) ;	
		
		double _halo = 0.36*(1.92/kpc/kpc*Myr*1.e28) * pow(beta[ip], -0.4) * pow(pp[ip]*A/(fabs(Z)*0.25), slope_halo) ;
		
		double connected = 0.5 * _disk * (1.0 - tanh( (fabs(z_vec[iz]) - 0.5*disk_size)/L_connecting  ) ) + 0.5*_halo * (1.0 + tanh( (fabs(z_vec[iz]) - 0.5*disk_size)/L_connecting  ) );
		
		
		thisspe[ipzr] = connected;
		if (ir == irsun && iz == 0) //izsun) 
		  thissp[ip] = connected;
		
		
		ipzr += 1;
		
	      }
	    }
	    
	  }
	  
	  sp.push_back(thissp);
	  spectrum_extended.push_back(thisspe);
	} //end of nuclei loop
	
      }

      
      else{ //if (Diff_table = true)
	
	vector<double> r_vec = Coord->GetR();
	vector<double> z_vec = Coord->GetZ();
	unsigned int irsun = (unsigned int) ((in->robs-r_vec.front())/(r_vec.back()-r_vec.front())*(double)(r_vec.size()-1));
	unsigned int izsun = (unsigned int) ((in->zobs-z_vec.front())/(z_vec.back()-z_vec.front())*(double)(z_vec.size()-1));
	
	double L_connecting = 0.050; //kpc   //in->L_connecting;
	double disk_size = 1.;   //kpc   //in->disk_size;
	
	datafileDiff_H.open(DiffFile_nameHalo.c_str());
	datafileDiff_ED.open(DiffFile_nameEDisc.c_str());
	
	if (!datafileDiff_ED){
	  std::cout << "Extended disc diffusion file not found!! --> Check the folder data in your setup for " << DiffFile_nameEDisc << "\n" << std::endl;
	  exit(-1);}
	
	if (!datafileDiff_H){
	  std::cout << "Halo diffusion file not found!! --> Check the folder data in your setup for " << DiffFile_nameHalo << "\n" << std::endl;
	  exit(-1);}
	
	std::cout << std::endl << "Reading diffusion coefficients: " << "\n    Extended Disc -> " << DiffFile_nameEDisc << "\n    Halo -> " << DiffFile_nameHalo << std::endl;
	
	
	double moment = 0.;
	double Dval, rlL; 	     	
	
	while (datafileDiff_H >> moment >> rlL >> Dval) {
	  
	  std::cout << "halo diffusion: " <<  moment << " " << Dval<< std::endl;
	  
	  DiffmomVector_H.push_back(moment); //   CHECK UNITS!!! The DragonEnergy is in GeV
	  Diff_vect_H.push_back(Dval/kpc/kpc*Myr); //*pow(10.,-28)); //Protons!
	  
	} // end of the while loop
	datafileDiff_H.close(); 
	
	
	moment = 0.; Dval = 0;	     
	
	while (datafileDiff_ED >> moment >> rlL >> Dval) {
	  
	  std::cout << moment << " " << Dval<< std::endl;
	  
	  DiffmomVector_ED.push_back(moment); //   CHECK UNITS!!! The DragonEnergy is in GeV
	  Diff_vect_ED.push_back(Dval/kpc/kpc*Myr); //*pow(10.,-28)); //Protons! 
	  
	} // end of the while loop
	datafileDiff_ED.close(); 
	
	std::cout << "...diffusion coefficients successfully read from files" << std::endl << std::endl;
	
	
	for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
	  int Z = -1000;  // nucleus charge
	  int A = -1000;  // nucleus mass
	  Utility::id_nuc(nuclei[nuloop], A, Z); 
	  //std::cout << A << "  " << Z << "  " << (float)fabs(Z)/(float)A << std::endl;
	  std::vector<double> thissp (pp.size(), 0.);
	  std::vector<double> thisspe (pp.size()*z_vec.size()*r_vec.size(), 0.);
	  int iezr = 0;
	  double AZ = (double)A/(double)fabs(Z);
	  //std::cout << A << "  " << Z << "  " << AZ << std::endl;
	  
	  for (int ie = 0; ie < pp.size(); ie++) {
	    double Dffinter_H = GetDInterp(pp[ie]*AZ, AZ, "Halo");     //D(p,A,Z) = D(p',A=1,Z=1) --> p = p' *A/fabs(Z) p of the nuclei is p'
	    //std::cout << "Diff in halo for nuclei " << nuclei[nuloop] <<  " is " << Dffinter_H << std::endl;
	    if (Dffinter_H > 100. || Dffinter_H < 0.)  
	      std::cout << "\n\nSomething strange happens for H";
	    
	    double Dffinter_ED = GetDInterp(pp[ie]*AZ, AZ, "EDisc");     //D(p,A,Z) = D(p',A=1,Z=1) --> p = p' *A/fabs(Z) p of the nuclei is p'
	    //std::cout << "Diff in extended disc for nuclei " << nuclei[nuloop] <<  " is " << Dffinter_ED << std::endl;
	    if (Dffinter_ED > 100. || Dffinter_ED < 0.)  
	      std::cout << "\n\nSomething strange happens for ED";
	    
	    for (int iz = 0; iz < z_vec.size(); iz++) {
	      for (int ir = 0; ir < r_vec.size(); ir++) {
		
		double diff_connected = 0.5 * in->normDisk * Dffinter_ED * (1.0 - tanh( (fabs(z_vec[iz]) - 0.5*disk_size)/L_connecting  ) ) + 0.5 * in->normHalo * Dffinter_H * (1.0 + tanh( (fabs(z_vec[iz]) - 0.5*disk_size)/L_connecting  ) );  
		
		
		if (diff_connected > 1000. | diff_connected < 0.) 
		  std::cout << "Something strange happens with diff_connected" << diff_connected  << std::endl;
		thisspe[iezr] = diff_connected;
		
		if (ir == irsun && iz == izsun){ 
		  thissp[ie] = diff_connected;		  
		}
   
		iezr +=1;
	      }
	    }
	  }
	  sp.push_back(thissp);
	  spectrum_extended.push_back(thisspe);
	}
	
      } //end of Diff_table
      
    }//End two zone
    
  } //End of else variable diffusion
  
  
  vector<double> r = Coord->GetR();
  vector<double> z = Coord->GetZ();
  
  dimr = r.size();
  dimz = z.size();
  
  
  //MW130705
  vector<double> dr_up = Coord->GetDeltaR_up();
  vector<double> dz_up = Coord->GetDeltaZ_up();
  vector<double> dr_down = Coord->GetDeltaR_down();
  vector<double> dz_down = Coord->GetDeltaZ_down();
  
  for (unsigned int i = 0; i < dimr; ++i) {
    double radius = r[i];
    double Deltar = Coord->GetDeltaR(i);
    
    for (unsigned int j = 0; j < dimz; ++j) {
      double zeta = z[j];
      double Deltaz = Coord->GetDeltaZ(j);
      
      dperp.push_back(GetProfile(radius, 0, zeta, SourceTerm));
      phi.push_back(0.5/Deltar*(dperp.back()/max(u,radius)+GetRDerivative(radius, 0, zeta, SourceTerm)));

    }
  }
  
  
  double sp_     = 0.;
  double D       = 0.;
  double D_rup   = 0.;
  double D_rdown = 0.;
  double D_zup   = 0.;
  double D_zdown = 0.;
  
  
  double sp_rup_   = 0.;
  double sp_rdown_ = 0;
  double sp_zup_   = 0.;
  double sp_zdown_ = 0;
  
  for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
    
    vector <double> r1 (pp.size()*dimr*dimz, 0.);
    vector <double> r2 (pp.size()*dimr*dimz, 0.);
    vector <double> r3 (pp.size()*dimr*dimz, 0.);
    
    vector <double> z1 (pp.size()*dimr*dimz, 0.);
    vector <double> z2 (pp.size()*dimr*dimz, 0.);
    vector <double> z3 (pp.size()*dimr*dimz, 0.);
    
    int ipzr = 0;
    for (unsigned int i = 0; i < dimr; ++i) {
      double dr_central = 0.5 * (dr_up[i]+dr_down[i]);
      
      for (unsigned int k = 0; k < dimz; ++k) {
	double dz_central = 0.5 * (dz_up[k]+dz_down[k]);
	
	//MW130624: Trying to gain speed at the cost of memory...
	//Tests gave me a speed up of 45% with no considerable memory loss!
	double indspat       = index(i,k);
	double indspat_rup   = index(i+1,k);
	double indspat_rdown = index(i-1,k);
	double indspat_zup   = index(i,k+1);
	double indspat_zdown = index(i,k-1);
	
	for  (int ip = 0; ip < pp.size(); ip++) {  
	  
	  if (in->VariableD == false)    {                  
	    sp_ = sp[nuloop][ip];
	    
	    D       = dperp[indspat]*sp_;
	    D_rup   = dperp[indspat_rup]*sp_;
	    D_rdown = dperp[indspat_rdown]*sp_;
	    D_zup   = dperp[indspat_zup]*sp_;
	    D_zdown = dperp[indspat_zdown]*sp_;
	    
	  }
	  
	  else if (in->VariableD == true)  {
	    
	    if (i==0 && k==0 && ip== 0)
	      cout << "Building variable delta CN coefficients... " << endl;
	    
	    sp_ = spectrum_extended[nuloop][index_pzr(ip,k,i)];
	    
	    if (i == dimr-1)
	      sp_rup_ = sp_;
	    else
	      sp_rup_ = spectrum_extended[nuloop][index_pzr(ip,k,i+1)];
	    
	    if (i == 0)
	      sp_rdown_ = sp_;
	    else
	      sp_rdown_ = spectrum_extended[nuloop][index_pzr(ip,k,i-1)];
	    
	    if (k == dimz-1)
	      sp_zup_ = sp_;
	    else
	      sp_zup_ = spectrum_extended[nuloop][index_pzr(ip,k+1,i)];
	    
	    if (k == 0)
	      sp_zdown_ = sp_;
	    else
	      sp_zdown_ = spectrum_extended[nuloop][index_pzr(ip,k-1,i)];
	    
	    D       = dperp[indspat]*sp_;
	    D_rup   = dperp[indspat_rup]*sp_rup_;
	    D_rdown = dperp[indspat_rdown]*sp_rdown_;
	    D_zup   = dperp[indspat_zup]*sp_zup_;
	    D_zdown = dperp[indspat_zdown]*sp_zdown_; 

	    
	  } 
	  
	  
	  r1[ipzr] = D/(dr_central*dr_down[i]) - (D_rup-D_rdown)/(4*dr_central*dr_central);
	  r2[ipzr] = D/(dr_central*dr_up[i]) + D/(dr_down[i]*dr_central);
	  r3[ipzr] = (D_rup-D_rdown)/(4*dr_central*dr_central) + D/(dr_up[i]*dr_central);
	  
	  if (in->VariableD == true && i > 0){	    
	    r1[ipzr] -= 0.5/dr_central * (D/max(u,r[i])); //this D/r term is due to the choice of the coordinates
	    r3[ipzr] += 0.5/dr_central * (D/max(u,r[i]));
	  }
	  
	  z1[ipzr] = D/(dz_central*dz_down[k]) - (D_zup-D_zdown)/(4*dz_central*dz_central);
	  z2[ipzr] = D/(dz_central*dz_up[k]) + D/(dz_down[k]*dz_central);
	  z3[ipzr] = (D_zup-D_zdown)/(4*dz_central*dz_central) + D/(dz_up[k]*dz_central);
  
	  
	  ipzr += 1;
	} //End energy loop
      } 
    }
    CNdiff_alpha1_r.push_back(r1);
    CNdiff_alpha2_r.push_back(r2);
    CNdiff_alpha3_r.push_back(r3);
    CNdiff_alpha1_z.push_back(z1);
    CNdiff_alpha2_z.push_back(z2);
    CNdiff_alpha3_z.push_back(z3); 
  }//End of nuloop
  
}



double TDiffusionCoefficient::GetDInterp(double men, double AZ, string zone) {
  
  std::vector<double> momvec, diffvec;
  
  if (zone == "Halo"){       momvec = DiffmomVector_H;    diffvec = Diff_vect_H;}
  else if (zone == "EDisc"){ momvec = DiffmomVector_ED;   diffvec = Diff_vect_ED;}
  else{ momvec = DiffmomVector;  diffvec = Diff_vect;}
  

  if (men/AZ - 0.01 > DRAGONmomentum.back() || men < DRAGONmomentum.front()) {
    std::cout << "Energy value over or below the energy of particles in DRAGON " << AZ << " " << men/AZ << " " << DRAGONmomentum.back() << std::endl;
    return 0.;
  }
  
  else if (men > momvec.back() || men < momvec.front()) { //Extrapolation   
    int je = 0, jee;
    double t;
    
    if (men > momvec.back()){ 
      je = momvec.size()-1; 	  jee = je;
      t = (men - momvec[jee])/(momvec[je]-momvec[je-1]);
    }
    
    else { // if (men < momvec.front()){
      je = 1; jee = 0;
      t = (men - momvec[jee])/(momvec[je]-momvec[je-1]);
      
    }
    
    double yH  = diffvec[jee] + t*(diffvec[je] - diffvec[je-1]);
    
    
    return max(yH, 0.);
    
  } //End of extrapolation block 
  
  int je = 0, jee;
  
  while (men > momvec[je])
    je++;
  
  if ((momvec[je] - men) > (men - momvec[je-1])) jee = je-1; 
  else
    jee = je;
  
  
  double Diffarr[momvec.size()], themom[momvec.size()];
  for(int dan=0; dan < momvec.size(); dan++){
    Diffarr[dan] = diffvec[dan];
    themom[dan] = momvec[dan]; 
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, momvec.size());
  
  gsl_spline_init (spline, themom, Diffarr, momvec.size()); //x axis must be always increasing!!
  double yH = gsl_spline_eval (spline, men, acc);
  
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
    
  if (diffvec[jee] == 0. && diffvec[jee-1] == 0.)
    yH = 0.;
  
  
  if (yH < 0.) std::cout << "bad gsl spline H" << yH << std::endl;
  
  yH = std::max(yH, 0.);
  return yH;
}




// 3D

TDiffusionCoefficient3D::TDiffusionCoefficient3D(TGrid* Coord, Input* in, TSource* SourceTerm, TBField* bf, TGeometry* geom, double A, double Z, std::vector<int> nuclei, bool El) : TDiffusionCoefficient(Coord, in, SourceTerm, bf, nuclei, El) {

  diffT = in->DiffT;
  
  vector<double> pp;
  vector<double> beta;
  if (El) {
    pp = Coord->GetMomentumEl();
    beta = Coord->GetBetaEl();
  }
  else {
    pp = Coord->GetMomentum();
    beta = Coord->GetBeta();
  }
  
  vector<double> x = Coord->GetX();
  vector<double> y = Coord->GetY();
  vector<double> z = Coord->GetZ();
  
  dimx = x.size();
  dimy = y.size();
  dimz = z.size();
  dimE = pp.size();
  
  //MW130624
  vector<double> dx_up = Coord->GetDeltaX_up();
  vector<double> dy_up = Coord->GetDeltaY_up();
  vector<double> dz_up = Coord->GetDeltaZ_up();
  vector<double> dx_down = Coord->GetDeltaX_down();
  vector<double> dy_down = Coord->GetDeltaY_down();
  vector<double> dz_down = Coord->GetDeltaZ_down();
  
  if (delta < 0){std::cout << "\nThis delta value is not allowed in the 3D model" << std::endl; exit(-55);}
  double smooth_factor = 0.04;
  if (diffT == Isotropic) {

    if (in->BreakDiff == false) cout << "3D building diff. coeff. with constant delta" << endl;
    else cout << "3D Diff. coeff. with break in delta from " << delta << " to " << in->delta_H << " at " << in-> rho_H << endl;
    for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
      std::vector <double> thissp (pp.size(), 0.);
	
	int Z = -1000;  // nucleus charge
	int A = -1000;  // nucleus mass
	Utility::id_nuc(nuclei[nuloop], A, Z);   
	
	if (A==0){ 
	  A=1;
	  Z=1;
	}

      for ( int i = 0; i < pp.size(); ++i ){

	if(in->D_ref_rig<=rho_H){   
	  if (in->BreakDiff == true)  // *** Smooth diffusion coefficient in rigidity, from Genolini et. al, PRL 119, 241101 (2017) *** // 
	    thissp[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta)/pow(1 + pow(pp[i]*A/(fabs(Z)*in->rho_H), in->delta_H/smooth_factor),smooth_factor);  
	    
	  else
	    thissp[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta);      
	    
	}

	if(in->D_ref_rig>rho_H){
	  if(pp[i]<=rho_H) thissp[i] = in->D0*pow(beta[i], in->etaT)*pow(rho_H/in->D_ref_rig, delta_H)*pow(pp[i]*A/(fabs(Z)*rho_H), delta);
	  if(pp[i]>rho_H)  thissp[i] = in->D0*pow(beta[i], in->etaT)*pow(pp[i]*A/(fabs(Z)*in->D_ref_rig), delta_H);
	}
	
	
      } // End energy loop
      sp.push_back(thissp);
    }
    
    for (unsigned int i = 0; i < dimx; ++i) {
      for (unsigned int j = 0; j < dimy; ++j) {
	for (unsigned int k = 0; k < dimz; ++k) {
	  double dperp_profile = GetProfile(x[i], y[j], z[k], SourceTerm);
	  
	  //mw, 130326, 130415
	  double spiral_factor_dperp = max( min( pow(geom->GetPattern(i,j,k), in->SA_diff), in->SA_cut_diff), 1./in->SA_cut_diff );
	  dperp.push_back( dperp_profile * pow( in->LB_diff, Coord->IsInLocalBubble(x[i],y[j],z[k]) ) * spiral_factor_dperp);
	}
      }
    }
    
     
    
    for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
      
      vector <double> x1 (pp.size()*dimx*dimy*dimz, 0.);
      vector <double> x2 (pp.size()*dimx*dimy*dimz, 0.);
      vector <double> x3 (pp.size()*dimx*dimy*dimz, 0.);
      
      vector <double> y1 (pp.size()*dimx*dimy*dimz, 0.);
      vector <double> y2 (pp.size()*dimx*dimy*dimz, 0.);
      vector <double> y3 (pp.size()*dimx*dimy*dimz, 0.);
      
      vector <double> z1 (pp.size()*dimx*dimy*dimz, 0.);
      vector <double> z2 (pp.size()*dimx*dimy*dimz, 0.);
      vector <double> z3 (pp.size()*dimx*dimy*dimz, 0.);

      int ipxyz = 0;      
      for (unsigned int i = 0; i < dimx; ++i) {
	double dx_central = 0.5 * (dx_up[i]+dx_down[i]);
	
	for (unsigned int j = 0; j < dimy; ++j) {
	  double dy_central = 0.5 * (dy_up[j]+dy_down[j]);
	  
	  for (unsigned int k = 0; k < dimz; ++k) {
	    double dz_central = 0.5 * (dz_up[k]+dz_down[k]);
	    
	    //MW130624: Trying to gain speed at the cost of memory...
	    //Tests gave me a speed up of 45% with no considerable memory loss!
	    double indspat       = index(i,j,k);
	    double indspat_xup   = index(i+1,j,k);
	    double indspat_xdown = index(i-1,j,k);
	    double indspat_yup   = index(i,j+1,k);
	    double indspat_ydown = index(i,j-1,k);
	    double indspat_zup   = index(i,j,k+1);
	    double indspat_zdown = index(i,j,k-1);
	    
	    for  (int ip = 0; ip < dimE; ip++) {
	      
	      double sp_ = sp[nuloop][ip];	      
	      double D       = dperp[indspat]*sp_;
	      double D_xup   = dperp[indspat_xup]*sp_;
	      double D_xdown = dperp[indspat_xdown]*sp_;
	      double D_yup   = dperp[indspat_yup]*sp_;
	      double D_ydown = dperp[indspat_ydown]*sp_;
	      double D_zup   = dperp[indspat_zup]*sp_;
	      double D_zdown = dperp[indspat_zdown]*sp_;
	      
	      x1[ipxyz] = D/(dx_central*dx_down[i]) - (D_xup-D_xdown)/(4*dx_central*dx_central);
	      x2[ipxyz] = D/(dx_central*dx_up[i]) + D/(dx_down[i]*dx_central);
	      x3[ipxyz] = (D_xup-D_xdown)/(4*dx_central*dx_central) + D/(dx_up[i]*dx_central);
	      
	      y1[ipxyz] = D/(dy_central*dy_down[j]) - (D_yup-D_ydown)/(4*dy_central*dy_central);
	      y2[ipxyz] = D/(dy_central*dy_up[j]) + D/(dy_down[j]*dy_central);
	      y3[ipxyz] = (D_yup-D_ydown)/(4*dy_central*dy_central) + D/(dy_up[j]*dy_central);
	      
	      z1[ipxyz] = D/(dz_central*dz_down[k]) - (D_zup-D_zdown)/(4*dz_central*dz_central);
	      z2[ipxyz] = D/(dz_central*dz_up[k]) + D/(dz_down[k]*dz_central);
	      z3[ipxyz] = (D_zup-D_zdown)/(4*dz_central*dz_central) + D/(dz_up[k]*dz_central);

	      ipxyz += 1;

	    }// End of Energy loop
	  }
	}
      }
      CNdiff_alpha1_x.push_back(x1);
      CNdiff_alpha2_x.push_back(x2);
      CNdiff_alpha3_x.push_back(x3);
      
      CNdiff_alpha1_y.push_back(y1);
      CNdiff_alpha2_y.push_back(y2);
      CNdiff_alpha3_y.push_back(y3);
      
      CNdiff_alpha1_z.push_back(z1);
      CNdiff_alpha2_z.push_back(z2);
      CNdiff_alpha3_z.push_back(z3);
    } //End of nuclei loop
    
  }//End of isotropic
  
  else if (diffT == Anisotropic) {
    
    cout << "In TDiffusionCoefficient3D, Anisotropic" << endl;
    for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
      std::vector<double> thissp (dimE, 0.);
      for (unsigned int ip = 0; ip < dimE; ip++){
	//MW: this still contains no break. 130725
	thissp[ip] =  (El) ? in->Dpar*pow(beta[ip], in->etaTpar)*pow(pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPar) : in->Dpar*pow(beta[ip], in->etaTpar)*pow(A*pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPar) ;
      }
      sp.push_back(thissp);
    }
    
    
    for (int i = 0; i < dimx; ++i) {
      for (int k = 0; k < dimy; ++k) {
	//double radius = sqrt(x[i]*x[i]+y[k]*y[k]);
	
	for (int j = 0; j < dimz; ++j) {
	  //  double zeta = z[j];
	  
	  double profDPerp = GetProfileDPerp(x[i], y[k], z[j], SourceTerm);
	  double profDPar = GetProfileDPar(x[i], y[k], z[j], SourceTerm);
	  
	  Dpar_profile.push_back(profDPar);
	  
	  //cout << "calling GetBregVersors(int,int,int) " << endl;	
	  std::vector<double> versors = bf->GetBregVersors(i,k,j);
	  //cout << versors[0] << " " << versors[1] << " " << versors[2] << endl;	
	  
	  std::vector<double> GradProfDPerp = GetGradProfileDPerp(x[i], y[k], z[j], SourceTerm);
	  
	  std::vector<double> GradProfDPar = GetGradProfileDPar(x[i], y[k], z[j], SourceTerm);
	  
	  std::map<char, std::vector<double> > GradVersors = GetGradVersors(x[i], y[k], z[j], bf);
	  
	  for (int nuloop = 0; nuloop < nuclei.size(); ++nuloop){  
	    for (unsigned int ip = 0; ip < dimE; ip++) {
	      
	      double Dpar = sp[nuloop][ip];
	      
	      double Dparprof = Dpar*profDPar;
	      dpar_total.push_back(Dparprof);
	      
	      double DPerp = (El) ? in->Dperp*pow(beta[ip], in->etaTperp)*pow(pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPerp) : in->Dperp*pow(beta[ip], in->etaTperp)*pow(A*pp[ip]/fabs(Z)/in->D_ref_rig, in->DeltaPerp);
	      
	      double Dperpprof = DPerp*profDPerp;
	      dperp_total.push_back(Dperpprof);
	      
	      double axx = (Dparprof - Dperpprof)*pow(versors[0],2) + Dperpprof;
	      double ayy = (Dparprof - Dperpprof)*pow(versors[1],2) + Dperpprof;
	      double azz = (Dparprof - Dperpprof)*pow(versors[2],2) + Dperpprof;
	      	      
	      double axy = 2.0*(Dparprof - Dperpprof)*versors[0]*versors[1];
	      double ayz = 2.0*(Dparprof - Dperpprof)*versors[2]*versors[1];
	      double axz = 2.0*(Dparprof - Dperpprof)*versors[0]*versors[2];
	      
	      alpha_xx.push_back(axx);
	      alpha_yy.push_back(ayy);
	      alpha_zz.push_back(azz);
	      
	      alpha_xy.push_back(axy);
	      alpha_yz.push_back(ayz);
	      alpha_xz.push_back(axz);
	      
	      double ux =
		(Dpar*GradProfDPar[0] - DPerp*GradProfDPerp[0])*pow(versors[0],2)
		+ DPerp*GradProfDPerp[0]
		+ (Dparprof - Dperpprof)*2.0*versors[0]*GradVersors['x'][0]    // dx alpha_xx
		+ (Dpar*GradProfDPar[1] - DPerp*GradProfDPerp[1])*versors[0] * versors[1]
		+ (Dparprof - Dperpprof)*versors[0]*GradVersors['y'][1]
		+ (Dparprof - Dperpprof)*versors[1]*GradVersors['y'][0]   // dy alpha_xy
		+ (Dpar*GradProfDPar[2] - DPerp*GradProfDPerp[2])*versors[2] * versors[0]
		+ (Dparprof - Dperpprof)*versors[2]*GradVersors['z'][0]
		+ (Dparprof - Dperpprof)*versors[0]*GradVersors['z'][2]   // dz alpha_xz
		;
	      
	      double uy =
		(Dpar*GradProfDPar[1] - DPerp*GradProfDPerp[1])*pow(versors[1],2)
		+ DPerp*GradProfDPerp[1]
		+ (Dparprof - Dperpprof)*2.0*versors[1]*GradVersors['y'][1]    // dy alpha_yy
		+ (Dpar*GradProfDPar[0] - DPerp*GradProfDPerp[0])*versors[0] * versors[1]
		+ (Dparprof - Dperpprof)*versors[0]*GradVersors['x'][1]
		+ (Dparprof - Dperpprof)*versors[1]*GradVersors['x'][0]   // dx alpha_xy
		+ (Dpar*GradProfDPar[2] - DPerp*GradProfDPerp[2])*versors[2] * versors[1]
		+ (Dparprof - Dperpprof)*versors[2]*GradVersors['z'][1]
		+ (Dparprof - Dperpprof)*versors[1]*GradVersors['z'][2]   // dz alpha_yz
		;
	      
	      double uz =
		(Dpar*GradProfDPar[2] - DPerp*GradProfDPerp[2])*pow(versors[2],2)
		+ DPerp*GradProfDPerp[2]
		+ (Dparprof - Dperpprof)*2.0*versors[2]*GradVersors['z'][2]    // dz alpha_zz
		+ (Dpar*GradProfDPar[1] - DPerp*GradProfDPerp[1])*versors[2] * versors[1]
		+ (Dparprof - Dperpprof)*versors[2]*GradVersors['y'][1]
		+ (Dparprof - Dperpprof)*versors[1]*GradVersors['y'][2]   // dy alpha_yz
		+ (Dpar*GradProfDPar[0] - DPerp*GradProfDPerp[0])*versors[2] * versors[0]
		+ (Dparprof - Dperpprof)*versors[2]*GradVersors['x'][0]
		+ (Dparprof - Dperpprof)*versors[0]*GradVersors['x'][2]   // dx alpha_xz
		;
	      
	      phix.push_back(ux);
	      phiy.push_back(uy);
	      phiz.push_back(uz);
	      
	    } // End of nuloop
	    
	  }
	}
      }
    }
    
  }
  
  return ;
}

double TDiffusionCoefficient3D::GetProfileDPar(double x, double y, double z, TSource* SourceTerm) {
  double radial = 1;
  double result = 1.;
  //    cout << " nrn_sn in GetProfile ->" << nrn_sn << endl;
  
  switch(set_profile) {
  case Constant :
    return 1.0;
    break;
    
  case Exp :
    return exp(fabs(z)/zt);
    break;
    
  case Qtau :
    radial = SourceTerm->GetSource(x,y,z);
    radial /= nrn_sn;
    return pow(radial, -index_radial);
    break;

  case Test:
    return exp(x/5.0);
    break; 
  default :
    return -1;
  }
}

double TDiffusionCoefficient3D::GetXDerivativeDPar(double x, double y, double z, TSource* SourceTerm) {
  
  double delta_x = 0.00002; //Coord->GetDeltaX();
  double xup = x + 0.5*delta_x;
  double xdo = x - 0.5*delta_x;
  return ( GetProfileDPar(xup,y,z, SourceTerm) - GetProfileDPar(xdo,y,z, SourceTerm) )/delta_x ;
}

double TDiffusionCoefficient3D::GetYDerivativeDPar(double x, double y, double z, TSource* SourceTerm) {
  
  double delta_y = 0.00002; //Coord->GetDeltaY();
  double yup = y + 0.5*delta_y;
  double ydo = y - 0.5*delta_y;
  return ( GetProfileDPar(x,yup,z, SourceTerm) - GetProfileDPar(x,ydo,z, SourceTerm) )/delta_y;
}

double TDiffusionCoefficient3D::GetZDerivativeDPar(double x, double y, double z, TSource* SourceTerm) {
  
  double delta_z = 0.00002; //Coord->GetDeltaZ();
  double zup = z + 0.5*delta_z;
  double zdo = z - 0.5*delta_z;
  return ( GetProfileDPar(x,y,zup, SourceTerm) - GetProfileDPar(x,y,zdo, SourceTerm) )/delta_z;
}

vector<double> TDiffusionCoefficient3D::GetXDerivativeBF(double x, double y, double z, TBField* SourceTerm) {
  
  double delta_x = 0.00002; //Coord->GetDeltaX();
  double xup = x + 0.5*delta_x;
  double xdo = x - 0.5*delta_x;
  vector<double> vers1 = SourceTerm->GetBregVersors(xup,y,z);
  vector<double> vers2 = SourceTerm->GetBregVersors(xdo,y,z);
  vector<double> final(3,0);
  final[0] = (vers1[0]-vers2[0])/delta_x;
  final[1] = (vers1[1]-vers2[1])/delta_x;
  final[2] = (vers1[2]-vers2[2])/delta_x;
  return final ;
}

vector<double> TDiffusionCoefficient3D::GetYDerivativeBF(double x, double y, double z, TBField* SourceTerm) {
  
  double delta_y = 0.00002; //Coord->GetDeltaY();
  double yup = y + 0.5*delta_y;
  double ydo = y - 0.5*delta_y;
  vector<double> vers1 = SourceTerm->GetBregVersors(x,yup,z);
  vector<double> vers2 = SourceTerm->GetBregVersors(x,ydo,z);
  vector<double> final(3,0);
  final[0] = (vers1[0]-vers2[0])/delta_y;
  final[1] = (vers1[1]-vers2[1])/delta_y;
  final[2] = (vers1[2]-vers2[2])/delta_y;
  return final ;
}

vector<double> TDiffusionCoefficient3D::GetZDerivativeBF(double x, double y, double z, TBField* SourceTerm) {
  
  double delta_z = 0.00002; //Coord->GetDeltaZ();
  double zup = z + 0.5*delta_z;
  double zdo = z - 0.5*delta_z;
  vector<double> vers1 = SourceTerm->GetBregVersors(x,y,zup);
  vector<double> vers2 = SourceTerm->GetBregVersors(x,y,zdo);
  vector<double> final(3,0);
  final[0] = (vers1[0]-vers2[0])/delta_z;
  final[1] = (vers1[1]-vers2[1])/delta_z;
  final[2] = (vers1[2]-vers2[2])/delta_z;
  return final ;
}


std::vector<double> TDiffusionCoefficient3D::GetGradProfileDPerp(double x, double y, double z, TSource* ST) {
  vector<double> prof(3,0);
  prof[0] = GetXDerivative(x,y,z,ST);
  prof[1] = GetYDerivative(x,y,z,ST);
  prof[2] = GetZDerivative(x,y,z,ST);
  return prof;
}

std::vector<double> TDiffusionCoefficient3D::GetGradProfileDPar(double x, double y, double z, TSource* ST) {
  vector<double> prof(3,0);
  prof[0] = GetXDerivativeDPar(x,y,z,ST);
  prof[1] = GetYDerivativeDPar(x,y,z,ST);
  prof[2] = GetZDerivativeDPar(x,y,z,ST);
  return prof;
  
}

std::map<char, std::vector<double> > TDiffusionCoefficient3D::GetGradVersors(double x, double y, double z, TBField* bf) {
  map<char, std::vector<double> > gradver;
  
  gradver['x'] = GetXDerivativeBF(x,y,z,bf);
  gradver['y'] = GetYDerivativeBF(x,y,z,bf);
  gradver['z']  = GetZDerivativeBF(x,y,z,bf);
  
  return gradver;
}
