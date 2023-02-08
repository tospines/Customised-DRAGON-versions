#include "xsec.h"
#include "kamae.h"


//namespace DRAGON {
//using namespace std;
//const double TSpallationNetwork::ethr = 5.e1;

  TSpallationNetwork::TSpallationNetwork(TGrid* coord, Input* in, std::vector<TXSecBase*> xsecmodel, std::vector<int>& nuclei) {
    
    std::cout << "\nXSecs option " << in->spallationxsec << std::endl;

    if (in->spallationxsec == Fluka){
      //spallXsecFile = "data/FlukaSpallInclusiveXsec_last_cumulative.dat";
      spallXsecFile = "data/Fluka_2020_DRAGON2.dat";
      //spallXsecFile = "data/Fluka_2020Renorm.dat";
      spallXsecFile_p = "data/Fluka2021_Proton_production.dat";   //"data/Fluka_Proton_production_new.dat";  
    }
    else if (in->spallationxsec == GalpropXSec){
      spallXsecFile = "data/galprop_eval_iso_cs_updated.dat";
      spallXsecFile_p = "";
    }
    else if (in->spallationxsec == Webber03){
      spallXsecFile = "data/webber_xsec_total.dat";
      spallXsecFile_p = "";}
    else if (in->spallationxsec == DRAGON2){
      spallXsecFile = "data/crxsecs_fragmentation_Evoli2019_cumulative_Complete_modified.dat";
      //spallXsecFile = "data/HXSec_Carmelo_Best.dat"; //--> Carmelo renormalized in all channels without scaling
      //spallXsecFile = "data/HXSec_Carmelo_Combined.dat";      
      //spallXsecFile = "data/HXSec_FLUKA_Combined.dat";      
      //spallXsecFile = "data/HXSec_GALPROP_Combined.dat";      
      //spallXsecFile = "data/HXSec_Galprop_even_scaled_Best.dat"; //--> Galprop renormalized in all channels scaling from 500MeV to 10GeV
      //spallXsecFile = "data/HXSec_Galprop_even_Best.dat"; //--> Galprop renormalized in all channels without scaling
      //spallXsecFile = "data/HXSec_Galprop_Best.dat"; //--> Galprop renormalized in the main channels
      //spallXsecFile = "data/Fluka_2020_DRAGON2.dat";
      spallXsecFile_p = "";}
    else {
      std::cout << "Wrong SpallationXSec option!! --> Taking the default one (DRAGON2 for nuclei, DRAGON treatment for Secondary protons) \n" << std::endl;
      spallXsecFile = "data/crxsecs_fragmentation_Evoli2019_cumulativeComplete_modified.dat";
      spallXsecFile_p = "";
    }

          
    DRAGONEnergyVector.clear();
    
    const double factor = Clight*1.e-27; 
    std::vector<double> beta = coord->GetBeta();
    DRAGONEnergyVector = coord->GetEk();
    const int dimEn = DRAGONEnergyVector.size();
    
    
    std::cout << std::endl << "Energy min " << *DRAGONEnergyVector.begin() << " Energy max " << DRAGONEnergyVector.back() << "\n" << "Energy entries " << DRAGONEnergyVector.size() << "\n" << std::endl;  
    
    if (in->spallationxsec == Fluka || in->spallationxsec == DRAGON2){
      
      std::cout << std::endl << "Reading spallation cross sections... " << std::endl;
      
      if(in->spallationxsec == Fluka) std::cout << " From Fluka table! \n" << std::endl;
      else if(in->spallationxsec == DRAGON2) std::cout << " From DRAGON2 table! \n" << std::endl;   
      else std::cout << " From default (DRAGON2 algorithm)" << std::endl;
      
      datafileSpall.open(spallXsecFile.c_str());
      if (!datafileSpall){
	std::cout << "SpallationXSec file not found!! --> Check the folder data in your setup for the file " << spallXsecFile << "\n" << std::endl;
	exit(-1);}
      std::cout << spallXsecFile << std::endl;
      
      
      std::string header;
      getline(datafileSpall, header); //skipping the header
      std::cout << "this is the header: " << header << std::endl;
      double parentUID_tableFormat;
      double daughterUID_tableFormat;
      int parentUID;
      int daughterUID;
      double kineticEnergy = 0.;
      double xsec;
      double xsec_He;
      std::string comment1;  	     
      
      int counter = 0;
      spallationKineticEnergyVector.clear();
      
      
      while (datafileSpall >> parentUID_tableFormat >> daughterUID_tableFormat >> kineticEnergy >> xsec >> xsec_He >> comment1) {
	
	//if (counter%100==0)
	//std::cout << parentUID_tableFormat << " " << daughterUID_tableFormat<< " " << kineticEnergy << " " << xsec << " " << xsec_He << " " << comment1.c_str() << std::endl;
	
	if (parentUID_tableFormat == 8.16 && daughterUID_tableFormat == 6.12) 
	  spallationKineticEnergyVector.push_back(kineticEnergy/1000.); //   * MeV --> CHECK UNITS!!! //The DragonEnergy is in GeV
	
	double fracPart, intPart=0;
	fracPart = modf(parentUID_tableFormat, &intPart);
	parentUID = int(intPart*1000. + fracPart*100.) ;  //1000*Z + A
	double fracPart2, intPart2=0;
	fracPart2 = modf(daughterUID_tableFormat, &intPart2);
	daughterUID = int(intPart2*1000. + fracPart2*100.) ;  //1000*Z + A
	
	// If nucleus is antiproton of leptons skip it. Skip also if mass(daughter) > mass(parent)
	if (daughterUID <= 1001 || intPart < intPart2) continue;    
	
	std::pair<int,int> couple1(parentUID, daughterUID);
	std::pair<double,double> couple2(xsec, xsec_He);
	spallationXsections[couple1].push_back(couple2);  // mbarn --> CHECK UNITS!!!                       
	
	//std::cout << couple1.first << " " << couple1.second << "    " << couple2.first << std::endl;
	counter++;
	
      } // end of the while loop
      datafileSpall.close(); 
      
      std::cout << "...spallation cross sections successfully read from file" << std::endl << std::endl;
      
      
    } // End if Carmelo or Fluka
    
    TGalpropXSec default_Galpropobject(coord, in);
    default_Galpropobject.set_sigma_cc();
    
    int ISS = -1;
    default_Galpropobject.sigtap_cc(ISS);

    
    
    if(in->spallationxsec == DRAGON2){ //Adding the He contribution to the DRAGON2 XSecs
      std::vector <double> AddedXSec_vec;
      
      for (int iloop = 0; iloop < nuclei.size()-1; ++iloop) {
	
	// If antiprotons or leptons, do not compute spallation. If protons, do not compute, because their spallation products will be computed elsewhere
	if (nuclei[iloop] <= 1001) continue;
        
	int iz = -1000;  // parent nucleus charge
	int ia = -1000;  // parent nucleus mass
	Utility::id_nuc(nuclei[iloop], ia, iz);
        
	for (int idaught = iloop+1; idaught < nuclei.size(); ++idaught) {
	  
	  int jz = -1000; // daughter nucleus charge
	  int ja = -1000; // daughter nucleus mass
	  Utility::id_nuc(nuclei[idaught], ja, jz);

	  // If nucleus is antiproton of leptons skip it. Skip also if mass(daughter) > mass(parent)
	  if (nuclei[idaught] < 1001 || ia < ja) continue;
	  
	  pair<int,int> couple(1000*iz+ia,1000*jz+ja);
	  if (jz > 2){ //> jz > 2 
	    
	    if (spallationXsections[couple].empty()){	      
	      for(unsigned int is = 0; is < spallationKineticEnergyVector.size(); ++is) {
		spallationXsections[couple].push_back(std::pair <double, double> (0., 0.));
	      }
	    }

	    else{

	      //AddedXSec_vec = default_Galpropobject.GetXSec(iz, ia, jz, ja, spallationXsections[couple], spallationKineticEnergyVector);
	      //std::cout << couple.first << " Added " << couple.second << "   " << AddedXSec_vec.back() << std::endl;
	      for(unsigned int is = 0; is < spallationKineticEnergyVector.size(); ++is) {
		//spallationXsections[couple][is].first = AddedXSec_vec[is];
		spallationXsections[couple][is].first = (spallationXsections[couple][is].first*(1.0 + He_abundance*default_Galpropobject.He_to_H_CS_ratio(spallationKineticEnergyVector[is],iz,ia,jz,ja)));

		//if (iz == 14 & ia == 32 & jz == 14 & ja == 30)
		//std::cout << spallationXsections[couple][is].first << " My final val " << spallationKineticEnergyVector[is] << " " << couple.first << " " << couple.second << std::endl;
	      }
	    }	 
	  }
	  else{ //For helium, tritium or deuterion creation we need Webber Spall XSecs
	    //spall[couple] = default_Galpropobject.GetXSec(iz, ia, jz, ja);
	    spall[couple] = TGalpropXSec(coord, in).GetXSec(iz, ia, jz, ja);
	    //for(unsigned int a = 0; a < dimEn; a++) 
	    //std::cout << spall[couple][a] << " Added " << iz << ia << "   " << jz << ja << std::endl;
	  }
	  
	}
	
      }

    } //End of He contribution in DRAGON2 XSecs

                 /////////////////////  PROTON XSECS  \\\\\\\\\\\\\\\\\\\

    int dimFlukaTables = 121; //160;  // Primary
    int dimprot = dimFlukaTables;
    int dimsec = 547; //160; // secondary
    double Ep[dimsec], Elept[dimsec], Eap[dimsec];
    double T[dimEn], ET[dimprot];  //dimEn is the DRAGON energy file -> energy of the produced protons
    double DBlog = 1./16.;   // Factor of change in FLUKA tables for primaries! 
    double DBlog2 = 1./64.;   // Factor of change in FLUKA tables for secondaries! 
    double E0prim = 1.e-3; // 1 MeV -> 1e-3 GeV    1e-2 for the past version
    double E0sec = 1.e-4; // 0.1 MeV -> 1e-4 GeV   1e-2 for the past version
    
    const double factorprot = factor*coord->GetDeltaE();  // DeltaE is the relative delta! dE/E
    const double factorelpos = 1.e3*factorprot;  // barn/GeV --- *10^3 ---> mbarn/GeV --- *factorprot ---> cm^3/s
    //std::cout << "Filling secondary protons table " << std::endl; 
    std::pair<int,int> coupleppr(1001, 1001);  // Secondary protons, from protons
    std::pair<int,int> couplepHe(2004, 1001);  // Secondary protons, from Helium
    
    if (in->spallationxsec != Fluka) {

      std::cout << "\nFilling secondary protons table"; 
      std::cout << " --> using DRAGON treatment of secondary protons" << std::endl;
      
      double xsec;
      spall[coupleppr] = vector<double>(dimEn, 0);

      for(unsigned int ip = 0; ip < dimEn; ++ip) {	
	PP_inel = 0.0;
	PA_inel = 0.0;
	aPP_non = 0.0;
	aPA_non = 0.0;
	aPP_ann = 0.0;
	aPA_ann = 0.0;
	default_Galpropobject.nucleon_cs(2, DRAGONEnergyVector[ip], 1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);  // CROSEC INEXS
	spall[coupleppr][ip] = in->SPALL*factorprot*(PP_inel + He_abundance*PA_inel);
       }
	std::cout << "\n Secondary protons OK!!\n" << std::endl;
    }
    
    else { //fluka      
      std::ifstream infile;
      infile.open(spallXsecFile_p.c_str());

      std::cout << "Opening secondary proton spallXsec file: " << spallXsecFile_p.c_str() << std::endl;
      if (!infile.is_open()){ std::cout << "problem opening Secproton file: "<< spallXsecFile_p << " !!!" << std::endl;
	exit(-1);}

      spall_apel[coupleppr] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
      spall_apel[couplepHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));

    
      double a,b,c,d;
      Matrix_p_pp   = std::vector<double>(dimsec*dimprot, 0.0);
      Matrix_p_pHe  = std::vector<double>(dimsec*dimprot, 0.0);
      Matrix_p_Hep  = std::vector<double>(dimsec*dimprot, 0.0);
      Matrix_p_HeHe = std::vector<double>(dimsec*dimprot, 0.0);

      for (int j = 0; j < dimprot; j++) {  //Eprim
	for (int i = 0; i < dimsec; i++) {    //Esec	  
	  infile >> a >> b >> c >> d;
	  
	  ind = j*dimsec + i;//index_matrix(i,j);
	  
	  Matrix_p_pp[ind] = a;
	  Matrix_p_pHe[ind] = b;
	  Matrix_p_Hep[ind] = c;
	  Matrix_p_HeHe[ind] = d;
	}
      }
      infile.close();
      std::cout << "Filling secondary protons table " << std::endl;       
      
      std::pair<int,int> coupleppr(1001, 1001);  // Secondary protons, from protons
      std::pair<int,int> couplepHe(2004, 1001);  // Secondary protons, from Helium
      
      std::vector <double> Eprimp, Esecp;
      for (int i = 0; i < dimprot; i++) { //Here we are filling at the beggining of the ApEl species all the energy vectors for FLUKA. If different
	  //tables were to have different dimensions or different energy configurations this is not correct anymore.
	//Ep[i] = pow(10., log10(E0sec) + double(i)*DBlog2); //+0.5*DBlog  The sum term is to position in the middle of the bin! But it's not needed since the interpolation procedure does not need it.
	ET[i] = pow(10., log10(E0prim) + double(i)*DBlog);
	Eprimp.push_back(ET[i]);
      }

      for (int i = 0; i < dimsec; i++) { Ep[i] = pow(10., log10(E0sec) + double(i)*DBlog2); //+0.5*DBlog      
	Esecp.push_back(Ep[i]);     } 

      
      double inx;
      std::vector <std::vector <double> > vec_ppF, vec_HepF;
      for (int e = 0; e < Eprimp.size(); e++){
	std::vector <double> vec_pp_F, vec_Hep_F;
	for (int u = 0; u < Esecp.size(); u++){
	  inx = e*Esecp.size() + u; // u*Eprimp.size() + e;  //This is due to the order of the Fluka table!!
	  vec_pp_F.push_back(Matrix_p_pp[inx] + He_abundance*Matrix_p_pHe[inx]);
	  vec_Hep_F.push_back(Matrix_p_Hep[inx] + He_abundance*Matrix_p_HeHe[inx]);
	}
	
	vec_ppF.push_back(vec_pp_F);
	vec_HepF.push_back(vec_Hep_F);     
      }

      std::vector < std::vector <double> > vec_pp_intF, vec_Hep_intF;
      for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
	std::vector <double> vec_pp_int_F, vec_Hep_int_F;
	for (int u = 0; u < Esecp.size(); u++){
	  int je = 0, jee;
	  double t;
	  if (DRAGONEnergyVector[e] > Eprimp.back()){  //Extrapolation
	    je = Eprimp.size()-1; 	  jee = je;
	    t = (DRAGONEnergyVector[e]  - Eprimp[jee])/(Eprimp[je]-Eprimp[je-1]);
	    if (DRAGONEnergyVector[e] > 7.*Eprimp.back())
	      t = (7.*Eprimp.back()  - Eprimp[jee])/(Eprimp[je]-Eprimp[je-1]);
	  }
	  else if ( DRAGONEnergyVector[e] < Eprimp.front()){
	    je = 1; jee = 0;
	    t =  (DRAGONEnergyVector[e] - Eprimp[jee])/(Eprimp[je]-Eprimp[je-1]);
	    if (Eprimp.front()/7. > DRAGONEnergyVector[e])
	      t = (Eprimp[0]/7. - Eprimp[jee])/(Eprimp[je]-Eprimp[je-1]);
	  }
	  
	  else{ //Interpolation
	    while (DRAGONEnergyVector[e] > Eprimp[je])
	      je++;
	    
	    if ((Eprimp[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprimp[je-1])) jee = je-1; 
	    else
	      jee = je;
	    t = (DRAGONEnergyVector[e]  - Eprimp[jee])/(Eprimp[je]-Eprimp[je-1]);
	  }
	  
	  vec_pp_int_F.push_back(std::max(vec_ppF[jee][u] + t*(vec_ppF[je][u] - vec_ppF[je-1][u]), 0.));
	  vec_Hep_int_F.push_back(std::max(vec_HepF[jee][e] + t*(vec_HepF[je][u] - vec_HepF[je-1][u]), 0.)); 
	}
	
	vec_pp_intF.push_back(vec_pp_int_F);
	vec_Hep_intF.push_back(vec_Hep_int_F);
      }


      for (int e = 0; e < dimEn; e++){
	for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
	  
	  int je = 0, jee;
	  double t;
	  if (DRAGONEnergyVector[u] > Esecp.back()){  //Extrapolation
	    je = Esecp.size()-1; 	  jee = je;
	    t =  (DRAGONEnergyVector[u]  - Esecp[jee])/(Esecp[je]-Esecp[je-1]);
	    if (DRAGONEnergyVector[u] > 7.*Esecp.back())
	      t = (7.*Esecp.back()  - Esecp[jee])/(Esecp[je]-Esecp[je-1]);
	  }
	  
	  else if ( DRAGONEnergyVector[u] < Esecp.front()){
	    je = 1; jee = 0;
	    t = (DRAGONEnergyVector[u] - Esecp[jee])/(Esecp[je]-Esecp[je-1]);
	    if (Esecp.front()/7. > DRAGONEnergyVector[u])
	      t = (Esecp[0]/7. - Esecp[jee])/(Esecp[je]-Esecp[je-1]);
	  }
	  
	  else{ //Interpolation
	    while (DRAGONEnergyVector[u] > Esecp[je])
	      je++;
	    
	    if ((Esecp[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esecp[je-1])) jee = je-1; 
	    else
	      jee = je;	
	    
	    t =  (DRAGONEnergyVector[u] - Esecp[jee])/(Esecp[je]-Esecp[je-1]);
	  }
	  
	  double ppF = std::max(vec_pp_intF[e][jee] + t*(vec_pp_intF[e][je] - vec_pp_intF[e][je-1]), 0.);
	  double HepF = std::max(vec_Hep_intF[e][jee] + t*(vec_Hep_intF[e][je] - vec_Hep_intF[e][je-1]), 0.);

	  spall_apel[coupleppr][u][e] = DRAGONEnergyVector[e]*factorelpos * ppF;
	  spall_apel[couplepHe][u][e] = 4.*DRAGONEnergyVector[e]*factorelpos * HepF; //Factor 4 because we integrate for energies (not E/nucleon) of the primary particle! see eqs. 3.10 and 3.11 of DRAGON2-2
	}
	
      }
      //for (int e = 0; e < dimEn; e++){
      //std::cout << spall_apel[coupleppr][e][100] << " cross sections, at " << DRAGONEnergyVector[e] << "GeV of secondary proton and " << DRAGONEnergyVector[100] << "GeV for primary!" << endl;}
      //exit(10);
      
    }  
    

    // If antiprotons and/or leptons are wanted in output, add them, and add also secondary protons
    
    if(in->prop_ap || in->prop_lep || in->prop_Adeuteron || in->prop_AHe3) {  //if(in->prop_ap || in->prop_lep || in->prop_Adeuteron || in->prop_AHe3 || in->prop_AHe4)
      vector<double> pp = coord->GetMomentum();
      size_t limit = 100000; 
      gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
                 
      if(in->prop_ap) {

	std::pair<int,int> coupleapap(-999, -999); // Tertiary antiprotons
	spall_apel[coupleapap] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	
	std::pair<int,int> coupleappr(1001,-999);  // Secondary antiprotons, from protons
	std::pair<int,int> coupleapHe(2004,-999);  // Secondary antiprotons, from Helium
	spall_apel[coupleappr] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	
	std::pair<int,int> coupleapC12(6012,-999);  // Secondary antiprotons, from C12
	std::pair<int,int> coupleapC13(6013,-999);  // Secondary antiprotons, from C13
	std::pair<int,int> coupleapO16(8016,-999);  // Secondary antiprotons, from O
	std::pair<int,int> coupleapN14(7014,-999);  // Secondary antiprotons, from N14
	std::pair<int,int> coupleapN15(7015,-999);  // Secondary antiprotons, from N15
	
	std::pair<int,int> coupleapLi6(3006,-999);  // Secondary antiprotons, from Li6
	std::pair<int,int> coupleapLi7(3007,-999);  // Secondary antiprotons, from Li7
	std::pair<int,int> coupleapBe7(4007,-999);  // Secondary antiprotons, from Be7
	std::pair<int,int> coupleapBe9(4009,-999);  // Secondary antiprotons, from Be9
	std::pair<int,int> coupleapBe10(4010,-999);  // Secondary antiprotons, from Be10
	std::pair<int,int> coupleapB10(5010,-999);  // Secondary antiprotons, from B10
	std::pair<int,int> coupleapB11(5011,-999);  // Secondary antiprotons, from B11
	
	std::pair<int,int> coupleapNe(10020,-999);  // Secondary antiprotons, from Ne
	std::pair<int,int> coupleapMg(12024,-999);  // Secondary antiprotons, from Mg
	std::pair<int,int> coupleapSi(14028,-999);  // Secondary antiprotons, from Si
	spall_apel[coupleapC12] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapC13] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));  
	spall_apel[coupleapN14] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));  
	spall_apel[coupleapN15] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));  
	spall_apel[coupleapO16] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));  
	spall_apel[coupleapNe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));  
	spall_apel[coupleapMg] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapSi] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));  
	
	spall_apel[coupleapLi6] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapLi7] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapBe7] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapBe9] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapBe10] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapB10] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleapB11] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));  
	
	
	if (in->apy == Winkler){
	  std::cout << "Reading sec antiprotons tables method Winkler... " << std::endl; 	
	  TSpallationNetwork::InitXSecWinkler(factorelpos, in);
	  std::cout << "\n Winkler secondary antiprotons OK!! \n" << std::endl;
	  std::cout << "Filling tert antiprotons table DRAGON procedure ... " << std::endl;
	  
	  spall[coupleapap] = vector<double>(dimEn, 0);
	  for(unsigned int ip = 0; ip < dimEn; ++ip) {
	    PP_inel = 0.0;
	    PA_inel = 0.0;
	    aPP_non = 0.0;
	    aPA_non = 0.0;
	    aPP_ann = 0.0;
	    aPA_ann = 0.0;
	    default_Galpropobject.nucleon_cs(2, DRAGONEnergyVector[ip], -1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);
	    spall[coupleapap][ip] = factorprot*( aPP_non + He_abundance*aPA_non );
	    //std::cout << "\n This: " << aPP_non << " " << aPA_non << " " << spall[coupleapap][ip] << std::endl;
	  }
	}
	
	else if (in->apy == GalpropFunction){ //Sec and Tertiary XS from Galprop
	  std::cout << "Filling sec antiprotons table method Galprop... " << std::endl; 	
	  for (unsigned int i = 0; i < dimEn; i++) {
	    for (unsigned int ii=i+1; ii < dimEn; ii++) { //ii is the progenitor energy!
	      
	      spall_apel[coupleappr][i][ii] = factorelpos*DRAGONEnergyVector[ii]*(default_Galpropobject.antiproton_cc1(w,limit,in->antiproton_cs, pp[i], pp[ii], 1, 1, 1, 1) * ( (!in->scaling) + (in->scaling)*(0.12 * pow(DRAGONEnergyVector[i], -1.67) + 1.78)) + (!in->scaling)*He_abundance*default_Galpropobject.antiproton_cc1(w,limit,in->antiproton_cs, pp[i], pp[ii], 1, 1, 2, 4)); 
	      
	      spall_apel[coupleapHe][i][ii] = factorelpos*DRAGONEnergyVector[ii]*4.0*(default_Galpropobject.antiproton_cc1(w,limit,in->antiproton_cs, pp[i], 4.0*pp[ii], 2, 4, 1, 1) * ( (!in->scaling) + (in->scaling)*(0.12 * pow(DRAGONEnergyVector[i], -1.67) + 1.78)) + (!in->scaling)*He_abundance*default_Galpropobject.antiproton_cc1(w,limit,in->antiproton_cs, pp[i], 4.0*pp[ii], 2, 4, 2, 4));
	      
	    }
	  }
	  
	  std::cout << "\n Secondary antiprotons OK!!\n " << std::endl;   
	  std::cout << "Filling tert antiprotons table DRAGON procedure ... " << std::endl;
	  spall[coupleapap] = vector<double>(dimEn, 0);
	  
	  for(unsigned int ip = 0; ip < dimEn; ++ip) {	   
	    PP_inel = 0.0;
	    PA_inel = 0.0;
	    aPP_non = 0.0;
	    aPA_non = 0.0;
	    aPP_ann = 0.0;
	    aPA_ann = 0.0;
	    default_Galpropobject.nucleon_cs(2, DRAGONEnergyVector[ip], -1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);
	    spall[coupleapap][ip] = factorprot*( aPP_non + He_abundance*aPA_non );
	    //std::cout << "\n This: " << aPP_non << " " << aPA_non << " " << spall[coupleapap][ip] << std::endl;

	    // Exploits approximation d\sigma/dEkin \propto 1/Ekin' [Tan & Ng '83]
	    //xset = factorprot*( aPP_non + He_abundance*aPA_non );
	  }    	  
	}
	
	else if (in->apy == FlukaAp){
	  spallXsecFile_ap = "data/Fluka2021_AProton_production.dat";  //"data/Fluka_Antiproton_production_new.dat";
	  spallXsecFile_Tap = "data/Fluka_TertiaryAntiproton_production_new.dat";
	
	  std::ifstream infile;
	  infile.open(spallXsecFile_ap.c_str());
	  std::cout << "Opening secondary antiproton spallXsec file " << spallXsecFile_ap.c_str() << std::endl;
	  if (!infile.is_open()){ std::cout << "problem opening SecAp Fluka file!!"<< std::endl;
	    exit(-1);}
	  double a,b,c,d;
	  Matrix_Ap_pp   = std::vector<double>(dimsec*dimprot, 0.0);
	  Matrix_Ap_pHe  = std::vector<double>(dimsec*dimprot, 0.0);
	  Matrix_Ap_Hep  = std::vector<double>(dimsec*dimprot, 0.0);
	  Matrix_Ap_HeHe = std::vector<double>(dimsec*dimprot, 0.0);


	  for (int i = 0; i < dimprot; i++) {
	    for (int j = 0; j < dimsec; j++) {
	      
	      infile >> a >> b >> c >> d;
	      
	      ind = i * dimsec + j;//index_matrix(i,j);
	      
	      Matrix_Ap_pp[ind] = a;
	      Matrix_Ap_pHe[ind] = b;
	      Matrix_Ap_Hep[ind] = c;
	      Matrix_Ap_HeHe[ind] = d;
	    }
	  }
	  
	  infile.close();
	  std::cout << "Filling sec antiprotons table " << std::endl; 	  
	  //std::ofstream outfileFluk;	  
	  //outfileFluk.open ("FlukaApXSec.txt");
	  //outfileFluk << "PrimEnergy(GeV)\t" << "Ap_Energy(GeV)\t" << "XSec_p+H (mbarn)\t" << "XSec_p+He (mbarn)\t" << "total_p+gas (mbarn)\t" << "XSec_He+H (mbarn)\t" << "XSec_He+He (mbarn)\t" << "total_He+gas (mbarn)\t\n";

	  std::vector <double> Eprim, Esec;
	  for (int j = 0; j < dimsec; j++){     Eap[j] = pow(10., log10(E0sec)+double(j)*DBlog2);  //+0.5*DBlog
	    Esec.push_back(Eap[j]);      }
          for (int ies = 0; ies < dimprot; ies++){   Eprim.push_back(ET[ies]);    }
	  

	  double inx;
	  std::vector <std::vector <double> > vec_papF, vec_HeapF;
	  for (int e = 0; e < Eprim.size(); e++){
	    std::vector <double> vec_pap_F, vec_Heap_F;
	    for (int u = 0; u < Esec.size(); u++){
	      inx = e*Esec.size() + u;  //This is due to the order of the Fluka table!! 
	      vec_pap_F.push_back(Matrix_Ap_pp[inx] + He_abundance*Matrix_Ap_pHe[inx]);
	      vec_Heap_F.push_back(Matrix_Ap_Hep[inx] + He_abundance*Matrix_Ap_HeHe[inx]);
	      
	    }
	    
	    vec_papF.push_back(vec_pap_F);
	    vec_HeapF.push_back(vec_Heap_F);     
	  }
	  
	  std::vector < std::vector <double> > vec_pap_intF, vec_Heap_intF;
	  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
	    std::vector <double> vec_pap_int_F, vec_Heap_int_F;
	    for (int u = 0; u < Esec.size(); u++){
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[e] > Eprim.back()){  //Extrapolation
		je = Eprim.size()-1; 	  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
		if (DRAGONEnergyVector[e] > 7.*Eprim.back())
		  t = (7.*Eprim.back()  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      else if ( DRAGONEnergyVector[e] < Eprim.front()){
		je = 1; jee = 0;
		t =  (DRAGONEnergyVector[e] - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
		if (Eprim.front()/7. > DRAGONEnergyVector[e])
		  t = (Eprim[0]/7. - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      
	      else{ //Interpolation
		while (DRAGONEnergyVector[e] > Eprim[je])
		  je++;
		
		if ((Eprim[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim[je-1])) jee = je-1; 
		else
		  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      
	      vec_pap_int_F.push_back(std::max(vec_papF[jee][u] + t*(vec_papF[je][u] - vec_papF[je-1][u]), 0.));
	      vec_Heap_int_F.push_back(std::max(vec_HeapF[jee][u] + t*(vec_HeapF[je][u] - vec_HeapF[je-1][u]), 0.));	
	    }
	    
	    vec_pap_intF.push_back(vec_pap_int_F);
	    vec_Heap_intF.push_back(vec_Heap_int_F);
	  }
	  
	  for (int e = 0; e < dimEn; e++){
	    
	    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
	      
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[u] > Esec.back()){  //Extrapolation
		je = Esec.size()-1; 	  jee = je;
		t =  (DRAGONEnergyVector[u]  - Esec[jee])/(Esec[je]-Esec[je-1]);
		if (DRAGONEnergyVector[u] > 7.*Esec.back())
		  t = (7.*Esec.back()  - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
	      
	      else if ( DRAGONEnergyVector[u] < Esec.front()){
		je = 1; jee = 0;
		t = (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
		if (Esec.front()/7. > DRAGONEnergyVector[u])
		  t = (Esec[0]/7. - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
      
	      else{ //Interpolation
		while (DRAGONEnergyVector[u] > Esec[je])
		  je++;
	
		if ((Esec[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec[je-1])) jee = je-1; 
		else
		  jee = je;	
		
		t =  (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
	     	      
	      double papF = std::max(vec_pap_intF[e][jee] + t*(vec_pap_intF[e][je] - vec_pap_intF[e][je-1]), 0.);
	      double HeapF = std::max(vec_Heap_intF[e][jee] + t*(vec_Heap_intF[e][je] - vec_Heap_intF[e][je-1]), 0.);

	      //std::cout << "XS: " <<papF << ", upper XS: " << vec_pap_intF[e][je] << ", lower XS: " << vec_pap_intF[e][je-1] << ", closest XS: " << vec_pap_intF[e][jee] << std::endl;

	      
	      spall_apel[coupleappr][u][e] = DRAGONEnergyVector[e]*factorelpos * papF;
	      spall_apel[coupleapHe][u][e] = 4.*DRAGONEnergyVector[e]*factorelpos * HeapF;

      	      //outfileFluk << DRAGONEnergyVector[j] << "\t\t" << DRAGONEnergyVector[i] << "\t\t" << cs_pp*1000. << "\t\t" << cs_pHe*1000. << "\t\t" << 1000.*(cs_pp + He_abundance*cs_pHe) << "\t\t" << 1000.*cs_Hep << "\t\t" << 1000.*cs_HeHe << "\t\t" << 1000.*(cs_Hep + He_abundance*cs_HeHe) << "\n";
	      
	    }
	    //outfileFluk << "\n\n\n";
	    //std::cout << " " <<  spall_apel[coupleapHe][e][60] << " " << spall_apel[coupleappr][e][70] << "  " <<  std::endl;  
	  }
	  
	  
	  //outfileFluk.close();
	  std::cout << "\n Secondary Antiprotons OK!!\n" << std::endl;	
	  
	  
	  //std::ifstream infile;
	  infile.open(spallXsecFile_Tap.c_str());
	  std::cout << "Opening tertiary antiproton spallXsec file " << spallXsecFile_Tap.c_str() << std::endl;
	  if (!infile.is_open()){ std::cout << "problem opening TertAp file!!"<< std::endl;
	    exit(-1);}
	  //a = 0;
	  //b = 0;
	  Matrix_3Ap_app   = std::vector<double>(160*160); //, 0.0); //dimatprot*dimprot
	  Matrix_3Ap_apHe  = std::vector<double>(160*160); //, 0.0);
	  
	  for (int i = 0; i < 160; i++) {
	    for (int j = 0; j < 160; j++) {
	      infile >> a >> b;
	      ind = i*dimprot + j;//index_matrix(i,j);
	      Matrix_3Ap_app.push_back(a);		    
	      Matrix_3Ap_apHe.push_back(b);
	    }
	  }
	  
	  infile.close();
	  std::cout << "Filling tert antiprotons table " << std::endl;   
	  
	  std::vector <std::vector <double> > vec_apapF;
	  for (int e = 0; e < 160; e++){
	    std::vector <double> vec_apap_F;
	    for (int u = 0; u < 160; u++){
	      inx = u*160 + e;  //This is due to the order of the Fluka table!! 
	      vec_apap_F.push_back(Matrix_3Ap_app[inx] + He_abundance*Matrix_3Ap_apHe[inx]);	      
	    }
	    
	    vec_apapF.push_back(vec_apap_F);
	  }
	  
	  std::vector <double> Eprim_ap, Esec_ap;
	  for (int j = 0; j < 160; j++){  
	    Esec_ap.push_back(pow(10., log10(1.e-2)+double(j)*(1./16.) ));
	    Eprim_ap.push_back(pow(10., log10(1.e-2)+double(j)*(1./16.) ));
	  }

	  
	  std::vector < std::vector <double> > vec_apap_intF;
	  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
	    std::vector <double> vec_apap_int_F;
	    for (int u = 0; u < Esec_ap.size(); u++){
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[e] > Eprim_ap.back()){  //Extrapolation
		je = Eprim_ap.size()-1; 	  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim_ap[jee])/(Eprim_ap[je]-Eprim_ap[je-1]);
		if (DRAGONEnergyVector[e] > 7.*Eprim_ap.back())
		  t = (7.*Eprim_ap.back()  - Eprim_ap[jee])/(Eprim_ap[je]-Eprim_ap[je-1]);
	      }
	      else if ( DRAGONEnergyVector[e] < Eprim_ap.front()){
		je = 1; jee = 0;
		t =  (DRAGONEnergyVector[e] - Eprim_ap[jee])/(Eprim_ap[je]-Eprim_ap[je-1]);
		if (Eprim_ap.front()/7. > DRAGONEnergyVector[e])
		  t = (Eprim_ap[0]/7. - Eprim_ap[jee])/(Eprim_ap[je]-Eprim_ap[je-1]);
	      }
	      
	      else{ //Interpolation
		while (DRAGONEnergyVector[e] > Eprim_ap[je])
		  je++;
		
		if ((Eprim_ap[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim_ap[je-1])) jee = je-1; 
		else
		  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim_ap[jee])/(Eprim_ap[je]-Eprim_ap[je-1]);
	      }
	      
	      vec_apap_int_F.push_back(std::max(vec_apapF[jee][u] + t*(vec_apapF[je][u] - vec_apapF[je-1][u]), 0.));
	    }
	    
	    vec_apap_intF.push_back(vec_apap_int_F);
	  }
	  
	  for (int e = 0; e < dimEn; e++){
	    
	    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
	      
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[u] > Esec_ap.back()){  //Extrapolation
		je = Esec_ap.size()-1; 	  jee = je;
		t =  (DRAGONEnergyVector[u]  - Esec_ap[jee])/(Esec_ap[je]-Esec_ap[je-1]);
		if (DRAGONEnergyVector[u] > 7.*Esec_ap.back())
		  t = (7.*Esec_ap.back()  - Esec_ap[jee])/(Esec_ap[je]-Esec_ap[je-1]);
	      }
	      
	      else if ( DRAGONEnergyVector[u] < Esec_ap.front()){
		je = 1; jee = 0;
		t = (DRAGONEnergyVector[u] - Esec_ap[jee])/(Esec_ap[je]-Esec_ap[je-1]);
		if (Esec_ap.front()/7. > DRAGONEnergyVector[u])
		  t = (Esec_ap[0]/7. - Esec_ap[jee])/(Esec_ap[je]-Esec_ap[je-1]);
	      }
	      
	      else{ //Interpolation
		while (DRAGONEnergyVector[u] > Esec_ap[je])
		  je++;
		
		if ((Esec_ap[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec_ap[je-1])) jee = je-1; 
		else
		  jee = je;	
		
		t =  (DRAGONEnergyVector[u] - Esec_ap[jee])/(Esec_ap[je]-Esec_ap[je-1]);
	      }
	     	      
	      double apapF = std::max(vec_apap_intF[e][jee] + t*(vec_apap_intF[e][je] - vec_apap_intF[e][je-1]), 0.);

	      spall_apel[coupleapap][u][e] = DRAGONEnergyVector[e]*factorelpos * apapF;

	    }
	    
	  }
	} //End tertAntip 
	std::cout << "Tertiary Antiprotons OK!!\n" << std::endl;	  
	
      } // ap
      
      
      if (in->prop_lep) {
	std::pair<int,int> coupleel(1001, -1000); // Electrons from protons
	std::pair<int,int> couplepos(1001, 1000); // Positrons from protons
	
	std::pair<int,int> coupleelHe(2004, -1000); // Electrons from He
	std::pair<int,int> coupleposHe(2004, 1000); // Positrons from He
	
	std::pair<int,int> coupleelC12(6012, -1000); // Electrons from Carbon
	std::pair<int,int> coupleposC12(6012, 1000); // Positrons from Carbon
	std::pair<int,int> coupleelC13(6013, -1000); // Electrons from Carbon13
	std::pair<int,int> coupleposC13(6013, 1000); // Positrons from Carbon13
	std::pair<int,int> coupleelO16(8016, -1000); // Electrons from Oxygen
	std::pair<int,int> coupleposO16(8016, 1000); // Positrons from Oxygen
	std::pair<int,int> coupleelN14(7014, -1000); // Electrons from Nitrogen14
	std::pair<int,int> coupleposN14(7014, 1000); // Positrons from Nitrogen14
	std::pair<int,int> coupleelN15(7015, -1000); // Electrons from Nitrogen15
	std::pair<int,int> coupleposN15(7015, 1000); // Positrons from Nitrogen15
	std::pair<int,int> coupleelNe(10020, -1000); // Electrons from Neon
	std::pair<int,int> coupleposNe(10020, 1000); // Positrons from Neon
	std::pair<int,int> coupleelMg(12024, -1000); // Electrons from Magnesium
	std::pair<int,int> coupleposMg(12024, 1000); // Positrons from Magnesium
	std::pair<int,int> coupleelSi(14028, -1000); // Electrons from Silicon
	std::pair<int,int> coupleposSi(14028, 1000); // Positrons from Silicon
	std::pair<int,int> coupleelB10(5010, -1000); // Electrons from B10
	std::pair<int,int> coupleposB10(5010, 1000); // Positrons from B10
	std::pair<int,int> coupleelB11(5011, -1000); // Electrons from B11
	std::pair<int,int> coupleposB11(5011, 1000); // Positrons from B11
	std::pair<int,int> coupleelBe10(4010, -1000); // Electrons from Be10
	std::pair<int,int> coupleposBe10(4010, 1000); // Positrons from Be10
	std::pair<int,int> coupleelBe9(4009, -1000); // Electrons from Be9
	std::pair<int,int> coupleposBe9(4009, 1000); // Positrons from Be9
	std::pair<int,int> coupleelBe7(4007, -1000); // Electrons from Be7
	std::pair<int,int> coupleposBe7(4007, 1000); // Positrons from Be7
	std::pair<int,int> coupleelLi7(3007, -1000); // Electrons from Nitrogen14
	std::pair<int,int> coupleposLi7(3007, 1000); // Positrons from Nitrogen14
	std::pair<int,int> coupleelLi6(3006, -1000); // Electrons from Nitrogen14
	std::pair<int,int> coupleposLi6(3006, 1000); // Positrons from Nitrogen14
	
	
	spall_apel[coupleel] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
        spall_apel[couplepos] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));

	spall_apel[coupleelHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
        spall_apel[coupleposHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));

	spall_apel[coupleelC12] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelC13] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelN14] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelN15] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelO16] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0)); 
	spall_apel[coupleelNe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelMg] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelSi] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelB10] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelB11] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelBe10] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelBe9] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelBe7] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelLi7] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleelLi6] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	

	spall_apel[coupleposC12] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposC13] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposN14] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposN15] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposO16] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposNe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposMg] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposSi] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposB10] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposB11] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposBe10] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposBe9] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposBe7] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposLi7] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	
	spall_apel[coupleposLi6] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));	

      
	if (in->ly == GalpropTable) {
	  std::cout << "\nGalprop leptons XSecs, loading..." << std::endl;
	  
	  spallXsecFile_pos = "data/ElTablefile";//Galprop_Positron_production_new.dat";
	  spallXsecFile_el = "data/PosTablefile";//Galprop_Electron_production_new.dat";

	  Nelectrons = 160;//401;
	  Nprotons = 160;//801;

	  // TSpallationNetwork::InitXSecGalprop(factorelpos);
	}//End Elpos Galprop
	
	else if(in->ly == Pohl) TSpallationNetwork::InitXSecPohl(factorelpos);
	else if (in->ly == Kamae){std::cout << "Kamae leptons " << std::endl;  TSpallationNetwork::InitXSecKamae(factorelpos); 	}
	if (in->ly == LeptUser){
	  std::cout << "Reading sec leptons cross sections from User defined table... " << std::endl; 	
	  TSpallationNetwork::InitLeptUser(factorelpos, in);
	  std::cout << "\n User-defined secondary leptons OK!! \n" << std::endl;
	}
	
       	else if (in->ly == FlukaLep) {
	  std::cout << "\nFluka leptons XSecs, loading..." << std::endl;

	  spallXsecFile_el = "data/Fluka2021_Electron_production.dat";  //"data/Fluka_Electron_production_new.dat";
	  spallXsecFile_pos = "data/Fluka2021_Positron_production.dat";  "data/Fluka_Positron_production_new.dat";

	  Nelectrons = dimsec; //547; //160;
	  Nprotons = dimprot; //121; //160;

	}
	else{
	  std::cout << " Taking the default Kamae cross sections \n" << std::endl;
	  TSpallationNetwork::InitXSecKamae(factorelpos);
	}

	if(in->ly == FlukaLep || in->ly == GalpropTable){
	  std::cout << "Opening electron production file " << spallXsecFile_el.c_str() << std::endl;

	  /*
	  Nelectrons = 547; //160; dimsec;
	  Nprotons = 121; //160;  dimprot;
	  DBlog = 1./16.;   // Factor of change in FLUKA tables for primaries! 
	  DBlog2 = 1./64.;   // Factor of change in FLUKA tables for secondaries!
	  E0prim = 1.e-3; //1.e-3; // 1 MeV -> 1e-3 GeV    1e-2 for the past version
	  E0sec = 1.e-4; //1.e-4; // 0.1 MeV -> 1e-4 GeV   1e-2 for the past version
	  */
	  //double ETl[Nprotons], Elept2[Nelectrons];
          //for (int j = 0; j < Nprotons; j++)    ETl[j] = pow(10, log10(E0prim) + double(j)*DBlog);
	  for (int j = 0; j < Nelectrons; j++)     Elept[j] = pow(10., log10(E0sec)+double(j)*DBlog2);  //+0.5*DBlog
	  
	  std::ifstream infile;
	  infile.open(spallXsecFile_el.c_str());
	  if (!infile.is_open()){ std::cout << "problem opening spall_el file!!"<< std::endl;
	    exit(-1);}
	  double a,b,c,d;
	  
	  Matrix_El_pp = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  Matrix_El_pHe = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  Matrix_El_Hep = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  Matrix_El_HeHe = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  
	  Matrix_Pos_pp = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  Matrix_Pos_pHe = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  Matrix_Pos_Hep = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  Matrix_Pos_HeHe = std::vector<double>(Nelectrons*Nprotons, 0.0);
	  
	  for (int j = 0; j < Nprotons; j++) {
	    for (int i = 0; i < Nelectrons; i++) {    
	      infile >> a >> b >> c >> d;
	      ind = j*Nelectrons + i; //index_matrix(i,j);
	      
	      Matrix_El_pp[ind] = a;
	      Matrix_El_pHe[ind] = b;
	      Matrix_El_Hep[ind] = c;
	      Matrix_El_HeHe[ind] = d;
	    }
	  }
	  
	  infile.close();
	  
	  infile.open(spallXsecFile_pos.c_str());
	  std::cout << "Opening positron production file " << spallXsecFile_pos.c_str() << std::endl;
	  if (!infile.is_open()){ std::cout << "problem opening spall_pos file!!"<< std::endl;
	    exit(-1);}
	  
	  for (int j = 0; j < Nprotons; j++) {
	    for (int i = 0; i < Nelectrons; i++) {
	      
	      infile >> a >> b >> c >> d;
	      ind = j*Nelectrons + i; //index_matrix(i,j);
	      Matrix_Pos_pp[ind] = a;
	      Matrix_Pos_pHe[ind] = b;
	      Matrix_Pos_Hep[ind] = c;
	      Matrix_Pos_HeHe[ind] = d;
	    }
	  }
	  infile.close();
	  
	  int dimlept = Nelectrons; //160; //Be careful with this, it can be different depending on the table you have!
	  int dimpr = Nprotons; 
	  //int i_lept;
	  
	  std::cout << "filling leptons... " << std::endl;
	  //std::ofstream outfileel;
	  //std::ofstream outfilepos;
	  
	  //outfileel.open ("FlukaElXSec.txt");
	  //outfileel << "PrimEnergy(GeV)\t" << "El_Energy(GeV)\t" << "XSec_p+H (barn)\t" << "XSec_p+He (barn)\t" << "total_p+gas (barn)\t" << "XSec_He+H (barn)\t" << "XSec_He+He (barn)\t" << "total_He+gas (barn)\t\n";
	  
	  //outfilepos.open ("FlukaPosXSec.txt");
	  //outfilepos << "PrimEnergy(GeV)\t" << "Pos_Energy(GeV)\t" << "XSec_p+H (barn)\t" << "XSec_p+He (barn)\t" << "total_p+gas (barn)\t" << "XSec_He+H (barn)\t" << "XSec_He+He (barn)\t" << "total_He+gas (barn)\t\n";
	  
	  std::vector <double> Eprim, Esec;
	  //Eprim.insert(Eprim.begin(), ET[0], ET[dimprot-1]);
	  //Esec.insert(Esec.begin(), Eap[0], Eap[dimprot-1]);
          for (int ies = 0; ies < dimpr; ies++)   Eprim.push_back(ET[ies]);
	  for (int ies = 0; ies < dimlept; ies++) Esec.push_back(Elept[ies]);
          
	  
	  double inx;
	  std::vector <std::vector <double> > vec_pel, vec_Heel;
	  for (int e = 0; e < Eprim.size(); e++){
	    std::vector <double> vec_pel_, vec_Heel_;
	    for (int u = 0; u < Esec.size(); u++){
	      inx = e*Esec.size() + u;  //This is due to the order of the Fluka table!!
	      vec_pel_.push_back(Matrix_El_pp[inx] + He_abundance*Matrix_El_pHe[inx]);
	      vec_Heel_.push_back(Matrix_El_Hep[inx] + He_abundance*Matrix_El_HeHe[inx]);
	    }	    
	    vec_pel.push_back(vec_pel_);
	    vec_Heel.push_back(vec_Heel_);     
	  }
	  
	  std::vector < std::vector <double> > vec_pel_int, vec_Heel_int;
	  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
	    std::vector <double> vec_pel_int_, vec_Heel_int_;
	    for (int u = 0; u < Esec.size(); u++){
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[e] > Eprim.back()){  //Extrapolation
		je = Eprim.size()-1; 	  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
		if (DRAGONEnergyVector[e] > 7.*Eprim.back())
		  t = (7.*Eprim.back()  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      else if ( DRAGONEnergyVector[e] < Eprim.front()){
		je = 1; jee = 0;
		t =  (DRAGONEnergyVector[e] - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
		if (Eprim.front()/7. > DRAGONEnergyVector[e])
		  t = (Eprim[0]/7. - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      
	      else{ //Interpolation
		while (DRAGONEnergyVector[e] > Eprim[je])
		  je++;
		
		if ((Eprim[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim[je-1])) jee = je-1; 
		else
		  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      
	      vec_pel_int_.push_back(std::max(vec_pel[jee][u] + t*(vec_pel[je][u] - vec_pel[je-1][u]), 0.));
	      vec_Heel_int_.push_back(std::max(vec_Heel[jee][u] + t*(vec_Heel[je][u] - vec_Heel[je-1][u]), 0.));	
	    }
	    
	    vec_pel_int.push_back(vec_pel_int_);
	    vec_Heel_int.push_back(vec_Heel_int_);
	  }
	  
	  for (int e = 0; e < dimEn; e++){
	    
	    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
	      
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[u] > Esec.back()){  //Extrapolation
		je = Esec.size()-1; 	  jee = je;
		t =  (DRAGONEnergyVector[u]  - Esec[jee])/(Esec[je]-Esec[je-1]);
		if (DRAGONEnergyVector[u] > 7.*Esec.back())
		  t = (7.*Esec.back()  - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
	      
	      else if ( DRAGONEnergyVector[u] < Esec.front()){
		je = 1; jee = 0;
		t = (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
		if (Esec.front()/7. > DRAGONEnergyVector[u])
		  t = (Esec[0]/7. - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
      
	      else{ //Interpolation
		while (DRAGONEnergyVector[u] > Esec[je])
		  je++;
	
		if ((Esec[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec[je-1])) jee = je-1; 
		else
		  jee = je;	
		
		t =  (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
	     	      
	      double pel = std::max(vec_pel_int[e][jee] + t*(vec_pel_int[e][je] - vec_pel_int[e][je-1]), 0.);
	      double Heel = std::max(vec_Heel_int[e][jee] + t*(vec_Heel_int[e][je] - vec_Heel_int[e][je-1]), 0.);

	      //std::cout << "XS: " <<papF << ", upper XS: " << vec_pap_intF[e][je] << ", lower XS: " << vec_pap_intF[e][je-1] << ", closest XS: " << vec_pap_intF[e][jee] << std::endl;

	      spall_apel[coupleel][u][e] = DRAGONEnergyVector[e]*factorelpos * pel;
	      spall_apel[coupleelHe][u][e] = 4.*DRAGONEnergyVector[e]*factorelpos * Heel;

	      
	      //outfileel << Epr << "\t\t" << Eel << "\t\t" << cs_pp << "\t\t" << cs_pHe << "\t\t" << cs_pp + He_abundance*cs_pHe << "\t\t" << cs_Hep << "\t\t" << cs_HeHe << "\t\t" << cs_Hep + He_abundance*cs_HeHe << "\n"; 
	    }
	    //outfileel << "\n\n\n";
	  }

	  //outfileel.close();
	  
	  
	  std::vector <std::vector <double> > vec_ppos, vec_Hepos;
	  for (int e = 0; e < Eprim.size(); e++){
	    std::vector <double> vec_ppos_, vec_Hepos_;
	    for (int u = 0; u < Esec.size(); u++){
	      inx = e*Esec.size() + u;  //This is due to the order of the Fluka table!!
	      vec_ppos_.push_back(Matrix_Pos_pp[inx] + He_abundance*Matrix_Pos_pHe[inx]);
	      vec_Hepos_.push_back(Matrix_Pos_Hep[inx] + He_abundance*Matrix_Pos_HeHe[inx]);
	    }	    
	    vec_ppos.push_back(vec_ppos_);
	    vec_Hepos.push_back(vec_Hepos_);     
	  }
	  
	  std::vector < std::vector <double> > vec_ppos_int, vec_Hepos_int;
	  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
	    std::vector <double> vec_ppos_int_, vec_Hepos_int_;
	    for (int u = 0; u < Esec.size(); u++){
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[e] > Eprim.back()){  //Extrapolation
		je = Eprim.size()-1; 	  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
		if (DRAGONEnergyVector[e] > 7.*Eprim.back())
		  t = (7.*Eprim.back()  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      else if ( DRAGONEnergyVector[e] < Eprim.front()){
		je = 1; jee = 0;
		t =  (DRAGONEnergyVector[e] - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
		if (Eprim.front()/7. > DRAGONEnergyVector[e])
		  t = (Eprim[0]/7. - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      
	      else{ //Interpolation
		while (DRAGONEnergyVector[e] > Eprim[je])
		  je++;
		
		if ((Eprim[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim[je-1])) jee = je-1; 
		else
		  jee = je;
		t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	      }
	      
	      vec_ppos_int_.push_back(std::max(vec_ppos[jee][u] + t*(vec_ppos[je][u] - vec_ppos[je-1][u]), 0.));
	      vec_Hepos_int_.push_back(std::max(vec_Hepos[jee][u] + t*(vec_Hepos[je][u] - vec_Hepos[je-1][u]), 0.));	
	    }
	    
	    vec_ppos_int.push_back(vec_ppos_int_);
	    vec_Hepos_int.push_back(vec_Hepos_int_);
	  }
	  
	  for (int e = 0; e < dimEn; e++){
	    
	    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
	      
	      int je = 0, jee;
	      double t;
	      if (DRAGONEnergyVector[u] > Esec.back()){  //Extrapolation
		je = Esec.size()-1; 	  jee = je;
		t =  (DRAGONEnergyVector[u]  - Esec[jee])/(Esec[je]-Esec[je-1]);
		if (DRAGONEnergyVector[u] > 7.*Esec.back())
		  t = (7.*Esec.back()  - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
	      
	      else if ( DRAGONEnergyVector[u] < Esec.front()){
		je = 1; jee = 0;
		t = (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
		if (Esec.front()/7. > DRAGONEnergyVector[u])
		  t = (Esec[0]/7. - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
	      
	      else{ //Interpolation
		while (DRAGONEnergyVector[u] > Esec[je])
		  je++;
		
		if ((Esec[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec[je-1])) jee = je-1; 
		else
		  jee = je;	
		
		t =  (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
	      }
	      
	      double ppos = std::max(vec_ppos_int[e][jee] + t*(vec_ppos_int[e][je] - vec_ppos_int[e][je-1]), 0.);
	      double Hepos = std::max(vec_Hepos_int[e][jee] + t*(vec_Hepos_int[e][je] - vec_Hepos_int[e][je-1]), 0.);
	      
	      spall_apel[couplepos][u][e] = DRAGONEnergyVector[e]*factorelpos * ppos;
	      spall_apel[coupleposHe][u][e] = 4.0*DRAGONEnergyVector[e]*factorelpos * Hepos;
	      
	      //outfilepos << Epr << "\t\t" << Eel << "\t\t" << cs_pp << "\t\t" << cs_pHe << "\t\t" << cs_pp + He_abundance*cs_pHe << "\t\t" << cs_Hep << "\t\t" << cs_HeHe << "\t\t" << cs_Hep + He_abundance*cs_HeHe << "\n";
	    }
	    
	    //outfilepos << "\n\n\n";
	    
	    //spall_apel[coupleel].push_back(elp_vec);
	    //spall_apel[coupleelHe].push_back(elHe_vec);
	    
	    //spall_apel[couplepos].push_back(posp_vec);
	    //spall_apel[coupleposHe].push_back(posHe_vec);	
	  }
	  //outfilepos.close();	  
	  
	}//End table reading 
	std::cout << "\n Leptons OK!! \n " << std::endl;
	
      }//el  
      
      if (in->prop_Adeuteron) {
	std::cout << "\nDedicated antideuteron propagation..." << std::endl;
	std::pair<int,int> coupleadad(-998,-998);  // Tertiary antideuterons
	std::pair<int,int> coupleadap(-999,-998);  // Secondary antideuterons, from antiprotons
	std::pair<int,int> coupleadp(1001,-998);  // Secondary antideuterons, from protons
	std::pair<int,int> coupleadHe(2004,-998);  // Secondary antideuterons, from Helium
	spall_apel[coupleadp] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleadap] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[coupleadHe] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	
	if (in->any == UserAN){ 
	  std::cout << "Filling sec antideuterons with User-defined table... " << std::endl;
	  Adtable_pp = in->Adtable_pp;
	  Adtable_pHe = in->Adtable_pHe;
	  Adtable_app = in->Adtable_app;
	  Adtable_apHe = in->Adtable_apHe;
	  Adtable_Hep =  in->Adtable_Hep;
	  Adtable_HeHe = in->Adtable_HeHe;
	}
	
	else {
	  if (in->any == !DRAGONAN) std::cout << "Wrong option for sec antideuterons... Taking DRAGON Anti-D cross sections!" << std::endl;
	  std::cout << "Reading sec anti-Deuteron tables method DRAGON... " << std::endl;
	  Adtable_pp = "data/DRAGON_dbar_pp.dat";
	  Adtable_pHe = "data/DRAGON_dbar_pHe.dat"; 
	  Adtable_app = "null";
	  Adtable_apHe = "null";
	  Adtable_Hep = "data/DRAGON_dbar_Hep.dat";
	  Adtable_HeHe = "data/DRAGON_dbar_HeHe.dat";
	}

	if (in->Scale_pAd){
	  std::cout << "scaling pp->aD table to approximate cross sections of the rest of interactions\n" << std::endl;
	  Adtable_pHe =  Adtable_app = Adtable_apHe = Adtable_Hep = Adtable_HeHe = Adtable_pp;
	}
	
	
	TSpallationNetwork::InitXSecDRAGONAN(factorelpos, default_Galpropobject, Adtable_pp, Adtable_pHe, Adtable_app, Adtable_apHe, Adtable_Hep, Adtable_HeHe, in->AN_DRAGON_IneXSec, in->Scale_pAd, "D");
	//TSpallationNetwork::InitXSecDRAGONAN(coord, in, factorelpos, Adtable_pp, Adtable_pHe, Adtable_app, Adtable_apHe, Adtable_Hep, Adtable_HeHe);
	std::cout << "DRAGON secondary anti-Deuterons OK!! \n" << std::endl;
	
	std::cout << "Filling tert anti-Deuterons table DRAGON procedure ... " << std::endl;

	spall[coupleadad] = vector<double>(dimEn, 0);
	for(unsigned int ip = 0; ip < dimEn; ++ip) {
	  PP_inel = 0.0;
	  PA_inel = 0.0;
	  aPP_non = 0.0;
	  aPA_non = 0.0;
	  aPP_ann = 0.0;
	  aPA_ann = 0.0;
	  default_Galpropobject.nucleon_cs(2, DRAGONEnergyVector[ip], -1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);
	  spall[coupleadad][ip] = 2*factorprot*( aPP_non + He_abundance*aPA_non );  //Inelastic ad+p = 2*ap+p 
	  //std::cout << "\n These are the tertiary components!: " << 2*aPP_non << " " << 2*aPA_non << " " << spall[coupleadad][ip] << std::endl;
	}

	std::cout << "Tertiary Anti-D OK!!\n" << std::endl;	  
      } //End of if prop_Adeuteron

      
      if (in->prop_AHe3) {
	std::cout << "\nDedicated anti-He3 propagation..." << std::endl;
	std::pair<int,int> couplea3a3(-997,-997);  // Tertiary antiHe3
	std::pair<int,int> couplea3ap(-999,-997);  // Secondary antiHe3, from antiprotons
	std::pair<int,int> couplea3p(1001,-997);  // Secondary antiHe3, from protons
	std::pair<int,int> couplea3He(2004,-997);  // Secondary antiHe3, from Helium
	spall_apel[couplea3p] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[couplea3ap] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	spall_apel[couplea3He] = vector< vector<double> >(dimEn, vector<double>(dimEn, 0.0));
	
	if (in->any == UserAN){ 
	  std::cout << "Filling sec anti-He3 with User-defined table... " << std::endl;
	  AHe3table_pp = in->AHe3table_pp;
	  AHe3table_pHe = in->AHe3table_pHe;
	  AHe3table_app = in->AHe3table_app;
	  AHe3table_apHe = in->AHe3table_apHe;
	  AHe3table_Hep =  in->AHe3table_Hep;
	  AHe3table_HeHe = in->AHe3table_HeHe;
	}
	
	else {
	  if (in->any == !DRAGONAN) std::cout << "Wrong option for sec anti-He3... Taking DRAGON Anti-He3 cross sections!" << std::endl;
	  std::cout << "Reading sec anti-He3 tables method DRAGON... " << std::endl;
	  AHe3table_pp = "data/DRAGON_hebar_pp.dat";
	  AHe3table_pHe = "data/DRAGON_hebar_pHe.dat"; 
	  AHe3table_app = "null";
	  AHe3table_apHe = "null";
	  AHe3table_Hep = "data/DRAGON_hebar_Hep.dat";
	  AHe3table_HeHe = "data/DRAGON_hebar_HeHe.dat";
	}
	
	if (in->Scale_pAHe3){
	  std::cout << "scaling pp->aHe3 table to approximate cross sections of the rest of interactions\n" << std::endl;
	  AHe3table_pHe =  AHe3table_app = AHe3table_apHe = AHe3table_Hep = AHe3table_HeHe = AHe3table_pp;
	}
	
	
	TSpallationNetwork::InitXSecDRAGONAN(factorelpos, default_Galpropobject, AHe3table_pp, AHe3table_pHe, AHe3table_app, AHe3table_apHe, AHe3table_Hep, AHe3table_HeHe, in->AN_DRAGON_IneXSec, in->Scale_pAHe3, "3He");
	std::cout << "DRAGON secondary anti-3He OK!! \n" << std::endl;
	
	std::cout << "Filling tert anti-3He table DRAGON procedure ... " << std::endl;
	
	spall[couplea3a3] = vector<double>(dimEn, 0);
	for(unsigned int ip = 0; ip < dimEn; ++ip) {
	  PP_inel = 0.0;
	  PA_inel = 0.0;
	  aPP_non = 0.0;
	  aPA_non = 0.0;
	  aPP_ann = 0.0;
	  aPA_ann = 0.0;

	  default_Galpropobject.nucleon_cs(2, DRAGONEnergyVector[ip], -1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);
	  spall[couplea3a3][ip] = 3*factorprot*( aPP_non + He_abundance*aPA_non );  // Approximation is that aHe3 + p = 3*ap + p
	  //std::cout << "\n These are the tertiary components!: " << 3*aPP_non << " " << 3*aPA_non << " " << spall[couplea3a3][ip] << std::endl;
	}

	std::cout << "Tertiary Anti-3He OK!!\n" << std::endl;	  
      } //End of if prop_AHe3

      
    } //End of triple option

    
    if (in->spallationxsec == Fluka || in->spallationxsec == DRAGON2){
      std::pair <double, double> xsecpair;
      //double xsecHint, xsecHeint;
      
      for(std::map<std::pair<int, int>, std::vector<std::pair<double, double> > >::iterator it = spallationXsections.begin(); it != spallationXsections.end(); ++it) {
	std::pair<int, int> couple1 = it->first; // (*it).first	
	if (couple1.second <= 1001) continue;  
	//std::cout <<  "Getting spallation cross sections from " << couple1.first << " to " << couple1.second << std::endl;

	for (int ie = 0; ie < DRAGONEnergyVector.size(); ie++) {
	  xsecpair = TSpallationNetwork::GetXSecPair(couple1.first, couple1.second, DRAGONEnergyVector[ie]);
	  //std::cout << "XSecH " << xsecpair.first <<  " XSecHe " << xsecpair.second << std::endl;
	  
	  if (xsecpair.first > 10000. || xsecpair.first < 0.)  
	    std::cout << "\n\n\n\nSomething strange happens with H when filling " << xsecpair.first << " nuclei: " << couple1.first << "  " << couple1.second << "  Ener: " << DRAGONEnergyVector[ie] << "  " << ie;
	  //if (xsecpair.second != 0. || xsecpair.second < 0.)
  	  if (xsecpair.second > 10000. || xsecpair.second < 0.)
	    std::cout << "\n\nSomething strange happens with He when filling";
	  
	  spallationXsectionsInterpolated[couple1].push_back(xsecpair);

	  if (in->spallationxsec == DRAGON2)
	    xsecpair.second = 0;
	  
	  spall[couple1].push_back(factor*beta[ie]*(spallationXsectionsInterpolated[couple1][ie].first + He_abundance*spallationXsectionsInterpolated[couple1][ie].second));
	  
	  //std::cout << couple1.first << " Added " << couple1.second << "   " << spall[couple1][ie] << std::endl;

	}  //End of Energy loop
	
      } // End of nuclei loop

    } //option Fluka and Carmelo


    else{ //If we have Webber or Galprop we do not have anything to read!
      if (in->spallationxsec == Webber03) std::cout << "\nWebber spallation XSecs\n" << std::endl;
      
      else if (in->spallationxsec == GalpropXSec){ default_Galpropobject.InitXSec();   std::cout << "\nGalprop spallation XSecs\n" << std::endl; }
      
      for (int iloop = 0; iloop < nuclei.size(); ++iloop) {
	
	// If antiprotons or leptons, do not compute spallation. If protons, do not compute, because their spallation products will be computed elsewhere
	if (nuclei[iloop] <= 1001) continue;
        
	int iz = -1000;  // parent nucleus charge
	int ia = -1000;  // parent nucleus mass
	Utility::id_nuc(nuclei[iloop], ia, iz);
	
	for (int idaught = iloop+1; idaught < nuclei.size(); ++idaught) {
	  
	  int jz = -1000; // daughter nucleus charge
	  int ja = -1000; // daughter nucleus mass
	  Utility::id_nuc(nuclei[idaught], ja, jz);
	  
	  std::pair <double, double> couple1 (nuclei[iloop], nuclei[idaught]);
	  //std::cout <<  "Getting spallations " << couple1.first << " to " << couple1.second << std::endl;
	  
	  // If nucleus is antiproton of leptons skip it. Skip also if mass(daughter) > mass(parent)
	  if (nuclei[idaught] < 1001 || ia < ja) continue;
          
	  spall[couple1] = default_Galpropobject.GetXSec(iz, ia, jz, ja);

	  //for (int g; g< dimEn; g++)
	    //std::cout << "\n\n" <<  spall[couple1][g] << " " << couple1.first << " " << couple1.second << std::endl;
	  //if (jz < 3){
	  //for(unsigned int a = 0; a < dimEn; a++) 
	  //    std::cout << spall[couple1][a] << " Added " << iz << ia << "   " << jz << ja << std::endl;
	  //}
	}
      }
    } //End Webber and Galprop option

    std::cout << "\n...General grid inter and extrapolation done\n\n" << std::endl << std::endl;
    //TSpallationNetwork::Retrieve_XSecs(in, nuclei, coord);
    //exit(1);
    
    return ;
  }
    
  
/*
std::vector<double> TGalpropXSec::GetXSec(int iz, int ia, int jz, int ja, std::vector<std::pair<double, double> > csmb_map, std::vector <double> spa_energy) {
    int IZ1,IA1,IZ3,IA3, K_electron = 0;
    int galdef_network_par=0; // temporary solution; value 1 caused problems in nuc_package
    double branching_ratio,t_half;

    const int diz=3;                // maximum delta Z considered
    const double factor = Clight*1.e-27;
   
    std::vector<double> spalla = std::vector<double>(spa_energy.size(), 0.0);

    // a loop over an intermediate state; final state must be as requested
    for (IZ1=(jz-diz>1) ? jz-diz: 1; IZ1<=iz && IZ1<=jz+diz; IZ1++) {
      for (IA1=(2*IZ1-4>ja) ? 2*IZ1-4: ja; IA1<ia && IA1<=2.5*IZ1+4.2; IA1++) {
	// channel selection procedure	

	if (IA1 < IZ1 || ia-IA1 < iz-IZ1) continue;

	t_half = t_half / year;
	
	// skip if long-lived intermediate state
	if(t_half>=t_half_limit                            // IMOS20010816
	   && 100*IZ1+IA1!=100*jz+ja && 100*IZ3+IA3!=100*jz+ja) continue;
	
	for(unsigned int ip = 0; ip < csmb_map.size(); ip++){
	  //if (iz == 14 & ia == 32 & jz == 14 & ja == 30){
	  //std::cout <<  iz << "0" << ia << "   " << jz << "0" << ja << " " <<  csmb_map[ip].first << "  " << ip << std::endl;
	  if (spalla[ip] == 0)
	    spalla[ip] += (csmb_map[ip].first*(1.0 + He_abundance*TGalpropXSec::He_to_H_CS_ratio(spa_energy[ip],iz,ia,IZ1,IA1)));
	  //std::cout << " " << spalla[ip] << "  " << spa_energy[ip] << " " << He_abundance*TGalpropXSec::He_to_H_CS_ratio(spa_energy[ip],iz,ia,IZ1,IA1) << "\n" << std::endl;
	  //}
	  
  	}		
      }   
    } //End nuclei loops
    
    return spalla;
}
*/
  
  std::vector<double> TGalpropXSec::GetXSec(int iz, int ia, int jz, int ja) {
    int IZ1,IA1,IZ3,IA3,kopt,info, K_electron =0;
    int galdef_network_par=0; // temporary solution; value 1 caused problems in nuc_package
    double branching_ratio,t_half;
    
    const int diz=3;                // maximum delta Z considered
    const double factor = Clight*1.e-27;
    
    kopt = cross_section_option; //= 12, defined in constants.h

    if (inn->spallationxsec == Webber03 || inn->spallationxsec == DRAGON2){  //In case we ask for Webber
      kopt = 1;
      //std::cout << "Here! kopt = " << kopt << "for jz ja: " << jz << " " << ja << std::endl;
      if (jz == 3)
	kopt = 2;
    }

    std::vector<double> spalla = std::vector<double>(energy.size(), 0.0);
    
    // a loop over an intermediate state; final state must be as requested
    for (IZ1=(jz-diz>1) ? jz-diz: 1; IZ1<=iz && IZ1<=jz+diz; IZ1++) {
      for (IA1=(2*IZ1-4>ja) ? 2*IZ1-4: ja; IA1<ia && IA1<=2.5*IZ1+4.2; IA1++) {

	// channel selection procedure
	if (IA1 < IZ1 || ia-IA1 < iz-IZ1) continue;
	
	// IMOS20010816 line below
	branching_ratio = TGalpropXSec::nucdata(galdef_network_par,IZ1,IA1,K_electron,jz,ja, &IZ3,&IA3,&t_half);

	if (branching_ratio == 0) continue;	
	t_half = t_half / year;
	
	// skip if long-lived intermediate state
	if(t_half>=t_half_limit                            // IMOS20010816
	   && 100*IZ1+IA1!=100*jz+ja && 100*IZ3+IA3!=100*jz+ja) continue;
	
	for(unsigned int ip = 0; ip < energy.size(); ip++){ 

	  //if (iz == 14 & ia == 32 & jz == 14 & ja == 30){
	  //std::cout << "\n\n" <<  iz << "0" << ia << "   " << jz << "0" << ja << "\n" << TGalpropXSec::isotope_cs(energy[ip]*1000.0,iz,ia,IZ1,IA1,kopt,&info) << std::endl;
	    spalla[ip] += (TGalpropXSec::isotope_cs(energy[ip]*1000.0,iz,ia,IZ1,IA1,kopt,&info)*branching_ratio*(1.0+ He_abundance*TGalpropXSec::He_to_H_CS_ratio(energy[ip],iz,ia,IZ1,IA1))*beta[ip]*factor);	    
	    //std::cout << "\n" << spalla[ip]/(factor*beta[ip]) << std::endl;
	    //}
	    //if (jz > 0)
	    //std::cout << "\n\n" <<  iz << "0" << ia << "   " << jz << "0" << ja << "\n" << TGalpropXSec::isotope_cs(energy[ip]*1000.0,iz,ia,IZ1,IA1,kopt,&info) << std::endl;
	}
	
      } //ja1 //jz1
    }  
    
    return spalla;
  }
  

  
  void TGalpropXSec::InitXSec(){
    
    int i,j,k,size;
    const int BufferSize=200;
    char readBuffer[BufferSize];
    
    std::ifstream data;

    data_filename.push_back("data/galprop_nucdata.dat");    
    data_filename.push_back("data/galprop_p_cs_fits.dat");
    data_filename.push_back("data/galprop_eval_iso_cs_updated.dat"); 
    
    for (int j = 0; j < 3 ; j++){
      std::cout << "Opening file " << data_filename[j].c_str() << std::endl;
      data.open(data_filename[j].c_str());                    // open file if exists
      if(data.fail())    {std::cout <<"Error opening file "<< data_filename[j] << std::endl;
	exit(1);  }
      
      while(!isspace(data.get()) && !data.eof())      // skip comments:
	data.getline(readBuffer,BufferSize,'\n');    // any symbol in 1st col. 
      
      for(i=0; i<3; data >> n_data[i++][j]);          // read array's sizes
      data.getline(readBuffer,BufferSize,'\n');       // skip the rest of line
      
      for(size=1, i=0; i<3; size*=n_data[i++][j]);    // allocate space
      
      if (j == 0) nucdat = new float[size];
      else if (j == 1) protdat = new float[size];
      else csdat = new float[size];
      
      for(k = 0; k < size && !data.eof();)            // read data loop
	{
	  while(!isspace(data.get()) && !data.eof())   // skip comments:
	    data.getline(readBuffer,BufferSize,'\n'); // any symbol in 1st col.
	  if (j == 0)
	    for(i=0; i < n_data[0][j]; i++) data >> *(nucdat+k++);
	  else if (j == 1)
	    for(i=0; i < n_data[0][j]; i++) data >> *(protdat+k++);
	  else
	    for(i=0; i < n_data[0][j]; i++) data >> *(csdat+k++);
	  
	  data.getline(readBuffer,BufferSize,'\n');    // skip the rest of line
	}
      
      data.close();
    }
    
  }
  

  double TGalpropXSec::isotope_cs(double emev,int iz,int ia,int izf,int iaf,int kopt,int* info) {
    int a1,a2,i,j,size, itable=0, info1;
    double e1,y,err2,xi1,xi2, f1=0., f2=0., T[11], a[3]={1.,1.,0.}, b[6];
    float *cs_data = csdat, *p_cs = protdat;
    double *tp=T;
    double ej;
    double CSmb = 0.0;
    
    e1 = emev;
    *info = kopt;

    // CHECK if user wants to use specific program (the value of "kopt")
    if(kopt == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);       // Webber's code IMOS20020502
    if(kopt == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);         // TS code       IMOS20020502
    CSmb = max(0.,CSmb);
    if(kopt == 1 || kopt == 2) return(CSmb);
    
    a1 = fnuc(iz, ia);
    a2 = fnuc(izf,iaf);
    
    // EVALUATED CROSS SECTIONS
    
    if(kopt == 12 || kopt == 22)
      {
	CSmb = eval_cs(emev,a1,a2,&info1);
	if (info1 > 0) return(std::max(0.,CSmb));
	kopt--;                                          // if evaluation doesn't exist,
      }                                                   // try other options
    
    // if user wants, use THE CROSS SECTION FITS

    if(kopt == 11 || kopt == 21)
      {      
	// special cases: Be, B - recursion calls
	if(izf != 0)
	  {
	    // A = 10
	    if(10 == iaf)  // B10 = B10 + C10 = a10 - Be10
	      {
		b[0] = isotope_cs(emev,iz,ia,0,iaf,21,&j);
		if(j == -21)
		  {                                          // B10
		    if(510 == a2) 
		      {
			b[0]-=isotope_cs(emev,iz,ia,4,iaf,21,&j);
			return(std::max(0.,b[0]));
		      }
		    if(5 < izf) return(0.);                // C10, =0
		  }
	      }
	    // A = 11
	    if(11 == iaf)  // B11 = a11 = Be11 + B11 + C11
	      {
		b[0] = isotope_cs(emev,iz,ia,0,iaf,21,&j);
		if(j == -21) 
		  {
		    if(511 == a2) return(std::max(0.,b[0]));     // B11
		    return(0.);                             // =0 for the rest
		  }
	      }
	  }
	// straight search in the table

	for(i=0; i<n_data[1][1]-1; i++, p_cs+=n_data[0][1]) // -1 fixes the reading error in the line below
	  if(a1 == inuc(*p_cs) && a2 == inuc(*(p_cs+1)))

	    {
	      for(p_cs+=2, j=0; j<6; b[j++]=*p_cs++);    // take the parameters
	      if(b[0] >= 0.)                             // if positive use fit
		{
		  *info=-kopt;
		  if(emev < b[5]) return(0);              // fitting function
		  b[0]*=(1.+sin(b[1]*pow(log10(emev),1.*b[2]))*exp(-b[3]*(emev-b[4])));

		  return(std::max(0.,b[0]));
		}
	      kopt = (int)(-b[0]+0.1);                   // negative b[0] gives kopt
	    }
	if(izf == 0) return(0.);
      }


    // CHECK if user wants to use specific program (the value of "kopt")
    if(kopt == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);       // Webber's code IMOS20020502
    if(kopt == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);         // TS code       IMOS20020502
    CSmb = max(0.,CSmb);
    if(kopt == 1 || kopt == 2) return(CSmb);
    
    // STARTING THE ALGHORITHM
    
    for(i=0; i<11; T[i++] = 0.);
    
    // CHECK the array: cs_data (is there a channel we are looking for ?)
    
    for(size=1, i=0; i<3; size*=n_data[i++][2]);
    for(tp = T, i=0; i<size; i+=n_data[0][2], tp = T, f1=0., f2=0.)
      {
	if(a1 != inuc(*(cs_data+i)))   continue;
	if(a2 != inuc(*(cs_data+i+1))) continue;

	// if there is such a channel then the LEAST-SQUARE FIT
	
	itable++;
	if(*(cs_data+i+4) < 0.) *(cs_data+i+4) *= -*(cs_data+i+3);  // calc.abs.err.
	err2 = pow(*(cs_data+i+4),2);                               // err^2
	
	y = *(cs_data+i+3);                                         // cs measured
	ej = *(cs_data+i+2);                                        // @ energy
	
	if(kopt/10 != 2) f1= wsigma_cc(iz,ia,izf,iaf,ej);             // Webber IMOS20020502
	if(kopt/10 != 1) f2= yieldx_cc(iz,ia,izf,iaf,*(cs_data+i+2)); // TS     IMOS20020502
	
	// calculations of the separate terms:
	*tp++ += f1*y /err2;       // Webber
	*tp++ += f1*f1/err2;
	*tp++ += f2*y /err2;       // TS
	*tp++ += f2*f2/err2;
	*tp++ += y    /err2;       // const cs
	*tp++ += 1.   /err2;
	
	// calculation of terms for the Xi2 estimates
	*tp++ += y*y    /err2;
	*tp++ += 2.*f1*y/err2;
	*tp++ += f1*f1/err2;
	*tp++ += 2.*f2*y/err2;
	*tp   += f2*f2/err2;
	
	// calculation of renormalization coefficients 
	for(j=0; j<3; j++) {
	  a[j]= (T[2*j+1] != 0.) ? T[2*j]/T[2*j+1]: a[j];
	}
      }
    
    if(kopt == 3 && a[2] != 0.) return(a[2]);                  // const cr.sect.
    if(kopt/10 == 1) CSmb = wsigma_cc(iz,ia,izf,iaf,emev);     // Webber code IMOS20020502
    if(kopt/10 == 2) CSmb = yieldx_cc(iz,ia,izf,iaf,e1);       // TS code     IMOS20020502
    if(kopt/10 == 1 || kopt/10 == 2) return(std::max(0.,CSmb*a[kopt/10-1]));

    // CHOOSE THE BEST APPROXIMATION (kopt = 0)
    if(itable < 2)                                         // no data or 1 pt.
      {	    
	*info = itable;
	CSmb = a[0]*wsigma_cc(iz,ia,izf,iaf,emev);          // use Webber code     IMOS20020502
	if(CSmb <= 0.) 
	  {                                                   // if W-code give 0,
	    CSmb = a[1]*yieldx_cc(iz,ia,izf,iaf,e1);         // take the TS approx. IMOS20020502
	    if(CSmb != 0. && itable == 1) *info = 2;
	  }
      } else                                                 // data exists
      {

	xi1= T[6] -a[0]*T[7] +a[0]*a[0]*T[8];               // Xi2 evaluation 1
	xi2= T[6] -a[1]*T[9] +a[1]*a[1]*T[10];              // Xi2 evaluation 2
	if(xi1 < xi2)
	  {
	    *info = 1;
	    CSmb = a[0]*wsigma_cc(iz,ia,izf,iaf,emev);       // renorm. Webber approx. IMOS20020502
	  } else
	  {
	    *info = 2;
	    CSmb = a[1]*yieldx_cc(iz,ia,izf,iaf,e1);         // renorm. TS approx.     IMOS20020502
	  }
      }
    return(std::max(0.,CSmb));
  }
  
  
  
  double TGalpropXSec::eval_cs(double emev,int za1,int za2,int* info) {
    int i,size;
    float *eval = csdat;
    double x[2]={-1.e10,1.e10},y[2]={0.,0.};
    
    // CHECK the array: eval (is there a channel we are looking for ?)
	
    for(size=1, i=0; i<3; size*=n_data[i++][2]);
      for(*info=0, i=0; i<size; i+=n_data[0][2])
      {
	if(za1 != inuc(*(eval+i)))   continue;
	if(za2 != inuc(*(eval+i+1))) continue;
	if(x[0] < *(eval+i+2) && *(eval+i+2) <= emev)   // find lower energy pt
	  {		    
	    x[0] = *(eval+i+2); 
	    y[0] = *(eval+i+3);		
	  }
	if(emev <= *(eval+i+2) && *(eval+i+2) < x[1])   // find higher energy pt 
	  {

	    x[1] = *(eval+i+2); 
	    y[1] = *(eval+i+3);

	  }
      }

    if(x[0]*x[1] < -1.e19) { *info = -1; return(0.);} // no evaluation found, return 0
    
    if(x[0] <   0.) { *info = 1; return(0.); }         // no lower grid pt, return 0
    if(x[1] > 9.e9) { *info = 2; return(y[0]); }       // no higher grid pt, extrapolate
    
    if(x[1]-x[0] == 0.) { *info = 3; return(y[1]); }   // emev falls exactly on the grid

    for(*info = 4, i=0; i<2; i++){ x[i] = log10(x[i]);}
    return(y[0]+(log10(emev)-x[0])*(y[1]-y[0])/(x[1]-x[0]));// interpolate
  }
  
  
  
  
  
  
  double TGalpropXSec::nucdata(int ksp,int iz,int ia,int K_electron,int izf,int iaf,int* izl,int* ial,double* To) {
    int i,j,k,l,m,n,iy,iz0,ia0,iz4,ia4,iz5,ia5,iw[121],
      nksp=ksp*n_data[0][0]*n_data[1][0];
    
    float w[2][121], *decay = nucdat+nksp;
    double b, xxx;
    
    // STABLE & LONG-LIVED ISOTOPES (numbers in the table are the proton numbers)
    // The long-lived isotopes listed in "longliv" table are included as stable;
    int stable[64][3] = {      // second index changes faster
      1,  0,  1,     1,  0,  1,     1,  0,  2,     2,  0,  2, //  A = 1- 4
      0,  0,  0,     3,  0,  3,     3,  0,  4,     0,  0,  0, //  A = 5- 8
      4,  0,  4,     4,  0,  5,     5,  0,  5,     6,  0,  6, //  A = 9-12
      6,  0,  6,     6,  0,  7,     7,  0,  7,     8,  0,  8, //  A =13-16
      8,  0,  8,     8,  0,  8,     9,  0,  9,    10,  0, 10, //  A =17-20
      10,  0, 10,    10,  0, 11,    11,  0, 11,    12,  0, 12, //  A =21-24
      12,  0, 12,    12,  0, 13,    13,  0, 13,    14,  0, 14, //  A =25-28
      14,  0, 14,    14,  0, 14,    15,  0, 15,    14,  0, 16, //  A =29-32
      16,  0, 16,    16,  0, 16,    17,  0, 17,    16, 17, 18, //  A =33-36
      17,  0, 18,    18,  0, 18,    18,  0, 19,    18, 19, 20, //  A =37-40
      19,  0, 20,    18,  0, 20,    20,  0, 20,    20,  0, 22, //  A =41-44
      21,  0, 21,    20,  0, 22,    22,  0, 22,    20,  0, 22, //  A =45-48
      22,  0, 23,    22, 23, 24,    23,  0, 24,    24,  0, 24, //  A =49-52
      24,  0, 25,    24, 25, 26,    25,  0, 26,    26,  0, 28, //  A =53-56
      26,  0, 27,    26,  0, 28,    27,  0, 28,    26, 27, 28, //  A =57-60
      28,  0, 28,    28,  0, 28,    28,  0, 29,    28,  0, 28  //  A =61-64
    };
    
    // LONG-LIVED ISOTOPES (>~1y): Zi.Ai T_1/2(y) Zf.Af - [][][0] half-life shown for naked nucleus
    int nll = 25;                                 // - [][][1] half-life shown for H2-like atoms
    float longliv[25][2][3] = {   // third index changes faster
      1.03,  12.33,    2.03,     //  3H (b-) 3He   100% [ToI]
      1.03,  12.33,    2.03,     // no EC
      
      4.07,  0.,       4.07,     // stable
      4.07,  0.1459,   3.07,     //  7Be(EC) 7Li   100% [ToI]
      
      4.10,  1.60e6,   5.10,     // 10Be(b-)10B    100% [ToI]
      4.10,  1.60e6,   5.10,     // no EC
      
      6.14,  5.73e3,   7.14,     // 14C (b-)14N    100% [ToI]
      6.14,  5.73e3,   7.14,     // no EC
      
      11.22,  4.80e3,  10.22,     // 22Na(b+)22Ne        [M98]
      11.22,  2.60e0,  10.22,     // 22Na(EC?)22Ne       [ToI] T1/2(Lab)=2.60e0 y
      
      13.26,  9.10e5,  12.26,     // 26Al(b+)26Mg        [M98]
      13.26,  4.075e6, 12.26,     // 26Al(EC)26Mg        [M98] T1/2(Lab)=7.4e5 y [ToI]
      
      14.32,  172.,    16.32,     // 32Si(2b-)32S   100% [ToI] Si-P -S 
      14.32,  172.,    16.32,     // no EC
      
      17.36,  3.07e5,  18.36,     // 36Cl(b-)36Ar        [ToI]
      17.36,  1.58e7,  16.36,     // 36Cl(EC)36S         [ToI] T1/2(Lab)=3.01e5 y
      
      18.37,  0.,      18.37,     // stable
      18.37,  0.1,     17.37,     // 37Ar(EC)37Cl   100% [ToI] T1/2(Lab)=35.04 d
      
      18.39,  2.69e2,  19.39,     // 39Ar(b-)39K    100% [ToI]
      18.39,  2.69e2,  19.39,     // no EC
      
      19.40,  1.43e9,  20.40,     // 40K (b-)40Ca   89.3%[ToI] T1/2(Lab)=1.277e9 y incl 10.7% ECb+
      19.40,  1.43e9,  20.40,     // no EC
      
      20.41,  0.,      20.41,     // stable
      20.41,  1.03e5,  19.41,     // 41Ca(EC)41K    100% [ToI]
      
      18.42,  32.9,    20.42,     // 42Ar(2b-)42Ca  100% [ToI] Ar-K -Ca
      18.42,  32.9,    20.42,     // no EC
      
      22.44,  0.,      22.44,     // stable
      22.44,  49.,     20.44,     // 44Ti(ECb+)44Ca 100% [ToI] Ti(EC)Sc(b+)Ca
      
      23.49,  0.,      23.49,     // stable
      23.49,  0.903,   22.49,     // 49V (EC)49Ti   100% [ToI] 
      
      24.51,  0.,      24.51,     // stable
      24.51,  0.076,   23.51,     // 51Cr(EC)51V   <100% [ToI] 
      
      25.53,  0.,      25.53,     // stable
      25.53,  3.74e6,  24.53,     // 53Mn(EC)53Cr   100% [ToI]
      
      25.54,  6.30e5,  26.54,     // 54Mn(b-)54Fe        [W98]
      25.54,  0.855,   24.54,     // 54Mn(EC)54Cr        [ToI] T1/2(Lab)=312.3 d
      
      26.55,  0.,      26.55,     // stable
      26.55,  2.73e0,  25.55,     // 55Fe(EC)55Mn   100% [ToI]
      
      28.56,  4.00e4,  26.56,     // 56Ni(2b+)56Fe <100% [F99] Ni-Co-Fe
      28.56,  0.1,     26.56,     // 56Ni(ECb+)56Fe      [ToI] T1/2(Lab)=~30 d Ni(EC)Co(b+)Fe
      
      27.57,  0.,      27.57,     // stable
      27.57,  0.744,   26.57,     // 57Co(EC)57Fe   100% [ToI]
      
      28.59,  0.,      28.59,     // stable
      28.59,  7.60e4,  27.59,     // 59Ni(EC)59Co  <100% [ToI] [B76]
      
      27.60,  5.27e0,  28.60,     // 60Co(b-)60Ni   100% [ToI]
      27.60,  5.27e0,  28.60,     // no EC
      
      26.60,  1.50e6,  27.60,     // 60Fe(b-)60Co   100% [ToI]
      26.60,  1.50e6,  27.60,     // no EC
        
      28.63,  1.00e2,  29.63,     // 63Ni(b-)63Cu   100% [ToI]
      28.63,  1.00e2,  29.63      // no EC
    };
    // K-capture nuclei - factor of 2 because only 1 electron
    for(i=0;i<nll;i++) if(longliv[i][0][1]!=longliv[i][1][1]) longliv[i][1][1] *=2.; 
    
    // BOUNDARY NUCLEI 
    // on the left side from the left boundary, the proton-emission is assumed;
    // on the right side from the right boundary, the neutron-emission is assumed.
    // ZB.AB; left boundary[][0]: Nn=0(1)28; right boundary[][1]: Np=1(1)28 
    int nb = 29;
    float boundary[29][2] = {  // second index changes faster
      1.01,  1.04,    3.04,  2.08,
      3.05,  3.11,    6.09,  4.14,
      6.10,  5.15,    8.13,  6.16,
      8.14,  7.21,   10.17,  8.22,
      12.20,  9.24,   12.21, 10.26,
      14.24, 11.30,   14.25, 12.30,
      15.27, 13.34,   16.29, 14.35,
      16.30, 15.38,   18.33, 16.40,
      20.36, 17.43,   20.37, 18.46,
      20.38, 19.48,   22.41, 20.51,
      22.42, 21.52,   24.45, 22.53,
      24.46, 23.55,   24.47, 24.59,
      25.49, 25.63,   26.51, 26.65,
      27.53, 27.69,   27.54, 28.69,
      28.56, 00.99
    };
    
    b = *To = *izl = *ial = 0;
    if(iz <= 0 || ia <= 0) return(0.);    // check against negative numbers,
    if(iz*ia > 1 && iz >= ia) return(0.); // non-existed nuclei,
    if(2864 < fnuc(iz,ia)) return(0.);    // Ni64 is the heaviest nucleus
    if(64 < ia) return(0.);               // A=64 is the maximal atomic number
    
    // CHECK FOR NUCLEI OUTSIDE THE BOUNDARIES (p/n decay)
    iz0 = iz;
    ia0 = ia;
    if(ia>inuc(boundary[iz-1][1]-iz)) ia0=inuc(boundary[iz-1][1]-iz); // n -decay
    if(29>ia-iz) if(ia>inuc(modf(boundary[ia-iz][0], &xxx)))          // p -decay
		   { 
		     iz0=(int)boundary[ia-iz][0];
		     ia0=inuc(boundary[ia-iz][0]-iz0);
		   }
    
    for(i=0; i<121; iw[i++]=-1)  for(j=0; j<2; w[j++][i]=0.);
    
    // SEARCH FOR A SPECIAL CASE (non beta decay)
    for(i=0; i<n_data[1][0]; i++)
      if(fnuc(iz0,ia0) == inuc(*(decay +i*n_data[0][0])))
	{
	  iw[0]  = i;            // if found, save the line number
	  w[1][0]= 1.00;         // assign 1 to the branching ratio
	  break;
	}
    
    // STANDARD CASE (beta decay & long-lived isotopes)
    if(iw[0] < 0)
      {
	iz5 = iz0;
	ia5 = ia0;
	// *** BETA DECAY ***
	if(iz0 > stable[ia0-1][2]) iz5 = stable[ia0-1][2];   // b+ decay
	if(iz0 < stable[ia0-1][0]) iz5 = stable[ia0-1][0];   // b- decay
	// *** LONG-LIVED ISOTOPES (>~1 y) ***
	for(i=0; i<nll; i++)
	  if(fnuc(iz5,ia5) == inuc(longliv[i][K_electron][0]))
	  {
	    *izl = iz5;
	    *ial = ia5;
	    *To = longliv[i][K_electron][1]*year;
	    if(!*To) *izl=*ial=0;
	    iz5 = (int) longliv[i][K_electron][2];
	    ia5 = inuc(longliv[i][K_electron][2]-iz5);
	    break;
	  }
	if(fnuc(izf,iaf)==fnuc(iz5,ia5) || fnuc(izf,iaf)==fnuc(*izl,*ial)) b = 1.;
	if(fnuc(iz0,ia0) == fnuc(*izl,*ial)) *izl = *ial = 0;
	return(b);
      }
    
    // DEVELOPING A NETWORK OF DECAYS
    for(l=-1, m=0, ia4=1, i=0; i<4; ia4 =(int) pow(3.,++i))
      {
	for(l+=ia4, iy=0, n=0; n<ia4; n++, m++)
	  {                                                      // check if there is
	    if(iw[m] < 0) continue;                             // a required channel
	    for(w[0][m]=0., k=2; k<8; k+=2)
	      {
		w[0][l+3*n+k/2]                                  // store sec.nuclei
		  =*(decay +iw[m]*n_data[0][0]+k-1);
		w[1][l+3*n+k/2]                                  // store branchings
		  =*(decay +iw[m]*n_data[0][0] +k)*w[1][m];
		for(j=0; j<iw[m]; j++)                           // check if sec.nucleus
		  if(*(decay +iw[m]*n_data[0][0] +k-1) // also develops a
		     == *(decay +j*n_data[0][0]))   // network of decays
		    {
		      iw[l+3*n+k/2] = j;                         // store such a nucleus
		      iy = l+3*n+k/2;
		    }  ///printf("%d %d %d %d %d %d\n",l,n,m,k,j,iy);
	      }
	  }
	if(iy == 0) break;
      }
    
    // CHECK FOR STABILITY OF THE FINAL NUCLEI
    for(k=0; k<=l+3*n; k++)
      {
	*To = *izl = *ial = 0;
	if(w[0][k] == 0.) continue;
	iz4 = (int) w[0][k];
	ia4 = inuc(w[0][k]-iz4);
	iz5 = iz4;
	ia5 = ia4;
	// *** BETA DECAY ***
	if(iz4 > stable[ia4-1][2]) iz5 = stable[ia4-1][2];   // b+ decay
	if(iz4 < stable[ia4-1][0]) iz5 = stable[ia4-1][0];   // b- decay
	// *** LONG-LIVED ISOTOPES (>~1 y) ***
	for(i=0; i<nll; i++)
	  {
	    if(fnuc(iz5,ia5) != inuc(longliv[i][K_electron][0])) continue;
	    *izl = iz5;
	    *ial = ia5;
	    *To = longliv[i][K_electron][1]*year;
	    if(!*To) *izl=*ial=0;
	    iz5 = (int) longliv[i][K_electron][2];
	    ia5 = inuc(longliv[i][K_electron][2]-iz5);
	    break;
	  }
	if(fnuc(izf,iaf) == fnuc(*izl,*ial) || fnuc(izf,iaf) == fnuc(iz5,ia5))
	  return(w[1][k]);
      }
    return(b);
  }

  
  
  //TPP
  std::vector<double> TSpallationNetwork::GetXSecTPP(std::vector<double> nu_vector) { // nu_vector -> vector of frequencies of ISRF
    
    //double photon_energy = 10.*1.e-9; //GeV
    double m_e           = 0.511e-3 ; //GeV
    double hPlanck = 4.135667e-15;
    
    std::vector<double> result(DRAGONEnergyVector.size()*DRAGONEnergyVector.size()*nu_vector.size());
    
    for (int k =0; k<DRAGONEnergyVector.size(); ++k) {
      for (int inu=0; inu < nu_vector.size(); ++inu) {
	for (int l = 0; l < DRAGONEnergyVector.size(); ++l) { // l -> energy index of parent electron
	  result[(k*nu_vector.size() + inu)*DRAGONEnergyVector.size() + l] = 0.;
	}
      }
    }
    
    for (int k =0; k<DRAGONEnergyVector.size(); ++k) { // k -> energy index of secondary particle
      
      for (int inu=0; inu < nu_vector.size(); ++inu) {
        
	double photon_energy = hPlanck*nu_vector[inu]*1.e-9; //GeV
	
	if (k == 0)
	  std::cout << "***" << std::endl;
	
	for (int l = 0; l < DRAGONEnergyVector.size(); ++l) { // l -> energy index of parent electron
          
	  int index = (k*nu_vector.size() + inu)*DRAGONEnergyVector.size() + l;
	  
	  double Average_TPP_positron_energy = 0.5 * sqrt(DRAGONEnergyVector[l] / photon_energy ) * m_e;
	  
	  double s = DRAGONEnergyVector[l] * hPlanck*nu_vector[inu]*1.e-9 / pow(m_e, 2.); //tutto in GeV; s: numero puro
	  
	  double correction = 0.;
	  if (s>4.) { 
	    correction = 0.1193662073189215 * exp(-3.839009766244388 - 0.1223186374016053*pow(log(-4. + s), 2) ) * pow((-4. + s),1.8194204144025334);
	    
	    if (s>79.)
	      correction = 0.13 * (-8.07407 + 3.11111 * log(0.0153186 * DRAGONEnergyVector[l] * 3.));
	    
	  }
	  
	  if (k == 0)
	    std::cout << " photon energy [eV] " << photon_energy*1.e9 << " energy of electron [GeV] -> " << DRAGONEnergyVector[l] << " s-> " << s << " correction-> " << correction << std::endl;
	  
          
	  double sigma_tot =  9.46 * 6.65 * 1.e-2 * (1/137.) * correction;// * temp; // c sigma_Thompson alpha_F; units: cm/Myr * cm^2
	  //double std_dev = 30;//Average_TPP_positron_energy;
	  double std_dev = sqrt(Average_TPP_positron_energy); // GeV
	  //double x = DRAGONEnergyVector[l]*photon_energy/(m_e*m_e);
	  //if (x>4)
	  result[index] = sigma_tot * DRAGONEnergyVector[l] * exp( - pow(DRAGONEnergyVector[k] - Average_TPP_positron_energy,2.0) / (2.*(pow(std_dev,2.0))) ) / (sqrt(2*3.14) * std_dev);
	  //			   cm^3/Myr   GeV          1/GeV
	}
      }
    }
    
    return result;
    
  }  

  void TSpallationNetwork::InitXSecPohl(double factorelpos) {
    
    std::pair<int,int> coupleel(1001, -1000); // Electrons from protons
    std::pair<int,int> couplepos(1001, 1000); // Positrons from protons
    std::pair<int,int> coupleelHe(2004, -1000); // Electrons from He
    std::pair<int,int> coupleposHe(2004, 1000); // Positrons from He
    
    InitDataTablesPohl();
    
    // Huang & Pohl tables, arXiv:0711.2528,  http://cherenkov.physics.iastate.edu/gamma-prod/
    
    //ofstream outfile("elposspectr_P.dat", ios::out);
    double Emin_table_lepton = 0.01;
    double Emax_table_lepton = 1.e8;
    
    /* Formula for secondary lepton KINETIC energy bin: 
       
       int i = int((log(Ek)-log(Emin_table_lepton))/(log(Emax_table_lepton)-log(Emin_table_lepton))*201.0 + 1.0);
       // but be careful with the names of indexes. Also below.
       */
    
    const int dimleptHP = 201;
    const int dimlept = dimleptHP-1;
    const double DBlogHP = log(Emax_table_lepton/Emin_table_lepton)/(double)dimleptHP;
    double EleptHP[dimleptHP];
    for (int i=0; i < dimleptHP; i++) EleptHP[i] = exp(log(0.01)+double(i)*DBlogHP);
    
    /* Formula for primary proton (He) TOTAL energy:   */
    const int dimETHP = 374;
    double ET[dimETHP];
    for (int j = 0; j < dimETHP; j++) {
      ET[j] = 1.24*pow(1.05,j)-mp;
      //cout << ET[j] << " " << ProdXsec[j] <<  endl;    
    }

    int Nelectrons = 374;
    double DBlogpr = (log(ET[Nelectrons-1])-log(ET[0]))/(double)Nelectrons;
    //cout << "ProdXsec.size() = " << ProdXsec.size() << " " << ProdXsec.size()/2 << endl;
    //cout << "Matrix.size() = " << Matrix_El_pp.size() << " " << Matrix_El_pp.size()/Nelectrons << endl;
    
    gsl_spline *spline_ProdXsec_pp = gsl_spline_alloc(gsl_interp_cspline, Nelectrons);
    gsl_spline *spline_ProdXsec_He = gsl_spline_alloc(gsl_interp_cspline, Nelectrons);
    
    gsl_spline_init(spline_ProdXsec_pp, ET, &(ProdXsec[0]), Nelectrons);
    gsl_spline_init(spline_ProdXsec_He, ET, &(ProdXsec[Nelectrons]), Nelectrons);
    
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline_El_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    gsl_spline *spline_El_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    gsl_spline *spline_Pos_pp = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    gsl_spline *spline_Pos_Hep = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    
    gsl_spline *spline_El_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    gsl_spline *spline_El_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    gsl_spline *spline_Pos_pp_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    gsl_spline *spline_Pos_Hep_up = gsl_spline_alloc(gsl_interp_cspline, dimlept+1);
    
    double Vectorlept_El_pp[dimlept+1];
    double Vectorlept_El_Hep[dimlept+1];
    double Vectorlept_Pos_pp[dimlept+1];
    double Vectorlept_Pos_Hep[dimlept+1];
    double Vectorlept_El_pp_up[dimlept+1];
    double Vectorlept_El_Hep_up[dimlept+1];
    double Vectorlept_Pos_pp_up[dimlept+1];
    double Vectorlept_Pos_Hep_up[dimlept+1];
    
    for (unsigned int j=0; j<dimEn; j++){
      
      //double Epr = DRAGONEnergyVector[j];
      //double momentum = sqrt(Epr*Epr + 2.0*mp*Epr);
      int j_pr = int(log(Epr/ET[0])/DBlogpr)-1;
      if (j_pr > Nelectrons-2) j_pr = Nelectrons-2;
      double u = (Epr-ET[j_pr])/(ET[j_pr+1]-ET[j_pr]);
      double interp_prodxsec;
      double interp_prodxsecHe;
      
      if (Epr < ET[0]) {
	interp_prodxsec = ProdXsec[0];
	interp_prodxsecHe = ProdXsec[Nelectrons];
      }
      else if (Epr > ET[Nelectrons-1]) {
	interp_prodxsec = ProdXsec[Nelectrons-1];
	interp_prodxsecHe = ProdXsec[2*Nelectrons-1];
      }
      else {
	interp_prodxsec = gsl_spline_eval(spline_ProdXsec_pp, Epr, acc); 
	interp_prodxsecHe = gsl_spline_eval(spline_ProdXsec_He, Epr, acc);
      }
      //cout << Epr << " " << interp_prodxsec << " " << interp_prodxsecHe << endl;
      
      for (int i = 0; i <= dimlept; i++) {
	//int index = index_matrix(i,j_pr);
	Vectorlept_El_pp[i] = ElppPohl[i][j_pr];//[index];
	Vectorlept_El_Hep[i] = ElHepPohl[i][j_pr];//Matrix_El_Hep[index];
	Vectorlept_Pos_pp[i] = PosppPohl[i][j_pr];//Matrix_Pos_pp[index];
	Vectorlept_Pos_Hep[i] = PosHepPohl[i][j_pr];//Matrix_Pos_Hep[index];
	//
	//if (j_pr <= Nelectrons)
	//index = index_matrix(i,j_pr+1);
	//else
	//    index = Nelectrons;
	Vectorlept_El_pp_up[i] = ElppPohl[i][j_pr+1];//[index];
	Vectorlept_El_Hep_up[i] = ElHepPohl[i][j_pr+1];//Matrix_El_Hep[index];
	Vectorlept_Pos_pp_up[i] = PosppPohl[i][j_pr+1];//Matrix_Pos_pp[index];
	Vectorlept_Pos_Hep_up[i] = PosHepPohl[i][j_pr+1];//Matrix_Pos_Hep[index];
      }
      
      gsl_spline_init(spline_El_pp, EleptHP, Vectorlept_El_pp, dimlept+1);
      gsl_spline_init(spline_El_Hep, EleptHP, Vectorlept_El_Hep, dimlept+1);
      gsl_spline_init(spline_Pos_pp, EleptHP, Vectorlept_Pos_pp, dimlept+1);
      gsl_spline_init(spline_Pos_Hep, EleptHP, Vectorlept_Pos_Hep, dimlept+1);
      
      gsl_spline_init(spline_El_pp_up, EleptHP, Vectorlept_El_pp_up, dimlept+1);
      gsl_spline_init(spline_El_Hep_up, EleptHP, Vectorlept_El_Hep_up, dimlept+1);
      gsl_spline_init(spline_Pos_pp_up, EleptHP, Vectorlept_Pos_pp_up, dimlept+1);
      gsl_spline_init(spline_Pos_Hep_up, EleptHP, Vectorlept_Pos_Hep_up, dimlept+1);
      
      for (unsigned int i = 0; i < dimEn; i++) {
	
	double Eel = std::min(1e8, DRAGONEnergyVector[i]);
	
	bool stopp = false;
	bool stopHe = false;
	if (Eel < EleptHP[0] || Eel > EleptHP[200] || Epr > ET[Nelectrons-1] || Epr < ET[0]) {
	  spall_apel[coupleel][i][j] = 0.0;
	  spall_apel[coupleelHe][i][j] = 0.0;
	  spall_apel[couplepos][i][j] = 0.0;
	  spall_apel[coupleposHe][i][j] = 0.0;
	  //   continue;
	  stopp = true;
	  stopHe = true;
	}
	
	if (Eel+MeleGeV > Epr + mp) {
	  spall_apel[coupleel][i][j] = 0.0;
	  spall_apel[couplepos][i][j] = 0.0;
	  stopp = true;
	  //                continue;
	}
	if (Eel+MeleGeV > 4.0*(Epr + mp)) {
	  spall_apel[coupleel][i][j] = 0.0;
	  spall_apel[coupleelHe][i][j] = 0.0;
	  spall_apel[couplepos][i][j] = 0.0;
	  spall_apel[coupleposHe][i][j] = 0.0;
	  stopp = true;
	  stopHe = true;
	}
	
	int i_lept = int(floor(log(Eel/EleptHP[0])/DBlogHP));
	if (i_lept > dimlept-1) i_lept = dimlept-1;
	
	
	
	//double valuefix = gsl_spline_eval(spline_El_pp, Eel, acc);
	//double valueup = gsl_spline_eval(spline_El_pp_up, Eel, acc);
	
	double t;
	
	if (!stopp) {
          
	  //t = (Eel-EleptHP[i_lept])/(EleptHP[i_lept+1]-EleptHP[i_lept]);
	  valuefix = gsl_spline_eval(spline_El_pp, Eel, acc);//Vectorlept_El_pp[i_lept]*(1.0-t) + Vectorlept_El_pp[i_lept+1]*t;
	  valueup = gsl_spline_eval(spline_El_pp_up, Eel, acc);//Vectorlept_El_pp_up[i_lept]*(1.0-t) + Vectorlept_El_pp_up[i_lept+1]*t;
	  
	  cs_pp = std::max(0.0,valuefix*(1-u) + valueup*u);
	  
	  //cout << Epr << " " << Eel << " " << cs_pp << endl;
	  
	  spall_apel[coupleel][i][j] = (Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_pp)*interp_prodxsec*factorelpos/1000.0;
	  
	  //valuefix = gsl_spline_eval(spline_Pos_pp, Eel, acc);
	  //valueup = gsl_spline_eval(spline_Pos_pp_up, Eel, acc);
	  valuefix = gsl_spline_eval(spline_Pos_pp, Eel, acc);//Vectorlept_Pos_pp[i_lept]*(1.0-t) + Vectorlept_Pos_pp[i_lept+1]*t;
	  valueup = gsl_spline_eval(spline_Pos_pp_up, Eel, acc);//Vectorlept_Pos_pp_up[i_lept]*(1.0-t) + Vectorlept_Pos_pp_up[i_lept+1]*t;
	  cs_pp = std::max(0.0,valuefix*(1-u) + valueup*u);
	  //cout << cs_pp << " ";
	  
	  spall_apel[couplepos][i][j] = 4.0*(Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_pp)*interp_prodxsec*factorelpos/1000.0;   // H
	  
	}
	if (!stopHe) {
	  //valuefix = gsl_spline_eval(spline_El_Hep, Eel, acc);
	  //valueup = gsl_spline_eval(spline_El_Hep_up, Eel, acc);
	  valuefix = gsl_spline_eval(spline_El_Hep, Eel, acc);//Vectorlept_El_Hep[i_lept]*(1.0-t) + Vectorlept_El_Hep[i_lept+1]*t;
	  valueup = gsl_spline_eval(spline_El_Hep_up, Eel, acc);//Vectorlept_El_Hep_up[i_lept]*(1.0-t) + Vectorlept_El_Hep_up[i_lept+1]*t;
	  cs_Hep = std::max(0.0,valuefix*(1-u) + valueup*u);
	  
	  //cout << cs_Hep << " ";
	  
	  spall_apel[coupleelHe][i][j] = (Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_Hep)*interp_prodxsecHe*factorelpos/1000.0;	
	  
          
	  //valuefix = gsl_spline_eval(spline_Pos_Hep, Eel, acc);
	  //valueup = gsl_spline_eval(spline_Pos_Hep_up, Eel, acc);
	  valuefix = gsl_spline_eval(spline_Pos_Hep, Eel, acc);//Vectorlept_Pos_Hep[i_lept]*(1.0-t) + Vectorlept_Pos_Hep[i_lept+1]*t;
	  valueup = gsl_spline_eval(spline_Pos_Hep_up, Eel, acc);//Vectorlept_Pos_Hep_up[i_lept]*(1.0-t) + Vectorlept_Pos_Hep_up[i_lept+1]*t;
	  cs_Hep = std::max(0.0,valuefix*(1-u) + valueup*u);
	  //cout << cs_Hep << endl;
	  
	  spall_apel[coupleposHe][i][j] = 4.0*(Epr)*sqrt(1.0-1.0/pow(1.0+Epr/mp,2))*(cs_Hep)*interp_prodxsecHe*factorelpos/1000.0;   // He
	}
	
      }
      //double sum = 0;
      //for (int i = 0; i < dimEn; i++) sum += spall_apel[coupleel][i][j];
      //cout << sum << endl;
    }
    
    //outfile.close();
    
  }
  
  
  void TSpallationNetwork::InitDataTablesPohl() {
    int Nelectrons = 374;
    
    std::ifstream infile;
    infile.open("espectra_eminus.decay.p.matrix.final.data");
    std::cout << "Opening Pohl file espectra_eminus.decay.p.matrix.final" << std::endl;
    if (!infile.is_open()){ std::cout << "problem opening Pohl espectra_eminus.decay.p.matrix.final.data!!"<< std::endl;
      exit(-1);}
    
    char* s = new char[20];
    std::vector<double> dnde(Nelectrons,0.0);

    for (int i = 0; i < 201; i++) {
      for (int j = 0; j < Nelectrons; j++) {
	infile.get(s,16);
	s[11] = 'E';
	dnde[j] = atof(s);
      }
      infile.get(*s);
      ElppPohl.push_back(dnde);
    }
    infile.close();
    
    
    
    infile.open("espectra_eminus.decay.he.matrix.final.data");
    std::cout << "Opening Pohl file espectra_eminus.decay.he.matrix.final" << std::endl;
    if (!infile.is_open()){ std::cout << "problem opening Pohl espectra_eminus.decay.he.matrix.final.data!!"<< std::endl;
      exit(-1);}

    dnde.clear();
    //std::vector<double> dnde(Nelectrons,0.0);
    for (int i = 0; i < 201; i++) {
      for (int j = 0; j < Nelectrons; j++) {
	infile.get(s,16);
	s[11] = 'E';
	dnde[j] = atof(s);
      }
      infile.get(*s);
      ElHepPohl.push_back(dnde);
    }
    infile.close();
    
    
    infile.open("espectra_eplus.decay.p.matrix.final.data");
    std::cout << "Opening Pohl file espectra_eplus.decay.p.matrix.final" << std::endl;
    if (!infile.is_open()){ std::cout << "problem opening Pohl espectra_eplus.decay.p.matrix.final.data!!"<< std::endl;
      exit(-1);}

    dnde.clear();
    //std::vector<double> dnde(Nelectrons,0.0);
    for (int i = 0; i < 201; i++) {
      for (int j = 0; j < Nelectrons; j++) {
	infile.get(s,16);
	s[11] = 'E';
	dnde[j] = atof(s);
      }
      infile.get(*s);
      PosppPohl.push_back(dnde);
    }
    infile.close();
    
    
    infile.open("espectra_eplus.decay.he.matrix.final.data");
    std::cout << "Opening Pohl file espectra_eplus.decay.he.matrix.final" << std::endl;
    if (!infile.is_open()){ std::cout << "problem opening Pohl espectra_eplus.decay.he.matrix.final.data!!"<< std::endl;
      exit(-1);}

    
    dnde.clear();
    //std::vector<double> dnde(Nelectrons,0.0);
    for (int i = 0; i < 201; i++) {
      for (int j = 0; j < Nelectrons; j++) {
	infile.get(s,16);
	s[11] = 'E';
	dnde[j] = atof(s);
      }
      infile.get(*s);
      PosHepPohl.push_back(dnde);
    }
    infile.close();
    
       
    infile.open("prodxsection.p.matrix.final.data");
    std::cout << "Opening Pohl file prodxsection.p.matrix.final" << std::endl;
    if (!infile.is_open()){ std::cout << "problem opening Pohl prodxsection.p.matrix.final.data!!"<< std::endl;
      exit(-1);}
    
    ProdXsec = std::vector<double>(2*Nelectrons, 0);
    for (int j = 0; j < Nelectrons; j++) {
      infile.get(s,16);
      s[11] = 'E';
      ProdXsec[j] = atof(s);
    }
    infile.close();
    
    
    infile.open("prodxsection.he.matrix.final.data");
    std::cout << "Opening Pohl file prodxsection.he.matrix.final" << std::endl;
    if (!infile.is_open()){ std::cout << "problem opening Pohl prodxsection.he.matrix.final.data!!"<< std::endl;
      exit(-1);}
    
    for (int j = 0; j < Nelectrons; j++) {
      infile.get(s,16);
      s[11] = 'E';
      ProdXsec[Nelectrons+j] = atof(s);
    }
    infile.close();
    
    delete [] s;
    return ;
  }  //End Pohl InitXSecs
  

void TSpallationNetwork::InitLeptUser(double factorelpos, Input* inn){  // Energies in the tables must be given in GeV and cross sections in mb/GeV!
  dimEn = DRAGONEnergyVector.size();
  double a,b,c;
  std::string d;
  
  std::vector<double> Matrix_El_pp;
  std::vector<double> Matrix_El_pHe;
  std::vector<double> Matrix_El_Hep;
  std::vector<double> Matrix_El_HeHe;
  
  std::ifstream infile;
  std::vector <double> Eprim, Esec;
  
  infile.open(inn->User_pp_El);
  std::cout << "Opening User pp -> e- file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User pp -> e- file!!"<< std::endl;
    exit(-1);}
  
  int fla = 0, counter = 0;
  //for (int i = 0; i < 1; i++) getline(infile, d);//skipping the header
  
  Eprim.push_back(0.);
  Esec.push_back(0.);
  while(getline(infile, d)){
    infile >> a >> b >> c;
    if (b == Esec.front())
      fla = 1;
    
    if (fla == 0) Esec.push_back(b);
    
    if (Eprim.back() != a)
      Eprim.push_back(a);
    
    if (counter == 0){
      Eprim.erase(Eprim.begin());
      Esec.erase(Esec.begin());
    }
    counter++;
    Matrix_El_pp.push_back(c);
  }
  infile.close();

    
  infile.open(inn->User_pHe_El);
  std::cout << "Opening User pHe -> e- file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User pHe -> e- file!!"<< std::endl;
    exit(-1);}
  
  while(getline(infile, d)){
    infile >> a >> b >> c;
    Matrix_El_pHe.push_back(c);}
  infile.close();
  
  infile.open(inn->User_Hep_El);
  std::cout << "Opening User Hep -> e- file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User Hep -> e- file!! " <<  std::endl;
    exit(-1);}
  
  while(getline(infile, d)){
    infile >> a >> b >> c;
    Matrix_El_Hep.push_back(c); }
  infile.close();
  
  
  infile.open(inn->User_HeHe_El);
  std::cout << "Opening User HeHe -> e- file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User HeHe -> e- file!!"<< std::endl;
    exit(-1);}
  
  while(getline(infile, d)){
    infile >> a >> b >> c;
    Matrix_El_HeHe.push_back(c); }
  infile.close();
  
  std::cout << "Filling User-defined e- tables \n" << std::endl;

    
  double inx;
  std::vector <std::vector <double> > vec_pEl, vec_HeEl;
  for (int e = 0; e < Eprim.size(); e++){
    std::vector <double> vec_pEl_, vec_HeEl_;
    for (int u = 0; u < Esec.size(); u++){
      inx = e*Esec.size() + u;
      vec_pEl_.push_back(Matrix_El_pp[inx] + He_abundance*Matrix_El_pHe[inx]);
      vec_HeEl_.push_back(Matrix_El_Hep[inx] + He_abundance*Matrix_El_HeHe[inx]);

    }

    vec_pEl.push_back(vec_pEl_);
    vec_HeEl.push_back(vec_HeEl_);     
  }


  std::vector < std::vector <double> > vec_pEl_int, vec_HeEl_int;
  
  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
    std::vector <double> vec_pEl_int_, vec_HeEl_int_;
    for (int u = 0; u < Esec.size(); u++){
      int je = 0, jee;
      double t;

      if (DRAGONEnergyVector[e] > Eprim.back()){  //Extrapolation
	je = Eprim.size()-1; 	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (DRAGONEnergyVector[e] > 7.*Eprim.back())
	  t = (7.*Eprim.back()  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      else if ( DRAGONEnergyVector[e] < Eprim.front()){
	je = 1; jee = 0;
	t =  (DRAGONEnergyVector[e] - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (Eprim.front()/7. > DRAGONEnergyVector[e])
	  t = (Eprim[0]/7. - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      else{ //Interpolation
	while (DRAGONEnergyVector[e] > Eprim[je])
	  je++;
	
	if ((Eprim[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim[je-1])) jee = je-1; 
	else
	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      vec_pEl_int_.push_back(std::max(vec_pEl[jee][u] + t*(vec_pEl[je][u] - vec_pEl[je-1][u]), 0.));
      vec_HeEl_int_.push_back(std::max(vec_HeEl[jee][u] + t*(vec_HeEl[je][u] - vec_HeEl[je-1][u]), 0.));	
    }
    
    vec_pEl_int.push_back(vec_pEl_int_);
    vec_HeEl_int.push_back(vec_HeEl_int_);
  }

  
  std::pair<int,int> coupleElpr(1001, -1000); // Electrons from protons
  std::pair<int,int> coupleElHe(2004, -1000); // Electrons from He
  
  std::pair<int,int> coupleElC12(6012, -1000); // Electrons from Carbon
  std::pair<int,int> coupleElC13(6013, -1000); // Electrons from Carbon13
  std::pair<int,int> coupleElO16(8016, -1000); // Electrons from Oxygen
  std::pair<int,int> coupleElN14(7014, -1000); // Electrons from Nitrogen14
  std::pair<int,int> coupleElN15(7015, -1000); // Electrons from Nitrogen15
  std::pair<int,int> coupleElNe(10020, -1000); // Electrons from Neon
  std::pair<int,int> coupleElMg(12024, -1000); // Electrons from Magnesium
  std::pair<int,int> coupleElSi(14028, -1000); // Electrons from Silicon
  std::pair<int,int> coupleElB10(5010, -1000); // Electrons from B10
  std::pair<int,int> coupleElB11(5011, -1000); // Electrons from B11
  std::pair<int,int> coupleElBe10(4010, -1000); // Electrons from Be10
  std::pair<int,int> coupleElBe9(4009, -1000); // Electrons from Be9
  std::pair<int,int> coupleElBe7(4007, -1000); // Electrons from Be7
  std::pair<int,int> coupleElLi7(3007, -1000); // Electrons from Nitrogen14
  std::pair<int,int> coupleElLi6(3006, -1000); // Electrons from Nitrogen14
	
  double El_A_Scaling = inn->El_A_Scal; //0.9; //2.0/3.0;
  for (int e = 0; e < dimEn; e++){
    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
      
      int je = 0, jee;
      double t;
      if (DRAGONEnergyVector[u] > Esec.back()){  //Extrapolation
	je = Esec.size()-1; 	  jee = je;
	t =  (DRAGONEnergyVector[u]  - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (DRAGONEnergyVector[u] > 7.*Esec.back())
	  t = (7.*Esec.back()  - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else if ( DRAGONEnergyVector[u] < Esec.front()){
	je = 1; jee = 0;
	t = (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (Esec.front()/7. > DRAGONEnergyVector[u])
	  t = (Esec[0]/7. - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else{ //Interpolation
	while (DRAGONEnergyVector[u] > Esec[je])
	  je++;
	
	if ((Esec[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec[je-1])) jee = je-1; 
	else
	  jee = je;	
	
	t =  (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      double pEl = std::max(vec_pEl_int[e][jee] + t*(vec_pEl_int[e][je] - vec_pEl_int[e][je-1]), 0.);
      double HeEl = std::max(vec_HeEl_int[e][jee] + t*(vec_HeEl_int[e][je] - vec_HeEl_int[e][je-1]), 0.);

      spall_apel[coupleElpr][u][e] = DRAGONEnergyVector[e]*factorelpos * pEl/1000.;
      spall_apel[coupleElHe][u][e] = 4.*DRAGONEnergyVector[e]*factorelpos * HeEl/1000.;
      if (El_A_Scaling > 0.) {
	  //Factor A* because we integrate for energies of the primary particle! see eqs. 3.10 and 3.11 of DRAGON2-2
	  spall_apel[coupleElC12][u][e] = 12.*DRAGONEnergyVector[e]*factorelpos * pow(12.0, El_A_Scaling)*pEl/1000.; 
	  spall_apel[coupleElC13][u][e] = 13.*DRAGONEnergyVector[e]*factorelpos * pow(13.0, El_A_Scaling)*pEl/1000.; 
	  spall_apel[coupleElN14][u][e] = 14.*DRAGONEnergyVector[e]*factorelpos * pow(14.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElN15][u][e] = 15.*DRAGONEnergyVector[e]*factorelpos * pow(15.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElO16][u][e] = 18.*DRAGONEnergyVector[e]*factorelpos * pow(18.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElNe][u][e] = 20.*DRAGONEnergyVector[e]*factorelpos * pow(20.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElMg][u][e] = 24.*DRAGONEnergyVector[e]*factorelpos * pow(24.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElSi][u][e] = 28.*DRAGONEnergyVector[e]*factorelpos * pow(28.0, El_A_Scaling)*pEl/1000.;
	  
	  spall_apel[coupleElB11][u][e] = 11.*DRAGONEnergyVector[e]*factorelpos * pow(11.0, El_A_Scaling)*pEl/1000.; 
	  spall_apel[coupleElB10][u][e] = 10.*DRAGONEnergyVector[e]*factorelpos * pow(10.0, El_A_Scaling)*pEl/1000.; 
	  spall_apel[coupleElBe10][u][e] = 10.*DRAGONEnergyVector[e]*factorelpos * pow(10.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElBe9][u][e] = 9.*DRAGONEnergyVector[e]*factorelpos * pow(9.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElBe7][u][e] = 7.*DRAGONEnergyVector[e]*factorelpos * pow(7.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElLi7][u][e] = 7.*DRAGONEnergyVector[e]*factorelpos * pow(7.0, El_A_Scaling)*pEl/1000.;
	  spall_apel[coupleElLi6][u][e] = 6.*DRAGONEnergyVector[e]*factorelpos * pow(6.0, El_A_Scaling)*pEl/1000.; 
	}
    }
    
  }
  
  // Positron part
  double e,f,g;
  std::string h;
  
  std::vector<double> Matrix_Pos_pp;
  std::vector<double> Matrix_Pos_pHe;
  std::vector<double> Matrix_Pos_Hep;
  std::vector<double> Matrix_Pos_HeHe;
    
  infile.open(inn->User_pp_Pos);
  std::cout << "Opening User pp -> e+ file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User pp -> e+ file!!"<< std::endl;
    exit(-1);}

  while(getline(infile, h)){
    infile >> e >> f >> g;
    Matrix_Pos_pp.push_back(g);}
  infile.close();
  
  
  infile.open(inn->User_pHe_Pos);
  std::cout << "Opening User pHe -> e+ file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User pHe -> e+ file!!"<< std::endl;
    exit(-1);}
  
  while(getline(infile, h)){
    infile >> e >> f >> g;
    Matrix_Pos_pHe.push_back(g); }
  infile.close();
  
  
  infile.open(inn->User_Hep_Pos);
  std::cout << "Opening User Hep -> e+ file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User Hep -> e+ file!!"<< std::endl;
    exit(-1);}
  
  while(getline(infile, h)){
    infile >> e >> f >> g;
    Matrix_Pos_Hep.push_back(g); }
  infile.close();
  
  
  infile.open(inn->User_HeHe_Pos);
  std::cout << "Opening User HeHe -> e+ file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening User HeHe -> e+ file!!"<< std::endl;
    exit(-1);}
  
  while(getline(infile, h)){
    infile >> e >> f >> g;
    Matrix_Pos_HeHe.push_back(g); }
  infile.close();
  
  std::cout << "Filling User-defined e+ tables \n" << std::endl;

  std::pair<int,int> couplePospr(1001, 1000); // Positrons from protons
  std::pair<int,int> couplePosHe(2004, 1000); // Positrons from He
 
  std::pair<int,int> couplePosC12(6012, 1000); // Positrons from Carbon
  std::pair<int,int> couplePosC13(6013, 1000); // Positrons from Carbon13
  std::pair<int,int> couplePosO16(8016, 1000); // Positrons from Oxygen
  std::pair<int,int> couplePosN14(7014, 1000); // Positrons from Nitrogen14
  std::pair<int,int> couplePosN15(7015, 1000); // Positrons from Nitrogen15
  std::pair<int,int> couplePosNe(10020, 1000); // Positrons from Neon
  std::pair<int,int> couplePosMg(12024, 1000); // Positrons from Magnesium
  std::pair<int,int> couplePosSi(14028, 1000); // Positrons from Silicon
  std::pair<int,int> couplePosB10(5010, 1000); // Positrons from B10
  std::pair<int,int> couplePosB11(5011, 1000); // Positrons from B11
  std::pair<int,int> couplePosBe10(4010, 1000); // Positrons from Be10
  std::pair<int,int> couplePosBe9(4009, 1000); // Positrons from Be9
  std::pair<int,int> couplePosBe7(4007, 1000); // Positrons from Be7
  std::pair<int,int> couplePosLi7(3007, 1000); // Positrons from Nitrogen14
  std::pair<int,int> couplePosLi6(3006, 1000); // Positrons from Nitrogen14

  
  inx = 0;
  std::vector <std::vector <double> > vec_pPos, vec_HePos;
  for (int e = 0; e < Eprim.size(); e++){
    std::vector <double> vec_pPos_, vec_HePos_;
    for (int u = 0; u < Esec.size(); u++){
      inx = e*Esec.size() + u;
      vec_pPos_.push_back(Matrix_Pos_pp[inx] + He_abundance*Matrix_Pos_pHe[inx]);
      vec_HePos_.push_back(Matrix_Pos_Hep[inx] + He_abundance*Matrix_Pos_HeHe[inx]);

    }

    vec_pPos.push_back(vec_pPos_);
    vec_HePos.push_back(vec_HePos_);     
  }

  std::vector < std::vector <double> > vec_pPos_int, vec_HePos_int;
  
  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
    std::vector <double> vec_pPos_int_, vec_HePos_int_;
    for (int u = 0; u < Esec.size(); u++){
      int je = 0, jee;
      double t;
      if (DRAGONEnergyVector[e] > Eprim.back()){  //Extrapolation
	je = Eprim.size()-1; 	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (DRAGONEnergyVector[e] > 7.*Eprim.back())
	  t = (7.*Eprim.back()  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      else if ( DRAGONEnergyVector[e] < Eprim.front()){
	je = 1; jee = 0;
	t =  (DRAGONEnergyVector[e] - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (Eprim.front()/7. > DRAGONEnergyVector[e])
	  t = (Eprim[0]/7. - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      
      else{ //Interpolation
	while (DRAGONEnergyVector[e] > Eprim[je])
	  je++;
	
	if ((Eprim[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim[je-1])) jee = je-1; 
	else
	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }

      vec_pPos_int_.push_back(std::max(vec_pPos[jee][u] + t*(vec_pPos[je][u] - vec_pPos[je-1][u]), 0.));
      vec_HePos_int_.push_back(std::max(vec_HePos[jee][u] + t*(vec_HePos[je][u] - vec_HePos[je-1][u]), 0.));	
    }
    
    vec_pPos_int.push_back(vec_pPos_int_);
    vec_HePos_int.push_back(vec_HePos_int_);
  }
  
	
  double Pos_A_Scaling = inn->Pos_A_Scal; //0.9; //2.0/3.0;

  for (int e = 0; e < dimEn; e++){
    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
      
      int je = 0, jee;
      double t;
      if (DRAGONEnergyVector[u] > Esec.back()){  //Extrapolation
	je = Esec.size()-1; 	  jee = je;
	t =  (DRAGONEnergyVector[u]  - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (DRAGONEnergyVector[u] > 7.*Esec.back())
	  t = (7.*Esec.back()  - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else if ( DRAGONEnergyVector[u] < Esec.front()){
	je = 1; jee = 0;
	t = (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (Esec.front()/7. > DRAGONEnergyVector[u])
	  t = (Esec[0]/7. - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else{ //Interpolation
	while (DRAGONEnergyVector[u] > Esec[je])
	  je++;
	
	if ((Esec[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec[je-1])) jee = je-1; 
	else
	  jee = je;	
	
	t =  (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      double pPos = std::max(vec_pPos_int[e][jee] + t*(vec_pPos_int[e][je] - vec_pPos_int[e][je-1]), 0.);
      double HePos = std::max(vec_HePos_int[e][jee] + t*(vec_HePos_int[e][je] - vec_HePos_int[e][je-1]), 0.);

      spall_apel[couplePospr][u][e] = DRAGONEnergyVector[e]*factorelpos * pPos/1000.;
      spall_apel[couplePosHe][u][e] = 4.*DRAGONEnergyVector[e]*factorelpos * HePos/1000.;

      if (Pos_A_Scaling > 0.){
	  //Factor A* because we integrate for energies of the primary particle! see eqs. 3.10 and 3.11 of DRAGON2-2
	  spall_apel[couplePosC12][u][e] = 12.*DRAGONEnergyVector[e]*factorelpos * pow(12.0, Pos_A_Scaling)*pPos/1000.; 
	  spall_apel[couplePosC13][u][e] = 13.*DRAGONEnergyVector[e]*factorelpos * pow(13.0, Pos_A_Scaling)*pPos/1000.; 
	  spall_apel[couplePosN14][u][e] = 14.*DRAGONEnergyVector[e]*factorelpos * pow(14.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosN15][u][e] = 15.*DRAGONEnergyVector[e]*factorelpos * pow(15.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosO16][u][e] = 18.*DRAGONEnergyVector[e]*factorelpos * pow(18.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosNe][u][e] = 20.*DRAGONEnergyVector[e]*factorelpos * pow(20.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosMg][u][e] = 24.*DRAGONEnergyVector[e]*factorelpos * pow(24.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosSi][u][e] = 28.*DRAGONEnergyVector[e]*factorelpos * pow(28.0, Pos_A_Scaling)*pPos/1000.;
	  
	  spall_apel[couplePosB11][u][e] = 11.*DRAGONEnergyVector[e]*factorelpos * pow(11.0, Pos_A_Scaling)*pPos/1000.; 
	  spall_apel[couplePosB10][u][e] = 10.*DRAGONEnergyVector[e]*factorelpos * pow(10.0, Pos_A_Scaling)*pPos/1000.; 
	  spall_apel[couplePosBe10][u][e] = 10.*DRAGONEnergyVector[e]*factorelpos * pow(10.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosBe9][u][e] = 9.*DRAGONEnergyVector[e]*factorelpos * pow(9.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosBe7][u][e] = 7.*DRAGONEnergyVector[e]*factorelpos * pow(7.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosLi7][u][e] = 7.*DRAGONEnergyVector[e]*factorelpos * pow(7.0, Pos_A_Scaling)*pPos/1000.;
	  spall_apel[couplePosLi6][u][e] = 6.*DRAGONEnergyVector[e]*factorelpos * pow(6.0, Pos_A_Scaling)*pPos/1000.; 
	}
    }

  }
} //End of InitUserLept




void TSpallationNetwork::InitXSecWinkler(double factorelpos, Input* inn){
  dimEn = DRAGONEnergyVector.size();
  double a,b,c;
  std::string d;
  
  std::vector<double> Matrix_Ap_pp;
  std::vector<double> Matrix_Ap_pHe;
  std::vector<double> Matrix_Ap_Hep;
  std::vector<double> Matrix_Ap_HeHe;
  
  std::ifstream infile;
  std::vector <double> Eprim, Esec;
  
  infile.open("data/Winkler-ap_pp.dat");
  //infile.open("data/ap-xsecs-tables/Feng2016_Beta_EPOS-ap_pp.dat");
  std::cout << "Opening Winkler pp -> ap file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening Winkler pp -> ap file!!"<< std::endl;
    exit(-1);}
  
  int fla = 0, counter = 0;
  for (int i = 0; i < 6; i++) getline(infile, d); //skipping the header
  
  Eprim.push_back(0.);
  Esec.push_back(0.);
  while(infile >> a >> b >> c){
    if (b == Esec.front())
      fla = 1;
    
    if (fla == 0) Esec.push_back(b);
    
    if (Eprim.back() != a)
      Eprim.push_back(a);
    
    if (counter == 0){
      Eprim.erase(Eprim.begin());
      Esec.erase(Esec.begin());
    }
    counter++;
    Matrix_Ap_pp.push_back(c);
  }
  infile.close();
  
  
  infile.open("data/Winkler-ap_pHe.dat");
  //infile.open("data/ap-xsecs-tables/Feng2016_Beta_EPOS-ap_pHe.dat");
  std::cout << "Opening Winkler pHe -> ap file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening Winkler pHe -> ap file!!"<< std::endl;
    exit(-1);}
  
  for (int i = 0; i < 6; i++) getline(infile, d); //skipping the header
  while(infile >> a >> b >> c){ Matrix_Ap_pHe.push_back(c); }
  infile.close();
  
  
  infile.open("data/Winkler-ap_Hep.dat");
  //infile.open("data/ap-xsecs-tables/Feng2016_Beta_EPOS-ap_Hep.dat");
  std::cout << "Opening Winkler Hep -> ap file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening Winkler Hep -> ap file!!"<< std::endl;
    exit(-1);}
  
  for (int i = 0; i < 6; i++) getline(infile, d); //skipping the header
  while(infile >> a >> b >> c){ Matrix_Ap_Hep.push_back(c); }
  infile.close();
  
  
  infile.open("data/Winkler-ap_HeHe.dat");
  //infile.open("data/ap-xsecs-tables/Feng2016_Beta_EPOS-ap_HeHe.dat");
  std::cout << "Opening Winkler HeHe -> ap file" << std::endl;
  if (!infile.is_open()){ std::cout << "problem opening Winkler HeHe -> ap file!!"<< std::endl;
    exit(-1);}
  
  for (int i = 0; i < 6; i++) getline(infile, d); //skipping the header
  while(infile >> a >> b >> c){ Matrix_Ap_HeHe.push_back(c); }
  infile.close();
  
  std::cout << "Filling Winkler ap tables" << std::endl;

  //std::ofstream outfileWin;
  
  //outfileWin.open ("WinklerApXSec.txt");
  //outfileWin << "PrimEnergy(GeV)\t" << "Ap_Energy(GeV)\t" << "XSec_p+H (mbarn)\t" << "XSec_p+He (mbarn)\t" << "total_p+gas (mbarn)\t" << "XSec_He+H (mbarn)\t" << "XSec_He+He (mbarn)\t" << "total_He+gas (mbarn)\t\n";

  double inx;
  std::vector <std::vector <double> > vec_pap, vec_Heap;
  for (int e = 0; e < Eprim.size(); e++){
    std::vector <double> vec_pap_, vec_Heap_;
    for (int u = 0; u < Esec.size(); u++){
      inx = e*Esec.size() + u;
      vec_pap_.push_back(Matrix_Ap_pp[inx] + He_abundance*Matrix_Ap_pHe[inx]);
      vec_Heap_.push_back(Matrix_Ap_Hep[inx] + He_abundance*Matrix_Ap_HeHe[inx]);

      //outfileWin << Eprim[e] << "\t\t" << Esec[u] << "\t\t" << Matrix_Ap_pp[inx] << "\t\t" << Matrix_Ap_pHe[inx] << "\t\t" << Matrix_Ap_pp[inx] + He_abundance*Matrix_Ap_pHe[inx] << "\t\t" << Matrix_Ap_Hep[inx] << "\t\t" << Matrix_Ap_HeHe[inx] << "\t\t" << Matrix_Ap_Hep[inx] + He_abundance*Matrix_Ap_HeHe[inx] << "\n";
    }

    //outfileWin << "\n\n\n";
    vec_pap.push_back(vec_pap_);
    vec_Heap.push_back(vec_Heap_);     
  }
  //outfileWin.close();
  std::vector < std::vector <double> > vec_pap_int, vec_Heap_int;
  
  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
    std::vector <double> vec_pap_int_, vec_Heap_int_;
    for (int u = 0; u < Esec.size(); u++){
      int je = 0, jee;
      double t;
      if (DRAGONEnergyVector[e] > Eprim.back()){  //Extrapolation
	je = Eprim.size()-1; 	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (DRAGONEnergyVector[e] > 7.*Eprim.back())
	  t = (7.*Eprim.back()  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      else if ( DRAGONEnergyVector[e] < Eprim.front()){
	je = 1; jee = 0;
	t =  (DRAGONEnergyVector[e] - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (Eprim.front()/7. > DRAGONEnergyVector[e])
	  t = (Eprim[0]/7. - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      
      else{ //Interpolation
	while (DRAGONEnergyVector[e] > Eprim[je])
	  je++;
	
	if ((Eprim[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim[je-1])) jee = je-1; 
	else
	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }

      vec_pap_int_.push_back(std::max(vec_pap[jee][u] + t*(vec_pap[je][u] - vec_pap[je-1][u]), 0.));
      vec_Heap_int_.push_back(std::max(vec_Heap[jee][u] + t*(vec_Heap[je][u] - vec_Heap[je-1][u]), 0.));	
    }
    
    vec_pap_int.push_back(vec_pap_int_);
    vec_Heap_int.push_back(vec_Heap_int_);
  }
  
  std::pair<int,int> coupleappr(1001,-999);  // Secondary antiprotons, from protons
  std::pair<int,int> coupleapHe(2004,-999);  // Secondary antiprotons, from Helium
  std::pair<int,int> coupleapLi6(3006,-999);  // 
  std::pair<int,int> coupleapLi7(3007,-999);  //
  std::pair<int,int> coupleapBe7(4007,-999);  //
  std::pair<int,int> coupleapBe9(4009,-999);  //
  std::pair<int,int> coupleapBe10(4010,-999);  //
  std::pair<int,int> coupleapB10(5010,-999);  // 
  std::pair<int,int> coupleapB11(5011,-999);  // 
  std::pair<int,int> coupleapC12(6012,-999);  // Secondary antiprotons, from C12
  std::pair<int,int> coupleapC13(6013,-999);  // Secondary antiprotons, from C13
  std::pair<int,int> coupleapC14(6014,-999);  // Secondary antiprotons, from C14
  std::pair<int,int> coupleapO16(8016,-999);  // Secondary antiprotons, from O
  std::pair<int,int> coupleapN14(7014,-999);  // Secondary antiprotons, from N14
  std::pair<int,int> coupleapN15(7015,-999);  // Secondary antiprotons, from N15
  std::pair<int,int> coupleapNe(10020,-999);  // Secondary antiprotons, from Ne
  std::pair<int,int> coupleapMg(12024,-999);  // Secondary antiprotons, from Mg
  std::pair<int,int> coupleapSi(14028,-999);  // Secondary antiprotons, from Si

  //std::ofstream oufileWin;
  //oufileWin.open ("WinklerApXSec_interp.txt");
  //oufileWin << "PrimEnergy(GeV)\t" << "Ap_Energy(GeV)\t" << "total_p+gas (mbarn)\t" << "total_He+gas (mbarn)\t\n";

  double Ap_XS_scal = inn->Ap_XS_Scal; //1.;
  double Ap_A_Scaling = inn->Ap_A_Scal; //0.9; //2.0/3.0;

  for (int e = 0; e < dimEn; e++){
    //std::vector <double> vec_pap_, vec_Heap_;
    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
      
      int je = 0, jee;
      double t;
      if (DRAGONEnergyVector[u] > Esec.back()){  //Extrapolation
	je = Esec.size()-1; 	  jee = je;
	t =  (DRAGONEnergyVector[u]  - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (DRAGONEnergyVector[u] > 7.*Esec.back())
	  t = (7.*Esec.back()  - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else if ( DRAGONEnergyVector[u] < Esec.front()){
	je = 1; jee = 0;
	t = (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (Esec.front()/7. > DRAGONEnergyVector[u])
	  t = (Esec[0]/7. - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else{ //Interpolation
	while (DRAGONEnergyVector[u] > Esec[je])
	  je++;
	
	if ((Esec[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec[je-1])) jee = je-1; 
	else
	  jee = je;	
	
	t =  (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      double pap = std::max(vec_pap_int[e][jee] + t*(vec_pap_int[e][je] - vec_pap_int[e][je-1]), 0.);
      double Heap = std::max(vec_Heap_int[e][jee] + t*(vec_Heap_int[e][je] - vec_Heap_int[e][je-1]), 0.);

      //oufileWin << DRAGONEnergyVector[e] << "\t\t" << DRAGONEnergyVector[u] << "\t\t" << pap << "\t\t" << Heap << "\n";

      spall_apel[coupleappr][u][e] = Ap_XS_scal*DRAGONEnergyVector[e]*factorelpos * pap/1000.;
      spall_apel[coupleapHe][u][e] = Ap_XS_scal*4.*DRAGONEnergyVector[e]*factorelpos * Heap/1000.;
      
      //Factor A* because we integrate for energies of the primary particle! see eqs. 3.10 and 3.11 of DRAGON2-2
      spall_apel[coupleapC12][u][e] = Ap_XS_scal*12.*DRAGONEnergyVector[e]*factorelpos * pow(12.0, Ap_A_Scaling)*pap/1000.; 
      spall_apel[coupleapC13][u][e] = Ap_XS_scal*13.*DRAGONEnergyVector[e]*factorelpos * pow(13.0, Ap_A_Scaling)*pap/1000.; 
      spall_apel[coupleapN14][u][e] = Ap_XS_scal*14.*DRAGONEnergyVector[e]*factorelpos * pow(14.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapN15][u][e] = Ap_XS_scal*15.*DRAGONEnergyVector[e]*factorelpos * pow(15.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapO16][u][e] = Ap_XS_scal*18.*DRAGONEnergyVector[e]*factorelpos * pow(18.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapNe][u][e] = Ap_XS_scal*20.*DRAGONEnergyVector[e]*factorelpos * pow(20.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapMg][u][e] = Ap_XS_scal*24.*DRAGONEnergyVector[e]*factorelpos * pow(24.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapSi][u][e] = Ap_XS_scal*28.*DRAGONEnergyVector[e]*factorelpos * pow(28.0, Ap_A_Scaling)*pap/1000.;

      spall_apel[coupleapB11][u][e] = Ap_XS_scal*11.*DRAGONEnergyVector[e]*factorelpos * pow(11.0, Ap_A_Scaling)*pap/1000.; 
      spall_apel[coupleapB10][u][e] = Ap_XS_scal*10.*DRAGONEnergyVector[e]*factorelpos * pow(10.0, Ap_A_Scaling)*pap/1000.; 
      spall_apel[coupleapBe10][u][e] = Ap_XS_scal*10.*DRAGONEnergyVector[e]*factorelpos * pow(10.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapBe9][u][e] = Ap_XS_scal*9.*DRAGONEnergyVector[e]*factorelpos * pow(9.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapBe7][u][e] = Ap_XS_scal*7.*DRAGONEnergyVector[e]*factorelpos * pow(7.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapLi7][u][e] = Ap_XS_scal*7.*DRAGONEnergyVector[e]*factorelpos * pow(7.0, Ap_A_Scaling)*pap/1000.;
      spall_apel[coupleapLi6][u][e] = Ap_XS_scal*6.*DRAGONEnergyVector[e]*factorelpos * pow(6.0, Ap_A_Scaling)*pap/1000.; 
      
    }

    //oufileWin << "\n\n\n";
    //std::cout << " " <<  spall_apel[coupleapHe][e][60] << " " << spall_apel[coupleappr][e][70] << "  " <<  std::endl;
  }
  //oufileWin.close();
} //End of InitWinkler




void TSpallationNetwork::InitXSecDRAGONAN(double factorelpos, TGalpropXSec GPobject, string Atablepp, string AtablepHe, string Atableapp, string AtableapHe, string AtableHep, string AtableHeHe, bool UseDRAGON_IneXSec, bool Scale_pp, string ANuc){  //TGrid* cord, Input* inp
  dimEn = DRAGONEnergyVector.size();
  double a,b,c;
  std::string d;
  
  std::vector<double> Matrix_AN_pp;
  std::vector<double> Matrix_AN_pHe;
  std::vector<double> Matrix_AN_app;
  std::vector<double> Matrix_AN_apHe;
  std::vector<double> Matrix_AN_Hep;
  std::vector<double> Matrix_AN_HeHe;
  
  std::ifstream infile;
  std::vector <double> Eprim, Esec;
  
  infile.open(Atablepp);
  if (ANuc == "D") std::cout << "Opening pp -> aD file: " << Atablepp << std::endl;
  else if (ANuc == "3He") std::cout << "Opening pp -> a3He file: " << Atablepp << std::endl;
  else { std::cout << "That is not a valid option for DRAGONAN function" << ANuc << std::endl;
    exit(-1);}
  if (!infile.is_open()){ std::cout << "problem opening pp file!!"<< std::endl;
    exit(-1);}
  
  int fla = 0, counter = 0;
  for (int i = 0; i < 7; i++) getline(infile, d); //skipping the header
  
  Eprim.push_back(0.);
  Esec.push_back(0.);

  // protons as projectile
  while(infile >> a >> b >> c){
    if (b == Esec.front())
      fla = 1;
    
    if (fla == 0) Esec.push_back(b);
    
    if (Eprim.back() != a)
      Eprim.push_back(a);
    
    if (counter == 0){
      Eprim.erase(Eprim.begin());
      Esec.erase(Esec.begin());
    }
    counter++;
    Matrix_AN_pp.push_back(c);
  }
  infile.close();
  
  infile.open(AtablepHe);
  if (ANuc == "D") std::cout << "Opening pHe -> aD file: " << AtablepHe <<  std::endl;
  else if (ANuc == "3He") std::cout << "Opening pHe -> a3He file: " << AtablepHe << std::endl;
  else { std::cout << "That is not a valid option for DRAGONAN function" << ANuc << std::endl;
    exit(-1);}
    
  if (!infile.is_open()){ std::cout << "problem opening pHe file!!"<< std::endl;
    exit(-1);}
  
  for (int i = 0; i < 7; i++) getline(infile, d); //skipping the header
  while(infile >> a >> b >> c){ Matrix_AN_pHe.push_back(c); }
  infile.close();

  
  // Antiprotons as projectile
  if (Atableapp != "null"){
    infile.open(Atableapp);
    if (ANuc == "D") std::cout << "Opening app -> aD file: " << Atableapp <<  std::endl;
    else if (ANuc == "3He") std::cout << "Opening app -> a3He file: " << Atableapp << std::endl;
    else { std::cout << "That is not a valid option for DRAGONAN function" << ANuc << std::endl;
      exit(-1);}
    if (!infile.is_open()){ std::cout << "problem opening app file!!"<< std::endl;
      exit(-1);}
  
    for (int i = 0; i < 7; i++) getline(infile, d); //skipping the header
    while(infile >> a >> b >> c){ Matrix_AN_app.push_back(c); }
    infile.close();
  }

  if (AtableapHe != "null"){
    infile.open(AtableapHe);
    if (ANuc == "D") std::cout << "Opening apHe -> aD file: " << AtableapHe <<  std::endl;
    else if (ANuc == "3He") std::cout << "Opening apHe -> a3He file: " << AtableapHe << std::endl;
    else { std::cout << "That is not a valid option for DRAGONAN function" << ANuc << std::endl;
      exit(-1);}
    
    if (!infile.is_open()){ std::cout << "problem opening apHe file!!"<< std::endl;
      exit(-1);}
    
    for (int i = 0; i < 7; i++) getline(infile, d); //skipping the header
    while(infile >> a >> b >> c){ Matrix_AN_apHe.push_back(c); }
    infile.close();
  }
  
  // Helium as projectile  
  infile.open(AtableHep);
  if (ANuc == "D") std::cout << "Opening Hep -> aD file: " << AtableHep <<  std::endl;
  else if (ANuc == "3He") std::cout << "Opening Hep -> a3He file: " << AtableHep << std::endl;
  else { std::cout << "That is not a valid option for DRAGONAN function" << ANuc << std::endl;
    exit(-1);}
  
  if (!infile.is_open()){ std::cout << "problem opening Hep file!!"<< std::endl;
    exit(-1);}
  
  for (int i = 0; i < 7; i++) getline(infile, d); //skipping the header
  while(infile >> a >> b >> c){ Matrix_AN_Hep.push_back(c); }
  infile.close();
  
  
  infile.open(AtableHeHe);
  if (ANuc == "D") std::cout << "Opening HeHe -> aD file: " << AtableHeHe <<  std::endl;
  else if (ANuc == "3He") std::cout << "Opening HeHe -> a3He file: " << AtableHeHe << std::endl;
  else { std::cout << "That is not a valid option for DRAGONAN function" << ANuc << std::endl;
    exit(-1);}
    
  if (!infile.is_open()){ std::cout << "problem opening HeHe file!!"<< std::endl;
    exit(-1);}
  
  for (int i = 0; i < 7; i++) getline(infile, d); //skipping the header
  while(infile >> a >> b >> c){ Matrix_AN_HeHe.push_back(c); }
  infile.close();

  if (UseDRAGON_IneXSec)
    std::cout << "\nMultiplicity tables done..." << std::endl;
  else
    std::cout << "\nAnti-D cross sections tables read..." << std::endl;
  //TGalpropXSec GPobject(cord, inp); 
  double inx;
  std::vector <std::vector <double> > vec_pan, vec_apan, vec_Hean ;
  for (int e = 0; e < Eprim.size(); e++){
    std::vector <double> vec_pan_, vec_apan_, vec_Hean_;
    for (int u = 0; u < Esec.size(); u++){
      inx = e*Esec.size() + u;
      PP_inel = 1.0;
      PA_inel = 1.0;
      aPP_non = 1.0;
      aPA_non = 1.0;
      aPP_ann = 1.0;
      aPA_ann = 1.0;
      
      if (UseDRAGON_IneXSec)
	GPobject.nucleon_cs(2, Eprim[e], 1, 2, 4, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann); //IneCrossSecs in mb!
      
      if (Scale_pp){
	vec_pan_.push_back(Matrix_AN_pp[inx]*PP_inel + pow(4, 0.8)*He_abundance*Matrix_AN_pp[inx]*PA_inel);
	vec_apan_.push_back(0); //nothing to do for antiprotons from protons
	//vec_apan_.push_back(Matrix_AN_app[inx]*(aPP_non + aPP_ann) + pow(4, 0.8)*He_abundance*Matrix_AN_apHe[inx]*(aPA_non + aPA_ann));	
      }
      
      else{
	vec_pan_.push_back(Matrix_AN_pp[inx]*PP_inel + He_abundance*Matrix_AN_pHe[inx]*PA_inel);
	if (Atableapp == "null")
	  vec_apan_.push_back(0);
	else{
	  if (AtableapHe != "null")
	    vec_apan_.push_back(Matrix_AN_app[inx]*(aPP_non + aPP_ann) + He_abundance*Matrix_AN_apHe[inx]*(aPA_non + aPA_ann));
	  else
	    vec_apan_.push_back(Matrix_AN_app[inx]*(aPP_non + aPP_ann));
	}
	
      }
      
      vec_Hean_.push_back(Matrix_AN_Hep[inx]*PA_inel + He_abundance*Matrix_AN_HeHe[inx]*4*PA_inel);
    }
    
    vec_pan.push_back(vec_pan_);
    vec_apan.push_back(vec_apan_);
    vec_Hean.push_back(vec_Hean_);     
  }

  std::vector < std::vector <double> > vec_pan_int, vec_apan_int, vec_Hean_int;
  
  for (int e = 0; e < dimEn; e++){ //Interpolation by nearest neighbor method -> for primary energies
    std::vector <double> vec_pan_int_, vec_apan_int_, vec_Hean_int_;
    for (int u = 0; u < Esec.size(); u++){
      int je = 0, jee;
      double t;
      if (DRAGONEnergyVector[e] > Eprim.back()){  //Extrapolation
	je = Eprim.size()-1; 	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (DRAGONEnergyVector[e] > 7.*Eprim.back())
	  t = (7.*Eprim.back()  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      else if ( DRAGONEnergyVector[e] < Eprim.front()){
	je = 1; jee = 0;
	t =  (DRAGONEnergyVector[e] - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
	if (Eprim.front()/7. > DRAGONEnergyVector[e])
	  t = (Eprim[0]/7. - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }
      
      else{ //Interpolation
	while (DRAGONEnergyVector[e] > Eprim[je])
	  je++;
	
	if ((Eprim[je] - DRAGONEnergyVector[e] ) > (DRAGONEnergyVector[e]  - Eprim[je-1])) jee = je-1; 
	else
	  jee = je;
	t = (DRAGONEnergyVector[e]  - Eprim[jee])/(Eprim[je]-Eprim[je-1]);
      }

      vec_pan_int_.push_back(std::max(vec_pan[jee][u] + t*(vec_pan[je][u] - vec_pan[je-1][u]), 0.));
      vec_apan_int_.push_back(std::max(vec_apan[jee][u] + t*(vec_apan[je][u] - vec_apan[je-1][u]), 0.));
      vec_Hean_int_.push_back(std::max(vec_Hean[jee][u] + t*(vec_Hean[je][u] - vec_Hean[je-1][u]), 0.));	
    }
    
    vec_pan_int.push_back(vec_pan_int_);
    vec_apan_int.push_back(vec_apan_int_);
    vec_Hean_int.push_back(vec_Hean_int_);
  }

  for (int e = 0; e < dimEn; e++){
    for (int u = 0; u < dimEn; u++){ //Interpolation by nearest neighbor method -> for secondary energies
      
      int je = 0, jee;
      double t;
      if (DRAGONEnergyVector[u] > Esec.back()){  //Extrapolation
	je = Esec.size()-1; 	  jee = je;
	t =  (DRAGONEnergyVector[u]  - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (DRAGONEnergyVector[u] > 7.*Esec.back())
	  t = (7.*Esec.back()  - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else if ( DRAGONEnergyVector[u] < Esec.front()){
	je = 1; jee = 0;
	t = (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
	if (Esec.front()/7. > DRAGONEnergyVector[u])
	  t = (Esec[0]/7. - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      else{ //Interpolation
	while (DRAGONEnergyVector[u] > Esec[je])
	  je++;
	
	if ((Esec[je] - DRAGONEnergyVector[u] ) > (DRAGONEnergyVector[u]  - Esec[je-1])) jee = je-1; 
	else
	  jee = je;	
	
	t =  (DRAGONEnergyVector[u] - Esec[jee])/(Esec[je]-Esec[je-1]);
      }
      
      double pan = std::max(vec_pan_int[e][jee] + t*(vec_pan_int[e][je] - vec_pan_int[e][je-1]), 0.);
      double apan = std::max(vec_apan_int[e][jee] + t*(vec_apan_int[e][je] - vec_apan_int[e][je-1]), 0.);
      double Hean = std::max(vec_Hean_int[e][jee] + t*(vec_Hean_int[e][je] - vec_Hean_int[e][je-1]), 0.);

      if (ANuc == "D") {
	std::pair<int,int> coupleadp(1001,-998);  // Secondary antiD, from protons
	std::pair<int,int> coupleadap(-999,-998);  // Secondary antiD, from antiprotons
	std::pair<int,int> coupleadHe(2004,-998);  // Secondary antiD, from Helium
	spall_apel[coupleadp][u][e] = DRAGONEnergyVector[e]*factorelpos * pan/1000.;
	spall_apel[coupleadap][u][e] = DRAGONEnergyVector[e]*factorelpos * apan/1000.;
	spall_apel[coupleadHe][u][e] = 4.*DRAGONEnergyVector[e]*factorelpos * Hean/1000.;
	if (Scale_pp)
	  spall_apel[coupleadHe][u][e] =  pow(4, 0.8)*4*spall_apel[coupleadp][u][e]; // If user wants only to use pp table and scale the He coliisions XSecs, we leave antid cross sections as the proton cross sections and scale the He cross sections by a 4^0.8      //spall_apel[coupleadHe][u][e] *=  pow(4, 0.8);
      }
      else if (ANuc == "3He") {
	std::pair<int,int> couplea3p(1001,-997);  // Secondary antiHe3, from protons
	std::pair<int,int> couplea3ap(-999,-997);  // Secondary antiHe3, from antiprotons
	std::pair<int,int> couplea3He(2004,-997);  // Secondary antiHe3, from Helium
	spall_apel[couplea3p][u][e] = DRAGONEnergyVector[e]*factorelpos * pan/1000.;
	spall_apel[couplea3ap][u][e] = DRAGONEnergyVector[e]*factorelpos * apan/1000.;
	spall_apel[couplea3He][u][e] = 4.*DRAGONEnergyVector[e]*factorelpos * Hean/1000.;
	if (Scale_pp)
	  spall_apel[couplea3He][u][e] =  pow(4, 0.8)*4*spall_apel[couplea3p][u][e]; 
      }
      
      else { std::cout << "That is not a valid option for DRAGONAN function" << ANuc << std::endl;
	exit(-1);}
            
    }  //End of primary energy loop
    //std::cout << " " <<  spall_apel[coupleadHe][e][60] << " " << spall_apel[coupleadp][e][70] << "  " <<  std::endl;
  } // End of secondary energy loop

  if (ANuc == "D") std::cout << "Filling aD cross sections tables... done" << std::endl;
  else if (ANuc == "3He") std::cout << "Filling aHe3 cross sections tables... done" << std::endl;
  
  
} //End of InitAntiD





void TSpallationNetwork::InitXSecKamae(double factorelpos) { 
    
  std::pair<int,int> coupleel(1001, -1000); // Electrons from protons
  std::pair<int,int> couplepos(1001, 1000); // Positrons from protons
  std::pair<int,int> coupleelHe(2004, -1000); // Electrons from He
  std::pair<int,int> coupleposHe(2004, 1000); // Positrons from He
  
  std::pair<int,int> coupleelC12(6012, -1000); // Electrons from Carbon
  std::pair<int,int> coupleposC12(6012, 1000); // Positrons from Carbon
  std::pair<int,int> coupleelC13(6013, -1000); // Electrons from Carbon13
  std::pair<int,int> coupleposC13(6013, 1000); // Positrons from Carbon13
  std::pair<int,int> coupleelO16(8016, -1000); // Electrons from Oxygen
  std::pair<int,int> coupleposO16(8016, 1000); // Positrons from Oxygen
  std::pair<int,int> coupleelN14(7014, -1000); // Electrons from Nitrogen14
  std::pair<int,int> coupleposN14(7014, 1000); // Positrons from Nitrogen14
  std::pair<int,int> coupleelN15(7015, -1000); // Electrons from Nitrogen15
  std::pair<int,int> coupleposN15(7015, 1000); // Positrons from Nitrogen15
  std::pair<int,int> coupleelNe(10020, -1000); // Electrons from Neon
  std::pair<int,int> coupleposNe(10020, 1000); // Positrons from Neon
  std::pair<int,int> coupleelMg(12024, -1000); // Electrons from Magnesium
  std::pair<int,int> coupleposMg(12024, 1000); // Positrons from Magnesium
  std::pair<int,int> coupleelSi(14028, -1000); // Electrons from Silicon
  std::pair<int,int> coupleposSi(14028, 1000); // Positrons from Silicon
  std::pair<int,int> coupleelB10(5010, -1000); // Electrons from B10
  std::pair<int,int> coupleposB10(5010, 1000); // Positrons from B10
  std::pair<int,int> coupleelB11(5011, -1000); // Electrons from B11
  std::pair<int,int> coupleposB11(5011, 1000); // Positrons from B11
  std::pair<int,int> coupleelBe10(4010, -1000); // Electrons from Be10
  std::pair<int,int> coupleposBe10(4010, 1000); // Positrons from Be10
  std::pair<int,int> coupleelBe9(4009, -1000); // Electrons from Be9
  std::pair<int,int> coupleposBe9(4009, 1000); // Positrons from Be9
  std::pair<int,int> coupleelBe7(4007, -1000); // Electrons from Be7
  std::pair<int,int> coupleposBe7(4007, 1000); // Positrons from Be7
  std::pair<int,int> coupleelLi7(3007, -1000); // Electrons from Nitrogen14
  std::pair<int,int> coupleposLi7(3007, 1000); // Positrons from Nitrogen14
  std::pair<int,int> coupleelLi6(3006, -1000); // Electrons from Nitrogen14
  std::pair<int,int> coupleposLi6(3006, 1000); // Positrons from Nitrogen14

  
  dimEn = DRAGONEnergyVector.size();
  double Eel;
  
  //std::ofstream oufileel;
  //std::ofstream oufilepos;

  //oufileel.open ("KamaeElXSec.txt");
  //oufileel << "PrimEnergy(GeV)\t" << "El_Energy(GeV)\t" << "XSec_p+H (barn)\t" << "XSec_p+He (barn)\t" << "total_p+gas (barn)\t" << "XSec_He+H (barn)\t" << "XSec_He+He (barn)\t" << "total_He+gas (mbarn)\t\n";

  //oufilepos.open ("KamaePosXSec.txt");
  //oufilepos << "PrimEnergy(GeV)\t" << "Pos_Energy(GeV)\t" << "XSec_p+H (barn)\t" << "XSec_p+He (barn)\t" << "total_p+gas (barn)\t" << "XSec_He+H (barn)\t" << "XSec_He+He (barn)\t" << "total_He+gas (barn)\t\n";

  double Pos_XS_scal = 1.;
  double Lep_scale = 0.9; //2.0/3.0;
  for (unsigned int j=0; j< dimEn; j++){
    double Epr = DRAGONEnergyVector[j];
    std::vector<double> elpxsec;
    std::vector<double> elHexsec;
    std::vector<double> pospxsec;
    std::vector<double> posHexsec;

    for (unsigned int i = 0; i < dimEn; i++) {
      
      //Eel = std::min(1.e5, DRAGONEnergyVector[i]);
      double Eel = DRAGONEnergyVector[i];

      
      cs_pp = KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);
      
      cs_pHe = pow(4.0, Lep_scale)*KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);
      
      cs_Hep = pow(4.0, Lep_scale)*KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);
      
      cs_HeHe = pow(4.0*4.0, Lep_scale)*KamaeYields::GetSigma(Eel, Epr, ID_ELECTRON);  
      
      spall_apel[coupleel][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;
      spall_apel[coupleelHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;
      
      spall_apel[coupleelC12][i][j] = 12.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(12.0, Lep_scale);
      spall_apel[coupleelC13][i][j] = 13.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(13.0, Lep_scale);
      spall_apel[coupleelN14][i][j] = 14.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(14.0, Lep_scale);
      spall_apel[coupleelN15][i][j] = 15.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(15.0, Lep_scale);
      spall_apel[coupleelO16][i][j] = 16.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(16.0, Lep_scale);
      spall_apel[coupleelNe][i][j] = 20.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(20.0, Lep_scale);
      spall_apel[coupleelMg][i][j] = 24.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(24.0, Lep_scale);
      spall_apel[coupleelSi][i][j] = 28.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(28.0, Lep_scale);
      spall_apel[coupleelB10][i][j] = 10.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(10.0, Lep_scale);
      spall_apel[coupleelB11][i][j] = 11.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(11.0, Lep_scale);
      spall_apel[coupleelBe10][i][j] = 10.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(10.0, Lep_scale);
      spall_apel[coupleelBe9][i][j] = 9.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(9.0, Lep_scale);
      spall_apel[coupleelBe7][i][j] = 7.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(7.0, Lep_scale);
      spall_apel[coupleelLi7][i][j] = 7.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(7.0, Lep_scale);
      spall_apel[coupleelLi6][i][j] = 6.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(6.0, Lep_scale);
										    
      
      //std::cout << "\n\n hehehe " << Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos << " and also: " <<  4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos << std::endl;
      
      cs_pp = Pos_XS_scal*KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);
      
      cs_pHe = Pos_XS_scal*pow(4.0, Lep_scale)*KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);
      
      cs_Hep = Pos_XS_scal*pow(4.0, Lep_scale)*KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);
      
      cs_HeHe = Pos_XS_scal*pow(4.0*4.0, Lep_scale)*KamaeYields::GetSigma(Eel, Epr, ID_POSITRON);

      spall_apel[couplepos][i][j] = Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos;   // H
      spall_apel[coupleposHe][i][j] = 4.0*Epr*(cs_Hep + He_abundance*cs_HeHe)*factorelpos;   // He
      
      spall_apel[coupleposC12][i][j] = 12.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(12.0, Lep_scale);
      spall_apel[coupleposC13][i][j] = 13.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(13.0, Lep_scale);
      spall_apel[coupleposN14][i][j] = 14.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(14.0, Lep_scale);
      spall_apel[coupleposN15][i][j] = 15.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(15.0, Lep_scale);
      spall_apel[coupleposO16][i][j] = 16.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(16.0, Lep_scale);
      spall_apel[coupleposNe][i][j] = 20.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(20.0, Lep_scale);
      spall_apel[coupleposMg][i][j] = 24.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(24.0, Lep_scale);
      spall_apel[coupleposSi][i][j] = 28.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(28.0, Lep_scale);
      spall_apel[coupleposB10][i][j] = 10.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(10.0, Lep_scale);
      spall_apel[coupleposB11][i][j] = 11.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(11.0, Lep_scale);
      spall_apel[coupleposBe10][i][j] = 10.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(10.0, Lep_scale);
      spall_apel[coupleposBe9][i][j] = 9.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(9.0, Lep_scale);
      spall_apel[coupleposBe7][i][j] = 7.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(7.0, Lep_scale);
      spall_apel[coupleposLi7][i][j] = 7.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(7.0, Lep_scale);
      spall_apel[coupleposLi6][i][j] = 6.*Epr*(cs_pp + He_abundance*cs_pHe)*factorelpos * pow(6.0, Lep_scale);
      
      //oufileel << Epr << "\t\t" << Eel << "\t\t" << cs_pp << "\t\t" << cs_pHe << "\t\t" << cs_pp + He_abundance*cs_pHe << "\t\t" << cs_Hep << "\t\t" << cs_HeHe << "\t\t" << cs_Hep + He_abundance*cs_HeHe << "\n";

      //oufilepos << Epr << "\t\t" << Eel << "\t\t" << cs_pp << "\t\t" << cs_pHe << "\t\t" << cs_pp + He_abundance*cs_pHe << "\t\t" << cs_Hep << "\t\t" << cs_HeHe << "\t\t" << cs_Hep + He_abundance*cs_HeHe << "\n";
    }

    //oufileel << "\n\n\n";
    //oufilepos << "\n\n\n";


  }
  
  //oufileel.close();
  //oufilepos.close();

} //End Kamae XSecs

  
  
  std::vector<double> TSpallationNetwork::GetXSecApEl(int ParentID, int DaughterID, int k){
    
    if (k >= DRAGONEnergyVector.size()) {
      std::cout << "The index you required is out of range!, max index is: " << DRAGONEnergyVector.size()-1 << std::endl;
      return std::vector<double>();
    }
    
    std::pair<int,int> input(ParentID, DaughterID);
    std::map< std::pair<int,int>, std::vector< std::vector<double> > >::iterator it = spall_apel.find(input);
    if (it != spall_apel.end()){
      //std::cout << "your index, " << k << " corresponds to a energy of the secondary particle of " << DRAGONEnergyVector[k] << " GeV" << std::endl;
      return (it->second)[k];//spall_apel[input][k];
    }
    else {
      //std::cout << "No channel in TSpallationNetwork. Parent = " << ParentID << " daughter = " << DaughterID << std::endl;
      return std::vector<double>();
    }
    
  }

  
  std::vector<double> TSpallationNetwork::GetXSec(int parentParticleID, int daughterParticleID) {    
    //std::cout << "Getting the total spallation XSec = (XsecH + fhe*XsecHe) from " << parentParticleID << " to " << daughterParticleID << std::endl;
    std::pair<int,int> inputParticlePair(parentParticleID, daughterParticleID);
    
    std::map< std::pair<int,int>, std::vector<double> >::iterator it = spall.find(inputParticlePair);
    if (it != spall.end()){
      //std::cout << "Good spallation cross section for this pair " << std::endl;
      return (*it).second;
    }
    else{
      //std::cout << "\n\n\n\nNo spallation cross section for this pair " << parentParticleID << " " << daughterParticleID << std::endl;
      return std::vector<double> ();
    }
  }
  
void TSpallationNetwork::Retrieve_XSecs(Input* in, std::vector<int>& nuclei, TGrid* Coord){
  
  std::vector<double> beta = Coord->GetBeta();
  const double factor = 1e-27*Clight;
  int dimEn = DRAGONEnergyVector.size();
  std::ofstream outfileX;

  if (He_abundance > 0.){
    if (in->spallationxsec == Webber03) outfileX.open ("TotalXSec_Webber.dat");
    if (in->spallationxsec == DRAGON2) outfileX.open ("TotalXSec_Carmelo.dat");
    if (in->spallationxsec == GalpropXSec) outfileX.open ("TotalXSec_Galprop.dat");
    //if (in->spallationxsec == Fluka) outfileX.open ("TotalXSec_Fluka.dat");
  }
  else{
    if (in->spallationxsec == Webber03) outfileX.open ("HXSec_Webber.dat");
    if (in->spallationxsec == DRAGON2) outfileX.open ("HXSec_Carmelo.dat");
    if (in->spallationxsec == GalpropXSec) outfileX.open ("HXSec_Galprop.dat");
    //if (in->spallationxsec == Fluka) outfileX.open ("HXSec_Fluka.dat");
  }
    
  outfileX << "Zi.Ai\t" << "Zo.Ao\t\t" <<"Energy(GeV)\t\t" << "XSec(mb)\n";
  for (int iloop = 0; iloop < nuclei.size(); ++iloop) {
    
    // If antiprotons or leptons, do not compute spallation. If protons, do not compute, because their spallation products will be computed elsewhere
    if (nuclei[iloop] <= 1001) continue;
    
    int iz = -1000;  // parent nucleus charge
    int ia = -1000;  // parent nucleus mass
    Utility::id_nuc(nuclei[iloop], ia, iz);
    
    for (int idaught = iloop+1; idaught < nuclei.size(); ++idaught) {
      
      int jz = -1000; // daughter nucleus charge
      int ja = -1000; // daughter nucleus mass
      Utility::id_nuc(nuclei[idaught], ja, jz);
      
      std::pair <double, double> couple (nuclei[iloop], nuclei[idaught]);
      //std::cout <<  "Getting spallations " << couple.first << " to " << couple.second <<  std::endl;
      
      if (nuclei[idaught] < 1001 || ia < ja) continue;
      for (unsigned int j = 0; j < dimEn; j++) {
	//if (iz == 14 & ia == 32 & jz == 14 & ja == 30){
	//std::cout <<  "Getting spallations " << spall[couple][j]/(beta[j]*factor) << "  " << spall[couple][j] << std::endl;
	//}
	outfileX << iz << "." << ia << "\t" << jz << "." << ja << "\t\t" << DRAGONEnergyVector[j] << "\t\t" << spall[couple][j]/(beta[j]*factor) << "\n";
      }
      outfileX << "\n\n";
    }
  }
  outfileX.close();
  
}

  
  std::vector<std::pair<double, double> > TSpallationNetwork::GetXSecPair(int parentParticleUID, int daughterParticleUID) {
    
    //std::cout <<  "Getting spallation cross sections from " << parentParticleUID << " to " << daughterParticleUID << std::endl;
    std::pair<int,int> inputParticlePair(parentParticleUID,daughterParticleUID);
    std::map< std::pair<int,int>, std::vector<std::pair<double, double> > >::iterator it = spallationXsectionsInterpolated.find(inputParticlePair);
    if (it != spallationXsectionsInterpolated.end()) 
      return (*it).second;
    else{
      //std::cout << "No spallation cross section for this pair " << std::endl;
      //std::pair<double,double> zeroPair(0.,0.);
      return std::vector<std::pair<double,double> >();
    }
    
  }
  
std::pair<double,double> TSpallationNetwork::GetXSecPair(int parentParticleUID, int daughterParticleUID, double en) {
  
  std::pair<int,int> inputParticlePair(parentParticleUID,daughterParticleUID);

  if (en > DRAGONEnergyVector.back() || en < DRAGONEnergyVector.front()) {
    std::pair<double,double> zeroPair(0.,0.);
    std::cout << "Energy value over or below the energy of particles in DRAGON" << std::endl;
    return zeroPair;
  }
  
  else if (en > spallationKineticEnergyVector.back() || en < spallationKineticEnergyVector.front()) { //Extrapolation   
    int je = 0, jee;
    double t;
    if (en > spallationKineticEnergyVector.back()){ 
      je = spallationKineticEnergyVector.size()-1; 	  jee = je;
      t = (en - spallationKineticEnergyVector[jee])/(spallationKineticEnergyVector[je]-spallationKineticEnergyVector[je-1]);
      
      if (en > 7.*spallationKineticEnergyVector.back())
	t = (7.*spallationKineticEnergyVector.back() - spallationKineticEnergyVector[jee])/(spallationKineticEnergyVector[je]-spallationKineticEnergyVector[je-1]);
    }
    
    else { // if (en < spallationKineticEnergyVector.front()){
      je = 1; jee = 0;
      t = (en - spallationKineticEnergyVector[jee])/(spallationKineticEnergyVector[je]-spallationKineticEnergyVector[je-1]);
      
      if (spallationKineticEnergyVector.front()/7. > en)
	t = (spallationKineticEnergyVector[0]/7. - spallationKineticEnergyVector[jee])/(spallationKineticEnergyVector[je]-spallationKineticEnergyVector[je-1]);
    }
    
    double yH  = spallationXsections[inputParticlePair][jee].first + t*(spallationXsections[inputParticlePair][je].first - spallationXsections[inputParticlePair][je-1].first);
    double yHe = spallationXsections[inputParticlePair][jee].second + t*(spallationXsections[inputParticlePair][je].second - spallationXsections[inputParticlePair][je-1].second);
    
    //std::cout << "Extrapolation! " << std::endl;
    
    std::pair <double, double> Xsecs(std::max(yH, 0.), std::max(yHe, 0.));
    return Xsecs;

  } //End of extrapolation block 
  
    //std::cout << "Getting spallation cross sections from " << parentParticleUID << " to " << daughterParticleUID << " --> ";
  std::map< std::pair<int,int>,  std::vector<std::pair<double, double> > >::iterator it = spallationXsections.find(inputParticlePair);
  if (it != spallationXsections.end()) {
    
      int je = 0, jee;
      
      while (en > spallationKineticEnergyVector[je])
	je++;
      
      if ((spallationKineticEnergyVector[je] - en) > (en - spallationKineticEnergyVector[je-1])) jee = je-1; 
      else
	jee = je;

      /*
	double t =  (en - spallationKineticEnergyVector[jee])/(spallationKineticEnergyVector[je]-spallationKineticEnergyVector[je-1]);
	
	double xsecHint  = spallationXsections[inputParticlePair][jee].first + t*(-spallationXsections[inputParticlePair][je].first + spallationXsections[inputParticlePair][je-1].first);
	double xsecHeint = spallationXsections[inputParticlePair][jee].second + t*(-spallationXsections[inputParticlePair][je].second + spallationXsections[inputParticlePair][je-1].second);
      */
      
      double XsecHarr[spallationKineticEnergyVector.size()], XsecHearr[spallationKineticEnergyVector.size()], theEn[spallationKineticEnergyVector.size()];
      //theEn[0] = 0.;
      
      for(int dan=0; dan < spallationKineticEnergyVector.size(); dan++){
	XsecHarr[dan] = spallationXsections[inputParticlePair][dan].first;
	//if (theEn[0] == 0.)
	theEn[dan] = spallationKineticEnergyVector[dan];
	XsecHearr[dan] = (spallationXsections[inputParticlePair][dan].second);
      }
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, spallationKineticEnergyVector.size());
      
      gsl_spline_init (spline, theEn, XsecHarr, spallationKineticEnergyVector.size()); //x axis must be always increasing!!
      double yH = gsl_spline_eval (spline, en, acc);
      
      
      gsl_spline_init (spline, theEn, XsecHearr, spallationKineticEnergyVector.size()); //x axis must be always increasing!!
      double yHe = gsl_spline_eval (spline, en, acc);
      
      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
      
      //std::cout << "\nClosest energy to "<< en << " is my vector energy = " <<  spallationKineticEnergyVector[jee] << " MeV with a Xsec value of: " <<spallationXsections[inputParticlePair][jee].first << " and the former was: " << spallationXsections[inputParticlePair][jee-1].first << std::endl;
      
      if (spallationXsections[inputParticlePair][jee].first == 0. && spallationXsections[inputParticlePair][jee-1].first == 0.)
	yH = 0.;
      if (spallationXsections[inputParticlePair][jee].second == 0. && spallationXsections[inputParticlePair][jee-1].second == 0.)
	yHe = 0.;	 
      
      //std::cout << "...interpolation done\n\n" << std::endl << std::endl;
      
      /*if (yH < 0.){std::cout << "bad gsl spline H" << yH << std::endl;
	//std::cout << "\nClosest energy to "<< en << " is my vector energy = " <<  spallationKineticEnergyVector[jee] << " MeV with a Xsec value of: " <<spallationXsections[inputParticlePair][jee].first << " and the former was: " << spallationXsections[inputParticlePair][jee-1].first << std::endl;
	yH = 0.;}
      
      if (yHe < 0.){ yHe = 0.;
	std::cout << "bad gsl spline He" << yHe << std::endl;}
      */
      
      yH = std::max(yH, 0.);
      yHe = std::max(yHe, 0.);
      std::pair <double, double> Xsecs(yH, yHe);
      return Xsecs;
    }
    
    else{
      std::cout << "No spallation cross section for this pair " << std::endl;
      std::pair<double,double> zeroPair(0.,0.);
      return zeroPair;
    }





  }
  
  std::pair<double,double> TSpallationNetwork::GetXSecPair(int parentParticleUID, int daughterParticleUID, int k) {
    //Daniele20.01.2018
    
    //std::cout << "Getting spallation cross sections from " << parentParticleUID << " to " << daughterParticleUID << "for the energy: " << DRAGONEnergyVector[k] << std::endl;
    std::pair<int,int> inputParticlePair(parentParticleUID,daughterParticleUID);
    std::map< std::pair<int,int>,   std::vector<std::pair<double, double> > >::iterator it = spallationXsectionsInterpolated.find(inputParticlePair);
    if (it != spallationXsectionsInterpolated.end()) 
      return (*it).second[k];
    else{
      std::cout << "No spallation cross section for this pair " << std::endl;
      std::pair<double,double> zeroPair(0.,0.);
      return zeroPair;
    }	
  }
  


  
//*********************************************************************************************************************************************
//*********************************************************************************************************************************************


  

// Inelastic cross sections
  TInelasticCrossSection::TInelasticCrossSection(TGrid* Coord, Input* in, int uid, int K_electron, std::vector<TXSecBase*> xsecmodel) { //double* coord, std::string* file, int uid, int K_electron){

    //std::vector<double> e_vec;
    e_vec = Coord->GetEk();
    const int dimEn = e_vec.size();
    vector<double> beta_vec = Coord->GetBetaEl();
    vector<double> gamma_vec = Coord->GetGammaEl();
    vector<double> beta = Coord->GetBeta();
    //std::vector<double> w_vec;

    const double factor = 1e-27*Clight;
    TGalpropXSec default_Ineobject(Coord, in);
  
   
    if (in->feedback > 0) 
    std::cout << "[debug] xsec test " << std::endl;    
        
    
    if (in->spallationxsec == Fluka){
      
      std::cout << "Fluka inelastic XSecs" << std::endl;
      //std::string inf = "data/Fluka2020_Inelastic.dat"; //FlukaIneXsec_new.dat";
      std::string inf = "data/Fluka2021_Inelastic.dat"; //FlukaIneXsec_new.dat";
      dataFileIne.open(inf.c_str());
      
      int projectile, target;//, ix = 0;
      std::vector<double> table_energy_vec;
      double energy_, xsec_value;
      
      std::vector<double> xsec_H, xsec_He;
      
      // reads the datafile
      while (dataFileIne >> projectile >> target >> energy_ >> xsec_value) {
	//std::cout << " \n" << projectile << "  " << target << "  " << energy_ << "  " << xsec_value << std::endl;
	if (projectile == uid) {
	  //std::cout << " \n" << projectile << std::endl;

	  //if (projectile == 5011 || projectile == 5010){ xsec_value = 1.1*xsec_value;}
	  //if (projectile == 3007 || projectile == 3006){ xsec_value = 1.1*xsec_value;}
	  //if (projectile == 4010 || projectile == 4009 || projectile == 4007){ xsec_value = 1.1*xsec_value;} 

	  if (target == 1001)
	    table_energy_vec.push_back(energy_); //The energy read is in GeV
	  
	  if (target == 1001)
	    xsec_H.push_back(xsec_value);  //Inelastic XSecs in mb
	  
	  else if(target == 2004)
	    xsec_He.push_back(xsec_value);
	  
	}
	else continue;
	
      } //End while loop to read inel
      
      if (xsec_H.empty() || xsec_He.empty())  { //If the provided table does not have IneXSec for this projectile...
	std::cout << "\nNo inelasctic XSec for this projectile from Fluka!! --> Taking the default one (CROSEC) \n" << std::endl;
	std::cout << "CROSEC inelastic XSecs" << std::endl;
	totalxsec = default_Ineobject.GetXSec(uid, K_electron);
      }
      
      
      else {
	std::vector<double> IneXSec; //Table of summed inelastic XSecs
	//std::pair<int, int> thecouple(projectile, target);
	
	for(int ip = 0; ip < table_energy_vec.size(); ip++)
	  IneXSec.push_back(xsec_H[ip] + He_abundance*xsec_He[ip]);  //We create the total vector in the file energy points
	
	
	for(int ip = 0; ip < e_vec.size(); ip++){
	  totalxsec.push_back(factor*beta[ip]*TInelasticCrossSection::GetTotalIne(table_energy_vec, IneXSec, e_vec[ip])); //interpolation to have vector at DRAGON energy points
	}
	
	std::cout << "\nFluka inelastic XSec OK!!\n" << std::endl;

      }

      
    } //End Fluka inelastic
    

    else if (in->spallationxsec == Webber03 || in->spallationxsec == DRAGON2){
      
      std::cout << "CROSEC inelastic XSecs" << std::endl;
      totalxsec = default_Ineobject.GetXSec(uid, K_electron);
      //std::cout << "\n\nWhat? " << totalxsec[3] << std::endl;
      
      /*
      int jz = -1000;  // parent nucleus charge
      int ja = -1000;  // parent nucleus mass
      Utility::id_nuc(uid, ja, jz);
      //if (outfileX.is_open() == false) {
      if (uid==4009){ //13027, -999
	outfileX.open("Ine_pBe_XS.dat");
	outfileX << "Zi.Ai\t" << "Energy(GeV)\t\t" << "XSec(mb)\n";
	for(int ip = 0; ip < e_vec.size(); ip++){      
	  //if (uid < 1001) continue;
	  //if (uid < 1001 || ja < jz) continue;
	  outfileX << jz << "." << ja << "\t\t" << e_vec[ip] << "\t\t" <<  totalxsec[ip]/beta[ip]/factor << "\n"; 
	}
	outfileX.close();
      }
      */
      
      
    } //End Webber or Carmelo inelastic
    
    
    else {
      
      if  (in->spallationxsec != GalpropXSec)
	std::cout << "Wrong SpallationXSec option!! --> Taking the default one (CROSEC) \n" << std::endl;
      
      //ineXsecFile = "GalpropIneXsec_new.dat";
      std::cout << "\nCROSEC inelastic XSecs" << std::endl;
      totalxsec = default_Ineobject.GetXSec(uid, K_electron);
	//std::cout << "\n\nWhat? " << totalxsec[3] << std::endl;

    }
        
    //Positron Annihilation
    
    if (uid == 1000) {
      
      std::cout << std::endl << "*** Positron annihilation cross section" << std::endl;
      
      for(unsigned int ip = 0; ip < e_vec.size(); ip++) {
	
	double gamma = gamma_vec[ip];
	
	// cross section from P. A. M. Dirac, Proc. Camb. Phil. Soc. 26, 361 (1930). sigma = pi r_e^2 /(gamma+1) * [�some function of gamma��] 
	double annihilation_cross_section = 249.46/(gamma+1) * ( (gamma*gamma + 4*gamma + 1)/(gamma*gamma -1) * log(gamma + sqrt(gamma*gamma-1)) - (gamma+3)/(sqrt(gamma*gamma-1)) );  //mbarn
	
	//std::cout << factor*beta_vec[ip]*annihilation_cross_section << std::endl;
	
	totalxsec[ip] += (factor*beta_vec[ip])*annihilation_cross_section;	
	
      }
    }
    
  } //End TInelasticCrossSection constructor

  
  double TInelasticCrossSection::GetTotalIne (std::vector<double> evec_table, std::vector<double> IneXSec_, double en) {
    
    //std::cout << "Interpolating IneXSecs from Fluka tables..." << std::endl;
    
    if (en > evec_table.back() || en < evec_table.front()) { //Extrapolation   
    int je = 0, jee;
    double t;
    if (en > evec_table.back()){ 
      je = evec_table.size()-1; 	  jee = je;
      t = (en - evec_table[jee])/(evec_table[je]-evec_table[je-1]);
      
      if (en > 7.*evec_table.back())
	t = (7.*evec_table.back() - evec_table[jee])/(evec_table[je]-evec_table[je-1]);
    }
    
    else { // if (en < spallationKineticEnergyVector.front()){
      je = 1; jee = 0;
      t = (en - evec_table[jee])/(evec_table[je]-evec_table[je-1]);
      
      if (evec_table.front()/7. > en)
	t = (evec_table[0]/7. - evec_table[jee])/(evec_table[je]-evec_table[je-1]);
    }
    
    double yH  = IneXSec_[jee] + t*(IneXSec_[je] - IneXSec_[je-1]);
    
    //std::cout << "Extrapolation! " << std::endl;
    return std::max(yH, 0.);
    
  } //End of extrapolation block 


    else{      //Interpolating...
      int je = 0, jee;
      
      while (en > evec_table[je])
	je++;
      
      if (( evec_table[je] - en) > (en - evec_table[je])) jee = je-1; 
      else
	jee = je;
      
      double IneXsecarr[evec_table.size()], theEn[evec_table.size()];  //gsl needs double[] as arguments, not vector...
      for(int dan=0; dan < evec_table.size(); dan++){
	IneXsecarr[dan] = IneXSec_[dan];
	theEn[dan] =  evec_table[dan];
      }
      
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, evec_table.size());
      
      gsl_spline_init (spline, theEn, IneXsecarr, evec_table.size()); //x axis must be always increasing!!
      double yH = gsl_spline_eval (spline, en, acc);
      
      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
      
      //std::cout << "\nClosest energy to "<< en << " is my vector energy = " <<  spallationKineticEnergyVector[jee] << " MeV with a Xsec value of: " <<spallationXsections[inputParticlePair][jee].first << " and the former was: " << spallationXsections[inputParticlePair][jee-1].first << std::endl;
      
      return std::max(yH, 0.);
    }  
  }//End GetTotalIne method 
  

std::vector<double> TGalpropXSec::GetXSec(int uid, int K_electron) { //modified
  int Z = -1000;
  int A = -1000;
  Utility::id_nuc(uid, A, Z);

  std::vector<double> result = std::vector<double>(energy.size(), 0.0);
  
  if (A == 0) return result;
  
  int target = 1;
  const bool protons = (A == 1 && Z == 1);
  const bool antiprotons = (A == 1 && Z == -1);
  const bool antideuteron = (A == 2 && Z == -1);
  //const bool antitritium = (A == 3 && Z == -1);
  const bool antiHe3 = (A == 3 && Z == -2);
  const bool antiHe4 = (A == 4 && Z == -2);

  if (protons) {
    A = 4;
    Z = 2;
  }
  /*
    if (antiHe3) {
    A = 4;
    Z = 2;
    target = -1;
  }
  */
    
  if (antiprotons || antiHe4 || antiHe3 || antideuteron) {
    A = 4;
    Z = 2;
    target = -1;
  }
  
  double PP_inel = 0.0;
  double PA_inel = 0.0;
  double aPP_non = 0.0;
  double aPA_non = 0.0;
  double aPP_ann = 0.0;
  double aPA_ann = 0.0;
  
  const double factor = 1.e-27*Clight;
  
  const double Hecontribution = He_abundance*TGalpropXSec::He_to_H_CS_totratio(A);
  
  for (unsigned int k = 0; k < energy.size(); ++k) {     
    TGalpropXSec::nucleon_cs(2, energy[k], target, Z, A, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &aPA_ann);  // Inelastic CROSEC
    //std::cout << "\n\n\n\n Here it is!!   " << PA_inel << std::endl;

    if (protons) result[k] = (PP_inel + He_abundance*PA_inel)*factor*beta[k];
    else if (antiprotons) result[k] = factor*(aPP_non + aPP_ann + He_abundance*(aPA_non + aPA_ann))*beta[k];
    else if (antideuteron) result[k] = pow(2, 0.8)*factor*(aPP_non + aPP_ann + He_abundance*(aPA_non + aPA_ann))*beta[k];
    else if (antiHe3) result[k] = pow(3, 0.8)*factor*(aPP_non + aPP_ann + He_abundance*(aPA_non + aPA_ann))*beta[k];
    //else if (antiHe4) result[k] = 
    else result[k] = PA_inel*factor*beta[k]*(1.0 + Hecontribution);
    
    //if (uid == -999){  std::cout << "\n\n\n This is Antiproton" << uid << "\n\n\n" << std::endl; }
    //if (uid == -998){  std::cout << "\n\n\n This is Anti-D" << uid << "\n\n\n" << std::endl; }
    //  result[k] *= 1.1;}
    //if (uid == 6012){  //std::cout << "\n\n\n This is what I wanted" << uid << "\n\n\n" << std::endl;
    //  result[k] *= 1.1;}
    
    

    /*
    //modified
    if (K_electron >= 0) {
    
    double attach_H=0., attach_He=0., strip_H=0., strip_He=0.;
    
    //modified
    TGalpropXSec::Kcapture_cs(energy[k], Z,1,&attach_H,&strip_H);
    TGalpropXSec::Kcapture_cs(energy[k], Z,2,&attach_He,&strip_He);
    
    //modified
    if (K_electron == 0) 
    result[k] += ( factor*beta[k]*(attach_H + He_abundance*attach_He) ); //The naked nucleus unstable for EC may attach an electron
    if (K_electron == 1)
    result[k] += ( factor*beta[k]*(strip_H + He_abundance*strip_He) );   //The "dressed" nucleus unstable for EC may lose an electron and become naked again 
    
    }
    */  
  }
  return result;
  
}
    
    
    
  
  void TGalpropXSec::Kcapture_cs(double Ek_GeV, int Zp, int Zt, double *attach, double *strip) {
    
    double gam = 1.+Ek_GeV/amu, sT =1.e27 * 8./3.*Pi*pow(Rele,2), //AWS20010829
      beta,Mbet,Nbet,a,fcor;
    
    beta = sqrt(1.-pow(gam,-2));
    Mbet = 4./3.+gam*(gam-2.)/(gam+1.)*(1.-log((1.+beta)/(1.-beta))/2./beta/pow(gam,2));
    Nbet = pow(beta,-3) *((-4.*gam +34. -63./gam +25./pow(gam,2) +8./pow(gam,3))/15.
			  -(gam-2.)*(gam-1.)/2./beta/pow(gam,3) *log((1.+beta)/(1.-beta)));
    a = Zp*ALPHAf;
    fcor = pow(a,2.*(sqrt(1.-a*a)-1.)) *exp(-2.*a/beta*acos(a)) *(1.+Pi*a*Nbet/Mbet);
    //    factor 1.202 accounts for contribution of states higher than 1s
    *attach = 1.202 *fcor *3./2.*sT*pow(1.*Zp,5)*Zt*pow(ALPHAf,4) *beta*gam *pow(gam-1.,-3)*Mbet;
    *strip = 3./2.*sT*pow(Zp*beta*ALPHAf,-2) 
      *0.285*(2.*log(2./0.219*beta*gam/Zp/ALPHAf)-pow(beta,2))*Zt*(Zt+1.);
    
  }
  
  
  double TGalpropXSec::He_to_H_CS_ratio(double E1,int IZI,int IAI,int IZF,int IAF) {
    double X,Y,a,b;
    double AMU,DELTA,FZI;
    double  E[]= {-999,  0.43,0.73,1.51};  // zeroth element dummy to conform to F77 original
    double d[]= {-999,  2.45,3.45,4.40};
    double am[]= {-999,  0.082,0.047,0.031};
    double Z[]= {-999,  6.,8.,26.};
    double F[]= {-999, -0.76,-0.41,+1.00};


    double CSratio    = 0.;
    if(IAI <= IAF) return 0;
    
    if(E1 < E[2])          // interpolation/extrapolation
      {
  	AMU  = FI(E1,E[1],E[2],am[1],am[2]);
  	DELTA= FI( log(E1), log(E[1]), log(E[2]),d[1],d[2]);  //log scale
      } else
      {
  	AMU  = FI(E1,E[2],E[3],am[2],am[3]);
  	DELTA= FI( log(E1), log(E[2]), log(E[3]),d[2],d[3]);  //log scale
      }    
    if(E1 < E[1])          // asymptotics E1->0
      {
  	AMU  = am[1];
  	DELTA=  d[1];
      }             
    if(E1 > E[3])          // asymptotics E1->inf
      {
  	AMU  = am[3];
  	DELTA=  d[3];
      }             
   
    if(IZI < Z[2]) FZI = FI( log(IZI*1.), log(Z[1]), log(Z[2]),F[1],F[2]); // inter-/extra-polation on log scale
    else           FZI = FI( log(IZI*1.), log(Z[2]), log(Z[3]),F[2],F[3]);   
    
    if(IZI-IZF <= IZI/2) CSratio = exp(AMU*pow(fabs(IZI-IZF-FZI*DELTA),1.43));
    else                 // if (IZI-IZF) > IZI/2, then another approximation (requires continuity in Z)
      {
  	X =   exp(AMU*pow(fabs(IZI/2-  FZI*DELTA), 1.43));      //value
  	Y = X-exp(AMU*pow(fabs(IZI/2-1-FZI*DELTA), 1.43));      // derivative in dZ
  	b = IZI/2*Y/X;
  	a = X/pow(IZI/2,b);
  	CSratio = a*pow(IZI-IZF,b);
      }
    return CSratio;
  }


void TGalpropXSec::nucleon_cs(int option, double Ek, int Zp, int Zt, int At, double *PP_inel,double *PA_inel,double *aPP_non,double *aPA_non,double *aPP_ann,double *aPA_ann) {
    
  *PP_inel= *PA_inel= *aPP_non= *aPA_non= *aPP_ann= *aPA_ann=0.;
  if(Ek <= 0.) return;
    
  double U,Cp,C1,s2,Z,A, aPP_inel=0.,aPA_inel=0., Emev=Ek, b0,Fcorr,rN,s0,p1,p2,p3,p4,p5,x,f1,f2;
  double PZ=fabs(1.*Zp),Em=1000.*Ek, TZ=Zt, TA=At;
  int ISS=2;
    
  //## Proton-Proton INELASTIC cross section, mb [TN83,p.234]
  if(Ek > 0.3)
    {
      U = log((Ek+mp)/200.);
      *PP_inel = 32.2*(1.+0.0273*U);
      if(U >= 0.) *PP_inel += 32.2*0.01*U*U;
      if(Ek < 3.) *PP_inel /= 1.+0.00262/pow(Ek,17.9+13.8*log(Ek)+4.41*pow(log(Ek),2));
    }

  //P-P parametrization from Kafexhiu, Aharonian, A.Taylor and Gabriela Vila (2014)
  double Eth = 0.2797; //GeV  = 2*mpi + mpi**2/2mp
  *PP_inel = std::max((30.7 - 0.96*log(Ek/Eth) + 0.18*pow(log(Ek/Eth),2))*pow((1-pow((Eth/Ek), 1.9)), 3), 0.); //mb 
  if(Zp*At == 1) return;
    
  //## Proton-Nucleus INELASTIC cross section, mb
  switch(option) {
  case 0: // [L83]
    C1 = (At == 4) ? 0.8 : 1.;                  // a correction for He
    if(At == 9) C1 = 1.+0.75*exp(-Ek/0.075);    // a correction for Be
    *PA_inel = C1 *45. *pow(TA,0.7) *(1.+0.016*sin(5.3-2.63*log(TA)));
    if(Ek < 5.) *PA_inel *= 1.-0.62 *exp(-Ek/0.2) *sin(10.9/pow(Ek*1.e3,0.28));
    if(At == 4) *PA_inel = (Ek > 0.01) ?        // pHe, my fit
		  111.*(1.-(1.-sin(9.72*pow(log10(Ek*1000.),0.319)-4.14))*exp(-3.84*(Ek-0.1))) : 0.;
    break;
            
  case 1: // Zt>5 [WA96], Zt<=5 [BP01]
    if(Zt>5) {
      b0 = 2.247-0.915*(1.-pow(TA,-1./3.));
      Fcorr = (1.-0.15*exp(-Emev))/(1.-0.0007*At); // high-energy correction
      rN = (At-Zt>1.5) ? log(TA-Zt) : 1.;
      s0 = Pi*10.*pow(1.36,2.)*Fcorr*rN*(1.+pow(TA,1./3.)-b0*(1.-pow(TA,-1./3.)));
      p1 = 8.-8./At-0.008*At;
      p2 = 2.*(1.17-2.7/At-0.0014*At);
      p3 = 0.8+18./At-0.002*At;
      p4 = 5.6-0.016*At;
      p5 = 1.37*(1.+1./At);
      x = log10(Emev);
      f1 = 1./(1.+exp(-p1*(x+p2))); // low-energy return to zero
      f2 = 1. +p3 *( 1. -1./(1.+exp(-p4*(x+p5))) ); // low-energy rise
      *PA_inel = f1*f2*s0;
    }
    break;
            
  case 2: // [BP01]
  default:
    if (Em<14.) Em=14.;
    if (Em>1.e6) Em=1.e6;
    *PA_inel = sighad_cc(ISS,PZ,PZ,TA,TZ,Em); // IMOS20020502
  }
  if(Zp*At >= 1) return;
    
  //## AntiProton-Proton ANNIHILATION cross section [TN83]
  if(Ek < 10.) *aPP_ann = 661.*(1.+0.0115/pow(Ek,0.774)-0.948*pow(Ek,0.0151)); // 0.1GeV<Ek<12GeV
  else
    {
      // assuming aPP_ann = aPP_tot -PP_tot (i.e., aPP_elast = PP_elast); (aPP_tot-PP_tot) from [PDG00]
      s2 = 2.*mp*(Ek+2*mp);                   // square of the total CMS energy
      *aPP_ann = 2*35.43/pow(s2,0.560);
    }
    
  //## AntiProton-Proton TOTAL INELASTIC cross section
  aPP_inel = *PP_inel + *aPP_ann;
  if(Ek <= 14.)
    { 
      aPP_inel = 24.7*(1.+0.584/pow(Ek,0.115)+0.856/pow(Ek,0.566));
      if(*aPP_ann > aPP_inel) *aPP_ann = aPP_inel;
    }
    
  //## AntiProton-Proton TOTAL INELASTIC NON-ANNIHILATION cross section
  *aPP_non = aPP_inel - *aPP_ann;
  if(*aPP_non < 0.) *aPP_non = 0.;
    
  //## AntiProton-NUCLEUS cross sections
  if(At > 1)
    {
      //## AntiProton-NUCLEUS TOTAL INELASTIC NON-ANNIHILATION cross section
      *aPA_non = *PA_inel;
        
      //## AntiProton-NUCLEUS ANNIHILATION cross section on 12C-nucleus [mb] (0.4<Pp<300) [MO97]
      A = At;
      Z = Zt;                               // Z = 0.59*pow(A,.927);  for Z > 2 nuclei
      if(At == 4) { Z = 2.; A = 3.30; }     // Modified to agree with HE p-He cs / imos
      *aPA_ann = pow(A,2./3.)               // Scaling to other nuclei
        //         *(48.2 +19./pow(Ek-0.02,0.55) 
        *(48.2 +19./pow(Ek,0.55)           // modified to agree w. He@<100 MeV / imos
          +(0.1-0.18/pow(Ek,1.2))*Z +0.0012/pow(Ek,1.5)*Z*Z)  - *aPA_non;
      if(*aPA_ann < 0.) *aPA_ann = 0.;
      if(*aPA_ann < *aPP_ann) *aPA_ann = *aPP_ann;
      if(At == 4 && Ek > 5.)  *aPA_ann = *aPP_ann;
      /*
      //## my fit to AntiProton-NUCLEUS total cross section on 12C-nucleus [mb] (0.4<Pp<300)
      double Pp =sqrt(pow(Ek+mp,2)-mp*mp); // GeV,kin. momentum per nucleon
      *aPA_ann = (Pp > 40.) ? 236.*(1.+6.9e-2*exp(-Pp/100.)) :
      209.7*pow(.542/Pp,.565)+29.6*cos(log10(Pp/9.29)*5.11)+257.9;
      *aPA_ann *= pow(At/12.,2./3.);                             // scaling to other nuclei
      */
    }
  return;
}



  double TGalpropXSec::antiproton_cc1(gsl_integration_workspace* w, size_t limit, int key, double Pap, double Pp1, int NZ1, int NA1, int NZ2, int NA2) {
    
    double Pp = Pp1/double(NA1);
    
    double s2 = 2.*mp*( sqrt( pow(Pp,2.) + pow(mp,2.) ) + mp );
    double s1 = sqrt(s2); // Center of Mass total energy
    
    if (s1 <= 4.*mp) 
      return 0.;
    
    double MX = 3.*mp;
    
    double gamma_CM = s1/(2.*mp);  // Center of Mass Lorentz factor
    double betagamma_CM = sqrt( pow(gamma_CM, 2.) - 1.  ); // Center of Mass beta*gamma
    double Eap = sqrt ( pow(Pap,2.) + mp*mp );
    
    if (NA2 > 1.) {
      Eap = Eap + 0.06;
      Pap = sqrt( Eap*Eap - mp*mp);
    }
    
    double Eap_MAX = (s2 - MX*MX + mp*mp)/(2.*s1);
    double Ek = sqrt(Pp*Pp + mp*mp)-mp; 
    
    double result = 0.;
    double cosX_MAX = -1.;
    
    if (betagamma_CM*Pap > 0)
      cosX_MAX = (gamma_CM*Eap - Eap_MAX)/(betagamma_CM*Pap); 
    
    if (cosX_MAX < -1.)
      cosX_MAX = -1.;
    
    double param[2] = {Pap, Pp};
    
    if (cosX_MAX <= 1.0) {
      
      //size_t limit = 1000;
      
      double AI = 0;
      double error = 0;
      gsl_function F;
      F.function = &(tan_ng);
      F.params = param;
      gsl_integration_qag(&F, 1.0, cosX_MAX, 0, 1e-4, limit, 6, w, &AI, &error);
      // dF/dEap - production spectrum vs. energy; factor 2 accounts for antineutrons
      result = -AI * 2.0 * Pap  * (2. * Pi) * 1.e-3 * Pap / Eap  ;   // conversion from mbarn c^3/GeV^2 to mbarn/GeV
      //   gsl_integration_workspace_free(w);
      
    }
    if (NA1 * NA2 == 1 || Ek <= 0.3) return result;
    
    // end of p p -> ap X computation
    
    double PP_inel = 0.;
    double PA_inel = 0.;
    double aPP_non = 0.;
    double aPA_non = 0.;
    double aPP_ann = 0.;
    double apA_ann = 0.;
    nucleon_cs(key, Ek, 1, NZ2, NA2, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &apA_ann);
    double CS2 = PP_inel;
    if (NA2 > 1) CS2 = PA_inel;
    double CS1 = PP_inel;
    if (NA1 > 1) {
      nucleon_cs(key, Ek, 1, NZ1, NA1, &PP_inel, &PA_inel, &aPP_non, &aPA_non, &aPP_ann, &apA_ann);
      CS1 = PA_inel;
    }
    
    //  multiplicity of anti-p in pA-,Ap-,AA-coll.; Gaisser&Schaeffer 1992,ApJ,394,174
    double correction = 1.2 * (double(NA1)*CS2 + double(NA2)*CS1) / (2.*PP_inel);
    return result*correction;
  }
  
  double tan_ng(double cosX, void* param) {
    
    // This routine computes the inclusive antiproton production cross section
    // result is E d^3sigma/dp^3 expressed in mbarn c^3/GeV^2
    // antiprotons come from the interaction between Cosmic Rays protons and
    // nuclei and interstellar gas
    
    double* param1 = (double*)param;
    double Pap = param1[0];
    double Pp = param1[1];
    
    
    double s2 = 2.*mp*( sqrt(Pp*Pp+mp*mp) + mp );
    double s1 = sqrt(s2); // Center of Mass total energy
    
    if (s1 < 4.*mp) // if CM energy < 4 GeV, no antiprotons are produced (reaction threshold)
      return 0;
    
    double gamma_CM = s1/(2.*mp);
    double betagamma_CM = sqrt( pow(gamma_CM, 2.) - 1.  ); // Center of Mass beta*gamma
    double Eap = sqrt(Pap*Pap + mp*mp);
    double sinX = 0.;
    if (1.-pow(cosX,2) > 0)
      sinX = sqrt(1.-pow(cosX,2));
    double pt = Pap*sinX; // anti-proton transverse momentum
    double MX = 3.*mp;
    
    double Xr = 2.*s1*( gamma_CM*Eap - betagamma_CM*Pap*cosX )/( s2 - MX*MX + mp*mp  )   ;  
    
    if (Xr > 1)
      return 0;
    
    //Cross section calculation. See Tan&Ng 1983
    
    double A = 0.465 * exp(-0.037*Xr) + 2.31 * exp(0.014*Xr); // c/GeV
    double B = 0.0302 * exp(-3.19*(Xr + 0.399))*pow(Xr+0.399, 8.39); 
    //(c/GeV)^2
    double f = (3.15-1.05e-4)*pow((1.-Xr), 7.9);
    if (Xr <= 0.5)
      f += 1.05e-4*exp(-10.1*Xr);
    double result = f*exp(-A*pt + B*pt*pt);
    
    if (s1 < 10.) {
      
      double Xt = ((gamma_CM*Eap-betagamma_CM*Pap*cosX)-mp)/( (s2 - MX*MX +
							       mp*mp) / (2.*s1)-mp); 
      double a = 0.306*exp(-0.12*Xt);
      double b = 0.0552*exp(2.72*Xt);
      double c = 0.758 - 0.68*Xt + 1.54*pow(Xt,2.);
      double d = 0.594*exp(2.87*Xt);
      double Q1 = s1 - 4.*mp;
      double correction = (1. - exp(-exp(c*Q1-d)*(1.-exp(-a*pow(Q1,b))) ));
      result /= correction;
      
    }
    
    return result;
    
  }
  
  
void TGalpropXSec::set_sigma_cc() { 
  int  cdr=99;
  set_sigma_(&cdr);
}

void TGalpropXSec::sigtap_cc(int ISS) { 
  sigtap2_(&ISS);
}

double TGalpropXSec::wsigma_cc(int IZ, int IA, int IZF, int IAF, double E) {
  return( wsigma_(&IZ,&IA,&IZF,&IAF,&E) );
}

double TGalpropXSec::yieldx_cc(int IZ, int IA, int IZF, int IAF, float E) {
  float CSmb;
  yieldx_(&IZ,&IA,&IZF,&IAF,&E,&CSmb);
  return( 1.*CSmb );
}

double TGalpropXSec::sighad_cc(int IS, double PA, double PZ, double TA, double TZ, double E) { 

  return( sighad_(&IS, &PA, &PZ, &TA, &TZ, &E) );
}

  
//} //End DRAGON namespace
