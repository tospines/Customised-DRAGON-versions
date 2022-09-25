#ifndef _DIFFUSION_H
#define _DIFFUSION_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include "utilities.h"
#include "constants.h"

#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>

using namespace std;

class TGrid;
class Input;
class TBField;
class TGeometry;
class TSource;

/**
 * @class TDiffusionCoefficient
 * @author Luca Maccione
 * @author Daniele Gaggero
 * @email luca.maccione@lmu.de
 * @email daniele.gaggero@sissa.it
 * @brief Position dependent diffusion coefficient.
 */
class TDiffusionCoefficient {
    
public:
    TDiffusionCoefficient() { }
    TDiffusionCoefficient(TGrid* /**< Geometry. */, Input*, TSource*, TBField*, std::vector<int>&, bool El = false /**< Electrons/positrons need special treatment*/);
    virtual ~TDiffusionCoefficient() { }
    /**< Destructor. */
    
    // 2D
    
    inline vector<double>& GetDiffusionCoefficient() { return dperp;  cout << "[MW-DEBUG] called GetDiffusionCoefficient(i)"<<endl; }
    inline double GetDiffusionCoefficient(int i /**< Linearized position index. */) { return dperp[i]; }  
    
    
    inline vector<vector<double> >& GetSpectrum() { return sp;  cout << "[MW-DEBUG] called GetSpectrum()"<<endl; }

    inline vector<double>& GetSpectrum(int part /**<index of the particle*/) { return sp[part];  cout << "[MW-DEBUG] called GetSpectrum()"<<endl; }
        
    inline double GetDiffusionCoefficient(int isp /**< Linearized position index. */, int ip /**< Energy index. */, int part /**<index of the particle*/) { return dperp[isp]*sp[part][ip];  cout << "[MW-DEBUG] called GetDiffusionCoefficient(isp,ip)"<<endl; }
    
    inline double GetDiffusionCoefficient(int ir /**< Radial index. */, int iz /**< Vertical index. */, int ip /**< Energy index. */, int part /**<code of the particle*/) { return dperp[index(ir,iz)]*sp[part][ip];  cout << "[MW-DEBUG] called GetDiffusionCoefficient(ir,iz,ip)"<<endl; }

    inline double GetSpectrum(int i /**< Energy index. */, int part /**<code of the particle*/) { return sp[part][i]; cout << "[MW-DEBUG] called GetSpectrum(i)"<<endl; }
    
    std::vector<double> DiffmomVector, DiffmomVector_H, DiffmomVector_ED; /**< momentum array of Diffusion table */
    std::vector<double> Diff_vect, Diff_vect_H, Diff_vect_ED; /**< diffusion vector of Diffusion table */
    std::vector<double> DRAGONmomentum;

    std::string DiffFile_name, DiffFile_nameHalo, DiffFile_nameEDisc;
    
    inline double GetDelta() { return delta; }
    /**< Retrieve delta. */
    inline double GetDelta_h() { return delta_H; }
    /**< Retrieve delta_H. */
    inline double Getrho_b() { return rho_H; }
    /**< Retrieve rho_H. */
   
    //SK130830
    inline const vector<double>& GetPhi() const { return phi; } 
    inline double GetPhi(int i /**< Linearized position index. */) { return phi[i]; }

    TGrid* GetCoord() { return _fCoordinates; cout << "[MW-DEBUG] my new sweet Diffusion GetCoord() method is called! yay!" << endl;}
    /**< Obtain the Geometry description. */
    
    // 3D
    inline double GetPhiX(int i) { return phix[i]; }
    inline double GetPhiY(int i) { return phiy[i]; }
    inline double GetPhiZ(int i) { return phiz[i]; }
    
    inline double GetPhiX(int ix,  int iy, int iz, int ip) { return phix[index(ix,iy,iz,ip)]; }
    inline double GetPhiY(int ix,  int iy, int iz, int ip) { return phiy[index(ix,iy,iz,ip)]; }
    inline double GetPhiZ(int ix,  int iy, int iz, int ip) { return phiz[index(ix,iy,iz,ip)]; }
    
    inline double GetPhiX(int isp, int ip) { return phix[isp*dimE+ip]; }
    inline double GetPhiY(int isp, int ip) { return phiy[isp*dimE+ip]; }
    inline double GetPhiZ(int isp, int ip) { return phiz[isp*dimE+ip]; }
    
    inline double GetAlpha_XX(int isp, int ip) { return alpha_xx[isp*dimE+ip]; }
    inline double GetAlpha_YY(int isp, int ip) { return alpha_yy[isp*dimE+ip]; }
    inline double GetAlpha_ZZ(int isp, int ip) { return alpha_zz[isp*dimE+ip]; }
    
    inline double GetDparTotal(int i)  { return dpar_total[i]; }
    inline double GetDperpTotal(int i) { return dperp_total[i]; }
    inline double GetAlpha_XY(int isp, int ip) { return alpha_xy[isp*dimE+ip]; }
    inline double GetAlpha_YZ(int isp, int ip) { return alpha_yz[isp*dimE+ip]; }
    inline double GetAlpha_XZ(int isp, int ip) { return alpha_xz[isp*dimE+ip]; }
    inline int GetDimR() const { return dimr; }
    inline int GetDimX() const { return dimx; }
    inline int GetDimY() const { return dimy; }
    inline int GetDimZ() const { return dimz; }
    inline int GetDimE() const { return dimE; }
    vector<double>& GetDPar() { return Dpar_profile; }

    //MW130624
    inline double GetCNdiff_alpha1_x(int indspat, int part) { return CNdiff_alpha1_x[part][indspat]; }
    inline double GetCNdiff_alpha2_x(int indspat, int part) { return CNdiff_alpha2_x[part][indspat]; }
    inline double GetCNdiff_alpha3_x(int indspat, int part) { return CNdiff_alpha3_x[part][indspat]; }
    inline double GetCNdiff_alpha1_y(int indspat, int part) { return CNdiff_alpha1_y[part][indspat]; }
    inline double GetCNdiff_alpha2_y(int indspat, int part) { return CNdiff_alpha2_y[part][indspat]; }
    inline double GetCNdiff_alpha3_y(int indspat, int part) { return CNdiff_alpha3_y[part][indspat]; }
    inline double GetCNdiff_alpha1_z(int indspat, int part) { return CNdiff_alpha1_z[part][indspat]; }
    inline double GetCNdiff_alpha2_z(int indspat, int part) { return CNdiff_alpha2_z[part][indspat]; }
    inline double GetCNdiff_alpha3_z(int indspat, int part) { return CNdiff_alpha3_z[part][indspat]; }
    //MW130705
    inline double GetCNdiff_alpha1_r(int indspat, int part) { return CNdiff_alpha1_r[part][indspat]; }
    inline double GetCNdiff_alpha2_r(int indspat, int part) { return CNdiff_alpha2_r[part][indspat]; }
    inline double GetCNdiff_alpha3_r(int indspat, int part) { return CNdiff_alpha3_r[part][indspat]; }
    

protected:
    
    TGrid* _fCoordinates; /**< Geometry and kinematics. */
    
    int dimr; /**< Radial dimension of simulation box. */
    int dimx; /**< Radial dimension of simulation box. */
    int dimy; /**< Radial dimension of simulation box. */
    int dimz; /**< Vertical dimension of simulation box. */
    int dimE;
    
    double zt; /**< Vertical scale of the diffusion coefficient. */
    DPerpType set_profile; /**< Spatial profile. \sa constants.h */
    DiffusionType diffT;
    double index_radial; /**< */
    double GetProfile(double, double, double, TSource*);
    /**< Compute the spatial profile. */
    double GetXDerivative(double, double, double, TSource*);
    double GetYDerivative(double, double, double, TSource*);
    double GetZDerivative(double, double, double, TSource*);

    double GetDInterp(double , double, string);
    std::ifstream datafileDiff, datafileDiff_H, datafileDiff_ED; 


    inline double GetRDerivative(double x, double y, double z, TSource* ST) { return GetXDerivative(x,y,z,ST); }
    //MW130708: fix boundaries in here to make diffusion.cc more readable!
    inline int index(int ir /**< Radial index. */, int iz /**< Vertical index. */)
    {
        if(ir<0) ir = 0;
        if(iz<0) iz = 0;
        if(ir>=dimr) ir = dimr - 1;
        if(iz>=dimz) iz = dimz - 1;
        
        return ir*dimz + iz;
    }
    inline int index(int ix, int iy, int iz) {
        if(ix<0) ix = 0;
        if(iy<0) iy = 0;
        if(iz<0) iz = 0;
        if(ix>=dimx) ix = dimx - 1;
        if(iy>=dimy) iy = dimy - 1;
        if(iz>=dimz) iz = dimz - 1;
        
        return (ix*dimy + iy)*dimz + iz;
    }
    
    inline int index_pzr(int ip, int iz, int ir) {
      
      if(ir<0) ir = 0;
      if(iz<0) iz = 0;
      if(ir>=dimr) ir = dimr - 1;
      if(iz>=dimz) iz = dimz - 1;      
      return (ip*dimz + iz)*dimr+ir;
    }
    


    inline int index(int ix /**< Radial index. */, int iy, int iz /**< Vertical index. */, int ip) { return index(ix,iy,iz)*dimE+ip; }
    /**
     * @fn inline int index(int ir, int iz)
     * @brief convert from matrix to linear representation.
     * @return ir*dimz+iz
     */
    
    // 2D
    vector<double> dperp;
    
//     vector<double> psi;
    vector<double> phi;
    
    vector<vector<double> > sp;
    vector<vector<double> > spectrum_extended;
    double delta; /**< Power-law index of the energy dependence of the diffusion coefficient below rho_H. */
    double delta_H, delta_L; /**< Power-law index of the energy dependence of the diffusion coefficient above rho_H. */
    double R_H, R_L, rho_H, rho_L, s_H, s_L, Deltadelta_H, Deltadelta_L;
    double nrn_sn;
    
    // 3D
    vector<double> dperp_total;
    vector<double> dpar_total;	

    vector<double> alpha_xx;
    vector<double> alpha_yy;
    vector<double> alpha_zz;
    
    vector<double> alpha_xy;
    vector<double> alpha_yz;
    vector<double> alpha_xz;
    
    vector<double> phix;
    vector<double> phiy;
    vector<double> phiz;

    vector<double> Dpar_profile;
    
    //MW130623: Trying to gain speed on the cost of memory...
    vector<vector<double> > CNdiff_alpha1_x;
    vector<vector<double> > CNdiff_alpha2_x;
    vector<vector<double> > CNdiff_alpha3_x;
    vector<vector<double> > CNdiff_alpha1_y;
    vector<vector<double> > CNdiff_alpha2_y;
    vector<vector<double> > CNdiff_alpha3_y;
    vector<vector<double> > CNdiff_alpha1_z;
    vector<vector<double> > CNdiff_alpha2_z;
    vector<vector<double> > CNdiff_alpha3_z;
    //MW130705
    vector<vector<double> > CNdiff_alpha1_r;
    vector<vector<double> > CNdiff_alpha2_r;
    vector<vector<double> > CNdiff_alpha3_r;

};

/**
 * @class TDiffusionCoefficient
 * @author Luca Maccione
 * @email luca.maccione@lmu.de
 * @brief Position dependent diffusion coefficient. Part of it independent of the nucleus. Needed coefficient will be restored in TCREvolutor.
 */
class TDiffusionCoefficient2D : public TDiffusionCoefficient {
    
public:
    TDiffusionCoefficient2D() {} /**< Default constructor. */
    TDiffusionCoefficient2D(TGrid* /**< Geometry. */, Input*, TSource*, TBField*, std::vector<int>&, bool El = false /**< Electrons/positrons need special treatment*/);
    /**
     * @fn TDiffusionCoefficient(TGrid*, Input*, bool El = false);
     * @brief Constructor given a geometry, and according to user's input.
     */
    
    ~TDiffusionCoefficient2D() { }
    /**< Destructor. */
    
};

/**
 * @class TDiffusionCoefficient3D
 * @author Luca Maccione
 * @author Daniele Gaggero
 * @email luca.maccione@lmu.de
 * @email daniele.gaggero@sissa.it
 * @brief Position dependent diffusion coefficient. Part of it independent of the nucleus. Needed coefficient will be restored in TCREvolutor.
 */
class TDiffusionCoefficient3D : public TDiffusionCoefficient {
    
public:
    TDiffusionCoefficient3D() {} /**< Default constructor. */
    TDiffusionCoefficient3D(TGrid* /**< Geometry. */, Input*, TSource*, TBField*, TGeometry*, double A, double Z, std::vector<int>, bool El = false /**< Electrons/positrons need special treatment*/);
    /**
     * @fn TDiffusionCoefficient(TGrid*, Input*, bool El = false);
     * @brief Constructor given a geometry, and according to user's input.
     */
    
    ~TDiffusionCoefficient3D() { }
    /**< Destructor. */
    
protected:
    
    double GetProfileDPar(double, double, double, TSource*);
    double GetXDerivativeDPar(double, double, double, TSource*);
    double GetYDerivativeDPar(double, double, double, TSource*);
    double GetZDerivativeDPar(double, double, double, TSource*);
    std::vector<double> GetXDerivativeBF(double, double, double, TBField*);
    std::vector<double> GetYDerivativeBF(double, double, double, TBField*);
    std::vector<double> GetZDerivativeBF(double, double, double, TBField*);

    inline double GetProfileDPerp(double x, double y, double z, TSource* ST) { return GetProfile(x, y, z, ST); }
    std::vector<double> GetGradProfileDPerp(double, double, double, TSource*);
    std::vector<double> GetGradProfileDPar(double, double, double, TSource*);
    std::map<char, std::vector<double> > GetGradVersors(double, double, double, TBField*);
};


#endif
