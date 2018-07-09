#include <XSFunctions/Utilities/funcType.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/functionMap.h>
#include "xsTypes.h"
#include <cmath>
#include <iomanip>

//#include <stlToCArrays.h>
#include <sstream>
#include <iostream>
#include <string>
#include <fitsio.h>
#include <fstream>

using namespace XSutility;
using namespace std;

#include "acx.h"
void acx_init(bool *ptr);
void get_abund(float **abZ);
void calc_cx(const RealArray& energyArray, RealArray& flux, 
	     const int model,  const Real fracH, const Real fracHe, 
	     const Real z, float *vabZ, double **ionpopZ);
void find_ion_bal(int Z, double Te, double Ne, double Ion_pop[]);

extern "C" 
void acxfunc(const RealArray& energyArray, 
	     const RealArray& params, 
	     int spectrumNumber, 
	     RealArray& flux, 
	     RealArray& fluxError, 
	     const string& initString) {
  
  /* energy : size Nflux+1 (has low and high bins, in keV)
     parameter: all the parameters
     spectrum: number of component being calculated.  Think I can ignore.
     flux: output flux
     fluxError: ignore
     init: initialization string - may hold location of input file.
  */

  std::ostringstream errMsg;
  static bool isFirst = true;
  static bool ZhasCX[MAXIMUM_Z+1];
  static double *ionpopZ[MAXIMUM_Z+1];
  static double *SumIonpopZ[MAXIMUM_Z+1];
  bool *ZhasCXptr;
  double ionpop[MAXIMUM_Z+1];
  float vabZ[MAXIMUM_Z+1];
  float *abZ;
  abZ = (float *) malloc((MAXIMUM_Z+1)*sizeof(float));
  get_abund(&abZ);

  const Real kT = params[0];
  const Real fracHe = params[1];
  const Real fracH = 1.0 - fracHe;
  const Real abund = params[2];
  const Real z = params[3];
  int  model = 1;
  int  swcx = (int) params[4];  
  /* swcx=1 for SWCX model, where each ion undergoes only one CX reaction */
  /* swcx=0 is used for all non-solar-wind cases, as the CX runs to completion
     and each ion is neutralized */

  for (int ii=0; ii<MAXIMUM_Z+1;ii++) vabZ[ii] = abund*abZ[ii];

  if (isFirst) {
    ZhasCXptr = (bool *) malloc((MAXIMUM_Z+1)*sizeof(bool));
    for (int iZ=1;iZ<MAXIMUM_Z;iZ++) ZhasCXptr[iZ] = false;

    acx_init(ZhasCXptr);
   
    for (int iZ=1;iZ<MAXIMUM_Z;iZ++) ZhasCX[iZ] = ZhasCXptr[iZ];

    /* Make plenty of space for the ionization balances */
    for (int iZ=1;iZ<MAXIMUM_Z;iZ++) {
      if (ZhasCX[iZ]) {
	ionpopZ[iZ] = (double *) malloc((iZ+1)*sizeof(double));
	SumIonpopZ[iZ] = (double *) malloc((iZ+1)*sizeof(double));
      } else {
	ionpopZ[iZ] = NULL;
	SumIonpopZ[iZ] = NULL;
      }
    }
  }

  if (((int) params[5] < 1)||((int) params[5] > 2*NUMCXMODELS+0.5)) {
    errMsg << "\nacx:\nModel " << 
      (int) params[5] << " must be between 1 and " << 2*NUMCXMODELS << "\n";
    throw YellowAlert(errMsg.str());
  } else {
    model = (int) params[5];
  }
  
  size_t N(energyArray.size());
  flux.resize(N-1);

  /* ************************** */
  /* First get all the ion pops */
  /* ************************** */
  for (int iZ=1; iZ<MAXIMUM_Z+1;iZ++) {
    if (ZhasCX[iZ]) {
      find_ion_bal(iZ, kT/KBOLTZ, 0.0, ionpop);
      for (int jj=0;jj<iZ+1;jj++) ionpopZ[iZ][jj]=ionpop[jj];
      if (swcx!=1) { /* Don't bother if SWCX case*/
	SumIonpopZ[iZ][iZ] = ionpopZ[iZ][iZ];
	for (int jj=iZ-1;jj>=0;jj--) 
	  SumIonpopZ[iZ][jj]=ionpop[jj] + SumIonpopZ[iZ][jj+1];
      }
    }
  }

  /* ************************** */
  /* Now do the CX calculation  */
  /* ************************** */

  if (swcx==1) {
    calc_cx(energyArray, flux, model, fracH, fracHe, z, vabZ, ionpopZ);
  } else {
    calc_cx(energyArray, flux, model, fracH, fracHe, z, vabZ, SumIonpopZ);
  }

  /* And we're done */

  isFirst = false;
  fluxError.resize(0);
  return;
}

extern "C" 
void vacxfunc(const RealArray& energyArray, 
	     const RealArray& params, 
	     int spectrumNumber, 
	     RealArray& flux, 
	     RealArray& fluxError, 
	     const string& initString) {
  
  /* energy : size Nflux+1 (has low and high bins, in keV)
     parameter: all the parameters
     spectrum: number of component being calculated.  Think I can ignore.
     flux: output flux
     fluxError: ignore
     init: initialization string - may hold location of input file.
  */

  std::ostringstream errMsg;
  static bool isFirst = true;
  static bool ZhasCX[MAXIMUM_Z+1];
  static double *ionpopZ[MAXIMUM_Z+1];
  static double *SumIonpopZ[MAXIMUM_Z+1];
  bool *ZhasCXptr;
  double ionpop[MAXIMUM_Z+1];
  float vabZ[MAXIMUM_Z+1];

  float *abZ;
  abZ = (float *) malloc((MAXIMUM_Z+1)*sizeof(float));
  get_abund(&abZ);

  const Real kT = params[0];
  const Real fracHe = params[1];
  const Real fracH = 1.0 - fracHe;
  for (int ii=0; ii<MAXIMUM_Z+1;ii++) vabZ[ii] = abZ[ii];
  vabZ[6] = abZ[6]*params[2];   // C
  vabZ[7] = abZ[7]*params[3];   // N
  vabZ[8] = abZ[8]*params[4];   // O
  vabZ[10] = abZ[10]*params[5]; // Ne
  vabZ[12] = abZ[12]*params[6]; // Mg
  vabZ[13] = abZ[13]*params[7]; // Al
  vabZ[14] = abZ[14]*params[8]; // Si
  vabZ[16] = abZ[16]*params[9]; // S
  vabZ[18] = abZ[18]*params[10]; // Ar
  vabZ[20] = abZ[20]*params[11];// Ca
  vabZ[26] = abZ[26]*params[12];// Fe
  vabZ[28] = abZ[28]*params[13];// Ni

  int  swcx = (int) params[15];  
  /* swcx=1 for SWCX model, where each ion undergoes only one CX reaction */
  /* swcx=0 is used for all non-solar-wind cases, as the CX runs to completion
     and each ion is neutralized */

  if (isFirst) {
    ZhasCXptr = (bool *) malloc((MAXIMUM_Z+1)*sizeof(bool));
    for (int iZ=1;iZ<MAXIMUM_Z+1;iZ++) ZhasCXptr[iZ] = false;

    acx_init(ZhasCXptr);
   
    for (int iZ=1;iZ<MAXIMUM_Z+1;iZ++) ZhasCX[iZ] = ZhasCXptr[iZ];

    /* Make plenty of space for the ionization balances */
    for (int iZ=1;iZ<MAXIMUM_Z+1;iZ++) {
      if (ZhasCX[iZ]) {
	ionpopZ[iZ] = (double *) malloc((iZ+1)*sizeof(double));
	SumIonpopZ[iZ] = (double *) malloc((iZ+1)*sizeof(double));
      } else {
	ionpopZ[iZ] = NULL;
	SumIonpopZ[iZ] = NULL;
      }
    }
  }
  
  const Real z = params[14];
  int  model=1;
  if (((int) params[16] < 1)||((int) params[16] > 2*NUMCXMODELS+0.5)) {
    errMsg << "\nvacx:\nModel #" << 
      (int) params[16] << " must be between 1 and " << 2*NUMCXMODELS << "\n";
    throw YellowAlert(errMsg.str());
  } else {
    model = (int) params[16];
  }

  size_t N(energyArray.size());
  flux.resize(N-1);

  /* ************************** */
  /* First get all the ion pops */
  /* ************************** */
  for (int iZ=1; iZ<MAXIMUM_Z+1;iZ++) {
    if (ZhasCX[iZ]) {
      find_ion_bal(iZ, kT/KBOLTZ, 0.0, ionpop);
      for (int jj=0;jj<iZ+1;jj++) ionpopZ[iZ][jj]=ionpop[jj];
      if (swcx!=1) {
	SumIonpopZ[iZ][iZ] = ionpopZ[iZ][iZ];
	for (int jj=iZ-1;jj>=0;jj--) 
	  SumIonpopZ[iZ][jj]=ionpop[jj] + SumIonpopZ[iZ][jj+1];
      }
    }
  }

  /* ************************** */
  /* Now do the CX calculation  */
  /* ************************** */

  if (swcx==1) {
    calc_cx(energyArray, flux, model, fracH, fracHe, z, vabZ, ionpopZ);
  } else {
    calc_cx(energyArray, flux, model, fracH, fracHe, z, vabZ, SumIonpopZ);
  }

  /* And we're done */

  isFirst = false;
  fluxError.resize(0);
  return;
}

