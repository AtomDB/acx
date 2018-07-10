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

extern "C" 
void acxion(const RealArray& energyArray, 
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
  bool *ZhasCXptr;
  double *ionpopZ[MAXIMUM_Z+1];

  float *abZ;
  abZ = (float *) malloc((MAXIMUM_Z+1)*sizeof(float));
  get_abund(&abZ);

  const Real fracHe = params[0];
  const Real fracH = 1.0 - fracHe;
  const int ZElement = (int) params[1];
  const int rmJIon = (int) params[2];
  const Real z = params[3];
  float vabZ[MAXIMUM_Z+1];
  int  model=1;

  for (int ii=0; ii<MAXIMUM_Z+1;ii++) vabZ[ii] = abZ[ii];

  if (isFirst) {

    ZhasCXptr = (bool *) malloc((MAXIMUM_Z+1)*sizeof(bool));
    for (int iZ=1;iZ<MAXIMUM_Z;iZ++) ZhasCXptr[iZ] = false;

    acx_init(ZhasCXptr);
   
    for (int iZ=1;iZ<MAXIMUM_Z;iZ++) ZhasCX[iZ] = ZhasCXptr[iZ];

    /* Make plenty of space for the ionization balances */
    for (int iZ=1;iZ<MAXIMUM_Z;iZ++) {
      if (ZhasCX[iZ]) {
	ionpopZ[iZ] = (double *) malloc((iZ+1)*sizeof(double));
      } else {
	ionpopZ[iZ] = NULL;
      }
    }
  }
  
  if (((int) params[4] < 1)||((int) params[4] > 2*NUMCXMODELS+0.5)) {
    errMsg << "\nacx:\nModel " << 
      (int) params[4] << " must be between 1 and " << 2*NUMCXMODELS << "\n";
    throw YellowAlert(errMsg.str());
  } else {
    model = (int) params[4];
  }

  size_t N(energyArray.size());
  flux.resize(N-1);

  /* ************************** */
  /* First get all the ion pops */
  /* ************************** */
  for (size_t iZ=1; iZ<MAXIMUM_Z+1;iZ++) {
    if (ZhasCX[iZ]) {
      for (size_t jj=0;jj<iZ+1;jj++) ionpopZ[iZ][jj]=0.0;
    }
  }
  
  ionpopZ[ZElement][rmJIon] = 1.0;  /* CHECK THIS */
  //  printf("check this.\n");

  /* ************************** */
  /* Now do the CX calculation  */
  /* ************************** */

  calc_cx(energyArray, flux, model, fracH, fracHe, z, vabZ, ionpopZ);

  /* And we're done */

  isFirst = false;
  fluxError.resize(0);
  return;
}

