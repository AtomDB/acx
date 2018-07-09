#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acx.h"
#include "dacx.h"
#include "messages.h"

/* BUGFIX 2015-06-26, Adam Foster. 
 * Replaced NUMCXMODELS with 2*NUMCXMODELS (line 62)
 * This is a fix so DACX can access models 9 and up
 * 
 */

void calc_cx(struct PARAMETERS *params, const int model,
	     const double fracH, const double fracHe, 
	     const double z, float *vabZ, double **ionpopZ);
void find_ion_bal(int Z, double Te, double Ne, double Ion_pop[]);
void acx_init(int *ptr);
void get_abund(float **abZ);
void get_input_data(int argc, char *argv[], Param_File_Type *p_input,
		    struct PARAMETERS *params, int *status);

int main(int argc, char *argv[]) {

  static double *ionpopZ[MAXIMUM_Z+1];
  int *ZhasCX;
  double *SumIonpopZ[MAXIMUM_Z+1];
  double ionpop[MAXIMUM_Z+1];
  float vabZ[MAXIMUM_Z+1];
  int status=0;
  int ii, jj, iZ;

  struct PARAMETERS pdata, *params;

  params = &pdata;

  get_input_data(argc, argv, NULL, params, &status);
  if (status != 0) {
    printf("dacx: Error reading input data.  Exiting without running.\n");
    return(1); /* error reading input */
  }

  float *abZ;
  abZ = (float *) malloc((MAXIMUM_Z+1)*sizeof(float));
  get_abund(&abZ);

  const double kT = params->kT;
  const double fracHe = params->FracHe;
  const double fracH = 1.0 - fracHe;

  for (ii=0; ii<MAXIMUM_Z+1;ii++) vabZ[ii] = abZ[ii];
  vabZ[6]  = params->Abundance*abZ[6]*params->AbC;   // C
  vabZ[7]  = params->Abundance*abZ[7]*params->AbN;   // N
  vabZ[8]  = params->Abundance*abZ[8]*params->AbO;   // O
  vabZ[10] = params->Abundance*abZ[10]*params->AbNe; // Ne
  vabZ[12] = params->Abundance*abZ[12]*params->AbMg; // Mg
  vabZ[13] = params->Abundance*abZ[13]*params->AbAl; // Al
  vabZ[14] = params->Abundance*abZ[14]*params->AbSi; // Si
  vabZ[16] = params->Abundance*abZ[16]*params->AbS; // S
  vabZ[18] = params->Abundance*abZ[18]*params->AbAr; // Ar
  vabZ[20] = params->Abundance*abZ[20]*params->AbCa;// Ca
  vabZ[26] = params->Abundance*abZ[26]*params->AbFe;// Fe
  vabZ[28] = params->Abundance*abZ[28]*params->AbNi;// Ni
  const double z = params->redshift;
  int model;

  if ((params->Model < 1)||(params->Model > 2*NUMCXMODELS)) {
    errmess("dacx","Model %d must be between 1 and 8.",params->Model);
  } else {
    model = params->Model;
  }

  ZhasCX = (int *) malloc((MAXIMUM_Z+1)*sizeof(int));
  for (iZ=1;iZ<MAXIMUM_Z+1;iZ++) ZhasCX[iZ] = FALSE;

  acx_init(ZhasCX);

  /* Make plenty of space for the ionization balances */
  for (iZ=1;iZ<MAXIMUM_Z+1;iZ++) {
    if (ZhasCX[iZ]) {
      ionpopZ[iZ]    = (double *) malloc((iZ+1)*sizeof(double));
      SumIonpopZ[iZ] = (double *) malloc((iZ+1)*sizeof(double));
    } else {
      ionpopZ[iZ] = NULL;
      SumIonpopZ[iZ] = NULL;
    }
  }  

  /* ************************** */
  /* First get all the ion pops */
  /* ************************** */
  for (iZ=1; iZ<MAXIMUM_Z+1;iZ++) {
    if (ZhasCX[iZ]) {
      find_ion_bal(iZ, kT/KBOLTZ, 0.0, ionpop);
      for (jj=0;jj<iZ+1;jj++) ionpopZ[iZ][jj]=ionpop[jj];
      SumIonpopZ[iZ][0] = ionpopZ[iZ][0];
      for (jj=1;jj<iZ+1;jj++) 
	SumIonpopZ[iZ][jj]=ionpop[jj] + SumIonpopZ[iZ][jj-1];
    }
  }

  /* ************************** */
  /* Now do the CX calculation  */
  /* ************************** */

  if (params->swcx) {
    calc_cx(params, model, fracH, fracHe, z, vabZ, ionpopZ);
  } else {
    calc_cx(params, model, fracH, fracHe, z, vabZ, SumIonpopZ);
  }    

  /* And we're done */
  return 0;
}
