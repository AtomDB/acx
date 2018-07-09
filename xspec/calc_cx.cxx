#include <XSUtil/FunctionUtils/funcType.h>
#include <XSUtil/Numerics/Integrate.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/FunctionUtils/FunctionUtility.h>
#include <XSUtil/FunctionUtils/xsFortran.h>
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

#include "acx.h"

using namespace XSutility;
using namespace std;

struct ION_BAL_TABLE *ionbal;

int Ncxh, NcxhOLD;
int Ncxhlines, NcxhlinesOLD;
struct CX_TABLE cxh, cxhOLD;

int Ncxhe, NcxheOLD;
int Ncxhelines, NcxhelinesOLD;
struct CX_TABLE cxhe, cxheOLD;

void read_ion_bal(const char *ionbal_file);
void find_ion_bal(int Z, double Te, double Ne, double Ion_pop[]);

double A_E_twoph(double A, double E0, double E);
double twophot_lambda(int Z, int rmJ);
void fits_read_cx_data(const char *file, int *Ncxlines, int *Ncx, 
		       struct CX_TABLE *cx);

void calc_cx(const RealArray& energyArray, RealArray& flux, 
	     const int model,  const Real fracH, const Real fracHe, 
	     const Real z, float *vabZ, double **ionpopZ) {
    
  std::ostringstream errMsg;
  int iCX;
  int iNeut;
  int Ncxl=0;
  int Ncxtot=0;
  int Actualmodel=0;
  Real fracNeutral=0.0;
  Real enz;
  bool done;
  struct CX_TABLE *cxm=NULL;

  size_t N(energyArray.size());

  /* ************************** */
  /* neutral loop: H, then He   */
  /* ************************** */

  for(iNeut=0;iNeut<2;iNeut++) {

    iCX = 0;
    done = false;

    if (iNeut==0) { /* H case */
      if (model < 9) {
	Actualmodel = model;
	cxm = &cxh;
	Ncxl = Ncxhlines;
	Ncxtot = Ncxh;
      } else {
	Actualmodel = model-8;
	cxm = &cxhOLD;
	Ncxl = NcxhlinesOLD;
	Ncxtot = NcxhOLD;
      }
      fracNeutral = fracH;
    }
    if (iNeut==1) { /* He case */
      if (model < 9) {
	Actualmodel = model;
	cxm = &cxhe;
	Ncxl = Ncxhelines;
	Ncxtot = Ncxhe;
      } else {
	Actualmodel = model - 8;
	cxm = &cxheOLD;
	Ncxl = NcxhelinesOLD;
	Ncxtot = NcxheOLD;
      }	
      fracNeutral = fracHe;
    }
    
    enz = cxm->energy[iCX]/(1.+z); /* Redshifted energy of CX photon */

    for (size_t ii=0; ii<N-1; ii++) {
      if (iNeut==0) flux[ii] = 0.0;
      /* Assuming energy-sorted list of CX lines, step through until we */ 
      /* find one that is inside our lower bound */
      while ((!done)&&(enz < energyArray[ii])) {  
	iCX++;
	enz = cxm->energy[iCX]/(1.+z);
	if (iCX == Ncxl) done = true;
      }
      /* iCX is now cued up to a 'real' photon */ 
      while ((!done)&&(enz < energyArray[ii+1])) { // calc emission
	/* Note that ionpopZ[cxm.Z[iCX]][0] is for the neutral state, rmJ=1 */
	/* So ionpopZ[cxm.Z[iCX]][cxm.rmJ[iCX]] is not for Z, rmJ but rather */
	/* for the next ion stage up, as it should be */
	/* printf("%d, %e\n",iCX,cxm->feps[iCX*NUMCXMODELS + model-1]); */
	flux[ii] += 1e5 * fracNeutral * vabZ[cxm->Z[iCX]] *
	  (ionpopZ[cxm->Z[iCX]][cxm->rmJ[iCX]]) * 
	  cxm->feps[iCX*NUMCXMODELS + Actualmodel-1];
	iCX++;
	if (iCX == Ncxl) {
	  done = true;
	} else {
	  enz = cxm->energy[iCX]/(1.+z);
	}
      }
    }
  
    /* Now we go through and add in two photon continua */ 
    while (iCX < Ncxtot) {
      
      double E0 = HC_IN_KEV_A/twophot_lambda(cxm->Z[iCX], cxm->rmJ[iCX]);
      
      if (ionpopZ[cxm->Z[iCX]][cxm->rmJ[iCX]] > ION_POP_MIN) {
	for (size_t iBin=0; iBin<N-1; iBin++) {
	  double E = 0.5*(energyArray[iBin] + energyArray[iBin+1]);
	  if (E < E0) {
	    double dE = energyArray[iBin+1] - energyArray[iBin];
	    flux[iBin] += 1e5 * fracNeutral * vabZ[cxm->Z[iCX]] * 
	      ionpopZ[cxm->Z[iCX]][cxm->rmJ[iCX]] * 
	      cxm->feps[iCX*NUMCXMODELS + Actualmodel-1] * 
	      (dE/E0) * A_E_twoph(1.0, E0, E); 
	  }
	}
      }
      iCX++;
    }
  }
}

void acx_init(bool *ZhasCX) {

  std::ostringstream errMsg;
  /* 
  string cxfile[NUMCXMODELS];
  string cxhefile[NUMCXMODELS];
  */
  string ibfile, cxfile, cxhefile, cxfileOLD, cxhefileOLD;
  int Ncx, Ncxlines;
  
  ibfile      = DIRECTORY "/data/v2.0.2_ionbal.fits";
  cxfile      = DIRECTORY "/data/AtomDB_CX_H_v1.0.0.fits";
  cxhefile    = DIRECTORY "/data/AtomDB_CX_He_v1.0.0.fits";
  cxfileOLD   = DIRECTORY "/data/AtomDB_CX_H_v0.5.1.fits";
  cxhefileOLD = DIRECTORY "/data/AtomDB_CX_He_v0.5.1.fits";

  try {
    ifstream file_exists(ibfile.c_str());
      if (!file_exists) {
	errMsg << "\nacx:\nError finding ionization balance file.\n";
	throw YellowAlert(errMsg.str());
      }
      read_ion_bal(ibfile.c_str());
  }
  catch (...) {
    return;
  }

  try {
    ifstream file_exists(cxfile.c_str());
    if (!file_exists) {
      errMsg << "\nacx:\nError finding charge exchange file" << cxfile << "\n";
      throw YellowAlert(errMsg.str());
    }
    fits_read_cx_data(cxfile.c_str(), &Ncxlines, &Ncx, &cxh);
    Ncxh = Ncx;
    Ncxhlines = Ncxlines;
    for (int iCX=0;iCX<Ncx;iCX++) ZhasCX[(cxh.Z[iCX])] = true;
  }
  catch (...) {
    return;
  }

  try {
    ifstream file_exists(cxhefile.c_str());
    if (!file_exists) {
      errMsg << "\nacx:\nError finding charge exchange file" <<cxhefile<<"\n";
      throw YellowAlert(errMsg.str());
    }
    fits_read_cx_data(cxhefile.c_str(), &Ncxlines, &Ncx, &cxhe);
    Ncxhe = Ncx;
    Ncxhelines = Ncxlines;
    for (int iCX=0;iCX<Ncxhe;iCX++) ZhasCX[(cxhe.Z[iCX])] = true;
  }
  catch (...) {
    return;
  }

  try {
    ifstream file_exists(cxfileOLD.c_str());
    if (!file_exists) {
      errMsg << "\nacx:\nError finding charge exchange file" << 
	cxfileOLD << "\n";
      throw YellowAlert(errMsg.str());
    }
    fits_read_cx_data(cxfileOLD.c_str(), &Ncxlines, &Ncx, &cxhOLD);
    NcxhOLD = Ncx;
    NcxhlinesOLD = Ncxlines;
    for (int iCX=0;iCX<Ncx;iCX++) ZhasCX[(cxhOLD.Z[iCX])] = true;
  }
  catch (...) {
    return;
  }

  try {
    ifstream file_exists(cxhefileOLD.c_str());
    if (!file_exists) {
      errMsg << "\nacx:\nError finding charge exchange file" << 
	cxhefileOLD <<"\n";
      throw YellowAlert(errMsg.str());
    }
    fits_read_cx_data(cxhefileOLD.c_str(), &Ncxlines, &Ncx, &cxheOLD);
    NcxheOLD = Ncx;
    NcxhelinesOLD = Ncxlines;
    for (int iCX=0;iCX<Ncx;iCX++) ZhasCX[(cxheOLD.Z[iCX])] = true;
  }
  catch (...) {
    return;
  }

}

void fits_read_cx_data(const char *file, int *Ncxlines, int *Ncx, 
		       struct CX_TABLE *cx) {

  fitsfile *fptr;
  int hdunum, casesen;
  
  int colnum_Z, colnum_rmJ, colnum_UL, colnum_LL, colnum_Lambda, colnum_Prob;
  char column[MAXSTRLEN],comment[MAXSTRLEN];
  int status, hdutype, anynul=0, nulval=0;
  int nrows;


  std::ostringstream errMsg;
  status = 0;

  if (fits_open_file(&fptr, file, READONLY, &status)) {
    errMsg << status << "\nFailed to open " << file << "\n";
    throw YellowAlert(errMsg.str());
  }

  hdunum = 2;
  if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) {
    errMsg << status << "\nread_ion_bal\nError finding second HDU\n";
    throw YellowAlert(errMsg.str());
  }
    
  casesen = FALSE;
  
  fits_get_colname(fptr, casesen, "ELEMENT", column,
		   &colnum_Z, &status);
  fits_get_colname(fptr, casesen, "ION", column,
		   &colnum_rmJ, &status);
  fits_get_colname(fptr, casesen, "UPPERLEV", column,
		   &colnum_UL, &status);
  fits_get_colname(fptr, casesen, "LOWERLEV", column,
		   &colnum_LL, &status);
  fits_get_colname(fptr, casesen, "WAVELEN", column,
		   &colnum_Lambda, &status);
  fits_get_colname(fptr, casesen, "LINEPROB", column,
		   &colnum_Prob, &status);

  if (status) {
    errMsg << status << "\nread_cx_data\nError getting column names\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_read_key(fptr, TINT, "NAXIS2", &nrows, comment, &status)) {
    errMsg << status << "\nread_cx_data\nError getting number of rows\n";
    throw YellowAlert(errMsg.str());
  }

  *Ncx = nrows;
  *Ncxlines = nrows;

  cx->Z      = (int *) malloc(nrows*sizeof(int));
  cx->rmJ    = (int *) malloc(nrows*sizeof(int));
  cx->ul     = (int *) malloc(nrows*sizeof(int));
  cx->ll     = (int *) malloc(nrows*sizeof(int));
  cx->lambda = (float *) malloc(nrows*sizeof(float));
  cx->energy = (float *) malloc(nrows*sizeof(float));
  cx->feps   = (float *) malloc(nrows*NUMCXMODELS*sizeof(float));

  if (fits_read_col(fptr, TINT, colnum_Z, 1, 1, nrows, 0,
		    cx->Z, &anynul, &status)) {
    errMsg << status << "\nread_cx_data\nError reading element Z value\n";
    throw YellowAlert(errMsg.str());
  }

  if (fits_read_col(fptr, TINT, colnum_rmJ, 1, 1, nrows, &nulval,
		    cx->rmJ, &anynul, &status)) {
    errMsg << status << "\nread_cx_data\nError reading ion rmJ value\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_read_col(fptr, TINT, colnum_UL, 1, 1, nrows, &nulval,
		    cx->ul, &anynul, &status)) {
    errMsg << status << "\nread_cx_data\nError reading upper level value\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_read_col(fptr, TINT, colnum_LL, 1, 1, nrows, &nulval,
		    cx->ll, &anynul, &status)) {
    errMsg << status << "\nread_cx_data\nError reading lower level value\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_read_col(fptr, TFLOAT, colnum_Lambda, 1, 1, nrows, &nulval,
		    cx->lambda, &anynul, &status)) {
    errMsg << status << "\nread_cx_data\nError reading wavelength value\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_read_col(fptr, TFLOAT, colnum_Prob, 1, 1, nrows*NUMCXMODELS,&nulval,
		    cx->feps, &anynul, &status)) {
    errMsg << status << "\nread_cx_data\nError reading line prob value\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_close_file(fptr, &status)) {
    errMsg << status << "\nread_cx_data\nError closing fits file\n";
    throw YellowAlert(errMsg.str());
  }
  /* 
  printf("%d %d %d %d  %e :  %e %e %e %e %e %e %e %e\n",
	 cx->Z[1], cx->rmJ[1], cx->ul[1], cx->ll[1], 
	 cx->lambda[1], 
	 cx->feps[8+0],cx->feps[8+1],cx->feps[8+2],cx->feps[8+3],
	 cx->feps[8+4],cx->feps[8+5],cx->feps[8+6],cx->feps[8+7]);
  */
  /* Now work out actual number of lines. */
  for (int ii=0;ii<nrows;ii++) {
    if (cx->lambda[ii] > 0.0) { /* Not a two-photon transition */
      cx->energy[ii] = HC_IN_KEV_A/cx->lambda[ii];
    } else {
      cx->energy[ii] = TWOPHOT;
      (*Ncxlines)--; /* Only counts REAL lines */
      /* two photon transitions are counted separately. */
    }
  }
}

void find_ion_bal(int Z, double Te, double Ne, double Ion_pop[]) {

  int i,iIonState, NT, ND;
  int iTlow=-1, iThigh=-1, iNlow=-1, iNhigh=-1;
  int indx;
  double w_Tlow=0., w_Thigh=0., w_Nlow=0., w_Nhigh=0.;
  double dT, dN;
  double a1,a2,a3,a4;
  std::ostringstream errMsg;
  
  NT = ionbal->Ntemp;
  ND = ionbal->Ndens;
  
  /* printf("find_ion_bal: NT = %d, ND = %d\n",NT, ND); */
  
  /* Now, find where we are...*/
  for (i=0;i<NT;i++) {  /* Assumes T is monotonic increasing */
    if (((ionbal->T)[i] > Te)&&(iTlow==-1)) {
      if (i==0) {
	errMsg << "\nfind_ion_bal\nTemp = " << Te <<
	  " out of range for table, minimum value = " << 
	  (ionbal->T)[i] << "\n";
	throw YellowAlert(errMsg.str());
      } else {
	iTlow = i-1;
	iThigh = i;
	dT = (ionbal->T)[iThigh] - (ionbal->T)[iTlow];
	w_Tlow = 1 - (Te - (ionbal->T)[iTlow])/dT;
	w_Thigh= 1 - w_Tlow;
      }
    }
  }
  if (iTlow==-1) {
    	errMsg << "\nfind_ion_bal\nTemp = " << Te <<
	  " out of range for table\n";
	throw YellowAlert(errMsg.str());
  }

  if (ND > 1) {
    for (i=0;i<ND;i++) {  /* Assumes N is monotonic increasing */
      if (((ionbal->dens)[i] > Ne)&&(iNlow==-1)) {
	if (i==0) {
	  errMsg << "\nfind_ion_bal\nDensity = " << Ne <<
	    " out of range for table, minimum value = " << 
	    (ionbal->dens)[i] << "\n";
	  throw YellowAlert(errMsg.str());
	} else {
	  iNlow = i-1;
	  iNhigh = i;
	  dN = (ionbal->dens)[iThigh] - (ionbal->dens)[iNlow];
	  w_Nlow = 1 - (Ne - (ionbal->dens)[iNlow])/dN;
	  w_Nhigh= 1 - w_Nlow;
	}
      }
    }
  } else {
    iNlow = 0;
    iNhigh = 0;
    w_Nlow = 1.;
    w_Nhigh = 0.;
  }
  
  for(iIonState=0;iIonState<(Z+1);iIonState++) {
    indx = ionbal->indx[Z] + iIonState;
    a1 = (ionbal->tableX[iTlow*ND +iNlow])[indx];
    a2 = (ionbal->tableX[iThigh*ND+iNlow])[indx];
    a3 = (ionbal->tableX[iTlow*ND +iNhigh])[indx];
    a4 = (ionbal->tableX[iThigh*ND+iNhigh])[indx];
    if ((a1>0)&&(a2>0)&&(a3>0)&(a4>0)) {
      a1 = log10(a1); a2 = log10(a2); a3 = log10(a3); a4 = log10(a4);
      Ion_pop[iIonState] = pow(10.,(a1*w_Tlow*w_Nlow + a2*w_Thigh*w_Nlow +
				    a3*w_Tlow*w_Nhigh+ a4*w_Thigh*w_Nhigh));
    } else {
      Ion_pop[iIonState] = a1*w_Tlow*w_Nlow + a2*w_Thigh*w_Nlow +
	a3*w_Tlow*w_Nhigh+ a4*w_Thigh*w_Nhigh;
    }
    
    if (Ion_pop[iIonState] == 0) {
      Ion_pop[iIonState] = 1/BLOWUP;
    }
  }
}

/* *********************************************************** */
/* read_ion_bal: function to read ionization balance from file */
/* *********************************************************** */

void read_ion_bal(const char *ionbal_file) {
  
  fitsfile *fptr;
  int hdunum, casesen;
  
  int colnum_temp, colnum_dens, colnum_Z, colnum_X;
  char column[MAXSTRLEN],comment[MAXSTRLEN];
  int loc_indx[MAXIMUM_Z+1];
  int status, hdutype, anynul;
  int nrows, NZ, NX, Ntemp, Ndens;
  int Tlinear, Nlinear;
  int row;
  int iTemp=0, iDens=0, iZ, i;
  
  ionbal = (struct ION_BAL_TABLE *) malloc(sizeof(struct ION_BAL_TABLE));

  std::ostringstream errMsg;
  status = 0;

  if (fits_open_file(&fptr, ionbal_file, READONLY, &status)) {
    errMsg << status << "\nFailed to open " << ionbal_file << "\n";
    throw YellowAlert(errMsg.str());
  }

  hdunum = 2;
  if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) {
    errMsg << status << "\nread_ion_bal\nError finding second HDU\n";
    throw YellowAlert(errMsg.str());
  }
    
  casesen = FALSE;
  
  fits_get_colname(fptr, casesen, "Temperature", column,
		   &colnum_temp, &status);
  fits_get_colname(fptr, casesen, "Density", column,
		   &colnum_dens, &status);
  fits_get_colname(fptr, casesen, "Z_ELEMENT", column,
		   &colnum_Z, &status);
  fits_get_colname(fptr, casesen, "X_IONPOP", column,
		   &colnum_X, &status);

  if (status) {
    errMsg << status << "\nread_ion_bal\nError getting column names\n";
    throw YellowAlert(errMsg.str());
  }
  
  if (fits_read_key(fptr, TINT, "NAXIS2", &nrows, comment, &status)) {
    errMsg << status << "\nread_ion_bal\nError getting number of rows\n";
    throw YellowAlert(errMsg.str());
  }
  if (fits_read_key(fptr, TINT, "T_NUMBER", &Ntemp, comment, &status)) {
    errMsg << status << "\nread_ion_bal\nError getting number of temps\n";
    throw YellowAlert(errMsg.str());
  }
  if (fits_read_key(fptr, TLOGICAL, "T_LINEAR", &Tlinear, comment, &status)) {
    errMsg << status << "\nread_ion_bal\nError getting temperature logical\n";
    throw YellowAlert(errMsg.str());
  }
  if (fits_read_key(fptr, TINT, "N_NUMBER", &Ndens, comment, &status)) {
    errMsg << status << "\nread_ion_bal\nError getting number of densities\n";
    throw YellowAlert(errMsg.str());
  }
  if (fits_read_key(fptr, TLOGICAL, "N_LINEAR", &Nlinear, comment, &status)) {
    errMsg << status << "\nread_ion_bal\nError getting density logical\n";
    throw YellowAlert(errMsg.str());
  }
  if (fits_read_key(fptr, TINT, "N_ELEMEN", &NZ, comment, &status)) {
    errMsg << status << "\nread_ion_bal\nError getting number of elements\n";
    throw YellowAlert(errMsg.str());
  }
  if (fits_read_key(fptr, TINT, "N_IONS", &NX, comment, &status)) {
    errMsg << status << "\nread_ion_bal\nError getting number of ions\n";
    throw YellowAlert(errMsg.str());
  }

  ionbal->NZ = NZ;
  ionbal->Nrows = nrows;
  ionbal->Ntemp = Ntemp;
  ionbal->Ndens = Ndens;
  ionbal->Ntime = 1; /* reading an ionization balance implies only one time. */
  ionbal->T = (float *) malloc(ionbal->Ntemp*sizeof(float));
  ionbal->dens = (float *) malloc(ionbal->Ndens*sizeof(float));
  ionbal->time = (float *) malloc(ionbal->Ntime*sizeof(float));

  ionbal->time[0] = 0.0; /* Defined as equilibrium run */

  if (Ntemp*Ndens != nrows) {
    errMsg << status << "\nread_ion_bal\nNumber of rows" << nrows << 
      "!= Number of temperatures (" << Ntemp <<") * Number of densities (" <<
      Ndens <<")\n";
    throw YellowAlert(errMsg.str());
  }

  /* ionbal->tableX = (float **) malloc(Ntemp*Ndens*sizeof(float)); */
  ionbal->tableX = (float **) malloc(nrows*sizeof(float *)); 

  if (fits_read_col(fptr, TINT, colnum_Z, 1, 1, NZ, NULL,
		    ionbal->Z, &anynul, &status)) {
    errMsg << status << "\nread_ion_bal\nError reading element Z value\n";
    throw YellowAlert(errMsg.str());
  }
  /* Create the local index variable into the X vector */
  loc_indx[0] = 0;
  for (iZ=1;iZ<NZ;iZ++) {
    loc_indx[iZ] = loc_indx[iZ-1] + ionbal->Z[iZ-1] + 1;
  }

  /* Set up the real index variable to blow up if asked for an illegal Z */
  /* value, but otherwise to be correct. */

  for (iZ=0;iZ<MAXIMUM_Z+1;iZ++) ionbal->indx[iZ] = -1000;
  for (iZ=0;iZ<NZ;iZ++) {
    ionbal->indx[ionbal->Z[iZ]] = loc_indx[iZ];
    if (iZ == 0) {
	ionbal->MaxIndx = ionbal->Z[iZ] + 1;
    } else {
	ionbal->MaxIndx = ionbal->indx[ionbal->Z[iZ]] + ionbal->Z[iZ] + 1;
    }
  }
  ionbal->X = (float *) malloc(ionbal->MaxIndx*sizeof(float));

  for (iTemp=0;iTemp<Ntemp;iTemp++) {
    for (iDens=0;iDens<Ndens;iDens++) {
      row = iTemp*Ndens + iDens + 1;
      
      if (iDens==0) {
	if (fits_read_col(fptr, TFLOAT, colnum_temp, row, 1, 1, NULL,
			  &(ionbal->T[iTemp]), &anynul, &status)) {
	  errMsg << status << "\nread_ion_bal\nError reading ion pop values\n";
	  throw YellowAlert(errMsg.str());
	}
      }
      
      if (iTemp==0) {
	if (fits_read_col(fptr, TFLOAT, colnum_dens, row, 1, 1, NULL,
			  &(ionbal->dens[iDens]), &anynul, &status)) {
	  errMsg << status << "\nread_ion_bal\nError reading ion pop values\n";
	  throw YellowAlert(errMsg.str());
	}
      }
      
      (ionbal->tableX)[row-1] = (float *) malloc(NX*sizeof(float));
      if (fits_read_col(fptr, TFLOAT, colnum_X, row, 1, NX, NULL,
			(ionbal->tableX)[row-1], &anynul, &status)) {
	errMsg << status << "\nread_ion_bal\nError reading ion pop values\n";
	throw YellowAlert(errMsg.str());
      }
      for(i=0;i<NX;i++) {
	if ((ionbal->tableX)[row-1][i] == 0) {
	  (ionbal->tableX)[row-1][i] = ION_POP_SET_MIN;
	}
      }
    }      
  }

  if (fits_close_file(fptr, &status)) {
    errMsg << status << "\nread_ion_bal\nError closing fits file\n";
    throw YellowAlert(errMsg.str());
  }

}
