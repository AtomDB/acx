#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "fitsio.h"
#include "dacx.h"
#include "acx.h"

struct ION_BAL_TABLE *ionbal;

int Ncxh, NcxhOLD;
int Ncxhlines, NcxhlinesOLD;
struct CX_TABLE cxh, cxhOLD;

int Ncxhe, NcxheOLD;
int Ncxhelines, NcxhelinesOLD;
struct CX_TABLE cxhe, cxheOLD;

void fits_read_cx_data(const char *file, int *Ncxlines, int *Ncx, 
		       struct CX_TABLE *cx);
void read_ion_bal(const char *ionbal_file);
void find_ion_bal(int Z, double Te, double Ne, double Ion_pop[]);

void calc_cx(struct PARAMETERS *params, const int model, 
	     const double fracH, const double fracHe, 
	     const double z, float *vabZ, double **ionpopZ) {
    
  int iCX;
  int iNeut;
  int Ncxl=0;
  int Ncxtot=0;
  int Actualmodel=0;
  double fracNeutral=0.0;
  double enz;
  double lineflux;
  double value;
  int done;
  struct CX_TABLE *cxm=NULL;
  FILE *fp=NULL;

  const double Emin = HC_IN_KEV_A/params->LambdaMax;
  const double Emax = HC_IN_KEV_A/params->LambdaMin;
  
  /* ************************** */
  /* Open output file if needed */
  /* ************************** */
  
  if (strncmp("STDOUT.dat",params->OutputFileName,MAXSTRLEN) != 0) {
    fp = fopen(params->OutputFileName,"w");
  }
  
  /* ************************** */
  /* neutral loop: H, then He   */
  /* ************************** */

  for(iNeut=0;iNeut<2;iNeut++) {

    iCX = 0;
    done = FALSE;

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

    /* Assuming energy-sorted list of CX lines, step through until we */ 
    /* find one that is inside our lower bound */
    while ((!done)&&(enz < Emin)) {  
      iCX++;
      enz = cxm->energy[iCX]/(1.+z);
      if (iCX == Ncxtot) done = TRUE;
    }
    
    /* iCX is now cued up to a photon of interest */ 
    while ((!done)&&(enz < Emax)) { // calc emission
      /* Note that ionpopZ[cxm.Z[iCX]][0] is for the neutral state, rmJ=1 */
      /* So ionpopZ[cxm.Z[iCX]][cxm.rmJ[iCX]] is not for Z, rmJ but rather */
      /* for the next ion stage up, as it should be */
	
      lineflux = 1e5 * params->Norm * fracNeutral * vabZ[cxm->Z[iCX]] *
	(ionpopZ[cxm->Z[iCX]][cxm->rmJ[iCX]]) * 
	cxm->feps[iCX*NUMCXMODELS + Actualmodel-1];

      if (lineflux > params->MinEmiss) {
	if (iCX < Ncxl) {
	  value = enz;
	  if (params->Wavelength) value = HC_IN_KEV_A/value;
	  if (fp == NULL) {
	    printf("%2d %3d %3d %5d %5d %f %e\n",iNeut, 
		   cxm->Z[iCX], cxm->rmJ[iCX], cxm->ul[iCX], cxm->ll[iCX], 
		   value, lineflux);
	  } else {
	    fprintf(fp, "%2d %3d %3d %5d %5d %f %e\n",iNeut, 
		   cxm->Z[iCX], cxm->rmJ[iCX], cxm->ul[iCX], cxm->ll[iCX], 
		   value, lineflux);
	  }
	} else { /* Two photon */
	  if (fp == NULL) {
	    printf("%2d %3d %3d %5d %5d 0.0 %e\n",iNeut, 
		   cxm->Z[iCX], cxm->rmJ[iCX], cxm->ul[iCX], cxm->ll[iCX], 
		   lineflux);
	  } else {
	    fprintf(fp, "%2d %3d %3d %5d %5d 0.0 %e\n",iNeut, 
		    cxm->Z[iCX], cxm->rmJ[iCX], cxm->ul[iCX], cxm->ll[iCX], 
		    lineflux);
	  }
	}
      }
      iCX++;
      if (iCX == Ncxtot) {
	done = TRUE;
      } else {
	enz = cxm->energy[iCX]/(1.+z);
      }
    }
  }

  if (fp != NULL) fclose(fp);

}

void acx_init(int *ZhasCX) {

  int ii, iCX;
  char *ibfile, *cxfile, *cxhefile, *cxfileOLD, *cxhefileOLD;
  int Ncx, Ncxlines;

  ibfile      = (char *) malloc(MAXSTRLEN*sizeof(char));
  cxfile      = (char *) malloc(MAXSTRLEN*sizeof(char));
  cxhefile    = (char *) malloc(MAXSTRLEN*sizeof(char));
  cxfileOLD   = (char *) malloc(MAXSTRLEN*sizeof(char));
  cxhefileOLD = (char *) malloc(MAXSTRLEN*sizeof(char));

  sprintf(ibfile,"%s%s",DIRECTORY,"/data/v2.0.2_ionbal.fits");
  sprintf(cxfile,"%s%s",DIRECTORY,"/data/AtomDB_CX_H_v1.0.0.fits");
  sprintf(cxhefile,"%s%s",DIRECTORY,"/data/AtomDB_CX_He_v1.0.0.fits");
  sprintf(cxfileOLD,"%s%s",DIRECTORY,"/data/AtomDB_CX_H_v0.5.1.fits");
  sprintf(cxhefileOLD,"%s%s",DIRECTORY,"/data/AtomDB_CX_He_v0.5.1.fits");

  read_ion_bal(ibfile);

  fits_read_cx_data(cxfile, &Ncxlines, &Ncx, &cxh);
  Ncxh = Ncx;
  Ncxhlines = Ncxlines;

  fits_read_cx_data(cxhefile, &Ncxlines, &Ncx, &cxhe);
  Ncxhe = Ncx;
  Ncxhelines = Ncxlines;
  
  fits_read_cx_data(cxfileOLD, &Ncxlines, &Ncx, &cxhOLD);
  NcxhOLD = Ncx;
  NcxhlinesOLD = Ncxlines;

  fits_read_cx_data(cxhefileOLD, &Ncxlines, &Ncx, &cxheOLD);
  NcxheOLD = Ncx;
  NcxhelinesOLD = Ncxlines;
  
  for (iCX=0;iCX<Ncxh;iCX++) ZhasCX[cxh.Z[iCX]] = TRUE;
  for (iCX=0;iCX<Ncxhe;iCX++) ZhasCX[cxhe.Z[iCX]] = TRUE;
  for (iCX=0;iCX<NcxhOLD;iCX++) ZhasCX[cxhOLD.Z[iCX]] = TRUE;
  for (iCX=0;iCX<NcxheOLD;iCX++) ZhasCX[cxheOLD.Z[iCX]] = TRUE;

}

void fits_read_cx_data(const char *file, int *Ncxlines, int *Ncx, 
		       struct CX_TABLE *cx) {

  fitsfile *fptr;
  int hdunum, casesen;
  
  int ii;
  int colnum_Z, colnum_rmJ, colnum_UL, colnum_LL, colnum_Lambda, colnum_Prob;
  char column[MAXSTRLEN],comment[MAXSTRLEN];
  int status, hdutype, anynul=0, nulval=0;
  int nrows;

  status = 0;

  if (fits_open_file(&fptr, file, READONLY, &status)) {
    message("fits_read_cx_data","Failed to open %s",file);
  }

  hdunum = 2;
  if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) {
    message("fits_read_cx_data","Error finding second HDU");
  }
    
  casesen = FALSE;
  
  fits_get_colname(fptr, casesen, "ELEMENT", column, &colnum_Z, &status);
  fits_get_colname(fptr, casesen, "ION", column,     &colnum_rmJ, &status);
  fits_get_colname(fptr, casesen, "UPPERLEV", column,&colnum_UL, &status);
  fits_get_colname(fptr, casesen, "LOWERLEV", column,&colnum_LL, &status);
  fits_get_colname(fptr, casesen, "WAVELEN", column, &colnum_Lambda, &status);
  fits_get_colname(fptr, casesen, "LINEPROB", column,&colnum_Prob, &status);

  if (status) {
    message("fits_read_cx_data","Error getting column names");
  }
  
  if (fits_read_key(fptr, TINT, "NAXIS2", &nrows, comment, &status)) {
    message("fits_read_cx_data","Error getting number of rows");
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
    message("fits_read_cx_data","Error reading element Z value");
  }

  if (fits_read_col(fptr, TINT, colnum_rmJ, 1, 1, nrows, &nulval,
		    cx->rmJ, &anynul, &status)) {
    message("fits_read_cx_data","Error reading ion rmJ value");
  }
  
  if (fits_read_col(fptr, TINT, colnum_UL, 1, 1, nrows, &nulval,
		    cx->ul, &anynul, &status)) {
    message("fits_read_cx_data","Error reading upper level value");
  }
  
  if (fits_read_col(fptr, TINT, colnum_LL, 1, 1, nrows, &nulval,
		    cx->ll, &anynul, &status)) {
    message("fits_read_cx_data","Error reading lower level value");
  }
  
  if (fits_read_col(fptr, TFLOAT, colnum_Lambda, 1, 1, nrows, &nulval,
		    cx->lambda, &anynul, &status)) {
    message("fits_read_cx_data","Error reading wavelength value");
  }
  
  if (fits_read_col(fptr, TFLOAT, colnum_Prob, 1, 1, nrows*NUMCXMODELS,&nulval,
		    cx->feps, &anynul, &status)) {
    message("fits_read_cx_data","Error reading line prob value");
  }
  
  if (fits_close_file(fptr, &status)) {
    message("fits_read_cx_data","Error closing fits file");
  }
  /* 
  printf("%d %d %d %d  %e :  %e %e %e %e %e %e %e %e\n",
	 cx->Z[1], cx->rmJ[1], cx->ul[1], cx->ll[1], 
	 cx->lambda[1], 
	 cx->feps[8+0],cx->feps[8+1],cx->feps[8+2],cx->feps[8+3],
	 cx->feps[8+4],cx->feps[8+5],cx->feps[8+6],cx->feps[8+7]);
  */
  /* Now work out actual number of lines. */
  for (ii=0;ii<nrows;ii++) {
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
  
  NT = ionbal->Ntemp;
  ND = ionbal->Ndens;
  
  /* printf("find_ion_bal: NT = %d, ND = %d\n",NT, ND); */
  
  /* Now, find where we are...*/
  for (i=0;i<NT;i++) {  /* Assumes T is monotonic increasing */
    if (((ionbal->T)[i] > Te)&&(iTlow==-1)) {
      if (i==0) {
	errmess("find_ion_bal","Te=%e out of range for table, min val=%e.",
		Te, (ionbal->T)[i]);
      } else {
	iTlow = i-1;
	iThigh = i;
	dT = (ionbal->T)[iThigh] - (ionbal->T)[iTlow];
	w_Tlow = 1 - (Te - (ionbal->T)[iTlow])/dT;
	w_Thigh= 1 - w_Tlow;
      }
    }
  }
  if (iTlow==-1) errmess("find_ion_bal","Te=%e out of range of table.",Te);
  
  if (ND > 1) {
    for (i=0;i<ND;i++) {  /* Assumes N is monotonic increasing */
      if (((ionbal->dens)[i] > Ne)&&(iNlow==-1)) {
	if (i==0) {
	  errmess("find_ion_bal","Ne=%e out of range for table, min val=%e",
		  Ne,(ionbal->dens)[i]);
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

  status = 0;

  if (fits_open_file(&fptr, ionbal_file, READONLY, &status)) {
    message("read_ion_bal","Error opening fits file %s",ionbal_file);
    /* fits_error(status); */
  }

  hdunum = 2;
  if (fits_movabs_hdu(fptr, hdunum, &hdutype, &status)) {
    message("read_ion_bal","Error finding second HDU");
    /* fits_error(status); */
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
    message("read_ion_bal","Error getting column names");
    /* fits_error(status); */
  }
  
  if (fits_read_key(fptr, TINT, "NAXIS2", &nrows, comment, &status)) {
    message("read_ion_bal","Error getting number of rows.");
    /* fits_error(status); */
  }
  if (fits_read_key(fptr, TINT, "T_NUMBER", &Ntemp, comment, &status)) {
    message("read_ion_bal","Error getting number of temperatures.");
    /* fits_error(status); */
  }
  if (fits_read_key(fptr, TLOGICAL, "T_LINEAR", &Tlinear, comment, &status)) {
    message("read_ion_bal","Error getting temperature logical.");
    /* fits_error(status); */
  }
  if (fits_read_key(fptr, TINT, "N_NUMBER", &Ndens, comment, &status)) {
    message("read_ion_bal","Error getting number of densities.");
    /* fits_error(status);*/
  }
  if (fits_read_key(fptr, TLOGICAL, "N_LINEAR", &Nlinear, comment, &status)) {
    message("read_ion_bal","Error getting density logical.");
    /* fits_error(status); */
  }
  if (fits_read_key(fptr, TINT, "N_ELEMEN", &NZ, comment, &status)) {
    message("read_ion_bal","Error getting number of elements.");
    /* fits_error(status); */
  }
  if (fits_read_key(fptr, TINT, "N_IONS", &NX, comment, &status)) {
    message("read_ion_bal","Error getting number of ions.");
    /* fits_error(status); */
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
    errmess("read_ion_bal",
	    "Number of rows %d != Number of temperatures (%d) * Number of densities (%d)",
	    nrows,Ntemp,Ndens);
  }

  /* ionbal->tableX = (float **) malloc(Ntemp*Ndens*sizeof(float)); */
  ionbal->tableX = (float **) malloc(nrows*sizeof(float *)); 

  if (fits_read_col(fptr, TINT, colnum_Z, 1, 1, NZ, NULL,
		    ionbal->Z, &anynul, &status)) {
    message("read_ion_bal","Error reading element Z values.");
    /* fits_error(status);*/
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
	  message("read_ion_bal","Error reading ion population values.");
	  /* fits_error(status); */
	}
      }
      
      if (iTemp==0) {
	if (fits_read_col(fptr, TFLOAT, colnum_dens, row, 1, 1, NULL,
			  &(ionbal->dens[iDens]), &anynul, &status)) {
	  message("read_ion_bal","Error reading ion population values.");
	  /* fits_error(status); */
	}
      }
      
      (ionbal->tableX)[row-1] = (float *) malloc(NX*sizeof(float));
      if (fits_read_col(fptr, TFLOAT, colnum_X, row, 1, NX, NULL,
			(ionbal->tableX)[row-1], &anynul, &status)) {
	message("read_ion_bal","Error reading ion population values.");
	/* fits_error(status);*/
      }
      for(i=0;i<NX;i++) {
	if ((ionbal->tableX)[row-1][i] == 0) {
	  (ionbal->tableX)[row-1][i] = ION_POP_SET_MIN;
	}
      }
    }      
  }

  if (fits_close_file(fptr, &status)) {
    fitsmess("read_ion_bal",status,"Error closing fits file.");
  }

}
