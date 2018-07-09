#define MAXIMUM_Z 100
#define MAXSTRLEN 1024
#define ION_POP_SET_MIN 1e-20
#define ION_POP_MIN 1e-10
#define BLOWUP 1e30
#ifndef KBOLTZ
#define KBOLTZ 8.617e-8
#endif
/* TWOPHOT is a bad idea, but we'll do it for now. */
#define TWOPHOT -1.0
#define HC_IN_KEV_A 12.398
#define ERG_KEV 1.60219e-9

#define TWOPH_MAXZ 28

/* #define DIRECTORY "/data/plato/rsmith/pse/projects/CX/acx_v0.3.1/data" */

#define NUMCXMODELS 8

struct ION_BAL_TABLE {
  int Nrows;    /* Total number of rows */
  int Ntemp;    /* Number of temperatures */
  int Ndens;    /* Number of densities */
  int Ntime;    /* Number of timesteps */
  float *T;     /* temperature (K) array */
  float *dens;  /* density (cm^-3) array */
  float *time;  /* time (s) array */
  int NZ;       /* actual number of elements in table */
  int MaxIndx;  /* maximum number of values in X */
  int Z[MAXIMUM_Z+1];    /* protons (Z) in i-th element */
  int indx[MAXIMUM_Z+1]; /* Index into X array */
  float *X;    /* ionization balance (unitless, normed to unity) */
  float **tableX; /* Table of ionization balances (unitless, normed to unity)*/
};

struct CX_TABLE {
  int Nrows;    /* Total number of rows */
  int *Z;    /* Element */
  int *rmJ;  /* ion, roman notation -- stripped = Z+1 */
  int *ul;   /* Upper level */
  int *ll;   /* Lower level */
  float *lambda;    /* Energy of line (in keV) */
  float *energy;    /* Energy of line (in keV) */
  float *feps;  /* scaled strength - 'fake' epsilon */
};
