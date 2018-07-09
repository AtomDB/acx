#ifndef dacx_h
#define dacx_h

#ifdef __cplusplus
extern "C" {
#endif

#define VERSION "0.1"
  
#define MAXLINELENGTH 132
#ifndef MAXSTRLEN
#define MAXSTRLEN 1024
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef PI
#define PI 3.1415926535
#endif
#ifndef M_E
#define M_E 2.7182818284590452354
#endif
#ifndef DEG2RAD
#define DEG2RAD 0.017453293
#endif
  
#define BOLTZ 11605.
#ifndef KBOLTZ
#define KBOLTZ 8.617385e-8  /* kB = [keV/K] */
#endif
#define RYDBERG 0.013605804
#define ERG_KEV 1.60219e-9
#define HC_IN_KEV_A 12.398
#define HC_IN_ERG_A 1.9848e-8   /* = hc (erg A) = 12400 * 1.602e-12 */

/* ****************************************************** */
/* PARAMETERS contains values that are set at the outset  */
/* of the run, and do not change afterward.               */
/* ****************************************************** */

struct PARAMETERS {       /* structure to collect all the plasma parameters */
  char OutputFileName[MAXSTRLEN];
  int Wavelength;
  double LambdaMin;
  double LambdaMax;
  double MinEmiss;
  double kT;
  double FracHe;
  double Abundance;
  int Model;
  double Norm;
  double redshift;
  int swcx;
  double AbC;
  double AbN;
  double AbO;
  double AbNe;
  double AbMg;
  double AbAl;
  double AbSi;
  double AbS;
  double AbAr;
  double AbCa;
  double AbFe;
  double AbNi;
  int clobber;
};

typedef struct Param_Table {
  char * name;
  unsigned int type;
  void * value;
} Param_Table_Type;

typedef void Param_File_Type;

void dacx_par_init(int argc, char *argv[], Param_File_Type *p_input,
		   struct PARAMETERS *params);

#ifdef __cplusplus
}
#endif

#endif
