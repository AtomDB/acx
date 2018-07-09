#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ape/ape_error.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_trad.h"

#include "dacx.h"
#include "messages.h"

enum {
  PF_UNKNOWN_TYPE,
  PF_BOOLEAN_TYPE,
  PF_DOUBLE_TYPE,
  PF_FILE_TYPE,
  PF_FLOAT_TYPE,
  PF_INTEGER_TYPE,
  PF_STRING_TYPE
};

#include <assert.h>

#define TRUE 1
#define FALSE 0

static char *OutputFileName;
static int Wavelength;
static double LambdaMin;
static double LambdaMax;
static double MinEmiss;
static double kT;
static double FracHe;
static double Abundance;
static int Model;
static double Norm;
static double redshift;
static int swcx;
static double AbC;
static double AbN;
static double AbO;
static double AbNe;
static double AbMg;
static double AbAl;
static double AbSi;
static double AbS;
static double AbAr;
static double AbCa;
static double AbFe;
static double AbNi;
static int clobber;

static Param_Table_Type Control_Parm_Table [] = {
  {"OutputFileName", PF_FILE_TYPE,    &OutputFileName},
  {"Wavelength",     PF_BOOLEAN_TYPE, &Wavelength},
  {"LambdaMin",      PF_DOUBLE_TYPE,  &LambdaMin},     
  {"LambdaMax",      PF_DOUBLE_TYPE,  &LambdaMax},     
  {"MinEmiss",       PF_DOUBLE_TYPE,  &MinEmiss},      
  {"kT",             PF_DOUBLE_TYPE,  &kT},            
  {"FracHe",         PF_DOUBLE_TYPE,  &FracHe},
  {"Abundance",      PF_DOUBLE_TYPE,  &Abundance},     
  {"Model",          PF_INTEGER_TYPE, &Model},         
  {"Norm",           PF_DOUBLE_TYPE,  &Norm},     
  {"redshift",       PF_DOUBLE_TYPE,  &redshift},      
  {"swcx",           PF_BOOLEAN_TYPE, &swcx},      
  {"C",              PF_DOUBLE_TYPE,  &AbC},             
  {"N",              PF_DOUBLE_TYPE,  &AbN},             
  {"O",              PF_DOUBLE_TYPE,  &AbO},             
  {"Ne",             PF_DOUBLE_TYPE,  &AbNe},            
  {"Mg",             PF_DOUBLE_TYPE,  &AbMg},            
  {"Al",             PF_DOUBLE_TYPE,  &AbAl},            
  {"Si",             PF_DOUBLE_TYPE,  &AbSi},            
  {"S",              PF_DOUBLE_TYPE,  &AbS},             
  {"Ar",             PF_DOUBLE_TYPE,  &AbAr},            
  {"Ca",             PF_DOUBLE_TYPE,  &AbCa},            
  {"Fe",             PF_DOUBLE_TYPE,  &AbFe},            
  {"Ni",             PF_DOUBLE_TYPE,  &AbNi},            
  {"clobber",        PF_BOOLEAN_TYPE, &clobber},
  {NULL, 0, NULL} };

int strip_strcmp(const char *s1, const char *s2);

void expand_env(char *string);

static int get_parameters(Param_Table_Type *params);

void get_input_data(int argc, char *argv[], Param_File_Type *p_input,
		    struct PARAMETERS *params, int *status) {

  char tmpstr[MAXSTRLEN];
  char *value;

  /* Read input parameters and populate params and debug structures. */
  dacx_par_init(argc, argv, p_input, params);

  strncpy(tmpstr,OutputFileName,MAXSTRLEN);
  expand_env(tmpstr);
  sprintf(params->OutputFileName,"%s.dat",tmpstr);
  params->Wavelength=Wavelength;
  if (Wavelength) {
    params->LambdaMin =LambdaMin;
    params->LambdaMax =LambdaMax;
    if (LambdaMin <= 0) params->LambdaMin = 0.1; /* Reset to a small number */
    if (LambdaMax <= 0) params->LambdaMax = 1e8;
  } else {
    params->LambdaMin = 10.0;  /* Energy Defaults */
    params->LambdaMin = 100.0; /* Energy Defaults */
    if (LambdaMax > 0) params->LambdaMin =HC_IN_KEV_A/LambdaMax;
    if (LambdaMin > 0) params->LambdaMax =HC_IN_KEV_A/LambdaMin;
  }
  params->MinEmiss  =MinEmiss;
  params->kT        =kT;
  params->FracHe    =FracHe;
  params->Abundance =Abundance;
  params->Model     =Model;
  params->Norm      =Norm;
  params->redshift  =redshift;
  params->swcx      =swcx;
  params->AbC       =AbC;
  params->AbN       =AbN;
  params->AbO       =AbO;
  params->AbNe      =AbNe;
  params->AbMg      =AbMg;
  params->AbAl      =AbAl;
  params->AbSi      =AbSi;
  params->AbS       =AbS;
  params->AbAr      =AbAr;
  params->AbCa      =AbCa;
  params->AbFe      =AbFe;
  params->AbNi      =AbNi;      
  params->clobber = clobber;

}

int strip_strcmp(const char *s1, const char *s2) {
  
  /* Strips any whitespace from string s1, and compares it to s2. */

  char strip_string[MAXSTRLEN];

  /*  strcpy(strip_string,s1); */
  sscanf(s1,"%s",strip_string);
  return strcmp(strip_string,s2);
}

static int get_parameters(Param_Table_Type *params) {

  int status=eOK;

  for (; NULL!=params->name && eOK == status; ++params) {
    /* Handle strings a little carefully; ape allocates just enough space. */

    char * tmp_str = 0;
    if (0 != tmp_str) fprintf(stderr, "You're wrong!!!!\n");

    switch (params->type) {
      case PF_BOOLEAN_TYPE: {
        int *value = (int *) params->value;
        char bool_value='\0';
        status = ape_trad_query_bool(params->name, &bool_value);
        *value = bool_value;
        break;
      }
      case PF_DOUBLE_TYPE:
        status = ape_trad_query_double(params->name, (double *) params->value);
        break;
      case PF_FILE_TYPE:
        status = ape_trad_query_file_name(params->name, &tmp_str);
        break;
      case PF_FLOAT_TYPE:
        status = ape_trad_query_float(params->name, (float *) params->value);
        break;
      case PF_INTEGER_TYPE:
        status = ape_trad_query_int(params->name, (int *) params->value);
        break;
      case PF_STRING_TYPE:
        /* Try to get string, matching the case of enumerated parameters. */
        status = ape_trad_query_string_case(params->name, &tmp_str, eEnumCase);
        /* If the parameter wasn't enumerated, that's OK, just ignore the error code. */
        if (eRangeNoEnum == status) status = eOK;
        break;
      default:
        status = -1;
    }
    /* If tmp_str was used, copy it over to the actual parameter value. */
    if (eOK == status && 0 != tmp_str) {
      char * real_str = calloc(MAXSTRLEN, sizeof(char));
      if (0 == real_str) status = eDynAllocFailed;
      if (eOK == status) {
        strncpy(real_str, tmp_str, MAXSTRLEN - 1);
        *(char **) params->value = real_str;
      }
      free(tmp_str);
    }
  }
  return status;
}


void dacx_par_init(int argc, char *argv[], Param_File_Type *p_input,
		      struct PARAMETERS *params) {

  int status=0;
  int ii=0;

  /* ************************** */
  /* Check the input parameters */
  /* ************************** */
  
  if (FALSE) { /* Skip this, not needed. */
    if (p_input == NULL) { 
      /* Make a copy of the arguments, including name of apec par file and
	 NULL one-past-last pointer. */
      char **new_argv = (char **) calloc(argc + 2, sizeof(char *));
      
      if (NULL == new_argv) {
	errmess ("dacx_par_init", 
		 "Unable to allocate array for command line arguments.\n");
      }
      
      for (ii = 0; ii < argc; ++ii) {
	new_argv[ii] = argv[ii];
      }
      
      if (eOK != (status = ape_trad_init (argc, new_argv))) {
	free (new_argv);
	errmess ("dacx_par_init", "Error opening parameter file.\n");
      } else {
	free (new_argv);
      }
    }
  }

  if (eOK != (status = ape_trad_init (argc, argv))) {
    errmess ("dacx_par_init", "Error opening parameter file.\n");
  }

  if (eOK == status) {
    if (eOK != (status = get_parameters (Control_Parm_Table))) {
      ape_trad_close (0);
      errmess("dacx_par_init", 
	      "Error getting parameters (Ape error code %d).\n", status);
    } else {
      ape_trad_close (0);
    }
  }

}

