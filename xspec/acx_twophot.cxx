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

/* **************************************************************** */
/* **************************************************************** */
/*                     Utility functions below here                 */
/* **************************************************************** */
/* **************************************************************** */


double A_E_twoph(double A, double E0, double E) {

  /* From Nussbaumer & Schmutz, 1984, A+A, 138,495 */
  /* Z is the element, and E is the energy of the bin, in keV */
  /* y is unitless, and is equal to nu/nu0 = lambda0/lambda, where */
  /* lambda0 = 1215.7 A for hydrogen--the base wavelength of the 2s->1s */
  /* transition.  This fit is accurate to better than 0.6% for */
  /* 0.01 < y < 0.99 */

  /* The A_norm is the A value for neutral hydrogen for this transition. */
  /* For other transitions, we renormalize to the appropriate A value. */

  /* This routine is used for BOTH hydrogenic and He-like two-photon */
  /* distributions.  This is justified using the result of */
  /* Derevianko & Johnson, 1997, Phys Rev A, 56, 1288 who show in */
  /* Figures 5 and 2 of that paper that the difference is everywhere */
  /* less than 10% between these two for Z=6-28 -- it is about 5% or so. */

  double C     = 202.0;  /* s^-1 */
  double alpha = 0.88;
  double beta  = 1.53;
  double gamma = 0.80;
  double A_norm= 8.2249; /* in s^-1 */
  double x, y, z;
  double result;
  
  y = E/E0;

  x = y*(1-y);
  z = pow(4*x,gamma);
  result = C * ( x*(1-z) + alpha*pow(x,beta)*z );
  
  result *= (A/A_norm);  /* Also need R_Z/R_H, but even for Z=26, this is */
                         /* only 1.0005, so we'll ignore it. */
  return result;  /* in same units as input 'A' value */
}

double twophot_lambda(int Z, int rmJ) {

  /* He-like 2-photon transitions (1s2s 1S0  ---> 1s2_1S0) From Drake */
  /* G.W Drake, 1982, Can.J.Phys. Vol 66, 1988 */

  double L_D82[TWOPH_MAXZ+1] = {0.0, 
				0.0, 601.40824, 0.0, 0.0, 0.0,
			  40.732571, 29.076157, 21.793820, 0.0, 13.545514, 
			  11.076951, 9.2261891, 7.8029406, 6.6849856, 0.0, 
			  5.0644832, 0.0, 3.9681272, 0.0, 3.1920069, 
			  0.0, 0.0, 0.0, 0.0, 0.0,
				1.8593994, 0.0, 1.59403};

  /* H-like 2-photon transitions (2s 1S0 -> 1s 1S0) */
  /* Primarily from Erickson 1977, JPCRD, 6, 831 */

  double L_E77[TWOPH_MAXZ+1] = {0.0, 
			  0.0, 0.0, 0.0, 0.0, 0.0,
			  33.7393, 24.7843, 18.9723, 14.9874, 12.1373,
			  10.0266, 8.42440, 7.17620, 6.18564, 5.38524,
			  4.73259, 4.18938, 3.73635, 3.35087, 3.02374,
			  2.74031, 2.49554, 2.28199, 2.09458, 1.92920,
				1.78332, 1.65187, 1.53565};
  
  double result = -1.0;
  
  if ((Z < 1)||(Z > TWOPH_MAXZ)) return result;

  if (Z == rmJ) result = L_E77[Z];
  if (Z == rmJ+1) result = L_D82[Z];
   
  return result;
  
}

