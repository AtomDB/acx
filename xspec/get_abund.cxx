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

using namespace XSutility;
using namespace std;

void get_abund(float **abZ) {

  float AG89[29] = {0.0, 12.0, 10.99, 1.160, 1.150, 2.6,
		    8.56, 8.050, 8.930, 4.56, 8.090,
		    6.330, 7.58, 6.47, 7.55, 5.45,
		    7.20, 5.5, 6.56, 5.12, 6.36, 
		    3.1, 4.99, 4.0, 5.67, 5.39,
		    7.67, 4.92, 6.25};

  for (size_t ii=1;ii<29; ii++) (*abZ)[ii] = pow(10, AG89[ii] - 12.0);

}
