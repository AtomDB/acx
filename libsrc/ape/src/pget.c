/** \file pget.c
    \brief
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_binary.h"
#include "ape/ape_error.h"

int main(int argc, char ** argv) {
  return eOK != ape_binary_pget(argc - 1, argv + 1) ? 1 : 0;
}

/*
 * $Log: pget.c,v $
 * Revision 1.1.1.1  2013/05/22 21:11:09  rsmith
 * Initial import
 *
 * Revision 1.1.1.1  2007/07/26 17:18:18  peachey
 * Import Ape version 2.1.1.
 *
 * Revision 1.3  2006/06/08 02:25:52  peachey
 * Rationalize pget, plist and pquery2 binaries to handle par files
 * the same way, or to use binary names and PFILES to find the parameter file(s).
 *
 * Revision 1.2  2006/05/19 17:44:15  peachey
 * Replace body of pget with call to ape_binary_pget.
 *
*/
