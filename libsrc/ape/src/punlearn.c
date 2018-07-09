/** \file punlearn.c
    \brief
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_binary.h"
#include "ape/ape_error.h"

int main(int argc, char ** argv) {
  return eOK != ape_binary_punlearn(argc - 1, argv + 1) ? 1 : 0;
}

/*
 * $Log: punlearn.c,v $
 * Revision 1.1.1.1  2013/05/22 21:11:09  rsmith
 * Initial import
 *
 * Revision 1.1.1.1  2007/07/26 17:18:18  peachey
 * Import Ape version 2.1.1.
 *
 * Revision 1.1  2006/05/23 19:35:31  peachey
 * Add punlearn executable.
 *
 */
