/** \file ape_test.c
    \brief Implementation of unit test functions. Module level tests are implemented elsewhere, in order to
    be able to access internal data.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_error.h"
#include "ape/ape_msg.h"
#include "ape/ape_test.h"

#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

static int s_status = 0;

int ape_test(void) {
  /* Handle for output file. */
  FILE * out_file = 0;

  /* Enable debug mode. */
  ape_msg_debug_enable(1);

  /* Open output file to contain all unit test messages for comparison to further tests. */
  out_file = fopen("ape_test-log", "w");
  if (0 == out_file) ape_test_failed("Cannot open file \"%s\" for output.\n", "ape_test-log");

  /* Reset streams to go to this file. */
  ape_msg_set_err_stream(out_file);
  ape_msg_set_out_stream(out_file);

  /* Test i/o facility. */
  ape_io_test();

  /* Test pil-compatible interface immediately after test file was set up. */
  ape_pil_test();

  /* Test list facility. */
  ape_list_test();

  /* Test message facility. */
  ape_msg_test();

  /* Test parameter facility. */
  ape_par_test();

  /* Test parameter group facility. */
  ape_par_group_test();

  /* Test internal utilities. */
  ape_util_test();

  /* Test traditional interface. */
  ape_trad_test();

  /* Test binary support interface. */
  ape_binary_test();

  /* Redirect errors back to stderr for final messages. */
  ape_msg_set_out_stream(stdout);
  ape_msg_set_err_stream(stderr);

  /* Close log file. */
  fclose(out_file);

  /* Issue final status message. */
  if (eOK == s_status) {
    ape_msg_info(0, "APE Unit Test Passed!\n");
  } else {
    ape_msg_error("APE Unit Test FAILED!\n");
  }
  return s_status;
}

void ape_test_cmp_double(const char * hint, double result, double expected_result, int status, int expected_status) {
  double tolerance = DBL_EPSILON * 10.;
  const char * join = ", ";
  if (0 == hint) {
    hint = "";
    join = "";
  }
  if (result != expected_result && (
    (0. != result && (fabs((result - expected_result) / result) > tolerance)) ||
    (0. != expected_result && (fabs((result - expected_result) / expected_result) > tolerance))
  )) {
    ape_test_failed("%s%slong result was %g, not %g, as expected.\n", hint, join, result, expected_result);
  }
  if (status != expected_status)
    ape_test_failed("%s%sstatus was %d, not %d, as expected.\n", hint, join, status, expected_status);
}

void ape_test_cmp_long(const char * hint, long result, long expected_result, int status, int expected_status) {
  const char * join = ", ";
  if (0 == hint) {
    hint = "";
    join = "";
  }
  if (result != expected_result) {
    ape_test_failed("%s%slong result was %ld, not %ld, as expected.\n", hint, join, result, expected_result);
  }
  if (status != expected_status)
    ape_test_failed("%s%sstatus was %d, not %d, as expected.\n", hint, join, status, expected_status);
}

void ape_test_cmp_ptr(const char * hint, void * result, void * expected_result, int status, int expected_status) {
  const char * join = ", ";
  if (0 == hint) {
    hint = "";
    join = "";
  }
  if (result != expected_result) {
    ape_test_failed("%s%spointer result was %p, not %p, as expected.\n", hint, join, result, expected_result);
  }
  if (status != expected_status)
    ape_test_failed("%s%sstatus was %d, not %d, as expected.\n", hint, join, status, expected_status);
}

void ape_test_cmp_string(const char * hint, const char * result, const char * expected_result, int status, int expected_status) {
  const char * join = ", ";
  if (0 == hint) {
    hint = "";
    join = "";
  }
  if (result != expected_result) {
    if (result == 0)
      ape_test_failed("%s%sstring result was 0, not \"%s\", as expected.\n", hint, join, expected_result);
    else if (expected_result == 0)
      ape_test_failed("%s%sstring result was \"%s\", not 0, as expected.\n", hint, join, result);
    else if (0 != strcmp(result, expected_result))
      ape_test_failed("%s%sstring result was \"%s\", not \"%s\", as expected.\n", hint, join, result, expected_result);
  }
  if (status != expected_status)
    ape_test_failed("%s%sstatus was %d, not %d, as expected.\n", hint, join, status, expected_status);
}

void ape_test_cmp_string_array(const char * hint, char ** result, const char ** expected_result,
  int status, int expected_status) {
  const char * join = ", ";
  if (0 == hint) {
    hint = "";
    join = "";
  }
  if ((void *) result != (void *) expected_result) {
    if (result == 0)
      ape_test_failed("%s%sstring array result was unexpectedly 0.\n", hint, join);
    else if (expected_result == 0)
      ape_test_failed("%s%sstring array result was unexpectedly not 0.\n", hint, join);
    else {
      for (; 0 != *result && 0 != *expected_result; ++result, ++expected_result) {
        /* Compare individual strings using ape_test_cmp_string, but don't use actual status;
           report that just once after all the strings are compared. */
        ape_test_cmp_string(hint, *result, *expected_result, expected_status, expected_status);
      }
      if (0 != *result) {
        ape_test_failed("%s%sstring array result had more elements than expected. The extra elements are:\n", hint, join);
        for (; 0 != *result; ++result) ape_test_failed("%s%s\t\"%s\"\n", hint, join, *result);
      }
      if (0 != *expected_result) {
        ape_test_failed("%s%sstring array result had fewer elements than expected. The missing elements are:\n", hint, join);
        for (; 0 != *expected_result; ++expected_result) ape_test_failed("%s%s\t\"%s\"\n", hint, join, *expected_result);
      }
    }
  }
  if (status != expected_status)
    ape_test_failed("%s%sstatus was %d, not %d, as expected.\n", hint, join, status, expected_status);
}

void ape_test_cmp_ulong(const char * hint, unsigned long result, unsigned long expected_result,
  int status, int expected_status) {
  const char * join = ", ";
  if (0 == hint) {
    hint = "";
    join = "";
  }
  if (result != expected_result) {
    ape_test_failed("%s%slong result was %lu, not %lu, as expected.\n", hint, join, result, expected_result);
  }
  if (status != expected_status)
    ape_test_failed("%s%sstatus was %d, not %d, as expected.\n", hint, join, status, expected_status);
}

void ape_test_set_status(int status) {
  /* Set the given status provided an error has not already been flagged. */
  if (eOK == s_status) s_status = status;
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_test.c,v $
 * Revision 1.1.1.1  2013/05/22 21:11:09  rsmith
 * Initial import
 *
 * Revision 1.1.1.1  2007/07/26 17:18:18  peachey
 * Import Ape version 2.1.1.
 *
 * Revision 1.12  2006/06/06 13:28:35  peachey
 * Add ape_pil_test.
 *
 * Revision 1.11  2006/05/19 17:31:30  peachey
 * Add ape_binary_test, and redirect output stream of unit test to same
 * log file where error goes.
 *
 * Revision 1.10  2006/05/01 19:12:33  peachey
 * Make ape_test_cmp_string_array give more specific error messages.
 *
 * Revision 1.9  2006/04/26 14:31:10  peachey
 * Add ape_test_cmp_string_array for comparing arrays of strings.
 *
 * Revision 1.8  2006/04/20 03:43:24  peachey
 * Add ape_test_cmp_ptr for comparing pointers in tests.
 *
 * Revision 1.7  2006/04/14 14:15:23  peachey
 * Add ape_test_cmp_ulong for comparing unsigned results.
 *
 * Revision 1.6  2006/04/13 18:42:25  peachey
 * Add ape_test_cmp_* family of functions to facilitate smaller and more
 * consistent tests.
 *
 * Revision 1.5  2006/04/12 17:57:00  peachey
 * Add final success/failure message.
 *
 * Revision 1.4  2006/04/12 14:26:37  peachey
 * Test traditional interface.
 *
 * Revision 1.3  2006/04/12 01:46:38  peachey
 * Close the log file at the end of the test.
 *
 * Revision 1.2  2006/04/10 21:12:18  peachey
 * Add test for ape_util module.
 *
 * Revision 1.1.1.1  2006/04/05 13:45:19  peachey
 * Initial import of All-purpose Parameter Environment (APE).
 *
*/
