/** \file ape_trad.c
    \brief Implementation of traditional XPI/PIL-compliant interface.
    \author James Peachey, HEASARC/EUD/GSFC.
*/
#include "ape/ape_error.h"
#include "ape/ape_io.h"
#include "ape/ape_list.h"
#include "ape/ape_msg.h"
#include "ape/ape_par.h"
#include "ape/ape_test.h"
#include "ape/ape_trad.h"
#include "ape/ape_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
static const char s_dir_separator[] = { '/', '\\' };
#define QC ";"
#define QS "|"
#define QD "\\"
#define QD2 "/"
#else
static const char s_dir_separator[] = { '/' };
#define QC ":"
#define QS ";"
#define QD "/"
#define QD2 "/"
#endif

enum { eSystem, eLocal, eMerged, eCurrent };

static ApeParFile * s_default_par_file[] = { 0, 0, 0, 0 };

static void ape_trad_destroy(void) {
  ape_io_destroy_file(s_default_par_file[eCurrent]); s_default_par_file[eCurrent] = 0;
  ape_io_destroy_file(s_default_par_file[eMerged]); s_default_par_file[eMerged] = 0;
  ape_io_destroy_file(s_default_par_file[eLocal]); s_default_par_file[eLocal] = 0;
  ape_io_destroy_file(s_default_par_file[eSystem]); s_default_par_file[eSystem] = 0;
}

static void ape_trad_atexit(void) {
  ape_trad_close(0);
}

static int find_par(const char * par_name, ApePar ** par) {
  int status = eOK;
  ApeParFile * par_file = s_default_par_file[eCurrent];
  ApeListIterator itor = 0;

  if (0 != par) *par = 0;
  else status = eNullPointer;

  if (eOK == status && 0 == par_name) status = eNullPointer;

  if (eOK == status && 0 != par_file) {
    /* Get the named parameter. */
    status = ape_io_find_par(par_name, par_file, &itor);
  } else {
    status = eNoParLoaded;
  }

  if (eOK == status) {
    *par = (ApePar *) ape_list_get(itor);
  }

  return status;
}

/* Internal utility which takes the full path to the executable as an argument, and looks for the
   last directory separator (think slash). Everything after the last separator is considered the
   name of the tool. */
static int get_pfile_name(const char * full_exec_name, char ** pfile_name) {
  size_t num_sep = sizeof(s_dir_separator) / sizeof(const char);
  size_t idx = 0;
  size_t last_idx = 0;
  const char * last_sep = 0;
  char * exec_name = 0;
  int status = eOK;

  /* Iterate over valid separators. */
  for (idx = 0; idx != num_sep; ++idx) {
    /* Look for right-most separator. */
    const char * sep = strrchr(full_exec_name, s_dir_separator[idx]);

    /* In case multiple separators are found, pick the right-most one. */
    if (sep > last_sep) {
      last_idx = idx;
      last_sep = sep;
    }
  }

  /* If a separator was found, start after it, otherwise just use the whole file name. */
  if (0 != last_sep) last_sep += sizeof(s_dir_separator[last_idx]) / sizeof(const char);
  else last_sep = full_exec_name;

  /* Make a copy of the file name part of the full name. */
  status = ape_util_copy_string(last_sep, &exec_name);
  if (0 == status) {
    /* Terminate the copy at the last period, in order to truncate any extension (e.g. .exe) before adding .par. */
    char * extension = strrchr(exec_name, '.');
    if (0 != extension) *extension = '\0';

    /* Add the .par in lieu of the extension. */
    status = ape_util_cat_string(exec_name, ".par", pfile_name);
  }
  /* Clean up. */
  free(exec_name);
  return status;
}

int ape_trad_init(int argc, char ** argv) {
  int status = eOK;
  char * pfile_name = 0;

  /* Destroy whatever parameter objects are currently open, saving them first. */
  ape_trad_close(1);

  /* Register cleanup function. */
  ape_util_atexit(&ape_trad_atexit);

  /* Check command line arguments for problems. */
  if (0 == argv) {
    status = eNullPointer;
    ape_msg_debug("ape_trad_init was called with argv == 0.\n");
  } else if (0 == argv[0]) {
    status = eNullPointer;
    ape_msg_debug("ape_trad_init was called with argv[0] == 0.\n");
  } else if (0 == *argv[0]) {
    status = eInvalidArgument;
    ape_msg_debug("ape_trad_init was called with argv[0] == \"\".\n");
  }
  if (0 >= argc) {
    status = eInvalidArgument;
    ape_msg_debug("ape_trad_init was called with argc <= 0.\n");
  }

  if (eOK == status) {
    /* Get the name of the parameter file, based on this executable name. */
    status = get_pfile_name(argv[0], &pfile_name);
  }

  if (eOK == status) {
    ApeList * loc_path = 0;
    ApeList * sys_path = 0;

    /* Get the system par file path, if any. */
    status = ape_io_get_sys_path(&sys_path);
    if (eOK == status) {
      /* Ignore the status because it's OK if this fails. */
      ape_io_read_file_path(pfile_name, sys_path, &s_default_par_file[eSystem]);
    }

    /* Get the local par file path, if any. */
    status = ape_io_get_loc_path(&loc_path);
    if (eOK == status) {
      /* Try to open a parameter file on the local path. */
      status = ape_io_read_file_path(pfile_name, loc_path, &s_default_par_file[eLocal]);
      if (eOK != status) {
        ApeListIterator itor = ape_list_begin(loc_path);
        if (itor != ape_list_end(loc_path)) {
          char * full_file_name = 0;

          /* Local path contains at least one directory, so clone the system file. */
          status = ape_io_clone_file(s_default_par_file[eSystem], s_default_par_file + eLocal);

          if (eOK == status) {
            const char * dir_name = (const char *) ape_list_get(itor);
            /* Get the name of the local parameter file. */
            status = ape_util_append_file_name(dir_name, pfile_name, &full_file_name);
          }
 
          if (eOK == status) {
            /* Ignore the status because it's OK if this fails. */
            ape_io_set_file_name(s_default_par_file[eLocal], full_file_name);
          }
          free(full_file_name); full_file_name = 0;
        }
      }
    }
    /* Some amount of error above may be tolerated. If something unrecoverable happened it will be detected below. */
    status = eOK;

    if (0 == s_default_par_file[eLocal] && 0 == s_default_par_file[eSystem]) {
      /* No luck searching the paths: try the parameter file name by itself. */
      status = ape_io_read(pfile_name, &s_default_par_file[eLocal]);
    }
  }

  if (eOK == status) {
    /* Merge the system and local parameter files to create a file which has the best set of parameter values. */
    status = ape_io_merge_par_files(s_default_par_file[eSystem], s_default_par_file[eLocal], &s_default_par_file[eMerged]);
  }

  if (eOK == status) {
    /* Clone the merged parameter file to create the "current" parameter file. */
    status = ape_io_clone_file(s_default_par_file[eMerged], &s_default_par_file[eCurrent]);
  }

  if (eOK == status) {
    /* Modify the merged parameter file using the command line arguments. */
    status = ape_io_apply_command_line(s_default_par_file[eCurrent], argc - 1, argv + 1);
  }

  if (eOK == status) {
    /* Check the content of the parameter file to ensure its format is valid. */
    status = ape_io_check_file_format(s_default_par_file[eCurrent], 0);
  }

  /* Clean up. */
  free(pfile_name); pfile_name = 0;

  return status;
}

int ape_trad_close(int save_flag) {
  if (0 != save_flag) ape_trad_save();
  ape_trad_destroy();
  return eOK;
}

int ape_trad_get_current(ApeParFile ** par_file) {
  int status = eOK;

  if (0 != par_file) {
    *par_file = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    if (0 != s_default_par_file[eCurrent]) {
      *par_file = s_default_par_file[eCurrent];
    } else {
      status = eUninitialized;
    }
  }

  return status;
}

int ape_trad_get_sys_pfile(ApeParFile ** par_file) {
  int status = eOK;

  if (0 != par_file) {
    *par_file = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    if (0 != s_default_par_file[eSystem]) {
      *par_file = s_default_par_file[eSystem];
    } else {
      status = eUninitialized;
    }
  }

  return status;
}

int ape_trad_get_par_names(char *** par_name) {
  int status = eOK;
  ApeList * par_cont = 0;
  ApeListIterator itor = 0;
  ApeListIterator end = 0;
  size_t num_par = 0;
  size_t idx = 0;

  if (0 != par_name) {
    *par_name = 0;
  } else {
    status = eNullPointer;
  }

  if (eOK == status && 0 == s_default_par_file[eCurrent]) {
    status = eUninitialized;
  }

  if (eOK == status) {
    status = ape_io_get_par_cont(s_default_par_file[eCurrent], &par_cont);
  }

  if (eOK == status) {
    num_par = ape_list_get_size(par_cont);
  }

  if (eOK == status) {
    *par_name = (char **) calloc(num_par + 1, sizeof(char *));
  }

  if (eOK == status) {
    end = ape_list_end(par_cont);
  }

  for (itor = ape_list_begin(par_cont); itor != end; itor = ape_list_next(itor)) {
    char * name = 0;
    ApePar * par = (ApePar *) ape_list_get(itor);
    int local_status = ape_par_get_field(par, eName, &name);
    if (eOK == local_status) {
      if (0 == name || '\0' == *name) {
        free(name); name = 0;
      } else {
        (*par_name)[idx] = name;
        ++idx;
      }
    } else {
      free(name); name = 0;
      status = eOK == status ? local_status : status;
    }
  }

  return status;   
}

int ape_trad_get_type(const char * par_name, char * par_type) {
  int status = eOK;
  ApePar * par = 0;

  if (0 != par_type) *par_type = 0;
  else status = eNullPointer;

  if (eOK == status) status = find_par(par_name, &par);

  if (eOK == status) status = ape_par_get_type(par, par_type);

  return status;
}

int ape_trad_get_value_string(const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  if (0 != value) *value = 0;
  else status = eNullPointer;

  if (eOK == status) status = find_par(par_name, &par);

  if (eOK == status) status = ape_par_get_string(par, value);

  return status;
}

int ape_trad_get_value_string_case(const char * par_name, char ** value, char case_code) {
  int status = eOK;
  ApePar * par = 0;

  if (0 != value) *value = 0;
  else status = eNullPointer;

  if (eOK == status) status = find_par(par_name, &par);

  if (eOK == status) status = ape_par_get_string_case(par, value, case_code);

  return status;
}

static int query_value(const char * par_name, ApePar ** par) {
  int status = eOK;
  ApeListIterator itor = 0;
  ApeParFile * par_file = s_default_par_file[eCurrent];
  char default_mode[APE_PAR_MODE_CODE_LEN] = "";

  /* Check arguments. */
  if (0 != par) {
    *par = '\0';
  } else if (eOK == status) {
    status = eNullPointer;
  }
  if (eOK == status && 0 == par_name) status = eNullPointer;
  
  if (eOK == status && 0 != par_file) {
    /* Find the named parameter in the merged (working) parameter file. */
    status = ape_io_find_par(par_name, par_file, &itor);
  } else {
    status = eNoParLoaded;
  }

  if (eOK == status) {
    /* Get the default mode from the parameter file. */
    status = ape_io_get_default_mode(par_file, default_mode);

    if (eOK != status) {
      /* No default mode defined, so use ql. */
      status = eOK;
      strcpy(default_mode, "ql");
    }
  }

  if (eOK == status) {
    *par = (ApePar *) ape_list_get(itor);

    /* Query for the parameter. */
    status = ape_par_query(*par, default_mode);

    /* Ignore status; try to save the current parameters here. */
    ape_trad_save();
  }

  return status;
}

int ape_trad_query_bool(const char * par_name, char * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_bool(par, value);
  }

  return status;
}

int ape_trad_query_double(const char * par_name, double * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_double(par, value);
  }

  return status;
}

int ape_trad_query_float(const char * par_name, float * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_float(par, value);
  }

  return status;
}

int ape_trad_query_file_name(const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_file_name(par, value);
  }

  return status;
}

int ape_trad_query_int(const char * par_name, int * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_int(par, value);
  }

  return status;
}

int ape_trad_query_long(const char * par_name, long * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_long(par, value);
  }

  return status;
}

int ape_trad_query_string(const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_string(par, value);
  }

  return status;
}

int ape_trad_query_string_case(const char * par_name, char ** value, char case_code) {
  int status = eOK;
  ApePar * par = 0;

  /* Check output parameter. */
  if (0 != value) {
    *value = '\0';
  } else {
    status = eNullPointer;
  }

  if (eOK == status) {
    /* Find parameter, prompt for it and return pointer to it. */
    status = query_value(par_name, &par);
  }

  if (eOK == status) {
    status = ape_par_get_string_case(par, value, case_code);
  }

  return status;
}

int ape_trad_get_bool(const char * par_name, char * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_bool(par, value);

  /* Check parameter for validity, include the value. */
  if (eOK == status) status = ape_par_check(par, 1);

  return status;
}

int ape_trad_get_double(const char * par_name, double * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_double(par, value);

  /* Check parameter for validity, include the value. */
  if (eOK == status) status = ape_par_check(par, 1);

  return status;
}

int ape_trad_get_float(const char * par_name, float * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_float(par, value);

  /* Check parameter for validity, include the value. */
  if (eOK == status) status = ape_par_check(par, 1);

  return status;
}

int ape_trad_get_file_name(const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_file_name(par, value);

  /* Check parameter for validity, include the value. */
  if (eOK == status) status = ape_par_check(par, 1);

  return status;
}

int ape_trad_get_long(const char * par_name, long * value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_long(par, value);

  /* Check parameter for validity, include the value. */
  if (eOK == status) status = ape_par_check(par, 1);

  return status;
}

int ape_trad_get_string(const char * par_name, char ** value) {
  int status = eOK;
  ApePar * par = 0;

  /* Find the given parameter. */
  status = find_par(par_name, &par);

  /* Get the value from the parameter. */
  if (eOK == status) status = ape_par_get_string(par, value);

  /* Check parameter for validity, include the value. */
  if (eOK == status) status = ape_par_check(par, 1);

  return status;
}

int ape_trad_set_bool(const char * par_name, char value) {
  ApePar * par = 0;
  int status = find_par(par_name, &par);

  if (eOK == status) {
    if (0 != value) status = ape_par_set_value_string(par, "yes");
    else status = ape_par_set_value_string(par, "no");
  }

  return status;
}

int ape_trad_set_double(const char * par_name, double value) {
  ApePar * par = 0;
  int status = find_par(par_name, &par);

  if (eOK == status) {
    char buf[128] = "";
    sprintf(buf, "%1.15g", value);
    status = ape_par_set_value_string(par, buf);
  }

  return status;
}

int ape_trad_set_float(const char * par_name, float value) {
  ApePar * par = 0;
  int status = find_par(par_name, &par);

  if (eOK == status) {
    char buf[128] = "";
    sprintf(buf, "%1.7g", value);
    status = ape_par_set_value_string(par, buf);
  }

  return status;
}

int ape_trad_set_file_name(const char * par_name, const char * value) {
  ApePar * par = 0;
  int status = find_par(par_name, &par);

  if (eOK == status) {
    status = ape_par_set_value_string(par, value);
  }

  return status;
}

int ape_trad_set_long(const char * par_name, long value) {
  ApePar * par = 0;
  int status = find_par(par_name, &par);

  if (eOK == status) {
    char buf[128] = "";
    sprintf(buf, "%ld", value);
    status = ape_par_set_value_string(par, buf);
  }

  return status;
}

int ape_trad_set_string(const char * par_name, const char * value) {
  ApePar * par = 0;
  int status = find_par(par_name, &par);

  if (eOK == status) {
    status = ape_par_set_value_string(par, value);
  }

  return status;
}

int ape_trad_save(void) {
  int status = eOK;
  ApeParFile * tmp_file = 0;

  /* Make a temporary copy of current parameters. */
  status = ape_io_clone_file(s_default_par_file[eCurrent], &tmp_file);

  if (eOK == status) {
    /* Reset unlearned parameters in the copy to their pre-command-line state. */
    status = ape_io_revert_unlearned(tmp_file, s_default_par_file[eMerged]);
  }

  if (eOK == status) {
    /* Write temporary copy. */
    status = ape_io_write(tmp_file, 0);
  }

  /* Clean up. */
  ape_io_destroy_file(tmp_file); tmp_file = 0;

  return status;
}

int ape_trad_find_par(const char * par_name, ApePar ** par) {
  return find_par(par_name, par);
}

void ape_trad_test(void) {
  int status = eOK;

  /* Test get_pfile_name. */
  /* TODO: test more thoroughly. */
  { char * pfile_name = 0;
    const char * expected_pfile_name = "binary.par";
    status = get_pfile_name("c:/some/path/to\\a/binary.exe", &pfile_name);
    if (0 == pfile_name) {
      ape_test_failed("ape_trad_test: get_pfile_name(long path) returned pfile_name == 0, not \"%s\" as expected.\n",
        expected_pfile_name);
    } else if (0 != strcmp(pfile_name, expected_pfile_name)) {
      ape_test_failed("ape_trad_test: get_pfile_name(long path) returned \"%s\", not \"%s\" as expected.\n",
        pfile_name, expected_pfile_name);
    }
    free(pfile_name); pfile_name = 0;
  }
  { char * pfile_name = 0;
    const char * expected_pfile_name = "binary.par";
    status = get_pfile_name("binary.exe", &pfile_name);
    if (0 == pfile_name) {
      ape_test_failed("ape_trad_test: get_pfile_name(short path) returned pfile_name == 0, not \"%s\" as expected.\n",
        expected_pfile_name);
    } else if (0 != strcmp(pfile_name, expected_pfile_name)) {
      ape_test_failed("ape_trad_test: get_pfile_name(short path) returned \"%s\", not \"%s\" as expected.\n",
        pfile_name, expected_pfile_name);
    }
    free(pfile_name); pfile_name = 0;
  }

  /* Attempt to initialize a task using invalid inputs. */
  /* Null argc and argv. */
  { int argc = 0;
    char ** argv = 0;
    status = ape_trad_init(argc, argv);
    if (eInvalidArgument != status) {
      ape_test_failed("ape_trad_init(0, 0) returned status %d, not %d as expected.\n", status, eInvalidArgument);
    }
  }
  /* Null argv. */
  { int argc = 1;
    char ** argv = 0;
    status = ape_trad_init(argc, argv);
    if (eNullPointer != status) {
      ape_test_failed("ape_trad_init(1, 0) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
  }
  /* Null argv[0]. */
  { int argc = 1;
    char * argv[] = { 0 };
    status = ape_trad_init(argc, argv);
    if (eNullPointer != status) {
      ape_test_failed("ape_trad_init(1, argv == { 0 }) returned status %d, not %d as expected.\n", status, eNullPointer);
    }
  }
  /* Empty argv[0]. */
  { int argc = 1;
    char * argv[] = { "" };
    status = ape_trad_init(argc, argv);
    if (eInvalidArgument != status) {
      ape_test_failed("ape_trad_init(1, argv == \"\") returned status %d, not %d as expected.\n", status, eInvalidArgument);
    }
  }

  /* Non-existent parameter file. */
  { int argc = 1;
    char * argv[] = { "non-existent-task" };
    status = ape_trad_init(argc, argv);
    if (eFileReadError != status) {
      ape_test_failed("ape_trad_init(1, argv == \"%s\") returned status %d, not %d as expected.\n",
        argv[0], status, eInvalidArgument);
    }
  }

  /* Existing parameter file. */
  { char * argv[] = {
/* Note this #ifdef appears backwards, but this is deliberate, to cross-up the parameter files to ensure
   that parameter files written on Windows can be read on Unix and vice versa. */
#ifdef WIN32
      "ape_test_unix",
#else
      "ape_test_win32",
#endif
      "sa",
      "sal",
      "sq",
      "sql",
      "bh",
      "=",
      "yes",
      "fh=.",
      0
    };
    const char * name[] = { "sa", "sal", "sq", "sql", "bh", "fh" };
    const char * expected[] = { "sa", "sal", "sq", "sql", "yes", "." };
    int argc = sizeof(argv) / sizeof(char *) - 1;
    int idx = 0;
    char * pfiles_orig = 0;
    status = ape_io_get_pfiles(&pfiles_orig);
    if (eOK == status) {
      status = ape_io_set_pfiles(QS".");
    }
    if (eOK == status) {
      status = ape_trad_init(argc, argv);
      if (eOK == status) {
        for (idx = 0; idx != sizeof(name) / sizeof(char *); ++idx) {
          ApeListIterator itor = 0;
          status = ape_io_find_par(name[idx], s_default_par_file[eCurrent], &itor);
          if (eOK == status) {
            ApePar * par = (ApePar *) ape_list_get(itor);
            char * value = 0;
            char buf[128];
            status = ape_par_get_field(par, eValue, &value);
            sprintf(buf, "after ape_trad_init(8 args) parameter %s", name[idx]);
            ape_test_cmp_string(buf, value, expected[idx], status, eOK);
            free(value); value = 0;
          }
        }
        /* Tests for ape_trad_query_bool. */
        { char value = '\0';
          status = ape_trad_query_bool("bh", &value);
          ape_test_cmp_long("ape_trad_query_bool(\"bh\")", value, '\1', status, eOK);
        }
        { char value = '\1';
          status = ape_trad_query_bool("nobool", &value);
          ape_test_cmp_long("ape_trad_query_bool(\"nobool\")", value, '\0', status, eParNotFound);
        }
        { char value = '\1';
          status = ape_trad_query_bool("sa", &value);
          ape_test_cmp_long("ape_trad_query_bool(\"sa\")", value, '\0', status, eConversionError);
        }
        /* Tests for ape_trad_query_double and ape_trad_set_double. */
        { double value = 1.;
          char msg[128] = "ape_trad_query_double(\"dh\", &value)";
          status = ape_trad_query_double("dh", &value);
          ape_test_cmp_double(msg, value, 0., status, eUndefinedValue);

          value = 1.23456789012345e200;
          sprintf(msg, "ape_trad_set_double(\"dh\", %1.15g)", value);
          status = ape_trad_set_double("dh", value);
          if (eOK == status) {
            char * string_value = 0;
            char expected_string_value[64] = "";
            double expected_value = value;
            sprintf(expected_string_value, "%1.15g", expected_value);
            status = ape_trad_get_string("dh", &string_value);
            ape_test_cmp_string("ape_trad_get_string(\"dh\")", string_value, expected_string_value, status, eOK);
            free(string_value); string_value = 0;
            value = 1.;
            status = ape_trad_get_double("dh", &value);
            ape_test_cmp_double("ape_trad_get_double(\"dh\")", value, expected_value, status, eOK);
          } else {
            ape_test_failed("%s returned status %d, not %d, as expected.\n", msg, status, eOK);
          }
        }
        /* Tests for ape_trad_query_long and ape_trad_set_long. */
        { long value = 0;
          status = ape_trad_query_long("ih", &value);
          ape_test_cmp_long("ape_trad_query_long(\"ih\")", value, 0, status, eUndefinedValue);
        }
        { long value = 1000000001l;
          status = ape_trad_set_long("ih", value);
          ape_test_cmp_long("ape_trad_set_long(\"ih\", 1000000001l)", value, value, status, eOK);
        }
        if (eOK == status) {
          char * value = 0;
          status = ape_trad_get_string("ih", &value);
          ape_test_cmp_string("ape_trad_get_string(\"ih\")", value, "1000000001", status, eOK);
          free(value); value = 0;
        }
        /* Test ape_trad_get_par_names. */
        { const char * expected[] = {
            "sa", "sal", "sh", "shl", "sq",
            "sql", "sh_comment", "bh", "dh", "fh",
            "frh", "fwh", "ih", "rh", "bhmin",
            "dhmin", "fhmin", "frhmin", "fwhmin", "ihmin",
            "rhmin", "shmin", "bhenum", "dhenum", "fhenum",
            "frhenum", "fwhenum", "ihenum", "rhenum", "shenum",
            "bhmax", "dhmax", "fhmax", "frhmax", "fwhmax",
            "ihmax", "rhmax", "shmax", "bhrange", "dhrange",
            "fhrange", "frhrange", "fwhrange", "ihrange", "rhrange",
            "shrange", "bhvalid", "dhvalid", "fhvalid", "frhvalid",
            "fwhvalid", "ihvalid", "rhvalid", "shvalid", "bhlow",
            "dhlow", "fhlow", "frhlow", "fwhlow", "ihlow",
            "rhlow", "shlow", "bhhigh", "dhhigh", "fhhigh",
            "frhhigh", "fwhhigh", "ihhigh", "rhhigh", "shhigh",
            "dhinvalid", "ihinvalid", "rhinvalid", "ihindef", "ihinf",
            "ihinfinity", "ihnan", "ihnone", "ihundef", "ihundefined",
            "rhindef", "rhinf", "rhinfinity", "rhnan", "rhnone",
            "rhundef", "rhundefined", "mode",
            0
          };
          char ** par_name = 0;
          status = ape_trad_get_par_names(&par_name);
          ape_test_cmp_string_array("ape_trad_par_get_names", par_name, expected, status, eOK);
          ape_util_free_string_array(par_name); par_name = 0;
        }
        /* Test ape_trad_get_type. */
        { const char * par_name = "frh";
          char par_type[APE_PAR_TYPE_CODE_LEN] = "qZX";
          const char * expected_type = "fr";
          int expected_status = eOK;
          status = ape_trad_get_type(par_name, par_type);
          ape_test_cmp_string("ape_trad_get_type(\"frh\")", par_type, expected_type, status, expected_status);

          par_name = "non-existent-par";
          expected_type = "";
          expected_status = eParNotFound;
          status = ape_trad_get_type(par_name, par_type);
          ape_test_cmp_string("ape_trad_get_type(\"non-existent-par\")", par_type, expected_type, status, expected_status);
        }
      } else {
        ape_test_failed("ape_trad_init(argv == \"%s\", ...) returned status %d, not %d as expected.\n",
          argv[0], status, eOK);
      }
      ape_trad_close(0);
      ape_io_set_pfiles(pfiles_orig);
    } else {
      ape_test_failed("Could not set up to test ape_trad_init (status is %d.)\n", status);
    }
    free(pfiles_orig); pfiles_orig = 0;
  }
}

#ifdef __cplusplus
}
#endif

/*
 * $Log: ape_trad.c,v $
 * Revision 1.1.1.1  2013/05/22 21:11:09  rsmith
 * Initial import
 *
 * Revision 1.2  2007/11/02 15:22:39  peachey
 * Import Ape v2r2.
 *
 * Revision 1.37  2007/05/15 20:52:48  peachey
 * Additional speed optimizations:
 * o Reduce number of calls to ape_list_end by calling it once and storing
 *   the returned iterator for use in loop logic when it is safe to do so.
 * o Make ape_list_destroy more efficient by iterating through the list
 *   back to front once and freeing the elements instead of repeatedly
 *   calling ape_list_remove_entry.
 *
 * Revision 1.36  2007/01/25 17:49:57  peachey
 * Fixed bug in ape_par_get_string_case that caused it not to behave
 * correctly when called for a non-enumerated parameter. Added functions
 * ape_trad_prompt_string_case and ape_trad_get_string_case for getting
 * parameters with control over case of returned string.
 *
 * Revision 1.35  2006/11/30 20:38:35  peachey
 * Add support for direct getting of integer values.
 *
 * Revision 1.34  2006/11/30 16:40:27  peachey
 * Fix bug in test code: pfiles needs to be set correctly for Windows test.
 *
 * Revision 1.33  2006/11/30 16:28:47  peachey
 * Fix bug in test code: pfiles needs to be set correctly for Windows test.
 *
 * Revision 1.32  2006/11/24 19:47:28  peachey
 * Add ape_trad_query_float, ape_trad_get_float and ape_trad_set_float.
 *
 * Revision 1.31  2006/08/25 19:45:58  peachey
 * Add tests for undefined and infinite numerical parameters.
 *
 * Revision 1.30  2006/06/06 13:30:21  peachey
 * Add ape_trad_find_par, for getting underlying ape parameter object.
 *
 * Revision 1.29  2006/06/05 18:59:34  peachey
 * Add apd_trad_get_type for getting the type code of a parameter.
 *
 * Revision 1.28  2006/05/31 03:53:00  peachey
 * Add body of ape_trad_query_file_name.
 *
 * Revision 1.27  2006/05/31 03:21:44  peachey
 * Rationalize ape_par_get_* family of functions, and use them to implement
 * ape_trad_get_* and ape_trad_query_*.
 *
 * Revision 1.26  2006/05/31 01:36:51  peachey
 * Rename ape_par_get_string to ape_par_get_field and ape_par_get_string_array to
 * ape_par_get_field_array.
 *
 * Revision 1.25  2006/05/26 15:22:11  peachey
 * Add ape_trad_set* functions complete ape_trad_query* functions, for completeness
 * and to support PIL compatibility.
 *
 * Revision 1.24  2006/05/23 19:31:40  peachey
 * Add ape_trad_get_sys_pfile, for getting the system parameter file.
 *
 * Revision 1.23  2006/05/23 16:27:25  peachey
 * Update call to ape_io_write to add new force_write argument.
 *
 * Revision 1.22  2006/05/22 17:35:40  peachey
 * Add ape_trad_query_string, for prompting for and getting string parameters.
 *
 * Revision 1.21  2006/05/18 14:00:50  peachey
 * Rename ape_trad_get_bool ape_trad_query_bool. Add ape_trad_get_string
 * which only gets string, without prompting.
 *
 * Revision 1.20  2006/05/18 03:18:54  peachey
 * Correct handling of situation in which local pfiles is not empty
 * but no file exists in that location yet.
 *
 * Revision 1.19  2006/05/17 03:15:30  peachey
 * Add ape_trad_get_current, for getting the current parameter file object.
 *
 * Revision 1.18  2006/05/17 02:20:24  peachey
 * Add ape_trad_get_par_names, for getting all named parameters in a file.
 *
 * Revision 1.17  2006/05/16 19:47:22  peachey
 * Set pfiles to ;. to prevent local file from being modified before
 * testing ape_trad_init. Consolidate ape_trad_init tests.
 *
 * Revision 1.16  2006/05/16 17:22:27  peachey
 * Update ape_trad tests to use latest, greatest ape_test.par.
 *
 * Revision 1.15  2006/05/16 15:13:21  peachey
 * Add note reminding of need to write better test for ape_trad_init.
 *
 * Revision 1.14  2006/05/16 15:09:46  peachey
 * Add test for ape_trad_init in which the long and complex ape_test.par is opened.
 *
 * Revision 1.13  2006/05/13 18:45:56  peachey
 * Add ape_trad_save and use it to save parameter file after each prompt.
 *
 * Revision 1.12  2006/05/12 17:33:21  peachey
 * Add fourth default parameter file and apply command line arguments to that;
 * preserve pre-command line version as 'merged' version.
 *
 * Revision 1.11  2006/05/12 03:32:09  peachey
 * Add ape_trad_get_bool, and infrastructure for all the ape_trad_get_... family.
 *
 * Revision 1.10  2006/05/10 16:48:50  peachey
 * Generalize format tests to allow them to be run with or without
 * value checking, to let them be used both in contexts where the whole parameter
 * needs to be checked (checking user input) and contexts in which a bad value
 * is not a show-stopper, because the user or client may yet correct the problem.
 *
 * Revision 1.9  2006/05/09 19:32:08  peachey
 * In ape_trad_init, check the merged parameter file using ape_io_check_file_format.
 *
 * Revision 1.8  2006/05/06 00:43:16  peachey
 * Add most of implemention and tests for ape_trad_init, which is
 * modeled after PILInit.
 *
 * Revision 1.7  2006/04/28 01:24:04  peachey
 * Add a couple casts to get the constness right when calling free.
 * Enable al mode, and remove eFileOpenError, using eFileReadError in lieu thereof.
 * The latter change was to make error messages consistent between Windows and
 * Unix, because although the two both fail to open a directory as if it were
 * a file, they fail at different points in the process of opening the file.
 *
 * Revision 1.6  2006/04/19 16:04:20  peachey
 * Add rational shutdown/cleanup/atexit stuff. Public function
 * ape_trad_close performs all shutdown/cleanup. Static function ape_trad_destroy
 * does all deallocations. Static function ape_trad_atexit is registered with
 * ape_util_atexit, and calls ape_trad_close. System and local par files are now
 * kept in static storage by ape_trad_init.
 *
 * Revision 1.5  2006/04/19 03:07:09  peachey
 * Uncomment final effort to open local parameter file, and correct
 * the expected status.
 *
 * Revision 1.4  2006/04/15 01:56:29  peachey
 * Revamp ape_trad_init using recent developments in ape_io.c to get local
 * and system paths, and from them to load local and system parameter
 * files.
 *
 * Revision 1.3  2006/04/13 18:44:44  peachey
 * Continue to chip away at ape_trad_init, adding get_pfile_name helper function.
 *
 * Revision 1.2  2006/04/12 14:29:29  peachey
 * Add forgotten header file.
 *
 * Revision 1.1  2006/04/12 14:18:14  peachey
 * Add beginning of traditional facility for duplicating PIL/XPI behavior.
 *
*/
