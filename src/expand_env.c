#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "messages.h"
#include "dacx.h"

void expand_env(char *filename) {

  /* Now, at this point, filename could be a fully-qualified filename.*/
  /* Or it could have an environment variable in it. */
  /* We need to distinguish between these, and do some checks. */

  int i;
  char *p1, *p2, EnvName[MAXSTRLEN], EnvStr[MAXSTRLEN];
  char *result; 

  result = (char *) malloc(MAXSTRLEN*sizeof(char));
  if (strchr(filename,'$') != NULL) { /* Environment var */
    p1 = strchr(filename,'$');
    i=0;
	while (0 != isalnum(p1[i+1])) {
      EnvName[i] = p1[i+1];
      i++;
    }
    EnvName[i] = '\0';
    p2 = getenv(EnvName);
    if (p2 == NULL) { 
      errmess("expand_env","No environment variable $%s",EnvName);
    } else {
      char * first_slash = 0;
      for (i=0;i<strlen(p2)+1;i++) EnvStr[i] = p2[i];
      first_slash = strchr(filename,'/');
	  if (0 == first_slash) first_slash = strchr(filename, '\\');
	  if (0 != first_slash) {
        strcat(EnvStr,first_slash);
        for (i=0;i<strlen(EnvStr)+1;i++) filename[i] = EnvStr[i];
      }
	}
  }
#ifdef WIN32
  for (i=0;i<strlen(filename);++i) {
    /* Rationalize path on Windows to contain all backslashes. */
    if ('/' == filename[i]) filename[i] = '\\';
  }
#endif
}
