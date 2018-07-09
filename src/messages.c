#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "messages.h"

static void (*s_apec_abort_func)(void) = &abort;
static void (*s_apec_exit_func)(int) = &exit;

void message(const char *routine, const char *fmt, ...) {
  
  va_list ap;

  va_start(ap, fmt);  /* Point ap to first message. */
  
#ifdef MESSAGES_DEBUG
  printf("%s: ",routine);
  vprintf(fmt,ap);
  va_end(ap);
  printf("\n");
  fflush(stdout);
#endif

}


void errmess(const char *routine, const char *fmt, ...) {

  va_list ap;

  va_start(ap, fmt);
  printf("Error reported in: %s\n",routine);
  vprintf(fmt,ap);
  printf("\n Bombing out. \n");
  fflush(stdout);
  s_apec_abort_func();
}

void quitmess(const char *routine, const char *fmt, ...) {

  va_list ap;

  va_start(ap, fmt);
  printf("Error reported in: %s\n",routine);
  vprintf(fmt,ap);
  printf("\n Exiting. \n");
  fflush(stdout);
  s_apec_exit_func(1);
}


void memmess(const char *routine) {

  printf("Memory allocation error in: %s\n",routine);
  printf("\n Bombing out. \n");
  fflush(stdout);
  s_apec_abort_func();
}

void fitsmess(const char *routine, int status, const char *fmt, ...) {
  
  va_list ap;

  va_start(ap, fmt);  /* Point ap to first message. */
  
  printf("FITS error number #%d\n", status);
  printf("%s: ",routine);
  vprintf(fmt,ap);
  va_end(ap);
  printf("\n");
  fflush(stdout);
}


char * malloc_or_die(size_t size, char *routine) {

  char *ptr;

  ptr = (char *) malloc(size);
  if (ptr == NULL) {
    memmess(routine);
  }
  return ptr;
}

void apec_set_abort(void (*func)(void)) {
  s_apec_abort_func = func;
}

void apec_set_exit(void (*func)(int)) {
  s_apec_exit_func = func;
}

