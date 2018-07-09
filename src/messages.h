#define MESSAGES_DEBUG
void message(const char *routine, const char *fmt, ...);
void errmess(const char *routine, const char *fmt, ...);
void quitmess(const char *routine, const char *fmt, ...);
void memmess(const char *routine);
void fitsmess(const char *routine, int status, const char *fmt, ...);
char * malloc_or_die(size_t size, char *routine);

#define malloc_safe(d) (malloc_or_die(d , __FILE__))
void apec_set_abort(void (*func)(void));
void apec_set_exit(void (*func)(int));
