#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void error(const char *format, ...) {
  va_list args;
  va_start(args, format);
  fprintf(stderr, "error: ");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
  va_end(args);
  exit(EXIT_FAILURE);
}

void warning(const char *format, ...) {
  va_list args;
  va_start(args, format);
  fprintf(stderr, "warning: ");
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
  va_end(args);
}
