#include <stdarg.h>
#include "ghl_io.h"

void ghl_info(
      const char *format,
      ...) {

  printf("(GRHayL) ");

  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

void ghl_Warn_Error(
      const char *type,
      const int exit_code,
      const char *filename,
      const int line,
      const char *funcname,
      const char *format,
      ...) {

  fprintf(stderr, "(GRHayL) %s in file: %s, line: %d, function: %s\n", type, filename, line, funcname);
  fprintf(stderr, "(GRHayL) %s message: ", type);

  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);

  switch (exit_code) {
  case ghl_io_success:
    return;
    break;
  case ghl_error_abort:
    abort();
    break;
  default:
    exit(exit_code);
    break;
  }
}
