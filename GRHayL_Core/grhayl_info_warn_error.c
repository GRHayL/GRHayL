#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void grhayl_info(const char *format, ...) {

  printf("(GRHayL) ");
  
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

void grhayl_Warn_Error(
      const char *type,
      const int exit_code,
      const char *filename,
      const int line,
      const char *funcname,
      const char *format,
      va_list args) {

  fprintf(stderr, "(GRHayL) %s in file: %s, line: %d, function: %s\n", type, filename, line, funcname);
  fprintf(stderr, "(GRHayL) %s message: ", type);

  vfprintf(stderr, format, args);
  va_end(args);

  if( exit_code )
    exit(exit_code);
}

void grhayl_Warn(
      const char *filename,
      const int line,
      const char *funcname,
      const char *format,
      ...) {
  va_list args;
  va_start(args, format);
  grhayl_Warn_Error("Warning", 0, filename, line, funcname, format, args);
}

void grhayl_Error(
      const int exit_code,
      const char *filename,
      const int line,
      const char *funcname,
      const char *format,
      ...) {
  va_list args;
  va_start(args, format);
  grhayl_Warn_Error("Error", exit_code, filename, line, funcname, format, args);
}
