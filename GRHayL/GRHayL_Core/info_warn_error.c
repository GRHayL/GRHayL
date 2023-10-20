#include "ghl_io.h"
#include <stdarg.h>

/// @brief Print information message: (GRHayL) <message>
/// @param format Input format string
/// @param ... Additional parameters for \p format string
void ghl_info(const char *format, ...) {

  printf("(GRHayL) ");

  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

/// @brief Print warning or error message.
/// @param type Type of message, generally "Error" or "Warning"
/// @param exit_code Exit code; if \p ghl_success returns with no error
/// @param filename Name of the file the warning/error was found
/// @param line Line in file where the warning/error was found
/// @param funcname Name of the function where the warning/error was found
/// @param format Format string for the warning/error message
/// @param ... Additional parameters for \p format string
void ghl_Warn_Error(
    const char *type,
    const int exit_code,
    const char *filename,
    const int line,
    const char *funcname,
    const char *format,
    ...) {

  fprintf(
      stderr, "(GRHayL) %s in file: %s, line: %d, function: %s\n", type, filename, line, funcname);
  fprintf(stderr, "(GRHayL) %s message: ", type);

  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);

  switch(exit_code) {
    case ghl_success:
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
