#ifndef GRHAYL_IO_H_
#define GRHAYL_IO_H_

#include <stdio.h>
#include <stdlib.h>

void grhayl_info(const char *format, ...);

void grhayl_Warn_Error(
      const char *type,
      const int exit_code,
      const char *filename,
      const int line,
      const char *funcname,
      const char *format,
      ...);

typedef enum  {
  grhayl_error_abort=-1, grhayl_success=0
} grhayl_error_keys;

#define grhayl_warn(format, ...) \
  grhayl_Warn_Error("Warning", grhayl_success, __FILE__, __LINE__, __func__, format __VA_OPT__(,) __VA_ARGS__)

#define grhayl_error(format, ...) \
  grhayl_Warn_Error("Error", 1, __FILE__, __LINE__, __func__, format __VA_OPT__(,) __VA_ARGS__)

#define grhayl_abort(format, ...) \
  grhayl_Warn_Error("Error", grhayl_error_abort, __FILE__, __LINE__, __func__, format __VA_OPT__(,) __VA_ARGS__)

#define grhayl_Error(exit_code, format, ...) \
  grhayl_Warn_Error("Error", exit_code, __FILE__, __LINE__, __func__, format __VA_OPT__(,) __VA_ARGS__)

#endif // GRHAYL_IO_H_
