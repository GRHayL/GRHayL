#ifndef GRHAYL_IO_H_
#define GRHAYL_IO_H_

void grhayl_info(const char *format, ...);

void grhayl_Warn(
      const char *filename,
      const int line,
      const char *funcname,
      const char *format,
      ...);

void grhayl_Error(
      const int exit_code,
      const char *filename,
      const int line,
      const char *funcname,
      const char *format,
      ...);

#define grhayl_warn(format, ...) grhayl_Warn(__FILE__, __LINE__, __func__, format __VA_OPT__(,) __VA_ARGS__)
#define grhayl_error(exit_code, format, ...) grhayl_Error(exit_code, __FILE__, __LINE__, __func__, format __VA_OPT__(,) __VA_ARGS__)

#endif // GRHAYL_IO_H_
