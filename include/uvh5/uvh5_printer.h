/* printer.h
 */
#ifndef _PRINTER_H
#define _PRINTER_H

#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>

// #define PRINTER_VERBOSE

#ifdef __cplusplus
extern "C" {
#endif

/* Call this to log an error message */
void uvh5_print_error(const char *name, const char *msg, ...);

/* Call this to log an warning message */
void uvh5_print_warn(const char *name, const char *msg, ...);

/* Call this to log an informational message */
void uvh5_print_info(const char *name, const char *msg, ...);

/* Call this to log verbose message */
void uvh5_print_verbose(const char *name, const char *msg, ...);

#ifdef __cplusplus
}
#endif

#endif // _PRINTER_H
