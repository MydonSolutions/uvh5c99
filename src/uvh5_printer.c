/* printer.c
 */
#include "uvh5/uvh5_printer.h"

/* For now just put it all to stderr.
 * Maybe do something clever like a stack in the future?
 */
void UVH5print_error(const char *name, const char *msg, ...) {
    fprintf(stderr, "Error (%s)", name);
    if(msg) {
        va_list ap;
        va_start(ap, msg);
        fprintf(stderr, ": ");
        vfprintf(stderr, msg, ap);
        va_end(ap);
    }
    if(errno) {
        fprintf(stderr, " [%s]", strerror(errno));
    }
    fprintf(stderr, "\n");
    fflush(stderr);
}

void UVH5print_warn(const char *name, const char *msg, ...) {
    fprintf(stderr, "Warning (%s)", name);
    if(msg) {
        va_list ap;
        va_start(ap, msg);
        fprintf(stderr, ": ");
        vfprintf(stderr, msg, ap);
        va_end(ap);
    }
    fprintf(stderr, "\n");
    fflush(stderr);
}

void UVH5print_info(const char *name, const char *msg, ...) {
    fprintf(stdout, "Info (%s)", name);
    if(msg) {
        va_list ap;
        va_start(ap, msg);
        fprintf(stdout, ": ");
        vfprintf(stdout, msg, ap);
        va_end(ap);
    }
    fprintf(stdout, "\n");
    fflush(stdout);
}

void UVH5print_verbose(const char *name, const char *msg, ...) {
#ifdef PRINTER_VERBOSE
    fprintf(stdout, "Verbage (%s)", name);
    if(msg) {
        va_list ap;
        va_start(ap, msg);
        fprintf(stdout, ": ");
        vfprintf(stdout, msg, ap);
        va_end(ap);
    }
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
}