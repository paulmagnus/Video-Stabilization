#ifdef OUTPUT
#define OUTPUT
#include <stdarg.h>
#include <stdio.h>

// comment to remove print statements
#define PRINT_PROGRESS_REPORT

void print(const char* format, ...) {
#ifdef PRINT_PROGRESS_REPORT
	va_list argptr;
	va_start(argptr, format);
	vfprintf(stderr, format, argptr);
	va_end(argptr);
#endif
}
#endif