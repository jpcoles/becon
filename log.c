#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include "log.h"

void Log(const char *fmt, ...)
{
    char s[100];
    va_list ap;
    time_t t;
    struct tm *tmp;
    t = time(NULL);
    tmp = localtime(&t);

    strftime(s, sizeof(s), "%F %T", tmp);

    va_start(ap, fmt);
    fprintf(stderr, "%15s  |  ", s);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
}
