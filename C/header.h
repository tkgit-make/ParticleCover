#ifndef HEADER_H
#define HEADER_H

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) < (Y) ? (Y) : (X))

static inline int floatCompare(const void* a, const void* b) {
    float diff = *(const float*)a - *(const float*)b;
    return (diff < 0) ? -1 : (diff > 0); //-1 if less than, 1 if greater than, 0 if equal
}

#endif