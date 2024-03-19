#include "header.h"

int floatCompare(const void* a, const void* b) {
    float diff = *(const float*)a - *(const float*)b;
    return (diff < 0) ? -1 : (diff > 0);
}
