#include "header.h"

int floatCompare(const void *a, const void *b)
{
    float diff = *(const float *)a - *(const float *)b;
    if (diff < 0)
    {
        return -1;
    }
    else if (diff > 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
