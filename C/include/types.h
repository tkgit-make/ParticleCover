#ifndef __TYPES_H__
#define __TYPES_H__

typedef struct
{
    int layer;
    int radius;
    double angle;
    double z;
} DataPoint_s;

typedef struct
{
    unsigned int num_points;
    DataPoint_s *points;
} DataPointArr_s;

#endif