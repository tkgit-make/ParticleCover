#ifndef __TYPES_H__
#define __TYPES_H__

#include "constants.h"

typedef struct {
  int layer;
  int radius;
  double angle;
  double z;
} Point_s;

typedef struct {
  unsigned int num_points;
  Point_s points[MAX_NUM_POINTS];
} PointArr_s;

#endif