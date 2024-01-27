#include "debug.h"
#include <stdio.h>

void print_DataPointArr(PointArr_s *point_arr) {
  for (int i = 0; i < point_arr->num_points; i++) {
    printf("layer: %d, radius: %d, angle: %lf, z: %lf\n",
           point_arr->points[i].layer, point_arr->points[i].radius,
           point_arr->points[i].angle, point_arr->points[i].z);
  }
}

void print_constants() {
  printf("CONSTANTS:\n");

  // max number of points
  printf("MAX_NUM_POINTS: %d\n", MAX_NUM_POINTS);

  // max line length
  printf("MAX_LINE_LENGTH: %d\n", MAX_LINE_LENGTH);

  return;
}