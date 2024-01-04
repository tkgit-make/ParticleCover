#include "debug.h"
#include <stdio.h>

void print_DataPointArr(DataPointArr_s *data_points_arr) {
  for (int i = 0; i < data_points_arr->num_points; i++) {
    printf("layer: %d, radius: %d, angle: %lf, z: %lf\n",
           data_points_arr->points[i].layer, data_points_arr->points[i].radius,
           data_points_arr->points[i].angle, data_points_arr->points[i].z);
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