#include "parser.h"
#include "constants.h"
#include "types.h"

#include <stdio.h>
#include <stdlib.h>

DataPointArr_s init_DataPointArr() {
  DataPointArr_s data_point_arr;
  data_point_arr.num_points = 0;

  return data_point_arr;
}

void read_points(char *file_path, DataPointArr_s *data_points_arr_out) {
  FILE *file = fopen(file_path, "r");

  if (file == NULL) {
    printf("[ERR] Error opening file!\n");
    exit(EXIT_FAILURE);
  }

  /**
   * read only the first line.
   * line format: (layer,radius,angle,z),(layer,radius,angle,z),...
   * example: (0,0,0,0),(0,0,0,0),(0,0,0,0)
   */
  // loop through each character in the file,
  // until reaching the end of the line
  // data format: (layer,radius,angle,z), separated by commas
  DataPoint_s temp_point;

  while (fscanf(file, "(%d,%d,%lf,%lf)", &temp_point.layer, &temp_point.radius,
                &temp_point.angle, &temp_point.z) == 4) {
    // increment the number of points
    data_points_arr_out->num_points++;

    // add the new point to the array
    data_points_arr_out->points[data_points_arr_out->num_points - 1] =
        temp_point;

    // skip the comma if next character is a comma
    if (fgetc(file) == ',') {
      continue;
    }

    // if reach next line, break
    if (fgetc(file) == '\n') {
      break;
    }
  }

  fclose(file);

  return;
}