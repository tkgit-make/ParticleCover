#include "debug.h"
#include "parser.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  // must be single argument
  if (argc != 2) {
    // usage: ./<argv[0]> <file_path>
    printf("usage: %s <file_path>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  print_constants();

  // read the points from the file
  char *file_path = argv[1];
  PointArr_s data_points_arr = init_PointArr();
  read_points(file_path, &data_points_arr);

  // print the points
  print_DataPointArr(&data_points_arr);

  return 0;
}