#include <stdio.h>
#include <stdlib.h>
#include "parser.h"

int main()
{
    char *file_path = "/Users/crinstaniev/Research/aifpga/ParticleCover/C/data/wedgeData_v3_128.txt";

    DataPointArr_s *data_points_arr = malloc(sizeof(DataPointArr_s));
    data_points_arr->num_points = 0;
    data_points_arr->points = NULL;

    read_points(file_path, data_points_arr);

    // DEBUG: print the points
    for (int i = 0; i < data_points_arr->num_points; i++)
    {
        printf(
            "layer: %d, radius: %d, angle: %lf, z: %lf\n",
            data_points_arr->points[i].layer,
            data_points_arr->points[i].radius,
            data_points_arr->points[i].angle,
            data_points_arr->points[i].z);
    }

    return 0;
}