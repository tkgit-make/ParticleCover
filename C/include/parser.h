#ifndef __PARSER_H__
#define __PARSER_H__

#include "types.h"

PointArr_s init_PointArr();

void read_points(char *file_path, PointArr_s *points_arr_out);

#endif