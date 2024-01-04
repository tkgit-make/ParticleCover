#ifndef __PARSER_H__
#define __PARSER_H__

#include "types.h"

DataPointArr_s init_DataPointArr();

void read_points(char *file_path, DataPointArr_s *data_points_arr_out);

#endif