#include "header.h"

//c doesn't support overloading
void initDataSetBase(DataSet* ds) {
    ds->total_points = 0;
    memset(ds->n_points, 0, sizeof(ds->n_points)); //initialize n_points with 0
}

void initDataSetExtra(DataSet* ds, Environment* envI) {
    ds->total_points = 0;
    //lhs is pointer, consistent with rhs
    ds->env = envI;
    memset(ds->n_points, 0, sizeof(ds->n_points)); 

    for (int i = 0; i < MAX_LAYERS; i++) {
        for (int j = 0; j < MAX_POINTS_FOR_DATASET; j++) {
            ds->array[i][j] = (Point){0, 0.0f, 0.0f, 0.0f}; //initialize with default values
        }
    }

    //redundant
    //for (int i = 0; i < MAX_LAYERS; i++) {
    //    ds->n_points[i] = 0;
    //}


    //continue writing importData and addBoundaryPoint
}
