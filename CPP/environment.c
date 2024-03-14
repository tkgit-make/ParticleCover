#include <stdlib.h>
#include "header.h"
#ifndef VECTOR_C
	#include "vector.c"
#endif

CREATE_VECTOR_OF_T(float)

typedef struct Environment {
    float top_layer_lim;
    float beam_axis_lim;
    int num_layers;
    floatVector* radii;
    floatVector* parallelogramSlopes;
    floatVector* radii_leverArm;
    floatVector* trapezoid_edges;
    float boundaryPoint_offset;
} Environment;

void initEnvironment(Environment* env, float top_layer_limI, float beam_axis_limI, int num_layersI, float* radiiArray) {
    if (top_layer_limI < beam_axis_limI) {
        printf("The top layer limits cannot be smaller than the bottom layer limits.");
        exit(0);
    }
    env->top_layer_lim = top_layer_limI;
    env->beam_axis_lim = beam_axis_limI;
    env->num_layers = num_layersI;

    env->radii = VectorOf_float(10, 5); //placeholder values
    env->parallelogramSlopes = VectorOf_float(64, 64);
    env->radii_leverArm = VectorOf_float(64, 64);
    env->trapezoid_edges = VectorOf_float(64, 84);

    //copy radii values into the vector. equivalent to radii = radiiI in CPP
    for (int i = 0; i < num_layersI; i++) {
        float* newRadiiItem = floatVector_newitem(env->radii);
        *newRadiiItem = radiiArray[i];
    }

    qsort(env->radii->data, env->radii->usedCount, sizeof(float), floatCompare);

    //calculate parallelogram slopes
    for (int i = 0; i < num_layersI - 1; i++) {
        float* newSlopeItem = floatVector_newitem(env->parallelogramSlopes);
        float radiiLast = *(floatVector_getitem(env->radii, num_layersI - 1));
        float radiiFirst = *(floatVector_getitem(env->radii, 0));
        float radiiCurrent = *(floatVector_getitem(env->radii, i));
        *newSlopeItem = (radiiFirst - radiiCurrent) / (radiiLast - radiiCurrent);
    }
}