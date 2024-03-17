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

void initEnvironment(Environment* env, float top_layer_limI, float beam_axis_limI, int num_layersI, floatVector* radiiI) {
    if (top_layer_limI < beam_axis_limI) {
        printf("The top layer limits cannot be smaller than the bottom layer limits.");
        exit(0);
    }
    env->top_layer_lim = top_layer_limI;
    env->beam_axis_lim = beam_axis_limI;

    if (((Vector*)radiiI)->usedCount != num_layersI) {
        printf("The radii do not match the number of layers.");
        exit(1);
    }

    env->num_layers = num_layersI;

    //env->radii = VectorOf_float(64, 64); //placeholder values
    env->parallelogramSlopes = VectorOf_float(64, 64);
    env->radii_leverArm = VectorOf_float(64, 64);
    env->trapezoid_edges = VectorOf_float(64, 64);

    //pointing env->radii to the same floatVector radiiI points to.
    env->radii = radiiI;
    //copy radii values into the vector. equivalent to radii = radiiI in CPP
    /*
    for (int i = 0; i < num_layersI; i++) {
        float* newRadiiItem = floatVector_newitem(env->radii);
        *newRadiiItem = radiiI[i];
    }
    */
    //qsort rearranges the elements within the block of memory
    qsort(((Vector*)env->radii)->data, ((Vector*)env->radii)->usedCount, sizeof(float), floatCompare);

    //populating parallelogramSlopes vector
    //env->num_layers is the same as ((Vector*)radiiI)->usedCount
    float radiiLast = *(floatVector_getitem(env->radii, env->num_layers - 1));
    for (int i = 0; i < env->num_layers - 1; i++) {
        float radiiFirst = *(floatVector_getitem(env->radii, 0));
        float radiiCurrent = *(floatVector_getitem(env->radii, i));
        float currentVal = (radiiFirst - radiiCurrent) / (radiiLast - radiiCurrent);
        //method returns pointer, must dereference (accessing value) to assign value
        //Without the *, newSlopeItem is just a memory address.
        float* newSlopeItem = floatVector_newitem(env->parallelogramSlopes);
        *newSlopeItem = currentVal;
    }

    //calculate radii lever arm
    for (int i = 0; i < ((Vector*)env->parallelogramSlopes)->usedCount; i++) {
        float slope = *(floatVector_getitem(env->parallelogramSlopes, i));
        float* newLeverArmItem = floatVector_newitem(env->radii_leverArm);
        *newLeverArmItem = 1 - slope;
    }

    env->boundaryPoint_offset = 0;

    //calculate trapezoid edges
    for (int i = 0; i < env->num_layers; i++) {
        float radiiCurrent = *(floatVector_getitem(env->radii, i));
        float currentVal = radiiCurrent * (env->top_layer_lim - env->beam_axis_lim) / radiiLast + env->beam_axis_lim;

        float* newTrapezoidEdge = floatVector_newitem(env->trapezoid_edges);
        *newTrapezoidEdge = currentVal;
    }


}