#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

typedef struct Point {
    int layer_num;
    float radius;
    float phi;
    float z;
} Point;

void initPoint(Point *point, int layerNum, float rad, float ph, float zVal) {
    point->layer_num = layerNum;
    point->radius = rad;
    point->phi = ph;
    point->z = zVal;
}

typedef struct Environment {
    float top_layer_lim;
    float beam_axis_lim;
    int num_layers;
    //vectors
    float* radii;
    float* parallelogramSlopes;
    float* radii_leverArm;
    float* trapezoid_edges;
    //

    float boundaryPoint_offset;
} Environment;

int cmpfunc (const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

void initEnvironment(Environment *env, float top_layer_limI, float beam_axis_limI, int num_layersI, float *radiiI) {
    if (top_layer_limI < beam_axis_limI) {
        printf("The top layer limits cannot be smaller than the bottom layer limits.");
        exit(0);
    }
    env->top_layer_lim = top_layer_limI;
    env->beam_axis_lim = beam_axis_limI;

    /**
    if (radiiI.size() != num_layersI) {
        printf("The radii do not match the number of layers.");
        exit(1);
    }
    */

    env->num_layers = num_layersI;
    env->radii = (float*)calloc(num_layersI, sizeof(float));

    qsort(env->radii, env->num_layers, sizeof(float), cmpfunc);

    for (int i = 0; i < num_layersI -1; i++) {
        //env->radii[i] = radiiI[i];
        float currentVal = (env->radii[0] - env->radii[i]) / (radii[radii.size() - 1] - radii[i]);
        //float currentVal = (radii[0] - radii[i]) / (radii[radii.size() - 1] - radii[i]);
        //parallelogramSlopes.push_back(currentVal);
    }

}
