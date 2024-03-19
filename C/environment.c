#include "header.h"

void initEnvironment(Environment* env, float top_layer_limI, float beam_axis_limI, int num_layersI, float* radiiI) {
    if (top_layer_limI < beam_axis_limI) {
        printf("The top layer limits cannot be smaller than the bottom layer limits.");
        exit(0);
    }
    env->top_layer_lim = top_layer_limI;
    env->beam_axis_lim = beam_axis_limI;

    /* //not checking size, assume parameters are compatible with each other
    if (radiiI.size() != num_layersI) {
        printf("The radii do not match the number of layers.");
        exit(1);
    }
    */
   if (MAX_LAYERS > num_layersI) {
    printf("Arrays aren't big enough");
    exit(2);
   }

   env->num_layers = num_layersI;

    //copy radii values
    for (int i = 0; i < env->num_layers; i++) {
        env->radii[i] = radiiI[i];
    }

    //qsort rearranges the elements within the block of memory
    qsort(env->radii, env->num_layers, sizeof(float), floatCompare);

    //populating parallelogramSlopes vector
    //bounds are correct
    float radiiLast = env->radii[num_layersI - 1]; 
    for (int i = 0; i < env->num_layers - 1; i++) {
        float radiiFirst = env->radii[0];
        float radiiCurrent = env->radii[i];
        env->parallelogramSlopes[i] = (radiiFirst - radiiCurrent) / (radiiLast - radiiCurrent);
    }

    //calculate radii lever arm
    //bounds are correct
    for (int i = 0; i < env->num_layers - 1; i++) {
        env->radii_leverArm[i] = 1 - env->parallelogramSlopes[i];
    }

    env->boundaryPoint_offset = 0;

    //calculate trapezoid edges
    //bounds are correct
    for (int i = 0; i < env->num_layers; i++) {
        float radiiCurrent = env->radii[i];
        env->trapezoid_edges[i] = radiiCurrent * (env->top_layer_lim - env->beam_axis_lim) / radiiLast + env->beam_axis_lim;
    }


}