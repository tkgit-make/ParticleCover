#include <stdio.h>
#include <stdlib.h>

struct Point { // represents a point
    int layer_num; 
    int radius;
    float phi;
    float z;
};

struct Environment { // environment of the collider
    float top_layer_lim;
    float beam_axis_lim;
    int num_layers;
    float* radii;
};

// create new Environment
struct Environment *env_new(float top_layer_lim, float beam_axis_lim, int num_layers, int radii_size) {
    struct Environment *env = (Environment *)malloc(sizeof(Environment));
    env->top_layer_lim = top_layer_lim;
    env->beam_axis_lim = beam_axis_lim;
    env->num_layers = num_layers;

    float* r = (float *)malloc(radii_size * sizeof(float));
    env->radii = r;

    return env;
}

struct Patch {
    struct Environment *env;
    float** wedgeSuperPoint;
};

// create new Patch
struct Patch *patch_new(struct Environment *env, int pts_per_layer) {
    struct Patch *patch = (Patch *)malloc(sizeof(Patch));
    patch->env = env;

    float** wedgeSuperPoint = (float **)malloc(env->num_layers * sizeof(float *));
    for (int i = 0; i < (env->num_layers); i++)
        wedgeSuperPoint[i] = (float*)malloc(pts_per_layer * sizeof(float));
    patch->wedgeSuperPoint = wedgeSuperPoint;

    return patch;
}