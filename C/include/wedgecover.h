#ifndef WEDGECOVER_H
#define WEDGECOVER_H

const int numLayers = 5;
const int numPtsPerSuperpoint = 16;
const int maxPatches = 100;

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
    float radii[numLayers];
};

struct Patch { // patch
    struct Environment env;
    float Superpoints[numLayers][numPtsPerSuperpoint];
};

struct Cover { // cover, consisting of patches
    struct Patch patches[maxPatches];
};

#endif