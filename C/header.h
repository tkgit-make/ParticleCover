// gdb bin/ProcessInput
// run < CPP/wedgeData_v3_128.txt
// bt
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>


//make conversion ratio to micron macro

#ifdef MAIN_C
#define EXTERN
#else
#define EXTERN extern
#endif

#define min(X, Y) ((X) < (Y) ? (X) : (Y))
#define max(X, Y) ((X) < (Y) ? (Y) : (X))

#define index_type int //change to unsigned int once it is verified that there are no errors

#define CLOSEST 11
#define ABOVE 12
#define BELOW 13
#define MAKE_PATCHES_SHADOW_QUILT_FROM_EDGES 33

#define MAX_LAYERS 5
#define MAX_POINTS_IN_EVENT 512
#define MAX_POINTS_PER_LAYER 256    // max size of vector of points "vect" in CPP. equivalent to MAX_POINTS_PER_DATASET
#define MAX_POINTS_FOR_DATASET MAX_POINTS_PER_LAYER    // max size of vector of points "vect" in CPP

#define MAX_POINTS_IN_LINE MAX_LAYERS // a point on the line is calculated for each layer in the environment.
#define MAX_POINTS_IN_SUPERPOINT 32
#define MAX_SUPERPOINTS_IN_PATCH 5
#define MAX_PARALLELOGRAMS_PER_PATCH MAX_LAYERS - 1 // layer 1 is a vertical ribbon, the other 4 layers are sloping, so each intersects with layer 1 to make a parallelogram
#define MAX_PATCHES 32                              // upper bound, 14-18 average.
// #define MAX_LINES __ //only used in visualization
#define MAX_SUPERPOINTS_IN_COVER (MAX_PATCHES * MAX_SUPERPOINTS_IN_PATCH)

#ifdef MAIN_C
    const index_type CONVERSION_FACTOR = 1;
    const float top_layer_lim = 50 * CONVERSION_FACTOR;
    const float beam_axis_lim = 15 * CONVERSION_FACTOR;
    const index_type num_layers = 5;
    const index_type radii[MAX_LAYERS] = {
        5 * CONVERSION_FACTOR, 
        10 * CONVERSION_FACTOR, 
        15 * CONVERSION_FACTOR, 
        20 * CONVERSION_FACTOR, 
        25 * CONVERSION_FACTOR
    };
    const float parallelogramSlopes[MAX_LAYERS-1] = {
        0,
        (radii[0] - radii[1]) / (radii[4] - radii[1]), 
        (radii[0] - radii[2]) / (radii[4] - radii[2]), 
        (radii[0] - radii[3]) / (radii[4] - radii[3])
    };
    const float radii_leverArm[MAX_LAYERS-1] = {
        1 - parallelogramSlopes[0], 
        1 - parallelogramSlopes[1], 
        1 - parallelogramSlopes[2], 
        1 - parallelogramSlopes[3]
    };
    const float trapezoid_edges[MAX_LAYERS] = {
        radii[0] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        radii[1] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        radii[2] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        radii[3] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR, 
        radii[4] * (top_layer_lim - beam_axis_lim) / radii[4] + beam_axis_lim + 0.0001 * CONVERSION_FACTOR
    };
#else
    extern index_type radii[MAX_LAYERS];
    extern float parallelogramSlopes[MAX_LAYERS-1];
    extern float radii_leverArm[MAX_LAYERS-1];
    extern float trapezoid_edges[MAX_LAYERS];
    extern index_type CONVERSION_FACTOR;
    extern float top_layer_lim;
    extern float beam_axis_lim;
    extern index_type num_layers;
#endif


typedef struct
{
    index_type layer_num;
    index_type radius;
    float phi;
    float z;
} Point;

typedef struct
{
    Point array[MAX_LAYERS][MAX_POINTS_FOR_DATASET]; // 2D array of points
    int n_points[MAX_LAYERS];                        // number of points in each layer of the array
    //index_type total_points; //not used
    float boundaryPoint_offset;
} DataSet;

typedef struct
{
    index_type layer_num;
    float pSlope;

    float shadow_bottomL_jR;
    float shadow_bottomR_jR;
    float shadow_bottomL_jL;
    float shadow_bottomR_jL;

    float z1_min;
    float z1_max;
} Parallelogram;

typedef struct
{
    Point points[MAX_POINTS_IN_SUPERPOINT];
    float z_values[MAX_POINTS_IN_SUPERPOINT];
    index_type point_count;
    float min;
    float max;
} wedgeSuperPoint;

typedef struct
{
    int end_layer;
    int left_end_layer;
    int right_end_layer;
    float left_end_lambdaZ;
    float right_end_lambdaZ;
    float apexZ0;

    float shadow_fromTopToInnermost_topL_jL;
    float shadow_fromTopToInnermost_topL_jR;
    float shadow_fromTopToInnermost_topR_jL;
    float shadow_fromTopToInnermost_topR_jR;

    float a_corner[2];
    float b_corner[2];
    float c_corner[2];
    float d_corner[2];

    wedgeSuperPoint superpoints[MAX_SUPERPOINTS_IN_PATCH]; //changed to direct assignment as opposed to pointer
    index_type superpoint_count;

    bool flatBottom;
    bool flatTop;

    bool squareAcceptance;
    bool triangleAcceptance;

    Parallelogram parallelograms[MAX_PARALLELOGRAMS_PER_PATCH];
    index_type parallelogram_count;
} wedgePatch;

extern int Point_load(Point *p);
extern void importData();
extern void addBoundaryPoint(float offset);
extern void initWedgeSuperPoint(wedgeSuperPoint *wsp, Point *points, int pointCount);
extern int areWedgeSuperPointsEqual(wedgeSuperPoint *wsp1, wedgeSuperPoint *wsp2);
extern void initParallelogram(Parallelogram *pg, int layer_numI, float z1_minI, float z1_maxI, float shadow_bottomL_jRI, float shadow_bottomR_jRI, float shadow_bottomL_jLI, float shadow_bottomR_jLI, float pSlopeI);
extern void wedgePatch_init(wedgePatch *wp, wedgeSuperPoint *superpointsI, int superpoint_count, float apexZ0I);
extern float straightLineProjectorFromLayerIJtoK(float z_i, float z_j, int i, int j, int k);
extern float straightLineProjector(float z_top, float z_j, int j);
extern void getParallelograms(wedgePatch *wp);
extern void getShadows(wedgePatch *wp, float zTopMin, float zTopMax);
extern void get_acceptanceCorners(wedgePatch *wp);
extern void get_end_layer(wedgePatch *wp);
extern void initWedgeCover();
extern int comparePoints(const void *a, const void *b);
extern void add_patch(wedgePatch *curr_patch);
extern void delete_patch(int index);
extern index_type get_index_from_z(int layer, float z_value);
extern void solve(float apexZ0, int ppl, int nlines, bool leftRight);
extern void makePatches_ShadowQuilt_fromEdges(float apexZ0, int stop, int ppl, bool leftRight);
extern void makePatch_alignedToLine(float apexZ0, float z_top, int ppl, bool leftRight, bool float_middleLayers_ppl);
extern void wedge_test(float apexZ0, int ppl, int wedges[]);
extern int floatCompare(const void *a, const void *b);

EXTERN DataSet Gdata;
EXTERN wedgePatch patches[MAX_PATCHES];
EXTERN index_type n_patches;
