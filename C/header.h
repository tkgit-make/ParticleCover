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

#define MAX_EVENTS_TO_READ 16000
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

//constant from environment that have been pulled out of structure
#define num_layers 5
#define top_layer_lim 50
#define beam_axis_lim 15
#ifdef MAIN_C
const float radii[MAX_LAYERS] = {5, 10, 15, 20, 25};
const float parallelogramSlopes[MAX_LAYERS-1] = {0, -0.333333, -1, -3};
#else
extern float radii[MAX_LAYERS];
extern float parallelogramSlopes[MAX_LAYERS-1];
#endif


typedef struct
{
    index_type layer_num;
    float radius;
    float phi;
    float z;
} Point;

typedef struct
{
    float radii_leverArm[MAX_LAYERS-1];
    float trapezoid_edges[MAX_LAYERS];
} Environment;

typedef struct
{
    Point points[MAX_POINTS_IN_EVENT];
    index_type count;
    Environment env;
} Event;

typedef struct
{
    Environment *env;
    Point array[MAX_LAYERS][MAX_POINTS_FOR_DATASET]; // 2D array of points
    int n_points[MAX_LAYERS];                        // number of points in each layer of the array
    //index_type total_points; //not used
    float boundaryPoint_offset;
} DataSet;

typedef struct
{
    Environment *env;
    float slope;
    float points[MAX_POINTS_IN_LINE];
    index_type num_points; // number of points in the line
} Line;

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
    Environment* env;
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

    // Parallelogram_v1* parallelograms_v1[MAX_PARALLELOGRAMS_PER_PATCH];
    // int parallelogram_v1_count;
} wedgePatch;

typedef struct
{
    index_type n_patches;
    wedgePatch patches[MAX_PATCHES];
    Environment *env;
    DataSet *data;
    // Line fitting_lines[MAX_LINES]; //only used in visualization
    //wedgeSuperPoint *superPoints[MAX_SUPERPOINTS_IN_COVER]; //not used
    //wedgePatch *all_patches[MAX_PATCHES]; //not needed anymore
    bool real_patch_list[MAX_PATCHES];
} wedgeCover;

extern int Point_load(Point *p, float *radius);
// extern void Point_init(Point* p, int layerNum, float rad, float ph, float zVal);
// extern void Event_init(Event* e, Environment* envI, Point* pointsArray, int numPoints);
extern index_type Event_load(Event *e);
extern void initEnvironment(Environment *env, float top_layer_limI, float beam_axis_limI, int num_layersI, float *radiiI);
extern void initDataSetBase(DataSet *ds);
extern void initDataSetExtra(DataSet *ds, Environment *envI);
extern void importData(DataSet *ds, Point *data_array, int data_array_size);
extern void addBoundaryPoint(DataSet *ds, float offset);
extern void initLine(Line *line, Environment *envI, float start, float slopeI);
extern void initWedgeSuperPoint(wedgeSuperPoint *wsp, Point *points, int pointCount);
extern int areWedgeSuperPointsEqual(wedgeSuperPoint *wsp1, wedgeSuperPoint *wsp2);
extern void initParallelogram(Parallelogram *pg, int layer_numI, float z1_minI, float z1_maxI, float shadow_bottomL_jRI, float shadow_bottomR_jRI, float shadow_bottomL_jLI, float shadow_bottomR_jLI, float pSlopeI);
extern void wedgePatch_init(wedgePatch *wp, Environment *envI, wedgeSuperPoint *superpointsI, int superpoint_count, float apexZ0I);
extern float straightLineProjectorFromLayerIJtoK(wedgePatch *wp, float z_i, float z_j, int i, int j, int k);
extern float straightLineProjector(float z_top, float z_j, int j, Environment *env);
extern void getParallelograms(wedgePatch *wp);
extern void getShadows(wedgePatch *wp, float zTopMin, float zTopMax);
extern void get_acceptanceCorners(wedgePatch *wp);
extern void get_end_layer(wedgePatch *wp);
extern void initWedgeCover(wedgeCover *wc, Environment *envI, DataSet *dataI);
extern int comparePoints(const void *a, const void *b);
extern void add_patch(wedgeCover *cover, wedgePatch *curr_patch);
extern void delete_patch(wedgeCover *cover, int index);
extern index_type get_index_from_z(DataSet *data, int layer, float z_value, int alignment);
extern void solve(wedgeCover *cover, int lining, float apexZ0, int ppl, int nlines, bool leftRight);
extern void makePatches_ShadowQuilt_fromEdges(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight);
extern void makePatch_alignedToLine(wedgeCover *cover, float apexZ0, float z_top, int ppl, bool leftRight, bool float_middleLayers_ppl);
extern void wedge_test(float apexZ0, float z0_spacing, int ppl, float z0_luminousRegion, int wedges[], int wedge_count, int lines, float top_layer_cutoff, float accept_cutoff, bool leftRightAlign, bool uniform_N_points, bool show_acceptance_of_cover, bool movie, bool savefig, int figSizeScale, int movieFigSizeScale);

extern int floatCompare(const void *a, const void *b);

//cm to microns x10000 for points.