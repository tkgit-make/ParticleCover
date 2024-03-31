
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) < (Y) ? (Y) : (X))

#define MAX_LAYERS 5
#define MAX_POINTS_IN_EVENT 512
#define MAX_POINTS_FOR_DATASET 512 //max size of vector of points "vect" in CPP
#define MAX_POINTS_PER_LAYER 256 //not yet used
#define MAX_POINTS_IN_LINE MAX_LAYERS //a point on the line is calculated for each layer in the environment.
#define MAX_POINTS_IN_WEDGESUPERPOINT 32
#define MAX_SUPERPOINTS_IN_PATCH 5
#define MAX_PARALLELOGRAMS_PER_PATCH MAX_SUPERPOINTS_IN_PATCH //not sure. could it be 1 parallelogram per superpoint?

#ifdef MAIN_C
	#define EXTERN
#else
	#define EXTERN extern
#endif

typedef struct
{
	int layer_num;
	float radius;
	float phi;
	float z;
} Point;

typedef struct
{
	Point points[MAX_POINTS_IN_EVENT];
	int count;
} Event;

typedef struct {
    float top_layer_lim;
    float beam_axis_lim;
    int num_layers;
    float radii[MAX_LAYERS];
    float parallelogramSlopes[MAX_LAYERS];
    float radii_leverArm[MAX_LAYERS];
    float trapezoid_edges[MAX_LAYERS];
    float boundaryPoint_offset; 
} Environment;

typedef struct {
    Environment* env;
    Point array[MAX_LAYERS][MAX_POINTS_FOR_DATASET]; //2D array of points
    int n_points[MAX_LAYERS]; //number of points in each layer of the array
    int total_points;
    float boundaryPoint_offset;
} DataSet;


typedef struct {
    Environment* env;
    float slope;
    float points[MAX_POINTS_IN_LINE]; 
    int num_points; //number of points in the line
} Line;

typedef struct {
    Environment* env;
    float start;
    float slope_ll;
    float slope_ul;
} LineGenerator;

typedef struct {
    int layer_num;
    float pSlope;

    float shadow_bottomL_jR;
    float shadow_bottomR_jR;
    float shadow_bottomL_jL;
    float shadow_bottomR_jL;

    float z1_min;
    float z1_max;
} Parallelogram;

typedef struct {
    int layer_num;
    float pSlope;

    float shadow_topR_jL;
    float shadow_topR_jR;
    float shadow_topL_jL;
    float shadow_topL_jR;

    float top_layer_zmin;
    float top_layer_zmax;
} Parallelogram_v1;

typedef struct {
    Point points[MAX_POINTS_IN_WEDGESUPERPOINT];
    float z_values[MAX_POINTS_IN_WEDGESUPERPOINT];
    int point_count;
    float min;
    float max;
} wedgeSuperPoint;

typedef struct {
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

    wedgeSuperPoint* superpoints[MAX_SUPERPOINTS_IN_PATCH]; //array of pointers
    int superpoint_count;

    bool flatBottom;
    bool flatTop;

    bool squareAcceptance;
    bool triangleAcceptance;

    Parallelogram* parallelograms[MAX_PARALLELOGRAMS_PER_PATCH];
    int parallelogram_count;

    Parallelogram_v1* parallelograms_v1[MAX_PARALLELOGRAMS_PER_PATCH];
    int parallelogram_v1_count;
} wedgePatch;


extern int Point_load(Point* p);
extern int Event_load(Event* e);
extern void initEnvironment(Environment* env, float top_layer_limI, float beam_axis_limI, int num_layersI, float* radiiI);
extern void initDataSetBase(DataSet* ds);
extern void initDataSetExtra(DataSet* ds, Environment* envI);
extern void importData(DataSet* ds, Point* data_array, int data_array_size);
extern void addBoundaryPoint(DataSet* ds, float offset);
extern void initLine(Line* line, Environment* envI, float start, float slopeI);
extern void initLineGenerator(LineGenerator* lg, Environment* envI, float startI);
extern void generateEvenGrid(LineGenerator* lg, Line* lines, int n);
extern void initWedgeSuperPoint(wedgeSuperPoint* wsp, Point* points, int pointCount);
extern int areWedgeSuperPointsEqual(wedgeSuperPoint* wsp1, wedgeSuperPoint* wsp2);
extern void initParallelogram(Parallelogram* pg, int layer_numI, float z1_minI, float z1_maxI, float shadow_bottomL_jRI, float shadow_bottomR_jRI, float shadow_bottomL_jLI, float shadow_bottomR_jLI, float pSlopeI);
extern void init_parallelogram_v1(Parallelogram_v1 *pg, int layer_numI, float top_layer_zminI, float top_layer_zmaxI, float shadow_topR_jLI, float shadow_topR_jRI, float pSlopeI);
extern void wedgePatch_init(wedgePatch* wp, Environment* envI, wedgeSuperPoint* superpointsI, int superpoint_count, float apexZ0I);
extern float straightLineProjectorFromLayerIJtoK(wedgePatch* wp, float z_i, float z_j, int i, int j, int k);
extern float straightLineProjector(float z_top, float z_j, int j, Environment* env);
extern void getParallelograms(wedgePatch* wp);
extern void getParallelograms_v1(wedgePatch* wp);
extern void getShadows(wedgePatch* wp, float zTopMin, float zTopMax);
extern void get_acceptanceCorners(wedgePatch* wp);
extern void get_end_layer(wedgePatch* wp);

extern int floatCompare(const void* a, const void* b);

EXTERN Event G_event;