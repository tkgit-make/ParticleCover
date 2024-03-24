
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
} parallelogram_v1;

typedef struct {
    Point points[MAX_POINTS_IN_WEDGESUPERPOINT];
    float z_values[MAX_POINTS_IN_WEDGESUPERPOINT];
    int point_count;
    float min;
    float max;
} wedgeSuperPoint;

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

extern int floatCompare(const void* a, const void* b);

EXTERN Event G_event;