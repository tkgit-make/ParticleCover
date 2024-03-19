
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) < (Y) ? (Y) : (X))

#define MAX_LAYERS 5
#define MAX_POINTS_IN_EVENT 512
#define MAX_POINTS_FOR_DATASET 512 //max size of vector of points "vect" in CPP
#define MAX_POINTS_PER_LAYER 256 //not yet used

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
    int n_points[MAX_LAYERS]; 
    int total_points;
    float boundaryPoint_offset;
} DataSet;

extern int Point_load(Point* p);
extern int Event_load(Event* e);
extern void initEnvironment(Environment* env, float top_layer_limI, float beam_axis_limI, int num_layersI, float* radiiI);
extern void initDataSetBase(DataSet* ds);
extern void initDataSetExtra(DataSet* ds, Environment* envI);

extern int floatCompare(const void* a, const void* b);

EXTERN Event G_event;