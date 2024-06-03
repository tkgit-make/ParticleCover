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


#define KEEP_DELETED_PATCHES true
//make conversion ratio to micron macro

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

//constant from environment that have been pulled out of structure
#define num_layers 5
#define top_layer_lim 50
#define beam_axis_lim 15


const float radii[MAX_LAYERS] = {5, 10, 15, 20, 25};
const float parallelogramSlopes[MAX_LAYERS-1] = {0, -0.333333, -1, -3};
const float radii_leverArm[MAX_LAYERS-1] = {1, 1.333333, 2, 4};
const float trapezoid_edges[MAX_LAYERS] = {22.0001, 29.0001, 36.0001, 43.0001, 50.0001};

typedef struct
{
    index_type layer_num;
    float radius;
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

typedef struct
{
    index_type n_patches; //make global
    wedgePatch patches[MAX_PATCHES];
    //DataSet *data; //make global
    #if KEEP_DELETED_PATCHES == true
        bool real_patch_list[MAX_PATCHES]; //only needed for video
    #endif

} wedgeCover;

int Point_load(Point *p);
void importData(DataSet *ds);
void addBoundaryPoint(DataSet *ds, float offset);
void initWedgeSuperPoint(wedgeSuperPoint *wsp, Point *points, int pointCount);
int areWedgeSuperPointsEqual(wedgeSuperPoint *wsp1, wedgeSuperPoint *wsp2);
void initParallelogram(Parallelogram *pg, int layer_numI, float z1_minI, float z1_maxI, float shadow_bottomL_jRI, float shadow_bottomR_jRI, float shadow_bottomL_jLI, float shadow_bottomR_jLI, float pSlopeI);
void wedgePatch_init(wedgePatch *wp, wedgeSuperPoint *superpointsI, int superpoint_count, float apexZ0I);
float straightLineProjectorFromLayerIJtoK(wedgePatch *wp, float z_i, float z_j, int i, int j, int k);
float straightLineProjector(float z_top, float z_j, int j);
void getParallelograms(wedgePatch *wp);
void getShadows(wedgePatch *wp, float zTopMin, float zTopMax);
void get_acceptanceCorners(wedgePatch *wp);
void get_end_layer(wedgePatch *wp);
void initWedgeCover(wedgeCover *wc);
int comparePoints(const void *a, const void *b);
void add_patch(wedgeCover *cover, wedgePatch *curr_patch);
void delete_patch(wedgeCover *cover, int index);
index_type get_index_from_z(int layer, float z_value);
void solve(wedgeCover *cover, float apexZ0, int ppl, int nlines, bool leftRight);
void makePatches_ShadowQuilt_fromEdges(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight);
float solveNextColumn(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight, bool fix42, float saved_apexZ0); 
void solveNextPatchPair(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight, bool fix42, float &saved_apexZ0, int &nPatchesInColumn, float &c_corner, float &projectionOfCornerToBeam, float &z_top_min, float &z_top_max, float &complementary_apexZ0);
void solveComplmentaryPatch(wedgeCover *cover, float &previous_white_space_height, int ppl, bool fix42, int nPatchesAtOriginal, float &previous_z_top_min, float complementary_apexZ0, float &white_space_height, index_type &lastPatchIndex, float original_c, float original_d, float &complementary_a, float &complementary_b, index_type &current_z_top_index, int &counter, int &counterUpshift, float &z_top_min, bool &repeat_patch, bool &repeat_original);
void makePatch_alignedToLine(wedgeCover *cover, float apexZ0, float z_top, int &ppl, bool leftRight, bool float_middleLayers_ppl);
void wedge_test(float apexZ0, float z0_spacing, int ppl, float z0_luminousRegion, int wedges[], int wedge_count, int lines, float top_layer_cutoff, float accept_cutoff);

int floatCompare(const void *a, const void *b);

DataSet Gdata;

int floatCompare(const void *a, const void *b)
{
    float diff = *(const float *)a - *(const float *)b;
    if (diff < 0)
    {
        return -1;
    }
    else if (diff > 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


/*
void Point_init(Point* p, int layerNum, float rad, float ph, float zVal) {
    p->layer_num = layerNum;
    p->radius = rad;
    p->phi = ph;
    p->z = zVal;
}
*/

int Point_load(Point *p)
{
    if (scanf("(%d,%f,%f,%f)", &p->layer_num, &p->radius, &p->phi, &p->z) == 4)
    {
        return 1;            // successful load
    }

    return 0; // failed to load
}

int comparePoints(const void *a, const void *b)
{

    const Point *pointA = (const Point *)a;
    const Point *pointB = (const Point *)b;
    
    //the below line is all that is needed. we perform additional checks to guarentee a unique ordering of points in debugging.
    //return (a_z < b_z) ? -1 : 1;

    if (pointA->z < pointB->z) return -1;
    if (pointA->z > pointB->z) return 1;

    if (pointA->layer_num < pointB->layer_num) return -1;
    if (pointA->layer_num > pointB->layer_num) return 1;

    if (pointA->phi < pointB->phi) return -1;
    if (pointA->phi > pointB->phi) return 1;

    return 0;
}


void adjustPointPositionFront(Point *array, int n_points, int start_index) {
    // move the point at start_index to its correct position to maintain sorted order
    Point toInsert = array[start_index];
    int j = start_index;
    //by checking if j < n_points-2, we are not going to the last index, thus, we will never have a situation where the end boundary point gets shifted before we position it
    //we check n_points-2 instead of n_points-1 because we have j+1 logic, j is the baseline and the comparison is with the next index, so we need j+1 to be not the end, but 1 away from it.
    //this is a valid cutoff because the z is the primary (first [and only in the case of the implementation, non-debugging comparator]) comparison in the comparator, and the trapezoid edges are always positive integers, so -x < x when x is a positive integer. 
    //it cannot be 0, so there is no possible equality as well, which could affect the debugging comparator.
    while (j < n_points - 2 && comparePoints(&array[j + 1], &toInsert) < 0) { // once we find one element does not need to be moved, we can stop, because the array is monotonic because it is sorted
        array[j] = array[j + 1]; // shift elements left, the other element(s) should come before the boundary point
        j++;
    }
    array[j] = toInsert; // place the element at its correct position
}

void adjustPointPositionBack(Point *array, int n_points, int start_index) {
    // move the point at start_index to its correct position to maintain sorted order
    Point toInsert = array[start_index];
    int j = start_index;
    //similarly, j > 1 ensures it doesn't reach the first index [j will end at 1 after checking if 2 should swap with 1], which while it wouldn't throw off the front position such that the adjustFront method doesn't work because that has already been called,
    //it is beneficial not to check the first index because it is a pointless computation. we can guarentee it will not shift.
    while (j > 1 && comparePoints(&array[j - 1], &toInsert) > 0) { // once we find one element does not need to be moved, we can stop, because the array is monotonic because it is sorted
        array[j] = array[j - 1]; // shift elements right, the other element(s) should come after the boundary point
        j--;
    } 
    array[j] = toInsert; // place the element at its correct position
}

void importData(DataSet *ds)
{
    // initDataSet line. The global DataSet is reused, so we just need to reset the number of points, or set it if this is the first time. 0 across all layers.
    // n_points being set to 0 when we reuse the DataSet will stop it from accessing any information from a past wedge.
    memset(ds->n_points, 0, sizeof(ds->n_points));

    index_type n = 0;
    char ch = ',';

    // read points until a non-comma is encountered or maximum points are read
    while ((ch == ',') && (n < MAX_POINTS_IN_EVENT))
    {
        Point p;
        if (Point_load(&p) < 1)
            break;
        
        index_type layer = p.layer_num - 1;
        ds->array[layer][ds->n_points[layer]+1] = p; //+1 leaves blank spot for the first boundary point
        ds->n_points[layer]++; //here n_points is not counting the blank spot at index 0. 
        
        n++;
        scanf("%c", &ch);

    }
       
    // iterating over the layers in DataSet
    for (index_type i = 0; i < num_layers; i++)
    {
        //sorts the points in the ith layer
        qsort(&ds->array[i][1], ds->n_points[i], sizeof(Point), comparePoints);
    }
}

void addBoundaryPoint(DataSet *ds, float offset)
{
    ds->boundaryPoint_offset = offset;

    for (index_type i = 0; i < num_layers; i++) {
        //adding two boundary points in each layer
        // inserting at the beginning
        ds->array[i][0].layer_num = i + 1;
        ds->array[i][0].radius = (i + 1) * 5;
        //is the phi for the boundary points used (answer: no)? so, instead of sorting in importData, we could wait and add boundary points, and then sort, without any shifting of boundary points needed. MlogM vs NlogN + 2N, where M = N+2
        ds->array[i][0].phi = ds->array[i][1].phi; // getting the first phi in the array sorted by z
        ds->array[i][0].z = -1 * ((trapezoid_edges[i]) - offset) - offset; //trapezoid edges is constant and initialized with the offset added. to preserve the original statement, we do it like this

        // appending at the end
        index_type lastIndex = ds->n_points[i] + 1; // after shifting, there's one more point
        ds->array[i][lastIndex].layer_num = i + 1;
        ds->array[i][lastIndex].radius = (i + 1) * 5;
        ds->array[i][lastIndex].phi = ds->array[i][1].phi; // getting the first phi in the array sorted by z
        ds->array[i][lastIndex].z = trapezoid_edges[i]; //here we want x.0001

        //now factors in the addition of both boundary points because n_points previously was counting true point additions, and did not count the blank index 0.
        ds->n_points[i] += 2;    

        // adjusting positions using insertion sort techniques as opposed to sorting the entire array. 
        // we have the guarentee from importData that the array was sorted
        // assigned points to indices first to avoid risk of comparing uninitialized "blank" points.
        // as opposed to full sorting algorithms like mergesort, each call here is O(N) and has the potential to escape much earlier. 
        adjustPointPositionFront(ds->array[i], ds->n_points[i], 0); // adjust the start boundary
        adjustPointPositionBack(ds->array[i], ds->n_points[i], lastIndex); // adjust the end boundary
    }

}


void initParallelogram(Parallelogram *pg, int layer_numI, float z1_minI, float z1_maxI,
                       float shadow_bottomL_jRI, float shadow_bottomR_jRI,
                       float shadow_bottomL_jLI, float shadow_bottomR_jLI,
                       float pSlopeI)
{
    pg->layer_num = layer_numI;
    pg->pSlope = pSlopeI;

    pg->shadow_bottomL_jR = shadow_bottomL_jRI;
    pg->shadow_bottomR_jR = shadow_bottomR_jRI;

    pg->shadow_bottomL_jL = shadow_bottomL_jLI;
    pg->shadow_bottomR_jL = shadow_bottomR_jLI;

    pg->z1_min = z1_minI;
    pg->z1_max = z1_maxI;
}


void initWedgeSuperPoint(wedgeSuperPoint *wsp, Point *points, int pointCount)
{
    wsp->point_count = pointCount;
    wsp->min = FLT_MAX;
    wsp->max = -FLT_MAX;

    // more efficient approach
    // instead of making a temp array and then transferring contents, add the values directly
    // simultaneously, you can determine the min and max, as opposed to doing it after the fact
    for (int i = 0; i < pointCount; i++)
    {
        wsp->points[i] = points[i];
        wsp->z_values[i] = points[i].z;

        if (points[i].z < wsp->min)
            wsp->min = points[i].z;
        if (points[i].z > wsp->max)
            wsp->max = points[i].z;
    }
}

// operator overloading not allowed in C, write separate method to check equality
int areWedgeSuperPointsEqual(wedgeSuperPoint *wsp1, wedgeSuperPoint *wsp2)
{
    //return (wsp1->min == wsp2->min) && (wsp1->max == wsp2->max);
    const float tolerance = 0.0001;
    return (fabs(wsp1->min - wsp2->min) < tolerance) && (fabs(wsp1->max - wsp2->max) < tolerance);
}

void getParallelograms(wedgePatch *wp)
{
    float z1_min = max(wp->superpoints[0].min, -1 * trapezoid_edges[0]);
    float z1_max = min(wp->superpoints[0].max, trapezoid_edges[0]);

    if (z1_min > z1_max)
    {
        z1_min = trapezoid_edges[0] + 1;
        z1_max = z1_min;
    }

    int previous_count = wp->parallelogram_count;

    // the code now handles the case below
    // if (! wp->parallelogram_count <= wp->superpoint_count - 1 ) {
    //     printf("Instead of assigning a temp array, we are overwriting the first superpoint_count-1 elements in the parallelogam array. If the current number of elements in the array is greater than superpoint_count-1, then we will have remaining elements that need to be deleted to replicate the functionality correctly");
    //     //exit(8);
    // }
    wp->parallelogram_count = 0; // we want to start at index 0 regardless and overwrite any old elements in the array to replicate the functionality of assigning a temp array.
    for (int i = 1; i < wp->superpoint_count; i++)
    {
        int j = i + 1;

        float z_j_min = wp->superpoints[i].min;
        float z_j_max = wp->superpoints[i].max;

        float a = straightLineProjectorFromLayerIJtoK(wp, z1_min, z_j_max, 1, j, num_layers);
        float b = straightLineProjectorFromLayerIJtoK(wp, z1_max, z_j_max, 1, j, num_layers);
        float c = straightLineProjectorFromLayerIJtoK(wp, z1_min, z_j_min, 1, j, num_layers);
        float d = straightLineProjectorFromLayerIJtoK(wp, z1_max, z_j_min, 1, j, num_layers);

        float pSlope = (j != num_layers) ? parallelogramSlopes[j - 1] : INT_MAX;

        // directly assign the values to the array
        if (wp->parallelogram_count < MAX_PARALLELOGRAMS_PER_PATCH)
        {
            Parallelogram *p = &wp->parallelograms[wp->parallelogram_count++]; // making a pointer to the address of first empty element in the array
            p->layer_num = j;                                                  // then dereferencing and assigning values to the properties
            p->pSlope = pSlope;
            p->shadow_bottomL_jR = a;
            p->shadow_bottomR_jR = b;
            p->shadow_bottomL_jL = c;
            p->shadow_bottomR_jL = d;
            p->z1_min = z1_min;
            p->z1_max = z1_max;
        }
    }
}

void wedgePatch_init(wedgePatch *wp, wedgeSuperPoint *superpointsI, int superpoint_count, float apexZ0I)
{
    wp->end_layer = -1; 
    wp->left_end_layer = -1;
    wp->right_end_layer = -1;
    wp->left_end_lambdaZ = 0;
    wp->right_end_lambdaZ = 0;
    wp->apexZ0 = apexZ0I;

    wp->shadow_fromTopToInnermost_topL_jL = 0;
    wp->shadow_fromTopToInnermost_topL_jR = 0;
    wp->shadow_fromTopToInnermost_topR_jL = 0;
    wp->shadow_fromTopToInnermost_topR_jR = 0;

    for (size_t i = 0; i < superpoint_count; i++)
    {                                          // size_t objects should only be non-negative and are more performant than ints
        wp->superpoints[i] = superpointsI[i]; // wp->superpoints is an array of pointers. Making the elements point to the elements in superpointsI.
    }
    wp->superpoint_count = superpoint_count;

    getParallelograms(wp);
    // getParallelograms_v1(wp);
    get_acceptanceCorners(wp);
    get_end_layer(wp);
}

float straightLineProjectorFromLayerIJtoK(wedgePatch *wp, float z_i, float z_j, int i, int j, int k)
{
    float radius_i = 0;
    float radius_j = 0;
    float radius_k = 0;

    if (i == 0)
    {
        radius_i = 0;
    }
    else
    {
        radius_i = radii[i - 1]; 
    }
    if (j == 0)
    {
        radius_j = 0;
    }
    else
    {
        radius_j = radii[j - 1];
    }
    if (k == 0)
    {
        radius_k = 0;
    }
    else
    {
        radius_k = radii[k - 1];
    }

    float radii_leverArm = (radius_k - radius_i) / (radius_j - radius_i);

    return z_i + (z_j - z_i) * radii_leverArm;
}

float straightLineProjector(float z_top, float z_j, int j)
{
    float temp = radii_leverArm[j - 1];
    return z_top - (z_top - z_j) * temp;
}

void getShadows(wedgePatch *wp, float zTopMin, float zTopMax)
{
    float zTop_min;
    float zTop_max;
    if (num_layers - 1 < 0)
    {
        zTop_min = zTopMin;
        zTop_max = zTopMax;
    }
    else
    {
        zTop_min = max(zTopMin, -trapezoid_edges[num_layers - 1]);
        zTop_max = min(zTopMax, trapezoid_edges[num_layers - 1]);
    }

    float topL_jL[MAX_SUPERPOINTS_IN_PATCH - 1];
    float topL_jR[MAX_SUPERPOINTS_IN_PATCH - 1];
    float topR_jL[MAX_SUPERPOINTS_IN_PATCH - 1];
    float topR_jR[MAX_SUPERPOINTS_IN_PATCH - 1];

    for (int i = 0; i < wp->superpoint_count - 1; ++i)
    {
        int j = i + 1;
        float z_j_min = wp->superpoints[i].min;
        float z_j_max = wp->superpoints[i].max;

        topL_jL[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_min, z_j_min, num_layers, j, 1);
        topL_jR[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_min, z_j_max, num_layers, j, 1);
        topR_jL[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_max, z_j_min, num_layers, j, 1);
        topR_jR[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_max, z_j_max, num_layers, j, 1);
    }

    wp->shadow_fromTopToInnermost_topL_jL = topL_jL[0];
    wp->shadow_fromTopToInnermost_topL_jR = topL_jR[0];
    wp->shadow_fromTopToInnermost_topR_jL = topR_jL[0];
    wp->shadow_fromTopToInnermost_topR_jR = topR_jR[0];

    // finding max in each of the respective arrays and saving to designated instance variables
    for (int i = 1; i < wp->superpoint_count - 1; ++i)
    {
        if (topL_jL[i] > wp->shadow_fromTopToInnermost_topL_jL)
        {
            wp->shadow_fromTopToInnermost_topL_jL = topL_jL[i];
        }
        if (topL_jR[i] < wp->shadow_fromTopToInnermost_topL_jR)
        {
            wp->shadow_fromTopToInnermost_topL_jR = topL_jR[i];
        }
        if (topR_jL[i] > wp->shadow_fromTopToInnermost_topR_jL)
        {
            wp->shadow_fromTopToInnermost_topR_jL = topR_jL[i];
        }
        if (topR_jR[i] < wp->shadow_fromTopToInnermost_topR_jR)
        {
            wp->shadow_fromTopToInnermost_topR_jR = topR_jR[i];
        }
    }
}

void get_acceptanceCorners(wedgePatch *wp)
{
    wp->squareAcceptance = true;
    wp->flatTop = true;
    wp->flatBottom = true;
    wp->triangleAcceptance = false;

    float a_corner_min = FLT_MAX;
    float b_corner_min = FLT_MAX;
    float c_corner_max = -FLT_MAX;
    float d_corner_max = -FLT_MAX;

    // getting min or max corners in all parallelograms
    for (int i = 0; i < wp->parallelogram_count; ++i)
    {
        Parallelogram *pg = &wp->parallelograms[i];
        if (pg->shadow_bottomL_jR < a_corner_min)
        {
            a_corner_min = pg->shadow_bottomL_jR;
        }
        if (pg->shadow_bottomR_jR < b_corner_min)
        {
            b_corner_min = pg->shadow_bottomR_jR;
        }
        if (pg->shadow_bottomL_jL > c_corner_max)
        {
            c_corner_max = pg->shadow_bottomL_jL;
        }
        if (pg->shadow_bottomR_jL > d_corner_max)
        {
            d_corner_max = pg->shadow_bottomR_jL;
        }
    }

    // assigning to the size-2 corner arrays
    wp->a_corner[0] = wp->parallelograms[0].z1_min;
    wp->a_corner[1] = a_corner_min;
    wp->b_corner[0] = wp->parallelograms[0].z1_max;
    wp->b_corner[1] = b_corner_min;
    wp->c_corner[0] = wp->parallelograms[0].z1_min;
    wp->c_corner[1] = c_corner_max;
    wp->d_corner[0] = wp->parallelograms[0].z1_max;
    wp->d_corner[1] = d_corner_max;

    // the nth element of shadow_bottom is the same as the nth element in the corner lists in CPP
    if (a_corner_min != wp->parallelograms[num_layers - 2].shadow_bottomL_jR)
    {
        wp->squareAcceptance = false;
        wp->flatTop = false;
    }
    if (b_corner_min != wp->parallelograms[num_layers - 2].shadow_bottomR_jR)
    {
        wp->squareAcceptance = false;
        wp->flatTop = false;
    }
    if (c_corner_max != wp->parallelograms[num_layers - 2].shadow_bottomL_jL)
    {
        wp->squareAcceptance = false;
        wp->flatBottom = false;
    }
    if (d_corner_max != wp->parallelograms[num_layers - 2].shadow_bottomR_jL)
    {
        wp->squareAcceptance = false;
        wp->flatBottom = false;
    }

    // adjusting corners for triangle acceptance
    if (wp->c_corner[1] > wp->a_corner[1])
    {
        wp->triangleAcceptance = true;
        wp->c_corner[1] = wp->b_corner[1];
        wp->a_corner[1] = wp->b_corner[1];
    }

    if (wp->b_corner[1] < wp->d_corner[1])
    {
        wp->triangleAcceptance = true;
        wp->b_corner[1] = wp->c_corner[1];
        wp->d_corner[1] = wp->c_corner[1];
    }
}

void get_end_layer(wedgePatch *wp)
{
    // naming counterintuitive
    float lambdaZLeftMax = -1 * INT_MAX + 2;
    float lambdaZRightMin = INT_MAX - 2;
    // assigning -1 to start because they should only hold non-negative integers
    wp->left_end_layer = -1;
    wp->right_end_layer = -1;

    // combined two independent loops
    for (int i = 0; i < num_layers; i++)
    {
        float lambdaZ_left = (wp->superpoints[i].min - wp->apexZ0) / radii[i];
        float lambdaZ_right = (wp->superpoints[i].max - wp->apexZ0) / radii[i];

        if (lambdaZ_left > lambdaZLeftMax)
        {
            wp->left_end_layer = i;
            lambdaZLeftMax = lambdaZ_left;
        }

        if (lambdaZ_right < lambdaZRightMin)
        {
            wp->right_end_layer = i;
            lambdaZRightMin = lambdaZ_right;
        }
    }
    // no need to find min/max of an array, we already have that data
    wp->left_end_lambdaZ = lambdaZLeftMax;
    wp->right_end_lambdaZ = lambdaZRightMin;
}


void initWedgeCover(wedgeCover *wc)
{
    wc->n_patches = 0;
    //wc->data = dataI;
}

void add_patch(wedgeCover *cover, wedgePatch *curr_patch)
{
    if (cover->n_patches == 0)
    {
        cover->patches[0] = *curr_patch; //copy the patch directly
        //cover->all_patches[0] = curr_patch;
        #if KEEP_DELETED_PATCHES == true
            cover->real_patch_list[0] = true;
        #endif
        cover->n_patches = 1;
    }
    else
    {
        wedgePatch *prev_patch = &(cover->patches[cover->n_patches - 1]);
        bool different = false;

        for (index_type i = 0; i < prev_patch->superpoint_count; i++)
        {
            if ((prev_patch->superpoints[i].min != curr_patch->superpoints[i].min) ||
                (prev_patch->superpoints[i].max != curr_patch->superpoints[i].max))
            {
                different = true;
                break;
            }
        }

        // if the min and max are the same for each superpoint, don't add a patch
        if (different)
        {
            if (cover->n_patches < MAX_PATCHES)
            {
                cover->patches[cover->n_patches] = *curr_patch;
                //cover->all_patches[cover->n_patches] = curr_patch;
                #if KEEP_DELETED_PATCHES == true
                    cover->real_patch_list[cover->n_patches] = true;
                #endif
                cover->n_patches += 1;
            }
        }
    }
}

void delete_patch(wedgeCover *cover, int index)
{
    if (index < 0 || index >= cover->n_patches)
    {
        return;
    }
    #if KEEP_DELETED_PATCHES == true
        cover->real_patch_list[index] = false;
    #endif
    for (index_type i = index; i < cover->n_patches - 1; i++)
    {
        cover->patches[i] = cover->patches[i + 1];
        #if KEEP_DELETED_PATCHES == true
            cover->real_patch_list[i] = cover->real_patch_list[i + 1];
        #endif
    }

    // resetting the last elements
    memset(&cover->patches[cover->n_patches - 1], 0, sizeof(wedgePatch));
    #if KEEP_DELETED_PATCHES == true
        cover->real_patch_list[cover->n_patches - 1] = false;
    #endif

    cover->n_patches -= 1;
}

// can't provide default parameters
index_type get_index_from_z(int layer, float z_value)
{
    // c doesn't support string comparison directly, using integer comparison for effiency
    // CLOSEST = 11, ABOVE = 12, BELOW = 13
    float minVal = 1000000;
    index_type index = 0;

    for (index_type i = 0; i < Gdata.n_points[layer]; i++)
    {
        float diff = fabs(Gdata.array[layer][i].z - z_value); // absolute difference
        if (diff < minVal)
        {
            minVal = diff;
            index = i;
        }
    }

    //alignment always equals closest so we can just return index
    return index;
}

// not implementing the logic corresponding to show = true (that would involve Line Generators, etc)
// Line Generator and its accompanying methods have been coded, but we are not going to implement show=true case here as main method passes with show=false.
// lining is always MAKE_PATCHES_SHADOW_QUILT_FROM_EDGES. assumes this in code
void solve(wedgeCover *cover, float apexZ0, int ppl, int nlines, bool leftRight)
{
    for (index_type i = 0; i < num_layers; i++)
    {
        bool foundIdentical = false;
        bool firstTime = true;

        while (foundIdentical || firstTime)
        {
            foundIdentical = false;
            for (index_type x = 0; x < Gdata.n_points[i] - 1; x++)
            {
                if (Gdata.array[i][x].z == Gdata.array[i][x + 1].z)
                {
                    Gdata.array[i][x + 1].z += 0.00001;
                    foundIdentical = true;
                }
            }

            firstTime = false;
            if (foundIdentical)
            {
                qsort(Gdata.array[i], Gdata.n_points[i], sizeof(Point), comparePoints); // Chip will ultimately have sorted data coming through, not needed for synthesis
            }
        }
    }
    makePatches_ShadowQuilt_fromEdges(cover, apexZ0, 1, ppl, leftRight);
}

void makePatches_ShadowQuilt_fromEdges(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight) // TOP-LEVEL FUNCTION FOR VITIS
{
    bool fix42 = true;
    apexZ0 = trapezoid_edges[0];
    float saved_apexZ0;

    while (apexZ0 > -1 * trapezoid_edges[0]) //consider how this works when we are expanding instead of retracting the trapezoid_edges
    {
        apexZ0 = solveNextColumn(cover, apexZ0, stop, ppl, leftRight, fix42, saved_apexZ0); 
        saved_apexZ0 = apexZ0; 
    }
}

float solveNextColumn(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight, bool fix42, float saved_apexZ0)
{
    float z_top_min = -1 * top_layer_lim;

    float complementary_apexZ0 = 0;
    index_type first_row_count = 0;
    float c_corner = LONG_MAX;

    float z_top_max = top_layer_lim;

    if (cover->n_patches > 0)
    {
        z_top_max = min(z_top_max, straightLineProjectorFromLayerIJtoK(&cover->patches[cover->n_patches - 1], -1 * beam_axis_lim, apexZ0, 0, 1, num_layers // includes passing a pointer to the last patch
                        ));
    }

    index_type nPatchesInColumn = 0;
    float projectionOfCornerToBeam = 0;
    //returnArray[6] = {nPatchesInColumn, c_corner, projectionOfCornerToBeam, z_top_min, z_top_max, complementary_apexZ0}
    //remove nPatchesInColumn once debugging finishes
    while((c_corner > -1 * trapezoid_edges[num_layers - 1]) && (nPatchesInColumn<100000000) && (projectionOfCornerToBeam < beam_axis_lim))
    {
        solveNextPatchPair(cover, apexZ0, stop, ppl, leftRight, fix42, saved_apexZ0, nPatchesInColumn, c_corner, projectionOfCornerToBeam, z_top_min, z_top_max, complementary_apexZ0); 
    }

    apexZ0 = cover->patches[cover->n_patches - 1].c_corner[0];
    apexZ0 = saved_apexZ0;
    printf("'=======================================================  z1_Align: %f\n", apexZ0);

    return apexZ0; 
}

void solveNextPatchPair(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight, bool fix42, float &saved_apexZ0, int &nPatchesInColumn, float &c_corner, float &projectionOfCornerToBeam, float &z_top_min, float &z_top_max, float &complementary_apexZ0)
{
    nPatchesInColumn++;
    printf("%f %d %f %d\n", apexZ0, ppl, z_top_max, leftRight);

    makePatch_alignedToLine(cover, apexZ0, z_top_max, ppl, false, false);

    index_type lastPatchIndex = cover->n_patches - 1;

    printf("top layer from %f to %f z_top_max: %f\n",
            cover->patches[lastPatchIndex].superpoints[num_layers - 1].max,
            cover->patches[lastPatchIndex].superpoints[num_layers - 1].min,
            z_top_max);
    printf("original: [%f, %f] for patch %d\n",
            cover->patches[lastPatchIndex].a_corner[0],
            cover->patches[lastPatchIndex].a_corner[1],
            cover->n_patches);
    printf("original: [%f, %f]\n",
            cover->patches[lastPatchIndex].b_corner[0],
            cover->patches[lastPatchIndex].b_corner[1]);

    printf("original: [%f, %f]\n",
            cover->patches[lastPatchIndex].c_corner[0],
            cover->patches[lastPatchIndex].c_corner[1]);

    printf("original: [%f, %f]\n",
            cover->patches[lastPatchIndex].d_corner[0],
            cover->patches[lastPatchIndex].d_corner[1]);

    for (index_type i = 1; i < cover->patches[lastPatchIndex].superpoint_count - 1; i++)
    {
        index_type j = i + 1;
        printf("%d superpoint: %f %f shadowTop from L1Max: %f %f from L1 min: %f %f\n",
                j,
                cover->patches[lastPatchIndex].superpoints[i].min,
                cover->patches[lastPatchIndex].superpoints[i].max,
                straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].max, cover->patches[lastPatchIndex].superpoints[i].min, 1, j, num_layers),
                straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].max, cover->patches[lastPatchIndex].superpoints[i].max, 1, j, num_layers),
                straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].min, cover->patches[lastPatchIndex].superpoints[i].min, 1, j, num_layers),
                straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].min, cover->patches[lastPatchIndex].superpoints[i].max, 1, j, num_layers));
    }

    float original_c = cover->patches[lastPatchIndex].c_corner[1];
    float original_d = cover->patches[lastPatchIndex].d_corner[1];

    c_corner = original_c;

    bool repeat_patch = false;
    bool repeat_original = false;

    // code written assuming number of layers is 5.
    /*
    if (cover->n_patches > 2) {
        int thirdLastPatchIndex = lastPatchIndex - 2;
        repeat_original = (cover->patches[lastPatchIndex].superpoints[num_layers - 1] == cover->patches[thirdLastPatchIndex].superpoints[num_layers - 1]) &&
                (cover->patches[lastPatchIndex].superpoints[0] == cover->patches[thirdLastPatchIndex].superpoints[0]) &&
                (cover->patches[lastPatchIndex].superpoints[1] == cover->patches[thirdLastPatchIndex].superpoints[1]) &&
                (cover->patches[lastPatchIndex].superpoints[2] == cover->patches[thirdLastPatchIndex].superpoints[2]) &&
                (cover->patches[lastPatchIndex].superpoints[3] == cover->patches[thirdLastPatchIndex].superpoints[3]);
    }
    */
    // dynamic version below
    if (cover->n_patches > 2)
    {
        index_type thirdLastPatchIndex = lastPatchIndex - 2;
        repeat_original = true; // assume they are the same initially
        for (index_type i = 0; i < MAX_SUPERPOINTS_IN_PATCH; i++)
        { // iterating over the first (five) superpoints
            if (!areWedgeSuperPointsEqual(&cover->patches[lastPatchIndex].superpoints[i], &cover->patches[thirdLastPatchIndex].superpoints[i]))
            {
                repeat_original = false; // if any pair of superpoints don't match, set to false
                break;                   // no need to check further if a mismatch is found
            }
        }
    }

    float seed_apexZ0 = apexZ0;
    wedgePatch *lastPatch = &cover->patches[lastPatchIndex];
    projectionOfCornerToBeam = straightLineProjectorFromLayerIJtoK(lastPatch, lastPatch->c_corner[1], lastPatch->c_corner[0], num_layers, 1, 0);
    bool squarePatch_alternate1 = (lastPatch->a_corner[1] > z_top_max) && (lastPatch->b_corner[1] > z_top_max) && (lastPatch->flatBottom);
    bool squarePatch_alternate2 = (lastPatch->a_corner[1] > z_top_max) && (lastPatch->flatBottom);

    bool notChoppedPatch = (lastPatch->squareAcceptance) || squarePatch_alternate1 || squarePatch_alternate2;
    bool madeComplementaryPatch = false;

    int nPatchesAtOriginal = cover->n_patches;

    printf("squareAcceptance: %d triangleAcceptance: %d projectionOfCornerToBeam: %f notChoppedPatch %d\n",
            lastPatch->squareAcceptance, lastPatch->triangleAcceptance,projectionOfCornerToBeam, notChoppedPatch);

    if (!notChoppedPatch && (lastPatch->c_corner[1] > -1 * trapezoid_edges[num_layers - 1]) && ((projectionOfCornerToBeam < beam_axis_lim)))
    {
        complementary_apexZ0 = lastPatch->superpoints[0].min;
        if (lastPatch->triangleAcceptance && !repeat_original)
        {
            z_top_min = lastPatch->d_corner[1];
        }
        else
        {
            printf("z_top_min before: %f superpoints[self.env.num_layers-1].min: %f\n", z_top_min, lastPatch->superpoints[num_layers - 1].min);
            z_top_min = max(-1 * top_layer_lim, lastPatch->superpoints[num_layers - 1].min);
        }

        makePatch_alignedToLine(cover, complementary_apexZ0, z_top_min, ppl, true, false);
        // updating the last patch index because makePatch_alignedToLine will add more patches to the patches array. Should revisit after writing method
        // makePatch_alignedToLine will call the add patch method, so we must get a new last patch index.
        lastPatchIndex = cover->n_patches - 1;

        madeComplementaryPatch = true;
        printf("complementary: [%f, %f] for z_top_min: %f\n", cover->patches[lastPatchIndex].a_corner[0], cover->patches[lastPatchIndex].a_corner[1], z_top_min);
        printf("complementary: [%f, %f] for patch %d\n", cover->patches[lastPatchIndex].b_corner[0], cover->patches[lastPatchIndex].b_corner[1], cover->n_patches);
        printf("complementary: [%f, %f]\n", cover->patches[lastPatchIndex].c_corner[0], cover->patches[lastPatchIndex].c_corner[1]);
        printf("complementary: [%f, %f]\n", cover->patches[lastPatchIndex].d_corner[0], cover->patches[lastPatchIndex].d_corner[1]);

        float complementary_a = cover->patches[lastPatchIndex].a_corner[1];
        float complementary_b = cover->patches[lastPatchIndex].b_corner[1];

        float white_space_height = max(original_c - complementary_a, original_d - complementary_b);
        float previous_white_space_height = -1;
        int counter = 0;
        int counterUpshift = 0;
        index_type current_z_top_index = -1;
        float previous_z_top_min = -999;

        while (!(white_space_height <= 0.0000005 && (previous_white_space_height >= 0)) && (fabs(white_space_height) > 0.000005) &&
                ((cover->patches[lastPatchIndex].c_corner[1] > -1 * trapezoid_edges[num_layers - 1]) ||
                (white_space_height > 0.000005)) &&
                (current_z_top_index < (int)(Gdata.n_points[num_layers - 1])) &&
                !(repeat_patch) && !(repeat_original))
        {
            solveComplmentaryPatch(cover, previous_white_space_height, ppl, fix42, nPatchesAtOriginal, previous_z_top_min, complementary_apexZ0, white_space_height, lastPatchIndex, original_c, original_d, complementary_a, complementary_b, current_z_top_index, counter, counterUpshift, z_top_min, repeat_patch, repeat_original);
        }
    }

    lastPatchIndex = cover->n_patches - 1; // just to keep fresh in case we use it
   c_corner = cover->patches[lastPatchIndex].c_corner[1];

   projectionOfCornerToBeam = straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex],c_corner, cover->patches[lastPatchIndex].c_corner[0], num_layers, 1, 0);

    saved_apexZ0 = cover->patches[lastPatchIndex].c_corner[0];

    if (madeComplementaryPatch) // Create separate function for this
    {
        int secondLastPatchIndex = lastPatchIndex - 1;

        // modifying patches, not adding patches, so index variables do not need to be updated.
        getShadows(&cover->patches[lastPatchIndex],z_top_min, z_top_max);
        getShadows(&cover->patches[secondLastPatchIndex],z_top_min, z_top_max);

        float original_topR_jL = cover->patches[secondLastPatchIndex].shadow_fromTopToInnermost_topR_jL;
        bool originalPartialTop = (original_topR_jL > complementary_apexZ0) && (original_topR_jL < apexZ0) &&
                                    (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[secondLastPatchIndex], original_topR_jL, z_top_max, 1, num_layers, 0)) < 20 * beam_axis_lim);

        float original_topL_jL = cover->patches[secondLastPatchIndex].shadow_fromTopToInnermost_topL_jL;
        
        bool originalPartialBottom = (original_topL_jL > complementary_apexZ0) && ((original_topL_jL - apexZ0) < -0.0001) &&
                                        (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[secondLastPatchIndex], original_topL_jL,z_top_min, 1, num_layers, 0)) < 20 * beam_axis_lim);                
        
        float complementary_topR_jR = cover->patches[lastPatchIndex].shadow_fromTopToInnermost_topR_jR;
        
        bool complementaryPartialTop = (complementary_topR_jR > complementary_apexZ0) && (complementary_topR_jR < apexZ0) &&
                                        (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_topR_jR, z_top_max, 1, num_layers, 0)) < 20 * beam_axis_lim);

        float complementary_topL_jR = cover->patches[lastPatchIndex].shadow_fromTopToInnermost_topL_jR;
        
        bool complementaryPartialBottom = (complementary_topL_jR > complementary_apexZ0) && ((complementary_topL_jR - apexZ0) < -0.0001) &&
                                            (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_topL_jR,z_top_min, 1, num_layers, 0)) < 20 * beam_axis_lim);

        float horizontalShiftTop = original_topR_jL - complementary_topR_jR;
        float horizontalShiftBottom = original_topL_jL - complementary_topL_jR;

        float complementary_topR_jL = cover->patches[lastPatchIndex].shadow_fromTopToInnermost_topR_jL;
        float complementary_topL_jL = cover->patches[lastPatchIndex].shadow_fromTopToInnermost_topL_jL;
        float original_topR_jR = cover->patches[secondLastPatchIndex].shadow_fromTopToInnermost_topR_jR;
        float original_topL_jR = cover->patches[secondLastPatchIndex].shadow_fromTopToInnermost_topL_jR;

        float horizontalOverlapTop = max(complementary_topR_jL - original_topR_jL, complementary_topR_jR - original_topR_jR);
        float horizontalOverlapBottom = max(complementary_topL_jL - original_topL_jL, complementary_topL_jR - original_topL_jR);

        horizontalOverlapTop = -1;
        horizontalOverlapBottom = -1;
        float newGapTop = -0.000001;
        float newGapBottom = -0.000001;

        bool makeHorizontallyShiftedPatch = false;
        float shifted_Align = apexZ0;
        bool doShiftedPatch = true;

        float newZtop = 0;

        float z0_original_bCorner = straightLineProjectorFromLayerIJtoK(&cover->patches[secondLastPatchIndex], apexZ0, z_top_max, 1, num_layers, 0);
        float z0_complementary_cCorner = straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_apexZ0,z_top_min, 1, num_layers, 0);
        bool shiftOriginal = true;

        if (z0_original_bCorner < 0)
        {
            shiftOriginal = false;
            shifted_Align = complementary_apexZ0;
        }

        if (z0_complementary_cCorner > 0)
        {
            shiftOriginal = true;
            shifted_Align = apexZ0;
        }

        //if (horizontalShiftTop > 0 || horizontalShiftBottom > 0)
        if (horizontalShiftTop > 0.000001 || horizontalShiftBottom > 0) // NOTE THAT horizontalShiftTop > 0.000001 is a "hack" to avoid infinite loop from Wedge 42 in this condition and the next
        {
            printf("originalPartialTop: %d complementaryPartialTop: %d originalPartialBottom: %d complementaryPartialBottom: %d %f %f %f %f horizontalOverlapTop: %f horizontalOverlapBottom: %f\n",
                    originalPartialTop, complementaryPartialTop, originalPartialBottom, complementaryPartialBottom,
                    original_topR_jL, original_topL_jL, complementary_topR_jR, complementary_topL_jR,
                    horizontalOverlapTop, horizontalOverlapBottom);
        }

        while ((((horizontalShiftTop > 0.000001) && originalPartialTop && complementaryPartialTop) || ((horizontalShiftBottom > 0.000001) && originalPartialBottom && complementaryPartialBottom)) && doShiftedPatch && (horizontalOverlapTop <= 0) && (horizontalOverlapBottom <= 0) && ((newGapTop < 0) || (newGapBottom < 0)))
        {
            printf("horizontalShifts: %f %f shifted_Align: %f\n", horizontalShiftTop, horizontalShiftBottom, shifted_Align);

            newZtop = z_top_max;

            if (shiftOriginal)
            {
                shifted_Align -= max(horizontalShiftTop, horizontalShiftBottom);
            }
            else
            {
                shifted_Align += max(horizontalShiftTop, horizontalShiftBottom);
                newZtop = z_top_min;
            }

            if (makeHorizontallyShiftedPatch)
            {
                delete_patch(cover, cover->n_patches - 1);
                // decrement n_patches is handled by delete_patch
            }

            makePatch_alignedToLine(cover, shifted_Align, newZtop, ppl, !shiftOriginal, false);

            getShadows(&cover->patches[cover->n_patches - 1], z_top_min, z_top_max);

            if (shiftOriginal)
            {
                original_topR_jL = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topR_jL;
                original_topL_jL = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topL_jL;
                original_topR_jR = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topR_jR;
                original_topL_jR = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topL_jR;
            }
            else
            {
                complementary_topR_jR = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topR_jR;
                complementary_topL_jR = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topL_jR;
                complementary_topR_jL = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topR_jL;
                complementary_topL_jL = cover->patches[cover->n_patches - 1].shadow_fromTopToInnermost_topL_jL;
            }

            horizontalShiftTop = original_topR_jL - complementary_topR_jR;
            horizontalShiftBottom = original_topL_jL - complementary_topL_jR;

            if (shiftOriginal && straightLineProjectorFromLayerIJtoK(&cover->patches[cover->n_patches - 1], original_topR_jR, z_top_max, 1, num_layers, 0) < beam_axis_lim)
            {
                horizontalOverlapTop = max(complementary_topR_jL - original_topR_jL, complementary_topR_jR - original_topR_jR);
                horizontalOverlapBottom = max(complementary_topL_jL - original_topL_jL, complementary_topL_jR - original_topL_jR);
                printf(" horizontalOverlapTop: %f horizontalOverlapBottom: %f\n", horizontalOverlapTop, horizontalOverlapBottom);
            }

            printf("original_topR_jL: %f complementary_topR_jR %f original_topL_jL %f complementary_topL_jR %f shiftOriginal %d\n",
                    original_topR_jL, complementary_topR_jR, original_topL_jL, complementary_topL_jR, shiftOriginal);

            makeHorizontallyShiftedPatch = true;

            printf("updated_horizontalShifts: %f %f shifted_Align: %f\n", horizontalShiftTop, horizontalShiftBottom, shifted_Align);
        }
        if (makeHorizontallyShiftedPatch)
        {
            if ((straightLineProjectorFromLayerIJtoK(&cover->patches[cover->n_patches - 1], shifted_Align, newZtop, 1, num_layers, 0) > beam_axis_lim) && shiftOriginal)
            {
                if (cover->n_patches > 2)
                {
                    delete_patch(cover, cover->n_patches - 3);
                }
            }
        }
    }

    z_top_max = c_corner;

    printf("+++++++++++++++++++++++ c_corner: %f\n", c_corner);
}

void solveComplmentaryPatch(wedgeCover *cover, float &previous_white_space_height, int ppl, bool fix42, int nPatchesAtOriginal, float &previous_z_top_min, float complementary_apexZ0, float &white_space_height, index_type &lastPatchIndex, float original_c, float original_d, float &complementary_a, float &complementary_b, index_type &current_z_top_index, int &counter, int &counterUpshift, float &z_top_min, bool &repeat_patch, bool &repeat_original)
{
    printf("\n");
    if (cover->n_patches > 2)
    {
        index_type secondLastPatchIndex = lastPatchIndex - 1;
        printf("original c: %f %f || original d: %f %f\n",
                original_c, cover->patches[secondLastPatchIndex].c_corner[1],
                original_d, cover->patches[secondLastPatchIndex].d_corner[1]);
    }
    printf("complementary_a: %f %f || complementary_b: %f %f\n",
            complementary_a, cover->patches[lastPatchIndex].a_corner[1],
            complementary_b, cover->patches[lastPatchIndex].b_corner[1]);

    current_z_top_index = get_index_from_z(num_layers - 1,z_top_min); 
    printf("current white_space_height: %f\n", white_space_height);
    printf("counter: %d counterUpshift: %d\n", counter, counterUpshift);
    printf("orig_ztop: %d orig_z_top_min: %f\n", current_z_top_index,z_top_min);

    index_type current_z_i_index[MAX_LAYERS];
    index_type new_z_i_index[MAX_LAYERS];

    for (index_type i = 0; i < num_layers; i++)
    {
        current_z_i_index[i] = get_index_from_z(i, straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_apexZ0,z_top_min, 1, num_layers, i + 1));
    }

    if (z_top_min == previous_z_top_min)
    {
        current_z_top_index += 1;
        for (index_type i = 0; i < num_layers; i++)
        {
            new_z_i_index[i] = current_z_i_index[i] + 1;
        }
    }

    previous_z_top_min = z_top_min;

    if (white_space_height < 0)
    {
        counter += 1;
        current_z_top_index -= 1;
        for (index_type i = 0; i < num_layers; i++)
        {
            new_z_i_index[i] = current_z_i_index[i] - 1;
        }
    }
    else
    {
        counterUpshift += 1;
        current_z_top_index += 1;
        for (index_type i = 0; i < num_layers; i++)
        {
            new_z_i_index[i] = current_z_i_index[i] + 1;
        }
    }

    int x = Gdata.n_points[num_layers - 1] - 1;
    current_z_top_index = min(current_z_top_index, Gdata.n_points[num_layers - 1] - 1); // n_points is an array of the sizes of each element of 'array'

    for (index_type i = 0; i < num_layers; i++)
    {
        new_z_i_index[i] = min(new_z_i_index[i], (float)Gdata.n_points[i] - 1);
    }

    for (index_type i = 0; i < num_layers; i++)
    { 
        new_z_i_index[i] = max(new_z_i_index[i], 0.0f);
    }
    float new_z_i[MAX_LAYERS];

    for (index_type i = 0; i < num_layers; i++)
    {
        new_z_i[i] = Gdata.array[i][new_z_i_index[i]].z;
    }

    float new_z_i_atTop[MAX_LAYERS - 1]; // note: the size is MAX_LAYERS - 1 because the loop starts from 1
    for (index_type i = 1; i < num_layers; i++)
    {
        new_z_i_atTop[i - 1] = straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex],
                                                                    complementary_apexZ0,
                                                                    new_z_i[i],
                                                                    1,
                                                                    i + 1,
                                                                    num_layers);
    }

    index_type layerWithSmallestShift = 0;
    float layerSMin = FLT_MAX;

    for (index_type i = 0; i < num_layers - 1; i++)
    {
        if (fabs(new_z_i_atTop[i] - previous_z_top_min) < layerSMin)
        { // fabs is for floats. abs is only int
            layerSMin = fabs(new_z_i_atTop[i] - previous_z_top_min);
            layerWithSmallestShift = i;
        }
    }

    layerWithSmallestShift += 1;

    for (index_type i = 0; i < num_layers - 1; i++)
    {
        printf("%u new_z_i_atTop: %f shift_i_ztop: %f layerWithSmallestShift: %u\n",
                i + 1, new_z_i_atTop[i], new_z_i_atTop[i] - previous_z_top_min, layerWithSmallestShift + 1);
    }

    z_top_min = Gdata.array[num_layers - 1][current_z_top_index].z;
    z_top_min = new_z_i_atTop[layerWithSmallestShift - 1];

    if (fabs(z_top_min - previous_z_top_min) < 0.000001)
    {
        z_top_min = Gdata.array[num_layers - 1][current_z_top_index].z;
    }

    if (fabs(z_top_min - previous_z_top_min) < 0.000001)
    {
        z_top_min = Gdata.array[num_layers - 2][current_z_top_index].z;
    }

    if (fabs(z_top_min - previous_z_top_min) < 0.000001)
    {
        z_top_min = Gdata.array[num_layers - 3][current_z_top_index].z;
    }

    if (((z_top_min - previous_z_top_min) * white_space_height) < 0)
    {
        z_top_min = new_z_i_atTop[num_layers - 2];
    }

    printf(" new_def_z_top_min_diff: %f\n", z_top_min - Gdata.array[num_layers - 1][current_z_top_index].z);

    printf(" new_ztop_index: %d new_z_i_index: ", current_z_top_index);
    for (index_type i = 0; i < num_layers; i++)
    {
        printf("%u ", new_z_i_index[i]);
    }
    printf("new_z_top_min: %f shift_ztop: %f\n", z_top_min, z_top_min - previous_z_top_min);

    int nPatchesAtComplementary = cover->n_patches;
    lastPatchIndex = cover->n_patches - 1; // this may have already been updated at the end of the last call, but just to be sure
    if (nPatchesAtComplementary > nPatchesAtOriginal)
    {
        printf("deleted complementary: [%f, %f] for patch %d\n",
                cover->patches[lastPatchIndex].a_corner[0],
                cover->patches[lastPatchIndex].a_corner[1],
                cover->n_patches);
        printf("deleted complementary: [%f, %f]\n",
                cover->patches[lastPatchIndex].b_corner[0],
                cover->patches[lastPatchIndex].b_corner[1]);
        printf("deleted complementary: [%f, %f]\n",
                cover->patches[lastPatchIndex].c_corner[0],
                cover->patches[lastPatchIndex].c_corner[1]);
        printf("deleted complementary: [%f, %f]\n",
                cover->patches[lastPatchIndex].d_corner[0],
                cover->patches[lastPatchIndex].d_corner[1]);

        // Call delete_patch to remove the last patch
        delete_patch(cover, lastPatchIndex);
        // no need to manually decrement n_patches, delete_patch will handle it
    }
    lastPatchIndex = cover->n_patches - 1; // lastPatchIndex has changed because of the delete patch
    // it may be not needed to update lastPatchIndex, but for now, I did it, so it wouldn't be forgotten later.

    // call makePatch_alignedToLine to add a new patch based on the complementary apex and top z values.
    makePatch_alignedToLine(cover, complementary_apexZ0, z_top_min, ppl, true, false);
    // update the lastPatchIndex to point to the newly added patch.
    lastPatchIndex = cover->n_patches - 1;

    // retrieve the a and b corner values from the latest patch.
    complementary_a = cover->patches[lastPatchIndex].a_corner[1];
    complementary_b = cover->patches[lastPatchIndex].b_corner[1];

    // update the previous white space height for the next iteration.
    previous_white_space_height = white_space_height;
    // calculate the new white space height based on the original and complementary corners.
    white_space_height = max(original_c - complementary_a, original_d - complementary_b);

    printf("complementary_a: %f %f || complementary_b: %f %f new z_top_min: %f\n",
            complementary_a, cover->patches[lastPatchIndex].a_corner[1],
            complementary_b, cover->patches[lastPatchIndex].b_corner[1],z_top_min);
    printf("new white_space_height: %f\n", white_space_height);
    printf("adjusted complementary: [%f, %f] for z_top_min: %f\n",
            cover->patches[lastPatchIndex].a_corner[0], cover->patches[lastPatchIndex].a_corner[1],z_top_min);
    printf("adjusted complementary: [%f, %f] for patch %d\n",
            cover->patches[lastPatchIndex].b_corner[0], cover->patches[lastPatchIndex].b_corner[1], cover->n_patches);
    printf("adjusted complementary: [%f, %f]\n",
            cover->patches[lastPatchIndex].c_corner[0], cover->patches[lastPatchIndex].c_corner[1]);
    printf("adjusted complementary: [%f, %f]\n",
            cover->patches[lastPatchIndex].d_corner[0], cover->patches[lastPatchIndex].d_corner[1]);

    if ((cover->n_patches > 3) && fix42)
    {
        index_type lastPatchIdx = cover->n_patches - 1;
        index_type thirdLastPatchIdx = lastPatchIdx - 2;

        // checking if the superpoints of the last and third last patches are the same
        repeat_patch = true;
        // turned this into a for loop, dynamic. if ((patches[patches.size() - 1].superpoints[env.num_layers - 1] == patches[patches.size() - 3].superpoints[env.num_layers - 1]) && (patches[patches.size() - 1].superpoints[0] == patches[patches.size() - 3].superpoints[0]) && (patches[patches.size() - 1].superpoints[1] == patches[patches.size() - 3].superpoints[1]) && (patches[patches.size() - 1].superpoints[2] == patches[patches.size() - 3].superpoints[2]) && (patches[patches.size() - 1].superpoints[3] == patches[patches.size() - 3].superpoints[3]))
        // that code checked 0 to 4
        for (index_type i = 0; i < num_layers; i++)
        {
            if (!areWedgeSuperPointsEqual(&cover->patches[lastPatchIdx].superpoints[i], &cover->patches[thirdLastPatchIdx].superpoints[i]))
            {
                repeat_patch = false;
                break;
            }
        }

        if (repeat_patch)
        {
            printf("%f %f repeat_patch: %d\n",
                    cover->patches[lastPatchIdx].superpoints[num_layers - 1].min,
                    cover->patches[lastPatchIdx].superpoints[num_layers - 1].max,
                    repeat_patch);

            delete_patch(cover, lastPatchIdx);

            current_z_top_index -= 1;

            z_top_min = Gdata.array[num_layers - 1][current_z_top_index].z;
            z_top_min = new_z_i_atTop[layerWithSmallestShift - 1];

            makePatch_alignedToLine(cover, complementary_apexZ0, z_top_min, ppl, true, false);
        }
    }
}

void makePatch_alignedToLine(wedgeCover *cover, float apexZ0, float z_top, int &ppl, bool leftRight, bool float_middleLayers_ppl)
{
    wedgeSuperPoint init_patch[MAX_LAYERS]; // correct
    int original_ppl = ppl;
    float alignmentAccuracy = 0.00001;
    // Point row_data[MAX_LAYERS][MAX_POINTS_FOR_DATASET];
    index_type init_patch_size = 0;

    for (index_type i = 0; i < num_layers; i++)
    {
        float y = radii[i];
        float row_list[MAX_POINTS_PER_LAYER];
        int row_list_size = 0;

        for (index_type j = 0; j < Gdata.n_points[i]; j++)
        {
            row_list[row_list_size++] = Gdata.array[i][j].z;
        }

        float r_max = radii[num_layers - 1];
        float projectionToRow = (z_top - apexZ0) * (y - radii[0]) / (r_max - radii[0]) + apexZ0;

        int start_index = 0;
        float start_value = 1000000;

        for (index_type j = 0; j < row_list_size; j++)
        {
            if (fabs(row_list[j] - projectionToRow) < fabs(start_value))
            {
                start_index = j;
                start_value = row_list[j] - projectionToRow;
            }
        }

        int left_bound = 0;
        float lbVal = INT_MAX;
        int right_bound = 0;
        float rbVal = INT_MAX;

        for (index_type j = 0; j < row_list_size; j++)
        {
            if (fabs((row_list[j] + trapezoid_edges[i])) < lbVal)
            {
                left_bound = j;
                lbVal = fabs((row_list[j] + trapezoid_edges[i]));
            }

            if (fabs((row_list[j] - trapezoid_edges[i])) < rbVal)
            {
                right_bound = j;
                rbVal = fabs((row_list[j] - trapezoid_edges[i]));
            }
        }

        if (float_middleLayers_ppl && i != 0 && i != num_layers - 1)
        {
            ppl = original_ppl * 2 - 1;
        }
        else
        {
            ppl = original_ppl;
        }

        Point temp[MAX_POINTS_PER_LAYER]; // check
        int temp_size = 0;

        if (leftRight)
        {
            if (start_index != 0 && start_value > alignmentAccuracy)
            {
                start_index -= 1;
            }
            // making and adding a new vector that is a subset of "row_data" or array, going from right+1-ppl to right+1?
            if ((start_index + ppl) > (right_bound + 1))
            {
                for (index_type j = right_bound + 1 - ppl; j <= right_bound; j++)
                {
                    temp[temp_size++] = Gdata.array[i][j];
                }
                // similarly
            }
            else
            {
                for (index_type j = start_index; j < start_index + ppl; j++)
                {
                    temp[temp_size++] = Gdata.array[i][j];
                }
            }
        }
        else
        {
            if (start_index != row_list_size - 1)
            {
                printf("row %d start_index %d start_value %f z: %f\n", i + 1, start_index, start_value, row_list[start_index]);
                if (start_value < -1 * alignmentAccuracy)
                {
                    start_index += 1;
                    start_value = row_list[start_index] - projectionToRow;
                    printf("row %d updated start_index %d start_value %f z: %f\n", i + 1, start_index, start_value, row_list[start_index]);
                }
            }
            // similarly adding subset of 'array' which represents row_data
            if ((start_index - ppl + 1) < left_bound)
            {
                for (index_type j = left_bound; j < left_bound + ppl; j++)
                {
                    temp[temp_size++] = Gdata.array[i][j];
                }
                // similarly
            }
            else
            {
                for (index_type j = start_index - ppl + 1; j <= start_index; j++)
                {
                    temp[temp_size++] = Gdata.array[i][j];
                }
            }
        }
        // passing in address to an uninitialized WedgeSuperPoint structure in the init_patch array with the points from temp to initialize it.
        initWedgeSuperPoint(&init_patch[init_patch_size++], temp, temp_size);
    }

    // once all points are added to patch new_patch, add the entire patch to the cover (first init it)
    wedgePatch new_patch;
    //new_patch will disappear from memory once makePatch_alignedToLine terminates, so we don't want wedgePatch_init to point superpoints to it. 
    //init_patch will also disappear for the same scope reasons
    wedgePatch_init(&new_patch, init_patch, init_patch_size, apexZ0);
    //indeed, add_patch is working fine as it is copying the values over: cover->patches[cover->n_patches] = *curr_patch;
    //doesn't matter how wedgePatch_init works since we're dereferencing the patch to store by value in an array belonging to cover.
    add_patch(cover, &new_patch);
}


void wedge_test(float apexZ0, float z0_spacing, int ppl, float z0_luminousRegion, int wedges[], int wedge_count, int lines, float top_layer_cutoff, float accept_cutoff)
{
    int numEventsLoaded = 0;

    FILE *myfile;
    myfile = fopen("cppForSynthesis/cppOutput.txt", "w"); 

    if (myfile == NULL)
    {
        printf("Error opening file");
        return;
    }

    for (index_type z = 0; z < wedges[1]; z++)
    { 
        if(z<wedges[0]) continue;
        printf("wedge %d\n", z); //main print
        fprintf(myfile, "wedge %d\n", z); //file to diff

        importData(&Gdata);
        
        addBoundaryPoint(&Gdata, 0.0001); // with default param

        wedgeCover cover;
        initWedgeCover(&cover);

        solve(&cover, apexZ0, ppl, 100, false); // solve modifies cover. false is from the left right align (previously a parameter in wedge test)

        for (int i = 0; i < (&cover)->n_patches; i++)
        {
            fprintf(myfile, "Patch \n");
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topL_jL * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topL_jR * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topR_jL * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topR_jR * 10000));

            for (int j = 0; j < cover.patches[i].superpoint_count; j++)
            {
                fprintf(myfile, "Superpoint \n");
                for (int r = 0; r < cover.patches[i].superpoints[j].point_count; r++)
                {
                    Point currentPt = cover.patches[i].superpoints[j].points[r];
                    fprintf(myfile, "%d %.4f %d %.4f\n",
                            currentPt.layer_num,
                            currentPt.phi,
                            (int)currentPt.radius,
                            currentPt.z);
                }
            }
        }
        for (int i = 0; i < cover.n_patches; i++)
        {
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(cover.patches[i].a_corner[0] * 10000),
                    lround(cover.patches[i].a_corner[1] * 10000));
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(cover.patches[i].b_corner[0] * 10000),
                    lround(cover.patches[i].b_corner[1] * 10000));
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(cover.patches[i].c_corner[0] * 10000),
                    lround(cover.patches[i].c_corner[1] * 10000));
            fprintf(myfile, "[%ld, %ld]\n",
                    lround(cover.patches[i].d_corner[0] * 10000),
                    lround(cover.patches[i].d_corner[1] * 10000));
            fprintf(myfile, "\n");
        }
        // instead of making an array of all events and passing them in, we only need access to them individually, so we will loop through and process as we create them.
    }

    fclose(myfile);
}

int main() // Not the top-level function, so you can do any FILE I/O or other non-synthesized actions here
{
    int wedgesToTest[] = {0, 1};

    wedge_test(0, 0.025, 16, 15.0, wedgesToTest, 2, 1000, 50, 15.0);

    //printf("Size of float in C: %zu bytes\n", sizeof(float));
    
    return 0;
}
