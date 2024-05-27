#include "header.h"

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


