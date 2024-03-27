#include "header.h"

void wedgePatch_init(wedgePatch* wp, Environment* envI, wedgeSuperPoint* superpointsI, int superpoint_count, float apexZ0I) {
    wp->env = envI; //accessing values in wp and changing the pointer of env to point to where envI points to
    wp->end_layer = -1; //end_layer not a pointer, accessing its value.
    wp->left_end_layer = -1;
    wp->right_end_layer = -1;
    wp->left_end_lambdaZ = 0;
    wp->right_end_lambdaZ = 0;
    wp->apexZ0 = apexZ0I;

    wp->shadow_fromTopToInnermost_topL_jL = 0;
    wp->shadow_fromTopToInnermost_topL_jR = 0;
    wp->shadow_fromTopToInnermost_topR_jL = 0;
    wp->shadow_fromTopToInnermost_topR_jR = 0;

    if (superpoint_count != wp->env->num_layers) {
        printf("The patch layers does not match environment layers.");
        exit(7);
    }

    for (size_t i = 0; i < superpoint_count; i++) { //size_t objects should only be non-negative and are more performant than ints
        wp->superpoints[i] = &superpointsI[i]; //wp->superpoints is an array of pointers. Making the elements point to the elements in superpointsI.
    }
    wp->superpoint_count = superpoint_count;

    getParallelograms(wp);
    getParallelograms_v1(wp);
    get_acceptanceCorners(wp);
    get_end_layer(wp);
}

float straightLineProjectorFromLayerIJtoK(wedgePatch* wp, float z_i, float z_j, int i, int j, int k) {
    float radius_i = 0;
    float radius_j = 0;
    float radius_k = 0;

    if (i == 0) 
    {
        radius_i = 0;
    }
    else 
    {
        radius_i = wp->env->radii[i - 1]; //[] directly accessing value
    }
    if (j == 0) 
    {
        radius_j = 0;
    } 
    else 
    {
        radius_j = wp->env->radii[j - 1];
    }
    if (k == 0) 
    {
        radius_k = 0;
    } 
    else 
    {
        radius_k = wp->env->radii[k - 1];
    }

    float radii_leverArm = (radius_k - radius_i) / (radius_j - radius_i);

    return z_i + (z_j - z_i) * radii_leverArm;
}

float straightLineProjector(float z_top, float z_j, int j, Environment* env) {
    float radii_leverArm = env->radii_leverArm[j - 1];
    return z_top - (z_top - z_j) * radii_leverArm;
}

void getParallelograms(wedgePatch* wp) {
    float z1_min = max(wp->superpoints[0]->min, -1 * wp->env->trapezoid_edges[0]);
    float z1_max = min(wp->superpoints[0]->max, wp->env->trapezoid_edges[0]);

    if (z1_min > z1_max) {
        z1_min = wp->env->trapezoid_edges[0] + 1;
        z1_max = z1_min;
    }

    int previous_count = wp->parallelogram_count;

    //the code now handles the case below
    //if (! wp->parallelogram_count <= wp->superpoint_count - 1 ) {
    //    printf("Instead of assigning a temp array, we are overwriting the first superpoint_count-1 elements in the parallelogam array. If the current number of elements in the array is greater than superpoint_count-1, then we will have remaining elements that need to be deleted to replicate the functionality correctly");
    //    //exit(8);
    //}
    wp->parallelogram_count = 0; //we want to start at index 0 regardless and overwrite any old elements in the array to replicate the functionality of assigning a temp array.
    for (int i = 1; i < wp->superpoint_count; i++) {
        int j = i + 1;

        float z_j_min = wp->superpoints[i]->min;
        float z_j_max = wp->superpoints[i]->max;

        float a = straightLineProjectorFromLayerIJtoK(wp, z1_min, z_j_max, 1, j, wp->env->num_layers);
        float b = straightLineProjectorFromLayerIJtoK(wp, z1_max, z_j_max, 1, j, wp->env->num_layers);
        float c = straightLineProjectorFromLayerIJtoK(wp, z1_min, z_j_min, 1, j, wp->env->num_layers);
        float d = straightLineProjectorFromLayerIJtoK(wp, z1_max, z_j_min, 1, j, wp->env->num_layers);

        float pSlope = (j != wp->env->num_layers) ? wp->env->parallelogramSlopes[j - 1] : INT_MAX;

        //directly assign the values to the array
        if (wp->parallelogram_count < MAX_PARALLELOGRAMS_PER_PATCH) {
            Parallelogram* p = &wp->parallelograms[wp->parallelogram_count++]; //making a pointer to the address of first empty element in the array
            p->layer_num = j; //then dereferencing and assigning values to the properties
            p->pSlope = pSlope;
            p->shadow_bottomL_jR = a;
            p->shadow_bottomR_jR = b;
            p->shadow_bottomL_jL = c;
            p->shadow_bottomR_jL = d;
            p->z1_min = z1_min;
            p->z1_max = z1_max;
        }
    }
    //clear out any remaining old parallelograms if the new count is less than the old
    //filling a block of memory with a certain value: 0. applies to every field within the parallelogram structure. for pointers, this essentially means they point to null.
    for (int i = wp->parallelogram_count; i < previous_count && i < MAX_PARALLELOGRAMS_PER_PATCH; i++) {
        memset(&wp->parallelograms[i], 0, sizeof(Parallelogram));
    }
}
