#include "header.h"

void wedgePatch_init(wedgePatch *wp, Environment *envI, wedgeSuperPoint *superpointsI, int superpoint_count, float apexZ0I)
{
    wp->env = *envI;    // accessing values in wp and changing the pointer of env to point to where envI points to
    wp->end_layer = -1; // end_layer not a pointer, accessing its value.
    wp->left_end_layer = -1;
    wp->right_end_layer = -1;
    wp->left_end_lambdaZ = 0;
    wp->right_end_lambdaZ = 0;
    wp->apexZ0 = apexZ0I;

    wp->shadow_fromTopToInnermost_topL_jL = 0;
    wp->shadow_fromTopToInnermost_topL_jR = 0;
    wp->shadow_fromTopToInnermost_topR_jL = 0;
    wp->shadow_fromTopToInnermost_topR_jR = 0;

    if (superpoint_count != wp->env.num_layers)
    {
        printf("The patch layers does not match environment layers.");
        exit(7);
    }

    for (size_t i = 0; i < superpoint_count; i++)
    {                                          // size_t objects should only be non-negative and are more performant than ints
        wp->superpoints[i] = &superpointsI[i]; // wp->superpoints is an array of pointers. Making the elements point to the elements in superpointsI.
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
        radius_i = wp->env.radii[i - 1]; //[] directly accessing value
    }
    if (j == 0)
    {
        radius_j = 0;
    }
    else
    {
        radius_j = wp->env.radii[j - 1];
    }
    if (k == 0)
    {
        radius_k = 0;
    }
    else
    {
        radius_k = wp->env.radii[k - 1];
    }

    float radii_leverArm = (radius_k - radius_i) / (radius_j - radius_i);

    return z_i + (z_j - z_i) * radii_leverArm;
}

float straightLineProjector(float z_top, float z_j, int j, Environment *env)
{
    float radii_leverArm = env->radii_leverArm[j - 1];
    return z_top - (z_top - z_j) * radii_leverArm;
}

void getParallelograms(wedgePatch *wp)
{
    float z1_min = max(wp->superpoints[0]->min, -1 * wp->env.trapezoid_edges[0]);
    float z1_max = min(wp->superpoints[0]->max, wp->env.trapezoid_edges[0]);

    if (z1_min > z1_max)
    {
        z1_min = wp->env.trapezoid_edges[0] + 1;
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

        float z_j_min = wp->superpoints[i]->min;
        float z_j_max = wp->superpoints[i]->max;

        float a = straightLineProjectorFromLayerIJtoK(wp, z1_min, z_j_max, 1, j, wp->env.num_layers);
        float b = straightLineProjectorFromLayerIJtoK(wp, z1_max, z_j_max, 1, j, wp->env.num_layers);
        float c = straightLineProjectorFromLayerIJtoK(wp, z1_min, z_j_min, 1, j, wp->env.num_layers);
        float d = straightLineProjectorFromLayerIJtoK(wp, z1_max, z_j_min, 1, j, wp->env.num_layers);

        float pSlope = (j != wp->env.num_layers) ? wp->env.parallelogramSlopes[j - 1] : INT_MAX;

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
    // clear out any remaining old parallelograms if the new count is less than the old
    // filling a block of memory with a certain value: 0. applies to every field within the parallelogram structure. for pointers, this essentially means they point to null.
    for (int i = wp->parallelogram_count; i < previous_count && i < MAX_PARALLELOGRAMS_PER_PATCH; i++)
    {
        memset(&wp->parallelograms[i], 0, sizeof(Parallelogram));
    }
}

void getShadows(wedgePatch *wp, float zTopMin, float zTopMax)
{
    float zTop_min;
    float zTop_max;
    if (wp->env.num_layers - 1 < 0)
    {
        zTop_min = zTopMin;
        zTop_max = zTopMax;
    }
    else
    {
        zTop_min = max(zTopMin, -wp->env.trapezoid_edges[wp->env.num_layers - 1]);
        zTop_max = min(zTopMax, wp->env.trapezoid_edges[wp->env.num_layers - 1]);
    }

    float topL_jL[MAX_SUPERPOINTS_IN_PATCH - 1];
    float topL_jR[MAX_SUPERPOINTS_IN_PATCH - 1];
    float topR_jL[MAX_SUPERPOINTS_IN_PATCH - 1];
    float topR_jR[MAX_SUPERPOINTS_IN_PATCH - 1];

    for (int i = 0; i < wp->superpoint_count - 1; ++i)
    {
        int j = i + 1;
        float z_j_min = wp->superpoints[i]->min;
        float z_j_max = wp->superpoints[i]->max;

        topL_jL[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_min, z_j_min, wp->env.num_layers, j, 1);
        topL_jR[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_min, z_j_max, wp->env.num_layers, j, 1);
        topR_jL[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_max, z_j_min, wp->env.num_layers, j, 1);
        topR_jR[i] = straightLineProjectorFromLayerIJtoK(wp, zTop_max, z_j_max, wp->env.num_layers, j, 1);
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

/*
void getParallelograms_v1(wedgePatch* wp) {
    float top_layer_zmin = max(wp->superpoints[wp->superpoint_count - 1]->min, -1 * wp->env->top_layer_lim);
    float top_layer_zmax = min(wp->superpoints[wp->superpoint_count - 1]->max, wp->env->top_layer_lim);

    if (top_layer_zmin > top_layer_zmax) {
        top_layer_zmin = wp->env->top_layer_lim + 1;
        top_layer_zmax = top_layer_zmin;
    }

    wp->parallelogram_v1_count = 0;

    for (int i = 0; i < wp->superpoint_count - 1; i++) {
        int j = i + 1;

        float z_j_min = wp->superpoints[i]->min;
        float z_j_max = wp->superpoints[i]->max;

        float a = straightLineProjector(top_layer_zmax, z_j_min, j, wp->env); //don't need to pass in the whole wp structure
        float b = straightLineProjector(top_layer_zmax, z_j_max, j, wp->env);

        float pSlope = wp->env->parallelogramSlopes[j - 1];

        if (wp->parallelogram_v1_count < MAX_PARALLELOGRAMS_PER_PATCH) {
            Parallelogram_v1* p = &wp->parallelograms_v1[wp->parallelogram_v1_count++]; //making a pointer to the address of the first empty element in the array
            p->layer_num = j; //then dereferencing and assigning values to the properties
            p->pSlope = pSlope;
            p->shadow_topR_jL = a;
            p->shadow_topR_jR = b;
            p->top_layer_zmin = top_layer_zmin;
            p->top_layer_zmax = top_layer_zmax;
        }
    }
}
*/

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
    if (a_corner_min != wp->parallelograms[wp->env.num_layers - 2].shadow_bottomL_jR)
    {
        wp->squareAcceptance = false;
        wp->flatTop = false;
    }
    if (b_corner_min != wp->parallelograms[wp->env.num_layers - 2].shadow_bottomR_jR)
    {
        wp->squareAcceptance = false;
        wp->flatTop = false;
    }
    if (c_corner_max != wp->parallelograms[wp->env.num_layers - 2].shadow_bottomL_jL)
    {
        wp->squareAcceptance = false;
        wp->flatBottom = false;
    }
    if (d_corner_max != wp->parallelograms[wp->env.num_layers - 2].shadow_bottomR_jL)
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
    for (int i = 0; i < wp->env.num_layers; i++)
    {
        float lambdaZ_left = (wp->superpoints[i]->min - wp->apexZ0) / wp->env.radii[i];
        float lambdaZ_right = (wp->superpoints[i]->max - wp->apexZ0) / wp->env.radii[i];

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
