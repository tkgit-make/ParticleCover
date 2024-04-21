#include "header.h"

void initWedgeCover(wedgeCover *wc, Environment *envI, DataSet *dataI)
{
    wc->n_patches = 0;
    wc->env = envI;
    wc->data = dataI;
    /*
    for (index_type i = 0; i < MAX_PATCHES; i++) {
        wc -> all_patches[i] = NULL;
        wc -> real_patch_list[i] = false;
    }
    for (index_type i = 0; i < MAX_SUPERPOINTS_IN_COVER; i++) {
        wc -> superPoints[i] = NULL;
    }
    */
}

void add_patch(wedgeCover *cover, wedgePatch *curr_patch)
{
    if (cover->n_patches == 0)
    {
        cover->patches[0] = *curr_patch; //copy the patch directly
        //cover->all_patches[0] = curr_patch;
        cover->real_patch_list[0] = true;
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
                cover->real_patch_list[cover->n_patches] = true;
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

    cover->real_patch_list[index] = false;

    for (index_type i = index; i < cover->n_patches - 1; i++)
    {
        cover->patches[i] = cover->patches[i + 1];
        cover->real_patch_list[i] = cover->real_patch_list[i + 1];
    }

    // resetting the last elements
    memset(&cover->patches[cover->n_patches - 1], 0, sizeof(wedgePatch));
    cover->real_patch_list[cover->n_patches - 1] = false;

    cover->n_patches -= 1;
}

// can't provide default parameters
index_type get_index_from_z(DataSet *data, int layer, float z_value, int alignment)
{
    // c doesn't support string comparison directly, using integer comparison for effiency
    // CLOSEST = 11, ABOVE = 12, BELOW = 13
    float minVal = 1000000;
    index_type index = 0;

    for (index_type i = 0; i < data->n_points[layer]; i++)
    {
        float diff = fabs(data->array[layer][i].z - z_value); // absolute difference
        if (diff < minVal)
        {
            minVal = diff;
            index = i;
        }
    }

    if (alignment == CLOSEST)
    {
        return index;
    }

    if (alignment == ABOVE)
    {
        if (data->array[layer][index].z > z_value)
        {
            return index;
        }
        // potential bounds issue
        return index + 1;
    }

    if (alignment == BELOW)
    {
        if (data->array[layer][index].z < z_value)
        {
            return index;
        }
        return index - 1;
    }

    return index;
}

// not implementing the logic corresponding to show = true (that would involve Line Generators, etc)
// Line Generator and its accompanying methods have been coded, but we are not going to implement show=true case here as main method passes with show=false.
// MAKE_PATCHES_SHADOW_QUILT_FROM_EDGES defined as 33
void solve(wedgeCover *cover, int lining, float apexZ0, int ppl, int nlines, bool leftRight)
{
    for (index_type i = 0; i < cover->env->num_layers; i++)
    {
        bool foundIdentical = false;
        bool firstTime = true;

        while (foundIdentical || firstTime)
        {
            foundIdentical = false;
            for (index_type x = 0; x < cover->data->n_points[i] - 1; x++)
            {
                if (cover->data->array[i][x].z == cover->data->array[i][x + 1].z)
                {
                    cover->data->array[i][x + 1].z += 0.00001;
                    foundIdentical = true;
                }
            }

            firstTime = false;
            if (foundIdentical)
            {
                qsort(cover->data->array[i], cover->data->n_points[i], sizeof(Point), comparePoints);
            }
        }
    }

    if (lining == MAKE_PATCHES_SHADOW_QUILT_FROM_EDGES)
    {
        makePatches_ShadowQuilt_fromEdges(cover, apexZ0, 1, ppl, leftRight);
    }
}

void makePatches_ShadowQuilt_fromEdges(wedgeCover *cover, float apexZ0, int stop, int ppl, bool leftRight)
{
    bool fix42 = true;
    apexZ0 = cover->env->trapezoid_edges[0];
    float saved_apexZ0;

    while (apexZ0 > -1 * cover->env->trapezoid_edges[0])
    {
        float z_top_min = -1 * cover->env->top_layer_lim;

        float complementary_apexZ0 = 0;
        index_type first_row_count = 0;
        float c_corner = LONG_MAX;

        float z_top_max = cover->env->top_layer_lim + cover->env->boundaryPoint_offset;

        if (cover->n_patches > 0)
        {
            z_top_max = min(z_top_max, straightLineProjectorFromLayerIJtoK(&cover->patches[cover->n_patches - 1], -1 * cover->env->beam_axis_lim, apexZ0, 0, 1, cover->env->num_layers // includes passing a pointer to the last patch
                                                                           ));
        }

        index_type nPatchesInColumn = 0;
        float projectionOfCornerToBeam = 0;

        while ((c_corner > -1 * cover->env->trapezoid_edges[cover->env->num_layers - 1]) && (projectionOfCornerToBeam < cover->env->beam_axis_lim))
        {
            nPatchesInColumn++;
            printf("%f %d %f %d\n", apexZ0, ppl, z_top_max, leftRight);
            // important assumption: assuming 'makePatch_alignedToLine' updates the 'patches' array in 'cover'
            // will need to revisit parameters when we write this method
            makePatch_alignedToLine(cover, apexZ0, z_top_max, ppl, false, false);

            index_type lastPatchIndex = cover->n_patches - 1;

            printf("top layer from %f to %f z_top_max: %f\n",
                   cover->patches[lastPatchIndex].superpoints[cover->env->num_layers - 1].max,
                   cover->patches[lastPatchIndex].superpoints[cover->env->num_layers - 1].min,
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
                       straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].max, cover->patches[lastPatchIndex].superpoints[i].min, 1, j, cover->env->num_layers),
                       straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].max, cover->patches[lastPatchIndex].superpoints[i].max, 1, j, cover->env->num_layers),
                       straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].min, cover->patches[lastPatchIndex].superpoints[i].min, 1, j, cover->env->num_layers),
                       straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], cover->patches[lastPatchIndex].superpoints[0].min, cover->patches[lastPatchIndex].superpoints[i].max, 1, j, cover->env->num_layers));
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
                repeat_original = (cover->patches[lastPatchIndex].superpoints[cover->env->num_layers - 1] == cover->patches[thirdLastPatchIndex].superpoints[cover->env->num_layers - 1]) &&
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
            projectionOfCornerToBeam = straightLineProjectorFromLayerIJtoK(lastPatch, lastPatch->c_corner[1], lastPatch->c_corner[0], cover->env->num_layers, 1, 0);
            bool squarePatch_alternate1 = (lastPatch->a_corner[1] > z_top_max) && (lastPatch->b_corner[1] > z_top_max) && (lastPatch->flatBottom);
            bool squarePatch_alternate2 = (lastPatch->a_corner[1] > z_top_max) && (lastPatch->flatBottom);

            bool notChoppedPatch = (lastPatch->squareAcceptance) || squarePatch_alternate1 || squarePatch_alternate2;
            bool madeComplementaryPatch = false;

            int nPatchesAtOriginal = cover->n_patches;

            printf("squareAcceptance: %d triangleAcceptance: %d projectionOfCornerToBeam: %f notChoppedPatch %d\n",
                   lastPatch->squareAcceptance, lastPatch->triangleAcceptance, projectionOfCornerToBeam, notChoppedPatch);

            if (!notChoppedPatch && (lastPatch->c_corner[1] > -1 * cover->env->trapezoid_edges[cover->env->num_layers - 1]) && (projectionOfCornerToBeam < cover->env->beam_axis_lim))
            {
                complementary_apexZ0 = lastPatch->superpoints[0].min;
                if (lastPatch->triangleAcceptance && !repeat_original)
                {
                    z_top_min = lastPatch->d_corner[1];
                }
                else
                {
                    printf("z_top_min before: %f superpoints[self.env.num_layers-1].min: %f\n", z_top_min, lastPatch->superpoints[cover->env->num_layers - 1].min);
                    z_top_min = max(-1 * cover->env->top_layer_lim, lastPatch->superpoints[cover->env->num_layers - 1].min);
                }
                // will need to revisit parameters when we write this method
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
                double previous_z_top_min = -999;

                while (!(white_space_height <= 0 && (previous_white_space_height >= 0)) && (fabs(white_space_height) > 0.000001) &&
                       ((cover->patches[lastPatchIndex].c_corner[1] > -1 * cover->env->trapezoid_edges[cover->env->num_layers - 1]) ||
                        (white_space_height > 0)) &&
                       (current_z_top_index < (int)(cover->data->n_points[cover->env->num_layers - 1]) &&
                        !(repeat_patch) && !(repeat_original)))
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

                    current_z_top_index = get_index_from_z(cover->data, cover->env->num_layers - 1, z_top_min, CLOSEST); // CLOSEST was default param in C++
                    printf("current white_space_height: %f\n", white_space_height);
                    printf("counter: %d counterUpshift: %d\n", counter, counterUpshift);
                    printf("orig_ztop: %d orig_z_top_min: %f\n", current_z_top_index, z_top_min);

                    index_type current_z_i_index[MAX_LAYERS];
                    index_type new_z_i_index[MAX_LAYERS];
                    // initalizing with 0s, revisit if this is needed
                    // memset(current_z_i_index, 0, sizeof(current_z_i_index));
                    // memset(new_z_i_index, 0, sizeof(new_z_i_index));

                    for (index_type i = 0; i < cover->env->num_layers; i++)
                    {
                        current_z_i_index[i] = get_index_from_z(cover->data, i, straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_apexZ0, z_top_min, 1, cover->env->num_layers, i + 1), CLOSEST);
                    }

                    if (z_top_min == previous_z_top_min)
                    {
                        current_z_top_index += 1;
                        for (index_type i = 0; i < cover->env->num_layers; i++)
                        {
                            new_z_i_index[i] = current_z_i_index[i] + 1;
                        }
                    }

                    previous_z_top_min = z_top_min;

                    if (white_space_height < 0)
                    {
                        counter += 1;
                        current_z_top_index -= 1;
                        for (index_type i = 0; i < cover->env->num_layers; i++)
                        {
                            new_z_i_index[i] = current_z_i_index[i] - 1;
                        }
                    }
                    else
                    {
                        counterUpshift += 1;
                        current_z_top_index += 1;
                        for (index_type i = 0; i < cover->env->num_layers; i++)
                        {
                            new_z_i_index[i] = current_z_i_index[i] + 1;
                        }
                    }

                    int x = cover->data->n_points[cover->env->num_layers - 1] - 1;
                    current_z_top_index = min(current_z_top_index, cover->data->n_points[cover->env->num_layers - 1] - 1); // n_points is an array of the sizes of each element of 'array'

                    for (index_type i = 0; i < cover->env->num_layers; i++)
                    {
                        new_z_i_index[i] = min(new_z_i_index[i], (float)cover->data->n_points[i] - 1);
                    }

                    for (index_type i = 0; i < cover->env->num_layers; i++)
                    { // replaced loop variables with macro index_type. made index_type unsigned int.
                        new_z_i_index[i] = max(new_z_i_index[i], 0.0f);
                    }
                    float new_z_i[MAX_LAYERS]; // unsigned integer.

                    for (index_type i = 0; i < cover->env->num_layers; i++)
                    {
                        new_z_i[i] = cover->data->array[i][new_z_i_index[i]].z;
                    }

                    float new_z_i_atTop[MAX_LAYERS - 1]; // note: the size is MAX_LAYERS - 1 because the loop starts from 1
                    for (index_type i = 1; i < cover->env->num_layers; i++)
                    {
                        new_z_i_atTop[i - 1] = straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex],
                                                                                   complementary_apexZ0,
                                                                                   new_z_i[i],
                                                                                   1,
                                                                                   i + 1,
                                                                                   cover->env->num_layers);
                    }

                    index_type layerWithSmallestShift = 0;
                    float layerSMin = FLT_MAX;

                    for (index_type i = 0; i < cover->env->num_layers - 1; i++)
                    {
                        if (fabs(new_z_i_atTop[i] - previous_z_top_min) < layerSMin)
                        { // fabs is for floats. abs is only int
                            layerSMin = fabs(new_z_i_atTop[i] - previous_z_top_min);
                            layerWithSmallestShift = i;
                        }
                    }

                    layerWithSmallestShift += 1;

                    for (index_type i = 0; i < cover->env->num_layers - 1; i++)
                    {
                        printf("%u new_z_i_atTop: %f shift_i_ztop: %f layerWithSmallestShift: %u\n",
                               i + 1, new_z_i_atTop[i], new_z_i_atTop[i] - previous_z_top_min, layerWithSmallestShift + 1);
                    }

                    z_top_min = cover->data->array[cover->env->num_layers - 1][current_z_top_index].z;
                    z_top_min = new_z_i_atTop[layerWithSmallestShift - 1];

                    if (fabs(z_top_min - previous_z_top_min) < 0.000001)
                    {
                        z_top_min = cover->data->array[cover->env->num_layers - 1][current_z_top_index].z;
                    }

                    if (fabs(z_top_min - previous_z_top_min) < 0.000001)
                    {
                        z_top_min = cover->data->array[cover->env->num_layers - 2][current_z_top_index].z;
                    }

                    if (fabs(z_top_min - previous_z_top_min) < 0.000001)
                    {
                        z_top_min = cover->data->array[cover->env->num_layers - 3][current_z_top_index].z;
                    }

                    if (((z_top_min - previous_z_top_min) * white_space_height) < 0)
                    {
                        z_top_min = new_z_i_atTop[cover->env->num_layers - 2];
                    }

                    printf(" new_def_z_top_min_diff: %f\n", z_top_min - cover->data->array[cover->env->num_layers - 1][current_z_top_index].z);

                    printf(" new_ztop_index: %d new_z_i_index: ", current_z_top_index);
                    for (index_type i = 0; i < cover->env->num_layers; i++)
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
                           complementary_b, cover->patches[lastPatchIndex].b_corner[1], z_top_min);
                    printf("new white_space_height: %f\n", white_space_height);
                    printf("adjusted complementary: [%f, %f] for z_top_min: %f\n",
                           cover->patches[lastPatchIndex].a_corner[0], cover->patches[lastPatchIndex].a_corner[1], z_top_min);
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
                        bool repeat_patch = true;
                        // turned this into a for loop, dynamic. if ((patches[patches.size() - 1].superpoints[env.num_layers - 1] == patches[patches.size() - 3].superpoints[env.num_layers - 1]) && (patches[patches.size() - 1].superpoints[0] == patches[patches.size() - 3].superpoints[0]) && (patches[patches.size() - 1].superpoints[1] == patches[patches.size() - 3].superpoints[1]) && (patches[patches.size() - 1].superpoints[2] == patches[patches.size() - 3].superpoints[2]) && (patches[patches.size() - 1].superpoints[3] == patches[patches.size() - 3].superpoints[3]))
                        // that code checked 0 to 4
                        for (index_type i = 0; i < cover->env->num_layers; i++)
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
                                   cover->patches[lastPatchIdx].superpoints[cover->env->num_layers - 1].min,
                                   cover->patches[lastPatchIdx].superpoints[cover->env->num_layers - 1].max,
                                   repeat_patch);

                            delete_patch(cover, lastPatchIdx);

                            current_z_top_index -= 1;

                            z_top_min = cover->data->array[cover->env->num_layers - 1][current_z_top_index].z;
                            z_top_min = new_z_i_atTop[layerWithSmallestShift - 1];

                            makePatch_alignedToLine(cover, complementary_apexZ0, z_top_min, ppl, true, false);
                        }
                    }
                }
            }

            lastPatchIndex = cover->n_patches - 1; // just to keep fresh in case we use it
            c_corner = cover->patches[lastPatchIndex].c_corner[1];

            projectionOfCornerToBeam = straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], c_corner, cover->patches[lastPatchIndex].c_corner[0], cover->env->num_layers, 1, 0);

            saved_apexZ0 = cover->patches[lastPatchIndex].c_corner[0];

            if (madeComplementaryPatch)
            {
                int secondLastPatchIndex = lastPatchIndex - 1;

                // modifying patches, not adding patches, so index variables do not need to be updated.
                getShadows(&cover->patches[lastPatchIndex], z_top_min, z_top_max);
                getShadows(&cover->patches[secondLastPatchIndex], z_top_min, z_top_max);

                float original_topR_jL = cover->patches[secondLastPatchIndex].shadow_fromTopToInnermost_topR_jL;
                bool originalPartialTop = (original_topR_jL > complementary_apexZ0) && (original_topR_jL < apexZ0) &&
                                          (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[secondLastPatchIndex], original_topR_jL, z_top_max, 1, cover->env->num_layers, 0)) < 20 * cover->env->beam_axis_lim);

                float original_topL_jL = cover->patches[secondLastPatchIndex].shadow_fromTopToInnermost_topL_jL;
                bool originalPartialBottom = (original_topL_jL > complementary_apexZ0) && (original_topL_jL < apexZ0) &&
                                             (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[secondLastPatchIndex], original_topL_jL, z_top_min, 1, cover->env->num_layers, 0)) < 20 * cover->env->beam_axis_lim);

                float complementary_topR_jR = cover->patches[lastPatchIndex].shadow_fromTopToInnermost_topR_jR;
                bool complementaryPartialTop = (complementary_topR_jR > complementary_apexZ0) && (complementary_topR_jR < apexZ0) &&
                                               (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_topR_jR, z_top_max, 1, cover->env->num_layers, 0)) < 20 * cover->env->beam_axis_lim);

                float complementary_topL_jR = cover->patches[lastPatchIndex].shadow_fromTopToInnermost_topL_jR;
                bool complementaryPartialBottom = (complementary_topL_jR > complementary_apexZ0) && (complementary_topL_jR < apexZ0) &&
                                                  (fabs(straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_topL_jR, z_top_min, 1, cover->env->num_layers, 0)) < 20 * cover->env->beam_axis_lim);

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

                float z0_original_bCorner = straightLineProjectorFromLayerIJtoK(&cover->patches[secondLastPatchIndex], apexZ0, z_top_max, 1, cover->env->num_layers, 0);
                float z0_complementary_cCorner = straightLineProjectorFromLayerIJtoK(&cover->patches[lastPatchIndex], complementary_apexZ0, z_top_min, 1, cover->env->num_layers, 0);
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

                if (horizontalShiftTop > 0 || horizontalShiftBottom > 0)
                {
                    printf("originalPartialTop: %d complementaryPartialTop: %d originalPartialBottom: %d complementaryPartialBottom: %d %f %f %f %f horizontalOverlapTop: %f horizontalOverlapBottom: %f\n",
                           originalPartialTop, complementaryPartialTop, originalPartialBottom, complementaryPartialBottom,
                           original_topR_jL, original_topL_jL, complementary_topR_jR, complementary_topL_jR,
                           horizontalOverlapTop, horizontalOverlapBottom);
                }

                while (((horizontalShiftTop > 0 && originalPartialTop && complementaryPartialTop) || (horizontalShiftBottom > 0 && originalPartialBottom && complementaryPartialBottom)) && doShiftedPatch && (horizontalOverlapTop <= 0) && (horizontalOverlapBottom <= 0) && (newGapTop < 0 || newGapBottom < 0))
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

                    if (shiftOriginal && straightLineProjectorFromLayerIJtoK(&cover->patches[cover->n_patches - 1], original_topR_jR, z_top_max, 1, cover->env->num_layers, 0) < cover->env->beam_axis_lim)
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
                    if ((straightLineProjectorFromLayerIJtoK(&cover->patches[cover->n_patches - 1], shifted_Align, newZtop, 1, cover->env->num_layers, 0) > cover->env->beam_axis_lim) && shiftOriginal)
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

        apexZ0 = cover->patches[cover->n_patches - 1].c_corner[0];
        apexZ0 = saved_apexZ0;
        printf("'=======================================================  z1_Align: %f\n", apexZ0);
    }
}

void makePatch_alignedToLine(wedgeCover *cover, float apexZ0, float z_top, int ppl, bool leftRight, bool float_middleLayers_ppl)
{
    wedgeSuperPoint init_patch[MAX_LAYERS]; // correct
    int original_ppl = ppl;
    float alignmentAccuracy = 0.00001;
    // Point row_data[MAX_LAYERS][MAX_POINTS_FOR_DATASET];
    index_type init_patch_size = 0;

    for (index_type i = 0; i < cover->env->num_layers; i++)
    {
        float y = cover->env->radii[i];
        float row_list[MAX_POINTS_PER_LAYER];
        int row_list_size = 0;

        for (index_type j = 0; j < cover->data->n_points[i]; j++)
        {
            row_list[row_list_size++] = cover->data->array[i][j].z;
        }

        float r_max = cover->env->radii[cover->env->num_layers - 1];
        float projectionToRow = (z_top - apexZ0) * (y - cover->env->radii[0]) / (r_max - cover->env->radii[0]) + apexZ0;

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
            if (fabs((row_list[j] + cover->env->trapezoid_edges[i] + cover->env->boundaryPoint_offset)) < lbVal)
            {
                left_bound = j;
                lbVal = fabs((row_list[j] + cover->env->trapezoid_edges[i] + cover->env->boundaryPoint_offset));
            }

            if (fabs((row_list[j] - cover->env->trapezoid_edges[i] - cover->env->boundaryPoint_offset)) < rbVal)
            {
                right_bound = j;
                rbVal = fabs((row_list[j] - cover->env->trapezoid_edges[i] - cover->env->boundaryPoint_offset));
            }
        }

        if (float_middleLayers_ppl && i != 0 && i != cover->env->num_layers - 1)
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
                    temp[temp_size++] = cover->data->array[i][j];
                }
                // similarly
            }
            else
            {
                for (index_type j = start_index; j < start_index + ppl; j++)
                {
                    temp[temp_size++] = cover->data->array[i][j];
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
                    temp[temp_size++] = cover->data->array[i][j];
                }
                // similarly
            }
            else
            {
                for (index_type j = start_index - ppl + 1; j <= start_index; j++)
                {
                    temp[temp_size++] = cover->data->array[i][j];
                }
            }
        }
        // passing in address to an uninitialized WedgeSuperPoint structure in the init_patch array with the points from temp to initialize it.
        initWedgeSuperPoint(&init_patch[init_patch_size++], temp, temp_size);
    }

    // once all points are added to init_patch, add the entire patch to the cover (first init it)
    wedgePatch new_patch;
    //init_patch will disappear from memory once makePatch_alignedToLine terminates, so we don't want wedgePatch_init to point superpoints to it. 
    wedgePatch_init(&new_patch, cover->env, init_patch, init_patch_size, apexZ0);
    add_patch(cover, &new_patch);
}
