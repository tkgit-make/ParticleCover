#include "header.h"

void initWedgeCover()
{
    n_patches = 0;
    //wc->data = dataI;
}

void add_patch(wedgePatch *curr_patch)
{
    if (n_patches == 0)
    {
        patches[0] = *curr_patch; //copy the patch directly
        //all_patches[0] = curr_patch;
        n_patches = 1;
    }
    else
    {
        wedgePatch *prev_patch = &(patches[n_patches - 1]);
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
            if (n_patches < MAX_PATCHES)
            {
                patches[n_patches] = *curr_patch;
                //all_patches[n_patches] = curr_patch;
                n_patches += 1;
            }
        }
    }
}

void delete_patch(int index)
{
    if (index < 0 || index >= n_patches)
    {
        return;
    }

    for (index_type i = index; i < n_patches - 1; i++)
    {
        patches[i] = patches[i + 1];
    }

    // resetting the last elements
    memset(&patches[n_patches - 1], 0, sizeof(wedgePatch));

    n_patches -= 1;
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
void solve(float apexZ0, int ppl, int nlines, bool leftRight)
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
                qsort(Gdata.array[i], Gdata.n_points[i], sizeof(Point), comparePoints);
            }
        }
    }
    makePatches_ShadowQuilt_fromEdges(apexZ0, 1, ppl, leftRight);
}

void makePatches_ShadowQuilt_fromEdges(float apexZ0, int stop, int ppl, bool leftRight)
{
    bool fix42 = true;
    apexZ0 = trapezoid_edges[0];
    float saved_apexZ0;

    while (apexZ0 > -1 * trapezoid_edges[0]) //consider how this works when we are expanding instead of retracting the trapezoid_edges
    {
        float z_top_min = -1 * top_layer_lim;

        float complementary_apexZ0 = 0;
        index_type first_row_count = 0;
        float c_corner = LONG_MAX;

        float z_top_max = top_layer_lim;

        if (n_patches > 0)
        {
            z_top_max = min(z_top_max, straightLineProjectorFromLayerIJtoK(&patches[n_patches - 1], -1 * beam_axis_lim, apexZ0, 0, 1, num_layers // includes passing a pointer to the last patch
                            ));
        }

        index_type nPatchesInColumn = 0;
        float projectionOfCornerToBeam = 0;
        
        //remove nPatchesInColumn once debugging finishes
        while((c_corner > -1 * trapezoid_edges[num_layers - 1]) && (nPatchesInColumn<100000000) && (projectionOfCornerToBeam < beam_axis_lim))
        {
            nPatchesInColumn++;
            printf("%f %d %f %d\n", apexZ0, ppl, z_top_max, leftRight);

            makePatch_alignedToLine(apexZ0, z_top_max, ppl, false, false);

            index_type lastPatchIndex = n_patches - 1;

            printf("top layer from %f to %f z_top_max: %f\n",
                   patches[lastPatchIndex].superpoints[num_layers - 1].max,
                   patches[lastPatchIndex].superpoints[num_layers - 1].min,
                   z_top_max);
            printf("original: [%f, %f] for patch %d\n",
                   patches[lastPatchIndex].a_corner[0],
                   patches[lastPatchIndex].a_corner[1],
                   n_patches);
            printf("original: [%f, %f]\n",
                   patches[lastPatchIndex].b_corner[0],
                   patches[lastPatchIndex].b_corner[1]);

            printf("original: [%f, %f]\n",
                   patches[lastPatchIndex].c_corner[0],
                   patches[lastPatchIndex].c_corner[1]);

            printf("original: [%f, %f]\n",
                   patches[lastPatchIndex].d_corner[0],
                   patches[lastPatchIndex].d_corner[1]);

            for (index_type i = 1; i < patches[lastPatchIndex].superpoint_count - 1; i++)
            {
                index_type j = i + 1;
                printf("%d superpoint: %f %f shadowTop from L1Max: %f %f from L1 min: %f %f\n",
                       j,
                       patches[lastPatchIndex].superpoints[i].min,
                       patches[lastPatchIndex].superpoints[i].max,
                       straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], patches[lastPatchIndex].superpoints[0].max, patches[lastPatchIndex].superpoints[i].min, 1, j, num_layers),
                       straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], patches[lastPatchIndex].superpoints[0].max, patches[lastPatchIndex].superpoints[i].max, 1, j, num_layers),
                       straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], patches[lastPatchIndex].superpoints[0].min, patches[lastPatchIndex].superpoints[i].min, 1, j, num_layers),
                       straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], patches[lastPatchIndex].superpoints[0].min, patches[lastPatchIndex].superpoints[i].max, 1, j, num_layers));
            }

            float original_c = patches[lastPatchIndex].c_corner[1];
            float original_d = patches[lastPatchIndex].d_corner[1];

            c_corner = original_c;

            bool repeat_patch = false;
            bool repeat_original = false;

            // code written assuming number of layers is 5.
            /*
            if (n_patches > 2) {
                int thirdLastPatchIndex = lastPatchIndex - 2;
                repeat_original = (patches[lastPatchIndex].superpoints[num_layers - 1] == patches[thirdLastPatchIndex].superpoints[num_layers - 1]) &&
                        (patches[lastPatchIndex].superpoints[0] == patches[thirdLastPatchIndex].superpoints[0]) &&
                        (patches[lastPatchIndex].superpoints[1] == patches[thirdLastPatchIndex].superpoints[1]) &&
                        (patches[lastPatchIndex].superpoints[2] == patches[thirdLastPatchIndex].superpoints[2]) &&
                        (patches[lastPatchIndex].superpoints[3] == patches[thirdLastPatchIndex].superpoints[3]);
            }
            */
            // dynamic version below
            if (n_patches > 2)
            {
                index_type thirdLastPatchIndex = lastPatchIndex - 2;
                repeat_original = true; // assume they are the same initially
                for (index_type i = 0; i < MAX_SUPERPOINTS_IN_PATCH; i++)
                { // iterating over the first (five) superpoints
                    if (!areWedgeSuperPointsEqual(&patches[lastPatchIndex].superpoints[i], &patches[thirdLastPatchIndex].superpoints[i]))
                    {
                        repeat_original = false; // if any pair of superpoints don't match, set to false
                        break;                   // no need to check further if a mismatch is found
                    }
                }
            }

            float seed_apexZ0 = apexZ0;
            wedgePatch *lastPatch = &patches[lastPatchIndex];
            projectionOfCornerToBeam = straightLineProjectorFromLayerIJtoK(lastPatch, lastPatch->c_corner[1], lastPatch->c_corner[0], num_layers, 1, 0);
            bool squarePatch_alternate1 = (lastPatch->a_corner[1] > z_top_max) && (lastPatch->b_corner[1] > z_top_max) && (lastPatch->flatBottom);
            bool squarePatch_alternate2 = (lastPatch->a_corner[1] > z_top_max) && (lastPatch->flatBottom);

            bool notChoppedPatch = (lastPatch->squareAcceptance) || squarePatch_alternate1 || squarePatch_alternate2;
            bool madeComplementaryPatch = false;

            int nPatchesAtOriginal = n_patches;

            printf("squareAcceptance: %d triangleAcceptance: %d projectionOfCornerToBeam: %f notChoppedPatch %d\n",
                   lastPatch->squareAcceptance, lastPatch->triangleAcceptance, projectionOfCornerToBeam, notChoppedPatch);

            if (!notChoppedPatch && (lastPatch->c_corner[1] > -1 * trapezoid_edges[num_layers - 1]) && (projectionOfCornerToBeam < beam_axis_lim))
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

                makePatch_alignedToLine(complementary_apexZ0, z_top_min, ppl, true, false);
                // updating the last patch index because makePatch_alignedToLine will add more patches to the patches array. Should revisit after writing method
                // makePatch_alignedToLine will call the add patch method, so we must get a new last patch index.
                lastPatchIndex = n_patches - 1;

                madeComplementaryPatch = true;
                printf("complementary: [%f, %f] for z_top_min: %f\n", patches[lastPatchIndex].a_corner[0], patches[lastPatchIndex].a_corner[1], z_top_min);
                printf("complementary: [%f, %f] for patch %d\n", patches[lastPatchIndex].b_corner[0], patches[lastPatchIndex].b_corner[1], n_patches);
                printf("complementary: [%f, %f]\n", patches[lastPatchIndex].c_corner[0], patches[lastPatchIndex].c_corner[1]);
                printf("complementary: [%f, %f]\n", patches[lastPatchIndex].d_corner[0], patches[lastPatchIndex].d_corner[1]);

                float complementary_a = patches[lastPatchIndex].a_corner[1];
                float complementary_b = patches[lastPatchIndex].b_corner[1];

                float white_space_height = max(original_c - complementary_a, original_d - complementary_b);
                float previous_white_space_height = -1;
                int counter = 0;
                int counterUpshift = 0;
                index_type current_z_top_index = -1;
                double previous_z_top_min = -999;

                while (!(white_space_height <= 0.0000005 && (previous_white_space_height >= 0)) && (fabs(white_space_height) > 0.000005) &&
                       ((patches[lastPatchIndex].c_corner[1] > -1 * trapezoid_edges[num_layers - 1]) ||
                        (white_space_height > 0.000005)) &&
                       (current_z_top_index < (int)(Gdata.n_points[num_layers - 1])) &&
                        !(repeat_patch) && !(repeat_original))
                {
                    printf("\n");
                    if (n_patches > 2)
                    {
                        index_type secondLastPatchIndex = lastPatchIndex - 1;
                        printf("original c: %f %f || original d: %f %f\n",
                               original_c, patches[secondLastPatchIndex].c_corner[1],
                               original_d, patches[secondLastPatchIndex].d_corner[1]);
                    }
                    printf("complementary_a: %f %f || complementary_b: %f %f\n",
                           complementary_a, patches[lastPatchIndex].a_corner[1],
                           complementary_b, patches[lastPatchIndex].b_corner[1]);

                    current_z_top_index = get_index_from_z(num_layers - 1, z_top_min); 
                    printf("current white_space_height: %f\n", white_space_height);
                    printf("counter: %d counterUpshift: %d\n", counter, counterUpshift);
                    printf("orig_ztop: %d orig_z_top_min: %f\n", current_z_top_index, z_top_min);

                    index_type current_z_i_index[MAX_LAYERS];
                    index_type new_z_i_index[MAX_LAYERS];

                    for (index_type i = 0; i < num_layers; i++)
                    {
                        current_z_i_index[i] = get_index_from_z(i, straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], complementary_apexZ0, z_top_min, 1, num_layers, i + 1));
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
                        new_z_i_atTop[i - 1] = straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex],
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

                    int nPatchesAtComplementary = n_patches;
                    lastPatchIndex = n_patches - 1; // this may have already been updated at the end of the last call, but just to be sure
                    if (nPatchesAtComplementary > nPatchesAtOriginal)
                    {
                        printf("deleted complementary: [%f, %f] for patch %d\n",
                               patches[lastPatchIndex].a_corner[0],
                               patches[lastPatchIndex].a_corner[1],
                               n_patches);
                        printf("deleted complementary: [%f, %f]\n",
                               patches[lastPatchIndex].b_corner[0],
                               patches[lastPatchIndex].b_corner[1]);
                        printf("deleted complementary: [%f, %f]\n",
                               patches[lastPatchIndex].c_corner[0],
                               patches[lastPatchIndex].c_corner[1]);
                        printf("deleted complementary: [%f, %f]\n",
                               patches[lastPatchIndex].d_corner[0],
                               patches[lastPatchIndex].d_corner[1]);

                        delete_patch(lastPatchIndex);
                        // no need to manually decrement n_patches, delete_patch will handle it
                    }
                    lastPatchIndex = n_patches - 1; // lastPatchIndex has changed because of the delete patch
                    // it may be not needed to update lastPatchIndex, but for now, I did it, so it wouldn't be forgotten later.

                    // call makePatch_alignedToLine to add a new patch based on the complementary apex and top z values.
                    makePatch_alignedToLine(complementary_apexZ0, z_top_min, ppl, true, false);
                    // update the lastPatchIndex to point to the newly added patch.
                    lastPatchIndex = n_patches - 1;

                    // retrieve the a and b corner values from the latest patch.
                    complementary_a = patches[lastPatchIndex].a_corner[1];
                    complementary_b = patches[lastPatchIndex].b_corner[1];

                    // update the previous white space height for the next iteration.
                    previous_white_space_height = white_space_height;
                    // calculate the new white space height based on the original and complementary corners.
                    white_space_height = max(original_c - complementary_a, original_d - complementary_b);

                    printf("complementary_a: %f %f || complementary_b: %f %f new z_top_min: %f\n",
                           complementary_a, patches[lastPatchIndex].a_corner[1],
                           complementary_b, patches[lastPatchIndex].b_corner[1], z_top_min);
                    printf("new white_space_height: %f\n", white_space_height);
                    printf("adjusted complementary: [%f, %f] for z_top_min: %f\n",
                           patches[lastPatchIndex].a_corner[0], patches[lastPatchIndex].a_corner[1], z_top_min);
                    printf("adjusted complementary: [%f, %f] for patch %d\n",
                           patches[lastPatchIndex].b_corner[0], patches[lastPatchIndex].b_corner[1], n_patches);
                    printf("adjusted complementary: [%f, %f]\n",
                           patches[lastPatchIndex].c_corner[0], patches[lastPatchIndex].c_corner[1]);
                    printf("adjusted complementary: [%f, %f]\n",
                           patches[lastPatchIndex].d_corner[0], patches[lastPatchIndex].d_corner[1]);

                    if ((n_patches > 3) && fix42)
                    {
                        index_type lastPatchIdx = n_patches - 1;
                        index_type thirdLastPatchIdx = lastPatchIdx - 2;

                        // checking if the superpoints of the last and third last patches are the same
                        bool repeat_patch = true;
                        // turned this into a for loop, dynamic. if ((patches[patches.size() - 1].superpoints[env.num_layers - 1] == patches[patches.size() - 3].superpoints[env.num_layers - 1]) && (patches[patches.size() - 1].superpoints[0] == patches[patches.size() - 3].superpoints[0]) && (patches[patches.size() - 1].superpoints[1] == patches[patches.size() - 3].superpoints[1]) && (patches[patches.size() - 1].superpoints[2] == patches[patches.size() - 3].superpoints[2]) && (patches[patches.size() - 1].superpoints[3] == patches[patches.size() - 3].superpoints[3]))
                        // that code checked 0 to 4
                        for (index_type i = 0; i < num_layers; i++)
                        {
                            if (!areWedgeSuperPointsEqual(&patches[lastPatchIdx].superpoints[i], &patches[thirdLastPatchIdx].superpoints[i]))
                            {
                                repeat_patch = false;
                                break;
                            }
                        }

                        if (repeat_patch)
                        {
                            printf("%f %f repeat_patch: %d\n",
                                   patches[lastPatchIdx].superpoints[num_layers - 1].min,
                                   patches[lastPatchIdx].superpoints[num_layers - 1].max,
                                   repeat_patch);

                            delete_patch(lastPatchIdx);

                            current_z_top_index -= 1;

                            z_top_min = Gdata.array[num_layers - 1][current_z_top_index].z;
                            z_top_min = new_z_i_atTop[layerWithSmallestShift - 1];

                            makePatch_alignedToLine(complementary_apexZ0, z_top_min, ppl, true, false);
                        }
                    }
                }
            }

            lastPatchIndex = n_patches - 1; // just to keep fresh in case we use it
            c_corner = patches[lastPatchIndex].c_corner[1];

            projectionOfCornerToBeam = straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], c_corner, patches[lastPatchIndex].c_corner[0], num_layers, 1, 0);

            saved_apexZ0 = patches[lastPatchIndex].c_corner[0];

            if (madeComplementaryPatch)
            {
                int secondLastPatchIndex = lastPatchIndex - 1;

                // modifying patches, not adding patches, so index variables do not need to be updated.
                getShadows(&patches[lastPatchIndex], z_top_min, z_top_max);
                getShadows(&patches[secondLastPatchIndex], z_top_min, z_top_max);

                float original_topR_jL = patches[secondLastPatchIndex].shadow_fromTopToInnermost_topR_jL;
                bool originalPartialTop = (original_topR_jL > complementary_apexZ0) && (original_topR_jL < apexZ0) &&
                                          (fabs(straightLineProjectorFromLayerIJtoK(&patches[secondLastPatchIndex], original_topR_jL, z_top_max, 1, num_layers, 0)) < 20 * beam_axis_lim);

                float original_topL_jL = patches[secondLastPatchIndex].shadow_fromTopToInnermost_topL_jL;
                
                bool originalPartialBottom = (original_topL_jL > complementary_apexZ0) && ((original_topL_jL - apexZ0) < -0.0001) &&
                                             (fabs(straightLineProjectorFromLayerIJtoK(&patches[secondLastPatchIndex], original_topL_jL, z_top_min, 1, num_layers, 0)) < 20 * beam_axis_lim);                
                
                float complementary_topR_jR = patches[lastPatchIndex].shadow_fromTopToInnermost_topR_jR;
                
                bool complementaryPartialTop = (complementary_topR_jR > complementary_apexZ0) && (complementary_topR_jR < apexZ0) &&
                                               (fabs(straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], complementary_topR_jR, z_top_max, 1, num_layers, 0)) < 20 * beam_axis_lim);

                float complementary_topL_jR = patches[lastPatchIndex].shadow_fromTopToInnermost_topL_jR;
                
                bool complementaryPartialBottom = (complementary_topL_jR > complementary_apexZ0) && ((complementary_topL_jR - apexZ0) < -0.0001) &&
                                                  (fabs(straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], complementary_topL_jR, z_top_min, 1, num_layers, 0)) < 20 * beam_axis_lim);

                float horizontalShiftTop = original_topR_jL - complementary_topR_jR;
                float horizontalShiftBottom = original_topL_jL - complementary_topL_jR;

                float complementary_topR_jL = patches[lastPatchIndex].shadow_fromTopToInnermost_topR_jL;
                float complementary_topL_jL = patches[lastPatchIndex].shadow_fromTopToInnermost_topL_jL;
                float original_topR_jR = patches[secondLastPatchIndex].shadow_fromTopToInnermost_topR_jR;
                float original_topL_jR = patches[secondLastPatchIndex].shadow_fromTopToInnermost_topL_jR;

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

                float z0_original_bCorner = straightLineProjectorFromLayerIJtoK(&patches[secondLastPatchIndex], apexZ0, z_top_max, 1, num_layers, 0);
                float z0_complementary_cCorner = straightLineProjectorFromLayerIJtoK(&patches[lastPatchIndex], complementary_apexZ0, z_top_min, 1, num_layers, 0);
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
                        delete_patch(n_patches - 1);
                        // decrement n_patches is handled by delete_patch
                    }

                    makePatch_alignedToLine(shifted_Align, newZtop, ppl, !shiftOriginal, false);

                    getShadows(&patches[n_patches - 1], z_top_min, z_top_max);

                    if (shiftOriginal)
                    {
                        original_topR_jL = patches[n_patches - 1].shadow_fromTopToInnermost_topR_jL;
                        original_topL_jL = patches[n_patches - 1].shadow_fromTopToInnermost_topL_jL;
                        original_topR_jR = patches[n_patches - 1].shadow_fromTopToInnermost_topR_jR;
                        original_topL_jR = patches[n_patches - 1].shadow_fromTopToInnermost_topL_jR;
                    }
                    else
                    {
                        complementary_topR_jR = patches[n_patches - 1].shadow_fromTopToInnermost_topR_jR;
                        complementary_topL_jR = patches[n_patches - 1].shadow_fromTopToInnermost_topL_jR;
                        complementary_topR_jL = patches[n_patches - 1].shadow_fromTopToInnermost_topR_jL;
                        complementary_topL_jL = patches[n_patches - 1].shadow_fromTopToInnermost_topL_jL;
                    }

                    horizontalShiftTop = original_topR_jL - complementary_topR_jR;
                    horizontalShiftBottom = original_topL_jL - complementary_topL_jR;

                    if (shiftOriginal && straightLineProjectorFromLayerIJtoK(&patches[n_patches - 1], original_topR_jR, z_top_max, 1, num_layers, 0) < beam_axis_lim)
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
                    if ((straightLineProjectorFromLayerIJtoK(&patches[n_patches - 1], shifted_Align, newZtop, 1, num_layers, 0) > beam_axis_lim) && shiftOriginal)
                    {
                        if (n_patches > 2)
                        {
                            delete_patch(n_patches - 3);
                        }
                    }
                }
            }

            z_top_max = c_corner;

            printf("+++++++++++++++++++++++ c_corner: %f\n", c_corner);
        }

        apexZ0 = patches[n_patches - 1].c_corner[0];
        apexZ0 = saved_apexZ0;
        printf("'=======================================================  z1_Align: %f\n", apexZ0);
    }
}

void makePatch_alignedToLine(float apexZ0, float z_top, int ppl, bool leftRight, bool float_middleLayers_ppl)
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

        Point temp[MAX_POINTS_PER_LAYER];
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
    //indeed, add_patch is working fine as it is copying the values over: patches[n_patches] = *curr_patch;
    //doesn't matter how wedgePatch_init works since we're dereferencing the patch to store by value in an array belonging to cover.
    add_patch(&new_patch);
}
