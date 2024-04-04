#include "header.h"

void wedgeCover_init(wedgeCover* wc, Environment* envI, DataSet* dataI) {
    wc->n_patches = 0;
    wc->env = envI;
    wc->data = dataI;
    for (int i = 0; i < MAX_PATCHES; i++) {
        wc->all_patches[i] = NULL;
        wc->real_patch_list[i] = false;
    }
    for (int i = 0; i < MAX_SUPERPOINTS_IN_COVER; i++) {
        wc->superPoints[i] = NULL;
    }
}

void add_patch(wedgeCover* cover, wedgePatch* curr_patch) {
    if (cover->n_patches == 0) {
        cover->patches[0] = *curr_patch; 
        cover->all_patches[0] = curr_patch; 
        cover->real_patch_list[0] = true; 
        cover->n_patches = 1; 
    } 
    else 
    {
        wedgePatch* prev_patch = &(cover->patches[cover->n_patches - 1]);
        bool different = false;

        for (int i = 0; i < MAX_SUPERPOINTS_IN_PATCH; i++) {
            if ((prev_patch->superpoints[i]->min != curr_patch->superpoints[i]->min) || 
                (prev_patch->superpoints[i]->max != curr_patch->superpoints[i]->max)) {
                different = true;
                break;
            }
        }

        if (different) {
            if (cover->n_patches < MAX_PATCHES) {
                cover->patches[cover->n_patches] = *curr_patch;
                cover->all_patches[cover->n_patches] = curr_patch;
                cover->real_patch_list[cover->n_patches] = true;
                cover->n_patches += 1;
            }
        }
    }
}

void delete_patch(wedgeCover* cover, int index) {
    if (index < 0 || index >= cover->n_patches) {
        return;
    }

    cover->real_patch_list[index] = false;

    for (int i = index; i < cover->n_patches - 1; i++) {
        cover->patches[i] = cover->patches[i + 1];
        cover->real_patch_list[i] = cover->real_patch_list[i + 1];
    }

    //resetting the last elements
    memset(&cover->patches[cover->n_patches - 1], 0, sizeof(wedgePatch));
    cover->real_patch_list[cover->n_patches - 1] = false;

    cover->n_patches -= 1;
}

//can't provide default parameters 
int get_index_from_z(DataSet* data, int layer, float z_value, int alignment) {
    //c doesn't support string comparison directly, using integer comparison for effiency
    //CLOSEST = 11, ABOVE = 12, BELOW = 13
    float minVal = 1000000;
    int index = 0;

    for (int i = 0; i < data->n_points[layer]; i++) {
        float diff = fabs(data->array[layer][i].z - z_value); //absolute difference
        if (diff < minVal) {
            minVal = diff;
            index = i;
        }
    }

    if (alignment == CLOSEST) {
        return index;
    }

    if (alignment == ABOVE) {
        if (data->array[layer][index].z > z_value) {
            return index;
        }
        //potential bounds issue
        return index+1;
    }

    if (alignment == BELOW) {
        if (data->array[layer][index].z < z_value) {
            return index;
        }
        return index-1;
    }

    return index;
}

//not implementing the logic corresponding to show = true (that would involve Line Generators, etc)
//Line Generator and its accompanying methods have been coded, but we are not going to implement show=true case here as main method passes with show=false.
//MAKE_PATCHES_SHADOW_QUILT_FROM_EDGES defined as 33
void solve(wedgeCover* cover, int lining, float apexZ0, int ppl, int nlines, bool leftRight) {
    for (int i = 0; i < cover->env->num_layers; i++) {
        bool foundIdentical = false;
        bool firstTime = true;

        while (foundIdentical || firstTime) {
            foundIdentical = false;
            for (int x = 0; x < cover->data->n_points[i] - 1; x++) {
                if (cover->data->array[i][x].z == cover->data->array[i][x + 1].z) {
                    cover->data->array[i][x + 1].z += 0.00001;
                    foundIdentical = true;
                }
            }

            firstTime = false;
            if (foundIdentical) {
                qsort(cover->data->array[i], cover->data->n_points[i], sizeof(Point), comparePoints);
            }
        }
    }

    if (lining == MAKE_PATCHES_SHADOW_QUILT_FROM_EDGES) {
        makePatches_ShadowQuilt_fromEdges(cover, apexZ0, 1, ppl, leftRight);
    }
}
