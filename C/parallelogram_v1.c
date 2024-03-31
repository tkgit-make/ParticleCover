#include "header.h"

/*
void init_parallelogram_v1(Parallelogram_v1 *pg, int layer_numI, float top_layer_zminI, float top_layer_zmaxI, float shadow_topR_jLI, float shadow_topR_jRI, float pSlopeI) {
    pg->layer_num = layer_numI;
    pg->pSlope = pSlopeI;

    pg->shadow_topR_jL = shadow_topR_jLI;
    pg->shadow_topR_jR = shadow_topR_jRI;

    float delta_ztop = top_layer_zmaxI - top_layer_zminI;
    float delta_z0 = delta_ztop * pSlopeI;

    pg->shadow_topL_jL = shadow_topR_jLI - delta_z0;
    pg->shadow_topL_jR = shadow_topR_jRI - delta_z0;

    pg->top_layer_zmin = top_layer_zminI;
    pg->top_layer_zmax = top_layer_zmaxI;

    if ((pg->top_layer_zmax > 100) || (pg->top_layer_zmin < -100)) {
        printf("z1_min or z1_max is out of expected range of +-100");
    }

}
*/
