#include "header.h"

void initParallelogram(Parallelogram* pg, int layer_numI, float z1_minI, float z1_maxI, 
                       float shadow_bottomL_jRI, float shadow_bottomR_jRI, 
                       float shadow_bottomL_jLI, float shadow_bottomR_jLI, 
                       float pSlopeI) {
    pg->layer_num = layer_numI;
    pg->pSlope = pSlopeI;

    pg->shadow_bottomL_jR = shadow_bottomL_jRI;
    pg->shadow_bottomR_jR = shadow_bottomR_jRI;

    pg->shadow_bottomL_jL = shadow_bottomL_jLI;
    pg->shadow_bottomR_jL = shadow_bottomR_jLI;

    pg->z1_min = z1_minI;
    pg->z1_max = z1_maxI;

    if ((pg->z1_min > 22.0) || (pg->z1_max < -22.0)) {
        printf("z1_min or z1_max is out of expected range of +-22");
    }
}
