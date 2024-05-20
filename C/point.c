#include "header.h"

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

    printf("Failed to load point.\n");
    return 0; // failed to load
}

int comparePoints(const void *a, const void *b)
{
    float a_z = ((const Point *)a)->z;
    float b_z = ((const Point *)b)->z;
    return (a_z < b_z) ? -1 : 1; //turnary equivalent treating point equality in the second case
    /*
    if (a_z < b_z)
        return -1;
    if (a_z > b_z)
        return 1;
    return 1;
    */
}
