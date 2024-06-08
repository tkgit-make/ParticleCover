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
    float phi_temp, z_temp;
    index_type radius_temp;
    index_type layer_num_temp;

    if (scanf("(%d,%d,%f,%f)", &layer_num_temp, &radius_temp, &phi_temp, &z_temp) == 4)
    {
        p->layer_num = layer_num_temp;
        p->phi = phi_temp * CONVERSION_FACTOR;
        p->z = z_temp * CONVERSION_FACTOR;
        p->radius = radius_temp * CONVERSION_FACTOR;
        //p->radius = layer_num_temp*5*CONVERSION_FACTOR;

        return 1; // successful load
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
