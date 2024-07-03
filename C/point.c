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
    
    if (scanf("(%d,%d,%f,%f)", &p->layer_num, &p->radius, &p->phi, &p->z) == 4)
    {
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
