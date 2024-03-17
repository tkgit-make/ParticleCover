#include "header.h"

#define EVENT_C
#ifndef POINT_C
	#include "point.c"
#endif

//an event only has a list of points
typedef struct
{
	Vector* list_of_Points;
} Event;

void Event_init(Event* e)
{
	e->list_of_Points = VectorOf_Point(512, 256); //initial size and linear growth factor
}

void Event_addpoint(Event* e, Point pt)
{
	Point* newpt = PointVector_newitem(e->list_of_Points);
	*newpt = pt; //dereferencing newpt to access the memory slot for the new element, then copying the fields of the pt structure into this new Point instance at the end of the vector.
}

Point* Event_getpoint(Event* e, int n)
{
	return PointVector_getitem(e->list_of_Points, n);
}

CREATE_VECTOR_OF_T(Event) //vector alias and type-specific methods
