#include "header.h"

void Event_init(Event* e) {
    e->count = 0;  //no points added yet
}

int Event_load(Event* e)
{
	int n = 0;
	char ch=',';

	while ((ch == ',') && (n < MAX_POINTS_IN_EVENT)) {
		if (Point_load(&e->points[n]) < 1) break;
		n++;

		scanf("%c", &ch);
		// event ends with a '\n' new line (continues with a ',')
	}

	AssertWithMessage(!((ch==',') && (n == MAX_POINTS_IN_EVENT)),"Reached maximum number of points");

	e->count = n;
	return n>0;
}

/*
void Event_addpoint(Event* e, Point pt) {
    if (e->count >= MAX_POINTS_IN_EVENT) {
        printf("Failure to add, array full");
        return;
    }
    e->points[e->count] = pt;  //dereferencing and assigning value
    e->count++;  //deferencing and incrementing
}

Point* Event_getpoint(Event* e, int index) {
    if (index < 0 || index >= e->count) {
        printf("Out of bounds or not yet assigned");
        return NULL;
    }
    return &(e->points[index]); //return pointer to the point
}
*/