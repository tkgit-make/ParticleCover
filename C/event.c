#include "header.h"

void Event_init(Event* e, Environment* envI, Point* pointsArray, int numPoints) {
    index_type old_count = e->count;
	e->count = 0;
	index_type i = 0;
    for (; i < numPoints && i < MAX_POINTS_IN_EVENT; i++) {
        e->points[i] = pointsArray[i]; 
        e->count++;
    }

	for (; i < old_count; i++) {
        //set each remaining field to zero if there were more indices filled before than what we are now. clearing out any leftovers
        e->points[i].layer_num = 0;
        e->points[i].radius = 0.0;
        e->points[i].phi = 0.0;
        e->points[i].z = 0.0;
    }

}

index_type Event_load(Event* e)
{
	index_type n = 0;
	char ch=',';

	while ((ch == ',') && (n < MAX_POINTS_IN_EVENT)) {
		if (Point_load(&e->points[n]) < 1) break;
		n++;

		scanf("%c", &ch);
		//if(scanf("%c", &ch) < 1) break; //error reading, don't use in conjunction with previous line
		// event ends with a '\n' new line (continues with a ',')
	}

	//printf("Filled event [in event.c]");

	e->count = n;
	return n>0;
}