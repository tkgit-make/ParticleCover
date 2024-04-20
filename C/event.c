#include "header.h"

void Event_init(Event* e, Environment* envI, Point* pointsArray, int numPoints) {
    index_type old_count = e->count;
	e->env = envI;
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

bool isUniqueRadius(float radius, float *uniqueRadii, int count) {
    for (int i = 0; i < count; i++) {
        if (uniqueRadii[i] == radius) {
            return false;  // not unique
        }
    }
    return true;  // unique
}

index_type Event_load(Event* e) {
    index_type n = 0;
    char ch = ',';
    float uniqueRadii[MAX_POINTS_IN_EVENT];  // Assuming no more than MAX_POINTS_IN_EVENT unique radii
    int numUniqueRadii = 0;
    float radius;

    while ((ch == ',') && (n < MAX_POINTS_IN_EVENT)) {
        if (Point_load(&e->points[n], &radius) < 1) break;  // Modified Point_load to return radius via a pointer argument
        if (isUniqueRadius(radius, uniqueRadii, numUniqueRadii)) {
            uniqueRadii[numUniqueRadii++] = radius;  // Add to unique radii if it's new
        }
        n++;
        scanf("%c", &ch);
    }

    e->count = n;

    // Optionally, you can initialize the environment here with the unique radii and numUniqueRadii
    initEnvironment(&e->env, 100.0, 15.0, numUniqueRadii, uniqueRadii);  // Assuming certain values for top and beam limits

    return n > 0;
}