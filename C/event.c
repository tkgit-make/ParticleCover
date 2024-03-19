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

	printf("Reached maximum number of points");

	e->count = n;
	return n>0;
}
