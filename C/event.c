#include "header.h"

index_type Event_load(Event *e)
{
    index_type n = 0;
    char ch = ',';

    while ((ch == ',') && (n < MAX_POINTS_IN_EVENT))
    {
        if (Point_load(&e->points[n]) < 1)
            break;
        n++;
        scanf("%c", &ch);
    }

    e->count = n;

    return n > 0;
}