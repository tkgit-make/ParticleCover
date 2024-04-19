#define MAIN_C
#include "header.h"

//main event loop.
//only global memory and minimal stack is utilized during the event loop

typedef struct
{
	int count_events;
	int count_points;
} EventStats;

EXTERN EventStats G_stats; 

void ProcessEvent()
{
	//printf("New event loaded, %d Points\n",G_event.count);
	G_stats.count_events++;
	G_stats.count_points += G_event.count;
}

int main()
{
	//event loop
	while (Event_load(&G_event) > 0)
	{
		ProcessEvent();
	}

	printf("Events: %d, Points: %d\n",G_stats.count_events,G_stats.count_points);
}