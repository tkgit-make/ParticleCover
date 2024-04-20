#define MAIN_C
#include "header.h"

typedef struct {
    int count_events;
    int count_points;
} EventStats;

EXTERN EventStats G_stats; 

Event events[MAX_EVENTS_TO_READ];

void ProcessEvent(Event event) {
    //printf("Event loaded, %d Points\n", event.count);
    G_stats.count_events++;
    G_stats.count_points += event.count;
}

int main() {
    int numEventsLoaded = 0;

    G_stats.count_events = 0;
    G_stats.count_points = 0;

    while (numEventsLoaded < MAX_EVENTS_TO_READ && Event_load(&events[numEventsLoaded]) > 0) {
        numEventsLoaded++;
    }

    for (int i = 0; i < numEventsLoaded; i++) {
        ProcessEvent(events[i]);
    }

    printf("Events: %d, Points: %d\n", G_stats.count_events, G_stats.count_points);

	assist(events);
    return 0;
}

int assist(Event* events) {
	const char* lining = "makePatches_Projective_center";
    float apexZ0 = 0.0;
    float z0_spacing = 0.5;
    int ppl = 16;
    float z0_luminousRegion = 15.0;
    int wedges[] = {0, 128};  
    int wedge_count = 2; //array size
    int lines = 1000;
    const char* v = "v3";
    float top_layer_cutoff = 50.0;
    float accept_cutoff = 10.0;
    bool leftRightAlign = true;
    bool uniform_N_points = false;
    const char* acceptance_method = "Analytic";
    bool show_acceptance_of_cover = false;
    bool movie = false;
    bool savefig = false;
    int figSizeScale = 6;
    int movieFigSizeScale = 3;
	wedge_test(events, lining, apexZ0, z0_spacing, ppl, z0_luminousRegion, wedges, wedge_count, lines, v, top_layer_cutoff, accept_cutoff, leftRightAlign, uniform_N_points, acceptance_method, show_acceptance_of_cover, movie, savefig, figSizeScale, movieFigSizeScale);
	return 0;
}