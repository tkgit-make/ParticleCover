#include "header.h"

void wedge_test(Event* events, const char* lining, float apexZ0, float z0_spacing, int ppl, float z0_luminousRegion, int wedges[], int wedge_count, int lines, const char* v, float top_layer_cutoff, float accept_cutoff, bool leftRightAlign, bool uniform_N_points, const char* acceptance_method, bool show_acceptance_of_cover, bool movie, bool savefig, int figSizeScale, int movieFigSizeScale) {
    //set default values for wedges if not given anything usable (would be unusual)
    if (wedge_count == 0) {
        static int default_wedges[] = {0, 128};
        wedges = default_wedges;
        wedge_count = 2; 
    }
    accept_cutoff = z0_luminousRegion;

    bool showZimperfect = false;

    if (wedges[1] - wedges[0] == 1) {
        showZimperfect = true;
    }

    if (wedges[1] - wedges[0] > 50) {
        show_acceptance_of_cover = false;
        z0_spacing = 0.2;
    }

    int num_covers[128]; //change
    int num_all_patches[128]; //change

    // Check if uniform points are needed
    if (uniform_N_points) {
        wedges[0] = 0;
        wedges[1] = 1;
    }
    //PRF variable not needed. data_string not needed
    /*
    //not needed
    #define MAX_Z_INNER_LAYER_SIZE 100 //move later
    #define MAX_Z0_ARRAY_SIZE 100  //move later

    float zInnerLayer[MAX_Z_INNER_LAYER_SIZE];
    float z0Array[MAX_Z0_ARRAY_SIZE];
    
    int zInnerLayerCount = 0;
    for (float i = -22; i < 22 + z0_spacing && zInnerLayerCount < MAX_Z_INNER_LAYER_SIZE; i += z0_spacing) {
        zInnerLayer[zInnerLayerCount++] = i;
    }

    int z0ArrayCount = 0;
    for (float i = -z0_luminousRegion; i < z0_luminousRegion + z0_spacing && z0ArrayCount < MAX_Z0_ARRAY_SIZE; i += z0_spacing) {
        z0Array[z0ArrayCount++] = i;
    }
    */
   //mean_list, vect not needed
   //z0Imperfect, z0OverEfficiency not needed

   //instead of getting all_data from readFile, we will pass it in as a parameter from main, which will have an array of events already from reading the data.

    int ik = 0;

    //added
    int numEventsLoaded = 0;

    G_stats.count_events = 0;
    G_stats.count_points = 0;

	while (Event_load(&G_event) > 0)
	{
			ProcessEvent();
            //add more code. instead of making an array of all events and passing them in, we only need access to them individually, so we will loop through and process as we create them.
	}

    printf("Events: %d, Points: %d\n", G_stats.count_events, G_stats.count_points);

}