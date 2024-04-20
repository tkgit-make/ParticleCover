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

    //int ik = 0;

    //added
    int numEventsLoaded = 0;

    G_stats.count_events = 0;
    G_stats.count_points = 0;

    FILE *myfile;
    myfile = fopen("C/cOutput.txt", "w"); //writing

    if (myfile == NULL) {
        printf("Error opening file");
        return -1; 
    }

    int k = wedges[0];
	while (/*Event_load(&G_event) > 0 &&*/ k < wedges[1])
	{
        Event_load(&G_event);
		ProcessEvent();
        printf(myfile, "wedge: %d\n", k);
        fprintf(myfile, "wedge: %d\n", k);

        Environment env = G_event.env; 
        Point* points = G_event.points;
        int num_points = G_event.count;

        env.top_layer_lim = top_layer_cutoff;
        env.beam_axis_lim = z0_luminousRegion;

        DataSet data;
        initDataSet(&data, env); 

        //no need to implement logic if show_acceptance_of_cover == true
        if (!uniform_N_points) {
        importData(&data, points, num_points);
        } else {
            //(1) will not run for our input. (2) else does nothing even if ran
        }

        addBoundaryPoint(&data, 0.0001); //with default param

        wedgeCover cover;
        initWedgeCover(&cover, &data);

        solve(&cover, lining, apexZ0, ppl, 100, leftRightAlign); //solve modifies cover

        //num_covers and num_all_patches not needed

        //revisit rounding
        for (int i = 0; i < (&cover)->n_patches; i++) {
            fprintf(myfile, "Patch\n");
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topL_jL * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topL_jR * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topR_jL * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topR_jR * 10000));

            for (int j = 0; j < cover.patches[i].superpoint_count; j++) {
                fprintf(myfile, "Superpoint\n");
                for (int r = 0; r < cover.patches[i].superpoints[j]->point_count; r++) {
                    Point currentPt = cover.patches[i].superpoints[j]->points[r];
                    fprintf(myfile, "%d %.4f %d %.4f\n", 
                            currentPt.layer_num, 
                            currentPt.phi, 
                            (int)currentPt.radius, 
                            currentPt.z);
                }
            }
        }
        //revisit rounding
        for (int i = 0; i < cover.n_patches; i++) {
            fprintf(myfile, "[%ld, %ld]\n", 
                    lround(cover.patches[i].a_corner[0] * 10000), 
                    lround(cover.patches[i].a_corner[1] * 10000));
            fprintf(myfile, "[%ld, %ld]\n", 
                    lround(cover.patches[i].b_corner[0] * 10000), 
                    lround(cover.patches[i].b_corner[1] * 10000));
            fprintf(myfile, "[%ld, %ld]\n", 
                    lround(cover.patches[i].c_corner[0] * 10000), 
                    lround(cover.patches[i].c_corner[1] * 10000));
            fprintf(myfile, "[%ld, %ld]\n", 
                    lround(cover.patches[i].d_corner[0] * 10000), 
                    lround(cover.patches[i].d_corner[1] * 10000));
            fprintf(myfile, "\n");
        }
        //instead of making an array of all events and passing them in, we only need access to them individually, so we will loop through and process as we create them.
        k++;
    }  

    fclose(myfile);

    printf("Events: %d, Points: %d\n", G_stats.count_events, G_stats.count_points);

}