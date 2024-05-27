#include "header.h"

void wedge_test(float apexZ0, float z0_spacing, int ppl, float z0_luminousRegion, int wedges[], int wedge_count, int lines, float top_layer_cutoff, float accept_cutoff)
{
    int numEventsLoaded = 0;

    FILE *myfile;
    myfile = fopen("C/cOutput.txt", "w"); 

    if (myfile == NULL)
    {
        printf("Error opening file");
        return;
    }

    for (index_type z = 0; z < wedges[1]; z++)
    {
        //Event event;
        //Event_load(&event);
        if(z<wedges[0]) continue;
        printf("wedge %d\n", z); //main print
        fprintf(myfile, "wedge %d\n", z); //file to diff

        //Point *points = event.points;

        //DataSet Gdata;
        initDataSet(&Gdata);
        
        //uniform_N_points, previously a parameter, is false, so we importData
        importData(&Gdata);
        
        addBoundaryPoint(&Gdata, 0.0001); // with default param

        wedgeCover cover;
        initWedgeCover(&cover);

        solve(&cover, apexZ0, ppl, 100, false); // solve modifies cover. false is from the left right align (previously a parameter in wedge test)

        for (int i = 0; i < (&cover)->n_patches; i++)
        {
            fprintf(myfile, "Patch \n");
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topL_jL * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topL_jR * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topR_jL * 10000));
            fprintf(myfile, "%ld\n", lround(cover.patches[i].shadow_fromTopToInnermost_topR_jR * 10000));

            for (int j = 0; j < cover.patches[i].superpoint_count; j++)
            {
                fprintf(myfile, "Superpoint \n");
                for (int r = 0; r < cover.patches[i].superpoints[j].point_count; r++)
                {
                    Point currentPt = cover.patches[i].superpoints[j].points[r];
                    fprintf(myfile, "%d %.4f %d %.4f\n",
                            currentPt.layer_num,
                            currentPt.phi,
                            (int)currentPt.radius,
                            currentPt.z);
                }
            }
        }
        for (int i = 0; i < cover.n_patches; i++)
        {
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
        // instead of making an array of all events and passing them in, we only need access to them individually, so we will loop through and process as we create them.
    }

    fclose(myfile);
}
