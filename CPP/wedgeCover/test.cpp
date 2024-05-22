//
// Created by Tanish Kumar on 1/8/24.
//

#include<iostream>
#include "makePatch_alignedtoLine.cpp"

using namespace std;

class Tester
{
public:
    void wedge_test(string str = "makePatches_Projective_center", double apexZ0 = 0, double z0_spacing = 0.5, double z0_luminousRegion = 15.0, vector<int> wedges = { 0, 128 }, int lines = 1000, string v = "v3", double top_layer_cutoff = 50.0, double accept_cutoff = 10.0, bool leftRightAlign = true, bool uniform_N_points = false, string acceptance_method = "Analytic", bool show_acceptance_of_cover = false, bool movie = false, bool savefig = false, int figSizeScale = 6, int movieFigSizeScale = 3)
    {
        accept_cutoff = z0_luminousRegion;

        bool showZimperfect = false;

        if(wedges[1] - wedges[0] == 1)
        {
            showZimperfect = true;
        }

        if(wedges[1] - wedges[0] > 50)
        {
            show_acceptance_of_cover = false;
            z0_spacing = 0.2;
        }

        vector<int> num_covers;
        vector<int> num_all_patches;

        vector< vector<int> > PRF;

        string data_string = v + " events";

        if(uniform_N_points != false)
        {
            data_string = "Uniform 1 points";
            wedges = {0, 1};
        }

        vector<double> zInnerLayer;

        for(double i = -22; i < 22 + z0_spacing; i += z0_spacing)
        {
            zInnerLayer.push_back(i);
        }

        vector<double> z0Array;

        for(double i = -1 * z0_luminousRegion; i < z0_luminousRegion + z0_spacing; i += z0_spacing)
        {
            z0Array.push_back(i);
        }

        vector< vector<int> > mean_list;

        for(int i = 0; i < (wedges[1] - wedges[0]); i++)
        {
            vector<int> vect(z0Array.size(), 0);
            mean_list.push_back(vect);
        }

        vector<double> z0Imperfect; //check data types
        vector<double> z0OverEfficiency;

        vector<Event> all_data = FileReader::readFile("wedgeData_" + v + "_128.txt", wedges[1]);

        int ik = 0;

        for(int k = wedges[0]; k < wedges[1]; k++)
        {
            cout << "wedge: " + k << endl;

            Environment env = all_data[k].env;
            vector<Point> points = all_data[k].list_of_Points;

            env = Environment(top_layer_cutoff, z0_luminousRegion);
            DataSet data(env);

            if(show_acceptance_of_cover)
            {
                cout << "unimplemented" << endl; //fix
                throw "show_acceptance_of_cover unimplemented";
            }
            if(uniform_N_points == false)
            {
                data.importData(points);
            }
            else
            {
                //write generateUniform
                vector<bool> vect(5, true);
                //data.generateUniform(vect);
            }

            data.addBoundaryPoint();

            wedgeCover cover(env, data);

            //cover.solve(apexZ0 = apexZ0, lining = lining, ppl = ppl, leftRight = leftRightAlign, show = False);

            //num_covers.append(cover.n_patches);
            //num_all_patches.append(len(cover.all_patches));

            ik++;
        }
    }
};