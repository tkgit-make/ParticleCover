#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include <regex>
#include <sstream>
#include <climits>

using namespace std;

class Point
{
public:
    int layer_num;
    double radius;
    double phi;
    double z;

    Point(int layerNum, double rad, double ph, double zVal)
    {
        layer_num = layerNum;
        radius = rad;
        phi = ph;
        z = zVal;
    }
};

class Environment
{
public:
    double top_layer_lim;
    double beam_axis_lim;
    int num_layers;
    vector<double> radii;
    vector<double> parallelogramSlopes;
    vector<double> radii_leverArm;
    double boundaryPoint_offset;
    vector<double> trapezoid_edges;

    Environment(double top_layer_limI = 100.0, double beam_axis_limI = 15.0, int num_layersI = 5, vector<double> radiiI = { 5.0, 10.0, 15.0, 20.0, 25.0 })
    {
        if(top_layer_limI < beam_axis_limI)
        {
            throw "The top layer limits cannot be smaller than the bottom layer limits.";
        }

        top_layer_lim = top_layer_limI;
        beam_axis_lim = beam_axis_limI;

        if(radiiI.size() != num_layersI)
        {
            throw "The radii do not match the number of layers.";
        }

        num_layers = num_layersI;
        radii = radiiI;

        sort(radii.begin(), radii.end());

        for(int i = 0; i < radii.size() - 1; i++)
        {
            double currentVal = (radii[0] - radii[i]) / (radii[radii.size() - 1] - radii[i]);
            parallelogramSlopes.push_back(currentVal);
        }

        for(int i = 0; i < parallelogramSlopes.size(); i++)
        {
            radii_leverArm.push_back(1 - parallelogramSlopes[i]);
        }

        boundaryPoint_offset = 0;

        for(int i = 0; i < radii.size(); i++)
        {
            double currentVal = radii[i] * (top_layer_lim - beam_axis_lim) / radii[radii.size() - 1] + beam_axis_lim;
            trapezoid_edges.push_back(currentVal);
        }
    }
};

class DataSet
{
public:
    Environment * env;
    vector< vector<Point> > array;
    vector<int> n_points;
    int total_points;
    double boundaryPoint_offset;

    DataSet()
    {
        total_points = 0;
    }

    DataSet(Environment& envI)
    {
        env = &envI;

        for(int i = 0; i < 5; i++)
        {
            vector<Point> vect;
            array.push_back(vect);
        }

        for(int i = 0; i < env->num_layers; i++)
        {
            n_points.push_back(0);
        }

        total_points = 0;
    }

    void importData(vector<Point> data_array)
    {
        total_points = data_array.size();

        for(int i = 0; i < data_array.size(); i++)
        {
            array[data_array[i].layer_num - 1].push_back(data_array[i]);
        }

        int ln = 0;

        for(int i = 0; i < array.size(); i++)
        {
            sort(array[i].begin(), array[i].end(), [](Point &a, Point &b){ return a.z < b.z; });
            n_points[i] = array[i].size();
            ln++;
        }
    }

    void addBoundaryPoint(double offset = 0.0001)
    {
        boundaryPoint_offset = offset;

        for(int i = 0; i < env->trapezoid_edges.size(); i++)
        {
            double phi0 = array[i][0].phi;

            array[i].insert(array[i].begin(), Point(i + 1, (i + 1) * 5, phi0, -1 * env->trapezoid_edges[i] - offset));
            array[i].push_back(Point(i + 1, (i + 1) * 5, phi0, env->trapezoid_edges[i] + offset));

            n_points[i] += 2;
        }

        total_points = array.size();

        int ln = 0;

        for(int i = 0; i < array.size(); i++)
        {
            sort(array[i].begin(), array[i].end(), [](Point &a, Point &b){ return a.z < b.z; });
            n_points[ln] = array[i].size();
            ln++;
        }

        for(int i = 0; i < env->trapezoid_edges.size(); i++)
        {
            env->trapezoid_edges[i] += offset;
        }
    }
};

class wedgeSuperPoint
{
public:
    vector<Point> points;
    vector<double> z_values;
    double min;
    double max;

    wedgeSuperPoint(vector<Point> pointsI)
    {
        if(pointsI.size() != 16)
        {
            if((pointsI.size() != 32) && (pointsI.size() != 31))
            {
                throw "This patch does not have 16 or 32/31 points in each layer";
            }
        }

        vector<double> z_list;

        for(int i = 0; i < pointsI.size(); i++)
        {
            z_list.push_back(pointsI[i].z);
        }

        z_values = z_list;
        points = pointsI;
        max = *max_element(z_list.begin(), z_list.end());
        min = *min_element(z_list.begin(), z_list.end());
    }
    /*
    bool contains(double p)
    {
        try
        {
            p = p.z;
        }
        catch
        {
            p = p;
        }

        if((min <= p) && (max >= p))
        {
            return true;
        }

        return false;
    } */

    bool eq(wedgeSuperPoint other)
    {
        if((min == other.min) && (max == other.max))
        {
            return true;
        }

        return false;
    }
};

class wedgePatch
{
public:
    Environment env;
    int end_layer;
    int left_end_layer;
    int right_end_layer;
    int left_end_lambdaZ;
    int right_end_lambdaZ;
    double apexZ0;

    int shadow_fromTopToInnermost_topL_jL;
    int shadow_fromTopToInnermost_topL_jR;
    int shadow_fromTopToInnermost_topR_jL;
    int shadow_fromTopToInnermost_topR_jR;

    vector<wedgeSuperPoint> superpoints;

    wedgePatch(Environment envI, vector<wedgeSuperPoint> superpointsI, double apexZ0I)
    {
        env = envI;
        end_layer = -1;
        left_end_layer = -1;
        right_end_layer = -1;
        left_end_lambdaZ = NULL;
        right_end_lambdaZ = NULL;
        apexZ0 = apexZ0I;

        shadow_fromTopToInnermost_topL_jL = NULL;
        shadow_fromTopToInnermost_topL_jR = NULL;
        shadow_fromTopToInnermost_topR_jL = NULL;
        shadow_fromTopToInnermost_topR_jR = NULL;

        if(superpointsI.size() != env.num_layers)
        {
            throw "The patch layers does not match environment layers. ";
        }

        superpoints = superpointsI;

        /*

        self.getParallelograms()
        self.getParallelograms_v1()
        self.get_acceptanceCorners()
        self.get_end_layer()

        */
    }
};

class wedgeCover
{
public:
	int n_patches; 
	vector<wedgePatch> patches;
	Environment env; 
	DataSet data; 
	vector<int> fitting_lines; //Check the data types for all of the 4 lists here
    vector<wedgeSuperPoint> superPoints;
    vector<wedgePatch> all_patches;
    vector<bool> real_patch_list; 

	wedgeCover(Environment envI, DataSet dataI)
	{
		n_patches = 0; 
		env = envI; 
		data = dataI; 
	}

	void add_patch(wedgePatch curr_patch)
	{
		if(n_patches == 0)
		{
			patches.push_back(curr_patch); 
			all_patches.push_back(curr_patch); 
			real_patch_list.push_back(true); 
			n_patches += 1; 
		}
		else
		{
			wedgePatch prev_patch = patches.back(); 
			vector<wedgeSuperPoint> prev_sp = prev_patch.superpoints;
			vector<wedgeSuperPoint> curr_sp = curr_patch.superpoints; 

			for(int i = 0; i < prev_sp.size(); i++)
			{
				if((prev_sp[i].min != curr_sp[i].min) || (prev_sp[i].max != curr_sp[i].max))
				{
					patches.push_back(curr_patch);
					all_patches.push_back(curr_patch);
					real_patch_list.push_back(true); 
					n_patches += 1; 
					break; 
				}
			}
		}
	}

	void makePatch_alignedToLine(double apexZ0 = 0, double z_top = -50, int ppl = 16, bool leftRight = true, bool double_middleLayers_ppl = false)
	{
		vector<wedgeSuperPoint> init_patch; 
		int original_ppl = ppl; 
		double alignmentAccuracy = 0.00001; 

		vector< vector<Point> > row_data = data.array;

		for(int i = 0; i < env.num_layers; i++)
		{
			double y = env.radii[i]; 

			vector<double> row_list; 

			for(int j = 0; j < row_data[i].size(); j++)
			{
				row_list.push_back(row_data[i][j].z); 
			} 

			double r_max = env.radii.back(); 

			double projectionToRow = (z_top - apexZ0) * (y - env.radii[0]) / (r_max - env.radii[0]) + apexZ0;

			int start_index = 0;
			double start_value = INT_MAX;

			for(int j = 0; j < row_list.size(); j++)
			{
				if(abs(row_list[j] - projectionToRow) < abs(start_value))
				{
					start_index = j; 
					start_value = row_list[j] - projectionToRow;
				}
			}

			int left_bound = 0;
            double lbVal = INT_MAX;
			int right_bound = 0;
            double rbVal = INT_MAX;

            for(int j = 0; j < row_list.size(); j++)
			{
				if(abs((row_list[j] + env.trapezoid_edges[i] + env.boundaryPoint_offset)) < lbVal)
				{
					left_bound = j;
                    lbVal = abs((row_list[j] + env.trapezoid_edges[i] + env.boundaryPoint_offset));
				}

				if(abs((row_list[j] - env.trapezoid_edges[i] - env.boundaryPoint_offset)) < rbVal)
				{
					right_bound = j;
                    rbVal = abs((row_list[j] - env.trapezoid_edges[i] - env.boundaryPoint_offset));
				}
			}

			if ((double_middleLayers_ppl == true) && (i != 0) && (i != env.num_layers - 1))
			{
				ppl = original_ppl * 2 - 1;
			}
			else
			{
				ppl = original_ppl; 
			}

			if (leftRight == true)
			{
				if(start_index != 0)
				{
					if(start_value > alignmentAccuracy)
					{
						start_index -= 1; 
					}
				}

				if((start_index + ppl) > (right_bound + 1))
				{
					init_patch.push_back(wedgeSuperPoint({row_data[i].begin() + right_bound + 1 - ppl, row_data[i].begin() + right_bound + 1})); 
				}
				else
				{
					init_patch.push_back(wedgeSuperPoint({row_data[i].begin() + start_index, row_data[i].begin() + start_index + ppl})); 
				}
			}
			else
			{
				if(start_index != (row_list.size() - 1))
				{
					cout << "row " + to_string(i + 1) + " start_index " + to_string(start_index) + " start_value " + to_string(start_value) + " z: " + to_string(row_list[start_index]) << endl;
                    if(start_value < -1 * alignmentAccuracy)
					{
						start_index += 1; 
						start_value = row_list[start_index] - projectionToRow;
						cout << "row " + to_string(i + 1) + " updated start_index " + to_string(start_index) + " start_value " + to_string(start_value) + " z: " + to_string(row_list[start_index]) << endl;
					}
				}

				if((start_index - ppl + 1) < left_bound)
				{
					init_patch.push_back(wedgeSuperPoint({row_data[i].begin() + left_bound, row_data[i].begin() + left_bound + ppl})); 
				}
				else
				{
					init_patch.push_back(wedgeSuperPoint({row_data[i].begin() + start_index + 1 - ppl, row_data[i].begin() + start_index + 1})); 
				}
			}
		}

        add_patch(wedgePatch(env, init_patch, apexZ0 = apexZ0));
    }

    void tester()
    {
        double apexZ0 = 7.75751;
        int ppl = 16;
        double z_top = -8.883753333333333;
        bool leftRight = false;

        makePatch_alignedToLine(apexZ0, z_top, ppl, leftRight);
    }
};

class Event
{
public:

    Environment env;
    vector<Point> list_of_Points;

    Event(Environment envI = NULL, vector<Point> listP = {})
    {
        env = envI;
        list_of_Points = listP;
    }
};

class FileReader
{
public:
    static vector<string> splitString(string str, string splitter = "),(")
    {
        vector<string> result;
        string currentStr = "";
        for (int i = 0; i < str.size(); i++)
        {
            bool flag = true;
            for (int j = 0; j < splitter.size(); j++)
            {
                if (str[i + j] != splitter[j]) flag = false;
            }
            if (flag) {
                if (currentStr.size() > 0) {
                    result.push_back(currentStr);
                    currentStr = "";
                    i += splitter.size() - 1;
                }
            }
            else
            {
                currentStr += str[i];
            }
        }
        result.push_back(currentStr);
        return result;
    }

    static vector<Event> readFile(string filepath, int stop = 128, bool performance = false)
    {
        ifstream currentFile;
        vector<Event> events;
        currentFile.open(filepath);
        string line;
        if(currentFile.is_open())
        {
            int line_index = 0;

            while(getline(currentFile, line))
            {
                line = regex_replace(line, regex("(^[ ]+)|([ ]+$)"),"");
                if(!line.empty())
                {
                    line = line.substr(1, line.size() - 2);

                    vector<string> tuples = splitString(line);
                    vector< vector<string> > finalTuples;

                    for(int i = 0; i < tuples.size(); i++)
                    {
                        vector<string> temp = splitString(tuples[i], ",");
                        finalTuples.push_back(temp);
                    }

                    vector<Point> list_of_Points;

                    for(int i = 0; i < finalTuples.size(); i++)
                    {
                        vector<string> ct = finalTuples[i];
                        Point temp = Point(stoi(ct[0]), stof(ct[1]), stof(ct[2]), stof(ct[3]));
                        list_of_Points.push_back(temp);
                    }

                    vector<double> radii;

                    for(int i = 0; i < list_of_Points.size(); i++)
                    {
                        radii.push_back(list_of_Points[i].radius);
                    }

                    sort(radii.begin(), radii.end());

                    int num_layers = radii.size();

                    Environment env = Environment(100.0, 15.0, num_layers, radii);

                    Event evt(env, list_of_Points);

                    events.push_back(evt);

                    line_index++;

                    if(line_index == stop)
                    {
                        break;
                    }
                }
            }

            currentFile.close();
        }
        else
        {
            cout << "Error opening file." << endl;
        }

        return events;
    }
};

int main()
{
	string filepath = "wedgeData_v3_128.txt";

    vector<Event> events = FileReader::readFile(filepath);

    Environment env = events[0].env;
    vector<Point> points = events[0].list_of_Points;

    env = Environment(50);
    DataSet ds(env);

    ds.importData(points);
    ds.addBoundaryPoint();

    wedgeCover cov(env, ds);
    cov.tester();
};
