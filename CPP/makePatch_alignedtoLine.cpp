#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include <regex>
#include <sstream>
#include <set>

using namespace std;

class Point
{
public:
    int layer_num;
    float radius;
    float phi;
    float z;

    Point(int layerNum, float rad, float ph, float zVal)
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
    float top_layer_lim;
    float beam_axis_lim;
    int num_layers;
    vector<float> radii;
    vector<float> parallelogramSlopes;
    vector<float> radii_leverArm;
    float boundaryPoint_offset;
    vector<float> trapezoid_edges;

    Environment(float top_layer_limI = 100.0, float beam_axis_limI = 15.0, int num_layersI = 5, vector<float> radiiI = { 5.0, 10.0, 15.0, 20.0, 25.0 })
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
            float currentVal = (radii[0] - radii[i]) / (radii[radii.size() - 1] - radii[i]);
            parallelogramSlopes.push_back(currentVal);
        }

        for(int i = 0; i < parallelogramSlopes.size(); i++)
        {
            radii_leverArm.push_back(1 - parallelogramSlopes[i]);
        }

        boundaryPoint_offset = 0;

        for(int i = 0; i < radii.size(); i++)
        {
            float currentVal = radii[i] * (top_layer_lim - beam_axis_lim) / radii[radii.size() - 1] + beam_axis_lim;
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
    float boundaryPoint_offset;

    DataSet()
    {
        total_points = 0;
    }

    DataSet(Environment& envI)
    {
        env = &envI;
        total_points = 0;

        for(int i = 0; i < 5; i++)
        {
            vector<Point> vect;
            array.push_back(vect);
        }

        for(int i = 0; i < env->num_layers; i++)
        {
            n_points.push_back(0);
        }

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

    void addBoundaryPoint(float offset = 0.0001)
    {
        boundaryPoint_offset = offset;

        for(int i = 0; i < env->trapezoid_edges.size(); i++)
        {
            float phi0 = array[i][0].phi;

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

class Event
{
public:

    Environment env;
    vector<Point> list_of_Points;

    Event(Environment envI, vector<Point> listP = {})
    {
        env = envI;
        list_of_Points = listP;
    }
};

class Line
{
public:
    Environment env;
    float slope;
    vector<float> points;

    Line(Environment envI, float start, float slopeI)
    {
        env = envI;
        slope = slopeI;

        points.push_back(start);

        for(int i = 0; i < env.radii.size(); i++)
        {
            points.push_back((1 / slope) * env.radii[i] + start);
        }
    }
};

class LineGenerator
{
public:
    Environment env;
    float start;
    float slope_ll;
    float slope_ul;

    LineGenerator(Environment envI, float startI)
    {
        env = envI;
        start = startI;

        if(start < (-env.beam_axis_lim) || start > env.beam_axis_lim)
        {
            throw "Start points is out of range.";
        }

        float max_height = env.radii[env.radii.size() - 1];

        slope_ll = max_height / (-1 * env.top_layer_lim - start);
        slope_ul = max_height / (env.top_layer_lim - start);
    }

    vector<Line> generateEvenGrid(int n = 100)
    {
        float Rcoor = env.radii[env.radii.size() - 1];

        float stepVal = (env.top_layer_lim) * 2 / (n - 1);

        vector<float> Zcoor;

        for(float i = 0; i < 100; i++)
        {
            Zcoor.push_back((-1 * env.top_layer_lim) + i * stepVal);
        }

        vector<Line> lines;

        for(int i = 0; i < Zcoor.size(); i++)
        {
            lines.push_back(Line(env, start, Rcoor / (Zcoor[i] - start)));
        }

        return lines;
    }
};

class parallelogram
{
public:
    int layer_num;
    float pSlope;

    float shadow_bottomL_jR;
    float shadow_bottomR_jR;
    float shadow_bottomL_jL;
    float shadow_bottomR_jL;

    float z1_min;
    float z1_max;

    parallelogram(int layer_numI, float z1_minI, float z1_maxI, float shadow_bottomL_jRI, float shadow_bottomR_jRI, float shadow_bottomL_jLI, float shadow_bottomR_jLI, float pSlopeI)
    {
        layer_num = layer_numI;
        pSlope = pSlopeI;

        shadow_bottomL_jR = shadow_bottomL_jRI;
        shadow_bottomR_jR = shadow_bottomR_jRI;

        shadow_bottomL_jL = shadow_bottomL_jLI;
        shadow_bottomR_jL = shadow_bottomR_jLI;

        z1_min = z1_minI;
        z1_max = z1_maxI;

        if((z1_min > 22.0) || (z1_max < -22.0))
        {
            //cout << top_layer_zmin << " " << top_layer_zmax << endl;
        }
    }
};

class parallelogram_v1
{
public:
    int layer_num;
    float pSlope;

    float shadow_topR_jL;
    float shadow_topR_jR;

    float shadow_topL_jL;
    float shadow_topL_jR;

    float top_layer_zmin;
    float top_layer_zmax;

    parallelogram_v1(int layer_numI, float top_layer_zminI, float top_layer_zmaxI, float shadow_topR_jLI, float shadow_topR_jRI, float pSlopeI)
    {
        layer_num = layer_numI;
        pSlope = pSlopeI;

        shadow_topR_jL = shadow_topR_jLI;
        shadow_topR_jR = shadow_topR_jRI;

        float delta_ztop = top_layer_zmax - top_layer_zmin;
        float delta_z0 = delta_ztop * pSlope;

        shadow_topL_jL = shadow_topR_jLI - delta_z0;
        shadow_topL_jR = shadow_topR_jRI - delta_z0;

        top_layer_zmin = top_layer_zminI;
        top_layer_zmax = top_layer_zmaxI;

        if(top_layer_zmax > 100 or top_layer_zmin < -100)
        {
            cout << top_layer_zmin << " " << top_layer_zmax << endl;
        }
    }
};

class wedgeSuperPoint
{
public:
    vector<Point> points;
    vector<float> z_values;
    float min;
    float max;

    wedgeSuperPoint(vector<Point> pointsI)
    {
        if(pointsI.size() != 16)
        {
            if((pointsI.size() != 32) && (pointsI.size() != 31))
            {
                throw "This patch does not have 16 or 32/31 points in each layer";
            }
        }

        vector<float> z_list;

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
    bool contains(float p)
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

    bool operator==(wedgeSuperPoint other)
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
    float left_end_lambdaZ;
    float right_end_lambdaZ;
    float apexZ0;

    float shadow_fromTopToInnermost_topL_jL;
    float shadow_fromTopToInnermost_topL_jR;
    float shadow_fromTopToInnermost_topR_jL;
    float shadow_fromTopToInnermost_topR_jR;

    vector<float> a_corner;
    vector<float> b_corner;
    vector<float> c_corner;
    vector<float> d_corner;

    vector<wedgeSuperPoint> superpoints;

    bool flatBottom;
    bool flatTop;

    bool squareAcceptance;
    bool triangleAcceptance;

    vector<parallelogram> parallelograms;
    vector<parallelogram_v1> parallelograms_v1;

    wedgePatch(Environment envI, vector<wedgeSuperPoint> superpointsI, float apexZ0I)
    {
        env = envI;
        end_layer = -1;
        left_end_layer = -1;
        right_end_layer = -1;
        left_end_lambdaZ = 0;
        right_end_lambdaZ = 0;
        apexZ0 = apexZ0I;

        shadow_fromTopToInnermost_topL_jL = 0;
        shadow_fromTopToInnermost_topL_jR = 0;
        shadow_fromTopToInnermost_topR_jL = 0;
        shadow_fromTopToInnermost_topR_jR = 0;

        if(superpointsI.size() != env.num_layers)
        {
            throw "The patch layers does not match environment layers. ";
        }

        superpoints = superpointsI;

        getParallelograms();
        getParallelograms_v1();
        get_acceptanceCorners();
        get_end_layer();

    }

    float straightLineProjectorFromLayerIJtoK(float z_i, float z_j, float i, float j, float k)
    {
        float radius_i = 0;
        float radius_j = 0;
        float radius_k = 0;

        if(i == 0)
        {
            radius_i = 0;
        }
        else
        {
            radius_i = env.radii[i - 1];
        }
        if(j == 0)
        {
            radius_j = 0;
        }
        else
        {
            radius_j = env.radii[j - 1];
        }
        if(k == 0)
        {
            radius_k = 0;
        }
        else
        {
            radius_k = env.radii[k - 1];
        }

        float radii_leverArm = (radius_k - radius_i) / (radius_j - radius_i);

        return z_i + (z_j - z_i) * radii_leverArm;
    }

    void getParallelograms()
    {
        vector<parallelogram> parallelogramsI;

        float z1_min = max(superpoints[0].min, -1 * env.trapezoid_edges[0]);
        float z1_max = min(superpoints[0].max, env.trapezoid_edges[0]);

        if (z1_min > z1_max)
        {
            z1_min = env.trapezoid_edges[0] + 1;
            z1_max = z1_min;
        }

        for(int i = 1; i < superpoints.size(); i++)
        {
            int j = i + 1;

            float z_j_min = superpoints[i].min;
            float z_j_max = superpoints[i].max;

            float a = straightLineProjectorFromLayerIJtoK(z1_min, z_j_max, 1, j, env.num_layers);
            float b = straightLineProjectorFromLayerIJtoK(z1_max, z_j_max, 1, j, env.num_layers);
            float c = straightLineProjectorFromLayerIJtoK(z1_min, z_j_min, 1, j, env.num_layers);
            float d = straightLineProjectorFromLayerIJtoK(z1_max, z_j_min, 1, j, env.num_layers);

            float pSlope;

            if (j != env.num_layers)
            {
                pSlope = env.parallelogramSlopes[j - 1];
            }
            else
            {
                pSlope = INT_MAX;
            }

            parallelogram Parallelogram = parallelogram(j, z1_min, z1_max, a, b, c, d, pSlope);
            parallelogramsI.push_back(Parallelogram);
        }

        parallelograms = parallelogramsI;
    }

    void getShadows(float zTopMin, float zTopMax)
    {
        float zTop_min = max(zTopMin, -env.trapezoid_edges[env.num_layers-1]);
        float zTop_max = min(zTopMax, env.trapezoid_edges[env.num_layers-1]);
        vector<float> topL_jL;
        vector<float> topL_jR;
        vector<float> topR_jL;
        vector<float> topR_jR;

        for(int i = 0; i < superpoints.size() - 1; i++)
        {
            int j = i + 1;
            float z_j_min = superpoints[i].min;
            float z_j_max = superpoints[i].max;

            topL_jL.push_back(straightLineProjectorFromLayerIJtoK(zTop_min, z_j_min, env.num_layers, j, 1));
            topL_jR.push_back(straightLineProjectorFromLayerIJtoK(zTop_min, z_j_max, env.num_layers, j, 1));
            topR_jL.push_back(straightLineProjectorFromLayerIJtoK(zTop_max, z_j_min, env.num_layers, j, 1));
            topR_jR.push_back(straightLineProjectorFromLayerIJtoK(zTop_max, z_j_max, env.num_layers, j, 1));
        }

        shadow_fromTopToInnermost_topL_jL = *max_element(topL_jL.begin(), topL_jL.end());
        shadow_fromTopToInnermost_topL_jR = *min_element(topL_jR.begin(), topL_jR.end());
        shadow_fromTopToInnermost_topR_jL = *max_element(topR_jL.begin(), topR_jL.end());
        shadow_fromTopToInnermost_topR_jR = *min_element(topR_jR.begin(), topR_jR.end());
    }

    void getParallelograms_v1()
    {
        vector<parallelogram_v1> parallelogramsI;

        float top_layer_zmin = max(superpoints[superpoints.size() - 1].min, -1 * env.top_layer_lim);
        float top_layer_zmax = min(superpoints[superpoints.size() - 1].max, env.top_layer_lim);

        if (top_layer_zmin > top_layer_zmax)
        {
            top_layer_zmin = env.top_layer_lim + 1;
            top_layer_zmax = top_layer_zmin;
        }

        for(int i = 0; i < superpoints.size() - 1; i++)
        {
            int j = i + 1;

            float z_j_min = superpoints[i].min;
            float z_j_max = superpoints[i].max;

            float a = straightLineProjector(top_layer_zmax, z_j_min, j);
            float b = straightLineProjector(top_layer_zmax, z_j_max, j);

            float pSlope = env.parallelogramSlopes[j - 1];

            parallelogram_v1 Parallelogram = parallelogram_v1(j, top_layer_zmin, top_layer_zmax, a, b, pSlope);
            parallelogramsI.push_back(Parallelogram);
        }

        parallelograms_v1 = parallelogramsI;
    }

    float straightLineProjector(float z_top, float z_j, int j)
    {
        float radii_leverArm = env.radii_leverArm[j - 1];
        return z_top - (z_top - z_j) * radii_leverArm;
    }

    void get_acceptanceCorners()
    {
        squareAcceptance = true;
        flatTop = true;
        flatBottom = true;
        triangleAcceptance = false;

        vector<float> a_corner_list;
        vector<float> b_corner_list;
        vector<float> c_corner_list;
        vector<float> d_corner_list;

        for(int i = 0; i < parallelograms.size(); i++)
        {
            a_corner_list.push_back(parallelograms[i].shadow_bottomL_jR);
            b_corner_list.push_back(parallelograms[i].shadow_bottomR_jR);
            c_corner_list.push_back(parallelograms[i].shadow_bottomL_jL);
            d_corner_list.push_back(parallelograms[i].shadow_bottomR_jL);
        }

        a_corner.push_back(parallelograms[0].z1_min);
        a_corner.push_back(*min_element(a_corner_list.begin(), a_corner_list.end()));
        b_corner.push_back(parallelograms[0].z1_max);
        b_corner.push_back(*min_element(b_corner_list.begin(), b_corner_list.end()));
        c_corner.push_back(parallelograms[0].z1_min);
        c_corner.push_back(*max_element(c_corner_list.begin(), c_corner_list.end()));
        d_corner.push_back(parallelograms[0].z1_max);
        d_corner.push_back(*max_element(d_corner_list.begin(), d_corner_list.end()));

        if (*min_element(a_corner_list.begin(), a_corner_list.end()) != a_corner_list[env.num_layers - 2])
        {
            squareAcceptance = false;
            flatTop = false;
        }
        if (*min_element(b_corner_list.begin(), b_corner_list.end()) != b_corner_list[env.num_layers - 2])
        {
            squareAcceptance = false;
            flatTop = false;
        }
        if (*max_element(c_corner_list.begin(), c_corner_list.end()) != c_corner_list[env.num_layers - 2])
        {
            squareAcceptance = false;
            flatBottom = false;
        }
        if (*max_element(d_corner_list.begin(), d_corner_list.end()) != d_corner_list[env.num_layers - 2])
        {
            squareAcceptance = false;
            flatBottom = false;
        }


        if (c_corner[1] > a_corner[1])
        {
            triangleAcceptance = true;
            c_corner[1] = b_corner[1];
            a_corner[1] = b_corner[1];
        }

        if (b_corner[1] < d_corner[1])
        {
            triangleAcceptance = true;
            b_corner[1] = c_corner[1];
            d_corner[1] = c_corner[1];
        }
    }

    void get_end_layer()
    {
        vector<float> lambdaZ_left_list;
        vector<float> lambdaZ_right_list;

        for(int i = 0; i < env.num_layers; i++)
        {
            lambdaZ_left_list.push_back((superpoints[i].min - apexZ0) / env.radii[i]);
            lambdaZ_right_list.push_back((superpoints[i].max - apexZ0)/ env.radii[i]);
        }

        float lambdaZLeftMax = -1 * INT_MAX + 2;
        float lambdaZRightMin = INT_MAX - 2;

        for(int i = 0; i < lambdaZ_left_list.size(); i++)
        {
            if(lambdaZ_left_list[i] > lambdaZLeftMax)
            {
                left_end_layer = i;
                lambdaZLeftMax = lambdaZ_left_list[i];
            }
        }

        for(int i = 0; i < lambdaZ_right_list.size(); i++)
        {
            if(lambdaZ_right_list[i] < lambdaZRightMin)
            {
                right_end_layer = i;
                lambdaZRightMin = lambdaZ_right_list[i];
            }
        }

        left_end_lambdaZ = *max_element(lambdaZ_left_list.begin(),lambdaZ_left_list.end());
        right_end_lambdaZ = *min_element(lambdaZ_right_list.begin(),lambdaZ_right_list.end());
    }
};

class wedgeCover
{
public:
	int n_patches; 
	vector<wedgePatch> patches;
	Environment env; 
	DataSet * data;
	vector<Line> fitting_lines; //Check the data types for all of the 4 lists here
    vector<wedgeSuperPoint> superPoints;
    vector<wedgePatch> all_patches;
    vector<bool> real_patch_list; 

	wedgeCover(Environment envI, DataSet& dataI)
	{
		n_patches = 0; 
		env = envI; 
		data = &dataI;
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

    int get_index_from_z(int layer, float z_value, string alignment = "closest")
    {
        vector<float> layer_data;

        for(int i = 0; i < data->array[layer].size(); i++)
        {
            layer_data.push_back(data->array[layer][i].z);
        }
        double minVal = 1000000;
        int index = 0;

        for(int i = 0; i < layer_data.size(); i++)
        {
            if(abs(layer_data[i] - z_value) < abs(minVal))
            {
                minVal = abs(layer_data[i] - z_value);
                index = i;
            }
        }

        if (alignment == "closest")
        {
            return index;
        }

        if (alignment == "above")
        {
            if (layer_data[index] > z_value)
            {
                return index;
            }
            else
            {
                return index + 1;
            }
        }

        if (alignment == "below")
        {
            if (layer_data[index] < z_value)
            {
                return index;
            }
            else
            {
                return index - 1;
            }
        }
        return 0; //need return type
    }

    void delete_patch(int index)
    {
        patches.erase(patches.begin() + index);
        real_patch_list[index] = false;
    }

    void solve(string lining = "makePatches_Projective", float apexZ0 = 0, int ppl = 16, int nlines = 100, bool leftRight = true, bool show = true)
    {
        for(int i = 0; i < env.num_layers; i++)
        {
            bool foundIdentical = false;
            bool firstTime = true;

            while(foundIdentical || firstTime)
            {
                foundIdentical = false;
                for(int x = 0; x < data->array[i].size() - 1; x++)
                {
                    if(data->array[i][x].z == data->array[i][x + 1].z)
                    {
                        data->array[i][x + 1].z += 0.00001;
                        foundIdentical = true;
                    }
                }

                firstTime = false;
                if(foundIdentical)
                {
                    sort(data->array[i].begin(), data->array[i].end(), [](Point &a, Point &b){ return a.z < b.z; });
                }
            }
        }

        if(show)
        {
            vector<Line> fitting_linesI;

            //implement check if apexZ0 is a single number;

            LineGenerator lGen = LineGenerator(env, apexZ0);
            vector<Line> temp = lGen.generateEvenGrid(nlines);
            fitting_linesI.insert(fitting_linesI.end(), temp.begin(), temp.end() );

            fitting_lines = fitting_linesI;
        }

        if(lining == "makePatches_ShadowQuilt_fromEdges")
        {
            try
            {
                makePatches_ShadowQuilt_fromEdges(apexZ0 = apexZ0, 1,ppl = ppl, leftRight = leftRight);
            }
            catch (string str)
            {
                makePatches_ShadowQuilt_fromEdges(apexZ0 = apexZ0, 1,ppl = ppl, leftRight = leftRight);
            }

            return;
        }
    }

    void makePatches_ShadowQuilt_fromEdges(float apexZ0 = 0, int stop = 1, int ppl = 16, bool leftRight = true)
    {
        bool fix42 = true;
        apexZ0 = env.trapezoid_edges[0];
        float saved_apexZ0;

        while(apexZ0 > -1 * env.trapezoid_edges[0])
        {
            float z_top_min = -1 * env.top_layer_lim;

            float complementary_apexZ0 = 0;
            int first_row_count = 0;
            float c_corner = LONG_MAX;

            float z_top_max = env.top_layer_lim + env.boundaryPoint_offset;

            if(patches.size() > 0)
            {
                z_top_max = min(z_top_max, patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(-1 * env.beam_axis_lim, apexZ0, 0, 1, env.num_layers));
            }

            int nPatchesInColumn = 0;
            float projectionOfCornerToBeam = 0;

            while((c_corner > -1 * env.trapezoid_edges[env.num_layers - 1]) && (projectionOfCornerToBeam < env.beam_axis_lim))
            {
                nPatchesInColumn++;
                cout << apexZ0 << " " << ppl << " " << z_top_max << " " << leftRight << endl;
                makePatch_alignedToLine(apexZ0, z_top_max, ppl = ppl, false);
                cout << "top layer from " << patches[patches.size() - 1].superpoints[env.num_layers - 1].max << " to " << patches[patches.size() - 1].superpoints[env.num_layers - 1].min << " z_top_max: " << z_top_max << endl;
                cout << "original: [" << patches[patches.size() - 1].a_corner[0] << ", " << patches[patches.size() - 1].a_corner[1] << "] for patch " << patches.size() << endl;
                cout << "original: [" << patches[patches.size() - 1].b_corner[0] << ", " << patches[patches.size() - 1].b_corner[1] << "]" << endl;
                cout << "original: [" << patches[patches.size() - 1].c_corner[0] << ", " << patches[patches.size() - 1].c_corner[1] << "]" << endl;
                cout << "original: [" << patches[patches.size() - 1].d_corner[0] << ", " << patches[patches.size() - 1].d_corner[1] << "]" << endl;

                for(int i = 1; i < patches[patches.size() - 1].superpoints.size() - 1; i++)
                {
                    int j = i + 1;
                    cout << j << " superpoint: " << patches[patches.size() - 1].superpoints[i].min << " " << patches[patches.size() - 1].superpoints[i].max <<
                            " shadowTop from L1Max: " << patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(patches[patches.size() - 1].superpoints[0].max, patches[patches.size() - 1].superpoints[i].min, 1, j, env.num_layers) << " " <<
                            patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(patches[patches.size() - 1].superpoints[0].max, patches[patches.size() - 1].superpoints[i].max, 1, j, env.num_layers) <<
                            " from L1 min: " << patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(patches[patches.size() - 1].superpoints[0].min, patches[patches.size() - 1].superpoints[i].min, 1, j, env.num_layers) << " " <<
                            patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(patches[patches.size() - 1].superpoints[0].min, patches[patches.size() - 1].superpoints[i].max, 1, j, env.num_layers) << endl;
                }

                float original_c = patches[patches.size() - 1].c_corner[1];
                float original_d = patches[patches.size() - 1].d_corner[1];

                c_corner = original_c;

                bool repeat_patch = false;
                bool repeat_original = false;

                if(patches.size() > 2)
                {
                    repeat_original = (patches[patches.size() - 1].superpoints[env.num_layers - 1] == patches[patches.size() - 3].superpoints[env.num_layers - 1]) &&
                            (patches[patches.size() - 1].superpoints[0] == patches[patches.size() - 3].superpoints[0]) &&
                            (patches[patches.size() - 1].superpoints[1] == patches[patches.size() - 3].superpoints[1]) &&
                            (patches[patches.size() - 1].superpoints[2] == patches[patches.size() - 3].superpoints[2]) &&
                            (patches[patches.size() - 1].superpoints[3] == patches[patches.size() - 3].superpoints[3]);
                }

                float seed_apexZ0 = apexZ0;
                projectionOfCornerToBeam = patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(patches[patches.size() - 1].c_corner[1], patches[patches.size() - 1].c_corner[0], env.num_layers, 1, 0);
                bool squarePatch_alternate1 = (patches[patches.size() - 1].a_corner[1] > z_top_max) && (patches[patches.size() - 1].b_corner[1] > z_top_max) && (patches[patches.size() - 1].flatBottom);
                bool squarePatch_alternate2 = (patches[patches.size() - 1].a_corner[1] > z_top_max) && (patches[patches.size() - 1].flatBottom);

                bool notChoppedPatch = (patches[patches.size() - 1].squareAcceptance) || squarePatch_alternate1 || squarePatch_alternate2;
                bool madeComplementaryPatch = false;

                int nPatchesAtOriginal = patches.size();

                cout << "squareAcceptance: " << patches[patches.size() - 1].squareAcceptance << " triangleAcceptance: " << patches[patches.size() - 1].triangleAcceptance << " projectionOfCornerToBeam: " << projectionOfCornerToBeam << " notChoppedPatch " << notChoppedPatch << endl;
                if (!(notChoppedPatch) && (patches[patches.size() - 1].c_corner[1] > -1 * env.trapezoid_edges[env.num_layers - 1]) && (projectionOfCornerToBeam < env.beam_axis_lim))
                {
                    complementary_apexZ0 = patches[patches.size() - 1].superpoints[0].min;
                    if ((patches[patches.size() - 1].triangleAcceptance) && !(repeat_original))
                    {
                        z_top_min = patches[patches.size() - 1].d_corner[1];
                    }
                    else
                    {
                        cout << "z_top_min before: " << z_top_min << " superpoints[self.env.num_layers-1].min: " << patches[patches.size() - 1].superpoints[env.num_layers - 1].min << endl;
                        z_top_min = max(-1 * env.top_layer_lim, patches[patches.size() - 1].superpoints[env.num_layers - 1].min);
                    }

                    makePatch_alignedToLine(complementary_apexZ0, z_top_min, ppl, true);
                    madeComplementaryPatch = true;
                    cout << "complementary: [" << patches[patches.size() - 1].a_corner[0] << ", " << patches[patches.size() - 1].a_corner[1] << "] for z_top_min: " << z_top_min << endl;
                    cout << "complementary: [" << patches[patches.size() - 1].b_corner[0] << ", " << patches[patches.size() - 1].b_corner[1] << "] for patch  " << patches.size() << endl;
                    cout << "complementary: [" << patches[patches.size() - 1].c_corner[0] << ", " << patches[patches.size() - 1].c_corner[1] << endl;
                    cout << "complementary: [" << patches[patches.size() - 1].d_corner[0] << ", " << patches[patches.size() - 1].d_corner[1] << endl;

                    float complementary_a = patches[patches.size() - 1].a_corner[1];
                    float complementary_b = patches[patches.size() - 1].b_corner[1];

                    float white_space_height = max(original_c - complementary_a, original_d - complementary_b);
                    float previous_white_space_height = -1;
                    int counter = 0;
                    int counterUpshift = 0;
                    int current_z_top_index = -1;
                    double previous_z_top_min = -999;

                    while (!(white_space_height <= 0 && (previous_white_space_height >= 0)) && (abs(white_space_height) > 0.000001) && ((patches[patches.size() - 1].c_corner[1] > -1 * env.trapezoid_edges[env.num_layers - 1]) || (white_space_height > 0)) && (current_z_top_index < (int) (data->array[env.num_layers - 1].size() - 1)) && !(repeat_patch) && !(repeat_original))
                    {
                        cout << endl;
                        if(patches.size() > 2)
                        {
                            cout << "original c: " << original_c << " " << patches[patches.size() - 2].c_corner[1] << " || original d: " << original_d << " " << patches[patches.size() - 2].d_corner[1] << endl;
                        }
                        cout << "complementary_a: " << complementary_a << " " << patches[patches.size() - 1].a_corner[1] << " || complementary_b: " << complementary_b << " " << patches[patches.size() - 1].b_corner[1] << endl;
                        current_z_top_index = get_index_from_z(env.num_layers - 1, z_top_min);
                        cout << "current white_space_height: " << white_space_height << endl;
                        cout << "counter: " << counter << " counterUpshift: " << counterUpshift << endl;
                        cout << "orig_ztop: " << current_z_top_index << " orig_z_top_min: " << z_top_min << endl;

                        vector<float> current_z_i_index;
                        vector<float> new_z_i_index;

                        for(int i = 0; i < env.num_layers; i++)
                        {
                            current_z_i_index.push_back(get_index_from_z(i, patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(complementary_apexZ0, z_top_min, 1, env.num_layers, i + 1)));
                        }

                        if (z_top_min == previous_z_top_min)
                        {
                            current_z_top_index += 1;
                            for(int i = 0; i < current_z_i_index.size(); i++)
                            {
                                new_z_i_index.push_back(current_z_i_index[i] + 1);
                            }
                        }

                        previous_z_top_min = z_top_min;

                        if (white_space_height < 0)
                        {
                            counter += 1;
                            current_z_top_index -= 1;
                            for(int i = 0; i < current_z_i_index.size(); i++)
                            {
                                new_z_i_index.push_back(current_z_i_index[i] - 1);
                            }
                        }
                        else
                        {
                            counterUpshift += 1;
                            current_z_top_index += 1;
                            for(int i = 0; i < current_z_i_index.size(); i++)
                            {
                                new_z_i_index.push_back(current_z_i_index[i] + 1);
                            }
                        }

                        int x = (int) data->array[env.num_layers - 1].size() - 1;
                        current_z_top_index = min(current_z_top_index,(int) data->array[env.num_layers - 1].size() - 1);

                        for(int i = 0; i < new_z_i_index.size(); i++)
                        {
                            new_z_i_index[i] = min(new_z_i_index[i], (float) data->array[i].size() - 1);
                        }

                        for(int i = 0; i < new_z_i_index.size(); i++)
                        {
                            new_z_i_index[i] = max(new_z_i_index[i], (float) 0.0);
                        }

                        vector<float> new_z_i;

                        for(int i = 0; i < env.num_layers; i++)
                        {
                            new_z_i.push_back(data->array[i][new_z_i_index[i]].z);
                        }

                        vector<float> new_z_i_atTop;

                        for(int i = 1; i < env.num_layers; i++)
                        {
                            new_z_i_atTop.push_back(patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(complementary_apexZ0,new_z_i[i],1,i + 1,env.num_layers));
                        }

                        int layerWithSmallestShift  = 0;
                        float layerSMin = INT_MAX;

                        for(int i = 0; i < new_z_i_atTop.size(); i++)
                        {
                            if(abs(new_z_i_atTop[i] - previous_z_top_min) < layerSMin)
                            {
                                layerSMin = abs(new_z_i_atTop[i] - previous_z_top_min);
                                layerWithSmallestShift = i;
                            }
                        }

                        layerWithSmallestShift += 1;

                        for(int i = 0; i < env.num_layers - 1; i++)
                        {
                            cout << i + 1 << " new_z_i_atTop: " << new_z_i_atTop[i] << " shift_i_ztop: " << new_z_i_atTop[i] - previous_z_top_min << " layerWithSmallestShift: " << layerWithSmallestShift << endl;
                        }

                        z_top_min = data->array[env.num_layers - 1][current_z_top_index].z;
                        z_top_min = new_z_i_atTop[layerWithSmallestShift - 1];

                        if (abs(z_top_min-previous_z_top_min) < 0.000001)
                        {
                            z_top_min = data->array[env.num_layers - 1][current_z_top_index].z;
                        }

                        if (abs(z_top_min-previous_z_top_min) < 0.000001)
                        {
                            z_top_min = data->array[env.num_layers - 2][current_z_top_index].z;
                        }

                        if (abs(z_top_min-previous_z_top_min) < 0.000001)
                        {
                            z_top_min = data->array[env.num_layers - 3][current_z_top_index].z;
                        }

                        if (((z_top_min - previous_z_top_min) * (white_space_height)) < 0)
                        {
                            z_top_min = new_z_i_atTop[env.num_layers - 2];
                        }

                        cout << " new_def_z_top_min_diff: " << z_top_min - data->array[env.num_layers - 1][current_z_top_index].z << endl;

                        cout << " new_ztop_index: " << current_z_top_index << " new_z_i_index: " << new_z_i_index[0] << " " << new_z_i_index[1] << " " << new_z_i_index[2] << " " << new_z_i_index[3] << " " << new_z_i_index[4] << " new_z_top_min: " << z_top_min << " shift_ztop: " << z_top_min - previous_z_top_min << endl;

                        int nPatchesAtComplementary = patches.size();

                        if (nPatchesAtComplementary > nPatchesAtOriginal)
                        {
                            cout << "deleted complementary: " << patches[patches.size() - 1].a_corner[0] << " " << patches[patches.size() - 1].a_corner[1]  << " for patch" << patches.size() << endl;
                            cout << "deleted complementary: " << patches[patches.size() - 1].b_corner[0] << " " << patches[patches.size() - 1].b_corner[1] << endl;
                            cout << "deleted complementary: " << patches[patches.size() - 1].c_corner[0] << " " << patches[patches.size() - 1].c_corner[1] << endl;
                            cout << "deleted complementary: " << patches[patches.size() - 1].d_corner[0] << " " << patches[patches.size() - 1].d_corner[1] << endl;

                            delete_patch(patches.size() - 1); //write delete patch
                            n_patches -= 1;
                        }

                        makePatch_alignedToLine(complementary_apexZ0, z_top_min, ppl, true);

                        complementary_a = patches[patches.size() - 1].a_corner[1];
                        complementary_b = patches[patches.size() - 1].b_corner[1];

                        previous_white_space_height = white_space_height;

                        white_space_height = max(original_c - complementary_a, original_d - complementary_b);

                        cout << "complementary_a:" << complementary_a << " " << patches[patches.size() - 1].a_corner[1] << " || complementary_b:" << complementary_b << " " << patches[patches.size() - 1].b_corner[1] << " new z_top_min: " << z_top_min << endl;
                        cout << "new white_space_height: " << white_space_height << endl;
                        cout << "adjusted complementary: " << patches[patches.size() - 1].a_corner[0] << " " << patches[patches.size() - 1].a_corner[1] << " for z_top_min:" << z_top_min << endl;
                        cout << "adjusted complementary: " << patches[patches.size() - 1].b_corner[0] << " " << patches[patches.size() - 1].b_corner[1] << "for patch " << patches.size() << endl;
                        cout << "adjusted complementary: " << patches[patches.size() - 1].c_corner[0] << " " << patches[patches.size() - 1].c_corner[1] << endl;
                        cout << "adjusted complementary: " << patches[patches.size() - 1].d_corner[0] << " " << patches[patches.size() - 1].d_corner[1] << endl;

                        if ((n_patches > 3) && fix42)
                        {
                            if ((patches[patches.size() - 1].superpoints[env.num_layers - 1] == patches[patches.size() - 3].superpoints[env.num_layers - 1]) && (patches[patches.size() - 1].superpoints[0] == patches[patches.size() - 3].superpoints[0]) && (patches[patches.size() - 1].superpoints[1] == patches[patches.size() - 3].superpoints[1]) && (patches[patches.size() - 1].superpoints[2] == patches[patches.size() - 3].superpoints[2]) && (patches[patches.size() - 1].superpoints[3] == patches[patches.size() - 3].superpoints[3]))
                            {
                                repeat_patch = true;
                                cout << patches[patches.size() - 1].superpoints[env.num_layers - 1].min << " " << patches[patches.size() - 1].superpoints[env.num_layers - 1].max << " repeat_patch: " << repeat_patch << endl;
                                delete_patch(patches.size() - 1);
                                n_patches -= 1;
                                current_z_top_index -= 1;
                                z_top_min = data->array[env.num_layers - 1][current_z_top_index].z;
                                z_top_min = new_z_i_atTop[layerWithSmallestShift - 1];
                                makePatch_alignedToLine( complementary_apexZ0, z_top_min,  ppl, true);
                            }
                        }


                    }
                }

                c_corner = patches[patches.size() - 1].c_corner[1];

                projectionOfCornerToBeam = patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(c_corner,patches[patches.size() - 1].c_corner[0],env.num_layers,1,0);

                saved_apexZ0 = patches[patches.size() - 1].c_corner[0];

                if (madeComplementaryPatch)
                {
                    patches[patches.size() - 1].getShadows(z_top_min,z_top_max);
                    patches[patches.size() - 2].getShadows(z_top_min,z_top_max);

                    float original_topR_jL = patches[patches.size() - 2].shadow_fromTopToInnermost_topR_jL;
                    bool originalPartialTop = (original_topR_jL > complementary_apexZ0) && (original_topR_jL < apexZ0) && (abs(patches[patches.size() - 2].straightLineProjectorFromLayerIJtoK(original_topR_jL,z_top_max,1,env.num_layers,0)) < 20 * env.beam_axis_lim);
                    float original_topL_jL = patches[patches.size() - 2].shadow_fromTopToInnermost_topL_jL;
                    bool originalPartialBottom = (original_topL_jL > complementary_apexZ0) && (original_topL_jL < apexZ0) && (abs(patches[patches.size() - 2].straightLineProjectorFromLayerIJtoK(original_topL_jL,z_top_min,1,env.num_layers,0)) < 20 * env.beam_axis_lim);
                    float complementary_topR_jR = patches[patches.size() - 1].shadow_fromTopToInnermost_topR_jR;
                    bool complementaryPartialTop = (complementary_topR_jR > complementary_apexZ0) && (complementary_topR_jR < apexZ0) && (abs(patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(complementary_topR_jR,z_top_max,1,env.num_layers,0)) < 20 * env.beam_axis_lim);
                    float complementary_topL_jR = patches[patches.size() - 1].shadow_fromTopToInnermost_topL_jR;
                    bool complementaryPartialBottom = (complementary_topL_jR > complementary_apexZ0) and (complementary_topL_jR < apexZ0) and (abs(patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(complementary_topL_jR,z_top_min,1,env.num_layers,0)) < 20 * env.beam_axis_lim);

                    float horizontalShiftTop = original_topR_jL - complementary_topR_jR;
                    float horizontalShiftBottom = original_topL_jL - complementary_topL_jR;

                    float complementary_topR_jL = patches[patches.size() - 1].shadow_fromTopToInnermost_topR_jL;
                    float complementary_topL_jL = patches[patches.size() - 1].shadow_fromTopToInnermost_topL_jL;
                    float original_topR_jR = patches[patches.size() - 2].shadow_fromTopToInnermost_topR_jR;
                    float original_topL_jR = patches[patches.size() - 2].shadow_fromTopToInnermost_topL_jR;

                    float horizontalOverlapTop = max(complementary_topR_jL - original_topR_jL, complementary_topR_jR - original_topR_jR);
                    float horizontalOverlapBottom = max(complementary_topL_jL - original_topL_jL, complementary_topL_jR - original_topL_jR);

                    horizontalOverlapTop = -1;
                    horizontalOverlapBottom = -1;
                    float newGapTop = -0.000001;
                    float newGapBottom = -0.000001;

                    bool makeHorizontallyShiftedPatch = false;
                    float shifted_Align = apexZ0;
                    bool doShiftedPatch = true;

                    float newZtop = 0;

                    float z0_original_bCorner = patches[patches.size() - 2].straightLineProjectorFromLayerIJtoK(apexZ0,z_top_max,1,env.num_layers,0);
                    float z0_complementary_cCorner = patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(complementary_apexZ0,z_top_min,1,env.num_layers,0);
                    bool shiftOriginal = true;

                    if (z0_original_bCorner < 0)
                    {
                        shiftOriginal = false;
                        shifted_Align = complementary_apexZ0;
                    }

                    if (z0_complementary_cCorner > 0)
                    {
                        shiftOriginal = true;
                        shifted_Align = apexZ0;
                    }

                    if (horizontalShiftTop > 0 or horizontalShiftBottom > 0)
                    {
                        cout << "originalPartialTop: " << originalPartialTop << " complementaryPartialTop: " << complementaryPartialTop << " originalPartialBottom: " << originalPartialBottom << " complementaryPartialBottom: " << complementaryPartialBottom << " " << original_topR_jL << " " << original_topL_jL << " " << complementary_topR_jR << " " << complementary_topL_jR << " h orizontalOverlapTop: " << horizontalOverlapTop << " horizontalOverlapBottom: " << horizontalOverlapBottom << endl;
                    }

                    while (((horizontalShiftTop > 0 && originalPartialTop && complementaryPartialTop) || (horizontalShiftBottom > 0 && originalPartialBottom && complementaryPartialBottom)) && doShiftedPatch && (horizontalOverlapTop <= 0) && (horizontalOverlapBottom <= 0) && (newGapTop < 0 || newGapBottom < 0))
                    {
                        cout << "horizontalShifts: " << horizontalShiftTop << " " << horizontalShiftBottom << " shifted_Align: " << shifted_Align << endl;

                        newZtop = z_top_max;

                        if (shiftOriginal)
                        {
                            shifted_Align -= max(horizontalShiftTop, horizontalShiftBottom);
                        }
                        else
                        {
                            shifted_Align += max(horizontalShiftTop,horizontalShiftBottom);
                            newZtop = z_top_min;
                        }

                        if (makeHorizontallyShiftedPatch)
                        {
                            delete_patch(patches.size() - 1);
                            n_patches -= 1;
                        }

                        makePatch_alignedToLine(shifted_Align, newZtop, ppl, (not shiftOriginal));

                        patches[patches.size() - 1].getShadows(z_top_min,z_top_max);

                        if (shiftOriginal)
                        {
                            original_topR_jL = patches[patches.size() - 1].shadow_fromTopToInnermost_topR_jL;
                            original_topL_jL = patches[patches.size() - 1].shadow_fromTopToInnermost_topL_jL;
                            original_topR_jR = patches[patches.size() - 1].shadow_fromTopToInnermost_topR_jR;
                            original_topL_jR = patches[patches.size() - 1].shadow_fromTopToInnermost_topL_jR;
                        }
                        else
                        {
                            complementary_topR_jR = patches[patches.size() - 1].shadow_fromTopToInnermost_topR_jR;
                            complementary_topL_jR = patches[patches.size() - 1].shadow_fromTopToInnermost_topL_jR;
                            complementary_topR_jL = patches[patches.size() - 1].shadow_fromTopToInnermost_topR_jL;
                            complementary_topL_jL = patches[patches.size() - 1].shadow_fromTopToInnermost_topL_jL;
                        }

                        horizontalShiftTop = original_topR_jL - complementary_topR_jR;
                        horizontalShiftBottom = original_topL_jL - complementary_topL_jR;

                        if (shiftOriginal && patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(original_topR_jR,z_top_max,1,env.num_layers,0) < env.beam_axis_lim)
                        {
                            horizontalOverlapTop = max(complementary_topR_jL - original_topR_jL, complementary_topR_jR - original_topR_jR);
                            horizontalOverlapBottom = max(complementary_topL_jL - original_topL_jL, complementary_topL_jR - original_topL_jR);
                            cout << " horizontalOverlapTop:" << " " << horizontalOverlapTop << " horizontalOverlapBottom: " << horizontalOverlapBottom << endl;
                        }

                        cout << "original_topR_jL: " << original_topR_jL << " complementary_topR_jR " << complementary_topR_jR << " original_topL_jL " << original_topL_jL << " complementary_topL_jR " << complementary_topL_jR << " shiftOriginal " << shiftOriginal << endl;

                        makeHorizontallyShiftedPatch = true;

                        cout << "updated_horizontalShifts: " << horizontalShiftTop << " " << horizontalShiftBottom << " shifted_Align: " << shifted_Align << endl;
                    }

                    if (makeHorizontallyShiftedPatch)
                    {
                        if (((patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(shifted_Align,newZtop,1,env.num_layers,0) > env.beam_axis_lim)) and shiftOriginal)
                        {
                            if (patches.size() > 2)
                            {
                                delete_patch(patches.size() - 3);
                                n_patches -= 1;
                            }
                        }
                    }
                }

                z_top_max = c_corner;

                cout << "+++++++++++++++++++++++ c_corner: " << c_corner << endl;
            }

            apexZ0 = patches[patches.size() - 1].c_corner[0];
            apexZ0 = saved_apexZ0;
            cout << "'=======================================================  z1_Align: " << apexZ0 << endl;

        }
    }

	void makePatch_alignedToLine(float apexZ0 = 0, float z_top = -50, int ppl = 16, bool leftRight = true, bool float_middleLayers_ppl = false)
	{
		vector<wedgeSuperPoint> init_patch; 
		int original_ppl = ppl; 
		float alignmentAccuracy = 0.00001; 

		vector< vector<Point> > row_data = data->array;

		for(int i = 0; i < env.num_layers; i++)
		{
			float y = env.radii[i]; 

			vector<float> row_list; 

			for(int j = 0; j < row_data[i].size(); j++)
			{
				row_list.push_back(row_data[i][j].z); 
			} 

			float r_max = env.radii.back(); 

			float projectionToRow = (z_top - apexZ0) * (y - env.radii[0]) / (r_max - env.radii[0]) + apexZ0;

			int start_index = 0;
			float start_value = 1000000;

			for(int j = 0; j < row_list.size(); j++)
			{
				if(abs(row_list[j] - projectionToRow) < abs(start_value))
				{
					start_index = j; 
					start_value = row_list[j] - projectionToRow;
				}
			}

			int left_bound = 0;
            float lbVal = INT_MAX;
			int right_bound = 0;
            float rbVal = INT_MAX;

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

			if ((float_middleLayers_ppl == true) && (i != 0) && (i != env.num_layers - 1))
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
        float apexZ0 = 0;
        int ppl = 32;
        float z_top = 1;
        bool leftRight = false;

        makePatch_alignedToLine(apexZ0, z_top, ppl, leftRight);
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

                    set<float> rawRadii;

                    for(int i = 0; i < list_of_Points.size(); i++)
                    {
                        rawRadii.insert(list_of_Points[i].radius);
                    }

                    vector<float> radii(rawRadii.begin(), rawRadii.end());

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

class Tester
{
public:
    void wedge_test(string lining = "makePatches_Projective_center", float apexZ0 = 0, float z0_spacing = 0.5, int ppl = 16, float z0_luminousRegion = 15.0, vector<int> wedges = { 0, 128 }, int lines = 1000, string v = "v3", float top_layer_cutoff = 50.0, float accept_cutoff = 10.0, bool leftRightAlign = true, bool uniform_N_points = false, string acceptance_method = "Analytic", bool show_acceptance_of_cover = false, bool movie = false, bool savefig = false, int figSizeScale = 6, int movieFigSizeScale = 3)
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

        vector<float> zInnerLayer;

        for(float i = -22; i < 22 + z0_spacing; i += z0_spacing)
        {
            zInnerLayer.push_back(i);
        }

        vector<float> z0Array;

        for(float i = -1 * z0_luminousRegion; i < z0_luminousRegion + z0_spacing; i += z0_spacing)
        {
            z0Array.push_back(i);
        }

        vector< vector<int> > mean_list;

        for(int i = 0; i < (wedges[1] - wedges[0]); i++)
        {
            vector<int> vect(z0Array.size(), 0);
            mean_list.push_back(vect);
        }

        vector<float> z0Imperfect; //check data types
        vector<float> z0OverEfficiency;

        vector<Event> all_data = FileReader::readFile("wedgeData_" + v + "_128.txt", wedges[1]);

        int ik = 0;

        for(int k = wedges[0]; k < wedges[1]; k++)
        {
            cout << "wedge: " << k << endl;

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

            cover.solve(lining, apexZ0, ppl, 100, leftRightAlign, false);

            num_covers.push_back(cover.n_patches);
            num_all_patches.push_back(cover.all_patches.size());

            ofstream myfile;
            myfile.open ("cppOutput.txt", ios::out | ios::trunc);

            for(int i = 0; i < cover.patches.size(); i++)
            {
                myfile << "Patch " << endl;
                for(int j = 0; j < cover.patches[i].superpoints.size(); j++)
                {
                    myfile << "Superpoint " << endl;
                    for(int r = 0; r < cover.patches[i].superpoints[j].points.size(); r++)
                    {
                        Point currentPt = cover.patches[i].superpoints[j].points[r];
                        myfile << currentPt.layer_num << " " << currentPt.phi << " " << currentPt.radius << " " << currentPt.z << endl;
                    }
                }
            }

            myfile << fixed;
            myfile.precision(5);

            for(int i = 0; i < cover.patches.size(); i++)
            {
                myfile << "['" << cover.patches[i].a_corner[0] << "', '" << cover.patches[i].a_corner[1] << "']" << endl;
                myfile << "['" << cover.patches[i].b_corner[0] << "', '" << cover.patches[i].b_corner[1] << "']" << endl;
                myfile << "['" << cover.patches[i].c_corner[0] << "', '" << cover.patches[i].c_corner[1] << "']" << endl;
                myfile << "['" << cover.patches[i].d_corner[0] << "', '" << cover.patches[i].d_corner[1] << "']" << endl;
                myfile << endl;
            }

            myfile.close();

            ik++;
        }
    }
};

int main()
{
    /*
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
    */

    string filepath = "wedgeData_v3_128.txt";
    vector<Event> events = FileReader::readFile(filepath);
    Environment env = events[0].env;
    vector<Point> points = events[0].list_of_Points;

    DataSet ds = DataSet(env);

    ds.importData(points);
    ds.addBoundaryPoint();

    wedgeCover cov = wedgeCover(env, ds);
    Tester test;
    test.wedge_test("makePatches_ShadowQuilt_fromEdges", 0, 0.5, 16, 15.0, {129, 130}, 1000, "v3", 50, 15.0, false, false, "Analytic", false, true, false, 6, 3);
};
