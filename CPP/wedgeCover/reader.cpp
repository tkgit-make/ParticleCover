//
// Created by Tanish Kumar on 1/9/24.
//

#include <string>
#include <iostream>
#include <regex>

#include "data_structs.cpp"

using namespace std;

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