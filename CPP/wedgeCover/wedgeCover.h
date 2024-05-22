//
// Created by Tanish Kumar on 1/8/24.
//


using namespace std;
#include <vector>

class Point {
public:
    Point();
};

class Environment {
public:
    Environment();
};

namespace N
{
    #ifndef WEDGECOVER_WEDGECOVER_H
    #define WEDGECOVER_WEDGECOVER_H
    #endif //WEDGECOVER_WEDGECOVER_H



    class DataSet {
    public:
        DataSet();

        DataSet(Environment env);

        void importData(vector<Point> data_array);

        void addBoundaryPoint(double offset = 0.0001);
    };

    class Event {
    public:
        Event();
    };



}