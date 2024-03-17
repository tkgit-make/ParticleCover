#ifndef POINT_C
    #include "point.c"
#endif 

#ifndef ENVIRONMENT_C
    #include "environment.c"
#endif 


CREATE_VECTOR_OF_T(int)
//CREATE_VECTOR_OF_T(Point)
CREATE_VECTOR_OF_T(PointVector)


typedef struct {
    Environment* env;
    PointVectorVector* array;  // A vector of PointVector
    intVector* n_points;
    int total_points;
    float boundaryPoint_offset;
} DataSet;

//c doesn't support overloading
void initDataSetBase(DataSet* ds) {
    ds->total_points = 0;
}

void initDataSetExtra(DataSet* ds, Environment* envI) {
    ds->total_points = 0;
    //lhs is pointer, consistent with rhs
    ds->env = envI;
    ds->array = VectorOf_PointVector(5, 64); //loops runs 5 times
    ds->n_points = VectorOf_int(envI->num_layers, 64); //initial size of num_layers based on the for loop below

    //adding 5 empty PointVectors to the vector
    for (int i = 0; i < 5; i++) {
        PointVector* vect = VectorOf_Point(64, 64);
        //points to the address where the pointer to the new PointVector will be stored in the PointVectorVector
        PointVector** next = PointVectorVector_newitem(ds->array);
        //dereferencing only once, changing the PointVector pointer to point to the next element. 
        *next = vect;
    }

    for (int i = 0; i < ds->env->num_layers; i++) {
        int* newItem = intVector_newitem(ds->n_points);
        *newItem = 0;
    }

    //continue writing importData and addBoundaryPoint
}
