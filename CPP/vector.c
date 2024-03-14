#define VECTOR_C
#include "header.h"

typedef struct
{
	int objSize;
	int allocCount;
	int growBy;
	int usedCount;
	void *data;
} Vector;

Vector* Vector_create(int objSize, int intialSize, int growBy)
{
	Vector* v = (Vector*)calloc(sizeof(Vector), 1);
	v->objSize = objSize;
	v->allocCount = intialSize;
	v->growBy = growBy;
	v->usedCount = 0;
	v->data = malloc(v->allocCount * v->objSize); //allocates enough memory to store allocCount elements, each of size objSize
	return v; //returning a reference (pointer) to v
}

void* Vector_getitem(Vector *v, int n)
{
	return (n < v->usedCount) ? v->data + (n * v->objSize) : 0; 
	/* longhand
	if (n < v->usedCount) {
		return v->data + (n * v->objSize);
		//data is a pointer to the first element in the vector
		//n*objSize is the number of bytes to skip
		//pointer arithmetic 
	} else {
		return 0;
	}
	*/
}

//newitem at end of vector
void* Vector_newitem(Vector *v)
{
	//if the vector is full
	if (v->usedCount == v->allocCount)
	{
		v->allocCount = v->allocCount + v->growBy; //linear growth
		v->data = realloc(v->data, v->allocCount * v->objSize); //accommodate the new size
		assert(v->data); //check if data == null. this would imply failure
	}

	return Vector_getitem(v, v->usedCount++);
	//post-increment operator. Calls (v, v->usedCount) and then increments usedCount
	//returns the pointer to the desired spot in memory, and outside, the caller must set the value at that location
}

// Create vector support functions around a structural type
//macros for type-specific vector operations for any data type T
#define CREATE_VECTOR_OF_T(T)                                                                                                \
	typedef Vector *T##Vector; /*e.g. int vector, alias for pointer to a vector that holds integers*/                                                                                          \
	T##Vector VectorOf_##T(int initialsize, int growby) { return (T##Vector)Vector_create(sizeof(T), initialsize, growby); } /*e.g. VectorOf_int method*/\
	T *T##Vector_newitem(T##Vector v) { return (T *)Vector_newitem(v); } /*e.g. call intVector_newitem method with parameter of an intVector*/                                                   \
	T *T##Vector_getitem(T##Vector v, int n) { return (T *)Vector_getitem(v, n); } /*e.g. call intVector_getItem method with parameters*/
