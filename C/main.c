#define MAIN_C
#include "header.h"

int main()
{
    int wedgesToTest[] = {0, 6400};

    wedge_test(0, 16, wedgesToTest);

    //printf("Size of float in C: %zu bytes\n", sizeof(float));
    
    return 0;
}
