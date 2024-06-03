// #define MAIN_C
#include "header.h"

int main()
{
    int wedgesToTest[] = {0, 1};

    wedge_test(0, 0.025, 16, 15.0, wedgesToTest, 2, 1000, 50, 15.0);

    //printf("Size of float in C: %zu bytes\n", sizeof(float));
    
    return 0;
}
