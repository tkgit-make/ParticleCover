#define MAIN_C
#include "header.h"

int main()
{
    int wedgesToTest[] = {0, 128};

    wedge_test(0, 0.5, 16, 15.0, wedgesToTest, 2, 1000, 50, 15.0);

    return 0;
}
