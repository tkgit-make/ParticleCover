#define MAIN_C
#include "header.h"

int main()
{
    int wedgesToTest[] = {0, 1};
    int wedge_count = 2;

    wedge_test(0, 0.5, 16, 15.0, wedgesToTest, wedge_count, 1000, "v3", 50, 15.0, false, false, "Analytic", false, false, false, 6, 3);

    return 0;
}
