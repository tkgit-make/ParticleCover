#define MAIN_C
#include "header.h"

void ProcessEvent()
{
	printf("New event loaded, %d Points\n",G_event.count);
	G_stats.count_events++;
	G_stats.count_points += G_event.count;
}

int main() {
    int wedgesToTest[] = {24, 25}; 
    int wedge_count = 2;

    wedge_test("makePatches_ShadowQuilt_fromEdges", 0, 0.5, 16, 15.0, wedgesToTest, wedge_count, 1000, "v3", 50, 15.0, false, false, "Analytic", false, false, false, 6, 3);

    return 0;
}
