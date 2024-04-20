#include "header.h"

void initLine(Line *line, Environment *envI, float start, float slopeI)
{
    line->env = envI; // line->env points to envI
    line->slope = slopeI;

    // adding the start point
    line->points[0] = start;
    line->num_points = 1;

    for (index_type i = 0; i < envI->num_layers; i++)
    {
        if (line->num_points < MAX_POINTS_IN_LINE)
        {
            line->points[line->num_points++] = (1.0f / slopeI) * envI->radii[i] + start;
        }
    }
}