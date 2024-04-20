#include "header.h"

void initLineGenerator(LineGenerator *lg, Environment *envI, float startI)
{
    lg->env = envI;
    lg->start = startI;

    if (startI < (-envI->beam_axis_lim) || startI > envI->beam_axis_lim)
    {
        printf("Start point is out of range.");
        exit(4);
    }

    float max_height = envI->radii[envI->num_layers - 1];
    lg->slope_ll = max_height / (-1 * envI->top_layer_lim - startI);
    lg->slope_ul = max_height / (envI->top_layer_lim - startI);
}

// editing lines through the pointer as opposed to returning a vector of lines.
void generateEvenGrid(LineGenerator *lg, Line *lines, int n)
{
    if (n > MAX_POINTS_IN_LINE)
    {
        printf("Number of lines exceeds maximum allowed.");
        exit(5);
    }

    float Rcoor = lg->env->radii[lg->env->num_layers - 1];
    float stepVal = (lg->env->top_layer_lim) * 2 / (n - 1);

    for (index_type i = 0; i < n; i++)
    {
        float z = (-1 * lg->env->top_layer_lim) + i * stepVal; // an entry in Zcoor in C++
        float slope = Rcoor / (z - lg->start);
        initLine(&lines[i], lg->env, lg->start, slope); // equivalent to lines.push_back(Line(env, start, Rcoor / (Zcoor[i] - start)));.
    }
}