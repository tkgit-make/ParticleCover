# Debugging Changes Made in C++ Version

While translating the algorithm from Python to C++, a number of bugs arose. These can primarily be attributed to the differences in how C++ stores and manipulates floats. Small errors in numerical values often lead to loop conditions evaluating to true instead of false or vice versa. 
The following sections outline the major changes. 

## Text File Outputs
The following code block (lines 1779 - 1814) is used to output key information about patches to [cppOutput.txt](wedgeCover%2Fcmake-build-debug%2FcppOutput.txt). The two main pieces of information are the a, b, c, and d corners of the patches, and the points within the superpoints of all patches in each wedge. 
```c++
for(int i = 0; i < cover.patches.size(); i++)
{
    myfile << "Patch " << endl;
    //myfile << fixed;
    //myfile.precision(4);
    myfile << round(cover.patches[i].shadow_fromTopToInnermost_topL_jL * 10000) << endl;
    myfile << round(cover.patches[i].shadow_fromTopToInnermost_topL_jR * 10000) << endl;
    myfile << round(cover.patches[i].shadow_fromTopToInnermost_topR_jL * 10000) << endl;
    myfile << round(cover.patches[i].shadow_fromTopToInnermost_topR_jR * 10000) << endl;

    for(int j = 0; j < cover.patches[i].superpoints.size(); j++)
    {
        myfile << "Superpoint " << endl;
        for(int r = 0; r < cover.patches[i].superpoints[j].points.size(); r++)
        {
            //myfile << fixed;
            //myfile.precision(4);
            Point currentPt = cover.patches[i].superpoints[j].points[r];
            myfile << currentPt.layer_num << " " <<  currentPt.phi << " " << int(currentPt.radius);
            //myfile << fixed;
            //myfile.precision(4);
            myfile << " " << currentPt.z << endl;
        }
    }
}
//myfile << fixed;
//myfile.precision(4);

for(int i = 0; i < cover.patches.size(); i++)
{
    myfile << "[" << round(cover.patches[i].a_corner[0] * 10000) << ", " << round(cover.patches[i].a_corner[1] * 10000) << "]" << endl;
    myfile << "[" << round(cover.patches[i].b_corner[0] * 10000) << ", " << round(cover.patches[i].b_corner[1] * 10000) << "]" << endl;
    myfile << "[" << round(cover.patches[i].c_corner[0] * 10000) << ", " << round(cover.patches[i].c_corner[1] * 10000) << "]" << endl;
    myfile << "[" << round(cover.patches[i].d_corner[0] * 10000) << ", " << round(cover.patches[i].d_corner[1] * 10000) << "]" << endl;
    myfile << endl;
}
```

The following code block (lines 1890 - 1909) is used to output the acceptance percentages for each wedge to [accOutputC.txt](wedgeCover%2Fcmake-build-debug%2FaccOutputC.txt). 
```c++
ofstream accFile;
accFile.open ("accOutputC.txt", ios::out | ios::trunc);

for (int a = 0; a < mean_list.size(); a++)
{
    vector<float> curVector = mean_list[a];
    float average = reduce(curVector.begin(), curVector.end(), 0.0) / curVector.size();
    accFile << "wedge " << a << endl;
    if (average == 100.0)
    {
        accFile << "100.0" << endl;

    }
    else
    {
        accFile << average << endl;
    }
}

accFile.close();
```

## Loop Condition Changes
### Line 1106

Original: 
```c++
while (!(white_space_height <= 0 && (previous_white_space_height >= 0)) && (abs(white_space_height) > 0.000001) && ((patches[patches.size() - 1].c_corner[1] > -1 * env.trapezoid_edges[env.num_layers - 1]) || (white_space_height > 0)) && (current_z_top_index < (int) (data->array[env.num_layers - 1].size() - 1)) && !(repeat_patch) && !(repeat_original))

```

Final: 
```c++
while (!(white_space_height <= 0.0000005 && (previous_white_space_height >= 0)) && (abs(white_space_height) > 0.000005) && ((patches[patches.size() - 1].c_corner[1] > -1 * env.trapezoid_edges[env.num_layers - 1]) || (white_space_height > 0.000005)) && (current_z_top_index < (int) (data->array[env.num_layers - 1].size() - 1)) && !(repeat_patch) && !(repeat_original)) 
```

### Line 1292 

Original:
```c++
bool originalPartialBottom = (original_topL_jL > complementary_apexZ0) && (original_topL_jL < apexZ0) && (abs(patches[patches.size() - 2].straightLineProjectorFromLayerIJtoK(original_topL_jL,z_top_min,1,env.num_layers,0)) < 20 * env.beam_axis_lim);
```

Final:
```c++
bool originalPartialBottom = (original_topL_jL > complementary_apexZ0) && ((original_topL_jL - apexZ0) < -0.0001) && (abs(patches[patches.size() - 2].straightLineProjectorFromLayerIJtoK(original_topL_jL,z_top_min,1,env.num_layers,0)) < 20 * env.beam_axis_lim);
```

### Line 1297

Original:
```c++
bool complementaryPartialBottom = (complementary_topL_jR > complementary_apexZ0) and (complementary_topL_jR < apexZ0) and (abs(self.patches[-1].straightLineProjectorFromLayerIJtoK(complementary_topL_jR,z_top_min,1,self.env.num_layers,0))<20*self.env.beam_axis_lim);
```
Final:
```c++
bool complementaryPartialBottom = (complementary_topL_jR > complementary_apexZ0) && ((complementary_topL_jR - apexZ0) < -0.0001) && (abs(patches[patches.size() - 1].straightLineProjectorFromLayerIJtoK(complementary_topL_jR,z_top_min,1,env.num_layers,0)) < 20 * env.beam_axis_lim);

```

### Line 1338 

Original: 
```c++
if (horizontalShiftTop > 0 or horizontalShiftBottom > 0)
```

Final: 
```c++
if (horizontalShiftTop > 0.000001 or horizontalShiftBottom > 0)
```

### Line 1344: 

Original: 
```c++
while (((horizontalShiftTop > 0 && originalPartialTop && complementaryPartialTop) || (horizontalShiftBottom > 0 && originalPartialBottom && complementaryPartialBottom)) && doShiftedPatch && (horizontalOverlapTop <= 0) && (horizontalOverlapBottom <= 0) && (newGapTop < 0 || newGapBottom < 0))

```

Final: 
```c++
while (((horizontalShiftTop > 0.000001 && originalPartialTop && complementaryPartialTop) || (horizontalShiftBottom > 0.000001 && originalPartialBottom && complementaryPartialBottom)) && doShiftedPatch && (horizontalOverlapTop <= 0) && (horizontalOverlapBottom <= 0) && (newGapTop < 0 || newGapBottom < 0))

```

