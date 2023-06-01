The `michelle` method for the `Cover` class generates a cover for a given set of data $D_{ij}$ where $i$ is the layer and $j$ indexes the points in the layer. The goal was to minimize the amount of patches. This assumes 150 points but the code will work for more or less points. The algorithm proceeds as follows:

1. The first patch is contrasted from a superpoint of the first 16 points of each layer. $$S_i = [D_{i0}, D_{i1}, ..., D_{i16}]$$
2. The `mk_repeated` method is run, which is the main part of the algorithm. This generates all the other patches.

3. `mk_repeated` loops through the rightmost point (index of 15) in each layer and finds which has the largest slope. The index of the layer is stored in `min_index`. 

4. Next, the layers are looped through to find the point in each layer that has the closest slope to $D_{min\_index,j}$. Let the index of this point be $j = m$. There are three cases. 
	1. $m$ is the leftmost point of a given layer ($m<1$), which happens the first few patches due to the radius offset. In this case, the first 16 points is the superpoint for that layer. $$S_i = [D_{i0}, D_{i1}, ..., D_{i16}]$$
	2. There are not enough points left for 16 new points ($m >= 150-16$.) The superpoint is then the rightmost 16 points in the layer. $$S_i = [D_{i(150-16)}, D_{i(150-15)}, ..., D_{i150}]$$
	3. The closest index is neither of the above cases. The superpoint starts at $m-1$ to ensure overlap then picks the next 15 points.$$S_i = [D_{i(m-1)}, D_{i(m)}, ..., D_{i(m+15)}]$$
	These superpoints are added to the list `patch_ingredients`.

5. Checks to see if the new patch equals the last patch added. This happens when Case 2 is true for all five layers. Once that happens, the algorithm terminates. Otherwise, `mk_repeated` is run again. 

After the algorithm terminates, the patches are stored in `cover.patches`.