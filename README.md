
# ParticlePartition

## Environment and Data Objects
Before we do anything, we must import the data.py module and create an environment. An environment object contains information about the kind of collider that the data will be generated in. It can be initialized with a constructor, and it has the following attributes 
```
env = Environment() 
print(env.top_layer_lim)        # a float representing the upper and lower limits of the top layer (default 100cm)
print(env.bottom_layer_lim)     # a float representing the upper and lower limits of layer 0 in m (default 15cm) 
print(env.layers)               # a integer representing the number of layers, excluding layer 0 (default 5) 
print(env.radii)                # a float representing the distance between consecutive layers in m (default 5.0cm)
```

Once we have generated the environment, we can generate a dataset by constructing the `Dataset` object with the environment and number of points per layer as its argument. 
``` 
data = DataSet(env, n_points = 150) 
```
To access the data, we just need to call `data.array`, which gives us a $L \times N$ matrix $A$, where ($L$ is the number of layers and $N$ is the number of points in each layer). Furthermore, for each $i = 1, \ldots, L$, row $i$ is already sorted in its points from least to greatest. Plotting it usig the `plot` method should give a shape that approximately looks like an inverted symmetric trapezoid.

## Using Realistic Versions of Collider Data

The `DataSet` object assumes that the data points are generated from uniform distribution. Hence it may be possible that such data may not holistically represent the data actually generated from the collider. To use realistic versions of collider data, we can simply read in the data file given in .txt format using `readFile` function in reader.py module.

`readFile` reads in the given data file line by line, where each line represents a spacepoint collection that can be used as an input data for our algorithms. As every line represents tuples separated by a comma, `readFile` stores the information in every tuple as an instance of class `SpacePoint` and appends all of the `SpacePoint` objects into an instance of class `SpacePtCollection`. In the end, `readFile` function would return a list of `SpacePtCollection` objects.

Then, one can use the converter function `convertToDataset` in converter.py module to convert each `SpacePtCollection` object into an instance of class `WedgeData`, which can be found in wedgedata.py module. `WedgeData` is a subclass of class `DataSet`, thus inheriting its attributes. Yet `WedgeData` class differs from `DataSet` in that the number of points for each layer vary instead of being the same. Hence the attribute `n_points` is a list with the dimensions $L \times 1$, with each element representing the number of points in the corresponding layer. Moreover, the attribute `array` is a nested list of `SpacePoint` objects, with the shape being $(L, )$. The shape being described in this way can be ascribed to the fact that each layer may have different number of points, so `array` may not be a $L \times N$ matrix, which is the case for our `DataSet` object. Essentially, `convertToDataset` function returns a `WedgeData` object that can be directly inputted into our algorithms.

Here is an example of how `readFile` and `convertToDataSet` functions can be used to read in a wedge data file and be converted to be applicable for use with our algorithms:
```
from reader import readFile                 # Import readFile function
from converter import convertToDataset      # Import convertToDataset function

wedges = readFile("Wedge_Data.txt")         # readFile reads Wedge_Data.txt, stores into a list of SpacePtCollection objects

firstWedge = convertToDataset(wedges[0])    # Converts the first element in wedges into a WedgeData object
firstWedge.plot(show_lines = True)          # Plots the converted first wedge data in a r vs. z plot
```
## Running Tests on Methods

The file `test_modules.py` contains the function `wedge_test` that generates the plots found in the LaTex document. All that is required to generate the plots is running the function on a python file that has the wedgeData text files in the same directory. The files should be named with the format `wedgeData_{VERSION}_128.txt`. Below is a list of the arguments in the function.
1. lining (str, optional): solving method, default is solveS
2. solve_at (int, optional): the z values the patches are being made in accordance to
3. z0 (num or list, optional): array of z0 values we are testing over, default is range of -15 to 15 with 0.5 spacing
4. n (int, optional): ppl (point per patch per layer), default is 16
5. wedges (list, optional): which wedges, enter in list of [starting wedge, ending wedge]. Default is [0, 128]
6. lines (int, optional): how many line to test acceptance with at each z0 value
7. savefig (bool, optional): `True` to save figure
8. v (str, optional): version of data, ensure data file is in directory as "wedgeData_{v}_128.txt"

## Constructing wedgeData Covers
In order to accommodate for the different data structure and avoid confusion, I made a new module for solving with `wedgeData` class called `wedgeCover`. It operates exactly the same as the `cover` class. There are currently four solving methods: `solveS`, `solveS_reverse`, `solveS_center2`, and `solveQ`. The easiest way to solve for a plot and obtain a visualization is to use the methods `solve` from `wedgeCover`.

The arguments for the `solve` method are
1. lining (str): solving method
2. z0 (num or list): where on the z axis the patches are being solved with respect to. This can be a single number or a list.
3. n (int): points per patch per layer
4. nlines (int): number of lines used in visualization.
5. show (bool): prevents line generation for visualization. This should be set to `False` when patches are solved for at multiple z0's.

After solving, the cover can be visualized with the `plot` method. An example of how to solve and visualize a cover follows.
```
env = Environment()                             #init environment
events = readFile('wedgeData_v3_128.txt', 128)  #read file that has wedge data
wedge1 = convertToDataset(events[0])            #convert into wedgeData format
cover = wedgeCover(env, wedge1)                 #init wedgeCover class
cover.solve('solveS', z0 = 0, show = True)      #solve for cover with respect to z0 = 0 
cover.plot()                                    #plots cover
```

## Point and Line Objects 
A line can be characterized by two things: a point $(x_0, y_0)$ on the line and its slope $m$. It is clear that $y_0$, the height, should be $0$, but the $x_0$ may vary (by default we set $x_0 = 0$). Therefore, a line should have two parameters: the $x_0$ value of its originating point and the slope $m$. We can construct it by calling: 
``` 
line = Line(env, start, slope)
``` 
where `env` is the embedding environment, `start` is $x_0$, and `slope` is $m$. It would be nice to add a check within this constructor to determine whether the inputted slope is viable (i.e. within the range of the maximum and minimum slope), but I thought this might take away computational speed, so I did not implement this. 

Ultimately, the line object is described by the `env.layer` = $5$ points of the line that lies on each layer. This should also be stored in the `points` attribute. 

A `LineGenerator` object simply generates lines within an environment and that has a start value equal to what is inputted. 
```
lg = LineGenerator(env, start) 
``` 

We can either choose to generate `n` lines with equal spacing across the entire environment (like a grid over angles) by calling `generateGridLines(n)` or generate them according to a uniform distribution (over the angle space $\theta = \arctan(m)$) by calling `generateRandomLines(n)`. They both return a Python list of `Line` objects. 


## Constructing Covers

We define three things:
 - A **cover** $C$ is a finite family of patches $C = \{P_i\}_{i=1}^n$. 
 - A **patch** $P$ is a list of $L$ superpoints, with the 1st superpoint element representing the superpoint of the 1st layer. 
 - A **superpoint** $S$ is a collection of 16 consecutive points in one layer. 
For further developers, I would recommend only adding extra attributes and methods to the `Patch` class. The cover itself should essentially be a collection of patches, no more, and a superpoint is something I've used for my personal algorithm. 

#### 1. Superpoints
A superpoint is characterized by the smallest interval of a layer that contains all 16 points. By abuse of notation, we can mathematically express it either a list of 16 numbers
$$S = [s_1, s_2, \ldots, s_{16}] \text{ with } s_1 < s_2 < \ldots < s_{16}$$
or as the closed interval between the minimum and maximum values 
 $$S = [s_1, s_{16}]$$ 
If we would like to see if a float value $p$ is contained within a superpoint, then we can run the `contains(p)` method, which returns a boolean describing whether $p \in [s_1, s_{16}]$. This is the basic functionality of the superpoint. 

#### 2. Patches
A patch is simply a collection of 5 superpoints. They should all be from different layers, but this check (for computational reasons) have not been implemented. Given superpoints $S_1, \ldots, S_5$, a patch $P$ is mathematically expressed as 
 $$P = (S_1, S_2, S_3, S_4, S_5)$$
Note that this is implemented as a tuple, since there should be no further mutation of this collection. We initialize a `Patch` object by inputting a tuple of superpoints objects and the underlying environment. 
```
# assume the sp1,...,sp5 have been initialized 
sp_array = (sp1, sp2, sp3, sp4, sp5) 
patch = Patch(env, sp_array) 
``` 
This constructor actually will check that the `sp_array` contains `env.layers` superpoints. I figured that this check should be important at least for now since this is a common error. The `contains` method is very important, since this allows the user to determine if a patch contains a line. 
```
patch.contains(line) 
# returns true if patch contains line 
``` 
Remember previously that a line can be characterized by the `env.layers`= $5$ points lying on each layer. We can describe a line as 
$$l = [l_1, l_2, l_3, l_4, l_5]$$ 
with $l_i$ is the float value on the $i$th layer. We can also characterize a patch as a 5-tuple of superpoints, which are essentially closed intervals. 
$$P = ([\min{S_1}, \max{S_1}], [\min{S_2}, \max{S_2}], [\min{S_3}, \max{S_3}], [\min{S_4}, \max{S_4}], [\min{S_5}, \max{S_5}])$$ 
Therefore, the contains method $\mathcal{C}$ determines whether $l \in P$ by determining whether $l_i \in [\min{S_i}, \max{S_i}]$ for $i \in [5]$. That is
$$\mathcal{C}(P, l) \coloneqq \begin{cases} \text{True} & \text{ if } l_i \in [\min{S_i}, \max{S_i}] \text{ for } i \in [5] \\ \text{False} & \text{ if else } \end{cases}$$
This function will be clearly useful for constructing a cover and performance testing it. 


#### 3. Covers

A cover $C$ is essentially a list of patches. Therefore, given $n$ patches, $P_1, \ldots, P_n$, a cover is mathematically described as 
$$C = [P_1, \ldots, P_n]$$
Now it is initialized in a bit of a different way. It takes in an `Environment` object and a `DataSet` object, and constructs an empty cover first, without any patches inside. It has 5 attributes 
```
cover = Cover(env, data)
cover.n_patches             # number of patches = 0
cover.env                   # the embedding environment object
cover.data                  # the inputted DataSet object
cover.patches               # an empty list of patches
cover.superPoints           # a list of lists of superpoints
``` 
In the constructor, we do the following for each layer $i$, which contains $n$ points each. 
- We initialize a list representing the superpoints in layer $i$. 
- A superpoint is constructed by taking the first 16 points $1, \ldots, 16$ and added to this list. 
- We take points $16, \ldots, 31$ to make a second superpoint and add this to the list. Note that there is an overlap of one point between adjacent superpoints. 
- We continue this until there are less than 16 points left for the final superpoint. 
- The final superpoint is constructed by taking the last 16 points and added to the list. 
For $n = 150$, this should create $10$ superpoints for each layer, and if we have 5 layers, `cover.superPoints` should be a list of 5 lists, each containing 10 superpoints. 

#### Attaining Estimations of the Minimal Cover 
The `solve` method of the `Cover` object is where the main algorithm resides. Now that we have a list of list of superpoints upon initialization, which we will label
$$[[S_{1,1}, S_{1, 10}], \ldots, [S_{5,1}, S_{5, 10}]]$$
where $S_{i, j}$ represents the $j$ th superpoint in the $i$ th layer. The algorithm is run as such: 
1. We construct a `LineGenerator` object and have it generate 100 equally spaced lines in the environment. These list of lines are generated "from left to right," meaning that the line that resides in the leftmost portion is the first element of the list, and the rightmost in the last element. 
2. For the first line $l = [l_1, l_2, l_3, l_4, l_5]$, we look at $l_i$ and find which superpoint in the $i$th list of `cover.superPoints` it is contained in. After $i = 1, \ldots, 5$, we should have a collection of $5$ superpoints, which are stored in the `patch_ingredients` variable. 
3. We construct the patch $P$ from the elements of `patch_ingredients`, which is guaranteed by construction to contain line $l$ and add it to `cover.patches` (and increment `cover.n_patches` by $1$). 
4. We move onto the next line and do the same. Note that for each line, we do not need to iterate through all the superpoints. This is because we're going from left to right, and so given that line $l_k$ is contained in the patch 
$$(S_{1, 3}, S_{2, 3}, S_{3, 4}, S_{4, 2}, S_{5, 1})$$
the next line $l_{k+1}$ must be contained in the superpoint that comes at least after those of line $l_k$. Therefore, the corresponding patches must be of form 
$$(S_{1, j_1}, S_{2, j_2}, S_{3, j_3}, S_{4, j_4}, S_{5, j_5})$$
where $j_1 \geq 3, j_2 \geq 3, j_3 \geq 4, j_4 \geq 2, j_5 \geq 1$. This will save a lot of computational time by design. 

5. We repeat the above steps for all 100 lines. To detect repeated patches, unfortunately we cannot use hashsets since patches are not hashable objects (implementing this would be very nice). Therefore, for every generated patch $P_k$ for line $l_k$, we just compare it with the latest patch that is stored in `cover.n_patches`. Just this one comparison is sufficient since again, since the patches are generated from left to right, and therefore, $l_k \not\in P_{k-1} \implies l_k \not\in P_1, \ldots, P_{k-2}$. 

6. After all repetitions we have our desired cover stored in the `cover.patches` list. 

### Additional Solving Algorithms
The `solveS` method for the `Cover` class generates a cover for a given set of data $D_{ij}$ where $i$ is the layer and $j$ indexes the points in the layer. The goal was to minimize the amount of patches. This assumes 150 points but the code will work for more or less points. The algorithm proceeds as follows:

1. The first patch is contrasted from a superpoint of the first 16 points of each layer. $$S_i = [D_{i,1}, D_{i,2}, ..., D_{i,16}]$$

2. The `S_loop15` method is run, which is the main part of the algorithm. This generates all the other patches.

3. `S_loop15` method loops through the rightmost point in each layer. It takes the point and 're-scales' it based on the r-value. Let the re-scaled value be $D_{i, j}'$ where $$D'{i, j} = \frac{D_{i, j}}{5i}$$ The minimum of $D_{i, j}'$ for $i \in [1, 5]$ is $i_{min}$. $i_{min}$ is stored in `min_index`. 
    
4. Next re-scale each layer to find the index of the point closest in value to $D_{i_{min},j}'$. Let the index of this point be $j = m$. There are three cases.
    
    1. $m$ is the leftmost point of a given layer ($m<1$), which happens the first few patches due to the radius offset. In this case, the first 16 points is the superpoint for that layer. $$S_i = [D_{i,0}, D_{i,1}, ..., D_{i,16}]$$

    2. There are not enough points left for 16 new points ($m >= 150-16$.) The superpoint is then the rightmost 16 points in the layer. $$S_i = [D_{i,(150-16)}, D_{i,(150-15)}, ..., D_{i,150}]$$

    3. The closest index is neither of the above cases. The superpoint starts at $m-1$ to ensure overlap then picks the next 15 points. $$S_i = [D_{i,(m-1)}, D_{i,(m)}, ..., D_{i,(m+15)}]$$
    
    These superpoints are added to the list `patch_ingredients`.
    
5. Checks to see if the new patch equals the last patch added. This happens when Case 2 is true for all five layers. Once that happens, the algorithm terminates. Otherwise, `S_repeated` is run again.

After the algorithm terminates, the patches are stored in `cover.patches`.

Using the same logic, there are several methods that vary on this premise.
1. `solveS`: Starts with the left most 16 points and works its way right to generate the cover using `S_repeated`.
2. `solveS_reverse`: Starts with the right most 16 points and works its way left to generate the cover using `S_repeated_reverse`.
3. `solveS_center1`: Starts with the center 16 points in each layer then uses `S_repeated` on the right half of data and `S_repeated_reverse` on left half of data.
4. `solveS_center2` Starts with a superpoint containing the 8 points to the left and right of 0 in each layer then uses `S_repeated` on the right half of data and `S_repeated_reverse` on left half of data.
5. `SolveQ`: Starts at Q1 and Q3 by z value and runs `solveS_center2` starting at those points and ending at the center
