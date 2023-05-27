# ParticlePartition

## Generating Data 
Before we do anything, we must import the data.py module and create an environment. An environment object contains information about the kind of collider that the data will be generated in. It can be initialized with a constructor, and it has the following attributes 
```
env = Environment() 
print(env.top_layer_lim)        # a float representing the upper and lower limits of the top layer (default 1.0m)
print(env.bottom_layer_lim)     # a float representing the upper and lower limits of layer 0 in m (default 0.15m) 
print(env.layers)               # a integer representing the number of layers, excluding layer 0 (default 5) 
print(env.radii)                # a float representing the distance between consecutive layers in m (default 5.0m)
```

Once we have generated the environment, we can generate a dataset by constructing the `Dataset` object with the environment and number of points per layer as its argument. 
``` 
data = DataSet(env, n_points = 150) 
```
To access the data, we just need to call `data.array`, which gives us a $L \times N$ matrix $A$, where ($L$ is the number of layers and $N$ is the number of points in each layer). Furthermore, for each $i = 1, \ldots, L$, row $i$ is already sorted in its points from least to greatest. Plotting it usig the `plot` method should give a shape that approximately looks like an inverted symmetric trapezoid. 

## Constructing Covers

We define three things: 
 - A **cover** $C$ is a finite family of patches $C = \{P_i\}_{i=1}^n$. 
 - A **patch** $P$ is a list of $L$ superpoints, with the 1st superpoint element representing the superpoint of the 1st layer. 
 - A **superpoint** $S$ is a collection of 16 consecutive points in one layer. 
For further developers, I would recommend only adding extra attributes and methods to the `Patch` class. The cover itself should essentially be a collection of patches, no more, and a superpoint is something I've used for my personal algorithm. 
