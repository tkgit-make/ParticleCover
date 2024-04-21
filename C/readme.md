# Required Particle Cover Algorithm Structures & Implementation Status

## Class Definitions

### Point
- Constructor: `Point(int layerNum, float rad, float ph, float zVal)` ✅

### Environment
- Constructor: `Environment(float top_layer_limI, float beam_axis_limI, int num_layersI, vector<float> radiiI)` ✅

### DataSet
- Constructors:
  - `DataSet()` ✅
  - `DataSet(Environment& envI)` ✅
- Methods:
  - `importData(vector<Point> data_array)` ✅
  - `addBoundaryPoint(float offset)` ✅

### Event
- Constructor: `Event(Environment envI, vector<Point> listP)` ✅

### Line
- Constructor: `Line(Environment envI, float start, float slopeI)` ✅

### Parallelogram
- Constructor: `parallelogram(int layer_numI, float z1_minI, float z1_maxI, ... )` ✅

### wedgeSuperPoint
- Constructor: `wedgeSuperPoint(vector<Point> pointsI)` ✅
- Operators:
  - `operator==` ✅

### wedgePatch
- Constructor: `wedgePatch(Environment envI, vector<wedgeSuperPoint> superpointsI, float apexZ0I)` ✅
- Methods:
  - `straightLineProjectorFromLayerIJtoK(...)` ✅
  - `getParallelograms()` ✅
  - `getShadows(float zTopMin, float zTopMax)` ✅
  - `straightLineProjector(float z_top, float z_j, int j)` ✅
  - `get_acceptanceCorners()` ✅
  - `get_end_layer()` ✅

### wedgeCover
- Constructor: `wedgeCover(Environment envI, DataSet& dataI)` ✅
- Methods:
  - `add_patch(wedgePatch curr_patch)` ✅
  - `get_index_from_z(int layer, float z_value, string alignment)` ✅
  - `delete_patch(int index)` ✅
  - `solve(...)` ✅
  - `makePatches_ShadowQuilt_fromEdges(...)` ✅
  - `makePatch_alignedToLine(...)` ✅

### FileReader
- Methods:
  - `splitString(string str, string splitter)` ✅ (equivalent)
  - `readFile(string filepath, int stop, bool performance)` ✅ (equivalent)

### Tester
- Method: `wedge_test(...)` ✅

## Main Function

- Method: `int main()` ✅

### Not Used
> ### LineGenerator
> - Constructor: `LineGenerator(Environment envI, float startI)` ✅
> - Method: `generateEvenGrid(int n)` ✅ 
> 
> ### wedgePatch
> - Method: `getParallelograms_v1()` ✅
> 
> ### Parallelogram_v1
> - Constructor: `parallelogram_v1(int layer_numI, float top_layer_zminI, float top_layer_zmaxI, ...)` ✅
> 
> ### wedgeCover
> - Method: `tester()`



