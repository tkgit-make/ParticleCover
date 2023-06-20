from wedgedata import WedgeData
from spaceptcollection import SpacePtCollection
from data import Environment

def convertToDataset(data:SpacePtCollection):
    layer_radius = data.spacePoints[1][0].radius
    n_points = [len(data.spacePoints[i+1]) for i in range(0, data.num_layers)]

    # Create a list of data.num_layers * list of respective SpacePoints in the layer
    pts = []
    for i in range(0, data.num_layers):
        pts.append(sorted(data.spacePoints[i+1], key = lambda pt: pt.z))

    env = Environment(top_layer_lim=100, bottom_layer_lim=15, layers=data.num_layers, radii=layer_radius)
    wedgeData = WedgeData(env, n_points, pts)

    return wedgeData

