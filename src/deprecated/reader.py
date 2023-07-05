from src.deprecated.spaceptcollection import *
from src.deprecated.wedgedata import WedgeData
from src.coverers.data_structs import Environment

import ast
import time

def readFile(filepath, stop = 128):
    events = []

    start = time.time()
    with open(filepath) as f:
        idx = 0
        for line in f:
            if line.strip():
                event = SpacePtCollection(idx)
                tuples = ast.literal_eval(line)
                num_layer = 0

                for tuple in tuples:
                    newSpacePoint = SpacePoint(tuple[0], tuple[1], tuple[2], tuple[3])
                    event.appendPoint(tuple[0], newSpacePoint)
                    num_layer = max(num_layer, tuple[0])
                
                event.retrieveNumLayer(num_layer)
                events.append(event)
                idx += 1
            if idx == stop:
                return events
    print(f"Time Taken to Read File : {time.time() - start}s")
    return events

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


