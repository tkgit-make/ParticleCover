from src.coverers.data_structs import Environment, Point
import time


def readFile(filepath, stop:int = 128, performance:bool=False): 
    start = time.time()
    
    # Create a list of events, with each event consisting of an environment and its list of Points 
    events = [] 
    with open(filepath) as f: 
        line_index = 0
        for line in f: 
            if line.strip():
                
                tuples = line.strip()[1:-1].split("),(")
                tuples = [tup.split(",") for tup in tuples]
                
                list_of_Points = [Point(int(tupl[0]), int(tupl[1]), float(tupl[2]), float(tupl[3])) for tupl in tuples]
                radii = sorted(list(set([point.radius for point in list_of_Points])))
                num_layers = len(radii)
                
                env = Environment(top_layer_lim= 100.0, 
                                  bottom_layer_lim = 15.0, 
                                  num_layers=num_layers, 
                                  radii=radii
                                )
                events.append((env, list_of_Points))
                
            line_index += 1 
            if line_index == stop: 
                break 
    
    if performance == True:
        print(f"Time Taken to Read File : {time.time() - start}s")
        
    return events 

