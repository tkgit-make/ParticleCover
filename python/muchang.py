from src.readers.reader import * 
from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.coverers.wedgecover import *
from src.testers.test import *
from src.debug import * 
import matplotlib.pyplot as plt 
import time 
import ast


wedge_test(lining= "makePatches_Projective_Rightleft", apexZ0 = 10, wedges = [56, 57], top_layer_cutoff = 50, show_acceptance_of_cover=True)
