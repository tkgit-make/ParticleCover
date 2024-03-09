from  src.testers.test import *

#wedge_test(movie = True, savefig=True)
#wedge_test()
#wedge_test(movie = True)
#wedge_test(movie = True, savefig=True)
wedge_test(lining = 'makePatches_ShadowQuilt_fromEdges',apexZ0=0,top_layer_cutoff=50, wedges=[244,245], z0_spacing=0.05, leftRightAlign=False, show_acceptance_of_cover=True, accept_cutoff=15, movie = True, movieFigSizeScale=3)