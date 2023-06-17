from utils import *
from prettytable import PrettyTable
from decimal import Decimal

import matplotlib.pyplot as plt
import numpy as np
import math

class SpacePoint:
    def __init__(self, layer_num, radius, phi, z):
        self.layer_num = layer_num
        self.radius = radius
        self.phi = phi
        self.z = z

class EventData:
    def __init__(self, event_num):
        self.event_num = event_num
        self.spacePoints = dict()

    def appendPoint(self, layer_num, SpacePoint:SpacePoint):
        appendToDict(self.spacePoints, layer_num, SpacePoint)
    
    def retrieveNumLayer(self, num_layers):
        self.num_layers = num_layers

    def locateWedgeBound(self, nPhiSlices, p, B):
        self.wedgeBound = np.ndarray(shape = (self.num_layers, nPhiSlices, 2))
        self.wedgeCenters = np.ndarray(shape = (nPhiSlices, 4))
        wedgeRadius = (100 * p) / (0.3 * B) # in our case, wedgeRadius = 1667 cm

        firstPt = (0, 0)
        r = self.spacePoints[1][0].radius
        for i in range(0, self.num_layers):
            for j in range(0, nPhiSlices):
                if i == 0:
                    angle_wrt_org = 0 + j * (2 * math.pi) / nPhiSlices
                    secondPt = (r * math.cos(angle_wrt_org), r * math.sin(angle_wrt_org))
                    [C_x1, C_y1], [C_x2, C_y2] = findCenter(firstPt[0], firstPt[1], secondPt[0], secondPt[1], wedgeRadius)
                    
                    C1_angle_wrt_org = math.atan2(C_y1, C_x1)
                    C1_angle_wrt_org = convertNegRadian(C1_angle_wrt_org)

                    C2_angle_wrt_org = math.atan2(C_y2, C_x2)
                    C2_angle_wrt_org = convertNegRadian(C2_angle_wrt_org)

                    C1_C2_angle_diff = C1_angle_wrt_org - C2_angle_wrt_org
                    C1_C2_angle_diff = convertNegRadian(C1_C2_angle_diff)

                    C2_C1_angle_diff = C2_angle_wrt_org - C1_angle_wrt_org
                    C2_C1_angle_diff = convertNegRadian(C2_C1_angle_diff)

                    if C1_C2_angle_diff < C2_C1_angle_diff:
                        self.wedgeCenters[j][0] = C_x2
                        self.wedgeCenters[j][1] = C_y2
                        self.wedgeCenters[j-1 if j > 0 else nPhiSlices - 1][2] = C_x1
                        self.wedgeCenters[j-1 if j > 0 else nPhiSlices - 1][3] = C_y1
                    else:
                        self.wedgeCenters[j][0] = C_x1
                        self.wedgeCenters[j][1] = C_y1
                        self.wedgeCenters[j-1 if j > 0 else nPhiSlices - 1][2] = C_x2
                        self.wedgeCenters[j-1 if j > 0 else nPhiSlices - 1][3] = C_y2
                    
                    self.wedgeBound[i][j][0] = angle_wrt_org
                    self.wedgeBound[i][j][1] = angle_wrt_org + (2 * math.pi) / nPhiSlices
                else:
                    curRWedgeCenter = (self.wedgeCenters[j][0], self.wedgeCenters[j][1])
                    curLWedgeCenter = (self.wedgeCenters[j][2], self.wedgeCenters[j][3])
                    firstLayerRBound = self.wedgeBound[0][j][0]
                    firstLayerLBound = self.wedgeBound[0][j][1]
                    r1 = self.spacePoints[i+1][0].radius

                    (i_x1, i_y1), (i_x2, i_y2) = getIntersection(firstPt[0], firstPt[1], r1, 
                                                                curRWedgeCenter[0], curRWedgeCenter[1], wedgeRadius)
                    self.wedgeBound[i][j][0] = determineWhichIntersection(i_x1, i_y1, i_x2, i_y2, firstLayerRBound)
                    (i_x1, i_y1), (i_x2, i_y2) = getIntersection(firstPt[0], firstPt[1], r1,
                                                                curLWedgeCenter[0], curLWedgeCenter[1], wedgeRadius)
                    self.wedgeBound[i][j][1] = determineWhichIntersection(i_x1, i_y1, i_x2, i_y2, firstLayerLBound)
        
    def calculateWedgeOverlap(self, nPhiSlices):
        myTable = PrettyTable(["Layer", "Overlap %"])

        for i in range(0, self.num_layers):
            overlap = 0.000000000
            for j in range(0, nPhiSlices):
                curEndAngle = self.wedgeBound[i][j][1]
                nextStartAngle = self.wedgeBound[i][j+1][0] if j < nPhiSlices - 1 else self.wedgeBound[i][0][0]

                if nextStartAngle == 0:
                    nextStartAngle = 2 * math.pi
                else:
                    nextStartAngle = nextStartAngle if j < nPhiSlices - 1 else nextStartAngle - 2 * math.pi

                nextWedge = j + 1 if j < nPhiSlices - 1 else 0

                if not math.isclose(curEndAngle, nextStartAngle):
                    overlap += (curEndAngle - nextStartAngle)
                
            overlap_percentage = overlap * 100 / (2 * math.pi)
            rounded_percentage = round(Decimal(str(overlap_percentage)), 9)
            myTable.add_row([str(i+1), str(rounded_percentage)])

        print(myTable)
    
    def produceWedgeData_(self, nPhiSlices):
        wedges = dict()
        for i in range(0, self.num_layers):
            ptsPerLayer = np.array(self.spacePoints[i + 1])
            for pt in ptsPerLayer:
                angle_wrt_org = convertNegRadian(pt.phi)
                respectiveWedges = np.where(angle_wrt_org < self.wedgeBound[:][:][1] | 
                                            angle_wrt_org >= self.wedgeBound[:][:][0])
                print(respectiveWedges)

    """
    Produces `nPhiSlices` wedge data with straight line boundaries.
    """
    
    def produceWedgeData(self, nPhiSlices):
        wedges = dict()
        for ptsPerLayer in self.spacePoints.values():
            ptsPerLayer = np.array(ptsPerLayer)
            for pt in ptsPerLayer:
                angle_wrt_org = math.degrees(pt.phi) % 360
                sector_num = int(angle_wrt_org / (360 / nPhiSlices)) + 1
                appendToDict(wedges, sector_num, (pt.layer_num, pt.radius, pt.phi, pt.z))
        
        filename = f'output/{nPhiSlices}_wedges_event_{self.event_num}.txt'
        with open(filename, 'w') as f:
            for i in range(1, nPhiSlices + 1):
                line_to_write = str(wedges[i]).replace('[', '').replace(']', '')
                f.write(line_to_write + '\n')

    """
    Methods designed to print or return number of SpacePoints per layer in an event.
    """

    def printNumPts(self):
        prt_string = ""
        total_pts = 0
        for i in range(1, self.num_layers + 1):
            if i != self.num_layers:
                prt_string = prt_string + str(len(self.spacePoints[i])) + ", "
            else:
                prt_string += str(len(self.spacePoints[i]))

            total_pts += len(self.spacePoints[i])
        # print(f"Pts per Layer: {prt_string}, Total: {total_pts} \n")

        return total_pts
    
    def returnPtsPerLayer(self):
        ptsPerLayer = np.zeros(self.num_layers)
        for i in range(1, self.num_layers + 1):
            ptsPerLayer[i - 1] = len(self.spacePoints[i])
        
        return ptsPerLayer
    
    def returnUniquePoints(self):
        uniquePoints = set()
        for i in range(1, self.num_layers + 1):
            pts_in_layer = self.spacePoints[i]
            for j in range(0, len(pts_in_layer)):
                pt = pts_in_layer[j]
                uniquePoints.add((pt.layer_num, pt.radius, pt.phi, pt.z))
        return uniquePoints

    """
    Methods below are designed for visualization.
    """

    def plotCylindrical(self):
        for i in range(1, self.num_layers + 1):
            z_arr = [x.z for x in self.spacePoints[i]]
            r_arr = self.spacePoints[i][0].radius * np.ones(len(z_arr)) 
            plt.scatter(z_arr, r_arr, s = 2, c = "b")
        
        plt.yticks(np.arange(0, self.spacePoints[self.num_layers][0].radius + 1, self.num_layers))
        plt.title(f'Distribution of All SpacePoints in Event {self.event_num}')
        plt.ylabel('Radius')
        plt.xlabel('z-coordinate')
        plt.show()
    
    def plotCartesian(self):
        fig, ax = plt.subplots()
        ax.set_box_aspect(1)
        for i in range(1, self.num_layers + 1):
            x_arr = [(x.radius * math.cos(x.phi)) for x in self.spacePoints[i]]
            y_arr = [(x.radius * math.sin(x.phi)) for x in self.spacePoints[i]]
            ax.scatter(x_arr, y_arr, s = 2, c = "b")

        plt.title(f'Distribution of All SpacePoints in Event {self.event_num}')
        plt.ylabel('y')
        plt.xlabel('x')
        plt.show()
    
    def nSlicePoints(self, nPhiSlices):
        fig, ax = plt.subplots()
        ax.set_box_aspect(1)

        plt.title(f'Distribution of All SpacePoints in Event {self.event_num}, Sliced in {nPhiSlices} Sectors')
        plt.ylabel('y')
        plt.xlabel('x')

        for i in range(1, self.num_layers + 1):
            for point in self.spacePoints[i]:
                x = point.radius * math.cos(point.phi)
                y = point.radius * math.sin(point.phi)

                angle_wrt_org = math.degrees(point.phi) % 360
                sector_num = int(angle_wrt_org / (360 / nPhiSlices))

                ax.plot(x, y, marker = 'o', markersize = 2, color = colorList[sector_num % 128])
        
        plt.show()