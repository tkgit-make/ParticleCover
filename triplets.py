import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import *
from test_modules import *
import math
import cv2 
import os 
import glob
import matplotlib
import ast
import time
from reader import *
from converter import *
from wedgecover import *
import matplotlib as mpl

def create_triplets(version = 'v2', random = False, anchor_layer = 1):
    if (anchor_layer < 1) or (anchor_layer > 3):
        raise('Anchor layer must not be the top or bottom layer')
    env = Environment()
    events = readFile(f'wedgeData_{version}_128.txt', 128)
    tripletsz = []
    tripletsphi = []
    hitres = 0.005 #50 microns
    for i in range(len(events)):
        wedge1 = convertToDataset(events[i])
        cover = wedgeCover(env, wedge1)
        cover.solveS_center2()
        for patch1 in cover.patches:
            
            if random == True:
                layer1 = np.random.rand(16) * hitres
                layer2 = np.random.rand(16) * hitres 
                layer3 = np.random.rand(16) * hitres
                for point2 in layer2:
                    for point1 in layer1:
                        for point3 in layer3:
                            tripletsz.append([point1, point2, point3])
                            tripletsphi.append([point1/15, point2/20, point3/25])
            
            else:
                layer1 = patch1.superpoints[anchor_layer - 1]
                layer2 = patch1.superpoints[anchor_layer]
                layer3 = patch1.superpoints[anchor_layer + 1]
                
                for point2 in layer2.points:
                    for point1 in layer1.points:
                        for point3 in layer3.points:
                            tripletsz.append([point1.z, point2.z, point3.z])
                            phi1 = point1.phi
                            phi2 = point2.phi
                            phi3 = point3.phi
                            if ((point1.phi - point2.phi) > np.pi):
                                phi1 = point1.phi-(2*np.pi)
                            if ((point2.phi - point1.phi) > np.pi):
                                phi2 =  point2.phi-(2*np.pi)
                            if ((point1.phi - point3.phi) > np.pi):
                                phi1 = point1.phi-(2*np.pi)
                            if ((point3.phi - point1.phi) > np.pi):
                                phi3 =  point3.phi-(2*np.pi)
                            if ((point2.phi - point3.phi) > np.pi):
                                phi2 = point2.phi-(2*np.pi)
                            if ((point3.phi - point2.phi) > np.pi):
                                phi3 =  point3.phi-(2*np.pi)
                            tripletsphi.append([phi1, phi2, phi3])
    return np.array([np.array(tripletsz), np.array(tripletsphi)])

def second_deriv(triplets):
    term1 = triplets[:,:,2]/(10-15)/(5-15)
    term2 = triplets[:,:,1]/(15-10)/(5-10)
    term3 = triplets[:,:,0]/(15-5)/(10-5)
    return 2*(term1+term2+term3)

def ratio(data, z_cutoff = -3, phi_cutoff = -4):
    z = np.log10(np.abs(second_deriv(data)))[0]
    phi = np.log10(np.abs(second_deriv(data)))[1]
    z_ratio = 0
    phi_ratio = 0
    both_ratio = 0 
    for i in range(len(z)):
        if (z[i] <= z_cutoff):
            z_ratio += 1
        if(phi[i] <= phi_cutoff):
            phi_ratio += 1
        if (z[i] <= z_cutoff)&(phi[i] <= phi_cutoff):
            both_ratio += 1
    z_ratio = z_ratio/len(z)
    phi_ratio = phi_ratio/len(z)
    both_ratio = both_ratio/len(z)
    return [z_ratio, phi_ratio, both_ratio]

triplets_rand = create_triplets(random = True, anchor_layer = 3)
triplets_v2 = create_triplets('v2', anchor_layer = 3)
triplets_v3 = create_triplets('v3', anchor_layer = 3)

log_second_deriv_rand = np.log10(np.abs(second_deriv(triplets_rand)))
log_second_deriv_v2 = np.log10(np.abs(second_deriv(triplets_v2)))
log_second_deriv_v3 = np.log10(np.abs(second_deriv(triplets_v3)))

rand_ratios = ratio(triplets_rand)
v2_ratios = ratio(triplets_v2)
v3_ratios = ratio(triplets_v3)

bins = [np.arange(-11, 1, .3), np.arange(-12, -1, 0.3)]
plt.figure(figsize = (15, 10))

plt.subplot(2, 3, 1)
plt.hist(log_second_deriv_rand[0], edgecolor='black', rwidth=0.8, bins = bins[0], label = f'Ratio: {np.round(rand_ratios[0], 5)}')
plt.title('Random z')
plt.ylabel('Count')
plt.xlabel(r'$log_{10}$(|z"|)')
plt.legend(loc = 'upper left')
plt.yscale('log')
plt.axvline(-3, ls = '--', color = 'r')

plt.subplot(2, 3, 2)
plt.hist(log_second_deriv_v2[0], edgecolor='black', rwidth=0.8, bins =bins[0], label = f'Ratio: {np.round(v2_ratios[0], 5)}')
plt.title('v2 z')
plt.ylabel('Count')
plt.xlabel(r'$log_{10}$(|z"|)')
plt.legend(loc = 'upper left')
plt.yscale('log')
plt.axvline(-3, ls = '--', color = 'r')

plt.subplot(2, 3, 3)
plt.hist(log_second_deriv_v3[0], edgecolor='black', rwidth=0.8, bins = bins[0], label = f'Ratio: {np.round(v3_ratios[0], 5)}')
plt.title('v3 z')
plt.ylabel('Count')
plt.xlabel(r'$log_{10}$( |z"| )')
plt.legend(loc = 'upper left')
plt.yscale('log')
plt.axvline(-3, ls = '--', color = 'r')

plt.subplot(2, 3, 4)
plt.hist(log_second_deriv_rand[1], edgecolor='black', rwidth=0.8, bins = bins[1], label = f'Ratio: {np.round(rand_ratios[1], 5)}')
plt.title('Random phi')
plt.ylabel('Count')
plt.xlabel(r'$log_{10}$(|$\phi$"|)')
plt.legend(loc = 'upper left')
plt.yscale('log')
plt.axvline(-4, ls = '--', color = 'r')

plt.subplot(2, 3, 5)
plt.hist(log_second_deriv_v2[1], edgecolor='black', rwidth=0.8, bins = bins[1], label = f'Ratio: {np.round(v2_ratios[1], 5)}')
plt.title('v2 phi')
plt.ylabel('Count')
plt.xlabel(r'$log_{10}$(|$\phi$" |)')
plt.legend(loc = 'upper left')
plt.yscale('log')
plt.axvline(-4, ls = '--', color = 'r')

plt.subplot(2, 3, 6)
plt.hist(log_second_deriv_v3[1], edgecolor='black', rwidth=0.8, bins = bins[1], label = f'Ratio: {np.round(v3_ratios[1], 5)}')
plt.title('v3 phi')
plt.ylabel('Count')
plt.xlabel(r'$log_{10}$(|$\phi$" |)')
plt.legend(loc = 'upper left')
plt.yscale('log')
plt.axvline(-4, ls = '--', color = 'r')

ranges = [[-11, 1], [-12, -1]]
vbound = (-4 - ranges[1][0])/(ranges[1][1]-ranges[1][0])
hbound = (-3 - ranges[0][0])/(ranges[0][1]-ranges[0][0])
plt.figure(figsize = (17, 5))
plt.subplot(1, 3, 1)
plt.hist2d(log_second_deriv_rand[0],log_second_deriv_rand[1],range = ranges, bins = 25, norm=mpl.colors.LogNorm())
plt.title('Random')
plt.xlabel(r'$log_{10}$(|z"|)')
plt.ylabel(r'$log_{10}$(|$\phi$" |)')
plt.axvline(-3, 0, vbound, ls = '--', color = 'r')
plt.axhline(-4, 0, hbound, ls = '--', color = 'r')
plt.text(ranges[0][0]+0.5, ranges[1][1]-0.7, f'Ratio: {np.round(rand_ratios[2], 5)}', bbox=dict(boxstyle="square", fc = 'None'))

plt.subplot(1, 3, 2)
plt.hist2d(log_second_deriv_v2[0],log_second_deriv_v2[1],range = ranges, bins = 25, norm=mpl.colors.LogNorm())
plt.title('v2')
plt.xlabel(r'$log_{10}$(|z"|)')
plt.ylabel(r'$log_{10}$(|$\phi$" |)')
plt.axvline(-3, 0, vbound, ls = '--', color = 'r')
plt.axhline(-4, 0, hbound, ls = '--', color = 'r')
plt.text(ranges[0][0]+0.5, ranges[1][1]-0.7, f'Ratio: {np.round(v2_ratios[2], 7)}', bbox=dict(boxstyle="square", fc = 'None'))

plt.subplot(1, 3, 3)
plt.hist2d(log_second_deriv_v3[0],log_second_deriv_v3[1], range = ranges, bins = 25, norm=mpl.colors.LogNorm())
plt.title('v3')
plt.xlabel(r'$log_{10}$(|z"|)')
plt.ylabel(r'$log_{10}$(|$\phi$" |)')
plt.colorbar()
plt.axvline(-3, 0, vbound, ls = '--', color = 'r')
plt.axhline(-4, 0, hbound, ls = '--', color = 'r')
plt.text(ranges[0][0]+0.5, ranges[1][1]-0.7, f'Ratio: {np.round(v3_ratios[2], 7)}', bbox=dict(boxstyle="square", fc = 'None'))

plt.show()