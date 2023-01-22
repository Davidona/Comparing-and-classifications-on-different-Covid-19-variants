# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 18:55:09 2022

@author: david
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.text import Annotation
import matplotlib.patches as mpatches
from utils import createFolder


class Kmedoids:
    max_iterations=100
    def __init__(self, max_iterations):
        self.max_iterations = max_iterations
    def generate_kmedoids(self, distance_matrix, k):
        # choose k random medoids as the initial clusters
        medoids_idx = np.random.choice(len(distance_matrix), size=k, replace=False)
        medoids = distance_matrix[medoids_idx]
        
        for i in range(self.max_iterations):
            # assign each point to the cluster with the nearest medoid
            clusters = [[] for _ in range(k)]
            for j, point in enumerate(distance_matrix):
                distances = np.linalg.norm(point - medoids, axis=1)
                clusters[np.argmin(distances)].append(j)
            
            # update the medoids to be the mean of the points in each cluster
            new_medoids = np.zeros((k, distance_matrix.shape[1]))
            for j, cluster in enumerate(clusters):
                new_medoids[j] = np.mean(distance_matrix[cluster], axis=0)
            
            # check for convergence
            if np.array_equal(new_medoids, medoids):
                break
            
            medoids = new_medoids
        
        return clusters
    
    
    
        
        
    
    
    def plot_kmedoids(self, distance_matrix, clusters, title, pointNames, filePath):
        pointOffsets = {}
        patchList=[]
        # create a scatter plot of the points, colored by cluster
        colors = plt.cm.rainbow(np.linspace(0, 1, len(clusters)))
        for i, cluster in enumerate(clusters):
            x = distance_matrix[cluster, 0]
            y = distance_matrix[cluster, 1]
            c =colors[i].reshape(1,-1)
            plt.scatter(x, y, c=c)
    
            # add labels to the points
            for j, point in enumerate(cluster):
                data_key = mpatches.Circle(xy = 0.5,color = colors[i], label=str(point)+""+pointNames[point])
                patchList.append(data_key)
                xmin, xmax, ymin, ymax = plt.axis()
                if y[j] in pointOffsets.keys():
                    pointOffsets[y[j]] = pointOffsets[y[j]]+10
                else:
                    pointOffsets[y[j]] = 10
                offset = (0, pointOffsets[y[j]])  # add 10 pixels of vertical spacing between labels
                if y[j] > ymax/2:
                    offset=(0,offset[1]*-1)
                annotation = Annotation(str(point), xy=(x[j], y[j]), xytext=offset, textcoords='offset points', ha='center', va='center', bbox=dict(facecolor='none', edgecolor='none'))
                plt.gca().add_artist(annotation)
        plt.legend(handles=patchList, loc='upper center',bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=3)        
        plt.title(title)
        createFolder(filePath+'/Clustering')
        plt.savefig(filePath+'/Clustering/'+title.split(':')[0]+'.png',dpi=300,bbox_inches='tight')           
        plt.show()
        
def getYmaxpixels():
    ymin, ymax = plt.ylim()
    pixels_per_inch = plt.gcf().get_dpi()
    units_per_inch = plt.gca().yaxis.get_scale()
    if units_per_inch == 'linear':
        units_per_inch = 1.0
    pixels_per_unit = pixels_per_inch / float(units_per_inch)
    ymax_pixels = ymax * pixels_per_unit
    return ymax_pixels
