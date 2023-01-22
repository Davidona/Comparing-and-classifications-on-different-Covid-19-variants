# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 15:36:26 2023

@author: david
"""
import pickle
from AFFunctions import *
from utils import *
from multiprocessing import Pool
from fastdtw import fastdtw
import numpy as np
import matplotlib.pyplot as plt
from kmedois import *


class Genome:
    name = ''
    sequence = ''
    GK = None
    
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
    def divideIntoKmers(self, windowSize):
        self.GK = CreateListOfSlidingWindow(self.sequence, windowSize, int(windowSize / 2))
    
class GenomeManipulation:    
    virusName = "coronavirus"
    windowSize = 0
    def __init__(self, windowSize):
        self.windowSize = windowSize
        
    def importGenomes(self, filePath): # This function imports genomes and creates a lists of imposters and covid

        imposterList = []  # list to hold all imposterGenomes
        covidList = []  # list to hold all covids
    
        f = open(filePath, "r")
        lines = f.readlines()
        a = ""
        currentVirus = ""
        for line in lines:
            if line[0] == ">":
                if a != "":
                    sequence = a.replace(" ", "")
                    name = currentVirus.split(" ")[0]
                    virus = Genome(name, sequence)
                    virus.divideIntoKmers(self.windowSize)
                    if self.virusName in currentVirus:  # if genome contains coronavirus in the name 
                        covidList.append(virus)
                    else:
                        imposterList.append(virus)
                    a = ""
                currentVirus = line
            else:
                a += line.replace("\n", "")
        sequence = a
        name = currentVirus.split(" ")[0]
        virus = Genome(name, sequence)
        virus.divideIntoKmers(self.windowSize)
        if self.virusName in currentVirus:
            covidList.append(virus)
        else:
            imposterList.append(virus)
        return covidList, imposterList
    
    
    def CreateImpostersCombinationNames(self, imposters): #this function justs combine all imposter names to make pairs
        combinationNames = []  # list to hold imposter pairs names
        for i in range(len(imposters) - 1):
            for j in range(i + 1, len(imposters)):
                combinationNames.append(
                    "combination:" + imposters[i].name + "-" + imposters[j].name
                )
        return combinationNames
    
    
    def compareAgainstpairs(self, imposterList, GK, dicOfImposters): # this function used to compare covid against pairs of imposters
        print('compareAgainstpairs')
        VList = []
        for i in range(len(imposterList) - 1):
            for j in range(i + 1, len(imposterList)):
                v = []
                for kmer in GK:
                    # if the covid genome is closer to the first imposter we add 0 to the vector otherwise we add 1
                    if dicOfImposters[i][kmer] < dicOfImposters[j][kmer]: 
                        v.append(0) 
                    else:
                        v.append(1)
                VList.append(v)
        return VList
    
    
    def calculateDistances(self, covidList, imposterList, pb3, pb3txt):
        distanceCalc = DistanceCalculator();
        print('calculateDistances')
        listofVlists = []
        count = 0
        allDicOfImposters = []
        for covidGenome in covidList:
            dicOfImposters = {}
            dic = dict.fromkeys(covidGenome.GK, 0)
            pb3["value"] = count / len(covidList) * 100 # this line is to edit the progress bar in GUI
            pb3txt["text"] = pb3["value"], "%" #this line is to edit the progress percentage in GUI
            count += 1
            for i in range(len(imposterList)):
                dicWithVal = dic.copy()
                for kmer in list(dicWithVal.keys()):
                    # kmerLis = [kmer]*len(impKmerList)
                    Sum = 0
                    for impKmer in imposterList[i].GK:
                        Sum += distanceCalc.CanberraAlignmentFree(kmer, impKmer)
                    Sum /= len(imposterList[i].GK)
                    dicWithVal[kmer] = Sum
                dicOfImposters[i] = dicWithVal
            allDicOfImposters.append(dicOfImposters)
            listofVlists.append(self.compareAgainstpairs(imposterList, covidGenome.GK, dicOfImposters))
        return listofVlists, allDicOfImposters
    
    
    def generateGenomeDistanceMatrix(self, 
        covidList, imposterList, allDicOfImposters):
        print('generateGenomeDistanceMatrix')
        listofVlists = []
        for covidGenome in range(len(covidList)):
            VList = []
            for i in range(len(imposterList)):
                v = []
                for kmer in covidList[covidGenome].GK:
                    v.append(allDicOfImposters[covidGenome][i][kmer])
                VList.append(v)
            listofVlists.append(VList)
        return listofVlists
    
    
    def plotAndSaveDistanceMatrices(self, listofVlists, covidList, filePath):
        print('plotAndSaveDistanceMatrices')
        covidnames = extract_attribute(covidList, 'name')
        for i in range(len(listofVlists)):
            plt.plot(np.array(listofVlists[i]))
            createFolder(filePath + "/SimilarityMatrices")
            plt.savefig(
                filePath + "/SimilarityMatrices/" + covidnames[i].replace(">", "") + ".png",
                dpi=300,
            )
            plt.title(covidnames[i])
            plt.show()
            plt.imshow(np.array(listofVlists[i]))
            plt.title(covidnames[i])
            plt.colorbar()
            plt.savefig(
                filePath
                + "/SimilarityMatrices/"
                + covidnames[i].replace(">", "")
                + "Matrix.png",
                dpi=300,
            )
            plt.show()
    
    
    # this function create a 2D matrix that counts how many genome i was clustered with genome j and returns it as a matrix
    def createClustersMatrix(self, listofAllLists, clustersList): 
        print('createClustersMatrix')
        mat = [[0] * len(listofAllLists) for i in range(len(listofAllLists))]
    
        for i in range(len(listofAllLists)):
            for clusters in clustersList:
                for cluster in clusters:
                    if i in cluster:
                        for j in cluster:
                            mat[i][j] += 1
        return mat
    
    
    def clusteringData(self, listofAllLists, covidList, combinationNames, filePath):
        covidNames = extract_attribute(covidList, 'name')
        kmedoids = Kmedoids(100)
        print('clusteringData')
        count = 0
        clustersList = []
        for combIndex in range(len(listofAllLists[0])):
            lis = []
            distanceList = [[0] * len(listofAllLists) for i in range(len(listofAllLists))]
            for i in listofAllLists:
                lis.append(i[combIndex])
            averagelis = []
            # if not checkCurrentRatio(lis, (len(listofAllLists)/2)-1, 0.9):
            # continue
            for i in lis:
                averagelis.append(Moving_Average(i))
            for i in range(len(averagelis)):
                for j in range(len(averagelis)):
                    distance, path = fastdtw(averagelis[i], averagelis[j])
                    distanceList[i][j] = distance
            count += 1
            distanceMatrx = np.array(distanceList)
    
            cl = kmedoids.generate_kmedoids(distanceMatrx, 2)
    
            kmedoids.plot_kmedoids(
                distanceMatrx,
                cl,
                "Index "+ str((count - 1)) + " " + combinationNames[count - 1],
                covidNames,
                filePath,
            )
            clustersList.append(cl)
    
        print(count)
    
        mat = self.createClustersMatrix(listofAllLists, clustersList)
    
        #normalize matrix:
        mat = np.array(mat)
        Mnorm = mat / count
        Mnorm = 1 - Mnorm
        
        cl = kmedoids.generate_kmedoids(Mnorm, 2)
    
        kmedoids.plot_kmedoids(Mnorm, cl, "Final Clustering", covidNames, filePath)
        
    
