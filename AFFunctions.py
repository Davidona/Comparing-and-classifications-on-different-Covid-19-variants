from collections import Counter
from scipy.spatial import distance
import numpy as np
import utils as ut

class DistanceCalculator:
    countsOfX = None
    countsOfY = None
    
    def EuclideanDistance(a, b): 
        NParrayA = np.array(a)
        NParrayB = np.array(b)
        return np.linalg.norm(NParrayA-NParrayB)
    
    # an alignment free function using the Euclidean method
    def EuclideanAlignmentFree(self, x, y, k = 3):
        self.prepareXandY(x, y, k)
        dis = self.EuclideanDistance(list(self.countsOfX.values()), list(self.countsOfY.values()))
        return dis
        
    
    # an alignment free function using the Canberra method
    def CanberraAlignmentFree(self, x, y, k = 3):
        self.prepareXandY(x, y, k)
        dis =  distance.canberra(list(self.countsOfX.values()), list(self.countsOfY.values()))
        return dis

    def prepareXandY(self, x, y, k = 3):
        KmerFromX = ut.CreateListOfSlidingWindow(x, k)
        KmerFromY = ut.CreateListOfSlidingWindow(y, k)
        setsWithZeros = dict.fromkeys(set(KmerFromX).union(set(KmerFromY)),0)
        
        self.countsOfX = Counter(KmerFromX)
        self.countsOfX = setsWithZeros | self.countsOfX # requires python 3.9! createse union between two sets
        self.countsOfY = Counter(KmerFromY)
        self.countsOfY = setsWithZeros | self.countsOfY # requires python 3.9! createse union between two sets