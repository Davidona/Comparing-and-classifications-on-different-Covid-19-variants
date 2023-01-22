  
import numpy as np
import os

def CreateListOfSlidingWindow(l, k,jump=1):
    return [l[i:i+k] for i in range(0,len(l)-k+1, jump)]   #l is for list, k is for silding window size
# Program to calculate moving average using numpy


def Moving_Average(arr, window_size = 3):
    i = 0
    # Initialize an empty list to store moving averages
    moving_averages = []
      
    # Loop through the array t o
    #consider every window of size 3
    while i < len(arr) - window_size + 1:
      
        # Calculate the average of current window
        window_average = round(np.sum(arr[
          i:i+window_size]) / window_size, 2)
          
        # Store the average of current
        # window in moving average list
        moving_averages.append(window_average)
          
        # Shift window to right by one position
        i += 1
    return moving_averages

def createFolder(filePath):
    if not os.path.exists(filePath):
      
    # if the filePath directory is not present 
    # then create it.
        os.makedirs(filePath)

def extract_attribute(obj_list, attr):
    return [getattr(obj, attr) for obj in obj_list]
