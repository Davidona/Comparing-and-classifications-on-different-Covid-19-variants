# -- coding: utf-8 --
"""
Created on Mon Dec 26 20:42:26 2022

@author: Hodaya

"""

# Import Module
from tkinter import *
from tkinter import ttk,messagebox
from tkinter.ttk import Progressbar
from tkinter import filedialog
import time
from GenomeManipulation import *
import pickle
import numpy as np
import matplotlib.pyplot as plt
import threading
import os


def upKey():
    result = kMersEntry.get()
    x = int(result)
    x = x + 1
    kMersEntry.delete(0, END)
    kMersEntry.insert(END, x)


def downKey():
    result = kMersEntry.get()
    x = int(result)
    x = x - 1
    kMersEntry.delete(0, END)
    if x < 0:
        kMersEntry.delete(0, END)
        kMersEntry.insert(END, 0)
    else:
        kMersEntry.insert(END, x)

def checkInput():
    if ResultField.get() == '':
        messagebox.showerror(title="Error", message="Please add destination folder")
        return False
    
    if covidGenomField.get() == '':
        messagebox.showerror(title="Error", message="Please add genomes file path")
        return False
    
    if int(kMersEntry.get()) < 1:
        messagebox.showerror(title="Error", message="K-mer must be above 0")
        return False
    if int(kMersEntry.get()) > 2000:
        messagebox.showerror(title="Error", message="K-mer must be below 2000")
        return False
    if not (os.path.exists(ResultField.get()) and os.path.exists(ResultField.get())):
        messagebox.showerror(title="Error", message="path does not exist")
        return False
    return True
    
    
    
def exitFunc():
    root.destroy()


def browseGenomes():
    filename = filedialog.askopenfilename(
        initialdir="/",
        title="Select a File",
        filetypes=(("Text files", "*.txt*"), ("all files", "*.*")),
    )
    covidGenomField.delete(0, "end")
    covidGenomField.insert(0, filename)


def start_thread(event):
    global Agenda_thread
    calculate_thread = threading.Thread(target=analayze)
    calculate_thread.daemon = True
    calculate_thread.start()


def getFolderPath():
    folder_selected = filedialog.askdirectory()
    ResultField.delete(0, "end")
    ResultField.insert(0, folder_selected)


def analayze():
    if checkInput():
        windowSize = int(kMersEntry.get())
        destination = ResultField.get()
        gm = GenomeManipulation(windowSize)
        covidList, imposterList, = gm.importGenomes(covidGenomField.get())
        combinationNames = gm.CreateImpostersCombinationNames(imposterList)

        listofVlists, allDicOfImposters = gm.calculateDistances(
            covidList, imposterList, pb3, pb3txt)
        
        DistanceMatrices = gm.generateGenomeDistanceMatrix(covidList, imposterList, allDicOfImposters)
        gm.plotAndSaveDistanceMatrices(DistanceMatrices, covidList, destination)
        gm.clusteringData(listofVlists, covidList, combinationNames, destination)
        pb3["value"] = 100 # this line is to edit the progress bar in GUI
        pb3txt["text"] = pb3["value"], "%" #this line is to edit the progress percentage in GUI
        messagebox.showinfo(title="Job finished", message="Job finished please check Destination folder")



if __name__ == "__main__":
    # create root window
    root = Tk()

    # root window title and dimension
    root.title("Genome Classification")
    # Set geometry(widthxheight)
    root.geometry("600x500")

    # Init fields
    covidGenomField = Entry(root, width=40)
    covidGenomField.place(x=100, y=59)

    ResultField = Entry(root, width=40)
    ResultField.place(x=100, y=159)

    ##### K-Mers

    kMersLabel = Label(root, width=20, text="K-Mers")
    kMersLabel.place(x=170, y=230)
    kMersEntry = Entry(root, width=5)
    kMersEntry.insert(0, "0")
    kMersEntry.place(x=270, y=230)


    kMersUp = Button(
        root, text="▲", fg="black", font=("arial", "5", "bold italic"), command=upKey
    )
    kMersDown = Button(
        root, text="▼", fg="black", font=("arial", "5", "bold italic"), command=downKey
    )
    kMersUp.place(x=300, y=230)
    kMersDown.place(x=300, y=240)

    ##Functions

    ######
    # General buttons
    closeButton = Button(
        root,
        text="Close",
        fg="black",
        font=("arial", "10", "bold italic"),
        command=exitFunc,
    )
    closeButton.place(x=520, y=440)

    ### Progress bars

    ### First bar

    # styleFirst = ttk.Style()
    # styleFirst.theme_use('clam')
    # styleFirst.configure("red.Horizontal.TProgressbar", foreground='white', background='orange')

    btn = Button(
        root,
        text="Browse Genomes",
        fg="black",
        font=("arial", "10", "bold italic"),
        command=browseGenomes,
    )
    btn.place(x=350, y=50)

    resultBtn = Button(
        root,
        text="Select Results Folder",
        fg="black",
        font=("arial", "10", "bold italic"),
        command=getFolderPath,
    )
    resultBtn.place(x=350, y=150)

    ### Second bar

    ## - ## - ## - ##

    ###Third bar
    pb3 = Progressbar(root, orient=HORIZONTAL, length=100, mode="determinate")

    pb3.place(x=200, y=350)

    pb3txt = Label(
        root,
        text="0%",
        # bg = '#04b424',
        fg="#04b424",
    )

    pb3txt.place(x=300, y=350)

    analyzeBtn = Button(
        root,
        text="Analayze",
        fg="black",
        font=("arial", "10", "bold italic"),
        command=lambda: start_thread(None),
    )
    analyzeBtn.place(x=220, y=310)
    #####

    # Execute Tkinter
    root.mainloop()
