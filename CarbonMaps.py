import os.path
import tkfilebrowser
from SaveToHDF5 import load_from_hdf5
import numpy as np
import PySimpleGUI as sg
import matplotlib
import matplotlib.pyplot as plt

def carbon_maps():
    # tells matplotlib to use TkAgg as the backend GUI integration
    matplotlib.use("TkAgg")
    progdir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    dirsave = open(progdir + "\\" + "lastdirsave.txt", "r")
    startdirbros = dirsave.read()
    dirsave.close()

    try:
        file_chosen = tkfilebrowser.askopenfilename(title="Pick a .hdf5 file with the processed data",
                                                    okbuttontext="Select", initialdir=startdirbros, filetypes=[("HDF5", "*.hdf5")])
    except:
        file_chosen = ""

    startdirbros = file_chosen.rsplit("\\", 1)[0]
    dirsave = open(progdir + "\\" + "lastdirbros.txt", "w")
    dirsave.write(startdirbros)
    dirsave.close()

    if file_chosen != "":
        data = load_from_hdf5(file_chosen)
        snew = carbon_map_creator(data)
        if snew == "Error":
            return


def carbon_map_creator(snew, *args):
    energy = snew.eVenergy
    stack = snew.spectr
    subdim = np.ceil(np.sqrt(len(energy)))

    if len(snew.eVenergy) < 2:
        sg.popup("Too few images for this mapping routine")
        return "Error"

    test = energy[(energy < 319)*(energy > 277)]
    if len(test) == 0:
        sg.popup("This is not the carbon edge")
        return "Error"

    if len(args) == 0:
        spthresh = 0.35
        figsav = 0
        nofig = 0
    elif len(args) == 1:
        spthresh = args(0)
        figsav = 0
        nofig = 0
    elif len(args) == 2:
        spthresh = args(0)
        figsav = 0
        nofig = 1
    elif len(args) == 3:
        spthresh = args(0)
        figsav = 1
        nofig = 0
        rootdir = args(1)
        sample = args(2)

    if spthresh > 1:
        spthresh = spthresh/100

    if nofig == 0 and len(energy) > 6:
        plt.figure()

    return snew
