import PySimpleGUI as sg
import tkfilebrowser as tk
import os.path
from SaveToHDF5 import append_dir_label_to_hdf5
from SaveToHDF5 import load_from_hdf5
import numpy as np
from CarbonMaps import carbon_map_creator
from SootCarboxSizeHist import soot_carbox_size_hist

def dir_label_struct():
    progdir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    rawstackdir = []
    rawstackname = []
    data_list = []

    layout = [
        [
            sg.Text("sp2"),
            sg.Push(),
            sg.Input(default_text="0", justification="right", size=10, key="-SP2 BOX-")
        ],
        [
            sg.Text("Dataset Name"),
            sg.Push(),
            sg.Input(justification="right", size=10, key="-DATASET BOX-")
        ],
        [
            sg.Text("Cropped Image Type"),
            sg.Push(),
            sg.Combo(["RGB", "SVD"], default_value="RGB", readonly=True, key="-TYPE BOX-")
        ],
        [
            sg.Button(button_text="GO!", key="-GO BUTT-")
        ]
    ]

    window = sg.Window("Particle Analysis Python Edition V1.0", layout)

    while True:
        event, values = window.read()
        if event == "Exit" or event == sg.WIN_CLOSED:
            break

        if event == "-GO BUTT-":
            sp2 = values["-SP2 BOX-"]
            dataset = values["-DATASET BOX-"]
            imagetype = values["-TYPE BOX-"]

            window.close()
            break

    dirsave = open(progdir + "\\" + "lastdirbros.txt", "r")
    startdirbros = dirsave.read()
    dirsave.close()

    # setup formatting for left side of GUI
    filelayout = [
        [
            sg.Push(),
            sg.Text("Choose Directories to Process"),
            sg.Button("Browse", key="-BROS BUTT-")

        ],
        [
            sg.Listbox(values=[], enable_events=True, expand_x=True, size=(1, 10), key="-FILE LIST-")
        ],
        [
            sg.Button("Process", key="-PROC BUTT-")
        ],
    ]

    # generate a window with the above layout
    filewindow = sg.Window("Particle Analysis Python Edition V1.0", filelayout)

    # event loop for GUI
    while True:
        event, values = filewindow.read()
        # end program if window is closed
        if event == "Exit" or event == sg.WIN_CLOSED:
            break

        if event == "-BROS BUTT-":
            dirsave = open(progdir + "\\" + "lastdirbros.txt", "r")
            startdirbros = dirsave.read()
            dirsave.close()
            try:
                qcdir = list(
                    tk.askopendirnames(title="Choose Directories to Process", okbuttontext="Add", initialdir=startdirbros))
            except:
                qcdir = []

            rawstackdir.extend(qcdir)
            for i in range(len(rawstackdir)):
                rawstackname.extend([rawstackdir[i].rsplit("\\", 1)[1]])

            if len(rawstackname) > 0:
                filewindow["-FILE LIST-"].update(rawstackname)
                startdirbros = rawstackdir[-1].rsplit("\\", 1)[0]
                dirsave = open(progdir + "\\" + "lastdirbros.txt", "w")
                dirsave.write(startdirbros)
                dirsave.close()
                rawstackname = []

        if event == "-PROC BUTT-":
            if not rawstackdir:
                sg.popup("No Folders Selected")
            else:
                while True:
                    dirsave = open(progdir + "\\" + "lastdirsave.txt", "r")
                    startdirsave = dirsave.read()
                    dirsave.close()
                    try:
                        mappath = tk.askopendirname(title="Chose a place to save output", okbuttontext="Save", initialdir=startdirsave)
                    except:
                        mappath = ""

                    if len(mappath) > 0:
                        print('Processing...')
                        sampleid = rawstackname.copy()

                        dirlabels = dir_label_maps_struct(qcdir, dataset, 0, mappath, sp2, imagetype)

                        for i in range(len(data_list)):
                            append_dir_label_to_hdf5(data_list[i], mappath, rawstackdir[i].rsplit("\\", 1)[1])
                        startdirsave = mappath
                        dirsave = open(progdir + "\\" + "lastdirsave.txt", "w")
                        dirsave.write(startdirsave)
                        dirsave.close()
                        filewindow.close()
                        print('Complete!')

    filewindow.close()
    window.close()

    return

def dir_label_maps_struct(indir, sample, saveflg, savepath, sp2, impar):
    # Modified from DirLabelMapsStruct.m a matlab program created by Ryan Moffet 7/7/2016
    # Modified version by TJJ 6/27/2022

    labelcnt = np.zeros((4, len(indir)))
    stackcnt = 0
    mapcnt = 0

    # initialize values
    partsize = [] # particle sizes
    label = [] # particle labels( as strings)
    cmpsiz = [] # areas of components
    sootcarbox = []
    totalcarbon = [] # OD of total carbon(OD(320) - OD(278))
    carbox = []
    sp2 = []
    sootdistcent = [] # relative distance of the soot inclusion from the particle center
    sootdistcentinscribed = []
    sootecc = [] # eccentricity of the soot inclusion
    sootmaj = [] # major axis of soot inclusion
    sootmin = [] # minor axis of the soot inclusion
    sootcvex = [] # convexity of the soot inclusion
    sootarea = []
    croppedparts = [] # Cropped RGB Images of particles
    imageprops = [] # image properties: [Xvalue, Yvalue,  # of X pixels,# of Y pixels]
    partdirs = []
    partsn = []

    # begin loop over sample directories
    for j in range(len(indir)):
        directory = os.listdir(indir[j])
        loopctr = 1
        totpix = 0

        # loop over stacks
        for l in range(len(directory)):
            ind = directory[l].find(".hdf5")
            # if the directory has a hdf5 file...
            if ind != -1:
                file_name = indir[j] + "\\" + directory[l]
                snew, dirdata = load_from_hdf5(file_name)
                print(file_name)
                # Error Checking
                if len(snew.eVenergy) < 3:
                    print("Too few images for CarbonMaps")
                    continue
                test = snew.eVenergy[(snew.eVenergy < 325) & (snew.eVenergy > 275)]
                if len(test) == 0:
                    print("This is not a carbon edge")
                    continue
                if np.amax(test) < 315:
                    print("No post edge, stack skipped")
                    continue
                if len(test) > 8:
                    stackcnt = stackcnt + 1
                else:
                    mapcnt = mapcnt + 1
                print("# of Stacks = " + str(stackcnt) + ", # of Maps = " + str(mapcnt))

                # Run DiffMaps and label particles
                if saveflg == 1:
                    sinp = carbon_map_creator(snew, sp2, savepath, sample[j])
                    sinp = soot_carbox_size_hist(sinp, 1)
                else:
                    sinp = carbon_map_creator(snew, sp2, 1)
                    sinp = soot_carbox_size_hist(sinp, 1)
                if sinp.size() == 0:
                    print('no particles identified, stack skipped')
                    continue

                # count particle classes
                newcount = label_count(sinp)
                labelcnt[:, j] = labelcnt[:, j] + newcount

                # find soot inclusion distance from the center


    return dirlabels

def label_count(sin):
    occnt = 0
    oceccnt = 0
    ocincnt = 0
    ocecincnt = 0

    for i in len(sin.partlabel):
        ecidx = sin.partlabel[i].find('EC')
        inidx = sin.partlabl[i].find('In')
        kidx = sin.partlabel[i].find('K')
        if sin.partlabel[i] == 'OC':
            occnt = occnt + 1
        elif ecidx != -1 and (inidx != -1 or kidx != -1):
            ocecincnt = ocecincnt + 1
        elif ecidx != -1 and (inidx == -1 or kidx == -1):
            oceccnt = oceccnt + 1
        elif ecidx == -1 and (inidx != -1 or kidx != -1):
            ocincnt = ocincnt + 1

    out = [occnt, oceccnt, ocecincnt, ocincnt]

    return out

def dist_cent(sin):
    nx, ny = np.shape(sin.labelmat)
    totimage = np.nansum(sin.spectr, axis=2)
    totimage[totimage == np.inf] = 0
    sootmap = sin.maps[2, : , :]
    mask = np.zeros((nx, ny))
    rad, cs, cy =

    return sout

def max_inscribed_circle(contourimage, display)
    # get contour
    sz = contourimage.shape()
    y, x = np.nonzero(contourimage == 255)
    y = y[0]
    x = x[0]
    contour =

    return r, cx, cy
