import PySimpleGUI as sg
import os.path
import numpy as np
from LoadStackRaw import LoadStackRaw
from DeglitchStack import deglitch_stack
from AlignStack import align_stack
from ODStack import od_stack
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def sing_stack_proc_gui(datafolder, plotfig, varargin=None):
    if varargin is None:
        varargin = []
    bade = []
    os.chdir(datafolder)
    foldstruct = os.listdir(datafolder)
    print("Processing ", foldstruct[0])
    if len(varargin):
        evals = ""
        vlayout = [
            [sg.Text("Enter energies to remove separated by a space:")],
            [sg.InputText(enable_events=True, key="-E ENTRY-")],
            [sg.Button("Ok"), sg.Button("Cancel")]
        ]
        vwindow = sg.Window('Sample', vlayout)
        while True:
            vevent, vvalues = vwindow.read()
            if vevent == sg.WIN_CLOSED or vevent == "Cancel":
                break
            if vevent == "-E ENTRY-":
                evals = vvalues["-E ENTRY-"]
            if vevent == "Ok" and len(evals):
                bade = [int(i) for i in evals.split()]
                vwindow.close()
        vwindow.close()
    else:
        bade = []

    numobj = len(foldstruct)

    s = LoadStackRaw(datafolder)
    sdim = s.spectr.shape
    if sdim[0] > len(s.eVenergy):
        s.spectr = np.delete(s.spectr, np.s_[len(s.eVenergy):sdim[0]], 0)

    for i in range(numobj):
        bidx = foldstruct[i].find("hdr")
        if bidx != -1:
            s.particle = foldstruct[i].removesuffix(".hdr")

    s = deglitch_stack(s, bade)
    s = align_stack(s)

    if len(s.eVenergy) < 10:
        snew = od_stack(s, 'map', plotfig)
    else:
        snew = od_stack(s, 'O', plotfig)

    return snew

def stack_movie(stack):
    fig2 = plt.figure()
    ims = []
    for i in range(len(stack.eVenergy)):
        im = plt.imshow(update_frame(i, stack), animated=True)
        ims.append([im])
    ani = animation.ArtistAnimation(fig2, ims, interval=50, blit=True,
                                    repeat_delay=1000)
    return ani

def update_frame(iev, stack):
    imdat = stack.spectr[int(iev)].copy()
    return imdat
