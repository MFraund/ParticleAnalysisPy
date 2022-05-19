# This is a test script: It is a Python 3.9 conversion for the Particle Analysis Matlab program
# This file is where the basic GUI will be developed and will call the functions that do the actual work

import PySimpleGUI as sg
import os.path
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
import matplotlib.pyplot as plt
from SingStackProcGUI import sing_stack_proc_gui
from SingStackProcGUI import stack_movie
import tkfilebrowser
from SaveToHDF5 import save_to_hdf5

# function that will draw plots onto PySimpleGUI's Canvas
def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=1)

    return figure_canvas_agg


# tells matplotlib to use TkAgg as the backend GUI integration
matplotlib.use("TkAgg")

# setup formatting for left side of GUI
layout = [
    [
        sg.Push(),
        sg.Text("Choose Directories to Process"),
        sg.Button("Browse", key="-BROS BUTT-")

    ],
    [
        sg.Listbox(
            values=[], enable_events=True, size=(60, 20), key="-FILE LIST-"
        )
    ],
    [
        sg.Button("Process", key="-PROC BUTT-")
    ],
]

figure_quest_form = [
    [
        sg.Text("Do you want to supress the figures?"),
        sg.Button("Yes"),
        sg.Button("No"),
        sg.Checkbox("Include animation:", key="-CHECK BOX-")
    ]
]

SaveDataLay = [
    [
        sg.Text("Would you like to save your data?"),
        sg.Button("Yes"),
        sg.Button("No")
    ]
]

# generate a window with the above layout
window = sg.Window("Particle Analysis Python Edition V1.0", layout)

rawStackDir = []
rawStackName = []
data_list = []
ani = []
dirlen2 = 0
plotFig = 0
aniflag = 0
startdirbros = os.getcwd()
startdirsave = os.getcwd()
progdir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

# event loop for GUI
while True:
    event, values = window.read()
    # end program if window is closed
    if event == "Exit" or event == sg.WIN_CLOSED:
        break

    if event == "-BROS BUTT-":
        dirsave = open(progdir + "\\" + "lastdirbros.txt", "r")
        startdirbros = dirsave.read()
        dirsave.close()
        try:
            file_list = list(tkfilebrowser.askopendirnames(title="Choose Directories", okbuttontext="Add", initialdir=startdirbros))
        except:
            file_list = []

        rawStackDir.extend(file_list)
        for i in range(len(rawStackDir)):
            rawStackName.extend([rawStackDir[i].rsplit("\\", 1)[1]])

        if len(rawStackName) > 0:
            window["-FILE LIST-"].update(rawStackName)
            startdirbros = rawStackDir[-1].rsplit("\\", 1)[0]
            dirsave = open(progdir + "\\" + "lastdirbros.txt", "w")
            dirsave.write(startdirbros)
            dirsave.close()
            rawStackName = []

    if event == "-PROC BUTT-":
        questWin = sg.Window("Particle Analysis Python Edition V1.0", figure_quest_form)

        while True:
            questEvent, questVal = questWin.read()

            if questEvent == "Yes":
                plotFig = 0
                questWin.close()
            if questEvent == "No":
                plotFig = 1
                questWin.close()
            if questWin["-CHECK BOX-"]:
                aniflag = 1
            if not questWin["-CHECK BOX-"]:
                aniflag = 0

            if questEvent == "Exit" or questEvent == sg.WIN_CLOSED:
                break

        questWin.close()

        for i in rawStackDir:
            data_list.append(sing_stack_proc_gui(i, plotFig))
            if aniflag == 1 and plotFig == 1:
                ani.append(stack_movie(data_list[-1]))
                plt.show(block=False)

        saveWin = sg.Window("Particle Analysis Python Edition V1.0", SaveDataLay)
        while True:
            saveEvent, saveVal = saveWin.read()
            dirsave = open(progdir + "\\" + "lastdirsave.txt", "r")
            startdirsave = dirsave.read()
            dirsave.close()

            if saveEvent == "Yes":
                try:
                    file_save = tkfilebrowser.askopendirname(title="Chose directory to save to", okbuttontext="Save",
                                                             initialdir=startdirsave)
                except:
                    file_save = ""

                if len(file_save) > 0:
                    for i in range(len(data_list)):
                        save_to_hdf5(data_list[i], file_save, rawStackDir[i].rsplit("\\", 1)[1])
                    startdirsave = file_save
                    dirsave = open(progdir + "\\" + "lastdirsave.txt", "w")
                    dirsave.write(startdirsave)
                    dirsave.close()
                    saveWin.close()
                    window.close()
            if saveEvent == "No":
                saveWin.close()
                window.close()
            if saveEvent == "Exit" or saveEvent == sg.WIN_CLOSED:
                break

window.close()
