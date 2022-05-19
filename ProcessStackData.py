# This is a test script: It is a Python 3.9 conversion for the Particle Analysis Matlab program
# This file is where the basic GUI will be developed and will call the functions that do the actual work

import PySimpleGUI as sg
import os.path
import matplotlib
import matplotlib.pyplot as plt
from SingStackProcGUI import sing_stack_proc_gui
from SingStackProcGUI import stack_movie
import tkfilebrowser
from SaveToHDF5 import save_to_hdf5


def process_stack_data():    
    
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
                values=[], enable_events=True, expand_x=True, size=(1, 10), key="-FILE LIST-"
            )
        ],
        [
            sg.Button("Process", key="-PROC BUTT-")
        ],
    ]
    
    figure_quest_form = [
            [sg.Text("Do you want to supress the figures?")],
            [sg.Button("Yes"),
             sg.Button("No")],
            [sg.Checkbox("Include animation", key="-CHECK BOX-")]
    ]
    
    savedatalay = [
        [
            sg.Text("Would you like to save your data?"),
            sg.Button("Yes"),
            sg.Button("No")
        ]
    ]
    
    # generate a window with the above layout
    window = sg.Window("Particle Analysis Python Edition V1.0", layout)
    
    rawstackdir = []
    rawstackname = []
    data_list = []
    ani = []
    plotfig = 0
    aniflag = 0
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
                file_list = list(
                    tkfilebrowser.askopendirnames(title="Choose Directories", okbuttontext="Add", initialdir=startdirbros))
            except:
                file_list = []
    
            rawstackdir.extend(file_list)
            for i in range(len(rawstackdir)):
                rawstackname.extend([rawstackdir[i].rsplit("\\", 1)[1]])
    
            if len(rawstackname) > 0:
                window["-FILE LIST-"].update(rawstackname)
                startdirbros = rawstackdir[-1].rsplit("\\", 1)[0]
                dirsave = open(progdir + "\\" + "lastdirbros.txt", "w")
                dirsave.write(startdirbros)
                dirsave.close()
                rawstackname = []
    
        if event == "-PROC BUTT-":
            if not rawstackdir:
                sg.popup("No Folders Selected")
            else:

                questwin = sg.Window("Particle Analysis Python Edition V1.0", figure_quest_form)

                while True:
                    questevent, questval = questwin.read()
                    if questevent == "Exit" or questevent == sg.WIN_CLOSED:
                        break
                    if questevent == "Yes":
                        plotfig = 0
                        questwin.close()
                    if questevent == "No":
                        plotfig = 1
                        questwin.close()
                    if questval["-CHECK BOX-"]:
                        aniflag = 1
                    if not questval["-CHECK BOX-"]:
                        aniflag = 0

                questwin.close()

                for i in rawstackdir:
                    data_list.append(sing_stack_proc_gui(i, plotfig))
                    if aniflag == 1 and plotfig == 1:
                        ani.append(stack_movie(data_list[-1]))
                        plt.show(block=False)

                savewin = sg.Window("Particle Analysis Python Edition V1.0", savedatalay)
                while True:
                    saveevent, saveval = savewin.read()
                    dirsave = open(progdir + "\\" + "lastdirsave.txt", "r")
                    startdirsave = dirsave.read()
                    dirsave.close()

                    if saveevent == "Yes":
                        try:
                            file_save = tkfilebrowser.askopendirname(title="Chose directory to save to", okbuttontext="Save",
                                                                     initialdir=startdirsave)
                        except:
                            file_save = ""

                        if len(file_save) > 0:
                            for i in range(len(data_list)):
                                save_to_hdf5(data_list[i], file_save, rawstackdir[i].rsplit("\\", 1)[1])
                            startdirsave = file_save
                            dirsave = open(progdir + "\\" + "lastdirsave.txt", "w")
                            dirsave.write(startdirsave)
                            dirsave.close()
                            savewin.close()
                            window.close()
                    if saveevent == "No":
                        savewin.close()
                        window.close()
                    if saveevent == "Exit" or saveevent == sg.WIN_CLOSED:
                        break
    
    window.close()
