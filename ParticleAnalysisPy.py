# This is a test script: It is a Python 3.9 conversion for the Particle Analysis Matlab program
# This file is where the basic GUI will be developed and will call the functions that do the actual work

import PySimpleGUI as sg
from ProcessStackData import process_stack_data
from CarbonMaps import carbon_maps

layout = [
        [sg.Text("Particle Analysis Py", justification="center", font=(None, 15, ["bold", "italic"]))],
        [sg.Frame("Raw Data Processing", [
                [sg.Button("Process Stack Data", expand_x=True)],
                [sg.Button("Stack From Images", expand_x=True)],
                [sg.Button("Deglitch Stack", expand_x=True)]
        ], element_justification="center", expand_x=True, title_location=sg.TITLE_LOCATION_TOP)],

        [sg.Frame("Processed Data Analysis", [
                [sg.Button("DirLabelStruct", expand_x=True)],
                [sg.Button("StackLab", expand_x=True)],
                [sg.Button("CarbonMaps", expand_x=True)],
                [sg.Button("FullSpecCMaps", expand_x=True)]
        ], element_justification="center", expand_x=True, title_location=sg.TITLE_LOCATION_TOP)],

        [sg.Frame("Processed Data Analysis", [
                [sg.Button("Explore Particles", expand_x=True)],
                [sg.Button("Remove Bad Particles", expand_x=True)]
        ], element_justification="center", expand_x=True, title_location=sg.TITLE_LOCATION_TOP)]
]

window = sg.Window("Particle Analysis Python Edition V1.0", layout, element_justification= "center")

while True:
    event, values = window.read()
    # end program if window is closed
    if event == "Exit" or event == sg.WIN_CLOSED:
        break

    if event == "Process Stack Data":
        process_stack_data()

    if event == "Stack From Images":
        sg.popup("Work in Progress")

    if event == "Deglitch Stack":
        sg.popup("Work in Progress")

    if event == "DirLabelStruct":
        sg.popup("Work in Progress")

    if event == "StackLab":
        sg.popup("Work in Progress")

    if event == "CarbonMaps":
        carbon_maps()

    if event == "FullSpecCMaps":
        sg.popup("Work in Progress")

    if event == "Explore Particles":
        sg.popup("Work in Progress")

    if event == "Remove Bad Particles":
        sg.popup("Work in Progress")

window.close()
