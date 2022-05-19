import PySimpleGUI as sg
from ProcessStackData import process_stack_data

layout = [
        [sg.Text("Particle Analysis Py", justification="center", font=(None, 15, ["bold", "italic"]))],
        [sg.Frame("Raw Data Processing", [
                [sg.Button("Process Stack Data", expand_x=True)],
                [sg.Button("Stack From Images", expand_x=True)],
                [sg.Button("Process Stack Data", expand_x=True)]
        ], element_justification="center", expand_x=True, title_location=sg.TITLE_LOCATION_TOP)],

        [sg.Frame("Processed Data Analysis", [
                [sg.Button("DirLabelStruct", expand_x=True)],
                [sg.Button("StackLab", expand_x=True)],
                [sg.Button("CarbonMaps", expand_x=True)],
                [sg.Button("FullSpecCMaps", expand_x=True)]
        ], element_justification="center", expand_x=True, title_location=sg.TITLE_LOCATION_TOP)],

        [sg.Frame("Processed Data Analysis", [
                [sg.Button("ExploreParticles", expand_x=True)],
                [sg.Button("RemoveBadParticles", expand_x=True)]
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

window.close()
