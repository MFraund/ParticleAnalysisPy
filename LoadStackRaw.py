import os.path
import numpy as np

class LoadStackRaw:
    def __init__(self, filedir):
        os.chdir(filedir)
        filestruct = os.listdir(filedir)
        self.particle = "" # for future use

        firstrun = 0

        for i in range(len(filestruct)):
            stridx = filestruct[i].find("xim")
            hdridx = filestruct[i].find("hdr")

            if stridx != -1:
                if firstrun == 0:
                    self.spectr = np.flipud(np.loadtxt(filestruct[i], dtype=int))
                    firstrun = 1
                else:
                    self.spectr = np.dstack((self.spectr, np.flipud(np.loadtxt(filestruct[i], dtype=int))))

            elif hdridx != -1:
                self.eVenergy, self.Xvalue, self.Yvalue = read_hdr(filestruct[i])

        # corrects for np.dstack's incorrect ordering of dimensions
        self.spectr = np.swapaxes(self.spectr, 0, 2)
        self.spectr = np.swapaxes(self.spectr, 1, 2)

        if self.spectr.shape[0] < self.eVenergy.size:
            self.eVenergy = np.delete(self.eVenergy, np.s_[self.spectr.shape[0]:self.eVenergy.size], 0)

        self.izero = np.atleast_2d(np.zeros((len(self.eVenergy), 2))) # for future use


def read_hdr(file):
    filestream = open(file, "r")
    cnt = 1

    evenergy = []
    xvalue = []
    yvalue = []

    # with open(file, "r") as filestream:
    #     for line in filestream:
    line = filestream.readline()
    while line != "":
        energypos = line.find("StackAxis")
        pos3 = line.find("; XRange =")

        if energypos != -1:
            line = filestream.readline()
            pos1 = line.find("(")
            pos2 = line.find(")")
            energyvec = np.array(line[pos1 + 1:pos2].split(', ')).astype(np.float)
            evenergy = np.atleast_2d(energyvec[1:len(energyvec)]).conj().T
        elif pos3 != -1:
            pos4 = line.find("; YRange =")
            pos5 = line.find("; XStep =")
            xvalue = float(line[pos3 + 10:pos4])
            yvalue = float(line[pos4 + 10:pos5])

        cnt = cnt + 1
        line = filestream.readline()

    filestream.close()

    return evenergy, xvalue, yvalue
