import numpy as np
import h5py
import os
from datetime import datetime

def save_to_hdf5(data, directory, filename):
    fullfilename = directory + "\\" + filename + ".hdf5"
    if os.path.exists(fullfilename):
        os.remove(fullfilename)

    file = h5py.File(fullfilename, 'w')

    file.create_group("Mantis")
    exchange_grp = file.create_group("exchange")
    file.create_dataset("implements", data="information:exchange:spectromicroscopy:Mantis".encode("ascii"))
    info_grp = file.create_group("information")
    file.create_group("spectromicroscopy")
    file.create_dataset("version", data="1.0.".encode("ascii"))

    adjust_data = np.swapaxes(data.spectr, 0, 2)
    adjust_data = np.swapaxes(adjust_data, 0, 1)
    data_set = exchange_grp.create_dataset("data", data=adjust_data)
    data_set.attrs["axes"] = "x:y:energy".encode("ascii")
    data_set.attrs["signal"] = 1

    energy_set = exchange_grp.create_dataset("energy", data=data.eVenergy)
    energy_set.attrs["units"] = "ev".encode("ascii")

    exchange_grp.create_dataset("x", data=data.Xvalue)
    exchange_grp.create_dataset("y", data=data.Yvalue)

    info_grp.create_dataset("comment", data="Converted in Mantis".encode("ascii"))
    info_grp.create_dataset("file_creation_datetime", data=datetime.now().strftime("%Y-%m-%dT%H:%M").encode("ascii"))

    file.create_dataset("particle", data=data.particle.encode("ascii"))

def load_from_hdf5(filename):
    file = h5py.File(filename, 'r')

    spectr = file["exchange/data"]
    spectr = np.swapaxes(spectr, 0, 2)
    spectr = np.swapaxes(spectr, 1, 2)
    evenergy = np.array(file["exchange/energy"])
    xvalue = float(np.array(file["exchange/x"]))
    yvalue = float(np.array(file["exchange/y"]))
    particle = file["particle"][()].decode("utf-8")

    data = QuickDataStruct(spectr, evenergy, xvalue, yvalue, particle)

    return data

class QuickDataStruct:
    def __init__(self, spectr, evenergy, xvalue, yvalue, particle):
        self.spectr = spectr
        self.eVenergy = evenergy
        self.Xvalue = xvalue
        self.Yvalue = yvalue
        self.particle = particle

