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
    file.create_dataset("implements", data="information:exchange:spectromicroscopy:Mantis".encode("ascii"),)
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
