import numpy as np
import h5py
import os
from datetime import datetime

def save_spectr_to_hdf5(data, directory, filename):

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

    # stuff unique to particle analysis py
    file.create_dataset("Particle_Analysis_Py", data="particle_analysis_py_ver_1.0".encode("ascii"))

    try:
        file.create_dataset("particle", data=data.particle.encode("ascii"))
    except:
        pass

def append_carbon_to_hdf5(data, fullfilename):

    file = h5py.File(fullfilename, 'a')

    carbon_grp = file.create_group("Carbon_Map")
    carbon_grp.create_dataset("totc", data=data.totc)
    carbon_grp.create_dataset("sp2", data=data.sp2)
    carbon_grp.create_dataset("soot_eccentricity", data=data.sooteccentricity)
    carbon_grp.create_dataset("soot_major_axis_length", data=data.sootmajoraxislength)
    carbon_grp.create_dataset("soot_minor_axis_length", data=data.sootminoraxislength)
    carbon_grp.create_dataset("soot_convex_area", data=data.sootconvexarea)
    carbon_grp.create_dataset("soot_area", data=data.sootarea)

    carbon_grp.create_dataset("label_mat", data=data.labelmat)
    carbon_grp.create_dataset("part_label", data=data.partlabel)
    carbon_grp.create_dataset("part_sn", data=data.partsn)
    carbon_grp.create_dataset("binmap", data=data.binmap)
    carbon_grp.create_dataset("comp_size", data=data.compsize)
    carbon_grp.create_dataset("part_dirs", data=data.partdirs)

    carbon_grp.create_dataset("size", data=data.size)
    carbon_grp.create_dataset("rgb_comp_map", data=data.rgbcompmap)
    carbon_grp.create_dataset("maps", data=data.maps)
    carbon_grp.create_dataset("bin_comp_map", data=data.bincompmap)

def append_dir_label_to_hdf5(data, fullfilename):

    file = h5py.File(fullfilename, 'a')

    dir_label_s = file.create_group("DirLabelS")
    dir_label_s.create_dataset("label_cnt", data=data.labelcnt)
    dir_label_s.create_dataset("part_size", data=data.partsize)
    dir_label_s.create_dataset("label", data=data.label)
    dir_label_s.create_dataset("cmp_siz", data=data.cmpsiz)
    dir_label_s.create_dataset("soot_carbox", data=data.sootcarbox)
    dir_label_s.create_dataset("total_carbon", data=data.totalcarbon)
    dir_label_s.create_dataset("carbox", data=data.carbox)
    dir_label_s.create_dataset("sp2", data=data.sp2)
    dir_label_s.create_dataset("out_rad", data=data.outrad)
    dir_label_s.create_dataset("rad_std", data=data.radstd)
    dir_label_s.create_dataset("sing_rad", data=data.singrad)
    dir_label_s.create_dataset("soot_dist_cent", data=data.sootdistcent)
    dir_label_s.create_dataset("soot_dist_cent_inscribed", data=data.sootdistcentinscribed)
    dir_label_s.create_dataset("soot_ecc", data=data.sootecc)
    dir_label_s.create_dataset("soot_maj", data=data.sootmaj)
    dir_label_s.create_dataset("soot_min", data=data.sootmin)
    dir_label_s.create_dataset("soot_cvex", data=data.sootcvex)
    dir_label_s.create_dataset("soot_area", data=data.sootarea)
    dir_label_s.create_dataset("cropped_props", data=data.croppedparts)
    dir_label_s.create_dataset("image_props", data=data.imageprops)
    dir_label_s.create_dataset("parts_dirs", data=data.partsdirs)
    dir_label_s.create_dataset("part_sn", data=data.partsn)


def load_from_hdf5(filename):

    file = h5py.File(filename, 'r')

    spectr = np.array(file["exchange/data"])
    spectr = np.swapaxes(spectr, 0, 1)
    spectr = np.swapaxes(spectr, 0, 2)
    evenergy = np.array(file["exchange/energy"])
    xvalue = float(np.array(file["exchange/x"]))
    yvalue = float(np.array(file["exchange/y"]))

    try:
        particle = file["particle"][()].decode("utf-8")
    except:
        particle = ''

    # initialize values for carbon mapping
    totc = []
    sp2 = []
    ecc = []
    major = []
    minor = []
    cvex = []
    area = []
    labelmat = []
    partlabel = []
    partsn = []
    binmap = []
    compsize = []
    partdirs = []
    size = []
    rgbcompmap = []
    maps = []
    bincompmap = []
    labelcnt = []
    partsize = []
    label = []
    cmpsiz = []
    sootcarbox = []
    totalcarbon = []
    carbox = []
    outrad = []
    radstd = []
    singrad = []
    sootdistcent = []
    sootdistcentinscribed = []
    sootecc = []
    sootmaj = []
    sootmin = []
    sootcvex = []
    sootarea = []
    croppedparts = []
    imageprops = []
    partsdirs = []
    avsootcarb = []
    avtotc = []
    avcarbox = []
    avsp2 = []

    data = QuickDataStruct(spectr, evenergy, xvalue, yvalue, particle , totc, sp2, ecc, major, minor, cvex, area,
                           labelmat, partlabel, partsn, binmap, compsize, partdirs, size, rgbcompmap, maps, bincompmap,
                           labelcnt, partsize, label, cmpsiz, sootcarbox, totalcarbon, carbox, outrad, radstd, singrad,
                           sootdistcent, sootdistcentinscribed, sootecc, sootmaj, sootmin, sootcvex, sootarea, croppedparts,
                           imageprops, partsdirs, avsootcarb, avtotc, avcarbox, avsp2)

    return data

class QuickDataStruct:
    def __init__(self, spectr, evenergy, xvalue, yvalue, particle, totc, sp2, ecc, major, minor, cvex, area,
                 labelmat, partlabel, partsn, binmap, compsize, partdirs, size, rgbcompmap, maps, bincompmap,
                 labelcnt, partsize, label, cmpsiz, sootcarbox, totalcarbon, carbox, outrad, radstd, singrad,
                 sootdistcent, sootdistcentinscribed, sootecc, sootmaj, sootmin, sootcvex, sootarea, croppedparts,
                 imageprops, partsdirs, avsootcarb, avtotc, avcarbox, avsp2):

        self.spectr = spectr
        self.eVenergy = evenergy
        self.Xvalue = xvalue
        self.Yvalue = yvalue
        self.particle = particle
        self.totc = totc
        self.sp2 = sp2
        self.sooteccentricity = ecc
        self.sootmajoraxislength = major
        self.sootminoraxislength = minor
        self.sootconvexarea = cvex
        self.sootarea = area
        self.labelmat = labelmat
        self.partlabel = partlabel
        self.partsn = partsn
        self.binmap = binmap
        self.compsize = compsize
        self.partdirs = partdirs
        self.size = size
        self.size = rgbcompmap
        self.maps = maps
        self.bincompmap = bincompmap
        self.labelcnt = labelcnt
        self.partsize = partsize
        self.label = label
        self.cmpsiz = cmpsiz
        self.sootcarbox = sootcarbox
        self.totalcarbon = totalcarbon
        self.carbox = carbox
        self.outrad = outrad
        self.radstd = radstd
        self.singrad = singrad
        self.sootdistcent = sootdistcent
        self.sootdistcentinscribed = sootdistcentinscribed
        self.sootecc = sootecc
        self.sootmaj = sootmaj
        self.sootmin = sootmin
        self.sootcvex = sootcvex
        self.sootarea = sootarea
        self.croppedparts = croppedparts
        self.imageprops = imageprops
        self.partsdirs = partsdirs
        self.avsootcarb = avsootcarb
        self.avtotc = avtotc
        self.avcarbox = avcarbox
        self.avsp2 = avsp2

# below is for future expansion
#     totc = np.array(file["Carbon_Map/totc"])
#     sp2 = np.array(file["Carbon_Map/sp2"])
#     ecc = file["Carbon_Map/soot_eccentricity"]
#     major = file["Carbon_Map/soot_major_axis_length"]
#     minor = file["Carbon_Map/soot_minor_axis_length"]
#     cvex = file["Carbon_Map/soot_convex_area"]
#     area = file["Carbon_Map/soot_area"]
#
#     labelmat = file["Carbon_Map/label_mat"]
#     partlabel = file["Carbon_Map/part_label"]
#     partsn = file["Carbon_Map/part_sn"]
#     binmap = file["Carbon_Map/binmap"]
#     compsize = file["Carbon_Map/comp_size"]
#     partdirs = file["Carbon_Map/part_dirs"]
#
#     size = file["Carbon_Map/size"]
#     rgbcompmap = file["Carbon_Map/rgb_comp_map"]
#     maps = file["Carbon_Map/maps"]
#     bincompmap = file["Carbon_Maps/bin_comp_map"]


#     labelcnt = file["dir_label_s/label_cnt"]
#     partsize = file["dir_label_s/part_size"]
#     label = file["dir_label_s/label"]
#     cmpsiz = file["dir_label_s/cmp_siz"]
#     sootcarbox = file["dir_label_s/soot_carbox"]
#     totalcarbon = file["dir_label_s/total_carbon"]
#     carbox = file["dir_label_s/carbox"]
#     sp2_dir = file["dir_label_s/sp2"]
#     outrad = file["dir_label_s/out_rad"]
#     radstd = file["dir_label_s/rad_std"]
#     singrad = file["dir_label_s/sing_rad"]
#     sootdistcent = file["dir_label_s/soot_dist_cent"]
#     sootdistcentinscribed = file["dir_label_s/soot_dist_cent_inscribed"]
#     sootecc = file["dir_label_s/soot_ecc"]
#     sootmaj = file["dir_label_s/soot_maj"]
#     sootmin = file["dir_label_s/soot_min"]
#     sootcvex = file["dir_label_s/soot_cvex"]
#     sootarea = file["dir_label_s/soot_area"]
#     croppedparts = file["dir_label_s/cropped_props"]
#     imageprops = file["dir_label_s/image_props"]
#     partsdirs = file["dir_label_s/parts_dirs"]
#     partsn_dir = file["dir_label_s/part_sn"]
