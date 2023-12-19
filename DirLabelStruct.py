import PySimpleGUI as sg
import tkfilebrowser as tk
import os.path
# from SaveToHDF5 import append_dir_label_to_hdf5
# from SaveToHDF5 import load_from_hdf5
import numpy as np
from skimage.measure import find_contours
from skimage.measure import regionprops
from CarbonMaps import carbon_map_creator
from SootCarboxSizeHist import soot_carbox_size_hist
from scipy.ndimage import distance_transform_edt
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
import matplotlib
import matplotlib.pyplot as plt
import re
from copy import deepcopy
import pickle as pkl
from Initiate_Struct import initiate_carbon_struct

def dir_label_struct():
    # Adapted to Python from Matlab by TJ Jones 2022
    matplotlib.use("TkAgg")

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
                    if event == "Exit" or event == sg.WIN_CLOSED:
                        break
                    dirsave = open(progdir + "\\" + "lastdirsave.txt", "r")
                    startdirsave = dirsave.read()
                    dirsave.close()
                    try:
                        mappath = tk.askopendirname(title="Chose a place to save output", okbuttontext="Save", initialdir=startdirsave)
                        if mappath == '':
                            break
                    except:
                        mappath = ""
                        break

                    if len(mappath) > 0:
                        print('Processing...')

                        dirlabels = dir_label_maps_struct(rawstackdir, dataset, 0, mappath, sp2, imagetype)

                        savefile = mappath + '\\' + dataset + '.pkl'
                        # append_dir_label_to_hdf5(dirlabels, savefile)
                        with open(savefile, 'w+b') as openedfile:
                            pkl.dump(dirlabels, openedfile)

                        startdirsave = mappath
                        dirsave = open(progdir + "\\" + "lastdirsave.txt", "w")
                        dirsave.write(startdirsave)
                        dirsave.close()
                        filewindow.close()
                        window.close()
                        print('Complete!')
                        break

    filewindow.close()
    window.close()

    return

def dir_label_maps_struct(indir, sample, saveflg, savepath, sp2, impar):
    # Modified from DirLabelMapsStruct.m a matlab program created by Ryan Moffet 7/7/2016
    # Adapted to Python from Matlab by TJ Jones 2022

    sp2 = int(sp2)
    labelcnt = np.zeros((4, len(indir)))
    stackcnt = 0
    mapcnt = 0

    # initialize values
    filenamelst = [[]]
    labelcntlst = [[]]
    partsize = [[]] # particle sizes
    label = [[]] # particle labels( as strings)
    cmpsiz = [[]] # areas of components
    sootcarbox = [[]]
    totalcarbon = [[]] # OD of total carbon(OD(320) - OD(278))
    carbox = [[]]
    asp2 = [[]]
    sootdistcent = [[]] # relative distance of the soot inclusion from the particle center
    sootdistcentinscribed = [[]]
    sootecc = [[]] # eccentricity of the soot inclusion
    sootmaj = [[]] # major axis of soot inclusion
    sootmin = [[]] # minor axis of the soot inclusion
    sootcvex = [[]] # convexity of the soot inclusion
    sootarea = [[]]
    croppedparts = [[]]# Cropped RGB Images of particles
    imageprops = [[]] # image properties: [Xvalue, Yvalue,  # of X pixels,# of Y pixels]
    partdirs = [[]]
    partsn = [[]]

    avgradclass = [[]]
    avgstdclass = [[]]
    npart = [[]]
    singradscans = [[]]

    outrad = [[]]
    radstd = [[]]
    singrad = [[]]

    # begin loop over sample directories
    for j in range(len(indir)):
        directory = os.listdir(indir[j])

        # loop over stacks
        for l in range(len(directory)):
            ind = directory[l].find(".pkl")
            # if the directory has a pickle file...
            if ind != -1:
                file_name = indir[j] + "\\" + directory[l]
                snew = initiate_carbon_struct(file_name)

                if snew is None:
                    break
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
                    sinp = carbon_map_creator(snew, file_name, sp2, savepath, sample[j])
                    sinp = soot_carbox_size_hist(sinp, 1)
                else:
                    sinp = carbon_map_creator(snew, file_name, sp2, 1)
                    sinp = soot_carbox_size_hist(sinp, 1)
                if len(sinp.size) == 0:
                    print('no particles identified, stack skipped')
                    continue

                # count particle classes
                newcount = label_count(sinp)
                labelcnt[:, j] = labelcnt[:, j] + newcount

                # find soot inclusion distance from the center
                sinp = dist_cent(sinp)
                
                # Do radial scans
                radtemp, stdtemp, ntemp, gradtemp = map_radial_scan_spline(sinp, 1, 0)
                avgradclass.append(radtemp)
                avgstdclass.append(stdtemp)
                npart.append(ntemp)
                singradscans.append(gradtemp)
                
                # Crop particle comp maps
                sinp = crop_part(sinp, 0, impar)
                plt.close(plt.gcf())
                for i in range(0, len(sinp.croppedparts)):
                    cpartsiz = np.shape(sinp.croppedparts[0])
                    if np.sum(cpartsiz[0:1]) == 0:
                        print(str(indir[j]) + '\\' + str(directory[1]) + ' IS GIVING EMPTY RGBs!!')
                    
                # Give Labels and size for plotting chemical size distributions
                aclass, siz = chem_siz(sinp)
                # Get image size:
                imsiz = np.shape(sinp.labelmat)

                # Append data from previous files
                filenamelst.append(file_name.split('\\')[-1].split('.')[0])
                labelcntlst.append(labelcnt)
                partsize.append(siz)
                label.append(aclass)
                cmpsiz.append(sinp.compsize) # area of the different components
                sootcarbox.append(sinp.avsootcarb)
                totalcarbon.append(sinp.avtotc)  # height of total carbon/particle
                carbox.append(sinp.avcarbox)  # this is the height of the carbox peak/particle
                asp2.append(sinp.avsp2)  # This is the height of the Sp2 peak/particle
                sootdistcent.append(sinp.sootdistcent)
                sootdistcentinscribed.append(sinp.sootdistcentinscribed)
                sootecc.append(sinp.sooteccentricity)
                sootmaj.append(sinp.sootmajoraxislength)
                sootmin.append(sinp.sootminoraxislength)
                sootcvex.append(sinp.sootconvexarea)
                sootarea.append(sinp.sootarea)
                croppedparts.append(sinp.croppedparts)

                outradtemp, radstdtemp, singradtemp = sum_rad(avgradclass, avgstdclass, npart, singradscans)
                outrad.append(outradtemp)
                radstd.append(radstdtemp)
                singrad.append(singradtemp)

                imageprops.append(sinp.imageprops)
                partdirs.append(sinp.partdirs)
                partsn.append(sinp.partsn)

    # Assign output to data structure
    DirLabels = DirLabelsStruct()
    DirLabels.filenames = filenamelst
    DirLabels.labelcnt = labelcntlst
    DirLabels.partsize = partsize
    DirLabels.label = label
    DirLabels.cmpsiz = cmpsiz
    DirLabels.sootcarbox = sootcarbox
    DirLabels.totalcarbon = totalcarbon
    DirLabels.carbox = carbox
    DirLabels.sp2 = asp2
    DirLabels.outrad = outrad # radial scans: for each particle, a matrix of radial scans of [pre/
    DirLabels.radstd = radstd
    DirLabels.singrad = singrad
    DirLabels.sootdistcent = sootdistcent
    DirLabels.sootdistcentinscribed = sootdistcentinscribed
    DirLabels.sootecc = sootecc
    DirLabels.sootmaj = sootmaj
    DirLabels.sootmin = sootmin
    DirLabels.sootcvex = sootcvex
    DirLabels.sootarea = sootarea
    DirLabels.croppedparts = croppedparts # Extract data by calling like this: DirLabels.CroppedParts{1}{1}
    DirLabels.imageprops = imageprops
    DirLabels.partdirs = partdirs
    DirLabels.partsn = partsn

    return DirLabels

def label_count(sin):
    # Adapted to Python from Matlab by TJ Jones 2022

    occnt = 0
    oceccnt = 0
    ocincnt = 0
    ocecincnt = 0

    for i in range(len(sin.partlabel)):
        ecidx = sin.partlabel[i].find('EC')
        inidx = sin.partlabel[i].find('In')
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
    # Adapted to Python from Matlab by TJ Jones 2022

    nx, ny = np.shape(sin.labelmat)
    totimage = np.nansum(sin.spectr, axis=0)
    totimage[totimage == np.inf] = 0
    sootmap = sin.maps[2, : , :]
    mask = np.zeros((nx, ny))
    rad, cs, cy = max_inscribed_circle(sin, 0)

    xtmpsoot, ytmpsoot = np.nonzero(sin.bincompmap[2] > 0)
    sootlinidx = np.ravel_multi_index((xtmpsoot, ytmpsoot), dims=np.shape(sootmap))

    labelnum = 1
    sootcentlindix = []
    linidx = [[]] * np.max(sin.labelmat)
    cent = [[]] * np.max(sin.labelmat)
    for i in range(0, np.max(sin.labelmat)):
        xtmp, ytmp = np.nonzero(sin.labelmat == labelnum)
        labelnum = labelnum + 1
        linidx[i] = np.ravel_multi_index((xtmp, ytmp), dims=np.shape(sin.labelmat))
        xhat = np.round(np.nansum(xtmp * totimage.flat[linidx[i]]) / np.nansum(totimage.flat[linidx[i]]))
        yhat = np.round(np.nansum(ytmp * totimage.flat[linidx[i]]) / np.nansum(totimage.flat[linidx[i]]))
        cent[i] = np.array((yhat, xhat))
        mask[np.abs(int(xhat)), np.abs(int(yhat))] = 1
        sootinclinidx = np.intersect1d(sootlinidx, linidx[i])

        if len(sootinclinidx) != 0:
            subsootx, subsooty = np.unravel_index(sootinclinidx, np.shape(sootmap))
            sootxhat = np.round(np.nansum(subsootx * sootmap.flat[sootinclinidx]) / np.nansum(sootmap.flat[sootinclinidx]))
            sootyhat = np.round(np.nansum(subsooty * sootmap.flat[sootinclinidx]) / np.nansum(sootmap.flat[sootinclinidx]))
            sootcentlindix.append(np.ravel_multi_index((int(sootxhat), int(sootyhat)), dims=np.shape(sootmap)))

    d1 = distance_transform_edt(np.logical_not(mask))

    mask2 = np.zeros((nx, ny))
    mask2[sin.labelmat > 0] = 1
    if np.array_equal(np.shape(d1), np.shape(mask2)):
        d2 = d1 * mask2 * np.mean([sin.Xvalue / nx, sin.Yvalue / ny])
        sout = deepcopy(sin)
        sout.disttocent = d2.copy()
        d3 = d2.copy()
        d4 = d1.copy()

        # loop over all particles in LabelMat
        for i in range(0, len(linidx)):
            d3.flat[linidx[i]] = d3.flat[linidx[i]] / np.amax(d3.flat[linidx[i]])
            d4.flat[linidx[i]] = d4.flat[linidx[i]] / rad[i]

        # mark out relative distance of soot inclusion from the center.
        if len(sootcentlindix) == 0:
            sout.sootdistcent = []
            sout.sootdistcentinscribed = []
        else:
            sout.sootdistcent = d3.flat[sootcentlindix]
            sout.sootdistcentinscribed = d4.flat[sootcentlindix]

    else:
        sout = deepcopy(sin)
        sout.disttocent = np.zeros(np.shape(sin.spectr[0, : , :]))
        sout.sootdistcent = []
        sout.sootdistcentinscribed = []

    return sout

def map_radial_scan_spline(sin, comp, fig):
    # plots radial scans of components comp in STXM stack data structure
    # Sin. Must run CarbonMaps.m first.

    # COMPONENTS:
    # comp=1 -- Total Carbon
    # comp=2 -- Pre/Post ratio
    # comp=3 -- sp2
    # Adapted to Python from Matlab by TJ Jones 2022
    
    # Compute single particle radial scans.
    if comp == 3:
        sp2idx = np.nonzero((sin.eVenergy > 284.5) & (sin.eVenergy < 285.5))
        preidx = np.nonzero((sin.eVenergy > 277) & (sin.eVenergy < 283))
        if len(sin.eVenergy) > 10:
            sp2idx = np.round(np.mean(sp2idx))
            preidx = preidx[0]

    bins = np.arange(0, 1, 0.1)
    distspec = np.zeros((len(bins), len(sin.size)))
    labelnum = 1
    for i in range(0, np.amax(sin.labelmat)):
        xind, yind = np.nonzero(sin.labelmat == labelnum)
        labelnum = labelnum + 1
        linidx = np.ravel_multi_index((xind, yind), dims=np.shape(sin.labelmat))
        if comp == 3:
            cmpmap = sin.spectr[sp2idx, :, :] - sin.spectr[preidx, :, :]
            cmpmapnois = np.std(cmpmap[sin.binmap == 0])
            cmpmap[cmpmap < 3 * cmpmapnois] = 0
        else:
            cmpmap = sin.maps[comp, :, :] - np.min(sin.maps[comp, :, :])

        rawpeak = cmpmap.flat[linidx] / np.amax(cmpmap.flat[linidx])
        rawdist = sin.disttocent.flat[linidx] / np.amax(sin.disttocent.flat[linidx])
        avg = []
        for j in range(0, len(bins) - 1):
            tidx = np.nonzero((rawdist >= bins[j]) & (rawdist < bins[j + 1]))
            if (len(tidx) != 0) or (tidx != np.nan):
                avg.append([np.nanmean(rawdist[tidx]), np.nanmean(rawpeak[tidx])])

        avg = np.array(avg)
        if len(avg) != 0:
            if np.amax(np.shape(avg)) > 2:
                f = interp1d(avg[:, 0], avg[:, 1], kind='linear', fill_value='extrapolate')
                distspec[:, i] = f(bins)
                avg = []
            else:
                avg = []
                continue
        else:
            continue

        if fig == 1:
            plt.figure()
            plt.plot(rawdist, rawpeak, 'k.', ms=10)
            plt.plot(bins, distspec[:, i], 'r-', lw=3)
            plt.title('Particle # ' + str(i) + ', ' + str(sin.partlabel[i]))
            plt.legend('Raw Data', 'Boxcar Average', 'Linear Interpolation')
            plt.xlabel('Rel. Dist. from Center')
            plt.ylabel('COOH (Norm. Abs.')

    # Now average single particle radial scans by particle type
    aclass, siz = chem_siz(sin)
    classstr = ['OC', 'ECOCIn', 'ECOC', 'InOC', 'NoID']
    avgradclass = np.zeros((len(bins), 5))
    stdradclass = np.zeros((len(bins), 5))
    npart = np.zeros(5)
    singradscans = [None] * 5 # closest python has to a blank cell array

    aclassnum = 1
    for i in range(0, int(np.amax(aclass))):
        cidx = np.array(np.nonzero(aclass == aclassnum)).flatten()
        aclassnum = aclassnum + 1
        npart[i] = len(cidx)
        avgradclass[:, i] = np.nanmean(distspec[:, cidx], axis=1).flatten()
        singradscans[i] = distspec[:, cidx]
        stdradclass[:, i] = np.std(np.abs(distspec[:, cidx]), ddof=1, axis=1).flatten()

        if fig == 1:
            plt.figure()
            plt.errorbar(bins, avgradclass[:, i], stdradclass[:, i], 'k.-')
            plt.title('Average Radial Scan for ' + str(classstr[i]))
            plt.xlabel('Rel. Dist. from Center')
            plt.ylabel('COOH (Norm. Abs.')

    return avgradclass, stdradclass, npart, singradscans

def chem_siz(sin):
    siz = sin.size
    aclass = np.zeros(len(sin.partlabel))

    for i in range(0, len(sin.partlabel)):
        noididx = 'NoID' == sin.partlabel[i]
        ocidx = 'OC' == sin.partlabel[i]
        ocaidx = 'OCsp2' == sin.partlabel[i]
        ecidx = [_.start() for _ in re.finditer('EC', sin.partlabel[i])]
        inidx = [_.start() for _ in re.finditer('In', sin.partlabel[i])]
        kidx = [_.start() for _ in re.finditer('K', sin.partlabel[i])]

        if (ocidx is True) or (ocaidx is True): # OC
            aclass[i] = 1
        elif (len(ecidx) != 0) and ((len(inidx) != 0) or (len(kidx) != 0)): # ECOCIn
            aclass[i] = 2
        elif (len(ecidx) != 0) and ((len(inidx) == 0) or (len(kidx) == 0)): # ECOC
            aclass[i] = 3
        elif (len(ecidx) == 0) and ((len(inidx) != 0) or (len(kidx) != 0)): # InOC
            aclass[i] = 4
        elif noididx: # NoID
            aclass[i] = 5

    return aclass, siz

def max_inscribed_circle(sin, display=1):
    #
    #
    #     by Tolga Birdal
    #
    #
    #     Maximum Inscribed Circle
    #
    #     Or in other words, "largest inner circle" , "maximum empty circle" etc.
    #
    #     This is a very common problem in computational geometry, and it is not
    #     simple to solve efficiently.
    #     Addressing 2D image/contour processing, I couldn't find a good
    #     implementation on the web. Generally, the reasonable way of solving
    #     this problem is to make use of Voronoi Diagrams, which are generally
    #     O(nlogn).
    #
    #     After analyzing the problem a bit, I noticed that the problem can
    #     easily and approximately be solved using well-known distance transform.
    #
    #     Here is how:
    #
    #     The computational aim can be written as:
    #     (x, y) maximizes r = min_{i} r_{i}
    #     where r_i = ||(x_i, y_i) ? (x, y)|| and d_i = r_i ? r
    #
    #     (x_i, y_i): Pairs data points
    #     (x, y), r : Pair, scalar circle centre and radius
    #
    #     In non-mathematical terms:
    #
    #     1.  The center of the maximum inscribed circle will lie inside the
    #         polygon
    #     2.  The center of such a circle will be furthest from any point on the
    #         edges of the polygon.
    #
    #     So we seek for the point that lies inside the polygon and has maximal
    #     distance to the closest edge. This is exactly the maximum pixel of the
    #     distance transform lying inside the contour.
    #
    #     Notice that this approach is completely in-precise. It is only
    #     pixel-precise and never subpixel accurate. However, unlike
    #     optimization approaches, it does guarantee a global convergence. In
    #     the case of ambiguity, any of the solutions will be valid.
    #
    #     To detect the points inside the region, inpolygon remains very slow.
    #     So, I make use of the great code of Darren Engwirda, here. As well as
    #     being contained in this package, it can also be downloaded from:
    #     http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test
    #
    #     Here are other implemnatations, which are more accurate, but much slower
    #     than my approach (only slower in Matalb of course!)
    #
    #     Using
    #     http://www.mathworks.com/matlabcentral/fileexchange/2794-cvoronoi
    #
    #     Using "Largest pixel that doesn't cross any point" approach:
    #     http://www.mathworks.com/matlabcentral/newsreader/view_thread/283296
    #
    #     Here is a sample call:
    #
    #       I=imread('hand_contour.png');
    #       [R cx cy]=max_inscribed_circle(I)
    #
    #   Cheers,
    # Adapted to Python from Matlab by TJ Jones 2022

    b = find_contours(sin.binmap, 0)
    sz = np.shape(sin.binmap)
    #rad = np.zeros((np.max(np.shape(b)), 1))
    rad = np.zeros((len(b), 1))

    #for i in range(0, np.max(np.shape(b))):
    for i in range(0, len(b)):
        # get the contour
        contourimage = np.zeros(sz)
        linb = np.ravel_multi_index((b[i][:, 0].astype(int), b[i][:, 1].astype(int)), dims=sz)
        contourimage.flat[linb] = 1
        # python has no equiv to bwtraceboundary from matlab. But simply doing b[i] seems to get the same result...
        contour = b[i]
        x = contour[:, 1]
        y = contour[:, 0]

        # find the maximum inscribed circle:
        # The point that has the maximum distance inside the given contour is the
        # center. The distance of to the closest edge (tangent edge) is the radius.
        bw = distance_transform_edt(np.logical_not(contourimage)) # python's distance transform uses inverted logic compared to matlab's
        mx, my = np.meshgrid(range(0, sz[1]), range(0, sz[0]))
        vin, von = in_poly(np.array([mx.flatten(), my.flatten()]).T, np.array([x, y]).T)
        ind = np.ravel_multi_index((my.flat[vin], mx.flat[vin]), dims=sz)
        r = np.amax(bw.flat[ind])
        rind = np.argmax(bw.flat[ind])
        # if multiple solutions, just take first.
        if np.shape(r) != () and np.shape(rind) != (): # shape returns () if it isn't an array
            rad[i] = r[0]
            rind = rind[0]
        else:
            rad[i] = r
            rind = rind
        cy, cx = np.unravel_index(ind[rind], sz)

        # display result
        if display:
            # bw[cy, cx] = r + 20 # to emphasize the center
            # plt.figure()
            # plt.imshow(bw)
            plt.plot(x, y, color='red', linewidth=2)
            theta = np.append(np.linspace(0, 2*np.pi), 0)
            plt.plot(np.cos(theta) * rad[i] + cx, np.sin(theta) * rad[i] + cy, color='green', linewidth=2)
            plt.show()

    return rad, cx, cy

def in_poly(p, node, edge = [], tol = 1.0e-12):
    #  INPOLY: Point-in-polygon testing.
    #
    # Determine whether a series of points lie within the bounds of a polygon
    # in the 2D plane. General non-convex, multiply-connected polygonal
    # regions can be handled.
    #
    # SHORT SYNTAX:
    #
    #   in = inpoly(p,node);
    #
    #   p   : The points to be tested as an Nx2 array [x1 y1; x2 y2; etc].
    #   node: The vertices of the polygon as an Mx2 array [X1 Y1; X2 Y2; etc].
    #         The standard syntax assumes that the vertices are specified in
    #         consecutive order.
    #
    #   in  : An Nx1 logical array with IN(i) = TRUE if P(i,:) lies within the
    #         region.
    #
    # LONG SYNTAX:
    #
    #  [in,on] = inpoly(p,node,edge);
    #
    #  edge: An Mx2 array of polygon edges, specified as connections between
    #        the vertices in NODE: [n1 n2; n3 n4; etc]. The vertices in NODE
    #        do not need to be specified in connsecutive order when using the
    #        extended syntax.
    #
    #  on  : An Nx1 logical array with ON(i) = TRUE if P(i,:) lies on a
    #        polygon edge. (A tolerance is used to deal with numerical
    #        precision, so that points within a distance of
    #        eps^0.8*norm(node(:),inf) from a polygon edge are considered "on"
    #        the edge.
    #
    # EXAMPLE:
    #
    #   polydemo;       # Will run a few examples
    #
    # See also INPOLYGON

    # The algorithm is based on the crossing number test, which counts the
    # number of times a line that extends from each point past the right-most
    # region of the polygon intersects with a polygon edge. Points with odd
    # counts are inside. A simple implementation of this method requires each
    # wall intersection be checked for each point, resulting in an O(N*M)
    # operation count.
    #
    # This implementation does better in 2 ways:
    #
    #   1. The test points are sorted by y-value and a binary search is used to
    #      find the first point in the list that has a chance of intersecting
    #      with a given wall. The sorted list is also used to determine when we
    #      have reached the last point in the list that has a chance of
    #      intersection. This means that in general only a small portion of
    #      points are checked for each wall, rather than the whole set.
    #
    #   2. The intersection test is simplified by first checking against the
    #      bounding box for a given wall segment. Checking against the bbox is
    #      an inexpensive alternative to the full intersection test and allows
    #      us to take a number of shortcuts, minimising the number of times the
    #      full test needs to be done.
    #
    #   Based on MATLAB Code by Darren Engwirda: 2005-2007
    # Adapted to Python from Matlab by TJ Jones 2022

    # Error Checking
    #####################################################################################

    # p = np.array(p)
    # node = np.array(node).T
    nnode = np.shape(node)[0]
    if len(edge) == 0:
        edge1 = np.array(range(1, nnode))[:, None]
        edge2 = np.array(range(2, nnode+1))[:, None]
        edge3 = np.array((nnode, 1))[None, :]
        edge = np.append(edge1, edge2, 1)
        edge = np.append(edge, edge3, 0)
        edge = edge - 1

    if np.shape(p)[1] != 2:
        print('P must be an Nx2 array.')
        cn = []
        on = []
        return cn, on
    if np.shape(node)[1] != 2:
        print('Node must be an Mx2 array.')
        cn = []
        on = []
        return cn, on
    if np.shape(edge)[1] != 2:
        print('EDGE must be an Mx2 array.')
        cn = []
        on = []
        return cn, on
    if (np.amax(edge) > nnode) or np.any(edge < 0):
        print('Invalid EDGE.')
        cn = []
        on = []
        return cn, on

    # Pre-Processing
    #########################################################################
    n = np.shape(p)[0]
    nc = np.shape(edge)[0]

    # Choose the direction with the biggest range as the "y-coordinate" for the
    # test. This should ensure that the sorting is done along the best
    # direction for long and skinny problems wrt either the x or y axes.
    dxy = np.amax(p, axis=0) - np.amin(p, axis=0)
    if dxy[0] > dxy[1]:
        # Flip co-ords if x rang eis bigger
        p = p[:, [1, 0]]
        node = node[:, [1, 0]]
    tol = tol * np.amin(dxy)

    # Sort test points by y-value
    y = np.sort(p[:, 1])
    i = np.searchsorted(y, p[:, 1])
    x = p[i, 0]

    # Main Loop
    ############################################################################
    # Because we're dealing with mod(cn,2) we don't have
    # to actually increment the crossing number, we can
    # just flip a logical at each intersection (faster!)
    cn = np.zeros((n, 1))
    on = cn.copy()

    for k in range(0, nc):

        # Nodes in current edge
        n1 = edge[k, 0]
        n2 = edge[k, 1]

        # Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
        #           - also get xmin = min(x1,x2), xmax = max(x1,x2)
        y1 = node[n1, 1]
        y2 = node[n2, 1]
        if y1 < y2:
            x1 = node[n1, 0]
            x2 = node[n2, 0]
        else:
            yt = y1
            y1 = y2
            y2 = yt
            x1 = node[n2, 0]
            x2 = node[n1, 0]

        if x1 > x2:
            xmin = x2
            xmax = x1
        else:
            xmin = x1
            xmax = x2

        # Binary search to find first point with y<=y1 cor current edge
        if y[0] >= y1:
            start = 0
        elif y[n - 1] < y1:
            start = n - 1
        else:
            lower = 0
            upper = n - 1
            for j in range(0, n):
                start = int(np.round(0.5 * (lower + upper)))
                if y[start] < y1:
                    lower = start
                elif y[start - 1] < y1:
                    break
                else:
                    upper = start

        # Loop through points
        for j in range(start, n):
            # Check the bounding-box for the edge before doing the intersection
            # test. Take shortcuts wherever possible!

            yy = y[j]
            if yy <= y2:
                xx = x[j]
                if xx <= xmin:
                    if xx <= xmax:

                        # Check if we're "on" the edge
                        on[j] = on[j] or (np.abs((y2-yy)*(x1-xx)-(y1-yy)*(x2-xx)) < tol)

                        # Do the actual intersection test
                        if (yy < y2) and ((y2-y1)*(xx-x1) < (yy-y1)*(x2-x1)):
                            cn[j] = not cn[j]

                # Deal with point exactly at vertices
                elif yy < y2:
                    # Has to cross edge
                    cn[j] = not cn[j]

            else:
                # Due to sorting, no point with > y
                # value need to be checked
                break

    # Re-index to undo the sorting
    cn[i] = np.logical_or(cn, on)
    on[i] = on

    return cn.astype(bool).flatten(), on.astype(bool).flatten()

def crop_part(snew, plotfig, atype):
    # need to run CarbonMaps.m first!
    # Input-------------------------------------------------------------------
    # Snew: the STXM data structure obtained from SingStackProc.m or ProcDir.m AND CarbonMaps.m
    
    # Output------------------------------------------------------------------
    # Snew.ImageProps - matrix containing [XLength,YLength,NXPix,NYPix] and
    #   rows corresponding to the number of particles
    # Snew.CroppedParts - cell array of cropped rgb images for later plotting
    #   by CropRGBPlot.m
    # 2016 RCMoffet
    
    # Crop Particles ---------------------------------------------------------
    # Portions of cropping code below was posted on the matlab forum by  Anton
    # Semechko: https://www.mathworks.com/matlabcentral/answers/28158#answer_36447
    # Adapted to Python from Matlab by TJ Jones 2022

    im = snew.binmap.copy()
    siz = np.shape(im)

    if 'SVD' == atype:
        ecsum = np.sum(snew.rgbcompmap[0, :, :])
        ocsum = np.sum(snew.rgbcompmap[1, :, :])
        insum = np.sum(snew.rgbcompmap[2, :, :])
        compvec = np.array([ecsum, ocsum, insum])
        compvec = compvec > 1
        if np.sum(compvec) > 1:
            timage = np.zeros((3, np.shape(snew.binmap)[0], np.shape(snew.binmap)[1]))
            aslice = np.mean(snew.spectr, 2)
            aslice[snew.binmap == 0] = 0
            aslice = aslice / np.amax(aslice) * 255
            timage[compvec > 0, :, :] = aslice
        elif np.sum(compvec) == 2:
            if (compvec[1] > 0) and (compvec[2] > 0) and (compvec[0] == 0):
                orgspec = component_spec_svd(snew, 'COOH', 0)
                inorgspec = component_spec_svd(snew, 'Inorg', 0)
                timage = stack_svd(snew, ['COOH', 'Inorg'], orgspec, inorgspec)
            elif (compvec[0] > 0) and (compvec[1] > 0) and (compvec[2] == 0):
                orgspec = component_spec_svd(snew, 'COOH', 0)
                scspec = component_spec_svd(snew, 'sp2', 0)
                timage = stack_svd(snew, ['sp2', 'COOH'], scspec, orgspec)
        elif np.sum(compvec) == 3:
            inorgspec = component_spec_svd(snew, 'Inorg', 0)
            orgspec = component_spec_svd(snew, 'COOH', 0)
            ecspec = component_spec_svd(snew, 'sp2', 0)
            timage = stack_svd(snew, ['sp2', 'COOH', 'Inorg'], ecspec, orgspec, inorgspec)
    else:
        timage = snew.rgbcompmap.copy()

    # Label the disconnected foreground regions (using 8 conned neighbourhood)
    l = snew.labelmat.copy()
    # get the bounding box around each object
    bb = np.array([region.bbox for region in regionprops(l)])
    # crop the individual objects and store them in a cell
    n = np.amax(l[:])
    rgbim = [[]] * n
    for i in range(0, n):
        # get the bb of the i-th object and offset by 2 pixels in all directions
        bb_i = np.ceil(bb[i])
        idx_x = np.array([bb_i[0] - 2, bb_i[0] + bb_i[2] + 2])
        idx_y = np.array([bb_i[1] - 2, bb_i[1] + bb_i[3] + 2])
        if idx_x[0] < 1:
            idx_x[0] = 1
        if idx_y[0] < 1:
            idx_y[0] = 1
        if idx_x[1] > siz[1]:
            idx_x[1] = siz[1]
        if idx_y[1] > siz[0]:
            idx_y[1] = siz[0]

        # crop the object and write the objcell
        im = (l == i)
        rgbim[i] = timage[:, int(idx_y[0]):int(idx_y[1]), int(idx_x[0]):int(idx_x[1])]

    # Prepare output
    improps = np.array([snew.Xvalue, snew.Yvalue, siz[1], siz[0]])
    xres = snew.Xvalue / siz[1]
    yres = snew.Yvalue / siz[0]
    improps = np.tile(improps, (n, 1))
    for i in range(0, len(rgbim)):
        subimsiz = np.shape(rgbim[i])
        improps[i, :] = np.array([subimsiz[1] * xres, subimsiz[0] * yres, subimsiz[1], subimsiz[0]])

    snew.imageprops = improps.copy()
    snew.croppedparts = rgbim.copy()

    # plot interactive figure if needed
    # if plotfig == 1:
        # crop_part_rgb_plot(snew.croppedparts, snew.imageprops, snew.partdirs, snew.partsn)

    return snew

def component_spec_svd(sin, comp, plot):
    # uses maps to isolate regions of interest to extract a spectrum.
    # Sin.DiffMaps have a m by n by j structure with m and n being the 
    # spatial dimensions and j being the compoenent maps
    #
    # j=1 -- COOH
    # j=2 -- inorganic
    # j=3 -- sp2
    # j=4 -- potassium
    # j=5 -- CO3
    #
    # specout is the output spectrum of the form [energies,OD]
    # RC Moffet, 2010
    # Adapted to Python from Matlab by TJ Jones 2022

    energy = sin.eVenergy.copy()
    spectr = sin.spectr.copy()

    # setup flags for isolating ROIs
    if 'COOH' == comp:
        cmp = 1
        exclude = 3
    elif 'Inorg' == comp:
        cmp = 2
        exclude = []
    elif 'K' == comp:
        cmp = 2
        exclude = []
    elif 'sp2' == comp:
        cmp = 3
        excldue = []
    elif 'co3' == comp:
        cmp = 5
        exclude = []

    cnt = 1

    # Intersect component maps and exclude overlapping pixels
    linidx = [[]] * len(sin.bincompmap) # alternate way to make empty 'cell' array
    for i in range(len(sin.bincompmap)):
        j, k = np.nonzero(sin.bincompmap[i] > 0)
        linidx[i] = np.ravel_multi_index((j, k), dims=np.shape(sin.bincompmap[i]))

    overlapidx = []
    cmpidx = []
    crap = np.zeros(np.shape(sin.bincompmap[0])) # ??

    if len(excldue) == 0:
        if len(cmp) == 1:
            cmpidx = linidx[cmp]
        else:
            for i in range(len(cmp)):
                cmpidx = np.union1d(cmpidx, linidx[i])
    else:
        for i in exclude:
            if (i == cmp) or (len(linidx[i]) == 0) or (i == 6):
                continue
            else:
                tidx = np.intersect1d(linidx[cmp], linidx[i])
                overlapidx = [overlapidx, tidx]
        cmpidx = np.setdiff1d(linidx[cmp], overlapidx)

    # Plot spectrum and binary map
    partspec = np.zeros(np.shape(energy))
    for j in range(len(energy)):
        junkmat = spectr[j, :, :]
        partspec[j] = np.mean(junkmat[cmpidx])

    amap = np.zeros(np.shape(sin.bincompmap[0]))
    amap[cmpidx] = 1
    if plot == 1:
        h = plt.figure()
        plt.subplot(121)
        plt.plot(energy, partspec, lw=3)
        plt.title(str(plt.gca) + ' Spectrum ' + str(comp), fontsize=28)
        plt.xlabel('Energy (eV)', fontsize=24)
        plt.ylabel('OD', fontsize=24)
        plt.setp(plt.gca, fontsize=18)
        plt.xlim((np.amin(energy), np.amax(energy)))
        if (len(partspec) == 0) or np.isnan(partspec):
            ymax = 1
        else:
            ymax = np.abs(np.amax(partspec))
        plt.ylim((0, ymax))

        plt.subplot(122)
        plt.imshow(amap)
        plt.colormaps('bone')
        plt.axis('image')
        plt.title(str(plt.gca) + ' Regions ' + str(comp))
        plt.setp(h, position=[1, 1, 1013, 639])

    # Return output
    # dime = np.shape(energy) # This was a fix for a bug in Matlab, it isn't an issue in Python
    # if dime[1] > 1:
    #     energy = np.transpose(energy)
    #     partspect = np.transpose(partspec)

    specout = [energy, partspec]
    npix = len(cmpidx)

    return specout, npix

def stack_svd(s, cmpstr, *args):
    # function svdM = stackSVD(S,varargin)
    #
    # Generation of component maps from STXM stack data using Singular Value Decomposition (SVD) analysis with reference spectra
    # R.C. Moffet, T.R. Henn February 2009
    #
    # Inputs
    # ------
    # S          Standard stack structure array containing the STXM data converted to optical density
    # varargin   arbitrary number of two-column vectors representing the chemical component's reference spectra.
    #            The first column of the reference spectra vectors contains the energy points (in eV), the second column
    #            holds the corresponding absorbance in optical density.
    #
    # Output
    # ------
    # svdM       m x n x k array containing the component maps. (m, n correspond to the number of stack pixels in y and x dimension,
    #            k is determined by the number of component reference spectra passed to the script as varargin.
    #            svdM(:,:,i) contains the component map corresponding to the ith reference spectrum in varargin
    #
    # Adapted to Python from Matlab by TJ Jones 2022

    # args contains the reference data for the # nofcps of components cpnts
    cps = args
    nofcps = len(args)

    # determine number of components
    sevmin = np.amin(s.eVenergy)
    sevmax = np.amax(s.eVenergy)

    # determine energy range of reference spectra
    cpsevmin = np.amin(cps[0][:, 0])
    cpsevmax = np.amax(cps[0][:, 0])

    for k in range(1, nofcps):
        if np.amin(cps[k][:, 0]) > cpsevmin:
            cpsevmin = np.amin(cps[k][:, 0])
        if np.amax(cps[k][: ,0]):
            cpsevmax = np.amax(cps[k][:, 0])

    # truncate experimental data if reference spectra energy range is smaller than experimental data
    if cpsevmin > sevmin:
        fstidx = np.nonzero(s.eVenergy < cpsevmin)[0]
        lstidx = np.nonzero(s.eVenergy < cpsevmin)[-1]

        # remove experimental data that is out of reference spectra data range
        s.eVenergy[fstidx:lstidx] = []
        s.spectr[fstidx:lstidx] = []

    if cpsevmax < sevmax:
        fstidx = np.nonzero(s.eVenergy < cpsevmax)[0]
        lstidx = np.nonzero(s.eVenergy < cpsevmax)[-1]

        # remove experimental data that is out of reference spectra data range
        s.eVenergy[fstidx:lstidx] = []
        s.spectr[fstidx:lstidx] = []

    # construction of the reference components coefficient matrix
    m = np.zeros((len(s.eVenergy), nofcps))

    for k in range(0, nofcps):
        ppoly = CubicSpline(cps[k][:, 0], cps[k][:, 1])
        m[:, k] = ppoly(s.eVenergy)

    # construction of the coefficient matrix pseudo inverse using SVD
    pseudoinvm = np.linalg.pinv(m)
    svdm = np.zeros((nofcps, np.shape(s.spectr, 0), np.shape(s.spectr, 1)))
    for y in range(0, np.shape(s.spectr, 0)):
        for x in range(0, np.shape(s.spectr, 1)):
            temp = s.spectr[:, y, x]
            svdm[:, y, x] = pseudoinvm * temp

    svdplot = plot_svd_rgb_new(svdm, s, cmpstr)

    return svdplot

def plot_svd_rgb_new(svdin, sin, cmpstr):
    # svdin is the svd matrix
    # Sin is the stack data
    # cmpstr is an array of strings that contain the different components
    # Adapted to Python from Matlab by TJ Jones 2022

    svdin[svdin < 0] = 0
    numcmp = len(cmpstr)
    idx = 0
    svdidx = 0
    dims = np.shape(svdin[0, :, :])
    colormat = np.zeros((3, dims[0], dims[1]))

    # red
    sootidx = np.nonzero(('sp2' == cmpstr) == 1)
    if len(sootidx) == 0:
        colormat[0, :, :] = np.zeros(np.shape(svdin[0, :, :]))
        idx = idx + 1
    else:
        maxcontrib = np.amax(svdin[svdidx, :, :])
        tmpmat = svdin[svdidx, :, :]
        tmpmat[tmpmat > (0.95 * maxcontrib)] = maxcontrib
        colormat[idx, :, :] = (255 / maxcontrib) * tmpmat
        idx = idx + 1
        svdidx = svdidx + 1

    # green
    orgidx = np.nonzero(('COOH' == cmpstr) == 1)
    if len(orgidx) == 0:
        colormat[1, :, :] = np.zeros(np.shape(svdin[0, :, :]))
        idx = idx + 1
    else:
        maxcontrib = np.amax(svdin[svdidx, :, :])
        tmpmat = svdin[svdidx, :, :]
        tmpmat[tmpmat > (0.95 * maxcontrib)] = maxcontrib
        colormat[idx, :, :] = (255 / maxcontrib) * tmpmat
        idx = idx + 1
        svdidx = svdidx + 1

    # blue
    inorgidx = np.nonzero(('Inorg' == cmpstr) == 1)
    if len(orgidx) == 0:
        colormat[2, :, :] = np.zeros(np.shape(svdin[0, :, :]))
        idx = idx + 1
    else:
        maxcontrib = np.amax(svdin[svdidx, :, :])
        tmpmat = svdin[svdidx, :, :] ^ 2
        tmpmat[tmpmat > (0.95 * maxcontrib)] = maxcontrib
        colormat[idx, :, :] = (255 / maxcontrib) * tmpmat
        idx = idx + 1
        svdidx = svdidx + 1

    cnt = 1

    ymax = len(svdin[0, :, 0])
    xmax = len(svdin[:, 0, 0])
    xticks = np.array(range(0, xmax, xmax / 6))
    yticks = np.array(range(0, ymax, ymax / 6))

    plt.figure()
    idx = np.nonzero(sin.labelmat == 0)
    c1 = colormat[0, :, :]
    c2 = colormat[1, :, :]
    c3 = colormat[2, :, :]

    c1[idx] = 0
    c2[idx] = 0
    c3[idx] = 0
    colormat1 = np.zeros((3, np.shape(c1)[0], np.shape(c1)[1]))
    if (numcmp == 3) or (numcmp == 4):
        colormat1[0, :, :] = c1
        colormat1[1, :, :] = c2
        colormat1[2, :, :] = c3
    elif numcmp == 2:
        colormat1[0, :, :] = c1
        colormat1[1, :, :] = c2
        colormat1[2, :, :] = c3

    plt.imshow(colormat1, extent=[0, sin.xvalue, 0, sin.Yvalue])
    plt.setp(plt.gca, 'TickDir', 'Out', 'FontSize', 14, 'FontName', 'Ariel')
    plt.axis('image')
    plt.show()

    return colormat1

def crop_part_rgb_plot(dirlabels, sample):
    # This function calls another function that uses STACKLab which will be converted to python later,
    # and it isn't called by the previous function. I still wrote it for future use.
    # Adapted to Python from Matlab by TJ Jones 2022

    rgbin = dirlabels.croppedparts[sample]
    improps = dirlabels.imageprops[sample]
    partsdirs = dirlabels.partdirs[sample]
    partsn = dirlabels.partsn[sample]
    fileloc = dirlabels.saveloc
    filename = dirlabels.filename

    # Plots cropped rgb images from carbons maps. When you click on a particle, carbon maps is run and STACKLab is opened.
    n = len(rgbin)
    msize = np.ceil(np.sqrt(n))
    figsiz = 1
    window = plt.get_current_fig_manager().window
    scrsz = [1, 1, window.winfo_screenwidth(), window.winfo_screenheight()] # simulates matlab's get(0, 'ScreenSize')
    plt.figure('units', 'normalized', 'position', [0, 0, figsiz * scrsz[3] / scrsz[3], figsiz * scrsz[2] / scrsz[3]])
    plt.setp(plt.gcf(), 'color', 'k')

    x = np.zeros((msize, msize))
    y = np.zeros((msize, msize))
    xsiz = np.zeros(n)
    ysiz = np.zeros(n)
    for i in range(0, n):
        xsiz[i] = improps[n, 0] / improps[n, 2]
        ysiz[i] = improps[n, 1] / improps[n, 3]
        imsiz = np.shape(rgbin[i][0, :, :])
        x[i] = imsiz[1] * xsiz[i]
        y[i] = imsiz[0] * ysiz[i]
    x = x.T
    y = y.T

    # Calculate maximum dimensions of the entire figure, initialize variables
    ymax = np.sum(np.amax(y, axis=1))
    xmax = np.amax(np.sum(x, axis=1))
    xymax = np.amax([ymax, xmax])
    x = x / xymax
    y = y / xymax
    normxpixsize = xsiz / xymax
    normypixsize = ysiz / xymax
    xmax = np.amax(x)
    ymax = np.amax(y, axis=1)

    rcnt = 0
    ccnt = 0
    ybot = 0
    xpos = 0
    ypos = ybot - ymax[rcnt]
    xmove = 0
    ctr = 0
    emptycnt = 0
    xt = x.T
    yt = y.T
    for i in range(0, n):
        if len(rgbin[i]) == 0:
            ctr = ctr + 1
            emptycnt = emptycnt + 1
            continue
        ydat = np.array(range(0, yt[ctr], normypixsize[ctr]))
        xdat = np.array(range(0, xt[ctr], normypixsize[ctr]))
        if ccnt <= msize:
            ccnt = ccnt + 1
            xpos = xpos + xmove
        else:
            ccnt = 1
            xpos = 0
            rcnt = rcnt + 1
            ypos = ypos - ymax[rcnt]
            ccnt = ccnt + 1
        plt.axis('position', [xpos, ypos, xt[ctr], yt[ctr]])
        plt.imshow(rgbin[i], extent=[0, xdat, 0, ydat])
        disp_data_file(partsdirs[i, :], partsn[i], [fileloc, filename]) # incomplete function
        plt.setp(plt.gca, 'XTick', [], 'YTick', [], 'visible', 'off')
        xmove = xt[ctr]
        ctr = ctr + 1

    return

def disp_data_file(src, evt, fileinfo):
    # This function requires another program that will be developed later (STACKLab)
    # It will not be written until that function has been made
    # Adapted to Python from Matlab by TJ Jones 2022

    print('Error: Incomplete program!')
    print('How did you get here?')

def sum_rad(radscan, std, nums, singscans):

    radout = np.zeros(np.shape(radscan[0]))
    stdout = np.zeros(np.shape(std[0]))
    num = np.zeros(len(nums[0]))
    # sing = [[]] * (len(num) - 1)
    sing = [np.array([])] * 5

    for i in range(0, len(radscan)):
        for j in range(0, len(num) - 1):
            if not np.all(np.shape(sing[j])):
                sing[j] = singscans[i][j]
            else:
                sing[j] = np.append(sing[j], singscans[i][j], axis=1)
            sing[j][sing[j] <= 0] = np.nan

    for i in range(0, len(num)):
        if (len(sing[i]) == 0) or (i == len(radout[0, :])):
            stdout[:, i] = np.zeros(len(radout[:, i]))
        else:
            radout[:, i] = np.nanmean(sing[i], axis=1)
            stdout[:, i] = np.nanstd(sing[i], axis=1)

    return radout, stdout, sing

# Data structure to store final data. Is kind of redundant with snew, might clean up later.
class DirLabelsStruct:
    def __init__(self):
        self.labelcnt = []
        self.partsize = []
        self.label = []
        self.cmpsiz = []
        self.sootcarbox = []
        self.totalcarbon = []
        self.carbox = []
        self.sp2 = []
        self.outrad = []
        self.radstd = []
        self.singrad = []
        self.sootdistcent = []
        self.sootdistcentinscribed = []
        self.sootecc = []
        self.sootmaj = []
        self.sootmin = []
        self.sootcvex = []
        self.sootarea = []
        self.croppedparts = []
        self.imageprops = []
        self.partdirs = []
        self.partsn = []
        
    