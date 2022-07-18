import os.path
import tkfilebrowser
from SaveToHDF5 import load_from_hdf5
from SaveToHDF5 import append_carbon_to_hdf5
import numpy as np
import PySimpleGUI as sg
import matplotlib
import matplotlib.pyplot as plt
from ODStack import mat2_gray
from scipy.ndimage.filters import median_filter
from skimage.measure import label
from skimage.measure import regionprops
from skimage.segmentation import clear_border
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects

def carbon_maps():
    # User interface for selecting file to be processed
    # TJJ 2022

    # tells matplotlib to use TkAgg as the backend GUI integration
    matplotlib.use("TkAgg")
    progdir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    # Allows file browser to open in most recent directory for ease of use
    dirsave = open(progdir + "\\" + "lastdirsave.txt", "r")
    startdirbros = dirsave.read()
    dirsave.close()

    try:
        file_chosen = tkfilebrowser.askopenfilename(title="Pick a .hdf5 file with the processed data",
                                                    okbuttontext="Select", initialdir=startdirbros, filetypes=[("HDF5", "*.hdf5")])
    except:
        file_chosen = ""

    # Save current directory, so it will reopen in the same location when program is run again
    startdirbros = file_chosen.rsplit("\\", 1)[0]
    dirsave = open(progdir + "\\" + "lastdirbros.txt", "w")
    dirsave.write(startdirbros)
    dirsave.close()

    if file_chosen != "":
        data = load_from_hdf5(file_chosen)
        snew = carbon_map_creator(data, file_chosen)

        save_question = sg.Window('Particle Analysis Python Edition V1.0',
                                  [[sg.T('Do you want save the Carbon Map data?')],
                                   [sg.Yes(s=10), sg.No(s=10)]], disable_close=True).read(close=True)

        if save_question[0] == 'Yes':
            append_carbon_to_hdf5(snew, file_chosen)

    return


def carbon_map_creator(snew, data_file, *args):
    # Based on CarbonMaps program by RCM and UOP 7/6/2016
    # Originally for Matlab, converted to Python 3 by TJJ 6/2022

    energy = snew.eVenergy
    stack = snew.spectr
    subdim = int(np.ceil(np.sqrt(len(energy))))

    if len(snew.eVenergy) < 2:
        sg.popup("Too few images for this mapping routine")
        return

    test = energy[(energy < 319) & (energy > 277)]
    if np.size(test) == 0:
        sg.popup("This is not the carbon edge")
        return

    if np.size(args) == 0:
        spthresh = 0.35
        figsav = 0
        nofig = 0
    elif np.size(args) == 1:
        spthresh = args[0]
        figsav = 0
        nofig = 0
    elif np.size(args) == 2:
        spthresh = args[0]
        figsav = 0
        nofig = 1
    elif np.size(args) == 3:
        spthresh = args[0]
        figsav = 1
        nofig = 0
        rootdir = args[1]
        sample = args[2]
    else:
        sg.popup("Invalid argument")
        return

    if spthresh > 1:
        spthresh = spthresh/100

    if nofig == 0 and len(energy) < 6:
        plt.figure(tight_layout = True)
        plt.title(snew.particle + " Single Energy Images")
        for i in range(len(energy)):
            plt.subplot(int(str(subdim)+str(subdim)+str(i+1)))
            plt.imshow(stack[i, :, :], extent=[0, snew.Xvalue, snew.Yvalue, 0])
            plt.axis('image')
            plt.xlabel("X (\u03BCm)")
            plt.ylabel("Y (\u03BCm)")
            plt.set_cmap('gray')
            plt.colorbar()
            plt.title('eV' + str(energy[i]))

    if figsav == 1:
        filename = str(rootdir) + str(sample) + str(snew.particle) + "_SingleEnergyIms"
        plt.savefig(filename, format = 'png')

    sp2idx = np.flatnonzero((284.5 < energy) & (energy < 285.6))
    preidx = np.flatnonzero((277 < energy) & (energy < 283))
    carboxidx = np.flatnonzero((288 < energy) & (energy < 289))
    postidx = np.flatnonzero((310 < energy) & (energy < 325))
    if len(snew.eVenergy) > 10:
        sp2idx = np.round(np.mean(sp2idx))
        preidx = preidx[0]
        carboxidx = np.round(np.mean(carboxidx))
        postidx = postidx[-1]
    else: # numpy method flatnonzero returns single element array if only one value is found, this fixes that.
        if np.size(sp2idx) != 0:
            sp2idx = sp2idx[0]
        if np.size(preidx) != 0:
            preidx = preidx[0]
        if np.size(carboxidx) != 0:
            carboxidx = carboxidx[0]
        if np.size(postidx) != 0:
            postidx = postidx[0]

    # Find Particles
    # meanim = np.sum(snew.spectr, axis=0)
    if np.size(carboxidx) == 0:
        if np.size(postidx) == 0:
            sg.popup("Energy isn't within usable ranges")
            return
        else:
            meanim = snew.spectr[int(postidx), : , :]
    else:
        meanim = snew.spectr[int(carboxidx), :, :]

    meanim[meanim > 0.2] = 0.2
    meanim[meanim < 0] = 0

    grayimage = mat2_gray(meanim)
    thresh = threshold_otsu(grayimage)
    binmap = np.zeros(grayimage.shape)
    binmap[grayimage >= thresh] = 1
    binmap = clear_border(binmap)

    binmap = median_filter(binmap, size=(3, 3))

    # Define Label Matrix
    labelmat = label(binmap, connectivity=2)

    # Filter noise that appears as small particles
    for i in range(np.amax(labelmat)):
        a1, b1 = np.nonzero(labelmat == i)
        linidx1 = np.ravel_multi_index((a1, b1), dims=np.shape(labelmat))
        if len(linidx1) < 7:
            labelmat.flat[linidx1] = 0

    labelmat[labelmat > 0] = 1
    labelmat = label(labelmat, connectivity=2)

    # Assign Particle Serial Numbers and Directories
    numpart = np.amax(labelmat)
    partzero = int(snew.particle[4:] + '0000')
    partsn = (np.array(range(numpart)) + 1) + partzero
    partdir = data_file
    partdirs = np.zeros((numpart, 2), dtype=object)
    partdirs[:, 0] = partdir.rsplit("\\", 1)[0]
    # searchstring = partdir + '\\*' + snew.particle + '*.hdf5'
    # for fname in glob.glob(searchstring):
    partdirs[:, 1] = partdir.split('\\')[-1]

    # Carbon Map
    carb = stack[postidx, :, :] - stack[preidx, :, :]
    if nofig == 0:
        plt.figure(tight_layout=True)
        plt.title(snew.particle + " Maps")
        plt.subplot(221)
        plt.imshow(carb, extent=[0, snew.Xvalue, snew.Yvalue, 0])
        plt.set_cmap('gray')
        plt.colorbar()
        plt.clim(0, np.amax(carb))
        plt.axis('image')
        plt.title('PostEdge-PreEdge')
        plt.xlabel('X (\u03BCm)')
        plt.ylabel('Y (\u03BCm)')

    carb1 = carb.copy()
    carb1[carb1 < 0] = 0
    carb1 = carb1 * binmap
    snew.totc = carb1.copy()
    carb1[carb1 != 0] = 1 # converts to binary (matlab equiv of remove_small_objects does this automatically)
    carb1 = remove_small_objects(carb1.astype(bool), min_size=3)
    carbmask = carb1.copy()
    carbmask[carbmask > 0] = 1
    carbmask = median_filter(carbmask, (3, 3))

    # Inorganic Map
    prepost = stack[int(preidx), :, :] / stack[int(postidx), :, :]
    prev = stack[int(preidx), :, :]
    noise = 3 * np.std(prev[binmap == 0])
    prev1 = prev.copy()
    prev1[prev < noise] = 0
    premask = prev1.copy()
    premask[premask > 0] = 1

    prepost = prepost * premask * binmap
    if nofig == 0:
        plt.subplot(222)
        plt.imshow(prepost, extent=[0, snew.Xvalue, snew.Yvalue, 0])
        plt.set_cmap('gray')
        plt.colorbar()
        plt.axis('image')
        plt.clim(0, 1.0)
        plt.xlabel('X (\u03BCm)')
        plt.ylabel('Y (\u03BCm)')
        plt.title('PreEdge/PostEdge')

    prepost[prepost < 0.5] = 0
    prepostmask = prepost.copy()
    prepostmask[prepost > 1] = 1
    carb1[carb1 != 0] = 1
    prepostmask = remove_small_objects(prepostmask.astype(bool), min_size=5)

    if np.size(sp2idx) != 0:
        doublecarb = stack[int(sp2idx), : , :] - stack[int(preidx), :, :]

        doubcarbnois = np.std(doublecarb[binmap == 0])
        doublecarb1 = doublecarb.copy()
        doublecarb2 = doublecarb / (stack[postidx, : , :] - stack[preidx, :, :]) * (0.4512 / 0.8656)
        doublecarb1[doublecarb1 < (3 * doubcarbnois)] = 0
        doublecarb1 = doublecarb1 * binmap
        spmask = doublecarb1.copy()
        spmask[doublecarb1 > 0] = 1

        sp2 = (stack[int(sp2idx), : , :] - stack[int(preidx), :, :]) / (stack[int(postidx), : , :] - stack[int(preidx), :, :]) * (0.4512 / 0.8656) * spmask
        if nofig == 0:
            plt.subplot(223)
            plt.imshow(sp2, extent=[0, snew.Xvalue, snew.Yvalue, 0])
            plt.set_cmap('gray')
            plt.clim(0, 1.0)
            plt.axis('image')
            plt.colorbar()
            plt.xlabel('X (\u03BCm)')
            plt.ylabel('Y (\u03BCm)')
            plt.title('sp^2 Map')

        sp2nothresh = doublecarb2.copy()
        sp2nothresh[sp2nothresh < 0] = np.NaN
        sp2nothresh[sp2nothresh > 1] = 1
        snew.sp2 = sp2nothresh.copy()
        sp2[sp2 < spthresh] = 0

        finsp2mask = sp2.copy()
        finsp2mask[sp2 > 0] = 1

        # Get rid of small few pixel regions
        finsp2mask = label(finsp2mask, connectivity=2)

        for i in range(1, np.amax(finsp2mask) + 1):
            a1, b1 = np.nonzero(finsp2mask == i)
            linidx1 = np.ravel_multi_index((a1, b1), dims=np.shape(finsp2mask))
            if len(linidx1) < 7:
                finsp2mask.flat[linidx1] = 0
        finsp2mask[finsp2mask > 0] = 1
        bw = np.zeros(finsp2mask.shape)
        bw[finsp2mask >= 0.5] = 1

        imstruct = regionprops(bw.astype(int))
        if len(imstruct) == 0:
            ecc = []
            major = []
            minor = []
            cvex = []
            area = []
        else:
            ecc = np.reshape(imstruct[0].eccentricity, np.shape(imstruct))
            major = np.reshape(imstruct[0].major_axis_length, np.shape(imstruct))
            minor = np.reshape(imstruct[0].minor_axis_length, np.shape(imstruct))
            cvex = np.reshape(imstruct[0].convex_area, np.shape(imstruct))
            area = np.reshape(imstruct[0].area, np.shape(imstruct))

        snew.sooteccentricity = ecc
        snew.sootmajoraxislength = major
        snew.sootminoraxislength = minor
        snew.sootconvexarea = cvex
        snew.sootarea = area

    else:
        sp2 = np.zeros(np.shape(binmap))
        finsp2mask = np.zeros(np.shape(binmap))
        # doublecarb = np.zeros(np.shape(binmap))

    # Combine Maps
    bincompmap = np.array((carbmask, prepostmask, finsp2mask))

    # This first loop creates masks for the individual components over the
    # entire field of view. Each component is then defined as a colored
    # component for visualization.
    colorvec = np.array(((0, 170, 0), (0, 255, 255), (255, 0, 0), (255, 170, 0), (255, 255, 255)))

    matsiz = np.shape(labelmat)
    rgbmat = np.zeros((3, matsiz[0], matsiz[1]))
    redmat = np.zeros(matsiz)
    gremat = np.zeros(matsiz)
    blumat = np.zeros(matsiz)

    l, m = np.nonzero(labelmat > 0)
    labidx = np.ravel_multi_index((l, m), dims=matsiz)
    ccmap = []
    for i in range(len(bincompmap)):
        j, k = np.nonzero(bincompmap[i] > 0)
        if np.size(j) != 0 or np.size(k) != 0:
            linidx = np.ravel_multi_index((j, k), dims=np.shape(bincompmap[i]))
            rejidx = np.setdiff1d(linidx, labidx)
            bincompmap[i].flat[rejidx] = 0
            linidx = np.ravel_multi_index((np.nonzero(bincompmap[i] > 0)), dims=np.shape(bincompmap[i]))
            redmat.flat[linidx] = colorvec[i, 0]
            gremat.flat[linidx] = colorvec[i, 1]
            blumat.flat[linidx] = colorvec[i, 2]
            trmat = np.zeros(np.shape(redmat))
            tgmat = np.zeros(np.shape(gremat))
            tbmat = np.zeros(np.shape(blumat))
            trmat.flat[linidx] = colorvec[i, 0]
            tgmat.flat[linidx] = colorvec[i, 1]
            tbmat.flat[linidx] = colorvec[i, 2]
            ccmap.append(np.dstack((trmat, tgmat, tbmat)))
        else:
            ccmap.append(np.zeros((3, matsiz[0], matsiz[1])))

    # This second loop assigns labels over individual particles defined previously
    labelstr = ['OC', 'In', 'EC']

    compsize = np.zeros((numpart, len(labelstr) + 1))
    partlabel = []
    for i in range(numpart):
        partlabel.append('')
        for j in range(len(labelstr)):
            a1, b1 = np.nonzero(labelmat == i)
            a2, b2 = np.nonzero(bincompmap[j] > 0)
            if (np.size(a1) != 0) and (np.size(b1) != 0) and (np.size(a2) != 0) and (np.size(b2) != 0):
                linidx1 = np.ravel_multi_index((a1, b1), dims = np.shape(labelmat))
                linidx2 = np.ravel_multi_index((a2, b2), dims=np.shape(bincompmap[j]))
                idxcom = np.intersect1d(linidx1, linidx2)
                if len(idxcom) > 3:
                    partlabel[i] = partlabel[i] + labelstr[j]
                    compsize[i , j] = len(idxcom)

        if np.size(partlabel[i]) == 0:
            partlabel[i] = 'NoID'

        compsize[i, j + 1] = len(linidx1)

    if np.size(partlabel) == 0:
        partlabel = 'NoID'
    else:
        partlabel = partlabel

    # Define Outputs
    snew.labelmat = labelmat
    snew.partlabel = partlabel
    snew.partsn = partsn.conj().T
    binmap = np.zeros(np.shape(labelmat))
    binmap[labelmat > 0] = 1
    snew.binmap = binmap
    snew = particle_size(snew)
    xsiz = snew.Xvalue/matsiz[0]
    ysiz = snew.Yvalue/matsiz[1]
    compsize = compsize * (xsiz * ysiz)
    snew.compsize = compsize
    snew.partdirs = partdirs

    rgbmat[0, :, :] = redmat
    rgbmat[1, :, :] = gremat
    rgbmat[2, :, :] = blumat
    snew.rgbcompmap = rgbmat

    temp = bincompmap.copy()
    for i in range(len(bincompmap)):
        temp[i][temp[i] > 1] = 1

    # xdat = np.arange(0, snew.Xvalue, xsiz)
    # ydat = np.arange(0, snew.Yvalue, ysiz)

    # Combined Masks
    if nofig == 0:
        plt.subplot(224)
        plt.imshow(np.uint8(rgbmat.transpose((1, 2, 0))), extent=[0, snew.Xvalue, snew.Yvalue, 0])
        # plt.xticks(xdat)
        # plt.yticks(ydat)
        plt.title('Red = sp2 > ' + str(spthresh) + '\nBlue = pre/post > 0.5\nGreen = Organic')
        plt.axis('image')
        plt.xlabel('X (\u03BCm)')
        plt.ylabel('Y (\u03BCm)')
        plt.show()
        if figsav == 1:
            filename = str(rootdir) + str(sample) + str(snew.particle) + '_Maps.png'
            plt.savefig(filename)

    xysiz = np.shape(carb)
    snew.maps = np.zeros((3, xysiz[0], xysiz[1]))
    snew.maps[0, :, :] = carb
    snew.maps[1, :, :] = prepost
    snew.maps[2, :, :] = sp2

    snew.bincompmap = temp

    return snew

def particle_size(s):
    matsiz = np.shape(s.labelmat)
    xsiz = s.Xvalue/matsiz[1]
    ysiz = s.Yvalue/matsiz[0]
    partsiz = []

    for i in range(np.amax(s.labelmat)):
        j, k = np.nonzero(s.labelmat == i)
        idx = np.ravel_multi_index((j, k), dims=np.shape(s.labelmat))
        partsiz.append(2 * np.sqrt((len(idx) * (xsiz * ysiz)) * np.pi))

    if len(partsiz) > 0:
        s.size = partsiz
    else:
        s.size = []

    return s
