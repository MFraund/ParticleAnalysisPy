import numpy as np
from scipy.ndimage.filters import median_filter
from skimage import exposure
from skimage.filters import threshold_otsu
import matplotlib.pyplot as plt


def od_stack(structin, method, plotflag):
    stack = structin.spectr
    evlength = len(structin.eVenergy)

    s = structin
    s.spectr = np.zeros((structin.spectr.shape[0], structin.spectr.shape[1], structin.spectr.shape[2]))

    if method == 'C':
        imagebuffer = np.mean(stack, 0)
        imagebuffer = median_filter(imagebuffer, (3, 3))
        grayimage = mat2_gray(imagebuffer)
        mask = np.zeros(grayimage.shape)
        mask[grayimage >= 0.85] = 1

    elif method == 'O':
        imagebuffer = np.mean(stack, 0)
        grayimage = mat2_gray(imagebuffer)
        grayimage = exposure.adjust_gamma(grayimage, 15)

        thresh = threshold_otsu(grayimage)
        mask = np.zeros(grayimage.shape)
        mask[grayimage >= thresh] = 1

    elif method == 'map':
        imagebuffer = np.mean(stack, 0)
        grayimage = mat2_gray(imagebuffer)
        grayimage = exposure.adjust_gamma(grayimage, 2)

        thresh = threshold_otsu(grayimage)
        mask = np.zeros(grayimage.shape)
        mask[grayimage >= thresh] = 1

    else:
        print("Error! No thresholding method defined! Input structure not converted!")
        return

    izero = np.atleast_2d(np.zeros((evlength, 2)))
    izero[:, 0] = s.eVenergy[:, 0]

    for cnt in range(evlength):
        buffer = stack[cnt]
        izero[cnt, 1] = np.mean(np.mean(buffer[mask == 1]))

    s.izero = izero

    for k in range(evlength):
        s.spectr[k] = -np.log(stack[k]/izero[k, 1])

    if plotflag == 1:
        xlim = [0, s.Xvalue]
        ylim = [s.Yvalue, 0]

        plt.figure(tight_layout = True)
        plt.subplot(221)
        plt.imshow(imagebuffer, extent=[xlim[0], xlim[1], ylim[0], ylim[1]])
        plt.axis('image')
        plt.colorbar()
        plt.set_cmap('gray')
        plt.title("Raw Intensity Stack Mean")
        plt.xlabel("X-Position (µm)")
        plt.ylabel("Y-Position (µm)")

        plt.subplot(222)
        plt.imshow(np.mean(s.spectr, 0), extent=[xlim[0], xlim[1], ylim[0], ylim[1]])
        plt.axis('image')
        plt.colorbar()
        plt.set_cmap('gray')
        plt.title("Optical Denity Stack Mean")
        plt.xlabel("X-Position (µm)")
        plt.ylabel("Y-Position (µm)")

        plt.subplot(223)
        plt.imshow(mask, extent=[xlim[0], xlim[1], ylim[0], ylim[1]])
        plt.axis('image')
        plt.colorbar()
        plt.title("Izero Region Mask")
        plt.xlabel("X-Position (µm)")
        plt.ylabel("Y-Position (µm)")

        plt.subplot(224)
        plt.plot(izero[:, 0], izero[:, 1])
        plt.title("Izero")
        plt.xlabel("Photon energy (eV)")
        plt.ylabel("Raw Counts")

        plt.show(block=False)

    return s


def mat2_gray(inmat):
    limits = [np.min(inmat), np.max(inmat)]
    delta = limits[1] - limits[0]
    outmat = (inmat - limits[0]) / delta
    return outmat
