import numpy as np
import matplotlib.pyplot as plt

def soot_carbox_size_hist(sin, ocopt):
    # ocopt is a flag for excluding soot
    # The follwing fields are appended to Sout:
    # sout.avsootcarb = avsootcarb -- single particle average carbon for soot??
    # sout.avtotc = actotc -- single particle total carbon
    # sout.avcarbox = avcarbox -- single particle carboxylic
    # sout.avsp2 = avsp2 -- average single particle sp2

    # run sootcarboxmap to map the carbox-c=c maps
    sin = soot_carbox_map(sin, 0)

    if len(sin.sootcarbox) == 0:
        sout = sin.copy()
        sout.avsootcarb = np.zeros(np.shape(sin.size))
        sout.avsootcarb[sout.avsootcarb == 0] = np.NaN
        sout.avtotc = np.zeros(np.shape(sin.size))
        sout.avtotc[sout.avtotc == 0] = np.NaN
        sout.avsp2 = np.zeros(np.shape(sin.size))
        sout.avsp2[sout.avsp2 == 0] = np.NaN
        sout.avcarbox = np.zeros(np.shape(sin.size))
        sout.avcarbox[sout.avcarbox == 0] = np.NaN
        return sout

    avsootcarb = np.zeros((1, np.amax(sin.labelmat)))
    avtotc = np.zeros((1, np.amax(sin.labelmat)))
    avsp2 = np.zeros((1, np.amax(sin.labelmat)))
    avcarbox = np.zeros((1, np.amax(sin.labelmat)))

    # Find regions having OC if OC is to be excluded from particle peak averages
    if ocopt == 1:
        ks, ls = np.nonzero(sin.maps[2, :, :] > 0)
        locs = np.ravel_multi_index((ks, ls), dims=np.shape(sin.labelmat))

    # average raw maps (totc, sp2, carbox, sootcarbox) over particle regions
    partspec = np.zeros((len(sin.spectr[:, 0, 0]) , np.max(sin.labelmat)))
    for i in range(np.max(sin.labelmat)):
        k, l = np.nonzero(sin.labelmat == 0)
        loc = np.ravel_multi_index((k, l), dims=np.shape(sin.labelmat))
        if ocopt == 1:
            loc = np.setdiff1d(loc, locs)
        for j in range(len(sin.spectr[:, 0, 0])):
            junkmat = sin.spectr[j, :, :]
            partspec[j , i] = np.mean(np.mean(junkmat[loc]))

        avsootcarb[i] = np.nanmean(np.nanmean(sin.sootcarbox[loc]))
        avtotc[i] = np.nanmean(np.nanmean(sin.totc[loc]))
        avsp2[i] = (partspec[2, i] - partspec[1, i]) / (partspec[-1, i] - partspec[1, i]) * (0.4512/0.8656)
        avcarbox[i] = np.nanmean(np.nanmean(sin.carbox[loc]))

    avsp2[avsp2 < 0] = 0
    avsp2[avsp2 > 1] = 1

    sout = sin.copy()
    sout.avsootcarb = avsootcarb
    sout.avtotc = avtotc
    sout.avcarbox = avcarbox
    sout.avsp2 = avsp2

    return

def soot_carbox_map(snew, plotflg):
    energy = snew.eVenergy.copy()
    stack = snew.spectr.copy()
    subdim = np.ceil(np.sqrt(len(energy)))

    sp2mapidx = np.flatnonzero((energy > 284.9) & (energy < 285.5))
    sp2mapidx = np.round(np.mean(sp2mapidx))
    carboxmapidx = np.flatnonzero((energy >= 288.5) & (energy <= 288.8))
    carboxmapidx = np.round(np.mean(carboxmapidx))
    if np.isnan(carboxmapidx) or np.isnan(sp2mapidx):
        print("Missing one of the energyies for carbox-sp2/TotC map")
        print("returning empty matrices")
        sout = snew.copy()
        sout.sootcarbox = []
        sout.totc = []
        return sout

    # plot selected images
    if plotflg == 1:
        plt.figure()
        plt.subplot(211)
        plt.imshow(stack[sp2mapidx, :, :], extent=[0, snew.Xvalue, 0, snew.Yvalue])
        plt.clim(0, 1.5)
        plt.set_cmap('gray')
        plt.axis('image')
        plt.colorbar()
        plt.title(str(energy[sp2mapidx]) + 'eV')

        plt.subplot(212)
        plt.imshow(stack[carboxmapidx, :, :], extent=[0, snew.Xvalue, 0, snew.Yvalue])
        plt.clim(0, 1.5)
        plt.set_cmap('gray')
        plt.axis('image')
        plt.colorbar()
        plt.title(str(energy[carboxmapidx]) + 'eV')

    # make map of Total Carbon
    sp2carbox = stack[carboxmapidx, :, :] - stack[sp2mapidx, :, :]
    totc = stack[-1, :, :] = stack[0, :, :]
    carbox = stack[carboxmapidx, :, :] - stack[0, :, :]
    binmap = snew.binmap.copy()

    # Mask Maps
    sp2carbox = sp2carbox * binmap
    totc = totc * binmap

    carbox = carbox * binmap
    totc[totc < 0] = 0
    sp2carbox[sp2carbox < 0] = 0
    carbox[carbox < 0] = 0

    if plotflg == 1:
        plt.figure()
        plt.imshow(sp2carbox, extent=[0, snew.Xvalue, 0, snew.Yvalue])
        plt.clim(0, 1.5)
        plt.set_cmap('jet')
        plt.axis('image')
        plt.colorbar()
        plt.title(str(energy[carboxmapidx]) + ' eV - ' + str(energy[sp2mapidx]) + ' eV/total C')

    sout = snew.copy()
    sout.sootcarbox = sp2carbox
    sout.carbox = carbox

    return sout
