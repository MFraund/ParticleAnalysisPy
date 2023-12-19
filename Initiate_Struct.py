import pickle as pkl

# Pull existing data structure from Process Stack Data and initiate values for Carbon Map data structure
def initiate_carbon_struct(filename):

    try:
        # Pull existing data structure
        with open(filename,'rb') as openedfile:
            prevdata = pkl.load(openedfile)
        structcheck = [i for i in dir(prevdata) if '__' not in i]
        if len(structcheck) > 6:
            data = prevdata
        else:
            spectr = prevdata.spectr
            evenergy = prevdata.eVenergy
            xvalue = prevdata.Xvalue
            yvalue = prevdata.Yvalue
            particle = prevdata.particle

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
    except:
        print(filename + ' is not compatible. Data not processed')
        return

    return data

# Complete Data Structure. Initiated ahead of time for simplicity of Matlab conversion.
# Python can't add new parameters to a class (structure) on the fly like Matlab.
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
