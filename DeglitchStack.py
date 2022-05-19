import numpy as np

def deglitch_stack(s, badenergy):
    numbad = len(badenergy)
    resolution = 0.1

    for counter in range(numbad):
        rmindex = np.delete(np.argwhere((s.eVenergy < (badenergy[counter] + resolution)) &
                                        (s.eVenergy > (badenergy[counter] - resolution))), 1, 1)

        if len(rmindex) != 0:
            s.eVenergy = np.delete(s.eVenergy, rmindex, 0)
            s.spectr = np.delete(s.spectr, rmindex, 0)

        else:
            print('Error: no energy match for input = ', badenergy[counter], ' eV! Image not removed!!')

    return s
