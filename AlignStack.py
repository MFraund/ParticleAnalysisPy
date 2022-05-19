import numpy as np
from DFTRegistration import dft_registration

def align_stack(stack):
    s = stack

    [emax, ymax, xmax] = s.spectr.shape

    xresolution = s.Xvalue/xmax
    yresoltuion = s.Yvalue/ymax

    center = np.ceil(emax/4*3)-1

    spectr = np.zeros((emax, ymax, xmax))

    shifts = np.zeros((emax, 4))

    for k in range(emax):
        shifts[k, :] = dft_registration(np.fft.fft2(s.spectr[int(center)]), np.fft.fft2(s.spectr[k]), 50)
        spectr[k] = ft_matrix_shift(s.spectr[k], -shifts[k, 2], -shifts[k, 3])

    shiftymax = np.ceil(np.max(shifts[:, 2]))
    shiftxmax = np.ceil(np.max(shifts[:, 3]))
    shiftymin = np.ceil(np.abs(np.min(shifts[:, 2])))
    shiftxmin = np.ceil(np.abs(np.min(shifts[:, 3])))

    shiftmatrix = np.zeros((int(emax), int(ymax - shiftymin - shiftymax), int(xmax - shiftxmax - shiftxmin)))

    shiftmatrix[:, :, :] = spectr[:, int(shiftymax):int(ymax - shiftymin), int(shiftxmax):int(xmax - shiftxmin)]

    s.spectr = abs(shiftmatrix)

    s.Xvalue = s.spectr.shape[2] * xresolution
    s.Yvalue = s.spectr.shape[1] * yresoltuion

    return s

def ft_matrix_shift(a, dy, dx):
    ny, nx = np.shape(a)
    rx = np.floor(nx / 2) + 1
    fx = ((np.array(range(nx)) + 1 - rx) / (nx / 2))
    ry = np.floor(ny / 2) + 1
    fy = ((np.array(range(ny)) + 1 - ry) / (ny / 2))

    px = np.atleast_2d(np.fft.ifftshift(np.exp(-1j * dx * np.pi * fx))).conj().T
    py = np.atleast_2d(np.fft.ifftshift(np.exp(-1j * dy * np.pi * fy))).conj().T

    yphase, xphase = np.meshgrid(py, px)

    yphase = yphase.T
    xphase = np.rot90(xphase)

    b = np.abs(np.fft.ifft2(np.fft.fft2(a)) * yphase * xphase)

    return b
