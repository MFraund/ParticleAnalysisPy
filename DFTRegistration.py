# Code translated from Matlab program of same name by Manuel Guizar
# code translation by Taylor Jones

import numpy as np

def dft_registration(buf1ft, buf2ft, usfac = 1):

    if usfac == 0:
        ccmax = np.sum(np.sum(buf1ft * buf2ft.conj()()))
        rfzero = np.sum(np.power(np.abs(buf1ft[:]), 2))
        rgzero = np.sum(np.power(np.abs(buf2ft[:]), 2))
        error = 1.0 - ccmax * ccmax.conj()/(rgzero * rfzero)
        error = np.sqrt(np.abs(error))[0][0]
        diffphase = np.arctan2(np.imag(ccmax), np.real(ccmax))
        output = [error, diffphase]

    elif usfac == 1:
        [m, n] = buf1ft.shape
        cc = np.fft.ifft2(buf1ft * buf2ft.conj())
        [max1, loc1] = [np.amax(cc, 0), np.argmax(cc, 0)]
        loc2 = np.argmax(max1, 0)
        rloc = loc1[loc2]
        cloc = loc2
        ccmax = cc[rloc, cloc]
        rfzero = np.sum(np.power(np.abs(buf1ft[:]), 2)) / (m*n)
        rgzero = np.sum(np.power(np.abs(buf2ft[:]), 2)) / (m*n)
        error = 1.0 - ccmax * ccmax.conj() / (rgzero[1, 1] * rfzero[1, 1])
        error = np.sqrt(np.abs(error))[0][0]
        diffphase = np.arctan2(np.imag(ccmax), np.real(ccmax))
        md2 = np.fix(m / 2)
        nd2 = np.fix(n / 2)
        if rloc > md2:
            row_shift = rloc - m
        else:
            row_shift = rloc

        if cloc > nd2:
            col_shift = cloc - n
        else:
            col_shift = cloc

        output = [error, diffphase, row_shift, col_shift]

    else:
        [m, n] = buf1ft.shape
        mlarge = m * 2
        nlarge = n * 2
        cc = np.zeros((mlarge, nlarge), dtype=complex)
        cc[int(m - np.fix(m / 2)) : int(m + 1 + np.fix((m-1) / 2)), int(n - np.fix(n / 2)) : int(n + 1 + np.fix((n - 1) / 2))] = \
            np.fft.fftshift(buf1ft) * np.fft.fftshift(buf2ft).conj()

        cc = np.fft.ifft2(np.fft.ifftshift(cc))
        [max1, loc1] = [np.amax(cc, 0), np.argmax(cc, 0)]
        loc2 = np.argmax(max1, 0)
        rloc = loc1[loc2]
        cloc = loc2
        ccmax = cc[rloc, cloc]

        [m, n] = cc.shape
        md2 = np.fix(m / 2)
        nd2 = np.fix(n / 2)
        if rloc > md2:
            row_shift = rloc - m
        else:
            row_shift = rloc

        if cloc > nd2:
            col_shift = cloc - n
        else:
            col_shift = cloc

        row_shift = row_shift / 2
        col_shift = col_shift / 2

        if usfac > 2:
            row_shift = np.rint(row_shift * usfac) / usfac
            col_shift = np.rint(col_shift * usfac) / usfac
            dftshift = np.fix(np.ceil(usfac * 1.5) / 2)

            cc = (dftups(buf2ft * buf1ft.conj(), np.ceil(usfac * 1.5), np.ceil(usfac * 1.5), usfac,
                         dftshift - row_shift * usfac, dftshift - col_shift * usfac)).conj() / (md2 * nd2 * np.power(usfac, 2))

            [max1, loc1] = [np.amax(cc, 0), np.argmax(cc, 0)]
            loc2 = np.argmax(max1, 0)
            rloc = loc1[loc2]
            cloc = loc2
            ccmax = cc[rloc, cloc]
            rg00 = dftups(buf1ft * buf1ft.conj(), 1, 1, usfac) / (md2 * nd2 * np.power(usfac, 2))
            rf00 = dftups(buf2ft * buf2ft.conj(), 1, 1, usfac) / (md2 * nd2 * np.power(usfac, 2))
            rloc = rloc - dftshift
            cloc = cloc - dftshift
            row_shift = row_shift + rloc / usfac
            col_shift = col_shift + cloc / usfac

        else:
            rg00 = np.sum(np.sum(buf1ft * buf1ft.conj())) / m / n
            rf00 = np.sum(np.sum(buf2ft * buf2ft.conj())) / m / n

        error = 1.0 - ccmax * ccmax.conj() / (rg00 * rf00)
        error = np.sqrt(np.abs(error))[0][0]
        diffphase = np.arctan2(np.imag(ccmax), np.real(ccmax))

        if md2 == 1:
            row_shift = 0

        if nd2 == 1:
            col_shift = 0

        output = [error, diffphase, row_shift, col_shift]

    return output

def dftups(din, nor = None, noc = None, usfac = 1, roff = 0, coff = 0):
    [nr, nc] = din.shape
    if noc is None:
        noc = nc
    if nor is None:
        nor = nr

    kernc = np.exp((-1j * 2 * np.pi / (nc * usfac)) *
                   (np.atleast_2d(np.fft.ifftshift(np.array(range(nc)))).conj().T - np.floor(nc / 2)) * (np.array(range(int(noc))) - coff))
    kernr = np.exp((-1j * 2 * np.pi / (nr * usfac)) *
                   (np.atleast_2d(np.array(range(int(nor)))).conj().T - roff) * (np.fft.ifftshift(np.array(range(nr))) - np.floor(nr / 2)))
    out = np.matmul(np.matmul(kernr, din), kernc)
    return out
