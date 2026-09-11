"""
Reproduces: "Projected versus fully atomistic MATBG/orderparams real space.ipynb"

Plots real-space order parameters (Inter/Intra sublattice, Inter/Intra valley)
from Hartree-Fock output files, reconstructed on the atomistic real-space
lattice using the coordinate files produced for a given twist angle index
`ntheta`.

Cells are delimited with "" so this file can be run cell-by-cell in
editors/IDEs that support Jupyter-style cells (VS Code, Spyder, etc.), just
like the original notebook.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import re
import os

# geometry routines

def rot(angle, g):
    return np.array([np.cos(angle) * g[0] - np.sin(angle) * g[1], np.cos(angle) * g[1] + np.sin(angle) * g[0]])


def _wigner_seitz_cell(ntheta, ndim, nlayers, RotateLayers, t1, t2, cs, sn, a1, a2):
    """Atom positions before lattice relaxation, FORGE Geometry.f90:WignerSeitzCell.

    Returns an array of shape (2, ndim): row 0 is x, row 1 is y, in the same
    atom order FORGE uses (and therefore the same order as every HF output
    file indexed by site i=1..ndim).
    """
    nrad = 3 * ntheta
    rMax = 3.0 * ntheta ** 2 + 3.0 * ntheta + 1.0
    sq3 = np.sqrt(3.0)

    n = np.arange(-nrad, nrad + 1)
    n1, n2 = np.meshgrid(n, n, indexing='ij')
    n1 = n1.ravel()
    n2 = n2.ravel()

    def within_cell(m1, m2):
        rTemp1 = np.abs(m1 * (3.0 * ntheta + 1) + m2 * (3.0 * ntheta + 2)) + 1e-6
        rTemp2 = np.abs(m1 - m2 * (3.0 * ntheta + 1)) + 1e-6
        rTemp3 = np.abs(m1 * (3.0 * ntheta + 2) + m2) + 1e-6
        return (rTemp1 < rMax) & (rTemp2 < rMax) & (rTemp3 < rMax)

    def lattice_points(m1, m2):
        mask = within_cell(m1, m2)
        return np.stack([m1[mask] * a1[0] + m2[mask] * a2[0],
                          m1[mask] * a1[1] + m2[mask] * a2[1]])

    blocks = []
    for nlayer in range(1, nlayers + 1):
        if RotateLayers[nlayer - 1] == -1:
            an1 = n1 + 1.0 / 3.0
            an2 = n2 - 2.0 / 3.0
            bn1 = n1 + 2.0 / 3.0
            bn2 = n2 - 1.0 / 3.0
            A_boundary = np.array([[(-t2[0] - t1[0]) / 3.0], [(-t2[1] - t1[1]) / 3.0]])
            B_boundary = np.array([[(t2[0] + t1[0]) / 3.0], [(t2[1] + t1[1]) / 3.0]])
        else:
            an1 = (n1 + 1.0 / 3.0) * (cs - sn / sq3) - 2 * (n2 - 2.0 / 3.0) * sn / sq3
            an2 = (n2 - 2.0 / 3.0) * (cs + sn / sq3) + 2 * (n1 + 1.0 / 3.0) * sn / sq3
            bn1 = (n1 + 2.0 / 3.0) * (cs - sn / sq3) - 2 * (n2 - 1.0 / 3.0) * sn / sq3
            bn2 = (n2 - 1.0 / 3.0) * (cs + sn / sq3) + 2 * (n1 + 2.0 / 3.0) * sn / sq3
            A_boundary = np.array([[(t1[0] + t2[0]) / 3.0], [(t1[1] + t2[1]) / 3.0]])
            B_boundary = np.array([[(-t1[0] - t2[0]) / 3.0], [(-t1[1] - t2[1]) / 3.0]])

        blocks.append(lattice_points(an1, an2))
        blocks.append(A_boundary)
        blocks.append(lattice_points(bn1, bn2))
        blocks.append(B_boundary)

    coords = np.concatenate(blocks, axis=1)
    assert coords.shape[1] == ndim
    return coords


def set_angle(ntheta):

    # graphene lattice vectors and layer-stacking parameters, FORGE Setup.f90
    a1 = np.array([0.5, np.sqrt(3.0) / 2])
    a2 = np.array([-0.5, np.sqrt(3.0) / 2])
    nlayers = 2
    RotateLayers = [-1, 1]

    ndim = 4 * (3 * ntheta ** 2 + 3 * ntheta + 1)

    # reciprocal and moire lattice vectors, twist-angle rotation, FORGE Main.f90
    aMoire = 3.0 * ntheta ** 2 + 3.0 * ntheta + 1.0
    g1 = (4.0 * np.pi / 3.0) / aMoire * ((3 * ntheta + 1) * a1 + a2)
    g2 = (4.0 * np.pi / 3.0) / aMoire * ((3 * ntheta + 2) * a2 - a1)  # = g1+g2 in FORGE's g12

    cs = 1.0 - 1.0 / (2.0 * aMoire)
    sn = np.sqrt(1.0 - cs ** 2)
    ang = -0.5 * np.arccos(cs)  # rotation that symmetrizes the cell around the x-axis

    g1 = rot(ang, g1)
    g2 = rot(ang, g2)

    t1 = ntheta * a1 + (ntheta + 1) * a2
    t2 = -(ntheta + 1) * a1 + (2 * ntheta + 1) * a2

    # Coords(:,1), Coords(:,2) before lattice relaxation, FORGE Geometry.f90/Main.f90
    coords = _wigner_seitz_cell(ntheta, ndim, nlayers, RotateLayers, t1, t2, cs, sn, a1, a2)
    coords = rot(ang, coords)
    t1 = rot(ang, t1)
    t2 = rot(ang, t2)

    x, y = coords[0], coords[1]

    vk1 = rot(-np.pi / 6, t1 / np.sqrt(3))
    vk2 = rot(-np.pi / 6 - np.pi / 3, t1 / np.sqrt(3))
    vk3 = rot(-np.pi / 6 - 2 * np.pi / 3, t1 / np.sqrt(3))
    vk4 = rot(-np.pi / 6 - 3 * np.pi / 3, t1 / np.sqrt(3))
    vk5 = rot(-np.pi / 6 - 4 * np.pi / 3, t1 / np.sqrt(3))
    vk6 = rot(-np.pi / 6 - 5 * np.pi / 3, t1 / np.sqrt(3))

    return ndim, t1, t2, g1, g2, x, y, vk1, vk2, vk3, vk4, vk5, vk6

def show_order_params(pairs, nrows=2):
    """Print the order parameters and draw them under the bottom panels.

    A label is one index per Pauli matrix, in the order sigma (sublattice),
    tau (valley), mu (layer), so '0z0' renders as sigma_0 tau_z mu_z. Values
    are laid out on `nrows` centred rows in the blank strip left at the
    bottom of the figure by subplots_adjust.
    """
    print('#################### order parameters ###################')
    for label, value in pairs:
        print(label, value)

    fig = plt.gcf()
    ncols = int(np.ceil(len(pairs) / nrows))
    for k, (label, value) in enumerate(pairs):
        s, t, m = label
        row, col = divmod(k, ncols)
        fig.text((col + .5) / ncols, .085 - .045 * row,
                 rf'$\langle \sigma_{s} \tau_{t} \mu_{m} \rangle = {value:.4f}$',
                 ha='center', va='center', fontsize=14)


# read filename

filename = sys.argv[1]

srch = re.search(r"(Inter|Intra)Sub(Inter|Intra)Val", filename)
if srch is None:
    raise ValueError(f"unexpected filename: {filename}")
sub, val = srch.groups()   # e.g. ("Inter", "Intra")

srch = re.search(r"-i(\d+)-", filename)
if srch is None:
    raise ValueError(f"no I-index in: {filename}")

# geometry
ntheta = int(srch.group(1))
ndim, t1, t2, g1, g2, x, y, vk1, vk2, vk3, vk4, vk5, vk6 = set_angle(ntheta)
ang = -1 / 2 * np.arccos(1 - 1 / (6 * ntheta ** 2 + 6 * ntheta + 2))


size = 36000 / (3 * ntheta ** 2 + 3 * ntheta + 1)


# Order Parameters
# Inter-sublattice inter-valley
if sub == "Inter" and val == "Inter":
        
    op = np.genfromtxt(filename, dtype=np.dtype(float))

    fkakpb1 = (op[:ndim // 4, 0] + 1j * op[:ndim // 4, 1]) * \
        np.exp(-1 * 1j * 8 / 3 * np.pi * (np.cos(ang) * x[:ndim // 4] + np.sin(ang) * y[:ndim // 4])) * \
        np.exp(0 * 1j * 8 / 3 * np.pi * (np.sin(ang) * (y)[:ndim // 4]))

    fkakpb2 = (op[ndim // 4:, 0] + 1j * op[ndim // 4:, 1]) * \
        np.exp(-1 * 1j * 8 / 3 * np.pi * (np.cos(ang) * x[ndim // 2:3 * ndim // 4] - np.sin(ang) * y[ndim // 2:3 * ndim // 4])) * \
        np.exp(-0 * 1j * 8 / 3 * np.pi * (np.sin(ang) * (y)[:ndim // 4]))

    fkbkpa1 = (op[:ndim // 4, 2] + 1j * op[:ndim // 4, 3]) * \
        np.exp(-1 * 1j * 8 / 3 * np.pi * (np.cos(ang) * x[:ndim // 4] + np.sin(ang) * y[:ndim // 4])) * \
        np.exp(0 * 1j * 8 / 3 * np.pi * (np.sin(ang) * (y)[:ndim // 4]))

    fkbkpa2 = (op[ndim // 4:, 2] + 1j * op[ndim // 4:, 3]) * \
        np.exp(-1 * 1j * 8 / 3 * np.pi * (np.cos(ang) * x[ndim // 2:3 * ndim // 4] - np.sin(ang) * y[ndim // 2:3 * ndim // 4])) * \
        np.exp(-0 * 1j * 8 / 3 * np.pi * (np.sin(ang) * (y)[:ndim // 4]))


    # bottom layer -> left half (columns 1-2)
    plt.figure(figsize=(20, 11))
    plt.subplots_adjust(left=.03, right=.97, top=.90, bottom=.12, wspace=.35, hspace=.15)
    plt.gcf().text(.5, .98, os.path.basename(filename), ha='center', va='top', fontsize=16)

    plt.subplot(2, 4, 1)
    plt.title("$Re(\\rho_{K'Bb, KAb})$", fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=np.real(fkakpb1), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 2)
    plt.title("$Im(\\rho_{K'Bb, KAb})$", fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=np.imag(fkakpb1), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 5)
    plt.title("$Re(\\rho_{K'Ab, KBb})$", fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=np.real(fkbkpa1), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 6)
    plt.title("$Im(\\rho_{K'Ab, KBb})$", fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=np.imag(fkbkpa1), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    # plt.show()

    # top layer -> right half (columns 3-4)

    plt.subplot(2, 4, 3)
    plt.title("$Re(\\rho_{K'Bt, KAt})$", fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=np.real(fkakpb2), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 4)
    plt.title("$Im(\\rho_{K'Bt, KAt})$", fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=np.imag(fkakpb2), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 7)
    plt.title("$Re(\\rho_{K'At, KBt})$", fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=np.real(fkbkpa2), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 8)
    plt.title("$Im(\\rho_{K'At, KBt})$", fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=np.imag(fkbkpa2), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    # plt.ylim([-t1[1] / 1.7, t1[1] / 1.7])
    # plt.xlim([-t1[1] / 1.7, t1[1] / 1.7])
    plt.axis('equal')
    plt.axis('off')

    # plt.show()

    show_order_params([
        ('xx0', 1 * np.real(np.sum(fkakpb1)) + 1 * np.real(np.sum(fkakpb2)) + 1 * np.real(np.sum(fkbkpa1)) + 1 * np.real(np.sum(fkbkpa2))),
        ('yy0', -1 * np.real(np.sum(fkakpb1)) - 1 * np.real(np.sum(fkakpb2)) + 1 * np.real(np.sum(fkbkpa1)) + 1 * np.real(np.sum(fkbkpa2))),
        ('xy0', 1 * np.imag(np.sum(fkakpb1)) + 1 * np.imag(np.sum(fkakpb2)) + 1 * np.imag(np.sum(fkbkpa1)) + 1 * np.imag(np.sum(fkbkpa2))),
        ('yx0', 1 * np.imag(np.sum(fkakpb1)) + 1 * np.imag(np.sum(fkakpb2)) - 1 * np.imag(np.sum(fkbkpa1)) - 1 * np.imag(np.sum(fkbkpa2))),
        ('xxz', -1 * np.real(np.sum(fkakpb1)) + 1 * np.real(np.sum(fkakpb2)) - 1 * np.real(np.sum(fkbkpa1)) + 1 * np.real(np.sum(fkbkpa2))),
        ('yyz', +1 * np.real(np.sum(fkakpb1)) - 1 * np.real(np.sum(fkakpb2)) - 1 * np.real(np.sum(fkbkpa1)) + 1 * np.real(np.sum(fkbkpa2))),
        ('xyz', -1 * np.imag(np.sum(fkakpb1)) + 1 * np.imag(np.sum(fkakpb2)) - 1 * np.imag(np.sum(fkbkpa1)) + 1 * np.imag(np.sum(fkbkpa2))),
        ('yxz', -1 * np.imag(np.sum(fkakpb1)) + 1 * np.imag(np.sum(fkakpb2)) + 1 * np.imag(np.sum(fkbkpa1)) - 1 * np.imag(np.sum(fkbkpa2))),
    ])
# Intra-sublattice inter-valley
if sub == "Intra" and val == "Inter":

    op = np.genfromtxt(filename, dtype=np.dtype(float))
    op = np.array(op, dtype=float)

    opn = np.zeros(ndim, dtype=complex)
    opn[:ndim // 2] = (op[:ndim // 2, 0] + 1j * op[:ndim // 2, 1]) * \
        np.exp(1 * 1j * 8 / 3 * np.pi * (np.cos(ang) * x[:ndim // 2] + np.sin(ang) * y[:ndim // 2])) * \
        np.exp(-0 * 1j * 8 / 3 * np.pi * (np.sin(ang) * (y)[:ndim // 2]))
    opn[ndim // 2:] = (op[ndim // 2:, 0] + 1j * op[ndim // 2:, 1]) * \
        np.exp(1 * 1j * 8 / 3 * np.pi * (np.cos(-ang) * x[ndim // 2:] + np.sin(-ang) * y[ndim // 2:])) * \
        np.exp(-0 * 1j * 8 / 3 * np.pi * (np.sin(-ang) * (y)[ndim // 2:]))
    op = opn

    # bottom layer -> left half (columns 1-2)
    plt.figure(figsize=(20, 11))
    plt.subplots_adjust(left=.03, right=.97, top=.90, bottom=.12, wspace=.35, hspace=.15)
    plt.gcf().text(.5, .98, os.path.basename(filename), ha='center', va='top', fontsize=16)

    plt.subplot(2, 4, 1)
    plt.title("$Re(\\rho_{KAb,K'Ab})$", fontsize=14)
    plt.scatter(x[1 * ndim // 4:2 * ndim // 4], y[1 * ndim // 4:2 * ndim // 4], c=np.real(op[1 * ndim // 4:2 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 2)
    plt.title("$Im(\\rho_{KAb,K'Ab})$", fontsize=14)
    plt.scatter(x[1 * ndim // 4:2 * ndim // 4], y[1 * ndim // 4:2 * ndim // 4], c=np.imag(op[1 * ndim // 4:2 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 5)
    plt.title("$Re(\\rho_{KBb,K'Bb})$", fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=np.real(op[0 * ndim // 4:1 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 6)
    plt.title("$Im(\\rho_{KBb,K'Bb})$", fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=np.imag(op[0 * ndim // 4:1 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')


    # top layer -> right half (columns 3-4)

    plt.subplot(2, 4, 3)
    plt.title("$Re(\\rho_{KAt,K'At})$", fontsize=14)
    plt.scatter(x[3 * ndim // 4:4 * ndim // 4], y[3 * ndim // 4:4 * ndim // 4], c=np.real(op[3 * ndim // 4:4 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 4)
    plt.title("$Im(\\rho_{KAt,K'At})$", fontsize=14)
    plt.scatter(x[3 * ndim // 4:4 * ndim // 4], y[3 * ndim // 4:4 * ndim // 4], c=np.imag(op[3 * ndim // 4:4 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 7)
    plt.title("$Re(\\rho_{KBt,K'Bt})$", fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=np.real(op[2 * ndim // 4:3 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 8)
    plt.title("$Im(\\rho_{KBt,K'Bt})$", fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=np.imag(op[2 * ndim // 4:3 * ndim // 4]), marker='h', s=size, cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    show_order_params([
        ('0x0', np.real(np.sum(op[0 * ndim // 4:4 * ndim // 4]))),
        ('0y0', -np.imag(np.sum(op[0 * ndim // 4:4 * ndim // 4]))),
        ('zx0', np.real(np.sum(op[1 * ndim // 4:2 * ndim // 4]) + np.sum(op[3 * ndim // 4:4 * ndim // 4]) - np.sum(op[0 * ndim // 4:1 * ndim // 4]) - np.sum(op[2 * ndim // 4:3 * ndim // 4]))),
        ('zy0', -np.imag(np.sum(op[1 * ndim // 4:2 * ndim // 4]) + np.sum(op[3 * ndim // 4:4 * ndim // 4]) - np.sum(op[0 * ndim // 4:1 * ndim // 4]) - np.sum(op[2 * ndim // 4:3 * ndim // 4]))),
        ('0xz', np.real(-np.sum(op[0 * ndim // 4:2 * ndim // 4]) + np.sum(op[2 * ndim // 4:4 * ndim // 4]))),
        ('0yz', -np.imag(-np.sum(op[0 * ndim // 4:2 * ndim // 4]) + np.sum(op[2 * ndim // 4:4 * ndim // 4]))),
        ('zxz', np.real(-np.sum(op[1 * ndim // 4:2 * ndim // 4]) + np.sum(op[3 * ndim // 4:4 * ndim // 4]) + np.sum(op[0 * ndim // 4:1 * ndim // 4]) - np.sum(op[2 * ndim // 4:3 * ndim // 4]))),
        ('zyz', -np.imag(-np.sum(op[1 * ndim // 4:2 * ndim // 4]) + np.sum(op[3 * ndim // 4:4 * ndim // 4]) + np.sum(op[0 * ndim // 4:1 * ndim // 4]) - np.sum(op[2 * ndim // 4:3 * ndim // 4]))),
    ])
# Intra-sublattice intra-valley
if sub == "Intra" and val == "Intra":

    opr = np.genfromtxt(filename, dtype=np.dtype(float))
    opr = np.array(opr, dtype=float)


    # bottom layer -> left half (columns 1-2)
    plt.figure(figsize=(20, 11))
    plt.subplots_adjust(left=.03, right=.97, top=.90, bottom=.12, wspace=.35, hspace=.15)
    plt.gcf().text(.5, .98, os.path.basename(filename), ha='center', va='top', fontsize=16)
    
    op = opr[:, 1]

    plt.subplot(2, 4, 1)
    plt.title("$density \ \\sum_\\eta \\rho_{\\eta A b, \\eta A b}$", fontsize=14)
    plt.scatter(x[1 * ndim // 4:2 * ndim // 4], y[1 * ndim // 4:2 * ndim // 4], c=op[1 * ndim // 4:2 * ndim // 4], s=size, marker='h', cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 2)
    plt.title("$density \ \\sum_\\eta \\rho_{\\eta B b, \\eta B b}$",  fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=op[0 * ndim // 4:1 * ndim // 4], s=size, marker='h', cmap='coolwarm')  # ,vmin=-.0001,vmax=.0001)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')
    
    op = opr[:, 0]

    plt.subplot(2, 4, 5)
    plt.title("$valley \ polarization \ \\sum_\\eta \\eta \\rho_{\\eta A b, \\eta A b}$",  fontsize=14)
    plt.scatter(x[1 * ndim // 4:2 * ndim // 4], y[1 * ndim // 4:2 * ndim // 4], c=op[1 * ndim // 4:2 * ndim // 4], s=size, marker='h', cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 6)
    plt.title("$valley \ polarization \ \\sum_\\eta \\eta \\rho_{\\eta B b, \\eta B b}$", fontsize=14)
    plt.scatter(x[0 * ndim // 4:1 * ndim // 4], y[0 * ndim // 4:1 * ndim // 4], c=op[0 * ndim // 4:1 * ndim // 4], s=size, marker='h', cmap='coolwarm')  # ,vmin=-.0001,vmax=.0001)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')


    # top layer -> right half (columns 3-4)

    op = opr[:, 1]

    plt.subplot(2, 4, 3)
    plt.title("$density \ \\sum_\\eta \\rho_{\\eta A t, \\eta A t}$",  fontsize=14)
    plt.scatter(x[3 * ndim // 4:4 * ndim // 4], y[3 * ndim // 4:4 * ndim // 4], c=op[3 * ndim // 4:4 * ndim // 4], s=size, marker='h', cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 4)
    plt.title("$density \ \\sum_\\eta \\rho_{\\eta B t, \\eta B t}$",  fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=op[2 * ndim // 4:3 * ndim // 4], s=size, marker='h', cmap='coolwarm')  # ,vmin=-.0001,vmax=.0001)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')
    
    op = opr[:, 0]

    plt.subplot(2, 4, 7)
    plt.title("$valley \ polarization \ \\sum_\\eta \\eta \\rho_{\\eta A t, \\eta A t}$",fontsize=14)
    plt.scatter(x[3 * ndim // 4:4 * ndim // 4], y[3 * ndim // 4:4 * ndim // 4], c=op[3 * ndim // 4:4 * ndim // 4], s=size, marker='h', cmap='coolwarm')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 8)
    plt.title("$valley \ polarization \ \\sum_\\eta \\eta \\rho_{\\eta B t, \\eta B t}$",fontsize=14)
    plt.scatter(x[2 * ndim // 4:3 * ndim // 4], y[2 * ndim // 4:3 * ndim // 4], c=op[2 * ndim // 4:3 * ndim // 4], s=size, marker='h', cmap='coolwarm')  # ,vmin=-.0001,vmax=.0001)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')


    show_order_params([
        ('0z0', .5 * np.sum(opr[0 * ndim // 4:4 * ndim // 4, 0])),
        ('0zz', .5 * (-np.sum(opr[0 * ndim // 4:2 * ndim // 4, 0]) + np.sum(opr[2 * ndim // 4:4 * ndim // 4, 0]))),
        ('zz0', .5 * (-np.sum(opr[0 * ndim // 4:1 * ndim // 4, 0]) + np.sum(opr[1 * ndim // 4:2 * ndim // 4, 0]) - np.sum(opr[2 * ndim // 4:3 * ndim // 4, 0]) + np.sum(opr[3 * ndim // 4:4 * ndim // 4, 0]))),
        ('zzz', .5 * (+np.sum(opr[0 * ndim // 4:1 * ndim // 4, 0]) - np.sum(opr[1 * ndim // 4:2 * ndim // 4, 0]) - np.sum(opr[2 * ndim // 4:3 * ndim // 4, 0]) + np.sum(opr[3 * ndim // 4:4 * ndim // 4, 0]))),
        ('000', .5 * np.sum(opr[0 * ndim // 4:4 * ndim // 4, 1])),
        ('00z', .5 * (-np.sum(opr[0 * ndim // 4:2 * ndim // 4, 1]) + np.sum(opr[2 * ndim // 4:4 * ndim // 4, 1]))),
        ('z00', .5 * (-np.sum(opr[0 * ndim // 4:1 * ndim // 4, 1]) + np.sum(opr[1 * ndim // 4:2 * ndim // 4, 1]) - np.sum(opr[2 * ndim // 4:3 * ndim // 4, 1]) + np.sum(opr[3 * ndim // 4:4 * ndim // 4, 1]))),
        ('z0z', .5 * (+np.sum(opr[0 * ndim // 4:1 * ndim // 4, 1]) - np.sum(opr[1 * ndim // 4:2 * ndim // 4, 1]) - np.sum(opr[2 * ndim // 4:3 * ndim // 4, 1]) + np.sum(opr[3 * ndim // 4:4 * ndim // 4, 1]))),
    ])
# Inter-sublattice intra-valley
if sub == "Inter" and val == "Intra":

    op = np.genfromtxt(filename, dtype=np.dtype(float))
    op = np.array(op, dtype=float)

    fkakb1 = op[:ndim // 4, 0] + 1j * op[:ndim // 4, 1]
    fkakb2 = op[ndim // 4:ndim // 2, 0] + 1j * op[ndim // 4:ndim // 2, 1]
    fkpakpb1 = op[ndim // 2:3 * ndim // 4, 0] + 1j * op[ndim // 2:3 * ndim // 4, 1]
    fkpakpb2 = op[3 * ndim // 4:, 0] + 1j * op[3 * ndim // 4:, 1]

    # bottom layer -> left half (columns 1-2)
    plt.figure(figsize=(20, 11))
    plt.subplots_adjust(left=.03, right=.97, top=.90, bottom=.12, wspace=.35, hspace=.15)
    plt.gcf().text(.5, .98, os.path.basename(filename), ha='center', va='top', fontsize=16)

    plt.subplot(2, 4, 1)
    plt.title("$| \\rho_{KBb, KAb} |$", fontsize=14)
    plt.scatter(x[:ndim // 4], y[:ndim // 4], c=np.abs(fkakb1), marker='h', s=size, cmap='gist_ncar')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 2)
    plt.title("$angle( \\rho_{KBb, KAb} )$", fontsize=14)
    plt.scatter(x[:ndim // 4], y[:ndim // 4], c=np.angle(fkakb1), marker='h', s=size, cmap='hsv', vmin=-np.pi, vmax=np.pi)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 5)
    plt.title("$| \\rho_{K'Bb, K'Ab} |$", fontsize=14)
    plt.scatter(x[:ndim // 4], y[:ndim // 4], c=np.abs(fkpakpb1), marker='h', s=size, cmap='gist_ncar')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 6)
    plt.title("$angle( \\rho_{K'Bb, K'Ab} )$", fontsize=14)
    plt.scatter(x[:ndim // 4], y[:ndim // 4], c=np.angle(fkpakpb1), marker='h', s=size, cmap='hsv',vmin=-np.pi,vmax=np.pi)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    # top layer -> right half (columns 3-4)

    plt.subplot(2, 4, 3)
    plt.title("$| \\rho_{KBt, KAt} |$", fontsize=14)
    plt.scatter(x[ndim//2:3*ndim // 4], y[ndim//2:3*ndim // 4], c=np.abs(fkakb2), marker='h', s=size, cmap='gist_ncar')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 4)
    plt.title("$angle() \\rho_{KBt, KAt} )$", fontsize=14)
    plt.scatter(x[ndim//2:3*ndim // 4], y[ndim//2:3*ndim // 4], c=np.angle(fkakb2), marker='h', s=size, cmap='hsv', vmin=-np.pi, vmax=np.pi)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 7)
    plt.title("$| \\rho_{K'Bt, K'At} |$", fontsize=14)
    plt.scatter(x[ndim//2:3*ndim // 4], y[ndim//2:3*ndim // 4], c=np.abs(fkpakpb2), marker='h', s=size, cmap='gist_ncar')
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')

    plt.subplot(2, 4, 8)
    plt.title("$angle(\\rho_{K'Bt, K'At} )$", fontsize=14)
    plt.scatter(x[ndim//2:3*ndim // 4], y[ndim//2:3*ndim // 4], c=np.angle(fkpakpb2), marker='h', s=size, cmap='hsv',vmin=-np.pi,vmax=np.pi)
    plt.colorbar()
    plt.plot([vk1[0], vk2[0], vk3[0], vk4[0], vk5[0], vk6[0], vk1[0]], [vk1[1], vk2[1], vk3[1], vk4[1], vk5[1], vk6[1], vk1[1]], c='k', alpha=1, linewidth=1, linestyle='-')
    plt.axis('equal')
    plt.axis('off')


    show_order_params([
        ('x00', np.real(np.sum(fkakb1)) + np.real(np.sum(fkpakpb1)) + np.real(np.sum(fkakb2)) + np.real(np.sum(fkpakpb2))),
        ('y00', np.imag(np.sum(fkakb1)) + np.imag(np.sum(fkpakpb1)) + np.imag(np.sum(fkakb2)) + np.imag(np.sum(fkpakpb2))),
        ('xz0', np.real(np.sum(fkakb1)) - np.real(np.sum(fkpakpb1)) + np.real(np.sum(fkakb2)) - np.real(np.sum(fkpakpb2))),
        ('yz0', np.imag(np.sum(fkakb1)) - np.imag(np.sum(fkpakpb1)) + np.imag(np.sum(fkakb2)) - np.imag(np.sum(fkpakpb2))),
        ('x0z', -np.real(np.sum(fkakb1)) - np.real(np.sum(fkpakpb1)) + np.real(np.sum(fkakb2)) + np.real(np.sum(fkpakpb2))),
        ('y0z', -np.imag(np.sum(fkakb1)) - np.imag(np.sum(fkpakpb1)) + np.imag(np.sum(fkakb2)) + np.imag(np.sum(fkpakpb2))),
        ('xzz', -np.real(np.sum(fkakb1)) + np.real(np.sum(fkpakpb1)) + np.real(np.sum(fkakb2)) - np.real(np.sum(fkpakpb2))),
        ('yzz', -np.imag(np.sum(fkakb1)) + np.imag(np.sum(fkpakpb1)) + np.imag(np.sum(fkakb2)) - np.imag(np.sum(fkpakpb2))),
    ])
plt.show()
