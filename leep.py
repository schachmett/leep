#!/usr/bin/env python3
"""
a: real
b: reciprocal vector
"""


import sys
from signal import signal, SIGINT
import re
import argparse
import copy

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


DEFAULTS = {
    "hex": [[1, 0], [np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3)]],
    "square": [[1, 0], [0, 1]],
    "rect": [[1, 0], [0, 2]],
    "centered": [[1, 0], [0.5, 0.3]],
    "monoclinic": [[1, 0], [0.2, 1.3]],

    "ru": [[2.71, 0], [2.71 * np.cos(2 * np.pi / 3), 2.71 * np.sin(2 * np.pi / 3)]],
    "ruo2-110": [[3.11, 0], [0, 6.38]],
    "ruo2-100": [[3.11, 0], [0, 4.49]],
    "vo2-100": [[0, -2.8514], [4.5546, 0]],
    "vo2-110": [[0, -2.8514], [4.5546 * np.sqrt(2), 0]]
}


def main():
    # parse arguments
    args = parse_args(sys.argv[1:])

    # initialize figure title to nothing (yet)
    title = ""

    # make a list of all phases that should be plotted
    phases = []

    if args.special:
        uc_ru = UC(*DEFAULTS["ru"], s=1)
        phases.append(uc_ru)
        if "110" in args.special[0].lower():
            if "ruo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["ruo2-110"], parent=uc_ru, color="b", s=0.6))
            if "vo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["vo2-110"], parent=uc_ru, color="b", s=0.6))
        if "100" in args.special[0].lower():
            if "ruo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["ruo2-100"], parent=uc_ru, color="b", s=0.6))
            if "vo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["vo2-100"], parent=uc_ru, color="b", s=0.6))
        for phase in phases[1:].copy():
            for color, theta in zip(("g", "r"), (2 * np.pi / 3, 4 * np.pi / 3)):
                phases.append(phase.get_reconstruction(theta=theta, color=color, label="rot_dom"))
        plot_phases(phases, title=args.special[0])
        plt.show()
        return

    if args.vectors is not None:
        # if the unit vectors are given in the command line, use them to make the UC
        uc = UC(args.vectors[:2], args.vectors[2:], label="substrate")
        title += "a1=" + str(args.vectors[:2]) + ", a2=" + str(args.vectors[2:])
    else:
        # else, use one of the bases defined in DEFAULTS
        uc = UC(*DEFAULTS[args.base[0]], label="substrate")
        title += args.base[0]
    # add the UC to the list of phases
    phases.append(uc)

    # calculate the list of rotation angles to use from the number of rotational domains
    thetas = np.linspace(0, 2 * np.pi, args.rotations[0] + 1)[1:-1]

    # parse the "matrix" argument into actual matrices
    matrices = list([matrix[:2], matrix[2:]] for matrix in args.matrix)

    for rec in args.woods + matrices:
        # make a reconstruction unit cell from rec, which is either a woods
        # notation or a matrix
        reconstruction = uc.get_reconstruction(rec, color="b", s=0.6, label="rec")
        title += " | " + str(rec) + " reconstr."
        # add the reconstruction to the list of phases
        phases.append(reconstruction)
        for theta in thetas:
            # if rotational domains are required, add them
            rot_dom = reconstruction.get_reconstruction(theta=theta, label="rot_dom")
            phases.append(rot_dom)

    if args.rotations[0] > 1:
        title += " | " + str(args.rotations[0]) + " rot. dom."

    # do the plotting
    plot_phases(phases, title=title)

    if args.output:
        # save an image if the argument is given
        plt.savefig(args.output)
    else:
        # show the plot
        plt.show()



def parse_args(arglist):
    """Parse arguments and define helptext."""
    helptext = (
        "leep.py visualizes reconstructions. If -b and -v are both"
        "not given, it assumes a square unit mesh"
    )
    parser = argparse.ArgumentParser(description=helptext)
    base_group = parser.add_mutually_exclusive_group()
    base_group.add_argument(
        "-b", "--base", nargs=1, default=["square"], type=str,
        help="Real space base. Either one of 'hex,square,rect,centered,"
             "monoclinic'.", metavar="BASE")
    base_group.add_argument(
        "-v", "--vectors", nargs=4, type=float,
        help="Explicit base vectors explicitly (instead of 'BASE').",
        metavar=("X1", "Y1", "X2", "Y2"))
    base_group.add_argument(
        "-s", "--special", nargs=1, type=str,
        help="Special cases instead of a bravais base or custom vectors.")

    parser.add_argument(
        "-w", "--woods", action="append", type=str, default=[],
        help="A reconstruction in wood's notation, e.g. 'c(2x2)'"
             "or 'sqrt3xsqrt3 R30' or '2x1'.")
    parser.add_argument(
        "-m", "--matrix", action="append", type=int, nargs=4, default=[],
        help="A reconstruction in matrix notation.",
        metavar=("A11", "A12", "A21", "A22"))
    parser.add_argument(
        "-r", "--rotations", type=int, nargs=1, default=[1],
        help="Number of rotational domains. Applies to all defined reconstructions.")
    parser.add_argument(
        "-o", "--output", type=str, help="Save image file.",
        metavar="fname")
    args = parser.parse_args(arglist)
    return args


def plot_phases(phases, title=""):
    # make 2 axes with a 1:1 aspect ratio
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for ax in axes.flatten():
        ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle(title)

    # plot each phases' base and pattern in reciprocal and in real space
    for phase in phases:
        if not "rot_dom" in phase.label:
            phase.plot_pattern(ax=axes[0])
            phase.plot_base(ax=axes[0])
        phase.plot_pattern(rec=True, ax=axes[1])
        phase.plot_base(rec=True, ax=axes[1])


class UC:
    """
    Unit Cell object. It is defined by its three base vectors a1, a2 and a3 which
    are combined in the _real_base attribute.
    If only two 2D-vectors are given, the third component is made 0 (see _make_vector_3d())
    and the third base vector a3 is set to (0, 0, 1) so that all 3D calculations can
    be used.
    The reciprocal base _rec_base is calculated as the matrix inverse of the real base.
    If a "parent" unit cell is given, its size can be used for the axes dimensions when plotting.
    The color and size of the pattern spots are given by color and s.
    """
    def __init__(self, a1, a2, a3=None, parent=None, color="k", s=1, label=""):
        if a3 is None:
            a3 = np.array([0, 0, 1])
        a1 = _make_vector_3d(a1)
        a2 = _make_vector_3d(a2)
        a3 = _make_vector_3d(a3)

        self._real_base = np.array([a1, a2, a3]).transpose()
        self._rec_base = np.linalg.inv(self.real_base).transpose() * 2 * np.pi

        self.label = label
        if parent is None:
            parent = self
        self.parent = parent
        self.color = color
        self.s = s

    def coords(self, uvw):
        """Transform a lattice vector into cartesian coordinates."""
        uvw = _make_vector_3d(uvw)
        return np.dot(self.real_base, uvw)

    def rec_coords(self, hkl):
        """Transform a reciprocal lattice vector into cartesian coordinates."""
        hkl = _make_vector_3d(hkl)
        return np.dot(self.rec_base, hkl)

    @property
    def size(self):
        """Size of the unit cell."""
        a1, a2, a3 = self.real_base.transpose()
        return np.dot(a1, np.cross(a2, a3))
        
    def d(self, hkl):
        """Lattice plane distance"""
        return 2 * np.pi / np.linalg.norm(self.rec_coords(hkl))

    @property
    def is_2d(self):
        """Is True only if the first 2 base vectors have 0 z component
        and the third base vector is the unit vector in z direction."""
        a1a2_component3 = (self.real_base[:2, 2] == 0).all()
        a3_2d = (self.real_base[2, :] == np.array([0, 0, 1])).all()
        return a1a2_component3 and a3_2d

    @staticmethod
    def _apply_matrix(base, M):
        """Apply a reconstruction matrix to a base"""
        return np.dot(base, M.transpose())

    @staticmethod
    def _apply_rotation(base, theta):
        """Rotate a base by theta (given in radians)"""
        if theta != 0:
            R = np.array([
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [0, 0, 1]
            ])
            return np.dot(R, base)
        return base

    def get_reconstruction(self, M="1x1", theta=0, **kwargs):
        """Returns a "child" unit cell with the base transformed by M.
        M can either be a reconstruction matrix or a woods notation.
        """
        # only works if the unit cell is 2D
        assert self.is_2d
        # if M is a string, this will make it into a transformation matrix.
        # also, 2x2 matrices are expanded into 3x3 matrices
        M, theta = _parse_woods(M, theta)
        # first the matrix, then the rotation is applied
        base = self._apply_matrix(self.real_base, M)
        base = self._apply_rotation(base, theta)

        # by default, the color and spot size are given to the child UC
        sub_kwargs = dict(color=self.color, s=self.s)
        # if color and s are given in kwargs, they override the default
        sub_kwargs.update(kwargs)

        cell = UC(*base.transpose(), parent=self.parent, **sub_kwargs)
        return cell

    def plot_base(self, rec=False, ax=None):
        """Plots the base vectors as arrows, either in reciprocal or real space."""
        assert self.is_2d
        if rec:
            base = self.rec_base
        else:
            base = self.real_base

        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect("equal")

        ax.annotate(
            "", xy=(base[0, 0], base[1, 0]), xytext=(0, 0),
            arrowprops=dict(arrowstyle="-|>, head_width=0.4, head_length=1", color=self.color))
        ax.annotate(
            "", xy=(base[0, 1], base[1, 1]), xytext=(0, 0),
            arrowprops=dict(arrowstyle="-|>, head_width=0.4, head_length=1", color=self.color))
        ax.annotate(
            "", xy=(base[0, 0] + base[0, 1], base[1, 0] + base[1, 1]), xytext=(base[0, 0], base[1, 0]),
            arrowprops=dict(arrowstyle="-", linestyle="--", color=self.color))
        ax.annotate(
            "", xy=(base[0, 0] + base[0, 1], base[1, 0] + base[1, 1]), xytext=(base[0, 1], base[1, 1]),
            arrowprops=dict(arrowstyle="-", linestyle="--", color=self.color))
        return ax

    def plot_pattern(self, rec=False, ax=None, ucs=None, d=None):
        """Plots the symmetry pattern generated by the unit cell.
        In reciprocal space, it will show the LEED pattern."""
        assert self.is_2d
        if ucs is None:
            if rec is True:
                ucs = 3
            else:
                ucs = 10
        
        # find the minimal lattice plane distance and define the base
        base = self.real_base
        if d is None:
            d = min(self.parent.d([1, 0]), self.parent.d([0, 1]))
        if rec:
            base = self.rec_base
            d = 2 * np.pi / d
        size = ucs * d

        # maximal indices that can be in the image
        i1, j1, _ = np.dot(np.linalg.inv(base), np.array([size, size, 0]))
        i2, j2, _ = np.dot(np.linalg.inv(base), np.array([-size, size, 0]))
        i_max, j_max = int(max(np.abs([i1, i2])) / 2), int(max(np.abs([j1, j2])) / 2)

        # minimal distances of atoms in x and y-direction, detect if atom is (barely) in image
        xs, ys = abs(max(*base[0, :2], key=abs)), abs(max(*base[1, :2], key=abs))
        def inside(p):
            if abs(p[0]) < size / 2 + xs and abs(p[1]) < size / 2 + ys:
                return True
            return False

        # find points to plot
        points = []
        for i in range(-i_max, i_max + 1):
            for j in range(-j_max, j_max + 1):
                y = np.dot(base, [i, j, 0])
                if inside(y):
                    points.append(y)
        points = np.unique(np.array(points), axis=0)

        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect("equal")
        ax.set_xlim(-size/2, size/2)
        ax.set_ylim(-size/2, size/2)

        if rec:
            ax.scatter(points[:, 0], points[:, 1], color=self.color, s=self.s**3 * 250 / ucs)
            ax.set_title("Reciprocal space")
        else:
            gauplot([(x, y) for x, y, _ in points], self.s * d, size=size*1.2, ax=ax, color=self.color)
            ax.set_title("Real space")
        return ax

    @property
    def real_base(self):
        return self._real_base
    @property
    def rec_base(self):
        return self._rec_base


def gauplot(centers, radius=5, size=20, ax=None, color="r", thresh=2):
    """adapted from https://stackoverflow.com/questions/10958835"""
    color = colors.to_rgb(color)
    cmap_mat = np.array([np.linspace(c, 1, 256) for c in color]).T
    cmap = colors.ListedColormap(cmap_mat)
    cmap.set_bad(alpha=0)
    
    nx, ny = 1000.,1000.
#    xgrid, ygrid = np.mgrid[-m:m:m*2/nx, -m:m:m*2/ny]
    xgrid, ygrid = np.mgrid[-size/2:size/2:size/nx, -size/2:size/2:size/ny]
    img = xgrid * 0 + np.nan
    for center in centers:      # center[1]: y axis is negative, flip sign
        atom = np.sqrt((xgrid - center[0])**2 + (ygrid + center[1])**2) / (radius / 2) * thresh
        img[atom < thresh] = np.exp(-.5 * atom**2)[atom < thresh]
    ax.imshow(img.T, cmap=cmap, extent=(-size/2, size/2, -size/2, size/2))


def _make_vector_3d(x, third=0):
    """Make sure that a vector is 3D. If it is 2D, make the third component zero."""
    x = np.array(x)
    if len(x) == 2:
        x = np.array([*x, third])
    assert len(x) == 3
    return x

def _parse_woods(M, theta=0):
    """Parse a woods notation string into a 3x3 transformation matrix."""
    if isinstance(M, str):
        a = re.findall("([\.sqrt0-9]*)\s*x", M)[0]
        a0 = float(re.findall("\d+\.*\d*", a)[0])
        if any(s in a for s in ("s", "r", "sqrt")):
            a = np.sqrt(a0)
        else:
            a = a0

        b = re.findall("x\s*([\.sqrt0-9]*)", M)[0]
        b0 = float(re.findall("\d+\.*\d*", b)[0])
        if any(s in b for s in ("s", "r", "sqrt")):
            b = np.sqrt(b0)
        else:
            b = b0

        t = re.findall("[Rr]\s*([-0-9]*)", M)
        if t:
            theta = np.pi * float(t[-1]) / 180

        if M.startswith("c"):
            if a != b:
                raise ValueError(f"Unknown reconstruction {M}")
            M = np.array([[a, 0], [1, 1]])
        else:
            M = np.array([[a, 0], [0, b]])

    M = np.pad(M, ((0, 1), (0, 1)), constant_values=0)
    M[2, 2] = 1
    return M, theta

def norm(x):
    """Return the euclidian length of a vector"""
    return np.sqrt(sum(xi**2 for xi in x))


def sigint_handler(signal_received, frame):
    print("\nExiting...")
    sys.exit()


if __name__ == "__main__":
    signal(SIGINT, sigint_handler)
    main()

