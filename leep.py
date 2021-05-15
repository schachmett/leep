#!/usr/bin/env python3
"""
a: real
b: reciprocal vector
"""


import sys
from signal import signal, SIGINT
import re
import argparse

import numpy as np
import matplotlib.pyplot as plt


DEFAULTS = {
    "hex": [[1, 0], [np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3)]],
    "square": [[1, 0], [0, 1]],
    "rect": [[1, 0], [0, 2]],
    "centered": [[1, 0], [0.5, 0.3]],
    "monoclinic": [[1, 0], [0.2, 1.3]],

    "ru": [[2.71, 0], [2.71 * np.cos(2 * np.pi / 3), 2.71 * np.sin(2 * np.pi / 3)]],
    "ruo2-110": [[3.11, 0], [0, 6.38]],
    "ruo2-100": [[3.11, 0], [0, 4.49]],
    "vo2-100": [[2.8514, 0], [0, 4.5546]],
    "vo2-110": [[2.8514, 0], [0, 4.5546 * np.sqrt(2)]]
}


def main():
    # parse arguments
    args = parse_args(sys.argv[1:])

    # initialize figure title to nothing (yet)
    title = ""

    # make a list of all phases that should be plotted
    phases = []

    if args.special:
        uc_ru = UC(*DEFAULTS["ru"], s=10)
        if "110" in args.special[0].lower():
            if "ruo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["ruo2-110"], parent=uc_ru, color="r", s=30))
            if "vo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["vo2-110"], parent=uc_ru, color="r", s=30))
        if "100" in args.special[0].lower():
            if "ruo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["ruo2-100"], parent=uc_ru, color="r", s=30))
            if "vo2" in args.special[0].lower():
                phases.append(UC(*DEFAULTS["vo2-100"], parent=uc_ru, color="r", s=30))
        for phase in phases.copy():
            for color, theta in zip(("g", "b"), (2 * np.pi / 3, 4 * np.pi / 3)):
                phases.append(phase.get_reconstruction(theta=theta, color=color))
        phases.append(uc_ru)
        plot_phases(phases, title=args.special[0])
        plt.show()
        return

    if args.vectors is not None:
        # if the unit vectors are given in the command line, use them to make the UC
        uc = UC(args.vectors[:2], args.vectors[2:])
        title += "a1=" + str(args.vectors[:2]) + ", a2=" + str(args.vectors[2:])
    else:
        # else, use one of the bases defined in DEFAULTS
        uc = UC(*DEFAULTS[args.base[0]])
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
        reconstruction = uc.get_reconstruction(rec, color="b", s=30)
        title += " | " + str(rec) + " reconstr."
        # add the reconstruction to the list of phases
        phases.append(reconstruction)
        for theta in thetas:
            # if rotational domains are required, add them
            phases.append(reconstruction.get_reconstruction(theta=theta))

    if args.rotations[0] > 1:
        title += " | " + str(args.rotations[0]) + " rot. dom."

    # do the plotting
    plot_phases(phases, title=title)

    # save an image if the argument is given
    if args.output:
        plt.savefig(args.output)

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
    # make 4 axes with a 1:1 aspect ratio
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    for ax in axes.flatten():
        ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle(title)

    # plot each phases' base and pattern in reciprocal and in real space
    for phase in phases:
        phase.plot_base(ax=axes[0, 0])
        phase.plot_base(rec=True, ax=axes[1, 0])
        phase.plot_pattern(ax=axes[0, 1])
        phase.plot_pattern(rec=True, ax=axes[1, 1])


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
    def __init__(self, a1, a2, a3=None, parent=None, color="k", s=10):
        if a3 is None:
            a3 = np.array([0, 0, 1])
        a1 = _make_vector_3d(a1)
        a2 = _make_vector_3d(a2)
        a3 = _make_vector_3d(a3)

        self._real_base = np.array([a1, a2, a3]).transpose()
        self._rec_base = np.linalg.inv(self.real_base).transpose() * 2 * np.pi

        if parent is None:
            parent = self
        self.parent = parent
        self.color = color
        self.s = s

    def coords(self, x):
        """Transform a lattice vector into cartesian coordinates."""
        x = _make_vector_3d(x)
        return np.dot(self.real_base, x)

    def rec_coords(self, y):
        """Transform a reciprocal lattice vector into cartesian coordinates."""
        y = _make_vector_3d(y)
        return np.dot(self.rec_base, y)

    @property
    def size(self):
        """Size of the unit cell."""
        a1, a2, a3 = self.real_base.transpose()
        return np.dot(a1, np.cross(a2, a3))

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
            ax.set_title("Reciprocal base vectors")
        else:
            base = self.real_base
            ax.set_title("Real base vectors")

        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect("equal")

        ax.scatter(base[0, :2], base[1, :2], color="w")
        ax.scatter([0], [0], color="k")

        ax.annotate(
            "", xy=(base[0, 0], base[1, 0]), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color=self.color))
        ax.annotate(
            "", xy=(base[0, 1], base[1, 1]), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color=self.color))
        return ax

    def plot_pattern(self, rec=False, ax=None, w=None, h=None):
        """Plots the symmetry pattern generated by the unit cell.
        In reciprocal space, it will show the LEED pattern."""
        assert self.is_2d
        if rec:
            base = self.rec_base
            ax.set_title("Reciprocal space")
            s = 500 / self.s
            if w is None:
                w = 6 * np.pi / np.sqrt(self.parent.size)
            if h is None:
                h = 6 * np.pi / np.sqrt(self.parent.size)
        else:
            base = self.real_base
            ax.set_title("Real space")
            s = self.s
            if w is None:
                w = 10 * np.sqrt(self.parent.size)
            if h is None:
                h = 10 * np.sqrt(self.parent.size)

        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect("equal")

        x1, x2, _ = base.transpose()

        i1, j1, _ = map(int, np.dot(np.linalg.inv(base), np.array([w, h, 0])))
        i2, j2, _ = map(int, np.dot(np.linalg.inv(base), np.array([-w, h, 0])))
        i_max = max(np.abs([i1, i2]))
        j_max = max(np.abs([j1, j2]))

        points = []
        for i in range(-i_max, i_max + 1):
            for j in range(-j_max, j_max + 1):
                y = i * x1 + j * x2
                if -w / 2 <= y[0] <= w / 2 and -h / 2 <= y[1] <= h / 2:
                    points.append(y)
        points = np.unique(np.array(points), axis=0)
        ax.scatter(points[:, 0], points[:, 1], color=self.color, s=s)
        ax.set_xlim(-w/2, w/2)
        ax.set_ylim(-h/2, h/2)
        return ax

    @property
    def real_base(self):
        return self._real_base
    @property
    def rec_base(self):
        return self._rec_base


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

