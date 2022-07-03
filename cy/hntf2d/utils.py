"""Miscellaneous utilities."""

import numpy
import math
import itertools
import functools


def poly2circle(polyline, collapse_intervals=[], closed_tol=0.0001):
    """Map a given poly-line to the unit circle.

    The poly-line must be given as a numpy array with shape (N, 2),
    where each line corresponds to a vertex of the poly-line;
    each pair of consecutive vertices constitutes a segment of the poly-line.

    Specifically, this function maps the center C of each segment to
    a point on the unit circle using t = l / L as the generalized coordinate,
    where:
    - l : length of the arc formed by the first vertex and the center point C,
    - L : total length of the poly-line.

    Additionally, zero or more intervals of the poly-line vertices
    can be specified such that the vertices each interval will be
    collapsed to the same point on the unit disk.

    The poly-line must be closed, i.e.,
    the first and last vertices must be equal.
    This is checked internally by comparing the norm of their difference to
    to the value of CLOSED_TOL.
    """
    # Check if polygon is closed.
    if numpy.linalg.norm(polyline[0] - polyline[-1]) > closed_tol:
        msg = "polygon is not closed"
        raise ValueError(msg)
    # A closed polygon is defined by at least 4 points.
    if len(polyline) < 4:
        msg = "closed polygon requires at least 4 points"
        raise ValueError(msg)
    if not collapse_intervals:
        lens = list(numpy.linalg.norm(
            polyline[1:, :] - polyline[:-1, :], axis=1))
        lens.insert(0, 0.0)
        clens = numpy.cumsum(lens)
        s = 0.5*(clens[:-1] + clens[1:])/clens[-1]  # interval (0, 1)
        t = 2.0 * numpy.pi * s                      # interval (0, 2*pi)
        ct, st = numpy.cos(t), numpy.sin(t)
        return numpy.array((ct, st)).T
    else:
        res = numpy.zeros((len(polyline)-1, 2))
        collapse_intervals = sorted(
            collapse_intervals,
            key=functools.cmp_to_key(
                lambda x, y: 1 if x[0] > y[0] else -1 if x[0] < y[0] else 0)
        )
        # Add virtual collapse intervals at the beginning and end,
        # if needed.
        if collapse_intervals[0][0] != 0:
            collapse_intervals.insert(0, (-1,0))
        if collapse_intervals[-1][1] != len(polyline)-1:
            collapse_intervals.append((len(polyline)-1,-1))
        print(collapse_intervals)
        rints = [(i[1], j[0]) for i, j in
                 zip(collapse_intervals[:-1], collapse_intervals[1:])]
        # Remove virtual collapse intervals from the beginning and end,
        # if added.
        if collapse_intervals[0][0] == -1:
            collapse_intervals.pop(0)
        if collapse_intervals[-1][1] == -1:
            collapse_intervals.pop(-1)
        # Extract regular parts.
        rparts = [polyline[ri[0]:ri[1]+1] for ri in rints]
        rlens = [
            list(numpy.linalg.norm(rp[1:]-rp[:-1], axis=1))
            for rp in rparts
        ]
        rlens.insert(0, [0.0])
        rclens = numpy.cumsum(tuple(itertools.chain(*rlens)))
        rs = 0.5*(rclens[:-1] + rclens[1:])/rclens[-1]  # interval (0, 1)
        rt = 2.0 * numpy.pi * rs                        # interval (0, 2*pi)
        # Split rt in slices, where each slice corresponds to a normal part.
        rt = numpy.split(rt, numpy.cumsum([len(rp)-1 for rp in rparts]))
        rt = [elt for elt in rt if len(elt)]
        for i, t in zip(rints, rt):
            res[i[0]:i[1]] = numpy.array([numpy.cos(t), numpy.sin(t)]).T
        ct = [(tp[-1]+tn[0])/2 for tp, tn in zip(rt[:-1], rt[1:])]
        if collapse_intervals[0][0] == 0:
            ct.insert(0, 0.0)
        if collapse_intervals[-1][1] == len(polyline)-1:
            ct.append(2.0*math.pi)
        for i, t in zip(collapse_intervals, ct):
            t = numpy.ones(i[1]-i[0]) * t
            res[i[0]:i[1]] = numpy.array([numpy.cos(t), numpy.sin(t)]).T
        return res


class Disk2PlaneMap(object):
    """Map the image of the harmonic navigation transformation to the plane.

    Given a harmonic navigation transformation mapping R^2 to the unit disk,
    an instance of this class maps its image (i.e., the unit disk),
    back to R^2, leaving any punctures as points but, potentially, displaced.
    """

    def __init__(self, tf, radius=1.0):
        """Constructor.

        Parameters
        ----------
        tf : HarmonicMap2D
            Transformation which maps an original domain to the 2D unit disk.

        """
        self._tf = tf
        self._radius = radius

    def map(self, *args):
        """Get the image of a given point of the original domain."""
        r = self._radius
        p = self._tf.map(*args)
        d = numpy.linalg.norm(p)
        q = d/(r-d)*p
        return q

    def jacob(self, *args):
        """Get the composition's Jacobian at a specified point."""
        r = self._radius
        p = self._tf.map(*args)
        d = numpy.linalg.norm(p)
        x, y = p.flatten()
        x2, y2 = x*x, y*y
        m = numpy.array([[-d*x2 - d*y2 + 2*r*x2 + r*y2, r*x*y],
                         [r*x*y, -d*x2 - d*y2 + r*x2 + 2*r*y2]])
        m = m/(d*(d-r)*(d-r))
        return m.dot(self._tf.jacob(*args))

    def puncts(self):
        """Get the finale location of the unit disk's punctures."""
        r = self._radius
        puncts = self._tf.puncts()
        res = numpy.zeros_like(puncts)
        for k, p in enumerate(puncts):
            d = numpy.linalg.norm(p)
            q = d/(r-d)*p
            res[k] = q
        return res
