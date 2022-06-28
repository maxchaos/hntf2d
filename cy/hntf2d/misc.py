"""Miscellaneous utilities."""

import numpy


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
