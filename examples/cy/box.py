import hntf2d.map
import hntf2d.utils
import numpy

hm = hntf2d.map.HarmonicMap2D()

# XXX: All boundaries must be closed curves

# Add a box-shaped outer boundary:
# XXX: Outer boundary must be counter-clock-wise oriented.
t = numpy.linspace(-1., +1., 101)[:100]
outer_boundary = numpy.zeros((400, 2))
# Outer boundary: line segment ((-1, -1) -> (1, -1))
outer_boundary[:100, 0] = t[:]
outer_boundary[:100, 1] = -1
# Outer boundary: line segment ((1, -1) -> (1, 1))
outer_boundary[100:200, 0] = 1
outer_boundary[100:200, 1] = t[:]
# Outer boundary: line segment ((1, 1) -> (-1, 1))
outer_boundary[200:300, 0] = t[::-1]
outer_boundary[200:300, 1] = 1
# Outer boundary: line segment ((-1, 1) -> (-1, -1))
outer_boundary[300:400, 0] = -1
outer_boundary[300:400, 1] = t[::-1]
# Outer boundary: first and last vertices must be equal
outer_boundary[-1, :] = outer_boundary[0, :]

hm.boundary_append(outer_boundary)

# Add a box-shaped inner boundary:
# XXX: Inner boundaries must be clock-wise oriented.
t = numpy.linspace(-1., +1., 21)[:20]
inner_boundary = numpy.zeros((80, 2))
# Inner boundary: line segment ((-0.1, -0.1) -> (-0.1, 0.1))
inner_boundary[:20, 0] = -0.1
inner_boundary[:20, 1] = t[:]
# Outer boundary: line segment ((-0.1, 0.1) -> (0.1, 0.1))
inner_boundary[20:40, 0] = t[:]
inner_boundary[20:40, 1] = 0.1
# Outer boundary: line segment ((0.1, 0.1) -> (0.1, -0.1))
inner_boundary[40:60, 0] = 0.1
inner_boundary[40:60, 1] = t[::-1]
# Outer boundary: line segment ((0.1, -0.1) -> (-0.1, -0.1))
inner_boundary[60:80, 0] = t[::-1]
inner_boundary[60:80, 1] = -0.1
# Outer boundary: first and last vertices must be equal
inner_boundary[-1, :] = inner_boundary[0, :]

hm.boundary_append(inner_boundary)

# Print the registered boundaries (1st it the outer, followed by inner ones).
hm.print_boundaries()

# Map the center of each segment belonging to the outer boundary to
# the unit circle:
uo, vo = hntf2d.utils.poly2circle(outer_boundary).T

hm.solve(uo, vo)

print(hm.map(0., 0.5))
print(hm.jacob(0., 0.5))
