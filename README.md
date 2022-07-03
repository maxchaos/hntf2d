
# Table of Contents

1.  [Synopsis](#sec/synopsis)
2.  [Installation](#org5ec5c9d)
    1.  [Shared Library](#org5837930)
    2.  [Python Wrapper through Cython](#org7337359)
3.  [Use](#orga221729)
    1.  [Cython wrapper](#org051a8d4)
4.  [Guix](#orgd4ec3b2)
5.  [External References](#sec/external_references)



<a id="sec/synopsis"></a>

# Synopsis

This software package provides an implementation of
a harmonic transformation that maps any compact, multiply connected planar
domain onto the unit disk;
each of the holes that the original domain may have are collapsed onto
distinct points inside its image.

This transformation has been employed for addressing the navigation problem
of a single planar robot operating in arbitrarily shaped workspaces in
[1.1](#orga38c052), [1.2](#org28b0dc0).

If you end up using this transformation and/or this software package,
please cite the aforementioned paper as follows:

    P. Vlantis, C. Vrohidis, C. P. Bechlioulis and K. J. Kyriakopoulos, "Robot Navigation in Complex Workspaces Using Harmonic Maps," 2018 IEEE International Conference on Robotics and Automation (ICRA), 2018, pp. 1726-1731, doi: 10.1109/ICRA.2018.8460695.


<a id="org5ec5c9d"></a>

# Installation

In general,
this package provides the following two distinct outputs:

-   a shared library implemented in C++,
-   a python wrapper of the aforementioned library.

All these outputs depend on the following libraries,
which must be installed on the system in order to build them:

-   [eigen3](#org3a0a4f9),
-   [tinyxml2](#org2df8873).

The recommended ways to build each of these outputs are described in
the following subsections.


<a id="org5837930"></a>

## Shared Library

In order to build the shared library from the contained C++ source code
using CMake, execute the following steps:

1.  Create sub-directory `build/cpp/` inside this package's root directory.

2.  Inside `build/cpp`, execute the following shell command:
    
        cmake -DCMAKE_INSTALL_PREFIX=prefix ../.. && make && make install
    
    where `prefix` is a placeholder for the path to the directory
    where the library and C++ headers will be installed at
    (defaults to `/usr/local/` on GNU/Linux systems).

After executing these steps,
the library and header files should be installed as follows:

    prefix/
    ├── include
    │   └── hntf2d
    │       └── ...
    └── lib
        └── libhntf2d.so


<a id="org7337359"></a>

## Python Wrapper through Cython

A python wrapper for the C++ implementation is also provided using Cython.

In order to build the python wrapper,
the following procedure should be executed:

1.  Download a customized version of the software package [eigency](#orgbce8c56)
    located at [2.4](#org305b335) by executing the following shell command inside
    this package's root directory:
    
        git submodule update --init --remote --recursive

2.  Install the following dependencies:
    
    -   `numpy`
    -   `cython`
    
    This can be done by executing the following shell commands:
    
        python3 -m pip install cython numpy

3.  Build and install `eigency`.
    This can be done by executing the following shell command inside
    the directory `cy/eigency`:
    
        python3 setup.py build_ext --inplace
        python3 setup.py install

4.  Finally, build and install the wrapper.
    This can be done  by executing the following shell command inside
    the directory `cy/hntf2d`:
    
        python3 setup.py build_ext --inplace
        python3 setup.py install


<a id="orga221729"></a>

# Use


<a id="org051a8d4"></a>

## Cython wrapper

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


<a id="orgd4ec3b2"></a>

# Guix

In order to develop this package in a Guix shell with
all the necessary dependencies,
the following manifest file can be used.

    (use-modules (gnu packages)
                 (guix packages)
                 ((gnu packages cmake) :prefix gnu:)
                 ((gnu packages llvm) :prefix gnu:)
                 ((gnu packages algebra) :prefix gnu:)
                 ((gnu packages python) :prefix gnu:)
                 ((gnu packages python-xyz) :prefix gnu:)
                 ((gnu packages xml) :prefix gnu:))
    
    (packages->manifest (list gnu:cmake
                              gnu:clang
                              gnu:python-scikit-build
                              gnu:eigen
                              gnu:tinyxml2
                              gnu:python
                              gnu:python-cython
                              gnu:python-numpy
                              gnu:python-setuptools-scm))

    guix shell -m manifest.scm


<a id="sec/external_references"></a>

# External References

1.  Papers:
    1.  <a id="orga38c052"></a> [scholar.google.com &#x2013; Robot navigation in complex workspaces using harmonic maps](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=R5c4qS8AAAAJ&citation_for_view=R5c4qS8AAAAJ:u-x6o8ySG0sC)
    2.  <a id="org28b0dc0"></a> [ieeexplore.ieee.org &#x2013; Robot Navigation in Complex Workspaces Using Harmonic Maps](https://ieeexplore.ieee.org/abstract/document/8460695)
2.  Software Packages:
    1.  <a id="org3a0a4f9"></a> [eigen.tuxfamily.org](https://eigen.tuxfamily.org/index.php?title=Main_Page)
    2.  <a id="org2df8873"></a> [leethomason.github.io &#x2013; tinyxml2](http://leethomason.github.io/tinyxml2/)
    3.  <a id="orgbce8c56"></a> [github.com &#x2013; wouterboomsma/eigency](https://github.com/wouterboomsma/eigency)
    4.  <a id="org305b335"></a> [github.com &#x2013; maxchaos/eigency](https://github.com/maxchaos/eigency)

