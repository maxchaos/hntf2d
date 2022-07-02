
cimport cython

import numpy
cimport numpy as cnumpy
cnumpy.import_array()

from eigency.core cimport *

from libcpp.vector cimport vector
from libcpp.string cimport string

# cdef extern from "Eigen/Dense" namespace "Eigen":
#     cdef cppclass MatrixXf:
#         MatrixXf() except +
#         MatrixXf(int x, int y) except +
#         # double& elem "operator()"(int row, int col)
#         float coeff(int, int)
#         float* data()
#         int rows()
#         int cols()
#     cdef cppclass VectorXf:
#         VectorXf() except +
#         double& element "operator()"(int row, int col)
#         float* data()
#     # cdef cppclass RowMajor:
#     #     pass
#     # cdef cppclass ColMajor:
#     #     pass
#     # cdef cppclass Matrix2f:
#     #     Matrix2f() except +
#     #     Matrix2f(float x, float y) except +

cdef extern from "hntf2d/map.h":
    cdef cppclass _HarmonicMap2D "HarmonicMap2D":
          _HarmonicMap2D()
          # vector[MatrixXf*] boundaries
          # vector[VectorXf*] ku, kv
          # MatrixXf* boundary(size_t index)
          # vector[MatrixXf*].iterator boundaries();
          MatrixXf boundary_get(int index) except +
          void boundary_set(int index, MatrixXf &mtx) except +
          void boundary_append(MatrixXf&) except +
          void print_boundaries()
          VectorXf ku_get(int index) except +
          void ku_set(int index, VectorXf &vec) except +
          VectorXf kv_get(int index) except +
          void kv_set(int index, VectorXf &vec) except +
          Matrix2f jacob(float x, float y)
          Vector2f map(float x, float y)
          MatrixXf mtx_geom() except +
          MatrixXf mtx_univ() except +
          MatrixXf mtx_aug() except +
          void solve(VectorXf uo, VectorXf vo) except +
          MatrixXf puncts() except +
          void load_xml(string &fn) except +
          void save_xml(string &fn) except +

cdef class HarmonicMap2D:

    cdef _HarmonicMap2D* _thisptr

    def __cinit__(self, boundaries=[]):
        self._thisptr = new _HarmonicMap2D()
        for bnd in boundaries:
            self.boundary_append(bnd)

    def __dealloc__(self):
        if self._thisptr != NULL:
            del self._thisptr

    # @property
    # def boundaries(self):
    #     return [ndarray_copy(elt)
    #             for elt in self._thisptr.boundaries()]
    # @boundaries.setter
    # def boundaries(self, val):
    #     pass

    def boundary_get(self, int index):
        # cdef MatrixXf *mtx = self._thisptr.boundaries.data()[index]
        cdef MatrixXf mtx = self._thisptr.boundary_get(index)
        # cdef int rows, cols
        # rows, cols = mtx[0].rows(), mtx[0].cols()
        # cdef float[::1,:] mem_view = <float[:rows:1,:cols]>(mtx[0].data())
        # dtype = 'float32'
        # cdef int itemsize = numpy.dtype(dtype).itemsize
        # return numpy.copy(mem_view, order="F")
        return ndarray_copy(mtx)

    def boundary_set(self, int index, cnumpy.ndarray mtx):
        """ """
        mtx = mtx.astype('float32')
        if not numpy.isfortran(mtx):
            mtx = numpy.asfortranarray(mtx)
        cdef cnumpy.ndarray[float, ndim=2] mtx_ = mtx
        self._thisptr.boundary_set(index, Map[MatrixXf](mtx_).toMatrix())

    def boundary_append(self, cnumpy.ndarray bnd):
        bnd = bnd.astype('float32')
        if not numpy.isfortran(bnd):
            bnd = numpy.asfortranarray(bnd)
        cdef cnumpy.ndarray[float, ndim=2] mtx = bnd
        self._thisptr.boundary_append(Map[MatrixXf](mtx).toMatrix())
        # self._thisptr.boundaries.push_back(Map[MatrixXf](mtx).toMatrixPtr())

    # def get_boundaries(self):
    #     res = [ndarray_copy(elt) for elt in self._thisptr.boundaries]
    #     return res

    def print_boundaries(self):
        self._thisptr.print_boundaries()

    def ku_get(self, int index):
        return ndarray_copy(self._thisptr.ku_get(index))

    def ku_set(self, int index, cnumpy.ndarray ku):
        ku = ku.astype('float32')
        cdef cnumpy.ndarray[float, ndim=1] vec = ku
        self._thisptr.ku_set(index, Map[VectorXf](vec).toMatrix())

    def kv_get(self, int index):
        return ndarray_copy(self._thisptr.kv_get(index))

    def kv_set(self, int index, cnumpy.ndarray kv):
        kv = kv.astype('float32')
        cdef cnumpy.ndarray[float, ndim=1] vec = kv
        self._thisptr.kv_set(index, Map[VectorXf](vec).toMatrix())

    def mtx_geom(self):
        return ndarray_copy(self._thisptr.mtx_geom())
        # self._thisptr.mtx_geom()

    def mtx_univ(self):
        return ndarray_copy(self._thisptr.mtx_univ())

    def mtx_aug(self):
        return ndarray_copy(self._thisptr.mtx_aug())

    def solve(self, cnumpy.ndarray uo, cnumpy.ndarray vo):
        if not (uo.ndim == 1 or uo.shape[0] == 1 or uo.shape[1] == 1):
            msg = "uo is not a vector"
            raise ValueError( msg )
        if not (vo.ndim == 1 or vo.shape[0] == 1 or vo.shape[1] == 1):
            msg = "uo is not a vector"
            raise ValueError( msg )
        uo = uo.astype('float32').flatten()
        vo = vo.astype('float32').flatten()
        self._thisptr.solve(
            Map[VectorXf](uo).toMatrix(), Map[VectorXf](vo).toMatrix()
        )

    # def map(self, float x, float y):
    #     return ndarray_copy(self._thisptr.map(x, y))
    def map(self, *args):
        cdef float x, y
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (numpy.ndarray, tuple, list)):
                x, y = numpy.array(arg).flatten()
        elif len(args) == 2:
            x, y = float(args[0]), float(args[1])
        return ndarray_copy(self._thisptr.map(x, y)).reshape(1,2)

    # def jacob(self, float x, float y):
    #     return ndarray_copy(self._thisptr.jacob(x, y))
    def jacob(self, *args):
        cdef float x, y
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (numpy.ndarray, tuple, list)):
                x, y = numpy.array(arg).flatten()
        elif len(args) == 2:
            x, y = float(args[0]), float(args[1])
        return ndarray_copy(self._thisptr.jacob(x, y))

    def puncts(self):
        cdef MatrixXf puncts = self._thisptr.puncts()
        if puncts.rows() > 0:
            return ndarray_copy(puncts)
        else:
            return None
        # return ndarray_copy(self._thisptr.puncts())

    def load_xml(self, fn):
        cdef fn_ = fn.encode("UTF-8")
        self._thisptr.load_xml(fn_)

    def save_xml(self, fn):
        cdef fn_ = fn.encode("UTF-8")
        self._thisptr.save_xml(fn_)
