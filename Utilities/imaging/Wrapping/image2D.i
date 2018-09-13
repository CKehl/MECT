%module Image2D

%{
#define SWIG_FILE_WITH_INIT
#include "../../imaging/src/image2D.hpp"
#include "../../imaging/src/image2D_utils.hpp"
#include "../../imaging/src/image.hpp"
#include "../../common/src/std_typedefs.h"
%}

%feature("autodoc", "1");

%include <std_iostream.i>
%include <std_sstream.i>
%include <std_vector.i>
%include <typemaps.i>
%include "numpy.i"
%init %{
	import_array();
%}

// ------------------------------------------------
// raw c array typemaps
// ------------------------------------------------
//%numpy_typemaps(int8_t, NPY_BYTE, int)
//%numpy_typemaps(uint8_t, NPY_UBYTE, int)
//%numpy_typemaps(int16_t, NPY_SHORT, int)
//%numpy_typemaps(uint16_t, NPY_USHORT, int)
//%numpy_typemaps(int32_t, NPY_INT, int)
//%numpy_typemaps(uint32_t, NPY_UINT, int)
//%numpy_typemaps(long, NPY_LONG, int)
//%numpy_typemaps(unsigned long, NPY_ULONG, int)
//%numpy_typemaps(int64_t, NPY_LONGLONG, int)
//%numpy_typemaps(uint64_t, NPY_ULONGLONG, int)
//%numpy_typemaps(float, NPY_FLOAT, int)
//%numpy_typemaps(double, NPY_DOUBLE, int)

// original source: https://gitlab.lam.fr/jclamber/unsio/blob/6969997e086c4041a00ad4d3e0b6000a481e7c21/py/swig/py_unsio.i
// below we change the type of dimension of numpy array.
// By default it's an "int", but we have std::vector object returning size()
// method as dimension, and which are of type "unsigned int"
// This may cause a problem with BIG simulations. Indeed a simulation with 1.1 billions particles
// store in a xyz positions array a size of 3.3 billions which overtake size of SIGNED int, which go
// only to 2^31 bits = 2.14 billons !!!!! (damn nasty bug....)
%numpy_typemaps(double       , NPY_DOUBLE , unsigned int)
%numpy_typemaps(float        , NPY_FLOAT  , unsigned int)
%numpy_typemaps(unsigned long, NPY_ULONG  , unsigned int)
%numpy_typemaps(long         , NPY_LONG   , unsigned int)
%numpy_typemaps(unsigned int , NPY_UINT   , unsigned int)
%numpy_typemaps(size_t       , NPY_ULONG   , unsigned int)
%numpy_typemaps(int          , NPY_INT    , unsigned int)
%numpy_typemaps(unsigned short ,NPY_USHORT, unsigned int)
%numpy_typemaps(short        , NPY_SHORT  , unsigned int)
%numpy_typemaps(double       , NPY_DOUBLE , int)
%numpy_typemaps(float        , NPY_FLOAT  , int)
%numpy_typemaps(unsigned long, NPY_ULONG  , int)
%numpy_typemaps(long         , NPY_LONG   , int)
%numpy_typemaps(unsigned int , NPY_UINT   , int)
%numpy_typemaps(size_t       , NPY_ULONG   , int)
%numpy_typemaps(int          , NPY_INT    , int)
%numpy_typemaps(unsigned short ,NPY_USHORT, int)
%numpy_typemaps(short        , NPY_SHORT  , int)
%numpy_typemaps(double       , NPY_DOUBLE , unsigned long)
%numpy_typemaps(float        , NPY_FLOAT  , unsigned long)
%numpy_typemaps(unsigned long, NPY_ULONG  , unsigned long)
%numpy_typemaps(long         , NPY_LONG   , unsigned long)
%numpy_typemaps(unsigned int , NPY_UINT   , unsigned long)
%numpy_typemaps(size_t       , NPY_ULONG   , unsigned long)
%numpy_typemaps(int          , NPY_INT    , unsigned long)
%numpy_typemaps(unsigned short ,NPY_USHORT, unsigned long)
%numpy_typemaps(short        , NPY_SHORT  , unsigned long)
%numpy_typemaps(double       , NPY_DOUBLE , long)
%numpy_typemaps(float        , NPY_FLOAT  , long)
%numpy_typemaps(unsigned long, NPY_ULONG  , long)
%numpy_typemaps(long         , NPY_LONG   , long)
%numpy_typemaps(unsigned int , NPY_UINT   , long)
%numpy_typemaps(size_t       , NPY_ULONG   , long)
%numpy_typemaps(int          , NPY_INT    , long)
%numpy_typemaps(unsigned short ,NPY_USHORT, long)
%numpy_typemaps(short        , NPY_SHORT  , long)

//%np_vector_typemaps(short, NPY_SHORT);
//%np_vector_typemaps(int, NPY_INT);
//%np_vector_typemaps(size_t, NPY_UINT);
//%np_vector_typemaps(unsigned int, NPY_UINT);
//%np_vector_typemaps(long, NPY_LONG);
//%np_vector_typemaps(float, NPY_FLOAT);
//%np_vector_typemaps(double, NPY_DOUBLE);

// apply numpy typemaps based on arg types + names
%apply (unsigned int* IN_ARRAY1, unsigned int DIM1) {(unsigned int* dimensions_in, unsigned int field_size)};
%apply (double* IN_ARRAY1, unsigned long DIM1) {(double* data_in, unsigned long field_size)};
%apply (float* IN_ARRAY1, unsigned long DIM1) {(float* data_in, unsigned long field_size)};
%apply (unsigned long* IN_ARRAY1, unsigned long DIM1) {(unsigned long* data_in, unsigned long field_size)};
%apply (long* IN_ARRAY1, unsigned long DIM1) {(long* data_in, unsigned long field_size)};
%apply (unsigned int* IN_ARRAY1, unsigned long DIM1) {(unsigned int* data_in, unsigned long field_size)};
%apply (int* IN_ARRAY1, unsigned long DIM1) {(int* data_in, unsigned long field_size)};
%apply (unsigned short* IN_ARRAY1, unsigned long DIM1) {(unsigned short* data_in, unsigned long field_size)};
%apply (short* IN_ARRAY1, unsigned long DIM1) {(short* data_in, unsigned long field_size)};
%apply (unsigned char* IN_ARRAY1, unsigned long DIM1) {(unsigned char* data_in, unsigned long field_size)};
%apply (char* IN_ARRAY1, unsigned long DIM1) {(char* data_in, unsigned long field_size)};

%apply (unsigned int* ARGOUT_ARRAY1, unsigned int DIM1) {(unsigned int* dimensions_out, unsigned int field_size)};
%apply (double* ARGOUT_ARRAY1, unsigned long DIM1) {(double* data_out, unsigned long field_size)};
%apply (float* ARGOUT_ARRAY1, unsigned long DIM1) {(float* data_out, unsigned long field_size)};
%apply (unsigned long* ARGOUT_ARRAY1, unsigned long DIM1) {(unsigned long* data_out, unsigned long field_size)};
%apply (long* ARGOUT_ARRAY1, unsigned long DIM1) {(long* data_out, unsigned long field_size)};
%apply (unsigned int* ARGOUT_ARRAY1, unsigned long DIM1) {(unsigned int* data_out, unsigned long field_size)};
%apply (int* ARGOUT_ARRAY1, unsigned long DIM1) {(int* data_out, unsigned long field_size)};
%apply (unsigned short* ARGOUT_ARRAY1, unsigned long DIM1) {(unsigned short* data_out, unsigned long field_size)};
%apply (short* ARGOUT_ARRAY1, unsigned long DIM1) {(short* data_out, unsigned long field_size)};
%apply (unsigned char* IN_ARRAY1, unsigned long DIM1) {(unsigned char* data_out, unsigned long field_size)};
%apply (char* ARGOUT_ARRAY1, unsigned long DIM1) {(char* data_out, unsigned long field_size)};

%rename(equals) operator ==;

//%rename(Image2D_double) Image2D<double>;
//%rename(Image2D_float)  Image2D<float>;
//%rename(Image2D_ulong)  Image2D<ulong>;
//%rename(Image2D_long)   Image2D<long>;
//%rename(Image2D_uint)   Image2D<uint>;
//%rename(Image2D_int)    Image2D<int>;
//%rename(Image2D_ushort) Image2D<ushort>;
//%rename(Image2D_short)  Image2D<short>;
//%rename(Image2D_ubyte)  Image2D<uchar>;
//%rename(Image2D_byte)   Image2D<char>;

%include "../../imaging/src/image2D.hpp"
%include "../../imaging/src/image2D_utils.hpp"
%include "../../imaging/src/image.hpp"
#%include "../../common/src/std_typedefs.h"
%import (module="std_typedefs") "../../common/Wrapping/std_typedefs.i"

//%extend Image2D {};
%template(Image2Ddouble) Image2D<double>;
%template(Image2Dfloat)  Image2D<float>;
%template(Image2Dulong)  Image2D<unsigned long>;
%template(Image2Dlong)   Image2D<long>;
%template(Image2Duint)   Image2D<unsigned int>;
%template(Image2Dint)    Image2D<int>;
%template(Image2Dushort) Image2D<unsigned short>;
%template(Image2Dshort)  Image2D<short>;
%template(Image2Dubyte)  Image2D<unsigned char>;
%template(Image2Dbyte)   Image2D<char>;
