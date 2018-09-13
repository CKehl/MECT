%module MhdImage

%{
#define SWIG_FILE_WITH_INIT
#include "../../imaging/src/image.hpp"
#include "../../imaging/src/image_utils.hpp"
#include "../../imaging/src/image2D.hpp"
#include "../../imaging/src/image2D_utils.hpp"
#include "../../imaging/src/image3D.hpp"
#include "../../imaging/src/image3D_utils.hpp"
#include "../../imaging/src/image4D.hpp"
#include "../../imaging/src/image4D_utils.hpp"
#include "../../imaging/src/mhdImage.hpp"
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
%numpy_typemaps(size_t       , NPY_ULONG  , unsigned int)
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
%numpy_typemaps(size_t       , NPY_UINT   , long)
%numpy_typemaps(int          , NPY_INT    , long)
%numpy_typemaps(unsigned short ,NPY_USHORT, long)
%numpy_typemaps(short        , NPY_SHORT  , long)

// apply numpy typemaps based on arg types + names
%apply (unsigned int* IN_ARRAY1, unsigned int DIM1) {(unsigned int* dimensions_in, unsigned int field_size)};
%apply (double* IN_ARRAY1, unsigned long DIM1) {(double* image_data_in, unsigned long field_size)};
%apply (float* IN_ARRAY1, unsigned long DIM1) {(float* image_data_in, unsigned long field_size)};
%apply (unsigned long* IN_ARRAY1, unsigned long DIM1) {(unsigned long* image_data_in, unsigned long field_size)};
%apply (long* IN_ARRAY1, unsigned long DIM1) {(long* image_data_in, unsigned long field_size)};
%apply (unsigned int* IN_ARRAY1, unsigned long DIM1) {(unsigned int* image_data_in, unsigned long field_size)};
%apply (int* IN_ARRAY1, unsigned long DIM1) {(int* image_data_in, unsigned long field_size)};
%apply (unsigned short* IN_ARRAY1, unsigned long DIM1) {(unsigned short* image_data_in, unsigned long field_size)};
%apply (short* IN_ARRAY1, unsigned long DIM1) {(short* image_data_in, unsigned long field_size)};
%apply (unsigned char* IN_ARRAY1, unsigned long DIM1) {(unsigned char* image_data_in, unsigned long field_size)};
%apply (char* IN_ARRAY1, unsigned long DIM1) {(char* image_data_in, unsigned long field_size)};

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

// ---- special addition needed for the FORTRAN-aligned, column-major order arrays we use here ... these only apply to the input types---- //
//%apply (unsigned int* IN_FARRAY1, int DIM1) {(unsigned int* dimensions_in, unsigned int field_size)};
//%apply (double* IN_FARRAY1, long DIM1) {(double* image_data_in_F, unsigned long field_size)};
//%apply (float* IN_FARRAY1, unsigned long DIM1) {(float* image_data_in_F, unsigned long _field_size)};
//%apply (unsigned long* IN_FARRAY1, unsigned long DIM1) {(unsigned long* image_data_in_F, unsigned long field_size)};
//%apply (long* IN_FARRAY1, unsigned long DIM1) {(long* image_data_in_F, unsigned long field_size)};
//%apply (unsigned int* IN_FARRAY1, unsigned long DIM1) {(unsigned int* image_data_in_F, unsigned long field_size)};
//%apply (int* IN_FARRAY1, unsigned long DIM1) {(int* image_data_in_F , unsigned long field_size)};
//%apply (unsigned short* IN_FARRAY1, unsigned long DIM1) {(unsigned short* image_data_in_F, unsigned long field_size)};
//%apply (short* IN_FARRAY1, unsigned long DIM1) {(short* image_data_in_F, unsigned long field_size)};
//%apply (unsigned char* IN_FARRAY1, unsigned long DIM1) {(unsigned char* image_data_in_F, unsigned long field_size)};
//%apply (char* IN_FARRAY1, unsigned long DIM1) {(char* image_data_in_F, unsigned long field_size)};

%rename(equals) operator ==;

%include "../../imaging/src/mhdImage.hpp"
%import (module="Image2D") "../../imaging/src/image2D.hpp"
%import (module="Image3D") "../../imaging/src/image3D.hpp"
%import (module="Image4D") "../../imaging/src/image4D.hpp"
%import (module="std_typedefs") "../../common/Wrapping/std_typedefs.i"

//%extend Image4D {};
%template(MhdImage2Ddouble) MhdImage2D<double>;
%template(MhdImage2Dfloat)  MhdImage2D<float>;
%template(MhdImage2Dulong)  MhdImage2D<unsigned long>;
%template(MhdImage2Dlong)   MhdImage2D<long>;
%template(MhdImage2Duint)   MhdImage2D<unsigned int>;
%template(MhdImage2Dint)    MhdImage2D<int>;
%template(MhdImage2Dushort) MhdImage2D<unsigned short>;
%template(MhdImage2Dshort)  MhdImage2D<short>;
%template(MhdImage2Dubyte)  MhdImage2D<unsigned char>;
%template(MhdImage2Dbyte)   MhdImage2D<char>;

%template(MhdImage3Ddouble) MhdImage2D<double>;
%template(MhdImage3Dfloat)  MhdImage3D<float>;
%template(MhdImage3Dulong)  MhdImage3D<unsigned long>;
%template(MhdImage3Dlong)   MhdImage3D<long>;
%template(MhdImage3Duint)   MhdImage3D<unsigned int>;
%template(MhdImage3Dint)    MhdImage3D<int>;
%template(MhdImage3Dushort) MhdImage3D<unsigned short>;
%template(MhdImage3Dshort)  MhdImage3D<short>;
%template(MhdImage3Dubyte)  MhdImage3D<unsigned char>;
%template(MhdImage3Dbyte)   MhdImage3D<char>;

%template(MhdImage4Ddouble) MhdImage4D<double>;
%template(MhdImage4Dfloat)  MhdImage4D<float>;
%template(MhdImage4Dulong)  MhdImage4D<unsigned long>;
%template(MhdImage4Dlong)   MhdImage4D<long>;
%template(MhdImage4Duint)   MhdImage4D<unsigned int>;
%template(MhdImage4Dint)    MhdImage4D<int>;
%template(MhdImage4Dushort) MhdImage4D<unsigned short>;
%template(MhdImage4Dshort)  MhdImage4D<short>;
%template(MhdImage4Dubyte)  MhdImage4D<unsigned char>;
%template(MhdImage4Dbyte)   MhdImage4D<char>;