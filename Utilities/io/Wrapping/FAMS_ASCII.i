%module fams_ascii

%{
#define SWIG_FILE_WITH_INIT
#include "../src/FAMS_ASCII.hpp"
#include "../../common/src/std_typedefs.h"
%}

%feature("autodoc", "1");

%include <std_iostream.i>
%include <std_sstream.i>
%include <std_string.i>
%include <typemaps.i>
%include "numpy.i"
%init %{
	import_array();
%}


// here comes the "meat" of the wrapping ...


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
%numpy_typemaps(long         , NPY_LONG    , unsigned int)
%numpy_typemaps(unsigned int , NPY_UINT   , unsigned int)
%numpy_typemaps(size_t       , NPY_UINT   , unsigned int)
%numpy_typemaps(int          , NPY_INT    , unsigned int)
%numpy_typemaps(short        , NPY_SHORT  , unsigned int)
%numpy_typemaps(double       , NPY_DOUBLE , int)
%numpy_typemaps(float        , NPY_FLOAT  , int)
%numpy_typemaps(long         , NPY_LONG    , int)
%numpy_typemaps(unsigned int , NPY_UINT   , int)
%numpy_typemaps(size_t       , NPY_UINT   , int)
%numpy_typemaps(int          , NPY_INT    , int)
%numpy_typemaps(short        , NPY_SHORT  , int)
%numpy_typemaps(double       , NPY_DOUBLE , long)
%numpy_typemaps(float        , NPY_FLOAT  , long)
%numpy_typemaps(unsigned long, NPY_ULONG  , long)
%numpy_typemaps(long         , NPY_LONG   , long)
%numpy_typemaps(unsigned int , NPY_UINT   , long)
%numpy_typemaps(size_t       , NPY_ULONG  , long)
%numpy_typemaps(int          , NPY_INT    , long)
%numpy_typemaps(unsigned short, NPY_USHORT, long)
%numpy_typemaps(short        , NPY_SHORT  , long)
%numpy_typemaps(unsigned char, NPY_UBYTE  , long)
%numpy_typemaps(char         , NPY_BYTE   , long)


// apply numpy typemaps based on arg types + names
%apply (unsigned int* IN_ARRAY1, int DIM1) {(unsigned int* numpy_input_dimensions, int field_size)};
%apply (unsigned short* IN_ARRAY1, long DIM1) {(unsigned short* numpy_input_data, long field_size)};
%apply (unsigned int* IN_ARRAY1, long DIM1) {(unsigned int* numpy_input_data, long field_size)};

%apply (unsigned int* ARGOUT_ARRAY1, int DIM1) {(unsigned int* numpy_output_dimensions, int field_size)};
%apply (unsigned short* ARGOUT_ARRAY1, long DIM1) {(unsigned short* numpy_output_data, long field_size)};
%apply (unsigned int* ARGOUT_ARRAY1, long DIM1) {(unsigned int* numpy_output_data, long field_size)};

// ---- special addition needed for the FORTRAN-aligned, column-major order arrays we use here ... these only apply to the input types---- //
%apply (unsigned int* IN_FARRAY1, int DIM1) {(unsigned int* numpy_input_dimensions, int field_size)};
%apply (unsigned short* IN_FARRAY1, long DIM1) {(unsigned short* numpy_input_data, long field_size)}
%apply (unsigned int* IN_FARRAY1, long DIM1) {(unsigned int* numpy_input_data, long field_size)}

%include "../../common/src/std_typedefs.h"
%include "../src/FAMS_ASCII.hpp"
