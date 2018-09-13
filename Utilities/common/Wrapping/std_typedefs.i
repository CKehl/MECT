%module std_typedefs

%{
#define SWIG_FILE_WITH_INIT
#include "../../common/src/std_typedefs.h"
%}

%feature("autodoc", "1");

%include <std_iostream.i>
%include <std_sstream.i>
%include <std_vector.i>
%include <typemaps.i>

%rename(equals) operator ==;

%include "../../common/src/std_typedefs.h"
