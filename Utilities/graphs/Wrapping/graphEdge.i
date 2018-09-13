%module graphEdge

%{
#define SWIG_FILE_WITH_INIT
#include "../../graphs/src/graphEdge.hpp"
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

%include "../../graphs/src/graphEdge.hpp"
%import (module="std_typedefs") "../../common/Wrapping/std_typedefs.i"

