%module GraphCut3D

%{
#define SWIG_FILE_WITH_INIT
#include "../../graphCut3D/src/GraphCut3D.hpp"
#include "../../../Utilities/graphs/src/graphDisjointForest.hpp"
#include "../../../Utilities/graphs/src/graphEdge.hpp"
#include "../../../Utilities/imaging/src/image2D.hpp"
#include "../../../Utilities/imaging/src/image2D_utils.hpp"
#include "../../../Utilities/imaging/src/image3D.hpp"
#include "../../../Utilities/imaging/src/image3D_utils.hpp"
#include "../../../Utilities/imaging/src/image.hpp"
#include "../../../Utilities/common/src/std_typedefs.h"
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

%include "../../graphCut3D/src/GraphCut3D.hpp"
%import (module="graphDisjointForest") "../../../Utilities/graphs/src/graphDisjointForest.hpp"
%import (module="graphEdge") "../../../Utilities/graphs/src/graphEdge.hpp"
%import (module="Image3D") "../../../Utilities/imaging/src/image2D.hpp"
%import (module="Image4D") "../../../Utilities/imaging/src/image3D.hpp"
%import (module="std_typedefs") "../../../Utilities/common/Wrapping/std_typedefs.i"

%template(GraphCut3Ddouble) GraphCut3D<double>;
%template(GraphCut3Dfloat)  GraphCut3D<float>;
%template(GraphCut3Dulong)  GraphCut3D<unsigned long>;
%template(GraphCut3Dlong)   GraphCut3D<long>;
%template(GraphCut3Duint)   GraphCut3D<unsigned int>;
%template(GraphCut3Dint)    GraphCut3D<int>;
%template(GraphCut3Dushort) GraphCut3D<unsigned short>;
%template(GraphCut3Dshort)  GraphCut3D<short>;
%template(GraphCut3Dubyte)  GraphCut3D<unsigned char>;
%template(GraphCut3Dbyte)   GraphCut3D<char>;
