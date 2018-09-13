%module GraphCut4D

%{
#define SWIG_FILE_WITH_INIT
#include "../../graphCut4D/src/GraphCut4D.hpp"
#include "../../../Utilities/graphs/src/graphDisjointForest.hpp"
#include "../../../Utilities/graphs/src/graphEdge.hpp"
#include "../../../Utilities/imaging/src/image3D.hpp"
#include "../../../Utilities/imaging/src/image3D_utils.hpp"
#include "../../../Utilities/imaging/src/image4D.hpp"
#include "../../../Utilities/imaging/src/image4D_utils.hpp"
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

%include "../../graphCut4D/src/GraphCut4D.hpp"
%import (module="graphDisjointForest") "../../../Utilities/graphs/src/graphDisjointForest.hpp"
%import (module="graphEdge") "../../../Utilities/graphs/src/graphEdge.hpp"
%import (module="Image3D") "../../../Utilities/imaging/src/image3D.hpp"
%import (module="Image4D") "../../../Utilities/imaging/src/image4D.hpp"
%import (module="std_typedefs") "../../../Utilities/common/Wrapping/std_typedefs.i"

//%extend Image2D {};
%template(GraphCut4Ddouble) GraphCut4D<double>;
%template(GraphCut4Dfloat)  GraphCut4D<float>;
%template(GraphCut4Dulong)  GraphCut4D<unsigned long>;
%template(GraphCut4Dlong)   GraphCut4D<long>;
%template(GraphCut4Duint)   GraphCut4D<unsigned int>;
%template(GraphCut4Dint)    GraphCut4D<int>;
%template(GraphCut4Dushort) GraphCut4D<unsigned short>;
%template(GraphCut4Dshort)  GraphCut4D<short>;
%template(GraphCut4Dubyte)  GraphCut4D<unsigned char>;
%template(GraphCut4Dbyte)   GraphCut4D<char>;
