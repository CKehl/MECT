%module volumeFileInterface

%{
#define SWIG_FILE_WITH_INIT
%}

%feature("autodoc", "1");
%include <std_iostream.i>
%include <std_sstream.i>
%include <std_string.i>
%import (module="ini") "INI.i"
%import (module="mhd") "MHD.i"
%import (module="volumeSegmentationDat") "volumeSegmentationDAT.i"
%import (module="fams_ascii") "FAMS_ASCII.i"
