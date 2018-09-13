/**
 * created on Mar 1 2018
 * author: Christian Kehl
 */
#ifndef FAMS_ASCII_HPP_
#define FAMS_ASCII_HPP_

#include <stdlib.h>
#include <cstring>
#include <iostream>
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif

class FamsFile_ASCII {
public:
	inline FamsFile_ASCII() : _data(NULL), _ndims(0), _dim(NULL), _filename(""), _datatype(USHORT), _dataorder(ROW_MAJOR) {}
	inline virtual ~FamsFile_ASCII() {
		if(_data!=NULL) {
			delete [] ((ushort*)_data);
		}
		if(_dim!=NULL)
			delete [] _dim;
	}
	
	void setFilename(std::string filename);
	inline void setDatatype(dtype datatype) { _datatype = datatype; }
	inline void setDataorder(dorder dataorder) { _dataorder = dataorder; }
	inline void setNumberOfDimensions(int ndims) { _ndims = ndims; }
	void setDimensions(uint* dimensions);
	void SetDimensions(unsigned int* numpy_input_dimensions, int field_size);
	void setData(void* data, dtype datatype, bool memcopy=true); // needs to be "something" unitary originally
	void setDataByUShort(unsigned short* numpy_input_data, long field_size);
	void setDataByUInt(unsigned int* numpy_input_data, long field_size);

	inline void* getData() { return _data; }
	void getDataAsUShort(unsigned short* numpy_output_data, long field_size);
	void getDataAsUInt(unsigned int* numpy_output_data, long field_size);
	inline uint* getDimensions() { return _dim; }
	void getDimensions(unsigned int* numpy_output_dimensions, int field_size);
	inline dtype getDatatype() { return _datatype; }
	inline dorder getDataorder() { return _dataorder; }
	inline int getNumberOfDimensions() { return _ndims; }
	inline std::string getFileName() { return _filename; }

	void readFile(void);
	void writeFile(void);

protected:
	inline void clearData(void) {
		if(_data==NULL)
			return;
		else {
			delete [] ((ushort*)_data);
		}
	}
	void* _data;
	int _ndims;				// auto-deducted when setting the dimensions
	unsigned int * _dim;
	std::string _filename;
	dtype _datatype;
	dorder _dataorder;
};



#endif /* FAMS_ASCII_HPP_ */
