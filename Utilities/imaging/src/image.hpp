/**
 * Created on 14 Jan 2018
 * author: christian
 */

#ifndef IMAGE_HPP_
#define IMAGE_HPP_

#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif

template<typename Scalar, uint ndims>
class Image {
public:
	inline Image() : _dim(NULL), _data(NULL), _datatype(NONE), _dataorder(UNDEF_ORDER) {

	}
	inline virtual ~Image() {

	}

	inline void setDimension(size_t index, uint value) {
		if(index<0)
			_dim[0] = value;
		if(index>(ndims-1))
			return;
		_dim[index]=value;
	}
	inline void setDimensions(uint* dim) {
		if(dim!=NULL)
			_dim=dim;
	}
	//template<uint ndims>
	inline void setDimensions(uint dimX, uint dimY, uint dimZ, uint dimW) {
		if(_dim!=NULL) {
			delete [] _dim;
		}
		_dim = new unsigned int[ndims];
		_dim[0] = dimX; _dim[1] = dimY;
		if(ndims>2)
			_dim[2]=dimZ;
		if(ndims>3)
			_dim[3]=dimW;
	}
	inline void setDatatype(dtype datatype) {
		_datatype = datatype;
	}
	inline void setDataorder(dorder dataorder) {
		_dataorder = dataorder;
	}
	inline uint* getDimensions() {
		return _dim;
	}
	//template<uint ndims>
	inline uint getDimension(size_t index) {
		if(index<0)
			return _dim[0];
		if(index>(ndims-1))
			return 0;
		return _dim[index];
	}
	inline virtual Scalar* getData() {
		return _data;
	}
	inline virtual dtype getDatatype() {
		return _datatype;
	}
	inline virtual dorder getDataorder() {
		return _dataorder;
	}
	//template<uint ndims>
	inline unsigned long getDataAddress(uint x, uint y, uint z, uint w) {
		if(getDataorder()==ROW_MAJOR)
			return getDataAddress_rowMajor(x,y,z,w);
		else if(getDataorder()==COLUMN_MAJOR)
			return getDataAddress_columnMajor(x,y,z,w);
		return 0;
	}
protected:
	//template<uint ndims>
	inline unsigned long getDataAddress_columnMajor(uint x, uint y, uint z, uint w) {
		uint dim3=0, dim4=0;
		if(ndims>2)
			dim3 = _dim[2];
		else
			z = 0;
		if(ndims>3)
			dim4 = _dim[3];
		else
			w = 0;
		return x+_dim[0]*(y+_dim[1]*(z+dim3*(w+dim4*(0))));
	}
	//template<uint ndims>
	inline unsigned long getDataAddress_rowMajor(uint x, uint y, uint z, uint w) {
		uint dim3=0, dim4=0;
		if(ndims>2)
			dim3 = _dim[2];
		else {
			z = 0;
			return y+_dim[1]*(x+_dim[0]*(0));
		}
		if(ndims>3)
			dim4 = _dim[3];
		else {
			w = 0;
			return z+dim3*(y+_dim[1]*(x+_dim[0]*(0)));
		}
		return w+dim4*(z+dim3*(y+_dim[1]*(x+_dim[0]*(0))));
	}

	inline void getIndicesFromAddress_columnMajor(unsigned long address, uint& x, uint& y, uint& z, uint& w) {
		uint dim3=0, dim4=0;
		if(ndims>2)
			dim3 = _dim[2];
		else {
			z = 0;
			y = address / _dim[0];
			x = address % _dim[0];
			return;
		}
		if(ndims>3)
			dim4 = _dim[3];
		else {
			w = 0;
			z = address / (_dim[0]*_dim[1]);
			unsigned long rest = address % (_dim[0]*_dim[1]);
			y = rest / _dim[0];
			x = rest % _dim[0];
			return;
		}
		w = address / (_dim[0]*_dim[1]*dim3);
		unsigned long rest = address % (_dim[0]*_dim[1]*dim3);
		z = rest / (_dim[0]*_dim[1]);
		rest = rest % (_dim[0]*_dim[1]);
		y = rest / _dim[0];
		x = rest % _dim[0];
		return;
	}
	inline void getIndicesFromAddress_rowMajor(unsigned long address, uint& x, uint& y, uint& z, uint& w) {
		uint dim3=0, dim4=0;
		if(ndims>2)
			dim3 = _dim[2];
		else {
			z = 0;
			x = address / _dim[1];
			y = address % _dim[1];
			return;
		}
		if(ndims>3)
			dim4 = _dim[3];
		else {
			w = 0;
			x = address / (dim3*_dim[1]);
			unsigned long rest = address % (dim3*_dim[1]);
			y = rest / dim3;
			z = rest % dim3;
			return;
		}
		x = address / (dim4*dim3*_dim[1]);
		unsigned long rest = address % (dim4*dim3*_dim[1]);
		y = rest / (dim4*dim3);
		rest = rest % (dim4*dim3);
		z = rest / dim4;
		w = rest % dim4;
		return;
	}

	uint* _dim;
	Scalar* _data;
	dtype _datatype;
	dorder _dataorder;
};

#include "image_utils.hpp"



#endif /* IMAGE_HPP_ */
