/**
 * Created on 14 Jan 2018
 * author: christian
 */

#ifndef IMAGE2D_HPP_
#define IMAGE2D_HPP_

#include "image.hpp"
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif
#include <iostream>
#include <string.h>

typedef enum _interpolation_mode_2D {
	NEAREST_NEIGHBOUR_2D=0,
	BILINEAR=1
} interpolationMode2D;

template<typename Scalar>
class Image2D : public Image<Scalar, 2> {
public:
	Image2D();
	/**
	 * Image2D shallow-copy constructor
	 */
	Image2D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	/**
	 * Image2D deep-copy constructor
	 */
	Image2D(uint dimX, uint dimY, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	Image2D(const Image2D& X);
	virtual ~Image2D();
	
	inline Image2D<Scalar>& operator=(Image2D<Scalar>& X) { swap(X); return *this; }
	inline Image2D<Scalar>& operator=(const Image2D<Scalar>& X) {
		swap(X);
		return *this;
	}
	bool operator==(Image2D& X);
	bool operator!=(Image2D& X);

	Scalar& operator[](unsigned long index);
	Scalar& operator()(uint coordX, uint coordY);
	inline friend std::ostream& operator<< (std::ostream& os, const Image2D& rhs) {
		if(rhs._dim!=NULL)
			os << "shape = ("<< rhs._dim[0] <<","<< rhs._dim[1] <<"); ";
		if(rhs._data!=NULL) {

		}
		os << "dtype = " << getDTypeString(rhs._datatype);
		os << "dorder = " << getDOrderString(rhs._dataorder);
		return os;
	}

	inline void setDimension(size_t index, uint value) {
		if(index<0)
			Image<Scalar,2>::_dim[0] = value;
		if(index>1)
			return;
		Image<Scalar,2>::_dim[index]=value;
	}
	inline void setDimensions(uint* dim, bool deepCopy=false) {
		if(dim!=NULL) {
			if(deepCopy) {
				if(Image<Scalar,2>::_dim!=NULL)
					delete [] Image<Scalar,2>::_dim;
				Image<Scalar,2>::_dim = new uint[2];
				memcpy(Image<Scalar,2>::_dim, dim, 2*sizeof(uint));
			} else
				Image<Scalar,2>::_dim=dim;
		}
	}
	inline void setDimensions(uint dimX, uint dimY) {
		if(Image<Scalar,2>::_dim!=NULL) {
			delete [] Image<Scalar,2>::_dim;
		}
		Image<Scalar,2>::_dim = new unsigned int[2];
		Image<Scalar,2>::_dim[0]=dimX; Image<Scalar,2>::_dim[1]=dimY;
	}
	inline void setData(Scalar* data) {
		setData(data, false);
	}
	inline void setData(Scalar* data, bool deepCopy) {
		if(data==NULL)
			return;
		clearData();
		if(deepCopy) {
			allocateData();
			copyData(data);
		} else {
			Image<Scalar,2>::_data=data;
		}
	}
	inline void setDatatype(dtype datatype) {
		Image<Scalar,2>::_datatype = datatype;
	}
	inline void setDataorder(dorder dataorder) {
		Image<Scalar,2>::_dataorder = dataorder;
	}
	inline uint* getDimensions() {
		return Image<Scalar,2>::_dim;
	}
	inline uint getDimension(size_t index) {
		if(index<0)
			return Image<Scalar,2>::_dim[0];
		if(index>1)
			return Image<Scalar,2>::_dim[1];
		return Image<Scalar,2>::_dim[index];
	}
	inline virtual Scalar* getData() {
		return Image<Scalar,2>::_data;
	}
	inline virtual dtype getDatatype() {
		return Image<Scalar,2>::_datatype;
	}
	inline virtual dorder getDataorder() {
		return Image<Scalar,2>::_dataorder;
	}

	// WRAPPING funcs - to be 'specialised'
	/*
	inline void SetData(Scalar* data_in, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1])) ? field_size : (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]));
		if(Image<Scalar,2>::_data!=NULL)
			delete [] Image<Scalar,2>::_data;
		Image<Scalar,2>::_data = new Scalar[Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]];
		memcpy(Image<Scalar,2>::_data, data_in, sizeof(Scalar)*fsize);
	}
	*/
	inline void SetData(Scalar* data_in, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1])) ? field_size : (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]));
		if(Image<Scalar,2>::_data!=NULL)
			delete [] Image<Scalar,2>::_data;
		Image<Scalar,2>::_data = new Scalar[Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]];
		memcpy(Image<Scalar,2>::_data, data_in, sizeof(Scalar)*fsize);
		Image<Scalar,2>::_dataorder = ROW_MAJOR;
	}
	inline void GetData(Scalar* data_out, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1])) ? field_size : (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]));
		memcpy(data_out, Image<Scalar,2>::_data, sizeof(Scalar)*fsize);
	}
	inline Scalar GetImageValue(uint coordX, uint coordY) {
		if(Image<Scalar,2>::_dim==NULL)
			return 0;
		if(Image<Scalar,2>::_data==NULL)
			return 0;
		unsigned long index = 0;
		if(Image<Scalar,2>::_dataorder == ROW_MAJOR) {
			index = getDataAddress_rowMajor(coordX, coordY);
		} else if(Image<Scalar,2>::_dataorder == COLUMN_MAJOR) {
			index = getDataAddress_columnMajor(coordX, coordY);
		} else {
			return 0;
		}
		if(index >= (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]))
			return 0;
		return Image<Scalar,2>::_data[index];
	}
	// WRAPPING funcs - to be 'specialised'

	Scalar interpolate(float coordX, float coordY, interpolationMode2D mode);
	Scalar& getImageValue(uint coordX, uint coordY);

	void setImageValue(Scalar value, uint coordX, uint coordY);
	void createImage(void);
	void createImage(uint* dim, dtype datatype);
	void clean(void);
	void pad(uint totalPadX, uint totalPadY);

	void clamp(const Scalar& lowValuePtr, const Scalar& highValuePtr);
	Scalar mean(void);
	Scalar variance(void);
	Scalar stdDev(void);
	Scalar max(void);
	Scalar min(void);


	inline void SetDimensions(unsigned int* dimensions_in, unsigned int field_size) {
		size_t fsize = (field_size<=2 ? field_size : 2);
		if(Image<Scalar,2>::_dim!=NULL)
			delete [] Image<Scalar,2>::_dim;
		Image<Scalar,2>::_dim = new unsigned int[2];
		memcpy(Image<Scalar,2>::_dim, dimensions_in, sizeof(unsigned int)*fsize);
	}
	inline void GetDimensions(unsigned int* dimensions_out, unsigned int field_size) {
		size_t fsize = (field_size<=2 ? field_size : 2);
		memcpy(dimensions_out, Image<Scalar,2>::_dim, sizeof(unsigned int)*fsize);
	}

	inline unsigned long getDataAddress(uint x, uint y) {
		if(Image<Scalar,2>::_dataorder==ROW_MAJOR)
			return getDataAddress_rowMajor(x,y);
		else if(Image<Scalar,2>::_dataorder==COLUMN_MAJOR)
			return getDataAddress_columnMajor(x,y);
		return 0;
	}

protected:
	typedef Image<Scalar,2> Base;
	void swap(Image2D<Scalar>& X);
	void swap(const Image2D<Scalar>& X);
	Scalar interpolateBilinear(float coordX, float coordY);
	Scalar interpolateNearestNeighbour(float coordX, float coordY);

	inline void clearData() {
		if(Image<Scalar,2>::_data!=NULL) {
			delete [] Image<Scalar,2>::_data;
		}
		Image<Scalar,2>::_data=NULL;
	}
	inline void allocateData() {
		unsigned long fsize=1;
		if(Image<Scalar,2>::_dim==NULL)
			return;
		for(size_t index=0; index<2; index++)
			fsize*=Image<Scalar,2>::_dim[index];
		if(Image<Scalar,2>::_data==NULL) {
			Image<Scalar,2>::_data=new Scalar[fsize];
			//bzero(Image<Scalar,2>::_data, sizeof(Scalar)*fsize);
			memset(Image<Scalar,2>::_data, 0, sizeof(Scalar)*fsize);
		}
	}
	inline void zeroData() {
		if(Image<Scalar,2>::_data!=NULL) {
			unsigned long fsize=1;
			if(Image<Scalar,2>::_dim==NULL)
				return;
			for(size_t index=0; index<2; index++)
				fsize*=Image<Scalar,2>::_dim[index];
			//bzero(Image<Scalar,2>::_data, sizeof(Scalar));
			memset(Image<Scalar,2>::_data, 0, sizeof(Scalar)*fsize);
		}
	}
	inline void copyData(Scalar* source) {
		unsigned long fsize=1;
		if(Image<Scalar,2>::_dim==NULL)
			return;
		for(size_t index=0; index<2; index++)
			fsize*=Image<Scalar,2>::_dim[index];
		if((Image<Scalar,2>::_data!=NULL) && (source!=NULL)) {
			memcpy(Image<Scalar,2>::_data, source, fsize*sizeof(Scalar));
		}
	}
	inline unsigned long getDataAddress_columnMajor(uint x, uint y) {
		//return z*_dim[1]*_dim[0]+x*_dim[1]+y;
		return x+Image<Scalar,2>::_dim[0]*(y+Image<Scalar,2>::_dim[1]*(0));
	}
	inline unsigned long getDataAddress_rowMajor(uint x, uint y) {
		//return y*_dim[0]*dim_[2]+x*_dim[2]+y;
		return y+Image<Scalar,2>::_dim[1]*(x+Image<Scalar,2>::_dim[0]*(0));
	}
	// ??? ??? recalculate before change ??? ???  judegment: need to switch //
	inline void getIndicesFromAddress_columnMajor(unsigned long address, uint& x, uint& y) {
		y = address / Base::_dim[0];
		x = address % Base::_dim[0];
	}
	inline void getIndicesFromAddress_rowMajor(unsigned long address, uint& x, uint& y) {
		x = address / Base::_dim[1];
		y = address % Base::_dim[1];
	}

	//uint* _dim;
	//Scalar* _data;
	//dtype _datatype;
	//dorder _dataorder;
};

inline unsigned long getDataAddress_columnMajor_field2D(uint x, uint dimX, uint y, uint dimY) {
	//return z*_dim[1]*_dim[0]+x*_dim[1]+y;
	return x+dimX*(y+dimY*(0));
}
inline unsigned long getDataAddress_rowMajor_field2D(uint x, uint dimX, uint y, uint dimY) {
	//return y*_dim[0]*dim_[2]+x*_dim[2]+y;
	return y+dimY*(x+dimX*(0));
}
inline unsigned long getDataAddress_field2D(uint x, uint dimX, uint y, uint dimY, dorder dataorder) {
	if(dataorder==ROW_MAJOR)
		return getDataAddress_rowMajor_field2D(x,dimX,y,dimY);
	else if(dataorder==COLUMN_MAJOR)
		return getDataAddress_columnMajor_field2D(x,dimX,y,dimY);
	return 0;
}

// ??? ??? recalculate before change ??? ???  judegment: need to switch //
inline void getIndicesFromAddress_columnMajor_field2D(unsigned long address, uint& x, uint dimX, uint& y, uint dimY) {
	y = address / dimX;
	x = address % dimX;
}
inline void getIndicesFromAddress_rowMajor_field2D(unsigned long address, uint& x, uint dimX, uint& y, uint dimY) {
	x = address / dimY;
	y = address % dimY;
}
inline void getIndicesFromAddress_field2D(unsigned long address,uint &x, uint dimX, uint &y, uint dimY, dorder dataorder) {
	if(dataorder==ROW_MAJOR)
		getIndicesFromAddress_rowMajor_field2D(address,x,dimX,y,dimY);
	else if(dataorder==COLUMN_MAJOR)
		getIndicesFromAddress_columnMajor_field2D(address,x,dimX,y,dimY);
}


/*
 *
 	 switch(_datatype) {
			case CHAR: {
				break;
			}
			case UCHAR: {
				break;
			}
			case SHORT: {
				break;
			}
			case USHORT: {
				break;
			}
			case INT: {
				break;
			}
			case UINT: {
				break;
			}
			case LONG: {
				break;
			}
			case ULONG: {
				break;
			}
			case FLOAT: {
				break;
			}
			case DOUBLE: {
				break;
			}
			default: {
				break;
			}
		}
 */

/*
template<>
inline void Image2D<char>::SetData(char* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new char[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(char)*fsize);
}
template<>
inline void Image2D<char>::GetData(char* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(char)*fsize);
}
template<>
inline char Image2D<char>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<char>;



template<>
inline void Image2D<uchar>::SetData(uchar* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new uchar[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(uchar)*fsize);
}
template<>
inline void Image2D<uchar>::GetData(uchar* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(uchar)*fsize);
}
template<>
inline uchar Image2D<uchar>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<uchar>;



template<>
inline void Image2D<short>::SetData(short* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new short[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(short)*fsize);
}
template<>
inline void Image2D<short>::GetData(short* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(short)*fsize);
}
template<>
inline short Image2D<short>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<short>;



template<>
inline void Image2D<ushort>::SetData(ushort* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new ushort[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(ushort)*fsize);
}
template<>
inline void Image2D<ushort>::GetData(ushort* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(ushort)*fsize);
}
template<>
inline ushort Image2D<ushort>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<ushort>;



template<>
inline void Image2D<int>::SetData(int* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new int[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(int)*fsize);
}
template<>
inline void Image2D<int>::GetData(int* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(int)*fsize);
}
template<>
inline int Image2D<int>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<int>;



template<>
inline void Image2D<uint>::SetData(uint* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new uint[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(uint)*fsize);
}
template<>
inline void Image2D<uint>::GetData(uint* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(uint)*fsize);
}
template<>
inline uint Image2D<uint>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<uint>;



template<>
inline void Image2D<long>::SetData(long* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new long[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(long)*fsize);
}
template<>
inline void Image2D<long>::GetData(long* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(long)*fsize);
}
template<>
inline long Image2D<long>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<long>;



template<>
inline void Image2D<ulong>::SetData(ulong* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new ulong[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(ulong)*fsize);
}
template<>
inline void Image2D<ulong>::GetData(ulong* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(ulong)*fsize);
}
template<>
inline ulong Image2D<ulong>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<ulong>;



template<>
inline void Image2D<float>::SetData(float* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new float[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(float)*fsize);
}
template<>
inline void Image2D<float>::GetData(float* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(float)*fsize);
}
template<>
inline float Image2D<float>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<float>;



template<>
inline void Image2D<double>::SetData(double* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	if(_data!=NULL)
		delete [] _data;
	_data = new double[_dim[0]*_dim[1]];
	memcpy(_data, data_in, sizeof(double)*fsize);
}
template<>
inline void Image2D<double>::GetData(double* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1])) ? field_size : (_dim[0]*_dim[1]));
	memcpy(data_out, _data, sizeof(double)*fsize);
}
template<>
inline double Image2D<double>::GetImageValue(uint coordX, uint coordY) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]))
		return 0;
	return _data[index];
}
template class Image2D<double>;
*/

/*
template class Image2D<char>;
template class Image2D<uchar>;
template class Image2D<short>;
template class Image2D<ushort>;
template class Image2D<int>;
template class Image2D<uint>;
template class Image2D<long>;
template class Image2D<ulong>;
template class Image2D<float>;
template class Image2D<double>;
*/

#include "image2D_utils.hpp"


#endif /* IMAGE2D_HPP_ */
