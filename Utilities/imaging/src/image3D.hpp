/**
 * Created on 18 Jan 2018
 * author: christian
 */

#ifndef IMAGE3D_HPP_
#define IMAGE3D_HPP_

#include "image.hpp"
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif
#include <iostream>
#include <string.h>
#include <cmath>
#include <omp.h>

typedef enum _interpolation_mode_3D {
	NEAREST_NEIGHBOUR_3D=0,
	TRILINEAR=1
} interpolationMode3D;

template<typename Scalar>
class Image3D : public Image<Scalar,3> {
public:
	Image3D();
	/**
	 * Image3D shallow-copy constructor
	 */
	Image3D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	/**
	 * Image3D deep-copy constructor
	 */
	Image3D(uint dimX, uint dimY, uint dimZ, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	Image3D(const Image3D& X);
	virtual ~Image3D();
	
	inline Image3D<Scalar>& operator=(Image3D<Scalar>& X) { swap(X); return *this; }
	inline Image3D<Scalar>& operator=(const Image3D<Scalar>& X) {
		swap(X);
		return *this;
	}
	bool operator==(Image3D& X);
	bool operator!=(Image3D& X);

	Scalar& operator[](unsigned long index);
	Scalar& operator()(uint coordX, uint coordY, uint coordZ);
	inline friend std::ostream& operator<< (std::ostream& os, const Image3D& rhs) {
		if(rhs._dim!=NULL)
			os << "shape = ("<< rhs._dim[0] <<","<< rhs._dim[1] <<","<< rhs._dim[2] <<"); ";
		if(rhs._data!=NULL) {

		}
		os << "dtype = " << getDTypeString(rhs._datatype);
		os << "dorder = " << getDOrderString(rhs._dataorder);
		return os;
	}

	inline void setDimension(size_t index, uint value) {
		if(index<0)
			Image<Scalar,3>::_dim[0] = value;
		if(index>2)
			return;
		Image<Scalar,3>::_dim[index]=value;
	}
	inline void setDimensions(uint* dim, bool deepCopy=false) {
		if(dim!=NULL) {
			if(deepCopy) {
				if(Image<Scalar,3>::_dim!=NULL)
					delete [] Image<Scalar,3>::_dim;
				Image<Scalar,3>::_dim = new uint[3];
				memcpy(Image<Scalar,3>::_dim, dim, 3*sizeof(uint));
			} else
				Image<Scalar,3>::_dim=dim;
		}
	}
	inline void setDimensions(uint dimX, uint dimY, uint dimZ) {
		if(Image<Scalar,3>::_dim!=NULL) {
			delete [] Image<Scalar,3>::_dim;
		}
		Image<Scalar,3>::_dim = new unsigned int[3];
		Image<Scalar,3>::_dim[0]=dimX; Image<Scalar,3>::_dim[1]=dimY; Image<Scalar,3>::_dim[2]=dimZ;
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
			Image<Scalar,3>::_data=data;
		}
	}
	inline void setDatatype(dtype datatype) {
		Image<Scalar,3>::_datatype = datatype;
	}
	inline void setDataorder(dorder dataorder) {
		Image<Scalar,3>::_dataorder = dataorder;
	}
	inline uint* getDimensions() {
		return Image<Scalar,3>::_dim;
	}
	inline uint getDimension(size_t index) {
		if(index<0)
			return Image<Scalar,3>::_dim[0];
		if(index>2)
			return Image<Scalar,3>::_dim[2];
		return Image<Scalar,3>::_dim[index];
	}
	inline virtual Scalar* getData() {
		return Image<Scalar,3>::_data;
	}
	inline virtual dtype getDatatype() {
		return Image<Scalar,3>::_datatype;
	}
	inline virtual dorder getDataorder() {
		return Image<Scalar,3>::_dataorder;
	}

	// WRAPPING funcs - to be 'specialised'
	/*
	inline void SetData(Scalar* data_in, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2])) ? field_size : (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]));
		if(Image<Scalar,3>::_data!=NULL)
			delete [] Image<Scalar,3>::_data;
		Image<Scalar,3>::_data = new Scalar[Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]];
		memcpy(Image<Scalar,3>::_data, data_in, sizeof(Scalar)*fsize);
	}
	*/
	inline void SetData(Scalar* data_in, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2])) ? field_size : (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]));
		if(Image<Scalar,3>::_data!=NULL)
			delete [] Image<Scalar,3>::_data;
		Image<Scalar,3>::_data = new Scalar[Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]];
		memcpy(Image<Scalar,3>::_data, data_in, sizeof(Scalar)*fsize);
		Image<Scalar,3>::_dataorder = ROW_MAJOR;
	}
	inline void GetData(Scalar* data_out, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2])) ? field_size : (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]));
		memcpy(data_out, Image<Scalar,3>::_data, sizeof(Scalar)*fsize);
	}
	inline Scalar GetImageValue(uint coordX, uint coordY, uint coordZ) {
		if(Image<Scalar,3>::_dim==NULL)
			return 0;
		if(Image<Scalar,3>::_data==NULL)
			return 0;
		unsigned long index = 0;
		if(Image<Scalar,3>::_dataorder == ROW_MAJOR) {
			index = getDataAddress_rowMajor(coordX, coordY, coordZ);
		} else if(Image<Scalar,3>::_dataorder == COLUMN_MAJOR) {
			index = getDataAddress_columnMajor(coordX, coordY, coordZ);
		} else {
			return 0;
		}
		if(index >= (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]))
			return 0;
		return Image<Scalar,3>::_data[index];
	}
	// WRAPPING funcs - to be 'specialised'

	Scalar interpolate(float coordX, float coordY, float coordZ, interpolationMode3D mode);
	Scalar& getImageValue(uint coordX, uint coordY, uint coordZ);

	void setImageValue(Scalar value, uint coordX, uint coordY, uint coordZ);
	void createImage(void);
	void createImage(uint* dim, dtype datatype);
	void clean(void);
	void pad(uint totalPadX, uint totalPadY, uint totalPadZ);
	void clip(size_t origin0, size_t origin1, size_t origin2, size_t dim0, size_t dim1, size_t dim2);
	inline void clipByBounds(size_t low_dim0, size_t high_dim0, size_t low_dim1, size_t high_dim1, size_t low_dim2, size_t high_dim2) {
		clip(low_dim0, low_dim1, low_dim2, 1+(high_dim0-low_dim0), 1+(high_dim1-low_dim1), 1+(high_dim2-low_dim2));
	}
	inline void clipByBounds(MinMaxUI bounds) {
		clipByBounds(bounds.min_x, bounds.max_x, bounds.min_y, bounds.max_y, bounds.min_z, bounds.max_z);
	}
	inline void crop(size_t origin0, size_t origin1, size_t origin2, size_t dim0, size_t dim1, size_t dim2) {
		clip(origin0, origin1, origin2, dim0, dim1, dim2);
	}
	inline void cropByBounds(size_t low_dim0, size_t high_dim0, size_t low_dim1, size_t high_dim1, size_t low_dim2, size_t high_dim2) {
		clipByBounds(low_dim0, high_dim0, low_dim1, high_dim1, low_dim2, high_dim2);
	}
	inline void cropByBounds(MinMaxUI bounds) {
		clipByBounds(bounds);
	}

	void clamp(const Scalar& lowValuePtr, const Scalar& highValuePtr);
	Scalar mean(void);
	Scalar variance(void);
	Scalar stdDev(void);
	Scalar max(void);
	Scalar min(void);


	inline void SetDimensions(unsigned int* dimensions_in, unsigned int field_size) {
		size_t fsize = (field_size<=3 ? field_size : 3);
		if(Image<Scalar,3>::_dim!=NULL)
			delete [] Image<Scalar,3>::_dim;
		Image<Scalar,3>::_dim = new unsigned int[fsize];
		memcpy(Image<Scalar,3>::_dim, dimensions_in, sizeof(unsigned int)*fsize);
	}
	inline void GetDimensions(unsigned int* dimensions_out, unsigned int field_size) {
		size_t fsize = (field_size<=3 ? field_size : 3);
		memcpy(dimensions_out, Image<Scalar,3>::_dim, sizeof(unsigned int)*fsize);
	}

	inline unsigned long getDataAddress(uint x, uint y, uint z) {
		if(Image<Scalar,3>::_dataorder==ROW_MAJOR)
			return getDataAddress_rowMajor(x,y,z);
		else if(Image<Scalar,3>::_dataorder==COLUMN_MAJOR)
			return getDataAddress_columnMajor(x,y,z);
		return 0;
	}

	inline void getIndicesFromAddress(ulong address, uint& x, uint& y, uint &z) {
		if(Base::_dataorder==ROW_MAJOR) {
			getIndicesFromAddress_rowMajor(address,x,y,z);
		} else if(Base::_dataorder==COLUMN_MAJOR) {
			getIndicesFromAddress_columnMajor(address,x,y,z);
		}
		return;
	}

protected:
	typedef Image<Scalar,3> Base;
	void swap(Image3D& X);
	void swap(const Image3D<Scalar>& X);
	Scalar interpolateTrilinear(float coordX, float coordY, float coordZ);
	Scalar interpolateNearestNeighbour(float coordX, float coordY, float coordZ);

	inline void clearData() {
		if(Image<Scalar,3>::_data!=NULL) {
			delete [] Image<Scalar,3>::_data;
		}
		Image<Scalar,3>::_data=NULL;
	}
	inline void allocateData() {
		unsigned long fsize=1;
		if(Image<Scalar,3>::_dim==NULL)
			return;
		for(size_t index=0; index<3; index++)
			fsize*=Image<Scalar,3>::_dim[index];
		if(Image<Scalar,3>::_data==NULL) {
			Image<Scalar,3>::_data=new Scalar[fsize];
			//bzero(Image<Scalar,3>::_data, sizeof(Scalar)*fsize);
			memset(Image<Scalar,3>::_data, 0, sizeof(Scalar)*fsize);
		}
	}
	inline void zeroData() {
		if(Image<Scalar,3>::_data!=NULL) {
			unsigned long fsize=1;
			if(Image<Scalar,3>::_dim==NULL)
				return;
			for(size_t index=0; index<3; index++)
				fsize*=Image<Scalar,3>::_dim[index];
			//bzero(Image<Scalar,3>::_data, sizeof(Scalar));
			memset(Image<Scalar,3>::_data, 0, sizeof(Scalar)*fsize);
		}
	}
	inline void copyData(Scalar* source) {
		unsigned long fsize=1;
		if(Image<Scalar,3>::_dim==NULL)
			return;
		for(size_t index=0; index<3; index++)
			fsize*=Image<Scalar,3>::_dim[index];
		if((Image<Scalar,3>::_data!=NULL) && (source!=NULL)) {
			memcpy(Image<Scalar,3>::_data, source, fsize*sizeof(Scalar));
		}
	}
	inline unsigned long getDataAddress_columnMajor(uint x, uint y, uint z) {
		//return z*_dim[1]*_dim[0]+x*_dim[1]+y;
		return x+Image<Scalar,3>::_dim[0]*(y+Image<Scalar,3>::_dim[1]*(z+Image<Scalar,3>::_dim[2]*(0)));
	}
	inline unsigned long getDataAddress_rowMajor(uint x, uint y, uint z) {
		//return y*_dim[0]*dim_[2]+x*_dim[2]+y;
		return z+Image<Scalar,3>::_dim[2]*(y+Image<Scalar,3>::_dim[1]*(x+Image<Scalar,3>::_dim[0]*(0)));
	}
	// ??? ??? recalculate before change ??? ???  judegment: need to switch //
	inline void getIndicesFromAddress_columnMajor(unsigned long address, uint& x, uint& y, uint& z) {
		z = address / (Base::_dim[0]*Base::_dim[1]);
		unsigned long rest = address % (Base::_dim[0]*Base::_dim[1]);
		y = rest / Base::_dim[0];
		x = rest % Base::_dim[0];
	}
	inline void getIndicesFromAddress_rowMajor(unsigned long address, uint& x, uint& y, uint& z) {
		x = address / (Base::_dim[2]*Base::_dim[1]);
		unsigned long rest = address % (Base::_dim[2]*Base::_dim[1]);
		y = rest / Base::_dim[2];
		z = rest % Base::_dim[2];
	}
};

inline unsigned long getDataAddress_columnMajor_field3D(uint x, uint dimX, uint y, uint dimY, uint z, uint dimZ) {
	//return z*_dim[1]*_dim[0]+x*_dim[1]+y;
	return x+dimX*(y+dimY*(z+dimZ*(0)));
}
inline unsigned long getDataAddress_rowMajor_field3D(uint x, uint dimX, uint y, uint dimY, uint z, uint dimZ) {
	//return y*_dim[0]*dim_[2]+x*_dim[2]+y;
	return z+dimZ*(y+dimY*(x+dimX*(0)));
}
inline unsigned long getDataAddress_field3D(uint x, uint dimX, uint y, uint dimY, uint z, uint dimZ, dorder dataorder) {
	if(dataorder==ROW_MAJOR)
		return getDataAddress_rowMajor_field3D(x,dimX,y,dimY,z,dimZ);
	else if(dataorder==COLUMN_MAJOR)
		return getDataAddress_columnMajor_field3D(x,dimX,y,dimY,z,dimZ);
	return 0;
}

// ??? ??? recalculate before change ??? ???  judegment: need to switch //
inline void getIndicesFromAddress_columnMajor_field3D(unsigned long address, uint& x, uint dimX, uint& y, uint dimY, uint& z, uint dimZ) {
	z = address / (dimX*dimY);
	unsigned long rest = address % (dimX*dimY);
	y = rest / dimX;
	x = rest % dimX;
}
inline void getIndicesFromAddress_rowMajor_field3D(unsigned long address, uint& x, uint dimX, uint& y, uint dimY, uint& z, uint dimZ) {
	x = address / (dimZ*dimY);
	unsigned long rest = address % (dimZ*dimY);
	y = rest / dimZ;
	z = rest % dimZ;
}
inline void getIndicesFromAddress_field3D(unsigned long address,uint &x, uint dimX, uint &y, uint dimY, uint &z, uint dimZ, dorder dataorder) {
	if(dataorder==ROW_MAJOR)
		getIndicesFromAddress_rowMajor_field3D(address,x,dimX,y,dimY,z,dimZ);
	else if(dataorder==COLUMN_MAJOR)
		getIndicesFromAddress_columnMajor_field3D(address,x,dimX,y,dimY,z,dimZ);
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
inline void Image3D<char>::SetData(char* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new char[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(char)*fsize);
}
template<>
inline void Image3D<char>::GetData(char* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(char)*fsize);
}
template<>
inline char Image3D<char>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<char>;



template<>
inline void Image3D<uchar>::SetData(uchar* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new uchar[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(uchar)*fsize);
}
template<>
inline void Image3D<uchar>::GetData(uchar* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(uchar)*fsize);
}
template<>
inline uchar Image3D<uchar>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<uchar>;



template<>
inline void Image3D<short>::SetData(short* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new short[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(short)*fsize);
}
template<>
inline void Image3D<short>::GetData(short* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(short)*fsize);
}
template<>
inline short Image3D<short>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<short>;



template<>
inline void Image3D<ushort>::SetData(ushort* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new ushort[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(ushort)*fsize);
}
template<>
inline void Image3D<ushort>::GetData(ushort* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(ushort)*fsize);
}
template<>
inline ushort Image3D<ushort>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<ushort>;



template<>
inline void Image3D<int>::SetData(int* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new int[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(int)*fsize);
}
template<>
inline void Image3D<int>::GetData(int* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(int)*fsize);
}
template<>
inline int Image3D<int>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<int>;



template<>
inline void Image3D<uint>::SetData(uint* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new uint[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(uint)*fsize);
}
template<>
inline void Image3D<uint>::GetData(uint* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(uint)*fsize);
}
template<>
inline uint Image3D<uint>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<uint>;



template<>
inline void Image3D<long>::SetData(long* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new long[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(long)*fsize);
}
template<>
inline void Image3D<long>::GetData(long* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(long)*fsize);
}
template<>
inline long Image3D<long>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<long>;



template<>
inline void Image3D<ulong>::SetData(ulong* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new ulong[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(ulong)*fsize);
}
template<>
inline void Image3D<ulong>::GetData(ulong* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(ulong)*fsize);
}
template<>
inline ulong Image3D<ulong>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<ulong>;



template<>
inline void Image3D<float>::SetData(float* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new float[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(float)*fsize);
}
template<>
inline void Image3D<float>::GetData(float* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(float)*fsize);
}
template<>
inline float Image3D<float>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<float>;



template<>
inline void Image3D<double>::SetData(double* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	if(_data!=NULL)
		delete [] _data;
	_data = new double[_dim[0]*_dim[1]*_dim[2]];
	memcpy(_data, data_in, sizeof(double)*fsize);
}
template<>
inline void Image3D<double>::GetData(double* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2])) ? field_size : (_dim[0]*_dim[1]*_dim[2]));
	memcpy(data_out, _data, sizeof(double)*fsize);
}
template<>
inline double Image3D<double>::GetImageValue(uint coordX, uint coordY, uint coordZ) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]))
		return 0;
	return _data[index];
}
template class Image3D<double>;
*/

/*
template class Image3D<char>;
template class Image3D<uchar>;
template class Image3D<short>;
template class Image3D<ushort>;
template class Image3D<int>;
template class Image3D<uint>;
template class Image3D<long>;
template class Image3D<ulong>;
template class Image3D<float>;
template class Image3D<double>;
*/


#include "image3D_utils.hpp"

#endif /* IMAGE3D_HPP_ */
