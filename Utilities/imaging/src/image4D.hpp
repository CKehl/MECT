/**
 * Created on 18 Jan 2018
 * author: christian
 */

#ifndef IMAGE4D_HPP_
#define IMAGE4D_HPP_

#include "image.hpp"
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif
#include <iostream>
#include <string.h>

typedef enum _interpolation_mode_4D {
	NEAREST_NEIGHBOUR_4D=0,
	QUADRILINEAR=1
} interpolationMode4D;

template<typename Scalar>
class Image4D : public Image<Scalar, 4> {
public:
	Image4D();
	/**
	 * Image4D shallow-copy constructor
	 */
	Image4D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	/**
	 * Image4D deep-copy constructor
	 */
	Image4D(uint dimX, uint dimY, uint dimZ, uint dimW, dtype datatype, dorder dataorder, Scalar* originalData=NULL);
	Image4D(const Image4D& X);
	virtual ~Image4D();
	
	inline Image4D<Scalar>& operator=(Image4D<Scalar>& X) { swap(X); return *this; }
	inline Image4D<Scalar>& operator=(const Image4D<Scalar>& X) {
		swap(X);
		return *this;
	}
	bool operator==(Image4D& X);
	bool operator!=(Image4D& X);

	Scalar& operator[](unsigned long index);
	Scalar& operator()(uint coordX, uint coordY, uint coordZ, uint coordW);
	inline friend std::ostream& operator<< (std::ostream& os, const Image4D& rhs) {
		if(rhs._dim!=NULL)
			os << "shape = ("<< rhs._dim[0] <<","<< rhs._dim[1] <<","<< rhs._dim[2] <<","<< rhs._dim[3] <<"); ";
		if(rhs._data!=NULL) {

		}
		os << "dtype = " << getDTypeString(rhs._datatype);
		os << "dorder = " << getDOrderString(rhs._dataorder);
		return os;
	}

	inline void setDimension(size_t index, uint value) {
		if(index<0)
			Image<Scalar,4>::_dim[0] = value;
		if(index>3)
			return;
		Image<Scalar,4>::_dim[index]=value;
	}
	inline void setDimensions(uint* dim, bool deepCopy=false) {
		if(dim!=NULL) {
			if(deepCopy) {
				if(Image<Scalar,4>::_dim!=NULL)
					delete [] Image<Scalar,4>::_dim;
				Image<Scalar,4>::_dim = new uint[4];
				memcpy(Image<Scalar,4>::_dim, dim, 4*sizeof(uint));
			} else
				Image<Scalar,4>::_dim=dim;
		}
	}
	inline void setDimensions(uint dimX, uint dimY, uint dimZ, uint dimW) {
		if(Image<Scalar,4>::_dim!=NULL) {
			delete [] Image<Scalar,4>::_dim;
		}
		Image<Scalar,4>::_dim = new unsigned int[4];
		Image<Scalar,4>::_dim[0]=dimX; Image<Scalar,4>::_dim[1]=dimY; Image<Scalar,4>::_dim[2]=dimZ; Image<Scalar,4>::_dim[3]=dimW;
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
			Image<Scalar,4>::_data=data;
		}
	}
	inline void setDatatype(dtype datatype) {
		Image<Scalar,4>::_datatype = datatype;
	}
	inline void setDataorder(dorder dataorder) {
		Image<Scalar,4>::_dataorder = dataorder;
	}
	inline uint* getDimensions() {
		return Image<Scalar,4>::_dim;
	}
	inline uint getDimension(size_t index) {
		if(index<0)
			return Image<Scalar,4>::_dim[0];
		if(index>3)
			return Image<Scalar,4>::_dim[3];
		return Image<Scalar,4>::_dim[index];
	}
	inline virtual Scalar* getData() {
		return Image<Scalar,4>::_data;
	}
	inline virtual dtype getDatatype() {
		return Image<Scalar,4>::_datatype;
	}
	inline virtual dorder getDataorder() {
		return Image<Scalar,4>::_dataorder;
	}

	// WRAPPING funcs - to be 'specialised'
	/*
	inline void SetData(Scalar* data_in, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3])) ? field_size : (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]));
		if(Image<Scalar,4>::_data!=NULL)
			delete [] Image<Scalar,4>::_data;
		Image<Scalar,4>::_data = new Scalar[Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]];
		memcpy(Image<Scalar,4>::_data, data_in, sizeof(Scalar)*fsize);
	}
	*/
	inline void SetData(Scalar* data_in, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3])) ? field_size : (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]));
		if(Image<Scalar,4>::_data!=NULL)
			delete [] Image<Scalar,4>::_data;
		Image<Scalar,4>::_data = new Scalar[Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]];
		memcpy(Image<Scalar,4>::_data, data_in, sizeof(Scalar)*fsize);
		Image<Scalar,4>::_dataorder = ROW_MAJOR;
	}
	inline void GetData(Scalar* data_out, unsigned long field_size) {
		size_t fsize = ((field_size <= (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3])) ? field_size : (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]));
		memcpy(data_out, Image<Scalar,4>::_data, sizeof(Scalar)*fsize);
	}
	inline Scalar GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
		if(Image<Scalar,4>::_dim==NULL)
			return 0;
		if(Image<Scalar,4>::_data==NULL)
			return 0;
		unsigned long index = 0;
		if(Image<Scalar,4>::_dataorder == ROW_MAJOR) {
			index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
		} else if(Image<Scalar,4>::_dataorder == COLUMN_MAJOR) {
			index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
		} else {
			return 0;
		}
		if(index >= (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]))
			return 0;
		return Image<Scalar,4>::_data[index];
	}
	// WRAPPING funcs - to be 'specialised'

	Scalar interpolate(float coordX, float coordY, float coordZ, float coordW, interpolationMode4D mode);
	Scalar& getImageValue(uint coordX, uint coordY, uint coordZ, uint coordW);

	void setImageValue(Scalar value, uint coordX, uint coordY, uint coordZ, uint coordW);
	void createImage(void);
	void createImage(uint* dim, dtype datatype);
	void clean(void);
	void pad(uint totalPadX, uint totalPadY, uint totalPadZ, uint totalPadW);

	void clamp(const Scalar& lowValuePtr, const Scalar& highValuePtr);
	Scalar mean(void);
	Scalar variance(void);
	Scalar stdDev(void);
	Scalar max(void);
	Scalar min(void);


	inline void SetDimensions(unsigned int* dimensions_in, unsigned int field_size) {
		size_t fsize = (field_size<=4 ? field_size : 4);
		if(Image<Scalar,4>::_dim!=NULL)
			delete [] Image<Scalar,4>::_dim;
		Image<Scalar,4>::_dim = new unsigned int[fsize];
		memcpy(Image<Scalar,4>::_dim, dimensions_in, sizeof(unsigned int)*fsize);
	}
	inline void GetDimensions(unsigned int* dimensions_out, unsigned int field_size) {
		size_t fsize = (field_size<=4 ? field_size : 4);
		memcpy(dimensions_out, Image<Scalar,4>::_dim, sizeof(unsigned int)*fsize);
	}

	inline unsigned long getDataAddress(uint x, uint y, uint z, uint w) {
		if(Image<Scalar,4>::_dataorder==ROW_MAJOR)
			return getDataAddress_rowMajor(x,y,z,w);
		else if(Image<Scalar,4>::_dataorder==COLUMN_MAJOR)
			return getDataAddress_columnMajor(x,y,z,w);
		return 0;
	}

	inline void getIndicesFromAddress(ulong address, uint& x, uint& y, uint &z, uint &w) {
		if(Base::_dataorder==ROW_MAJOR) {
			getIndicesFromAddress_rowMajor(address,x,y,z,w);
		} else if(Base::_dataorder==COLUMN_MAJOR) {
			getIndicesFromAddress_columnMajor(address,x,y,z,w);
		}
		return;
	}

protected:
	typedef Image<Scalar,4> Base;
	void swap(Image4D<Scalar>& X);
	void swap(const Image4D<Scalar>& X);
	Scalar interpolateQuadrilinear(float coordX, float coordY, float coordZ, float coordW);
	Scalar interpolateNearestNeighbour(float coordX, float coordY, float coordZ, float coordW);

	inline void clearData() {
		if(Image<Scalar,4>::_data!=NULL) {
			delete [] Image<Scalar,4>::_data;
		}
		Image<Scalar,4>::_data=NULL;
	}
	inline void allocateData() {
		unsigned long fsize=1;
		if(Image<Scalar,4>::_dim==NULL)
			return;
		for(size_t index=0; index<4; index++)
			fsize*=Image<Scalar,4>::_dim[index];
		if(Image<Scalar,4>::_data==NULL) {
			Image<Scalar,4>::_data=new Scalar[fsize];
			//bzero(Image<Scalar,4>::_data, sizeof(Scalar)*fsize);
			memset(Image<Scalar,4>::_data, 0, sizeof(Scalar)*fsize);
		}
	}
	inline void zeroData() {
		if(Image<Scalar,4>::_data!=NULL) {
			unsigned long fsize=1;
			if(Image<Scalar,4>::_dim==NULL)
				return;
			for(size_t index=0; index<4; index++)
				fsize*=Image<Scalar,4>::_dim[index];
			//bzero(Image<Scalar,4>::_data, sizeof(Scalar));
			memset(Image<Scalar,4>::_data, 0, sizeof(Scalar)*fsize);
		}
	}
	inline void copyData(Scalar* source) {
		unsigned long fsize=1;
		if(Image<Scalar,4>::_dim==NULL)
			return;
		for(size_t index=0; index<4; index++)
			fsize*=Image<Scalar,4>::_dim[index];
		if((Image<Scalar,4>::_data!=NULL) && (source!=NULL)) {
			memcpy(Image<Scalar,4>::_data, source, fsize*sizeof(Scalar));
		}
	}
	inline unsigned long getDataAddress_columnMajor(uint x, uint y, uint z, uint w) {
		//return z*_dim[1]*_dim[0]+x*_dim[1]+y;
		return x+Image<Scalar,4>::_dim[0]*(y+Image<Scalar,4>::_dim[1]*(z+Image<Scalar,4>::_dim[2]*(w+Image<Scalar,4>::_dim[3]*(0))));
	}
	inline unsigned long getDataAddress_rowMajor(uint x, uint y, uint z, uint w) {
		//return y*_dim[0]*dim_[2]+x*_dim[2]+z;
		return w+Image<Scalar,4>::_dim[3]*(z+Image<Scalar,4>::_dim[2]*(y+Image<Scalar,4>::_dim[1]*(x+Image<Scalar,4>::_dim[0]*(0))));
	}
	// ??? ??? recalculate before change ??? ???  judegment: need to switch //
	inline void getIndicesFromAddress_columnMajor(unsigned long address, uint& x, uint& y, uint& z, uint& w) {
		w = address / (Base::_dim[0]*Base::_dim[1]*Base::_dim[2]);
		unsigned long rest = address % (Base::_dim[0]*Base::_dim[1]*Base::_dim[2]);
		z = rest / (Base::_dim[0]*Base::_dim[1]);
		rest = rest % (Base::_dim[0]*Base::_dim[1]);
		y = rest / Base::_dim[0];
		x = rest % Base::_dim[0];
	}
	inline void getIndicesFromAddress_rowMajor(unsigned long address, uint& x, uint& y, uint& z, uint& w) {
		x = address / (Base::_dim[3]*Base::_dim[2]*Base::_dim[1]);
		unsigned long rest = address % (Base::_dim[3]*Base::_dim[2]*Base::_dim[1]);
		y = rest / (Base::_dim[3]*Base::_dim[2]);
		rest = rest % (Base::_dim[3]*Base::_dim[2]);
		z = rest / Base::_dim[3];
		w = rest % Base::_dim[3];
	}
};

inline unsigned long getDataAddress_columnMajor_field4D(uint x, uint dimX, uint y, uint dimY, uint z, uint dimZ, uint w, uint dimW) {
	//return z*_dim[1]*_dim[0]+x*_dim[1]+y;
	return x+dimX*(y+dimY*(z+dimZ*(w+dimW*(0))));
}
inline unsigned long getDataAddress_rowMajor_field4D(uint x, uint dimX, uint y, uint dimY, uint z, uint dimZ, uint w, uint dimW) {
	//return y*_dim[0]*dim_[2]+x*_dim[2]+y;
	return w+dimW*(z+dimZ*(y+dimY*(x+dimX*(0))));
}
inline unsigned long getDataAddress_field4D(uint x, uint dimX, uint y, uint dimY, uint z, uint dimZ, uint w, uint dimW, dorder dataorder) {
	if(dataorder==ROW_MAJOR)
		return getDataAddress_rowMajor_field4D(x,dimX,y,dimY,z,dimZ,w,dimW);
	else if(dataorder==COLUMN_MAJOR)
		return getDataAddress_columnMajor_field4D(x,dimX,y,dimY,z,dimZ,w,dimW);
	return 0;
}

// ??? ??? recalculate before change ??? ??? //
inline void getIndicesFromAddress_columnMajor_field4D(unsigned long address, uint& x, uint dimX, uint& y, uint dimY, uint& z, uint dimZ, uint& w, uint dimW) {
	w = address / (dimX*dimY*dimZ);
	unsigned long rest = address % (dimX*dimY*dimZ);
	z = rest / (dimX*dimY);
	rest = rest % (dimX*dimY);
	y = rest / dimX;
	x = rest % dimX;
}
inline void getIndicesFromAddress_rowMajor_field4D(unsigned long address, uint& x, uint dimX, uint& y, uint dimY, uint& z, uint dimZ, uint& w, uint dimW) {
	x = address / (dimW*dimZ*dimY);
	unsigned long rest = address % (dimW*dimZ*dimY);
	y = rest / (dimW*dimZ);
	rest = rest % (dimW*dimZ);
	z = rest / dimW;
	w = rest % dimW;
}
inline void getIndicesFromAddress_field4D(unsigned long address,uint &x, uint dimX, uint &y, uint dimY, uint &z, uint dimZ, uint &w, uint dimW, dorder dataorder) {
	if(dataorder==ROW_MAJOR)
		getIndicesFromAddress_rowMajor_field4D(address,x,dimX,y,dimY,z,dimZ,w,dimW);
	else if(dataorder==COLUMN_MAJOR)
		getIndicesFromAddress_columnMajor_field4D(address,x,dimX,y,dimY,z,dimZ,w,dimW);
}


/*
template<>
inline void Image4D<char>::SetData(char* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new char[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(char)*fsize);
}
template<>
inline void Image4D<char>::GetData(char* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(char)*fsize);
}
template<>
inline char Image4D<char>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<char>;



template<>
inline void Image4D<uchar>::SetData(uchar* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new uchar[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(uchar)*fsize);
}
template<>
inline void Image4D<uchar>::GetData(uchar* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(uchar)*fsize);
}
template<>
inline uchar Image4D<uchar>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<uchar>;



template<>
inline void Image4D<short>::SetData(short* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new short[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(short)*fsize);
}
template<>
inline void Image4D<short>::GetData(short* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(short)*fsize);
}
template<>
inline short Image4D<short>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<short>;



template<>
inline void Image4D<ushort>::SetData(ushort* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new ushort[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(ushort)*fsize);
}
template<>
inline void Image4D<ushort>::GetData(ushort* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(ushort)*fsize);
}
template<>
inline ushort Image4D<ushort>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<ushort>;



template<>
inline void Image4D<int>::SetData(int* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new int[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(int)*fsize);
}
template<>
inline void Image4D<int>::GetData(int* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(int)*fsize);
}
template<>
inline int Image4D<int>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<int>;



template<>
inline void Image4D<uint>::SetData(uint* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new uint[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(uint)*fsize);
}
template<>
inline void Image4D<uint>::GetData(uint* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(uint)*fsize);
}
template<>
inline uint Image4D<uint>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<uint>;



template<>
inline void Image4D<long>::SetData(long* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new long[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(long)*fsize);
}
template<>
inline void Image4D<long>::GetData(long* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(long)*fsize);
}
template<>
inline long Image4D<long>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<long>;



template<>
inline void Image4D<ulong>::SetData(ulong* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new ulong[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(ulong)*fsize);
}
template<>
inline void Image4D<ulong>::GetData(ulong* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(ulong)*fsize);
}
template<>
inline ulong Image4D<ulong>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<ulong>;



template<>
inline void Image4D<float>::SetData(float* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new float[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(float)*fsize);
}
template<>
inline void Image4D<float>::GetData(float* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(float)*fsize);
}
template<>
inline float Image4D<float>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<float>;



template<>
inline void Image4D<double>::SetData(double* data_in, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	if(_data!=NULL)
		delete [] _data;
	_data = new double[_dim[0]*_dim[1]*_dim[2]*_dim[3]];
	memcpy(_data, data_in, sizeof(double)*fsize);
}
template<>
inline void Image4D<double>::GetData(double* data_out, unsigned long field_size) {
	size_t fsize = ((field_size <= (_dim[0]*_dim[1]*_dim[2]*_dim[3])) ? field_size : (_dim[0]*_dim[1]*_dim[2]*_dim[3]));
	memcpy(data_out, _data, sizeof(double)*fsize);
}
template<>
inline double Image4D<double>::GetImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(_dim==NULL)
		return 0;
	if(_data==NULL)
		return 0;
	unsigned long index = 0;
	if(_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return 0;
	}
	if(index >= (_dim[0]*_dim[1]*_dim[2]*_dim[3]))
		return 0;
	return _data[index];
}
template class Image4D<double>;
*/

/*
template class Image4D<char>;
template class Image4D<uchar>;
template class Image4D<short>;
template class Image4D<ushort>;
template class Image4D<int>;
template class Image4D<uint>;
template class Image4D<long>;
template class Image4D<ulong>;
template class Image4D<float>;
template class Image4D<double>;
*/

#include "image4D_utils.hpp"



#endif /* IMAGE4D_HPP_ */
