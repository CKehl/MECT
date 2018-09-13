/**
 * Created on 14 Jan 2018
 * author: christian
 */
#include "image2D.hpp"
#include <cmath>
#include <float.h>
#include <stdint.h>
//#include <cassert>
#include <assert.h>

template<typename Scalar>
Image2D<Scalar>::Image2D() {
	Image<Scalar,2>::_dim = NULL;
	Image<Scalar,2>::_data = NULL;
	Image<Scalar,2>::_datatype =NONE;
	Image<Scalar,2>::_dataorder = UNDEF_ORDER;
}

template<typename Scalar>
Image2D<Scalar>::Image2D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData) {
	Image<Scalar,2>::_dim = dimensions;
	Image<Scalar,2>::_data = originalData;
	Image<Scalar,2>::_datatype = datatype;
	Image<Scalar,2>::_dataorder = dataorder;
}

template<typename Scalar>
Image2D<Scalar>::Image2D(uint dimX, uint dimY, dtype datatype, dorder dataorder, Scalar* originalData) {
	Image<Scalar,2>::_data = NULL;
	Image<Scalar,2>::_dim = new uint[2];
	Image<Scalar,2>::_dim[0] = dimX; Image<Scalar,2>::_dim[1] = dimY;
	Image<Scalar,2>::_datatype = datatype;
	Image<Scalar,2>::_dataorder = dataorder;
	allocateData();
	if(originalData!=NULL)
		copyData(originalData);
}

template<typename Scalar>
Image2D<Scalar>::Image2D(const Image2D& X) {
	Image<Scalar,2>::_dim = X._dim;
	Image<Scalar,2>::_data = X._data;
	Image<Scalar,2>::_datatype = X._datatype;
	Image<Scalar,2>::_dataorder = X._dataorder;
}

template<typename Scalar>
Image2D<Scalar>::~Image2D() {
	clearData();
	if(Image<Scalar,2>::_dim!=NULL)
		delete [] Image<Scalar,2>::_dim;
	Image<Scalar,2>::_dim=NULL;
}

template<typename Scalar>
bool Image2D<Scalar>::operator==(Image2D<Scalar>& X) {
	bool result = true;
	if(Image<Scalar,2>::_dim != X.Image<Scalar,2>::_dim) {
		if((Image<Scalar,2>::_dim!=NULL) && (X.Image<Scalar,2>::_dim!=NULL)) {
			for(size_t i=0; i<2; i++)
				result &= (Image<Scalar,2>::_dim[i]==X.Image<Scalar,2>::_dim[i]);
		} else {
			result &= false;
		}
	}
	result &= (Image<Scalar,2>::_datatype==X.Image<Scalar,2>::_datatype);
	result &= (Image<Scalar,2>::_dataorder==X.Image<Scalar,2>::_dataorder);
	if(Image<Scalar,2>::_data != X.Image<Scalar,2>::_data) {
		unsigned long fsize_this=1, fsize_X=1;
		if((Image<Scalar,2>::_dim!=NULL) && (X.Image<Scalar,2>::_dim!=NULL)) {
			for(size_t index=0; index<2; index++) {
				fsize_this*=Image<Scalar,2>::_dim[index];
				fsize_X*=X.Image<Scalar,2>::_dim[index];
			}
			if(fsize_this==fsize_X) {
				for(unsigned long i=0; i<fsize_this;i++) {
					result &= (Image<Scalar,2>::_data[i]==X.Image<Scalar,2>::_data[i]);
				}
			} else {
				result &= false;
			}
		} else {
			result &= false;
		}
	}
	return result;
}

template<typename Scalar>
bool Image2D<Scalar>::operator!=(Image2D<Scalar>& X) {
	bool result = true;
	if(Image<Scalar,2>::_dim != X.Image<Scalar,2>::_dim) {
		if((Image<Scalar,2>::_dim!=NULL) && (X.Image<Scalar,2>::_dim!=NULL)) {
			for(size_t i=0; i<2; i++)
				result &= (Image<Scalar,2>::_dim[i]==X.Image<Scalar,2>::_dim[i]);
		} else {
			result &= false;
		}
	}
	result &= (Image<Scalar,2>::_datatype==X.Image<Scalar,2>::_datatype);
	result &= (Image<Scalar,2>::_dataorder==X.Image<Scalar,2>::_dataorder);
	if(Image<Scalar,2>::_data != X.Image<Scalar,2>::_data) {
		unsigned long fsize_this=1, fsize_X=1;
		if((Image<Scalar,2>::_dim!=NULL) && (X.Image<Scalar,2>::_dim!=NULL)) {
			for(size_t index=0; index<2; index++) {
				fsize_this*=Image<Scalar,2>::_dim[index];
				fsize_X*=X.Image<Scalar,2>::_dim[index];
			}
			if(fsize_this==fsize_X) {
				for(unsigned long i=0; i<fsize_this;i++) {
					result &= (Image<Scalar,2>::_data[i]==X.Image<Scalar,2>::_data[i]);
				}
			} else {
				result &= false;
			}
		} else {
			result &= false;
		}
	}
	return !result;
}

template<typename Scalar>
void Image2D<Scalar>::setImageValue(Scalar value, uint coordX, uint coordY) {
	if(Image<Scalar,2>::_dim==NULL)
		return;
	if(Image<Scalar,2>::_data==NULL)
		return;
	unsigned long index = 0;
	if(Image<Scalar,2>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(Image<Scalar,2>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		return;
	}
	if(index >= (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]))
		return;
	Image<Scalar,2>::_data[index] = value;
}

template<typename Scalar>
void Image2D<Scalar>::createImage(void) {
	if(Image<Scalar,2>::_dim==NULL)
		return;
	clearData();
	if(Image<Scalar,2>::_datatype!=NONE)
		allocateData();
}

template<typename Scalar>
void Image2D<Scalar>::createImage(uint* dim, dtype datatype) {
	if((dim==NULL) || (datatype==NONE))
		return;
	clearData();
	if(Image<Scalar,2>::_dim!=NULL)
		delete [] Image<Scalar,2>::_dim;
	Image<Scalar,2>::_dim = new uint[2];
	memcpy(Image<Scalar,2>::_dim, dim, sizeof(uint)*2);
	Image<Scalar,2>::_datatype=datatype;
	allocateData();
}

template<typename Scalar>
void Image2D<Scalar>::clean(void) {
	zeroData();
}

template<typename Scalar>
void Image2D<Scalar>::pad(uint totalPadX, uint totalPadY) {
	unsigned long fsize_old=1;
	if(Image<Scalar,2>::_dim==NULL)
		return;
	if(Image<Scalar,2>::_data==NULL)
		return;
	if(Image<Scalar,2>::_datatype==NONE)
		return;
	if(Image<Scalar,2>::_dataorder==UNDEF_ORDER)
		return;
	for(size_t index=0; index<2; index++)
		fsize_old*=Image<Scalar,2>::_dim[index];
	uint* ndim = new uint[2];
	ndim[0]=Image<Scalar,2>::_dim[0]+totalPadX; ndim[1]=Image<Scalar,2>::_dim[1]+totalPadY;
	unsigned long fsize_new = ndim[0]*ndim[1];
	uint offsetX = uint(std::floor(totalPadX/2));
	uint offsetY = uint(std::floor(totalPadY/2));

	Scalar* tbuffer = new Scalar[fsize_new];
	ulong oldAddress=0, newAddress=0;
		for(uint y=0; y<Image<Scalar,2>::_dim[1]; y++) {
			for(uint x=0; x<Image<Scalar,2>::_dim[0]; x++) {
				oldAddress = getDataAddress(x,y);
				newAddress = getDataAddress_field2D(x+offsetX, ndim[0], y+offsetY, ndim[1], Image<Scalar,2>::_dataorder);
				tbuffer[newAddress] = Image<Scalar,2>::_data[oldAddress];
			}
		}
	delete [] Image<Scalar,2>::_data; Image<Scalar,2>::_data=NULL;
	memcpy(Image<Scalar,2>::_dim, ndim, sizeof(uint)*2);
	Image<Scalar,2>::_data = new Scalar[fsize_new];
	memcpy(Image<Scalar,2>::_data, tbuffer, sizeof(Scalar)*fsize_new);
}

template<typename Scalar>
void Image2D<Scalar>::swap(Image2D<Scalar>& X) {
	if((Image<Scalar,2>::_data!=NULL) && (X.Image<Scalar,2>::_data!=NULL)) {
		ulong fsize_this = Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1];
		ulong fsize_X = X.Image<Scalar,2>::_dim[0]*X.Image<Scalar,2>::_dim[1];
		Scalar* temp = new Scalar[fsize_this];
		memcpy(temp, Image<Scalar,2>::_data, sizeof(Scalar)*fsize_this);
		clearData();
		Image<Scalar,2>::_data = new Scalar[fsize_X];
		memcpy(Image<Scalar,2>::_data, X.Image<Scalar,2>::_data, sizeof(Scalar)*fsize_X);
		X.clearData();
		X.Image<Scalar,2>::_data = new Scalar[fsize_this];
		memcpy(X.Image<Scalar,2>::_data, temp, sizeof(Scalar)*fsize_this);
		delete [] temp;
	}
	uint tempDim[2];
	if((Image<Scalar,2>::_dim!=NULL) && (X.Image<Scalar,2>::_dim!=NULL)) {
		memcpy((uint*)tempDim, Image<Scalar,2>::_dim, sizeof(uint)*2);
		memcpy(Image<Scalar,2>::_dim, X.Image<Scalar,2>::_dim, sizeof(uint)*2);
		memcpy(X.Image<Scalar,2>::_dim, (uint*)tempDim, sizeof(uint)*2);
	}
	std::swap(Image<Scalar,2>::_datatype, X.Image<Scalar,2>::_datatype);
	std::swap(Image<Scalar,2>::_dataorder, X.Image<Scalar,2>::_dataorder);
}

template<typename Scalar>
void Image2D<Scalar>::swap(const Image2D<Scalar>& X) {

	if(X.Image<Scalar,2>::_dim!=NULL) {
		memcpy(Image<Scalar,2>::_dim, X.Image<Scalar,2>::_dim, sizeof(uint)*2);

		if(Image<Scalar,2>::_data!=NULL) {
			delete [] Image<Scalar,2>::_data;
			Image<Scalar,2>::_data = new Scalar[Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]];
		}
	} else {
		if(Image<Scalar,2>::_dim!=NULL)
			delete [] Image<Scalar,2>::_dim;
		Image<Scalar,2>::_dim=NULL;
	}

	if(X.Image<Scalar,2>::_data!=NULL) {
		memcpy(Image<Scalar,2>::_data, X.Image<Scalar,2>::_data, sizeof(Scalar)*Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]);
	} else {
		if(Image<Scalar,2>::_data!=NULL)
			delete [] Image<Scalar,2>::_data;
		Image<Scalar,2>::_data=NULL;
	}

	//std::swap(Image<Scalar,2>::_dim, X.Image<Scalar,2>::_dim);
	//std::swap(Image<Scalar,2>::_data, X.Image<Scalar,2>::_data);
	Image<Scalar,2>::_datatype = X.Image<Scalar,2>::_datatype;
	Image<Scalar,2>::_dataorder = X.Image<Scalar,2>::_dataorder;
}

template<typename Scalar>
Scalar& Image2D<Scalar>::operator[](unsigned long index) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(index < (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]));

	return Image<Scalar,2>::_data[index];
}

template<typename Scalar>
Scalar& Image2D<Scalar>::operator()(uint coordX, uint coordY) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(getDataorder()!=UNDEF_ORDER);
	unsigned long index = 0;
	if(Image<Scalar,2>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(Image<Scalar,2>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		index = 0;
	}
	assert(index < (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]));

	return Image<Scalar,2>::_data[index];
}

template<typename Scalar>
Scalar& Image2D<Scalar>::getImageValue(uint coordX, uint coordY) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(getDataorder()!=UNDEF_ORDER);
	unsigned long index = 0;
	if(Image<Scalar,2>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY);
	} else if(Image<Scalar,2>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY);
	} else {
		index = 0;
	}
	assert(index < (Image<Scalar,2>::_dim[0]*Image<Scalar,2>::_dim[1]));

	return Image<Scalar,2>::_data[index];
}

template<typename Scalar>
Scalar Image2D<Scalar>::interpolate(float coordX, float coordY, interpolationMode2D mode) {
	switch(mode) {
	case NEAREST_NEIGHBOUR_2D: return interpolateNearestNeighbour(coordX, coordY);
	case BILINEAR: return interpolateBilinear(coordX, coordY);
	}
	return Scalar(0);
}

template<typename Scalar>
Scalar Image2D<Scalar>::interpolateBilinear(float coordX, float coordY) {
	uint xf = uint(floor(coordX)); uint yf = uint(floor(coordY));
	float x_offset = (coordX - float(xf)) - 0.5; float y_offset = (coordY - float(yf)) - 0.5;
	float pixDistance = .0f, kd=.0f, ld=.0f;
	float maxdist = sqrt(3*1.5);
	//float interpCumulative=.0f;
	double result=0;
	Scalar kernel[3][3];
	uint kernel_k=0, kernel_l=0;
	for(int k=(int(xf)-1); k<(int(xf)+2); k++) {
		for(int l=(int(yf)-1); l<(int(yf)+2); l++) {
				//zero-based addresses
				kernel_k=k-(xf-1);
				kernel_l=l-(yf-1);
				if((k<0)||(uint(k)>(Image<Scalar,2>::_dim[0]-1)) || (l<0)||(uint(l)>(Image<Scalar,2>::_dim[1]-1))) {
					kernel[kernel_k][kernel_l]=0;
				} else {
					kernel[kernel_k][kernel_l]=Image<Scalar,2>::_data[getDataAddress(k,l)];
				}
		}
	}
	for(uint k=0;k<3;k++) {
		for(uint l=0;l<3;l++) {
				kd=float(k-1)+x_offset;
				ld=float(l-1)+y_offset;
				pixDistance=sqrt(kd*kd+ld*ld)/maxdist;
				result+=floor( (1.0f-pixDistance)*kernel[k][l] );
		}
	}
	return Scalar(floor(result/9.0));
}

template<typename Scalar>
Scalar Image2D<Scalar>::interpolateNearestNeighbour(float coordX, float coordY) {
	uint xr = uint(round(coordX)); uint yr = uint(round(coordY));
	return getImageValue(xr, yr);
}

template<typename Scalar>
void Image2D<Scalar>::clamp(const Scalar& lowValuePtr, const Scalar& highValuePtr) {
	if(Image<Scalar,2>::_dim==NULL)
		return;
	if(Image<Scalar,2>::_data==NULL)
		return;
	unsigned long fsize=1;
	for(size_t index=0; index<2; index++)
		fsize*=Image<Scalar,2>::_dim[index];
	for(ulong i=0; i<fsize; i++) {
		Image<Scalar,2>::_data[i]=clampByValue(Image<Scalar,2>::_data[i], lowValuePtr, highValuePtr);
	}
}

template<typename Scalar>
Scalar Image2D<Scalar>::mean(void) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	unsigned long fsize=1;
	for(size_t index=0; index<2; index++)
		fsize*=Image<Scalar,2>::_dim[index];
	double result=0;
	for(uint i=0; i<fsize; i++)
		result+=double(Image<Scalar,2>::_data[i]);
	result/=double(fsize);
	return Scalar(result);
}

template<typename Scalar>
Scalar Image2D<Scalar>::variance(void) {
	double result = 0.0;
    Scalar _mean = mean();
	unsigned long fsize=1;
	for(size_t index=0; index<2; index++)
		fsize*=Image<Scalar,2>::_dim[index];
    for(uint i=0; i<fsize; i++)
    	result += pow(double(Image<Scalar,2>::_data[i] - _mean),2);
    return Scalar(result/double(fsize));
}

template<typename Scalar>
Scalar Image2D<Scalar>::stdDev(void) {
	return Scalar(sqrt(variance()));
}

template<typename Scalar>
Scalar Image2D<Scalar>::max(void) {
	Scalar result = -Scalar(DBL_MAX);
	unsigned long fsize=1;
	for(size_t index=0; index<2; index++)
		fsize*=Image<Scalar,2>::_dim[index];
	for(uint i=0; i<fsize; i++)
		result = std::max(result, Image<Scalar,2>::_data[i]);
	return result;
}

template<typename Scalar>
Scalar Image2D<Scalar>::min(void) {
	Scalar result = Scalar(DBL_MAX);
	unsigned long fsize=1;
	for(size_t index=0; index<2; index++)
		fsize*=Image<Scalar,2>::_dim[index];
	for(uint i=0; i<fsize; i++)
		result = std::min(result, Image<Scalar,2>::_data[i]);
	return result;
}

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
