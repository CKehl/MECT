/**
 * Created on 18 Jan 2018
 * author: christian
 */
#include "image4D.hpp"
#include <cmath>
#include <float.h>
#include <stdint.h>
//#include <cassert>
#include <assert.h>

template<typename Scalar>
Image4D<Scalar>::Image4D() {
	Image<Scalar,4>::_dim = NULL;
	Image<Scalar,4>::_data = NULL;
	Image<Scalar,4>::_datatype = NONE;
	Image<Scalar,4>::_dataorder = UNDEF_ORDER;
}

template<typename Scalar>
Image4D<Scalar>::Image4D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData) {
	Image<Scalar,4>::_dim = dimensions;
	Image<Scalar,4>::_data = originalData;
	Image<Scalar,4>::_datatype = datatype;
	Image<Scalar,4>::_dataorder = dataorder;
	allocateData();
	if(originalData!=NULL)
		copyData(originalData);
}

template<typename Scalar>
Image4D<Scalar>::Image4D(uint dimX, uint dimY, uint dimZ, uint dimW, dtype datatype, dorder dataorder, Scalar* originalData) {
	Image<Scalar,4>::_data = NULL;
	Image<Scalar,4>::_dim = new uint[4];
	Image<Scalar,4>::_dim[0] = dimX; Image<Scalar,4>::_dim[1] = dimY; Image<Scalar,4>::_dim[2] = dimZ; Image<Scalar,4>::_dim[3] = dimW;
	Image<Scalar,4>::_datatype = datatype;
	Image<Scalar,4>::_dataorder = dataorder;
	allocateData();
	if(originalData!=NULL)
		copyData(originalData);
}

template<typename Scalar>
Image4D<Scalar>::Image4D(const Image4D<Scalar>& X) {
	Image<Scalar,4>::_dim = X._dim;
	Image<Scalar,4>::_data = X._data;
	Image<Scalar,4>::_datatype = X._datatype;
	Image<Scalar,4>::_dataorder = X._dataorder;
}

template<typename Scalar>
Image4D<Scalar>::~Image4D() {
	clearData();
	if(Image<Scalar,4>::_dim!=NULL)
		delete [] Image<Scalar,4>::_dim;
	Image<Scalar,4>::_dim=NULL;
}

template<typename Scalar>
bool Image4D<Scalar>::operator==(Image4D<Scalar>& X) {
	bool result = true;
	if(Image<Scalar,4>::_dim != X.Image<Scalar,4>::_dim) {
		if((Image<Scalar,4>::_dim!=NULL) && (X.Image<Scalar,4>::_dim!=NULL)) {
			for(size_t i=0; i<4; i++)
				result &= (Image<Scalar,4>::_dim[i]==X.Image<Scalar,4>::_dim[i]);
		} else {
			result &= false;
		}
	}
	result &= (Image<Scalar,4>::_datatype==X.Image<Scalar,4>::_datatype);
	result &= (Image<Scalar,4>::_dataorder==X.Image<Scalar,4>::_dataorder);
	if(Image<Scalar,4>::_data != X.Image<Scalar,4>::_data) {
		unsigned long fsize_this=1, fsize_X=1;
		if((Image<Scalar,4>::_dim!=NULL) && (X.Image<Scalar,4>::_dim!=NULL)) {
			for(size_t index=0; index<4; index++) {
				fsize_this*=Image<Scalar,4>::_dim[index];
				fsize_X*=X.Image<Scalar,4>::_dim[index];
			}
			if(fsize_this==fsize_X) {
				for(unsigned long i=0; i<fsize_this;i++) {
					result &= (Image<Scalar,4>::_data[i]==X.Image<Scalar,4>::_data[i]);
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
bool Image4D<Scalar>::operator!=(Image4D<Scalar>& X) {
	bool result = true;
	if(Image<Scalar,4>::_dim != X.Image<Scalar,4>::_dim) {
		if((Image<Scalar,4>::_dim!=NULL) && (X.Image<Scalar,4>::_dim!=NULL)) {
			for(size_t i=0; i<4; i++)
				result &= (Image<Scalar,4>::_dim[i]==X.Image<Scalar,4>::_dim[i]);
		} else {
			result &= false;
		}
	}
	result &= (Image<Scalar,4>::_datatype==X.Image<Scalar,4>::_datatype);
	result &= (Image<Scalar,4>::_dataorder==X.Image<Scalar,4>::_dataorder);
	if(Image<Scalar,4>::_data != X.Image<Scalar,4>::_data) {
		unsigned long fsize_this=1, fsize_X=1;
		if((Image<Scalar,4>::_dim!=NULL) && (X.Image<Scalar,4>::_dim!=NULL)) {
			for(size_t index=0; index<4; index++) {
				fsize_this*=Image<Scalar,4>::_dim[index];
				fsize_X*=X.Image<Scalar,4>::_dim[index];
			}
			if(fsize_this==fsize_X) {
				for(unsigned long i=0; i<fsize_this;i++) {
					result &= (Image<Scalar,4>::_data[i]==X.Image<Scalar,4>::_data[i]);
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
void Image4D<Scalar>::setImageValue(Scalar value, uint coordX, uint coordY, uint coordZ, uint coordW) {
	if(Image<Scalar,4>::_dim==NULL)
		return;
	if(Image<Scalar,4>::_data==NULL)
		return;
	unsigned long index = 0;
	if(Image<Scalar,4>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(Image<Scalar,4>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		return;
	}
	if(index >= (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]))
		return;
	Image<Scalar,4>::_data[index] = value;
}

template<typename Scalar>
void Image4D<Scalar>::createImage(void) {
	if(Image<Scalar,4>::_dim==NULL)
		return;
	clearData();
	if(Image<Scalar,4>::_datatype!=NONE)
		allocateData();
}

template<typename Scalar>
void Image4D<Scalar>::createImage(uint* dim, dtype datatype) {
	if((dim==NULL) || (datatype==NONE))
		return;
	clearData();
	if(Image<Scalar,4>::_dim!=NULL)
		delete [] Image<Scalar,4>::_dim;
	Image<Scalar,4>::_dim = new uint[4];
	memcpy(Image<Scalar,4>::_dim, dim, sizeof(uint)*4);
	Image<Scalar,4>::_datatype=datatype;
	allocateData();
}

template<typename Scalar>
void Image4D<Scalar>::clean(void) {
	zeroData();
}

template<typename Scalar>
void Image4D<Scalar>::pad(uint totalPadX, uint totalPadY, uint totalPadZ, uint totalPadW) {
	unsigned long fsize_old=1;
	if(Image<Scalar,4>::_dim==NULL)
		return;
	if(Image<Scalar,4>::_data==NULL)
		return;
	if(Image<Scalar,4>::_datatype==NONE)
		return;
	if(Image<Scalar,4>::_dataorder==UNDEF_ORDER)
		return;
	for(size_t index=0; index<4; index++)
		fsize_old*=Image<Scalar,4>::_dim[index];
	uint* ndim = new uint[4];
	ndim[0]=Image<Scalar,4>::_dim[0]+totalPadX; ndim[1]=Image<Scalar,4>::_dim[1]+totalPadY; ndim[2]=Image<Scalar,4>::_dim[2]+totalPadZ; ndim[3]=Image<Scalar,4>::_dim[3]+totalPadW;
	unsigned long fsize_new = ndim[0]*ndim[1]*ndim[2]*ndim[3];
	uint offsetX = uint(std::floor(totalPadX/2));
	uint offsetY = uint(std::floor(totalPadY/2));
	uint offsetZ = uint(std::floor(totalPadZ/2));
	uint offsetW = uint(std::floor(totalPadW/2));

	Scalar* tbuffer = new Scalar[fsize_new];
	ulong oldAddress=0, newAddress=0;
	for(uint w=0; w<Image<Scalar,4>::_dim[3]; w++) {
		for(uint z=0; z<Image<Scalar,4>::_dim[2]; z++) {
			for(uint y=0; y<Image<Scalar,4>::_dim[1]; y++) {
				for(uint x=0; x<Image<Scalar,4>::_dim[0]; x++) {
					oldAddress = getDataAddress(x,y,z,w);
					newAddress = getDataAddress_field4D(x+offsetX, ndim[0], y+offsetY, ndim[1], z+offsetZ, ndim[2], w+offsetW, ndim[3], Image<Scalar,4>::_dataorder);
					tbuffer[newAddress] = Image<Scalar,4>::_data[oldAddress];
				}
			}
		}
	}
	delete [] Image<Scalar,4>::_data; Image<Scalar,4>::_data=NULL;
	memcpy(Image<Scalar,4>::_dim, ndim, sizeof(uint)*4);
	Image<Scalar,4>::_data = new Scalar[fsize_new];
	memcpy(Image<Scalar,4>::_data, tbuffer, sizeof(Scalar)*fsize_new);
}

template<typename Scalar>
void Image4D<Scalar>::swap(Image4D<Scalar>& X) {
	if((Image<Scalar,4>::_data!=NULL) && (X.Image<Scalar,4>::_data!=NULL)) {
		ulong fsize_this = Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3];
		ulong fsize_X = X.Image<Scalar,4>::_dim[0]*X.Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3];
		Scalar* temp = new Scalar[fsize_this];
		memcpy(temp, Image<Scalar,4>::_data, sizeof(Scalar)*fsize_this);
		clearData();
		Image<Scalar,4>::_data = new Scalar[fsize_X];
		memcpy(Image<Scalar,4>::_data, X.Image<Scalar,4>::_data, sizeof(Scalar)*fsize_X);
		X.clearData();
		X.Image<Scalar,4>::_data = new Scalar[fsize_this];
		memcpy(X.Image<Scalar,4>::_data, temp, sizeof(Scalar)*fsize_this);
		delete [] temp;
	}
	uint tempDim[4];
	if((Image<Scalar,4>::_dim!=NULL) && (X.Image<Scalar,4>::_dim!=NULL)) {
		memcpy((uint*)tempDim, Image<Scalar,4>::_dim, sizeof(uint)*4);
		memcpy(Image<Scalar,4>::_dim, X.Image<Scalar,4>::_dim, sizeof(uint)*4);
		memcpy(X.Image<Scalar,4>::_dim, (uint*)tempDim, sizeof(uint)*4);
	}
	std::swap(Image<Scalar,4>::_datatype, X.Image<Scalar,4>::_datatype);
	std::swap(Image<Scalar,4>::_dataorder, X.Image<Scalar,4>::_dataorder);
}

template<typename Scalar>
void Image4D<Scalar>::swap(const Image4D<Scalar>& X) {

	if(X.Image<Scalar,4>::_dim!=NULL) {
		memcpy(Image<Scalar,4>::_dim, X.Image<Scalar,4>::_dim, sizeof(uint)*4);

		if(Image<Scalar,4>::_data!=NULL) {
			delete [] Image<Scalar,4>::_data;
			Image<Scalar,4>::_data = new Scalar[Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]];
		}
	} else {
		if(Image<Scalar,4>::_dim!=NULL)
			delete [] Image<Scalar,4>::_dim;
		Image<Scalar,4>::_dim=NULL;
	}

	if(X.Image<Scalar,4>::_data!=NULL) {
		memcpy(Image<Scalar,4>::_data, X.Image<Scalar,4>::_data, sizeof(Scalar)*Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]);
	} else {
		if(Image<Scalar,4>::_data!=NULL)
			delete [] Image<Scalar,4>::_data;
		Image<Scalar,4>::_data=NULL;
	}

	Image<Scalar,4>::_datatype = X.Image<Scalar,4>::_datatype;
	Image<Scalar,4>::_dataorder = X.Image<Scalar,4>::_dataorder;
}

template<typename Scalar>
Scalar& Image4D<Scalar>::operator[](unsigned long index) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(index < (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]));

	return Image<Scalar,4>::_data[index];
}

template<typename Scalar>
Scalar& Image4D<Scalar>::operator()(uint coordX, uint coordY, uint coordZ, uint coordW) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(getDataorder()!=UNDEF_ORDER);
	unsigned long index = 0;
	if(Image<Scalar,4>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(Image<Scalar,4>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		index = 0;
	}
	assert(index < (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]));

	return Image<Scalar,4>::_data[index];
}

template<typename Scalar>
Scalar& Image4D<Scalar>::getImageValue(uint coordX, uint coordY, uint coordZ, uint coordW) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(getDataorder()!=UNDEF_ORDER);
	unsigned long index = 0;
	if(Image<Scalar,4>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ, coordW);
	} else if(Image<Scalar,4>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ, coordW);
	} else {
		index = 0;
	}
	assert(index < (Image<Scalar,4>::_dim[0]*Image<Scalar,4>::_dim[1]*Image<Scalar,4>::_dim[2]*Image<Scalar,4>::_dim[3]));

	return Image<Scalar,4>::_data[index];
}

template<typename Scalar>
Scalar Image4D<Scalar>::interpolate(float coordX, float coordY, float coordZ, float coordW, interpolationMode4D mode) {
	switch(mode) {
	case NEAREST_NEIGHBOUR_4D: return interpolateNearestNeighbour(coordX, coordY, coordZ, coordW);
	case QUADRILINEAR: return interpolateQuadrilinear(coordX, coordY, coordZ, coordW);
	}
	return Scalar(0);
}

template<typename Scalar>
Scalar Image4D<Scalar>::interpolateQuadrilinear(float coordX, float coordY, float coordZ, float coordW) {
	uint xf = uint(floor(coordX)); uint yf = uint(floor(coordY)); uint zf = uint(floor(coordZ)); uint wf = uint(floor(coordW));
	float x_offset = (coordX - float(xf)) - 0.5; float y_offset = (coordY - float(yf)) - 0.5; float z_offset = (coordZ - float(zf)) - 0.5; float w_offset = (coordW - float(wf)) - 0.5;
	float pixDistance = .0f, kd=.0f, ld=.0f, md=.0f, nd=.0f;
	float maxdist = sqrt(3*1.5);
	//float interpCumulative=.0f;
	double result=0;
	Scalar kernel[3][3][3][3];
	uint kernel_k=0, kernel_l=0, kernel_m=0, kernel_n=0;
	for(int k=(int(xf)-1); k<(int(xf)+2); k++) {
		for(int l=(int(yf)-1); l<(int(yf)+2); l++) {
			for(int m=(int(zf)-1); m<(int(zf)+2); m++) {
				for(int n=(int(wf)-1); n<(int(wf)+2); n++) {
					//zero-based addresses
					kernel_k=k-(xf-1);
					kernel_l=l-(yf-1);
					kernel_m=m-(zf-1);
					kernel_n=n-(wf-1);
					if((k<0)||(uint(k)>(Image<Scalar,4>::_dim[0]-1)) || (l<0)||(uint(l)>(Image<Scalar,4>::_dim[1]-1)) || (m<0)||(uint(m)>(Image<Scalar,4>::_dim[2]-1)) || (n<0)||(uint(n)>(Image<Scalar,4>::_dim[3]-1))) {
						kernel[kernel_k][kernel_l][kernel_m][kernel_n]=0;
					} else {
						kernel[kernel_k][kernel_l][kernel_m][kernel_n]=Image<Scalar,4>::_data[getDataAddress(k,l,m,n)];
					}
				}
			}
		}
	}
	for(uint k=0;k<3;k++) {
		for(uint l=0;l<3;l++) {
			for(uint m=0;m<3;m++) {
				for(uint n=0;n<3;n++) {
					kd=float(k-1)+x_offset;
					ld=float(l-1)+y_offset;
					md=float(m-1)+z_offset;
					nd=float(n-1)+w_offset;
					pixDistance=sqrt(kd*kd+ld*ld+md*md+nd*nd)/maxdist;
					result+=floor( (1.0f-pixDistance)*kernel[k][l][m][n] );
				}
			}
		}
	}
	return Scalar(floor(result/81.0));
}

template<typename Scalar>
Scalar Image4D<Scalar>::interpolateNearestNeighbour(float coordX, float coordY, float coordZ, float coordW) {
	uint xr = uint(round(coordX)); uint yr = uint(round(coordY)); uint zr = uint(round(coordZ)); uint wr = uint(round(coordW));
	return getImageValue(xr, yr, zr, wr);
}

template<typename Scalar>
void Image4D<Scalar>::clamp(const Scalar& lowValuePtr, const Scalar& highValuePtr) {
	if(Image<Scalar,4>::_dim==NULL)
		return;
	if(Image<Scalar,4>::_data==NULL)
		return;
	unsigned long fsize=1;
	for(size_t index=0; index<4; index++)
		fsize*=Image<Scalar,4>::_dim[index];
	for(ulong i=0; i<fsize; i++) {
		Image<Scalar,4>::_data[i]=clampByValue(Image<Scalar,4>::_data[i], lowValuePtr, highValuePtr);
	}
}

template<typename Scalar>
Scalar Image4D<Scalar>::mean(void) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	unsigned long fsize=1;
	for(size_t index=0; index<4; index++)
		fsize*=Image<Scalar,4>::_dim[index];
	double result=0;
	for(uint i=0; i<fsize; i++)
		result+=double(Image<Scalar,4>::_data[i]);
	result/=double(fsize);
	return Scalar(result);
}

template<typename Scalar>
Scalar Image4D<Scalar>::variance(void) {
	double result = 0.0;
    Scalar _mean = mean();
	unsigned long fsize=1;
	for(size_t index=0; index<4; index++)
		fsize*=Image<Scalar,4>::_dim[index];
    for(uint i=0; i<fsize; i++)
    	result += pow(double(Image<Scalar,4>::_data[i] - _mean),2);
    return Scalar(result/double(fsize));
}

template<typename Scalar>
Scalar Image4D<Scalar>::stdDev(void) {
	return Scalar(sqrt(variance()));
}

template<typename Scalar>
Scalar Image4D<Scalar>::max(void) {
	Scalar result = -Scalar(DBL_MAX);
	unsigned long fsize=1;
	for(size_t index=0; index<4; index++)
		fsize*=Image<Scalar,4>::_dim[index];
	for(uint i=0; i<fsize; i++)
		result = std::max(result, Image<Scalar,4>::_data[i]);
	return result;
}

template<typename Scalar>
Scalar Image4D<Scalar>::min(void) {
	Scalar result = Scalar(DBL_MAX);
	unsigned long fsize=1;
	for(size_t index=0; index<4; index++)
		fsize*=Image<Scalar,4>::_dim[index];
	for(uint i=0; i<fsize; i++)
		result = std::min(result, Image<Scalar,4>::_data[i]);
	return result;
}

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


