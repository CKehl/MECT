/**
 * Created on 18 Jan 2018
 * author: christian
 */
#include "image3D.hpp"
#include <cmath>
#include <float.h>
#include <stdint.h>
#include <cassert>

#include <omp.h>

template<typename Scalar>
Image3D<Scalar>::Image3D() {
	Image<Scalar,3>::_dim = NULL;
	Image<Scalar,3>::_data = NULL;
	Image<Scalar,3>::_datatype = NONE;
	Image<Scalar,3>::_dataorder = UNDEF_ORDER;
}

template<typename Scalar>
Image3D<Scalar>::Image3D(uint* dimensions, dtype datatype, dorder dataorder, Scalar* originalData) {
	Image<Scalar,3>::_dim = dimensions;
	Image<Scalar,3>::_data = originalData;
	Image<Scalar,3>::_datatype = datatype;
	Image<Scalar,3>::_dataorder = dataorder;
}

template<typename Scalar>
Image3D<Scalar>::Image3D(uint dimX, uint dimY, uint dimZ, dtype datatype, dorder dataorder, Scalar* originalData) {
	Image<Scalar,3>::_data = NULL;
	Image<Scalar,3>::_dim = new uint[3];
	Image<Scalar,3>::_dim[0] = dimX; Image<Scalar,3>::_dim[1] = dimY; Image<Scalar,3>::_dim[2] = dimZ;
	Image<Scalar,3>::_datatype = datatype;
	Image<Scalar,3>::_dataorder = dataorder;
	allocateData();
	if(originalData!=NULL)
		copyData(originalData);
}

template<typename Scalar>
Image3D<Scalar>::Image3D(const Image3D<Scalar>& X) {
	Image<Scalar,3>::_dim = X._dim;
	Image<Scalar,3>::_data = X._data;
	Image<Scalar,3>::_datatype = X._datatype;
	Image<Scalar,3>::_dataorder = X._dataorder;
}

template<typename Scalar>
Image3D<Scalar>::~Image3D() {
	clearData();
	if(Image<Scalar,3>::_dim!=NULL)
		delete [] Image<Scalar,3>::_dim;
	Image<Scalar,3>::_dim=NULL;
}

template<typename Scalar>
bool Image3D<Scalar>::operator==(Image3D<Scalar>& X) {
	bool result = true;
	if(Image<Scalar,3>::_dim != X._dim) {
		if((Image<Scalar,3>::_dim!=NULL) && (X.Image<Scalar,3>::_dim!=NULL)) {
			for(size_t i=0; i<3; i++)
				result &= (Image<Scalar,3>::_dim[i]==X.Image<Scalar,3>::_dim[i]);
		} else {
			result &= false;
		}
	}
	result &= (Image<Scalar,3>::_datatype==X.Image<Scalar,3>::_datatype);
	result &= (Image<Scalar,3>::_dataorder==X.Image<Scalar,3>::_dataorder);
	if(Image<Scalar,3>::_data != X.Image<Scalar,3>::_data) {
		unsigned long fsize_this=1, fsize_X=1;
		if((Image<Scalar,3>::_dim!=NULL) && (X.Image<Scalar,3>::_dim!=NULL)) {
			for(size_t index=0; index<3; index++) {
				fsize_this*=Image<Scalar,3>::_dim[index];
				fsize_X*=X.Image<Scalar,3>::_dim[index];
			}
			if(fsize_this==fsize_X) {
				for(unsigned long i=0; i<fsize_this;i++) {
					result &= (Image<Scalar,3>::_data[i]==X.Image<Scalar,3>::_data[i]);
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
bool Image3D<Scalar>::operator!=(Image3D<Scalar>& X) {
	bool result = true;
	if(Image<Scalar,3>::_dim != X.Image<Scalar,3>::_dim) {
		if((Image<Scalar,3>::_dim!=NULL) && (X.Image<Scalar,3>::_dim!=NULL)) {
			for(size_t i=0; i<3; i++)
				result &= (Image<Scalar,3>::_dim[i]==X.Image<Scalar,3>::_dim[i]);
		} else {
			result &= false;
		}
	}
	result &= (Image<Scalar,3>::_datatype==X.Image<Scalar,3>::_datatype);
	result &= (Image<Scalar,3>::_dataorder==X._dataorder);
	if(Image<Scalar,3>::_data != X._data) {
		unsigned long fsize_this=1, fsize_X=1;
		if((Image<Scalar,3>::_dim!=NULL) && (X._dim!=NULL)) {
			for(size_t index=0; index<3; index++) {
				fsize_this*=Image<Scalar,3>::_dim[index];
				fsize_X*=X._dim[index];
			}
			if(fsize_this==fsize_X) {
				for(unsigned long i=0; i<fsize_this;i++) {
					result &= (Image<Scalar,3>::_data[i]==X._data[i]);
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
void Image3D<Scalar>::setImageValue(Scalar value, uint coordX, uint coordY, uint coordZ) {
	if(Image<Scalar,3>::_dim==NULL)
		return;
	if(Image<Scalar,3>::_data==NULL)
		return;
	unsigned long index = 0;
	if(Image<Scalar,3>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(Image<Scalar,3>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		return;
	}
	if(index >= (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]))
		return;
	Image<Scalar,3>::_data[index] = value;
}

template<typename Scalar>
void Image3D<Scalar>::createImage(void) {
	if(Image<Scalar,3>::_dim==NULL)
		return;
	clearData();
	if(Image<Scalar,3>::_datatype!=NONE)
		allocateData();
}

template<typename Scalar>
void Image3D<Scalar>::createImage(uint* dim, dtype datatype) {
	if((dim==NULL) || (datatype==NONE))
		return;
	clearData();
	if(Image<Scalar,3>::_dim!=NULL)
		delete [] Image<Scalar,3>::_dim;
	Image<Scalar,3>::_dim = new uint[3];
	memcpy(Image<Scalar,3>::_dim, dim, sizeof(uint)*3);
	Image<Scalar,3>::_datatype=datatype;
	allocateData();
}

template<typename Scalar>
void Image3D<Scalar>::clean(void) {
	zeroData();
}

template<typename Scalar>
void Image3D<Scalar>::pad(uint totalPadX, uint totalPadY, uint totalPadZ) {
	unsigned long fsize_old=1;
	if(Image<Scalar,3>::_dim==NULL)
		return;
	if(Image<Scalar,3>::_data==NULL)
		return;
	if(Image<Scalar,3>::_datatype==NONE)
		return;
	if(Image<Scalar,3>::_dataorder==UNDEF_ORDER)
		return;
	for(size_t index=0; index<3; index++)
		fsize_old*=Image<Scalar,3>::_dim[index];
	uint* ndim = new uint[3];
	ndim[0]=Image<Scalar,3>::_dim[0]+totalPadX; ndim[1]=Image<Scalar,3>::_dim[1]+totalPadY; ndim[2]=Image<Scalar,3>::_dim[2]+totalPadZ;
	unsigned long fsize_new = ndim[0]*ndim[1]*ndim[2];
	uint offsetX = uint(std::floor(totalPadX/2));
	uint offsetY = uint(std::floor(totalPadY/2));
	uint offsetZ = uint(std::floor(totalPadZ/2));

	Scalar* tbuffer = new Scalar[fsize_new];
	ulong oldAddress=0, newAddress=0;
	for(uint z=0; z<Image<Scalar,3>::_dim[2]; z++) {
		for(uint y=0; y<Image<Scalar,3>::_dim[1]; y++) {
			for(uint x=0; x<Image<Scalar,3>::_dim[0]; x++) {
				oldAddress = getDataAddress(x,y,z);
				newAddress = getDataAddress_field3D(x+offsetX, ndim[0], y+offsetY, ndim[1], z+offsetZ, ndim[2], Image<Scalar,3>::_dataorder);
				tbuffer[newAddress] = Image<Scalar,3>::_data[oldAddress];
			}
		}
	}
	delete [] Image<Scalar,3>::_data; Image<Scalar,3>::_data=NULL;
	memcpy(Image<Scalar,3>::_dim, ndim, sizeof(uint)*3);
	Image<Scalar,3>::_data = new Scalar[fsize_new];
	memcpy(Image<Scalar,3>::_data, tbuffer, sizeof(Scalar)*fsize_new);
}

template<typename Scalar>
void Image3D<Scalar>::clip(size_t origin0, size_t origin1, size_t origin2, size_t dim0, size_t dim1, size_t dim2) {
	unsigned long fsize_old=1;
	if(Image<Scalar,3>::_dim==NULL)
		return;
	if(Image<Scalar,3>::_data==NULL)
		return;
	if(Image<Scalar,3>::_datatype==NONE)
		return;
	if(Image<Scalar,3>::_dataorder==UNDEF_ORDER)
		return;
	for(size_t index=0; index<3; index++)
		fsize_old*=Image<Scalar,3>::_dim[index];
	uint* ndim = new uint[3];
	ndim[0]=dim0; ndim[1]=dim1; ndim[2]=dim2;
	unsigned long fsize_new = ndim[0]*ndim[1]*ndim[2];

	Scalar* tbuffer = new Scalar[fsize_new];
	ulong oldAddress=0, newAddress=0;
	for(uint z=0; z<dim2; z++) {
		for(uint y=0; y<dim1; y++) {
			for(uint x=0; x<dim0; x++) {
				oldAddress = getDataAddress(origin0+x,origin1+y,origin2+z);
				newAddress = getDataAddress_field3D(x, ndim[0], y, ndim[1], z, ndim[2], Image<Scalar,3>::_dataorder);
				tbuffer[newAddress] = Image<Scalar,3>::_data[oldAddress];
			}
		}
	}
	delete [] Image<Scalar,3>::_data; Image<Scalar,3>::_data=NULL;
	memcpy(Image<Scalar,3>::_dim, ndim, sizeof(uint)*3);
	Image<Scalar,3>::_data = new Scalar[fsize_new];
	memcpy(Image<Scalar,3>::_data, tbuffer, sizeof(Scalar)*fsize_new);
}

template<typename Scalar>
void Image3D<Scalar>::swap(Image3D<Scalar>& X) {
	if((Image<Scalar,3>::_data!=NULL) && (X.Image<Scalar,3>::_data!=NULL)) {
		ulong fsize_this = Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2];
		ulong fsize_X = X.Image<Scalar,3>::_dim[0]*X.Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2];
		Scalar* temp = new Scalar[fsize_this];
		memcpy(temp, Image<Scalar,3>::_data, sizeof(Scalar)*fsize_this);
		clearData();
		Image<Scalar,3>::_data = new Scalar[fsize_X];
		memcpy(Image<Scalar,3>::_data, X.Image<Scalar,3>::_data, sizeof(Scalar)*fsize_X);
		X.clearData();
		X.Image<Scalar,3>::_data = new Scalar[fsize_this];
		memcpy(X.Image<Scalar,3>::_data, temp, sizeof(Scalar)*fsize_this);
		delete [] temp;
	}
	uint tempDim[3];
	if((Image<Scalar,3>::_dim!=NULL) && (X.Image<Scalar,3>::_dim!=NULL)) {
		memcpy((uint*)tempDim, Image<Scalar,3>::_dim, sizeof(uint)*3);
		memcpy(Image<Scalar,3>::_dim, X.Image<Scalar,3>::_dim, sizeof(uint)*3);
		memcpy(X.Image<Scalar,3>::_dim, (uint*)tempDim, sizeof(uint)*3);
	}
	std::swap(Image<Scalar,3>::_datatype, X._datatype);
	std::swap(Image<Scalar,3>::_dataorder, X._dataorder);
}

template<typename Scalar>
void Image3D<Scalar>::swap(const Image3D<Scalar>& X) {

	if(X.Image<Scalar,3>::_dim!=NULL) {
		memcpy(Image<Scalar,3>::_dim, X.Image<Scalar,3>::_dim, sizeof(uint)*3);

		if(Image<Scalar,3>::_data!=NULL) {
			delete [] Image<Scalar,3>::_data;
			Image<Scalar,3>::_data = new Scalar[Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]];
		}
	} else {
		if(Image<Scalar,3>::_dim!=NULL)
			delete [] Image<Scalar,3>::_dim;
		Image<Scalar,3>::_dim=NULL;
	}

	if(X.Image<Scalar,3>::_data!=NULL) {
		memcpy(Image<Scalar,3>::_data, X.Image<Scalar,3>::_data, sizeof(Scalar)*Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]);
	} else {
		if(Image<Scalar,3>::_data!=NULL)
			delete [] Image<Scalar,3>::_data;
		Image<Scalar,3>::_data=NULL;
	}

	Image<Scalar,3>::_datatype = X.Image<Scalar,3>::_datatype;
	Image<Scalar,3>::_dataorder = X.Image<Scalar,3>::_dataorder;
}

template<typename Scalar>
Scalar& Image3D<Scalar>::operator[](unsigned long index) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(index < (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]));

	return Image<Scalar,3>::_data[index];
}

template<typename Scalar>
Scalar& Image3D<Scalar>::operator()(uint coordX, uint coordY, uint coordZ) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(getDataorder()!=UNDEF_ORDER);
	unsigned long index = 0;
	if(Image<Scalar,3>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(Image<Scalar,3>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		index = 0;
	}
	assert(index < (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]));

	return Image<Scalar,3>::_data[index];
}

template<typename Scalar>
Scalar& Image3D<Scalar>::getImageValue(uint coordX, uint coordY, uint coordZ) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	assert(getDataorder()!=UNDEF_ORDER);
	unsigned long index = 0;
	if(Image<Scalar,3>::_dataorder == ROW_MAJOR) {
		index = getDataAddress_rowMajor(coordX, coordY, coordZ);
	} else if(Image<Scalar,3>::_dataorder == COLUMN_MAJOR) {
		index = getDataAddress_columnMajor(coordX, coordY, coordZ);
	} else {
		index = 0;
	}
	assert(index < (Image<Scalar,3>::_dim[0]*Image<Scalar,3>::_dim[1]*Image<Scalar,3>::_dim[2]));

	return Image<Scalar,3>::_data[index];
}

template<typename Scalar>
Scalar Image3D<Scalar>::interpolate(float coordX, float coordY, float coordZ, interpolationMode3D mode) {
	switch(mode) {
	case NEAREST_NEIGHBOUR_3D: return interpolateNearestNeighbour(coordX, coordY, coordZ);
	case TRILINEAR: return interpolateTrilinear(coordX, coordY, coordZ);
	}
	return Scalar(0);
}

template<typename Scalar>
Scalar Image3D<Scalar>::interpolateTrilinear(float coordX, float coordY, float coordZ) {
	uint xf = uint(floor(coordX)); uint yf = uint(floor(coordY)); uint zf = uint(floor(coordZ));
	float x_offset = (coordX - float(xf)) - 0.5; float y_offset = (coordY - float(yf)) - 0.5; float z_offset = (coordZ - float(zf)) - 0.5;
	float pixDistance = .0f, kd=.0f, ld=.0f, md=.0f;
	float maxdist = sqrt(3*1.5);
	//float interpCumulative=.0f;
	double result=0;
	Scalar kernel[3][3][3];
	uint kernel_k=0, kernel_l=0, kernel_m=0;
	for(int k=(int(xf)-1); k<(int(xf)+2); k++) {
		for(int l=(int(yf)-1); l<(int(yf)+2); l++) {
			for(int m=(int(zf)-1); m<(int(zf)+2); m++) {
				//zero-based addresses
				kernel_k=k-(xf-1);
				kernel_l=l-(yf-1);
				kernel_m=m-(zf-1);
				if((k<0)||(uint(k)>(Image<Scalar,3>::_dim[0]-1)) || (l<0)||(uint(l)>(Image<Scalar,3>::_dim[1]-1)) || (m<0)||(uint(m)>(Image<Scalar,3>::_dim[2]-1))) {
					kernel[kernel_k][kernel_l][kernel_m]=0;
				} else {
					kernel[kernel_k][kernel_l][kernel_m]=Image<Scalar,3>::_data[getDataAddress(k,l,m)];
				}
			}
		}
	}
	for(uint k=0;k<3;k++) {
		for(uint l=0;l<3;l++) {
			for(uint m=0;m<3;m++) {
				kd=float(k-1)+x_offset;
				ld=float(l-1)+y_offset;
				md=float(m-1)+z_offset;
				pixDistance=sqrt(kd*kd+ld*ld+md*md)/maxdist;
				result+=floor( (1.0f-pixDistance)*kernel[k][l][m] );
			}
		}
	}
	return Scalar(floor(result/27.0));
}

template<typename Scalar>
Scalar Image3D<Scalar>::interpolateNearestNeighbour(float coordX, float coordY, float coordZ) {
	uint xr = uint(round(coordX)); uint yr = uint(round(coordY)); uint zr = uint(round(coordZ));
	return getImageValue(xr, yr, zr);
}

template<typename Scalar>
void Image3D<Scalar>::clamp(const Scalar& lowValuePtr, const Scalar& highValuePtr) {
	if(Image<Scalar,3>::_dim==NULL)
		return;
	if(Image<Scalar,3>::_data==NULL)
		return;
	unsigned long fsize=1;
	for(size_t index=0; index<3; index++)
		fsize*=Image<Scalar,3>::_dim[index];
	for(ulong i=0; i<fsize; i++) {
		Image<Scalar,3>::_data[fsize]=clampByValue(Image<Scalar,3>::_data[fsize], lowValuePtr, highValuePtr);
	}
}

template<typename Scalar>
Scalar Image3D<Scalar>::mean(void) {
	assert(getDimensions()!=NULL);
	assert(getData()!=NULL);
	unsigned long fsize=1;
	for(size_t index=0; index<3; index++)
		fsize*=Image<Scalar,3>::_dim[index];
	double result=0;
	for(uint i=0; i<fsize; i++)
		result+=double(Image<Scalar,3>::_data[i]);
	result/=double(fsize);
	return Scalar(result);
}

template<typename Scalar>
Scalar Image3D<Scalar>::variance(void) {
	double result = 0.0;
    Scalar _mean = mean();
	unsigned long fsize=1;
	for(size_t index=0; index<3; index++)
		fsize*=Image<Scalar,3>::_dim[index];
    for(uint i=0; i<fsize; i++)
    	result += pow(double(Image<Scalar,3>::_data[i] - _mean),2);
    return Scalar(result/double(fsize));
}

template<typename Scalar>
Scalar Image3D<Scalar>::stdDev(void) {
	return Scalar(sqrt(variance()));
}

template<typename Scalar>
Scalar Image3D<Scalar>::max(void) {
	Scalar result = -Scalar(DBL_MAX);
	unsigned long fsize=1;
	for(size_t index=0; index<3; index++)
		fsize*=Image<Scalar,3>::_dim[index];
	for(uint i=0; i<fsize; i++)
		result = std::max(result, Image<Scalar,3>::_data[i]);
	return result;
}

template<typename Scalar>
Scalar Image3D<Scalar>::min(void) {
	Scalar result = Scalar(DBL_MAX);
	unsigned long fsize=1;
	for(size_t index=0; index<3; index++)
		fsize*=Image<Scalar,3>::_dim[index];
	for(uint i=0; i<fsize; i++)
		result = std::min(result, Image<Scalar,3>::_data[i]);
	return result;
}

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




