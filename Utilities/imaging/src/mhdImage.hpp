/**
 * Created on 14 Jan 2018
 * author: christian
 */

#ifndef MHDIMAGE_HPP_
#define MHDIMAGE_HPP_

#include <stdlib.h>
#include <cstring>
#include <cassert>
#include <iostream>
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif
#include "image2D.hpp"
#include "image3D.hpp"
#include "image4D.hpp"

template<typename T>
class MhdImage2D {
public:
	MhdImage2D();
	virtual ~MhdImage2D();

	inline void createNewImage() {

	}
	inline void setFilename(std::string filename) {
		_filename = filename;
	}
	inline void setElementSpacing(float* spacing) {
		_spacing = spacing;
	}
	inline void setElementSpacing(float spacing0, float spacing1) {
		if(_spacing!=NULL)
			delete [] _spacing;
		_spacing = new float[2];
		_spacing[0] = spacing0;
		_spacing[1] = spacing1;
	}
	inline void setDatatype(dtype datatype) {
		if(_data==NULL)
			return;
		_data->setDatatype(datatype);
	}
	inline void setDataorder(dorder dataorder) {
		if(_data==NULL)
			return;
		_data->setDataorder(dataorder);
	}
	inline void setDimensions(unsigned int* dims) {
		if(_data==NULL)
			return;
		_data->setDimensions(dims, false);
	}
	inline void setDimensions(unsigned int dim0, unsigned int dim1) {
		if(_data==NULL)
			return;
		_data->setDimensions(dim0, dim1);
	}
	// comment to 'initialise' function: only does something if the image data is empty - otherwise, leaves the old data to avoid data leakage
	inline void initialise(uint dimX, uint dimY, dtype datatype, dorder dataorder) {
		initialise(dimX, dimY, datatype, dataorder, NULL);
	}
	inline void initialise(uint dimX, uint dimY, dtype datatype, dorder dataorder, T* originalData) {
		if(_data!=NULL)
			return;
		if(originalData==NULL) {
			_data = new Image2D<T>(dimX, dimY,datatype, dataorder);
			_data->createImage();
			_externalData = false;
		} else {
			_data = new Image2D<T>(dimX, dimY,datatype, dataorder, originalData);
			_externalData = true;
		}
	}
	inline void setImage(Image2D<T>* data) {
		_data = data;
		_externalData = true;
	}
	inline void setImageData(T* data) {
		if(_data==NULL)
			return;
		_data->setData(data);
		_externalData = true;
	}

	inline std::string getFilename(void) {
		return _filename;
	}
	inline float* getSpacing(void) {
		return _spacing;
	}
	inline Image2D<T>* getImage(void) {
		return _data;
	}

	// read/write functions //
	void read(void);
	void write(void);

	// Wrapping-exclusive functions //
	inline void setImageDataByNumpy(T* image_data_in, unsigned long field_size) {

	}
	inline void setImageDimsByNumpy(unsigned int* dims_in, unsigned int field_size) {

	}
	inline void getImageDataToNumpy(T* image_data_out, unsigned long field_size) {

	}
	inline void getImageDimsToNumpy(unsigned int* dims_out, unsigned int field_size) {

	}

protected:
	Image2D<T>* _data;
	float* _spacing;
	std::string _filename;
	bool _externalData;
};

/*
template class MhdImage2D<char>;
template class MhdImage2D<unsigned char>;
template class MhdImage2D<short>;
template class MhdImage2D<unsigned short>;
template class MhdImage2D<int>;
template class MhdImage2D<unsigned int>;
template class MhdImage2D<long>;
template class MhdImage2D<unsigned long>;
template class MhdImage2D<float>;
template class MhdImage2D<double>;
*/

template<typename T>
class MhdImage3D {
public:
	MhdImage3D();
	virtual ~MhdImage3D();

	// setter- and initializer functions
	inline void createNewImage() {

	}
	inline void setFilename(std::string filename) {
		_filename = filename;
	}
	inline void setElementSpacing(float* spacing) {
		if(_spacing!=NULL)
			delete [] _spacing;
		_spacing = new float[3];
		memcpy(_spacing, spacing, 3*sizeof(float));
	}
	inline void setElementSpacing(float spacing0, float spacing1, float spacing2) {
		if(_spacing!=NULL)
			delete [] _spacing;
		_spacing = new float[3];
		_spacing[0] = spacing0;
		_spacing[1] = spacing1;
		_spacing[2] = spacing2;
	}
	inline void setImageCenter(float* center) {
		if(_volumeCenter!=NULL)
			delete [] _volumeCenter;
		_volumeCenter = new float[3];
		memcpy(_volumeCenter, center, 3*sizeof(float));
	}
	inline void setImageCenter(float center0, float center1, float center2) {
		if(_volumeCenter!=NULL)
			delete [] _volumeCenter;
		_volumeCenter = new float[3];
		_volumeCenter[0] = center0;
		_volumeCenter[1] = center1;
		_volumeCenter[2] = center2;
	}
	inline void setDatatype(dtype datatype) {
		if(_data==NULL)
			return;
		_data->setDatatype(datatype);
	}
	inline void setDataorder(dorder dataorder) {
		if(_data==NULL)
			return;
		_data->setDataorder(dataorder);
	}
	inline void setDimensions(unsigned int* dims) {
		if(_data==NULL)
			return;
		_data->setDimensions(dims, false);
	}
	inline void setDimensions(unsigned int dim0, unsigned int dim1, unsigned int dim2) {
		if(_data==NULL)
			return;
		_data->setDimensions(dim0, dim1, dim2);
	}
	// comment to 'initialise' function: only does something if the image data is empty - otherwise, leaves the old data to avoid data leakage
	inline void initialise(uint dimX, uint dimY, uint dimZ, dtype datatype, dorder dataorder) {
		initialise(dimX, dimY, dimZ, datatype, dataorder, NULL);
	}
	inline void initialise(uint dimX, uint dimY, uint dimZ, dtype datatype, dorder dataorder, T* originalData) {
		if(_data!=NULL)
			return;
		if(originalData==NULL) {
			_data = new Image3D<T>(dimX, dimY, dimZ, datatype, dataorder);
			_data->createImage();
			_externalData = false;
		} else {
			_data = new Image3D<T>(dimX, dimY, dimZ, datatype, dataorder, originalData);
			_externalData = true;
		}
	}
	inline void setImage(Image3D<T>* data) {
		_data = data;
		_externalData = true;
	}
	inline void setImageData(T* data) {
		if(_data==NULL)
			return;
		_data->setData(data);
		_externalData = true;
	}

	// getter functions //
	inline std::string getFilename(void) {
		return _filename;
	}
	inline float* getSpacing(void) {
		return _spacing;
	}
	inline float* getImageCenter(void) {
		return _volumeCenter;
	}
	inline Image3D<T>* getImage(void) {
		return _data;
	}

	// read/write functions //
	void read(void);
	void write(void);

	// Wrapping-exclusive functions //
	inline void setImageDataByNumpy(T* image_data_in, unsigned long field_size) {

	}
	inline void setImageDimsByNumpy(unsigned int* dims_in, unsigned int field_size) {

	}
	inline void getImageDataToNumpy(T* image_data_out, unsigned long field_size) {

	}
	inline void getImageDimsToNumpy(unsigned int* dims_out, unsigned int field_size) {

	}

protected:
	Image3D<T>* _data;
	float* _spacing;
	float* _volumeCenter;	// always 3D coord; computed automatically
	std::string _filename;
	bool _externalData;
};

/*
template class MhdImage3D<char>;
template class MhdImage3D<unsigned char>;
template class MhdImage3D<short>;
template class MhdImage3D<unsigned short>;
template class MhdImage3D<int>;
template class MhdImage3D<unsigned int>;
template class MhdImage3D<long>;
template class MhdImage3D<unsigned long>;
template class MhdImage3D<float>;
template class MhdImage3D<double>;
*/

template<typename T>
class MhdImage4D {
public:
	MhdImage4D();
	virtual ~MhdImage4D();

	inline void createNewImage() {

	}
	inline void setFilename(std::string filename) {
		_filename = filename;
	}
	inline void setElementSpacing(float* spacing) {
		if(_spacing!=NULL)
			delete [] _spacing;
		_spacing = new float[4];
		memcpy(_spacing, spacing, 4*sizeof(float));
	}
	inline void setElementSpacing(float spacing0, float spacing1, float spacing2, float spacing3) {
		if(_spacing!=NULL)
			delete [] _spacing;
		_spacing = new float[4];
		_spacing[0] = spacing0;
		_spacing[1] = spacing1;
		_spacing[2] = spacing2;
		_spacing[3] = spacing3;
	}
	inline void setImageCenter(float* center) {
		if(_volumeCenter!=NULL)
			delete [] _volumeCenter;
		_volumeCenter = new float[3];
		memcpy(_volumeCenter, center, 3*sizeof(float));
	}
	inline void setImageCenter(float center0, float center1, float center2) {
		if(_volumeCenter!=NULL)
			delete [] _volumeCenter;
		_volumeCenter = new float[3];
		_volumeCenter[0] = center0;
		_volumeCenter[1] = center1;
		_volumeCenter[2] = center2;
	}
	inline void setDatatype(dtype datatype) {
		if(_data==NULL)
			return;
		_data->setDatatype(datatype);
	}
	inline void setDataorder(dorder dataorder) {
		if(_data==NULL)
			return;
		_data->setDataorder(dataorder);
	}
	inline void setDimensions(unsigned int* dims) {
		if(_data==NULL)
			return;
		_data->setDimensions(dims, false);
	}
	inline void setDimensions(unsigned int dim0, unsigned int dim1, unsigned int dim2, unsigned int dim3) {
		if(_data==NULL)
			return;
		_data->setDimensions(dim0, dim1, dim2, dim3);
	}
	// comment to 'initialise' function: only does something if the image data is empty - otherwise, leaves the old data to avoid data leakage
	inline void initialise(uint dimX, uint dimY, uint dimZ, uint dimW, dtype datatype, dorder dataorder) {
		initialise(dimX, dimY, dimZ, dimW, datatype, dataorder, NULL);
	}
	inline void initialise(uint dimX, uint dimY, uint dimZ, uint dimW, dtype datatype, dorder dataorder, T* originalData) {
		if(_data!=NULL)
			return;
		if(originalData==NULL) {
			_data = new Image4D<T>(dimX, dimY, dimZ, dimW, datatype, dataorder);
			_data->createImage();
			_externalData = false;
		} else {
			_data = new Image4D<T>(dimX, dimY, dimZ, dimW, datatype, dataorder, originalData);
			_externalData = true;
		}
	}
	inline void setImage(Image4D<T>* data) {
		_data = data;
		_externalData = true;
	}
	inline void setImageData(T* data) {
		if(_data==NULL)
			return;
		_data->setData(data);
		_externalData = true;
	}

	inline std::string getFilename(void) {
		return _filename;
	}
	inline float* getSpacing(void) {
		return _spacing;
	}
	inline float* getImageCenter(void) {
		return _volumeCenter;
	}
	inline Image4D<T>* getImage(void) {
		return _data;
	}

	// read/write functions //
	void read(void);
	void write(void);

	// Wrapping-exclusive functions //
	inline void setImageDataByNumpy(T* image_data_in, unsigned long field_size) {

	}
	inline void setImageDimsByNumpy(unsigned int* dims_in, unsigned int field_size) {

	}
	inline void getImageDataToNumpy(T* image_data_out, unsigned long field_size) {

	}
	inline void getImageDimsToNumpy(unsigned int* dims_out, unsigned int field_size) {

	}

protected:
	Image4D<T>* _data;
	float* _spacing;
	float* _volumeCenter;	// always 3D coord; computed automatically
	std::string _filename;
	bool _externalData;
};

/*
template class MhdImage4D<char>;
template class MhdImage4D<unsigned char>;
template class MhdImage4D<short>;
template class MhdImage4D<unsigned short>;
template class MhdImage4D<int>;
template class MhdImage4D<unsigned int>;
template class MhdImage4D<long>;
template class MhdImage4D<unsigned long>;
template class MhdImage4D<float>;
template class MhdImage4D<double>;
*/





#endif /* MHDIMAGE_HPP_ */
