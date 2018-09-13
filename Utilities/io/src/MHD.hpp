/**
 * created on Jan 9 2018
 * author: Christian Kehl
 */
#ifndef MHD_HPP_
#define MHD_HPP_

#include <stdlib.h>
#include <cstring>
#include <iostream>
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif

/*
example:
ObjectType = Image
NDims = 2
BinaryData = True
BinaryDataByteOrderMSB = False
CompressedData = False
TransformMatrix = 1 0 0 0 1 0 0 0 1	=> ONLY IN 3D
Offset = 0 0 0						=> ONLY IN 3D
CenterOfRotation = 0 0 0			=> ONLY IN 3D
ElementSpacing = 1 1
DimSize = 100 100
AnatomicalOrientation = ???
ElementType = MET_SHORT
ElementDataFile = segmented.raw
 */

/*
template<typename T>
void ARTReconstruct<T>::SaveVolume(const char * pcFilename)
{
	printf("Saving volume %s\n", pcFilename);
	std::ofstream kMhaFile(pcFilename, std::fstream::binary);
	kMhaFile << "ObjectType = Image\n";
	kMhaFile << "NDims = 3\n";
	kMhaFile << "DimSize = " << m_uiImageWidth << " " << m_uiImageHeight << " " << m_uiImageSlices << "\n";
	if ( sizeof(T) == 4 )
	{
		kMhaFile << "ElementType = MET_FLOAT\n";
	}
	else
	{
		kMhaFile << "ElementType = MET_DOUBLE\n";
	}
	kMhaFile << "ElementSize = " << m_dDelta << " " << m_dDelta << " " << m_dDelta << "\n";
	kMhaFile << "ElementByteOrderMSB = False\n";
	kMhaFile << "ElementDataFile = LOCAL\n";
	kMhaFile.write(reinterpret_cast<const char*>(GetOutput()), m_uiImageWidth*m_uiImageHeight*m_uiImageSlices*sizeof(T));
	kMhaFile.close();
}
 */

class MhdFile {
public:
	inline MhdFile() : _data(NULL), _ndims(0), _dim(NULL), _spacing(NULL), _volumeCenter(NULL), _filename(""), _datatype(NONE) {}
	inline virtual ~MhdFile() {
		if(_data!=NULL) {
			switch(_datatype) {
				case CHAR: {
					delete [] ((char*)_data);
					break;
				}
				case UCHAR: {
					delete [] ((uchar*)_data);
					break;
				}
				case SHORT: {
					delete [] ((short*)_data);
					break;
				}
				case USHORT: {
					delete [] ((ushort*)_data);
					break;
				}
				case INT: {
					delete [] ((int*)_data);
					break;
				}
				case UINT: {
					delete [] ((uint*)_data);
					break;
				}
				case LONG: {
					delete [] ((long*)_data);
					break;
				}
				case ULONG: {
					delete [] ((ulong*)_data);
					break;
				}
				case FLOAT: {
					delete [] ((float*)_data);
					break;
				}
				case DOUBLE: {
					delete [] ((double*)_data);
					break;
				}
				default: {
					break;
				}
			}
		}
		if(_dim!=NULL)
			delete [] _dim;
		if(_spacing!=NULL)
			delete [] _spacing;
		if(_volumeCenter!=NULL)
			delete [] _volumeCenter;
	}
	
	void setFilename(std::string filename);
	inline void setDatatype(dtype datatype) { _datatype = datatype; }
	inline void setNumberOfDimensions(int ndims) { _ndims = ndims; }
	void setDimensions(uint* dimensions);
	void SetDimensions(unsigned int* numpy_input_dimensions, int field_size);
	void setSpacing(float* spacing);
	void SetSpacing(float* numpy_input_spacing, int field_size);
	void SetSpacing(double* numpy_input_spacing, int field_size);
	void setData(void* data, dtype datatype, bool memcopy=true); // needs to be "something" unitary originally
	void setDataAsChar(char* numpy_input_data, long field_size);
	void setDataAsUChar(unsigned char* numpy_input_data, long field_size);
	void setDataAsShort(short* numpy_input_data, long field_size);
	void setDataAsUShort(unsigned short* numpy_input_data, long field_size);
	void setDataAsInt(int* numpy_input_data, long field_size);
	void setDataAsUInt(unsigned int* numpy_input_data, long field_size);
	void setDataAsLong(long* numpy_input_data, long field_size);
	void setDataAsULong(unsigned long* numpy_input_data, long field_size);
	void setDataAsFloat(float* numpy_input_data, long field_size);
	void setDataAsDouble(double* numpy_input_data, long field_size);
	//void setMinMax(void* minmax_structure);
	//void setMinMax(char minx, char miny, char minz, char maxx, char maxy, char maxz);
	//void setMinMax(uchar minx, uchar miny, uchar minz, uchar maxx, uchar maxy, uchar maxz);
	//void setMinMax(short minx, short miny, short minz, short maxx, short maxy, short maxz);
	//void setMinMax(ushort minx, ushort miny, ushort minz, ushort maxx, ushort maxy, ushort maxz);
	//void setMinMax(int minx, int miny, int minz, int maxx, int maxy, int maxz);
	//void setMinMax(uint minx, uint miny, uint minz, uint maxx, uint maxy, uint maxz);
	//void setMinMax(long minx, long miny, long minz, long maxx, long maxy, long maxz);
	//void setMinMax(ulong minx, ulong miny, ulong minz, ulong maxx, ulong maxy, ulong maxz);
	//void setMinMax(float minx, float miny, float minz, float maxx, float maxy, float maxz);
	//void setMinMax(double minx, double miny, double minz, double maxx, double maxy, double maxz);

	inline void* getData() { return _data; }
	void getDataAsChar(char* numpy_output_data, long field_size);
	void getDataAsUChar(unsigned char* numpy_output_data, long field_size);
	void getDataAsShort(short* numpy_output_data, long field_size);
	void getDataAsUShort(unsigned short* numpy_output_data, long field_size);
	void getDataAsInt(int* numpy_output_data, long field_size);
	void getDataAsUInt(unsigned int* numpy_output_data, long field_size);
	void getDataAsLong(long* numpy_output_data, long field_size);
	void getDataAsULong(unsigned long* numpy_output_data, long field_size);
	void getDataAsFloat(float* numpy_output_data, long field_size);
	void getDataAsDouble(double* numpy_output_data, long field_size);
	inline uint* getDimensions() { return _dim; }
	void GetDimensions(unsigned int* numpy_output_dimensions, int field_size);
	inline float* getSpacing() { return _spacing; }
	void GetSpacing(float* numpy_output_spacing, int field_size);
	inline dtype getDatatype() { return _datatype; }
	inline int getNumberOfDimensions() { return _ndims; }
	inline std::string getFileName() { return _filename; }

	void readFile(void);
	void writeFile(void);

protected:
	inline void clearData(void) {
		if(_data==NULL)
			return;
		else {
			switch(_datatype) {
				case CHAR: {
					delete [] ((char*)_data);
					break;
				}
				case UCHAR: {
					delete [] ((uchar*)_data);
					break;
				}
				case SHORT: {
					delete [] ((short*)_data);
					break;
				}
				case USHORT: {
					delete [] ((ushort*)_data);
					break;
				}
				case INT: {
					delete [] ((int*)_data);
					break;
				}
				case UINT: {
					delete [] ((uint*)_data);
					break;
				}
				case LONG: {
					delete [] ((long*)_data);
					break;
				}
				case ULONG: {
					delete [] ((ulong*)_data);
					break;
				}
				case FLOAT: {
					delete [] ((float*)_data);
					break;
				}
				case DOUBLE: {
					delete [] ((double*)_data);
					break;
				}
				default: {
					break;
				}
			}
		}
	}
	void* _data;
	int _ndims;				// auto-deducted when setting the dimensions
	unsigned int * _dim;
	float* _spacing;
	float* _volumeCenter;	// always 3D coord; computed automatically
	// void* _minmax; // will be a pointer to the appropriate MinMax<T> structure - not needed for MHD
	std::string _filename;
	dtype _datatype;
};




#endif /* MHD_HPP_ */
