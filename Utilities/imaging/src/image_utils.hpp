#ifndef IMAGE_UTILS_HPP_
#define IMAGE_UTILS_HPP_

#include <iostream>
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif

template<typename Scalar, uint ndims>
inline bool isValid(Image<Scalar, ndims>* data) {
	if(data==NULL) {
		std::cout << "Image pointer NULL." << std::endl;
		return false;
	}
	if(data->getDataorder() == UNDEF_ORDER) {
		std::cout << "Image data order undefined." << std::endl;
		return false;
	}
	if(data->getDatatype() == NONE) {
		std::cout << "Image data type undefined." << std::endl;
		return false;
	}
	if(data->getDimensions()==NULL) {
		std::cout << "Image dimensions undefined." << std::endl;
		return false;
	}
	if(data->getData()==NULL) {
		std::cout << "image data itself undefined." << std::endl;
		return false;
	}
	return true;
}

#endif /* IMAGE_UTILS_HPP_ */
