/**
 * Created on 25 Apr 2018
 * author: christian
 */

#ifndef IMAGE2D_UTILS_HPP_
#define IMAGE2D_UTILS_HPP_

#include "image.hpp"
#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif
#include <iostream>
#include <vector>
#include <string.h>
#include <cmath>
#include <omp.h>

typedef enum _NEIGHBOURHOOD2D {
	NBR_5   = 0,
	NBR_9  = 1,
	NBR_9W = 2,
} NEIGHBOURHOOD2D;







// dissimilarity measure between pixels
template<typename Scalar>
inline Scalar SSD_imageStack(Image2D<Scalar>& data, int x1, int x2) {
	float result = 0;
	for(unsigned int c = 0; c < data.getDimension(1); c++) {
		result += (data.getImageValue(x1, c)-data.getImageValue(x2, c))*(data.getImageValue(x1, c)-data.getImageValue(x2, c));
	}
	return std::sqrt(result);
}

#endif /* IMAGE2D_UTILS_HPP_ */
