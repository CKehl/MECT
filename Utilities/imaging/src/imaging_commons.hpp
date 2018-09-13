/**
 * Created on 01 May 2018
 * author: christian
 */

#ifndef IMAGING_COMMONS_HPP_
#define IMAGING_COMMONS_HPP_

#include "image.hpp"
#include "image2D.hpp"
#include "image3D.hpp"
#include "image4D.hpp"
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

double EbMultix_calibrated[] = {20.726,21.768,22.809,23.851,24.893,25.934,26.976,28.017,29.059,30.101,31.142,32.184,33.226,34.267,
		35.309,36.35,37.392,38.434,39.475,40.517,41.559,42.6,43.642,44.684,45.725,46.767,47.808,48.85,49.892,50.933,51.975,53.017,
		54.058,55.1,56.141,57.183,58.225,59.266,60.308,61.35,62.391,63.433,64.474,65.516,66.558,67.599,68.641,69.683,70.724,71.766,
		72.807,73.849,74.891,75.932,76.974,78.016,79.057,80.099,81.14,82.182,83.224,84.265,85.307,86.349,87.39,88.432,89.473,90.515,
		91.557,92.598,93.64,94.682,95.723,96.765,97.806,98.848,99.89,100.931,101.973,103.015,104.056,105.098,106.139,107.181,108.223,
		109.264,110.306,111.348,112.389,113.431,114.472,115.514,116.556,117.597,118.639,119.681,120.722,121.764,122.805,123.847,124.889,
		125.93,126.972,128.014,129.055,130.097,131.138,132.18,133.222,134.263,135.305,136.347,137.388,138.43,139.471,140.513,141.555,142.596,
		143.638,144.68,145.721,146.763,147.804,148.846,149.888,150.929,151.971,153.013};


// dissimilarity measure between pixels
template<typename Scalar>
inline std::vector< Image3D<Scalar>* > toArray(Image4D<Scalar>& data) {
	uint arrayLength = data.getDimension(3);
	uint dim0 = data.getDimension(0), dim1 = data.getDimension(1), dim2 = data.getDimension(2);

	std::vector< Image3D<Scalar>* > result;
	for(uint i=0; i<arrayLength; i++) {
		Image3D<Scalar>* container = new Image3D<Scalar>();
		container->setDimensions(dim0,dim1,dim2);
		container->setDataorder(data.getDataorder());
		container->setDatatype(data.getDatatype());
		container->createImage();
		result.push_back(container);
	}

	#pragma omp parallel for shared(result,data)
	for(unsigned long i=0; i<(dim0*dim1*dim2*arrayLength); i++) {
		uint channel=0; uint z=0; uint y=0; uint x=0;
		getIndicesFromAddress_field4D(i,x,dim0,y,dim1,z,dim2,channel,arrayLength,data.getDataorder());
		result.at(channel)->setImageValue(data.getImageValue(x,y,z,channel),x,y,z);
	}

	return result;
}

template<typename Scalar>
inline std::vector< Image2D<Scalar>* > toArray(Image3D<Scalar>& data) {
	uint arrayLength = data.getDimension(2);
	uint dim0 = data.getDimension(0), dim1 = data.getDimension(1);

	std::vector< Image2D<Scalar>* > result;
	for(uint i=0; i<arrayLength; i++) {
		Image2D<Scalar>* container = new Image2D<Scalar>();
		container->setDimensions(dim0,dim1);
		container->setDataorder(data.getDataorder());
		container->setDatatype(data.getDatatype());
		container->createImage();
		result.push_back(container);
	}

	#pragma omp parallel for shared(result,data)
	for(unsigned long i=0; i<(dim0*dim1*arrayLength); i++) {
		uint channel=0; uint y=0; uint x=0;
		getIndicesFromAddress_field3D(i,x,dim0,y,dim1,channel,arrayLength,data.getDataorder());
		result.at(channel)->setImageValue(data.getImageValue(x,y,channel),x,y);
	}

	return result;
}

#endif /* IMAGING_COMMONS_HPP_ */
