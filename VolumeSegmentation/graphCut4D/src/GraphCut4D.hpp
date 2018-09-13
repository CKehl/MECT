/*
 * GraphCut4D.h
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#ifndef GRAPHCUT4D_H_
#define GRAPHCUT4D_H_

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <vector>
#ifndef _CXX_COMPILE_
#include "../../../Utilities/common/src/std_typedefs.h"
#include "../../../Utilities/imaging/src/image3D.hpp"
#include "../../../Utilities/imaging/src/image4D.hpp"
#include "../../../Utilities/graphs/src/graphDisjointForest.hpp"
#include "../../../Utilities/graphs/src/graphEdge.hpp"
#else
#include "Utilities/common/src/std_typedefs.h"
#include "Utilities/imaging/src/image3D.hpp"
#include "Utilities/imaging/src/image4D.hpp"
#include "Utilities/graphs/src/graphDisjointForest.hpp"
#include "Utilities/graphs/src/graphEdge.hpp"
#endif

/* make filters */
#define MAKE_FILTER(name, fun)                                \
std::vector<float> make_ ## name (float sigma) {       	      \
  sigma = std::max(sigma, 0.01F);			                  \
  int len = (int)ceil(sigma * WIDTH) + 1;                     \
  std::vector<float> mask(len);                               \
  for (int i = 0; i < len; i++) {                             \
    mask[i] = fun;                                            \
  }                                                           \
  return mask;                                                \
}

void normalizeMask(std::vector<float> &mask);

typedef enum _graph_edge_weight {
	SSD               = 0,
	SAD               = 1,
	ENTROPY_LOSS      = 2,
	INFORMATION_LOSS  = 3,
	CLUSTER_NORM2     = 4,
	CLUSTER_NORM1     = 5,
	SLIC_NORM1        = 6,
	SLIC_NORM2        = 7,
	HISTOGRAM_SSD     = 8,
	HISTOGRAM_SAD     = 9,
	VORONOI_NORM1     = 10,
	VORONOI_NORM2     = 11
} graph_edge_weight;

template<typename Scalar>
class GraphCut4D {
public:
	GraphCut4D();
	virtual ~GraphCut4D();

	/*
	 * Segment an image
	 *
	 * Returns a color image representing the segmentation.
	 *
	 * image: image to segment.
	 * output: segment map result
	 * sigma: to smooth the image.
	 * c: constant for treshold function.
	 * min_size: minimum component size (enforced by post-processing stage).
	 * @return num_ccs: number of connected components in the segmentation.
	 */

	int segmentImageStack_7(Image4D<Scalar>& image, Image3D<short>& output, float sigma, float constantTreshold, int min_size);
	int segmentImageStack_27(Image4D<Scalar>& image, Image3D<short>& output, float sigma, float constantTreshold, int min_size);
	int segmentImageStack_27weighted(Image4D<Scalar>& image, Image3D<short>& output, float sigma, float constantTreshold, int min_size);

private:
	/* convolve src with mask.  dst is flipped! */
	void convolve_x(Image4D<Scalar> *src, Image4D<Scalar> *dst, std::vector<float> &mask);
	void convolve_y(Image4D<Scalar> *src, Image4D<Scalar> *dst, std::vector<float> &mask);
	void convolve_z(Image4D<Scalar> *src, Image4D<Scalar> *dst, std::vector<float> &mask);
	void smooth(Image4D<Scalar>& input, Image4D<Scalar>& result, float sigma);
	Image4D<Scalar>* smooth(Image4D<Scalar>& input, float sigma);

	graph_edge_weight _weightMethod;
};

/*
template class GraphCut4D<char>;
template class GraphCut4D<unsigned char>;
template class GraphCut4D<short>;
template class GraphCut4D<unsigned short>;
template class GraphCut4D<int>;
template class GraphCut4D<unsigned int>;
template class GraphCut4D<long>;
template class GraphCut4D<unsigned long>;
template class GraphCut4D<float>;
template class GraphCut4D<double>;
*/

#endif /* GRAPHCUT4D_H_ */
