/*
 * graphEdge.h
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#ifndef GRAPHEDGE_H_
#define GRAPHEDGE_H_

#ifndef _CXX_COMPILE_
#include "../../common/src/std_typedefs.h"
#else
#include "Utilities/common/src/std_typedefs.h"
#endif

class graphEdge {
public:
	graphEdge();
	graphEdge(const graphEdge& X);
	virtual ~graphEdge();

	inline bool operator<(const graphEdge &X) {
		return this->w < X.w;
	}
	inline bool operator==(const graphEdge &X) {
		return (((this->a==X.a) && (this->b==X.b)) || ((this->a==X.b) && (this->b==X.a)));
	}
	inline bool operator!=(const graphEdge &X) {
		return (((this->a!=X.a) && (this->b!=X.b)) && ((this->a!=X.b) && (this->b!=X.a)));
	}

	unsigned long a;
	unsigned long b;
	float w;
};

#endif /* GRAPHEDGE_H_ */
