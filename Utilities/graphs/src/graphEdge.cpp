/*
 * graphEdge.cpp
 *
 *  Created on: Apr 5, 2018
 *      Author: christian
 */

#include "graphEdge.hpp"

graphEdge::graphEdge() : a(0), b(0), w(1.0) {
	// TODO Auto-generated constructor stub

}

graphEdge::graphEdge(const graphEdge& X) : a(X.a), b(X.b), w(X.w) {

}

graphEdge::~graphEdge() {
	// TODO Auto-generated destructor stub
}
