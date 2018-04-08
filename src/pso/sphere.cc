#include "sphere.h"
#include <armadillo>

using namespace arma;

namespace pso {

Sphere::Sphere() {
	// TODO Auto-generated constructor stub

}

Sphere::~Sphere() {
	// TODO Auto-generated destructor stub
}

double Sphere::CostFunction(vec x) {
	vec result = x % x;
	return sum(result);
}

} /* namespace pso */
