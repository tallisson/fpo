#ifndef SPHERE_H_
#define SPHERE_H_

#include "problem.h"

namespace pso {

class Sphere: public Problem {
public:
	Sphere();
	virtual ~Sphere();

	virtual double CostFunction(vec x);
};

} /* namespace pso */

#endif /* SPHERE_H_ */
