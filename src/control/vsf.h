/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef VSF_H_
#define VSF_H_

#include "../control/v-control.h"
#include "../model/graph.h"
#include "../load/jacobian.h"

#include <armadillo>

using namespace arma;
using namespace load;

namespace control
{

class Vsf : public VControl
{
public:
	Vsf();
	virtual ~Vsf();

	virtual bool DoControl (mat jac, Graph* graph);
	static const double LIMIAR = 0.01;
};

}

#endif /* VSF_H_ */
