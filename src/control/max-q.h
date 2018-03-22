/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef MAX_Q_H_
#define MAX_Q_H_

#include "../control/v-control.h"
#include "../model/graph.h"
#include "../load/jacobian.h"
#include <armadillo>

using namespace arma;

namespace control {

class MaxQ: public VControl {
public:
	MaxQ();
	virtual ~MaxQ();

	virtual bool DoControl(mat jqv, Graph* graph);
	static const double LIMIAR = 0.01;
};

}

#endif /* MaxQ_H_ */
