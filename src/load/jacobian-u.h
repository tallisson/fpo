/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef JACOBIAN_U_H_
#define JACOBIAN_U_H_

#include <armadillo>

#include "../model/graph.h"

using namespace arma;
using namespace model;

namespace load {

class JacobianU {
public:
	JacobianU();
	virtual ~JacobianU();

	mat CalcJac(Graph* graph);
	mat GetJac(void) const;
private:
	mat m_jac;

	void CalcDPk(Graph* graph);
	void CalcDQk(Graph* graph);
};

} /* namespace load */

#endif /* JACOBIAN_U_H_ */
