/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef V_CONTROL_H_
#define V_CONTROL_H_

#include "../model/graph.h"

#include <armadillo>

using namespace model;

namespace control {

class VControl {
public:
	VControl();
	virtual ~VControl();

	Bus* MaxDsv(Graph* graph);
	int MaxV(Graph* graph, arma::vec vt, Bus* modBus);
	virtual bool DoControl(arma::mat jqv, Graph* graph) = 0;

};

}

#endif /* V_CONTROL_H_ */
