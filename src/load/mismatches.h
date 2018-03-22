/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */


#ifndef MISMATCHES_H_
#define MISMATCHES_H_

#include "../model/graph.h"

#include <armadillo>

using namespace arma;
using namespace model;

namespace load {

class Mismatch {
private:
	vec m_mis;

	void CalcPkB(Graph* graph);
	void CalcQkB(Graph* graph);
public:
	Mismatch();
	virtual ~Mismatch();

	vec CalcMismatches(Graph* graph);

	vec GetMis(void);
};

}

#endif  // MISMATCHES_H_
