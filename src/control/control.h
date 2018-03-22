/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef CONTROL_H_
#define CONTROL_H_

#include "../model/graph.h"

using namespace model;

namespace control {

class Control {
public:
	virtual bool DoControl(Graph* graph) = 0;
	virtual bool DoRestore(Graph* graph) = 0;
	Control();
	virtual ~Control();
};

}

#endif /* CONTROL_H_ */
