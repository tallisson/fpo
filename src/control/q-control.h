/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef QCONTROL_H_
#define QCONTROL_H_

#include "control.h"
#include "../model/graph.h"

using namespace model;

namespace control {

class QControl: public Control {
private:
	void UpdateOrd(Graph* graph);

public:
	QControl();
	virtual ~QControl();

	virtual bool DoControl(Graph* graph);
	virtual bool DoRestore(Graph* graph);
};

}

#endif /* QCONTROL_H_ */
