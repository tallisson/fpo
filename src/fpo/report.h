/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef REPORT_H_
#define REPORT_H_

#include "../model/graph.h"

#include <string>
#include <stdio.h>

using namespace model;
using namespace std;

namespace fpo {

class Report {
public:
	Report();
	virtual ~Report();

	void Output(Graph* graph);
};

} /* namespace fpo */

#endif /* REPORT_H_ */
