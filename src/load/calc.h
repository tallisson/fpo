/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef CALC_H_
#define CALC_H_

#include "../model/graph.h"
#include "load-flow.h"
#include <armadillo>
#include <string>
#include <sstream>

using namespace arma;
using namespace model;
using namespace std;

namespace load {

struct data {
	double m_c;
	int m_conv;
	string print(void) {
		ostringstream os;
		os << "c = " << m_c << ", conv = " << m_conv;
		return os.str();
	}
};
typedef data Data_t;

struct aresta {
	aresta(int k, int m) : m_k(k), m_m(m){}
	int m_k;
	int m_m;
};
typedef aresta Aresta_t;
class Calc {
public:
	Calc();
	virtual ~Calc();

	vec GrafX(Graph* graph, double w);
	vec GradLU(Graph* graph, vec lamb);
	vec GradFU(Graph* graph);
	double Fitness(Graph* graph, double w);
	Data_t Busca(LoadFlow* lf, Graph* graph, double w, double c, double erro, vec dLdu);

	void IsVerbose(bool verbose);
private:
	bool m_verbose;
};

}

#endif /* CALC_H_ */
