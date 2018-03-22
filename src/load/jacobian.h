/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef JACOBIAN_H_
#define JACOBIAN_H_

#include "../model/graph.h"

#include <armadillo>

using namespace arma;
using namespace model;

namespace load {

class Jacobian {
public:
	Jacobian();
	virtual ~Jacobian();

	mat GetMatrix() const;
	void SetMatrix(int numLines, int numCols);

	void SetJ1(int numRows, int numCols);
	mat GetJ1() const;

	void SetJ2(int numRows, int numCols);
	mat GetJ2() const;

	void SetJ3(int numRows, int numCols);
	mat GetJ3() const;

	void SetJ4(int numRows, int numCols);
	mat GetJ4() const;

	mat CalcJac(Graph* graph);
	vec SolveSys(vec b);

	void Zeros(void);
	void Zeros(int numRows, int numCols);
	void Zeros(Graph* graph);

	mat GetJqv(void);

private:
	mat m_matrix;
	mat m_j1;
	mat m_j2;
	mat m_j3;
	mat m_j4;
	mat m_jqv;

	void CalcDPk(Graph* graph);
	void CalcDQk(Graph* graph);
	void CalcJDQ(Graph* graph);
	void CalcJqv(void);
};

}

#endif  // JACOBIAN_H_
