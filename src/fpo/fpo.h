/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef SRC_FPO_H_
#define SRC_FPO_H_

#include "../model/graph.h"
#include "../load/load-flow.h"

#include "../load/calc.h"

#include <armadillo>
#include <string>

using namespace model;
using namespace arma;
using namespace load;

namespace fpo {

class Fpo {
public:
	Fpo();
	virtual ~Fpo();

	void SetErro(double erro);
	double GetErro(void) const;

	void SetW(double w);
	double GetW(void) const;

	void SetMaxIter(double maxIter);
	double GetMaxIter(void) const;

	void SetRo(double ro);
	double GetRo(void) const;

	void SetVerbose(bool verbose);
	bool GetVerbose(void) const;

	void SetU(vec u);
	vec GetU(void) const;

	void Execute(std::string cdf);

private:
	double m_erro;
	// Parâmetro de Penalidade
	double m_w;
	double m_maxIter;
	// Fator de ajuste de w. Se ro=1, w não é ajustado.
	double m_ro;
	double m_co;

	int m_k;
	bool m_verbose;

	// Vetor de Variáveis de Controle
	vec m_u;
	vec m_x;
	vec m_foMod;
	vec m_perdas;

	LoadFlow* lf;

	Calc* m_calc;
};

}

#endif /* SRC_FPO_H_ */
