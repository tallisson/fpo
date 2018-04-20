/*
 * fpo-reactive.h
 *
 *  Created on: 17 de abr de 2018
 *      Author: pc-thiago
 */

#ifndef FPO_FPO_REACTIVE_H_
#define FPO_FPO_REACTIVE_H_

#include "../model/graph.h"
#include "../load/load-flow.h"

#include "../load/calc.h"

#include <armadillo>
#include <string>

using namespace arma;
using namespace std;

namespace fpo {

class FpoReactive {
public:
	FpoReactive();
	virtual ~FpoReactive();

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

	void Report(LoadFlow* lf, int numIt);

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

	LoadFlow* m_lf;

	Calc* m_calc;
};

} /* namespace load */

#endif /* SRC_FPO_FPO_REACTIVE_H_ */
