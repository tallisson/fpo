/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef LOAD_FLOW_H_
#define LOAD_FLOW_H_

#include "../model/graph.h"
#include "jacobian.h"
#include "mismatches.h"
#include "../utils/utils.h"
#include "../control/control.h"
#include "../control/q-control.h"
#include "../control/max-q.h"
#include "../control/vsf.h"

#include <string>
#include <sstream>

using namespace model;
using namespace utils;
using namespace control;


namespace load {

struct dataLoss {
	vec m_pkm;
	vec m_pmk;
	vec m_qkm;
	vec m_qmk;
	vec m_l;
	double m_total;

	dataLoss(vec pkm, vec pmk, vec qkm, vec qmk, vec l, double total) {
		m_pkm = pkm;
		m_pmk = pmk;
		m_qkm = qkm;
		m_qmk = qmk;
		m_l = l;
		m_total = total;
	}
};
typedef dataLoss DataLoss_t;

class LoadFlow {
private:
	Graph* m_graph;
	Control* m_qControl;
	Mismatch* m_mismatches;
	Jacobian* m_jac;
	VControl* m_vControl;

	Sts_t m_sts;

	double m_error;
	double m_totalL;
	uint32_t m_iter;
	uint32_t m_maxIter;

	vec m_losses;
	vec m_b;

	void InitX0(void);

	void InitJ(void);

	bool m_verbose;
	bool m_hasQControl;
	std::string m_dir;
	std::string m_file;

	std::ostringstream os;

	int m_conv;

	bool m_hasInit;
public:
	LoadFlow(void);
	virtual ~LoadFlow(void);

	Graph* GetGraph(void) const;

	Control* GetQControl(void) const;
	void SetQControl(Control* qControl);

	void SetVControl(VControl* vControl);

	void Prepare(std::string cdf);

	int Execute();

	void SetError(double error);

	void SetDir(std::string dir);

	void SetFile(std::string file);

	void SetQControl(bool hasQControl);

	void Reset(void);

	void SetVerbose(bool verbose);
	bool GetVerbose(void) const;

	mat ExecuteJ(void);

	int GetConv(void) const;

	Jacobian* GetJac(void) const;

	void SetInit(bool init);

	DataLoss_t CalcLosses(void);
};

}

#endif /* LOAD_FLOW */
