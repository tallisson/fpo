/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 * Author: Rafael <enzasampaiof@hotmail.com>
 */

#ifndef BUS_H_
#define BUS_H_

#include "branch.h"

#include <armadillo>
#include <string>
#include <vector>

using namespace arma;

namespace model {

/*
 * Struct to store buses variables
 */
struct bus {
public:
	int m_nex;
	int m_nin;
	std::string m_nome;
	double m_area;
	double m_zona;
	int m_tipo;
	double m_v;
	double m_ang;
	double m_pc;
	double m_qc;
	double m_pg;
	double m_qg;
	double m_base_kV;
	double m_vg_o;
	double m_qgmax;
	double m_qgmin;
	double m_vmax;
	double m_vmin;
	double m_gsh;
	double m_bsh;
	double m_crt;
	int m_ctrl_rem;
	int m_ordPV;
	int m_posPV;
	int m_ordPQ;
	int m_posPQ;
	int m_ord;
};
typedef bus DBus_t;

class Bus {
public:
	enum Type {
		SLACK = 3, GENERATION = 2, LOAD = 0, LOSS_CONTROL_REACT = 4, NONE = -1
	};

	enum TapType {
		IMP = 0, TAP = 1
	};

	enum Violation {
		NORMAL = 0, MIN_VOLTAGE_VIOLATION = 1, MAX_VOLTAGE_VIOLATION = 2
	};

	static const double MIN_VOLTAGE_GR = 0.95;
	static const double MIN_VOLTAGE_IEEE = 1.10;
	static const double MIN_VOLTAGE_ONS = 0.95;
	static const double MAX_VOLTAGE_ONS = 1.05;
	static const double MAX_VOLTAGE_IEEE = 0.90;
	static const double INCREMENT_STEP = 0.01;

	Bus();
	virtual ~Bus();

	double GetACalc();
	double GetVCalc();
	double GetPCalc();
	double GetQCalc();

	void SetACalc(double aCalc);
	void SetVCalc(double vCalc);
	void SetPCalc(double pCalc);
	void SetQCalc(double qCalc);

	DBus_t GetBus(void) const;
	void SetBus(DBus_t bus);

	Type GetType(void);
	void SetType(enum Type);

	void AddBranch(Branch* branch, Bus* neighbor);

	std::vector<Branch*> GetBranches() const;
	std::vector<Bus*> GetNeighbors() const;

	double CalcPg(void);
	double CalcQg(void);

	void SetTap(TapType tap);
	TapType GetTap(void) const;

	void Print(void);

	void SetCrt(double crt);
	double GetCrt(void) const;

	double CalcDsv(void);
	double GetDsv(void);

	void SetVPreBusca(double v);
	double GetVPreBusca(void) const;
	void SetAPreBusca(double a);
	double GetAPreBusca(void) const;

	Bus::Violation GetStatus(void);

	bool IsControlled(void);

	void SetOrdG(int ordG);
	int GetOrdG(void) const;

	void AddControl(double value);
	vec GetControl(void) const;
	void ClearControl(void);
	double CalcPG(void);
private:
	double m_aCalc;
	double m_vCalc;
	double m_qgCalc;
	double m_pgCalc;
	double m_dsv;
	double m_ordG;
	double m_vPreBusca;
	double m_aPreBusca;

	TapType m_tap;

	DBus_t m_bus;
	std::vector<Branch*> m_branches;
	std::vector<Bus*> m_neighbors;

	Type m_type;
	// @ToDo
	double m_crt;
	std::vector<double> m_control;

	Violation m_status;

	bool m_isControlled;
	int m_numControl;
};
}
#endif /* BUS_H_ */
