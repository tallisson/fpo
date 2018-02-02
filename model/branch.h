/*
 * Copyright (c) 2016 UFPI (Universidade Federal do Piauí)
 * Author: Thiago Allisson <allissonribeiro02@gmail.com>
 * Author: Enza Rafaela <enzasampaiof@hotmail.com>
 */

#ifndef BRANCH_H_
#define BRANCH_H_

namespace model {

/*
 * Struct to store branches variables
 */
struct branch {
public:
	int m_ni;
	int m_nf;
	int m_area;
	int m_zona;
	int m_circ;
	int m_tipo;
	int m_ordtap;
	int m_postap;
	double m_r;
	double m_x;
	double m_bsh;
	double m_line_rating1;
	double m_line_rating2;
	double m_line_rating3;
	double m_t_ctrl;
	double m_side;
	double m_tap;
	double m_def;
	double m_tapmin;
	double m_tapmax;
	double m_passo;
	double m_ctrl_min;
	double m_ctrl_max;
	double m_g;
	double m_b;
	double m_size;
};
typedef branch DBranch_t;

class Branch {
private:
	DBranch_t m_branch;

	double m_p_km;
	double m_p_mk;

	double m_q_km;
	double m_q_mk;

	double m_p_km_L;
	double m_p_mk_L;
	double m_q_km_L;
	double m_q_mk_L;

	double m_l;
public:
	Branch();
	virtual ~Branch();

	void SetBranch(DBranch_t branch);
	DBranch_t GetBranch() const;

	double CalcPkmL(double vK, double vM, double aK, double aM);
	double CalcPmkL(double vK, double vM, double aK, double aM);

	double CalcQkmL(double vK, double vM, double aK, double aM);
	double CalcQmkL(double vK, double vM, double aK, double aM);

	double CalcL(double vK, double vM, double aK, double aM);

	void Print(void);
};

}

#endif /* BRANCH_H_ */
