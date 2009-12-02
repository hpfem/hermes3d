// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _ADAPT_H1_H_
#define _ADAPT_H1_H_

/// hp-Adaptivity module for H1 space
///
/// TODO
///
///
/// @ingroup hp-adaptivity
class H1Adapt {
public:
	/// Initializes the class. 'num' is the number of mesh-space pairs to be adapted.
	/// After 'num', exactly that many space pointers must follow.
	H1Adapt(int num, ...);
	~H1Adapt();

	/// Type-safe version of calc_error_n() for one solution.
	double calc_error(Solution *sln, Solution *rsln) {
		if (num != 1) EXIT("Wrong number of solutions.");
		return calc_error_n(1, sln, rsln);
	}

	/// Type-safe version of calc_error_n() for two solutions.
	double calc_error_2(Solution *sln1,  Solution *sln2, Solution *rsln1, Solution *rsln2) {
		if (num != 2) EXIT("Wrong number of solutions.");
		return calc_error_n(2, sln1, sln2, rsln1, rsln2);
	}

	/// Calculates the error of the solution. 'n' must be the same
	/// as 'num' in the constructor. After that, n coarse solution
	/// pointers are passed, followed by n fine solution pointers.
	double calc_error_n(int n, ...);


	typedef scalar (*biform_t)(ScalarFunction *fu, ScalarFunction *fv, RefMap *ru, RefMap *rv);

	/// Selects elements to refine (based on results from calc_error() or calc_energy_error())
	/// and performs their optimal hp-refinement.
	void adapt(double thr, bool h_only = false, int strat = 0);

	/// Get the number of elements refined in the last adaptivity iteration
	int get_num_refined_elements() { return reft_elems; }

protected:
	// spaces & solutions
	static const int NUM = 10;

	int num;
	Space *spaces[NUM];
	Solution *sln[NUM];
	Solution *rsln[NUM];

	// element error arrays
	double *errors[NUM];
	double  norms[NUM]; // ?
	bool    have_errors;
	double  total_err;
	int2 *esort;
	int   nact;
	int reft_elems;

	/// Used by adapt(). Can be utilized in specialized adaptivity
	/// procedures, for which adapt() is not sufficient.
	void get_optimal_refinement(Mesh *mesh, Element *e, order3_t order, Solution *rsln, Shapeset *ss, int &split, order3_t p[8], bool aniso = true, bool h_adapt = false);
	double get_projection_error(Element *e, int split, int son, order3_t order, Solution *rsln, Shapeset *ss);
	int get_dof_count(int split, order3_t order[]);

	struct ProjKey {
		int split;			// transformation index
		int son;
		order3_t order;			// element order

		ProjKey(int t, int s, order3_t o) {
			split = t;
			son = s;
			order = o;
		}
	};
	Map<ProjKey, double> proj_err;				// cache for projection errors
};

#endif
