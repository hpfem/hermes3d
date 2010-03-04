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

#ifndef _LINPROBLEM_H_
#define _LINPROBLEM_H_

#include <common/array.h>
#include "weakform.h"

class Space;
class Matrix;
class Vector;
class FacePos;

class LinProblem {
public:
	LinProblem(WeakForm *wf);
	virtual ~LinProblem();
	void free();

	void set_spaces(int n, ...);

	// @return true if successful, otherwise false
	bool assemble(Matrix *matrix, Vector *rhs = NULL);

	double get_time() { return time; }

protected:
	WeakForm *wf;

	int ndofs;					/// number of DOFs
	int *sp_seq;				/// sequence numbers of spaces
	Space **spaces;
	double time;				/// time of the assembling (in secs)

	scalar **matrix_buffer;		/// buffer for holding square matrix (during assembling)
	int matrix_buffer_dim;		/// dimension of the matrix held by 'matrix_buffer'
	inline scalar **get_matrix_buffer(int n);

	void create(Matrix *matrix, Vector *rhs = NULL);

	// pre-transforming and fn. caching
	struct fn_key_t {
		int index;
		int order;
		int sub_idx;
		int ss_id;			// shapeset id

		fn_key_t(int index, int order, int sub_idx, int ss_id = -1) {
			this->index = index;
			this->order = order;
			this->sub_idx = sub_idx;
			this->ss_id = ss_id;
		}
	};

	struct FnCache {
		Array<double *> jwt;			// jacobian x weight
		Array<geom_t<double> > e;		// geometries
		Map<fn_key_t, sfn_t*> fn;		// shape functions
		Map<fn_key_t, mfn_t*> ext;		// external functions

		~FnCache() {
			free();
		}

		void free() {
			for (Word_t i = jwt.first(); i != INVALID_IDX; i = jwt.next(i))
				delete [] jwt[i];
			jwt.remove_all();
			for (Word_t i = e.first(); i != INVALID_IDX; i = e.next(i))
				free_geom(&e[i]);
			e.remove_all();
			for (Word_t i = fn.first(); i != INVALID_IDX; i = fn.next(i))
				free_fn(fn[i]);
			fn.remove_all();
			for (Word_t i = ext.first(); i != INVALID_IDX; i = ext.next(i))
				delete ext[i];
			ext.remove_all();
		}
	} fn_cache;

	scalar eval_form(WeakForm::BiFormVol *bf, ShapeFunction *fu, ShapeFunction *fv, RefMap *ru, RefMap *rv);
	scalar eval_form(WeakForm::LiFormVol *lf, ShapeFunction *fv, RefMap *rv);
	scalar eval_form(WeakForm::BiFormSurf *bf, ShapeFunction *fu, ShapeFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp);
	scalar eval_form(WeakForm::LiFormSurf *lf, ShapeFunction *fv, RefMap *rv, FacePos *fp);

	sfn_t *get_fn(ShapeFunction *fu, int order, RefMap *rm, const int np, const QuadPt3D *pt);
	sfn_t *get_fn(ShapeFunction *fu, int order, RefMap *rm, int iface, const int np, const QuadPt3D *pt);

	void init_ext_fns(user_data_t<ord_t> &fake_ud, std::vector<MeshFunction *> &ext);
	void init_ext_fns(user_data_t<scalar> &ud, std::vector<MeshFunction *> &ext, int order, RefMap *rm, const int np, const QuadPt3D *pt);

};

#endif /* _LINPROBLEM_H_ */
