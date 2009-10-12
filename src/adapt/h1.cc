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

#include "../h3dconfig.h"
#include "../common.h"
#include "../solution.h"
#include "../refmap.h"
#include "../quad.h"
#include "../matrix.h"
#include "../traverse.h"
#include "../norm.h"
#include "h1.h"
#include "h1proj.h"
#include <common/timer.h>
#include <common/callstack.h>

// enable this to save refinements into a file
#undef SAVE_REFTS

#define PRINTF(...)
//#define PRINTF			printf

#ifndef COMPLEX

//

inline double int_h1_error(Solution *fu, Solution *fv, RefMap *ru, RefMap *rv) {
	_F_
	Quad3D *quad = get_quadrature(fv->get_active_element()->get_mode());

	order3_t o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	o.limit();
	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);

	fu->precalculate(np, pt, FN_DEFAULT);
	fv->precalculate(np, pt, FN_DEFAULT);

	scalar *uval = fu->get_fn_values();
	scalar *vval = fv->get_fn_values();
	scalar *dudx, *dudy, *dudz, *dvdx, *dvdy, *dvdz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudz);
	fv->get_dx_dy_dz_values(dvdx, dvdy, dvdz);

	double *jac = ru->get_jacobian(np, pt, true);
	double result = 0.0;
	for (int i = 0; i < np; i++)
		result += jac[i] * (sqr(uval[i] - vval[i]) + sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]) + sqr(dudz[i] - dvdz[i]));
	delete [] jac;
	return result;
}

inline double int_h1_norm(Solution *fu, RefMap *ru) {
	_F_
	Quad3D *quad = get_quadrature(fu->get_active_element()->get_mode());

	order3_t o = fu->get_fn_order() + ru->get_inv_ref_order();
	o.limit();

	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);

	fu->precalculate(np, pt, FN_DEFAULT);

	scalar *uval = fu->get_fn_values();
	scalar *dudx, *dudy, *dudz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudz);

	double *jac = ru->get_jacobian(np, pt, true);
	double result = 0.0;
	for (int i = 0; i < np; i++)
		result += jac[i] * (sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i]) + sqr(dudz[i]));
	delete [] jac;
	return result;
}

// H1 adapt ///////////////////////////////////////////////////////////////////////////////////////

H1Adapt::H1Adapt(int num, ...) {
	_F_
	this->num = num;

	va_list ap;
	va_start(ap, num);
	for (int i = 0; i < num; i++)
		spaces[i] = va_arg(ap, Space *);
	va_end(ap);

	memset(errors, 0, sizeof(errors));
	esort = NULL;
	have_errors = false;
}

H1Adapt::~H1Adapt() {
	_F_
	for (int i = 0; i < num; i++)
		if (errors[i] != NULL)
			delete [] errors[i];

	if (esort != NULL)
		delete [] esort;
}

double H1Adapt::get_projection_error(Element *e, int split, int son, order3_t order, Solution *rsln, Shapeset *ss) {
	_F_
	ProjKey key(split, son, order);
	double err;
	if (proj_err.lookup(key, err))
		return err;
	else {
		H1Projection proj(rsln, e, ss);
		err = proj.get_error(split, son, order);
		proj_err.set(key, err);
		return err;
	}


}

//// optimal refinement search /////////////////////////////////////////////////////////////////////

static inline int ndofs_elem(order3_t order) {
	assert(order.type == MODE_HEXAHEDRON);
	return (order.x + 1) * (order.y + 1) * (order.z + 1);
}

static inline int ndofs_bubble(order3_t order) {
	return (order.x - 1) * (order.y - 1) * (order.z - 1);
}

static inline int ndofs_face(int face, order3_t order1, order3_t order2) {
	order2_t forder[] = { order1.get_face_order(face), order2.get_face_order(face) };
	return (
		(std::min(forder[0].x, forder[1].x) - 1) *
		(std::min(forder[0].y, forder[1].y) - 1));
}

static inline int ndofs_edge(int edge, order3_t o) {
	return o.get_edge_order(edge) - 1;
}

static inline int ndofs_edge(int edge, order3_t o1, order3_t o2) {
	return std::min(o1.get_edge_order(edge), o2.get_edge_order(edge)) - 1;
}

static inline int ndofs_edge(int edge, order3_t o1, order3_t o2, order3_t o3, order3_t o4) {
	return
		std::min(
			std::min(o1.get_edge_order(edge), o2.get_edge_order(edge)),
			std::min(o3.get_edge_order(edge), o4.get_edge_order(edge))
		) - 1;
}

int H1Adapt::get_dof_count(int split, order3_t order[]) {
	_F_
	int dofs = 0;
	switch (split) {
		case REFT_HEX_NONE:
			dofs = ndofs_elem(order[0]);
			break;

		case REFT_HEX_X:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]);
			dofs -= ndofs_face(0, order[0], order[1]);			// face
			dofs -= 2 * ndofs_edge(3, order[0], order[1]); 		// edge
			dofs -= 2 * ndofs_edge(4, order[0], order[1]);
			dofs -= 4;											// vertex
			break;

		case REFT_HEX_Y:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]);
			dofs -= ndofs_face(2, order[0], order[1]);
			dofs -= 2 * ndofs_edge(0, order[0], order[1]); 		// edge
			dofs -= 2 * ndofs_edge(5, order[0], order[1]);
			dofs -= 4;											// vertex
			break;

		case REFT_HEX_Z:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]);
			dofs -= ndofs_face(4, order[0], order[1]);
			dofs -= 2 * ndofs_edge(0, order[0], order[1]); 		// edge
			dofs -= 2 * ndofs_edge(1, order[0], order[1]);
			dofs -= 4;											// vertex
			break;

		case REFT_HEX_XY:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]) + ndofs_elem(order[2]) + ndofs_elem(order[3]);
			dofs -= ndofs_face(1, order[0], order[1]);			// faces
			dofs -= ndofs_face(3, order[1], order[2]);
			dofs -= ndofs_face(0, order[2], order[3]);
			dofs -= ndofs_face(2, order[3], order[0]);
			dofs -= 5 * ndofs_edge(4, order[0]);				// edge
			dofs -= 2 * ndofs_edge(0, order[0], order[3]) + 2 * ndofs_edge(0, order[1], order[2]);
			dofs -= 2 * ndofs_edge(1, order[0], order[1]) + 2 * ndofs_edge(1, order[2], order[3]);
			dofs -= 10;											// vertex
			break;

		case REFT_HEX_XZ:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]) + ndofs_elem(order[2]) + ndofs_elem(order[3]);
			dofs -= ndofs_face(1, order[0], order[1]);			// faces
			dofs -= ndofs_face(5, order[1], order[2]);
			dofs -= ndofs_face(0, order[2], order[3]);
			dofs -= ndofs_face(4, order[3], order[0]);
			dofs -= 5 * ndofs_edge(1, order[0]);				// edge
			dofs -= 2 * ndofs_edge(0, order[0], order[3]) + 2 * ndofs_edge(0, order[1], order[2]);
			dofs -= 2 * ndofs_edge(4, order[0], order[1]) + 2 * ndofs_edge(4, order[2], order[3]);
			dofs -= 10;											// vertex
			break;

		case REFT_HEX_YZ:
			dofs = ndofs_elem(order[0]) + ndofs_elem(order[1]) + ndofs_elem(order[2]) + ndofs_elem(order[3]);
			dofs -= ndofs_face(3, order[0], order[1]);			// faces
			dofs -= ndofs_face(5, order[1], order[2]);
			dofs -= ndofs_face(2, order[2], order[3]);
			dofs -= ndofs_face(4, order[3], order[0]);
			dofs -= 5 * ndofs_edge(0, order[0]);				// edge
			dofs -= 2 * ndofs_edge(1, order[0], order[3]) + 2 * ndofs_edge(1, order[1], order[2]);
			dofs -= 2 * ndofs_edge(4, order[0], order[1]) + 2 * ndofs_edge(4, order[2], order[3]);
			dofs -= 10;											// vertex
			break;

		case REFT_HEX_XYZ:
			for (int i = 0; i < 8; i++) dofs += ndofs_elem(order[i]);
			dofs -= 15;			// vertex fns
			dofs -= ndofs_edge(0, order[0], order[4]) + ndofs_edge(0, order[3], order[7]) + ndofs_edge(0, order[0], order[3], order[4], order[7]);
			dofs -= ndofs_edge(0, order[1], order[5]) + ndofs_edge(0, order[2], order[6]) + ndofs_edge(0, order[1], order[2], order[5], order[6]);

			dofs -= ndofs_edge(1, order[0], order[4]) + ndofs_edge(1, order[1], order[5]) + ndofs_edge(1, order[0], order[1], order[4], order[5]);
			dofs -= ndofs_edge(1, order[3], order[7]) + ndofs_edge(1, order[2], order[6]) + ndofs_edge(1, order[2], order[3], order[6], order[7]);

			dofs -= ndofs_edge(4, order[0], order[1]) + ndofs_edge(4, order[2], order[3]) + ndofs_edge(4, order[0], order[1], order[2], order[3]);
			dofs -= ndofs_edge(4, order[4], order[5]) + ndofs_edge(4, order[6], order[7]) + ndofs_edge(4, order[4], order[5], order[6], order[7]);

			dofs -= ndofs_edge(4, order[0], order[3]) + ndofs_edge(4, order[4], order[7]);
			dofs -= ndofs_edge(4, order[1], order[2]) + ndofs_edge(4, order[5], order[6]);

			dofs -= ndofs_edge(1, order[0], order[1]) + ndofs_edge(1, order[2], order[3]);
			dofs -= ndofs_edge(1, order[4], order[5]) + ndofs_edge(1, order[6], order[7]);

			dofs -= ndofs_edge(0, order[1], order[2]) + ndofs_edge(0, order[0], order[3]);
			dofs -= ndofs_edge(0, order[5], order[6]) + ndofs_edge(0, order[4], order[7]);


			dofs -= ndofs_face(1, order[0], order[1]) + ndofs_face(3, order[1], order[2]);
			dofs -= ndofs_face(0, order[2], order[3]) + ndofs_face(2, order[0], order[3]);

			dofs -= ndofs_face(5, order[0], order[4]) + ndofs_face(5, order[1], order[5]);
			dofs -= ndofs_face(5, order[2], order[6]) + ndofs_face(5, order[3], order[7]);

			dofs -= ndofs_face(1, order[4], order[5]) + ndofs_face(3, order[5], order[6]);
			dofs -= ndofs_face(0, order[6], order[7]) + ndofs_face(2, order[7], order[4]);
			break;

		default: assert(false);
	}

	return dofs;
}

void H1Adapt::get_optimal_refinement(Mesh *mesh, Element *e, order3_t order, Solution *rsln, Shapeset *ss, int &split, order3_t p[8], bool aniso, bool h_adapt) {
	_F_
	int i, k, n = 0;
	const int MAX_CAND = 5000;

	int max_order = MAX_ELEMENT_ORDER;					// FIXME: will be diferent for curv. elements

	struct Cand {
		double error;
		int dofs, split;
		order3_t p[Hex::NUM_SONS];				// polynomial degree

		Cand() {
			error = 0.0;
			dofs = 0;
			split = -1;
			memset(p, 0, sizeof(p));
		}
	};
	Cand cand[MAX_CAND];

#define MAKE_P_CAND(q) { \
    assert(n < MAX_CAND);   \
    cand[n].split = REFT_HEX_NONE; \
    cand[n].p[1] = cand[n].p[2] = cand[n].p[3] = cand[n].p[4] = cand[n].p[5] = cand[n].p[6] = cand[n].p[7] = 0; \
    cand[n].p[0] = (q); \
    n++; }

#define MAKE_HP_CAND(q0, q1, q2, q3, q4, q5, q6, q7) { \
    assert(n < MAX_CAND);  \
    cand[n].split = REFT_HEX_XYZ; \
    cand[n].p[0] = (q0); \
    cand[n].p[1] = (q1); \
    cand[n].p[2] = (q2); \
    cand[n].p[3] = (q3); \
    cand[n].p[4] = (q4); \
    cand[n].p[5] = (q5); \
    cand[n].p[6] = (q6); \
    cand[n].p[7] = (q7); \
    n++; }

#define MAKE_ANI2_CAND(s, q0, q1) { \
	if (mesh->can_refine_element(e->id, s)) {\
		assert(n < MAX_CAND);  \
		cand[n].split = s; \
		cand[n].p[2] = cand[n].p[3] = cand[n].p[4] = cand[n].p[5] = cand[n].p[6] = cand[n].p[7] = 0; \
		cand[n].p[0] = (q0); \
		cand[n].p[1] = (q1); \
		n++; }}

#define MAKE_ANI4_CAND(s, q0, q1, q2, q3) { \
	if (mesh->can_refine_element(e->id, s)) {\
		assert(n < MAX_CAND);  \
		cand[n].split = s; \
		cand[n].p[4] = cand[n].p[5] = cand[n].p[6] = cand[n].p[7] = 0; \
		cand[n].p[0] = (q0); \
		cand[n].p[1] = (q1); \
		cand[n].p[2] = (q2); \
		cand[n].p[3] = (q3); \
		n++; }}

	// prepare p-candidates
	int dord[3] = { order.x, order.y, order.z };

	{
		int q[] = { std::min(max_order, dord[0] + 1), std::min(max_order, dord[1] + 1), std::min(max_order, dord[2] + 1) };

		MAKE_P_CAND(order);
		MAKE_P_CAND(order3_t(q[0], q[1], q[2]));
	}
	// prepare hp-candidates
	{
		order3_t pp[] = {
			order3_t(dord[0], dord[1], dord[2]),
			order3_t((dord[0] + 1) / 2, (dord[1] + 1) / 2, (dord[2] + 1) / 2),
			order3_t(std::min(((dord[0] + 1) / 2) + 1, max_order), std::min(((dord[1] + 1) / 2) + 1, max_order), std::min(((dord[2] + 1) / 2) + 1, max_order))
		};

		for (int q0 = 1; q0 < 3; q0++)
			for (int q1 = 1; q1 < 3; q1++)
				for (int q2 = 1; q2 < 3; q2++)
					for (int q3 = 1; q3 < 3; q3++)
						for (int q4 = 1; q4 < 3; q4++)
							for (int q5 = 1; q5 < 3; q5++)
								for (int q6 = 1; q6 < 3; q6++)
									for (int q7 = 1; q7 < 3; q7++)
										MAKE_HP_CAND(pp[q0], pp[q1], pp[q2], pp[q3], pp[q4], pp[q5], pp[q6], pp[q7]);
	}

#ifdef DEBUG_PRINT
	const char *split_str[] = {
		"NONE", "X   ", "Y   ", "Z   ", "XY  ", "XZ  ", "YZ  ", "XYZ "
	};
#endif

	// calculate their errors
	double avg = 0.0;
	double dev = 0.0;
	for (i = k = 0; i < n; i++) {
		Cand *c = cand + i;

		c->error = 0.0;
		switch (c->split) {
			case REFT_HEX_NONE:
				c->error += get_projection_error(e, c->split, -1, c->p[0], rsln, ss);
				break;

			case REFT_HEX_XYZ:
				for (int j = 0; j < 8; j++)
					c->error += get_projection_error(e, c->split, j, c->p[j], rsln, ss);
				break;

			case REFT_HEX_X:
				c->error += get_projection_error(e, c->split, 20, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 21, c->p[1], rsln, ss);
				break;

			case REFT_HEX_Y:
				c->error += get_projection_error(e, c->split, 22, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 23, c->p[1], rsln, ss);
				break;

			case REFT_HEX_Z:
				c->error += get_projection_error(e, c->split, 24, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 25, c->p[1], rsln, ss);
				break;

			case REFT_HEX_XY:
				c->error += get_projection_error(e, c->split,  8, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split,  9, c->p[1], rsln, ss);
				c->error += get_projection_error(e, c->split, 10, c->p[2], rsln, ss);
				c->error += get_projection_error(e, c->split, 11, c->p[3], rsln, ss);
				break;

			case REFT_HEX_XZ:
				c->error += get_projection_error(e, c->split, 12, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 13, c->p[1], rsln, ss);
				c->error += get_projection_error(e, c->split, 14, c->p[2], rsln, ss);
				c->error += get_projection_error(e, c->split, 15, c->p[3], rsln, ss);
				break;

			case REFT_HEX_YZ:
				c->error += get_projection_error(e, c->split, 16, c->p[0], rsln, ss);
				c->error += get_projection_error(e, c->split, 17, c->p[1], rsln, ss);
				c->error += get_projection_error(e, c->split, 18, c->p[2], rsln, ss);
				c->error += get_projection_error(e, c->split, 19, c->p[3], rsln, ss);
				break;

			default:
				EXIT(ERR_NOT_IMPLEMENTED);
				break;
		}
		c->error = sqrt(c->error);
		c->dofs = get_dof_count(c->split, c->p);

		if (!i || c->error <= cand[0].error) {
			avg += log10(c->error);
			dev += sqr(log10(c->error));
			k++;
		}
	}
	avg /= k; // mean
	dev /= k; // second moment
	dev = sqrt(dev - sqr(avg)); // deviation is square root of variance

	// select an above-average candidate with the steepest error decrease
	int imax = 0;
	double score, maxscore = 0.0;
	for (i = 1; i < n; i++) {
		if ((log10(cand[i].error) < avg + dev) && (cand[i].dofs > cand[0].dofs)) {
			score = (log10(cand[0].error) - log10(cand[i].error)) / (cand[i].dofs - cand[0].dofs);

			if (score > maxscore) {
				maxscore = score;
				imax = i;
			}
		}
	}

	// return result
	split = cand[imax].split;
	memcpy(p, cand[imax].p, Hex::NUM_SONS * sizeof(order3_t));

#ifdef DEBUG_PRINT
	printf(": best cand: #%d, split = %s", imax, split_str[cand[imax].split]);
	for (int i = 0; i < 8; i++)
		printf(", (%d, %d, %d)", p[i].x, p[i].y, p[i].z);
	printf(" | order = (%d, %d, %d)", order.x, order.y, order.z);
	printf("\n");
#endif
}

//// adapt /////////////////////////////////////////////////////////////////////////////////////////

void H1Adapt::adapt(double thr, bool h_only, int strat) {
	_F_
	if (!have_errors)
		EXIT("Element errors have to be calculated first, see calc_error().");

	Mesh *mesh[NUM];
	for (int j = 0; j < num; j++) {
		mesh[j] = spaces[j]->get_mesh();
		rsln[j]->enable_transform(false);
	}

#ifdef SAVE_REFTS
	// save the reft to the file
	FILE *file = fopen("adapt", "a");
	fprintf(file, "--\n");
#endif

	double err0 = 1000.0;
	double processed_error = 0.0;
	int i = 0;
	for (i = 0; i < nact; i++) {
		int comp = esort[i][1];
		int id = esort[i][0];
		double err = errors[comp][id];

		// first refinement strategy:
		// refine elements until prescribed amount of error is processed
		// if more elements have similar error refine all to keep the mesh symmetric
		if ((strat == 0) && (processed_error > sqrt(thr) * total_err) && fabs((err - err0) / err0) > 1e-3)
			break;

		// second refinement strategy:
		// refine all elements whose error is bigger than some portion of maximal error
		if ((strat == 1) && (err < thr * errors[esort[0][1]][esort[0][0]]))
			break;

		assert(mesh[comp]->elements.exists(id));
		Element *e = mesh[comp]->elements[id];
#ifdef DEBUG_PRINT
		printf("  - element #%d", id);
#endif

		int split = 0;
		order3_t p[8];													// polynomial order of sons
		for (int k = 0; k < 8; k++) p[k] = order3_t(0, 0, 0);
		order3_t cur_order = spaces[comp]->get_element_order(id);

		if (h_only) {
			p[0] = p[1] = p[2] = p[3] = p[4] = p[5] = p[6] = p[7] = cur_order;
			split = REFT_HEX_XYZ;
#ifdef DEBUG_PRINT
			printf("\n");			// new-line
#endif
		}
		else
			get_optimal_refinement(mesh[comp], e, cur_order, rsln[comp], spaces[comp]->get_shapeset(), split, p);

#ifdef SAVE_REFTS
		fprintf(file, "%ld %d %d %d %d %d %d %d %d %d\n", e->id, split, p[0].get_idx(), p[1].get_idx(), p[2].get_idx(), p[3].get_idx(),
				p[4].get_idx(), p[5].get_idx(), p[6].get_idx(), p[7].get_idx());
#endif

		switch (split) {
			case REFT_HEX_NONE:
				spaces[comp]->set_element_order(id, p[0]);
				break;

			case REFT_HEX_XYZ:
				mesh[comp]->refine_element(id, REFT_HEX_XYZ);
				for (int j = 0; j < Hex::NUM_SONS; j++)			// FIXME: hex specific
					spaces[comp]->set_element_order(e->get_son(j), p[j]);
				break;

			case REFT_HEX_X:
			case REFT_HEX_Y:
			case REFT_HEX_Z:
				mesh[comp]->refine_element(id, split);
				for (int j = 0; j < 2; j++)
					spaces[comp]->set_element_order(e->get_son(j), p[j]);
				break;

			case REFT_HEX_XY:
			case REFT_HEX_XZ:
			case REFT_HEX_YZ:
				mesh[comp]->refine_element(id, split);
				for (int j = 0; j < 4; j++)
					spaces[comp]->set_element_order(e->get_son(j), p[j]);
				break;

			default: assert(false);
		}

		err0 = err;
		processed_error += err;

		proj_err.remove_all();
	}

	for (int j = 0; j < num; j++)
		rsln[j]->enable_transform(true);

	have_errors = false;

	reft_elems = i;

#ifdef SAVE_REFTS
	fclose(file);
#endif
}

//// Unrefinements /////////////////////////////////////////////////////////////////////////////////

// TODO: unrefts

//// error calculation /////////////////////////////////////////////////////////////////////////////

static double **cmp_err;
static int compare(const void* p1, const void* p2) {
	const int2 (*e1) = ((const int2 *) p1);
	const int2 (*e2) = ((const int2 *) p2);
	return cmp_err[(*e1)[1]][(*e1)[0]] < cmp_err[(*e2)[1]][(*e2)[0]] ? 1 : -1;
}

double H1Adapt::calc_error_n(int n, ...) {
	_F_
	int i, j, k;

	if (n != num) EXIT("Wrong number of solutions.");

	va_list ap;
	va_start(ap, n);
	for (i = 0; i < n; i++) {
		sln[i] = va_arg(ap, Solution *);
		sln[i]->enable_transform(true);
	}
	for (i = 0; i < n; i++) {
		rsln[i] = va_arg(ap, Solution *);
		rsln[i]->enable_transform(true);
	}
	va_end(ap);

	nact = 0;
	for (j = 0; j < num; j++)
		nact += sln[j]->get_mesh()->get_num_active_elements();
	if (esort != NULL) delete [] esort;
	esort = new int2[nact];
	MEM_CHECK(esort);

	double total_error = 0.0, total_norm = 0.0;

	double norms[n];
	memset(norms, 0, n * sizeof(double));
	for (i = 0; i < n; i++)
		norms[i] = sqr(h1_norm(rsln[i]));

	for (j = k = 0; j < num; j++) {
		Mesh *cmesh = sln[j]->get_mesh();
		Mesh *fmesh = rsln[j]->get_mesh();

		int max = cmesh->get_max_element_id() + 1;
		if (errors[j] != NULL)
			delete [] errors[j];
		errors[j] = new double[max];
		MEM_CHECK(errors[j]);
		memset(errors[j], 0, sizeof(double) * max);

		FOR_ALL_ACTIVE_ELEMENTS(eid, cmesh) {
			Element *e = cmesh->elements[eid];
			EMode3D mode = e->get_mode();
			assert(mode == MODE_HEXAHEDRON);

			sln[j]->set_active_element(e);

			// prepare transformations (hex specific)
			int ns = 8;									// number of sons
			int trf[] = { 0, 1, 2, 3, 4, 5, 6, 7 };		// transformations corresponding to sons

			for (i = 0; i < ns; i++) {
				sln[j]->push_transform(trf[i]);

				assert(fmesh->elements.exists(e->id));
				Element *fe = fmesh->elements[e->id];
				assert(fe != NULL);
				Word_t son_idx = fe->get_son(i);
				if (fe->active || son_idx == INVALID_IDX) EXIT("Bad reference solution.");
				assert(fmesh->elements.exists(son_idx));
				Element *son = fmesh->elements[son_idx];
				assert(son != NULL);
				if (!son->active) EXIT("Bad reference solution (son not active).");
				assert(son->get_mode() == mode);

				rsln[j]->set_active_element(son);

				RefMap *crm = sln[j]->get_refmap();
				RefMap *frm = rsln[j]->get_refmap();

		        double err = int_h1_error(sln[j], rsln[j], crm, frm);
		        total_norm += int_h1_norm (rsln[j], frm);

		        errors[j][e->id] += err / norms[j];
		        total_error += err;

				sln[j]->pop_transform();
			}

			esort[k][0] = e->id;
			esort[k++][1] = j;
		}

		FOR_ALL_INACTIVE_ELEMENTS(eid, cmesh) {
			Element *e = cmesh->elements[eid];
			errors[j][e->id] = -1.0;
		}
	}

	assert(k == nact);
	cmp_err = errors;
	qsort(esort, nact, sizeof(int2), compare);

#ifdef DEBUG_PRINT
	Mesh *cmesh = sln[0]->get_mesh();
	for (int i = 0; i < std::min(10, (int) cmesh->elements.count()); i++) {
		printf("  - elem #% 3d: err = % e\n", esort[i][0], errors[0][esort[i][0]]);
	}
#endif

	have_errors = true;
	total_err = total_error / total_norm;

	return sqrt(total_error / total_norm);
}

#endif
