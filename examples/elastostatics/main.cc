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

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

const double E  = 200e9; 		// steel: 200GPa
const double nu = 0.3;
const double f  = 1e4;   		// 10^4 N

const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

// usage info

void usage() {
	printf("Usage\n");
	printf("\n");
	printf("  elastostatics <mesh-file>\n");
	printf("\n");
}

// integrals

template<typename f_t, typename res_t>
res_t int_a_dx_b_dy_c_dz(double a, double b, double c, int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e) {
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (a * u->dx[i] * v->dx[i] + b * u->dy[i] * v->dy[i] + c * u->dz[i] * v->dz[i]);
	return res;
}

template<typename f_t, typename res_t>
res_t int_a_dudx_dvdy_b_dudy_dvdx(double a, double b, int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e) {
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (a * u->dx[i] * v->dy[i] + b * u->dy[i] * v->dx[i]);
	return res;
}

template<typename f_t, typename res_t>
res_t int_a_dudx_dvdz_b_dudz_dvdx(double a, double b, int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e) {
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (a * u->dx[i] * v->dz[i] + b * u->dz[i] * v->dx[i]);
	return res;
}

template<typename f_t, typename res_t>
res_t int_a_dudy_dvdz_b_dudz_dvdy(double a, double b, int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e) {
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (a * u->dy[i] * v->dz[i] + b * u->dz[i] * v->dy[i]);
	return res;
}

// BCs

EBCType bc_types_x(int marker) {
	return BC_NATURAL;
}

EBCType bc_types_y(int marker) {
	return BC_NATURAL;
}

EBCType bc_types_z(int marker) {
	return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL;
}

// 1. equation

template<typename f_t, typename res_t>
res_t bilinear_form_0_0(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_a_dx_b_dy_c_dz<f_t, res_t>(lambda + 2*mu, mu, mu, n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form_0_1(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_a_dudx_dvdy_b_dudy_dvdx<f_t, res_t>(lambda, mu, n, wt, v, u, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form_0_2(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_a_dudx_dvdz_b_dudz_dvdx<f_t, res_t>(lambda, mu, n, wt, v, u, e);
}

template<typename f_t, typename res_t>
res_t surf_linear_form_0(int n, double *wt, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return 0.0;
}

// 2. equation

template<typename f_t, typename res_t>
res_t bilinear_form_1_1(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_a_dx_b_dy_c_dz<f_t, res_t>(mu, lambda + 2*mu, mu, n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form_1_2(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_a_dudy_dvdz_b_dudz_dvdy<f_t, res_t>(lambda, mu, n, wt, v, u, e);
}

template<typename f_t, typename res_t>
res_t surf_linear_form_1(int n, double *wt, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return 0.0;
}

// 3. equation

template<typename f_t, typename res_t>
res_t bilinear_form_2_2(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_a_dx_b_dy_c_dz<f_t, res_t>(mu, mu, lambda + 2*mu, n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t surf_linear_form_2(int n, double *wt, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (f * v->fn[i]);
	return res;
}

//

void out_fn(MeshFunction *x, MeshFunction *y, MeshFunction *z, const char *name) {
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s.vtk", name);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(x, y, z, name);
		fclose(ofile);
	}
	else {
		error("Can not not open '%s' for writing.", of_name);
	}
}


// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif

	if (argc < 2) {
		usage();
		return 0;
	}

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mloader;
	if (!mloader.load(argv[1], &mesh)) die("Loading mesh file '%s'\n", argv[1]);

	mesh.refine_all_elements(REFT_HEX_XYZ);
	mesh.refine_all_elements(REFT_HEX_XYZ);
	mesh.refine_all_elements(REFT_HEX_XYZ);

	H1ShapesetLobattoHex shapeset;

	printf("* Setting the space up\n");
	H1Space xdisp(&mesh, &shapeset);
	xdisp.set_bc_types(bc_types_x);
	xdisp.set_uniform_order(order3_t(2, 2, 2));

	H1Space ydisp(&mesh, &shapeset);
	ydisp.set_bc_types(bc_types_y);
	ydisp.set_uniform_order(order3_t(2, 2, 2));

	H1Space zdisp(&mesh, &shapeset);
	zdisp.set_bc_types(bc_types_z);
	zdisp.set_uniform_order(order3_t(2, 2, 2));

	int ndofs = 0;
	ndofs += xdisp.assign_dofs(ndofs);
	ndofs += ydisp.assign_dofs(ndofs);
	ndofs += zdisp.assign_dofs(ndofs);
	printf("  - Number of DOFs: %d\n", ndofs);

	// weak formulation
	WeakForm wf(3);
	wf.add_biform(0, 0, bilinear_form_0_0<double, scalar>, bilinear_form_0_0<ord_t, ord_t>, SYM);
	wf.add_biform(0, 1, bilinear_form_0_1<double, scalar>, bilinear_form_0_1<ord_t, ord_t>, SYM);
	wf.add_biform(0, 2, bilinear_form_0_2<double, scalar>, bilinear_form_0_2<ord_t, ord_t>, SYM);
	wf.add_liform_surf(0, surf_linear_form_0<double, scalar>, surf_linear_form_0<ord_t, ord_t>);

	wf.add_biform(1, 1, bilinear_form_1_1<double, scalar>, bilinear_form_1_1<ord_t, ord_t>, SYM);
	wf.add_biform(1, 2, bilinear_form_1_2<double, scalar>, bilinear_form_1_2<ord_t, ord_t>, SYM);
	wf.add_liform_surf(1, surf_linear_form_1<double, scalar>, surf_linear_form_1<ord_t, ord_t>);

	wf.add_biform(2, 2, bilinear_form_2_2<double, scalar>, bilinear_form_2_2<ord_t, ord_t>, SYM);
	wf.add_liform_surf(2, surf_linear_form_2<double, scalar>, surf_linear_form_2<ord_t, ord_t>, 5);

	LinProblem lp(&wf);
	lp.set_spaces(3, &xdisp, &ydisp, &zdisp);

#if defined WITH_UMFPACK
	UMFPackMatrix mat;
	UMFPackVector rhs;
	UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PARDISO
	PardisoMatrix mat;
	PardisoVector rhs;
	PardisoLinearSolver solver(&mat, &rhs);
#elif defined WITH_PETSC
	PetscMatrix mat;
	PetscVector rhs;
	PetscLinearSolver solver(&mat, &rhs);
#elif defined WITH_MUMPS
	MumpsMatrix mat;
	MumpsVector rhs;
	MumpsSolver solver(&mat, &rhs);
#endif

	// assemble stiffness matrix
	printf("  - assembling... "); fflush(stdout);
	Timer tmr_assemble;
	tmr_assemble.start();
	bool assembled = lp.assemble(&mat, &rhs);
	tmr_assemble.stop();
	if (assembled)
		printf("done in %s (%lf secs)\n", tmr_assemble.get_human_time(), tmr_assemble.get_seconds());
	else
		die("failed!");

	// solve the stiffness matrix
	printf("  - solving... "); fflush(stdout);
	Timer tmr_solve;
	tmr_solve.start();
	bool solved = solver.solve();
	tmr_solve.stop();

	if (solved) {
		printf("done in %s (%lf secs)\n", tmr_solve.get_human_time(), tmr_solve.get_seconds());
	}
	else {
		printf("failed\n");
		return -1;
	}

	double *s = solver.get_solution();
	Solution xsln(&mesh), ysln(&mesh), zsln(&mesh);
	xsln.set_fe_solution(&xdisp, s);
	ysln.set_fe_solution(&ydisp, s);
	zsln.set_fe_solution(&zdisp, s);

	printf("  - output... "); fflush(stdout);
	out_fn(&xsln, &ysln, &zsln, "disp");
	printf("done\n");

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	return 0;
}
