// This file is part of Hermes3D
//
// Copyright (c) 2010 hp-FEM group at the University of Nevada, Reno (UNR).
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
#include <float.h>
#include <getopt.h>

// problem constants
const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;

// commnad line arguments
bool do_output = true;				// generate output files (if true)
char *mesh_file_name = NULL;		// the name of the mesh file

// usage info

void usage() {
	printf("Usage\n");
	printf("\n");
	printf("  bessel <mesh-file>\n");
	printf("\n");
	printf("Options:\n");
	printf("  --no-output         - do not generate output files\n");
	printf("\n");
}

bool process_cmd_line(int argc, char **argv)
{
	static struct option long_options[] = {
		{ "no-output", no_argument, (int *) &do_output, false },
		{ 0, 0, 0, 0 }
	};

	// getopt_long stores the option index here.
	int option_index = 0;
	int c;
	while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
		switch (c) {
            case 0:
				break;

			case '?':
				// getopt_long already printed an error message
				break;

			default:
				return false;
		}
	}

	if (optind < argc) {
		mesh_file_name = argv[optind++];
		return true;
	}
	else
		return false;
}

// Bessel function of the first kind, order n, defined in bessel.cpp
double jv(double n, double x);

static void exact_sol_val(double x, double y, double z, scalar &e0, scalar &e1)
{
	double t1 = x*x;
	double t2 = y*y;
	double t4 = sqrt(t1+t2);
	double t5 = jv(-1.0/3.0,t4);
	double t6 = 1/t4;
	double t7 = jv(2.0/3.0,t4);
	double t11 = (t5-2.0/3.0*t6*t7)*t6;
	double t12 = atan2(y,x);
	if (t12 < 0) t12 += 2.0*M_PI;
	double t13 = 2.0/3.0*t12;
	double t14 = cos(t13);
	double t17 = sin(t13);
	double t18 = t7*t17;
	double t20 = 1/t1;
	double t23 = 1/(1.0+t2*t20);
	e0 = t11*y*t14-2.0/3.0*t18/x*t23;
	e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
}

static void exact_sol(double x, double y, double z, scalar &e0, scalar &e1, scalar &e1dx, scalar &e0dy)
{
	exact_sol_val(x, y, z, e0, e1);

	double t1 = x*x;
	double t2 = y*y;
	double t3 = t1+t2;
	double t4 = sqrt(t3);
	double t5 = jv(2.0/3.0,t4);
	double t6 = 1/t4;
	double t7 = jv(-1.0/3.0,t4);
	double t11 = (-t5-t6*t7/3.0)*t6;
	double t14 = 1/t4/t3;
	double t15 = t14*t5;
	double t21 = t7-2.0/3.0*t6*t5;
	double t22 = 1/t3*t21;
	double t27 = atan2(y,x);
	if (t27 < 0) t27 += 2.0*M_PI;
	double t28 = 2.0/3.0*t27;
	double t29 = cos(t28);
	double t32 = t21*t14;
	double t35 = t21*t6;
	double t36 = t35*t29;
	double t39 = sin(t28);
	double t41 = 1/t1;
	double t43 = 1.0+t2*t41;
	double t44 = 1/t43;
	double t47 = 4.0/3.0*t35/x*t39*y*t44;
	double t48 = t5*t29;
	double t49 = t1*t1;
	double t52 = t43*t43;
	double t53 = 1/t52;
	double t57 = t5*t39;
	double t59 = 1/t1/x;
	e1dx =-(t11*x+2.0/3.0*t15*x-2.0/3.0*t22*x)
		*t6*x*t29+t32*t1*t29-t36-t47+4.0/9.0*t48*t2/t49*t53+4.0/3.0*t57*y*t59*t44-4.0/3.0*t57*t2*y/t49/x*t53;
	e0dy = (t11*y+2.0/3.0*t15*y-2.0/3.0*t22*y)*t6*y*t29-t32*t2*t29+t36-t47-4.0/9.0*t48*t41*t53+4.0/3.0*t57*t59*t53*y;
}

// exact solution
scalar3 &exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz)
{
	static scalar3 ex;

	ex[0] = ex[1] = ex[2] = 0;
	exact_sol(x, y, z, ex[0], ex[1], dx[1], dy[0]);
	return ex;
}


// BCs

EBCType bc_types(int marker)
{
	if (marker == 1 || marker == 6)
		return BC_ESSENTIAL; // perfect conductor
	else
		return BC_NATURAL; // impedance
}

template<typename f_t, typename res_t>
res_t biform(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *ext)
{
	return 1.0/mu_r * hcurl_int_curl_u_curl_v<f_t, res_t>(n, wt, u, v, e)
		- sqr(kappa) * hcurl_int_u_v<f_t, res_t>(n, wt, u, v, e);
}

ord_t biform_surf_ord(int n, double *wt, fn_t<ord_t> *u, fn_t<ord_t> *v, geom_t<ord_t> *e, user_data_t<ord_t> *ext)
{
	return ord_t(v->fn[0].get_max_order());
}

scalar biform_surf(int n, double *wt, fn_t<double> *u, fn_t<double> *v, geom_t<double> *e, user_data_t<scalar> *ext)
{
	// j * kappa * E_T * F_T
	// E_T = nu x E x nu  (nu is outer normal)
	std::complex<double> ii = std::complex<double>(0.0, 1.0);
	scalar result = 0;
	for (int i = 0; i < n; i++) {
		scalar uu[3] = { u->fn0[i], u->fn1[i], u->fn2[i] };
		scalar tpu[3];
		calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], uu, tpu);

		scalar vv[3] = { v->fn0[i], v->fn1[i], v->fn2[i] };
		scalar tpv[3];
		calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], vv, tpv);

		result += wt[i] * (uu[0] * vv[0] + uu[1] * vv[1] + uu[2] * vv[2]);
	}

	return ii * (-kappa) * result;
}

scalar liform_surf(int n, double *wt, fn_t<double> *v, geom_t<double> *e, user_data_t<scalar> *ext)
{
	std::complex<double> ii = std::complex<double>(0.0, 1.0);
	scalar result = 0;
	for (int i = 0; i < n; i++) {
		scalar dx[3], dy[3], dz[3];
		scalar3 ev = exact_solution(e->x[i], e->y[i], e->z[i], dx, dy, dz);

		scalar curl_e[3];
		calc_curl(dx, dy, dz, curl_e);
		scalar tpe[3];
		calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], ev, tpe);

		scalar g[3] = {
			(e->nz[i] * curl_e[1] - e->ny[i] * curl_e[2]) - ii * kappa * tpe[0],
			(e->nx[i] * curl_e[2] - e->nz[i] * curl_e[0]) - ii * kappa * tpe[1],
			(e->ny[i] * curl_e[0] - e->nx[i] * curl_e[1]) - ii * kappa * tpe[2],
		};

		// tpv is tangencial projection of v (test function)
		scalar vv[3] = { v->fn0[i], v->fn1[i], v->fn2[i] };
		scalar tpv[3];
		calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], vv, tpv);

		result += wt[i] * (g[0] * tpv[0] + g[1] * tpv[1] + g[2] * tpv[2]);
	}

	return result;
}

// maximal polynomial order to integrate surface linear form
ord_t liform_surf_ord(int n, double *wt, fn_t<ord_t> *v, geom_t<ord_t> *e, user_data_t<ord_t> *ext)
{
	return ord_t(v->fn[0].get_max_order());
}

void out_mesh(Mesh *mesh, const char *name)
{
	char fname[1024];
	sprintf(fname, "%s.vtk", name);
	FILE *f = fopen(fname, "w");
	if (f != NULL) {
		VtkOutputEngine vtk(f);
		vtk.out(mesh);
		fclose(f);
	}
	else
		error("Could not open file '%s' for writing.", fname);
}

void out_fn(MeshFunction *fn, const char *name)
{
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s.vtk", name);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(fn, name);
		fclose(ofile);
	}
	else {
		error("Can not not open '%s' for writing.", of_name);
	}
}


// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif

	if (!process_cmd_line(argc, argv)) {
		usage();
		return 0;
	}

	printf("* Loading mesh '%s'\n", mesh_file_name);
	Mesh mesh;
	Mesh3DReader mloader;
	if (!mloader.load(mesh_file_name, &mesh))
		die("Loading mesh file '%s'\n", mesh_file_name);

	mesh.refine_all_elements(REFT_HEX_XYZ);
	mesh.refine_all_elements(REFT_HEX_XY);
	mesh.refine_all_elements(REFT_HEX_XY);
	mesh.refine_all_elements(REFT_HEX_XY);

	HcurlShapesetLobattoHex shapeset;

	printf("* Setting the space up\n");
	HcurlSpace sp(&mesh, &shapeset);
	sp.set_bc_types(bc_types);
	sp.set_uniform_order(order3_t(1, 1, 1));

	int ndofs = sp.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	// weak formulation
	WeakForm wf(1);
	wf.add_biform(0, 0, biform<double, scalar>, biform<ord_t, ord_t>, SYM);
	wf.add_biform_surf(0, 0, biform_surf, biform_surf_ord);
	wf.add_liform_surf(0, liform_surf, liform_surf_ord);

	LinProblem lp(&wf);
	lp.set_spaces(1, &sp);

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
	lp.assemble(&mat, &rhs);
	tmr_assemble.stop();
	printf("done in %s (%lf secs)\n", tmr_assemble.get_human_time(), tmr_assemble.get_seconds());

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

	std::complex<double> *s = solver.get_solution();
	Solution sln(&mesh);
	sln.set_fe_solution(&sp, s);

	if (do_output) {
		printf("  - output... "); fflush(stdout);
		out_fn(&sln, "solution");
		out_mesh(&mesh, "mesh");
		printf("done\n");
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	return 0;
}
