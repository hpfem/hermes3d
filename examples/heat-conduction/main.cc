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
#include <getopt.h>
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

double tm = 0.0;							// "global" time
const double tau = 0.05;					// time step


// error should be smaller than this epsilon
#define EPS								10e-10F

// commnad line arguments
bool do_output = true;				// generate output files (if true)
char *mesh_file_name = NULL;		// the name of the mesh file

// usage info

void usage() {
	printf("Usage\n");
	printf("\n");
	printf("  heat-conduction [options] <mesh-file>\n");
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

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = -2 * sin(tm) * x * (1 - y*y) * (1 - z*z);
	dy = -2 * sin(tm) * (1 - x*x) * y * (1 - z*z);
	dz = -2 * sin(tm) * (1 - x*x) * (1 - y*y) * z;

	return sin(tm) * (1 - x*x) * (1 - y*y) * (1 - z*z);
}

double exact_solution_prev(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = -2 * sin(tm - tau) * x * (1 - y*y) * (1 - z*z);
	dy = -2 * sin(tm - tau) * (1 - x*x) * y * (1 - z*z);
	dz = -2 * sin(tm - tau) * (1 - x*x) * (1 - y*y) * z;

	return sin(tm - tau) * (1 - x*x) * (1 - y*y) * (1 - z*z);
}

//

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return
		int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e) +
		int_u_v<f_t, res_t>(n, wt, u, v, e) / tau;
}

template<typename T>
T f(T x, T y, T z) {
	T ddxx = -2 * sin(tm) * (1 - y*y) * (1 - z*z);
	T ddyy = -2 * sin(tm) * (1 - x*x) * (1 - z*z);
	T ddzz = -2 * sin(tm) * (1 - x*x) * (1 - y*y);
	T dt = cos(tm) * (1 - x*x) * (1 - y*y) * (1 - z*z);

	return dt - (ddxx + ddyy + ddzz);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return
		int_F_v<f_t, res_t>(n, wt, f, u, e) +
		int_u_v<f_t, res_t>(n, wt, data->ext + 0, u, e) / tau;
}

//

void out_fn(MeshFunction *x, const char *name, int i) {
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s-%d.vtk", name, i);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(x, name);
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

	if (!process_cmd_line(argc, argv)) {
		usage();
		return 0;
	}

	printf("* Loading mesh '%s'\n", mesh_file_name);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(mesh_file_name, &mesh)) die("Loading mesh file '%s'\n", mesh_file_name);

	mesh.refine_all_elements(REFT_HEX_XYZ);
	mesh.refine_all_elements(REFT_HEX_XYZ);

	H1ShapesetLobattoHex shapeset;

	Solution uprev(&mesh);

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_uniform_order(order3_t(2, 2, 2));

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	// initial condition is zero
	uprev.set_zero();

	// discretization
	WeakForm wf(1);
	wf.add_biform(0, 0, bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
	wf.add_liform(0, linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY, 1, &uprev);

	LinProblem lp(&wf);
	lp.set_spaces(1, &space);

	Solution sln(&mesh);

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

	// main loop
	int n_iter = 2 * M_PI / tau;
	tm = tau;
	for (int i = 0; i < n_iter; i++) 	{
		printf("\n-- Iteration %d ---------------------\n", i);

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
			break;
		}

		double *s = solver.get_solution();
		sln.set_fe_solution(&space, s);

		ExactSolution esln(&mesh, exact_solution);

		if (do_output) {
			printf("  - output... "); fflush(stdout);
			out_fn(&sln, "temp", i);
			out_fn(&esln, "ex", i);
			printf("done\n");
		}

		// check our solution
		// norm
		double h1_err_norm = h1_error(&sln, &esln);
		printf("  - H1 error norm:      % le\n", h1_err_norm);

		// next step
		uprev = sln;
		tm += tau;
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	return 0;
}
