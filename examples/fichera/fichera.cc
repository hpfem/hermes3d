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

//
// fichera.cc
//
// usage: $0 <mesh file>
//
//

#include "config.h"
#ifdef USE_PETSC
#include <petsc.h>
#endif
#ifdef USE_UMFPACK
#include <umfpack.h>
#endif

#include <hermes3d.h>

// error should be smaller than this epsilon
#define EPS								10e-14F

double fnc(double x, double y, double z) {
	return pow(x*x + y*y + z*z, .25);
}

double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 0.5 * x * pow(x*x + y*y + z*z, -.75);
	dy = 0.5 * y * pow(x*x + y*y + z*z, -.75);
	dz = 0.5 * z * pow(x*x + y*y + z*z, -.75);

	return fnc(x, y, z);
}

// weak formulation

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

double bc_values(int marker, double x, double y, double z) {
	return fnc(x, y, z);
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T dfnc(T x, T y, T z) {
	return -0.75 * pow(x*x + y*y + z*z, -0.75);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, dfnc, u, e);
}


// helpers ////////////////////////////////////////////////////////////////////////////////////////

void out_orders(Space *space, const char *name, int iter) {
#ifdef OUTPUT_DIR
	char fname[1024];
	sprintf(fname, "%s/iter-%d-%s.gmsh", OUTPUT_DIR, iter, name);
	FILE *f = fopen(fname, "w");
	if (f != NULL) {
		GmshOutputEngine gmsh(f);
		gmsh.out_orders(space, name);
		fclose(f);
	}
	else
		error("Could not open file '%s' for writing.", fname);
#endif
}

void out_fn(MeshFunction *fn, const char *name, int iter) {
#ifdef OUTPUT_DIR
	char fname[1024];
	sprintf(fname, "%s/iter-%d-%s.gmsh", OUTPUT_DIR, iter, name);
	FILE *f = fopen(fname, "w");
	if (f != NULL) {
		GmshOutputEngine gmsh(f);
		gmsh.out(fn, name);
		fclose(f);
	}
	else
		error("Could not open file '%s' for writing.", fname);
#endif
}

//

const double TOLERANCE = 0.001;		// error tolerance in percent
const double THRESHOLD = 0.3;		// error threshold for element refinement

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
	PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL); // disable PETSc error handler
#endif

	if (argc < 2) die("Not enough parameters");

	// load the inital mesh
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) die("Unable to load mesh file '%s'\n", args[1]);

	H1ShapesetLobattoHex shapeset;

#if defined WITH_UMFPACK
	UMFPackMatrix mat;
	UMFPackVector rhs;
	UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PETSC
	PetscMatrix mat;
	PetscVector rhs;
	PetscLinearSolver solver(&mat, &rhs);
#elif defined WITH_MUMPS
	MumpsMatrix mat;
	MumpsVector rhs;
	MumpsSolver solver(&mat, &rhs);
#endif

	// Graphs
	GnuplotGraph graph;
	graph.set_captions("", "Degrees of Freedom", "Error [%]");
	graph.set_log_y();
	graph.add_row("Total error", "k", "-", "O");

	// problem setup
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_bc_values(bc_values);
	space.set_uniform_order(order3_t(2, 2, 2));

	WeakForm wf(1);
	wf.add_biform(0, 0, bilinear_form<double, double>, bilinear_form<ord_t, ord_t>, SYM, ANY, 0);
	wf.add_liform(0, linear_form<double, double>, linear_form<ord_t, ord_t>, ANY, 0);

	// ADAPT loop
	LinProblem lp(&wf);
	lp.set_spaces(1, &space);

	int iter = 0;
	bool done = false;
	do {
		printf("\n=== Iter #%d ================================================================\n", iter);

		Timer t_assemble;
		Timer t_solver;

		printf("\nSolution\n");

		// assign DOFs
		int ndofs = space.assign_dofs();
		printf("  - Number of DOFs: %d\n", ndofs);
		// assemble stiffness matrix and rhs
		printf("  - Assembling... "); fflush(stdout);
		t_assemble.start();
		lp.assemble(&mat, &rhs);
		t_assemble.stop();
		printf("done in %s (%lf secs)\n", t_assemble.get_human_time(), t_assemble.get_seconds());

		// solve the system
		printf("  - Solving... "); fflush(stdout);
		t_solver.start();
		bool solved = solver.solve();
		t_solver.stop();
		if (solved)
			printf("done in %s (%lf secs)\n", t_solver.get_human_time(), t_solver.get_seconds());
		else {
			printf("Failed\n");
			break;
		}

		// construct a solution
		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

		// output the orders and the solution
		out_orders(&space, "order", iter);
		out_fn(&sln, "sln", iter);

		// reference solution
		printf("Reference solution\n");

#if defined WITH_UMFPACK
		UMFPackLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PETSC
		PetscLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_MUMPS
		MumpsSolver rsolver(&mat, &rhs);
#endif

		// construct the mesh for reference solution
		Mesh rmesh;
		rmesh.copy(mesh);
		rmesh.refine_all_elements(REFT_HEX_XYZ);
		// setup space for the reference solution
		Space *rspace = space.dup(&rmesh);
		rspace->copy_orders(space, 1);

		LinProblem rlp(&wf);
		rlp.set_spaces(1, rspace);

		// assign DOFs
		int rndofs = rspace->assign_dofs();
		printf("  - Number of DOFs: %d\n", rndofs);

		// assemble stiffness matric and rhs
		printf("  - Assembling... "); fflush(stdout);
		t_assemble.start();
		rlp.assemble(&mat, &rhs);
		t_assemble.stop();
		printf("done in %s (%lf secs)\n", t_assemble.get_human_time(), t_assemble.get_seconds());

		// solve the system
		printf("  - Solving... "); fflush(stdout);
		t_solver.start();
		bool rsolved = rsolver.solve();
		t_solver.stop();
		if (rsolved)
			printf("done in %s (%lf secs)\n", t_solver.get_human_time(), t_solver.get_seconds());
		else {
			printf("failed\n");
			break;
		}

		// construct the reference solution
		Solution rsln(&rmesh);
		rsln.set_fe_solution(rspace, rsolver.get_solution());

		// calculate the error estimate
		double err = h1_error(&sln, &rsln);
		printf("  - H1 error: % lf\n", err * 100);
#ifdef OUTPUT_DIR
		// save it to the graph
		graph.add_value(0, ndofs, err * 100);
		graph.save(OUTPUT_DIR"/conv.gp");
#endif

		// do the hp-adaptivity
		printf("Adaptivity\n");
		printf("  - tolerance: "); fflush(stdout);
		H1Adapt hp(1, &space);
		Timer t_err;
		t_err.start();
		double tol = hp.calc_error(&sln, &rsln) * 100;		// calc error estimates on elements
		t_err.stop();
		printf("% lf\n", tol);

		if (tol < TOLERANCE) {
			// we are within the tolerance, so we can stop
			printf("\nDone\n");
			break;
		}

		printf("  - adapting... "); fflush(stdout);
		Timer t_adapt;
		t_adapt.start();
		hp.adapt(THRESHOLD);								// run the adaptivity algorithm
		t_adapt.stop();
		printf("done in %lf secs (refined %d element(s))\n", t_adapt.get_seconds(), hp.get_num_refined_elements());

		delete rspace;										// clean-up

		// next iteration
		iter++;

		mat.free();
		rhs.free();
	} while (!done);

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	return res;
}
