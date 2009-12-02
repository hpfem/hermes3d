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

// Solving a simple heat equation to demonstrate how to use CUBIT with Hermes3D
//
// Use mesh file from 'meshes/exodusII/cylynder2.e. Material IDs corresponds
// to elements markers, sideset IDs correspond to face (BC) markers
//

EBCType bc_types(int marker)
{
	if (marker == 1) return BC_ESSENTIAL;
	else return BC_NATURAL;
}

double bc_values(int marker, double x, double y, double z)
{
	return 10;
}

template<typename f_t, typename res_t>
res_t bilinear_form1(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                    user_data_t<res_t> *data)
{
	return 10 * int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form2(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                    user_data_t<res_t> *data)
{
	return 0.5 * int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T rhs(T x, T y, T z)
{
	return 4.0;
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	return -int_F_v<f_t, res_t>(n, wt, rhs, u, e);
}

//

void out_fn(MeshFunction *x, const char *name)
{
#ifdef OUTPUT_DIR
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s/%s.vtk", OUTPUT_DIR, name);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(x, name);
		fclose(ofile);
	}
	else {
		error("Can not not open '%s' for writing.", of_name);
	}
#endif
}

void out_bc(Mesh *mesh, const char *name)
{
#ifdef OUTPUT_DIR
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s/%s.vtk", OUTPUT_DIR, name);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out_bc(mesh, name);
		fclose(ofile);
	}
	else {
		error("Can not not open '%s' for writing.", of_name);
	}
#endif
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif

	if (argc < 2) die("Not enough parameters");

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	ExodusIIReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) die("Loading mesh file '%s'\n", argv[1]);

	H1ShapesetLobattoHex shapeset;

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_bc_values(bc_values);
	space.set_uniform_order(order3_t(1, 1, 1));

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	WeakForm wf(1);
	wf.add_biform(0, 0, bilinear_form1<double, scalar>, bilinear_form1<ord_t, ord_t>, SYM, 1);
	wf.add_biform(0, 0, bilinear_form2<double, scalar>, bilinear_form2<ord_t, ord_t>, SYM, 2);
	wf.add_liform(0, linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY);

	LinProblem lp(&wf);
	lp.set_spaces(1, &space);

#if defined WITH_UMFPACK
	UMFPackMatrix mat;
	UMFPackVector rhs;
	UMFPackLinearSolver solver(mat, rhs);
#elif defined WITH_PARDISO
	PardisoMatrix mat;
	PardisoVector rhs;
	PardisoLinearSolver solver(mat, rhs);
#elif defined WITH_PETSC
	PetscMatrix mat;
	PetscVector rhs;
	PetscLinearSolver solver(mat, rhs);
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
		double *s = solver.get_solution();

		Solution sln(&mesh);
		sln.set_fe_solution(&space, s);

#ifdef OUTPUT_DIR
		printf("  - output... "); fflush(stdout);
		out_bc(&mesh, "bc");
		out_fn(&sln, "temp");
		printf("done\n");
#endif
	}
	else {
		printf("failed\n");
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	return 0;
}
