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
#include "aztecoo.h"
#include <common/callstack.h>

#define AZTECOO_NOT_COMPILED "hermes3d was not built with AztecOO support."

// AztecOO solver //////////////////////////////////////////////////////////////////////////////////

AztecOOSolver::AztecOOSolver(EpetraMatrix &m, EpetraVector &rhs)
	: Solver(), m(m), rhs(rhs)
{
	_F_
#ifdef HAVE_AZTECOO
	// set default values
	max_iters = 10000;
	tolerance = 10e-8;
	pc = NULL;
#else
	die(AZTECOO_NOT_COMPILED);
#endif
}

AztecOOSolver::~AztecOOSolver()
{
	_F_
#ifdef HAVE_AZTECOO
#endif
}

void AztecOOSolver::set_solver(const char *name)
{
	_F_
#ifdef HAVE_AZTECOO
	int az_solver;
	if (strcasecmp(name, "gmres") == 0) az_solver = AZ_gmres;
	else if (strcasecmp(name, "cg") == 0) az_solver = AZ_cg;
	else if (strcasecmp(name, "cgs") == 0) az_solver = AZ_cgs;
	else if (strcasecmp(name, "tfqmr") == 0) az_solver = AZ_tfqmr;
	else if (strcasecmp(name, "bicgstab") == 0) az_solver = AZ_bicgstab;
	else az_solver = AZ_gmres;

	aztec.SetAztecOption(AZ_solver, az_solver);
#endif
}

void AztecOOSolver::set_precond(const char *name)
{
	_F_
#ifdef HAVE_AZTECOO
	int az_precond;
	if (strcasecmp(name, "none") == 0) az_precond = AZ_none;
	else if (strcasecmp(name, "jacobi") == 0) az_precond = AZ_Jacobi;
	else if (strcasecmp(name, "neumann") == 0) az_precond = AZ_Neumann;
	else if (strcasecmp(name, "least-squares") == 0) az_precond = AZ_ls;
	else az_precond = AZ_none;

	aztec.SetAztecOption(AZ_precond, az_precond);
#endif
}

bool AztecOOSolver::solve()
{
	_F_
#ifdef HAVE_AZTECOO
	assert(m.size == rhs.size);

	// no output
	aztec.SetAztecOption(AZ_output, AZ_none);	// AZ_all | AZ_warnings | AZ_last | AZ_summary

	// setup the problem
	aztec.SetUserMatrix(m.mat);
	aztec.SetRHS(rhs.vec);
	Epetra_Vector x(*rhs.std_map);
	aztec.SetLHS(&x);

	if (pc != NULL) {
		Epetra_Operator *op = pc->get_obj();
		assert(op != NULL);		// can work only with Epetra_Operators
		aztec.SetPrecOperator(op);
	}

	// solve it
	aztec.Iterate(max_iters, tolerance);

	delete [] sln;
	sln = new scalar[m.size];
	MEM_CHECK(sln);
	memset(sln, 0, m.size * sizeof(scalar));

	// copy the solution into sln vector
	for (int i = 0; i < m.size; i++) sln[i] = x[i];

	return true;
#else
	return false;
#endif
}

int AztecOOSolver::get_num_iters()
{
	_F_
#ifdef HAVE_AZTECOO
	return aztec.NumIters();
#else
	return -1;
#endif
}

double AztecOOSolver::get_residual()
{
	_F_
#ifdef HAVE_AZTECOO
	return aztec.TrueResidual();
#else
	return -1.0;
#endif
}

