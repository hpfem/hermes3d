#include "config.h"
#ifdef USE_PETSC
#include <petsc.h>
#endif
#ifdef USE_UMFPACK
#include <umfpack.h>
#endif
#include <getopt.h>
#include <hermes3d.h>

//  This is another example that allows you to compare h- and hp-adaptivity from the point of view
//  of both CPU time requirements and discrete problem size, look at the quality of the a-posteriori
//  error estimator used by Hermes (exact error is provided), etc. 
//  The problem is made harder for adaptive algorithms by increasing the parameter SLOPE.
//
//  PDE: -Laplace u = f.
//
//  Known exact solution, see functions fn() and fndd().
//
//  Domain: unit square (0, 0, 1)x(0, 1, 0)x(1, 0, 0), see the file hexahedron.mesh3d.
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:
//  
//
//  usage: $0 <mesh file>
//
//

const double ERR_STOP  = 1;			// Stopping criterion for adaptivity (rel. error tolerance between the
						// fine mesh and coarse mesh solution in percent).
const double THRESHOLD = 0.3;			// error threshold for element refinement of the adapt(...) function 
						// (default) STRATEGY = 0 ... refine elements elements until sqrt(THRESHOLD) 
						// times total error is processed. If more elements have similar errors, 
						// refine all to keep the mesh symmetric.
						// STRATEGY = 1 ... refine all elements whose error is larger
						// than THRESHOLD times maximum element error.
const int P_INIT = 2;				// Initial polynomial degree of all mesh elements.
const int NDOF_STOP = 100000;			// Adaptivity process stops when the number of degrees of freedom grows
						// over this limit. This is to prevent h-adaptivity to go on forever.
bool do_output = true;				// generate output files (if true); Cmd line arguments

// Error should be smaller than this epsilon. 
#define EPS	10e-14F

// Problem parameters. 
double SLOPE = 200.0;				// slope of the layer inside the domain


// Exact solution.
#include "exact_solution.cpp"

// Boundary condition types
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return fn(x, y, z);
}

// Weak forms.
#include "forms.cpp"

// Output the order of polynomial
void out_orders(Space *space, const char *name, int iter)
{
  char fname[1024];
  sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out_orders(space, name);
    fclose(f);
  }
  else
    error("Could not open file '%s' for writing.", fname);
}

// Output the solution
void out_fn(MeshFunction *fn, const char *name, int iter)
{
  char fname[1024];
  sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out(fn, name);
    fclose(f);
  }
  else
    error("Could not open file '%s' for writing.", fname);
}

/***********************************************************************************
 * main program                                                                    *
************************************************************************************/
int main(int argc, char **args) 
{
  int res = ERR_SUCCESS;						// error code for debug purpose

#ifdef WITH_PETSC
  PetscInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
  PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);		// disable PETSc error handler
#endif

  // Load the inital mesh.
  Mesh mesh;
  Mesh3DReader mesh_loader;
  mesh_loader.load("hexahedron.mesh3d", &mesh))				// hexahedron used

  //Initialize the shapeset and the cache.
  H1ShapesetLobattoHex shapeset;

  //Matrix solver.
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

  // Graphs of DOF convergence.
  GnuplotGraph graph;
  graph.set_captions("", "Degrees of Freedom", "Error [%]");
  graph.set_log_y();
  graph.add_row("Total error", "k", "-", "O");

  // Create H1 space to setup the problem.
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(order3_t(P_INIT, P_INIT, P_INIT));

  // Initialize the week formulation. 
  WeakForm wf(1);
  wf.add_biform(0, 0, biform<double, double>, biform<ord_t, ord_t>, SYM, ANY, 0);
  wf.add_liform(0, liform<double, double>, liform<ord_t, ord_t>, ANY, 0);

  // Initialize the coarse mesh problem.
  LinProblem lp(&wf);
  lp.set_spaces(1, &space);

  // Adaptivity loop.
  int iter = 0;  bool done = false;
  do {
    printf("\n---- Adaptivity step %d\n", iter);

    printf("\nSolving on coarse mesh\n");

    // Procedures for coarse mesh problem.
    // Assign DOFs.
    int ndof = space.assign_dofs();
    printf("  - Number of DOFs: %d\n", ndof);

    // Assemble stiffness matrix and rhs.
    printf("  - Assembling... "); fflush(stdout);
    if (lp.assemble(&mat, &rhs))
      printf("done in %lf secs\n", lp.get_time());
    else
      die("failed!");

    // Solve the system.
    printf("  - Solving... "); fflush(stdout);
    bool solved = solver.solve();
    if (solved)
      printf("done in %lf secs\n", solver.get_time());
    else 
    {
      printf("Failed\n");
      break;
    }

    // Construct a solution.
    Solution sln(&mesh);
    sln.set_fe_solution(&space, solver.get_solution());

    if (do_output) 
    {
    // Output the orders and the solution.
    out_orders(&space, "order", iter);
    out_fn(&sln, "sln", iter);
    }

    // Procedures for fine mesh problem.
    // Reference(refined) solution.
    printf("Solving on fine mesh\n");

    //Matrix solver.
#if defined WITH_UMFPACK
    UMFPackLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PETSC
    PetscLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_MUMPS
    MumpsSolver rsolver(&mat, &rhs);
#endif

    // Construct the refined mesh for reference(refined) solution.
    Mesh rmesh;
    rmesh.copy(mesh);
    rmesh.refine_all_elements(REFT_HEX_XYZ);		// REFT_HEX_XYZ is hexahedron refiment type

    // Setup space for the reference(refined) solution.
    Space *rspace = space.dup(&rmesh);
    rspace->copy_orders(space, 1);

    // Initialize the mesh problem for reference(refined) solution.
    LinProblem rlp(&wf);
    rlp.set_spaces(1, rspace);

    // Assign DOF.
    int rndof = rspace->assign_dofs();
    printf("  - Number of DOFs: %d\n", rndof);

    // Assemble stiffness matric and rhs.
    printf("  - Assembling... "); fflush(stdout);
    if (rlp.assemble(&mat, &rhs))
      printf("done in %lf secs\n", rlp.get_time());
    else
      die("failed!");

    // Solve the system.
    printf("  - Solving... "); fflush(stdout);
    bool rsolved = rsolver.solve();
    if (rsolved)
      printf("done in %lf secs\n", rsolver.get_time());
    else 
    {
      printf("failed\n");
      break;
    }

    // Construct the reference(refined) solution.
    Solution rsln(&rmesh);
    rsln.set_fe_solution(rspace, rsolver.get_solution());

    // Compare coarse and fine mesh. 
    // Calculate the error estimate wrt. refined mesh solution. 
    double err = h1_error(&sln, &rsln);
    printf("  - H1 error: % lf\n", err * 100);

    // Save it to the graph.
    graph.add_value(0, ndof, err * 100);
    if (do_output)
      graph.save("conv.gp");

    // Do the hp-adaptivity.
    printf("Adaptivity\n");
    printf("  - calculating error: "); fflush(stdout);
    H1Adapt hp(1, &space);
    double err_est = hp.calc_error(&sln, &rsln) * 100;		// calc error estimates on elements
    printf("% lf\n", err_est);

    if (err_est < ERR_STOP) 
    {
      printf("\nDone\n");
      break;
    }

    printf("  - adapting... "); fflush(stdout);
    hp.adapt(THRESHOLD);					// run the adaptivity algorithm
    printf("done in %lf secs (refined %d element(s))\n", hp.get_adapt_time(), hp.get_num_refined_elements());

    if (rndof >= NDOF_STOP) 
    {
      printf("\nDone\n");
      break;
    }

    delete rspace;						// clean-up

    // Next iteration.
    iter++;

    mat.free();
    rhs.free();
  } while (!done);

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	return res;
}
