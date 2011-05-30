#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example shows how to solve a simple PDE that describes stationary 
// heat transfer in an object consisting of two materials (aluminum and 
// copper). The object is heated by constant volumetric heat sources
// generated by a DC electric current. The temperature on the boundary 
// is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format 
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const int P_MAG_INIT = 2;                             // Uniform polynomial degree of mesh elements.
const int P_TEMP_INIT = 2;
const int P_ELAST_INIT = 2;
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double A_INIT = 0.0;
const double TEMP_INIT = 20.0;
const double DK_INIT = 0.0;

const double TIME_STEP = 1.0;
const double TIME_FINAL = 10.0;

const double frequency = 50;

std::string *str_marker;

scalar maxB, max_el_cond, min_el_cond;

//when true, parameters are taken as functions, otherwise constanst from tables
const bool USE_NONLINEARITIES = true;


#include "tables.cpp"

// Weak forms.
#include "definitions.cpp"

using namespace RefinementSelectors;
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows


void adapt_mesh(Hermes::vector<Space*> spaces, WeakForm *wf)
{
    // Initialize refinement selector.
    H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

    Solution sln1, sln2, ref_sln1, ref_sln2;
    Hermes::vector<Solution*> solutions = Hermes::vector<Solution*>(&sln1, &sln2);
    Hermes::vector<Solution*> ref_solutions = Hermes::vector<Solution*>(&ref_sln1, &ref_sln2);

    // Initialize views.
    OrderView  oview("Polynomial orders", new WinGeom(420, 0, 400, 600));

    // DOF and CPU convergence graphs initialization.
    SimpleGraph graph_dof, graph_cpu;

    // Time measurement.
    TimePeriod cpu_time;
    cpu_time.tick();

    // Adaptivity loop:
    int as = 1;
    bool done = false;
    do
    {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      Hermes::vector<Space*> *ref_spaces = Space::construct_refined_spaces(spaces, 1);
      int ndof_ref = Space::get_num_dofs(*ref_spaces);

      // Initialize matrix solver.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Assemble reference problem.
      info("Solving on reference mesh.");
      DiscreteProblem* dp = new DiscreteProblem(wf, *ref_spaces, true);

      // Time measurement.
      cpu_time.tick();

      dp->assemble(matrix, rhs);
      solver->solve();
      Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, ref_solutions);

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(spaces, ref_solutions, solutions, matrix_solver);

      // Time measurement.
      cpu_time.tick();

      // Skip visualization time.
      cpu_time.tick(HERMES_SKIP);

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(spaces);
      bool solutions_for_adapt = true;
      // In the following function, the Boolean parameter "solutions_for_adapt" determines whether
      // the calculated errors are intended for use with adaptivity (this may not be the case, for example,
      // when error wrt. an exact solution is calculated). The default value is solutions_for_adapt = true,
      // The last parameter "error_flags" determine whether the total and element errors are treated as
      // absolute or relative. Its default value is error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL.
      // In subsequent examples and benchmarks, these two parameters will be often used with
      // their default values, and thus they will not be present in the code explicitly.
      //      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt,
      //                           HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;
      double err_est_rel = adaptivity->calc_err_est(solutions, ref_solutions);

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
        Space::get_num_dofs(spaces), Space::get_num_dofs(*ref_spaces), err_est_rel);

      // Time measurement.
      cpu_time.tick();

      // View the coarse mesh solution and polynomial orders.
      if (HERMES_VISUALIZATION) {
        oview.show(spaces[0]);
      }

      // Add entry to DOF and CPU convergence graphs.
      graph_dof.add_values(Space::get_num_dofs(spaces), err_est_rel);
      graph_dof.save("conv_dof_est.dat");
      graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
      graph_cpu.save("conv_cpu_est.dat");

      // If err_est too large, adapt the mesh.
      if (err_est_rel < ERR_STOP) done = true;
      else
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        // Increase the counter of performed adaptivity steps.
        if (done == false)  as++;
      }
      if (Space::get_num_dofs(spaces) >= NDOF_STOP) done = true;

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      //if(done == false) delete ref_spaces[0]->get_mesh();
      delete dp;

    }
    while (done == false);

}

int main(int argc, char* argv[])
{
    str_marker = new std::string[NUM_EDGES + NUM_LABELS];
    char str[5];
    for(int i = 0; i < NUM_EDGES + NUM_LABELS; i++){
        sprintf(str, "%d", i);
        str_marker[i].assign(str);
    }

    // Instantiate a class with global functions.
    Hermes2D hermes2d;

    // Initialize tables from the file tables.cpp
    initTables();

    // Load the mesh.
    Mesh mesh_mag, mesh_temp, mesh_elast;
    H2DReader mloader;
    mloader.load("kotouc_mesh_mag.mesh", &mesh_mag);
    mloader.load("kotouc_mesh_temp.mesh", &mesh_temp);
    mloader.load("kotouc_mesh_elast.mesh", &mesh_elast);

//    MeshView mv_mag;
//    mv_mag.show(&mesh_mag);
//    MeshView mv_temp;
//    mv_temp.show(&mesh_temp);
//    MeshView mv_elast;
//    mv_elast.show(&mesh_elast);
//
//    mv_elast.wait();

    // Perform initial mesh refinements (optional).
     for (int i=0; i < INIT_REF_NUM; i++){
         mesh_temp.refine_all_elements();
         mesh_mag.refine_all_elements();
         mesh_elast.refine_all_elements();
    }

     Solution sln_mag_real(&mesh_mag, A_INIT);
     Solution sln_mag_imag(&mesh_mag, A_INIT);;
     Solution sln_temp(&mesh_temp, TEMP_INIT);
     DoNothingFilter filter_temp(&sln_temp);

    // Initialize the weak formulation.
    WeakFormMagnetic wf(2);
    wf.registerForms(magneticLabels, &sln_mag_real, &sln_mag_imag, &filter_temp);

    // Initialize boundary conditions.
    EssentialBCs bcs_mag;
    for(int i = 0; i < NUM_EDGES; i++){
        if(magneticEdge[i].type == PhysicFieldBC_Magnetic_VectorPotential){
            /// TODO pridavam jen value_real? cas=0?
            DefaultEssentialBCConst *bc = new DefaultEssentialBCConst(str_marker[i], magneticEdge[i].value_real);
            bcs_mag.add_boundary_condition(bc);
        }
    }

    // Create an H1 space with default shapeset.
    H1Space space_mag_real(&mesh_mag, &bcs_mag, P_MAG_INIT);
    H1Space space_mag_imag(&mesh_mag, &bcs_mag, P_MAG_INIT);
    // ndof
    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&space_mag_real, &space_mag_imag));
    std::cout << "ndof: " << ndof << std::endl;

    //adaptivity ?????
    adapt_mesh(Hermes::vector<Space*>(&space_mag_real, &space_mag_imag), &wf);


    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Initialize the FE problem.
    DiscreteProblem dp(&wf, Hermes::vector<Space *>(&space_mag_real, &space_mag_imag));

    WjFilter wjfilter(&sln_mag_real, &sln_mag_imag);
    MagneticVectorPotentialFilter afilter(&sln_mag_real);

    // Visualize the solution.
    ScalarView view_a("Ar - real", new WinGeom(0, 0, 440, 750));
    ScalarView view_wj("wj", new WinGeom(450, 0, 440, 750));
    //view_wj.set_min_max_range(-0.000001, 0.000001);
 //   View::wait();


//****************** TEMPERATURE **********************************************

    WeakFormTemp wf_temp(TIME_STEP);
    wf_temp.registerForms(heatLabels, &sln_temp, &wjfilter);

    double current_time = 0;

    EssentialBCs bcs_temp;

    for(int i = 0; i < NUM_EDGES; i++){
        if(heatEdge[i].type == PhysicFieldBC_Heat_Temperature){
            DefaultEssentialBCConst *bc = new DefaultEssentialBCConst(str_marker[i], heatEdge[i].temperature);
            bcs_temp.add_boundary_condition(bc);
        }
    }

    // Create an H1 space with default shapeset.
    H1Space space_temp(&mesh_temp, &bcs_temp, P_TEMP_INIT);
    int ndof_temp = space_temp.get_num_dofs();
    info("temperature ndof = %d", ndof_temp);

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp_temp(&wf_temp, &space_temp, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix_temp = create_matrix(matrix_solver);
    Vector* rhs_temp = create_vector(matrix_solver);
    Solver* solver_temp = create_linear_solver(matrix_solver, matrix_temp, rhs_temp);
    solver_temp->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

    // Initialize views.
    ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));

//    //******************** Elasticita *******************************
    Solution sln_elast_r(&mesh_elast, DK_INIT);
    Solution sln_elast_z(&mesh_elast, DK_INIT);
    WeakFormElast wf_elast;
    wf_elast.register_forms(elasticityLabels, &sln_temp);

    EssentialBCs bcs_elast_r, bcs_elast_z;

    for(int i = 1; i < NUM_EDGES; i++){
        if(elasticityEdge[i].typeX == PhysicFieldBC_Elasticity_Fixed){
            DefaultEssentialBCConst *bc = new DefaultEssentialBCConst(str_marker[i], 0.);
            bcs_elast_r.add_boundary_condition(bc);
        }
        if(elasticityEdge[i].typeY == PhysicFieldBC_Elasticity_Fixed){
            DefaultEssentialBCConst *bc = new DefaultEssentialBCConst(str_marker[i], 0.);
            bcs_elast_z.add_boundary_condition(bc);
        }
    }

    H1Space space_elast_r(&mesh_elast, &bcs_elast_r, P_ELAST_INIT);
    H1Space space_elast_z(&mesh_elast, &bcs_elast_z, P_ELAST_INIT);
    int ndof_elast = Space::get_num_dofs(Hermes::vector<Space*>(&space_elast_r, &space_elast_z));
    info("elasticity ndof = %d", ndof_elast);

    // Initialize the FE problem.
    is_linear = true;
    DiscreteProblem dp_elast(&wf_elast, Hermes::vector<Space*>(&space_elast_r, &space_elast_z), is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix_elast = create_matrix(matrix_solver);
    Vector* rhs_elast = create_vector(matrix_solver);
    Solver* solver_elast = create_linear_solver(matrix_solver, matrix_elast, rhs_elast);
    solver_elast->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

    // Initialize views.
    ScalarView view_elast_r("Elasticity r", new WinGeom(400, 0, 450, 600));
    view_elast_r.show(&sln_elast_r);
    ScalarView view_elast_z("Elasticity z", new WinGeom(800, 0, 450, 600));
    view_elast_z.show(&sln_elast_z);
    // Visualize the solution.
//    WinGeom* stress_win_geom = new WinGeom(0, 0, 800, 400);
//    ScalarView stress_view("Von Mises stress [Pa]", stress_win_geom);
//    VonMisesFilter stress_filter(Hermes::vector<MeshFunction*>(&sln_elast_r, &sln_elast_z), elasticityLabel[0].lambda(), elasticityLabel[0].mu());
//    stress_view.show_mesh(false);
    ScalarView disp_view("Displacement", new WinGeom(0, 0, 800, 400));
    DisplacementFilter disp_filter(Hermes::vector<MeshFunction*>(&sln_elast_r, &sln_elast_z));

    GnuplotGraph temp_graph_time("Temperature/time", "time", "temperature");
    temp_graph_time.add_row("inner","k", "-");
    temp_graph_time.add_row("middle","k", "--");
    temp_graph_time.add_row("outer","k", ":");

    GnuplotGraph deformation_graph_time("Deformation/time", "time", "radial deformation");
    deformation_graph_time.add_row("inner","k", "-");
    deformation_graph_time.add_row("middle","k", "--");
    deformation_graph_time.add_row("outer","k", ":");

    // Time stepping:
    int ts = 1;
    do
    {
      info("---- Time step %d, time %3.5f s", ts, current_time);
      maxB = -10000; max_el_cond = -1000; min_el_cond = 1e10;

      if(ts == 2)
          wf.push_previous_temperature(&sln_temp);

      info("Assembling the magnetic stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);

      if (solver->solve())
          Solution::vector_to_solutions(solver->get_solution(),
                     Hermes::vector<Space *>(&space_mag_real, &space_mag_imag),
                     Hermes::vector<Solution *>(&sln_mag_real, &sln_mag_imag));
      else
          error ("Matrix solver failed.\n");

      info("max B %lf, min el cond %lf, max el cond %lf",maxB, min_el_cond, max_el_cond);

      view_a.show(&afilter, HERMES_EPS_NORMAL);
      view_wj.show(&wjfilter, HERMES_EPS_NORMAL);

      info("Assembling the temperature stiffness matrix and right-hand side vector.");
      dp_temp.assemble(matrix_temp, rhs_temp);

      info("Solving the temperature matrix problem.");
      if(solver_temp->solve())
          Solution::vector_to_solution(solver_temp->get_solution(), &space_temp, &sln_temp);
      else error ("Matrix solver failed.\n");

      // Visualize the solution.
      char title[100];
      sprintf(title, "Temperature, Time %3.2f s", current_time);
      Tview.set_title(title);
      Tview.show(&sln_temp, HERMES_EPS_NORMAL, H2D_FN_VAL_0);
    //  Tview.wait();

      info("Assembling the elastic stiffness matrix and right-hand side vector.");
      dp_elast.assemble(matrix_elast, rhs_elast);

      info("Solving the elasticity matrix problem.");
      if(solver_elast->solve())
          Solution::vector_to_solutions(solver_elast->get_solution(),
                Hermes::vector<Space*>(&space_elast_r, &space_elast_z), Hermes::vector<Solution*>(&sln_elast_r, &sln_elast_z));
      else error ("Matrix solver failed.\n");

      // Visualize the solution.
      sprintf(title, "Time %3.2f s, elast r", current_time);
      view_elast_r.set_title(title);
      view_elast_r.show(&sln_elast_r);

      sprintf(title, "Time %3.2f s, elast z", current_time);
      view_elast_z.set_title(title);
      view_elast_z.show(&sln_elast_z);

//      stress_view.show(&stress_filter, HERMES_EPS_LOW, H2D_FN_VAL_0, &sln_elast_r, &sln_elast_z, 1.5e5);
//      disp_view.show(&disp_filter);//, HERMES_EPS_LOW, H2D_FN_VAL_0, &sln_elast_r, &sln_elast_z, 1.5e5);
      disp_view.show(&disp_filter, HERMES_EPS_HIGH, H2D_FN_VAL_0, &sln_elast_r, &sln_elast_z, 5e2);

      temp_graph_time.add_values(0, current_time, sln_temp.get_pt_value(0.06, 0));
      temp_graph_time.add_values(1, current_time, sln_temp.get_pt_value(0.245, 0));
      temp_graph_time.add_values(2, current_time, sln_temp.get_pt_value(0.43, 0));
      temp_graph_time.save("results/temp_time.gnu");

      deformation_graph_time.add_values(0, current_time, sln_elast_r.get_pt_value(0.06, 0));
      deformation_graph_time.add_values(1, current_time, sln_elast_r.get_pt_value(0.245, 0));
      deformation_graph_time.add_values(2, current_time, sln_elast_r.get_pt_value(0.43, 0));
      deformation_graph_time.save("results/deformation_time.gnu");

      GnuplotGraph loses_graph_axis("Joule loses on axis", "r", "temperature");
      loses_graph_axis.add_row();

      GnuplotGraph deformation_graph_axis("Radial deformation on axis", "r", "temperature");
      deformation_graph_axis.add_row();

      GnuplotGraph temp_graph_axis("Temperature on axis", "r", "temperature");
      temp_graph_axis.add_row();

      for (double rr = 0.06; rr <= 0.43; rr+=0.005){
          temp_graph_axis.add_values(0, rr, sln_temp.get_pt_value(rr, 0));
          //loses_graph_axis.add_values(0, rr, wjfilter.get_pt_value(rr, 0, H2D_FN_VAL_0));
          deformation_graph_axis.add_values(0, rr, sln_elast_r.get_pt_value(rr, 0));
      }
      char filename[100];
      sprintf(filename, "results/temp_axis_%3.1lf.gnu", current_time);
      temp_graph_axis.save(filename);
      //loses_graph_axis.save_numbered("results/loses_axis.gnu", ts);
      sprintf(filename, "results/deformation_axis_%3.1lf.gnu", current_time);
      deformation_graph_axis.save(filename);


      // Increase current time and time step counter.
      current_time += TIME_STEP;
      ts++;
    }
    while (current_time < TIME_FINAL);

    Tview.wait();

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;

    return 0;
}
