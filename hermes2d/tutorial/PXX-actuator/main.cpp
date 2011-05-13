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
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double A_INIT = 0.0;
const double TEMP_INIT = 20.0;
const double DK_INIT = 0.0;

const double TIME_STEP = 0.1;
const double TIME_FINAL = 20.;

const double frequency = 5000;

std::string *str_marker;

#include "tables.cpp"

// Weak forms.
#include "definitions.cpp"

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
    mloader.load("act_mesh_mag.mesh", &mesh_mag);
    mloader.load("act_mesh_temp.mesh", &mesh_temp);
    mloader.load("act_mesh_elast.mesh", &mesh_elast);

    MeshView mv_mag;
    mv_mag.show(&mesh_mag);
    MeshView mv_temp;
    mv_temp.show(&mesh_temp);
    MeshView mv_elast;
    mv_elast.show(&mesh_elast);

    mv_elast.wait();

    // Perform initial mesh refinements (optional).
     for (int i=0; i < INIT_REF_NUM; i++){
         mesh_temp.refine_all_elements();
         mesh_mag.refine_all_elements();
   //      mesh_elast.refine_all_elements();
    }
    // Initialize the weak formulation.
    WeakFormMagnetic wf(2);
    wf.registerForms(magneticLabels);

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

    Solution *sln_mag_real = new Solution();
    sln_mag_real->set_const(&mesh_mag, A_INIT);
    Solution *sln_mag_imag = new Solution();
    sln_mag_imag->set_const(&mesh_mag, A_INIT);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Initialize the FE problem.
    DiscreteProblem dp(&wf, Hermes::vector<Space *>(&space_mag_real, &space_mag_imag));
    dp.assemble(matrix, rhs);

    if (solver->solve())
        Solution::vector_to_solutions(solver->get_solution(),
                                      Hermes::vector<Space *>(&space_mag_real, &space_mag_imag),
                                      Hermes::vector<Solution *>(sln_mag_real, sln_mag_imag));
    else
        error ("Matrix solver failed.\n");

    WjFilter wjfilter(sln_mag_real, sln_mag_imag);
    MagneticVectorPotentialFilter afilter(sln_mag_real);

    // Visualize the solution.
    ScalarView view_a("Ar - real", new WinGeom(0, 0, 440, 750));
    view_a.show(&afilter, HERMES_EPS_NORMAL);
    ScalarView view_wj("wj", new WinGeom(450, 0, 440, 750));
    //view_wj.set_min_max_range(-0.000001, 0.000001);
    view_wj.show(&wjfilter, HERMES_EPS_NORMAL);
 //   View::wait();


//****************** TEMPERATURE **********************************************

//    Solution* sln_temp = new Solution();
//    sln_temp->set_const(&mesh_temp, TEMP_INIT);
    Solution sln_temp(&mesh_temp, TEMP_INIT);
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
    Tview.show(&sln_temp);

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

    // Time stepping:
    int ts = 1;
    do
    {
      info("---- Time step %d, time %3.5f s", ts, current_time);

      info("Assembling the temperature stiffness matrix and right-hand side vector.");
      dp_temp.assemble(matrix_temp, rhs_temp);
      FILE *matfile;
      matfile = fopen("matice.txt", "w");
      matrix_temp->dump(matfile, "matrix");
      rhs_temp->dump(matfile, "rhs");

      info("Solving the temperature matrix problem.");
      if(solver_temp->solve())
          Solution::vector_to_solution(solver_temp->get_solution(), &space_temp, &sln_temp);
      else error ("Matrix solver failed.\n");

      // Visualize the solution.
      char title[100];
      sprintf(title, "Time %3.2f s", current_time);
      Tview.set_title(title);
      Tview.show(&sln_temp);
    //  Tview.wait();

      info("Assembling the elastic stiffness matrix and right-hand side vector.");
      dp_elast.assemble(matrix_elast, rhs_elast);
//      FILE *matfile;
//      matfile = fopen("matice.txt", "w");
//      matrix_temp->dump(matfile, "matrix");
//      rhs_temp->dump(matfile, "rhs");

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
