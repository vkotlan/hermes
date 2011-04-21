#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

#include <iostream>

// A massive hollow conductor is heated by induction and cooled by water running inside.
// We will model this problem using linear thermoelasticity equations, where the x-displacement,
// y-displacement, and the temperature will be approximated on individual meshes equipped
// with mutually independent adaptivity mechanisms. Use MULTI = true to use multimesh,
// MULTI = false for single-mesh (all solution components on the samemesh).
//
// PDE: Linear thermoelasticity.
//
// BC: u_1 = u_2 = 0 on Gamma_1
//     du_1/dn = du_2/dn = 0 elsewhere
//     temp = TEMP_INNER on Gamma_4
//     negative heat flux with HEAT_FLUX_OUTER elsewhere.

const int P_INIT_TEMP = 3;
const int P_INIT_MAG = 2;
const int P_INIT_ELAST = 2;

using namespace RefinementSelectors;

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
// If false, projection of the fine mesh solution on the coarse mesh is used.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
// See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
// cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 3.0;            // Stopping criterion for adaptivity (rel. error tolerance between the
// fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.

// Global time variable.
double TIME = 0;

// Problem parameters.
const double A_INIT = 0.0;
const double T_INIT = 20.0;
const double DK_INIT = 0.0;

// Weak forms.
#include "forms.cpp"

class DisplacementFilter : public Filter
{
public:
    DisplacementFilter(Hermes::vector<MeshFunction*> solutions)
        : Filter(solutions)
    {

    }

    inline void precalculate(int order, int mask)
    {
        Quad2D* quad = quads[cur_quad];
        int np = quad->get_num_points(order);
        Node* node = new_node(H2D_FN_VAL_0, np);

        sln[0]->set_quad_order(order, H2D_FN_VAL);
        sln[1]->set_quad_order(order, H2D_FN_VAL);

        scalar *uval = sln[0]->get_fn_values();
        scalar *vval = sln[1]->get_fn_values();
        update_refmap();
        double *x = refmap->get_phys_x(order);

        for (int i = 0; i < np; i++)
        {
            // node->values[0][0][i] = sqrt(sqr(uval[i]) + sqr(vval[i]));
            node->values[0][0][i] = uval[i];
        }

        // remove the old node and attach the new one
        replace_cur_node(node);
    }

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); }
};

class WjFilter : public Filter
{
public:
    WjFilter(MeshFunction* slnr, MeshFunction* slni)
    {
        num = 2;

        sln[0] = slnr;
        sln[1] = slni;

        init();
    }

    inline void precalculate(int order, int mask)
    {
        Quad2D* quad = quads[cur_quad];
        int np = quad->get_num_points(order);
        Node* node = new_node(H2D_FN_DEFAULT, np);

        double *dudx1, *dudy1, *dudx2, *dudy2;
        double *value1, *value2;

        sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[0]->get_dx_dy_values(dudx1, dudy1);
        value1 = sln[0]->get_fn_values();

        sln[1]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[1]->get_dx_dy_values(dudx2, dudy2);
        value2 = sln[1]->get_fn_values();

        update_refmap();

        double *x = refmap->get_phys_x(order);
        double *y = refmap->get_phys_y(order);
        Element *e = refmap->get_active_element();

        for (int i = 0; i < np; i++)
        {
            node->values[0][0][i] = (magneticLabel[e->marker].conductivity > 0.0) ?
                        0.5 / magneticLabel[e->marker].conductivity * (
                            sqr(2 * M_PI * frequency * magneticLabel[e->marker].conductivity * value2[i]) +
                            sqr(2 * M_PI * frequency * magneticLabel[e->marker].conductivity * value1[i]))
                      :
                        0.0;
        }

        replace_cur_node(node);
    }

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); }
};

void set_magnetic_edge(Hermes::vector<MagneticEdge *> magneticEdge, PhysicFieldBC type, double value_real, double value_imag)
{
    for (int i = 0; i<magneticEdge.size(); i++)
    {
        magneticEdge[i]->type = type;
        magneticEdge[i]->value_real = value_real;
        magneticEdge[i]->value_imag = value_imag;
    }
}

void set_magnetic_label(Hermes::vector<MagneticLabel *> magneticLabel, double current_density_real, double current_density_imag, double permeability, double conductivity,
                        double remanence, double remanence_angle, double velocity_x, double velocity_y, double velocity_angular)
{
    for (int i = 0; i<magneticLabel.size(); i++)
    {
        magneticLabel[i]->current_density_real = current_density_real;
        magneticLabel[i]->current_density_imag = current_density_imag;
        magneticLabel[i]->permeability = permeability;
        magneticLabel[i]->conductivity = conductivity;
        magneticLabel[i]->remanence = remanence;
        magneticLabel[i]->remanence_angle = remanence_angle;
        magneticLabel[i]->velocity_x = velocity_x;
        magneticLabel[i]->velocity_y = velocity_y;
        magneticLabel[i]->velocity_angular = velocity_angular;
    }
}

void set_heat_edge(Hermes::vector<HeatEdge *> heatEdge, PhysicFieldBC type, double temperature, double heatFlux, double h, double externalTemperature)
{
    for (int i = 0; i<heatEdge.size(); i++)
    {
        heatEdge[i]->type = type;
        heatEdge[i]->temperature = temperature;
        heatEdge[i]->heatFlux = heatFlux;
        heatEdge[i]->h = h;
        heatEdge[i]->externalTemperature = externalTemperature;
    }
}

void set_heat_label(Hermes::vector<HeatLabel *> heatLabel, double thermal_conductivity, double volume_heat, double density, double specific_heat)
{
    for (int i = 0; i<heatLabel.size(); i++)
    {
        heatLabel[i]->thermal_conductivity = thermal_conductivity;
        heatLabel[i]->volume_heat = volume_heat;
        heatLabel[i]->density = density;
        heatLabel[i]->specific_heat = specific_heat;
    }
}

void set_elasticity_edge(Hermes::vector<ElasticityEdge *> elasticityEdge, PhysicFieldBC typeX, PhysicFieldBC typeY, double forceX, double forceY)
{
    for (int i = 0; i<elasticityEdge.size(); i++)
    {
        elasticityEdge[i]->typeX = typeX;
        elasticityEdge[i]->typeY = typeY;
        elasticityEdge[i]->forceX = forceX;
        elasticityEdge[i]->forceY = forceY;
    }
}

void set_elasticity_label(Hermes::vector<ElasticityLabel *> elasticityLabel, double young_modulus, double poisson_ratio, double forceX, double forceY, double thermal_expansion)
{
    for (int i = 0; i<elasticityLabel.size(); i++)
    {
        elasticityLabel[i]->young_modulus = young_modulus;
        elasticityLabel[i]->poisson_ratio = poisson_ratio;
        elasticityLabel[i]->forceX = forceX;
        elasticityLabel[i]->forceY = forceY;
        elasticityLabel[i]->thermal_expansion = thermal_expansion;
    }
}

void initTables()
{
    // thermal conductivity
    double temp_thermal_conductivity[] = { 0, 100, 200, 300, 400, 500, 600, 2000 };
    double data_thermal_conductivity[] = { 54, 49.4, 46.9, 43.9, 40.2, 37.2, 33.5, 23.0 };
    thermal_conductivity_fe.add(temp_thermal_conductivity, data_thermal_conductivity, 8);

    // electric conductivity
    double temp_electric_conductivity[] = { 0, 20, 100, 200, 400, 500, 600, 2000 };
    double data_electric_conductivity[] = { 5.5e6, 5.02e6, 4.08e6, 3.14e6, 1.96e6, 1.65e6, 1.31e6, 0.9e6 };
    electric_conductivity_fe.add(temp_electric_conductivity, data_electric_conductivity, 8);

    // relative permeabilty
    double temp_relative_mag_permeability_temp[] = { 0, 100, 200, 300, 400, 500, 600, 700, 800, 10000 };
    double data_relative_mag_permeability_temp[] = { 2600, 2600, 2598, 2595, 2500, 2100, 1490, 750, 2, 1 };
    relative_mag_permeability_temp.add(temp_relative_mag_permeability_temp, data_relative_mag_permeability_temp, 10);

    double temp_relative_mag_permeability[] = { 1.2566370614359172e-07, 0.2402, 0.8654, 1.1106, 1.2458, 1.331, 1.5, 1.6, 1.683, 1.741, 1.78, 1.905, 2.025, 2.085, 2.13, 2.165, 2.28, 2.485, 2.585 };
    double data_relative_mag_permeability[] = { 1000.0, 1202.1703563104797, 2165.6082979831167, 1852.803771466027, 1558.7675165399623, 1332.2970393415892, 750.25900263307039, 400.01242373080828, 280.53809093387082, 217.63175928604289, 178.01671402763208, 95.25295840089872, 50.624981898320513, 34.750351479349241, 26.624990479857455, 21.650128290457886, 11.399995923769952, 6.2124977786334066, 5.1700046482110427 };
    relative_mag_permeability.add(temp_relative_mag_permeability, data_relative_mag_permeability, 10);
}

int main(int argc, char* argv[])
{
    initTables();

    // harmonic elmag. field
    magneticEdge = new MagneticEdge[42];
    set_magnetic_edge(Hermes::vector<MagneticEdge *>(&magneticEdge[0]), PhysicFieldBC_None, 0.0, 0.0);
    set_magnetic_edge(Hermes::vector<MagneticEdge *>(&magneticEdge[16], &magneticEdge[17], &magneticEdge[18],
                                            &magneticEdge[39], &magneticEdge[40], &magneticEdge[41]), PhysicFieldBC_Magnetic_VectorPotential, 0.0, 0.0);

    magneticLabel = new MagneticLabel[7];
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[0]), 0.0, 0.0, 700.0, 1e6, 0.0, 0.0, 0.0, 0.0, 0.0);
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[1]), 6e7, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[2]), 6e7, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[3]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[4]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[5]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[6]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    // heat
//    heatEdge = new HeatEdge[42];
//    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[0]), PhysicFieldBC_None, 0.0, 0.0, 0.0, 0.0);
//    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[3], &heatEdge[2]), PhysicFieldBC_Heat_Temperature, 20.0, 0.0, 0.0, 0.0);
//    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[10], &heatEdge[12], &heatEdge[13], &heatEdge[14], &heatEdge[15], &heatEdge[1]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 10.0, 20.0);
//    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[4], &heatEdge[5], &heatEdge[6], &heatEdge[7], &heatEdge[8], &heatEdge[9]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 20.0, 20.0);
//    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[35], &heatEdge[36], &heatEdge[37], &heatEdge[38]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 0.0, 0.0);
//    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[40]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 0.0, 0.0);
//    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[11]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 10.0, 20.0);

//    heatLabel = new HeatLabel[7];
//    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[0]), 19.9, 1.0, 8670.0, 450);
//    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[6]), 1e-6, 0.0, 0.0, 0.0);

    // thermoelasticity
//    elasticityEdge = new ElasticityEdge[42];
//    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[0]), PhysicFieldBC_None, PhysicFieldBC_None, 0.0, 0.0);
//    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[2], &elasticityEdge[3]), PhysicFieldBC_Elasticity_Fixed, PhysicFieldBC_Elasticity_Fixed, 0.0, 0.0);
//    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[11], &elasticityEdge[12], &elasticityEdge[13], &elasticityEdge[14], &elasticityEdge[15], &elasticityEdge[1]), PhysicFieldBC_Elasticity_Free, PhysicFieldBC_Elasticity_Free, 0.0, 0.0);
//    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[4], &elasticityEdge[5], &elasticityEdge[6], &elasticityEdge[7], &elasticityEdge[8], &elasticityEdge[9], &elasticityEdge[10]), PhysicFieldBC_Elasticity_Free, PhysicFieldBC_Elasticity_Free, 0.0, 0.0);
//    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[35], &elasticityEdge[36], &elasticityEdge[37]), PhysicFieldBC_Elasticity_Free, PhysicFieldBC_Elasticity_Free, 0.0, 0.0);
//    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[38]), PhysicFieldBC_Elasticity_Fixed, PhysicFieldBC_Elasticity_Fixed, 0.0, 0.0);
//    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[40]), PhysicFieldBC_Elasticity_Fixed, PhysicFieldBC_Elasticity_Free, 0.0, 0.0);

//    elasticityLabel = new ElasticityLabel[7];
//    set_elasticity_label(Hermes::vector<ElasticityLabel *>(&elasticityLabel[0]), 200e9, 0.28, 0.0, 0.0, 9.7e-6);
//    set_elasticity_label(Hermes::vector<ElasticityLabel *>(&elasticityLabel[6]), 200e9, 0.25, 0.0, 0.0, 13e-6);

    // Time measurement.
    TimePeriod cpu_time;
    cpu_time.tick();

    // Load the mesh.
//    Mesh mesh_mag, mesh_temp, mesh_elastDrzak, mesh_elastStopka;
    Mesh mesh_mag;
    H2DReader mloader;
    mloader.load("mesh_mag.mesh", &mesh_mag);

//    mloader.load("mesh_temp.mesh", &mesh_temp);
//    mloader.load("mesh_elast_drzak.mesh", &mesh_elastDrzak);
//    mloader.load("mesh_elast_stopka.mesh", &mesh_elastStopka);
//    mesh_temp.refine_all_elements(0);
    mesh_mag.refine_all_elements();
//    mesh_elastDrzak.refine_all_elements(0);
//    mesh_elastStopka.refine_all_elements(0);

    // Create H1 spaces with default shapesets.
    H1Space space_mag_real(&mesh_mag, magnetic_bc_types, magnetic_bc_values_real, P_INIT_MAG);
    H1Space space_mag_imag(&mesh_mag, magnetic_bc_types, magnetic_bc_values_imag, P_INIT_MAG);
//    H1Space space_temp(&mesh_temp, heat_bc_types, heat_bc_values, P_INIT_TEMP);
//    H1Space space_elast_drzak_r(&mesh_elastDrzak, elasticity_bc_types_r, elasticity_bc_values_r, P_INIT_ELAST);
//    H1Space space_elast_drzak_z(&mesh_elastDrzak, elasticity_bc_types_z, elasticity_bc_values_z, P_INIT_ELAST);
//    H1Space space_elast_stopka_r(&mesh_elastStopka, elasticity_bc_types_r, elasticity_bc_values_r, P_INIT_ELAST);
//    H1Space space_elast_stopka_z(&mesh_elastStopka, elasticity_bc_types_z, elasticity_bc_values_z, P_INIT_ELAST);

    // ndof
//    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&space_mag_real,
//                                                  &space_mag_imag,
//                                                  &space_temp,
//                                                  &space_elast_drzak_r,
//                                                  &space_elast_drzak_z,
//                                                  &space_elast_stopka_r,
//                                                  &space_elast_stopka_z));

    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&space_mag_real,
                                                      &space_mag_imag));

    std::cout << "ndofs: " << ndof << std::endl;

    // Set initial condition.
    Solution *sln_mag_real = new Solution();
    sln_mag_real->set_const(&mesh_mag, A_INIT);
    Solution *sln_mag_imag = new Solution();
    sln_mag_imag->set_const(&mesh_mag, A_INIT);
//    Solution *sln_temp = new Solution();
//    sln_temp->set_const(&mesh_temp, T_INIT);
//    Solution *sln_elast_drzak_r = new Solution();
//    sln_elast_drzak_r->set_const(&mesh_elastDrzak, DK_INIT);
//    Solution *sln_elast_drzak_z = new Solution();
//    sln_elast_drzak_z->set_const(&mesh_elastDrzak, DK_INIT);
//    Solution *sln_elast_stopka_r = new Solution();
//    sln_elast_stopka_r->set_const(&mesh_elastStopka, DK_INIT);
//    Solution *sln_elast_stopka_z = new Solution();
//    sln_elast_stopka_z->set_const(&mesh_elastStopka, DK_INIT);

    WjFilter wjfilter(sln_mag_real, sln_mag_imag);

    // Initialize views.
    char title[100];
//    ScalarView TView("Temperature", new WinGeom(450, 0, 300, 900));
    ScalarView TView("Magnetic field", new WinGeom(450, 0, 300, 900));
//    sprintf(title, "Time %3.5f", TIME);
//    // Tview.set_min_max_range(16.0348, 16.0698);
//    TView.set_title(title);

//    ScalarView WjView("Joule losses", new WinGeom(800, 0, 400, 900));
//    ScalarView DisplacementDrzakView("Von Mises Stress", new WinGeom(1250, 0, 300, 900));
//    DisplacementDrzakView.set_min_max_range(-1e-6, 30e-6);
//    ScalarView DisplacementStopkaView("Von Mises Stress", new WinGeom(1650, 0, 200, 900));
//    DisplacementStopkaView.set_min_max_range(-1e-6, 5e-6);
//    ScalarView DispRView("Displacement r", new WinGeom(1250, 0, 300, 900));
//    ScalarView DispZView("Displacement z", new WinGeom(1600, 0, 300, 900));
//    // Tview.set_min_max_range(16.0348, 16.0698);

    int N = 20;
    double yy_total[N], yy2_total[N];
    for (int k = 0; k < N; k++)
    {
        yy_total[k] = 0.0;
        yy2_total[k] = 0.0;
    }

    // Time stepping:
    int nsteps = (int)(FINAL_TIME/timeStep + 0.5);
    bool rhsonly = false;
    for(int ts = 1; ts <= nsteps; ts++)
    {
        info("---- Time step %d, time %3.5f", ts, TIME);

        if (ts == 29)
        {
            int N = 30;
            for (int k = 0; k < N; k++)
            {
                double xx = 0.004;
                double yy = 0.0 + k * 0.0135/(double) N;

                //tempBreak.add(yy, sln_temp->get_pt_value(xx, yy));
            }

//            set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[11], &heatEdge[35], &heatEdge[36], &heatEdge[37], &heatEdge[38]), PhysicFieldBC_Heat_Temperature, 20.0, 0.0, 0.0, 0.0);

//            set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[0]), 19.9, 0.0, 1000.0, 1e6);
//            set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[6]), 1e-6, 0.0, 0.0, 0.0);
        }
        if (ts == 30)
        {
//            set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[35], &heatEdge[36], &heatEdge[37]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 10.0, 20.0);
//            set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[8], &heatEdge[9], &heatEdge[10], &heatEdge[37], &heatEdge[38]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 150.0, 10.0);
//            set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[11]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 0.0, 0.0);

//            set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[0]), 19.9, 0.0, 8670.0, 450);
//            set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[6]), 60, 0.0, 7800.0, 450.0);
        }

        // Initialize the weak formulation.
        WeakForm wf(2);
        // Magnetic field
        wf.add_matrix_form(0, 0, magnetic_matrix_form_linear_real_real, magnetic_matrix_form_linear_ord, HERMES_NONSYM, HERMES_ANY); //, Hermes::vector<MeshFunction *>(sln_temp));
        wf.add_matrix_form(0, 1, callback(magnetic_matrix_form_linear_real_imag));
        wf.add_matrix_form(1, 0, callback(magnetic_matrix_form_linear_imag_real));
        wf.add_matrix_form(1, 1, magnetic_matrix_form_linear_imag_imag, magnetic_matrix_form_linear_ord, HERMES_NONSYM, HERMES_ANY); //, Hermes::vector<MeshFunction *>(sln_temp));
        wf.add_vector_form(0, callback(magnetic_vector_form_linear_real));
        wf.add_vector_form(1, callback(magnetic_vector_form_linear_imag));
        wf.add_vector_form_surf(0, callback(magnetic_vector_form_linear_surf_real));
        wf.add_vector_form_surf(1, callback(magnetic_vector_form_linear_surf_imag));

//        // Temperature field
//        wf.add_matrix_form(2, 2, heat_matrix_form_linear, heat_matrix_form_linear_ord, HERMES_SYM, HERMES_ANY);
//        wf.add_vector_form(2, heat_vector_form_linear, heat_vector_form_linear_ord, HERMES_ANY, Hermes::vector<MeshFunction *>(&wjfilter, sln_temp));
//        wf.add_matrix_form_surf(2, 2, callback(heat_matrix_form_linear_surf));
//        wf.add_vector_form_surf(2, callback(heat_vector_form_linear_surf));

//        // Thermoelasticity - drzak
//        wf.add_matrix_form(3, 3, callback(elasticity_matrix_form_r_r));
//        wf.add_matrix_form(3, 4, callback(elasticity_matrix_form_r_z), HERMES_SYM);
//        wf.add_matrix_form(4, 4, callback(elasticity_matrix_form_z_z));
//        wf.add_matrix_form(3, 2, callback(elasticity_matrix_form_r_T));
//        wf.add_matrix_form(4, 2, callback(elasticity_matrix_form_z_T));
//        wf.add_vector_form(3, elasticity_vector_form_r, heat_vector_form_linear_ord, HERMES_ANY, Hermes::vector<MeshFunction *>(sln_temp));
//        wf.add_vector_form(4, elasticity_vector_form_z, heat_vector_form_linear_ord, HERMES_ANY, Hermes::vector<MeshFunction *>(sln_temp));

//        // Thermoelasticity - stopka
//        wf.add_matrix_form(5, 5, callback(elasticity_matrix_form_r_r));
//        wf.add_matrix_form(5, 6, callback(elasticity_matrix_form_r_z), HERMES_SYM);
//        wf.add_matrix_form(6, 6, callback(elasticity_matrix_form_z_z));
//        wf.add_matrix_form(5, 2, callback(elasticity_matrix_form_r_T));
//        wf.add_matrix_form(6, 2, callback(elasticity_matrix_form_z_T));
//        wf.add_vector_form(5, elasticity_vector_form_r, heat_vector_form_linear_ord, HERMES_ANY, Hermes::vector<MeshFunction *>(sln_temp));
//        wf.add_vector_form(6, elasticity_vector_form_z, heat_vector_form_linear_ord, HERMES_ANY, Hermes::vector<MeshFunction *>(sln_temp));

        // Set up the solver, matrix, and rhs according to the solver selection.
        SparseMatrix* matrix = create_matrix(SOLVER_UMFPACK);
        Vector* rhs = create_vector(SOLVER_UMFPACK);
        Solver* solver = create_linear_solver(SOLVER_UMFPACK, matrix, rhs);

        // Initialize the reference FE problem.
        bool is_linear = true;
//        DiscreteProblem fep(&wf, Hermes::vector<Space *>(&space_mag_real,
//                                                &space_mag_imag,
//                                                &space_temp,
//                                                &space_elast_drzak_r,
//                                                &space_elast_drzak_z,
//                                                &space_elast_stopka_r,
//                                                &space_elast_stopka_z), is_linear);

        DiscreteProblem fep(&wf, Hermes::vector<Space *>(&space_mag_real,
                                                &space_mag_imag), is_linear);


        // Assemble stiffness matrix and rhs.
        fep.assemble(matrix, rhs, rhsonly);
        // rhsonly = true;

        // Solve the matrix problem.
        if (!solver->solve()) error ("Matrix solver failed.\n");
//        Solution::vector_to_solutions(solver->get_solution(),
//                                      Hermes::vector<Space *>(&space_mag_real,
//                                                     &space_mag_imag,
//                                                     &space_temp,
//                                                     &space_elast_drzak_r,
//                                                     &space_elast_drzak_z,
//                                                     &space_elast_stopka_r,
//                                                     &space_elast_stopka_z),
//                                      Hermes::vector<Solution *>(sln_mag_real,
//                                                        sln_mag_imag,
//                                                        sln_temp,
//                                                        sln_elast_drzak_r,
//                                                        sln_elast_drzak_z,
//                                                        sln_elast_stopka_r,
//                                                        sln_elast_stopka_z));

        Solution::vector_to_solutions(solver->get_solution(),
                                      Hermes::vector<Space *>(&space_mag_real,
                                                     &space_mag_imag),
                                      Hermes::vector<Solution *>(sln_mag_real,
                                                        sln_mag_imag));


        delete matrix;
        delete rhs;
        delete solver;

        /*
        bool done = false;
        do
        {
            // Construct globally refined reference mesh and setup reference space.
            Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(Hermes::vector<Space *>(&space_mag_real,
                                                                                 &space_mag_imag,
                                                                                 &space_temp,
                                                                                 &space_elast_kotva_r,
                                                                                 &space_elast_kotva_z,
                                                                                 &space_elast_uchyt_r,
                                                                                 &space_elast_uchyt_z));

            // Initialize the reference FE problem.
            bool is_linear = true;
            FeProblem fep(&wf, *ref_spaces, is_linear);

            // Set up the solver, matrix, and rhs according to the solver selection.
            SparseMatrix* matrix = create_matrix(SOLVER_UMFPACK);
            Vector* rhs = create_vector(SOLVER_UMFPACK);
            Solver* solver = create_linear_solver(SOLVER_UMFPACK, matrix, rhs);

            // Assemble stiffness matrix and rhs.
            fep.assemble(matrix, rhs, rhsonly);
            // rhsonly = true;

            // Solve the matrix problem.
            if (!solver->solve()) error ("Matrix solver failed.\n");

            Solution sln_mag_real_ref, sln_mag_imag_ref, sln_temp_ref, sln_elast_kotva_r_ref, sln_elast_kotva_z_ref, sln_elast_uchyt_r_ref, sln_elast_uchyt_z_ref;

            Solution::vector_to_solutions(solver->get_solution(), *ref_spaces,
                                          Hermes::vector<Solution *>(&sln_mag_real_ref,
                                                            &sln_mag_imag_ref,
                                                            &sln_temp_ref,
                                                            &sln_elast_kotva_r_ref,
                                                            &sln_elast_kotva_z_ref,
                                                            &sln_elast_uchyt_r_ref,
                                                            &sln_elast_uchyt_z_ref));

            project_global(Hermes::vector<Space *>(&space_mag_real,
                                          &space_mag_imag,
                                          &space_temp,
                                          &space_elast_kotva_r,
                                          &space_elast_kotva_z,
                                          &space_elast_uchyt_r,
                                          &space_elast_uchyt_z),
                           Hermes::vector<Solution *>(&sln_mag_real_ref,
                                             &sln_mag_imag_ref, &sln_temp_ref,
                                             &sln_elast_kotva_r_ref,
                                             &sln_elast_kotva_z_ref,
                                             &sln_elast_uchyt_r_ref,
                                             &sln_elast_uchyt_z_ref),
                           Hermes::vector<Solution *>(sln_mag_real,
                                             sln_mag_imag,
                                             sln_temp,
                                             sln_elast_kotva_r,
                                             sln_elast_kotva_z,
                                             sln_elast_uchyt_r,
                                             sln_elast_uchyt_z));


            // Calculate element errors.
            info("Calculating error estimate and exact error.");
            Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&space_mag_real,
                                                         &space_mag_imag,
                                                         &space_temp,
                                                         &space_elast_kotva_r,
                                                         &space_elast_kotva_z,
                                                         &space_elast_uchyt_r,
                                                         &space_elast_uchyt_z),
                                          Hermes::vector<ProjNormType>(HERMES_H1_NORM,
                                                              HERMES_H1_NORM,
                                                              HERMES_H1_NORM,
                                                              HERMES_H1_NORM,
                                                              HERMES_H1_NORM,
                                                              HERMES_H1_NORM,
                                                              HERMES_H1_NORM));

            adaptivity->set_solutions(Hermes::vector<Solution *>(sln_mag_real,
                                                        sln_mag_imag,
                                                        sln_temp,
                                                        sln_elast_kotva_r,
                                                        sln_elast_kotva_z,
                                                        sln_elast_uchyt_r,
                                                        sln_elast_uchyt_z),
                                      Hermes::vector<Solution *>(&sln_mag_real_ref,
                                                        &sln_mag_imag_ref,
                                                        &sln_temp_ref,
                                                        &sln_elast_kotva_r_ref,
                                                        &sln_elast_kotva_z_ref,
                                                        &sln_elast_uchyt_r_ref,
                                                        &sln_elast_uchyt_z_ref));

            // Calculate error estimate for each solution component and the total error estimate.
            Hermes::vector<double> err_est_rel;
            double err_est_rel_total = adaptivity->calc_err_est(err_est_rel,
                                       HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;

            // Initialize refinement selector.
            H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

            // If err_est too large, adapt the mesh.
            if (err_est_rel_total < ERR_STOP)
                done = true;
            else
            {
                info("Adapting coarse mesh.");
                done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector),
                                         THRESHOLD, STRATEGY, MESH_REGULARITY);
            }

            delete matrix;
            delete rhs;
            delete solver;
        }
        while (!done);
        */

        WjFilter wjfilter(sln_mag_real, sln_mag_imag);
//        DisplacementFilter dispDrzak(Hermes::vector<MeshFunction *>(sln_elast_drzak_r, sln_elast_drzak_z));
//        DisplacementFilter dispStopka(Hermes::vector<MeshFunction *>(sln_elast_stopka_r, sln_elast_stopka_z));

        // Update the time variable.
        TIME += timeStep;

        // Visualize the solution.
//        sprintf(title, "Time %3.2f", TIME);
        sprintf(title, "Mag_real", TIME);

        TView.set_title(title);
//        TView.show(sln_temp);
        TView.show(sln_mag_real);

        TView.set_min_max_range(20, 300);

        // WjView.show(&wjfilter);

        //DisplacementDrzakView.show(&dispDrzak, HERMES_EPS_HIGH, H2D_FN_VAL_0, sln_elast_drzak_r, sln_elast_drzak_z, 5e1);
        //DisplacementStopkaView.show(&dispStopka, HERMES_EPS_HIGH, H2D_FN_VAL_0, sln_elast_stopka_r, sln_elast_stopka_z, 3e1);
        // VonMisesKView.show(&stress, HERMES_EPS_HIGH);
        // VonMisesUchytView.show(&stressUchyt, HERMES_EPS_HIGH);
        // DispRView.show(sln_elast_kotva_r);
        // DispZView.show(sln_elast_kotva_z);

        /*
        if (ts % 1 == 0)
        {
            FILE *f1, *f2;
            char fn1[50];
            sprintf(fn1, "chart_%d.dat", ts);
            char fn2[50];
            sprintf(fn2, "chart_temp_%d.dat", ts);

            f1 = fopen(fn1, "w");
            f2 = fopen(fn2, "w");

            for (int k = 0; k < N; k++)
            {
                double xx = 0.004;
                double yy = 0.0 + k * 0.0135/(double) N;

                double u = sln_elast_r->get_pt_value(xx, yy);
                double T = sln_temp->get_pt_value(xx, yy);

                yy_total[k] += u;

                printf("%3.5e\t%3.5e\t%3.5e\t%3.5e\n", yy, T, u, yy_total[k]);

                fprintf(f1, "%3.5e\t\t%3.5e\n", yy, u); // yy_total[k]
                fprintf(f2, "%3.5e\t\t%3.5e\n", yy, T);
            }
            fclose(f1);
            fclose(f2);
            printf("\n");
        }
        */
    }

    // Show the Von Mises stress on the reference mesh.
    /*
    WinGeom *winVonMisesStress = new WinGeom(0, 0, 450, 350);
    ScalarView sview("Von Mises Stress", winVonMisesStress);

    VonMisesFilter ref_stress(Hermes::vector<MeshFunction*>(xdisp_sln, ydisp_sln), lambda, mu);
    sview.show_mesh(false);
    sview.show(&ref_stress, H2D_EPS_HIGH);
    */
    // Wait for all views to be closed.
    View::wait();
    return 0;
};
