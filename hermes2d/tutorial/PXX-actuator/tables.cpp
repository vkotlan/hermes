#include <hermes2d.h>
#include "data_table.h"

const int NUM_EDGES = 61;
const int NUM_LABELS = 19;

DataTable thermal_conductivity_fe;
DataTable electric_conductivity_fe;
DataTable relative_mag_permeability_temp;
DataTable relative_mag_permeability;

DataTable tempBreak;

enum PhysicFieldBC
{
    PhysicFieldBC_Undefined,
    PhysicFieldBC_None,
    PhysicFieldBC_General_Value,
    PhysicFieldBC_General_Derivative,
    PhysicFieldBC_Electrostatic_Potential,
    PhysicFieldBC_Electrostatic_SurfaceCharge,
    PhysicFieldBC_Magnetic_VectorPotential,
    PhysicFieldBC_Magnetic_SurfaceCurrent,
    PhysicFieldBC_Heat_Temperature,
    PhysicFieldBC_Heat_Flux,
    PhysicFieldBC_Current_Potential,
    PhysicFieldBC_Current_InwardCurrentFlow,
    PhysicFieldBC_Elasticity_Fixed,
    PhysicFieldBC_Elasticity_Free,
    PhysicFieldBC_Flow_Velocity,
    PhysicFieldBC_Flow_Pressure,
    PhysicFieldBC_Flow_Outlet,
    PhysicFieldBC_Flow_Wall,
};

double permeability_function(double B, double T)
{
    double koef = 1.25350233e-28 * pow(T, 10) - 9.38795536e-25 * pow(T, 9) + 2.95485549e-21 * pow(T, 8) - 5.06321208e-18 * pow(T, 7) +
            5.11236483e-15 * pow(T, 6) - 3.08097605e-12 * pow(T, 5) + 1.08123707e-09 * pow(T, 4) - 2.12158135e-07 * pow(T, 3) +
            2.11535604e-05 * pow(T, 2) - 8.14585198e-04 * pow(T, 1) + 1.00001549e+00;

    double perm = relative_mag_permeability.value(B);

    perm = perm * koef;
    if (perm < 1.0) perm = 1.0;

    return perm;
}

// ************************************************************************************************************

// elmag. field
struct MagneticEdge
{
    PhysicFieldBC type;
    double value_real;
    double value_imag;
};

struct MagneticLabel
{
    double current_density_real;
    double current_density_imag;
    double permeability;
    double conductivity;
    double remanence;
    double remanence_angle;
    double velocity_x;
    double velocity_y;
    double velocity_angular;
};

MagneticEdge *magneticEdge;
MagneticLabel *magneticLabel;
Hermes::vector<int> magneticLabels;

//BCType magnetic_bc_types(int marker)
//{
//    switch (magneticEdge[marker].type)
//    {
//    case PhysicFieldBC_None:
//        return BC_NONE;
//    case PhysicFieldBC_Magnetic_VectorPotential:
//        return BC_ESSENTIAL;
//    case PhysicFieldBC_Magnetic_SurfaceCurrent:
//        return BC_NATURAL;
//    }
//}

scalar magnetic_bc_values_real(int marker, double x, double y)
{
    return  0.0;
    info("vnd = %2.5f", magneticEdge[marker].value_real);
    return magneticEdge[marker].value_real;
}

scalar magnetic_bc_values_imag(int marker, double x, double y)
{
    return magneticEdge[marker].value_imag;
}


// ************************************************************************************************************

// temperature
struct HeatEdge
{
    PhysicFieldBC type;
    double temperature;
    double heatFlux;
    double h;
    double externalTemperature;
};

struct HeatLabel
{
    double thermal_conductivity;
    double volume_heat;
    double density;
    double specific_heat;
};

HeatEdge *heatEdge;
HeatLabel *heatLabel;
Hermes::vector<int> heatLabels;
Hermes::vector<int> zelezoLabels;
Hermes::vector<int> vodiveZelezoMosazLabels;

//BCType heat_bc_types(int marker)
//{
//    // info("heat_bc_types");
//
//    switch (heatEdge[marker].type)
//    {
//    case PhysicFieldBC_None:
//        return BC_NONE;
//    case PhysicFieldBC_Heat_Temperature:
//        return BC_ESSENTIAL;
//    case PhysicFieldBC_Heat_Flux:
//        return BC_NATURAL;
//    }
//}

scalar heat_bc_values(int marker, double x, double y)
{
    // info("heat_bc_values");

    switch (heatEdge[marker].type)
    {
    case PhysicFieldBC_Heat_Temperature:
        if (marker == 11)
        {
            // break time
            return (tempBreak.value(y) + 20.0) * 60.0 / (60.0 + 19.9);
        }
        else
        {
            return heatEdge[marker].temperature;
        }
    case PhysicFieldBC_Heat_Flux:
        return heatEdge[marker].heatFlux;
    }
    assert(0);
}


// ************************************************************************************************************

struct ElasticityEdge
{
public:
    PhysicFieldBC typeX;
    PhysicFieldBC typeY;
    double forceX;
    double forceY;
};

struct ElasticityLabel
{
    double young_modulus;
    double poisson_ratio;
    double forceX;
    double forceY;
    double thermal_expansion;

    // Lame constant
    inline double lambda() { return (young_modulus * poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0*poisson_ratio)); }
    inline double mu() { return young_modulus / (2.0*(1.0 + poisson_ratio)); }
};

ElasticityEdge *elasticityEdge;
ElasticityLabel *elasticityLabel;
Hermes::vector<int> elasticityLabels;

//BCType elasticity_bc_types_r(int marker)
//{
//    switch (elasticityEdge[marker].typeX)
//    {
//    case PhysicFieldBC_None:
//        return BC_NONE;
//        break;
//    case PhysicFieldBC_Elasticity_Fixed:
//        return BC_ESSENTIAL;
//        break;
//    case PhysicFieldBC_Elasticity_Free:
//        return BC_NATURAL;
//        break;
//    }
//}

//BCType elasticity_bc_types_z(int marker)
//{
//    switch (elasticityEdge[marker].typeY)
//    {
//    case PhysicFieldBC_None:
//        return BC_NONE;
//        break;
//    case PhysicFieldBC_Elasticity_Fixed:
//        return BC_ESSENTIAL;
//        break;
//    case PhysicFieldBC_Elasticity_Free:
//        return BC_NATURAL;
//        break;
//    }
//}

scalar elasticity_bc_values_r(int marker, double x, double y)
{
    switch (elasticityEdge[marker].typeX)
    {
    case PhysicFieldBC_None:
        return 0;
        break;
    case PhysicFieldBC_Elasticity_Fixed:
        return 0;
        break;
    case PhysicFieldBC_Elasticity_Free:
        return elasticityEdge[marker].forceX;
        break;
    }
    assert(0);
}

scalar elasticity_bc_values_z(int marker, double x, double y)
{
    switch (elasticityEdge[marker].typeY)
    {
    case PhysicFieldBC_None:
        return 0;
        break;
    case PhysicFieldBC_Elasticity_Fixed:
        return 0;
        break;
    case PhysicFieldBC_Elasticity_Free:
        return elasticityEdge[marker].forceY;
        break;
    }
    assert(0);
}

void set_magnetic_edge(Hermes::vector<MagneticEdge *> magneticEdge, PhysicFieldBC type, double value_real, double value_imag)
{
    for (unsigned int i = 0; i<magneticEdge.size(); i++)
    {
        magneticEdge[i]->type = type;
        magneticEdge[i]->value_real = value_real;
        magneticEdge[i]->value_imag = value_imag;
    }
}

void set_magnetic_label(Hermes::vector<MagneticLabel *> magneticLabel, double current_density_real, double current_density_imag, double permeability, double conductivity,
                        double remanence, double remanence_angle, double velocity_x, double velocity_y, double velocity_angular)
{
    for (unsigned int i = 0; i<magneticLabel.size(); i++)
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
    for (unsigned int i = 0; i<heatEdge.size(); i++)
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
    for (unsigned int i = 0; i<heatLabel.size(); i++)
    {
        heatLabel[i]->thermal_conductivity = thermal_conductivity;
        heatLabel[i]->volume_heat = volume_heat;
        heatLabel[i]->density = density;
        heatLabel[i]->specific_heat = specific_heat;
    }
}

void set_elasticity_edge(Hermes::vector<ElasticityEdge *> elasticityEdge, PhysicFieldBC typeX, PhysicFieldBC typeY, double forceX, double forceY)
{
    for (unsigned int i = 0; i<elasticityEdge.size(); i++)
    {
        elasticityEdge[i]->typeX = typeX;
        elasticityEdge[i]->typeY = typeY;
        elasticityEdge[i]->forceX = forceX;
        elasticityEdge[i]->forceY = forceY;
    }
}

void set_elasticity_label(Hermes::vector<ElasticityLabel *> elasticityLabel, double young_modulus, double poisson_ratio, double forceX, double forceY, double thermal_expansion)
{
    for (unsigned int i = 0; i<elasticityLabel.size(); i++)
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

    // harmonic elmag. field
    magneticEdge = new MagneticEdge[NUM_EDGES];
    set_magnetic_edge(Hermes::vector<MagneticEdge *>(&magneticEdge[13], &magneticEdge[81], &magneticEdge[54], &magneticEdge[55], &magneticEdge[56]), PhysicFieldBC_Magnetic_VectorPotential, 0.0, 0.0); //antisymetrie (osa rotace) nevim jakej je to typ, mam ted zase rozsekanej hermes
    set_magnetic_edge(Hermes::vector<MagneticEdge *>(&magneticEdge[24], &magneticEdge[25], &magneticEdge[46],
                                            &magneticEdge[47], &magneticEdge[48], &magneticEdge[49]), PhysicFieldBC_Magnetic_VectorPotential, 0.0, 0.0); //nulovy potencial, melo by byt spravne

    //nelinearni elektricka vodivost pouze u vodiveho zeleza a mosazi
    magneticLabel = new MagneticLabel[NUM_LABELS];
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[1], &magneticLabel[2], &magneticLabel[3]), 0.0, 0.0, 750.0, 6.6e6, 0.0, 0.0, 0.0, 0.0, 0.0); //vodive zelezo - jadro,ukotveni,kostra spojky
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[4], &magneticLabel[5]), 0.0, 0.0, 750.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //el.nevodive zelezo - plast
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[9], &magneticLabel[10]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //azbest
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[8]), 0.0, 0.0, 1.0, 1.6e6, 0.0, 0.0, 0.0, 0.0, 0.0); //mosaz
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[13], &magneticLabel[14], &magneticLabel[15],
                                                       &magneticLabel[16], &magneticLabel[17]), 1.0e6, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //civky zapnute
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[18]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //civka spojky
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[0], &magneticLabel[6], &magneticLabel[7], &magneticLabel[11], &magneticLabel[12]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //vzduch

    magneticLabels = Hermes::vector<int>(1,2,3,4,5,9,10,8,13,14,15,16,17,18,0,6,7,11,12);
    zelezoLabels = Hermes::vector<int>(1,2,3,4,5);

    vodiveZelezoMosazLabels = Hermes::vector<int>(1,2,3,8);

    // heat
    heatEdge = new HeatEdge[NUM_EDGES];
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[13], &heatEdge[54], &heatEdge[55], &heatEdge[56]), PhysicFieldBC_None, 0.0, 0.0, 0.0, 0.0); //osa rotace ???
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[24], &heatEdge[25], &heatEdge[47]), PhysicFieldBC_Heat_Temperature, 20.0, 0.0, 0.0, 0.0); //spodni hrana ukotveni ma pevnou teplotu - kvuli srovnani s Bohousem
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[45], &heatEdge[46], &heatEdge[1], &heatEdge[2], &heatEdge[3],
                                             &heatEdge[59], &heatEdge[60], &heatEdge[11], &heatEdge[12], &heatEdge[39], &heatEdge[40], &heatEdge[41],
                                             &heatEdge[42], &heatEdge[43]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 20.0, 20.0); //hranice


    //nelinearni teplotni vodivost pro obe zeleza
    heatLabel = new HeatLabel[NUM_LABELS];
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[1], &heatLabel[2], &heatLabel[3]), 50.0, 1.0, 7850.0, 469); //Fe_vod
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[4], &heatLabel[5]), 50.0, 0.0, 7850.0, 469); //Fe_nevod
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[13], &heatLabel[14], &heatLabel[15], &heatLabel[16], &heatLabel[17], &heatLabel[18]), 306.1, 0.0, 8700.0, 385);//med
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[6], &heatLabel[7], &heatLabel[11], &heatLabel[12]), 3e-3, 0.0, 1.29, 1);//vzduch
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[8]), 116, 1.0, 8400.0, 380);//mosaz
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[9], &heatLabel[10]), 0.12, 0.0, 2100.0, 220); //azbest

    heatLabels = Hermes::vector<int>(1,2,3,4,5,13,14,15,16,17,18,6,7,11,12,8,9,10);

    // thermoelasticity
    elasticityEdge = new ElasticityEdge[NUM_EDGES];
    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[1], &elasticityEdge[2], &elasticityEdge[57], &elasticityEdge[58]), PhysicFieldBC_Elasticity_Free, PhysicFieldBC_Elasticity_Free, 0.0, 0.0);
    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[14], &elasticityEdge[32]), PhysicFieldBC_Elasticity_Fixed, PhysicFieldBC_Elasticity_Fixed, 0.0, 0.0); //fixace v obou smerech
    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[13]), PhysicFieldBC_Elasticity_Fixed, PhysicFieldBC_Elasticity_Free, 0.0, 0.0);

    elasticityLabel = new ElasticityLabel[NUM_LABELS];
    set_elasticity_label(Hermes::vector<ElasticityLabel *>(&elasticityLabel[8]), 1e11, 0.25, 0.0, 0.0, 2.5e-5);

   elasticityLabels = Hermes::vector<int>(8);
}

