#include <hermes2d.h>
#include "data_table.h"

const int NUM_EDGES = 162;
const int NUM_LABELS = 68;

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

double permeability(double B, double T)
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
    // thermal conductivity STEEL CSN 12040
    double temp_thermal_conductivity[] = { 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
    double data_thermal_conductivity[] = { 47.9, 49.3, 48.6, 46.6, 44.1, 41.4, 38.0, 35.0, 32.4, 31.0, 30.0 };
    thermal_conductivity_fe.add(temp_thermal_conductivity, data_thermal_conductivity, 11);

    // electric conductivity STEEL CSN 12040
    double temp_electric_conductivity[] = { 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
    double data_electric_conductivity[] = { 5.70e6, 4.45e6, 3.36e6, 2.60e6, 2.05e6, 1.62e6, 1.32e6, 1.12e6, 1.01e6, 9.10e5, 8.90e5 };
    electric_conductivity_fe.add(temp_electric_conductivity, data_electric_conductivity, 8);

    // relative permeabilty temp
    double temp_relative_mag_permeability_temp[] = { 0, 100, 200, 300, 400, 500, 600, 700, 800, 10000 };
    double data_relative_mag_permeability_temp[] = { 2600, 2600, 2598, 2595, 2500, 2100, 1490, 750, 2, 1 };
    relative_mag_permeability_temp.add(temp_relative_mag_permeability_temp, data_relative_mag_permeability_temp, 10);


    // mi-B charackteristic
    double temp_relative_mag_permeability[] = { 1.2566370614359172e-07, 0.2402, 0.8654, 1.1106, 1.2458, 1.331, 1.5, 1.6, 1.683, 1.741, 1.78, 1.905, 2.025, 2.085, 2.13, 2.165, 2.28, 2.485, 2.585 };
    double data_relative_mag_permeability[] = { 1000.0, 1202.1703563104797, 2165.6082979831167, 1852.803771466027, 1558.7675165399623, 1332.2970393415892, 750.25900263307039, 400.01242373080828, 280.53809093387082, 217.63175928604289, 178.01671402763208, 95.25295840089872, 50.624981898320513, 34.750351479349241, 26.624990479857455, 21.650128290457886, 11.399995923769952, 6.2124977786334066, 5.1700046482110427 };
    relative_mag_permeability.add(temp_relative_mag_permeability, data_relative_mag_permeability, 10);

    // harmonic elmag. field
    magneticEdge = new MagneticEdge[NUM_EDGES];
    set_magnetic_edge(Hermes::vector<MagneticEdge *>(&magneticEdge[150], &magneticEdge[155], &magneticEdge[144]), PhysicFieldBC_Magnetic_VectorPotential, 0.0, 0.0); //antisymetrie (osa rotace) nevim jakej je to typ, mam ted zase rozsekanej hermes
    set_magnetic_edge(Hermes::vector<MagneticEdge *>(&magneticEdge[143]), PhysicFieldBC_Magnetic_VectorPotential, 0.0, 0.0); //nulovy potencial, melo by byt spravne

    magneticLabel = new MagneticLabel[NUM_LABELS];
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[5], &magneticLabel[67]), 0.0, 0.0, 2000.0, 5e6, 0.0, 0.0, 0.0, 0.0, 0.0); //vodive zelezo - kotouc
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[0]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //teflon
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[3], &magneticLabel[1], &magneticLabel[2], &magneticLabel[4]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //vata
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[7], &magneticLabel[8], &magneticLabel[9], &magneticLabel[10], &magneticLabel[11], &magneticLabel[12],
                                                       &magneticLabel[13], &magneticLabel[14], &magneticLabel[15], &magneticLabel[16], &magneticLabel[17], &magneticLabel[18],
                                                       &magneticLabel[19]), 4.2441e6, 0.0, 1.0, 5.7e7, 0.0, 0.0, 0.0, 0.0, 0.0); //vodice
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[20], &magneticLabel[21], &magneticLabel[22], &magneticLabel[23], &magneticLabel[24],
                                                       &magneticLabel[25], &magneticLabel[26], &magneticLabel[27], &magneticLabel[28], &magneticLabel[29], &magneticLabel[30],
                                                       &magneticLabel[31], &magneticLabel[36], &magneticLabel[35], &magneticLabel[34], &magneticLabel[33],
                                                       &magneticLabel[32]), 4.2441e6, 0.0, 1.0, 5.7e7, 0.0, 0.0, 0.0, 0.0, 0.0); //vodice pokracovani
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[6]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //vzduch
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[37], &magneticLabel[38], &magneticLabel[39], &magneticLabel[40], &magneticLabel[41], &magneticLabel[42], &magneticLabel[43],
                                                       &magneticLabel[44], &magneticLabel[45], &magneticLabel[46], &magneticLabel[47], &magneticLabel[48], &magneticLabel[49], &magneticLabel[50],
                                                       &magneticLabel[51]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //chladici voda
    set_magnetic_label(Hermes::vector<MagneticLabel *>(&magneticLabel[66], &magneticLabel[65], &magneticLabel[64], &magneticLabel[63], &magneticLabel[62], &magneticLabel[61],
                                                       &magneticLabel[60], &magneticLabel[59], &magneticLabel[58], &magneticLabel[57], &magneticLabel[56], &magneticLabel[55], &magneticLabel[54],
                                                       &magneticLabel[53], &magneticLabel[52]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //chladici voda - pokracovani

    magneticLabels = Hermes::vector<int>();
    for (unsigned int i = 0; i<68; i++)
    {
        magneticLabels.push_back(i);
    }

    // heat
    heatEdge = new HeatEdge[NUM_EDGES];
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[155]), PhysicFieldBC_None, 0.0, 0.0, 0.0, 0.0); //osa rotace ???

    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[151], &heatEdge[132]), PhysicFieldBC_None, 0.0, 0.0, 0.0, 0.0); //spodni a horni hrana teflonove trubky - dlouha trubka

    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[129], &heatEdge[137], &heatEdge[136], &heatEdge[134], &heatEdge[135],
                                             &heatEdge[138], &heatEdge[139], &heatEdge[142], &heatEdge[14], &heatEdge[8], &heatEdge[126], &heatEdge[128],
                                             &heatEdge[127]), PhysicFieldBC_Heat_Flux, 0.0, 0.0, 20.0, 30.0); //hranice

    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[101], &heatEdge[102], &heatEdge[77], &heatEdge[78], &heatEdge[91], &heatEdge[92], &heatEdge[93], &heatEdge[94], &heatEdge[95], &heatEdge[96],
                                             &heatEdge[86]), PhysicFieldBC_Heat_Temperature, 50.0, 0.0, 0.0, 0.0); //vnitrky vodicu - pevna teplot adana chladici vodou
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[85], &heatEdge[97], &heatEdge[98], &heatEdge[105], &heatEdge[116], &heatEdge[113], &heatEdge[112], &heatEdge[108], &heatEdge[106],
                                             &heatEdge[117]), PhysicFieldBC_Heat_Temperature, 50.0, 0.0, 0.0, 0.0); //vnitrky vodicu - pevna teplot adana chladici vodou
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[118], &heatEdge[121], &heatEdge[122], &heatEdge[125], &heatEdge[71], &heatEdge[35], &heatEdge[32], &heatEdge[148], &heatEdge[147],
                                             &heatEdge[67]), PhysicFieldBC_Heat_Temperature, 50.0, 0.0, 0.0, 0.0); //vnitrky vodicu - pevna teplot adana chladici vodou
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[65], &heatEdge[64], &heatEdge[63], &heatEdge[60], &heatEdge[59], &heatEdge[56], &heatEdge[55], &heatEdge[52], &heatEdge[51],
                                             &heatEdge[22]), PhysicFieldBC_Heat_Temperature, 50.0, 0.0, 0.0, 0.0); //vnitrky vodicu - pevna teplot adana chladici vodou
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[24], &heatEdge[30], &heatEdge[33], &heatEdge[42], &heatEdge[21], &heatEdge[69], &heatEdge[68], &heatEdge[31], &heatEdge[34],
                                             &heatEdge[50]), PhysicFieldBC_Heat_Temperature, 50.0, 0.0, 0.0, 0.0); //vnitrky vodicu - pevna teplot adana chladici vodou
    set_heat_edge(Hermes::vector<HeatEdge *>(&heatEdge[49], &heatEdge[48], &heatEdge[47], &heatEdge[46], &heatEdge[45], &heatEdge[12], &heatEdge[11], &heatEdge[74],
                                             &heatEdge[73]), PhysicFieldBC_Heat_Temperature, 50.0, 0.0, 0.0, 0.0); //vnitrky vodicu - pevna teplot adana chladici vodou


    heatLabel = new HeatLabel[NUM_LABELS];
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[5], &heatLabel[67]), 50.0, 1.0, 7620.0, 550); //vodive zelezo - kotouc
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[0]), 0.24, 0.0, 2220.0, 1050); //teflon
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[3], &heatLabel[1], &heatLabel[2], &heatLabel[4]), 0.04, 0.0, 72.0, 670); //vata
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[7], &heatLabel[8], &heatLabel[9], &heatLabel[10], &heatLabel[11], &heatLabel[12],
                                               &heatLabel[13], &heatLabel[14], &heatLabel[15], &heatLabel[16], &heatLabel[17], &heatLabel[18]), 395.0, 1.0, 8960.0, 390); //vodice
    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[19], &heatLabel[20], &heatLabel[21], &heatLabel[22], &heatLabel[23], &heatLabel[24],
                                               &heatLabel[25], &heatLabel[26], &heatLabel[27], &heatLabel[28], &heatLabel[29], &heatLabel[30],
                                               &heatLabel[31], &heatLabel[36], &heatLabel[35], &heatLabel[34], &heatLabel[33],
                                               &heatLabel[32]), 395.0, 1.0, 8960.0, 390); //vodice
//    set_heat_label(Hermes::vector<HeatLabel *>(&heatLabel[37], &heatLabel[38], &heatLabel[39], &heatLabel[40], &heatLabel[41], &heatLabel[42], &heatLabel[43],
//                                                       &heatLabel[44], &heatLabel[45], &heatLabel[46], &heatLabel[47], &heatLabel[48], &heatLabel[49], &heatLabel[50],
//                                                       &heatLabel[51], &heatLabel[66], &heatLabel[65], &heatLabel[64], &heatLabel[63], &heatLabel[62], &heatLabel[61],
//                                                       &heatLabel[60], &heatLabel[59], &heatLabel[58], &heatLabel[57], &heatLabel[56], &heatLabel[55], &heatLabel[54],
//                                                       &heatLabel[53], &heatLabel[52]), 50.0, 1.0, 7850.0, 469); //chladici voda - neni v teplotnim modelu

    heatLabels = Hermes::vector<int>(0,1,2,3,4,5);
    for (unsigned int i=7; i<37; i++)
    {
        heatLabels.push_back(i);
    }
    heatLabels.push_back(67);


    // thermoelasticity
    elasticityEdge = new ElasticityEdge[NUM_EDGES];
    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[152]), PhysicFieldBC_Elasticity_Free, PhysicFieldBC_Elasticity_Fixed, 0.0, 0.0);//stred kotouce - fixovan v Z
    set_elasticity_edge(Hermes::vector<ElasticityEdge *>(&elasticityEdge[4], &elasticityEdge[154], &elasticityEdge[158], &elasticityEdge[153], &elasticityEdge[159], &elasticityEdge[5],
                                                         &elasticityEdge[6], &elasticityEdge[3], &elasticityEdge[2], &elasticityEdge[161], &elasticityEdge[157], &elasticityEdge[160],
                                                         &elasticityEdge[156], &elasticityEdge[1]), PhysicFieldBC_Elasticity_Free, PhysicFieldBC_Elasticity_Free, 0.0, 0.0);//free

    elasticityLabel = new ElasticityLabel[NUM_LABELS];
    set_elasticity_label(Hermes::vector<ElasticityLabel *>(&elasticityLabel[5],&elasticityLabel[67]), 1e11, 0.3, 0.0, 0.0, 1.25e-5);

   elasticityLabels = Hermes::vector<int>(5,67);
}

