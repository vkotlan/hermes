#include "data_table.h"

#define EPS0 8.854e-12
#define MU0 4*M_PI*1e-7

const double frequency = 5000.0;

const double FINAL_TIME = 70.0; // Length of time interval in seconds.
double timeStep = 0.5;

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

    double perm = relative_mag_permeability.value(T);

    perm = perm * koef;
    if (perm < 1.0) perm = 1.0;

    return perm;
}

template<typename Real, typename Scalar>
Scalar int_x_u_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->val[i] * v->val[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_v(int n, double *wt, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (v->val[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_u_dvdx_over_x(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (v->dx[i] * u->val[i]) / e->x[i];
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_grad_u_grad_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_dudx_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_dudy_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->dy[i] * v->dy[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_dudx_dvdy(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->dx[i] * v->dy[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_dudy_dvdx(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (v->dx[i] * u->dy[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_dudx_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->dx[i] * v->val[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar int_x_dudy_v(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * e->x[i] * (u->dy[i] * v->val[i]);
    return result;
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

BCType magnetic_bc_types(int marker)
{
    switch (magneticEdge[marker].type)
    {
    case PhysicFieldBC_None:
        return BC_NONE;
    case PhysicFieldBC_Magnetic_VectorPotential:
        return BC_ESSENTIAL;
    case PhysicFieldBC_Magnetic_SurfaceCurrent:
        return BC_NATURAL;
    }
}

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

template<typename Real, typename Scalar>
Scalar magnetic_vector_form_linear_surf_real(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double K = 0.0;

    if (magneticEdge[e->edge_marker].type == PhysicFieldBC_Magnetic_SurfaceCurrent)
        K = magneticEdge[e->edge_marker].value_real;

    return K * 2 * M_PI * int_x_v<Real, Scalar>(n, wt, v, e);
}

template<typename Real, typename Scalar>
Scalar magnetic_vector_form_linear_surf_imag(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double K = 0.0;

    if (magneticEdge[e->edge_marker].type == PhysicFieldBC_Magnetic_SurfaceCurrent)
        K = magneticEdge[e->edge_marker].value_imag;


    return K * 2 * M_PI * int_x_v<Real, Scalar>(n, wt, v, e);
}

template<typename Real, typename Scalar>
Scalar magnetic_matrix_form_linear_real_real(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        double perm = magneticLabel[e->elem_marker].permeability * MU0;
        if (e->elem_marker == 0)
            // perm = permeability(sqrt(sqr(u->dx[i]) + sqr(u->dy[i])), ext->fn[0]->val[i]) * MU0;
            perm = permeability(sqrt(sqr(u->dx[i]) + sqr(u->dy[i])), 20.0) * MU0;

        result +=  wt[i] * 1.0 / perm * (((v->dx[i] * u->val[i]) / e->x[i]) + (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
    }
    return result;
}

template<typename Real, typename Scalar>
Scalar magnetic_matrix_form_linear_real_imag(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return - 2 * M_PI * frequency * magneticLabel[e->elem_marker].conductivity * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar magnetic_matrix_form_linear_imag_real(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return + 2 * M_PI * frequency * magneticLabel[e->elem_marker].conductivity * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar magnetic_matrix_form_linear_imag_imag(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        double perm = magneticLabel[e->elem_marker].permeability * MU0;
        if (e->elem_marker == 0)
            // perm = permeability(sqrt(sqr(u->dx[i]) + sqr(u->dy[i])), ext->fn[0]->val[i]) * MU0;
            perm = permeability(sqrt(sqr(u->dx[i]) + sqr(u->dy[i])), 20.0) * MU0;

        result +=  wt[i] * 1.0 / perm * (((v->dx[i] * u->val[i]) / e->x[i]) + (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
    }
    return result;
}

template<typename Real, typename Scalar>
Scalar magnetic_vector_form_linear_real(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return magneticLabel[e->elem_marker].current_density_real * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar magnetic_vector_form_linear_imag(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return magneticLabel[e->elem_marker].current_density_imag * int_v<Real, Scalar>(n, wt, v);
}

Ord magnetic_matrix_form_linear_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
    return u->val[0] * v->val[0] * e->x[0] * e->x[0];
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

BCType heat_bc_types(int marker)
{
    // info("heat_bc_types");

    switch (heatEdge[marker].type)
    {
    case PhysicFieldBC_None:
        return BC_NONE;
    case PhysicFieldBC_Heat_Temperature:
        return BC_ESSENTIAL;
    case PhysicFieldBC_Heat_Flux:
        return BC_NATURAL;
    }
}

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
}

// linear forms
template<typename Real, typename Scalar>
Scalar heat_matrix_form_linear_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    double h = 0.0;

    if (heatEdge[e->edge_marker].type == PhysicFieldBC_Heat_Flux)
        h = heatEdge[e->edge_marker].h;

    // info("heat_matrix_form_linear_surf, h = %3.5f", h);

    return h * 2 * M_PI * int_x_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar heat_vector_form_linear_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    if (heatEdge[e->edge_marker].type == PhysicFieldBC_None)
        return 0.0;

    double q = 0.0;
    double h = 0.0;
    double Text = 0.0;

    if (heatEdge[e->edge_marker].type == PhysicFieldBC_Heat_Flux)
    {
        q = heatEdge[e->edge_marker].heatFlux;
        h = heatEdge[e->edge_marker].h;
        Text = heatEdge[e->edge_marker].externalTemperature;
    }

    // info("heat_vector_form_linear_surf, h = %3.5f", h);

    return (q + Text * h) * 2 * M_PI * int_x_v<Real, Scalar>(n, wt, v, e);
}

// Integration order for the thermoelasticity_matrix form
Ord heat_matrix_form_linear_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
    return u->val[0] * v->val[0] * e->x[0] * e->x[0];
}

// Integration order for the volumetric linear form
Ord heat_vector_form_linear_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
    return v->val[0] * e->x[0] * e->x[0];
}

template<typename Real, typename Scalar>
double heat_matrix_form_linear(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    return heatLabel[e->elem_marker].thermal_conductivity * 2 * M_PI * int_x_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e)
            + heatLabel[e->elem_marker].density * heatLabel[e->elem_marker].specific_heat * 2 * M_PI * int_x_u_v<Real, Scalar>(n, wt, u, v, e) / timeStep;
}

template<typename Real, typename Scalar>
double heat_vector_form_linear(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
    {
        result += ext->fn[0]->val[i] * heatLabel[e->elem_marker].volume_heat * 2 * M_PI * wt[i] * e->x[i] * (v->val[i])
                + heatLabel[e->elem_marker].density * heatLabel[e->elem_marker].specific_heat * 2 * M_PI * wt[i] * e->x[i] * (ext->fn[1]->val[i] * v->val[i]) / timeStep;
    }
    return result;
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

BCType elasticity_bc_types_r(int marker)
{
    switch (elasticityEdge[marker].typeX)
    {
    case PhysicFieldBC_None:
        return BC_NONE;
        break;
    case PhysicFieldBC_Elasticity_Fixed:
        return BC_ESSENTIAL;
        break;
    case PhysicFieldBC_Elasticity_Free:
        return BC_NATURAL;
        break;
    }
}

BCType elasticity_bc_types_z(int marker)
{
    switch (elasticityEdge[marker].typeY)
    {
    case PhysicFieldBC_None:
        return BC_NONE;
        break;
    case PhysicFieldBC_Elasticity_Fixed:
        return BC_ESSENTIAL;
        break;
    case PhysicFieldBC_Elasticity_Free:
        return BC_NATURAL;
        break;
    }
}

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
}

template<typename Real, typename Scalar>
Scalar elasticity_matrix_form_r_r(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (e->x[i] * elasticityLabel[e->elem_marker].lambda() * (u->dx[i] * v->dx[i] +
                                                                            u->val[i]/e->x[i] * v->dx[i] +
                                                                            u->dx[i] * v->val[i]/e->x[i] +
                                                                            1/sqr(e->x[i]) * u->val[i] * v->val[i]) +
                           e->x[i] * elasticityLabel[e->elem_marker].mu() * (2 * u->dx[i] * v->dx[i] +
                                                                        2 * 1/sqr(e->x[i]) * u->val[i] * v->val[i] +
                                                                        u->dy[i] * v->dy[i]));
    return result;
}

template<typename Real, typename Scalar>
Scalar elasticity_matrix_form_r_z(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (e->x[i] * elasticityLabel[e->elem_marker].lambda() * (u->dy[i] * v->dx[i] +
                                                                            u->dy[i] * v->val[i]/e->x[i]) +
                           e->x[i] * elasticityLabel[e->elem_marker].mu() * u->dx[i] * v->dy[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar elasticity_matrix_form_z_z(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += wt[i] * (e->x[i] * elasticityLabel[e->elem_marker].lambda() * (u->dy[i] * v->dy[i]) +
                           e->x[i] * elasticityLabel[e->elem_marker].mu() * (u->dx[i] * v->dx[i] +
                                                                        2 * u->dy[i] * v->dy[i]));
    return result;
}

template<typename Real, typename Scalar>
Scalar elasticity_vector_form_r(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += - wt[i] * (3 * elasticityLabel[e->elem_marker].lambda() + 2 * elasticityLabel[e->elem_marker].mu()) * elasticityLabel[e->elem_marker].thermal_expansion *
                (e->x[i] * (ext->fn[0]->val[i] - 20.0) * v->dx[i]); // e->x[i] *
    return result;
}

template<typename Real, typename Scalar>
Scalar elasticity_vector_form_z(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += - wt[i] * (3 * elasticityLabel[e->elem_marker].lambda() + 2 * elasticityLabel[e->elem_marker].mu()) * elasticityLabel[e->elem_marker].thermal_expansion *
                (e->x[i] * (ext->fn[0]->val[i] - 20.0) * v->dy[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar elasticity_matrix_form_r_T(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += - wt[i] * (3*elasticityLabel[e->elem_marker].lambda() + 2*elasticityLabel[e->elem_marker].mu()) * elasticityLabel[e->elem_marker].thermal_expansion *
                  (e->x[i] * u->dx[i] * v->val[i]);
    return result;
}

template<typename Real, typename Scalar>
Scalar elasticity_matrix_form_z_T(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
    Scalar result = 0;
    for (int i = 0; i < n; i++)
        result += - wt[i] * (3*elasticityLabel[e->elem_marker].lambda() + 2*elasticityLabel[e->elem_marker].mu()) * elasticityLabel[e->elem_marker].thermal_expansion *
                (e->x[i] * u->dy[i] * v->val[i]);
    return result;
}
