#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"
#include "weakform_library/maxwell.h"
#include "function/filter.h"

// zero
#define EPS_ZERO 1e-10

// physical constants
#define EPS0 8.854e-12
#define MU0 4*M_PI*1e-7

const double frequency = 500;

using namespace WeakFormsH1;
using namespace WeakFormsH1::SurfaceMatrixForms;
using namespace WeakFormsH1::SurfaceVectorForms;
using namespace WeakFormsH1::VolumetricMatrixForms;
using namespace WeakFormsH1::VolumetricVectorForms;

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

        if(nodes->present(order)) {
          assert(nodes->get(order) == cur_node);
          ::free(nodes->get(order));
        }
        nodes->add(node, order);
        cur_node = node;
    }

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); }
};

class MagneticVectorPotentialFilter : public Filter
{
public:
    MagneticVectorPotentialFilter(Hermes::vector<MeshFunction*> solutions)
        : Filter(solutions)
    {

    }

    inline void precalculate(int order, int mask)
    {
        Quad2D* quad = quads[cur_quad];
        int np = quad->get_num_points(order);
        Node* node = new_node(H2D_FN_VAL_0, np);

        sln[0]->set_quad_order(order, H2D_FN_VAL);

        scalar *u = sln[0]->get_fn_values();
        update_refmap();
        double *x = refmap->get_phys_x(order);

        for (int i = 0; i < np; i++)
        {
            node->values[0][0][i] = - x[i] * u[i];
        }

        if(nodes->present(order)) {
          assert(nodes->get(order) == cur_node);
          ::free(nodes->get(order));
        }
        nodes->add(node, order);
        cur_node = node;
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

        std::string user_marker = mesh->get_element_markers_conversion().get_user_marker(e->marker);
        int marker = atoi(user_marker.c_str());
        for (int i = 0; i < np; i++)
        {
            if(magneticLabel[marker].conductivity > 0.0){
//                printf("jednickovy marker %d, cond %lf\n", marker, magneticLabel[marker].conductivity);
                node->values[0][0][i] = 0.5 / magneticLabel[marker].conductivity * (
                    sqr(2 * M_PI * frequency * magneticLabel[marker].conductivity * value2[i]) +
                    sqr(2 * M_PI * frequency * magneticLabel[marker].conductivity * value1[i]));
//                node->values[0][0][i] = 1.0;
            }
            else
                node->values[0][0][i] = 0.0;

        }

        if(nodes->present(order)) {
          assert(nodes->get(order) == cur_node);
          ::free(nodes->get(order));
        }
        nodes->add(node, order);
        cur_node = node;
    }

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); }
};

class CustomVectorFormTimeDep : public WeakForm::VectorFormVol
{
public:
    CustomVectorFormTimeDep(int i, scalar coeff, Solution *solution, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i), coeff(coeff), gt(gt)
    {
        ext.push_back(solution);
    }
    CustomVectorFormTimeDep(int i, std::string area, scalar coeff, Solution *solution, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt)
    {
        ext.push_back(solution);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const {
        if (gt == HERMES_PLANAR)
            return coeff * int_u_v<double, scalar>(n, wt, ext->fn[0], v);
        else if (gt == HERMES_AXISYM_X)
            return coeff * int_y_u_v<double, scalar>(n, wt, ext->fn[0], v, e);
        else
            return coeff * int_x_u_v<double, scalar>(n, wt, ext->fn[0], v, e);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
        if (gt == HERMES_PLANAR)
            return int_u_v<Ord, Ord>(n, wt, ext->fn[0], v);
        else if (gt == HERMES_AXISYM_X)
            return int_y_u_v<Ord, Ord>(n, wt, ext->fn[0], v, e);
        else
            return int_x_u_v<Ord, Ord>(n, wt, ext->fn[0], v, e);
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::VectorFormVol* clone() {
        return new CustomVectorFormTimeDep(*this);
    }

private:
    scalar coeff;
    GeomType gt;
};

class DefaultLinearThermoelasticityX : public WeakForm::VectorFormVol
{
public:
    DefaultLinearThermoelasticityX(int i, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

    DefaultLinearThermoelasticityX(int i, int j, std::string area, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, area), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0.0;
        if (gt == HERMES_PLANAR)
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dx[i];
        else if (gt == HERMES_AXISYM_X)
            for (int i = 0; i < n; i++)
                result += wt[i] * e->x[i] * v->dy[i];
        else
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dx[i];

        return (3*lambda + 2*mu) * alpha * (temp - temp_ref) * result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR)
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dx[i];
        else if (gt == HERMES_AXISYM_X)
            for (int i = 0; i < n; i++)
                result += wt[i] * e->x[i] * v->dy[i];
        else
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dx[i];

        return result;
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::VectorFormVol* clone() {
        return new DefaultLinearThermoelasticityX(*this);
    }

private:
    scalar lambda, mu, alpha, temp, temp_ref, vel_ang;
    GeomType gt;
};

class DefaultLinearThermoelasticityY : public WeakForm::VectorFormVol
{
public:
    DefaultLinearThermoelasticityY(int i, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

    DefaultLinearThermoelasticityY(int i, std::string area, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, area), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const {
        scalar result = 0.0;
        if (gt == HERMES_PLANAR)
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dy[i];
        else if (gt == HERMES_AXISYM_X)
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dx[i];
        else
            for (int i = 0; i < n; i++)
                result += wt[i] * e->x[i] * v->dy[i];

        return (3*lambda + 2*mu) * alpha * (temp - temp_ref) * result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
        Ord result = 0;
        if (gt == HERMES_PLANAR)
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dy[i];
        else if (gt == HERMES_AXISYM_X)
            for (int i = 0; i < n; i++)
                result += wt[i] * v->dx[i];
        else
            for (int i = 0; i < n; i++)
                result += wt[i] * e->x[i] * v->dy[i];

        return result;
    }

    // This is to make the form usable in rk_time_step().
   virtual WeakForm::VectorFormVol* clone() {
        return new DefaultLinearThermoelasticityY(*this);
    }

private:
    scalar lambda, mu, alpha, temp, temp_ref, vel_ang;
    GeomType gt;
};

class WeakFormMagnetic : public WeakForm
{
public:
    WeakFormMagnetic(int neq) : WeakForm(neq) { }

    void registerForms()
    {
        // magnetic field
//        add_magnetic_material("0", 700.0, 1e6, 0.0);
//        add_magnetic_material("1", 1.0, 0.0, 6e7);
//        add_magnetic_material("2", 1.0, 0.0, 6e7);
//        add_magnetic_material("3", 1.0, 0.0, 0.0);
//        add_magnetic_material("4", 1.0, 0.0, 0.0);
//        add_magnetic_material("5", 1.0, 0.0, 0.0);
//        add_magnetic_material("6", 1.0, 0.0, 0.0);

        for(int i=0; i<NUM_LABELS; i++){
            char str[5];
            sprintf(str, "%d", i);
            //TODO only real part of current density
            add_magnetic_material(str, magneticLabel[i].permeability, magneticLabel[i].conductivity, magneticLabel[i].current_density_real);
        }
    }

    void add_magnetic_material(std::string marker, double permeability, double conductivity, double external_current_density)
    {
        // real part
        add_matrix_form(new WeakFormsMaxwell::VolumetricMatrixForms::DefaultLinearMagnetostatics(0, 0, marker, 1.0 / (permeability * MU0), HERMES_NONSYM, HERMES_AXISYM_Y));
        // imag part
        add_matrix_form(new WeakFormsMaxwell::VolumetricMatrixForms::DefaultLinearMagnetostatics(1, 1, marker, 1.0 / (permeability * MU0), HERMES_NONSYM, HERMES_AXISYM_Y));
        // conductivity
        if (fabs(conductivity) > EPS_ZERO)
        {
            add_matrix_form(new WeakFormsH1::VolumetricMatrixForms::DefaultLinearMass(0, 1, marker, - 2 * M_PI * frequency * conductivity, HERMES_NONSYM, HERMES_PLANAR));
            add_matrix_form(new WeakFormsH1::VolumetricMatrixForms::DefaultLinearMass(1, 0, marker,   2 * M_PI * frequency * conductivity, HERMES_NONSYM, HERMES_PLANAR));
        }
        // external current density
        if (fabs(external_current_density) > EPS_ZERO)
            add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormConst(0, marker, external_current_density, HERMES_PLANAR));
    }
};

const double TIME_CONSTANT = 10000.;

class WeakFormTemp : public WeakForm
{
public:
    WeakFormTemp(double time_step) : WeakForm(1), time_step(time_step) { };
    void registerForms(Solution *prev_time_sln)
    {
        // temperature field
        add_temperature_material("0", 1.0, 50.0, 1000., 1.0, prev_time_sln);
        add_temperature_material("6", 1.0, 50.0, 1000., 1.0, prev_time_sln);
    }
    void add_temperature_material(std::string marker, double thermal_conductivity, double volume_heat, double density, double specific_heat, Solution* prev_time_sln)
    {
        /// TODO vsude jsem smazal 2 * M_PI

        // Contribution of the time derivative term.
        /// TODO NONSYM
        printf("TIME STEP %g\n", time_step);
        add_matrix_form(new DefaultLinearMass(0, 0, marker, density * specific_heat / time_step, HERMES_NONSYM, HERMES_AXISYM_Y));

        // Contribution of the diffusion term.
        /// TODO NONSYM
        add_matrix_form(new DefaultLinearDiffusion(0, 0, marker, thermal_conductivity, HERMES_NONSYM, HERMES_AXISYM_Y));

        //Heat sources
        add_vector_form(new DefaultVectorFormConst(0, marker, volume_heat, HERMES_AXISYM_Y));

        // Right-hand side volumetric vector form.
        CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, marker, 2 * M_PI * density * specific_heat / time_step);
        vec_form_vol->ext.push_back(prev_time_sln);
        add_vector_form(vec_form_vol);

        //add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, bdy_air, alpha / (density * heatcap)));

        //add_vector_form_surf(new CustomVectorFormSurf(0, bdy_air, alpha, density, heatcap,
        //                     current_time_ptr, temp_init, t_final));

    }

private:
    double time_step;

    class CustomVectorFormSource : public WeakForm::VectorFormVol
    {
    public:
        CustomVectorFormSource(int i, std::string marker)
            : WeakForm::VectorFormVol(i, marker) { }

        template<typename Real, typename Scalar>
        Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
            Func<Real>* source = ext->fn[0];

        }
    };

    // This form is custom since it contains previous time-level solution.
    class CustomVectorFormVol : public WeakForm::VectorFormVol
    {
    public:
        CustomVectorFormVol(int i, std::string marker, double time_contrib_coef)
        : WeakForm::VectorFormVol(i, marker), time_contrib_coef(time_contrib_coef) { }

      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        Func<Real>* temp_prev_time = ext->fn[0];
       // printf("func: %p\n", temp_prev_time);

//        std::string user_marker = this-> wf->get_element_markers_conversion()-> .get_user_marker(e->marker);
//        int marker = atoi(user_marker.c_str());


        //double density = heatLabel[atoi(marker.c_str())].density;
        double density = 1000.;
        //printf("  density %lf\n", density);
       // if(atoi(marker.c_str()) == 0)
       //     density = 1000000000000000.;
        double specific_heat = 1.;
        return int_x_u_v<Real, Scalar>(n, wt, temp_prev_time, v, e)  * time_contrib_coef;
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        double contr = vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
        if (fabs(contr) > 1e-7)
            printf("contribution prev %g\n", contr);
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      double time_contrib_coef;
    };

    // This form is custom since it contains time-dependent exterior temperature.
    class CustomVectorFormSurf : public WeakForm::VectorFormSurf
    {
    public:
      CustomVectorFormSurf(int i, std::string area, double alpha, double density, double heatcap,
                                  double* current_time_ptr, double temp_init, double t_final)
        : WeakForm::VectorFormSurf(i, area), alpha(alpha), density(density), heatcap(heatcap), current_time_ptr(current_time_ptr),
                                   temp_init(temp_init), t_final(t_final) { }

     template<typename Real, typename Scalar>
      Scalar vector_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        return  0;//alpha / (density * heatcap) * temp_ext(*current_time_ptr + time_step) * int_v<Real>(n, wt, v);
      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
          return vector_form_surf<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
          return vector_form_surf<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }
      // Time-dependent exterior temperature.
      template<typename Real>
      Real temp_ext(Real t) const {
        return temp_init + 10. * sin(2*M_PI*t/t_final);
      }

      double alpha, density, heatcap, *current_time_ptr, temp_init, t_final;
    };


};

/*

class WeakFormElasticity : public WeakFormAgros
{
public:
    WeakFormElasticity() : WeakFormAgros(2) { }

    void registerForms()
    {
        // boundary conditions
        for (int i = 0; i<Util::scene()->edges.count(); i++)
        {
            SceneBoundaryElasticity *boundary = dynamic_cast<SceneBoundaryElasticity *>(Util::scene()->edges[i]->boundary);

            if (boundary && Util::scene()->edges[i]->boundary != Util::scene()->boundaries[0])
            {
                if (boundary->typeX == PhysicFieldBC_Elasticity_Free)
                {
                    if (fabs(boundary->forceX.number) > EPS_ZERO)
                        add_vector_form_surf(new WeakFormsH1::SurfaceVectorForms::DefaultVectorFormSurf(0,
                                                                                                        QString::number(i + 1).toStdString(),
                                                                                                        boundary->forceX.number,
                                                                                                        convertProblemType(Util::scene()->problemInfo()->problemType)));

                }
                if (boundary->typeY == PhysicFieldBC_Elasticity_Free)
                {
                    if (fabs(boundary->forceY.number) > EPS_ZERO)
                        add_vector_form_surf(new WeakFormsH1::SurfaceVectorForms::DefaultVectorFormSurf(1,
                                                                                                        QString::number(i + 1).toStdString(),
                                                                                                        boundary->forceY.number,
                                                                                                        convertProblemType(Util::scene()->problemInfo()->problemType)));
                }
            }
        }

        // materials (Default forms not implemented axisymmetric problems!)
        for (int i = 0; i<Util::scene()->labels.count(); i++)
        {
            SceneMaterialElasticity *material = dynamic_cast<SceneMaterialElasticity *>(Util::scene()->labels[i]->material);

            if (material && Util::scene()->labels[i]->material != Util::scene()->materials[0])
            {
                add_matrix_form(new WeakFormsElasticity::VolumetricMatrixForms::DefaultLinearXX(0, 0,
                                                                                                QString::number(i).toStdString(),
                                                                                                material->lambda(), material->mu(),
                                                                                                convertProblemType(Util::scene()->problemInfo()->problemType)));

                add_matrix_form(new WeakFormsElasticity::VolumetricMatrixForms::DefaultLinearXY(0, 1,
                                                                                                QString::number(i).toStdString(),
                                                                                                material->lambda(), material->mu(),
                                                                                                convertProblemType(Util::scene()->problemInfo()->problemType)));

                add_matrix_form(new WeakFormsElasticity::VolumetricMatrixForms::DefaultLinearYY(1, 1,
                                                                                                QString::number(i).toStdString(),
                                                                                                material->lambda(), material->mu(),
                                                                                                convertProblemType(Util::scene()->problemInfo()->problemType)));

                // inner forces
                if (fabs(material->forceX.number) > EPS_ZERO)
                    add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormConst(0,
                                                                                                   QString::number(i).toStdString(),
                                                                                                   material->forceX.number,
                                                                                                   convertProblemType(Util::scene()->problemInfo()->problemType)));
                if (fabs(material->forceY.number) > EPS_ZERO)
                    add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormConst(1,
                                                                                                   QString::number(i).toStdString(),
                                                                                                   material->forceY.number,
                                                                                                   convertProblemType(Util::scene()->problemInfo()->problemType)));

                // thermoelasticity
                if ((fabs(material->alpha.number) > EPS_ZERO) &&
                        (fabs(material->temp.number - material->temp_ref.number) > EPS_ZERO))
                    add_vector_form(new DefaultLinearThermoelasticityX(0, 0,
                                                                      QString::number(i).toStdString(),
                                                                      material->lambda(), material->mu(),
                                                                      material->alpha.number, material->temp.number, material->temp_ref.number,
                                                                      convertProblemType(Util::scene()->problemInfo()->problemType)));
                if ((fabs(material->alpha.number) > EPS_ZERO) &&
                        (fabs(material->temp.number - material->temp_ref.number) > EPS_ZERO))
                    add_vector_form(new DefaultLinearThermoelasticityY(1,
                                                                      QString::number(i).toStdString(),
                                                                      material->lambda(), material->mu(),
                                                                      material->alpha.number, material->temp.number, material->temp_ref.number,
                                                                      convertProblemType(Util::scene()->problemInfo()->problemType)));
            }
        }
    }

    class DefaultLinearThermoelasticityX : public WeakForm::VectorFormVol
    {
    public:
        DefaultLinearThermoelasticityX(int i, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
            : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

        DefaultLinearThermoelasticityX(int i, int j, std::string area, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
            : WeakForm::VectorFormVol(i, area), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                             Geom<double> *e, ExtData<scalar> *ext) const {
            scalar result = 0.0;
            if (gt == HERMES_PLANAR)
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dx[i];
            else if (gt == HERMES_AXISYM_X)
                for (int i = 0; i < n; i++)
                    result += wt[i] * e->x[i] * v->dy[i];
            else
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dx[i];

            return (3*lambda + 2*mu) * alpha * (temp - temp_ref) * result;
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const {
            Ord result = 0;
            if (gt == HERMES_PLANAR)
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dx[i];
            else if (gt == HERMES_AXISYM_X)
                for (int i = 0; i < n; i++)
                    result += wt[i] * e->x[i] * v->dy[i];
            else
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dx[i];

            return result;
        }

        // This is to make the form usable in rk_time_step().
        virtual WeakForm::VectorFormVol* clone() {
            return new DefaultLinearThermoelasticityX(*this);
        }

    private:
        scalar lambda, mu, alpha, temp, temp_ref, vel_ang;
        GeomType gt;
    };

    class DefaultLinearThermoelasticityY : public WeakForm::VectorFormVol
    {
    public:
        DefaultLinearThermoelasticityY(int i, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
            : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

        DefaultLinearThermoelasticityY(int i, std::string area, scalar lambda, scalar mu, scalar alpha, scalar temp, scalar temp_ref, GeomType gt = HERMES_PLANAR)
            : WeakForm::VectorFormVol(i, area), lambda(lambda), mu(mu), alpha(alpha), temp(temp), temp_ref(temp_ref), gt(gt) { }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                             Geom<double> *e, ExtData<scalar> *ext) const {
            scalar result = 0.0;
            if (gt == HERMES_PLANAR)
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dy[i];
            else if (gt == HERMES_AXISYM_X)
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dx[i];
            else
                for (int i = 0; i < n; i++)
                    result += wt[i] * e->x[i] * v->dy[i];

            return (3*lambda + 2*mu) * alpha * (temp - temp_ref) * result;
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                        Geom<Ord> *e, ExtData<Ord> *ext) const {
            Ord result = 0;
            if (gt == HERMES_PLANAR)
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dy[i];
            else if (gt == HERMES_AXISYM_X)
                for (int i = 0; i < n; i++)
                    result += wt[i] * v->dx[i];
            else
                for (int i = 0; i < n; i++)
                    result += wt[i] * e->x[i] * v->dy[i];

            return result;
        }

        // This is to make the form usable in rk_time_step().
       virtual WeakForm::VectorFormVol* clone() {
            return new DefaultLinearThermoelasticityY(*this);
        }

    private:
        scalar lambda, mu, alpha, temp, temp_ref, vel_ang;
        GeomType gt;
    };
};

class WeakFormHeat : public WeakFormAgros
{
public:
    WeakFormHeat() : WeakFormAgros() { }

    void registerForms()
    {
        // delete all forms
        delete_all();

        // boundary conditions
        for (int i = 0; i<Util::scene()->edges.count(); i++)
        {
            SceneBoundaryHeat *boundary = dynamic_cast<SceneBoundaryHeat *>(Util::scene()->edges[i]->boundary);

            if (boundary && Util::scene()->edges[i]->boundary != Util::scene()->boundaries[0])
            {
                if (boundary->type == PhysicFieldBC_Heat_Flux)
                {
                    // vector flux term
                    double flux = boundary->heatFlux.number + boundary->h.number * boundary->externalTemperature.number;

                    if (fabs(flux) > EPS_ZERO)
                        add_vector_form_surf(new WeakFormsH1::SurfaceVectorForms::DefaultVectorFormSurf(0,
                                                                                                        QString::number(i + 1).toStdString(),
                                                                                                        flux,
                                                                                                        convertProblemType(Util::scene()->problemInfo()->problemType)));

                    if (fabs(boundary->h.number) > EPS_ZERO)
                        add_matrix_form_surf(new WeakFormsH1::SurfaceMatrixForms::DefaultMatrixFormSurf(0, 0,
                                                                                                        QString::number(i + 1).toStdString(),
                                                                                                        boundary->h.number,
                                                                                                        convertProblemType(Util::scene()->problemInfo()->problemType)));
                }
            }
        }

        // materials
        for (int i = 0; i<Util::scene()->labels.count(); i++)
        {
            SceneMaterialHeat *material = dynamic_cast<SceneMaterialHeat *>(Util::scene()->labels[i]->material);

            if (material && Util::scene()->labels[i]->material != Util::scene()->materials[0])
            {
                add_matrix_form(new WeakFormsH1::VolumetricMatrixForms::DefaultLinearDiffusion(0, 0,
                                                                                               QString::number(i).toStdString(),
                                                                                               material->thermal_conductivity.number,
                                                                                               HERMES_SYM,
                                                                                               convertProblemType(Util::scene()->problemInfo()->problemType)));

                if (fabs(material->volume_heat.number) > EPS_ZERO)
                    add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormConst(0,
                                                                                                   QString::number(i).toStdString(),
                                                                                                   material->volume_heat.number,
                                                                                                   convertProblemType(Util::scene()->problemInfo()->problemType)));

                // transient analysis
                if (Util::scene()->problemInfo()->analysisType == AnalysisType_Transient)
                {
                    if ((fabs(material->density.number) > EPS_ZERO) && (fabs(material->specific_heat.number) > EPS_ZERO))
                    {
                        if (solution.size() > 0)
                        {
                            add_matrix_form(new WeakFormsH1::VolumetricMatrixForms::DefaultLinearMass(0, 0,
                                                                                                      QString::number(i).toStdString(),
                                                                                                      material->density.number * material->specific_heat.number / Util::scene()->problemInfo()->timeStep.number,
                                                                                                      HERMES_SYM,
                                                                                                      convertProblemType(Util::scene()->problemInfo()->problemType)));

                            add_vector_form(new CustomVectorFormTimeDep(0,
                                                                        QString::number(i).toStdString(),
                                                                        material->density.number * material->specific_heat.number / Util::scene()->problemInfo()->timeStep.number,
                                                                        solution[0],
                                                                        convertProblemType(Util::scene()->problemInfo()->problemType)));
                        }
                    }
                }
            }
        }
    }
};

class CustomVectorFormTimeDep : public WeakForm::VectorFormVol
{
public:
    CustomVectorFormTimeDep(int i, scalar coeff, Solution *solution, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i), coeff(coeff), gt(gt)
    {
        ext.push_back(solution);
    }
    CustomVectorFormTimeDep(int i, std::string area, scalar coeff, Solution *solution, GeomType gt = HERMES_PLANAR)
        : WeakForm::VectorFormVol(i, area), coeff(coeff), gt(gt)
    {
        ext.push_back(solution);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                         Geom<double> *e, ExtData<scalar> *ext) const {
        if (gt == HERMES_PLANAR)
            return coeff * int_u_v<double, scalar>(n, wt, ext->fn[0], v);
        else if (gt == HERMES_AXISYM_X)
            return coeff * int_y_u_v<double, scalar>(n, wt, ext->fn[0], v, e);
        else
            return coeff * int_x_u_v<double, scalar>(n, wt, ext->fn[0], v, e);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e,
                    ExtData<Ord> *ext) const {
        if (gt == HERMES_PLANAR)
            return int_u_v<Ord, Ord>(n, wt, ext->fn[0], v);
        else if (gt == HERMES_AXISYM_X)
            return int_y_u_v<Ord, Ord>(n, wt, ext->fn[0], v, e);
        else
            return int_x_u_v<Ord, Ord>(n, wt, ext->fn[0], v, e);
    }

    // This is to make the form usable in rk_time_step().
    virtual WeakForm::VectorFormVol* clone() {
        return new CustomVectorFormTimeDep(*this);
    }

private:
    scalar coeff;
    GeomType gt;
};

*/

