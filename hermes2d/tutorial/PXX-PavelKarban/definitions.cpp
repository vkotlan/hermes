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

        for (int i = 0; i < np; i++)
        {
            /*
            node->values[0][0][i] = (magneticLabel[e->marker].conductivity > 0.0) ?
                        0.5 / magneticLabel[e->marker].conductivity * (
                            sqr(2 * M_PI * frequency * magneticLabel[e->marker].conductivity * value2[i]) +
                            sqr(2 * M_PI * frequency * magneticLabel[e->marker].conductivity * value1[i]))
                      :
                        0.0;
            */
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
        add_magnetic_material("0", 700.0, 1e6, 0.0);
        add_magnetic_material("1", 1.0, 0.0, 6e7);
        add_magnetic_material("2", 1.0, 0.0, 6e7);
        add_magnetic_material("3", 1.0, 0.0, 0.0);
        add_magnetic_material("4", 1.0, 0.0, 0.0);
        add_magnetic_material("5", 1.0, 0.0, 0.0);
        add_magnetic_material("6", 1.0, 0.0, 0.0);
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

