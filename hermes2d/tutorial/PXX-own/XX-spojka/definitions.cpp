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
            node->values[0][0][i] = sqrt(sqr(uval[i]) + sqr(vval[i]));
            //node->values[0][0][i] = uval[i];
        }

        if(nodes->present(order)) {
          assert(nodes->get(order) == cur_node);
          ::free(nodes->get(order));
        }
        nodes->add(node, order);
        cur_node = node;
    }

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); return 0;}
};

class DoNothingFilter : public Filter
{
public:
    DoNothingFilter(MeshFunction* solution)
        : Filter()
    {
        num = 1;
        sln[0] = solution;
        init();
    }

    inline void precalculate(int order, int mask)
    {
        Quad2D* quad = quads[cur_quad];
        int np = quad->get_num_points(order);
        Node* node = new_node(H2D_FN_VAL_0, np);

        sln[0]->set_quad_order(order, H2D_FN_VAL);

        scalar *uval = sln[0]->get_fn_values();
        update_refmap();

        for (int i = 0; i < np; i++)
        {
            node->values[0][0][i] = uval[i];
        }

        if(nodes->present(order)) {
          assert(nodes->get(order) == cur_node);
          ::free(nodes->get(order));
        }
        nodes->add(node, order);
        cur_node = node;
    }

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); return 0;}
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

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); return 0; }
};

class WjFilter : public Filter
{
public:
    WjFilter(MeshFunction* slnr)
    {
        num = 1;

        sln[0] = slnr;


        init();
    }

    inline void precalculate(int order, int mask)
    {
        Quad2D* quad = quads[cur_quad];
        int np = quad->get_num_points(order);
        Node* node = new_node(H2D_FN_DEFAULT, np);

        double *dudx, *dudy;
        double *value;

        sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[0]->get_dx_dy_values(dudx, dudy);
        value = sln[0]->get_fn_values();

        update_refmap();

        double *x = refmap->get_phys_x(order);
        double *y = refmap->get_phys_y(order);
        Element *e = refmap->get_active_element();

        std::string user_marker = mesh->get_element_markers_conversion().get_user_marker(e->marker);
        int marker = atoi(user_marker.c_str());
        for (int i = 0; i < np; i++)
        {
            double vel_x = -magneticLabel[marker].velocity_angular * y[i];
            double vel_y = magneticLabel[marker].velocity_angular * x[i];

            if(magneticLabel[marker].conductivity > 0.0){
                node->values[0][0][i] = 1.0 / magneticLabel[marker].conductivity * (
                    sqr(std::abs(-magneticLabel[marker].conductivity * (vel_x * dudx[i] + vel_y * dudy[i]))));
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

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); return 0;}
};

class BFilter : public Filter
{
public:
    BFilter(MeshFunction* slnr)
    {
        num = 1;

        sln[0] = slnr;


        init();
    }

    inline void precalculate(int order, int mask)
    {
        Quad2D* quad = quads[cur_quad];
        int np = quad->get_num_points(order);
        Node* node = new_node(H2D_FN_DEFAULT, np);

        double *dudx1, *dudy1;

        sln[0]->set_quad_order(order, H2D_FN_VAL | H2D_FN_DX | H2D_FN_DY);
        sln[0]->get_dx_dy_values(dudx1, dudy1);

        update_refmap();

        double *x = refmap->get_phys_x(order);
        Element *e = refmap->get_active_element();

        for (int i = 0; i < np; i++)
        {
                node->values[0][0][i] = sqrt(sqr(dudx1[i]) + sqr(dudy1[i]));
        }

        if(nodes->present(order)) {
          assert(nodes->get(order) == cur_node);
          ::free(nodes->get(order));
        }
        nodes->add(node, order);
        cur_node = node;
    }

    inline scalar get_pt_value(double x, double y, int item) { error("Not implemented"); return 0;}
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


double prev_temp_set;  //TODO dat dovnitr, tady je to osklivy

// Nonlinear parameter, supposes that we consider only positive parameters !!!
const double NONLINEAR_PARAMETER = -1.;

class WeakFormMagnetic : public WeakForm
{
public:
    WeakFormMagnetic(int neq) : WeakForm(neq) { }

    void registerForms(Hermes::vector<int> labels, Solution *prev_mag_r_sln, Filter *prev_temp_sln)
    {
        for(std::vector<int>::iterator it = labels.begin(); it != labels.end(); ++it) {
            double perm = magneticLabel[*it].permeability;
            if (USE_NONLINEARITIES && (zelezoLabels.find_index(*it, false) != -1))
                perm = NONLINEAR_PARAMETER;

            double cond = magneticLabel[*it].conductivity;
            if (USE_NONLINEARITIES && (zelezoLabels.find_index(*it, false) != -1))
                cond = NONLINEAR_PARAMETER;

            double omega = magneticLabel[*it].velocity_angular;
            double reman = magneticLabel[*it].remanence;
            double reman_angle = magneticLabel[*it].remanence_angle;

            add_magnetic_material(str_marker[*it], perm, cond, magneticLabel[*it].current_density_real, prev_mag_r_sln, prev_temp_sln, reman, reman_angle, omega);
        }
        prev_temp_set = false;
    }

    void add_magnetic_material(std::string marker, double permeability, double conductivity, double external_current_density, Solution *prev_mag_r_sln,
                               Filter *prev_temp_filter, double remanence, double remanence_angle, double omega_velocity)
    {
        /// TODO PROC TO PADA, KDYZ DAM NASLEDUJICIM FORMAM HERMES_SYM ???????????
        /// TODO PROC TO PADA, KDYZ DAM NASLEDUJICIM FORMAM HERMES_SYM ???????????
        /// TODO PROC TO PADA, KDYZ DAM NASLEDUJICIM FORMAM HERMES_SYM ???????????

        // real part
        CustomMatrixFormVol *mat_form_real = new CustomMatrixFormVol(0, 0, marker, permeability);
        mat_form_real->ext.push_back(prev_mag_r_sln);
        // VK // mat_form_real->ext.push_back(prev_mag_i_sln);
        add_matrix_form(mat_form_real);

        CustomMatrixFormVolVel *mat_form_vel = new CustomMatrixFormVolVel(0, 0, marker, conductivity, omega_velocity);
        add_matrix_form(mat_form_vel);



// VK       // conductivity pouze u harmonického pole
//        if (fabs(conductivity) > EPS_ZERO)
//        {
////            add_matrix_form(new WeakFormsH1::VolumetricMatrixForms::DefaultLinearMass(0, 1, marker, - 2 * M_PI * frequency * conductivity, HERMES_NONSYM, HERMES_PLANAR));
////            add_matrix_form(new WeakFormsH1::VolumetricMatrixForms::DefaultLinearMass(1, 0, marker,   2 * M_PI * frequency * conductivity, HERMES_NONSYM, HERMES_PLANAR));
//            add_matrix_form(new CustomMatrixFormCond(0, 1, marker, - 2 * M_PI * frequency, conductivity));
//            add_matrix_form(new CustomMatrixFormCond(1, 0, marker,   2 * M_PI * frequency, conductivity));
// VK       }

        // external current density
        if (fabs(external_current_density) > EPS_ZERO)
            add_vector_form(new WeakFormsH1::VolumetricVectorForms::DefaultVectorFormConst(0, marker, external_current_density, HERMES_PLANAR));

        CustomVectorFormVolReman *vec_form_reman = new CustomVectorFormVolReman(0, marker, permeability, remanence, remanence_angle);
        vec_form_reman->ext.push_back(prev_mag_r_sln);
        add_vector_form(vec_form_reman);

    }

    void push_previous_temperature(Solution *prev_temp_sln)
    {
        for(std::vector<VectorFormVol *>::iterator it = vfvol.begin(); it != vfvol.end(); ++it) {
            (*it)->ext.push_back(prev_temp_sln);
        }

        for(std::vector<MatrixFormVol *>::iterator it = mfvol.begin(); it != mfvol.end(); ++it) {
            (*it)->ext.push_back(prev_temp_sln);
        }
        prev_temp_set = true;
    }

private:
    class CustomMatrixFormVol : public WeakForm::MatrixFormVol
    {
    public:
        CustomMatrixFormVol(int i, int j, std::string marker, double permeability_const)
            : WeakForm::MatrixFormVol(i, j, marker, HERMES_NONSYM), permeability_const(permeability_const){}

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
           // return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);

            Func<double>* sln_mag_r_prev = ext->fn[0];
           //VK // Func<double>* sln_mag_i_prev = ext->fn[1];
            Func<double>* sln_temp_prev = ext->fn[1];

            scalar result = 0;
            for (int i = 0; i < n; i++){
                scalar B = sqrt(sqr(sln_mag_r_prev->dx[i]) + sqr(sln_mag_r_prev->dy[i]));
                if(B>maxB) maxB = B;
                scalar T = (prev_temp_set) ? sln_temp_prev->val[i] : TEMP_INIT;
                scalar permeability = (permeability_const == NONLINEAR_PARAMETER) ? permeability_function(B,T) : permeability_const;

                result += wt[i] / (MU0 * permeability) * (
                            u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
            }
            return result;

        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            //return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
            return Ord(20);
        }

    private:
        double permeability_const;
    };

    // ******** Form for velocity (matrix volume) ***********
private:
    class CustomMatrixFormVolVel : public WeakForm::MatrixFormVol
    {
    public:
        CustomMatrixFormVolVel(int i, int j, std::string marker, double conductivity_const, double omega_velocity)
            : WeakForm::MatrixFormVol(i, j, marker, HERMES_NONSYM), conductivity_const(conductivity_const), omega_velocity(omega_velocity) {}

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
           // return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);

            Func<double>* sln_temp_prev = ext->fn[0];

            scalar result = 0;
            for (int i = 0; i < n; i++){
                scalar T = (prev_temp_set) ? sln_temp_prev->val[i] : TEMP_INIT;
                scalar conductivity = (conductivity_const == NONLINEAR_PARAMETER) ? electric_conductivity_fe.value(T) : conductivity_const;
                if(conductivity > max_el_cond) max_el_cond = conductivity;
                if(conductivity < min_el_cond) min_el_cond = conductivity;

                scalar vel_x = - omega_velocity * e->y[i];
                scalar vel_y = omega_velocity * e->x[i];

                result += wt[i] * conductivity * (
                            u->val[i] * vel_x * v->dx[i] + u->val[i] * vel_y * v->dy[i]);
            }
            return result;

        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            //return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
            return Ord(20);
        }

    private:
        double conductivity_const, omega_velocity;
    };


    // ******** Form for remanence of magnets (vector volume) ****************
private:
    class CustomVectorFormVolReman : public WeakForm::VectorFormVol
    {
    public:
        CustomVectorFormVolReman(int i, std::string marker, double permeability_const, double remanence, double remanence_angle)
            :WeakForm::VectorFormVol(i, marker), permeability_const(permeability_const), remanence(remanence), remanence_angle(remanence_angle) {}

        template<typename Real, typename Scalar>
        Scalar vector_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {

        }


        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {

            Func<double>* sln_mag_r_prev = ext->fn[0];
            Func<double>* sln_temp_prev = ext->fn[1];

            scalar result = 0;
            for (int i = 0; i < n; i++){
                scalar B = sqrt(sqr(sln_mag_r_prev->dx[i]) + sqr(sln_mag_r_prev->dy[i]));
                if(B>maxB) maxB = B;
                scalar T = (prev_temp_set) ? sln_temp_prev->val[i] : TEMP_INIT;
                scalar permeability = (permeability_const == NONLINEAR_PARAMETER) ? permeability_function(B,T) : permeability_const;

                scalar reman_x = remanence * cos(remanence_angle);
                scalar reman_y = remanence * sin(remanence_angle);

                result += - wt[i] / (MU0 * permeability) * (
                            reman_y * v->dx[i] - reman_x * v->dy[i]);
            }
            return result;

        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            //return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
            //return Ord(20);
            return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext); //nevim proc mas nekde navrat a nekde Ord(20)
        }

    private:
        double permeability_const, remanence, remanence_angle;
    };


};




class WeakFormTemp : public WeakForm
{
public:
    WeakFormTemp(double time_step) : WeakForm(1), time_step(time_step) { };
    void registerForms(Hermes::vector<int> labels, Solution *prev_time_sln, Filter *joule_loses)
    {
        for(std::vector<int>::iterator it = labels.begin(); it != labels.end(); ++it) {
            double cond = heatLabel[*it].thermal_conductivity;
            if (USE_NONLINEARITIES && (zelezoLabels.find_index(*it, false) != -1))
                cond = NONLINEAR_PARAMETER;

            add_temperature_material(str_marker[*it], cond, heatLabel[*it].volume_heat,
                                     heatLabel[*it].density, heatLabel[*it].specific_heat, prev_time_sln, joule_loses);
        }

        for(int i=0; i<NUM_EDGES; i++){
            if(heatEdge[i].type == PhysicFieldBC_Heat_Flux){
                add_temperature_edge(str_marker[i], heatEdge[i].heatFlux, heatEdge[i].h, heatEdge[i].externalTemperature);
            }
        }
    }

    void add_temperature_edge(std::string marker, double flux, double h, double ext_temp)
    {
        add_matrix_form_surf(new DefaultMatrixFormSurf(0, 0, marker, h, HERMES_PLANAR));
        add_vector_form_surf(new DefaultVectorFormSurf(0, marker, flux + h * ext_temp, HERMES_PLANAR));
    }

    void add_temperature_material(std::string marker, double thermal_conductivity, double volume_heat, double density, double specific_heat, Solution* prev_time_sln, Filter* joule_loses)
    {
        /// TODO vsude jsem smazal 2 * M_PI

        // Contribution of the time derivative term.
        /// TODO NONSYM
        add_matrix_form(new DefaultLinearMass(0, 0, marker, density * specific_heat / time_step, HERMES_SYM, HERMES_PLANAR));

        // Contribution of the diffusion term.
        /// TODO NONSYM
        //add_matrix_form(new DefaultLinearDiffusion(0, 0, marker, thermal_conductivity, HERMES_SYM, HERMES_AXISYM_Y));
        CustomNonlinearDiffusion* diffusion_form = new CustomNonlinearDiffusion(0, 0, marker, thermal_conductivity);
        diffusion_form->ext.push_back(prev_time_sln);
        add_matrix_form(diffusion_form);

        // Right-hand side volumetric vector form.
        CustomVectorFormVol* vec_form_vol = new CustomVectorFormVol(0, marker, density, specific_heat, volume_heat, time_step);
        vec_form_vol->ext.push_back(prev_time_sln);
        vec_form_vol->ext.push_back(joule_loses);
        add_vector_form(vec_form_vol);
    }

private:
    double time_step;

    class CustomNonlinearDiffusion : public WeakForm::MatrixFormVol
    {
    public:
        CustomNonlinearDiffusion(int i, int j, std::string marker, double thermal_conductivity_const)
            : WeakForm::MatrixFormVol(i, j, marker, HERMES_SYM), thermal_conductivity_const(thermal_conductivity_const) {};

        template<typename Real, typename Scalar>
        Scalar matrix_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
        }

        virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
            //return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
            scalar result = 0;
            Func<double>* temp_prev_time = ext->fn[0];

            for (int i = 0; i < n; i++){
                double temperature = temp_prev_time->val[i];
                double thermal_conductivity = (thermal_conductivity_const == NONLINEAR_PARAMETER) ? thermal_conductivity_fe.value(temperature) : thermal_conductivity_const;
                //result += wt[i] * thermal_conductivity * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
                result += wt[i] * thermal_conductivity * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
            }
            return result;
        }

        virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
            //return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
            return Ord(20);
        }
    private:
        double thermal_conductivity_const;
    };

    // This form is custom since it contains previous time-level solution.
    class CustomVectorFormVol : public WeakForm::VectorFormVol
    {
    public:
        CustomVectorFormVol(int i, std::string marker, double density, double specific_heat, double volume_heat, double time_step)
        : WeakForm::VectorFormVol(i, marker), density(density), specific_heat(specific_heat), volume_heat(volume_heat), time_step(time_step) { }

      template<typename Real, typename Scalar>
      Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
          Func<Real>* temp_prev_time = ext->fn[0];
          Func<Real>* joule_loses = ext->fn[1];

          /// TODO pridat 2 * M_PI ?????????
           /// JE to OK????
//          return int_x_u_v<Real, Scalar>(n, wt, temp_prev_time, v, e)  * density * specific_heat / time_step +
//                 int_x_u_v<Real, Scalar>(n, wt, joule_loses, v, e)  * volume_heat;

          return int_u_v<Real, Scalar>(n, wt, temp_prev_time, v)  * density * specific_heat / time_step +
                 int_u_v<Real, Scalar>(n, wt, joule_loses, v)  * volume_heat;

      }

      virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
        return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
      }

      virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
        return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
      }

      double density, specific_heat, volume_heat, time_step;
    };
};

class WeakFormElast : public WeakForm
{
public:
    WeakFormElast() : WeakForm(2) {}

    void register_forms(Hermes::vector<int> labels, Solution * temperature){
        for(std::vector<int>::iterator it = labels.begin(); it != labels.end(); ++it) {
            add_elasticity_material(str_marker[*it], elasticityLabel[*it].lambda(),  elasticityLabel[*it].mu(), elasticityLabel[*it].thermal_expansion, temperature);
        }
    }
    
    void add_elasticity_material(std::string marker, double lambda, double mu, double thermal_expansion, Solution* temperature){
        CustomMatrixFormXX* mat_form_x_x = new CustomMatrixFormXX(0, 0, marker, lambda, mu);
        mat_form_x_x->ext.push_back(temperature);
        add_matrix_form(mat_form_x_x);
        
        CustomMatrixFormXY* mat_form_x_y = new CustomMatrixFormXY(0, 1, marker, lambda, mu);
        add_matrix_form(mat_form_x_y);

        CustomMatrixFormYY* mat_form_y_y = new CustomMatrixFormYY(1, 1, marker, lambda, mu);
        add_matrix_form(mat_form_y_y);

        CustomVectorFormX* vec_form_x = new CustomVectorFormX(0, marker, lambda, mu, thermal_expansion);
        vec_form_x->ext.push_back(temperature);
        add_vector_form(vec_form_x);

        CustomVectorFormY* vec_form_y = new CustomVectorFormY(1, marker, lambda, mu, thermal_expansion);
        vec_form_y->ext.push_back(temperature);
        add_vector_form(vec_form_y);
    }

private:
    class CustomMatrixFormXX : public WeakForm::MatrixFormVol
    {
    public:
         CustomMatrixFormXX(int i, int j, std::string marker, double lambda, double mu)
             : WeakForm::MatrixFormVol(i, j, marker, HERMES_NONSYM), lambda(lambda), mu(mu) {}

         template<typename Real, typename Scalar>
         Scalar matrix_form_x_x(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
             Scalar result = 0;

             for (int i = 0; i < n; i++)
                 result += wt[i] * ((2.0 * mu + lambda) * (u->dx[i] * v->dx[i]) +
                                    mu * ( u->dy[i] * v->dy[i]));
             return result;
         }

         virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
           return matrix_form_x_x<double, scalar>(n, wt, u_ext, u, v, e, ext);
         }

         virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
           return matrix_form_x_x<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
         }

     private:
         double lambda, mu;
    };

    class CustomMatrixFormXY : public WeakForm::MatrixFormVol
    {
    public:
         CustomMatrixFormXY(int i, int j, std::string marker, double lambda, double mu)
             : WeakForm::MatrixFormVol(i, j, marker, HERMES_SYM), lambda(lambda), mu(mu) {}

         template<typename Real, typename Scalar>
         Scalar matrix_form_x_y(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
             Scalar result = 0;
             for (int i = 0; i < n; i++)
                 result += wt[i] * (lambda * (u->dx[i] * v->dy[i]) +
                                    mu * (u->dy[i] * v->dx[i]));
             return result;
         }

         virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
           return matrix_form_x_y<double, scalar>(n, wt, u_ext, u, v, e, ext);
         }

         virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
           return matrix_form_x_y<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
         }

     private:
         double lambda, mu;
    };

    class CustomMatrixFormYY : public WeakForm::MatrixFormVol
    {
    public:
         CustomMatrixFormYY(int i, int j, std::string marker, double lambda, double mu)
             : WeakForm::MatrixFormVol(i, j, marker, HERMES_NONSYM), lambda(lambda), mu(mu) {}

         template<typename Real, typename Scalar>
         Scalar matrix_form_y_y(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
             Scalar result = 0;
             for (int i = 0; i < n; i++)
                 result += wt[i] * ((2 * mu + lambda) * (u->dy[i] * v->dy[i]) +
                                    mu * (u->dx[i] * v->dx[i]));
             return result;
         }

         virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
           return matrix_form_y_y<double, scalar>(n, wt, u_ext, u, v, e, ext);
         }

         virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
           return matrix_form_y_y<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
         }

     private:
         double lambda, mu;
    };

    class CustomVectorFormX : public WeakForm::VectorFormVol
    {
    public:
         CustomVectorFormX(int i, std::string marker, double lambda, double mu, double expansion)
             : WeakForm::VectorFormVol(i, marker), lambda(lambda), mu(mu), expansion(expansion) {}

         template<typename Real, typename Scalar>
         Scalar vector_form_x(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
             Scalar result = 0;
             Func<Real>* temperature = ext->fn[0];
             for (int i = 0; i < n; i++)
                 result += wt[i] * (3*lambda + 2*mu) * expansion * (
                            // temperature->dx[i] * v->val[i]
                              (temperature->val[i] - TEMP_INIT) * v->dx[i] //  nejsem si jisty jestli ma byt i zbytek + v->val[i])
                             );
                 //result += wt[i] * (100000000) * e->x[i] * v->val[i];
                 //result += 0;
             return result;
         }

         virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
           return vector_form_x<double, scalar>(n, wt, u_ext, v, e, ext);
         }

         virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
           return vector_form_x<Ord, Ord>(n, wt, u_ext, v, e, ext);
         }

     private:
         double lambda, mu, expansion;
    };

    class CustomVectorFormY : public WeakForm::VectorFormVol
    {
    public:
         CustomVectorFormY(int i, std::string marker, double lambda, double mu, double expansion)
             : WeakForm::VectorFormVol(i, marker), lambda(lambda), mu(mu), expansion(expansion) {}

         template<typename Real, typename Scalar>
         Scalar vector_form_y(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
             Scalar result = 0;
             Func<Real>* temperature = ext->fn[0];
             for (int i = 0; i < n; i++)
                 result += wt[i] * (3*lambda + 2*mu) * expansion * (
                            // temperature->dy[i] * v->val[i]
                              (temperature->val[i] - TEMP_INIT) * v->dy[i]
                             );
                 //result += wt[i] * (-1000000000) * e->x[i] * v->val[i];
                 //result += 0;
             return result;
         }

         virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
           return vector_form_y<double, scalar>(n, wt, u_ext, v, e, ext);
         }

         virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
           return vector_form_y<Ord, Ord>(n, wt, u_ext, v, e, ext);
         }

     private:
         double lambda, mu, expansion;
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

