/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_548503771068811561);
void inv_err_fun(double *nom_x, double *true_x, double *out_2682236651435864745);
void H_mod_fun(double *state, double *out_388372280572216428);
void f_fun(double *state, double dt, double *out_6950273804263183812);
void F_fun(double *state, double dt, double *out_2419952130103293792);
void h_25(double *state, double *unused, double *out_5363960631214508521);
void H_25(double *state, double *unused, double *out_2239868391431655836);
void h_24(double *state, double *unused, double *out_2062856328637710990);
void H_24(double *state, double *unused, double *out_8930133330717296648);
void h_30(double *state, double *unused, double *out_680766082583151051);
void H_30(double *state, double *unused, double *out_4288408685475712954);
void h_26(double *state, double *unused, double *out_5928987183223659744);
void H_26(double *state, double *unused, double *out_1633925218784246052);
void h_27(double *state, double *unused, double *out_1644625139496102996);
void H_27(double *state, double *unused, double *out_8331498373242644936);
void h_29(double *state, double *unused, double *out_1949514412564825142);
void H_29(double *state, double *unused, double *out_4075708323965073868);
void h_28(double *state, double *unused, double *out_7141644414365002210);
void H_28(double *state, double *unused, double *out_1069462091680102462);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
