/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_741375720241655743);
void inv_err_fun(double *nom_x, double *true_x, double *out_8668550649948152875);
void H_mod_fun(double *state, double *out_7434977771084928310);
void f_fun(double *state, double dt, double *out_8209496335212906999);
void F_fun(double *state, double dt, double *out_6984683106181723232);
void h_25(double *state, double *unused, double *out_8861255887784400674);
void H_25(double *state, double *unused, double *out_3571542409095679594);
void h_24(double *state, double *unused, double *out_3063426930001039074);
void H_24(double *state, double *unused, double *out_435574090307051598);
void h_30(double *state, double *unused, double *out_4999572411355899022);
void H_30(double *state, double *unused, double *out_8006030188871679802);
void h_26(double *state, double *unused, double *out_8221264446518327647);
void H_26(double *state, double *unused, double *out_2168251107584336086);
void h_27(double *state, double *unused, double *out_6169214896276080926);
void H_27(double *state, double *unused, double *out_7329938489690128998);
void h_29(double *state, double *unused, double *out_6444408958560586815);
void H_29(double *state, double *unused, double *out_8580588802727215977);
void h_28(double *state, double *unused, double *out_7780434171282851314);
void H_28(double *state, double *unused, double *out_6594980122181637024);
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
