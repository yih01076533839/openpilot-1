
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_741375720241655743) {
   out_741375720241655743[0] = delta_x[0] + nom_x[0];
   out_741375720241655743[1] = delta_x[1] + nom_x[1];
   out_741375720241655743[2] = delta_x[2] + nom_x[2];
   out_741375720241655743[3] = delta_x[3] + nom_x[3];
   out_741375720241655743[4] = delta_x[4] + nom_x[4];
   out_741375720241655743[5] = delta_x[5] + nom_x[5];
   out_741375720241655743[6] = delta_x[6] + nom_x[6];
   out_741375720241655743[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8668550649948152875) {
   out_8668550649948152875[0] = -nom_x[0] + true_x[0];
   out_8668550649948152875[1] = -nom_x[1] + true_x[1];
   out_8668550649948152875[2] = -nom_x[2] + true_x[2];
   out_8668550649948152875[3] = -nom_x[3] + true_x[3];
   out_8668550649948152875[4] = -nom_x[4] + true_x[4];
   out_8668550649948152875[5] = -nom_x[5] + true_x[5];
   out_8668550649948152875[6] = -nom_x[6] + true_x[6];
   out_8668550649948152875[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_7434977771084928310) {
   out_7434977771084928310[0] = 1.0;
   out_7434977771084928310[1] = 0.0;
   out_7434977771084928310[2] = 0.0;
   out_7434977771084928310[3] = 0.0;
   out_7434977771084928310[4] = 0.0;
   out_7434977771084928310[5] = 0.0;
   out_7434977771084928310[6] = 0.0;
   out_7434977771084928310[7] = 0.0;
   out_7434977771084928310[8] = 0.0;
   out_7434977771084928310[9] = 1.0;
   out_7434977771084928310[10] = 0.0;
   out_7434977771084928310[11] = 0.0;
   out_7434977771084928310[12] = 0.0;
   out_7434977771084928310[13] = 0.0;
   out_7434977771084928310[14] = 0.0;
   out_7434977771084928310[15] = 0.0;
   out_7434977771084928310[16] = 0.0;
   out_7434977771084928310[17] = 0.0;
   out_7434977771084928310[18] = 1.0;
   out_7434977771084928310[19] = 0.0;
   out_7434977771084928310[20] = 0.0;
   out_7434977771084928310[21] = 0.0;
   out_7434977771084928310[22] = 0.0;
   out_7434977771084928310[23] = 0.0;
   out_7434977771084928310[24] = 0.0;
   out_7434977771084928310[25] = 0.0;
   out_7434977771084928310[26] = 0.0;
   out_7434977771084928310[27] = 1.0;
   out_7434977771084928310[28] = 0.0;
   out_7434977771084928310[29] = 0.0;
   out_7434977771084928310[30] = 0.0;
   out_7434977771084928310[31] = 0.0;
   out_7434977771084928310[32] = 0.0;
   out_7434977771084928310[33] = 0.0;
   out_7434977771084928310[34] = 0.0;
   out_7434977771084928310[35] = 0.0;
   out_7434977771084928310[36] = 1.0;
   out_7434977771084928310[37] = 0.0;
   out_7434977771084928310[38] = 0.0;
   out_7434977771084928310[39] = 0.0;
   out_7434977771084928310[40] = 0.0;
   out_7434977771084928310[41] = 0.0;
   out_7434977771084928310[42] = 0.0;
   out_7434977771084928310[43] = 0.0;
   out_7434977771084928310[44] = 0.0;
   out_7434977771084928310[45] = 1.0;
   out_7434977771084928310[46] = 0.0;
   out_7434977771084928310[47] = 0.0;
   out_7434977771084928310[48] = 0.0;
   out_7434977771084928310[49] = 0.0;
   out_7434977771084928310[50] = 0.0;
   out_7434977771084928310[51] = 0.0;
   out_7434977771084928310[52] = 0.0;
   out_7434977771084928310[53] = 0.0;
   out_7434977771084928310[54] = 1.0;
   out_7434977771084928310[55] = 0.0;
   out_7434977771084928310[56] = 0.0;
   out_7434977771084928310[57] = 0.0;
   out_7434977771084928310[58] = 0.0;
   out_7434977771084928310[59] = 0.0;
   out_7434977771084928310[60] = 0.0;
   out_7434977771084928310[61] = 0.0;
   out_7434977771084928310[62] = 0.0;
   out_7434977771084928310[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_8209496335212906999) {
   out_8209496335212906999[0] = state[0];
   out_8209496335212906999[1] = state[1];
   out_8209496335212906999[2] = state[2];
   out_8209496335212906999[3] = state[3];
   out_8209496335212906999[4] = state[4];
   out_8209496335212906999[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8209496335212906999[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8209496335212906999[7] = state[7];
}
void F_fun(double *state, double dt, double *out_6984683106181723232) {
   out_6984683106181723232[0] = 1;
   out_6984683106181723232[1] = 0;
   out_6984683106181723232[2] = 0;
   out_6984683106181723232[3] = 0;
   out_6984683106181723232[4] = 0;
   out_6984683106181723232[5] = 0;
   out_6984683106181723232[6] = 0;
   out_6984683106181723232[7] = 0;
   out_6984683106181723232[8] = 0;
   out_6984683106181723232[9] = 1;
   out_6984683106181723232[10] = 0;
   out_6984683106181723232[11] = 0;
   out_6984683106181723232[12] = 0;
   out_6984683106181723232[13] = 0;
   out_6984683106181723232[14] = 0;
   out_6984683106181723232[15] = 0;
   out_6984683106181723232[16] = 0;
   out_6984683106181723232[17] = 0;
   out_6984683106181723232[18] = 1;
   out_6984683106181723232[19] = 0;
   out_6984683106181723232[20] = 0;
   out_6984683106181723232[21] = 0;
   out_6984683106181723232[22] = 0;
   out_6984683106181723232[23] = 0;
   out_6984683106181723232[24] = 0;
   out_6984683106181723232[25] = 0;
   out_6984683106181723232[26] = 0;
   out_6984683106181723232[27] = 1;
   out_6984683106181723232[28] = 0;
   out_6984683106181723232[29] = 0;
   out_6984683106181723232[30] = 0;
   out_6984683106181723232[31] = 0;
   out_6984683106181723232[32] = 0;
   out_6984683106181723232[33] = 0;
   out_6984683106181723232[34] = 0;
   out_6984683106181723232[35] = 0;
   out_6984683106181723232[36] = 1;
   out_6984683106181723232[37] = 0;
   out_6984683106181723232[38] = 0;
   out_6984683106181723232[39] = 0;
   out_6984683106181723232[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6984683106181723232[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6984683106181723232[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6984683106181723232[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6984683106181723232[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6984683106181723232[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6984683106181723232[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6984683106181723232[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6984683106181723232[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6984683106181723232[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6984683106181723232[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6984683106181723232[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6984683106181723232[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6984683106181723232[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6984683106181723232[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6984683106181723232[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6984683106181723232[56] = 0;
   out_6984683106181723232[57] = 0;
   out_6984683106181723232[58] = 0;
   out_6984683106181723232[59] = 0;
   out_6984683106181723232[60] = 0;
   out_6984683106181723232[61] = 0;
   out_6984683106181723232[62] = 0;
   out_6984683106181723232[63] = 1;
}
void h_25(double *state, double *unused, double *out_8861255887784400674) {
   out_8861255887784400674[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3571542409095679594) {
   out_3571542409095679594[0] = 0;
   out_3571542409095679594[1] = 0;
   out_3571542409095679594[2] = 0;
   out_3571542409095679594[3] = 0;
   out_3571542409095679594[4] = 0;
   out_3571542409095679594[5] = 0;
   out_3571542409095679594[6] = 1;
   out_3571542409095679594[7] = 0;
}
void h_24(double *state, double *unused, double *out_3063426930001039074) {
   out_3063426930001039074[0] = state[4];
   out_3063426930001039074[1] = state[5];
}
void H_24(double *state, double *unused, double *out_435574090307051598) {
   out_435574090307051598[0] = 0;
   out_435574090307051598[1] = 0;
   out_435574090307051598[2] = 0;
   out_435574090307051598[3] = 0;
   out_435574090307051598[4] = 1;
   out_435574090307051598[5] = 0;
   out_435574090307051598[6] = 0;
   out_435574090307051598[7] = 0;
   out_435574090307051598[8] = 0;
   out_435574090307051598[9] = 0;
   out_435574090307051598[10] = 0;
   out_435574090307051598[11] = 0;
   out_435574090307051598[12] = 0;
   out_435574090307051598[13] = 1;
   out_435574090307051598[14] = 0;
   out_435574090307051598[15] = 0;
}
void h_30(double *state, double *unused, double *out_4999572411355899022) {
   out_4999572411355899022[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8006030188871679802) {
   out_8006030188871679802[0] = 0;
   out_8006030188871679802[1] = 0;
   out_8006030188871679802[2] = 0;
   out_8006030188871679802[3] = 0;
   out_8006030188871679802[4] = 1;
   out_8006030188871679802[5] = 0;
   out_8006030188871679802[6] = 0;
   out_8006030188871679802[7] = 0;
}
void h_26(double *state, double *unused, double *out_8221264446518327647) {
   out_8221264446518327647[0] = state[7];
}
void H_26(double *state, double *unused, double *out_2168251107584336086) {
   out_2168251107584336086[0] = 0;
   out_2168251107584336086[1] = 0;
   out_2168251107584336086[2] = 0;
   out_2168251107584336086[3] = 0;
   out_2168251107584336086[4] = 0;
   out_2168251107584336086[5] = 0;
   out_2168251107584336086[6] = 0;
   out_2168251107584336086[7] = 1;
}
void h_27(double *state, double *unused, double *out_6169214896276080926) {
   out_6169214896276080926[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7329938489690128998) {
   out_7329938489690128998[0] = 0;
   out_7329938489690128998[1] = 0;
   out_7329938489690128998[2] = 0;
   out_7329938489690128998[3] = 1;
   out_7329938489690128998[4] = 0;
   out_7329938489690128998[5] = 0;
   out_7329938489690128998[6] = 0;
   out_7329938489690128998[7] = 0;
}
void h_29(double *state, double *unused, double *out_6444408958560586815) {
   out_6444408958560586815[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8580588802727215977) {
   out_8580588802727215977[0] = 0;
   out_8580588802727215977[1] = 1;
   out_8580588802727215977[2] = 0;
   out_8580588802727215977[3] = 0;
   out_8580588802727215977[4] = 0;
   out_8580588802727215977[5] = 0;
   out_8580588802727215977[6] = 0;
   out_8580588802727215977[7] = 0;
}
void h_28(double *state, double *unused, double *out_7780434171282851314) {
   out_7780434171282851314[0] = state[5];
   out_7780434171282851314[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6594980122181637024) {
   out_6594980122181637024[0] = 0;
   out_6594980122181637024[1] = 0;
   out_6594980122181637024[2] = 0;
   out_6594980122181637024[3] = 0;
   out_6594980122181637024[4] = 0;
   out_6594980122181637024[5] = 1;
   out_6594980122181637024[6] = 0;
   out_6594980122181637024[7] = 0;
   out_6594980122181637024[8] = 0;
   out_6594980122181637024[9] = 0;
   out_6594980122181637024[10] = 0;
   out_6594980122181637024[11] = 0;
   out_6594980122181637024[12] = 0;
   out_6594980122181637024[13] = 0;
   out_6594980122181637024[14] = 1;
   out_6594980122181637024[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
