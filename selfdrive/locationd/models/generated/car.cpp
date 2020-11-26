
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
void err_fun(double *nom_x, double *delta_x, double *out_4793795122995228942) {
   out_4793795122995228942[0] = delta_x[0] + nom_x[0];
   out_4793795122995228942[1] = delta_x[1] + nom_x[1];
   out_4793795122995228942[2] = delta_x[2] + nom_x[2];
   out_4793795122995228942[3] = delta_x[3] + nom_x[3];
   out_4793795122995228942[4] = delta_x[4] + nom_x[4];
   out_4793795122995228942[5] = delta_x[5] + nom_x[5];
   out_4793795122995228942[6] = delta_x[6] + nom_x[6];
   out_4793795122995228942[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1664103941613909424) {
   out_1664103941613909424[0] = -nom_x[0] + true_x[0];
   out_1664103941613909424[1] = -nom_x[1] + true_x[1];
   out_1664103941613909424[2] = -nom_x[2] + true_x[2];
   out_1664103941613909424[3] = -nom_x[3] + true_x[3];
   out_1664103941613909424[4] = -nom_x[4] + true_x[4];
   out_1664103941613909424[5] = -nom_x[5] + true_x[5];
   out_1664103941613909424[6] = -nom_x[6] + true_x[6];
   out_1664103941613909424[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_212441384176172332) {
   out_212441384176172332[0] = 1.0;
   out_212441384176172332[1] = 0.0;
   out_212441384176172332[2] = 0.0;
   out_212441384176172332[3] = 0.0;
   out_212441384176172332[4] = 0.0;
   out_212441384176172332[5] = 0.0;
   out_212441384176172332[6] = 0.0;
   out_212441384176172332[7] = 0.0;
   out_212441384176172332[8] = 0.0;
   out_212441384176172332[9] = 1.0;
   out_212441384176172332[10] = 0.0;
   out_212441384176172332[11] = 0.0;
   out_212441384176172332[12] = 0.0;
   out_212441384176172332[13] = 0.0;
   out_212441384176172332[14] = 0.0;
   out_212441384176172332[15] = 0.0;
   out_212441384176172332[16] = 0.0;
   out_212441384176172332[17] = 0.0;
   out_212441384176172332[18] = 1.0;
   out_212441384176172332[19] = 0.0;
   out_212441384176172332[20] = 0.0;
   out_212441384176172332[21] = 0.0;
   out_212441384176172332[22] = 0.0;
   out_212441384176172332[23] = 0.0;
   out_212441384176172332[24] = 0.0;
   out_212441384176172332[25] = 0.0;
   out_212441384176172332[26] = 0.0;
   out_212441384176172332[27] = 1.0;
   out_212441384176172332[28] = 0.0;
   out_212441384176172332[29] = 0.0;
   out_212441384176172332[30] = 0.0;
   out_212441384176172332[31] = 0.0;
   out_212441384176172332[32] = 0.0;
   out_212441384176172332[33] = 0.0;
   out_212441384176172332[34] = 0.0;
   out_212441384176172332[35] = 0.0;
   out_212441384176172332[36] = 1.0;
   out_212441384176172332[37] = 0.0;
   out_212441384176172332[38] = 0.0;
   out_212441384176172332[39] = 0.0;
   out_212441384176172332[40] = 0.0;
   out_212441384176172332[41] = 0.0;
   out_212441384176172332[42] = 0.0;
   out_212441384176172332[43] = 0.0;
   out_212441384176172332[44] = 0.0;
   out_212441384176172332[45] = 1.0;
   out_212441384176172332[46] = 0.0;
   out_212441384176172332[47] = 0.0;
   out_212441384176172332[48] = 0.0;
   out_212441384176172332[49] = 0.0;
   out_212441384176172332[50] = 0.0;
   out_212441384176172332[51] = 0.0;
   out_212441384176172332[52] = 0.0;
   out_212441384176172332[53] = 0.0;
   out_212441384176172332[54] = 1.0;
   out_212441384176172332[55] = 0.0;
   out_212441384176172332[56] = 0.0;
   out_212441384176172332[57] = 0.0;
   out_212441384176172332[58] = 0.0;
   out_212441384176172332[59] = 0.0;
   out_212441384176172332[60] = 0.0;
   out_212441384176172332[61] = 0.0;
   out_212441384176172332[62] = 0.0;
   out_212441384176172332[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_3714352266407884559) {
   out_3714352266407884559[0] = state[0];
   out_3714352266407884559[1] = state[1];
   out_3714352266407884559[2] = state[2];
   out_3714352266407884559[3] = state[3];
   out_3714352266407884559[4] = state[4];
   out_3714352266407884559[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3714352266407884559[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3714352266407884559[7] = state[7];
}
void F_fun(double *state, double dt, double *out_4521413872731011395) {
   out_4521413872731011395[0] = 1;
   out_4521413872731011395[1] = 0;
   out_4521413872731011395[2] = 0;
   out_4521413872731011395[3] = 0;
   out_4521413872731011395[4] = 0;
   out_4521413872731011395[5] = 0;
   out_4521413872731011395[6] = 0;
   out_4521413872731011395[7] = 0;
   out_4521413872731011395[8] = 0;
   out_4521413872731011395[9] = 1;
   out_4521413872731011395[10] = 0;
   out_4521413872731011395[11] = 0;
   out_4521413872731011395[12] = 0;
   out_4521413872731011395[13] = 0;
   out_4521413872731011395[14] = 0;
   out_4521413872731011395[15] = 0;
   out_4521413872731011395[16] = 0;
   out_4521413872731011395[17] = 0;
   out_4521413872731011395[18] = 1;
   out_4521413872731011395[19] = 0;
   out_4521413872731011395[20] = 0;
   out_4521413872731011395[21] = 0;
   out_4521413872731011395[22] = 0;
   out_4521413872731011395[23] = 0;
   out_4521413872731011395[24] = 0;
   out_4521413872731011395[25] = 0;
   out_4521413872731011395[26] = 0;
   out_4521413872731011395[27] = 1;
   out_4521413872731011395[28] = 0;
   out_4521413872731011395[29] = 0;
   out_4521413872731011395[30] = 0;
   out_4521413872731011395[31] = 0;
   out_4521413872731011395[32] = 0;
   out_4521413872731011395[33] = 0;
   out_4521413872731011395[34] = 0;
   out_4521413872731011395[35] = 0;
   out_4521413872731011395[36] = 1;
   out_4521413872731011395[37] = 0;
   out_4521413872731011395[38] = 0;
   out_4521413872731011395[39] = 0;
   out_4521413872731011395[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4521413872731011395[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4521413872731011395[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4521413872731011395[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4521413872731011395[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4521413872731011395[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4521413872731011395[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4521413872731011395[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4521413872731011395[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4521413872731011395[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4521413872731011395[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4521413872731011395[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4521413872731011395[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4521413872731011395[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4521413872731011395[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4521413872731011395[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4521413872731011395[56] = 0;
   out_4521413872731011395[57] = 0;
   out_4521413872731011395[58] = 0;
   out_4521413872731011395[59] = 0;
   out_4521413872731011395[60] = 0;
   out_4521413872731011395[61] = 0;
   out_4521413872731011395[62] = 0;
   out_4521413872731011395[63] = 1;
}
void h_25(double *state, double *unused, double *out_9010470870943443566) {
   out_9010470870943443566[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1174950662950397075) {
   out_1174950662950397075[0] = 0;
   out_1174950662950397075[1] = 0;
   out_1174950662950397075[2] = 0;
   out_1174950662950397075[3] = 0;
   out_1174950662950397075[4] = 0;
   out_1174950662950397075[5] = 0;
   out_1174950662950397075[6] = 1;
   out_1174950662950397075[7] = 0;
}
void h_24(double *state, double *unused, double *out_2634279382811077010) {
   out_2634279382811077010[0] = state[4];
   out_2634279382811077010[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3293235470661685495) {
   out_3293235470661685495[0] = 0;
   out_3293235470661685495[1] = 0;
   out_3293235470661685495[2] = 0;
   out_3293235470661685495[3] = 0;
   out_3293235470661685495[4] = 1;
   out_3293235470661685495[5] = 0;
   out_3293235470661685495[6] = 0;
   out_3293235470661685495[7] = 0;
   out_3293235470661685495[8] = 0;
   out_3293235470661685495[9] = 0;
   out_3293235470661685495[10] = 0;
   out_3293235470661685495[11] = 0;
   out_3293235470661685495[12] = 0;
   out_3293235470661685495[13] = 1;
   out_3293235470661685495[14] = 0;
   out_3293235470661685495[15] = 0;
}
void h_30(double *state, double *unused, double *out_3695862822254458359) {
   out_3695862822254458359[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8438948247998786205) {
   out_8438948247998786205[0] = 0;
   out_8438948247998786205[1] = 0;
   out_8438948247998786205[2] = 0;
   out_8438948247998786205[3] = 0;
   out_8438948247998786205[4] = 1;
   out_8438948247998786205[5] = 0;
   out_8438948247998786205[6] = 0;
   out_8438948247998786205[7] = 0;
}
void h_26(double *state, double *unused, double *out_1838395483279779233) {
   out_1838395483279779233[0] = state[7];
}
void H_26(double *state, double *unused, double *out_4564842853729618605) {
   out_4564842853729618605[0] = 0;
   out_4564842853729618605[1] = 0;
   out_4564842853729618605[2] = 0;
   out_4564842853729618605[3] = 0;
   out_4564842853729618605[4] = 0;
   out_4564842853729618605[5] = 0;
   out_4564842853729618605[6] = 0;
   out_4564842853729618605[7] = 1;
}
void h_27(double *state, double *unused, double *out_8990603157829754913) {
   out_8990603157829754913[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8720213837874140099) {
   out_8720213837874140099[0] = 0;
   out_8720213837874140099[1] = 0;
   out_8720213837874140099[2] = 0;
   out_8720213837874140099[3] = 1;
   out_8720213837874140099[4] = 0;
   out_8720213837874140099[5] = 0;
   out_8720213837874140099[6] = 0;
   out_8720213837874140099[7] = 0;
}
void h_29(double *state, double *unused, double *out_9180946853595290814) {
   out_9180946853595290814[0] = state[1];
}
void H_29(double *state, double *unused, double *out_6183997056581933458) {
   out_6183997056581933458[0] = 0;
   out_6183997056581933458[1] = 1;
   out_6183997056581933458[2] = 0;
   out_6183997056581933458[3] = 0;
   out_6183997056581933458[4] = 0;
   out_6183997056581933458[5] = 0;
   out_6183997056581933458[6] = 0;
   out_6183997056581933458[7] = 0;
}
void h_28(double *state, double *unused, double *out_4490002131666605758) {
   out_4490002131666605758[0] = state[5];
   out_4490002131666605758[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3724597007574809371) {
   out_3724597007574809371[0] = 0;
   out_3724597007574809371[1] = 0;
   out_3724597007574809371[2] = 0;
   out_3724597007574809371[3] = 0;
   out_3724597007574809371[4] = 0;
   out_3724597007574809371[5] = 1;
   out_3724597007574809371[6] = 0;
   out_3724597007574809371[7] = 0;
   out_3724597007574809371[8] = 0;
   out_3724597007574809371[9] = 0;
   out_3724597007574809371[10] = 0;
   out_3724597007574809371[11] = 0;
   out_3724597007574809371[12] = 0;
   out_3724597007574809371[13] = 0;
   out_3724597007574809371[14] = 1;
   out_3724597007574809371[15] = 0;
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
