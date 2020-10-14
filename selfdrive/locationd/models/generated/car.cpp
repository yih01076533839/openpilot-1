
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
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_548503771068811561) {
   out_548503771068811561[0] = delta_x[0] + nom_x[0];
   out_548503771068811561[1] = delta_x[1] + nom_x[1];
   out_548503771068811561[2] = delta_x[2] + nom_x[2];
   out_548503771068811561[3] = delta_x[3] + nom_x[3];
   out_548503771068811561[4] = delta_x[4] + nom_x[4];
   out_548503771068811561[5] = delta_x[5] + nom_x[5];
   out_548503771068811561[6] = delta_x[6] + nom_x[6];
   out_548503771068811561[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2682236651435864745) {
   out_2682236651435864745[0] = -nom_x[0] + true_x[0];
   out_2682236651435864745[1] = -nom_x[1] + true_x[1];
   out_2682236651435864745[2] = -nom_x[2] + true_x[2];
   out_2682236651435864745[3] = -nom_x[3] + true_x[3];
   out_2682236651435864745[4] = -nom_x[4] + true_x[4];
   out_2682236651435864745[5] = -nom_x[5] + true_x[5];
   out_2682236651435864745[6] = -nom_x[6] + true_x[6];
   out_2682236651435864745[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_388372280572216428) {
   out_388372280572216428[0] = 1.0;
   out_388372280572216428[1] = 0.0;
   out_388372280572216428[2] = 0.0;
   out_388372280572216428[3] = 0.0;
   out_388372280572216428[4] = 0.0;
   out_388372280572216428[5] = 0.0;
   out_388372280572216428[6] = 0.0;
   out_388372280572216428[7] = 0.0;
   out_388372280572216428[8] = 0.0;
   out_388372280572216428[9] = 1.0;
   out_388372280572216428[10] = 0.0;
   out_388372280572216428[11] = 0.0;
   out_388372280572216428[12] = 0.0;
   out_388372280572216428[13] = 0.0;
   out_388372280572216428[14] = 0.0;
   out_388372280572216428[15] = 0.0;
   out_388372280572216428[16] = 0.0;
   out_388372280572216428[17] = 0.0;
   out_388372280572216428[18] = 1.0;
   out_388372280572216428[19] = 0.0;
   out_388372280572216428[20] = 0.0;
   out_388372280572216428[21] = 0.0;
   out_388372280572216428[22] = 0.0;
   out_388372280572216428[23] = 0.0;
   out_388372280572216428[24] = 0.0;
   out_388372280572216428[25] = 0.0;
   out_388372280572216428[26] = 0.0;
   out_388372280572216428[27] = 1.0;
   out_388372280572216428[28] = 0.0;
   out_388372280572216428[29] = 0.0;
   out_388372280572216428[30] = 0.0;
   out_388372280572216428[31] = 0.0;
   out_388372280572216428[32] = 0.0;
   out_388372280572216428[33] = 0.0;
   out_388372280572216428[34] = 0.0;
   out_388372280572216428[35] = 0.0;
   out_388372280572216428[36] = 1.0;
   out_388372280572216428[37] = 0.0;
   out_388372280572216428[38] = 0.0;
   out_388372280572216428[39] = 0.0;
   out_388372280572216428[40] = 0.0;
   out_388372280572216428[41] = 0.0;
   out_388372280572216428[42] = 0.0;
   out_388372280572216428[43] = 0.0;
   out_388372280572216428[44] = 0.0;
   out_388372280572216428[45] = 1.0;
   out_388372280572216428[46] = 0.0;
   out_388372280572216428[47] = 0.0;
   out_388372280572216428[48] = 0.0;
   out_388372280572216428[49] = 0.0;
   out_388372280572216428[50] = 0.0;
   out_388372280572216428[51] = 0.0;
   out_388372280572216428[52] = 0.0;
   out_388372280572216428[53] = 0.0;
   out_388372280572216428[54] = 1.0;
   out_388372280572216428[55] = 0.0;
   out_388372280572216428[56] = 0.0;
   out_388372280572216428[57] = 0.0;
   out_388372280572216428[58] = 0.0;
   out_388372280572216428[59] = 0.0;
   out_388372280572216428[60] = 0.0;
   out_388372280572216428[61] = 0.0;
   out_388372280572216428[62] = 0.0;
   out_388372280572216428[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6950273804263183812) {
   out_6950273804263183812[0] = state[0];
   out_6950273804263183812[1] = state[1];
   out_6950273804263183812[2] = state[2];
   out_6950273804263183812[3] = state[3];
   out_6950273804263183812[4] = state[4];
   out_6950273804263183812[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6950273804263183812[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6950273804263183812[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2419952130103293792) {
   out_2419952130103293792[0] = 1;
   out_2419952130103293792[1] = 0;
   out_2419952130103293792[2] = 0;
   out_2419952130103293792[3] = 0;
   out_2419952130103293792[4] = 0;
   out_2419952130103293792[5] = 0;
   out_2419952130103293792[6] = 0;
   out_2419952130103293792[7] = 0;
   out_2419952130103293792[8] = 0;
   out_2419952130103293792[9] = 1;
   out_2419952130103293792[10] = 0;
   out_2419952130103293792[11] = 0;
   out_2419952130103293792[12] = 0;
   out_2419952130103293792[13] = 0;
   out_2419952130103293792[14] = 0;
   out_2419952130103293792[15] = 0;
   out_2419952130103293792[16] = 0;
   out_2419952130103293792[17] = 0;
   out_2419952130103293792[18] = 1;
   out_2419952130103293792[19] = 0;
   out_2419952130103293792[20] = 0;
   out_2419952130103293792[21] = 0;
   out_2419952130103293792[22] = 0;
   out_2419952130103293792[23] = 0;
   out_2419952130103293792[24] = 0;
   out_2419952130103293792[25] = 0;
   out_2419952130103293792[26] = 0;
   out_2419952130103293792[27] = 1;
   out_2419952130103293792[28] = 0;
   out_2419952130103293792[29] = 0;
   out_2419952130103293792[30] = 0;
   out_2419952130103293792[31] = 0;
   out_2419952130103293792[32] = 0;
   out_2419952130103293792[33] = 0;
   out_2419952130103293792[34] = 0;
   out_2419952130103293792[35] = 0;
   out_2419952130103293792[36] = 1;
   out_2419952130103293792[37] = 0;
   out_2419952130103293792[38] = 0;
   out_2419952130103293792[39] = 0;
   out_2419952130103293792[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2419952130103293792[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2419952130103293792[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2419952130103293792[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2419952130103293792[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2419952130103293792[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2419952130103293792[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2419952130103293792[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2419952130103293792[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2419952130103293792[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2419952130103293792[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2419952130103293792[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2419952130103293792[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2419952130103293792[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2419952130103293792[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2419952130103293792[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2419952130103293792[56] = 0;
   out_2419952130103293792[57] = 0;
   out_2419952130103293792[58] = 0;
   out_2419952130103293792[59] = 0;
   out_2419952130103293792[60] = 0;
   out_2419952130103293792[61] = 0;
   out_2419952130103293792[62] = 0;
   out_2419952130103293792[63] = 1;
}
void h_25(double *state, double *unused, double *out_5363960631214508521) {
   out_5363960631214508521[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2239868391431655836) {
   out_2239868391431655836[0] = 0;
   out_2239868391431655836[1] = 0;
   out_2239868391431655836[2] = 0;
   out_2239868391431655836[3] = 0;
   out_2239868391431655836[4] = 0;
   out_2239868391431655836[5] = 0;
   out_2239868391431655836[6] = 1;
   out_2239868391431655836[7] = 0;
}
void h_24(double *state, double *unused, double *out_2062856328637710990) {
   out_2062856328637710990[0] = state[4];
   out_2062856328637710990[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8930133330717296648) {
   out_8930133330717296648[0] = 0;
   out_8930133330717296648[1] = 0;
   out_8930133330717296648[2] = 0;
   out_8930133330717296648[3] = 0;
   out_8930133330717296648[4] = 1;
   out_8930133330717296648[5] = 0;
   out_8930133330717296648[6] = 0;
   out_8930133330717296648[7] = 0;
   out_8930133330717296648[8] = 0;
   out_8930133330717296648[9] = 0;
   out_8930133330717296648[10] = 0;
   out_8930133330717296648[11] = 0;
   out_8930133330717296648[12] = 0;
   out_8930133330717296648[13] = 1;
   out_8930133330717296648[14] = 0;
   out_8930133330717296648[15] = 0;
}
void h_30(double *state, double *unused, double *out_680766082583151051) {
   out_680766082583151051[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4288408685475712954) {
   out_4288408685475712954[0] = 0;
   out_4288408685475712954[1] = 0;
   out_4288408685475712954[2] = 0;
   out_4288408685475712954[3] = 0;
   out_4288408685475712954[4] = 1;
   out_4288408685475712954[5] = 0;
   out_4288408685475712954[6] = 0;
   out_4288408685475712954[7] = 0;
}
void h_26(double *state, double *unused, double *out_5928987183223659744) {
   out_5928987183223659744[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1633925218784246052) {
   out_1633925218784246052[0] = 0;
   out_1633925218784246052[1] = 0;
   out_1633925218784246052[2] = 0;
   out_1633925218784246052[3] = 0;
   out_1633925218784246052[4] = 0;
   out_1633925218784246052[5] = 0;
   out_1633925218784246052[6] = 0;
   out_1633925218784246052[7] = 1;
}
void h_27(double *state, double *unused, double *out_1644625139496102996) {
   out_1644625139496102996[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8331498373242644936) {
   out_8331498373242644936[0] = 0;
   out_8331498373242644936[1] = 0;
   out_8331498373242644936[2] = 0;
   out_8331498373242644936[3] = 1;
   out_8331498373242644936[4] = 0;
   out_8331498373242644936[5] = 0;
   out_8331498373242644936[6] = 0;
   out_8331498373242644936[7] = 0;
}
void h_29(double *state, double *unused, double *out_1949514412564825142) {
   out_1949514412564825142[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4075708323965073868) {
   out_4075708323965073868[0] = 0;
   out_4075708323965073868[1] = 1;
   out_4075708323965073868[2] = 0;
   out_4075708323965073868[3] = 0;
   out_4075708323965073868[4] = 0;
   out_4075708323965073868[5] = 0;
   out_4075708323965073868[6] = 0;
   out_4075708323965073868[7] = 0;
}
void h_28(double *state, double *unused, double *out_7141644414365002210) {
   out_7141644414365002210[0] = state[5];
   out_7141644414365002210[1] = state[6];
}
void H_28(double *state, double *unused, double *out_1069462091680102462) {
   out_1069462091680102462[0] = 0;
   out_1069462091680102462[1] = 0;
   out_1069462091680102462[2] = 0;
   out_1069462091680102462[3] = 0;
   out_1069462091680102462[4] = 0;
   out_1069462091680102462[5] = 1;
   out_1069462091680102462[6] = 0;
   out_1069462091680102462[7] = 0;
   out_1069462091680102462[8] = 0;
   out_1069462091680102462[9] = 0;
   out_1069462091680102462[10] = 0;
   out_1069462091680102462[11] = 0;
   out_1069462091680102462[12] = 0;
   out_1069462091680102462[13] = 0;
   out_1069462091680102462[14] = 1;
   out_1069462091680102462[15] = 0;
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
