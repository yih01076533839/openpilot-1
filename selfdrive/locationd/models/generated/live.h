/******************************************************************************
 *                      Code generated with sympy 1.6.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6901136647582497383);
void inv_err_fun(double *nom_x, double *true_x, double *out_4696231310121475947);
void H_mod_fun(double *state, double *out_2775029769364007813);
void f_fun(double *state, double dt, double *out_3045977256224477729);
void F_fun(double *state, double dt, double *out_2367992463432873292);
void h_3(double *state, double *unused, double *out_4997303934406486452);
void H_3(double *state, double *unused, double *out_5248447944578785974);
void h_4(double *state, double *unused, double *out_510821180175175656);
void H_4(double *state, double *unused, double *out_5029204696515159069);
void h_9(double *state, double *unused, double *out_7911241815808842404);
void H_9(double *state, double *unused, double *out_1691658230340065029);
void h_10(double *state, double *unused, double *out_2553232778568391880);
void H_10(double *state, double *unused, double *out_5670948878238040113);
void h_12(double *state, double *unused, double *out_8699363372027606959);
void H_12(double *state, double *unused, double *out_1475708198208418628);
void h_31(double *state, double *unused, double *out_1162995685292145910);
void H_31(double *state, double *unused, double *out_1199381576935001193);
void h_32(double *state, double *unused, double *out_8292430729612212799);
void H_32(double *state, double *unused, double *out_4738965817137512777);
void h_13(double *state, double *unused, double *out_7901306498293793494);
void H_13(double *state, double *unused, double *out_2066657572675916836);
void h_14(double *state, double *unused, double *out_7911241815808842404);
void H_14(double *state, double *unused, double *out_1691658230340065029);
void h_19(double *state, double *unused, double *out_765182631658051776);
void H_19(double *state, double *unused, double *out_6124338081442496210);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);