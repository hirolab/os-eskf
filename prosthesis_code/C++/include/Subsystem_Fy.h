//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: Subsystem_Fy.h
//
// Code generated for Simulink model 'Subsystem_Fy'.
//
// Model version                  : 1.198
// Simulink Coder version         : 9.5 (R2021a) 14-Nov-2020
// C/C++ source code generated on : Mon Aug 30 09:00:55 2021
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Windows64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#ifndef RTW_HEADER_Subsystem_Fy_h_
#define RTW_HEADER_Subsystem_Fy_h_
#include <cmath>
#include <cstring>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "Subsystem_Fy_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"


// Macros for accessing real-time model data structure
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

// Class declaration for model Subsystem_Fy
class Subsystem_FyModelClass {
  // public data and function members
 public:
  // Block signals (default storage)
  struct B_Subsystem_Fy_T {
    real_T dv[1600];
    real_T b[1600];
    real_T b_m[1600];
    real_T b_c[1600];
    real_T c[1600];
    real_T G[1600];
    real_T c_y[1600];
    real_T dv1[1600];
    real_T dv2[1600];
    real_T dv3[1600];
    real_T G_k[1600];
    real_T c_c[1200];
    real_T c_b[1200];
    real_T H[520];
    real_T b_H_data[520];
    real_T b_Y_data[520];
    real_T b_B_data[520];
    real_T K_data[520];
    real_T y_data[520];
    int8_T iv[1600];
    real_T measurementCov_data[169];
    real_T b_data[169];
    real_T b_data_p[169];
    real_T b_A_data[169];
    real_T c_A_data[169];
    real_T b_A_data_c[169];
    real_T State[42];
    real_T dv4[42];
    real_T dx[40];
    creal_T eigVec[16];
    creal_T At[16];
    real_T R[16];
    real_T b_I[16];
    real_T a[16];
    real_T ym[13];
    real_T res_data[13];
    real_T yh[13];
    real_T tau_data[13];
    real_T work_data[13];
    real_T vn1_data[13];
    real_T vn2_data[13];
    real_T dv5[12];
    real_T R_f[9];
    real_T dcm[9];
    real_T b_dcm[9];
    real_T dcm_g[9];
    creal_T eigVal[4];
    creal_T beta1[4];
    creal_T work1[4];
    int32_T measurementCov_tmp_data[13];
    int32_T jpvt_data[13];
    int32_T ipiv_data[13];
    real_T State_g[4];
    real_T fquat[4];
    real_T squat[4];
    real_T varargin_1[4];
    real_T b_D[4];
    real_T work[4];
    real_T rworka[4];
    real_T dv6[3];
    real_T dv7[3];
    real_T s1[3];
    real_T b_dcm_m[3];
    real_T v[3];
    creal_T mm1_abs;
    creal_T mm2_abs;
    creal_T s;
    creal_T ctemp;
    creal_T ad22;
    creal_T shift;
    creal_T ascale;
    int32_T rscale[4];
    real_T axang_idx_0;
    real_T axang_idx_1;
    real_T axang_idx_2;
    real_T axang_idx_3;
    real_T q12_idx_0;
    real_T q12_idx_1;
    real_T q12_idx_2;
    real_T d;
    real_T d1;
    real_T bkj;
    real_T tol;
    real_T tol_tmp;
    real_T temp;
    real_T smax;
    real_T y;
    real_T t3;
    real_T t4;
    real_T t5;
    real_T t6;
    real_T t10;
    real_T t11;
    real_T t12;
    real_T t13;
    real_T t14;
    real_T t15;
    real_T t16;
    real_T t17;
    real_T t18;
    real_T t19;
    real_T t20;
    real_T t21;
    real_T t22;
    real_T t23;
    real_T t24;
    real_T t25;
    real_T t43;
    real_T t44;
    real_T t45;
    real_T t49;
    real_T t50;
    real_T t51;
    real_T t67;
    real_T t68;
    real_T t69;
    real_T t70;
    real_T t71;
    real_T t72;
    real_T t73;
    real_T t74;
    real_T t75;
    real_T t76;
    real_T t77;
    real_T t78;
    real_T t133;
    real_T t137;
    real_T t121;
    real_T t122;
    real_T t123;
    real_T t130;
    real_T t131;
    real_T t132;
    real_T t165;
    real_T t166;
    real_T t167;
    real_T t168;
    real_T t169;
    real_T t170;
    real_T t171;
    real_T t172;
    real_T t185;
    real_T t186;
    real_T t141;
    real_T t143;
    real_T t144;
    real_T t145;
    real_T t146;
    real_T t147;
    real_T t149;
    real_T t150;
    real_T t151;
    real_T t152;
    real_T t153;
    real_T t154;
    real_T t155;
    real_T t156;
    real_T t157;
    real_T t158;
    real_T t159;
    real_T t160;
    real_T t161;
    real_T t162;
    real_T t163;
    real_T t164;
    real_T t187;
    real_T t188;
    real_T t189;
    real_T t190;
    real_T t191;
    real_T t192;
    real_T t193;
    real_T t194;
    real_T t195;
    real_T t196;
    real_T t197;
    real_T t198;
    real_T t205;
    real_T t206;
    real_T t207;
    real_T t208;
    real_T t209;
    real_T t210;
    real_T t218;
    real_T t219;
    real_T t220;
    real_T t221;
    real_T t225;
    real_T t226;
    real_T t227;
    real_T t228;
    real_T t229;
    real_T t230;
    real_T t231;
    real_T t232;
    real_T t239;
    real_T t240;
    real_T t241;
    real_T t249;
    real_T t250;
    real_T t254;
    real_T t255;
    real_T t269;
    real_T t270;
    real_T t271;
    real_T t273;
    real_T t274;
    real_T t275;
    real_T t276;
    real_T t279;
    real_T t280;
    real_T t281;
    real_T t282;
    real_T t283;
    real_T t285;
    real_T t286;
    real_T t287;
    real_T t288;
    real_T t291;
    real_T t292;
    real_T t356;
    real_T t357;
    real_T t358;
    real_T t362;
    real_T t363;
    real_T t364;
    real_T t471;
    real_T t472;
    real_T t474;
    real_T t477;
    real_T t478;
    real_T t480;
    real_T t446;
    real_T t447;
    real_T t449;
    real_T t450;
    real_T t578;
    real_T t582;
    real_T t375;
    real_T t376;
    real_T t377;
    real_T t557;
    real_T t558;
    real_T t575;
    real_T t576;
    real_T t577;
    real_T t579;
    real_T t580;
    real_T t581;
    real_T t731;
    real_T t732;
    real_T t799;
    real_T t801;
    real_T t841;
    real_T t842;
    real_T t843;
    real_T t844;
    real_T t845;
    real_T t846;
    real_T t800;
    real_T t802;
    real_T t862;
    real_T t863;
    real_T t864;
    real_T t865;
    real_T t866;
    real_T t867;
    real_T t868;
    real_T t869;
    real_T t870;
    real_T t874;
    real_T t875;
    real_T t876;
    real_T t877;
    real_T t878;
    real_T t879;
    real_T t880;
    real_T t881;
    real_T t882;
    real_T t883;
    real_T t884;
    real_T t885;
    real_T t886;
    real_T t887;
    real_T t888;
    real_T t889;
    real_T t890;
    real_T t921;
    real_T t923;
    real_T t924;
    real_T t925;
    real_T t926;
    real_T t928;
    real_T t929;
    real_T t930;
    real_T t931;
    real_T t933;
    real_T t934;
    real_T t935;
    real_T t937;
    real_T t938;
    real_T t939;
    real_T t941;
    real_T t942;
    real_T t943;
    real_T t946;
    real_T t947;
    real_T t948;
    real_T t950;
    real_T t951;
    real_T t952;
    real_T t953;
    real_T t954;
    real_T t955;
    real_T t956;
    real_T t957;
    real_T t958;
    real_T t647;
    real_T t648;
    real_T t649;
    real_T t650;
    real_T t651;
    real_T t652;
    real_T t653;
    real_T t654;
    real_T t779;
    real_T t780;
    real_T t781;
    real_T t782;
    real_T t783;
    real_T t784;
    real_T t661;
    real_T t662;
    real_T t663;
    real_T t664;
    real_T t665;
    real_T t666;
    real_T t1079;
    real_T t1080;
    real_T t1082;
    real_T t1083;
    real_T t1094;
    real_T t1095;
    real_T t1096;
    real_T t828;
    real_T t959;
    real_T t960;
    real_T t961;
    real_T t962;
    real_T t963;
    real_T t964;
    real_T t965;
    real_T t966;
    real_T t967;
    real_T t968;
    real_T t969;
    real_T t970;
    real_T t835;
    real_T t836;
    real_T t837;
    real_T t838;
    real_T t839;
    real_T t840;
    real_T t1079_tmp;
    real_T t1079_tmp_n;
    real_T t1079_tmp_p;
    real_T t1079_tmp_l;
    real_T t1080_tmp;
    real_T t1080_tmp_j;
    real_T t1080_tmp_d;
    real_T t1080_tmp_g;
    real_T t1083_tmp;
    real_T c_x;
    real_T yaxis_idx_1;
    real_T mag_north_idx_2;
    real_T mag_north_idx_1;
    real_T yaxis_idx_2;
    real_T xaxis_idx_0;
    real_T mag_north_idx_0;
    real_T yaxis_idx_0;
    real_T dcm_tmp;
    real_T dcm_tmp_l;
    real_T dcm_tmp_d;
    real_T dcm_tmp_dy;
    real_T K12;
    real_T K13;
    real_T K14;
    real_T K23;
    real_T K24;
    real_T K34;
    real_T anrm;
    real_T b_absxk;
    real_T cfromc;
    real_T ctoc;
    real_T cto1;
    real_T mul;
    real_T stemp_im;
    real_T alpha1;
    real_T beta1_l;
    real_T xnorm_tmp;
    real_T tau_idx_1;
    real_T tau_idx_0;
    real_T xnorm_tmp_tmp;
    real_T tst;
    real_T htmp1;
    real_T htmp2;
    real_T ba;
    real_T aa;
    real_T h12;
    real_T h21s;
    real_T a__3;
    real_T a__4;
    real_T s_tmp;
    real_T p;
    real_T bcmax;
    real_T bcmis;
    real_T scale;
    real_T z;
    real_T tau;
    real_T a_o;
    real_T t2;
    real_T t3_b;
    real_T t4_n;
    real_T t5_b;
    real_T t6_l;
    real_T t7;
    real_T t8;
    real_T t9;
    real_T t12_h;
    real_T t13_b;
    real_T t14_d;
    real_T t15_e;
    real_T t16_b;
    real_T t17_j;
    real_T t18_f;
    real_T t19_a;
    real_T t20_j;
    real_T t21_j;
    real_T t22_o;
    real_T t23_n;
    real_T t28;
    real_T t30;
    real_T t31;
    real_T t32;
    real_T t33;
    real_T t45_i;
    real_T t58;
    real_T t59;
    real_T t60;
    real_T t61;
    real_T t62;
    real_T t63;
    real_T t64;
    real_T t65;
    real_T t66;
    real_T t67_o;
    real_T t68_n;
    real_T t69_m;
    real_T t70_c;
    real_T t71_m;
    real_T t72_m;
    real_T t73_j;
    real_T t74_h;
    real_T t75_c;
    real_T t82;
    real_T t83;
    real_T t84;
    real_T t88;
    real_T t89;
    real_T t90;
    real_T t94;
    real_T t100;
    real_T t101;
    real_T t76_c;
    real_T t81;
    real_T t98;
    real_T t106;
    real_T t110;
    real_T t111;
    real_T t112;
    real_T t113;
    real_T smax_p;
    real_T beta1_p;
    real_T c_a;
    real_T anorm;
    real_T scale_e;
    real_T ssq;
    real_T colscale;
    real_T absxk;
    real_T t;
    real_T ar;
    real_T ai;
    real_T shift_im;
    real_T eshift_re;
    real_T eshift_im;
    real_T shift_tmp;
    real_T scale_a;
    real_T g2;
    real_T f2s;
    real_T di;
    real_T x;
    real_T fs_re;
    real_T fs_im;
    real_T gs_re;
    real_T gs_im;
    real_T anorm_a;
    real_T ascale_i;
    real_T temp_l;
    real_T acoeff;
    real_T scale_o;
    real_T dmin;
    real_T f_y;
    real_T salpha_re;
    real_T salpha_im;
    real_T qin_idx_0;
    real_T qin_idx_1;
    real_T qin_idx_2;
    real_T qin_idx_3;
    real_T dcm_tmp_o;
    real_T t52;
    real_T t56;
    real_T t46;
    real_T t47;
    real_T t48;
    real_T t49_i;
    real_T t50_f;
    real_T t51_i;
    real_T t60_f;
    real_T t61_g;
    real_T t68_c;
    real_T scale_o3;
    real_T g2_l;
    real_T d_m;
    real_T f2s_m;
    real_T x_c;
    real_T fs_re_f;
    real_T fs_im_p;
    real_T gs_re_e;
    real_T xnorm;
    real_T scale_o4;
    real_T absxk_h;
    real_T t_l;
    real_T temp_h;
    real_T absxr;
    real_T temp_m;
    int32_T measurementCov_size[2];
    int32_T b_H_size[2];
    int32_T K_size[2];
    int32_T b_y_size[2];
    int32_T b_size[2];
    int32_T c_A_size[2];
    int32_T jpvt_size[2];
    int32_T b_A_size[2];
    int32_T ipiv_size[2];
    int32_T i;
    int32_T ch;
    int32_T loop_ub;
    int32_T i_m;
    int32_T measurementCov_tmp_size_idx_0;
    int32_T m;
    int32_T coffset;
    int32_T boffset;
    int32_T aoffset;
    int32_T i_h;
    int32_T b_i;
    int32_T b_m_c;
    int32_T c_aoffset;
    int32_T y_size_idx_0;
    int32_T rankA;
    int32_T minmn;
    int32_T maxmn;
    int32_T c_A;
    int32_T tau_size;
    int32_T b_A_size_idx_0;
    int32_T b_B_size_idx_0;
    int32_T b_Y_data_tmp;
    int32_T b_Y_data_tmp_tmp;
    int32_T b_Y_data_tmp_k;
    int32_T b_Y_data_tmp_p;
    int32_T n;
    int32_T jp;
    int32_T b_jAcol;
    int32_T b_jBcol;
    int32_T b_kBcol;
    int32_T i1;
    int32_T idx;
    int32_T k;
    int32_T j;
    int32_T jrow;
    int32_T jcol;
    int32_T c_j;
    int32_T d_i;
    int32_T stemp_re_tmp;
    int32_T stemp_re_tmp_p;
    int32_T knt;
    int32_T lastc;
    int32_T ix;
    int32_T iac;
    int32_T l;
    int32_T ia;
    int32_T jy;
    int32_T minmana;
    int32_T loop_ub_p;
    int32_T ma;
    int32_T minmn_a;
    int32_T ii;
    int32_T b_ix;
  };

  // Block states (default storage) for system '<Root>'
  struct DW_Subsystem_Fy_T {
    emxArray_real_T_1x13_Subsyste_T rowsf1;// '<S1>/MATLAB Function'
    real_T NominalState[42];           // '<S1>/MATLAB Function'
    real_T errorCov_[1600];            // '<S1>/MATLAB Function'
    boolean_T NominalState_not_empty;  // '<S1>/MATLAB Function'
    boolean_T rowsf1_not_empty;        // '<S1>/MATLAB Function'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU_Subsystem_Fy_T {
    real_T simu_a[3];                  // '<Root>/simu_a'
    real_T Fy;                         // '<Root>/Fy'
    real_T Fy_tol;                     // '<Root>/Fy_tol'
    real_T simu_w[3];                  // '<Root>/simu_w'
    real_T simu_m[3];                  // '<Root>/simu_m'
    real_T fimu_a[3];                  // '<Root>/fimu_a'
    real_T fimu_w[3];                  // '<Root>/fimu_w'
    real_T fimu_m[3];                  // '<Root>/fimu_m'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_Subsystem_Fy_T {
    real_T aavel[3];                   // '<Root>/aavel'
    real_T IE_DP_rad[2];               // '<Root>/IE_DP_rad'
    real_T IE_DP_deg[2];               // '<Root>/IE_DP_deg'
    real_T aavel_IE_DP_rad[2];         // '<Root>/aavel_IE_DP_rad'
    real_T aavel_IE_DP_deg[2];         // '<Root>/aavel_IE_DP_deg'
    real_T aangleh[3];                 // '<Root>/aangleh'
    real_T states[42];                 // '<Root>/states'
    real_T meas[19];                   // '<Root>/meas'
  };

  // Real-time Model Data Structure
  struct RT_MODEL_Subsystem_Fy_T {
    const char_T * volatile errorStatus;
  };

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  void terminate();

  // Constructor
  Subsystem_FyModelClass();

  // Destructor
  ~Subsystem_FyModelClass();

  // Root-level structure-based inputs set method

  // Root inports set method
  void setExternalInputs(const ExtU_Subsystem_Fy_T* pExtU_Subsystem_Fy_T)
  {
    Subsystem_Fy_U = *pExtU_Subsystem_Fy_T;
  }

  // Root-level structure-based outputs get method

  // Root outports get method
  const Subsystem_FyModelClass::ExtY_Subsystem_Fy_T & getExternalOutputs() const
  {
    return Subsystem_Fy_Y;
  }

  // Real-Time Model get method
  Subsystem_FyModelClass::RT_MODEL_Subsystem_Fy_T * getRTM();

  // private data and function members
 private:
  // Block signals
  B_Subsystem_Fy_T Subsystem_Fy_B;

  // Block states
  DW_Subsystem_Fy_T Subsystem_Fy_DW;

  // External inputs
  ExtU_Subsystem_Fy_T Subsystem_Fy_U;

  // External outputs
  ExtY_Subsystem_Fy_T Subsystem_Fy_Y;

  // Real-Time Model
  RT_MODEL_Subsystem_Fy_T Subsystem_Fy_M;

  // private member function(s) for subsystem '<Root>'
  boolean_T Subsystem_Fy_anyNonFinite(const real_T x[16]);
  real_T Subsystem_Fy_rt_hypotd_snf(real_T u0, real_T u1);
  void Subsystem_Fy_xzggbal(creal_T A[16], int32_T *ilo, int32_T *ihi, int32_T
    rscale[4]);
  void Subsystem_Fy_sqrt(creal_T *x);
  void Subsystem_Fy_xzlartg_k(const creal_T f, const creal_T g, real_T *cs,
    creal_T *sn);
  void Subsystem_Fy_xzlartg(const creal_T f, const creal_T g, real_T *cs,
    creal_T *sn, creal_T *r);
  void Subsystem_Fy_xzhgeqz(creal_T A[16], int32_T ilo, int32_T ihi, creal_T Z
    [16], int32_T *info, creal_T alpha1[4], creal_T beta1[4]);
  void Subsystem_Fy_xztgevc(const creal_T A[16], creal_T V[16]);
  real_T Subsystem_Fy_xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  void Subsystem_Fy_xzlarf(int32_T m, int32_T n, int32_T iv0, real_T tau, real_T
    C[16], int32_T ic0, real_T work[4]);
  real_T Subsystem_Fy_xnrm2_d(int32_T n, const real_T x[3]);
  real_T Subsystem_Fy_xzlarfg(int32_T n, real_T *alpha1, real_T x[3]);
  void Subsystem_Fy_xdlanv2(real_T *a, real_T *b, real_T *c, real_T *d, real_T
    *rt1r, real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *cs, real_T *sn);
  void Subsystem_Fy_xrot(int32_T n, real_T x[16], int32_T ix0, int32_T iy0,
    real_T c, real_T s);
  void Subsystem_Fy_xrot_m(real_T x[16], int32_T ix0, int32_T iy0, real_T c,
    real_T s);
  int32_T Subsystem_Fy_xhseqr(real_T h[16], real_T z[16]);
  void Subsystem__eigHermitianStandard(const real_T A[16], real_T V[16], real_T
    D[4]);
  void Subsystem_Fy_eig(const real_T A[16], creal_T V[16], creal_T D[4]);
  void Subsystem_Fy_rotm2quat(const real_T R[9], real_T quat[4]);
  void project_state_to_constrain_fcn(const real_T rf1[3], const real_T rfa[3],
    const real_T MAG_NORTH[3], const real_T rsa[3], const real_T rs2[3], const
    real_T s1[3], const real_T am1[3], const real_T mm1[3], const real_T am2[3],
    const real_T mm2[3], real_T xproj[42]);
  void Subsystem__estimate_measurement(const real_T x_nominal_prev[42], real_T
    yh[13], real_T H[520]);
  void Subsystem_Fy_xswap(int32_T n, real_T x_data[], int32_T ix0, int32_T iy0);
  real_T Subsystem_Fy_xnrm2_dg(int32_T n, const real_T x_data[], int32_T ix0);
  void Subsystem_Fy_xzlarf_f(int32_T m, int32_T n, int32_T iv0, real_T tau,
    real_T C_data[], int32_T ic0, int32_T ldc, real_T work_data[]);
  void Subsystem_Fy_qrpf(real_T A_data[], const int32_T A_size[2], int32_T m,
    int32_T n, real_T tau_data[], int32_T jpvt_data[]);
  void Subsystem_Fy_xzgeqp3(real_T A_data[], const int32_T A_size[2], int32_T m,
    int32_T n, real_T tau_data[], int32_T *tau_size, int32_T jpvt_data[],
    int32_T jpvt_size[2]);
  void Subsystem_Fy_xgetrf(int32_T m, int32_T n, real_T A_data[], const int32_T
    A_size[2], int32_T lda, int32_T ipiv_data[], int32_T ipiv_size[2], int32_T
    *info);
  void Subsystem_Fy_lusolve(const real_T A_data[], const int32_T A_size[2],
    const real_T B_data[], const int32_T B_size[2], real_T X_data[], int32_T
    X_size[2]);
  void Subsystem_Fy_mrdiv(const real_T A_data[], const int32_T A_size[2], const
    real_T B_data[], const int32_T B_size[2], real_T Y_data[], int32_T Y_size[2]);
  void Subsyst_eskf_cdprosthesis_reset(const real_T in1[42], const real_T in2[40],
    real_T Xt[42], real_T reset_dx[1600]);
  void Subsystem_Fy_mtimes(const real_T A_data[], const int32_T A_size[2], const
    real_T B_data[], const int32_T B_size[2], real_T C[1600]);
  void Subsystem_Fy_correct(const real_T get_NominalState[42], real_T errorCov_
    [1600], const real_T res_data[], const real_T measurementCov_data[], const
    int32_T measurementCov_size[2], const real_T H_data[], const int32_T H_size
    [2], real_T x_post[42]);
  void Subsystem_F_correct_measurement(const real_T errorCov_c[1600], const
    real_T measurementCov_[169], const real_T mm1[3], const real_T mm2[3],
    boolean_T is_stance, const real_T get_NominalState[42], real_T x_post[42],
    real_T errorCov_[1600]);
  void Subsystem_Fy_predict_eksf(const real_T processNoise[900], const real_T U
    [12], const real_T x_nominal_prev[42], const real_T errorstate_prev[1600],
    real_T NominalState_predict[42], real_T ErrorState_predict[1600]);
  void Subsystem_Fy_quatrotate_edit(const real_T q[4], const real_T r[3], real_T
    qout[3]);
};

//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Note that this particular code originates from a subsystem build,
//  and has its own system numbers different from the parent model.
//  Refer to the system hierarchy for this subsystem below, and use the
//  MATLAB hilite_system command to trace the generated code back
//  to the parent model.  For example,
//
//  hilite_system('eksf/Subsystem_Fy')    - opens subsystem eksf/Subsystem_Fy
//  hilite_system('eksf/Subsystem_Fy/Kp') - opens and selects block Kp
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'eksf'
//  '<S1>'   : 'eksf/Subsystem_Fy'
//  '<S2>'   : 'eksf/Subsystem_Fy/MATLAB Function'
//  '<S3>'   : 'eksf/Subsystem_Fy/MATLAB Function1'

#endif                                 // RTW_HEADER_Subsystem_Fy_h_

//
// File trailer for generated code.
//
// [EOF]
//
