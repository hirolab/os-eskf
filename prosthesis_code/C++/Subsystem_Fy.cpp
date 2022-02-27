//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: Subsystem_Fy.cpp
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
#include "Subsystem_Fy.h"
#include "Subsystem_Fy_private.h"

#include "Subsystem_Fy.h"
#include "Subsystem_Fy_private.h"
#include "rtwtypes.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"
#include "rtGetNaN.h"
#include "Subsystem_Fy_types.h"
#include "rt_nonfinite.cpp"
#include "rtGetInf.cpp"
#include "rtGetNaN.cpp"



// Function for MATLAB Function: '<S1>/MATLAB Function'
boolean_T Subsystem_FyModelClass::Subsystem_Fy_anyNonFinite(const real_T x[16])
{
  real_T x_0;
  int32_T k;
  boolean_T b_p;
  b_p = true;
  for (k = 0; k < 16; k++) {
    x_0 = x[k];
    if (b_p && (rtIsInf(x_0) || rtIsNaN(x_0))) {
      b_p = false;
    }
  }

  return !b_p;
}

real_T Subsystem_FyModelClass::Subsystem_Fy_rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  Subsystem_Fy_B.a_o = std::abs(u0);
  y = std::abs(u1);
  if (Subsystem_Fy_B.a_o < y) {
    Subsystem_Fy_B.a_o /= y;
    y *= std::sqrt(Subsystem_Fy_B.a_o * Subsystem_Fy_B.a_o + 1.0);
  } else if (Subsystem_Fy_B.a_o > y) {
    y /= Subsystem_Fy_B.a_o;
    y = std::sqrt(y * y + 1.0) * Subsystem_Fy_B.a_o;
  } else if (!rtIsNaN(y)) {
    y = Subsystem_Fy_B.a_o * 1.4142135623730951;
  }

  return y;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xzggbal(creal_T A[16], int32_T *ilo,
  int32_T *ihi, int32_T rscale[4])
{
  real_T atmp_im;
  real_T atmp_re;
  int32_T atmp_re_tmp_tmp;
  int32_T exitg1;
  int32_T exitg2;
  int32_T i;
  int32_T ii;
  int32_T j;
  int32_T jj;
  int32_T nzcount;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T found;
  rscale[0] = 1;
  rscale[1] = 1;
  rscale[2] = 1;
  rscale[3] = 1;
  *ilo = 1;
  *ihi = 4;
  do {
    exitg2 = 0;
    i = 0;
    j = 0;
    found = false;
    ii = *ihi;
    exitg3 = false;
    while ((!exitg3) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = *ihi;
      jj = 0;
      exitg4 = false;
      while ((!exitg4) && (jj <= *ihi - 1)) {
        atmp_re_tmp_tmp = ((jj << 2) + ii) - 1;
        if ((A[atmp_re_tmp_tmp].re != 0.0) || (A[atmp_re_tmp_tmp].im != 0.0) ||
            (jj + 1 == ii)) {
          if (nzcount == 0) {
            j = jj + 1;
            nzcount = 1;
            jj++;
          } else {
            nzcount = 2;
            exitg4 = true;
          }
        } else {
          jj++;
        }
      }

      if (nzcount < 2) {
        found = true;
        exitg3 = true;
      } else {
        ii--;
      }
    }

    if (!found) {
      exitg2 = 2;
    } else {
      if (i != *ihi) {
        atmp_re = A[i - 1].re;
        atmp_im = A[i - 1].im;
        A[i - 1] = A[*ihi - 1];
        A[*ihi - 1].re = atmp_re;
        A[*ihi - 1].im = atmp_im;
        atmp_re = A[i + 3].re;
        atmp_im = A[i + 3].im;
        A[i + 3] = A[*ihi + 3];
        A[*ihi + 3].re = atmp_re;
        A[*ihi + 3].im = atmp_im;
        atmp_re = A[i + 7].re;
        atmp_im = A[i + 7].im;
        A[i + 7] = A[*ihi + 7];
        A[*ihi + 7].re = atmp_re;
        A[*ihi + 7].im = atmp_im;
        atmp_re = A[i + 11].re;
        atmp_im = A[i + 11].im;
        A[i + 11] = A[*ihi + 11];
        A[*ihi + 11].re = atmp_re;
        A[*ihi + 11].im = atmp_im;
      }

      if (j != *ihi) {
        for (ii = 0; ii < *ihi; ii++) {
          i = ((j - 1) << 2) + ii;
          atmp_re = A[i].re;
          atmp_im = A[i].im;
          atmp_re_tmp_tmp = ((*ihi - 1) << 2) + ii;
          A[i] = A[atmp_re_tmp_tmp];
          A[atmp_re_tmp_tmp].re = atmp_re;
          A[atmp_re_tmp_tmp].im = atmp_im;
        }
      }

      rscale[*ihi - 1] = j;
      (*ihi)--;
      if (*ihi == 1) {
        rscale[0] = 1;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0);

  if (exitg2 == 1) {
  } else {
    do {
      exitg1 = 0;
      ii = 0;
      j = 0;
      found = false;
      i = *ilo;
      exitg3 = false;
      while ((!exitg3) && (i <= *ihi)) {
        nzcount = 0;
        ii = *ihi;
        j = i;
        jj = *ilo;
        exitg4 = false;
        while ((!exitg4) && (jj <= *ihi)) {
          atmp_re_tmp_tmp = (((i - 1) << 2) + jj) - 1;
          if ((A[atmp_re_tmp_tmp].re != 0.0) || (A[atmp_re_tmp_tmp].im != 0.0) ||
              (jj == i)) {
            if (nzcount == 0) {
              ii = jj;
              nzcount = 1;
              jj++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            jj++;
          }
        }

        if (nzcount < 2) {
          found = true;
          exitg3 = true;
        } else {
          i++;
        }
      }

      if (!found) {
        exitg1 = 1;
      } else {
        if (ii != *ilo) {
          for (nzcount = *ilo - 1; nzcount + 1 < 5; nzcount++) {
            atmp_re_tmp_tmp = nzcount << 2;
            i = (atmp_re_tmp_tmp + ii) - 1;
            atmp_re = A[i].re;
            atmp_im = A[i].im;
            atmp_re_tmp_tmp = (atmp_re_tmp_tmp + *ilo) - 1;
            A[i] = A[atmp_re_tmp_tmp];
            A[atmp_re_tmp_tmp].re = atmp_re;
            A[atmp_re_tmp_tmp].im = atmp_im;
          }
        }

        if (j != *ilo) {
          for (ii = 0; ii < *ihi; ii++) {
            i = ((j - 1) << 2) + ii;
            atmp_re = A[i].re;
            atmp_im = A[i].im;
            atmp_re_tmp_tmp = ((*ilo - 1) << 2) + ii;
            A[i] = A[atmp_re_tmp_tmp];
            A[atmp_re_tmp_tmp].re = atmp_re;
            A[atmp_re_tmp_tmp].im = atmp_im;
          }
        }

        rscale[*ilo - 1] = j;
        (*ilo)++;
        if (*ilo == *ihi) {
          rscale[*ilo - 1] = *ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_sqrt(creal_T *x)
{
  real_T absxi;
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      Subsystem_Fy_B.absxr = 0.0;
      absxi = std::sqrt(-x->re);
    } else {
      Subsystem_Fy_B.absxr = std::sqrt(x->re);
      absxi = 0.0;
    }
  } else if (x->re == 0.0) {
    if (x->im < 0.0) {
      Subsystem_Fy_B.absxr = std::sqrt(-x->im / 2.0);
      absxi = -Subsystem_Fy_B.absxr;
    } else {
      Subsystem_Fy_B.absxr = std::sqrt(x->im / 2.0);
      absxi = Subsystem_Fy_B.absxr;
    }
  } else if (rtIsNaN(x->re)) {
    Subsystem_Fy_B.absxr = x->re;
    absxi = x->re;
  } else if (rtIsNaN(x->im)) {
    Subsystem_Fy_B.absxr = x->im;
    absxi = x->im;
  } else if (rtIsInf(x->im)) {
    Subsystem_Fy_B.absxr = std::abs(x->im);
    absxi = x->im;
  } else if (rtIsInf(x->re)) {
    if (x->re < 0.0) {
      Subsystem_Fy_B.absxr = 0.0;
      absxi = x->im * -x->re;
    } else {
      Subsystem_Fy_B.absxr = x->re;
      absxi = 0.0;
    }
  } else {
    Subsystem_Fy_B.absxr = std::abs(x->re);
    absxi = std::abs(x->im);
    if ((Subsystem_Fy_B.absxr > 4.4942328371557893E+307) || (absxi >
         4.4942328371557893E+307)) {
      Subsystem_Fy_B.absxr *= 0.5;
      absxi = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.absxr, absxi * 0.5);
      if (absxi > Subsystem_Fy_B.absxr) {
        Subsystem_Fy_B.absxr = std::sqrt(Subsystem_Fy_B.absxr / absxi + 1.0) *
          std::sqrt(absxi);
      } else {
        Subsystem_Fy_B.absxr = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      Subsystem_Fy_B.absxr = std::sqrt((Subsystem_Fy_rt_hypotd_snf
        (Subsystem_Fy_B.absxr, absxi) + Subsystem_Fy_B.absxr) * 0.5);
    }

    if (x->re > 0.0) {
      absxi = x->im / Subsystem_Fy_B.absxr * 0.5;
    } else {
      if (x->im < 0.0) {
        absxi = -Subsystem_Fy_B.absxr;
      } else {
        absxi = Subsystem_Fy_B.absxr;
      }

      Subsystem_Fy_B.absxr = x->im / absxi * 0.5;
    }
  }

  x->re = Subsystem_Fy_B.absxr;
  x->im = absxi;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xzlartg_k(const creal_T f, const
  creal_T g, real_T *cs, creal_T *sn)
{
  real_T gs_im;
  boolean_T guard1 = false;
  Subsystem_Fy_B.d_m = std::abs(f.re);
  Subsystem_Fy_B.f2s_m = std::abs(f.im);
  Subsystem_Fy_B.scale_o3 = Subsystem_Fy_B.d_m;
  if (Subsystem_Fy_B.f2s_m > Subsystem_Fy_B.d_m) {
    Subsystem_Fy_B.scale_o3 = Subsystem_Fy_B.f2s_m;
  }

  Subsystem_Fy_B.gs_re_e = std::abs(g.re);
  gs_im = std::abs(g.im);
  if (gs_im > Subsystem_Fy_B.gs_re_e) {
    Subsystem_Fy_B.gs_re_e = gs_im;
  }

  if (Subsystem_Fy_B.gs_re_e > Subsystem_Fy_B.scale_o3) {
    Subsystem_Fy_B.scale_o3 = Subsystem_Fy_B.gs_re_e;
  }

  Subsystem_Fy_B.fs_re_f = f.re;
  Subsystem_Fy_B.fs_im_p = f.im;
  Subsystem_Fy_B.gs_re_e = g.re;
  gs_im = g.im;
  guard1 = false;
  if (Subsystem_Fy_B.scale_o3 >= 7.4428285367870146E+137) {
    do {
      Subsystem_Fy_B.fs_re_f *= 1.3435752215134178E-138;
      Subsystem_Fy_B.fs_im_p *= 1.3435752215134178E-138;
      Subsystem_Fy_B.gs_re_e *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      Subsystem_Fy_B.scale_o3 *= 1.3435752215134178E-138;
    } while (!(Subsystem_Fy_B.scale_o3 < 7.4428285367870146E+137));

    guard1 = true;
  } else if (Subsystem_Fy_B.scale_o3 <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        Subsystem_Fy_B.fs_re_f *= 7.4428285367870146E+137;
        Subsystem_Fy_B.fs_im_p *= 7.4428285367870146E+137;
        Subsystem_Fy_B.gs_re_e *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        Subsystem_Fy_B.scale_o3 *= 7.4428285367870146E+137;
      } while (!(Subsystem_Fy_B.scale_o3 > 1.3435752215134178E-138));

      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    Subsystem_Fy_B.scale_o3 = Subsystem_Fy_B.fs_re_f * Subsystem_Fy_B.fs_re_f +
      Subsystem_Fy_B.fs_im_p * Subsystem_Fy_B.fs_im_p;
    Subsystem_Fy_B.g2_l = Subsystem_Fy_B.gs_re_e * Subsystem_Fy_B.gs_re_e +
      gs_im * gs_im;
    Subsystem_Fy_B.x_c = Subsystem_Fy_B.g2_l;
    if (1.0 > Subsystem_Fy_B.g2_l) {
      Subsystem_Fy_B.x_c = 1.0;
    }

    if (Subsystem_Fy_B.scale_o3 <= Subsystem_Fy_B.x_c * 2.0041683600089728E-292)
    {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        Subsystem_Fy_B.d_m = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.gs_re_e,
          gs_im);
        sn->re = Subsystem_Fy_B.gs_re_e / Subsystem_Fy_B.d_m;
        sn->im = -gs_im / Subsystem_Fy_B.d_m;
      } else {
        Subsystem_Fy_B.scale_o3 = std::sqrt(Subsystem_Fy_B.g2_l);
        *cs = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.fs_re_f,
          Subsystem_Fy_B.fs_im_p) / Subsystem_Fy_B.scale_o3;
        if (Subsystem_Fy_B.f2s_m > Subsystem_Fy_B.d_m) {
          Subsystem_Fy_B.d_m = Subsystem_Fy_B.f2s_m;
        }

        if (Subsystem_Fy_B.d_m > 1.0) {
          Subsystem_Fy_B.d_m = Subsystem_Fy_rt_hypotd_snf(f.re, f.im);
          Subsystem_Fy_B.fs_re_f = f.re / Subsystem_Fy_B.d_m;
          Subsystem_Fy_B.fs_im_p = f.im / Subsystem_Fy_B.d_m;
        } else {
          Subsystem_Fy_B.fs_re_f = 7.4428285367870146E+137 * f.re;
          Subsystem_Fy_B.f2s_m = 7.4428285367870146E+137 * f.im;
          Subsystem_Fy_B.d_m = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.fs_re_f,
            Subsystem_Fy_B.f2s_m);
          Subsystem_Fy_B.fs_re_f /= Subsystem_Fy_B.d_m;
          Subsystem_Fy_B.fs_im_p = Subsystem_Fy_B.f2s_m / Subsystem_Fy_B.d_m;
        }

        Subsystem_Fy_B.gs_re_e /= Subsystem_Fy_B.scale_o3;
        gs_im = -gs_im / Subsystem_Fy_B.scale_o3;
        sn->re = Subsystem_Fy_B.fs_re_f * Subsystem_Fy_B.gs_re_e -
          Subsystem_Fy_B.fs_im_p * gs_im;
        sn->im = Subsystem_Fy_B.fs_re_f * gs_im + Subsystem_Fy_B.fs_im_p *
          Subsystem_Fy_B.gs_re_e;
      }
    } else {
      Subsystem_Fy_B.f2s_m = std::sqrt(Subsystem_Fy_B.g2_l /
        Subsystem_Fy_B.scale_o3 + 1.0);
      *cs = 1.0 / Subsystem_Fy_B.f2s_m;
      Subsystem_Fy_B.d_m = Subsystem_Fy_B.scale_o3 + Subsystem_Fy_B.g2_l;
      Subsystem_Fy_B.fs_re_f = Subsystem_Fy_B.f2s_m * Subsystem_Fy_B.fs_re_f /
        Subsystem_Fy_B.d_m;
      Subsystem_Fy_B.fs_im_p = Subsystem_Fy_B.f2s_m * Subsystem_Fy_B.fs_im_p /
        Subsystem_Fy_B.d_m;
      sn->re = Subsystem_Fy_B.fs_re_f * Subsystem_Fy_B.gs_re_e -
        Subsystem_Fy_B.fs_im_p * -gs_im;
      sn->im = Subsystem_Fy_B.fs_re_f * -gs_im + Subsystem_Fy_B.fs_im_p *
        Subsystem_Fy_B.gs_re_e;
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xzlartg(const creal_T f, const creal_T
  g, real_T *cs, creal_T *sn, creal_T *r)
{
  int32_T count;
  int32_T rescaledir;
  boolean_T guard1 = false;
  Subsystem_Fy_B.f2s = std::abs(f.re);
  Subsystem_Fy_B.di = std::abs(f.im);
  Subsystem_Fy_B.scale_a = Subsystem_Fy_B.f2s;
  if (Subsystem_Fy_B.di > Subsystem_Fy_B.f2s) {
    Subsystem_Fy_B.scale_a = Subsystem_Fy_B.di;
  }

  Subsystem_Fy_B.gs_re = std::abs(g.re);
  Subsystem_Fy_B.gs_im = std::abs(g.im);
  if (Subsystem_Fy_B.gs_im > Subsystem_Fy_B.gs_re) {
    Subsystem_Fy_B.gs_re = Subsystem_Fy_B.gs_im;
  }

  if (Subsystem_Fy_B.gs_re > Subsystem_Fy_B.scale_a) {
    Subsystem_Fy_B.scale_a = Subsystem_Fy_B.gs_re;
  }

  Subsystem_Fy_B.fs_re = f.re;
  Subsystem_Fy_B.fs_im = f.im;
  Subsystem_Fy_B.gs_re = g.re;
  Subsystem_Fy_B.gs_im = g.im;
  count = -1;
  rescaledir = 0;
  guard1 = false;
  if (Subsystem_Fy_B.scale_a >= 7.4428285367870146E+137) {
    do {
      count++;
      Subsystem_Fy_B.fs_re *= 1.3435752215134178E-138;
      Subsystem_Fy_B.fs_im *= 1.3435752215134178E-138;
      Subsystem_Fy_B.gs_re *= 1.3435752215134178E-138;
      Subsystem_Fy_B.gs_im *= 1.3435752215134178E-138;
      Subsystem_Fy_B.scale_a *= 1.3435752215134178E-138;
    } while (!(Subsystem_Fy_B.scale_a < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = true;
  } else if (Subsystem_Fy_B.scale_a <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        Subsystem_Fy_B.fs_re *= 7.4428285367870146E+137;
        Subsystem_Fy_B.fs_im *= 7.4428285367870146E+137;
        Subsystem_Fy_B.gs_re *= 7.4428285367870146E+137;
        Subsystem_Fy_B.gs_im *= 7.4428285367870146E+137;
        Subsystem_Fy_B.scale_a *= 7.4428285367870146E+137;
      } while (!(Subsystem_Fy_B.scale_a > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    Subsystem_Fy_B.scale_a = Subsystem_Fy_B.fs_re * Subsystem_Fy_B.fs_re +
      Subsystem_Fy_B.fs_im * Subsystem_Fy_B.fs_im;
    Subsystem_Fy_B.g2 = Subsystem_Fy_B.gs_re * Subsystem_Fy_B.gs_re +
      Subsystem_Fy_B.gs_im * Subsystem_Fy_B.gs_im;
    Subsystem_Fy_B.x = Subsystem_Fy_B.g2;
    if (1.0 > Subsystem_Fy_B.g2) {
      Subsystem_Fy_B.x = 1.0;
    }

    if (Subsystem_Fy_B.scale_a <= Subsystem_Fy_B.x * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = Subsystem_Fy_rt_hypotd_snf(g.re, g.im);
        r->im = 0.0;
        Subsystem_Fy_B.f2s = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.gs_re,
          Subsystem_Fy_B.gs_im);
        sn->re = Subsystem_Fy_B.gs_re / Subsystem_Fy_B.f2s;
        sn->im = -Subsystem_Fy_B.gs_im / Subsystem_Fy_B.f2s;
      } else {
        Subsystem_Fy_B.scale_a = std::sqrt(Subsystem_Fy_B.g2);
        *cs = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.fs_re,
          Subsystem_Fy_B.fs_im) / Subsystem_Fy_B.scale_a;
        if (Subsystem_Fy_B.di > Subsystem_Fy_B.f2s) {
          Subsystem_Fy_B.f2s = Subsystem_Fy_B.di;
        }

        if (Subsystem_Fy_B.f2s > 1.0) {
          Subsystem_Fy_B.f2s = Subsystem_Fy_rt_hypotd_snf(f.re, f.im);
          Subsystem_Fy_B.fs_re = f.re / Subsystem_Fy_B.f2s;
          Subsystem_Fy_B.fs_im = f.im / Subsystem_Fy_B.f2s;
        } else {
          Subsystem_Fy_B.fs_re = 7.4428285367870146E+137 * f.re;
          Subsystem_Fy_B.di = 7.4428285367870146E+137 * f.im;
          Subsystem_Fy_B.f2s = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.fs_re,
            Subsystem_Fy_B.di);
          Subsystem_Fy_B.fs_re /= Subsystem_Fy_B.f2s;
          Subsystem_Fy_B.fs_im = Subsystem_Fy_B.di / Subsystem_Fy_B.f2s;
        }

        Subsystem_Fy_B.gs_re /= Subsystem_Fy_B.scale_a;
        Subsystem_Fy_B.gs_im = -Subsystem_Fy_B.gs_im / Subsystem_Fy_B.scale_a;
        sn->re = Subsystem_Fy_B.fs_re * Subsystem_Fy_B.gs_re -
          Subsystem_Fy_B.fs_im * Subsystem_Fy_B.gs_im;
        sn->im = Subsystem_Fy_B.fs_re * Subsystem_Fy_B.gs_im +
          Subsystem_Fy_B.fs_im * Subsystem_Fy_B.gs_re;
        r->re = (sn->re * g.re - sn->im * g.im) + *cs * f.re;
        r->im = (sn->re * g.im + sn->im * g.re) + *cs * f.im;
      }
    } else {
      Subsystem_Fy_B.f2s = std::sqrt(Subsystem_Fy_B.g2 / Subsystem_Fy_B.scale_a
        + 1.0);
      r->re = Subsystem_Fy_B.f2s * Subsystem_Fy_B.fs_re;
      r->im = Subsystem_Fy_B.f2s * Subsystem_Fy_B.fs_im;
      *cs = 1.0 / Subsystem_Fy_B.f2s;
      Subsystem_Fy_B.f2s = Subsystem_Fy_B.scale_a + Subsystem_Fy_B.g2;
      Subsystem_Fy_B.fs_re = r->re / Subsystem_Fy_B.f2s;
      Subsystem_Fy_B.f2s = r->im / Subsystem_Fy_B.f2s;
      sn->re = Subsystem_Fy_B.fs_re * Subsystem_Fy_B.gs_re - Subsystem_Fy_B.f2s *
        -Subsystem_Fy_B.gs_im;
      sn->im = Subsystem_Fy_B.fs_re * -Subsystem_Fy_B.gs_im + Subsystem_Fy_B.f2s
        * Subsystem_Fy_B.gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else if (rescaledir < 0) {
        for (rescaledir = 0; rescaledir <= count; rescaledir++) {
          r->re *= 1.3435752215134178E-138;
          r->im *= 1.3435752215134178E-138;
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xzhgeqz(creal_T A[16], int32_T ilo,
  int32_T ihi, creal_T Z[16], int32_T *info, creal_T alpha1[4], creal_T beta1[4])
{
  int32_T absxk_tmp;
  int32_T col;
  int32_T ctemp_tmp;
  int32_T ctemp_tmp_tmp;
  int32_T exitg1;
  int32_T ifirst;
  int32_T iiter;
  int32_T ilastm1;
  int32_T j;
  int32_T jp1;
  int32_T nm1;
  int32_T shift_tmp;
  boolean_T exitg2;
  boolean_T failed;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  *info = 0;
  alpha1[0].re = 0.0;
  alpha1[0].im = 0.0;
  beta1[0].re = 1.0;
  beta1[0].im = 0.0;
  alpha1[1].re = 0.0;
  alpha1[1].im = 0.0;
  beta1[1].re = 1.0;
  beta1[1].im = 0.0;
  alpha1[2].re = 0.0;
  alpha1[2].im = 0.0;
  beta1[2].re = 1.0;
  beta1[2].im = 0.0;
  alpha1[3].re = 0.0;
  alpha1[3].im = 0.0;
  beta1[3].re = 1.0;
  beta1[3].im = 0.0;
  Subsystem_Fy_B.eshift_re = 0.0;
  Subsystem_Fy_B.eshift_im = 0.0;
  Subsystem_Fy_B.ctemp.re = 0.0;
  Subsystem_Fy_B.ctemp.im = 0.0;
  Subsystem_Fy_B.anorm = 0.0;
  if (ilo <= ihi) {
    Subsystem_Fy_B.scale_e = 3.3121686421112381E-170;
    Subsystem_Fy_B.ssq = 0.0;
    nm1 = ihi - ilo;
    for (ifirst = -1; ifirst < nm1; ifirst++) {
      Subsystem_Fy_B.colscale = 3.3121686421112381E-170;
      Subsystem_Fy_B.anorm = 0.0;
      col = ilo + ifirst;
      if (ifirst + 2 < nm1) {
        ilastm1 = ifirst + 2;
      } else {
        ilastm1 = nm1;
      }

      ilastm1 += ilo;
      for (iiter = ilo; iiter <= ilastm1; iiter++) {
        absxk_tmp = ((col << 2) + iiter) - 1;
        Subsystem_Fy_B.absxk = std::abs(A[absxk_tmp].re);
        if (Subsystem_Fy_B.absxk > Subsystem_Fy_B.colscale) {
          Subsystem_Fy_B.t = Subsystem_Fy_B.colscale / Subsystem_Fy_B.absxk;
          Subsystem_Fy_B.anorm = Subsystem_Fy_B.anorm * Subsystem_Fy_B.t *
            Subsystem_Fy_B.t + 1.0;
          Subsystem_Fy_B.colscale = Subsystem_Fy_B.absxk;
        } else {
          Subsystem_Fy_B.t = Subsystem_Fy_B.absxk / Subsystem_Fy_B.colscale;
          Subsystem_Fy_B.anorm += Subsystem_Fy_B.t * Subsystem_Fy_B.t;
        }

        Subsystem_Fy_B.absxk = std::abs(A[absxk_tmp].im);
        if (Subsystem_Fy_B.absxk > Subsystem_Fy_B.colscale) {
          Subsystem_Fy_B.t = Subsystem_Fy_B.colscale / Subsystem_Fy_B.absxk;
          Subsystem_Fy_B.anorm = Subsystem_Fy_B.anorm * Subsystem_Fy_B.t *
            Subsystem_Fy_B.t + 1.0;
          Subsystem_Fy_B.colscale = Subsystem_Fy_B.absxk;
        } else {
          Subsystem_Fy_B.t = Subsystem_Fy_B.absxk / Subsystem_Fy_B.colscale;
          Subsystem_Fy_B.anorm += Subsystem_Fy_B.t * Subsystem_Fy_B.t;
        }
      }

      if (Subsystem_Fy_B.scale_e >= Subsystem_Fy_B.colscale) {
        Subsystem_Fy_B.colscale /= Subsystem_Fy_B.scale_e;
        Subsystem_Fy_B.ssq += Subsystem_Fy_B.colscale * Subsystem_Fy_B.colscale *
          Subsystem_Fy_B.anorm;
      } else {
        Subsystem_Fy_B.scale_e /= Subsystem_Fy_B.colscale;
        Subsystem_Fy_B.ssq = Subsystem_Fy_B.scale_e * Subsystem_Fy_B.scale_e *
          Subsystem_Fy_B.ssq + Subsystem_Fy_B.anorm;
        Subsystem_Fy_B.scale_e = Subsystem_Fy_B.colscale;
      }
    }

    Subsystem_Fy_B.anorm = Subsystem_Fy_B.scale_e * std::sqrt(Subsystem_Fy_B.ssq);
  }

  Subsystem_Fy_B.scale_e = 2.2204460492503131E-16 * Subsystem_Fy_B.anorm;
  Subsystem_Fy_B.ssq = 2.2250738585072014E-308;
  if (Subsystem_Fy_B.scale_e > 2.2250738585072014E-308) {
    Subsystem_Fy_B.ssq = Subsystem_Fy_B.scale_e;
  }

  Subsystem_Fy_B.scale_e = 2.2250738585072014E-308;
  if (Subsystem_Fy_B.anorm > 2.2250738585072014E-308) {
    Subsystem_Fy_B.scale_e = Subsystem_Fy_B.anorm;
  }

  Subsystem_Fy_B.anorm = 1.0 / Subsystem_Fy_B.scale_e;
  failed = true;
  for (nm1 = ihi; nm1 + 1 < 5; nm1++) {
    alpha1[nm1] = A[(nm1 << 2) + nm1];
  }

  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    col = ilo;
    nm1 = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    absxk_tmp = 0;
    do {
      exitg1 = 0;
      if (absxk_tmp <= ((ihi - ilo) + 1) * 30 - 1) {
        if (nm1 + 1 == ilo) {
          goto60 = true;
        } else {
          jp1 = (ilastm1 << 2) + nm1;
          if (std::abs(A[jp1].re) + std::abs(A[jp1].im) <= Subsystem_Fy_B.ssq) {
            A[jp1].re = 0.0;
            A[jp1].im = 0.0;
            goto60 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                jp1 = ((j - 1) << 2) + j;
                if (std::abs(A[jp1].re) + std::abs(A[jp1].im) <=
                    Subsystem_Fy_B.ssq) {
                  A[jp1].re = 0.0;
                  A[jp1].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }

            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }
          }
        }

        if ((!goto60) && (!goto70)) {
          alpha1[0].re = (rtNaN);
          alpha1[0].im = 0.0;
          beta1[0].re = (rtNaN);
          beta1[0].im = 0.0;
          alpha1[1].re = (rtNaN);
          alpha1[1].im = 0.0;
          beta1[1].re = (rtNaN);
          beta1[1].im = 0.0;
          alpha1[2].re = (rtNaN);
          alpha1[2].im = 0.0;
          beta1[2].re = (rtNaN);
          beta1[2].im = 0.0;
          alpha1[3].re = (rtNaN);
          alpha1[3].im = 0.0;
          beta1[3].re = (rtNaN);
          beta1[3].im = 0.0;
          for (jp1 = 0; jp1 < 16; jp1++) {
            Z[jp1].re = (rtNaN);
            Z[jp1].im = 0.0;
          }

          *info = 1;
          exitg1 = 1;
        } else if (goto60) {
          goto60 = false;
          alpha1[nm1] = A[(nm1 << 2) + nm1];
          nm1 = ilastm1;
          ilastm1--;
          if (nm1 + 1 < ilo) {
            failed = false;
            guard2 = true;
            exitg1 = 1;
          } else {
            iiter = 0;
            Subsystem_Fy_B.eshift_re = 0.0;
            Subsystem_Fy_B.eshift_im = 0.0;
            absxk_tmp++;
          }
        } else {
          if (goto70) {
            goto70 = false;
            iiter++;
            if (iiter - iiter / 10 * 10 != 0) {
              j = (ilastm1 << 2) + ilastm1;
              Subsystem_Fy_B.ar = A[j].re * Subsystem_Fy_B.anorm;
              Subsystem_Fy_B.ai = A[j].im * Subsystem_Fy_B.anorm;
              if (Subsystem_Fy_B.ai == 0.0) {
                Subsystem_Fy_B.shift.re = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.shift.im = 0.0;
              } else if (Subsystem_Fy_B.ar == 0.0) {
                Subsystem_Fy_B.shift.re = 0.0;
                Subsystem_Fy_B.shift.im = Subsystem_Fy_B.ai / 0.5;
              } else {
                Subsystem_Fy_B.shift.re = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.shift.im = Subsystem_Fy_B.ai / 0.5;
              }

              j = (nm1 << 2) + nm1;
              Subsystem_Fy_B.ar = A[j].re * Subsystem_Fy_B.anorm;
              Subsystem_Fy_B.ai = A[j].im * Subsystem_Fy_B.anorm;
              if (Subsystem_Fy_B.ai == 0.0) {
                Subsystem_Fy_B.ad22.re = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.ad22.im = 0.0;
              } else if (Subsystem_Fy_B.ar == 0.0) {
                Subsystem_Fy_B.ad22.re = 0.0;
                Subsystem_Fy_B.ad22.im = Subsystem_Fy_B.ai / 0.5;
              } else {
                Subsystem_Fy_B.ad22.re = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.ad22.im = Subsystem_Fy_B.ai / 0.5;
              }

              Subsystem_Fy_B.absxk = (Subsystem_Fy_B.shift.re +
                Subsystem_Fy_B.ad22.re) * 0.5;
              Subsystem_Fy_B.t = (Subsystem_Fy_B.shift.im +
                                  Subsystem_Fy_B.ad22.im) * 0.5;
              j = (nm1 << 2) + ilastm1;
              Subsystem_Fy_B.ar = A[j].re * Subsystem_Fy_B.anorm;
              Subsystem_Fy_B.ai = A[j].im * Subsystem_Fy_B.anorm;
              if (Subsystem_Fy_B.ai == 0.0) {
                Subsystem_Fy_B.scale_e = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.colscale = 0.0;
              } else if (Subsystem_Fy_B.ar == 0.0) {
                Subsystem_Fy_B.scale_e = 0.0;
                Subsystem_Fy_B.colscale = Subsystem_Fy_B.ai / 0.5;
              } else {
                Subsystem_Fy_B.scale_e = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.colscale = Subsystem_Fy_B.ai / 0.5;
              }

              j = (ilastm1 << 2) + nm1;
              Subsystem_Fy_B.ar = A[j].re * Subsystem_Fy_B.anorm;
              Subsystem_Fy_B.ai = A[j].im * Subsystem_Fy_B.anorm;
              if (Subsystem_Fy_B.ai == 0.0) {
                Subsystem_Fy_B.ar /= 0.5;
                Subsystem_Fy_B.ai = 0.0;
              } else if (Subsystem_Fy_B.ar == 0.0) {
                Subsystem_Fy_B.ar = 0.0;
                Subsystem_Fy_B.ai /= 0.5;
              } else {
                Subsystem_Fy_B.ar /= 0.5;
                Subsystem_Fy_B.ai /= 0.5;
              }

              Subsystem_Fy_B.shift_im = Subsystem_Fy_B.shift.re *
                Subsystem_Fy_B.ad22.im + Subsystem_Fy_B.shift.im *
                Subsystem_Fy_B.ad22.re;
              Subsystem_Fy_B.shift.re = ((Subsystem_Fy_B.absxk *
                Subsystem_Fy_B.absxk - Subsystem_Fy_B.t * Subsystem_Fy_B.t) +
                (Subsystem_Fy_B.scale_e * Subsystem_Fy_B.ar -
                 Subsystem_Fy_B.colscale * Subsystem_Fy_B.ai)) -
                (Subsystem_Fy_B.shift.re * Subsystem_Fy_B.ad22.re -
                 Subsystem_Fy_B.shift.im * Subsystem_Fy_B.ad22.im);
              Subsystem_Fy_B.shift_tmp = Subsystem_Fy_B.absxk * Subsystem_Fy_B.t;
              Subsystem_Fy_B.shift.im = ((Subsystem_Fy_B.shift_tmp +
                Subsystem_Fy_B.shift_tmp) + (Subsystem_Fy_B.scale_e *
                Subsystem_Fy_B.ai + Subsystem_Fy_B.colscale * Subsystem_Fy_B.ar))
                - Subsystem_Fy_B.shift_im;
              Subsystem_Fy_sqrt(&Subsystem_Fy_B.shift);
              if ((Subsystem_Fy_B.absxk - Subsystem_Fy_B.ad22.re) *
                  Subsystem_Fy_B.shift.re + (Subsystem_Fy_B.t -
                   Subsystem_Fy_B.ad22.im) * Subsystem_Fy_B.shift.im <= 0.0) {
                Subsystem_Fy_B.shift.re += Subsystem_Fy_B.absxk;
                Subsystem_Fy_B.shift.im += Subsystem_Fy_B.t;
              } else {
                Subsystem_Fy_B.shift.re = Subsystem_Fy_B.absxk -
                  Subsystem_Fy_B.shift.re;
                Subsystem_Fy_B.shift.im = Subsystem_Fy_B.t -
                  Subsystem_Fy_B.shift.im;
              }
            } else {
              j = (ilastm1 << 2) + nm1;
              Subsystem_Fy_B.ar = A[j].re * Subsystem_Fy_B.anorm;
              Subsystem_Fy_B.ai = A[j].im * Subsystem_Fy_B.anorm;
              if (Subsystem_Fy_B.ai == 0.0) {
                Subsystem_Fy_B.scale_e = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.colscale = 0.0;
              } else if (Subsystem_Fy_B.ar == 0.0) {
                Subsystem_Fy_B.scale_e = 0.0;
                Subsystem_Fy_B.colscale = Subsystem_Fy_B.ai / 0.5;
              } else {
                Subsystem_Fy_B.scale_e = Subsystem_Fy_B.ar / 0.5;
                Subsystem_Fy_B.colscale = Subsystem_Fy_B.ai / 0.5;
              }

              Subsystem_Fy_B.eshift_re += Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.eshift_im += Subsystem_Fy_B.colscale;
              Subsystem_Fy_B.shift.re = Subsystem_Fy_B.eshift_re;
              Subsystem_Fy_B.shift.im = Subsystem_Fy_B.eshift_im;
            }

            j = ilastm1;
            jp1 = ilastm1 + 1;
            exitg2 = false;
            while ((!exitg2) && (j + 1 > ifirst)) {
              col = j + 1;
              ctemp_tmp_tmp = j << 2;
              ctemp_tmp = ctemp_tmp_tmp + j;
              Subsystem_Fy_B.ctemp.re = A[ctemp_tmp].re * Subsystem_Fy_B.anorm -
                Subsystem_Fy_B.shift.re * 0.5;
              Subsystem_Fy_B.ctemp.im = A[ctemp_tmp].im * Subsystem_Fy_B.anorm -
                Subsystem_Fy_B.shift.im * 0.5;
              Subsystem_Fy_B.scale_e = std::abs(Subsystem_Fy_B.ctemp.re) + std::
                abs(Subsystem_Fy_B.ctemp.im);
              jp1 += ctemp_tmp_tmp;
              Subsystem_Fy_B.colscale = (std::abs(A[jp1].re) + std::abs(A[jp1].
                im)) * Subsystem_Fy_B.anorm;
              Subsystem_Fy_B.absxk = Subsystem_Fy_B.scale_e;
              if (Subsystem_Fy_B.colscale > Subsystem_Fy_B.scale_e) {
                Subsystem_Fy_B.absxk = Subsystem_Fy_B.colscale;
              }

              if ((Subsystem_Fy_B.absxk < 1.0) && (Subsystem_Fy_B.absxk != 0.0))
              {
                Subsystem_Fy_B.scale_e /= Subsystem_Fy_B.absxk;
                Subsystem_Fy_B.colscale /= Subsystem_Fy_B.absxk;
              }

              jp1 = ((j - 1) << 2) + j;
              if ((std::abs(A[jp1].re) + std::abs(A[jp1].im)) *
                  Subsystem_Fy_B.colscale <= Subsystem_Fy_B.scale_e *
                  Subsystem_Fy_B.ssq) {
                goto90 = true;
                exitg2 = true;
              } else {
                jp1 = j;
                j--;
              }
            }

            if (!goto90) {
              col = ifirst;
              ctemp_tmp = (((ifirst - 1) << 2) + ifirst) - 1;
              Subsystem_Fy_B.ctemp.re = A[ctemp_tmp].re * Subsystem_Fy_B.anorm -
                Subsystem_Fy_B.shift.re * 0.5;
              Subsystem_Fy_B.ctemp.im = A[ctemp_tmp].im * Subsystem_Fy_B.anorm -
                Subsystem_Fy_B.shift.im * 0.5;
            }

            goto90 = false;
            j = ((col - 1) << 2) + col;
            Subsystem_Fy_B.ascale.re = A[j].re * Subsystem_Fy_B.anorm;
            Subsystem_Fy_B.ascale.im = A[j].im * Subsystem_Fy_B.anorm;
            Subsystem_Fy_xzlartg_k(Subsystem_Fy_B.ctemp, Subsystem_Fy_B.ascale,
              &Subsystem_Fy_B.scale_e, &Subsystem_Fy_B.ad22);
            j = col;
            jp1 = col - 2;
            while (j < nm1 + 1) {
              if (j > col) {
                Subsystem_Fy_xzlartg(A[(j + (jp1 << 2)) - 1], A[j + (jp1 << 2)],
                                     &Subsystem_Fy_B.scale_e,
                                     &Subsystem_Fy_B.ad22, &A[(j + (jp1 << 2)) -
                                     1]);
                jp1 = j + (jp1 << 2);
                A[jp1].re = 0.0;
                A[jp1].im = 0.0;
              }

              for (ctemp_tmp_tmp = j - 1; ctemp_tmp_tmp + 1 < 5; ctemp_tmp_tmp++)
              {
                jp1 = (ctemp_tmp_tmp << 2) + j;
                Subsystem_Fy_B.shift.re = A[jp1 - 1].re * Subsystem_Fy_B.scale_e
                  + (A[jp1].re * Subsystem_Fy_B.ad22.re - A[jp1].im *
                     Subsystem_Fy_B.ad22.im);
                Subsystem_Fy_B.shift.im = A[jp1 - 1].im * Subsystem_Fy_B.scale_e
                  + (A[jp1].im * Subsystem_Fy_B.ad22.re + A[jp1].re *
                     Subsystem_Fy_B.ad22.im);
                Subsystem_Fy_B.colscale = A[jp1 - 1].im;
                Subsystem_Fy_B.absxk = A[jp1 - 1].re;
                A[jp1].re = A[jp1].re * Subsystem_Fy_B.scale_e - (A[jp1 - 1].re *
                  Subsystem_Fy_B.ad22.re + A[jp1 - 1].im *
                  Subsystem_Fy_B.ad22.im);
                A[jp1].im = A[jp1].im * Subsystem_Fy_B.scale_e -
                  (Subsystem_Fy_B.ad22.re * Subsystem_Fy_B.colscale -
                   Subsystem_Fy_B.ad22.im * Subsystem_Fy_B.absxk);
                A[jp1 - 1] = Subsystem_Fy_B.shift;
              }

              Subsystem_Fy_B.ad22.re = -Subsystem_Fy_B.ad22.re;
              Subsystem_Fy_B.ad22.im = -Subsystem_Fy_B.ad22.im;
              ctemp_tmp_tmp = j;
              if (nm1 + 1 < j + 2) {
                ctemp_tmp_tmp = nm1 - 1;
              }

              for (ctemp_tmp = 0; ctemp_tmp < ctemp_tmp_tmp + 2; ctemp_tmp++) {
                jp1 = ((j - 1) << 2) + ctemp_tmp;
                shift_tmp = (j << 2) + ctemp_tmp;
                Subsystem_Fy_B.shift.re = (A[jp1].re * Subsystem_Fy_B.ad22.re -
                  A[jp1].im * Subsystem_Fy_B.ad22.im) + A[shift_tmp].re *
                  Subsystem_Fy_B.scale_e;
                Subsystem_Fy_B.shift.im = (A[jp1].im * Subsystem_Fy_B.ad22.re +
                  A[jp1].re * Subsystem_Fy_B.ad22.im) + A[shift_tmp].im *
                  Subsystem_Fy_B.scale_e;
                Subsystem_Fy_B.colscale = A[shift_tmp].im;
                Subsystem_Fy_B.absxk = A[shift_tmp].re;
                A[jp1].re = A[jp1].re * Subsystem_Fy_B.scale_e - (A[shift_tmp].
                  re * Subsystem_Fy_B.ad22.re + A[shift_tmp].im *
                  Subsystem_Fy_B.ad22.im);
                A[jp1].im = A[jp1].im * Subsystem_Fy_B.scale_e -
                  (Subsystem_Fy_B.ad22.re * Subsystem_Fy_B.colscale -
                   Subsystem_Fy_B.ad22.im * Subsystem_Fy_B.absxk);
                A[shift_tmp] = Subsystem_Fy_B.shift;
              }

              jp1 = (j - 1) << 2;
              shift_tmp = j << 2;
              Subsystem_Fy_B.shift.re = (Z[jp1].re * Subsystem_Fy_B.ad22.re -
                Z[jp1].im * Subsystem_Fy_B.ad22.im) + Z[shift_tmp].re *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.shift.im = (Z[jp1].im * Subsystem_Fy_B.ad22.re +
                Z[jp1].re * Subsystem_Fy_B.ad22.im) + Z[shift_tmp].im *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.colscale = Z[shift_tmp].im;
              Subsystem_Fy_B.absxk = Z[shift_tmp].re;
              Z[jp1].re = Z[jp1].re * Subsystem_Fy_B.scale_e - (Z[shift_tmp].re *
                Subsystem_Fy_B.ad22.re + Z[shift_tmp].im *
                Subsystem_Fy_B.ad22.im);
              Z[jp1].im = Z[jp1].im * Subsystem_Fy_B.scale_e -
                (Subsystem_Fy_B.ad22.re * Subsystem_Fy_B.colscale -
                 Subsystem_Fy_B.ad22.im * Subsystem_Fy_B.absxk);
              Z[shift_tmp] = Subsystem_Fy_B.shift;
              Subsystem_Fy_B.shift.re = (Z[jp1 + 1].re * Subsystem_Fy_B.ad22.re
                - Z[jp1 + 1].im * Subsystem_Fy_B.ad22.im) + Z[shift_tmp + 1].re *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.shift.im = (Z[jp1 + 1].im * Subsystem_Fy_B.ad22.re
                + Z[jp1 + 1].re * Subsystem_Fy_B.ad22.im) + Z[shift_tmp + 1].im *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.colscale = Z[shift_tmp + 1].im;
              Subsystem_Fy_B.absxk = Z[shift_tmp + 1].re;
              Z[jp1 + 1].re = Z[jp1 + 1].re * Subsystem_Fy_B.scale_e -
                (Z[shift_tmp + 1].re * Subsystem_Fy_B.ad22.re + Z[shift_tmp + 1]
                 .im * Subsystem_Fy_B.ad22.im);
              Z[jp1 + 1].im = Z[jp1 + 1].im * Subsystem_Fy_B.scale_e -
                (Subsystem_Fy_B.ad22.re * Subsystem_Fy_B.colscale -
                 Subsystem_Fy_B.ad22.im * Subsystem_Fy_B.absxk);
              Z[shift_tmp + 1] = Subsystem_Fy_B.shift;
              Subsystem_Fy_B.shift.re = (Z[jp1 + 2].re * Subsystem_Fy_B.ad22.re
                - Z[jp1 + 2].im * Subsystem_Fy_B.ad22.im) + Z[shift_tmp + 2].re *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.shift.im = (Z[jp1 + 2].im * Subsystem_Fy_B.ad22.re
                + Z[jp1 + 2].re * Subsystem_Fy_B.ad22.im) + Z[shift_tmp + 2].im *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.colscale = Z[shift_tmp + 2].im;
              Subsystem_Fy_B.absxk = Z[shift_tmp + 2].re;
              Z[jp1 + 2].re = Z[jp1 + 2].re * Subsystem_Fy_B.scale_e -
                (Z[shift_tmp + 2].re * Subsystem_Fy_B.ad22.re + Z[shift_tmp + 2]
                 .im * Subsystem_Fy_B.ad22.im);
              Z[jp1 + 2].im = Z[jp1 + 2].im * Subsystem_Fy_B.scale_e -
                (Subsystem_Fy_B.ad22.re * Subsystem_Fy_B.colscale -
                 Subsystem_Fy_B.ad22.im * Subsystem_Fy_B.absxk);
              Z[shift_tmp + 2] = Subsystem_Fy_B.shift;
              Subsystem_Fy_B.shift.re = (Z[jp1 + 3].re * Subsystem_Fy_B.ad22.re
                - Z[jp1 + 3].im * Subsystem_Fy_B.ad22.im) + Z[shift_tmp + 3].re *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.shift.im = (Z[jp1 + 3].im * Subsystem_Fy_B.ad22.re
                + Z[jp1 + 3].re * Subsystem_Fy_B.ad22.im) + Z[shift_tmp + 3].im *
                Subsystem_Fy_B.scale_e;
              Subsystem_Fy_B.colscale = Z[shift_tmp + 3].im;
              Subsystem_Fy_B.absxk = Z[shift_tmp + 3].re;
              Z[jp1 + 3].re = Z[jp1 + 3].re * Subsystem_Fy_B.scale_e -
                (Z[shift_tmp + 3].re * Subsystem_Fy_B.ad22.re + Z[shift_tmp + 3]
                 .im * Subsystem_Fy_B.ad22.im);
              Z[jp1 + 3].im = Z[jp1 + 3].im * Subsystem_Fy_B.scale_e -
                (Subsystem_Fy_B.ad22.re * Subsystem_Fy_B.colscale -
                 Subsystem_Fy_B.ad22.im * Subsystem_Fy_B.absxk);
              Z[shift_tmp + 3] = Subsystem_Fy_B.shift;
              jp1 = j - 1;
              j++;
            }
          }

          absxk_tmp++;
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }

  if (guard2) {
    if (failed) {
      *info = nm1 + 1;
      for (ifirst = 0; ifirst <= nm1; ifirst++) {
        alpha1[ifirst].re = (rtNaN);
        alpha1[ifirst].im = 0.0;
        beta1[ifirst].re = (rtNaN);
        beta1[ifirst].im = 0.0;
      }

      for (jp1 = 0; jp1 < 16; jp1++) {
        Z[jp1].re = (rtNaN);
        Z[jp1].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    for (nm1 = 0; nm1 <= ilo - 2; nm1++) {
      alpha1[nm1] = A[(nm1 << 2) + nm1];
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xztgevc(const creal_T A[16], creal_T
  V[16])
{
  real_T work2_idx_2_im;
  real_T work2_idx_3_im;
  real_T work2_idx_3_re;
  int32_T c_j;
  int32_T c_x_tmp;
  int32_T c_x_tmp_tmp;
  int32_T d_re_tmp;
  int32_T d_re_tmp_tmp;
  int32_T e_jr;
  int32_T i;
  int32_T j;
  boolean_T lscalea;
  boolean_T lscaleb;
  Subsystem_Fy_B.rworka[0] = 0.0;
  Subsystem_Fy_B.rworka[2] = 0.0;
  Subsystem_Fy_B.rworka[3] = 0.0;
  Subsystem_Fy_B.anorm_a = std::abs(A[0].re) + std::abs(A[0].im);
  Subsystem_Fy_B.rworka[1] = std::abs(A[4].re) + std::abs(A[4].im);
  Subsystem_Fy_B.ascale_i = (std::abs(A[5].re) + std::abs(A[5].im)) +
    Subsystem_Fy_B.rworka[1];
  if (Subsystem_Fy_B.ascale_i > Subsystem_Fy_B.anorm_a) {
    Subsystem_Fy_B.anorm_a = Subsystem_Fy_B.ascale_i;
  }

  for (i = 0; i < 2; i++) {
    Subsystem_Fy_B.rworka[2] += std::abs(A[i + 8].re) + std::abs(A[i + 8].im);
  }

  Subsystem_Fy_B.ascale_i = (std::abs(A[10].re) + std::abs(A[10].im)) +
    Subsystem_Fy_B.rworka[2];
  if (Subsystem_Fy_B.ascale_i > Subsystem_Fy_B.anorm_a) {
    Subsystem_Fy_B.anorm_a = Subsystem_Fy_B.ascale_i;
  }

  for (i = 0; i < 3; i++) {
    Subsystem_Fy_B.rworka[3] += std::abs(A[i + 12].re) + std::abs(A[i + 12].im);
  }

  Subsystem_Fy_B.ascale_i = (std::abs(A[15].re) + std::abs(A[15].im)) +
    Subsystem_Fy_B.rworka[3];
  if (Subsystem_Fy_B.ascale_i > Subsystem_Fy_B.anorm_a) {
    Subsystem_Fy_B.anorm_a = Subsystem_Fy_B.ascale_i;
  }

  Subsystem_Fy_B.ascale_i = Subsystem_Fy_B.anorm_a;
  if (2.2250738585072014E-308 > Subsystem_Fy_B.anorm_a) {
    Subsystem_Fy_B.ascale_i = 2.2250738585072014E-308;
  }

  Subsystem_Fy_B.ascale_i = 1.0 / Subsystem_Fy_B.ascale_i;
  for (i = 0; i < 4; i++) {
    c_x_tmp_tmp = (3 - i) << 2;
    c_x_tmp = (c_x_tmp_tmp - i) + 3;
    Subsystem_Fy_B.salpha_re = (std::abs(A[c_x_tmp].re) + std::abs(A[c_x_tmp].im))
      * Subsystem_Fy_B.ascale_i;
    if (1.0 > Subsystem_Fy_B.salpha_re) {
      Subsystem_Fy_B.salpha_re = 1.0;
    }

    Subsystem_Fy_B.temp_l = 1.0 / Subsystem_Fy_B.salpha_re;
    Subsystem_Fy_B.salpha_re = A[c_x_tmp].re * Subsystem_Fy_B.temp_l *
      Subsystem_Fy_B.ascale_i;
    Subsystem_Fy_B.salpha_im = A[c_x_tmp].im * Subsystem_Fy_B.temp_l *
      Subsystem_Fy_B.ascale_i;
    Subsystem_Fy_B.acoeff = Subsystem_Fy_B.temp_l * Subsystem_Fy_B.ascale_i;
    lscalea = ((Subsystem_Fy_B.temp_l >= 2.2250738585072014E-308) &&
               (Subsystem_Fy_B.acoeff < 4.0083367200179456E-292));
    Subsystem_Fy_B.dmin = std::abs(Subsystem_Fy_B.salpha_re) + std::abs
      (Subsystem_Fy_B.salpha_im);
    lscaleb = ((Subsystem_Fy_B.dmin >= 2.2250738585072014E-308) &&
               (Subsystem_Fy_B.dmin < 4.0083367200179456E-292));
    Subsystem_Fy_B.scale_o = 1.0;
    if (lscalea) {
      Subsystem_Fy_B.scale_o = Subsystem_Fy_B.anorm_a;
      if (2.4948003869184E+291 < Subsystem_Fy_B.anorm_a) {
        Subsystem_Fy_B.scale_o = 2.4948003869184E+291;
      }

      Subsystem_Fy_B.scale_o *= 4.0083367200179456E-292 / Subsystem_Fy_B.temp_l;
    }

    if (lscaleb) {
      work2_idx_2_im = 4.0083367200179456E-292 / Subsystem_Fy_B.dmin;
      if (work2_idx_2_im > Subsystem_Fy_B.scale_o) {
        Subsystem_Fy_B.scale_o = work2_idx_2_im;
      }
    }

    if (lscalea || lscaleb) {
      work2_idx_2_im = Subsystem_Fy_B.acoeff;
      if (1.0 > Subsystem_Fy_B.acoeff) {
        work2_idx_2_im = 1.0;
      }

      if (Subsystem_Fy_B.dmin > work2_idx_2_im) {
        work2_idx_2_im = Subsystem_Fy_B.dmin;
      }

      Subsystem_Fy_B.dmin = 1.0 / (2.2250738585072014E-308 * work2_idx_2_im);
      if (Subsystem_Fy_B.dmin < Subsystem_Fy_B.scale_o) {
        Subsystem_Fy_B.scale_o = Subsystem_Fy_B.dmin;
      }

      if (lscalea) {
        Subsystem_Fy_B.acoeff = Subsystem_Fy_B.scale_o * Subsystem_Fy_B.temp_l *
          Subsystem_Fy_B.ascale_i;
      } else {
        Subsystem_Fy_B.acoeff *= Subsystem_Fy_B.scale_o;
      }

      Subsystem_Fy_B.salpha_re *= Subsystem_Fy_B.scale_o;
      Subsystem_Fy_B.salpha_im *= Subsystem_Fy_B.scale_o;
    }

    std::memset(&Subsystem_Fy_B.work1[0], 0, sizeof(creal_T) << 2U);
    Subsystem_Fy_B.work1[3 - i].re = 1.0;
    Subsystem_Fy_B.work1[3 - i].im = 0.0;
    Subsystem_Fy_B.dmin = 2.2204460492503131E-16 * Subsystem_Fy_B.acoeff *
      Subsystem_Fy_B.anorm_a;
    Subsystem_Fy_B.temp_l = (std::abs(Subsystem_Fy_B.salpha_re) + std::abs
      (Subsystem_Fy_B.salpha_im)) * 2.2204460492503131E-16;
    if (Subsystem_Fy_B.temp_l > Subsystem_Fy_B.dmin) {
      Subsystem_Fy_B.dmin = Subsystem_Fy_B.temp_l;
    }

    if (2.2250738585072014E-308 > Subsystem_Fy_B.dmin) {
      Subsystem_Fy_B.dmin = 2.2250738585072014E-308;
    }

    for (c_x_tmp = 0; c_x_tmp <= 2 - i; c_x_tmp++) {
      d_re_tmp = c_x_tmp_tmp + c_x_tmp;
      Subsystem_Fy_B.work1[c_x_tmp].re = A[d_re_tmp].re * Subsystem_Fy_B.acoeff;
      Subsystem_Fy_B.work1[c_x_tmp].im = A[d_re_tmp].im * Subsystem_Fy_B.acoeff;
    }

    Subsystem_Fy_B.work1[3 - i].re = 1.0;
    Subsystem_Fy_B.work1[3 - i].im = 0.0;
    c_x_tmp = static_cast<int32_T>(((-1.0 - ((-static_cast<real_T>(i) + 4.0) -
      1.0)) + 1.0) / -1.0);
    for (c_j = 0; c_j < c_x_tmp; c_j++) {
      j = 2 - (i + c_j);
      d_re_tmp_tmp = j << 2;
      d_re_tmp = d_re_tmp_tmp + j;
      work2_idx_3_re = A[d_re_tmp].re * Subsystem_Fy_B.acoeff -
        Subsystem_Fy_B.salpha_re;
      Subsystem_Fy_B.scale_o = A[d_re_tmp].im * Subsystem_Fy_B.acoeff -
        Subsystem_Fy_B.salpha_im;
      if (std::abs(work2_idx_3_re) + std::abs(Subsystem_Fy_B.scale_o) <=
          Subsystem_Fy_B.dmin) {
        work2_idx_3_re = Subsystem_Fy_B.dmin;
        Subsystem_Fy_B.scale_o = 0.0;
      }

      work2_idx_2_im = std::abs(work2_idx_3_re);
      Subsystem_Fy_B.f_y = std::abs(Subsystem_Fy_B.scale_o);
      Subsystem_Fy_B.temp_l = work2_idx_2_im + Subsystem_Fy_B.f_y;
      if (Subsystem_Fy_B.temp_l < 1.0) {
        work2_idx_3_im = std::abs(Subsystem_Fy_B.work1[j].re) + std::abs
          (Subsystem_Fy_B.work1[j].im);
        if (work2_idx_3_im >= Subsystem_Fy_B.temp_l * 1.1235582092889474E+307) {
          Subsystem_Fy_B.temp_l = 1.0 / work2_idx_3_im;
          for (d_re_tmp = 0; d_re_tmp <= 3 - i; d_re_tmp++) {
            Subsystem_Fy_B.work1[d_re_tmp].re *= Subsystem_Fy_B.temp_l;
            Subsystem_Fy_B.work1[d_re_tmp].im *= Subsystem_Fy_B.temp_l;
          }
        }
      }

      Subsystem_Fy_B.temp_l = -Subsystem_Fy_B.work1[j].re;
      work2_idx_3_im = -Subsystem_Fy_B.work1[j].im;
      if (Subsystem_Fy_B.scale_o == 0.0) {
        if (work2_idx_3_im == 0.0) {
          Subsystem_Fy_B.work1[j].re = Subsystem_Fy_B.temp_l / work2_idx_3_re;
          Subsystem_Fy_B.work1[j].im = 0.0;
        } else if (Subsystem_Fy_B.temp_l == 0.0) {
          Subsystem_Fy_B.work1[j].re = 0.0;
          Subsystem_Fy_B.work1[j].im = work2_idx_3_im / work2_idx_3_re;
        } else {
          Subsystem_Fy_B.work1[j].re = Subsystem_Fy_B.temp_l / work2_idx_3_re;
          Subsystem_Fy_B.work1[j].im = work2_idx_3_im / work2_idx_3_re;
        }
      } else if (work2_idx_3_re == 0.0) {
        if (Subsystem_Fy_B.temp_l == 0.0) {
          Subsystem_Fy_B.work1[j].re = work2_idx_3_im / Subsystem_Fy_B.scale_o;
          Subsystem_Fy_B.work1[j].im = 0.0;
        } else if (work2_idx_3_im == 0.0) {
          Subsystem_Fy_B.work1[j].re = 0.0;
          Subsystem_Fy_B.work1[j].im = -(Subsystem_Fy_B.temp_l /
            Subsystem_Fy_B.scale_o);
        } else {
          Subsystem_Fy_B.work1[j].re = work2_idx_3_im / Subsystem_Fy_B.scale_o;
          Subsystem_Fy_B.work1[j].im = -(Subsystem_Fy_B.temp_l /
            Subsystem_Fy_B.scale_o);
        }
      } else if (work2_idx_2_im > Subsystem_Fy_B.f_y) {
        work2_idx_2_im = Subsystem_Fy_B.scale_o / work2_idx_3_re;
        Subsystem_Fy_B.scale_o = work2_idx_2_im * Subsystem_Fy_B.scale_o +
          work2_idx_3_re;
        Subsystem_Fy_B.work1[j].re = (work2_idx_2_im * work2_idx_3_im +
          Subsystem_Fy_B.temp_l) / Subsystem_Fy_B.scale_o;
        Subsystem_Fy_B.work1[j].im = (work2_idx_3_im - work2_idx_2_im *
          Subsystem_Fy_B.temp_l) / Subsystem_Fy_B.scale_o;
      } else if (Subsystem_Fy_B.f_y == work2_idx_2_im) {
        work2_idx_3_re = work2_idx_3_re > 0.0 ? 0.5 : -0.5;
        Subsystem_Fy_B.scale_o = Subsystem_Fy_B.scale_o > 0.0 ? 0.5 : -0.5;
        Subsystem_Fy_B.work1[j].re = (Subsystem_Fy_B.temp_l * work2_idx_3_re +
          work2_idx_3_im * Subsystem_Fy_B.scale_o) / work2_idx_2_im;
        Subsystem_Fy_B.work1[j].im = (work2_idx_3_im * work2_idx_3_re -
          Subsystem_Fy_B.temp_l * Subsystem_Fy_B.scale_o) / work2_idx_2_im;
      } else {
        work2_idx_2_im = work2_idx_3_re / Subsystem_Fy_B.scale_o;
        Subsystem_Fy_B.scale_o += work2_idx_2_im * work2_idx_3_re;
        Subsystem_Fy_B.work1[j].re = (work2_idx_2_im * Subsystem_Fy_B.temp_l +
          work2_idx_3_im) / Subsystem_Fy_B.scale_o;
        Subsystem_Fy_B.work1[j].im = (work2_idx_2_im * work2_idx_3_im -
          Subsystem_Fy_B.temp_l) / Subsystem_Fy_B.scale_o;
      }

      if (j + 1 > 1) {
        if (std::abs(Subsystem_Fy_B.work1[j].re) + std::abs
            (Subsystem_Fy_B.work1[j].im) > 1.0) {
          Subsystem_Fy_B.temp_l = 1.0 / (std::abs(Subsystem_Fy_B.work1[j].re) +
            std::abs(Subsystem_Fy_B.work1[j].im));
          if (Subsystem_Fy_B.acoeff * Subsystem_Fy_B.rworka[j] >=
              1.1235582092889474E+307 * Subsystem_Fy_B.temp_l) {
            for (d_re_tmp = 0; d_re_tmp <= 3 - i; d_re_tmp++) {
              Subsystem_Fy_B.work1[d_re_tmp].re *= Subsystem_Fy_B.temp_l;
              Subsystem_Fy_B.work1[d_re_tmp].im *= Subsystem_Fy_B.temp_l;
            }
          }
        }

        work2_idx_3_re = Subsystem_Fy_B.acoeff * Subsystem_Fy_B.work1[j].re;
        Subsystem_Fy_B.scale_o = Subsystem_Fy_B.acoeff * Subsystem_Fy_B.work1[j]
          .im;
        for (e_jr = 0; e_jr < j; e_jr++) {
          d_re_tmp = d_re_tmp_tmp + e_jr;
          Subsystem_Fy_B.work1[e_jr].re += A[d_re_tmp].re * work2_idx_3_re -
            A[d_re_tmp].im * Subsystem_Fy_B.scale_o;
          Subsystem_Fy_B.work1[e_jr].im += A[d_re_tmp].im * work2_idx_3_re +
            A[d_re_tmp].re * Subsystem_Fy_B.scale_o;
        }
      }
    }

    Subsystem_Fy_B.salpha_re = 0.0;
    Subsystem_Fy_B.salpha_im = 0.0;
    Subsystem_Fy_B.acoeff = 0.0;
    Subsystem_Fy_B.dmin = 0.0;
    Subsystem_Fy_B.scale_o = 0.0;
    work2_idx_2_im = 0.0;
    work2_idx_3_re = 0.0;
    work2_idx_3_im = 0.0;
    for (c_x_tmp = 0; c_x_tmp <= 3 - i; c_x_tmp++) {
      c_j = c_x_tmp << 2;
      Subsystem_Fy_B.salpha_re += V[c_j].re * Subsystem_Fy_B.work1[c_x_tmp].re -
        V[c_j].im * Subsystem_Fy_B.work1[c_x_tmp].im;
      Subsystem_Fy_B.salpha_im += V[c_j].re * Subsystem_Fy_B.work1[c_x_tmp].im +
        V[c_j].im * Subsystem_Fy_B.work1[c_x_tmp].re;
      Subsystem_Fy_B.acoeff += V[c_j + 1].re * Subsystem_Fy_B.work1[c_x_tmp].re
        - V[c_j + 1].im * Subsystem_Fy_B.work1[c_x_tmp].im;
      Subsystem_Fy_B.dmin += V[c_j + 1].re * Subsystem_Fy_B.work1[c_x_tmp].im +
        V[c_j + 1].im * Subsystem_Fy_B.work1[c_x_tmp].re;
      Subsystem_Fy_B.scale_o += V[c_j + 2].re * Subsystem_Fy_B.work1[c_x_tmp].re
        - V[c_j + 2].im * Subsystem_Fy_B.work1[c_x_tmp].im;
      work2_idx_2_im += V[c_j + 2].re * Subsystem_Fy_B.work1[c_x_tmp].im + V[c_j
        + 2].im * Subsystem_Fy_B.work1[c_x_tmp].re;
      work2_idx_3_re += V[c_j + 3].re * Subsystem_Fy_B.work1[c_x_tmp].re - V[c_j
        + 3].im * Subsystem_Fy_B.work1[c_x_tmp].im;
      work2_idx_3_im += V[c_j + 3].re * Subsystem_Fy_B.work1[c_x_tmp].im + V[c_j
        + 3].im * Subsystem_Fy_B.work1[c_x_tmp].re;
    }

    Subsystem_Fy_B.temp_l = std::abs(Subsystem_Fy_B.salpha_re) + std::abs
      (Subsystem_Fy_B.salpha_im);
    Subsystem_Fy_B.f_y = std::abs(Subsystem_Fy_B.acoeff) + std::abs
      (Subsystem_Fy_B.dmin);
    if (Subsystem_Fy_B.f_y > Subsystem_Fy_B.temp_l) {
      Subsystem_Fy_B.temp_l = Subsystem_Fy_B.f_y;
    }

    Subsystem_Fy_B.f_y = std::abs(Subsystem_Fy_B.scale_o) + std::abs
      (work2_idx_2_im);
    if (Subsystem_Fy_B.f_y > Subsystem_Fy_B.temp_l) {
      Subsystem_Fy_B.temp_l = Subsystem_Fy_B.f_y;
    }

    Subsystem_Fy_B.f_y = std::abs(work2_idx_3_re) + std::abs(work2_idx_3_im);
    if (Subsystem_Fy_B.f_y > Subsystem_Fy_B.temp_l) {
      Subsystem_Fy_B.temp_l = Subsystem_Fy_B.f_y;
    }

    if (Subsystem_Fy_B.temp_l > 2.2250738585072014E-308) {
      Subsystem_Fy_B.temp_l = 1.0 / Subsystem_Fy_B.temp_l;
      V[c_x_tmp_tmp].re = Subsystem_Fy_B.temp_l * Subsystem_Fy_B.salpha_re;
      V[c_x_tmp_tmp].im = Subsystem_Fy_B.temp_l * Subsystem_Fy_B.salpha_im;
      d_re_tmp = ((3 - i) << 2) + 1;
      V[d_re_tmp].re = Subsystem_Fy_B.temp_l * Subsystem_Fy_B.acoeff;
      V[d_re_tmp].im = Subsystem_Fy_B.temp_l * Subsystem_Fy_B.dmin;
      d_re_tmp = ((3 - i) << 2) + 2;
      V[d_re_tmp].re = Subsystem_Fy_B.temp_l * Subsystem_Fy_B.scale_o;
      V[d_re_tmp].im = Subsystem_Fy_B.temp_l * work2_idx_2_im;
      d_re_tmp = ((3 - i) << 2) + 3;
      V[d_re_tmp].re = Subsystem_Fy_B.temp_l * work2_idx_3_re;
      V[d_re_tmp].im = Subsystem_Fy_B.temp_l * work2_idx_3_im;
    } else {
      V[c_x_tmp_tmp].re = 0.0;
      V[c_x_tmp_tmp].im = 0.0;
      V[c_x_tmp_tmp + 1].re = 0.0;
      V[c_x_tmp_tmp + 1].im = 0.0;
      V[c_x_tmp_tmp + 2].re = 0.0;
      V[c_x_tmp_tmp + 2].im = 0.0;
      V[c_x_tmp_tmp + 3].re = 0.0;
      V[c_x_tmp_tmp + 3].im = 0.0;
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
real_T Subsystem_FyModelClass::Subsystem_Fy_xnrm2(int32_T n, const real_T x[16],
  int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      for (k = ix0; k <= ix0 + 1; k++) {
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xzlarf(int32_T m, int32_T n, int32_T
  iv0, real_T tau, real_T C[16], int32_T ic0, real_T work[4])
{
  real_T c;
  int32_T coltop;
  int32_T d;
  int32_T exitg1;
  int32_T ia;
  int32_T iac;
  int32_T ix;
  int32_T jy;
  int32_T lastc;
  int32_T lastv;
  boolean_T exitg2;
  if (tau != 0.0) {
    lastv = m;
    lastc = iv0 + m;
    while ((lastv > 0) && (C[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      coltop = (lastc << 2) + ic0;
      jy = coltop;
      do {
        exitg1 = 0;
        if (jy <= (coltop + lastv) - 1) {
          if (C[jy - 1] != 0.0) {
            exitg1 = 1;
          } else {
            jy++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    if (lastc + 1 != 0) {
      for (coltop = 0; coltop <= lastc; coltop++) {
        work[coltop] = 0.0;
      }

      coltop = 0;
      jy = (lastc << 2) + ic0;
      for (iac = ic0; iac <= jy; iac += 4) {
        ix = iv0;
        c = 0.0;
        d = (iac + lastv) - 1;
        for (ia = iac; ia <= d; ia++) {
          c += C[ia - 1] * C[ix - 1];
          ix++;
        }

        work[coltop] += c;
        coltop++;
      }
    }

    if (!(-tau == 0.0)) {
      coltop = ic0 - 1;
      jy = 0;
      for (iac = 0; iac <= lastc; iac++) {
        if (work[jy] != 0.0) {
          c = work[jy] * -tau;
          ix = iv0;
          d = lastv + coltop;
          for (ia = coltop; ia < d; ia++) {
            C[ia] += C[ix - 1] * c;
            ix++;
          }
        }

        jy++;
        coltop += 4;
      }
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
real_T Subsystem_FyModelClass::Subsystem_Fy_xnrm2_d(int32_T n, const real_T x[3])
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[1]);
    } else {
      Subsystem_Fy_B.scale_o4 = 3.3121686421112381E-170;
      Subsystem_Fy_B.absxk_h = std::abs(x[1]);
      if (Subsystem_Fy_B.absxk_h > 3.3121686421112381E-170) {
        y = 1.0;
        Subsystem_Fy_B.scale_o4 = Subsystem_Fy_B.absxk_h;
      } else {
        Subsystem_Fy_B.t_l = Subsystem_Fy_B.absxk_h / 3.3121686421112381E-170;
        y = Subsystem_Fy_B.t_l * Subsystem_Fy_B.t_l;
      }

      Subsystem_Fy_B.absxk_h = std::abs(x[2]);
      if (Subsystem_Fy_B.absxk_h > Subsystem_Fy_B.scale_o4) {
        Subsystem_Fy_B.t_l = Subsystem_Fy_B.scale_o4 / Subsystem_Fy_B.absxk_h;
        y = y * Subsystem_Fy_B.t_l * Subsystem_Fy_B.t_l + 1.0;
        Subsystem_Fy_B.scale_o4 = Subsystem_Fy_B.absxk_h;
      } else {
        Subsystem_Fy_B.t_l = Subsystem_Fy_B.absxk_h / Subsystem_Fy_B.scale_o4;
        y += Subsystem_Fy_B.t_l * Subsystem_Fy_B.t_l;
      }

      y = Subsystem_Fy_B.scale_o4 * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
real_T Subsystem_FyModelClass::Subsystem_Fy_xzlarfg(int32_T n, real_T *alpha1,
  real_T x[3])
{
  real_T tau;
  int32_T c_k;
  int32_T knt;
  tau = 0.0;
  if (n > 0) {
    Subsystem_Fy_B.xnorm = Subsystem_Fy_xnrm2_d(n - 1, x);
    if (Subsystem_Fy_B.xnorm != 0.0) {
      Subsystem_Fy_B.xnorm = Subsystem_Fy_rt_hypotd_snf(*alpha1,
        Subsystem_Fy_B.xnorm);
      if (*alpha1 >= 0.0) {
        Subsystem_Fy_B.xnorm = -Subsystem_Fy_B.xnorm;
      }

      if (std::abs(Subsystem_Fy_B.xnorm) < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          for (c_k = 1; c_k < n; c_k++) {
            x[c_k] *= 9.9792015476736E+291;
          }

          Subsystem_Fy_B.xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(Subsystem_Fy_B.xnorm) >= 1.0020841800044864E-292));

        Subsystem_Fy_B.xnorm = Subsystem_Fy_rt_hypotd_snf(*alpha1,
          Subsystem_Fy_xnrm2_d(n - 1, x));
        if (*alpha1 >= 0.0) {
          Subsystem_Fy_B.xnorm = -Subsystem_Fy_B.xnorm;
        }

        tau = (Subsystem_Fy_B.xnorm - *alpha1) / Subsystem_Fy_B.xnorm;
        *alpha1 = 1.0 / (*alpha1 - Subsystem_Fy_B.xnorm);
        for (c_k = 1; c_k < n; c_k++) {
          x[c_k] *= *alpha1;
        }

        for (c_k = 0; c_k <= knt; c_k++) {
          Subsystem_Fy_B.xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = Subsystem_Fy_B.xnorm;
      } else {
        tau = (Subsystem_Fy_B.xnorm - *alpha1) / Subsystem_Fy_B.xnorm;
        *alpha1 = 1.0 / (*alpha1 - Subsystem_Fy_B.xnorm);
        for (knt = 1; knt < n; knt++) {
          x[knt] *= *alpha1;
        }

        *alpha1 = Subsystem_Fy_B.xnorm;
      }
    }
  }

  return tau;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xdlanv2(real_T *a, real_T *b, real_T
  *c, real_T *d, real_T *rt1r, real_T *rt1i, real_T *rt2r, real_T *rt2i, real_T *
  cs, real_T *sn)
{
  int32_T b_0;
  int32_T c_0;
  boolean_T tmp;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    Subsystem_Fy_B.bcmax = *d;
    *d = *a;
    *a = Subsystem_Fy_B.bcmax;
    *b = -*c;
    *c = 0.0;
  } else {
    Subsystem_Fy_B.tau = *a - *d;
    if ((Subsystem_Fy_B.tau == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      Subsystem_Fy_B.p = 0.5 * Subsystem_Fy_B.tau;
      Subsystem_Fy_B.bcmis = std::abs(*b);
      Subsystem_Fy_B.z = std::abs(*c);
      tmp = rtIsNaN(Subsystem_Fy_B.z);
      if ((Subsystem_Fy_B.bcmis > Subsystem_Fy_B.z) || tmp) {
        Subsystem_Fy_B.bcmax = Subsystem_Fy_B.bcmis;
      } else {
        Subsystem_Fy_B.bcmax = Subsystem_Fy_B.z;
      }

      if ((Subsystem_Fy_B.bcmis < Subsystem_Fy_B.z) || tmp) {
        Subsystem_Fy_B.z = Subsystem_Fy_B.bcmis;
      }

      if (!(*b < 0.0)) {
        b_0 = 1;
      } else {
        b_0 = -1;
      }

      if (!(*c < 0.0)) {
        c_0 = 1;
      } else {
        c_0 = -1;
      }

      Subsystem_Fy_B.bcmis = Subsystem_Fy_B.z * static_cast<real_T>(b_0) *
        static_cast<real_T>(c_0);
      Subsystem_Fy_B.scale = std::abs(Subsystem_Fy_B.p);
      if ((!(Subsystem_Fy_B.scale > Subsystem_Fy_B.bcmax)) && (!rtIsNaN
           (Subsystem_Fy_B.bcmax))) {
        Subsystem_Fy_B.scale = Subsystem_Fy_B.bcmax;
      }

      Subsystem_Fy_B.z = Subsystem_Fy_B.p / Subsystem_Fy_B.scale *
        Subsystem_Fy_B.p + Subsystem_Fy_B.bcmax / Subsystem_Fy_B.scale *
        Subsystem_Fy_B.bcmis;
      if (Subsystem_Fy_B.z >= 8.8817841970012523E-16) {
        if (!(Subsystem_Fy_B.p < 0.0)) {
          Subsystem_Fy_B.tau = std::sqrt(Subsystem_Fy_B.scale) * std::sqrt
            (Subsystem_Fy_B.z);
        } else {
          Subsystem_Fy_B.tau = -(std::sqrt(Subsystem_Fy_B.scale) * std::sqrt
            (Subsystem_Fy_B.z));
        }

        Subsystem_Fy_B.z = Subsystem_Fy_B.p + Subsystem_Fy_B.tau;
        *a = *d + Subsystem_Fy_B.z;
        *d -= Subsystem_Fy_B.bcmax / Subsystem_Fy_B.z * Subsystem_Fy_B.bcmis;
        Subsystem_Fy_B.tau = Subsystem_Fy_rt_hypotd_snf(*c, Subsystem_Fy_B.z);
        *cs = Subsystem_Fy_B.z / Subsystem_Fy_B.tau;
        *sn = *c / Subsystem_Fy_B.tau;
        *b -= *c;
        *c = 0.0;
      } else {
        Subsystem_Fy_B.bcmax = *b + *c;
        Subsystem_Fy_B.tau = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.bcmax,
          Subsystem_Fy_B.tau);
        *cs = std::sqrt((std::abs(Subsystem_Fy_B.bcmax) / Subsystem_Fy_B.tau +
                         1.0) * 0.5);
        if (!(Subsystem_Fy_B.bcmax < 0.0)) {
          b_0 = 1;
        } else {
          b_0 = -1;
        }

        *sn = -(Subsystem_Fy_B.p / (Subsystem_Fy_B.tau * *cs)) *
          static_cast<real_T>(b_0);
        Subsystem_Fy_B.p = *a * *cs + *b * *sn;
        Subsystem_Fy_B.tau = -*a * *sn + *b * *cs;
        Subsystem_Fy_B.bcmax = *c * *cs + *d * *sn;
        Subsystem_Fy_B.bcmis = -*c * *sn + *d * *cs;
        *b = Subsystem_Fy_B.tau * *cs + Subsystem_Fy_B.bcmis * *sn;
        *c = -Subsystem_Fy_B.p * *sn + Subsystem_Fy_B.bcmax * *cs;
        Subsystem_Fy_B.bcmax = ((Subsystem_Fy_B.p * *cs + Subsystem_Fy_B.bcmax *
          *sn) + (-Subsystem_Fy_B.tau * *sn + Subsystem_Fy_B.bcmis * *cs)) * 0.5;
        *a = Subsystem_Fy_B.bcmax;
        *d = Subsystem_Fy_B.bcmax;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              Subsystem_Fy_B.z = std::sqrt(std::abs(*b));
              Subsystem_Fy_B.bcmis = std::sqrt(std::abs(*c));
              if (!(*c < 0.0)) {
                Subsystem_Fy_B.p = Subsystem_Fy_B.z * Subsystem_Fy_B.bcmis;
              } else {
                Subsystem_Fy_B.p = -(Subsystem_Fy_B.z * Subsystem_Fy_B.bcmis);
              }

              Subsystem_Fy_B.tau = 1.0 / std::sqrt(std::abs(*b + *c));
              *a = Subsystem_Fy_B.bcmax + Subsystem_Fy_B.p;
              *d = Subsystem_Fy_B.bcmax - Subsystem_Fy_B.p;
              *b -= *c;
              *c = 0.0;
              Subsystem_Fy_B.p = Subsystem_Fy_B.z * Subsystem_Fy_B.tau;
              Subsystem_Fy_B.tau *= Subsystem_Fy_B.bcmis;
              Subsystem_Fy_B.bcmax = *cs * Subsystem_Fy_B.p - *sn *
                Subsystem_Fy_B.tau;
              *sn = *cs * Subsystem_Fy_B.tau + *sn * Subsystem_Fy_B.p;
              *cs = Subsystem_Fy_B.bcmax;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            Subsystem_Fy_B.bcmax = *cs;
            *cs = -*sn;
            *sn = Subsystem_Fy_B.bcmax;
          }
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = std::sqrt(std::abs(*b)) * std::sqrt(std::abs(*c));
    *rt2i = -*rt1i;
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xrot(int32_T n, real_T x[16], int32_T
  ix0, int32_T iy0, real_T c, real_T s)
{
  int32_T ix;
  int32_T iy;
  int32_T k;
  if (n >= 1) {
    ix = ix0 - 1;
    iy = iy0 - 1;
    for (k = 0; k < n; k++) {
      Subsystem_Fy_B.temp_h = c * x[ix] + s * x[iy];
      x[iy] = c * x[iy] - s * x[ix];
      x[ix] = Subsystem_Fy_B.temp_h;
      iy++;
      ix++;
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xrot_m(real_T x[16], int32_T ix0,
  int32_T iy0, real_T c, real_T s)
{
  real_T temp_tmp;
  Subsystem_Fy_B.temp_m = x[iy0 - 1];
  temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = Subsystem_Fy_B.temp_m * c - temp_tmp * s;
  x[ix0 - 1] = temp_tmp * c + Subsystem_Fy_B.temp_m * s;
  Subsystem_Fy_B.temp_m = x[ix0] * c + x[iy0] * s;
  x[iy0] = x[iy0] * c - x[ix0] * s;
  x[ix0] = Subsystem_Fy_B.temp_m;
  Subsystem_Fy_B.temp_m = x[iy0 + 1];
  temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = Subsystem_Fy_B.temp_m * c - temp_tmp * s;
  x[ix0 + 1] = temp_tmp * c + Subsystem_Fy_B.temp_m * s;
  Subsystem_Fy_B.temp_m = x[iy0 + 2];
  temp_tmp = x[ix0 + 2];
  x[iy0 + 2] = Subsystem_Fy_B.temp_m * c - temp_tmp * s;
  x[ix0 + 2] = temp_tmp * c + Subsystem_Fy_B.temp_m * s;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
int32_T Subsystem_FyModelClass::Subsystem_Fy_xhseqr(real_T h[16], real_T z[16])
{
  int32_T L;
  int32_T c_j;
  int32_T hoffset;
  int32_T i;
  int32_T info;
  int32_T ix;
  int32_T k;
  int32_T m;
  int32_T m_tmp;
  int32_T nr;
  int32_T s_tmp;
  int32_T sum1_tmp;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T goto150;
  info = 0;
  Subsystem_Fy_B.v[0] = 0.0;
  Subsystem_Fy_B.v[1] = 0.0;
  Subsystem_Fy_B.v[2] = 0.0;
  h[2] = 0.0;
  h[3] = 0.0;
  h[7] = 0.0;
  i = 3;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    ix = 0;
    exitg2 = false;
    while ((!exitg2) && (ix < 301)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        s_tmp = ((k - 1) << 2) + k;
        Subsystem_Fy_B.htmp2 = std::abs(h[s_tmp]);
        if (Subsystem_Fy_B.htmp2 <= 4.0083367200179456E-292) {
          exitg3 = true;
        } else {
          m_tmp = (k << 2) + k;
          Subsystem_Fy_B.tst = std::abs(h[s_tmp - 1]) + std::abs(h[m_tmp]);
          if (Subsystem_Fy_B.tst == 0.0) {
            if (k - 1 >= 1) {
              Subsystem_Fy_B.tst = std::abs(h[(((k - 2) << 2) + k) - 1]);
            }

            if (k + 2 <= 4) {
              Subsystem_Fy_B.tst += std::abs(h[m_tmp + 1]);
            }
          }

          if (Subsystem_Fy_B.htmp2 <= 2.2204460492503131E-16 *
              Subsystem_Fy_B.tst) {
            Subsystem_Fy_B.htmp1 = std::abs(h[s_tmp]);
            Subsystem_Fy_B.htmp2 = std::abs(h[m_tmp - 1]);
            if (Subsystem_Fy_B.htmp1 > Subsystem_Fy_B.htmp2) {
              Subsystem_Fy_B.tst = Subsystem_Fy_B.htmp1;
              Subsystem_Fy_B.ba = Subsystem_Fy_B.htmp2;
            } else {
              Subsystem_Fy_B.tst = Subsystem_Fy_B.htmp2;
              Subsystem_Fy_B.ba = Subsystem_Fy_B.htmp1;
            }

            Subsystem_Fy_B.htmp2 = h[m_tmp];
            Subsystem_Fy_B.htmp1 = std::abs(Subsystem_Fy_B.htmp2);
            Subsystem_Fy_B.htmp2 = std::abs(h[s_tmp - 1] - Subsystem_Fy_B.htmp2);
            if (Subsystem_Fy_B.htmp1 > Subsystem_Fy_B.htmp2) {
              Subsystem_Fy_B.aa = Subsystem_Fy_B.htmp1;
              Subsystem_Fy_B.htmp1 = Subsystem_Fy_B.htmp2;
            } else {
              Subsystem_Fy_B.aa = Subsystem_Fy_B.htmp2;
            }

            Subsystem_Fy_B.htmp2 = Subsystem_Fy_B.aa + Subsystem_Fy_B.tst;
            Subsystem_Fy_B.htmp1 = Subsystem_Fy_B.aa / Subsystem_Fy_B.htmp2 *
              Subsystem_Fy_B.htmp1 * 2.2204460492503131E-16;
            if ((4.0083367200179456E-292 > Subsystem_Fy_B.htmp1) || rtIsNaN
                (Subsystem_Fy_B.htmp1)) {
              Subsystem_Fy_B.htmp1 = 4.0083367200179456E-292;
            }

            if (Subsystem_Fy_B.tst / Subsystem_Fy_B.htmp2 * Subsystem_Fy_B.ba <=
                Subsystem_Fy_B.htmp1) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }

      L = k + 1;
      if (k + 1 > 1) {
        h[k + ((k - 1) << 2)] = 0.0;
      }

      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        switch (ix) {
         case 10:
          s_tmp = (k << 2) + k;
          Subsystem_Fy_B.htmp2 = std::abs(h[(((k + 1) << 2) + k) + 2]) + std::
            abs(h[s_tmp + 1]);
          Subsystem_Fy_B.tst = h[s_tmp] + 0.75 * Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.h12 = -0.4375 * Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.aa = Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.htmp1 = Subsystem_Fy_B.tst;
          break;

         case 20:
          Subsystem_Fy_B.htmp2 = std::abs(h[(((i - 2) << 2) + i) - 1]) + std::
            abs(h[((i - 1) << 2) + i]);
          Subsystem_Fy_B.tst = h[(i << 2) + i] + 0.75 * Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.h12 = -0.4375 * Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.aa = Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.htmp1 = Subsystem_Fy_B.tst;
          break;

         default:
          Subsystem_Fy_B.tst = h[(((i - 1) << 2) + i) - 1];
          Subsystem_Fy_B.aa = h[((i - 1) << 2) + i];
          m_tmp = (i << 2) + i;
          Subsystem_Fy_B.h12 = h[m_tmp - 1];
          Subsystem_Fy_B.htmp1 = h[m_tmp];
          break;
        }

        Subsystem_Fy_B.htmp2 = ((std::abs(Subsystem_Fy_B.tst) + std::abs
          (Subsystem_Fy_B.h12)) + std::abs(Subsystem_Fy_B.aa)) + std::abs
          (Subsystem_Fy_B.htmp1);
        if (Subsystem_Fy_B.htmp2 == 0.0) {
          Subsystem_Fy_B.tst = 0.0;
          Subsystem_Fy_B.htmp1 = 0.0;
          Subsystem_Fy_B.ba = 0.0;
          Subsystem_Fy_B.aa = 0.0;
        } else {
          Subsystem_Fy_B.tst /= Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.htmp1 /= Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.ba = (Subsystem_Fy_B.tst + Subsystem_Fy_B.htmp1) / 2.0;
          Subsystem_Fy_B.tst = (Subsystem_Fy_B.tst - Subsystem_Fy_B.ba) *
            (Subsystem_Fy_B.htmp1 - Subsystem_Fy_B.ba) - Subsystem_Fy_B.h12 /
            Subsystem_Fy_B.htmp2 * (Subsystem_Fy_B.aa / Subsystem_Fy_B.htmp2);
          Subsystem_Fy_B.aa = std::sqrt(std::abs(Subsystem_Fy_B.tst));
          if (Subsystem_Fy_B.tst >= 0.0) {
            Subsystem_Fy_B.tst = Subsystem_Fy_B.ba * Subsystem_Fy_B.htmp2;
            Subsystem_Fy_B.ba = Subsystem_Fy_B.tst;
            Subsystem_Fy_B.htmp1 = Subsystem_Fy_B.aa * Subsystem_Fy_B.htmp2;
            Subsystem_Fy_B.aa = -Subsystem_Fy_B.htmp1;
          } else {
            Subsystem_Fy_B.tst = Subsystem_Fy_B.ba + Subsystem_Fy_B.aa;
            Subsystem_Fy_B.ba -= Subsystem_Fy_B.aa;
            if (std::abs(Subsystem_Fy_B.tst - Subsystem_Fy_B.htmp1) <= std::abs
                (Subsystem_Fy_B.ba - Subsystem_Fy_B.htmp1)) {
              Subsystem_Fy_B.tst *= Subsystem_Fy_B.htmp2;
              Subsystem_Fy_B.ba = Subsystem_Fy_B.tst;
            } else {
              Subsystem_Fy_B.ba *= Subsystem_Fy_B.htmp2;
              Subsystem_Fy_B.tst = Subsystem_Fy_B.ba;
            }

            Subsystem_Fy_B.htmp1 = 0.0;
            Subsystem_Fy_B.aa = 0.0;
          }
        }

        m = i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          s_tmp = ((m - 1) << 2) + m;
          Subsystem_Fy_B.h21s = h[s_tmp];
          Subsystem_Fy_B.h12 = h[s_tmp - 1];
          Subsystem_Fy_B.s_tmp = Subsystem_Fy_B.h12 - Subsystem_Fy_B.ba;
          Subsystem_Fy_B.htmp2 = (std::abs(Subsystem_Fy_B.s_tmp) + std::abs
            (Subsystem_Fy_B.aa)) + std::abs(Subsystem_Fy_B.h21s);
          Subsystem_Fy_B.h21s /= Subsystem_Fy_B.htmp2;
          s_tmp = (m << 2) + m;
          Subsystem_Fy_B.v[0] = (Subsystem_Fy_B.s_tmp / Subsystem_Fy_B.htmp2 *
            (Subsystem_Fy_B.h12 - Subsystem_Fy_B.tst) + h[s_tmp - 1] *
            Subsystem_Fy_B.h21s) - Subsystem_Fy_B.aa / Subsystem_Fy_B.htmp2 *
            Subsystem_Fy_B.htmp1;
          Subsystem_Fy_B.s_tmp = h[s_tmp];
          Subsystem_Fy_B.v[1] = (((Subsystem_Fy_B.h12 + Subsystem_Fy_B.s_tmp) -
            Subsystem_Fy_B.tst) - Subsystem_Fy_B.ba) * Subsystem_Fy_B.h21s;
          Subsystem_Fy_B.v[2] = h[s_tmp + 1] * Subsystem_Fy_B.h21s;
          Subsystem_Fy_B.htmp2 = (std::abs(Subsystem_Fy_B.v[0]) + std::abs
            (Subsystem_Fy_B.v[1])) + std::abs(Subsystem_Fy_B.v[2]);
          Subsystem_Fy_B.v[0] /= Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.v[1] /= Subsystem_Fy_B.htmp2;
          Subsystem_Fy_B.v[2] /= Subsystem_Fy_B.htmp2;
          if (k + 1 == m) {
            exitg3 = true;
          } else {
            s_tmp = ((m - 2) << 2) + m;
            if (std::abs(h[s_tmp - 1]) * (std::abs(Subsystem_Fy_B.v[1]) + std::
                 abs(Subsystem_Fy_B.v[2])) <= ((std::abs(h[s_tmp - 2]) + std::
                  abs(Subsystem_Fy_B.h12)) + std::abs(Subsystem_Fy_B.s_tmp)) *
                (2.2204460492503131E-16 * std::abs(Subsystem_Fy_B.v[0]))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }

        for (s_tmp = m; s_tmp <= i; s_tmp++) {
          nr = (i - s_tmp) + 2;
          if (3 < nr) {
            nr = 3;
          }

          if (s_tmp > m) {
            hoffset = ((s_tmp - 2) << 2) + s_tmp;
            for (m_tmp = 0; m_tmp < nr; m_tmp++) {
              Subsystem_Fy_B.v[m_tmp] = h[(m_tmp + hoffset) - 1];
            }
          }

          Subsystem_Fy_B.tst = Subsystem_Fy_B.v[0];
          Subsystem_Fy_B.htmp2 = Subsystem_Fy_xzlarfg(nr, &Subsystem_Fy_B.tst,
            Subsystem_Fy_B.v);
          Subsystem_Fy_B.v[0] = Subsystem_Fy_B.tst;
          if (s_tmp > m) {
            h[(s_tmp + ((s_tmp - 2) << 2)) - 1] = Subsystem_Fy_B.tst;
            h[s_tmp + ((s_tmp - 2) << 2)] = 0.0;
            if (s_tmp < i) {
              h[s_tmp + 1] = 0.0;
            }
          } else if (m > k + 1) {
            h[s_tmp - 1] *= 1.0 - Subsystem_Fy_B.htmp2;
          }

          Subsystem_Fy_B.tst = Subsystem_Fy_B.v[1];
          Subsystem_Fy_B.ba = Subsystem_Fy_B.htmp2 * Subsystem_Fy_B.v[1];
          switch (nr) {
           case 3:
            Subsystem_Fy_B.aa = Subsystem_Fy_B.v[2];
            Subsystem_Fy_B.h12 = Subsystem_Fy_B.htmp2 * Subsystem_Fy_B.v[2];
            for (hoffset = s_tmp - 1; hoffset + 1 < 5; hoffset++) {
              nr = (hoffset << 2) + s_tmp;
              Subsystem_Fy_B.htmp1 = (h[nr - 1] + h[nr] * Subsystem_Fy_B.tst) +
                h[nr + 1] * Subsystem_Fy_B.aa;
              h[nr - 1] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.htmp2;
              h[nr] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.ba;
              h[nr + 1] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.h12;
            }

            if (s_tmp + 3 < i + 1) {
              m_tmp = s_tmp + 3;
            } else {
              m_tmp = i + 1;
            }

            for (c_j = 0; c_j < m_tmp; c_j++) {
              nr = ((s_tmp - 1) << 2) + c_j;
              hoffset = (s_tmp << 2) + c_j;
              sum1_tmp = ((s_tmp + 1) << 2) + c_j;
              Subsystem_Fy_B.htmp1 = (h[nr] + h[hoffset] * Subsystem_Fy_B.tst) +
                h[sum1_tmp] * Subsystem_Fy_B.aa;
              h[nr] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.htmp2;
              h[hoffset] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.ba;
              h[sum1_tmp] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.h12;
            }

            for (m_tmp = 0; m_tmp < 4; m_tmp++) {
              nr = ((s_tmp - 1) << 2) + m_tmp;
              hoffset = (s_tmp << 2) + m_tmp;
              sum1_tmp = ((s_tmp + 1) << 2) + m_tmp;
              Subsystem_Fy_B.htmp1 = (z[nr] + z[hoffset] * Subsystem_Fy_B.tst) +
                z[sum1_tmp] * Subsystem_Fy_B.aa;
              z[nr] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.htmp2;
              z[hoffset] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.ba;
              z[sum1_tmp] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.h12;
            }
            break;

           case 2:
            for (hoffset = s_tmp - 1; hoffset + 1 < 5; hoffset++) {
              nr = (hoffset << 2) + s_tmp;
              Subsystem_Fy_B.aa = h[nr - 1];
              Subsystem_Fy_B.htmp1 = Subsystem_Fy_B.aa + h[nr] *
                Subsystem_Fy_B.tst;
              h[nr - 1] = Subsystem_Fy_B.aa - Subsystem_Fy_B.htmp1 *
                Subsystem_Fy_B.htmp2;
              h[nr] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.ba;
            }

            for (m_tmp = 0; m_tmp <= i; m_tmp++) {
              nr = ((s_tmp - 1) << 2) + m_tmp;
              hoffset = (s_tmp << 2) + m_tmp;
              Subsystem_Fy_B.htmp1 = h[nr] + h[hoffset] * Subsystem_Fy_B.tst;
              h[nr] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.htmp2;
              h[hoffset] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.ba;
            }

            for (m_tmp = 0; m_tmp < 4; m_tmp++) {
              nr = ((s_tmp - 1) << 2) + m_tmp;
              Subsystem_Fy_B.aa = z[nr];
              hoffset = (s_tmp << 2) + m_tmp;
              Subsystem_Fy_B.htmp1 = Subsystem_Fy_B.aa + z[hoffset] *
                Subsystem_Fy_B.tst;
              z[nr] = Subsystem_Fy_B.aa - Subsystem_Fy_B.htmp1 *
                Subsystem_Fy_B.htmp2;
              z[hoffset] -= Subsystem_Fy_B.htmp1 * Subsystem_Fy_B.ba;
            }
            break;
          }
        }

        ix++;
      }
    }

    if (!goto150) {
      info = i + 1;
      exitg1 = true;
    } else {
      if ((i + 1 != L) && (L == i)) {
        ix = ((i - 1) << 2) + i;
        Subsystem_Fy_B.ba = h[ix - 1];
        k = (i << 2) + i;
        Subsystem_Fy_B.htmp1 = h[k - 1];
        Subsystem_Fy_B.aa = h[ix];
        Subsystem_Fy_B.h12 = h[k];
        Subsystem_Fy_xdlanv2(&Subsystem_Fy_B.ba, &Subsystem_Fy_B.htmp1,
                             &Subsystem_Fy_B.aa, &Subsystem_Fy_B.h12,
                             &Subsystem_Fy_B.s_tmp, &Subsystem_Fy_B.h21s,
                             &Subsystem_Fy_B.a__3, &Subsystem_Fy_B.a__4,
                             &Subsystem_Fy_B.htmp2, &Subsystem_Fy_B.tst);
        h[ix - 1] = Subsystem_Fy_B.ba;
        h[k - 1] = Subsystem_Fy_B.htmp1;
        h[ix] = Subsystem_Fy_B.aa;
        h[k] = Subsystem_Fy_B.h12;
        if (4 > i + 1) {
          k = ((i + 1) << 2) + i;
          ix = k - 1;
          for (m = 0; m <= 2 - i; m++) {
            Subsystem_Fy_B.ba = Subsystem_Fy_B.htmp2 * h[ix] +
              Subsystem_Fy_B.tst * h[k];
            h[k] = Subsystem_Fy_B.htmp2 * h[k] - Subsystem_Fy_B.tst * h[ix];
            h[ix] = Subsystem_Fy_B.ba;
            k += 4;
            ix += 4;
          }
        }

        Subsystem_Fy_xrot(i - 1, h, ((i - 1) << 2) + 1, (i << 2) + 1,
                          Subsystem_Fy_B.htmp2, Subsystem_Fy_B.tst);
        Subsystem_Fy_xrot_m(z, ((i - 1) << 2) + 1, (i << 2) + 1,
                            Subsystem_Fy_B.htmp2, Subsystem_Fy_B.tst);
      }

      i = L - 2;
    }
  }

  h[3] = 0.0;
  return info;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem__eigHermitianStandard(const real_T A[16],
  real_T V[16], real_T D[4])
{
  int32_T exitg1;
  boolean_T exitg2;
  if (Subsystem_Fy_anyNonFinite(A)) {
    for (Subsystem_Fy_B.knt = 0; Subsystem_Fy_B.knt < 16; Subsystem_Fy_B.knt++)
    {
      V[Subsystem_Fy_B.knt] = (rtNaN);
    }

    Subsystem_Fy_B.knt = 2;
    while (Subsystem_Fy_B.knt < 5) {
      V[Subsystem_Fy_B.knt - 1] = 0.0;
      Subsystem_Fy_B.knt++;
    }

    Subsystem_Fy_B.knt = 3;
    while (Subsystem_Fy_B.knt < 5) {
      V[Subsystem_Fy_B.knt + 3] = 0.0;
      Subsystem_Fy_B.knt++;
    }

    V[11] = 0.0;
    for (Subsystem_Fy_B.knt = 0; Subsystem_Fy_B.knt < 16; Subsystem_Fy_B.knt++)
    {
      Subsystem_Fy_B.a[Subsystem_Fy_B.knt] = (rtNaN);
    }
  } else {
    std::memcpy(&Subsystem_Fy_B.a[0], &A[0], sizeof(real_T) << 4U);
    Subsystem_Fy_B.work[0] = 0.0;
    Subsystem_Fy_B.work[1] = 0.0;
    Subsystem_Fy_B.work[2] = 0.0;
    Subsystem_Fy_B.work[3] = 0.0;
    Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.a[1];
    Subsystem_Fy_B.tau_idx_0 = 0.0;
    Subsystem_Fy_B.xnorm_tmp = Subsystem_Fy_xnrm2(2, Subsystem_Fy_B.a, 3);
    if (Subsystem_Fy_B.xnorm_tmp != 0.0) {
      Subsystem_Fy_B.beta1_l = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.a[1],
        Subsystem_Fy_B.xnorm_tmp);
      if (Subsystem_Fy_B.a[1] >= 0.0) {
        Subsystem_Fy_B.beta1_l = -Subsystem_Fy_B.beta1_l;
      }

      if (std::abs(Subsystem_Fy_B.beta1_l) < 1.0020841800044864E-292) {
        Subsystem_Fy_B.knt = -1;
        do {
          Subsystem_Fy_B.knt++;
          Subsystem_Fy_B.lastc = 3;
          while (Subsystem_Fy_B.lastc <= 4) {
            Subsystem_Fy_B.a[Subsystem_Fy_B.lastc - 1] *= 9.9792015476736E+291;
            Subsystem_Fy_B.lastc++;
          }

          Subsystem_Fy_B.beta1_l *= 9.9792015476736E+291;
          Subsystem_Fy_B.alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(Subsystem_Fy_B.beta1_l) >= 1.0020841800044864E-292));

        Subsystem_Fy_B.beta1_l = Subsystem_Fy_rt_hypotd_snf
          (Subsystem_Fy_B.alpha1, Subsystem_Fy_xnrm2(2, Subsystem_Fy_B.a, 3));
        if (Subsystem_Fy_B.alpha1 >= 0.0) {
          Subsystem_Fy_B.beta1_l = -Subsystem_Fy_B.beta1_l;
        }

        Subsystem_Fy_B.tau_idx_0 = (Subsystem_Fy_B.beta1_l -
          Subsystem_Fy_B.alpha1) / Subsystem_Fy_B.beta1_l;
        Subsystem_Fy_B.alpha1 = 1.0 / (Subsystem_Fy_B.alpha1 -
          Subsystem_Fy_B.beta1_l);
        Subsystem_Fy_B.lastc = 3;
        while (Subsystem_Fy_B.lastc <= 4) {
          Subsystem_Fy_B.a[Subsystem_Fy_B.lastc - 1] *= Subsystem_Fy_B.alpha1;
          Subsystem_Fy_B.lastc++;
        }

        Subsystem_Fy_B.lastc = 0;
        while (Subsystem_Fy_B.lastc <= Subsystem_Fy_B.knt) {
          Subsystem_Fy_B.beta1_l *= 1.0020841800044864E-292;
          Subsystem_Fy_B.lastc++;
        }

        Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.beta1_l;
      } else {
        Subsystem_Fy_B.tau_idx_0 = (Subsystem_Fy_B.beta1_l - Subsystem_Fy_B.a[1])
          / Subsystem_Fy_B.beta1_l;
        Subsystem_Fy_B.alpha1 = 1.0 / (Subsystem_Fy_B.a[1] -
          Subsystem_Fy_B.beta1_l);
        Subsystem_Fy_B.knt = 3;
        while (Subsystem_Fy_B.knt <= 4) {
          Subsystem_Fy_B.a[Subsystem_Fy_B.knt - 1] *= Subsystem_Fy_B.alpha1;
          Subsystem_Fy_B.knt++;
        }

        Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.beta1_l;
      }
    }

    Subsystem_Fy_B.a[1] = 1.0;
    if (Subsystem_Fy_B.tau_idx_0 != 0.0) {
      Subsystem_Fy_B.knt = 2;
      Subsystem_Fy_B.lastc = 3;
      while ((Subsystem_Fy_B.knt + 1 > 0) &&
             (Subsystem_Fy_B.a[Subsystem_Fy_B.lastc] == 0.0)) {
        Subsystem_Fy_B.knt--;
        Subsystem_Fy_B.lastc--;
      }

      Subsystem_Fy_B.lastc = 4;
      exitg2 = false;
      while ((!exitg2) && (Subsystem_Fy_B.lastc > 0)) {
        Subsystem_Fy_B.ix = Subsystem_Fy_B.lastc + 4;
        do {
          exitg1 = 0;
          if (Subsystem_Fy_B.ix <= (Subsystem_Fy_B.lastc + (Subsystem_Fy_B.knt <<
                2)) + 4) {
            if (Subsystem_Fy_B.a[Subsystem_Fy_B.ix - 1] != 0.0) {
              exitg1 = 1;
            } else {
              Subsystem_Fy_B.ix += 4;
            }
          } else {
            Subsystem_Fy_B.lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      Subsystem_Fy_B.knt = -1;
      Subsystem_Fy_B.lastc = 0;
    }

    if (Subsystem_Fy_B.knt + 1 > 0) {
      if (Subsystem_Fy_B.lastc != 0) {
        Subsystem_Fy_B.ix = 0;
        while (Subsystem_Fy_B.ix <= Subsystem_Fy_B.lastc - 1) {
          Subsystem_Fy_B.work[Subsystem_Fy_B.ix] = 0.0;
          Subsystem_Fy_B.ix++;
        }

        Subsystem_Fy_B.ix = 1;
        Subsystem_Fy_B.jy = (Subsystem_Fy_B.knt << 2) + 5;
        Subsystem_Fy_B.iac = 5;
        while (Subsystem_Fy_B.iac <= Subsystem_Fy_B.jy) {
          Subsystem_Fy_B.b_ix = 0;
          Subsystem_Fy_B.l = (Subsystem_Fy_B.iac + Subsystem_Fy_B.lastc) - 1;
          Subsystem_Fy_B.ia = Subsystem_Fy_B.iac;
          while (Subsystem_Fy_B.ia <= Subsystem_Fy_B.l) {
            Subsystem_Fy_B.work[Subsystem_Fy_B.b_ix] +=
              Subsystem_Fy_B.a[Subsystem_Fy_B.ia - 1] *
              Subsystem_Fy_B.a[Subsystem_Fy_B.ix];
            Subsystem_Fy_B.b_ix++;
            Subsystem_Fy_B.ia++;
          }

          Subsystem_Fy_B.ix++;
          Subsystem_Fy_B.iac += 4;
        }
      }

      if (!(-Subsystem_Fy_B.tau_idx_0 == 0.0)) {
        Subsystem_Fy_B.ix = 4;
        Subsystem_Fy_B.jy = 1;
        Subsystem_Fy_B.iac = 0;
        while (Subsystem_Fy_B.iac <= Subsystem_Fy_B.knt) {
          if (Subsystem_Fy_B.a[Subsystem_Fy_B.jy] != 0.0) {
            Subsystem_Fy_B.beta1_l = Subsystem_Fy_B.a[Subsystem_Fy_B.jy] *
              -Subsystem_Fy_B.tau_idx_0;
            Subsystem_Fy_B.b_ix = 0;
            Subsystem_Fy_B.l = Subsystem_Fy_B.lastc + Subsystem_Fy_B.ix;
            Subsystem_Fy_B.ia = Subsystem_Fy_B.ix;
            while (Subsystem_Fy_B.ia + 1 <= Subsystem_Fy_B.l) {
              Subsystem_Fy_B.a[Subsystem_Fy_B.ia] +=
                Subsystem_Fy_B.work[Subsystem_Fy_B.b_ix] *
                Subsystem_Fy_B.beta1_l;
              Subsystem_Fy_B.b_ix++;
              Subsystem_Fy_B.ia++;
            }
          }

          Subsystem_Fy_B.jy++;
          Subsystem_Fy_B.ix += 4;
          Subsystem_Fy_B.iac++;
        }
      }
    }

    Subsystem_Fy_xzlarf(3, 3, 2, Subsystem_Fy_B.tau_idx_0, Subsystem_Fy_B.a, 6,
                        Subsystem_Fy_B.work);
    Subsystem_Fy_B.a[1] = Subsystem_Fy_B.alpha1;
    Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.a[6];
    Subsystem_Fy_B.tau_idx_1 = 0.0;
    Subsystem_Fy_B.xnorm_tmp = Subsystem_Fy_xnrm2(1, Subsystem_Fy_B.a, 8);
    if (Subsystem_Fy_B.xnorm_tmp != 0.0) {
      Subsystem_Fy_B.beta1_l = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.a[6],
        Subsystem_Fy_B.xnorm_tmp);
      if (Subsystem_Fy_B.a[6] >= 0.0) {
        Subsystem_Fy_B.beta1_l = -Subsystem_Fy_B.beta1_l;
      }

      if (std::abs(Subsystem_Fy_B.beta1_l) < 1.0020841800044864E-292) {
        Subsystem_Fy_B.knt = -1;
        do {
          Subsystem_Fy_B.knt++;
          Subsystem_Fy_B.a[7] *= 9.9792015476736E+291;
          Subsystem_Fy_B.beta1_l *= 9.9792015476736E+291;
          Subsystem_Fy_B.alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(Subsystem_Fy_B.beta1_l) >= 1.0020841800044864E-292));

        Subsystem_Fy_B.beta1_l = Subsystem_Fy_rt_hypotd_snf
          (Subsystem_Fy_B.alpha1, Subsystem_Fy_xnrm2(1, Subsystem_Fy_B.a, 8));
        if (Subsystem_Fy_B.alpha1 >= 0.0) {
          Subsystem_Fy_B.beta1_l = -Subsystem_Fy_B.beta1_l;
        }

        Subsystem_Fy_B.tau_idx_1 = (Subsystem_Fy_B.beta1_l -
          Subsystem_Fy_B.alpha1) / Subsystem_Fy_B.beta1_l;
        Subsystem_Fy_B.a[7] *= 1.0 / (Subsystem_Fy_B.alpha1 -
          Subsystem_Fy_B.beta1_l);
        Subsystem_Fy_B.lastc = 0;
        while (Subsystem_Fy_B.lastc <= Subsystem_Fy_B.knt) {
          Subsystem_Fy_B.beta1_l *= 1.0020841800044864E-292;
          Subsystem_Fy_B.lastc++;
        }

        Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.beta1_l;
      } else {
        Subsystem_Fy_B.tau_idx_1 = (Subsystem_Fy_B.beta1_l - Subsystem_Fy_B.a[6])
          / Subsystem_Fy_B.beta1_l;
        Subsystem_Fy_B.a[7] *= 1.0 / (Subsystem_Fy_B.a[6] -
          Subsystem_Fy_B.beta1_l);
        Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.beta1_l;
      }
    }

    Subsystem_Fy_B.a[6] = 1.0;
    if (Subsystem_Fy_B.tau_idx_1 != 0.0) {
      Subsystem_Fy_B.knt = 1;
      Subsystem_Fy_B.lastc = 7;
      while ((Subsystem_Fy_B.knt + 1 > 0) &&
             (Subsystem_Fy_B.a[Subsystem_Fy_B.lastc] == 0.0)) {
        Subsystem_Fy_B.knt--;
        Subsystem_Fy_B.lastc--;
      }

      Subsystem_Fy_B.lastc = 4;
      exitg2 = false;
      while ((!exitg2) && (Subsystem_Fy_B.lastc > 0)) {
        Subsystem_Fy_B.ix = Subsystem_Fy_B.lastc + 8;
        do {
          exitg1 = 0;
          if (Subsystem_Fy_B.ix <= (Subsystem_Fy_B.lastc + (Subsystem_Fy_B.knt <<
                2)) + 8) {
            if (Subsystem_Fy_B.a[Subsystem_Fy_B.ix - 1] != 0.0) {
              exitg1 = 1;
            } else {
              Subsystem_Fy_B.ix += 4;
            }
          } else {
            Subsystem_Fy_B.lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      Subsystem_Fy_B.knt = -1;
      Subsystem_Fy_B.lastc = 0;
    }

    if (Subsystem_Fy_B.knt + 1 > 0) {
      if (Subsystem_Fy_B.lastc != 0) {
        Subsystem_Fy_B.ix = 0;
        while (Subsystem_Fy_B.ix <= Subsystem_Fy_B.lastc - 1) {
          Subsystem_Fy_B.work[Subsystem_Fy_B.ix] = 0.0;
          Subsystem_Fy_B.ix++;
        }

        Subsystem_Fy_B.ix = 6;
        Subsystem_Fy_B.jy = (Subsystem_Fy_B.knt << 2) + 9;
        Subsystem_Fy_B.iac = 9;
        while (Subsystem_Fy_B.iac <= Subsystem_Fy_B.jy) {
          Subsystem_Fy_B.b_ix = 0;
          Subsystem_Fy_B.l = (Subsystem_Fy_B.iac + Subsystem_Fy_B.lastc) - 1;
          Subsystem_Fy_B.ia = Subsystem_Fy_B.iac;
          while (Subsystem_Fy_B.ia <= Subsystem_Fy_B.l) {
            Subsystem_Fy_B.work[Subsystem_Fy_B.b_ix] +=
              Subsystem_Fy_B.a[Subsystem_Fy_B.ia - 1] *
              Subsystem_Fy_B.a[Subsystem_Fy_B.ix];
            Subsystem_Fy_B.b_ix++;
            Subsystem_Fy_B.ia++;
          }

          Subsystem_Fy_B.ix++;
          Subsystem_Fy_B.iac += 4;
        }
      }

      if (!(-Subsystem_Fy_B.tau_idx_1 == 0.0)) {
        Subsystem_Fy_B.ix = 8;
        Subsystem_Fy_B.jy = 6;
        Subsystem_Fy_B.iac = 0;
        while (Subsystem_Fy_B.iac <= Subsystem_Fy_B.knt) {
          if (Subsystem_Fy_B.a[Subsystem_Fy_B.jy] != 0.0) {
            Subsystem_Fy_B.beta1_l = Subsystem_Fy_B.a[Subsystem_Fy_B.jy] *
              -Subsystem_Fy_B.tau_idx_1;
            Subsystem_Fy_B.b_ix = 0;
            Subsystem_Fy_B.l = Subsystem_Fy_B.lastc + Subsystem_Fy_B.ix;
            Subsystem_Fy_B.ia = Subsystem_Fy_B.ix;
            while (Subsystem_Fy_B.ia + 1 <= Subsystem_Fy_B.l) {
              Subsystem_Fy_B.a[Subsystem_Fy_B.ia] +=
                Subsystem_Fy_B.work[Subsystem_Fy_B.b_ix] *
                Subsystem_Fy_B.beta1_l;
              Subsystem_Fy_B.b_ix++;
              Subsystem_Fy_B.ia++;
            }
          }

          Subsystem_Fy_B.jy++;
          Subsystem_Fy_B.ix += 4;
          Subsystem_Fy_B.iac++;
        }
      }
    }

    Subsystem_Fy_xzlarf(2, 2, 7, Subsystem_Fy_B.tau_idx_1, Subsystem_Fy_B.a, 11,
                        Subsystem_Fy_B.work);
    Subsystem_Fy_B.a[6] = Subsystem_Fy_B.alpha1;
    Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.a[11];
    Subsystem_Fy_B.xnorm_tmp = 0.0;
    Subsystem_Fy_B.xnorm_tmp_tmp = Subsystem_Fy_xnrm2(0, Subsystem_Fy_B.a, 12);
    if (Subsystem_Fy_B.xnorm_tmp_tmp != 0.0) {
      Subsystem_Fy_B.beta1_l = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.a[11],
        Subsystem_Fy_B.xnorm_tmp_tmp);
      if (Subsystem_Fy_B.a[11] >= 0.0) {
        Subsystem_Fy_B.beta1_l = -Subsystem_Fy_B.beta1_l;
      }

      if (std::abs(Subsystem_Fy_B.beta1_l) < 1.0020841800044864E-292) {
        Subsystem_Fy_B.knt = -1;
        do {
          Subsystem_Fy_B.knt++;
          Subsystem_Fy_B.beta1_l *= 9.9792015476736E+291;
          Subsystem_Fy_B.alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(Subsystem_Fy_B.beta1_l) >= 1.0020841800044864E-292));

        Subsystem_Fy_B.beta1_l = Subsystem_Fy_rt_hypotd_snf
          (Subsystem_Fy_B.alpha1, Subsystem_Fy_B.xnorm_tmp_tmp);
        if (Subsystem_Fy_B.alpha1 >= 0.0) {
          Subsystem_Fy_B.beta1_l = -Subsystem_Fy_B.beta1_l;
        }

        Subsystem_Fy_B.xnorm_tmp = (Subsystem_Fy_B.beta1_l -
          Subsystem_Fy_B.alpha1) / Subsystem_Fy_B.beta1_l;
        Subsystem_Fy_B.lastc = 0;
        while (Subsystem_Fy_B.lastc <= Subsystem_Fy_B.knt) {
          Subsystem_Fy_B.beta1_l *= 1.0020841800044864E-292;
          Subsystem_Fy_B.lastc++;
        }

        Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.beta1_l;
      } else {
        Subsystem_Fy_B.xnorm_tmp = (Subsystem_Fy_B.beta1_l - Subsystem_Fy_B.a[11])
          / Subsystem_Fy_B.beta1_l;
        Subsystem_Fy_B.alpha1 = Subsystem_Fy_B.beta1_l;
      }
    }

    Subsystem_Fy_B.a[11] = 1.0;
    if (Subsystem_Fy_B.xnorm_tmp != 0.0) {
      Subsystem_Fy_B.knt = 0;
      Subsystem_Fy_B.lastc = 11;
      while ((Subsystem_Fy_B.knt + 1 > 0) &&
             (Subsystem_Fy_B.a[Subsystem_Fy_B.lastc] == 0.0)) {
        Subsystem_Fy_B.knt--;
        Subsystem_Fy_B.lastc--;
      }

      Subsystem_Fy_B.lastc = 4;
      exitg2 = false;
      while ((!exitg2) && (Subsystem_Fy_B.lastc > 0)) {
        Subsystem_Fy_B.ix = Subsystem_Fy_B.lastc + 12;
        do {
          exitg1 = 0;
          if (Subsystem_Fy_B.ix <= (Subsystem_Fy_B.lastc + (Subsystem_Fy_B.knt <<
                2)) + 12) {
            if (Subsystem_Fy_B.a[Subsystem_Fy_B.ix - 1] != 0.0) {
              exitg1 = 1;
            } else {
              Subsystem_Fy_B.ix += 4;
            }
          } else {
            Subsystem_Fy_B.lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      Subsystem_Fy_B.knt = -1;
      Subsystem_Fy_B.lastc = 0;
    }

    if (Subsystem_Fy_B.knt + 1 > 0) {
      if (Subsystem_Fy_B.lastc != 0) {
        Subsystem_Fy_B.ix = 0;
        while (Subsystem_Fy_B.ix <= Subsystem_Fy_B.lastc - 1) {
          Subsystem_Fy_B.work[Subsystem_Fy_B.ix] = 0.0;
          Subsystem_Fy_B.ix++;
        }

        Subsystem_Fy_B.ix = 11;
        Subsystem_Fy_B.jy = (Subsystem_Fy_B.knt << 2) + 13;
        Subsystem_Fy_B.iac = 13;
        while (Subsystem_Fy_B.iac <= Subsystem_Fy_B.jy) {
          Subsystem_Fy_B.b_ix = 0;
          Subsystem_Fy_B.l = (Subsystem_Fy_B.iac + Subsystem_Fy_B.lastc) - 1;
          Subsystem_Fy_B.ia = Subsystem_Fy_B.iac;
          while (Subsystem_Fy_B.ia <= Subsystem_Fy_B.l) {
            Subsystem_Fy_B.work[Subsystem_Fy_B.b_ix] +=
              Subsystem_Fy_B.a[Subsystem_Fy_B.ia - 1] *
              Subsystem_Fy_B.a[Subsystem_Fy_B.ix];
            Subsystem_Fy_B.b_ix++;
            Subsystem_Fy_B.ia++;
          }

          Subsystem_Fy_B.ix++;
          Subsystem_Fy_B.iac += 4;
        }
      }

      if (!(-Subsystem_Fy_B.xnorm_tmp == 0.0)) {
        Subsystem_Fy_B.ix = 12;
        Subsystem_Fy_B.jy = 11;
        Subsystem_Fy_B.iac = 0;
        while (Subsystem_Fy_B.iac <= Subsystem_Fy_B.knt) {
          if (Subsystem_Fy_B.a[Subsystem_Fy_B.jy] != 0.0) {
            Subsystem_Fy_B.beta1_l = Subsystem_Fy_B.a[Subsystem_Fy_B.jy] *
              -Subsystem_Fy_B.xnorm_tmp;
            Subsystem_Fy_B.b_ix = 0;
            Subsystem_Fy_B.l = Subsystem_Fy_B.lastc + Subsystem_Fy_B.ix;
            Subsystem_Fy_B.ia = Subsystem_Fy_B.ix;
            while (Subsystem_Fy_B.ia + 1 <= Subsystem_Fy_B.l) {
              Subsystem_Fy_B.a[Subsystem_Fy_B.ia] +=
                Subsystem_Fy_B.work[Subsystem_Fy_B.b_ix] *
                Subsystem_Fy_B.beta1_l;
              Subsystem_Fy_B.b_ix++;
              Subsystem_Fy_B.ia++;
            }
          }

          Subsystem_Fy_B.jy++;
          Subsystem_Fy_B.ix += 4;
          Subsystem_Fy_B.iac++;
        }
      }
    }

    Subsystem_Fy_xzlarf(1, 1, 12, Subsystem_Fy_B.xnorm_tmp, Subsystem_Fy_B.a, 16,
                        Subsystem_Fy_B.work);
    Subsystem_Fy_B.a[11] = Subsystem_Fy_B.alpha1;
    std::memcpy(&V[0], &Subsystem_Fy_B.a[0], sizeof(real_T) << 4U);
    Subsystem_Fy_B.knt = 0;
    while (Subsystem_Fy_B.knt <= 2) {
      V[Subsystem_Fy_B.knt + 12] = 0.0;
      Subsystem_Fy_B.knt++;
    }

    Subsystem_Fy_B.knt = 0;
    while (Subsystem_Fy_B.knt <= 1) {
      V[Subsystem_Fy_B.knt + 8] = 0.0;
      Subsystem_Fy_B.knt++;
    }

    Subsystem_Fy_B.knt = 1;
    while (Subsystem_Fy_B.knt + 3 < 5) {
      V[Subsystem_Fy_B.knt + 10] = V[Subsystem_Fy_B.knt + 6];
      Subsystem_Fy_B.knt++;
    }

    V[4] = 0.0;
    Subsystem_Fy_B.knt = 0;
    while (Subsystem_Fy_B.knt + 3 < 5) {
      V[Subsystem_Fy_B.knt + 6] = V[Subsystem_Fy_B.knt + 2];
      Subsystem_Fy_B.knt++;
    }

    Subsystem_Fy_B.work[0] = 0.0;
    V[1] = 0.0;
    Subsystem_Fy_B.work[1] = 0.0;
    V[2] = 0.0;
    Subsystem_Fy_B.work[2] = 0.0;
    V[3] = 0.0;
    Subsystem_Fy_B.work[3] = 0.0;
    V[0] = 1.0;
    V[15] = 1.0 - Subsystem_Fy_B.xnorm_tmp;
    Subsystem_Fy_B.knt = 0;
    while (Subsystem_Fy_B.knt <= 1) {
      V[14 - Subsystem_Fy_B.knt] = 0.0;
      Subsystem_Fy_B.knt++;
    }

    V[10] = 1.0;
    Subsystem_Fy_xzlarf(2, 1, 11, Subsystem_Fy_B.tau_idx_1, V, 15,
                        Subsystem_Fy_B.work);
    Subsystem_Fy_B.knt = 11;
    while (Subsystem_Fy_B.knt + 1 <= 12) {
      V[Subsystem_Fy_B.knt] *= -Subsystem_Fy_B.tau_idx_1;
      Subsystem_Fy_B.knt++;
    }

    V[10] = 1.0 - Subsystem_Fy_B.tau_idx_1;
    V[9] = 0.0;
    V[5] = 1.0;
    Subsystem_Fy_xzlarf(3, 2, 6, Subsystem_Fy_B.tau_idx_0, V, 10,
                        Subsystem_Fy_B.work);
    Subsystem_Fy_B.knt = 6;
    while (Subsystem_Fy_B.knt + 1 <= 8) {
      V[Subsystem_Fy_B.knt] *= -Subsystem_Fy_B.tau_idx_0;
      Subsystem_Fy_B.knt++;
    }

    V[5] = 1.0 - Subsystem_Fy_B.tau_idx_0;
    Subsystem_Fy_xhseqr(Subsystem_Fy_B.a, V);
  }

  D[0] = Subsystem_Fy_B.a[0];
  D[1] = Subsystem_Fy_B.a[5];
  D[2] = Subsystem_Fy_B.a[10];
  D[3] = Subsystem_Fy_B.a[15];
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_eig(const real_T A[16], creal_T V[16],
  creal_T D[4])
{
  int32_T exitg1;
  boolean_T exitg2;
  boolean_T guard1 = false;
  boolean_T ilascl;
  boolean_T notdone;
  if (Subsystem_Fy_anyNonFinite(A)) {
    for (Subsystem_Fy_B.stemp_re_tmp = 0; Subsystem_Fy_B.stemp_re_tmp < 16;
         Subsystem_Fy_B.stemp_re_tmp++) {
      V[Subsystem_Fy_B.stemp_re_tmp].re = (rtNaN);
      V[Subsystem_Fy_B.stemp_re_tmp].im = 0.0;
    }

    D[0].re = (rtNaN);
    D[0].im = 0.0;
    D[1].re = (rtNaN);
    D[1].im = 0.0;
    D[2].re = (rtNaN);
    D[2].im = 0.0;
    D[3].re = (rtNaN);
    D[3].im = 0.0;
  } else {
    ilascl = true;
    Subsystem_Fy_B.j = 0;
    exitg2 = false;
    while ((!exitg2) && (Subsystem_Fy_B.j < 4)) {
      Subsystem_Fy_B.d_i = 0;
      do {
        exitg1 = 0;
        if (Subsystem_Fy_B.d_i <= Subsystem_Fy_B.j) {
          if (!(A[(Subsystem_Fy_B.j << 2) + Subsystem_Fy_B.d_i] == A
                [(Subsystem_Fy_B.d_i << 2) + Subsystem_Fy_B.j])) {
            ilascl = false;
            exitg1 = 1;
          } else {
            Subsystem_Fy_B.d_i++;
          }
        } else {
          Subsystem_Fy_B.j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }

    if (ilascl) {
      Subsystem__eigHermitianStandard(A, Subsystem_Fy_B.b_I, Subsystem_Fy_B.b_D);
      for (Subsystem_Fy_B.stemp_re_tmp = 0; Subsystem_Fy_B.stemp_re_tmp < 16;
           Subsystem_Fy_B.stemp_re_tmp++) {
        V[Subsystem_Fy_B.stemp_re_tmp].re =
          Subsystem_Fy_B.b_I[Subsystem_Fy_B.stemp_re_tmp];
        V[Subsystem_Fy_B.stemp_re_tmp].im = 0.0;
      }

      D[0].re = Subsystem_Fy_B.b_D[0];
      D[0].im = 0.0;
      D[1].re = Subsystem_Fy_B.b_D[1];
      D[1].im = 0.0;
      D[2].re = Subsystem_Fy_B.b_D[2];
      D[2].im = 0.0;
      D[3].re = Subsystem_Fy_B.b_D[3];
      D[3].im = 0.0;
    } else {
      for (Subsystem_Fy_B.stemp_re_tmp = 0; Subsystem_Fy_B.stemp_re_tmp < 16;
           Subsystem_Fy_B.stemp_re_tmp++) {
        Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re =
          A[Subsystem_Fy_B.stemp_re_tmp];
        Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im = 0.0;
      }

      Subsystem_Fy_B.anrm = 0.0;
      Subsystem_Fy_B.j = 0;
      exitg2 = false;
      while ((!exitg2) && (Subsystem_Fy_B.j < 16)) {
        Subsystem_Fy_B.b_absxk = Subsystem_Fy_rt_hypotd_snf
          (Subsystem_Fy_B.At[Subsystem_Fy_B.j].re,
           Subsystem_Fy_B.At[Subsystem_Fy_B.j].im);
        if (rtIsNaN(Subsystem_Fy_B.b_absxk)) {
          Subsystem_Fy_B.anrm = (rtNaN);
          exitg2 = true;
        } else {
          if (Subsystem_Fy_B.b_absxk > Subsystem_Fy_B.anrm) {
            Subsystem_Fy_B.anrm = Subsystem_Fy_B.b_absxk;
          }

          Subsystem_Fy_B.j++;
        }
      }

      if (rtIsInf(Subsystem_Fy_B.anrm) || rtIsNaN(Subsystem_Fy_B.anrm)) {
        D[0].re = (rtNaN);
        D[0].im = 0.0;
        Subsystem_Fy_B.beta1[0].re = (rtNaN);
        Subsystem_Fy_B.beta1[0].im = 0.0;
        D[1].re = (rtNaN);
        D[1].im = 0.0;
        Subsystem_Fy_B.beta1[1].re = (rtNaN);
        Subsystem_Fy_B.beta1[1].im = 0.0;
        D[2].re = (rtNaN);
        D[2].im = 0.0;
        Subsystem_Fy_B.beta1[2].re = (rtNaN);
        Subsystem_Fy_B.beta1[2].im = 0.0;
        D[3].re = (rtNaN);
        D[3].im = 0.0;
        Subsystem_Fy_B.beta1[3].re = (rtNaN);
        Subsystem_Fy_B.beta1[3].im = 0.0;
        for (Subsystem_Fy_B.stemp_re_tmp = 0; Subsystem_Fy_B.stemp_re_tmp < 16;
             Subsystem_Fy_B.stemp_re_tmp++) {
          V[Subsystem_Fy_B.stemp_re_tmp].re = (rtNaN);
          V[Subsystem_Fy_B.stemp_re_tmp].im = 0.0;
        }
      } else {
        ilascl = false;
        Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.anrm;
        guard1 = false;
        if ((Subsystem_Fy_B.anrm > 0.0) && (Subsystem_Fy_B.anrm <
             6.7178761075670888E-139)) {
          Subsystem_Fy_B.b_absxk = 6.7178761075670888E-139;
          ilascl = true;
          guard1 = true;
        } else if (Subsystem_Fy_B.anrm > 1.4885657073574029E+138) {
          Subsystem_Fy_B.b_absxk = 1.4885657073574029E+138;
          ilascl = true;
          guard1 = true;
        }

        if (guard1) {
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.anrm;
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk;
          notdone = true;
          while (notdone) {
            Subsystem_Fy_B.stemp_im = Subsystem_Fy_B.cfromc *
              2.0041683600089728E-292;
            Subsystem_Fy_B.cto1 = Subsystem_Fy_B.ctoc / 4.9896007738368E+291;
            if ((Subsystem_Fy_B.stemp_im > Subsystem_Fy_B.ctoc) &&
                (Subsystem_Fy_B.ctoc != 0.0)) {
              Subsystem_Fy_B.mul = 2.0041683600089728E-292;
              Subsystem_Fy_B.cfromc = Subsystem_Fy_B.stemp_im;
            } else if (Subsystem_Fy_B.cto1 > Subsystem_Fy_B.cfromc) {
              Subsystem_Fy_B.mul = 4.9896007738368E+291;
              Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cto1;
            } else {
              Subsystem_Fy_B.mul = Subsystem_Fy_B.ctoc / Subsystem_Fy_B.cfromc;
              notdone = false;
            }

            for (Subsystem_Fy_B.stemp_re_tmp = 0; Subsystem_Fy_B.stemp_re_tmp <
                 16; Subsystem_Fy_B.stemp_re_tmp++) {
              Subsystem_Fy_B.s = Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp];
              Subsystem_Fy_B.s.re *= Subsystem_Fy_B.mul;
              Subsystem_Fy_B.s.im *= Subsystem_Fy_B.mul;
              Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp] = Subsystem_Fy_B.s;
            }
          }
        }

        Subsystem_Fy_xzggbal(Subsystem_Fy_B.At, &Subsystem_Fy_B.d_i,
                             &Subsystem_Fy_B.j, Subsystem_Fy_B.rscale);
        std::memset(&Subsystem_Fy_B.b_I[0], 0, sizeof(real_T) << 4U);
        Subsystem_Fy_B.b_I[0] = 1.0;
        Subsystem_Fy_B.b_I[5] = 1.0;
        Subsystem_Fy_B.b_I[10] = 1.0;
        Subsystem_Fy_B.b_I[15] = 1.0;
        for (Subsystem_Fy_B.stemp_re_tmp = 0; Subsystem_Fy_B.stemp_re_tmp < 16;
             Subsystem_Fy_B.stemp_re_tmp++) {
          V[Subsystem_Fy_B.stemp_re_tmp].re =
            Subsystem_Fy_B.b_I[Subsystem_Fy_B.stemp_re_tmp];
          V[Subsystem_Fy_B.stemp_re_tmp].im = 0.0;
        }

        if (Subsystem_Fy_B.j >= Subsystem_Fy_B.d_i + 2) {
          Subsystem_Fy_B.jcol = Subsystem_Fy_B.d_i - 1;
          while (Subsystem_Fy_B.jcol + 1 < Subsystem_Fy_B.j - 1) {
            Subsystem_Fy_B.jrow = Subsystem_Fy_B.j - 1;
            while (Subsystem_Fy_B.jrow + 1 > Subsystem_Fy_B.jcol + 2) {
              Subsystem_Fy_xzlartg(Subsystem_Fy_B.At[(Subsystem_Fy_B.jrow +
                (Subsystem_Fy_B.jcol << 2)) - 1],
                                   Subsystem_Fy_B.At[Subsystem_Fy_B.jrow +
                                   (Subsystem_Fy_B.jcol << 2)],
                                   &Subsystem_Fy_B.cfromc, &Subsystem_Fy_B.s,
                                   &Subsystem_Fy_B.At[(Subsystem_Fy_B.jrow +
                (Subsystem_Fy_B.jcol << 2)) - 1]);
              Subsystem_Fy_B.stemp_re_tmp = Subsystem_Fy_B.jrow +
                (Subsystem_Fy_B.jcol << 2);
              Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re = 0.0;
              Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im = 0.0;
              Subsystem_Fy_B.c_j = Subsystem_Fy_B.jcol + 1;
              while (Subsystem_Fy_B.c_j + 1 < 5) {
                Subsystem_Fy_B.stemp_re_tmp = (Subsystem_Fy_B.c_j << 2) +
                  Subsystem_Fy_B.jrow;
                Subsystem_Fy_B.ctoc =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].re *
                  Subsystem_Fy_B.cfromc +
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re *
                   Subsystem_Fy_B.s.re -
                   Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im *
                   Subsystem_Fy_B.s.im);
                Subsystem_Fy_B.stemp_im =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].im *
                  Subsystem_Fy_B.cfromc +
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im *
                   Subsystem_Fy_B.s.re +
                   Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re *
                   Subsystem_Fy_B.s.im);
                Subsystem_Fy_B.cto1 =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].re;
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re *
                  Subsystem_Fy_B.cfromc -
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].re *
                   Subsystem_Fy_B.s.re +
                   Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].im *
                   Subsystem_Fy_B.s.im);
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im *
                  Subsystem_Fy_B.cfromc -
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].im *
                   Subsystem_Fy_B.s.re - Subsystem_Fy_B.s.im *
                   Subsystem_Fy_B.cto1);
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].re =
                  Subsystem_Fy_B.ctoc;
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp - 1].im =
                  Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.c_j++;
              }

              Subsystem_Fy_B.s.re = -Subsystem_Fy_B.s.re;
              Subsystem_Fy_B.s.im = -Subsystem_Fy_B.s.im;
              Subsystem_Fy_B.c_j = 0;
              while (Subsystem_Fy_B.c_j + 1 <= Subsystem_Fy_B.j) {
                Subsystem_Fy_B.stemp_re_tmp = ((Subsystem_Fy_B.jrow - 1) << 2) +
                  Subsystem_Fy_B.c_j;
                Subsystem_Fy_B.stemp_re_tmp_p = (Subsystem_Fy_B.jrow << 2) +
                  Subsystem_Fy_B.c_j;
                Subsystem_Fy_B.ctoc =
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re *
                   Subsystem_Fy_B.s.re -
                   Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im *
                   Subsystem_Fy_B.s.im) +
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].re *
                  Subsystem_Fy_B.cfromc;
                Subsystem_Fy_B.stemp_im =
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im *
                   Subsystem_Fy_B.s.re +
                   Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re *
                   Subsystem_Fy_B.s.im) +
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].im *
                  Subsystem_Fy_B.cfromc;
                Subsystem_Fy_B.cto1 =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].re;
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].re *
                  Subsystem_Fy_B.cfromc -
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].re *
                   Subsystem_Fy_B.s.re +
                   Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].im *
                   Subsystem_Fy_B.s.im);
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im =
                  Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp].im *
                  Subsystem_Fy_B.cfromc -
                  (Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].im *
                   Subsystem_Fy_B.s.re - Subsystem_Fy_B.s.im *
                   Subsystem_Fy_B.cto1);
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].re =
                  Subsystem_Fy_B.ctoc;
                Subsystem_Fy_B.At[Subsystem_Fy_B.stemp_re_tmp_p].im =
                  Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.c_j++;
              }

              Subsystem_Fy_B.stemp_re_tmp = (Subsystem_Fy_B.jrow - 1) << 2;
              Subsystem_Fy_B.stemp_re_tmp_p = Subsystem_Fy_B.jrow << 2;
              Subsystem_Fy_B.ctoc = (V[Subsystem_Fy_B.stemp_re_tmp].re *
                Subsystem_Fy_B.s.re - V[Subsystem_Fy_B.stemp_re_tmp].im *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p].re *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.stemp_im = (V[Subsystem_Fy_B.stemp_re_tmp].im *
                Subsystem_Fy_B.s.re + V[Subsystem_Fy_B.stemp_re_tmp].re *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p].im *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.cto1 = V[Subsystem_Fy_B.stemp_re_tmp_p].re;
              V[Subsystem_Fy_B.stemp_re_tmp].re = V[Subsystem_Fy_B.stemp_re_tmp]
                .re * Subsystem_Fy_B.cfromc - (V[Subsystem_Fy_B.stemp_re_tmp_p].
                re * Subsystem_Fy_B.s.re + V[Subsystem_Fy_B.stemp_re_tmp_p].im *
                Subsystem_Fy_B.s.im);
              V[Subsystem_Fy_B.stemp_re_tmp].im = V[Subsystem_Fy_B.stemp_re_tmp]
                .im * Subsystem_Fy_B.cfromc - (V[Subsystem_Fy_B.stemp_re_tmp_p].
                im * Subsystem_Fy_B.s.re - Subsystem_Fy_B.s.im *
                Subsystem_Fy_B.cto1);
              V[Subsystem_Fy_B.stemp_re_tmp_p].re = Subsystem_Fy_B.ctoc;
              V[Subsystem_Fy_B.stemp_re_tmp_p].im = Subsystem_Fy_B.stemp_im;
              Subsystem_Fy_B.ctoc = (V[Subsystem_Fy_B.stemp_re_tmp + 1].re *
                Subsystem_Fy_B.s.re - V[Subsystem_Fy_B.stemp_re_tmp + 1].im *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p + 1].re *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.stemp_im = (V[Subsystem_Fy_B.stemp_re_tmp + 1].im *
                Subsystem_Fy_B.s.re + V[Subsystem_Fy_B.stemp_re_tmp + 1].re *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p + 1].im *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.cto1 = V[Subsystem_Fy_B.stemp_re_tmp_p + 1].re;
              V[Subsystem_Fy_B.stemp_re_tmp + 1].re =
                V[Subsystem_Fy_B.stemp_re_tmp + 1].re * Subsystem_Fy_B.cfromc -
                (V[Subsystem_Fy_B.stemp_re_tmp_p + 1].re * Subsystem_Fy_B.s.re +
                 V[Subsystem_Fy_B.stemp_re_tmp_p + 1].im * Subsystem_Fy_B.s.im);
              V[Subsystem_Fy_B.stemp_re_tmp + 1].im =
                V[Subsystem_Fy_B.stemp_re_tmp + 1].im * Subsystem_Fy_B.cfromc -
                (V[Subsystem_Fy_B.stemp_re_tmp_p + 1].im * Subsystem_Fy_B.s.re -
                 Subsystem_Fy_B.s.im * Subsystem_Fy_B.cto1);
              V[Subsystem_Fy_B.stemp_re_tmp_p + 1].re = Subsystem_Fy_B.ctoc;
              V[Subsystem_Fy_B.stemp_re_tmp_p + 1].im = Subsystem_Fy_B.stemp_im;
              Subsystem_Fy_B.ctoc = (V[Subsystem_Fy_B.stemp_re_tmp + 2].re *
                Subsystem_Fy_B.s.re - V[Subsystem_Fy_B.stemp_re_tmp + 2].im *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p + 2].re *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.stemp_im = (V[Subsystem_Fy_B.stemp_re_tmp + 2].im *
                Subsystem_Fy_B.s.re + V[Subsystem_Fy_B.stemp_re_tmp + 2].re *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p + 2].im *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.cto1 = V[Subsystem_Fy_B.stemp_re_tmp_p + 2].re;
              V[Subsystem_Fy_B.stemp_re_tmp + 2].re =
                V[Subsystem_Fy_B.stemp_re_tmp + 2].re * Subsystem_Fy_B.cfromc -
                (V[Subsystem_Fy_B.stemp_re_tmp_p + 2].re * Subsystem_Fy_B.s.re +
                 V[Subsystem_Fy_B.stemp_re_tmp_p + 2].im * Subsystem_Fy_B.s.im);
              V[Subsystem_Fy_B.stemp_re_tmp + 2].im =
                V[Subsystem_Fy_B.stemp_re_tmp + 2].im * Subsystem_Fy_B.cfromc -
                (V[Subsystem_Fy_B.stemp_re_tmp_p + 2].im * Subsystem_Fy_B.s.re -
                 Subsystem_Fy_B.s.im * Subsystem_Fy_B.cto1);
              V[Subsystem_Fy_B.stemp_re_tmp_p + 2].re = Subsystem_Fy_B.ctoc;
              V[Subsystem_Fy_B.stemp_re_tmp_p + 2].im = Subsystem_Fy_B.stemp_im;
              Subsystem_Fy_B.ctoc = (V[Subsystem_Fy_B.stemp_re_tmp + 3].re *
                Subsystem_Fy_B.s.re - V[Subsystem_Fy_B.stemp_re_tmp + 3].im *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p + 3].re *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.stemp_im = (V[Subsystem_Fy_B.stemp_re_tmp + 3].im *
                Subsystem_Fy_B.s.re + V[Subsystem_Fy_B.stemp_re_tmp + 3].re *
                Subsystem_Fy_B.s.im) + V[Subsystem_Fy_B.stemp_re_tmp_p + 3].im *
                Subsystem_Fy_B.cfromc;
              Subsystem_Fy_B.cto1 = V[Subsystem_Fy_B.stemp_re_tmp_p + 3].re;
              V[Subsystem_Fy_B.stemp_re_tmp + 3].re =
                V[Subsystem_Fy_B.stemp_re_tmp + 3].re * Subsystem_Fy_B.cfromc -
                (V[Subsystem_Fy_B.stemp_re_tmp_p + 3].re * Subsystem_Fy_B.s.re +
                 V[Subsystem_Fy_B.stemp_re_tmp_p + 3].im * Subsystem_Fy_B.s.im);
              V[Subsystem_Fy_B.stemp_re_tmp + 3].im =
                V[Subsystem_Fy_B.stemp_re_tmp + 3].im * Subsystem_Fy_B.cfromc -
                (V[Subsystem_Fy_B.stemp_re_tmp_p + 3].im * Subsystem_Fy_B.s.re -
                 Subsystem_Fy_B.s.im * Subsystem_Fy_B.cto1);
              V[Subsystem_Fy_B.stemp_re_tmp_p + 3].re = Subsystem_Fy_B.ctoc;
              V[Subsystem_Fy_B.stemp_re_tmp_p + 3].im = Subsystem_Fy_B.stemp_im;
              Subsystem_Fy_B.jrow--;
            }

            Subsystem_Fy_B.jcol++;
          }
        }

        Subsystem_Fy_xzhgeqz(Subsystem_Fy_B.At, Subsystem_Fy_B.d_i,
                             Subsystem_Fy_B.j, V, &Subsystem_Fy_B.jcol, D,
                             Subsystem_Fy_B.beta1);
        if (Subsystem_Fy_B.jcol == 0) {
          Subsystem_Fy_xztgevc(Subsystem_Fy_B.At, V);
          if (Subsystem_Fy_B.d_i > 1) {
            Subsystem_Fy_B.d_i -= 2;
            while (Subsystem_Fy_B.d_i + 1 >= 1) {
              Subsystem_Fy_B.jcol = Subsystem_Fy_B.rscale[Subsystem_Fy_B.d_i] -
                1;
              if (Subsystem_Fy_B.d_i + 1 !=
                  Subsystem_Fy_B.rscale[Subsystem_Fy_B.d_i]) {
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.d_i].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.d_i].im;
                V[Subsystem_Fy_B.d_i] = V[Subsystem_Fy_B.jcol];
                V[Subsystem_Fy_B.jcol].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol].im = Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.d_i + 4].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.d_i + 4].im;
                V[Subsystem_Fy_B.d_i + 4] = V[Subsystem_Fy_B.jcol + 4];
                V[Subsystem_Fy_B.jcol + 4].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol + 4].im = Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.d_i + 8].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.d_i + 8].im;
                V[Subsystem_Fy_B.d_i + 8] = V[Subsystem_Fy_B.jcol + 8];
                V[Subsystem_Fy_B.jcol + 8].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol + 8].im = Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.d_i + 12].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.d_i + 12].im;
                V[Subsystem_Fy_B.d_i + 12] = V[Subsystem_Fy_B.jcol + 12];
                V[Subsystem_Fy_B.jcol + 12].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol + 12].im = Subsystem_Fy_B.stemp_im;
              }

              Subsystem_Fy_B.d_i--;
            }
          }

          if (Subsystem_Fy_B.j < 4) {
            while (Subsystem_Fy_B.j + 1 < 5) {
              Subsystem_Fy_B.jcol = Subsystem_Fy_B.rscale[Subsystem_Fy_B.j] - 1;
              if (Subsystem_Fy_B.j + 1 != Subsystem_Fy_B.rscale[Subsystem_Fy_B.j])
              {
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.j].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.j].im;
                V[Subsystem_Fy_B.j] = V[Subsystem_Fy_B.jcol];
                V[Subsystem_Fy_B.jcol].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol].im = Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.j + 4].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.j + 4].im;
                V[Subsystem_Fy_B.j + 4] = V[Subsystem_Fy_B.jcol + 4];
                V[Subsystem_Fy_B.jcol + 4].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol + 4].im = Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.j + 8].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.j + 8].im;
                V[Subsystem_Fy_B.j + 8] = V[Subsystem_Fy_B.jcol + 8];
                V[Subsystem_Fy_B.jcol + 8].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol + 8].im = Subsystem_Fy_B.stemp_im;
                Subsystem_Fy_B.ctoc = V[Subsystem_Fy_B.j + 12].re;
                Subsystem_Fy_B.stemp_im = V[Subsystem_Fy_B.j + 12].im;
                V[Subsystem_Fy_B.j + 12] = V[Subsystem_Fy_B.jcol + 12];
                V[Subsystem_Fy_B.jcol + 12].re = Subsystem_Fy_B.ctoc;
                V[Subsystem_Fy_B.jcol + 12].im = Subsystem_Fy_B.stemp_im;
              }

              Subsystem_Fy_B.j++;
            }
          }

          for (Subsystem_Fy_B.j = 0; Subsystem_Fy_B.j < 4; Subsystem_Fy_B.j++) {
            Subsystem_Fy_B.d_i = Subsystem_Fy_B.j << 2;
            Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.d_i].re) + std::
              abs(V[Subsystem_Fy_B.d_i].im);
            Subsystem_Fy_B.ctoc = std::abs(V[Subsystem_Fy_B.d_i + 1].re) + std::
              abs(V[Subsystem_Fy_B.d_i + 1].im);
            if (Subsystem_Fy_B.ctoc > Subsystem_Fy_B.cfromc) {
              Subsystem_Fy_B.cfromc = Subsystem_Fy_B.ctoc;
            }

            Subsystem_Fy_B.ctoc = std::abs(V[Subsystem_Fy_B.d_i + 2].re) + std::
              abs(V[Subsystem_Fy_B.d_i + 2].im);
            if (Subsystem_Fy_B.ctoc > Subsystem_Fy_B.cfromc) {
              Subsystem_Fy_B.cfromc = Subsystem_Fy_B.ctoc;
            }

            Subsystem_Fy_B.ctoc = std::abs(V[Subsystem_Fy_B.d_i + 3].re) + std::
              abs(V[Subsystem_Fy_B.d_i + 3].im);
            if (Subsystem_Fy_B.ctoc > Subsystem_Fy_B.cfromc) {
              Subsystem_Fy_B.cfromc = Subsystem_Fy_B.ctoc;
            }

            if (Subsystem_Fy_B.cfromc >= 6.7178761075670888E-139) {
              Subsystem_Fy_B.cfromc = 1.0 / Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i].re *= Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i].im *= Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i + 1].re *= Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i + 1].im *= Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i + 2].re *= Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i + 2].im *= Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i + 3].re *= Subsystem_Fy_B.cfromc;
              V[Subsystem_Fy_B.d_i + 3].im *= Subsystem_Fy_B.cfromc;
            }
          }

          if (ilascl) {
            ilascl = true;
            while (ilascl) {
              Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk *
                2.0041683600089728E-292;
              Subsystem_Fy_B.ctoc = Subsystem_Fy_B.anrm / 4.9896007738368E+291;
              if ((Subsystem_Fy_B.cfromc > Subsystem_Fy_B.anrm) &&
                  (Subsystem_Fy_B.anrm != 0.0)) {
                Subsystem_Fy_B.stemp_im = 2.0041683600089728E-292;
                Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
              } else if (Subsystem_Fy_B.ctoc > Subsystem_Fy_B.b_absxk) {
                Subsystem_Fy_B.stemp_im = 4.9896007738368E+291;
                Subsystem_Fy_B.anrm = Subsystem_Fy_B.ctoc;
              } else {
                Subsystem_Fy_B.stemp_im = Subsystem_Fy_B.anrm /
                  Subsystem_Fy_B.b_absxk;
                ilascl = false;
              }

              D[0].re *= Subsystem_Fy_B.stemp_im;
              D[0].im *= Subsystem_Fy_B.stemp_im;
              D[1].re *= Subsystem_Fy_B.stemp_im;
              D[1].im *= Subsystem_Fy_B.stemp_im;
              D[2].re *= Subsystem_Fy_B.stemp_im;
              D[2].im *= Subsystem_Fy_B.stemp_im;
              D[3].re *= Subsystem_Fy_B.stemp_im;
              D[3].im *= Subsystem_Fy_B.stemp_im;
            }
          }
        }
      }

      Subsystem_Fy_B.anrm = 0.0;
      Subsystem_Fy_B.b_absxk = 3.3121686421112381E-170;
      Subsystem_Fy_B.j = 0;
      while (Subsystem_Fy_B.j + 1 <= 4) {
        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].re);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].im);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.j++;
      }

      Subsystem_Fy_B.anrm = Subsystem_Fy_B.b_absxk * std::sqrt
        (Subsystem_Fy_B.anrm);
      Subsystem_Fy_B.j = 0;
      while (Subsystem_Fy_B.j + 1 <= 4) {
        Subsystem_Fy_B.b_absxk = V[Subsystem_Fy_B.j].re;
        Subsystem_Fy_B.cfromc = V[Subsystem_Fy_B.j].im;
        if (Subsystem_Fy_B.cfromc == 0.0) {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = 0.0;
        } else if (Subsystem_Fy_B.b_absxk == 0.0) {
          V[Subsystem_Fy_B.j].re = 0.0;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        } else {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        }

        Subsystem_Fy_B.j++;
      }

      Subsystem_Fy_B.anrm = 0.0;
      Subsystem_Fy_B.b_absxk = 3.3121686421112381E-170;
      Subsystem_Fy_B.j = 4;
      while (Subsystem_Fy_B.j + 1 <= 8) {
        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].re);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].im);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.j++;
      }

      Subsystem_Fy_B.anrm = Subsystem_Fy_B.b_absxk * std::sqrt
        (Subsystem_Fy_B.anrm);
      Subsystem_Fy_B.j = 4;
      while (Subsystem_Fy_B.j + 1 <= 8) {
        Subsystem_Fy_B.b_absxk = V[Subsystem_Fy_B.j].re;
        Subsystem_Fy_B.cfromc = V[Subsystem_Fy_B.j].im;
        if (Subsystem_Fy_B.cfromc == 0.0) {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = 0.0;
        } else if (Subsystem_Fy_B.b_absxk == 0.0) {
          V[Subsystem_Fy_B.j].re = 0.0;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        } else {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        }

        Subsystem_Fy_B.j++;
      }

      Subsystem_Fy_B.anrm = 0.0;
      Subsystem_Fy_B.b_absxk = 3.3121686421112381E-170;
      Subsystem_Fy_B.j = 8;
      while (Subsystem_Fy_B.j + 1 <= 12) {
        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].re);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].im);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.j++;
      }

      Subsystem_Fy_B.anrm = Subsystem_Fy_B.b_absxk * std::sqrt
        (Subsystem_Fy_B.anrm);
      Subsystem_Fy_B.j = 8;
      while (Subsystem_Fy_B.j + 1 <= 12) {
        Subsystem_Fy_B.b_absxk = V[Subsystem_Fy_B.j].re;
        Subsystem_Fy_B.cfromc = V[Subsystem_Fy_B.j].im;
        if (Subsystem_Fy_B.cfromc == 0.0) {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = 0.0;
        } else if (Subsystem_Fy_B.b_absxk == 0.0) {
          V[Subsystem_Fy_B.j].re = 0.0;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        } else {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        }

        Subsystem_Fy_B.j++;
      }

      Subsystem_Fy_B.anrm = 0.0;
      Subsystem_Fy_B.b_absxk = 3.3121686421112381E-170;
      Subsystem_Fy_B.j = 12;
      while (Subsystem_Fy_B.j + 1 <= 16) {
        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].re);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.cfromc = std::abs(V[Subsystem_Fy_B.j].im);
        if (Subsystem_Fy_B.cfromc > Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.anrm = Subsystem_Fy_B.anrm * Subsystem_Fy_B.ctoc *
            Subsystem_Fy_B.ctoc + 1.0;
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.cfromc;
        } else {
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.anrm += Subsystem_Fy_B.ctoc * Subsystem_Fy_B.ctoc;
        }

        Subsystem_Fy_B.j++;
      }

      Subsystem_Fy_B.anrm = Subsystem_Fy_B.b_absxk * std::sqrt
        (Subsystem_Fy_B.anrm);
      Subsystem_Fy_B.j = 12;
      while (Subsystem_Fy_B.j + 1 <= 16) {
        Subsystem_Fy_B.b_absxk = V[Subsystem_Fy_B.j].re;
        Subsystem_Fy_B.cfromc = V[Subsystem_Fy_B.j].im;
        if (Subsystem_Fy_B.cfromc == 0.0) {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = 0.0;
        } else if (Subsystem_Fy_B.b_absxk == 0.0) {
          V[Subsystem_Fy_B.j].re = 0.0;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        } else {
          V[Subsystem_Fy_B.j].re = Subsystem_Fy_B.b_absxk / Subsystem_Fy_B.anrm;
          V[Subsystem_Fy_B.j].im = Subsystem_Fy_B.cfromc / Subsystem_Fy_B.anrm;
        }

        Subsystem_Fy_B.j++;
      }

      if (Subsystem_Fy_B.beta1[0].im == 0.0) {
        if (D[0].im == 0.0) {
          Subsystem_Fy_B.anrm = D[0].re / Subsystem_Fy_B.beta1[0].re;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[0].re == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = D[0].im / Subsystem_Fy_B.beta1[0].re;
        } else {
          Subsystem_Fy_B.anrm = D[0].re / Subsystem_Fy_B.beta1[0].re;
          Subsystem_Fy_B.b_absxk = D[0].im / Subsystem_Fy_B.beta1[0].re;
        }
      } else if (Subsystem_Fy_B.beta1[0].re == 0.0) {
        if (D[0].re == 0.0) {
          Subsystem_Fy_B.anrm = D[0].im / Subsystem_Fy_B.beta1[0].im;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[0].im == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = -(D[0].re / Subsystem_Fy_B.beta1[0].im);
        } else {
          Subsystem_Fy_B.anrm = D[0].im / Subsystem_Fy_B.beta1[0].im;
          Subsystem_Fy_B.b_absxk = -(D[0].re / Subsystem_Fy_B.beta1[0].im);
        }
      } else {
        Subsystem_Fy_B.b_absxk = std::abs(Subsystem_Fy_B.beta1[0].re);
        Subsystem_Fy_B.anrm = std::abs(Subsystem_Fy_B.beta1[0].im);
        if (Subsystem_Fy_B.b_absxk > Subsystem_Fy_B.anrm) {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[0].im /
            Subsystem_Fy_B.beta1[0].re;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [0].im + Subsystem_Fy_B.beta1[0].re;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[0].im + D[0].re) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (D[0].im - Subsystem_Fy_B.b_absxk * D[0].re) /
            Subsystem_Fy_B.cfromc;
        } else if (Subsystem_Fy_B.anrm == Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.beta1[0].re > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.beta1[0].im > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.anrm = (D[0].re * Subsystem_Fy_B.cfromc + D[0].im *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.b_absxk = (D[0].im * Subsystem_Fy_B.cfromc - D[0].re *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
        } else {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[0].re /
            Subsystem_Fy_B.beta1[0].im;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [0].re + Subsystem_Fy_B.beta1[0].im;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[0].re + D[0].im) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (Subsystem_Fy_B.b_absxk * D[0].im - D[0].re) /
            Subsystem_Fy_B.cfromc;
        }
      }

      D[0].re = Subsystem_Fy_B.anrm;
      D[0].im = Subsystem_Fy_B.b_absxk;
      if (Subsystem_Fy_B.beta1[1].im == 0.0) {
        if (D[1].im == 0.0) {
          Subsystem_Fy_B.anrm = D[1].re / Subsystem_Fy_B.beta1[1].re;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[1].re == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = D[1].im / Subsystem_Fy_B.beta1[1].re;
        } else {
          Subsystem_Fy_B.anrm = D[1].re / Subsystem_Fy_B.beta1[1].re;
          Subsystem_Fy_B.b_absxk = D[1].im / Subsystem_Fy_B.beta1[1].re;
        }
      } else if (Subsystem_Fy_B.beta1[1].re == 0.0) {
        if (D[1].re == 0.0) {
          Subsystem_Fy_B.anrm = D[1].im / Subsystem_Fy_B.beta1[1].im;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[1].im == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = -(D[1].re / Subsystem_Fy_B.beta1[1].im);
        } else {
          Subsystem_Fy_B.anrm = D[1].im / Subsystem_Fy_B.beta1[1].im;
          Subsystem_Fy_B.b_absxk = -(D[1].re / Subsystem_Fy_B.beta1[1].im);
        }
      } else {
        Subsystem_Fy_B.b_absxk = std::abs(Subsystem_Fy_B.beta1[1].re);
        Subsystem_Fy_B.anrm = std::abs(Subsystem_Fy_B.beta1[1].im);
        if (Subsystem_Fy_B.b_absxk > Subsystem_Fy_B.anrm) {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[1].im /
            Subsystem_Fy_B.beta1[1].re;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [1].im + Subsystem_Fy_B.beta1[1].re;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[1].im + D[1].re) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (D[1].im - Subsystem_Fy_B.b_absxk * D[1].re) /
            Subsystem_Fy_B.cfromc;
        } else if (Subsystem_Fy_B.anrm == Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.beta1[1].re > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.beta1[1].im > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.anrm = (D[1].re * Subsystem_Fy_B.cfromc + D[1].im *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.b_absxk = (D[1].im * Subsystem_Fy_B.cfromc - D[1].re *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
        } else {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[1].re /
            Subsystem_Fy_B.beta1[1].im;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [1].re + Subsystem_Fy_B.beta1[1].im;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[1].re + D[1].im) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (Subsystem_Fy_B.b_absxk * D[1].im - D[1].re) /
            Subsystem_Fy_B.cfromc;
        }
      }

      D[1].re = Subsystem_Fy_B.anrm;
      D[1].im = Subsystem_Fy_B.b_absxk;
      if (Subsystem_Fy_B.beta1[2].im == 0.0) {
        if (D[2].im == 0.0) {
          Subsystem_Fy_B.anrm = D[2].re / Subsystem_Fy_B.beta1[2].re;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[2].re == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = D[2].im / Subsystem_Fy_B.beta1[2].re;
        } else {
          Subsystem_Fy_B.anrm = D[2].re / Subsystem_Fy_B.beta1[2].re;
          Subsystem_Fy_B.b_absxk = D[2].im / Subsystem_Fy_B.beta1[2].re;
        }
      } else if (Subsystem_Fy_B.beta1[2].re == 0.0) {
        if (D[2].re == 0.0) {
          Subsystem_Fy_B.anrm = D[2].im / Subsystem_Fy_B.beta1[2].im;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[2].im == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = -(D[2].re / Subsystem_Fy_B.beta1[2].im);
        } else {
          Subsystem_Fy_B.anrm = D[2].im / Subsystem_Fy_B.beta1[2].im;
          Subsystem_Fy_B.b_absxk = -(D[2].re / Subsystem_Fy_B.beta1[2].im);
        }
      } else {
        Subsystem_Fy_B.b_absxk = std::abs(Subsystem_Fy_B.beta1[2].re);
        Subsystem_Fy_B.anrm = std::abs(Subsystem_Fy_B.beta1[2].im);
        if (Subsystem_Fy_B.b_absxk > Subsystem_Fy_B.anrm) {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[2].im /
            Subsystem_Fy_B.beta1[2].re;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [2].im + Subsystem_Fy_B.beta1[2].re;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[2].im + D[2].re) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (D[2].im - Subsystem_Fy_B.b_absxk * D[2].re) /
            Subsystem_Fy_B.cfromc;
        } else if (Subsystem_Fy_B.anrm == Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.beta1[2].re > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.beta1[2].im > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.anrm = (D[2].re * Subsystem_Fy_B.cfromc + D[2].im *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.b_absxk = (D[2].im * Subsystem_Fy_B.cfromc - D[2].re *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
        } else {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[2].re /
            Subsystem_Fy_B.beta1[2].im;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [2].re + Subsystem_Fy_B.beta1[2].im;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[2].re + D[2].im) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (Subsystem_Fy_B.b_absxk * D[2].im - D[2].re) /
            Subsystem_Fy_B.cfromc;
        }
      }

      D[2].re = Subsystem_Fy_B.anrm;
      D[2].im = Subsystem_Fy_B.b_absxk;
      if (Subsystem_Fy_B.beta1[3].im == 0.0) {
        if (D[3].im == 0.0) {
          Subsystem_Fy_B.anrm = D[3].re / Subsystem_Fy_B.beta1[3].re;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[3].re == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = D[3].im / Subsystem_Fy_B.beta1[3].re;
        } else {
          Subsystem_Fy_B.anrm = D[3].re / Subsystem_Fy_B.beta1[3].re;
          Subsystem_Fy_B.b_absxk = D[3].im / Subsystem_Fy_B.beta1[3].re;
        }
      } else if (Subsystem_Fy_B.beta1[3].re == 0.0) {
        if (D[3].re == 0.0) {
          Subsystem_Fy_B.anrm = D[3].im / Subsystem_Fy_B.beta1[3].im;
          Subsystem_Fy_B.b_absxk = 0.0;
        } else if (D[3].im == 0.0) {
          Subsystem_Fy_B.anrm = 0.0;
          Subsystem_Fy_B.b_absxk = -(D[3].re / Subsystem_Fy_B.beta1[3].im);
        } else {
          Subsystem_Fy_B.anrm = D[3].im / Subsystem_Fy_B.beta1[3].im;
          Subsystem_Fy_B.b_absxk = -(D[3].re / Subsystem_Fy_B.beta1[3].im);
        }
      } else {
        Subsystem_Fy_B.b_absxk = std::abs(Subsystem_Fy_B.beta1[3].re);
        Subsystem_Fy_B.anrm = std::abs(Subsystem_Fy_B.beta1[3].im);
        if (Subsystem_Fy_B.b_absxk > Subsystem_Fy_B.anrm) {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[3].im /
            Subsystem_Fy_B.beta1[3].re;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [3].im + Subsystem_Fy_B.beta1[3].re;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[3].im + D[3].re) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (D[3].im - Subsystem_Fy_B.b_absxk * D[3].re) /
            Subsystem_Fy_B.cfromc;
        } else if (Subsystem_Fy_B.anrm == Subsystem_Fy_B.b_absxk) {
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.beta1[3].re > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.ctoc = Subsystem_Fy_B.beta1[3].im > 0.0 ? 0.5 : -0.5;
          Subsystem_Fy_B.anrm = (D[3].re * Subsystem_Fy_B.cfromc + D[3].im *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
          Subsystem_Fy_B.b_absxk = (D[3].im * Subsystem_Fy_B.cfromc - D[3].re *
            Subsystem_Fy_B.ctoc) / Subsystem_Fy_B.b_absxk;
        } else {
          Subsystem_Fy_B.b_absxk = Subsystem_Fy_B.beta1[3].re /
            Subsystem_Fy_B.beta1[3].im;
          Subsystem_Fy_B.cfromc = Subsystem_Fy_B.b_absxk * Subsystem_Fy_B.beta1
            [3].re + Subsystem_Fy_B.beta1[3].im;
          Subsystem_Fy_B.anrm = (Subsystem_Fy_B.b_absxk * D[3].re + D[3].im) /
            Subsystem_Fy_B.cfromc;
          Subsystem_Fy_B.b_absxk = (Subsystem_Fy_B.b_absxk * D[3].im - D[3].re) /
            Subsystem_Fy_B.cfromc;
        }
      }

      D[3].re = Subsystem_Fy_B.anrm;
      D[3].im = Subsystem_Fy_B.b_absxk;
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_rotm2quat(const real_T R[9], real_T
  quat[4])
{
  boolean_T exitg1;
  Subsystem_Fy_B.K12 = R[1] + R[3];
  Subsystem_Fy_B.K13 = R[2] + R[6];
  Subsystem_Fy_B.K14 = R[5] - R[7];
  Subsystem_Fy_B.K23 = R[5] + R[7];
  Subsystem_Fy_B.K24 = R[6] - R[2];
  Subsystem_Fy_B.K34 = R[1] - R[3];
  Subsystem_Fy_B.R[0] = ((R[0] - R[4]) - R[8]) / 3.0;
  Subsystem_Fy_B.R[4] = Subsystem_Fy_B.K12 / 3.0;
  Subsystem_Fy_B.R[8] = Subsystem_Fy_B.K13 / 3.0;
  Subsystem_Fy_B.R[12] = Subsystem_Fy_B.K14 / 3.0;
  Subsystem_Fy_B.R[1] = Subsystem_Fy_B.K12 / 3.0;
  Subsystem_Fy_B.R[5] = ((R[4] - R[0]) - R[8]) / 3.0;
  Subsystem_Fy_B.R[9] = Subsystem_Fy_B.K23 / 3.0;
  Subsystem_Fy_B.R[13] = Subsystem_Fy_B.K24 / 3.0;
  Subsystem_Fy_B.R[2] = Subsystem_Fy_B.K13 / 3.0;
  Subsystem_Fy_B.R[6] = Subsystem_Fy_B.K23 / 3.0;
  Subsystem_Fy_B.R[10] = ((R[8] - R[0]) - R[4]) / 3.0;
  Subsystem_Fy_B.R[14] = Subsystem_Fy_B.K34 / 3.0;
  Subsystem_Fy_B.R[3] = Subsystem_Fy_B.K14 / 3.0;
  Subsystem_Fy_B.R[7] = Subsystem_Fy_B.K24 / 3.0;
  Subsystem_Fy_B.R[11] = Subsystem_Fy_B.K34 / 3.0;
  Subsystem_Fy_B.R[15] = ((R[0] + R[4]) + R[8]) / 3.0;
  Subsystem_Fy_eig(Subsystem_Fy_B.R, Subsystem_Fy_B.eigVec,
                   Subsystem_Fy_B.eigVal);
  Subsystem_Fy_B.varargin_1[0] = Subsystem_Fy_B.eigVal[0].re;
  Subsystem_Fy_B.varargin_1[1] = Subsystem_Fy_B.eigVal[1].re;
  Subsystem_Fy_B.varargin_1[2] = Subsystem_Fy_B.eigVal[2].re;
  Subsystem_Fy_B.varargin_1[3] = Subsystem_Fy_B.eigVal[3].re;
  if (!rtIsNaN(Subsystem_Fy_B.eigVal[0].re)) {
    Subsystem_Fy_B.idx = 1;
  } else {
    Subsystem_Fy_B.idx = 0;
    Subsystem_Fy_B.k = 2;
    exitg1 = false;
    while ((!exitg1) && (Subsystem_Fy_B.k < 5)) {
      if (!rtIsNaN(Subsystem_Fy_B.varargin_1[Subsystem_Fy_B.k - 1])) {
        Subsystem_Fy_B.idx = Subsystem_Fy_B.k;
        exitg1 = true;
      } else {
        Subsystem_Fy_B.k++;
      }
    }
  }

  if (Subsystem_Fy_B.idx != 0) {
    Subsystem_Fy_B.K12 = Subsystem_Fy_B.varargin_1[Subsystem_Fy_B.idx - 1];
    Subsystem_Fy_B.k = Subsystem_Fy_B.idx - 1;
    while (Subsystem_Fy_B.idx + 1 < 5) {
      if (Subsystem_Fy_B.K12 < Subsystem_Fy_B.varargin_1[Subsystem_Fy_B.idx]) {
        Subsystem_Fy_B.K12 = Subsystem_Fy_B.varargin_1[Subsystem_Fy_B.idx];
        Subsystem_Fy_B.k = Subsystem_Fy_B.idx;
      }

      Subsystem_Fy_B.idx++;
    }

    Subsystem_Fy_B.idx = Subsystem_Fy_B.k;
  }

  Subsystem_Fy_B.idx <<= 2;
  quat[0] = Subsystem_Fy_B.eigVec[Subsystem_Fy_B.idx + 3].re;
  quat[1] = Subsystem_Fy_B.eigVec[Subsystem_Fy_B.idx].re;
  quat[2] = Subsystem_Fy_B.eigVec[Subsystem_Fy_B.idx + 1].re;
  quat[3] = Subsystem_Fy_B.eigVec[Subsystem_Fy_B.idx + 2].re;
  if (quat[0] < 0.0) {
    quat[0] = -quat[0];
    quat[1] = -quat[1];
    quat[2] = -quat[2];
    quat[3] = -quat[3];
  }
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = std::atan2(static_cast<real_T>(u0_0), static_cast<real_T>(u1_0));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }

  return y;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::project_state_to_constrain_fcn(const real_T rf1[3],
  const real_T rfa[3], const real_T MAG_NORTH[3], const real_T rsa[3], const
  real_T rs2[3], const real_T s1[3], const real_T am1[3], const real_T mm1[3],
  const real_T am2[3], const real_T mm2[3], real_T xproj[42])
{
  Subsystem_Fy_B.yaxis_idx_2 = std::sqrt((am1[0] * am1[0] + am1[1] * am1[1]) +
    am1[2] * am1[2]);
  Subsystem_Fy_B.yaxis_idx_0 = am1[0] / Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.yaxis_idx_1 = am1[1] / Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.yaxis_idx_2 = am1[2] / Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.xaxis_idx_0 = std::sqrt((mm1[0] * mm1[0] + mm1[1] * mm1[1]) +
    mm1[2] * mm1[2]);
  Subsystem_Fy_B.mag_north_idx_0 = -mm1[0] / Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.mag_north_idx_1 = -mm1[1] / Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.mag_north_idx_2 = -mm1[2] / Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.xaxis_idx_0 = Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_2 - Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.yaxis_idx_2 - Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.mag_north_idx_1 = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_1 - Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.c_x = std::sqrt((Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0 + Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2) + Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.mag_north_idx_1);
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.xaxis_idx_0 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.R_f[0] = Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.R_f[1] = Subsystem_Fy_B.yaxis_idx_0;
  Subsystem_Fy_B.xaxis_idx_0 = Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.mag_north_idx_2 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.R_f[3] = Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.R_f[4] = Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.mag_north_idx_1 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.R_f[6] = Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.R_f[7] = Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.R_f[2] = Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.yaxis_idx_2 - Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.R_f[5] = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_0 - Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.R_f[8] = Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_1 - Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.yaxis_idx_2 = std::sqrt((am2[0] * am2[0] + am2[1] * am2[1]) +
    am2[2] * am2[2]);
  Subsystem_Fy_B.yaxis_idx_0 = am2[0] / Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.yaxis_idx_1 = am2[1] / Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.yaxis_idx_2 = am2[2] / Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.xaxis_idx_0 = std::sqrt((mm2[0] * mm2[0] + mm2[1] * mm2[1]) +
    mm2[2] * mm2[2]);
  Subsystem_Fy_B.mag_north_idx_0 = -mm2[0] / Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.mag_north_idx_1 = -mm2[1] / Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.mag_north_idx_2 = -mm2[2] / Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.xaxis_idx_0 = Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_2 - Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.yaxis_idx_2 - Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.mag_north_idx_1 = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_1 - Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.c_x = std::sqrt((Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0 + Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2) + Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.mag_north_idx_1);
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.xaxis_idx_0 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.xaxis_idx_0 = Subsystem_Fy_B.mag_north_idx_2 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.mag_north_idx_1 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_rotm2quat(Subsystem_Fy_B.R_f, Subsystem_Fy_B.fquat);
  Subsystem_Fy_B.R_f[0] = Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.R_f[1] = Subsystem_Fy_B.yaxis_idx_0;
  Subsystem_Fy_B.R_f[3] = Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.R_f[4] = Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.R_f[6] = Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.R_f[7] = Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.R_f[2] = Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_2 - Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.R_f[5] = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2 - Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.R_f[8] = Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.yaxis_idx_1 - Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_rotm2quat(Subsystem_Fy_B.R_f, Subsystem_Fy_B.squat);
  Subsystem_Fy_B.yaxis_idx_2 = ((Subsystem_Fy_B.squat[0] * Subsystem_Fy_B.squat
    [0] + Subsystem_Fy_B.squat[1] * Subsystem_Fy_B.squat[1]) +
    Subsystem_Fy_B.squat[2] * Subsystem_Fy_B.squat[2]) + Subsystem_Fy_B.squat[3]
    * Subsystem_Fy_B.squat[3];
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.squat[0] /
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.mag_north_idx_0 = -Subsystem_Fy_B.squat[1] /
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.mag_north_idx_1 = -Subsystem_Fy_B.squat[2] /
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.yaxis_idx_1 = -Subsystem_Fy_B.squat[3] /
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.yaxis_idx_2 = ((Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[0] - Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.fquat[1]) - Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.fquat[2]) - Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.fquat[3];
  Subsystem_Fy_B.yaxis_idx_0 = (Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[1] + Subsystem_Fy_B.fquat[0] *
    Subsystem_Fy_B.mag_north_idx_0) + (Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.fquat[3] - Subsystem_Fy_B.fquat[2] *
    Subsystem_Fy_B.yaxis_idx_1);
  Subsystem_Fy_B.xaxis_idx_0 = (Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[2] + Subsystem_Fy_B.fquat[0] *
    Subsystem_Fy_B.mag_north_idx_1) + (Subsystem_Fy_B.fquat[1] *
    Subsystem_Fy_B.yaxis_idx_1 - Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.fquat[3]);
  Subsystem_Fy_B.yaxis_idx_1 = (Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[3] + Subsystem_Fy_B.fquat[0] *
    Subsystem_Fy_B.yaxis_idx_1) + (Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.fquat[2] - Subsystem_Fy_B.fquat[1] *
    Subsystem_Fy_B.mag_north_idx_1);
  Subsystem_Fy_B.mag_north_idx_2 = std::sqrt(((Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_2 + Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_0) + Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0) + Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.yaxis_idx_1);
  Subsystem_Fy_B.yaxis_idx_2 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.yaxis_idx_0 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.xaxis_idx_0 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.yaxis_idx_1 /
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.yaxis_idx_1 = (Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0 + Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2) * 2.0;
  if (Subsystem_Fy_B.yaxis_idx_1 < -1.0) {
    Subsystem_Fy_B.yaxis_idx_1 = -1.0;
  }

  if (Subsystem_Fy_B.yaxis_idx_1 > 1.0) {
    Subsystem_Fy_B.yaxis_idx_1 = 1.0;
  }

  Subsystem_Fy_B.yaxis_idx_1 = std::asin(Subsystem_Fy_B.yaxis_idx_1) / 2.0;
  Subsystem_Fy_B.yaxis_idx_2 = rt_atan2d_snf((Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2 - Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_0) * -2.0, ((Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_2 - Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_0) + Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0) - Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2) / 2.0;
  Subsystem_Fy_B.mag_north_idx_2 = std::cos(Subsystem_Fy_B.yaxis_idx_1);
  Subsystem_Fy_B.yaxis_idx_1 = std::sin(Subsystem_Fy_B.yaxis_idx_1);
  Subsystem_Fy_B.mag_north_idx_0 = std::cos(Subsystem_Fy_B.yaxis_idx_2);
  Subsystem_Fy_B.mag_north_idx_1 = std::sin(Subsystem_Fy_B.yaxis_idx_2);
  Subsystem_Fy_B.yaxis_idx_2 = 0.9999761767522789 *
    Subsystem_Fy_B.mag_north_idx_2 * Subsystem_Fy_B.mag_north_idx_0 -
    -0.0069026029796791189 * Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_1;
  Subsystem_Fy_B.yaxis_idx_0 = 0.9999761767522789 *
    Subsystem_Fy_B.mag_north_idx_2 * Subsystem_Fy_B.mag_north_idx_1 +
    -0.0069026029796791189 * Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.xaxis_idx_0 = 0.9999761767522789 * Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_1 + -0.0069026029796791189 *
    Subsystem_Fy_B.mag_north_idx_2 * Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.yaxis_idx_1 = 0.9999761767522789 * Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.mag_north_idx_0 - -0.0069026029796791189 *
    Subsystem_Fy_B.mag_north_idx_2 * Subsystem_Fy_B.mag_north_idx_1;
  Subsystem_Fy_B.c_x = ((Subsystem_Fy_B.yaxis_idx_2 * Subsystem_Fy_B.yaxis_idx_2
    + Subsystem_Fy_B.yaxis_idx_0 * Subsystem_Fy_B.yaxis_idx_0) +
                        Subsystem_Fy_B.xaxis_idx_0 * Subsystem_Fy_B.xaxis_idx_0)
    + Subsystem_Fy_B.yaxis_idx_1 * Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.yaxis_idx_2 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.mag_north_idx_0 = -Subsystem_Fy_B.yaxis_idx_0 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.mag_north_idx_1 = -Subsystem_Fy_B.xaxis_idx_0 /
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.yaxis_idx_1 = -Subsystem_Fy_B.yaxis_idx_1 / Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.squat[0] = ((Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[0] - Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.fquat[1]) - Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.fquat[2]) - Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.fquat[3];
  Subsystem_Fy_B.squat[1] = (Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[1] + Subsystem_Fy_B.fquat[0] *
    Subsystem_Fy_B.mag_north_idx_0) + (Subsystem_Fy_B.mag_north_idx_1 *
    Subsystem_Fy_B.fquat[3] - Subsystem_Fy_B.fquat[2] *
    Subsystem_Fy_B.yaxis_idx_1);
  Subsystem_Fy_B.squat[2] = (Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[2] + Subsystem_Fy_B.fquat[0] *
    Subsystem_Fy_B.mag_north_idx_1) + (Subsystem_Fy_B.fquat[1] *
    Subsystem_Fy_B.yaxis_idx_1 - Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.fquat[3]);
  Subsystem_Fy_B.squat[3] = (Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.fquat[3] + Subsystem_Fy_B.fquat[0] *
    Subsystem_Fy_B.yaxis_idx_1) + (Subsystem_Fy_B.mag_north_idx_0 *
    Subsystem_Fy_B.fquat[2] - Subsystem_Fy_B.fquat[1] *
    Subsystem_Fy_B.mag_north_idx_1);
  Subsystem_Fy_B.yaxis_idx_1 = ((Subsystem_Fy_B.fquat[0] * Subsystem_Fy_B.fquat
    [0] + Subsystem_Fy_B.fquat[1] * Subsystem_Fy_B.fquat[1]) +
    Subsystem_Fy_B.fquat[2] * Subsystem_Fy_B.fquat[2]) + Subsystem_Fy_B.fquat[3]
    * Subsystem_Fy_B.fquat[3];
  Subsystem_Fy_B.yaxis_idx_2 = Subsystem_Fy_B.fquat[0] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.yaxis_idx_0 = -Subsystem_Fy_B.fquat[1] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.xaxis_idx_0 = -Subsystem_Fy_B.fquat[2] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.yaxis_idx_1 = -Subsystem_Fy_B.fquat[3] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.mag_north_idx_2 = std::sqrt(((Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_2 + Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_0) + Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0) + Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.yaxis_idx_1);
  Subsystem_Fy_B.yaxis_idx_2 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.yaxis_idx_0 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.xaxis_idx_0 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.yaxis_idx_1 /
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.yaxis_idx_1 = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_0;
  Subsystem_Fy_B.mag_north_idx_1 = Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.c_x = Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.dcm[0] = ((Subsystem_Fy_B.yaxis_idx_1 +
    Subsystem_Fy_B.mag_north_idx_0) - Subsystem_Fy_B.mag_north_idx_1) -
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.dcm_tmp = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.dcm_tmp_l = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.dcm[3] = (Subsystem_Fy_B.dcm_tmp + Subsystem_Fy_B.dcm_tmp_l) *
    2.0;
  Subsystem_Fy_B.dcm_tmp_d = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.dcm_tmp_dy = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.dcm[6] = (Subsystem_Fy_B.dcm_tmp_d - Subsystem_Fy_B.dcm_tmp_dy)
    * 2.0;
  Subsystem_Fy_B.dcm[1] = (Subsystem_Fy_B.dcm_tmp - Subsystem_Fy_B.dcm_tmp_l) *
    2.0;
  Subsystem_Fy_B.yaxis_idx_1 -= Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.dcm[4] = (Subsystem_Fy_B.yaxis_idx_1 +
    Subsystem_Fy_B.mag_north_idx_1) - Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.dcm_tmp = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_0;
  Subsystem_Fy_B.dcm[7] = (Subsystem_Fy_B.mag_north_idx_0 +
    Subsystem_Fy_B.dcm_tmp) * 2.0;
  Subsystem_Fy_B.dcm[2] = (Subsystem_Fy_B.dcm_tmp_d + Subsystem_Fy_B.dcm_tmp_dy)
    * 2.0;
  Subsystem_Fy_B.dcm[5] = (Subsystem_Fy_B.mag_north_idx_0 -
    Subsystem_Fy_B.dcm_tmp) * 2.0;
  Subsystem_Fy_B.dcm[8] = (Subsystem_Fy_B.yaxis_idx_1 -
    Subsystem_Fy_B.mag_north_idx_1) + Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.yaxis_idx_1 = ((Subsystem_Fy_B.squat[0] * Subsystem_Fy_B.squat
    [0] + Subsystem_Fy_B.squat[1] * Subsystem_Fy_B.squat[1]) +
    Subsystem_Fy_B.squat[2] * Subsystem_Fy_B.squat[2]) + Subsystem_Fy_B.squat[3]
    * Subsystem_Fy_B.squat[3];
  Subsystem_Fy_B.yaxis_idx_2 = Subsystem_Fy_B.squat[0] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.yaxis_idx_0 = -Subsystem_Fy_B.squat[1] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.xaxis_idx_0 = -Subsystem_Fy_B.squat[2] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.yaxis_idx_1 = -Subsystem_Fy_B.squat[3] /
    Subsystem_Fy_B.yaxis_idx_1;
  Subsystem_Fy_B.mag_north_idx_2 = std::sqrt(((Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_2 + Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_0) + Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0) + Subsystem_Fy_B.yaxis_idx_1 *
    Subsystem_Fy_B.yaxis_idx_1);
  Subsystem_Fy_B.yaxis_idx_2 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.yaxis_idx_0 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.xaxis_idx_0 /= Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.mag_north_idx_2 = Subsystem_Fy_B.yaxis_idx_1 /
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.yaxis_idx_1 = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_2;
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.yaxis_idx_0;
  Subsystem_Fy_B.mag_north_idx_1 = Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.c_x = Subsystem_Fy_B.mag_north_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.b_dcm[0] = ((Subsystem_Fy_B.yaxis_idx_1 +
    Subsystem_Fy_B.mag_north_idx_0) - Subsystem_Fy_B.mag_north_idx_1) -
    Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.dcm_tmp = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.dcm_tmp_l = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.b_dcm[3] = (Subsystem_Fy_B.dcm_tmp + Subsystem_Fy_B.dcm_tmp_l) *
    2.0;
  Subsystem_Fy_B.dcm_tmp_d = Subsystem_Fy_B.yaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.dcm_tmp_dy = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.xaxis_idx_0;
  Subsystem_Fy_B.b_dcm[6] = (Subsystem_Fy_B.dcm_tmp_d -
    Subsystem_Fy_B.dcm_tmp_dy) * 2.0;
  Subsystem_Fy_B.b_dcm[1] = (Subsystem_Fy_B.dcm_tmp - Subsystem_Fy_B.dcm_tmp_l) *
    2.0;
  Subsystem_Fy_B.yaxis_idx_1 -= Subsystem_Fy_B.mag_north_idx_0;
  Subsystem_Fy_B.b_dcm[4] = (Subsystem_Fy_B.yaxis_idx_1 +
    Subsystem_Fy_B.mag_north_idx_1) - Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.mag_north_idx_0 = Subsystem_Fy_B.xaxis_idx_0 *
    Subsystem_Fy_B.mag_north_idx_2;
  Subsystem_Fy_B.dcm_tmp = Subsystem_Fy_B.yaxis_idx_2 *
    Subsystem_Fy_B.yaxis_idx_0;
  Subsystem_Fy_B.b_dcm[7] = (Subsystem_Fy_B.mag_north_idx_0 +
    Subsystem_Fy_B.dcm_tmp) * 2.0;
  Subsystem_Fy_B.b_dcm[2] = (Subsystem_Fy_B.dcm_tmp_d +
    Subsystem_Fy_B.dcm_tmp_dy) * 2.0;
  Subsystem_Fy_B.b_dcm[5] = (Subsystem_Fy_B.mag_north_idx_0 -
    Subsystem_Fy_B.dcm_tmp) * 2.0;
  Subsystem_Fy_B.b_dcm[8] = (Subsystem_Fy_B.yaxis_idx_1 -
    Subsystem_Fy_B.mag_north_idx_1) + Subsystem_Fy_B.c_x;
  Subsystem_Fy_B.yaxis_idx_2 = std::sqrt((MAG_NORTH[0] * MAG_NORTH[0] +
    MAG_NORTH[1] * MAG_NORTH[1]) + MAG_NORTH[2] * MAG_NORTH[2]);
  Subsystem_Fy_B.yaxis_idx_0 = -rf1[0] + rfa[0];
  Subsystem_Fy_B.xaxis_idx_0 = -rsa[0] + rs2[0];
  Subsystem_Fy_B.yaxis_idx_1 = -rf1[1] + rfa[1];
  Subsystem_Fy_B.mag_north_idx_2 = -rsa[1] + rs2[1];
  Subsystem_Fy_B.mag_north_idx_0 = -rf1[2] + rfa[2];
  Subsystem_Fy_B.mag_north_idx_1 = -rsa[2] + rs2[2];
  for (Subsystem_Fy_B.i1 = 0; Subsystem_Fy_B.i1 < 3; Subsystem_Fy_B.i1++) {
    Subsystem_Fy_B.c_x = s1[Subsystem_Fy_B.i1];
    Subsystem_Fy_B.s1[Subsystem_Fy_B.i1] =
      ((Subsystem_Fy_B.dcm[Subsystem_Fy_B.i1 + 3] * Subsystem_Fy_B.yaxis_idx_1 +
        Subsystem_Fy_B.dcm[Subsystem_Fy_B.i1] * Subsystem_Fy_B.yaxis_idx_0) +
       Subsystem_Fy_B.dcm[Subsystem_Fy_B.i1 + 6] *
       Subsystem_Fy_B.mag_north_idx_0) + Subsystem_Fy_B.c_x;
    xproj[Subsystem_Fy_B.i1] = Subsystem_Fy_B.c_x;
    xproj[Subsystem_Fy_B.i1 + 3] = 0.0;
    Subsystem_Fy_B.b_dcm_m[Subsystem_Fy_B.i1] =
      Subsystem_Fy_B.b_dcm[Subsystem_Fy_B.i1 + 6] *
      Subsystem_Fy_B.mag_north_idx_1 + (Subsystem_Fy_B.b_dcm[Subsystem_Fy_B.i1 +
      3] * Subsystem_Fy_B.mag_north_idx_2 +
      Subsystem_Fy_B.b_dcm[Subsystem_Fy_B.i1] * Subsystem_Fy_B.xaxis_idx_0);
  }

  xproj[6] = Subsystem_Fy_B.fquat[0];
  xproj[7] = Subsystem_Fy_B.fquat[1];
  xproj[8] = Subsystem_Fy_B.fquat[2];
  xproj[9] = Subsystem_Fy_B.fquat[3];
  xproj[19] = Subsystem_Fy_B.s1[0] + Subsystem_Fy_B.b_dcm_m[0];
  xproj[22] = 0.0;
  xproj[20] = Subsystem_Fy_B.s1[1] + Subsystem_Fy_B.b_dcm_m[1];
  xproj[23] = 0.0;
  xproj[21] = Subsystem_Fy_B.s1[2] + Subsystem_Fy_B.b_dcm_m[2];
  xproj[24] = 0.0;
  xproj[25] = Subsystem_Fy_B.squat[0];
  xproj[26] = Subsystem_Fy_B.squat[1];
  xproj[27] = Subsystem_Fy_B.squat[2];
  xproj[28] = Subsystem_Fy_B.squat[3];
  for (Subsystem_Fy_B.i1 = 0; Subsystem_Fy_B.i1 < 9; Subsystem_Fy_B.i1++) {
    xproj[Subsystem_Fy_B.i1 + 10] = 0.0;
    xproj[Subsystem_Fy_B.i1 + 29] = 0.0;
  }

  // xproj[38] = -0.0138053155886834; 
  xproj[38] = 0.0055;  // updated Block_yaw
  xproj[39] = MAG_NORTH[0] / Subsystem_Fy_B.yaxis_idx_2;
  xproj[40] = MAG_NORTH[1] / Subsystem_Fy_B.yaxis_idx_2;
  xproj[41] = MAG_NORTH[2] / Subsystem_Fy_B.yaxis_idx_2;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem__estimate_measurement(const real_T
  x_nominal_prev[42], real_T yh[13], real_T H[520])
{
  real_T t114;
  real_T t115;
  real_T t122;
  real_T t123;
  real_T t124;
  real_T t128;
  real_T t129;
  Subsystem_Fy_B.t2 = x_nominal_prev[6] * x_nominal_prev[6];
  Subsystem_Fy_B.t3_b = x_nominal_prev[7] * x_nominal_prev[7];
  Subsystem_Fy_B.t4_n = x_nominal_prev[8] * x_nominal_prev[8];
  Subsystem_Fy_B.t5_b = x_nominal_prev[9] * x_nominal_prev[9];
  Subsystem_Fy_B.t6_l = x_nominal_prev[25] * x_nominal_prev[25];
  Subsystem_Fy_B.t7 = x_nominal_prev[26] * x_nominal_prev[26];
  Subsystem_Fy_B.t8 = x_nominal_prev[27] * x_nominal_prev[27];
  Subsystem_Fy_B.t9 = x_nominal_prev[28] * x_nominal_prev[28];
  Subsystem_Fy_B.t12_h = 1.0 / (((Subsystem_Fy_B.t6_l + Subsystem_Fy_B.t7) +
    Subsystem_Fy_B.t8) + Subsystem_Fy_B.t9);
  Subsystem_Fy_B.t13_b = 1.0 / (((Subsystem_Fy_B.t2 + Subsystem_Fy_B.t3_b) +
    Subsystem_Fy_B.t4_n) + Subsystem_Fy_B.t5_b);
  Subsystem_Fy_B.t14_d = x_nominal_prev[6] * x_nominal_prev[7] *
    Subsystem_Fy_B.t13_b * 2.0;
  Subsystem_Fy_B.t15_e = x_nominal_prev[6] * x_nominal_prev[8] *
    Subsystem_Fy_B.t13_b * 2.0;
  Subsystem_Fy_B.t16_b = x_nominal_prev[6] * x_nominal_prev[9] *
    Subsystem_Fy_B.t13_b * 2.0;
  Subsystem_Fy_B.t17_j = x_nominal_prev[7] * x_nominal_prev[8] *
    Subsystem_Fy_B.t13_b * 2.0;
  Subsystem_Fy_B.t18_f = x_nominal_prev[7] * x_nominal_prev[9] *
    Subsystem_Fy_B.t13_b * 2.0;
  Subsystem_Fy_B.t19_a = x_nominal_prev[8] * x_nominal_prev[9] *
    Subsystem_Fy_B.t13_b * 2.0;
  Subsystem_Fy_B.t20_j = Subsystem_Fy_B.t6_l * Subsystem_Fy_B.t12_h;
  Subsystem_Fy_B.t21_j = Subsystem_Fy_B.t7 * Subsystem_Fy_B.t12_h;
  Subsystem_Fy_B.t22_o = Subsystem_Fy_B.t8 * Subsystem_Fy_B.t12_h;
  Subsystem_Fy_B.t23_n = Subsystem_Fy_B.t9 * Subsystem_Fy_B.t12_h;
  Subsystem_Fy_B.t9 = x_nominal_prev[25] * x_nominal_prev[26] *
    Subsystem_Fy_B.t12_h * 2.0;
  Subsystem_Fy_B.t8 = x_nominal_prev[25] * x_nominal_prev[27] *
    Subsystem_Fy_B.t12_h * 2.0;
  Subsystem_Fy_B.t7 = x_nominal_prev[25] * x_nominal_prev[28] *
    Subsystem_Fy_B.t12_h * 2.0;
  Subsystem_Fy_B.t6_l = x_nominal_prev[26] * x_nominal_prev[27] *
    Subsystem_Fy_B.t12_h * 2.0;
  Subsystem_Fy_B.t28 = x_nominal_prev[26] * x_nominal_prev[28] *
    Subsystem_Fy_B.t12_h * 2.0;
  Subsystem_Fy_B.t12_h = x_nominal_prev[27] * x_nominal_prev[28] *
    Subsystem_Fy_B.t12_h * 2.0;
  Subsystem_Fy_B.t30 = Subsystem_Fy_B.t2 * Subsystem_Fy_B.t13_b;
  Subsystem_Fy_B.t31 = Subsystem_Fy_B.t3_b * Subsystem_Fy_B.t13_b;
  Subsystem_Fy_B.t32 = Subsystem_Fy_B.t4_n * Subsystem_Fy_B.t13_b;
  Subsystem_Fy_B.t33 = Subsystem_Fy_B.t5_b * Subsystem_Fy_B.t13_b;
  Subsystem_Fy_B.t5_b = Subsystem_Fy_B.t21_j * 2.0;
  Subsystem_Fy_B.t13_b = Subsystem_Fy_B.t22_o * 2.0;
  Subsystem_Fy_B.t4_n = Subsystem_Fy_B.t23_n * 2.0;
  Subsystem_Fy_B.t3_b = Subsystem_Fy_B.t31 * 2.0;
  Subsystem_Fy_B.t2 = Subsystem_Fy_B.t32 * 2.0;
  Subsystem_Fy_B.t45_i = Subsystem_Fy_B.t33 * 2.0;
  Subsystem_Fy_B.t58 = Subsystem_Fy_B.t9 + Subsystem_Fy_B.t12_h;
  Subsystem_Fy_B.t59 = Subsystem_Fy_B.t8 + Subsystem_Fy_B.t28;
  Subsystem_Fy_B.t60 = Subsystem_Fy_B.t7 + Subsystem_Fy_B.t6_l;
  Subsystem_Fy_B.t61 = Subsystem_Fy_B.t14_d + Subsystem_Fy_B.t19_a;
  Subsystem_Fy_B.t62 = Subsystem_Fy_B.t15_e + Subsystem_Fy_B.t18_f;
  Subsystem_Fy_B.t63 = Subsystem_Fy_B.t16_b + Subsystem_Fy_B.t17_j;
  Subsystem_Fy_B.t64 = Subsystem_Fy_B.t9 + -Subsystem_Fy_B.t12_h;
  Subsystem_Fy_B.t65 = Subsystem_Fy_B.t8 + -Subsystem_Fy_B.t28;
  Subsystem_Fy_B.t66 = Subsystem_Fy_B.t7 + -Subsystem_Fy_B.t6_l;
  Subsystem_Fy_B.t67_o = Subsystem_Fy_B.t14_d + -Subsystem_Fy_B.t19_a;
  Subsystem_Fy_B.t68_n = Subsystem_Fy_B.t15_e + -Subsystem_Fy_B.t18_f;
  Subsystem_Fy_B.t69_m = Subsystem_Fy_B.t16_b + -Subsystem_Fy_B.t17_j;
  Subsystem_Fy_B.t70_c = 0.04272878004382006 * Subsystem_Fy_B.t63;
  Subsystem_Fy_B.t71_m = 0.02402453088668155 * Subsystem_Fy_B.t61;
  Subsystem_Fy_B.t72_m = -0.0041177090717768 * Subsystem_Fy_B.t62;
  Subsystem_Fy_B.t73_j = -0.0053222135219039794 * Subsystem_Fy_B.t60;
  Subsystem_Fy_B.t74_h = -0.098224457113952918 * Subsystem_Fy_B.t58;
  Subsystem_Fy_B.t75_c = -0.02268340689026823 * Subsystem_Fy_B.t59;
  Subsystem_Fy_B.t82 = (Subsystem_Fy_B.t3_b + Subsystem_Fy_B.t2) - 1.0;
  Subsystem_Fy_B.t83 = (Subsystem_Fy_B.t3_b + Subsystem_Fy_B.t45_i) - 1.0;
  Subsystem_Fy_B.t84 = (Subsystem_Fy_B.t2 + Subsystem_Fy_B.t45_i) - 1.0;
  Subsystem_Fy_B.t88 = (Subsystem_Fy_B.t5_b + Subsystem_Fy_B.t13_b) - 1.0;
  Subsystem_Fy_B.t89 = (Subsystem_Fy_B.t5_b + Subsystem_Fy_B.t4_n) - 1.0;
  Subsystem_Fy_B.t90 = (Subsystem_Fy_B.t13_b + Subsystem_Fy_B.t4_n) - 1.0;
  Subsystem_Fy_B.t94 = Subsystem_Fy_B.t59 * Subsystem_Fy_B.t63;
  Subsystem_Fy_B.t100 = ((Subsystem_Fy_B.t30 + Subsystem_Fy_B.t31) +
    -Subsystem_Fy_B.t32) + -Subsystem_Fy_B.t33;
  Subsystem_Fy_B.t101 = ((Subsystem_Fy_B.t30 + Subsystem_Fy_B.t32) +
    -Subsystem_Fy_B.t31) + -Subsystem_Fy_B.t33;
  Subsystem_Fy_B.t32 = ((Subsystem_Fy_B.t30 + Subsystem_Fy_B.t33) +
                        -Subsystem_Fy_B.t31) + -Subsystem_Fy_B.t32;
  Subsystem_Fy_B.t30 = ((Subsystem_Fy_B.t20_j + Subsystem_Fy_B.t21_j) +
                        -Subsystem_Fy_B.t22_o) + -Subsystem_Fy_B.t23_n;
  Subsystem_Fy_B.t31 = ((Subsystem_Fy_B.t20_j + Subsystem_Fy_B.t22_o) +
                        -Subsystem_Fy_B.t21_j) + -Subsystem_Fy_B.t23_n;
  Subsystem_Fy_B.t20_j = ((Subsystem_Fy_B.t20_j + Subsystem_Fy_B.t23_n) +
    -Subsystem_Fy_B.t21_j) + -Subsystem_Fy_B.t22_o;
  Subsystem_Fy_B.t76_c = 0.04272878004382006 * Subsystem_Fy_B.t68_n;
  Subsystem_Fy_B.t21_j = 0.02402453088668155 * Subsystem_Fy_B.t69_m;
  Subsystem_Fy_B.t22_o = -0.0041177090717768 * Subsystem_Fy_B.t67_o;
  Subsystem_Fy_B.t23_n = -0.0053222135219039794 * Subsystem_Fy_B.t65;
  Subsystem_Fy_B.t33 = -0.098224457113952918 * Subsystem_Fy_B.t66;
  Subsystem_Fy_B.t81 = -0.02268340689026823 * Subsystem_Fy_B.t64;
  Subsystem_Fy_B.t98 = Subsystem_Fy_B.t59 * Subsystem_Fy_B.t68_n;
  Subsystem_Fy_B.t106 = Subsystem_Fy_B.t64 * Subsystem_Fy_B.t68_n;
  Subsystem_Fy_B.t110 = 0.04272878004382006 * Subsystem_Fy_B.t100;
  Subsystem_Fy_B.t111 = 0.02402453088668155 * Subsystem_Fy_B.t101;
  Subsystem_Fy_B.t112 = -0.0041177090717768 * Subsystem_Fy_B.t32;
  Subsystem_Fy_B.t113 = -0.0053222135219039794 * Subsystem_Fy_B.t30;
  t114 = -0.098224457113952918 * Subsystem_Fy_B.t31;
  t115 = -0.02268340689026823 * Subsystem_Fy_B.t20_j;
  t122 = (Subsystem_Fy_B.t63 * Subsystem_Fy_B.t64 + Subsystem_Fy_B.t59 *
          Subsystem_Fy_B.t84) + -(Subsystem_Fy_B.t68_n * Subsystem_Fy_B.t88);
  t123 = 1.0 / ((Subsystem_Fy_B.t60 * Subsystem_Fy_B.t63 + Subsystem_Fy_B.t65 *
                 Subsystem_Fy_B.t68_n) + Subsystem_Fy_B.t84 * Subsystem_Fy_B.t90);
  t124 = t123 * t123;
  t129 = (Subsystem_Fy_B.t63 * Subsystem_Fy_B.t65 + -(Subsystem_Fy_B.t60 *
           Subsystem_Fy_B.t68_n)) * t122 * t124;
  t128 = 1.0 / (t122 * t122 * t124 + 1.0);
  yh[0] = ((x_nominal_prev[40] * Subsystem_Fy_B.t63 + x_nominal_prev[16]) -
           x_nominal_prev[41] * Subsystem_Fy_B.t68_n) - x_nominal_prev[39] *
    Subsystem_Fy_B.t84;
  yh[1] = ((x_nominal_prev[41] * Subsystem_Fy_B.t61 + x_nominal_prev[17]) -
           x_nominal_prev[39] * Subsystem_Fy_B.t69_m) - x_nominal_prev[40] *
    Subsystem_Fy_B.t83;
  yh[2] = ((x_nominal_prev[39] * Subsystem_Fy_B.t62 + x_nominal_prev[18]) -
           x_nominal_prev[40] * Subsystem_Fy_B.t67_o) - x_nominal_prev[41] *
    Subsystem_Fy_B.t82;
  yh[3] = ((x_nominal_prev[40] * Subsystem_Fy_B.t60 + x_nominal_prev[35]) -
           x_nominal_prev[41] * Subsystem_Fy_B.t65) - x_nominal_prev[39] *
    Subsystem_Fy_B.t90;
  yh[4] = ((x_nominal_prev[41] * Subsystem_Fy_B.t58 + x_nominal_prev[36]) -
           x_nominal_prev[39] * Subsystem_Fy_B.t66) - x_nominal_prev[40] *
    Subsystem_Fy_B.t89;
  yh[5] = ((x_nominal_prev[39] * Subsystem_Fy_B.t59 + x_nominal_prev[37]) -
           x_nominal_prev[40] * Subsystem_Fy_B.t64) - x_nominal_prev[41] *
    Subsystem_Fy_B.t88;
  yh[6] = ((((((x_nominal_prev[0] - x_nominal_prev[19]) + Subsystem_Fy_B.t72_m)
              + Subsystem_Fy_B.t33) + -Subsystem_Fy_B.t75_c) +
            -Subsystem_Fy_B.t21_j) - 0.04272878004382006 * Subsystem_Fy_B.t84) +
    -0.0053222135219039794 * Subsystem_Fy_B.t90;
  yh[7] = ((((((x_nominal_prev[1] - x_nominal_prev[20]) + Subsystem_Fy_B.t70_c)
              + Subsystem_Fy_B.t81) + -Subsystem_Fy_B.t73_j) +
            -Subsystem_Fy_B.t22_o) - 0.02402453088668155 * Subsystem_Fy_B.t83) +
    -0.098224457113952918 * Subsystem_Fy_B.t89;
  yh[8] = ((((((x_nominal_prev[2] - x_nominal_prev[21]) + Subsystem_Fy_B.t71_m)
              + Subsystem_Fy_B.t23_n) + -Subsystem_Fy_B.t74_h) +
            -Subsystem_Fy_B.t76_c) - -0.0041177090717768 * Subsystem_Fy_B.t82) +
    -0.02268340689026823 * Subsystem_Fy_B.t88;
  yh[9] = std::atan(t122 * t123) + -x_nominal_prev[38];
  yh[10] = x_nominal_prev[3];
  yh[11] = x_nominal_prev[4];
  yh[12] = x_nominal_prev[5];
  H[0] = 0.0;
  H[1] = 0.0;
  H[2] = 0.0;
  H[3] = 0.0;
  H[4] = 0.0;
  H[5] = 0.0;
  H[6] = 1.0;
  std::memset(&H[7], 0, 13U * sizeof(real_T));
  H[20] = 1.0;
  std::memset(&H[21], 0, 13U * sizeof(real_T));
  H[34] = 1.0;
  std::memset(&H[35], 0, 14U * sizeof(real_T));
  H[49] = 1.0;
  std::memset(&H[50], 0, 13U * sizeof(real_T));
  H[63] = 1.0;
  std::memset(&H[64], 0, 13U * sizeof(real_T));
  H[77] = 1.0;
  H[78] = x_nominal_prev[41] * Subsystem_Fy_B.t63 + x_nominal_prev[40] *
    Subsystem_Fy_B.t68_n;
  H[79] = -x_nominal_prev[40] * Subsystem_Fy_B.t61 + x_nominal_prev[41] *
    Subsystem_Fy_B.t101;
  H[80] = -x_nominal_prev[41] * Subsystem_Fy_B.t67_o - x_nominal_prev[40] *
    Subsystem_Fy_B.t32;
  H[81] = 0.0;
  H[82] = 0.0;
  H[83] = 0.0;
  H[84] = 0.0;
  H[85] = (-Subsystem_Fy_B.t71_m + Subsystem_Fy_B.t76_c) - Subsystem_Fy_B.t112;
  H[86] = (Subsystem_Fy_B.t70_c + -Subsystem_Fy_B.t22_o) + Subsystem_Fy_B.t111;
  H[87] = ((Subsystem_Fy_B.t63 * Subsystem_Fy_B.t88 + Subsystem_Fy_B.t106) *
           t123 + t129) * t128;
  H[88] = 0.0;
  H[89] = 0.0;
  H[90] = 0.0;
  H[91] = -x_nominal_prev[39] * Subsystem_Fy_B.t68_n - x_nominal_prev[41] *
    Subsystem_Fy_B.t100;
  H[92] = x_nominal_prev[39] * Subsystem_Fy_B.t61 + x_nominal_prev[41] *
    Subsystem_Fy_B.t69_m;
  H[93] = -x_nominal_prev[41] * Subsystem_Fy_B.t62 + x_nominal_prev[39] *
    Subsystem_Fy_B.t32;
  H[94] = 0.0;
  H[95] = 0.0;
  H[96] = 0.0;
  H[97] = (Subsystem_Fy_B.t71_m + -Subsystem_Fy_B.t76_c) + Subsystem_Fy_B.t112;
  H[98] = 0.0;
  H[99] = (-Subsystem_Fy_B.t72_m + Subsystem_Fy_B.t21_j) - Subsystem_Fy_B.t110;
  Subsystem_Fy_B.t71_m = t122 * t124;
  H[100] = ((Subsystem_Fy_B.t98 - Subsystem_Fy_B.t88 * Subsystem_Fy_B.t100) *
            t123 - (Subsystem_Fy_B.t68_n * Subsystem_Fy_B.t90 +
                    Subsystem_Fy_B.t65 * Subsystem_Fy_B.t100) *
            Subsystem_Fy_B.t71_m) * t128;
  H[101] = 0.0;
  H[102] = 0.0;
  H[103] = 0.0;
  H[104] = -x_nominal_prev[39] * Subsystem_Fy_B.t63 + x_nominal_prev[40] *
    Subsystem_Fy_B.t100;
  H[105] = -x_nominal_prev[40] * Subsystem_Fy_B.t69_m - x_nominal_prev[39] *
    Subsystem_Fy_B.t101;
  H[106] = x_nominal_prev[40] * Subsystem_Fy_B.t62 + x_nominal_prev[39] *
    Subsystem_Fy_B.t67_o;
  H[107] = 0.0;
  H[108] = 0.0;
  H[109] = 0.0;
  H[110] = (-Subsystem_Fy_B.t70_c + Subsystem_Fy_B.t22_o) - Subsystem_Fy_B.t111;
  H[111] = (Subsystem_Fy_B.t72_m + -Subsystem_Fy_B.t21_j) + Subsystem_Fy_B.t110;
  H[112] = 0.0;
  H[113] = ((Subsystem_Fy_B.t64 * Subsystem_Fy_B.t100 + Subsystem_Fy_B.t94) *
            t123 - (Subsystem_Fy_B.t63 * Subsystem_Fy_B.t90 + Subsystem_Fy_B.t60
                    * Subsystem_Fy_B.t100) * Subsystem_Fy_B.t71_m) * t128;
  std::memset(&H[114], 0, 81U * sizeof(real_T));
  H[195] = 1.0;
  std::memset(&H[196], 0, 13U * sizeof(real_T));
  H[209] = 1.0;
  std::memset(&H[210], 0, 13U * sizeof(real_T));
  H[223] = 1.0;
  std::memset(&H[224], 0, sizeof(real_T) << 4U);
  H[240] = -1.0;
  std::memset(&H[241], 0, 13U * sizeof(real_T));
  H[254] = -1.0;
  std::memset(&H[255], 0, 13U * sizeof(real_T));
  H[268] = -1.0;
  std::memset(&H[269], 0, 46U * sizeof(real_T));
  H[315] = x_nominal_prev[41] * Subsystem_Fy_B.t60 + x_nominal_prev[40] *
    Subsystem_Fy_B.t65;
  H[316] = -x_nominal_prev[40] * Subsystem_Fy_B.t58 + x_nominal_prev[41] *
    Subsystem_Fy_B.t31;
  H[317] = -x_nominal_prev[41] * Subsystem_Fy_B.t64 - x_nominal_prev[40] *
    Subsystem_Fy_B.t20_j;
  H[318] = 0.0;
  H[319] = (Subsystem_Fy_B.t74_h - Subsystem_Fy_B.t23_n) + t115;
  H[320] = (Subsystem_Fy_B.t81 + -Subsystem_Fy_B.t73_j) - t114;
  H[321] = ((Subsystem_Fy_B.t106 - Subsystem_Fy_B.t63 * Subsystem_Fy_B.t20_j) *
            t123 + t129) * -t128;
  H[322] = 0.0;
  H[323] = 0.0;
  H[324] = 0.0;
  H[325] = 0.0;
  H[326] = 0.0;
  H[327] = 0.0;
  H[328] = -x_nominal_prev[39] * Subsystem_Fy_B.t65 - x_nominal_prev[41] *
    Subsystem_Fy_B.t30;
  H[329] = x_nominal_prev[39] * Subsystem_Fy_B.t58 + x_nominal_prev[41] *
    Subsystem_Fy_B.t66;
  H[330] = -x_nominal_prev[41] * Subsystem_Fy_B.t59 + x_nominal_prev[39] *
    Subsystem_Fy_B.t20_j;
  H[331] = (Subsystem_Fy_B.t23_n + -Subsystem_Fy_B.t74_h) - t115;
  H[332] = 0.0;
  H[333] = (Subsystem_Fy_B.t75_c - Subsystem_Fy_B.t33) + Subsystem_Fy_B.t113;
  H[334] = ((Subsystem_Fy_B.t98 - Subsystem_Fy_B.t84 * Subsystem_Fy_B.t20_j) *
            t123 + (Subsystem_Fy_B.t65 * Subsystem_Fy_B.t84 +
                    Subsystem_Fy_B.t68_n * Subsystem_Fy_B.t30) *
            Subsystem_Fy_B.t71_m) * -t128;
  H[335] = 0.0;
  H[336] = 0.0;
  H[337] = 0.0;
  H[338] = 0.0;
  H[339] = 0.0;
  H[340] = 0.0;
  H[341] = -x_nominal_prev[39] * Subsystem_Fy_B.t60 + x_nominal_prev[40] *
    Subsystem_Fy_B.t30;
  H[342] = -x_nominal_prev[40] * Subsystem_Fy_B.t66 - x_nominal_prev[39] *
    Subsystem_Fy_B.t31;
  H[343] = x_nominal_prev[40] * Subsystem_Fy_B.t59 + x_nominal_prev[39] *
    Subsystem_Fy_B.t64;
  H[344] = (Subsystem_Fy_B.t73_j - Subsystem_Fy_B.t81) + t114;
  H[345] = (Subsystem_Fy_B.t33 + -Subsystem_Fy_B.t75_c) - Subsystem_Fy_B.t113;
  H[346] = 0.0;
  H[347] = ((Subsystem_Fy_B.t94 - Subsystem_Fy_B.t64 * Subsystem_Fy_B.t84) *
            t123 + (Subsystem_Fy_B.t60 * Subsystem_Fy_B.t84 + Subsystem_Fy_B.t63
                    * Subsystem_Fy_B.t30) * Subsystem_Fy_B.t71_m) * -t128;
  std::memset(&H[348], 0, 84U * sizeof(real_T));
  H[432] = 1.0;
  std::memset(&H[433], 0, 13U * sizeof(real_T));
  H[446] = 1.0;
  std::memset(&H[447], 0, 13U * sizeof(real_T));
  H[460] = 1.0;
  std::memset(&H[461], 0, sizeof(real_T) << 4U);
  H[477] = -1.0;
  H[478] = 0.0;
  H[479] = 0.0;
  H[480] = 0.0;
  H[481] = (-Subsystem_Fy_B.t2 + -Subsystem_Fy_B.t45_i) + 1.0;
  H[482] = -Subsystem_Fy_B.t16_b + Subsystem_Fy_B.t17_j;
  H[483] = Subsystem_Fy_B.t62;
  H[484] = (-Subsystem_Fy_B.t13_b + -Subsystem_Fy_B.t4_n) + 1.0;
  H[485] = -Subsystem_Fy_B.t7 + Subsystem_Fy_B.t6_l;
  H[486] = Subsystem_Fy_B.t59;
  H[487] = 0.0;
  H[488] = 0.0;
  H[489] = 0.0;
  H[490] = 0.0;
  H[491] = 0.0;
  H[492] = 0.0;
  H[493] = 0.0;
  H[494] = Subsystem_Fy_B.t63;
  H[495] = (-Subsystem_Fy_B.t3_b + -Subsystem_Fy_B.t45_i) + 1.0;
  H[496] = -Subsystem_Fy_B.t14_d + Subsystem_Fy_B.t19_a;
  H[497] = Subsystem_Fy_B.t60;
  H[498] = (-Subsystem_Fy_B.t5_b + -Subsystem_Fy_B.t4_n) + 1.0;
  H[499] = -Subsystem_Fy_B.t9 + Subsystem_Fy_B.t12_h;
  H[500] = 0.0;
  H[501] = 0.0;
  H[502] = 0.0;
  H[503] = 0.0;
  H[504] = 0.0;
  H[505] = 0.0;
  H[506] = 0.0;
  H[507] = -Subsystem_Fy_B.t15_e + Subsystem_Fy_B.t18_f;
  H[508] = Subsystem_Fy_B.t61;
  H[509] = (-Subsystem_Fy_B.t3_b + -Subsystem_Fy_B.t2) + 1.0;
  H[510] = -Subsystem_Fy_B.t8 + Subsystem_Fy_B.t28;
  H[511] = Subsystem_Fy_B.t58;
  H[512] = (-Subsystem_Fy_B.t5_b + -Subsystem_Fy_B.t13_b) + 1.0;
  H[513] = 0.0;
  H[514] = 0.0;
  H[515] = 0.0;
  H[516] = 0.0;
  H[517] = 0.0;
  H[518] = 0.0;
  H[519] = 0.0;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xswap(int32_T n, real_T x_data[],
  int32_T ix0, int32_T iy0)
{
  real_T temp;
  int32_T ix;
  int32_T iy;
  int32_T k;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 0; k < n; k++) {
    temp = x_data[ix];
    x_data[ix] = x_data[iy];
    x_data[iy] = temp;
    ix++;
    iy++;
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
real_T Subsystem_FyModelClass::Subsystem_Fy_xnrm2_dg(int32_T n, const real_T
  x_data[], int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T k;
  int32_T kend;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x_data[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = std::abs(x_data[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xzlarf_f(int32_T m, int32_T n, int32_T
  iv0, real_T tau, real_T C_data[], int32_T ic0, int32_T ldc, real_T work_data[])
{
  int32_T coltop;
  int32_T d;
  int32_T exitg1;
  int32_T ia;
  int32_T iac;
  int32_T ix;
  int32_T jy;
  int32_T lastc;
  int32_T lastv;
  boolean_T exitg2;
  if (tau != 0.0) {
    lastv = m;
    lastc = iv0 + m;
    while ((lastv > 0) && (C_data[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      coltop = lastc * ldc + ic0;
      jy = coltop;
      do {
        exitg1 = 0;
        if (jy <= (coltop + lastv) - 1) {
          if (C_data[jy - 1] != 0.0) {
            exitg1 = 1;
          } else {
            jy++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = -1;
  }

  if (lastv > 0) {
    if (lastc + 1 != 0) {
      for (coltop = 0; coltop <= lastc; coltop++) {
        work_data[coltop] = 0.0;
      }

      coltop = 0;
      jy = ldc * lastc + ic0;
      iac = ic0;
      while (((ldc > 0) && (iac <= jy)) || ((ldc < 0) && (iac >= jy))) {
        ix = iv0;
        Subsystem_Fy_B.c_a = 0.0;
        d = (iac + lastv) - 1;
        for (ia = iac; ia <= d; ia++) {
          Subsystem_Fy_B.c_a += C_data[ia - 1] * C_data[ix - 1];
          ix++;
        }

        work_data[coltop] += Subsystem_Fy_B.c_a;
        coltop++;
        iac += ldc;
      }
    }

    if (!(-tau == 0.0)) {
      coltop = ic0 - 1;
      jy = 0;
      for (iac = 0; iac <= lastc; iac++) {
        if (work_data[jy] != 0.0) {
          Subsystem_Fy_B.c_a = work_data[jy] * -tau;
          ix = iv0;
          d = lastv + coltop;
          for (ia = coltop; ia < d; ia++) {
            C_data[ia] += C_data[ix - 1] * Subsystem_Fy_B.c_a;
            ix++;
          }
        }

        jy++;
        coltop += ldc;
      }
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_qrpf(real_T A_data[], const int32_T
  A_size[2], int32_T m, int32_T n, real_T tau_data[], int32_T jpvt_data[])
{
  int32_T b_j;
  int32_T ix;
  int32_T k;
  int32_T mmi;
  int32_T nmi;
  int32_T pvt;
  Subsystem_Fy_B.ma = A_size[0];
  if (m < n) {
    Subsystem_Fy_B.minmn_a = m;
  } else {
    Subsystem_Fy_B.minmn_a = n;
  }

  mmi = A_size[1];
  if (0 <= mmi - 1) {
    std::memset(&Subsystem_Fy_B.work_data[0], 0, mmi * sizeof(real_T));
  }

  mmi = A_size[1];
  if (0 <= mmi - 1) {
    std::memset(&Subsystem_Fy_B.vn1_data[0], 0, mmi * sizeof(real_T));
  }

  mmi = A_size[1];
  if (0 <= mmi - 1) {
    std::memset(&Subsystem_Fy_B.vn2_data[0], 0, mmi * sizeof(real_T));
  }

  for (b_j = 0; b_j < n; b_j++) {
    Subsystem_Fy_B.vn1_data[b_j] = Subsystem_Fy_xnrm2_dg(m, A_data, b_j *
      Subsystem_Fy_B.ma + 1);
    Subsystem_Fy_B.vn2_data[b_j] = Subsystem_Fy_B.vn1_data[b_j];
  }

  for (b_j = 0; b_j < Subsystem_Fy_B.minmn_a; b_j++) {
    Subsystem_Fy_B.ii = b_j * Subsystem_Fy_B.ma + b_j;
    nmi = n - b_j;
    mmi = m - b_j;
    if (nmi < 1) {
      pvt = -1;
    } else {
      pvt = 0;
      if (nmi > 1) {
        ix = b_j;
        Subsystem_Fy_B.smax_p = std::abs(Subsystem_Fy_B.vn1_data[b_j]);
        for (k = 2; k <= nmi; k++) {
          ix++;
          Subsystem_Fy_B.beta1_p = std::abs(Subsystem_Fy_B.vn1_data[ix]);
          if (Subsystem_Fy_B.beta1_p > Subsystem_Fy_B.smax_p) {
            pvt = k - 1;
            Subsystem_Fy_B.smax_p = Subsystem_Fy_B.beta1_p;
          }
        }
      }
    }

    pvt += b_j;
    if (pvt + 1 != b_j + 1) {
      Subsystem_Fy_xswap(m, A_data, pvt * Subsystem_Fy_B.ma + 1, b_j *
                         Subsystem_Fy_B.ma + 1);
      ix = jpvt_data[pvt];
      jpvt_data[pvt] = jpvt_data[b_j];
      jpvt_data[b_j] = ix;
      Subsystem_Fy_B.vn1_data[pvt] = Subsystem_Fy_B.vn1_data[b_j];
      Subsystem_Fy_B.vn2_data[pvt] = Subsystem_Fy_B.vn2_data[b_j];
    }

    if (b_j + 1 < m) {
      Subsystem_Fy_B.smax_p = A_data[Subsystem_Fy_B.ii];
      tau_data[b_j] = 0.0;
      if (mmi > 0) {
        Subsystem_Fy_B.beta1_p = Subsystem_Fy_xnrm2_dg(mmi - 1, A_data,
          Subsystem_Fy_B.ii + 2);
        if (Subsystem_Fy_B.beta1_p != 0.0) {
          Subsystem_Fy_B.beta1_p = Subsystem_Fy_rt_hypotd_snf
            (A_data[Subsystem_Fy_B.ii], Subsystem_Fy_B.beta1_p);
          if (A_data[Subsystem_Fy_B.ii] >= 0.0) {
            Subsystem_Fy_B.beta1_p = -Subsystem_Fy_B.beta1_p;
          }

          if (std::abs(Subsystem_Fy_B.beta1_p) < 1.0020841800044864E-292) {
            pvt = -1;
            ix = Subsystem_Fy_B.ii + mmi;
            do {
              pvt++;
              for (k = Subsystem_Fy_B.ii + 1; k < ix; k++) {
                A_data[k] *= 9.9792015476736E+291;
              }

              Subsystem_Fy_B.beta1_p *= 9.9792015476736E+291;
              Subsystem_Fy_B.smax_p *= 9.9792015476736E+291;
            } while (!(std::abs(Subsystem_Fy_B.beta1_p) >=
                       1.0020841800044864E-292));

            Subsystem_Fy_B.beta1_p = Subsystem_Fy_rt_hypotd_snf
              (Subsystem_Fy_B.smax_p, Subsystem_Fy_xnrm2_dg(mmi - 1, A_data,
                Subsystem_Fy_B.ii + 2));
            if (Subsystem_Fy_B.smax_p >= 0.0) {
              Subsystem_Fy_B.beta1_p = -Subsystem_Fy_B.beta1_p;
            }

            tau_data[b_j] = (Subsystem_Fy_B.beta1_p - Subsystem_Fy_B.smax_p) /
              Subsystem_Fy_B.beta1_p;
            Subsystem_Fy_B.smax_p = 1.0 / (Subsystem_Fy_B.smax_p -
              Subsystem_Fy_B.beta1_p);
            for (k = Subsystem_Fy_B.ii + 1; k < ix; k++) {
              A_data[k] *= Subsystem_Fy_B.smax_p;
            }

            for (ix = 0; ix <= pvt; ix++) {
              Subsystem_Fy_B.beta1_p *= 1.0020841800044864E-292;
            }

            Subsystem_Fy_B.smax_p = Subsystem_Fy_B.beta1_p;
          } else {
            tau_data[b_j] = (Subsystem_Fy_B.beta1_p - A_data[Subsystem_Fy_B.ii])
              / Subsystem_Fy_B.beta1_p;
            Subsystem_Fy_B.smax_p = 1.0 / (A_data[Subsystem_Fy_B.ii] -
              Subsystem_Fy_B.beta1_p);
            pvt = Subsystem_Fy_B.ii + mmi;
            for (ix = Subsystem_Fy_B.ii + 1; ix < pvt; ix++) {
              A_data[ix] *= Subsystem_Fy_B.smax_p;
            }

            Subsystem_Fy_B.smax_p = Subsystem_Fy_B.beta1_p;
          }
        }
      }

      A_data[Subsystem_Fy_B.ii] = Subsystem_Fy_B.smax_p;
    } else {
      tau_data[b_j] = 0.0;
    }

    if (b_j + 1 < n) {
      Subsystem_Fy_B.smax_p = A_data[Subsystem_Fy_B.ii];
      A_data[Subsystem_Fy_B.ii] = 1.0;
      Subsystem_Fy_xzlarf_f(mmi, nmi - 1, Subsystem_Fy_B.ii + 1, tau_data[b_j],
                            A_data, (Subsystem_Fy_B.ii + Subsystem_Fy_B.ma) + 1,
                            Subsystem_Fy_B.ma, Subsystem_Fy_B.work_data);
      A_data[Subsystem_Fy_B.ii] = Subsystem_Fy_B.smax_p;
    }

    Subsystem_Fy_B.ii = b_j + 1;
    while (Subsystem_Fy_B.ii + 1 <= n) {
      nmi = Subsystem_Fy_B.ii * Subsystem_Fy_B.ma + b_j;
      if (Subsystem_Fy_B.vn1_data[Subsystem_Fy_B.ii] != 0.0) {
        Subsystem_Fy_B.smax_p = std::abs(A_data[nmi]) /
          Subsystem_Fy_B.vn1_data[Subsystem_Fy_B.ii];
        Subsystem_Fy_B.smax_p = 1.0 - Subsystem_Fy_B.smax_p *
          Subsystem_Fy_B.smax_p;
        if (Subsystem_Fy_B.smax_p < 0.0) {
          Subsystem_Fy_B.smax_p = 0.0;
        }

        Subsystem_Fy_B.beta1_p = Subsystem_Fy_B.vn1_data[Subsystem_Fy_B.ii] /
          Subsystem_Fy_B.vn2_data[Subsystem_Fy_B.ii];
        Subsystem_Fy_B.beta1_p = Subsystem_Fy_B.beta1_p * Subsystem_Fy_B.beta1_p
          * Subsystem_Fy_B.smax_p;
        if (Subsystem_Fy_B.beta1_p <= 1.4901161193847656E-8) {
          if (b_j + 1 < m) {
            Subsystem_Fy_B.vn1_data[Subsystem_Fy_B.ii] = Subsystem_Fy_xnrm2_dg
              (mmi - 1, A_data, nmi + 2);
            Subsystem_Fy_B.vn2_data[Subsystem_Fy_B.ii] =
              Subsystem_Fy_B.vn1_data[Subsystem_Fy_B.ii];
          } else {
            Subsystem_Fy_B.vn1_data[Subsystem_Fy_B.ii] = 0.0;
            Subsystem_Fy_B.vn2_data[Subsystem_Fy_B.ii] = 0.0;
          }
        } else {
          Subsystem_Fy_B.vn1_data[Subsystem_Fy_B.ii] *= std::sqrt
            (Subsystem_Fy_B.smax_p);
        }
      }

      Subsystem_Fy_B.ii++;
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xzgeqp3(real_T A_data[], const int32_T
  A_size[2], int32_T m, int32_T n, real_T tau_data[], int32_T *tau_size, int32_T
  jpvt_data[], int32_T jpvt_size[2])
{
  if (A_size[0] < A_size[1]) {
    Subsystem_Fy_B.minmana = A_size[0];
  } else {
    Subsystem_Fy_B.minmana = A_size[1];
  }

  jpvt_size[0] = 1;
  jpvt_size[1] = A_size[1];
  Subsystem_Fy_B.loop_ub_p = A_size[1] - 1;
  if (0 <= Subsystem_Fy_B.loop_ub_p) {
    std::memset(&jpvt_data[0], 0, (Subsystem_Fy_B.loop_ub_p + 1) * sizeof
                (int32_T));
  }

  Subsystem_Fy_B.loop_ub_p = 0;
  while (Subsystem_Fy_B.loop_ub_p <= n - 1) {
    jpvt_data[Subsystem_Fy_B.loop_ub_p] = Subsystem_Fy_B.loop_ub_p + 1;
    Subsystem_Fy_B.loop_ub_p++;
  }

  *tau_size = Subsystem_Fy_B.minmana;
  if (0 <= Subsystem_Fy_B.minmana - 1) {
    std::memset(&tau_data[0], 0, Subsystem_Fy_B.minmana * sizeof(real_T));
  }

  Subsystem_Fy_qrpf(A_data, A_size, m, n, tau_data, jpvt_data);
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_xgetrf(int32_T m, int32_T n, real_T
  A_data[], const int32_T A_size[2], int32_T lda, int32_T ipiv_data[], int32_T
  ipiv_size[2], int32_T *info)
{
  int32_T b_c;
  int32_T b_ix;
  int32_T b_j;
  int32_T c_ix;
  int32_T d;
  int32_T e;
  int32_T ijA;
  int32_T ix;
  int32_T jy;
  int32_T mmj;
  int32_T yk;
  if (m < n) {
    d = m;
  } else {
    d = n;
  }

  ipiv_size[0] = 1;
  ipiv_size[1] = d;
  ipiv_data[0] = 1;
  yk = 1;
  for (mmj = 2; mmj <= d; mmj++) {
    yk++;
    ipiv_data[mmj - 1] = yk;
  }

  *info = 0;
  if (m - 1 < n) {
    d = m - 1;
  } else {
    d = n;
  }

  for (yk = 0; yk < d; yk++) {
    mmj = m - yk;
    b_c = (lda + 1) * yk;
    if (mmj < 1) {
      b_ix = -1;
    } else {
      b_ix = 0;
      if (mmj > 1) {
        ix = b_c;
        Subsystem_Fy_B.smax = std::abs(A_data[b_c]);
        for (jy = 1; jy < mmj; jy++) {
          ix++;
          Subsystem_Fy_B.y = std::abs(A_data[ix]);
          if (Subsystem_Fy_B.y > Subsystem_Fy_B.smax) {
            b_ix = jy;
            Subsystem_Fy_B.smax = Subsystem_Fy_B.y;
          }
        }
      }
    }

    if (A_data[b_c + b_ix] != 0.0) {
      if (b_ix != 0) {
        ix = yk + b_ix;
        ipiv_data[yk] = ix + 1;
        b_ix = yk;
        for (jy = 0; jy < n; jy++) {
          Subsystem_Fy_B.smax = A_data[b_ix];
          A_data[b_ix] = A_data[ix];
          A_data[ix] = Subsystem_Fy_B.smax;
          b_ix += lda;
          ix += lda;
        }
      }

      b_ix = b_c + mmj;
      for (ix = b_c + 1; ix < b_ix; ix++) {
        A_data[ix] /= A_data[b_c];
      }
    } else {
      *info = yk + 1;
    }

    b_ix = n - yk;
    ix = b_c + lda;
    jy = ix;
    for (b_j = 0; b_j <= b_ix - 2; b_j++) {
      Subsystem_Fy_B.smax = A_data[jy];
      if (A_data[jy] != 0.0) {
        c_ix = b_c + 1;
        e = mmj + ix;
        for (ijA = ix + 1; ijA < e; ijA++) {
          A_data[ijA] += A_data[c_ix] * -Subsystem_Fy_B.smax;
          c_ix++;
        }
      }

      jy += lda;
      ix += lda;
    }
  }

  if ((*info == 0) && (m <= n) && (!(A_data[(m + A_size[0] * (m - 1)) - 1] !=
        0.0))) {
    *info = m;
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_lusolve(const real_T A_data[], const
  int32_T A_size[2], const real_T B_data[], const int32_T B_size[2], real_T
  X_data[], int32_T X_size[2])
{
  int32_T b_A;
  int32_T b_k;
  int32_T d_i;
  int32_T loop_ub;
  Subsystem_Fy_B.n = A_size[1];
  Subsystem_Fy_B.b_A_size[0] = A_size[0];
  Subsystem_Fy_B.b_A_size[1] = A_size[1];
  loop_ub = A_size[0] * A_size[1];
  if (0 <= loop_ub - 1) {
    std::memcpy(&Subsystem_Fy_B.b_A_data_c[0], &A_data[0], loop_ub * sizeof
                (real_T));
  }

  Subsystem_Fy_xgetrf(A_size[1], A_size[1], Subsystem_Fy_B.b_A_data_c,
                      Subsystem_Fy_B.b_A_size, A_size[1],
                      Subsystem_Fy_B.ipiv_data, Subsystem_Fy_B.ipiv_size,
                      &Subsystem_Fy_B.jp);
  X_size[0] = 40;
  X_size[1] = B_size[1];
  loop_ub = 40 * B_size[1];
  if (0 <= loop_ub - 1) {
    std::memcpy(&X_data[0], &B_data[0], loop_ub * sizeof(real_T));
  }

  for (b_A = 0; b_A < Subsystem_Fy_B.n; b_A++) {
    Subsystem_Fy_B.b_jBcol = 40 * b_A;
    Subsystem_Fy_B.b_jAcol = Subsystem_Fy_B.n * b_A - 1;
    for (b_k = 0; b_k < b_A; b_k++) {
      Subsystem_Fy_B.b_kBcol = 40 * b_k;
      Subsystem_Fy_B.jp = (b_k + Subsystem_Fy_B.b_jAcol) + 1;
      if (Subsystem_Fy_B.b_A_data_c[Subsystem_Fy_B.jp] != 0.0) {
        for (d_i = 0; d_i < 40; d_i++) {
          loop_ub = d_i + Subsystem_Fy_B.b_jBcol;
          X_data[loop_ub] -= Subsystem_Fy_B.b_A_data_c[Subsystem_Fy_B.jp] *
            X_data[d_i + Subsystem_Fy_B.b_kBcol];
        }
      }
    }

    Subsystem_Fy_B.temp = 1.0 / Subsystem_Fy_B.b_A_data_c[(b_A +
      Subsystem_Fy_B.b_jAcol) + 1];
    for (loop_ub = 0; loop_ub < 40; loop_ub++) {
      Subsystem_Fy_B.jp = loop_ub + Subsystem_Fy_B.b_jBcol;
      X_data[Subsystem_Fy_B.jp] *= Subsystem_Fy_B.temp;
    }
  }

  for (b_A = A_size[1]; b_A > 0; b_A--) {
    Subsystem_Fy_B.b_jBcol = (b_A - 1) * 40;
    Subsystem_Fy_B.b_jAcol = (b_A - 1) * Subsystem_Fy_B.n - 1;
    for (b_k = b_A + 1; b_k <= Subsystem_Fy_B.n; b_k++) {
      Subsystem_Fy_B.b_kBcol = (b_k - 1) * 40;
      Subsystem_Fy_B.jp = b_k + Subsystem_Fy_B.b_jAcol;
      if (Subsystem_Fy_B.b_A_data_c[Subsystem_Fy_B.jp] != 0.0) {
        for (d_i = 0; d_i < 40; d_i++) {
          loop_ub = d_i + Subsystem_Fy_B.b_jBcol;
          X_data[loop_ub] -= Subsystem_Fy_B.b_A_data_c[Subsystem_Fy_B.jp] *
            X_data[d_i + Subsystem_Fy_B.b_kBcol];
        }
      }
    }
  }

  Subsystem_Fy_B.n = A_size[1] - 2;
  while (Subsystem_Fy_B.n + 1 > 0) {
    d_i = Subsystem_Fy_B.ipiv_data[Subsystem_Fy_B.n];
    if (Subsystem_Fy_B.n + 1 != d_i) {
      for (loop_ub = 0; loop_ub < 40; loop_ub++) {
        Subsystem_Fy_B.temp = X_data[loop_ub + 40 * Subsystem_Fy_B.n];
        Subsystem_Fy_B.jp = loop_ub + 40 * (d_i - 1);
        X_data[loop_ub + 40 * Subsystem_Fy_B.n] = X_data[Subsystem_Fy_B.jp];
        X_data[Subsystem_Fy_B.jp] = Subsystem_Fy_B.temp;
      }
    }

    Subsystem_Fy_B.n--;
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_mrdiv(const real_T A_data[], const
  int32_T A_size[2], const real_T B_data[], const int32_T B_size[2], real_T
  Y_data[], int32_T Y_size[2])
{
  if (B_size[0] == B_size[1]) {
    Subsystem_Fy_lusolve(B_data, B_size, A_data, A_size, Y_data, Y_size);
  } else {
    Subsystem_Fy_B.b_A_size_idx_0 = B_size[1];
    Subsystem_Fy_B.minmn = B_size[0];
    for (Subsystem_Fy_B.maxmn = 0; Subsystem_Fy_B.maxmn < Subsystem_Fy_B.minmn;
         Subsystem_Fy_B.maxmn++) {
      Subsystem_Fy_B.rankA = B_size[1];
      for (Subsystem_Fy_B.c_A = 0; Subsystem_Fy_B.c_A < Subsystem_Fy_B.rankA;
           Subsystem_Fy_B.c_A++) {
        Subsystem_Fy_B.b_A_data[Subsystem_Fy_B.c_A +
          Subsystem_Fy_B.b_A_size_idx_0 * Subsystem_Fy_B.maxmn] =
          B_data[Subsystem_Fy_B.maxmn + B_size[0] * Subsystem_Fy_B.c_A];
      }
    }

    Subsystem_Fy_B.c_A_size[0] = B_size[1];
    Subsystem_Fy_B.c_A_size[1] = B_size[0];
    Subsystem_Fy_B.minmn = B_size[0] * B_size[1];
    if (0 <= Subsystem_Fy_B.minmn - 1) {
      std::memcpy(&Subsystem_Fy_B.c_A_data[0], &Subsystem_Fy_B.b_A_data[0],
                  Subsystem_Fy_B.minmn * sizeof(real_T));
    }

    Subsystem_Fy_xzgeqp3(Subsystem_Fy_B.c_A_data, Subsystem_Fy_B.c_A_size,
                         B_size[1], B_size[0], Subsystem_Fy_B.tau_data,
                         &Subsystem_Fy_B.tau_size, Subsystem_Fy_B.jpvt_data,
                         Subsystem_Fy_B.jpvt_size);
    Subsystem_Fy_B.rankA = 0;
    if (Subsystem_Fy_B.c_A_size[0] < Subsystem_Fy_B.c_A_size[1]) {
      Subsystem_Fy_B.minmn = Subsystem_Fy_B.c_A_size[0];
      Subsystem_Fy_B.maxmn = Subsystem_Fy_B.c_A_size[1];
    } else {
      Subsystem_Fy_B.minmn = Subsystem_Fy_B.c_A_size[1];
      Subsystem_Fy_B.maxmn = Subsystem_Fy_B.c_A_size[0];
    }

    Subsystem_Fy_B.tol = 2.2204460492503131E-15 * static_cast<real_T>
      (Subsystem_Fy_B.maxmn) * std::abs(Subsystem_Fy_B.c_A_data[0]);
    while ((Subsystem_Fy_B.rankA < Subsystem_Fy_B.minmn) && (!(std::abs
             (Subsystem_Fy_B.c_A_data[Subsystem_Fy_B.rankA +
              Subsystem_Fy_B.c_A_size[0] * Subsystem_Fy_B.rankA]) <=
             Subsystem_Fy_B.tol))) {
      Subsystem_Fy_B.rankA++;
    }

    Subsystem_Fy_B.b_A_size_idx_0 = static_cast<int8_T>(Subsystem_Fy_B.c_A_size
      [1]);
    Subsystem_Fy_B.minmn = static_cast<int8_T>(Subsystem_Fy_B.c_A_size[1]) * 40
      - 1;
    if (0 <= Subsystem_Fy_B.minmn) {
      std::memset(&Subsystem_Fy_B.b_Y_data[0], 0, (Subsystem_Fy_B.minmn + 1) *
                  sizeof(real_T));
    }

    Subsystem_Fy_B.b_B_size_idx_0 = A_size[1];
    Subsystem_Fy_B.minmn = A_size[1];
    for (Subsystem_Fy_B.maxmn = 0; Subsystem_Fy_B.maxmn < 40;
         Subsystem_Fy_B.maxmn++) {
      for (Subsystem_Fy_B.c_A = 0; Subsystem_Fy_B.c_A < Subsystem_Fy_B.minmn;
           Subsystem_Fy_B.c_A++) {
        Subsystem_Fy_B.b_B_data[Subsystem_Fy_B.c_A +
          Subsystem_Fy_B.b_B_size_idx_0 * Subsystem_Fy_B.maxmn] =
          A_data[Subsystem_Fy_B.maxmn + 40 * Subsystem_Fy_B.c_A];
      }
    }

    Subsystem_Fy_B.minmn = Subsystem_Fy_B.c_A_size[0];
    if (Subsystem_Fy_B.c_A_size[0] < Subsystem_Fy_B.c_A_size[1]) {
      Subsystem_Fy_B.maxmn = Subsystem_Fy_B.c_A_size[0];
    } else {
      Subsystem_Fy_B.maxmn = Subsystem_Fy_B.c_A_size[1];
    }

    Subsystem_Fy_B.c_A = 0;
    while (Subsystem_Fy_B.c_A <= Subsystem_Fy_B.maxmn - 1) {
      if (Subsystem_Fy_B.tau_data[Subsystem_Fy_B.c_A] != 0.0) {
        for (Subsystem_Fy_B.b_Y_data_tmp_tmp = 0;
             Subsystem_Fy_B.b_Y_data_tmp_tmp < 40;
             Subsystem_Fy_B.b_Y_data_tmp_tmp++) {
          Subsystem_Fy_B.b_Y_data_tmp = Subsystem_Fy_B.b_B_size_idx_0 *
            Subsystem_Fy_B.b_Y_data_tmp_tmp;
          Subsystem_Fy_B.b_Y_data_tmp_p = Subsystem_Fy_B.c_A +
            Subsystem_Fy_B.b_Y_data_tmp;
          Subsystem_Fy_B.tol_tmp =
            Subsystem_Fy_B.b_B_data[Subsystem_Fy_B.b_Y_data_tmp_p];
          Subsystem_Fy_B.tol = Subsystem_Fy_B.tol_tmp;
          Subsystem_Fy_B.b_Y_data_tmp_k = Subsystem_Fy_B.c_A + 1;
          while (Subsystem_Fy_B.b_Y_data_tmp_k + 1 <= Subsystem_Fy_B.minmn) {
            Subsystem_Fy_B.tol +=
              Subsystem_Fy_B.c_A_data[Subsystem_Fy_B.b_Y_data_tmp_k +
              Subsystem_Fy_B.c_A_size[0] * Subsystem_Fy_B.c_A] *
              Subsystem_Fy_B.b_B_data[Subsystem_Fy_B.b_Y_data_tmp_k +
              Subsystem_Fy_B.b_Y_data_tmp];
            Subsystem_Fy_B.b_Y_data_tmp_k++;
          }

          Subsystem_Fy_B.tol *= Subsystem_Fy_B.tau_data[Subsystem_Fy_B.c_A];
          if (Subsystem_Fy_B.tol != 0.0) {
            Subsystem_Fy_B.b_B_data[Subsystem_Fy_B.b_Y_data_tmp_p] =
              Subsystem_Fy_B.tol_tmp - Subsystem_Fy_B.tol;
            Subsystem_Fy_B.b_Y_data_tmp_k = Subsystem_Fy_B.c_A + 1;
            while (Subsystem_Fy_B.b_Y_data_tmp_k + 1 <= Subsystem_Fy_B.minmn) {
              Subsystem_Fy_B.b_Y_data_tmp_p = Subsystem_Fy_B.b_Y_data_tmp_k +
                Subsystem_Fy_B.b_Y_data_tmp;
              Subsystem_Fy_B.b_B_data[Subsystem_Fy_B.b_Y_data_tmp_p] -=
                Subsystem_Fy_B.c_A_data[Subsystem_Fy_B.b_Y_data_tmp_k +
                Subsystem_Fy_B.c_A_size[0] * Subsystem_Fy_B.c_A] *
                Subsystem_Fy_B.tol;
              Subsystem_Fy_B.b_Y_data_tmp_k++;
            }
          }
        }
      }

      Subsystem_Fy_B.c_A++;
    }

    for (Subsystem_Fy_B.minmn = 0; Subsystem_Fy_B.minmn < 40;
         Subsystem_Fy_B.minmn++) {
      Subsystem_Fy_B.maxmn = 0;
      while (Subsystem_Fy_B.maxmn <= Subsystem_Fy_B.rankA - 1) {
        Subsystem_Fy_B.b_Y_data[(Subsystem_Fy_B.jpvt_data[Subsystem_Fy_B.maxmn]
          + Subsystem_Fy_B.b_A_size_idx_0 * Subsystem_Fy_B.minmn) - 1] =
          Subsystem_Fy_B.b_B_data[Subsystem_Fy_B.maxmn +
          Subsystem_Fy_B.b_B_size_idx_0 * Subsystem_Fy_B.minmn];
        Subsystem_Fy_B.maxmn++;
      }

      Subsystem_Fy_B.maxmn = Subsystem_Fy_B.rankA - 1;
      while (Subsystem_Fy_B.maxmn + 1 > 0) {
        Subsystem_Fy_B.b_Y_data_tmp_tmp = Subsystem_Fy_B.b_A_size_idx_0 *
          Subsystem_Fy_B.minmn;
        Subsystem_Fy_B.b_Y_data_tmp =
          (Subsystem_Fy_B.jpvt_data[Subsystem_Fy_B.maxmn] +
           Subsystem_Fy_B.b_Y_data_tmp_tmp) - 1;
        Subsystem_Fy_B.b_Y_data_tmp_k = Subsystem_Fy_B.c_A_size[0] *
          Subsystem_Fy_B.maxmn;
        Subsystem_Fy_B.b_Y_data[Subsystem_Fy_B.b_Y_data_tmp] /=
          Subsystem_Fy_B.c_A_data[Subsystem_Fy_B.maxmn +
          Subsystem_Fy_B.b_Y_data_tmp_k];
        Subsystem_Fy_B.c_A = 0;
        while (Subsystem_Fy_B.c_A <= Subsystem_Fy_B.maxmn - 1) {
          Subsystem_Fy_B.b_Y_data_tmp_p =
            (Subsystem_Fy_B.jpvt_data[Subsystem_Fy_B.c_A] +
             Subsystem_Fy_B.b_Y_data_tmp_tmp) - 1;
          Subsystem_Fy_B.b_Y_data[Subsystem_Fy_B.b_Y_data_tmp_p] -=
            Subsystem_Fy_B.b_Y_data[Subsystem_Fy_B.b_Y_data_tmp] *
            Subsystem_Fy_B.c_A_data[Subsystem_Fy_B.c_A +
            Subsystem_Fy_B.b_Y_data_tmp_k];
          Subsystem_Fy_B.c_A++;
        }

        Subsystem_Fy_B.maxmn--;
      }
    }

    Y_size[0] = 40;
    Y_size[1] = static_cast<int8_T>(Subsystem_Fy_B.c_A_size[1]);
    for (Subsystem_Fy_B.maxmn = 0; Subsystem_Fy_B.maxmn <
         Subsystem_Fy_B.b_A_size_idx_0; Subsystem_Fy_B.maxmn++) {
      for (Subsystem_Fy_B.c_A = 0; Subsystem_Fy_B.c_A < 40; Subsystem_Fy_B.c_A++)
      {
        Y_data[Subsystem_Fy_B.c_A + 40 * Subsystem_Fy_B.maxmn] =
          Subsystem_Fy_B.b_Y_data[Subsystem_Fy_B.maxmn +
          Subsystem_Fy_B.b_A_size_idx_0 * Subsystem_Fy_B.c_A];
      }
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsyst_eskf_cdprosthesis_reset(const real_T in1[42],
  const real_T in2[40], real_T Xt[42], real_T reset_dx[1600])
{
  real_T t69;
  real_T t70;
  real_T t71;
  real_T t72;
  real_T t73;
  real_T t76;
  real_T t77;
  Subsystem_Fy_B.t52 = ((in2[24] * in1[26] / 2.0 + -in1[25]) + in2[25] * in1[27]
                        / 2.0) + in2[26] * in1[28] / 2.0;
  Subsystem_Fy_B.t56 = ((in2[6] * in1[7] / 2.0 + -in1[6]) + in2[7] * in1[8] /
                        2.0) + in2[8] * in1[9] / 2.0;
  Subsystem_Fy_B.t46 = ((in2[6] * in1[6] / 2.0 + in1[7]) + in2[7] * in1[9] / 2.0)
    + -(in2[8] * in1[8] / 2.0);
  Subsystem_Fy_B.t47 = ((in1[6] * in2[7] / 2.0 + in1[8]) + in1[7] * in2[8] / 2.0)
    + -(in2[6] * in1[9] / 2.0);
  Subsystem_Fy_B.t48 = ((in2[6] * in1[8] / 2.0 + in1[9]) + in1[6] * in2[8] / 2.0)
    + -(in2[7] * in1[7] / 2.0);
  Subsystem_Fy_B.t49_i = ((in2[24] * in1[25] / 2.0 + in1[26]) + in2[25] * in1[28]
    / 2.0) + -(in2[26] * in1[27] / 2.0);
  Subsystem_Fy_B.t50_f = ((in2[25] * in1[25] / 2.0 + in1[27]) + in2[26] * in1[26]
    / 2.0) + -(in2[24] * in1[28] / 2.0);
  Subsystem_Fy_B.t51_i = ((in2[24] * in1[27] / 2.0 + in1[28]) + in1[25] * in2[26]
    / 2.0) + -(in2[25] * in1[26] / 2.0);
  Subsystem_Fy_B.t60_f = ((in2[6] * in2[6] / 4.0 + in2[7] * in2[7] / 4.0) + in2
    [8] * in2[8] / 4.0) + 1.0;
  Subsystem_Fy_B.t61_g = ((in2[24] * in2[24] / 4.0 + in2[25] * in2[25] / 4.0) +
    in2[26] * in2[26] / 4.0) + 1.0;
  Subsystem_Fy_B.t60_f = 1.0 / std::sqrt(Subsystem_Fy_B.t60_f *
    Subsystem_Fy_B.t60_f);
  Subsystem_Fy_B.t61_g = 1.0 / std::sqrt(Subsystem_Fy_B.t61_g *
    Subsystem_Fy_B.t61_g);
  Subsystem_Fy_B.t68_c = in2[6] * Subsystem_Fy_B.t60_f / 2.0;
  t69 = in2[7] * Subsystem_Fy_B.t60_f / 2.0;
  t70 = in2[8] * Subsystem_Fy_B.t60_f / 2.0;
  t71 = in2[24] * Subsystem_Fy_B.t61_g / 2.0;
  t72 = in2[25] * Subsystem_Fy_B.t61_g / 2.0;
  t73 = in2[26] * Subsystem_Fy_B.t61_g / 2.0;
  t76 = 1.0 / std::sqrt(((Subsystem_Fy_B.t49_i * Subsystem_Fy_B.t49_i +
    Subsystem_Fy_B.t50_f * Subsystem_Fy_B.t50_f) + Subsystem_Fy_B.t51_i *
    Subsystem_Fy_B.t51_i) + Subsystem_Fy_B.t52 * Subsystem_Fy_B.t52);
  t77 = 1.0 / std::sqrt(((Subsystem_Fy_B.t46 * Subsystem_Fy_B.t46 +
    Subsystem_Fy_B.t47 * Subsystem_Fy_B.t47) + Subsystem_Fy_B.t48 *
    Subsystem_Fy_B.t48) + Subsystem_Fy_B.t56 * Subsystem_Fy_B.t56);
  Xt[0] = in1[0] + in2[0];
  Xt[1] = in1[1] + in2[1];
  Xt[2] = in1[2] + in2[2];
  Xt[3] = in1[3] + in2[3];
  Xt[4] = in1[4] + in2[4];
  Xt[5] = in1[5] + in2[5];
  Xt[6] = -Subsystem_Fy_B.t56 * t77;
  Xt[7] = Subsystem_Fy_B.t46 * t77;
  Xt[8] = Subsystem_Fy_B.t47 * t77;
  Xt[9] = Subsystem_Fy_B.t48 * t77;
  Xt[10] = in2[9] + in1[10];
  Xt[11] = in2[10] + in1[11];
  Xt[12] = in2[11] + in1[12];
  Xt[13] = in2[12] + in1[13];
  Xt[14] = in2[13] + in1[14];
  Xt[15] = in2[14] + in1[15];
  Xt[16] = in2[15] + in1[16];
  Xt[17] = in2[16] + in1[17];
  Xt[18] = in2[17] + in1[18];
  Xt[19] = in2[18] + in1[19];
  Xt[20] = in2[19] + in1[20];
  Xt[21] = in2[20] + in1[21];
  Xt[22] = in2[21] + in1[22];
  Xt[23] = in2[22] + in1[23];
  Xt[24] = in2[23] + in1[24];
  Xt[25] = -Subsystem_Fy_B.t52 * t76;
  Xt[26] = Subsystem_Fy_B.t49_i * t76;
  Xt[27] = Subsystem_Fy_B.t50_f * t76;
  Xt[28] = Subsystem_Fy_B.t51_i * t76;
  Xt[29] = in2[27] + in1[29];
  Xt[30] = in2[28] + in1[30];
  Xt[31] = in2[29] + in1[31];
  Xt[32] = in2[30] + in1[32];
  Xt[33] = in2[31] + in1[33];
  Xt[34] = in2[32] + in1[34];
  Xt[35] = in2[33] + in1[35];
  Xt[36] = in2[34] + in1[36];
  Xt[37] = in2[35] + in1[37];
  Xt[38] = in2[36] + in1[38];
  Xt[39] = in2[37] + in1[39];
  Xt[40] = in2[38] + in1[40];
  Xt[41] = in2[39] + in1[41];
  reset_dx[0] = 1.0;
  std::memset(&reset_dx[1], 0, 40U * sizeof(real_T));
  reset_dx[41] = 1.0;
  std::memset(&reset_dx[42], 0, 40U * sizeof(real_T));
  reset_dx[82] = 1.0;
  std::memset(&reset_dx[83], 0, 40U * sizeof(real_T));
  reset_dx[123] = 1.0;
  std::memset(&reset_dx[124], 0, 40U * sizeof(real_T));
  reset_dx[164] = 1.0;
  std::memset(&reset_dx[165], 0, 40U * sizeof(real_T));
  reset_dx[205] = 1.0;
  std::memset(&reset_dx[206], 0, 40U * sizeof(real_T));
  reset_dx[246] = Subsystem_Fy_B.t60_f;
  reset_dx[247] = -t70;
  reset_dx[248] = t69;
  std::memset(&reset_dx[249], 0, 37U * sizeof(real_T));
  reset_dx[286] = t70;
  reset_dx[287] = Subsystem_Fy_B.t60_f;
  reset_dx[288] = -Subsystem_Fy_B.t68_c;
  std::memset(&reset_dx[289], 0, 37U * sizeof(real_T));
  reset_dx[326] = -t69;
  reset_dx[327] = Subsystem_Fy_B.t68_c;
  reset_dx[328] = Subsystem_Fy_B.t60_f;
  std::memset(&reset_dx[329], 0, 40U * sizeof(real_T));
  reset_dx[369] = 1.0;
  std::memset(&reset_dx[370], 0, 40U * sizeof(real_T));
  reset_dx[410] = 1.0;
  std::memset(&reset_dx[411], 0, 40U * sizeof(real_T));
  reset_dx[451] = 1.0;
  std::memset(&reset_dx[452], 0, 40U * sizeof(real_T));
  reset_dx[492] = 1.0;
  std::memset(&reset_dx[493], 0, 40U * sizeof(real_T));
  reset_dx[533] = 1.0;
  std::memset(&reset_dx[534], 0, 40U * sizeof(real_T));
  reset_dx[574] = 1.0;
  std::memset(&reset_dx[575], 0, 40U * sizeof(real_T));
  reset_dx[615] = 1.0;
  std::memset(&reset_dx[616], 0, 40U * sizeof(real_T));
  reset_dx[656] = 1.0;
  std::memset(&reset_dx[657], 0, 40U * sizeof(real_T));
  reset_dx[697] = 1.0;
  std::memset(&reset_dx[698], 0, 40U * sizeof(real_T));
  reset_dx[738] = 1.0;
  std::memset(&reset_dx[739], 0, 40U * sizeof(real_T));
  reset_dx[779] = 1.0;
  std::memset(&reset_dx[780], 0, 40U * sizeof(real_T));
  reset_dx[820] = 1.0;
  std::memset(&reset_dx[821], 0, 40U * sizeof(real_T));
  reset_dx[861] = 1.0;
  std::memset(&reset_dx[862], 0, 40U * sizeof(real_T));
  reset_dx[902] = 1.0;
  std::memset(&reset_dx[903], 0, 40U * sizeof(real_T));
  reset_dx[943] = 1.0;
  std::memset(&reset_dx[944], 0, 40U * sizeof(real_T));
  reset_dx[984] = Subsystem_Fy_B.t61_g;
  reset_dx[985] = -t73;
  reset_dx[986] = t72;
  std::memset(&reset_dx[987], 0, 37U * sizeof(real_T));
  reset_dx[1024] = t73;
  reset_dx[1025] = Subsystem_Fy_B.t61_g;
  reset_dx[1026] = -t71;
  std::memset(&reset_dx[1027], 0, 37U * sizeof(real_T));
  reset_dx[1064] = -t72;
  reset_dx[1065] = t71;
  reset_dx[1066] = Subsystem_Fy_B.t61_g;
  std::memset(&reset_dx[1067], 0, 40U * sizeof(real_T));
  reset_dx[1107] = 1.0;
  std::memset(&reset_dx[1108], 0, 40U * sizeof(real_T));
  reset_dx[1148] = 1.0;
  std::memset(&reset_dx[1149], 0, 40U * sizeof(real_T));
  reset_dx[1189] = 1.0;
  std::memset(&reset_dx[1190], 0, 40U * sizeof(real_T));
  reset_dx[1230] = 1.0;
  std::memset(&reset_dx[1231], 0, 40U * sizeof(real_T));
  reset_dx[1271] = 1.0;
  std::memset(&reset_dx[1272], 0, 40U * sizeof(real_T));
  reset_dx[1312] = 1.0;
  std::memset(&reset_dx[1313], 0, 40U * sizeof(real_T));
  reset_dx[1353] = 1.0;
  std::memset(&reset_dx[1354], 0, 40U * sizeof(real_T));
  reset_dx[1394] = 1.0;
  std::memset(&reset_dx[1395], 0, 40U * sizeof(real_T));
  reset_dx[1435] = 1.0;
  std::memset(&reset_dx[1436], 0, 40U * sizeof(real_T));
  reset_dx[1476] = 1.0;
  std::memset(&reset_dx[1477], 0, 40U * sizeof(real_T));
  reset_dx[1517] = 1.0;
  std::memset(&reset_dx[1518], 0, 40U * sizeof(real_T));
  reset_dx[1558] = 1.0;
  std::memset(&reset_dx[1559], 0, 40U * sizeof(real_T));
  reset_dx[1599] = 1.0;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_mtimes(const real_T A_data[], const
  int32_T A_size[2], const real_T B_data[], const int32_T B_size[2], real_T C
  [1600])
{
  real_T bkj;
  int32_T C_tmp;
  int32_T aoffset;
  int32_T b_i;
  int32_T boffset;
  int32_T coffset;
  int32_T i;
  int32_T j;
  for (j = 0; j < 40; j++) {
    coffset = j * 40;
    boffset = j * B_size[0];
    std::memset(&C[coffset], 0, 40U * sizeof(real_T));
    for (i = 0; i < A_size[1]; i++) {
      aoffset = i * 40;
      bkj = B_data[boffset + i];
      for (b_i = 0; b_i < 40; b_i++) {
        C_tmp = coffset + b_i;
        C[C_tmp] += A_data[aoffset + b_i] * bkj;
      }
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_correct(const real_T get_NominalState
  [42], real_T errorCov_[1600], const real_T res_data[], const real_T
  measurementCov_data[], const int32_T measurementCov_size[2], const real_T
  H_data[], const int32_T H_size[2], real_T x_post[42])
{
  static const int8_T b_I[1600] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  Subsystem_Fy_B.m = H_size[0];
  Subsystem_Fy_B.y_size_idx_0 = H_size[0];
  for (Subsystem_Fy_B.b_m_c = 0; Subsystem_Fy_B.b_m_c < 40; Subsystem_Fy_B.b_m_c
       ++) {
    Subsystem_Fy_B.coffset = Subsystem_Fy_B.b_m_c * Subsystem_Fy_B.m;
    Subsystem_Fy_B.boffset = Subsystem_Fy_B.b_m_c * 40;
    Subsystem_Fy_B.i_h = 0;
    while (Subsystem_Fy_B.i_h <= Subsystem_Fy_B.m - 1) {
      Subsystem_Fy_B.y_data[Subsystem_Fy_B.coffset + Subsystem_Fy_B.i_h] = 0.0;
      Subsystem_Fy_B.i_h++;
    }

    for (Subsystem_Fy_B.c_aoffset = 0; Subsystem_Fy_B.c_aoffset < 40;
         Subsystem_Fy_B.c_aoffset++) {
      Subsystem_Fy_B.aoffset = Subsystem_Fy_B.c_aoffset * H_size[0];
      Subsystem_Fy_B.bkj = errorCov_[Subsystem_Fy_B.boffset +
        Subsystem_Fy_B.c_aoffset];
      Subsystem_Fy_B.b_i = 1;
      while (Subsystem_Fy_B.b_i - 1 <= Subsystem_Fy_B.m - 1) {
        Subsystem_Fy_B.i_h = (Subsystem_Fy_B.coffset + Subsystem_Fy_B.b_i) - 1;
        Subsystem_Fy_B.y_data[Subsystem_Fy_B.i_h] += H_data
          [(Subsystem_Fy_B.aoffset + Subsystem_Fy_B.b_i) - 1] *
          Subsystem_Fy_B.bkj;
        Subsystem_Fy_B.b_i++;
      }
    }
  }

  Subsystem_Fy_B.m = 0;
  while (Subsystem_Fy_B.m <= H_size[0] - 1) {
    Subsystem_Fy_B.coffset = Subsystem_Fy_B.m * Subsystem_Fy_B.y_size_idx_0;
    Subsystem_Fy_B.i_h = 0;
    while (Subsystem_Fy_B.i_h <= Subsystem_Fy_B.y_size_idx_0 - 1) {
      Subsystem_Fy_B.b_data[Subsystem_Fy_B.coffset + Subsystem_Fy_B.i_h] = 0.0;
      Subsystem_Fy_B.i_h++;
    }

    for (Subsystem_Fy_B.boffset = 0; Subsystem_Fy_B.boffset < 40;
         Subsystem_Fy_B.boffset++) {
      Subsystem_Fy_B.c_aoffset = Subsystem_Fy_B.boffset *
        Subsystem_Fy_B.y_size_idx_0;
      Subsystem_Fy_B.bkj = H_data[Subsystem_Fy_B.boffset * H_size[0] +
        Subsystem_Fy_B.m];
      Subsystem_Fy_B.aoffset = 1;
      while (Subsystem_Fy_B.aoffset - 1 <= Subsystem_Fy_B.y_size_idx_0 - 1) {
        Subsystem_Fy_B.i_h = (Subsystem_Fy_B.coffset + Subsystem_Fy_B.aoffset) -
          1;
        Subsystem_Fy_B.b_data[Subsystem_Fy_B.i_h] += Subsystem_Fy_B.y_data
          [(Subsystem_Fy_B.c_aoffset + Subsystem_Fy_B.aoffset) - 1] *
          Subsystem_Fy_B.bkj;
        Subsystem_Fy_B.aoffset++;
      }
    }

    Subsystem_Fy_B.m++;
  }

  Subsystem_Fy_B.b_y_size[0] = 40;
  Subsystem_Fy_B.b_y_size[1] = H_size[0];
  Subsystem_Fy_B.b_m_c = 0;
  while (Subsystem_Fy_B.b_m_c <= H_size[0] - 1) {
    Subsystem_Fy_B.m = Subsystem_Fy_B.b_m_c * 40;
    std::memset(&Subsystem_Fy_B.y_data[Subsystem_Fy_B.m], 0, 40U * sizeof(real_T));
    for (Subsystem_Fy_B.coffset = 0; Subsystem_Fy_B.coffset < 40;
         Subsystem_Fy_B.coffset++) {
      Subsystem_Fy_B.c_aoffset = Subsystem_Fy_B.coffset * 40;
      Subsystem_Fy_B.bkj = H_data[Subsystem_Fy_B.coffset * H_size[0] +
        Subsystem_Fy_B.b_m_c];
      for (Subsystem_Fy_B.boffset = 0; Subsystem_Fy_B.boffset < 40;
           Subsystem_Fy_B.boffset++) {
        Subsystem_Fy_B.i_h = Subsystem_Fy_B.m + Subsystem_Fy_B.boffset;
        Subsystem_Fy_B.y_data[Subsystem_Fy_B.i_h] +=
          errorCov_[Subsystem_Fy_B.c_aoffset + Subsystem_Fy_B.boffset] *
          Subsystem_Fy_B.bkj;
      }
    }

    Subsystem_Fy_B.b_m_c++;
  }

  Subsystem_Fy_B.b_size[0] = H_size[0];
  Subsystem_Fy_B.b_size[1] = H_size[0];
  Subsystem_Fy_B.coffset = H_size[0] * H_size[0];
  for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < Subsystem_Fy_B.coffset;
       Subsystem_Fy_B.i_h++) {
    Subsystem_Fy_B.b_data_p[Subsystem_Fy_B.i_h] =
      Subsystem_Fy_B.b_data[Subsystem_Fy_B.i_h] +
      measurementCov_data[Subsystem_Fy_B.i_h];
  }

  Subsystem_Fy_mrdiv(Subsystem_Fy_B.y_data, Subsystem_Fy_B.b_y_size,
                     Subsystem_Fy_B.b_data_p, Subsystem_Fy_B.b_size,
                     Subsystem_Fy_B.K_data, Subsystem_Fy_B.K_size);
  std::memset(&Subsystem_Fy_B.dx[0], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b_m_c = 0;
  while (Subsystem_Fy_B.b_m_c <= Subsystem_Fy_B.K_size[1] - 1) {
    Subsystem_Fy_B.m = Subsystem_Fy_B.b_m_c * 40;
    for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 40; Subsystem_Fy_B.i_h++)
    {
      Subsystem_Fy_B.dx[Subsystem_Fy_B.i_h] +=
        Subsystem_Fy_B.K_data[Subsystem_Fy_B.m + Subsystem_Fy_B.i_h] *
        res_data[Subsystem_Fy_B.b_m_c];
    }

    Subsystem_Fy_B.b_m_c++;
  }

  Subsystem_Fy_B.b_y_size[1] = measurementCov_size[1];
  Subsystem_Fy_B.b_m_c = 0;
  while (Subsystem_Fy_B.b_m_c <= measurementCov_size[1] - 1) {
    Subsystem_Fy_B.m = Subsystem_Fy_B.b_m_c * 40;
    Subsystem_Fy_B.coffset = Subsystem_Fy_B.b_m_c * measurementCov_size[0];
    std::memset(&Subsystem_Fy_B.y_data[Subsystem_Fy_B.m], 0, 40U * sizeof(real_T));
    Subsystem_Fy_B.c_aoffset = 0;
    while (Subsystem_Fy_B.c_aoffset <= Subsystem_Fy_B.K_size[1] - 1) {
      Subsystem_Fy_B.aoffset = Subsystem_Fy_B.c_aoffset * 40;
      Subsystem_Fy_B.bkj = measurementCov_data[Subsystem_Fy_B.coffset +
        Subsystem_Fy_B.c_aoffset];
      for (Subsystem_Fy_B.boffset = 0; Subsystem_Fy_B.boffset < 40;
           Subsystem_Fy_B.boffset++) {
        Subsystem_Fy_B.i_h = Subsystem_Fy_B.m + Subsystem_Fy_B.boffset;
        Subsystem_Fy_B.y_data[Subsystem_Fy_B.i_h] +=
          Subsystem_Fy_B.K_data[Subsystem_Fy_B.aoffset + Subsystem_Fy_B.boffset]
          * Subsystem_Fy_B.bkj;
      }

      Subsystem_Fy_B.c_aoffset++;
    }

    Subsystem_Fy_B.b_m_c++;
  }

  for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 40; Subsystem_Fy_B.i_h++) {
    Subsystem_Fy_B.b_m_c = Subsystem_Fy_B.i_h * 40;
    std::memset(&Subsystem_Fy_B.c_y[Subsystem_Fy_B.b_m_c], 0, 40U * sizeof
                (real_T));
    Subsystem_Fy_B.coffset = 0;
    while (Subsystem_Fy_B.coffset <= Subsystem_Fy_B.b_y_size[1] - 1) {
      Subsystem_Fy_B.boffset = Subsystem_Fy_B.coffset * 40;
      Subsystem_Fy_B.bkj = Subsystem_Fy_B.K_data[Subsystem_Fy_B.coffset * 40 +
        Subsystem_Fy_B.i_h];
      for (Subsystem_Fy_B.m = 0; Subsystem_Fy_B.m < 40; Subsystem_Fy_B.m++) {
        Subsystem_Fy_B.c_aoffset = Subsystem_Fy_B.b_m_c + Subsystem_Fy_B.m;
        Subsystem_Fy_B.c_y[Subsystem_Fy_B.c_aoffset] +=
          Subsystem_Fy_B.y_data[Subsystem_Fy_B.boffset + Subsystem_Fy_B.m] *
          Subsystem_Fy_B.bkj;
      }

      Subsystem_Fy_B.coffset++;
    }
  }

  Subsyst_eskf_cdprosthesis_reset(get_NominalState, Subsystem_Fy_B.dx, x_post,
    Subsystem_Fy_B.G);
  Subsystem_Fy_mtimes(Subsystem_Fy_B.K_data, Subsystem_Fy_B.K_size, H_data,
                      H_size, Subsystem_Fy_B.dv1);
  for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 1600; Subsystem_Fy_B.i_h++)
  {
    Subsystem_Fy_B.b_m_c = b_I[Subsystem_Fy_B.i_h];
    Subsystem_Fy_B.G_k[Subsystem_Fy_B.i_h] = static_cast<real_T>
      (Subsystem_Fy_B.b_m_c) - Subsystem_Fy_B.dv1[Subsystem_Fy_B.i_h];
    Subsystem_Fy_B.iv[Subsystem_Fy_B.i_h] = static_cast<int8_T>
      (Subsystem_Fy_B.b_m_c);
  }

  for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 40; Subsystem_Fy_B.i_h++) {
    std::memset(&Subsystem_Fy_B.dv2[Subsystem_Fy_B.i_h * 40], 0, 40U * sizeof
                (real_T));
  }

  for (Subsystem_Fy_B.b_m_c = 0; Subsystem_Fy_B.b_m_c < 40; Subsystem_Fy_B.b_m_c
       ++) {
    for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 40; Subsystem_Fy_B.i_h++)
    {
      for (Subsystem_Fy_B.m = 0; Subsystem_Fy_B.m < 40; Subsystem_Fy_B.m++) {
        Subsystem_Fy_B.coffset = 40 * Subsystem_Fy_B.i_h + Subsystem_Fy_B.b_m_c;
        Subsystem_Fy_B.dv2[Subsystem_Fy_B.coffset] += Subsystem_Fy_B.G_k[40 *
          Subsystem_Fy_B.m + Subsystem_Fy_B.b_m_c] * errorCov_[40 *
          Subsystem_Fy_B.i_h + Subsystem_Fy_B.m];
      }

      Subsystem_Fy_B.m = 40 * Subsystem_Fy_B.b_m_c + Subsystem_Fy_B.i_h;
      Subsystem_Fy_B.dv3[Subsystem_Fy_B.b_m_c + 40 * Subsystem_Fy_B.i_h] =
        static_cast<real_T>(Subsystem_Fy_B.iv[Subsystem_Fy_B.m]) -
        Subsystem_Fy_B.dv1[Subsystem_Fy_B.m];
    }
  }

  for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 40; Subsystem_Fy_B.i_h++) {
    for (Subsystem_Fy_B.b_m_c = 0; Subsystem_Fy_B.b_m_c < 40;
         Subsystem_Fy_B.b_m_c++) {
      Subsystem_Fy_B.bkj = 0.0;
      for (Subsystem_Fy_B.m = 0; Subsystem_Fy_B.m < 40; Subsystem_Fy_B.m++) {
        Subsystem_Fy_B.bkj += Subsystem_Fy_B.dv2[40 * Subsystem_Fy_B.m +
          Subsystem_Fy_B.i_h] * Subsystem_Fy_B.dv3[40 * Subsystem_Fy_B.b_m_c +
          Subsystem_Fy_B.m];
      }

      Subsystem_Fy_B.m = 40 * Subsystem_Fy_B.b_m_c + Subsystem_Fy_B.i_h;
      Subsystem_Fy_B.dv1[Subsystem_Fy_B.m] = Subsystem_Fy_B.c_y[Subsystem_Fy_B.m]
        + Subsystem_Fy_B.bkj;
      Subsystem_Fy_B.G_k[Subsystem_Fy_B.b_m_c + 40 * Subsystem_Fy_B.i_h] = 0.0;
    }
  }

  for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 40; Subsystem_Fy_B.i_h++) {
    for (Subsystem_Fy_B.b_m_c = 0; Subsystem_Fy_B.b_m_c < 40;
         Subsystem_Fy_B.b_m_c++) {
      for (Subsystem_Fy_B.m = 0; Subsystem_Fy_B.m < 40; Subsystem_Fy_B.m++) {
        Subsystem_Fy_B.coffset = 40 * Subsystem_Fy_B.i_h + Subsystem_Fy_B.m;
        Subsystem_Fy_B.G_k[Subsystem_Fy_B.coffset] += Subsystem_Fy_B.G[40 *
          Subsystem_Fy_B.b_m_c + Subsystem_Fy_B.m] * Subsystem_Fy_B.dv1[40 *
          Subsystem_Fy_B.i_h + Subsystem_Fy_B.b_m_c];
      }

      errorCov_[Subsystem_Fy_B.i_h + 40 * Subsystem_Fy_B.b_m_c] = 0.0;
    }
  }

  for (Subsystem_Fy_B.i_h = 0; Subsystem_Fy_B.i_h < 40; Subsystem_Fy_B.i_h++) {
    for (Subsystem_Fy_B.b_m_c = 0; Subsystem_Fy_B.b_m_c < 40;
         Subsystem_Fy_B.b_m_c++) {
      for (Subsystem_Fy_B.m = 0; Subsystem_Fy_B.m < 40; Subsystem_Fy_B.m++) {
        Subsystem_Fy_B.coffset = 40 * Subsystem_Fy_B.b_m_c + Subsystem_Fy_B.i_h;
        errorCov_[Subsystem_Fy_B.coffset] += Subsystem_Fy_B.G_k[40 *
          Subsystem_Fy_B.m + Subsystem_Fy_B.i_h] * Subsystem_Fy_B.G[40 *
          Subsystem_Fy_B.m + Subsystem_Fy_B.b_m_c];
      }
    }
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_F_correct_measurement(const real_T
  errorCov_c[1600], const real_T measurementCov_[169], const real_T mm1[3],
  const real_T mm2[3], boolean_T is_stance, const real_T get_NominalState[42],
  real_T x_post[42], real_T errorCov_[1600])
{
  static const int8_T c[7] = { 1, 2, 3, 7, 8, 9, 10 };

  static const int8_T d[10] = { 1, 2, 3, 7, 8, 9, 10, 11, 12, 13 };

  Subsystem__estimate_measurement(get_NominalState, Subsystem_Fy_B.yh,
    Subsystem_Fy_B.H);
  Subsystem_Fy_B.ym[9] = 0.0;
  Subsystem_Fy_B.ym[0] = mm1[0];
  Subsystem_Fy_B.ym[3] = mm2[0];
  Subsystem_Fy_B.ym[6] = 0.0;
  Subsystem_Fy_B.ym[10] = 0.0;
  Subsystem_Fy_B.ym[1] = mm1[1];
  Subsystem_Fy_B.ym[4] = mm2[1];
  Subsystem_Fy_B.ym[7] = 0.0;
  Subsystem_Fy_B.ym[11] = 0.0;
  Subsystem_Fy_B.ym[2] = mm1[2];
  Subsystem_Fy_B.ym[5] = mm2[2];
  Subsystem_Fy_B.ym[8] = 0.0;
  Subsystem_Fy_B.ym[12] = 0.0;
  if (!Subsystem_Fy_DW.rowsf1_not_empty) {
    Subsystem_Fy_DW.rowsf1.size[0] = 1;
    Subsystem_Fy_DW.rowsf1.size[1] = 13;
    for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 13; Subsystem_Fy_B.ch++) {
      Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = static_cast<real_T>
        (Subsystem_Fy_B.ch) + 1.0;
    }

    Subsystem_Fy_DW.rowsf1_not_empty = true;
  }

  Subsystem_Fy_B.mm1_abs.re = (mm1[0] * mm1[0] + mm1[1] * mm1[1]) + mm1[2] *
    mm1[2];
  Subsystem_Fy_B.mm1_abs.im = 0.0;
  Subsystem_Fy_sqrt(&Subsystem_Fy_B.mm1_abs);
  Subsystem_Fy_B.mm2_abs.re = (mm2[0] * mm2[0] + mm2[1] * mm2[1]) + mm2[2] *
    mm2[2];
  Subsystem_Fy_B.mm2_abs.im = 0.0;
  Subsystem_Fy_sqrt(&Subsystem_Fy_B.mm2_abs);
  Subsystem_Fy_B.ch = 0;
  Subsystem_Fy_B.d = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.mm1_abs.re - 1.0,
    Subsystem_Fy_B.mm1_abs.im);
  if (Subsystem_Fy_B.d > 0.2) {
    Subsystem_Fy_B.d1 = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.mm2_abs.re -
      1.0, Subsystem_Fy_B.mm2_abs.im);
    if ((Subsystem_Fy_B.d1 > 0.2) && (!is_stance)) {
      Subsystem_Fy_DW.rowsf1.size[0] = 1;
      Subsystem_Fy_DW.rowsf1.size[1] = 4;
      Subsystem_Fy_DW.rowsf1.data[0] = 7.0;
      Subsystem_Fy_DW.rowsf1.data[1] = 8.0;
      Subsystem_Fy_DW.rowsf1.data[2] = 9.0;
      Subsystem_Fy_DW.rowsf1.data[3] = 10.0;
      Subsystem_Fy_DW.rowsf1_not_empty = true;
      Subsystem_Fy_B.ch = 1;
    }

    if ((Subsystem_Fy_B.d1 > 0.2) && is_stance) {
      Subsystem_Fy_DW.rowsf1.size[0] = 1;
      Subsystem_Fy_DW.rowsf1.size[1] = 7;
      for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 7; Subsystem_Fy_B.ch++) {
        Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = static_cast<real_T>
          (Subsystem_Fy_B.ch) + 7.0;
      }

      Subsystem_Fy_DW.rowsf1_not_empty = true;
      Subsystem_Fy_B.ch = 2;
    }

    if ((Subsystem_Fy_B.d1 < 0.2) && (!is_stance)) {
      Subsystem_Fy_DW.rowsf1.size[0] = 1;
      Subsystem_Fy_DW.rowsf1.size[1] = 7;
      for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 7; Subsystem_Fy_B.ch++) {
        Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = static_cast<real_T>
          (Subsystem_Fy_B.ch) + 4.0;
      }

      Subsystem_Fy_DW.rowsf1_not_empty = true;
      Subsystem_Fy_B.ch = 3;
    }

    if ((Subsystem_Fy_B.d1 < 0.2) && is_stance) {
      Subsystem_Fy_DW.rowsf1.size[0] = 1;
      Subsystem_Fy_DW.rowsf1.size[1] = 10;
      for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 10; Subsystem_Fy_B.ch++) {
        Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = static_cast<real_T>
          (Subsystem_Fy_B.ch) + 4.0;
      }

      Subsystem_Fy_DW.rowsf1_not_empty = true;
      Subsystem_Fy_B.ch = 4;
    }
  }

  if (Subsystem_Fy_B.d < 0.2) {
    Subsystem_Fy_B.d = Subsystem_Fy_rt_hypotd_snf(Subsystem_Fy_B.mm2_abs.re -
      1.0, Subsystem_Fy_B.mm2_abs.im);
    if ((Subsystem_Fy_B.d > 0.2) && (!is_stance)) {
      Subsystem_Fy_DW.rowsf1.size[0] = 1;
      Subsystem_Fy_DW.rowsf1.size[1] = 7;
      for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 7; Subsystem_Fy_B.ch++) {
        Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = c[Subsystem_Fy_B.ch];
      }

      Subsystem_Fy_DW.rowsf1_not_empty = true;
      Subsystem_Fy_B.ch = 5;
    }

    if ((Subsystem_Fy_B.d > 0.2) && is_stance) {
      Subsystem_Fy_DW.rowsf1.size[0] = 1;
      Subsystem_Fy_DW.rowsf1.size[1] = 10;
      for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 10; Subsystem_Fy_B.ch++) {
        Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = d[Subsystem_Fy_B.ch];
      }

      Subsystem_Fy_DW.rowsf1_not_empty = true;
      Subsystem_Fy_B.ch = 6;
    }

    if ((Subsystem_Fy_B.d < 0.2) && (!is_stance)) {
      Subsystem_Fy_DW.rowsf1.size[0] = 1;
      Subsystem_Fy_DW.rowsf1.size[1] = 10;
      for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 10; Subsystem_Fy_B.ch++) {
        Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = static_cast<real_T>
          (Subsystem_Fy_B.ch) + 1.0;
      }

      Subsystem_Fy_DW.rowsf1_not_empty = true;
      Subsystem_Fy_B.ch = 7;
    }
  }

  if (Subsystem_Fy_B.ch == 0) {
    Subsystem_Fy_DW.rowsf1.size[0] = 1;
    Subsystem_Fy_DW.rowsf1.size[1] = 13;
    for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 13; Subsystem_Fy_B.ch++) {
      Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch] = static_cast<real_T>
        (Subsystem_Fy_B.ch) + 1.0;
    }

    Subsystem_Fy_DW.rowsf1_not_empty = true;
    Subsystem_Fy_B.measurementCov_size[0] = 13;
    Subsystem_Fy_B.measurementCov_size[1] = 13;
    std::memcpy(&Subsystem_Fy_B.measurementCov_data[0], &measurementCov_[0],
                169U * sizeof(real_T));
    Subsystem_Fy_B.b_H_size[0] = 13;
    Subsystem_Fy_B.b_H_size[1] = 40;
    std::memcpy(&Subsystem_Fy_B.b_H_data[0], &Subsystem_Fy_B.H[0], 520U * sizeof
                (real_T));
    for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 13; Subsystem_Fy_B.ch++) {
      Subsystem_Fy_B.res_data[Subsystem_Fy_B.ch] =
        Subsystem_Fy_B.ym[Subsystem_Fy_B.ch] -
        Subsystem_Fy_B.yh[Subsystem_Fy_B.ch];
    }
  } else {
    Subsystem_Fy_B.measurementCov_tmp_size_idx_0 = Subsystem_Fy_DW.rowsf1.size[1];
    Subsystem_Fy_B.loop_ub = Subsystem_Fy_DW.rowsf1.size[1];
    for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < Subsystem_Fy_B.loop_ub;
         Subsystem_Fy_B.ch++) {
      Subsystem_Fy_B.measurementCov_tmp_data[Subsystem_Fy_B.ch] =
        static_cast<int32_T>(Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch]);
    }

    Subsystem_Fy_B.measurementCov_size[0] = Subsystem_Fy_DW.rowsf1.size[1];
    Subsystem_Fy_B.measurementCov_size[1] = Subsystem_Fy_DW.rowsf1.size[1];
    for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch <
         Subsystem_Fy_B.measurementCov_tmp_size_idx_0; Subsystem_Fy_B.ch++) {
      for (Subsystem_Fy_B.i_m = 0; Subsystem_Fy_B.i_m <
           Subsystem_Fy_B.measurementCov_tmp_size_idx_0; Subsystem_Fy_B.i_m++) {
        Subsystem_Fy_B.measurementCov_data[Subsystem_Fy_B.i_m +
          Subsystem_Fy_B.measurementCov_tmp_size_idx_0 * Subsystem_Fy_B.ch] =
          measurementCov_
          [((Subsystem_Fy_B.measurementCov_tmp_data[Subsystem_Fy_B.ch] - 1) * 13
            + Subsystem_Fy_B.measurementCov_tmp_data[Subsystem_Fy_B.i_m]) - 1];
      }
    }

    Subsystem_Fy_B.b_H_size[0] = Subsystem_Fy_DW.rowsf1.size[1];
    Subsystem_Fy_B.b_H_size[1] = 40;
    Subsystem_Fy_B.loop_ub = Subsystem_Fy_DW.rowsf1.size[1];
    for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < 40; Subsystem_Fy_B.ch++) {
      for (Subsystem_Fy_B.i_m = 0; Subsystem_Fy_B.i_m < Subsystem_Fy_B.loop_ub;
           Subsystem_Fy_B.i_m++) {
        Subsystem_Fy_B.b_H_data[Subsystem_Fy_B.i_m + Subsystem_Fy_B.b_H_size[0] *
          Subsystem_Fy_B.ch] = Subsystem_Fy_B.H[(13 * Subsystem_Fy_B.ch +
          static_cast<int32_T>(Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.i_m]))
          - 1];
      }
    }

    Subsystem_Fy_B.loop_ub = Subsystem_Fy_DW.rowsf1.size[1];
    for (Subsystem_Fy_B.ch = 0; Subsystem_Fy_B.ch < Subsystem_Fy_B.loop_ub;
         Subsystem_Fy_B.ch++) {
      Subsystem_Fy_B.measurementCov_tmp_size_idx_0 = static_cast<int32_T>
        (Subsystem_Fy_DW.rowsf1.data[Subsystem_Fy_B.ch]) - 1;
      Subsystem_Fy_B.res_data[Subsystem_Fy_B.ch] =
        Subsystem_Fy_B.ym[Subsystem_Fy_B.measurementCov_tmp_size_idx_0] -
        Subsystem_Fy_B.yh[Subsystem_Fy_B.measurementCov_tmp_size_idx_0];
    }
  }

  std::memcpy(&errorCov_[0], &errorCov_c[0], 1600U * sizeof(real_T));
  Subsystem_Fy_correct(get_NominalState, errorCov_, Subsystem_Fy_B.res_data,
                       Subsystem_Fy_B.measurementCov_data,
                       Subsystem_Fy_B.measurementCov_size,
                       Subsystem_Fy_B.b_H_data, Subsystem_Fy_B.b_H_size, x_post);
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T tmp;
  real_T tmp_0;
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = (rtNaN);
    } else {
      y = std::pow(u0, u1);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_predict_eksf(const real_T
  processNoise[900], const real_T U[12], const real_T x_nominal_prev[42], const
  real_T errorstate_prev[1600], real_T NominalState_predict[42], real_T
  ErrorState_predict[1600])
{
  real_T t1083_tmp;
  real_T t1083_tmp_0;
  real_T t1083_tmp_1;
  int32_T b_tmp;
  int32_T i;
  int32_T i_0;
  int32_T i_1;
  Subsystem_Fy_B.t3 = x_nominal_prev[6] * x_nominal_prev[6];
  Subsystem_Fy_B.t4 = x_nominal_prev[7] * x_nominal_prev[7];
  Subsystem_Fy_B.t5 = x_nominal_prev[8] * x_nominal_prev[8];
  Subsystem_Fy_B.t6 = x_nominal_prev[9] * x_nominal_prev[9];
  Subsystem_Fy_B.t10 = x_nominal_prev[25] * x_nominal_prev[25];
  Subsystem_Fy_B.t11 = x_nominal_prev[26] * x_nominal_prev[26];
  Subsystem_Fy_B.t12 = x_nominal_prev[27] * x_nominal_prev[27];
  Subsystem_Fy_B.t13 = x_nominal_prev[28] * x_nominal_prev[28];
  Subsystem_Fy_B.t195 = x_nominal_prev[6] * x_nominal_prev[7];
  Subsystem_Fy_B.t14 = Subsystem_Fy_B.t195 * 2.0;
  Subsystem_Fy_B.t194 = x_nominal_prev[6] * x_nominal_prev[8];
  Subsystem_Fy_B.t15 = Subsystem_Fy_B.t194 * 2.0;
  Subsystem_Fy_B.t254 = x_nominal_prev[6] * x_nominal_prev[9];
  Subsystem_Fy_B.t16 = Subsystem_Fy_B.t254 * 2.0;
  Subsystem_Fy_B.t193 = x_nominal_prev[7] * x_nominal_prev[8];
  Subsystem_Fy_B.t17 = Subsystem_Fy_B.t193 * 2.0;
  Subsystem_Fy_B.t225 = x_nominal_prev[7] * x_nominal_prev[9];
  Subsystem_Fy_B.t18 = Subsystem_Fy_B.t225 * 2.0;
  Subsystem_Fy_B.t218 = x_nominal_prev[8] * x_nominal_prev[9];
  Subsystem_Fy_B.t19 = Subsystem_Fy_B.t218 * 2.0;
  Subsystem_Fy_B.t240 = x_nominal_prev[25] * x_nominal_prev[26];
  Subsystem_Fy_B.t20 = Subsystem_Fy_B.t240 * 2.0;
  Subsystem_Fy_B.t239 = x_nominal_prev[25] * x_nominal_prev[27];
  Subsystem_Fy_B.t21 = Subsystem_Fy_B.t239 * 2.0;
  Subsystem_Fy_B.t250 = x_nominal_prev[25] * x_nominal_prev[28];
  Subsystem_Fy_B.t22 = Subsystem_Fy_B.t250 * 2.0;
  Subsystem_Fy_B.t249 = x_nominal_prev[26] * x_nominal_prev[27];
  Subsystem_Fy_B.t23 = Subsystem_Fy_B.t249 * 2.0;
  Subsystem_Fy_B.t209 = x_nominal_prev[26] * x_nominal_prev[28];
  Subsystem_Fy_B.t24 = Subsystem_Fy_B.t209 * 2.0;
  Subsystem_Fy_B.t210 = x_nominal_prev[27] * x_nominal_prev[28];
  Subsystem_Fy_B.t25 = Subsystem_Fy_B.t210 * 2.0;
  Subsystem_Fy_B.t43 = Subsystem_Fy_B.t4 * 2.0;
  Subsystem_Fy_B.t44 = Subsystem_Fy_B.t5 * 2.0;
  Subsystem_Fy_B.t45 = Subsystem_Fy_B.t6 * 2.0;
  Subsystem_Fy_B.t49 = Subsystem_Fy_B.t11 * 2.0;
  Subsystem_Fy_B.t50 = Subsystem_Fy_B.t12 * 2.0;
  Subsystem_Fy_B.t51 = Subsystem_Fy_B.t13 * 2.0;
  Subsystem_Fy_B.t67 = U[0] + -x_nominal_prev[10];
  Subsystem_Fy_B.t68 = U[1] + -x_nominal_prev[11];
  Subsystem_Fy_B.t69 = U[2] + -x_nominal_prev[12];
  Subsystem_Fy_B.t70 = U[3] + -x_nominal_prev[13];
  Subsystem_Fy_B.t71 = U[4] + -x_nominal_prev[14];
  Subsystem_Fy_B.t72 = U[5] + -x_nominal_prev[15];
  Subsystem_Fy_B.t73 = U[6] + -x_nominal_prev[29];
  Subsystem_Fy_B.t74 = U[7] + -x_nominal_prev[30];
  Subsystem_Fy_B.t75 = U[8] + -x_nominal_prev[31];
  Subsystem_Fy_B.t76 = U[9] + -x_nominal_prev[32];
  Subsystem_Fy_B.t77 = U[10] + -x_nominal_prev[33];
  Subsystem_Fy_B.t78 = U[11] + -x_nominal_prev[34];
  Subsystem_Fy_B.t133 = ((Subsystem_Fy_B.t3 + Subsystem_Fy_B.t4) +
    Subsystem_Fy_B.t5) + Subsystem_Fy_B.t6;
  Subsystem_Fy_B.t137 = ((Subsystem_Fy_B.t10 + Subsystem_Fy_B.t11) +
    Subsystem_Fy_B.t12) + Subsystem_Fy_B.t13;
  Subsystem_Fy_B.t121 = U[9] / 2.0 + -(x_nominal_prev[32] / 2.0);
  Subsystem_Fy_B.t122 = U[10] / 2.0 + -(x_nominal_prev[33] / 2.0);
  Subsystem_Fy_B.t123 = U[11] / 2.0 + -(x_nominal_prev[34] / 2.0);
  Subsystem_Fy_B.t130 = U[3] / 2.0 + -(x_nominal_prev[13] / 2.0);
  Subsystem_Fy_B.t131 = U[4] / 2.0 + -(x_nominal_prev[14] / 2.0);
  Subsystem_Fy_B.t132 = U[5] / 2.0 + -(x_nominal_prev[15] / 2.0);
  Subsystem_Fy_B.t165 = 1.0 / Subsystem_Fy_B.t137;
  Subsystem_Fy_B.t166 = (Subsystem_Fy_B.t16 + Subsystem_Fy_B.t17) *
    Subsystem_Fy_B.t67;
  Subsystem_Fy_B.t167 = (Subsystem_Fy_B.t14 + Subsystem_Fy_B.t19) *
    Subsystem_Fy_B.t68;
  Subsystem_Fy_B.t168 = (Subsystem_Fy_B.t15 + Subsystem_Fy_B.t18) *
    Subsystem_Fy_B.t69;
  Subsystem_Fy_B.t169 = 1.0 / Subsystem_Fy_B.t133;
  Subsystem_Fy_B.t170 = (Subsystem_Fy_B.t22 + Subsystem_Fy_B.t23) *
    Subsystem_Fy_B.t73;
  Subsystem_Fy_B.t171 = (Subsystem_Fy_B.t20 + Subsystem_Fy_B.t25) *
    Subsystem_Fy_B.t74;
  Subsystem_Fy_B.t172 = (Subsystem_Fy_B.t21 + Subsystem_Fy_B.t24) *
    Subsystem_Fy_B.t75;
  Subsystem_Fy_B.t185 = 1.0 / std::sqrt(Subsystem_Fy_B.t137);
  Subsystem_Fy_B.t186 = 1.0 / std::sqrt(Subsystem_Fy_B.t133);
  Subsystem_Fy_B.t141 = 0.0025 * x_nominal_prev[6] * Subsystem_Fy_B.t130;
  Subsystem_Fy_B.t133 = 0.0025 * x_nominal_prev[7] * Subsystem_Fy_B.t130;
  Subsystem_Fy_B.t143 = 0.0025 * x_nominal_prev[6] * Subsystem_Fy_B.t131;
  Subsystem_Fy_B.t144 = 0.0025 * x_nominal_prev[8] * Subsystem_Fy_B.t130;
  Subsystem_Fy_B.t145 = 0.0025 * x_nominal_prev[7] * Subsystem_Fy_B.t131;
  Subsystem_Fy_B.t146 = 0.0025 * x_nominal_prev[9] * Subsystem_Fy_B.t130;
  Subsystem_Fy_B.t147 = 0.0025 * x_nominal_prev[6] * Subsystem_Fy_B.t132;
  Subsystem_Fy_B.t137 = 0.0025 * x_nominal_prev[8] * Subsystem_Fy_B.t131;
  Subsystem_Fy_B.t149 = 0.0025 * x_nominal_prev[7] * Subsystem_Fy_B.t132;
  Subsystem_Fy_B.t150 = 0.0025 * x_nominal_prev[9] * Subsystem_Fy_B.t131;
  Subsystem_Fy_B.t151 = 0.0025 * x_nominal_prev[8] * Subsystem_Fy_B.t132;
  Subsystem_Fy_B.t152 = 0.0025 * x_nominal_prev[9] * Subsystem_Fy_B.t132;
  Subsystem_Fy_B.t153 = 0.0025 * x_nominal_prev[25] * Subsystem_Fy_B.t121;
  Subsystem_Fy_B.t154 = 0.0025 * x_nominal_prev[26] * Subsystem_Fy_B.t121;
  Subsystem_Fy_B.t155 = 0.0025 * x_nominal_prev[25] * Subsystem_Fy_B.t122;
  Subsystem_Fy_B.t156 = 0.0025 * x_nominal_prev[27] * Subsystem_Fy_B.t121;
  Subsystem_Fy_B.t157 = 0.0025 * x_nominal_prev[26] * Subsystem_Fy_B.t122;
  Subsystem_Fy_B.t158 = 0.0025 * x_nominal_prev[28] * Subsystem_Fy_B.t121;
  Subsystem_Fy_B.t159 = 0.0025 * x_nominal_prev[25] * Subsystem_Fy_B.t123;
  Subsystem_Fy_B.t160 = 0.0025 * x_nominal_prev[27] * Subsystem_Fy_B.t122;
  Subsystem_Fy_B.t161 = 0.0025 * x_nominal_prev[26] * Subsystem_Fy_B.t123;
  Subsystem_Fy_B.t162 = 0.0025 * x_nominal_prev[28] * Subsystem_Fy_B.t122;
  Subsystem_Fy_B.t163 = 0.0025 * x_nominal_prev[27] * Subsystem_Fy_B.t123;
  Subsystem_Fy_B.t164 = 0.0025 * x_nominal_prev[28] * Subsystem_Fy_B.t123;
  Subsystem_Fy_B.t123 = (Subsystem_Fy_B.t15 + -Subsystem_Fy_B.t18) *
    Subsystem_Fy_B.t67;
  Subsystem_Fy_B.t122 = (Subsystem_Fy_B.t16 + -Subsystem_Fy_B.t17) *
    Subsystem_Fy_B.t68;
  Subsystem_Fy_B.t121 = (Subsystem_Fy_B.t14 + -Subsystem_Fy_B.t19) *
    Subsystem_Fy_B.t69;
  Subsystem_Fy_B.t132 = (Subsystem_Fy_B.t21 + -Subsystem_Fy_B.t24) *
    Subsystem_Fy_B.t73;
  Subsystem_Fy_B.t131 = (Subsystem_Fy_B.t22 + -Subsystem_Fy_B.t23) *
    Subsystem_Fy_B.t74;
  Subsystem_Fy_B.t130 = (Subsystem_Fy_B.t20 + -Subsystem_Fy_B.t25) *
    Subsystem_Fy_B.t75;
  Subsystem_Fy_B.t187 = Subsystem_Fy_B.t195 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t188 = Subsystem_Fy_B.t194 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t189 = Subsystem_Fy_B.t254 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t190 = Subsystem_Fy_B.t193 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t191 = Subsystem_Fy_B.t225 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t192 = Subsystem_Fy_B.t218 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t193 = Subsystem_Fy_B.t240 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t194 = Subsystem_Fy_B.t239 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t195 = Subsystem_Fy_B.t250 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t196 = Subsystem_Fy_B.t249 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t197 = Subsystem_Fy_B.t209 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t198 = Subsystem_Fy_B.t210 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t14 *= Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t15 *= Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t16 *= Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t205 = Subsystem_Fy_B.t10 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t206 = Subsystem_Fy_B.t11 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t207 = Subsystem_Fy_B.t12 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t208 = Subsystem_Fy_B.t13 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t209 = Subsystem_Fy_B.t20 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t210 = Subsystem_Fy_B.t21 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t10 = Subsystem_Fy_B.t22 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t22 = ((Subsystem_Fy_B.t44 + Subsystem_Fy_B.t45) - 1.0) *
    Subsystem_Fy_B.t67;
  Subsystem_Fy_B.t21 = ((Subsystem_Fy_B.t43 + Subsystem_Fy_B.t45) - 1.0) *
    Subsystem_Fy_B.t68;
  Subsystem_Fy_B.t20 = ((Subsystem_Fy_B.t43 + Subsystem_Fy_B.t44) - 1.0) *
    Subsystem_Fy_B.t69;
  Subsystem_Fy_B.t218 = x_nominal_prev[6] * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t219 = x_nominal_prev[7] * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t220 = x_nominal_prev[8] * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t221 = x_nominal_prev[9] * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t13 = ((Subsystem_Fy_B.t50 + Subsystem_Fy_B.t51) - 1.0) *
    Subsystem_Fy_B.t73;
  Subsystem_Fy_B.t12 = ((Subsystem_Fy_B.t49 + Subsystem_Fy_B.t51) - 1.0) *
    Subsystem_Fy_B.t74;
  Subsystem_Fy_B.t11 = ((Subsystem_Fy_B.t49 + Subsystem_Fy_B.t50) - 1.0) *
    Subsystem_Fy_B.t75;
  Subsystem_Fy_B.t225 = x_nominal_prev[25] * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t226 = x_nominal_prev[26] * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t227 = x_nominal_prev[27] * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t228 = x_nominal_prev[28] * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t229 = Subsystem_Fy_B.t3 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t230 = Subsystem_Fy_B.t4 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t231 = Subsystem_Fy_B.t5 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t232 = Subsystem_Fy_B.t6 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t239 = Subsystem_Fy_B.t49 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t240 = Subsystem_Fy_B.t50 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t241 = Subsystem_Fy_B.t51 * Subsystem_Fy_B.t165;
  Subsystem_Fy_B.t3 = Subsystem_Fy_B.t43 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t249 = Subsystem_Fy_B.t44 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t250 = Subsystem_Fy_B.t45 * Subsystem_Fy_B.t169;
  Subsystem_Fy_B.t254 = -x_nominal_prev[6] * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t255 = -x_nominal_prev[25] * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t45 = (x_nominal_prev[7] * Subsystem_Fy_B.t70 / 2.0 +
                        x_nominal_prev[8] * Subsystem_Fy_B.t71 / 2.0) +
    x_nominal_prev[9] * Subsystem_Fy_B.t72 / 2.0;
  Subsystem_Fy_B.t44 = (x_nominal_prev[26] * Subsystem_Fy_B.t76 / 2.0 +
                        x_nominal_prev[27] * Subsystem_Fy_B.t77 / 2.0) +
    x_nominal_prev[28] * Subsystem_Fy_B.t78 / 2.0;
  Subsystem_Fy_B.t269 = Subsystem_Fy_B.t153 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t270 = Subsystem_Fy_B.t154 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t271 = Subsystem_Fy_B.t155 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t273 = Subsystem_Fy_B.t157 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t274 = Subsystem_Fy_B.t158 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t275 = Subsystem_Fy_B.t159 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t276 = Subsystem_Fy_B.t160 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t279 = Subsystem_Fy_B.t163 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t280 = Subsystem_Fy_B.t164 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t281 = Subsystem_Fy_B.t141 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t282 = Subsystem_Fy_B.t133 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t283 = Subsystem_Fy_B.t143 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t285 = Subsystem_Fy_B.t145 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t286 = Subsystem_Fy_B.t146 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t287 = Subsystem_Fy_B.t147 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t288 = Subsystem_Fy_B.t137 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t291 = Subsystem_Fy_B.t151 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t292 = Subsystem_Fy_B.t152 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t43 = (x_nominal_prev[7] * Subsystem_Fy_B.t71 / 2.0 +
                        x_nominal_prev[6] * Subsystem_Fy_B.t72 / 2.0) +
    -(x_nominal_prev[8] * Subsystem_Fy_B.t70 / 2.0);
  Subsystem_Fy_B.t51 = (x_nominal_prev[6] * Subsystem_Fy_B.t71 / 2.0 +
                        x_nominal_prev[9] * Subsystem_Fy_B.t70 / 2.0) +
    -(x_nominal_prev[7] * Subsystem_Fy_B.t72 / 2.0);
  Subsystem_Fy_B.t50 = (x_nominal_prev[6] * Subsystem_Fy_B.t70 / 2.0 +
                        x_nominal_prev[8] * Subsystem_Fy_B.t72 / 2.0) +
    -(x_nominal_prev[9] * Subsystem_Fy_B.t71 / 2.0);
  Subsystem_Fy_B.t49 = (x_nominal_prev[26] * Subsystem_Fy_B.t77 / 2.0 +
                        x_nominal_prev[25] * Subsystem_Fy_B.t78 / 2.0) +
    -(x_nominal_prev[27] * Subsystem_Fy_B.t76 / 2.0);
  Subsystem_Fy_B.t6 = (x_nominal_prev[25] * Subsystem_Fy_B.t77 / 2.0 +
                       x_nominal_prev[28] * Subsystem_Fy_B.t76 / 2.0) +
    -(x_nominal_prev[26] * Subsystem_Fy_B.t78 / 2.0);
  Subsystem_Fy_B.t5 = (x_nominal_prev[25] * Subsystem_Fy_B.t76 / 2.0 +
                       x_nominal_prev[27] * Subsystem_Fy_B.t78 / 2.0) +
    -(x_nominal_prev[28] * Subsystem_Fy_B.t77 / 2.0);
  Subsystem_Fy_B.t356 = Subsystem_Fy_B.t25 * Subsystem_Fy_B.t165 +
    Subsystem_Fy_B.t209;
  Subsystem_Fy_B.t357 = Subsystem_Fy_B.t24 * Subsystem_Fy_B.t165 +
    Subsystem_Fy_B.t210;
  Subsystem_Fy_B.t358 = Subsystem_Fy_B.t23 * Subsystem_Fy_B.t165 +
    Subsystem_Fy_B.t10;
  Subsystem_Fy_B.t362 = Subsystem_Fy_B.t19 * Subsystem_Fy_B.t169 +
    Subsystem_Fy_B.t14;
  Subsystem_Fy_B.t363 = Subsystem_Fy_B.t18 * Subsystem_Fy_B.t169 +
    Subsystem_Fy_B.t15;
  Subsystem_Fy_B.t364 = Subsystem_Fy_B.t17 * Subsystem_Fy_B.t169 +
    Subsystem_Fy_B.t16;
  Subsystem_Fy_B.t24 = Subsystem_Fy_B.t192 * 2.0 - Subsystem_Fy_B.t14;
  Subsystem_Fy_B.t17 = Subsystem_Fy_B.t24 * -0.0025;
  Subsystem_Fy_B.t25 = Subsystem_Fy_B.t191 * 2.0 - Subsystem_Fy_B.t15;
  Subsystem_Fy_B.t169 = Subsystem_Fy_B.t25 * -0.0025;
  Subsystem_Fy_B.t16 = Subsystem_Fy_B.t190 * 2.0 - Subsystem_Fy_B.t16;
  Subsystem_Fy_B.t18 = Subsystem_Fy_B.t16 * -0.0025;
  Subsystem_Fy_B.t15 = Subsystem_Fy_B.t198 * 2.0 - Subsystem_Fy_B.t209;
  Subsystem_Fy_B.t19 = Subsystem_Fy_B.t15 * -0.0025;
  Subsystem_Fy_B.t14 = Subsystem_Fy_B.t197 * 2.0 - Subsystem_Fy_B.t210;
  Subsystem_Fy_B.t23 = Subsystem_Fy_B.t14 * -0.0025;
  Subsystem_Fy_B.t10 = Subsystem_Fy_B.t196 * 2.0 - Subsystem_Fy_B.t10;
  Subsystem_Fy_B.t165 = Subsystem_Fy_B.t10 * -0.0025;
  Subsystem_Fy_B.t133 = ((-x_nominal_prev[6] + Subsystem_Fy_B.t133) +
    Subsystem_Fy_B.t137) + Subsystem_Fy_B.t152;
  Subsystem_Fy_B.t137 = ((-x_nominal_prev[25] + Subsystem_Fy_B.t154) +
    Subsystem_Fy_B.t160) + Subsystem_Fy_B.t164;
  Subsystem_Fy_B.t471 = Subsystem_Fy_B.t14 * Subsystem_Fy_B.t73;
  Subsystem_Fy_B.t472 = Subsystem_Fy_B.t10 * Subsystem_Fy_B.t74;
  Subsystem_Fy_B.t474 = Subsystem_Fy_B.t15 * Subsystem_Fy_B.t75;
  Subsystem_Fy_B.t477 = Subsystem_Fy_B.t25 * Subsystem_Fy_B.t67;
  Subsystem_Fy_B.t478 = Subsystem_Fy_B.t16 * Subsystem_Fy_B.t68;
  Subsystem_Fy_B.t480 = Subsystem_Fy_B.t24 * Subsystem_Fy_B.t69;
  Subsystem_Fy_B.t154 = (Subsystem_Fy_B.t193 + Subsystem_Fy_B.t198) * 6.25E-6;
  Subsystem_Fy_B.t152 = (Subsystem_Fy_B.t194 + Subsystem_Fy_B.t197) * 6.25E-6;
  Subsystem_Fy_B.t160 = (Subsystem_Fy_B.t195 + Subsystem_Fy_B.t196) * 6.25E-6;
  Subsystem_Fy_B.t164 = 0.0025 * Subsystem_Fy_B.t362;
  Subsystem_Fy_B.t24 = 0.0025 * Subsystem_Fy_B.t363;
  Subsystem_Fy_B.t25 = 0.0025 * Subsystem_Fy_B.t364;
  Subsystem_Fy_B.t16 = (Subsystem_Fy_B.t187 + Subsystem_Fy_B.t192) * 6.25E-6;
  Subsystem_Fy_B.t15 = (Subsystem_Fy_B.t188 + Subsystem_Fy_B.t191) * 6.25E-6;
  Subsystem_Fy_B.t10 = (Subsystem_Fy_B.t189 + Subsystem_Fy_B.t190) * 6.25E-6;
  Subsystem_Fy_B.t14 = 0.0025 * Subsystem_Fy_B.t356;
  Subsystem_Fy_B.t210 = 0.0025 * Subsystem_Fy_B.t357;
  Subsystem_Fy_B.t209 = 0.0025 * Subsystem_Fy_B.t358;
  Subsystem_Fy_B.t4 = ((Subsystem_Fy_B.t3 + Subsystem_Fy_B.t249) - 1.0) * 0.0025;
  Subsystem_Fy_B.t3 = ((Subsystem_Fy_B.t3 + Subsystem_Fy_B.t250) - 1.0) * 0.0025;
  Subsystem_Fy_B.t249 = ((Subsystem_Fy_B.t249 + Subsystem_Fy_B.t250) - 1.0) *
    0.0025;
  Subsystem_Fy_B.t250 = ((Subsystem_Fy_B.t239 + Subsystem_Fy_B.t240) - 1.0) *
    0.0025;
  Subsystem_Fy_B.t239 = ((Subsystem_Fy_B.t239 + Subsystem_Fy_B.t241) - 1.0) *
    0.0025;
  Subsystem_Fy_B.t240 = ((Subsystem_Fy_B.t240 + Subsystem_Fy_B.t241) - 1.0) *
    0.0025;
  Subsystem_Fy_B.t358 *= Subsystem_Fy_B.t73;
  Subsystem_Fy_B.t446 = Subsystem_Fy_B.t74 * Subsystem_Fy_B.t356;
  Subsystem_Fy_B.t447 = Subsystem_Fy_B.t75 * Subsystem_Fy_B.t357;
  Subsystem_Fy_B.t364 *= Subsystem_Fy_B.t67;
  Subsystem_Fy_B.t449 = Subsystem_Fy_B.t68 * Subsystem_Fy_B.t362;
  Subsystem_Fy_B.t450 = Subsystem_Fy_B.t69 * Subsystem_Fy_B.t363;
  Subsystem_Fy_B.t141 = ((x_nominal_prev[7] + Subsystem_Fy_B.t141) +
    Subsystem_Fy_B.t151) + -Subsystem_Fy_B.t150;
  Subsystem_Fy_B.t143 = ((x_nominal_prev[8] + Subsystem_Fy_B.t143) +
    Subsystem_Fy_B.t146) + -Subsystem_Fy_B.t149;
  Subsystem_Fy_B.t151 = ((x_nominal_prev[9] + Subsystem_Fy_B.t145) +
    Subsystem_Fy_B.t147) + -Subsystem_Fy_B.t144;
  Subsystem_Fy_B.t145 = ((Subsystem_Fy_B.t230 + Subsystem_Fy_B.t231) - 0.5) *
    6.25E-6;
  Subsystem_Fy_B.t147 = ((Subsystem_Fy_B.t230 + Subsystem_Fy_B.t232) - 0.5) *
    6.25E-6;
  Subsystem_Fy_B.t146 = ((Subsystem_Fy_B.t231 + Subsystem_Fy_B.t232) - 0.5) *
    6.25E-6;
  Subsystem_Fy_B.t153 = ((x_nominal_prev[26] + Subsystem_Fy_B.t153) +
    Subsystem_Fy_B.t163) + -Subsystem_Fy_B.t162;
  Subsystem_Fy_B.t155 = ((x_nominal_prev[27] + Subsystem_Fy_B.t155) +
    Subsystem_Fy_B.t158) + -Subsystem_Fy_B.t161;
  Subsystem_Fy_B.t157 = ((x_nominal_prev[28] + Subsystem_Fy_B.t157) +
    Subsystem_Fy_B.t159) + -Subsystem_Fy_B.t156;
  Subsystem_Fy_B.t159 = ((Subsystem_Fy_B.t206 + Subsystem_Fy_B.t207) - 0.5) *
    6.25E-6;
  Subsystem_Fy_B.t158 = ((Subsystem_Fy_B.t206 + Subsystem_Fy_B.t208) - 0.5) *
    6.25E-6;
  Subsystem_Fy_B.t163 = ((Subsystem_Fy_B.t207 + Subsystem_Fy_B.t208) - 0.5) *
    6.25E-6;
  Subsystem_Fy_B.t578 = ((Subsystem_Fy_B.t255 + Subsystem_Fy_B.t270) +
    Subsystem_Fy_B.t276) + Subsystem_Fy_B.t280;
  Subsystem_Fy_B.t582 = ((Subsystem_Fy_B.t254 + Subsystem_Fy_B.t282) +
    Subsystem_Fy_B.t288) + Subsystem_Fy_B.t292;
  Subsystem_Fy_B.t362 = -(0.0025 * Subsystem_Fy_B.t44 / 2.0) + x_nominal_prev[25];
  Subsystem_Fy_B.t357 = 0.0025 * Subsystem_Fy_B.t5 / 2.0 + x_nominal_prev[26];
  Subsystem_Fy_B.t356 = 0.0025 * Subsystem_Fy_B.t6 / 2.0 + x_nominal_prev[27];
  Subsystem_Fy_B.t241 = 0.0025 * Subsystem_Fy_B.t49 / 2.0 + x_nominal_prev[28];
  Subsystem_Fy_B.t375 = -(0.0025 * Subsystem_Fy_B.t45 / 2.0) + x_nominal_prev[6];
  Subsystem_Fy_B.t376 = 0.0025 * Subsystem_Fy_B.t50 / 2.0 + x_nominal_prev[7];
  Subsystem_Fy_B.t377 = 0.0025 * Subsystem_Fy_B.t51 / 2.0 + x_nominal_prev[8];
  Subsystem_Fy_B.t363 = 0.0025 * Subsystem_Fy_B.t43 / 2.0 + x_nominal_prev[9];
  Subsystem_Fy_B.t193 = (Subsystem_Fy_B.t193 + -Subsystem_Fy_B.t198) * 6.25E-6;
  Subsystem_Fy_B.t194 = (Subsystem_Fy_B.t194 + -Subsystem_Fy_B.t197) * 6.25E-6;
  Subsystem_Fy_B.t195 = (Subsystem_Fy_B.t195 + -Subsystem_Fy_B.t196) * 6.25E-6;
  Subsystem_Fy_B.t187 = (Subsystem_Fy_B.t187 + -Subsystem_Fy_B.t192) * 6.25E-6;
  Subsystem_Fy_B.t196 = (Subsystem_Fy_B.t188 + -Subsystem_Fy_B.t191) * 6.25E-6;
  Subsystem_Fy_B.t197 = (Subsystem_Fy_B.t189 + -Subsystem_Fy_B.t190) * 6.25E-6;
  Subsystem_Fy_B.t557 = (((Subsystem_Fy_B.t229 + Subsystem_Fy_B.t230) +
    -Subsystem_Fy_B.t231) + -Subsystem_Fy_B.t232) * Subsystem_Fy_B.t67;
  Subsystem_Fy_B.t558 = (((Subsystem_Fy_B.t229 + Subsystem_Fy_B.t231) +
    -Subsystem_Fy_B.t230) + -Subsystem_Fy_B.t232) * Subsystem_Fy_B.t68;
  Subsystem_Fy_B.t231 = (((Subsystem_Fy_B.t229 + Subsystem_Fy_B.t232) +
    -Subsystem_Fy_B.t230) + -Subsystem_Fy_B.t231) * Subsystem_Fy_B.t69;
  Subsystem_Fy_B.t229 = (((Subsystem_Fy_B.t205 + Subsystem_Fy_B.t206) +
    -Subsystem_Fy_B.t207) + -Subsystem_Fy_B.t208) * Subsystem_Fy_B.t73;
  Subsystem_Fy_B.t230 = (((Subsystem_Fy_B.t205 + Subsystem_Fy_B.t207) +
    -Subsystem_Fy_B.t206) + -Subsystem_Fy_B.t208) * Subsystem_Fy_B.t74;
  Subsystem_Fy_B.t232 = (((Subsystem_Fy_B.t205 + Subsystem_Fy_B.t208) +
    -Subsystem_Fy_B.t206) + -Subsystem_Fy_B.t207) * Subsystem_Fy_B.t75;
  Subsystem_Fy_B.t575 = ((Subsystem_Fy_B.t226 + Subsystem_Fy_B.t269) +
    Subsystem_Fy_B.t279) + -Subsystem_Fy_B.t162 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t576 = ((Subsystem_Fy_B.t227 + Subsystem_Fy_B.t271) +
    Subsystem_Fy_B.t274) + -Subsystem_Fy_B.t161 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t577 = ((Subsystem_Fy_B.t228 + Subsystem_Fy_B.t273) +
    Subsystem_Fy_B.t275) + -Subsystem_Fy_B.t156 * Subsystem_Fy_B.t185;
  Subsystem_Fy_B.t579 = ((Subsystem_Fy_B.t219 + Subsystem_Fy_B.t281) +
    Subsystem_Fy_B.t291) + -Subsystem_Fy_B.t150 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t580 = ((Subsystem_Fy_B.t220 + Subsystem_Fy_B.t283) +
    Subsystem_Fy_B.t286) + -Subsystem_Fy_B.t149 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t581 = ((Subsystem_Fy_B.t221 + Subsystem_Fy_B.t285) +
    Subsystem_Fy_B.t287) + -Subsystem_Fy_B.t144 * Subsystem_Fy_B.t186;
  Subsystem_Fy_B.t270 = ((-(Subsystem_Fy_B.t225 / 2.0) + Subsystem_Fy_B.t270 /
    2.0) + Subsystem_Fy_B.t276 / 2.0) + Subsystem_Fy_B.t280 / 2.0;
  Subsystem_Fy_B.t292 = ((-(Subsystem_Fy_B.t218 / 2.0) + Subsystem_Fy_B.t282 /
    2.0) + Subsystem_Fy_B.t288 / 2.0) + Subsystem_Fy_B.t292 / 2.0;
  Subsystem_Fy_B.t208 = Subsystem_Fy_B.t375 * Subsystem_Fy_B.t376 * 2.0;
  Subsystem_Fy_B.t206 = Subsystem_Fy_B.t375 * Subsystem_Fy_B.t377 * 2.0;
  Subsystem_Fy_B.t282 = Subsystem_Fy_B.t375 * Subsystem_Fy_B.t363 * 2.0;
  Subsystem_Fy_B.t288 = Subsystem_Fy_B.t376 * Subsystem_Fy_B.t377 * 2.0;
  Subsystem_Fy_B.t207 = Subsystem_Fy_B.t376 * Subsystem_Fy_B.t363 * 2.0;
  Subsystem_Fy_B.t198 = Subsystem_Fy_B.t377 * Subsystem_Fy_B.t363 * 2.0;
  Subsystem_Fy_B.t188 = Subsystem_Fy_B.t362 * Subsystem_Fy_B.t357 * 2.0;
  Subsystem_Fy_B.t189 = Subsystem_Fy_B.t362 * Subsystem_Fy_B.t356 * 2.0;
  Subsystem_Fy_B.t190 = Subsystem_Fy_B.t362 * Subsystem_Fy_B.t241 * 2.0;
  Subsystem_Fy_B.t191 = Subsystem_Fy_B.t357 * Subsystem_Fy_B.t356 * 2.0;
  Subsystem_Fy_B.t192 = Subsystem_Fy_B.t357 * Subsystem_Fy_B.t241 * 2.0;
  Subsystem_Fy_B.t205 = Subsystem_Fy_B.t356 * Subsystem_Fy_B.t241 * 2.0;
  Subsystem_Fy_B.t162 = ((Subsystem_Fy_B.t226 / 2.0 + Subsystem_Fy_B.t269 / 2.0)
    + Subsystem_Fy_B.t279 / 2.0) + -(Subsystem_Fy_B.t162 * Subsystem_Fy_B.t185 /
    2.0);
  Subsystem_Fy_B.t161 = ((Subsystem_Fy_B.t227 / 2.0 + Subsystem_Fy_B.t271 / 2.0)
    + Subsystem_Fy_B.t274 / 2.0) + -(Subsystem_Fy_B.t161 * Subsystem_Fy_B.t185 /
    2.0);
  Subsystem_Fy_B.t156 = ((Subsystem_Fy_B.t228 / 2.0 + Subsystem_Fy_B.t273 / 2.0)
    + Subsystem_Fy_B.t275 / 2.0) + -(Subsystem_Fy_B.t156 * Subsystem_Fy_B.t185 /
    2.0);
  Subsystem_Fy_B.t281 = ((Subsystem_Fy_B.t219 / 2.0 + Subsystem_Fy_B.t281 / 2.0)
    + Subsystem_Fy_B.t291 / 2.0) + -(Subsystem_Fy_B.t150 * Subsystem_Fy_B.t186 /
    2.0);
  Subsystem_Fy_B.t283 = ((Subsystem_Fy_B.t220 / 2.0 + Subsystem_Fy_B.t283 / 2.0)
    + Subsystem_Fy_B.t286 / 2.0) + -(Subsystem_Fy_B.t149 * Subsystem_Fy_B.t186 /
    2.0);
  Subsystem_Fy_B.t269 = ((Subsystem_Fy_B.t221 / 2.0 + Subsystem_Fy_B.t285 / 2.0)
    + Subsystem_Fy_B.t287 / 2.0) + -(Subsystem_Fy_B.t144 * Subsystem_Fy_B.t186 /
    2.0);
  Subsystem_Fy_B.t186 = Subsystem_Fy_B.t357 * Subsystem_Fy_B.t357 * 2.0;
  Subsystem_Fy_B.t144 = Subsystem_Fy_B.t356 * Subsystem_Fy_B.t356 * 2.0;
  Subsystem_Fy_B.t185 = Subsystem_Fy_B.t241 * Subsystem_Fy_B.t241 * 2.0;
  Subsystem_Fy_B.t149 = Subsystem_Fy_B.t376 * Subsystem_Fy_B.t376 * 2.0;
  Subsystem_Fy_B.t150 = Subsystem_Fy_B.t377 * Subsystem_Fy_B.t377 * 2.0;
  Subsystem_Fy_B.t285 = Subsystem_Fy_B.t363 * Subsystem_Fy_B.t363 * 2.0;
  Subsystem_Fy_B.t287 = ((Subsystem_Fy_B.t358 + Subsystem_Fy_B.t474) +
    Subsystem_Fy_B.t230) * 0.0025;
  Subsystem_Fy_B.t286 = ((Subsystem_Fy_B.t446 + Subsystem_Fy_B.t471) +
    Subsystem_Fy_B.t232) * 0.0025;
  Subsystem_Fy_B.t291 = ((Subsystem_Fy_B.t447 + Subsystem_Fy_B.t472) +
    Subsystem_Fy_B.t229) * 0.0025;
  Subsystem_Fy_B.t273 = ((Subsystem_Fy_B.t450 + Subsystem_Fy_B.t478) +
    Subsystem_Fy_B.t557) * 0.0025;
  Subsystem_Fy_B.t275 = ((Subsystem_Fy_B.t364 + Subsystem_Fy_B.t480) +
    Subsystem_Fy_B.t558) * 0.0025;
  Subsystem_Fy_B.t271 = ((Subsystem_Fy_B.t449 + Subsystem_Fy_B.t477) +
    Subsystem_Fy_B.t231) * 0.0025;
  Subsystem_Fy_B.t279 = 1.0 / std::sqrt(((Subsystem_Fy_B.t141 *
    Subsystem_Fy_B.t141 + Subsystem_Fy_B.t143 * Subsystem_Fy_B.t143) +
    Subsystem_Fy_B.t151 * Subsystem_Fy_B.t151) + Subsystem_Fy_B.t133 *
    Subsystem_Fy_B.t133);
  Subsystem_Fy_B.t274 = 1.0 / std::sqrt(((Subsystem_Fy_B.t153 *
    Subsystem_Fy_B.t153 + Subsystem_Fy_B.t155 * Subsystem_Fy_B.t155) +
    Subsystem_Fy_B.t157 * Subsystem_Fy_B.t157) + Subsystem_Fy_B.t137 *
    Subsystem_Fy_B.t137);
  Subsystem_Fy_B.t478 = ((Subsystem_Fy_B.t450 / 2.0 + Subsystem_Fy_B.t478 / 2.0)
    + Subsystem_Fy_B.t557 / 2.0) * 6.25E-6;
  Subsystem_Fy_B.t480 = ((Subsystem_Fy_B.t364 / 2.0 + Subsystem_Fy_B.t480 / 2.0)
    + Subsystem_Fy_B.t558 / 2.0) * 6.25E-6;
  Subsystem_Fy_B.t477 = ((Subsystem_Fy_B.t449 / 2.0 + Subsystem_Fy_B.t477 / 2.0)
    + Subsystem_Fy_B.t231 / 2.0) * 6.25E-6;
  Subsystem_Fy_B.t449 = ((Subsystem_Fy_B.t446 / 2.0 + Subsystem_Fy_B.t471 / 2.0)
    + Subsystem_Fy_B.t232 / 2.0) * 6.25E-6;
  Subsystem_Fy_B.t731 = ((Subsystem_Fy_B.t358 / 2.0 + Subsystem_Fy_B.t474 / 2.0)
    + Subsystem_Fy_B.t230 / 2.0) * 6.25E-6;
  Subsystem_Fy_B.t732 = ((Subsystem_Fy_B.t447 / 2.0 + Subsystem_Fy_B.t472 / 2.0)
    + Subsystem_Fy_B.t229 / 2.0) * 6.25E-6;
  Subsystem_Fy_B.t799 = 1.0 / std::sqrt(((Subsystem_Fy_B.t579 *
    Subsystem_Fy_B.t579 + Subsystem_Fy_B.t580 * Subsystem_Fy_B.t580) +
    Subsystem_Fy_B.t581 * Subsystem_Fy_B.t581) + Subsystem_Fy_B.t582 *
    Subsystem_Fy_B.t582);
  Subsystem_Fy_B.t801 = 1.0 / std::sqrt(((Subsystem_Fy_B.t575 *
    Subsystem_Fy_B.t575 + Subsystem_Fy_B.t576 * Subsystem_Fy_B.t576) +
    Subsystem_Fy_B.t577 * Subsystem_Fy_B.t577) + Subsystem_Fy_B.t578 *
    Subsystem_Fy_B.t578);
  Subsystem_Fy_B.t841 = ((0.0025 * Subsystem_Fy_B.t218 * Subsystem_Fy_B.t579 +
    0.0025 * Subsystem_Fy_B.t219 * Subsystem_Fy_B.t582) + 0.0025 *
    Subsystem_Fy_B.t221 * Subsystem_Fy_B.t580) + -(0.0025 * Subsystem_Fy_B.t220 *
    Subsystem_Fy_B.t581);
  Subsystem_Fy_B.t842 = ((0.0025 * Subsystem_Fy_B.t218 * Subsystem_Fy_B.t580 +
    0.0025 * Subsystem_Fy_B.t219 * Subsystem_Fy_B.t581) + 0.0025 *
    Subsystem_Fy_B.t220 * Subsystem_Fy_B.t582) + -(0.0025 * Subsystem_Fy_B.t221 *
    Subsystem_Fy_B.t579);
  Subsystem_Fy_B.t843 = ((0.0025 * Subsystem_Fy_B.t218 * Subsystem_Fy_B.t581 +
    0.0025 * Subsystem_Fy_B.t220 * Subsystem_Fy_B.t579) + 0.0025 *
    Subsystem_Fy_B.t221 * Subsystem_Fy_B.t582) + -(0.0025 * Subsystem_Fy_B.t219 *
    Subsystem_Fy_B.t580);
  Subsystem_Fy_B.t844 = ((0.0025 * Subsystem_Fy_B.t225 * Subsystem_Fy_B.t575 +
    0.0025 * Subsystem_Fy_B.t226 * Subsystem_Fy_B.t578) + 0.0025 *
    Subsystem_Fy_B.t228 * Subsystem_Fy_B.t576) + -(0.0025 * Subsystem_Fy_B.t227 *
    Subsystem_Fy_B.t577);
  Subsystem_Fy_B.t845 = ((0.0025 * Subsystem_Fy_B.t225 * Subsystem_Fy_B.t576 +
    0.0025 * Subsystem_Fy_B.t226 * Subsystem_Fy_B.t577) + 0.0025 *
    Subsystem_Fy_B.t227 * Subsystem_Fy_B.t578) + -(0.0025 * Subsystem_Fy_B.t228 *
    Subsystem_Fy_B.t575);
  Subsystem_Fy_B.t846 = ((0.0025 * Subsystem_Fy_B.t225 * Subsystem_Fy_B.t577 +
    0.0025 * Subsystem_Fy_B.t227 * Subsystem_Fy_B.t575) + 0.0025 *
    Subsystem_Fy_B.t228 * Subsystem_Fy_B.t578) + -(0.0025 * Subsystem_Fy_B.t226 *
    Subsystem_Fy_B.t576);
  Subsystem_Fy_B.t800 = rt_powd_snf(Subsystem_Fy_B.t799, 3.0);
  Subsystem_Fy_B.t802 = rt_powd_snf(Subsystem_Fy_B.t801, 3.0);
  Subsystem_Fy_B.t862 = 0.0025 * Subsystem_Fy_B.t226 * Subsystem_Fy_B.t153 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t863 = 0.0025 * Subsystem_Fy_B.t226 * Subsystem_Fy_B.t155 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t864 = 0.0025 * Subsystem_Fy_B.t226 * Subsystem_Fy_B.t157 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t865 = 0.0025 * Subsystem_Fy_B.t227 * Subsystem_Fy_B.t153 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t866 = 0.0025 * Subsystem_Fy_B.t227 * Subsystem_Fy_B.t155 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t867 = 0.0025 * Subsystem_Fy_B.t227 * Subsystem_Fy_B.t157 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t868 = 0.0025 * Subsystem_Fy_B.t228 * Subsystem_Fy_B.t153 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t869 = 0.0025 * Subsystem_Fy_B.t228 * Subsystem_Fy_B.t155 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t870 = 0.0025 * Subsystem_Fy_B.t228 * Subsystem_Fy_B.t157 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t874 = 0.0025 * Subsystem_Fy_B.t219 * Subsystem_Fy_B.t141 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t875 = 0.0025 * Subsystem_Fy_B.t219 * Subsystem_Fy_B.t143 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t876 = 0.0025 * Subsystem_Fy_B.t219 * Subsystem_Fy_B.t151 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t877 = 0.0025 * Subsystem_Fy_B.t220 * Subsystem_Fy_B.t141 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t878 = 0.0025 * Subsystem_Fy_B.t220 * Subsystem_Fy_B.t143 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t879 = 0.0025 * Subsystem_Fy_B.t220 * Subsystem_Fy_B.t151 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t880 = 0.0025 * Subsystem_Fy_B.t221 * Subsystem_Fy_B.t141 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t881 = 0.0025 * Subsystem_Fy_B.t221 * Subsystem_Fy_B.t143 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t882 = 0.0025 * Subsystem_Fy_B.t221 * Subsystem_Fy_B.t151 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t883 = 0.0025 * Subsystem_Fy_B.t218 * Subsystem_Fy_B.t133 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t884 = 0.0025 * Subsystem_Fy_B.t219 * Subsystem_Fy_B.t133 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t885 = 0.0025 * Subsystem_Fy_B.t220 * Subsystem_Fy_B.t133 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t886 = 0.0025 * Subsystem_Fy_B.t221 * Subsystem_Fy_B.t133 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799;
  Subsystem_Fy_B.t887 = 0.0025 * Subsystem_Fy_B.t225 * Subsystem_Fy_B.t137 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t888 = 0.0025 * Subsystem_Fy_B.t226 * Subsystem_Fy_B.t137 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t889 = 0.0025 * Subsystem_Fy_B.t227 * Subsystem_Fy_B.t137 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t890 = 0.0025 * Subsystem_Fy_B.t228 * Subsystem_Fy_B.t137 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801;
  Subsystem_Fy_B.t921 = Subsystem_Fy_B.t141 * Subsystem_Fy_B.t292 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t923 = Subsystem_Fy_B.t141 * Subsystem_Fy_B.t283 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t924 = Subsystem_Fy_B.t141 * Subsystem_Fy_B.t269 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t925 = Subsystem_Fy_B.t143 * Subsystem_Fy_B.t292 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t926 = Subsystem_Fy_B.t143 * Subsystem_Fy_B.t281 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t928 = Subsystem_Fy_B.t143 * Subsystem_Fy_B.t269 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t929 = Subsystem_Fy_B.t151 * Subsystem_Fy_B.t292 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t930 = Subsystem_Fy_B.t151 * Subsystem_Fy_B.t281 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t931 = Subsystem_Fy_B.t151 * Subsystem_Fy_B.t283 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t933 = Subsystem_Fy_B.t153 * Subsystem_Fy_B.t270 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t934 = Subsystem_Fy_B.t155 * Subsystem_Fy_B.t270 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t935 = Subsystem_Fy_B.t157 * Subsystem_Fy_B.t270 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t937 = Subsystem_Fy_B.t155 * Subsystem_Fy_B.t162 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t938 = Subsystem_Fy_B.t157 * Subsystem_Fy_B.t162 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t939 = Subsystem_Fy_B.t153 * Subsystem_Fy_B.t161 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t941 = Subsystem_Fy_B.t157 * Subsystem_Fy_B.t161 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t942 = Subsystem_Fy_B.t153 * Subsystem_Fy_B.t156 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t943 = Subsystem_Fy_B.t155 * Subsystem_Fy_B.t156 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t946 = Subsystem_Fy_B.t133 * Subsystem_Fy_B.t281 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t947 = Subsystem_Fy_B.t133 * Subsystem_Fy_B.t283 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t948 = Subsystem_Fy_B.t133 * Subsystem_Fy_B.t269 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.t950 = Subsystem_Fy_B.t137 * Subsystem_Fy_B.t162 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t951 = Subsystem_Fy_B.t137 * Subsystem_Fy_B.t161 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t952 = Subsystem_Fy_B.t137 * Subsystem_Fy_B.t156 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.t953 = ((Subsystem_Fy_B.t575 * Subsystem_Fy_B.t270 * 2.0 +
    Subsystem_Fy_B.t576 * Subsystem_Fy_B.t156 * 2.0) + -(Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t161 * 2.0)) + -(Subsystem_Fy_B.t578 * Subsystem_Fy_B.t162 *
    2.0);
  Subsystem_Fy_B.t954 = ((Subsystem_Fy_B.t575 * Subsystem_Fy_B.t161 * 2.0 +
    Subsystem_Fy_B.t577 * Subsystem_Fy_B.t270 * 2.0) + -(Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t162 * 2.0)) + -(Subsystem_Fy_B.t578 * Subsystem_Fy_B.t156 *
    2.0);
  Subsystem_Fy_B.t955 = ((Subsystem_Fy_B.t575 * Subsystem_Fy_B.t156 * 2.0 +
    Subsystem_Fy_B.t578 * Subsystem_Fy_B.t161 * 2.0) + -(Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t270 * 2.0)) + -(Subsystem_Fy_B.t577 * Subsystem_Fy_B.t162 *
    2.0);
  Subsystem_Fy_B.t956 = ((Subsystem_Fy_B.t579 * Subsystem_Fy_B.t292 * 2.0 +
    Subsystem_Fy_B.t580 * Subsystem_Fy_B.t269 * 2.0) + -(Subsystem_Fy_B.t581 *
    Subsystem_Fy_B.t283 * 2.0)) + -(Subsystem_Fy_B.t582 * Subsystem_Fy_B.t281 *
    2.0);
  Subsystem_Fy_B.t957 = ((Subsystem_Fy_B.t579 * Subsystem_Fy_B.t283 * 2.0 +
    Subsystem_Fy_B.t581 * Subsystem_Fy_B.t292 * 2.0) + -(Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t281 * 2.0)) + -(Subsystem_Fy_B.t582 * Subsystem_Fy_B.t269 *
    2.0);
  Subsystem_Fy_B.t958 = ((Subsystem_Fy_B.t579 * Subsystem_Fy_B.t269 * 2.0 +
    Subsystem_Fy_B.t582 * Subsystem_Fy_B.t283 * 2.0) + -(Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t292 * 2.0)) + -(Subsystem_Fy_B.t581 * Subsystem_Fy_B.t281 *
    2.0);
  Subsystem_Fy_B.t228 = Subsystem_Fy_B.t76 * Subsystem_Fy_B.t357;
  Subsystem_Fy_B.t227 = Subsystem_Fy_B.t77 * Subsystem_Fy_B.t356;
  Subsystem_Fy_B.t226 = Subsystem_Fy_B.t78 * Subsystem_Fy_B.t241;
  Subsystem_Fy_B.t647 = -(((Subsystem_Fy_B.t228 / 2.0 + Subsystem_Fy_B.t227 /
    2.0) + Subsystem_Fy_B.t226 / 2.0) * 0.0025 / 2.0) + x_nominal_prev[25];
  Subsystem_Fy_B.t221 = Subsystem_Fy_B.t76 * Subsystem_Fy_B.t362;
  Subsystem_Fy_B.t220 = Subsystem_Fy_B.t77 * Subsystem_Fy_B.t241;
  Subsystem_Fy_B.t219 = Subsystem_Fy_B.t78 * Subsystem_Fy_B.t356;
  Subsystem_Fy_B.t648 = ((Subsystem_Fy_B.t221 / 2.0 + Subsystem_Fy_B.t219 / 2.0)
    + -(Subsystem_Fy_B.t220 / 2.0)) * 0.0025 / 2.0 + x_nominal_prev[26];
  Subsystem_Fy_B.t472 = Subsystem_Fy_B.t77 * Subsystem_Fy_B.t362;
  Subsystem_Fy_B.t241 *= Subsystem_Fy_B.t76;
  Subsystem_Fy_B.t471 = Subsystem_Fy_B.t78 * Subsystem_Fy_B.t357;
  Subsystem_Fy_B.t649 = ((Subsystem_Fy_B.t241 / 2.0 + Subsystem_Fy_B.t472 / 2.0)
    + -(Subsystem_Fy_B.t471 / 2.0)) * 0.0025 / 2.0 + x_nominal_prev[27];
  Subsystem_Fy_B.t356 *= Subsystem_Fy_B.t76;
  Subsystem_Fy_B.t357 *= Subsystem_Fy_B.t77;
  Subsystem_Fy_B.t362 *= Subsystem_Fy_B.t78;
  Subsystem_Fy_B.t650 = ((Subsystem_Fy_B.t357 / 2.0 + Subsystem_Fy_B.t362 / 2.0)
    + -(Subsystem_Fy_B.t356 / 2.0)) * 0.0025 / 2.0 + x_nominal_prev[28];
  Subsystem_Fy_B.t474 = Subsystem_Fy_B.t70 * Subsystem_Fy_B.t376;
  Subsystem_Fy_B.t358 = Subsystem_Fy_B.t71 * Subsystem_Fy_B.t377;
  Subsystem_Fy_B.t447 = Subsystem_Fy_B.t72 * Subsystem_Fy_B.t363;
  Subsystem_Fy_B.t651 = -(((Subsystem_Fy_B.t474 / 2.0 + Subsystem_Fy_B.t358 /
    2.0) + Subsystem_Fy_B.t447 / 2.0) * 0.0025 / 2.0) + x_nominal_prev[6];
  Subsystem_Fy_B.t446 = Subsystem_Fy_B.t70 * Subsystem_Fy_B.t375;
  Subsystem_Fy_B.t276 = Subsystem_Fy_B.t71 * Subsystem_Fy_B.t363;
  Subsystem_Fy_B.t280 = Subsystem_Fy_B.t72 * Subsystem_Fy_B.t377;
  Subsystem_Fy_B.t652 = ((Subsystem_Fy_B.t446 / 2.0 + Subsystem_Fy_B.t280 / 2.0)
    + -(Subsystem_Fy_B.t276 / 2.0)) * 0.0025 / 2.0 + x_nominal_prev[7];
  Subsystem_Fy_B.t229 = Subsystem_Fy_B.t71 * Subsystem_Fy_B.t375;
  Subsystem_Fy_B.t363 *= Subsystem_Fy_B.t70;
  Subsystem_Fy_B.t364 = Subsystem_Fy_B.t72 * Subsystem_Fy_B.t376;
  Subsystem_Fy_B.t653 = ((Subsystem_Fy_B.t363 / 2.0 + Subsystem_Fy_B.t229 / 2.0)
    + -(Subsystem_Fy_B.t364 / 2.0)) * 0.0025 / 2.0 + x_nominal_prev[8];
  Subsystem_Fy_B.t377 *= Subsystem_Fy_B.t70;
  Subsystem_Fy_B.t376 *= Subsystem_Fy_B.t71;
  Subsystem_Fy_B.t375 *= Subsystem_Fy_B.t72;
  Subsystem_Fy_B.t654 = ((Subsystem_Fy_B.t376 / 2.0 + Subsystem_Fy_B.t375 / 2.0)
    + -(Subsystem_Fy_B.t377 / 2.0)) * 0.0025 / 2.0 + x_nominal_prev[9];
  Subsystem_Fy_B.t230 = Subsystem_Fy_B.t647 * Subsystem_Fy_B.t648 * 2.0;
  Subsystem_Fy_B.t231 = Subsystem_Fy_B.t647 * Subsystem_Fy_B.t649 * 2.0;
  Subsystem_Fy_B.t232 = Subsystem_Fy_B.t647 * Subsystem_Fy_B.t650 * 2.0;
  Subsystem_Fy_B.t450 = Subsystem_Fy_B.t648 * Subsystem_Fy_B.t649 * 2.0;
  Subsystem_Fy_B.t558 = Subsystem_Fy_B.t648 * Subsystem_Fy_B.t650 * 2.0;
  Subsystem_Fy_B.t557 = Subsystem_Fy_B.t649 * Subsystem_Fy_B.t650 * 2.0;
  Subsystem_Fy_B.t779 = Subsystem_Fy_B.t651 * Subsystem_Fy_B.t652 * 2.0;
  Subsystem_Fy_B.t780 = Subsystem_Fy_B.t651 * Subsystem_Fy_B.t653 * 2.0;
  Subsystem_Fy_B.t781 = Subsystem_Fy_B.t651 * Subsystem_Fy_B.t654 * 2.0;
  Subsystem_Fy_B.t782 = Subsystem_Fy_B.t652 * Subsystem_Fy_B.t653 * 2.0;
  Subsystem_Fy_B.t783 = Subsystem_Fy_B.t652 * Subsystem_Fy_B.t654 * 2.0;
  Subsystem_Fy_B.t784 = Subsystem_Fy_B.t653 * Subsystem_Fy_B.t654 * 2.0;
  Subsystem_Fy_B.t661 = Subsystem_Fy_B.t652 * Subsystem_Fy_B.t652 * 2.0;
  Subsystem_Fy_B.t662 = Subsystem_Fy_B.t653 * Subsystem_Fy_B.t653 * 2.0;
  Subsystem_Fy_B.t663 = Subsystem_Fy_B.t654 * Subsystem_Fy_B.t654 * 2.0;
  Subsystem_Fy_B.t664 = Subsystem_Fy_B.t648 * Subsystem_Fy_B.t648 * 2.0;
  Subsystem_Fy_B.t665 = Subsystem_Fy_B.t649 * Subsystem_Fy_B.t649 * 2.0;
  Subsystem_Fy_B.t666 = Subsystem_Fy_B.t650 * Subsystem_Fy_B.t650 * 2.0;
  Subsystem_Fy_B.t1079_tmp = Subsystem_Fy_B.t155 * Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1079_tmp_n = Subsystem_Fy_B.t153 * Subsystem_Fy_B.t578 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1079_tmp_p = Subsystem_Fy_B.t157 * Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1079_tmp_l = Subsystem_Fy_B.t137 * Subsystem_Fy_B.t575 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1079 = (((Subsystem_Fy_B.t1079_tmp * Subsystem_Fy_B.t844 +
    (((Subsystem_Fy_B.t866 + Subsystem_Fy_B.t870) + Subsystem_Fy_B.t887) +
     -Subsystem_Fy_B.t862)) + Subsystem_Fy_B.t1079_tmp_n * Subsystem_Fy_B.t844)
    + -(Subsystem_Fy_B.t1079_tmp_p * Subsystem_Fy_B.t844)) +
    -(Subsystem_Fy_B.t1079_tmp_l * Subsystem_Fy_B.t844);
  Subsystem_Fy_B.t1080_tmp = Subsystem_Fy_B.t157 * Subsystem_Fy_B.t575 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1080_tmp_j = Subsystem_Fy_B.t155 * Subsystem_Fy_B.t578 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1080_tmp_d = Subsystem_Fy_B.t153 * Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1080_tmp_g = Subsystem_Fy_B.t137 * Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t1080 = (((Subsystem_Fy_B.t1080_tmp * Subsystem_Fy_B.t845 +
    (((Subsystem_Fy_B.t862 + Subsystem_Fy_B.t870) + Subsystem_Fy_B.t887) +
     -Subsystem_Fy_B.t866)) + Subsystem_Fy_B.t1080_tmp_j * Subsystem_Fy_B.t845)
    + -(Subsystem_Fy_B.t1080_tmp_d * Subsystem_Fy_B.t845)) +
    -(Subsystem_Fy_B.t1080_tmp_g * Subsystem_Fy_B.t845);
  Subsystem_Fy_B.t576 = Subsystem_Fy_B.t153 * Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t578 = Subsystem_Fy_B.t157 * Subsystem_Fy_B.t578 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t575 = Subsystem_Fy_B.t155 * Subsystem_Fy_B.t575 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t577 = Subsystem_Fy_B.t137 * Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t802;
  Subsystem_Fy_B.t862 = (((Subsystem_Fy_B.t576 * Subsystem_Fy_B.t846 +
    (((Subsystem_Fy_B.t862 + Subsystem_Fy_B.t866) + Subsystem_Fy_B.t887) +
     -Subsystem_Fy_B.t870)) + Subsystem_Fy_B.t578 * Subsystem_Fy_B.t846) +
    -(Subsystem_Fy_B.t575 * Subsystem_Fy_B.t846)) + -(Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t846);
  Subsystem_Fy_B.t866 = Subsystem_Fy_B.t141 * Subsystem_Fy_B.t582 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t870 = Subsystem_Fy_B.t143 * Subsystem_Fy_B.t581 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t802 = Subsystem_Fy_B.t151 * Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t887 = Subsystem_Fy_B.t133 * Subsystem_Fy_B.t579 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t1082 = (((Subsystem_Fy_B.t866 * Subsystem_Fy_B.t841 +
    (((Subsystem_Fy_B.t878 + Subsystem_Fy_B.t882) + Subsystem_Fy_B.t883) +
     -Subsystem_Fy_B.t874)) + Subsystem_Fy_B.t870 * Subsystem_Fy_B.t841) +
    -(Subsystem_Fy_B.t802 * Subsystem_Fy_B.t841)) + -(Subsystem_Fy_B.t887 *
    Subsystem_Fy_B.t841);
  Subsystem_Fy_B.t1083_tmp = Subsystem_Fy_B.t141 * Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  t1083_tmp = Subsystem_Fy_B.t151 * Subsystem_Fy_B.t582 * Subsystem_Fy_B.t279 *
    Subsystem_Fy_B.t800;
  t1083_tmp_0 = Subsystem_Fy_B.t143 * Subsystem_Fy_B.t579 * Subsystem_Fy_B.t279 *
    Subsystem_Fy_B.t800;
  t1083_tmp_1 = Subsystem_Fy_B.t133 * Subsystem_Fy_B.t581 * Subsystem_Fy_B.t279 *
    Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t1083 = (((Subsystem_Fy_B.t1083_tmp * Subsystem_Fy_B.t843 +
    (((Subsystem_Fy_B.t874 + Subsystem_Fy_B.t878) + Subsystem_Fy_B.t883) +
     -Subsystem_Fy_B.t882)) + t1083_tmp * Subsystem_Fy_B.t843) + -(t1083_tmp_0 *
    Subsystem_Fy_B.t843)) + -(t1083_tmp_1 * Subsystem_Fy_B.t843);
  Subsystem_Fy_B.t582 = Subsystem_Fy_B.t143 * Subsystem_Fy_B.t582 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t579 = Subsystem_Fy_B.t151 * Subsystem_Fy_B.t579 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t581 = Subsystem_Fy_B.t141 * Subsystem_Fy_B.t581 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t580 = Subsystem_Fy_B.t133 * Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t800;
  Subsystem_Fy_B.t874 = (((Subsystem_Fy_B.t582 * Subsystem_Fy_B.t842 +
    (((Subsystem_Fy_B.t874 + Subsystem_Fy_B.t882) + Subsystem_Fy_B.t883) +
     -Subsystem_Fy_B.t878)) + Subsystem_Fy_B.t579 * Subsystem_Fy_B.t842) +
    -(Subsystem_Fy_B.t581 * Subsystem_Fy_B.t842)) + -(Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t842);
  Subsystem_Fy_B.t878 = ((((((0.0025 * Subsystem_Fy_B.t225 * Subsystem_Fy_B.t157
    * Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 + -Subsystem_Fy_B.t863) +
    -Subsystem_Fy_B.t865) + -Subsystem_Fy_B.t890) + Subsystem_Fy_B.t1079_tmp *
    Subsystem_Fy_B.t845) + Subsystem_Fy_B.t1079_tmp_n * Subsystem_Fy_B.t845) +
    -(Subsystem_Fy_B.t1079_tmp_p * Subsystem_Fy_B.t845)) +
    -(Subsystem_Fy_B.t1079_tmp_l * Subsystem_Fy_B.t845);
  Subsystem_Fy_B.t882 = ((((((0.0025 * Subsystem_Fy_B.t255 * Subsystem_Fy_B.t155
    * Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 + Subsystem_Fy_B.t889) +
    -Subsystem_Fy_B.t864) + -Subsystem_Fy_B.t868) + Subsystem_Fy_B.t1079_tmp *
    Subsystem_Fy_B.t846) + Subsystem_Fy_B.t1079_tmp_n * Subsystem_Fy_B.t846) +
    -(Subsystem_Fy_B.t1079_tmp_p * Subsystem_Fy_B.t846)) +
    -(Subsystem_Fy_B.t1079_tmp_l * Subsystem_Fy_B.t846);
  Subsystem_Fy_B.t863 = ((((((0.0025 * Subsystem_Fy_B.t255 * Subsystem_Fy_B.t157
    * Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 + Subsystem_Fy_B.t890) +
    -Subsystem_Fy_B.t863) + -Subsystem_Fy_B.t865) + Subsystem_Fy_B.t1080_tmp *
    Subsystem_Fy_B.t844) + Subsystem_Fy_B.t1080_tmp_j * Subsystem_Fy_B.t844) +
    -(Subsystem_Fy_B.t1080_tmp_d * Subsystem_Fy_B.t844)) +
    -(Subsystem_Fy_B.t1080_tmp_g * Subsystem_Fy_B.t844);
  Subsystem_Fy_B.t846 = ((((((0.0025 * Subsystem_Fy_B.t225 * Subsystem_Fy_B.t153
    * Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 + -Subsystem_Fy_B.t867) +
    -Subsystem_Fy_B.t869) + -Subsystem_Fy_B.t888) + Subsystem_Fy_B.t1080_tmp *
    Subsystem_Fy_B.t846) + Subsystem_Fy_B.t1080_tmp_j * Subsystem_Fy_B.t846) +
    -(Subsystem_Fy_B.t1080_tmp_d * Subsystem_Fy_B.t846)) +
    -(Subsystem_Fy_B.t1080_tmp_g * Subsystem_Fy_B.t846);
  Subsystem_Fy_B.t225 = ((((((0.0025 * Subsystem_Fy_B.t225 * Subsystem_Fy_B.t155
    * Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 + -Subsystem_Fy_B.t864) +
    -Subsystem_Fy_B.t868) + -Subsystem_Fy_B.t889) + Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t844) + Subsystem_Fy_B.t578 * Subsystem_Fy_B.t844) +
    -(Subsystem_Fy_B.t575 * Subsystem_Fy_B.t844)) + -(Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t844);
  Subsystem_Fy_B.t888 = ((((((0.0025 * Subsystem_Fy_B.t255 * Subsystem_Fy_B.t153
    * Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 + Subsystem_Fy_B.t888) +
    -Subsystem_Fy_B.t867) + -Subsystem_Fy_B.t869) + Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t845) + Subsystem_Fy_B.t578 * Subsystem_Fy_B.t845) +
    -(Subsystem_Fy_B.t575 * Subsystem_Fy_B.t845)) + -(Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t845);
  Subsystem_Fy_B.t883 = ((((((0.0025 * Subsystem_Fy_B.t218 * Subsystem_Fy_B.t143
    * Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 + -Subsystem_Fy_B.t876) +
    -Subsystem_Fy_B.t880) + -Subsystem_Fy_B.t885) + Subsystem_Fy_B.t1083_tmp *
    Subsystem_Fy_B.t841) + t1083_tmp * Subsystem_Fy_B.t841) + -(t1083_tmp_0 *
    Subsystem_Fy_B.t841)) + -(t1083_tmp_1 * Subsystem_Fy_B.t841);
  Subsystem_Fy_B.t889 = ((((((0.0025 * Subsystem_Fy_B.t254 * Subsystem_Fy_B.t141
    * Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 + Subsystem_Fy_B.t884) +
    -Subsystem_Fy_B.t879) + -Subsystem_Fy_B.t881) + Subsystem_Fy_B.t1083_tmp *
    Subsystem_Fy_B.t842) + t1083_tmp * Subsystem_Fy_B.t842) + -(t1083_tmp_0 *
    Subsystem_Fy_B.t842)) + -(t1083_tmp_1 * Subsystem_Fy_B.t842);
  Subsystem_Fy_B.t800 = ((((((0.0025 * Subsystem_Fy_B.t218 * Subsystem_Fy_B.t151
    * Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 + -Subsystem_Fy_B.t875) +
    -Subsystem_Fy_B.t877) + -Subsystem_Fy_B.t886) + Subsystem_Fy_B.t866 *
    Subsystem_Fy_B.t842) + Subsystem_Fy_B.t870 * Subsystem_Fy_B.t842) +
    -(Subsystem_Fy_B.t802 * Subsystem_Fy_B.t842)) + -(Subsystem_Fy_B.t887 *
    Subsystem_Fy_B.t842);
  Subsystem_Fy_B.t1094 = ((((((0.0025 * Subsystem_Fy_B.t254 *
    Subsystem_Fy_B.t143 * Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 +
    Subsystem_Fy_B.t885) + -Subsystem_Fy_B.t876) + -Subsystem_Fy_B.t880) +
    Subsystem_Fy_B.t866 * Subsystem_Fy_B.t843) + Subsystem_Fy_B.t870 *
    Subsystem_Fy_B.t843) + -(Subsystem_Fy_B.t802 * Subsystem_Fy_B.t843)) +
    -(Subsystem_Fy_B.t887 * Subsystem_Fy_B.t843);
  Subsystem_Fy_B.t1095 = ((((((0.0025 * Subsystem_Fy_B.t254 *
    Subsystem_Fy_B.t151 * Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 +
    Subsystem_Fy_B.t886) + -Subsystem_Fy_B.t875) + -Subsystem_Fy_B.t877) +
    Subsystem_Fy_B.t582 * Subsystem_Fy_B.t841) + Subsystem_Fy_B.t579 *
    Subsystem_Fy_B.t841) + -(Subsystem_Fy_B.t581 * Subsystem_Fy_B.t841)) +
    -(Subsystem_Fy_B.t580 * Subsystem_Fy_B.t841);
  Subsystem_Fy_B.t1096 = ((((((0.0025 * Subsystem_Fy_B.t218 *
    Subsystem_Fy_B.t141 * Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 +
    -Subsystem_Fy_B.t879) + -Subsystem_Fy_B.t881) + -Subsystem_Fy_B.t884) +
    Subsystem_Fy_B.t582 * Subsystem_Fy_B.t843) + Subsystem_Fy_B.t579 *
    Subsystem_Fy_B.t843) + -(Subsystem_Fy_B.t581 * Subsystem_Fy_B.t843)) +
    -(Subsystem_Fy_B.t580 * Subsystem_Fy_B.t843);
  Subsystem_Fy_B.t218 = Subsystem_Fy_B.t76 * Subsystem_Fy_B.t647;
  Subsystem_Fy_B.t254 = Subsystem_Fy_B.t77 * Subsystem_Fy_B.t650;
  Subsystem_Fy_B.t255 = Subsystem_Fy_B.t78 * Subsystem_Fy_B.t649;
  Subsystem_Fy_B.t843 = ((Subsystem_Fy_B.t218 / 2.0 + Subsystem_Fy_B.t255 / 2.0)
    + -(Subsystem_Fy_B.t254 / 2.0)) * 0.0025 + x_nominal_prev[26];
  Subsystem_Fy_B.t879 = Subsystem_Fy_B.t77 * Subsystem_Fy_B.t647;
  Subsystem_Fy_B.t881 = Subsystem_Fy_B.t76 * Subsystem_Fy_B.t650;
  Subsystem_Fy_B.t841 = Subsystem_Fy_B.t78 * Subsystem_Fy_B.t648;
  Subsystem_Fy_B.t875 = ((Subsystem_Fy_B.t879 / 2.0 + Subsystem_Fy_B.t881 / 2.0)
    + -(Subsystem_Fy_B.t841 / 2.0)) * 0.0025 + x_nominal_prev[27];
  Subsystem_Fy_B.t877 = Subsystem_Fy_B.t76 * Subsystem_Fy_B.t649;
  Subsystem_Fy_B.t884 = Subsystem_Fy_B.t77 * Subsystem_Fy_B.t648;
  Subsystem_Fy_B.t647 *= Subsystem_Fy_B.t78;
  Subsystem_Fy_B.t876 = ((Subsystem_Fy_B.t884 / 2.0 + Subsystem_Fy_B.t647 / 2.0)
    + -(Subsystem_Fy_B.t877 / 2.0)) * 0.0025 + x_nominal_prev[28];
  Subsystem_Fy_B.t880 = Subsystem_Fy_B.t70 * Subsystem_Fy_B.t651;
  Subsystem_Fy_B.t842 = Subsystem_Fy_B.t71 * Subsystem_Fy_B.t654;
  Subsystem_Fy_B.t885 = Subsystem_Fy_B.t72 * Subsystem_Fy_B.t653;
  Subsystem_Fy_B.t886 = ((Subsystem_Fy_B.t880 / 2.0 + Subsystem_Fy_B.t885 / 2.0)
    + -(Subsystem_Fy_B.t842 / 2.0)) * 0.0025 + x_nominal_prev[7];
  Subsystem_Fy_B.t845 = Subsystem_Fy_B.t71 * Subsystem_Fy_B.t651;
  Subsystem_Fy_B.t867 = Subsystem_Fy_B.t70 * Subsystem_Fy_B.t654;
  Subsystem_Fy_B.t869 = Subsystem_Fy_B.t72 * Subsystem_Fy_B.t652;
  Subsystem_Fy_B.t844 = ((Subsystem_Fy_B.t845 / 2.0 + Subsystem_Fy_B.t867 / 2.0)
    + -(Subsystem_Fy_B.t869 / 2.0)) * 0.0025 + x_nominal_prev[8];
  Subsystem_Fy_B.t864 = Subsystem_Fy_B.t70 * Subsystem_Fy_B.t653;
  Subsystem_Fy_B.t868 = Subsystem_Fy_B.t71 * Subsystem_Fy_B.t652;
  Subsystem_Fy_B.t651 *= Subsystem_Fy_B.t72;
  Subsystem_Fy_B.t865 = ((Subsystem_Fy_B.t868 / 2.0 + Subsystem_Fy_B.t651 / 2.0)
    + -(Subsystem_Fy_B.t864 / 2.0)) * 0.0025 + x_nominal_prev[9];
  Subsystem_Fy_B.t648 *= Subsystem_Fy_B.t76;
  Subsystem_Fy_B.t649 *= Subsystem_Fy_B.t77;
  Subsystem_Fy_B.t650 *= Subsystem_Fy_B.t78;
  Subsystem_Fy_B.t890 = -(((Subsystem_Fy_B.t648 / 2.0 + Subsystem_Fy_B.t649 /
    2.0) + Subsystem_Fy_B.t650 / 2.0) * 0.0025) + x_nominal_prev[25];
  Subsystem_Fy_B.t652 *= Subsystem_Fy_B.t70;
  Subsystem_Fy_B.t653 *= Subsystem_Fy_B.t71;
  Subsystem_Fy_B.t654 *= Subsystem_Fy_B.t72;
  Subsystem_Fy_B.t828 = -(((Subsystem_Fy_B.t652 / 2.0 + Subsystem_Fy_B.t653 /
    2.0) + Subsystem_Fy_B.t654 / 2.0) * 0.0025) + x_nominal_prev[6];
  Subsystem_Fy_B.t959 = Subsystem_Fy_B.t843 * Subsystem_Fy_B.t875 * 2.0;
  Subsystem_Fy_B.t960 = Subsystem_Fy_B.t843 * Subsystem_Fy_B.t876 * 2.0;
  Subsystem_Fy_B.t961 = Subsystem_Fy_B.t875 * Subsystem_Fy_B.t876 * 2.0;
  Subsystem_Fy_B.t962 = Subsystem_Fy_B.t886 * Subsystem_Fy_B.t844 * 2.0;
  Subsystem_Fy_B.t963 = Subsystem_Fy_B.t886 * Subsystem_Fy_B.t865 * 2.0;
  Subsystem_Fy_B.t964 = Subsystem_Fy_B.t844 * Subsystem_Fy_B.t865 * 2.0;
  Subsystem_Fy_B.t965 = Subsystem_Fy_B.t843 * Subsystem_Fy_B.t890 * 2.0;
  Subsystem_Fy_B.t966 = Subsystem_Fy_B.t875 * Subsystem_Fy_B.t890 * 2.0;
  Subsystem_Fy_B.t967 = Subsystem_Fy_B.t876 * Subsystem_Fy_B.t890 * 2.0;
  Subsystem_Fy_B.t968 = Subsystem_Fy_B.t886 * Subsystem_Fy_B.t828 * 2.0;
  Subsystem_Fy_B.t969 = Subsystem_Fy_B.t844 * Subsystem_Fy_B.t828 * 2.0;
  Subsystem_Fy_B.t970 = Subsystem_Fy_B.t865 * Subsystem_Fy_B.t828 * 2.0;
  Subsystem_Fy_B.t835 = Subsystem_Fy_B.t843 * Subsystem_Fy_B.t843 * 2.0;
  Subsystem_Fy_B.t836 = Subsystem_Fy_B.t875 * Subsystem_Fy_B.t875 * 2.0;
  Subsystem_Fy_B.t837 = Subsystem_Fy_B.t876 * Subsystem_Fy_B.t876 * 2.0;
  Subsystem_Fy_B.t838 = Subsystem_Fy_B.t886 * Subsystem_Fy_B.t886 * 2.0;
  Subsystem_Fy_B.t839 = Subsystem_Fy_B.t844 * Subsystem_Fy_B.t844 * 2.0;
  Subsystem_Fy_B.t840 = Subsystem_Fy_B.t865 * Subsystem_Fy_B.t865 * 2.0;
  Subsystem_Fy_B.b[0] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[41] = 1.0;
  std::memset(&Subsystem_Fy_B.b[42], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[82] = 1.0;
  std::memset(&Subsystem_Fy_B.b[83], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.b[120] = 0.0025;
  Subsystem_Fy_B.b[121] = 0.0;
  Subsystem_Fy_B.b[122] = 0.0;
  Subsystem_Fy_B.b[123] = 1.0;
  std::memset(&Subsystem_Fy_B.b[124], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.b[161] = 0.0025;
  Subsystem_Fy_B.b[162] = 0.0;
  Subsystem_Fy_B.b[163] = 0.0;
  Subsystem_Fy_B.b[164] = 1.0;
  std::memset(&Subsystem_Fy_B.b[165], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.b[202] = 0.0025;
  Subsystem_Fy_B.b[203] = 0.0;
  Subsystem_Fy_B.b[204] = 0.0;
  Subsystem_Fy_B.b[205] = 1.0;
  std::memset(&Subsystem_Fy_B.b[206], 0, 35U * sizeof(real_T));
  Subsystem_Fy_B.b[241] = -Subsystem_Fy_B.t477;
  Subsystem_Fy_B.b[242] = Subsystem_Fy_B.t480;
  Subsystem_Fy_B.b[243] = 0.0;
  Subsystem_Fy_B.b[244] = -Subsystem_Fy_B.t271;
  Subsystem_Fy_B.b[245] = Subsystem_Fy_B.t275;
  Subsystem_Fy_B.t133 = ((Subsystem_Fy_B.t141 * Subsystem_Fy_B.t281 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0 + Subsystem_Fy_B.t143 *
    Subsystem_Fy_B.t283 * Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0) +
    Subsystem_Fy_B.t151 * Subsystem_Fy_B.t269 * Subsystem_Fy_B.t279 *
    Subsystem_Fy_B.t799 * 2.0) + Subsystem_Fy_B.t133 * Subsystem_Fy_B.t292 *
    Subsystem_Fy_B.t279 * Subsystem_Fy_B.t799 * 2.0;
  Subsystem_Fy_B.b[246] = (((Subsystem_Fy_B.t866 * Subsystem_Fy_B.t956 +
    Subsystem_Fy_B.t133) + Subsystem_Fy_B.t870 * Subsystem_Fy_B.t956) -
    Subsystem_Fy_B.t802 * Subsystem_Fy_B.t956) - Subsystem_Fy_B.t887 *
    Subsystem_Fy_B.t956;
  Subsystem_Fy_B.b[247] = ((((((-Subsystem_Fy_B.t923 + Subsystem_Fy_B.t926) -
    Subsystem_Fy_B.t929) + Subsystem_Fy_B.t948) - Subsystem_Fy_B.t581 *
    Subsystem_Fy_B.t956) + Subsystem_Fy_B.t579 * Subsystem_Fy_B.t956) +
    Subsystem_Fy_B.t582 * Subsystem_Fy_B.t956) - Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t956;
  Subsystem_Fy_B.b[248] = (((Subsystem_Fy_B.t1083_tmp * Subsystem_Fy_B.t956 +
    (((-Subsystem_Fy_B.t924 + Subsystem_Fy_B.t925) + Subsystem_Fy_B.t930) -
     Subsystem_Fy_B.t947)) - t1083_tmp_0 * Subsystem_Fy_B.t956) + t1083_tmp *
    Subsystem_Fy_B.t956) - t1083_tmp_1 * Subsystem_Fy_B.t956;
  std::memset(&Subsystem_Fy_B.b[249], 0, 31U * sizeof(real_T));
  Subsystem_Fy_B.b[280] = Subsystem_Fy_B.t477;
  Subsystem_Fy_B.b[281] = 0.0;
  Subsystem_Fy_B.b[282] = -Subsystem_Fy_B.t478;
  Subsystem_Fy_B.b[283] = Subsystem_Fy_B.t271;
  Subsystem_Fy_B.b[284] = 0.0;
  Subsystem_Fy_B.b[285] = -Subsystem_Fy_B.t273;
  Subsystem_Fy_B.b[286] = ((((((Subsystem_Fy_B.t923 - Subsystem_Fy_B.t926) +
    Subsystem_Fy_B.t929) - Subsystem_Fy_B.t948) - Subsystem_Fy_B.t866 *
    Subsystem_Fy_B.t958) - Subsystem_Fy_B.t870 * Subsystem_Fy_B.t958) +
    Subsystem_Fy_B.t802 * Subsystem_Fy_B.t958) + Subsystem_Fy_B.t887 *
    Subsystem_Fy_B.t958;
  Subsystem_Fy_B.b[287] = (((Subsystem_Fy_B.t581 * Subsystem_Fy_B.t958 +
    Subsystem_Fy_B.t133) - Subsystem_Fy_B.t579 * Subsystem_Fy_B.t958) -
    Subsystem_Fy_B.t582 * Subsystem_Fy_B.t958) + Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t958;
  Subsystem_Fy_B.b[288] = ((((((-Subsystem_Fy_B.t921 - Subsystem_Fy_B.t928) +
    Subsystem_Fy_B.t931) + Subsystem_Fy_B.t946) - Subsystem_Fy_B.t1083_tmp *
    Subsystem_Fy_B.t958) + t1083_tmp_0 * Subsystem_Fy_B.t958) - t1083_tmp *
    Subsystem_Fy_B.t958) + t1083_tmp_1 * Subsystem_Fy_B.t958;
  std::memset(&Subsystem_Fy_B.b[289], 0, 31U * sizeof(real_T));
  Subsystem_Fy_B.b[320] = -Subsystem_Fy_B.t480;
  Subsystem_Fy_B.b[321] = Subsystem_Fy_B.t478;
  Subsystem_Fy_B.b[322] = 0.0;
  Subsystem_Fy_B.b[323] = -Subsystem_Fy_B.t275;
  Subsystem_Fy_B.b[324] = Subsystem_Fy_B.t273;
  Subsystem_Fy_B.b[325] = 0.0;
  Subsystem_Fy_B.b[326] = (((Subsystem_Fy_B.t866 * Subsystem_Fy_B.t957 +
    (((Subsystem_Fy_B.t924 - Subsystem_Fy_B.t925) - Subsystem_Fy_B.t930) +
     Subsystem_Fy_B.t947)) + Subsystem_Fy_B.t870 * Subsystem_Fy_B.t957) -
    Subsystem_Fy_B.t802 * Subsystem_Fy_B.t957) - Subsystem_Fy_B.t887 *
    Subsystem_Fy_B.t957;
  Subsystem_Fy_B.b[327] = ((((((Subsystem_Fy_B.t921 + Subsystem_Fy_B.t928) -
    Subsystem_Fy_B.t931) - Subsystem_Fy_B.t946) - Subsystem_Fy_B.t581 *
    Subsystem_Fy_B.t957) + Subsystem_Fy_B.t579 * Subsystem_Fy_B.t957) +
    Subsystem_Fy_B.t582 * Subsystem_Fy_B.t957) - Subsystem_Fy_B.t580 *
    Subsystem_Fy_B.t957;
  Subsystem_Fy_B.b[328] = (((Subsystem_Fy_B.t1083_tmp * Subsystem_Fy_B.t957 +
    Subsystem_Fy_B.t133) - t1083_tmp_0 * Subsystem_Fy_B.t957) + t1083_tmp *
    Subsystem_Fy_B.t957) - t1083_tmp_1 * Subsystem_Fy_B.t957;
  std::memset(&Subsystem_Fy_B.b[329], 0, 31U * sizeof(real_T));
  Subsystem_Fy_B.b[360] = Subsystem_Fy_B.t146;
  Subsystem_Fy_B.b[361] = -Subsystem_Fy_B.t10;
  Subsystem_Fy_B.b[362] = Subsystem_Fy_B.t196;
  Subsystem_Fy_B.b[363] = Subsystem_Fy_B.t249;
  Subsystem_Fy_B.b[364] = -Subsystem_Fy_B.t25;
  Subsystem_Fy_B.b[365] = Subsystem_Fy_B.t169;
  Subsystem_Fy_B.b[366] = 0.0;
  Subsystem_Fy_B.b[367] = 0.0;
  Subsystem_Fy_B.b[368] = 0.0;
  Subsystem_Fy_B.b[369] = 1.0;
  std::memset(&Subsystem_Fy_B.b[370], 0, 30U * sizeof(real_T));
  Subsystem_Fy_B.b[400] = Subsystem_Fy_B.t197;
  Subsystem_Fy_B.b[401] = Subsystem_Fy_B.t147;
  Subsystem_Fy_B.b[402] = -Subsystem_Fy_B.t16;
  Subsystem_Fy_B.b[403] = Subsystem_Fy_B.t18;
  Subsystem_Fy_B.b[404] = Subsystem_Fy_B.t3;
  Subsystem_Fy_B.b[405] = -Subsystem_Fy_B.t164;
  Subsystem_Fy_B.b[406] = 0.0;
  Subsystem_Fy_B.b[407] = 0.0;
  Subsystem_Fy_B.b[408] = 0.0;
  Subsystem_Fy_B.b[409] = 0.0;
  Subsystem_Fy_B.b[410] = 1.0;
  std::memset(&Subsystem_Fy_B.b[411], 0, 29U * sizeof(real_T));
  Subsystem_Fy_B.b[440] = -Subsystem_Fy_B.t15;
  Subsystem_Fy_B.b[441] = Subsystem_Fy_B.t187;
  Subsystem_Fy_B.b[442] = Subsystem_Fy_B.t145;
  Subsystem_Fy_B.b[443] = -Subsystem_Fy_B.t24;
  Subsystem_Fy_B.b[444] = Subsystem_Fy_B.t17;
  Subsystem_Fy_B.b[445] = Subsystem_Fy_B.t4;
  Subsystem_Fy_B.b[446] = 0.0;
  Subsystem_Fy_B.b[447] = 0.0;
  Subsystem_Fy_B.b[448] = 0.0;
  Subsystem_Fy_B.b[449] = 0.0;
  Subsystem_Fy_B.b[450] = 0.0;
  Subsystem_Fy_B.b[451] = 1.0;
  std::memset(&Subsystem_Fy_B.b[452], 0, 34U * sizeof(real_T));
  Subsystem_Fy_B.b[486] = Subsystem_Fy_B.t1082;
  Subsystem_Fy_B.b[487] = Subsystem_Fy_B.t1095;
  Subsystem_Fy_B.b[488] = Subsystem_Fy_B.t883;
  Subsystem_Fy_B.b[489] = 0.0;
  Subsystem_Fy_B.b[490] = 0.0;
  Subsystem_Fy_B.b[491] = 0.0;
  Subsystem_Fy_B.b[492] = 1.0;
  std::memset(&Subsystem_Fy_B.b[493], 0, 33U * sizeof(real_T));
  Subsystem_Fy_B.b[526] = Subsystem_Fy_B.t800;
  Subsystem_Fy_B.b[527] = Subsystem_Fy_B.t874;
  Subsystem_Fy_B.b[528] = Subsystem_Fy_B.t889;
  Subsystem_Fy_B.b[529] = 0.0;
  Subsystem_Fy_B.b[530] = 0.0;
  Subsystem_Fy_B.b[531] = 0.0;
  Subsystem_Fy_B.b[532] = 0.0;
  Subsystem_Fy_B.b[533] = 1.0;
  std::memset(&Subsystem_Fy_B.b[534], 0, sizeof(real_T) << 5U);
  Subsystem_Fy_B.b[566] = Subsystem_Fy_B.t1094;
  Subsystem_Fy_B.b[567] = Subsystem_Fy_B.t1096;
  Subsystem_Fy_B.b[568] = Subsystem_Fy_B.t1083;
  Subsystem_Fy_B.b[569] = 0.0;
  Subsystem_Fy_B.b[570] = 0.0;
  Subsystem_Fy_B.b[571] = 0.0;
  Subsystem_Fy_B.b[572] = 0.0;
  Subsystem_Fy_B.b[573] = 0.0;
  Subsystem_Fy_B.b[574] = 1.0;
  std::memset(&Subsystem_Fy_B.b[575], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[615] = 1.0;
  std::memset(&Subsystem_Fy_B.b[616], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[656] = 1.0;
  std::memset(&Subsystem_Fy_B.b[657], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[697] = 1.0;
  std::memset(&Subsystem_Fy_B.b[698], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[738] = 1.0;
  std::memset(&Subsystem_Fy_B.b[739], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[779] = 1.0;
  std::memset(&Subsystem_Fy_B.b[780], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[820] = 1.0;
  std::memset(&Subsystem_Fy_B.b[821], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.b[858] = 0.0025;
  Subsystem_Fy_B.b[859] = 0.0;
  Subsystem_Fy_B.b[860] = 0.0;
  Subsystem_Fy_B.b[861] = 1.0;
  std::memset(&Subsystem_Fy_B.b[862], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.b[899] = 0.0025;
  Subsystem_Fy_B.b[900] = 0.0;
  Subsystem_Fy_B.b[901] = 0.0;
  Subsystem_Fy_B.b[902] = 1.0;
  std::memset(&Subsystem_Fy_B.b[903], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.b[940] = 0.0025;
  Subsystem_Fy_B.b[941] = 0.0;
  Subsystem_Fy_B.b[942] = 0.0;
  Subsystem_Fy_B.b[943] = 1.0;
  std::memset(&Subsystem_Fy_B.b[944], 0, 35U * sizeof(real_T));
  Subsystem_Fy_B.b[979] = -Subsystem_Fy_B.t449;
  Subsystem_Fy_B.b[980] = Subsystem_Fy_B.t731;
  Subsystem_Fy_B.b[981] = 0.0;
  Subsystem_Fy_B.b[982] = -Subsystem_Fy_B.t286;
  Subsystem_Fy_B.b[983] = Subsystem_Fy_B.t287;
  Subsystem_Fy_B.t133 = ((Subsystem_Fy_B.t153 * Subsystem_Fy_B.t162 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0 + Subsystem_Fy_B.t155 *
    Subsystem_Fy_B.t161 * Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0) +
    Subsystem_Fy_B.t157 * Subsystem_Fy_B.t156 * Subsystem_Fy_B.t274 *
    Subsystem_Fy_B.t801 * 2.0) + Subsystem_Fy_B.t137 * Subsystem_Fy_B.t270 *
    Subsystem_Fy_B.t274 * Subsystem_Fy_B.t801 * 2.0;
  Subsystem_Fy_B.b[984] = (((Subsystem_Fy_B.t1079_tmp_n * Subsystem_Fy_B.t953 +
    Subsystem_Fy_B.t133) + Subsystem_Fy_B.t1079_tmp * Subsystem_Fy_B.t953) -
    Subsystem_Fy_B.t1079_tmp_p * Subsystem_Fy_B.t953) -
    Subsystem_Fy_B.t1079_tmp_l * Subsystem_Fy_B.t953;
  Subsystem_Fy_B.b[985] = ((((((-Subsystem_Fy_B.t935 + Subsystem_Fy_B.t937) -
    Subsystem_Fy_B.t939) + Subsystem_Fy_B.t952) - Subsystem_Fy_B.t1080_tmp_d *
    Subsystem_Fy_B.t953) + Subsystem_Fy_B.t1080_tmp * Subsystem_Fy_B.t953) +
    Subsystem_Fy_B.t1080_tmp_j * Subsystem_Fy_B.t953) -
    Subsystem_Fy_B.t1080_tmp_g * Subsystem_Fy_B.t953;
  Subsystem_Fy_B.b[986] = (((Subsystem_Fy_B.t576 * Subsystem_Fy_B.t953 +
    (((Subsystem_Fy_B.t934 + Subsystem_Fy_B.t938) - Subsystem_Fy_B.t942) -
     Subsystem_Fy_B.t951)) - Subsystem_Fy_B.t575 * Subsystem_Fy_B.t953) +
    Subsystem_Fy_B.t578 * Subsystem_Fy_B.t953) - Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t953;
  std::memset(&Subsystem_Fy_B.b[987], 0, 31U * sizeof(real_T));
  Subsystem_Fy_B.b[1018] = Subsystem_Fy_B.t449;
  Subsystem_Fy_B.b[1019] = 0.0;
  Subsystem_Fy_B.b[1020] = -Subsystem_Fy_B.t732;
  Subsystem_Fy_B.b[1021] = Subsystem_Fy_B.t286;
  Subsystem_Fy_B.b[1022] = 0.0;
  Subsystem_Fy_B.b[1023] = -Subsystem_Fy_B.t291;
  Subsystem_Fy_B.b[1024] = ((((((Subsystem_Fy_B.t935 - Subsystem_Fy_B.t937) +
    Subsystem_Fy_B.t939) - Subsystem_Fy_B.t952) - Subsystem_Fy_B.t1079_tmp_n *
    Subsystem_Fy_B.t955) - Subsystem_Fy_B.t1079_tmp * Subsystem_Fy_B.t955) +
    Subsystem_Fy_B.t1079_tmp_p * Subsystem_Fy_B.t955) +
    Subsystem_Fy_B.t1079_tmp_l * Subsystem_Fy_B.t955;
  Subsystem_Fy_B.b[1025] = (((Subsystem_Fy_B.t1080_tmp_d * Subsystem_Fy_B.t955 +
    Subsystem_Fy_B.t133) - Subsystem_Fy_B.t1080_tmp * Subsystem_Fy_B.t955) -
    Subsystem_Fy_B.t1080_tmp_j * Subsystem_Fy_B.t955) +
    Subsystem_Fy_B.t1080_tmp_g * Subsystem_Fy_B.t955;
  Subsystem_Fy_B.b[1026] = ((((((-Subsystem_Fy_B.t933 + Subsystem_Fy_B.t941) -
    Subsystem_Fy_B.t943) + Subsystem_Fy_B.t950) - Subsystem_Fy_B.t576 *
    Subsystem_Fy_B.t955) + Subsystem_Fy_B.t575 * Subsystem_Fy_B.t955) -
    Subsystem_Fy_B.t578 * Subsystem_Fy_B.t955) + Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t955;
  std::memset(&Subsystem_Fy_B.b[1027], 0, 31U * sizeof(real_T));
  Subsystem_Fy_B.b[1058] = -Subsystem_Fy_B.t731;
  Subsystem_Fy_B.b[1059] = Subsystem_Fy_B.t732;
  Subsystem_Fy_B.b[1060] = 0.0;
  Subsystem_Fy_B.b[1061] = -Subsystem_Fy_B.t287;
  Subsystem_Fy_B.b[1062] = Subsystem_Fy_B.t291;
  Subsystem_Fy_B.b[1063] = 0.0;
  Subsystem_Fy_B.b[1064] = (((Subsystem_Fy_B.t1079_tmp_n * Subsystem_Fy_B.t954 +
    (((-Subsystem_Fy_B.t934 - Subsystem_Fy_B.t938) + Subsystem_Fy_B.t942) +
     Subsystem_Fy_B.t951)) + Subsystem_Fy_B.t1079_tmp * Subsystem_Fy_B.t954) -
    Subsystem_Fy_B.t1079_tmp_p * Subsystem_Fy_B.t954) -
    Subsystem_Fy_B.t1079_tmp_l * Subsystem_Fy_B.t954;
  Subsystem_Fy_B.b[1065] = ((((((Subsystem_Fy_B.t933 - Subsystem_Fy_B.t941) +
    Subsystem_Fy_B.t943) - Subsystem_Fy_B.t950) - Subsystem_Fy_B.t1080_tmp_d *
    Subsystem_Fy_B.t954) + Subsystem_Fy_B.t1080_tmp * Subsystem_Fy_B.t954) +
    Subsystem_Fy_B.t1080_tmp_j * Subsystem_Fy_B.t954) -
    Subsystem_Fy_B.t1080_tmp_g * Subsystem_Fy_B.t954;
  Subsystem_Fy_B.b[1066] = (((Subsystem_Fy_B.t576 * Subsystem_Fy_B.t954 +
    Subsystem_Fy_B.t133) - Subsystem_Fy_B.t575 * Subsystem_Fy_B.t954) +
    Subsystem_Fy_B.t578 * Subsystem_Fy_B.t954) - Subsystem_Fy_B.t577 *
    Subsystem_Fy_B.t954;
  std::memset(&Subsystem_Fy_B.b[1067], 0, 31U * sizeof(real_T));
  Subsystem_Fy_B.b[1098] = Subsystem_Fy_B.t163;
  Subsystem_Fy_B.b[1099] = -Subsystem_Fy_B.t160;
  Subsystem_Fy_B.b[1100] = Subsystem_Fy_B.t194;
  Subsystem_Fy_B.b[1101] = Subsystem_Fy_B.t240;
  Subsystem_Fy_B.b[1102] = -Subsystem_Fy_B.t209;
  Subsystem_Fy_B.b[1103] = Subsystem_Fy_B.t23;
  Subsystem_Fy_B.b[1104] = 0.0;
  Subsystem_Fy_B.b[1105] = 0.0;
  Subsystem_Fy_B.b[1106] = 0.0;
  Subsystem_Fy_B.b[1107] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1108], 0, 30U * sizeof(real_T));
  Subsystem_Fy_B.b[1138] = Subsystem_Fy_B.t195;
  Subsystem_Fy_B.b[1139] = Subsystem_Fy_B.t158;
  Subsystem_Fy_B.b[1140] = -Subsystem_Fy_B.t154;
  Subsystem_Fy_B.b[1141] = Subsystem_Fy_B.t165;
  Subsystem_Fy_B.b[1142] = Subsystem_Fy_B.t239;
  Subsystem_Fy_B.b[1143] = -Subsystem_Fy_B.t14;
  Subsystem_Fy_B.b[1144] = 0.0;
  Subsystem_Fy_B.b[1145] = 0.0;
  Subsystem_Fy_B.b[1146] = 0.0;
  Subsystem_Fy_B.b[1147] = 0.0;
  Subsystem_Fy_B.b[1148] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1149], 0, 29U * sizeof(real_T));
  Subsystem_Fy_B.b[1178] = -Subsystem_Fy_B.t152;
  Subsystem_Fy_B.b[1179] = Subsystem_Fy_B.t193;
  Subsystem_Fy_B.b[1180] = Subsystem_Fy_B.t159;
  Subsystem_Fy_B.b[1181] = -Subsystem_Fy_B.t210;
  Subsystem_Fy_B.b[1182] = Subsystem_Fy_B.t19;
  Subsystem_Fy_B.b[1183] = Subsystem_Fy_B.t250;
  Subsystem_Fy_B.b[1184] = 0.0;
  Subsystem_Fy_B.b[1185] = 0.0;
  Subsystem_Fy_B.b[1186] = 0.0;
  Subsystem_Fy_B.b[1187] = 0.0;
  Subsystem_Fy_B.b[1188] = 0.0;
  Subsystem_Fy_B.b[1189] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1190], 0, 34U * sizeof(real_T));
  Subsystem_Fy_B.b[1224] = Subsystem_Fy_B.t1079;
  Subsystem_Fy_B.b[1225] = Subsystem_Fy_B.t863;
  Subsystem_Fy_B.b[1226] = Subsystem_Fy_B.t225;
  Subsystem_Fy_B.b[1227] = 0.0;
  Subsystem_Fy_B.b[1228] = 0.0;
  Subsystem_Fy_B.b[1229] = 0.0;
  Subsystem_Fy_B.b[1230] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1231], 0, 33U * sizeof(real_T));
  Subsystem_Fy_B.b[1264] = Subsystem_Fy_B.t878;
  Subsystem_Fy_B.b[1265] = Subsystem_Fy_B.t1080;
  Subsystem_Fy_B.b[1266] = Subsystem_Fy_B.t888;
  Subsystem_Fy_B.b[1267] = 0.0;
  Subsystem_Fy_B.b[1268] = 0.0;
  Subsystem_Fy_B.b[1269] = 0.0;
  Subsystem_Fy_B.b[1270] = 0.0;
  Subsystem_Fy_B.b[1271] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1272], 0, sizeof(real_T) << 5U);
  Subsystem_Fy_B.b[1304] = Subsystem_Fy_B.t882;
  Subsystem_Fy_B.b[1305] = Subsystem_Fy_B.t846;
  Subsystem_Fy_B.b[1306] = Subsystem_Fy_B.t862;
  Subsystem_Fy_B.b[1307] = 0.0;
  Subsystem_Fy_B.b[1308] = 0.0;
  Subsystem_Fy_B.b[1309] = 0.0;
  Subsystem_Fy_B.b[1310] = 0.0;
  Subsystem_Fy_B.b[1311] = 0.0;
  Subsystem_Fy_B.b[1312] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1313], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[1353] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1354], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[1394] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1395], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[1435] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1436], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[1476] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1477], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[1517] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1518], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[1558] = 1.0;
  std::memset(&Subsystem_Fy_B.b[1559], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.b[1599] = 1.0;
  Subsystem_Fy_B.c_c[0] = Subsystem_Fy_B.t146;
  Subsystem_Fy_B.c_c[1] = -Subsystem_Fy_B.t10;
  Subsystem_Fy_B.c_c[2] = Subsystem_Fy_B.t196;
  Subsystem_Fy_B.c_c[3] = Subsystem_Fy_B.t249;
  Subsystem_Fy_B.c_c[4] = -Subsystem_Fy_B.t25;
  Subsystem_Fy_B.c_c[5] = Subsystem_Fy_B.t169;
  std::memset(&Subsystem_Fy_B.c_c[6], 0, 34U * sizeof(real_T));
  Subsystem_Fy_B.c_c[40] = Subsystem_Fy_B.t197;
  Subsystem_Fy_B.c_c[41] = Subsystem_Fy_B.t147;
  Subsystem_Fy_B.c_c[42] = -Subsystem_Fy_B.t16;
  Subsystem_Fy_B.c_c[43] = Subsystem_Fy_B.t18;
  Subsystem_Fy_B.c_c[44] = Subsystem_Fy_B.t3;
  Subsystem_Fy_B.c_c[45] = -Subsystem_Fy_B.t164;
  std::memset(&Subsystem_Fy_B.c_c[46], 0, 34U * sizeof(real_T));
  Subsystem_Fy_B.c_c[80] = -Subsystem_Fy_B.t15;
  Subsystem_Fy_B.c_c[81] = Subsystem_Fy_B.t187;
  Subsystem_Fy_B.c_c[82] = Subsystem_Fy_B.t145;
  Subsystem_Fy_B.c_c[83] = -Subsystem_Fy_B.t24;
  Subsystem_Fy_B.c_c[84] = Subsystem_Fy_B.t17;
  Subsystem_Fy_B.c_c[85] = Subsystem_Fy_B.t4;
  std::memset(&Subsystem_Fy_B.c_c[86], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[126] = Subsystem_Fy_B.t1082;
  Subsystem_Fy_B.c_c[127] = Subsystem_Fy_B.t1095;
  Subsystem_Fy_B.c_c[128] = Subsystem_Fy_B.t883;
  std::memset(&Subsystem_Fy_B.c_c[129], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.c_c[166] = Subsystem_Fy_B.t800;
  Subsystem_Fy_B.c_c[167] = Subsystem_Fy_B.t874;
  Subsystem_Fy_B.c_c[168] = Subsystem_Fy_B.t889;
  std::memset(&Subsystem_Fy_B.c_c[169], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.c_c[206] = Subsystem_Fy_B.t1094;
  Subsystem_Fy_B.c_c[207] = Subsystem_Fy_B.t1096;
  Subsystem_Fy_B.c_c[208] = Subsystem_Fy_B.t1083;
  std::memset(&Subsystem_Fy_B.c_c[209], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[249] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[250], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[290] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[291], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[331] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[332], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[372] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[373], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[413] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[414], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[454] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[455], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[495] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[496], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[536] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[537], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[577] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[578], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[618] = Subsystem_Fy_B.t163;
  Subsystem_Fy_B.c_c[619] = -Subsystem_Fy_B.t160;
  Subsystem_Fy_B.c_c[620] = Subsystem_Fy_B.t194;
  Subsystem_Fy_B.c_c[621] = Subsystem_Fy_B.t240;
  Subsystem_Fy_B.c_c[622] = -Subsystem_Fy_B.t209;
  Subsystem_Fy_B.c_c[623] = Subsystem_Fy_B.t23;
  std::memset(&Subsystem_Fy_B.c_c[624], 0, 34U * sizeof(real_T));
  Subsystem_Fy_B.c_c[658] = Subsystem_Fy_B.t195;
  Subsystem_Fy_B.c_c[659] = Subsystem_Fy_B.t158;
  Subsystem_Fy_B.c_c[660] = -Subsystem_Fy_B.t154;
  Subsystem_Fy_B.c_c[661] = Subsystem_Fy_B.t165;
  Subsystem_Fy_B.c_c[662] = Subsystem_Fy_B.t239;
  Subsystem_Fy_B.c_c[663] = -Subsystem_Fy_B.t14;
  std::memset(&Subsystem_Fy_B.c_c[664], 0, 34U * sizeof(real_T));
  Subsystem_Fy_B.c_c[698] = -Subsystem_Fy_B.t152;
  Subsystem_Fy_B.c_c[699] = Subsystem_Fy_B.t193;
  Subsystem_Fy_B.c_c[700] = Subsystem_Fy_B.t159;
  Subsystem_Fy_B.c_c[701] = -Subsystem_Fy_B.t210;
  Subsystem_Fy_B.c_c[702] = Subsystem_Fy_B.t19;
  Subsystem_Fy_B.c_c[703] = Subsystem_Fy_B.t250;
  std::memset(&Subsystem_Fy_B.c_c[704], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[744] = Subsystem_Fy_B.t1079;
  Subsystem_Fy_B.c_c[745] = Subsystem_Fy_B.t863;
  Subsystem_Fy_B.c_c[746] = Subsystem_Fy_B.t225;
  std::memset(&Subsystem_Fy_B.c_c[747], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.c_c[784] = Subsystem_Fy_B.t878;
  Subsystem_Fy_B.c_c[785] = Subsystem_Fy_B.t1080;
  Subsystem_Fy_B.c_c[786] = Subsystem_Fy_B.t888;
  std::memset(&Subsystem_Fy_B.c_c[787], 0, 37U * sizeof(real_T));
  Subsystem_Fy_B.c_c[824] = Subsystem_Fy_B.t882;
  Subsystem_Fy_B.c_c[825] = Subsystem_Fy_B.t846;
  Subsystem_Fy_B.c_c[826] = Subsystem_Fy_B.t862;
  std::memset(&Subsystem_Fy_B.c_c[827], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[867] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[868], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[908] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[909], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[949] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[950], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[990] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[991], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[1031] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[1032], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[1072] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[1073], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[1113] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[1114], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[1154] = 0.0025;
  std::memset(&Subsystem_Fy_B.c_c[1155], 0, 40U * sizeof(real_T));
  Subsystem_Fy_B.c_c[1195] = 0.0025;
  Subsystem_Fy_B.c_c[1196] = 0.0;
  Subsystem_Fy_B.c_c[1197] = 0.0;
  Subsystem_Fy_B.c_c[1198] = 0.0;
  Subsystem_Fy_B.c_c[1199] = 0.0;
  Subsystem_Fy_B.t17 = Subsystem_Fy_B.t67 * ((Subsystem_Fy_B.t150 +
    Subsystem_Fy_B.t285) - 1.0);
  Subsystem_Fy_B.t169 = Subsystem_Fy_B.t69 * (Subsystem_Fy_B.t206 +
    Subsystem_Fy_B.t207);
  Subsystem_Fy_B.t18 = Subsystem_Fy_B.t68 * (Subsystem_Fy_B.t282 +
    -Subsystem_Fy_B.t288);
  Subsystem_Fy_B.t19 = Subsystem_Fy_B.t67 * ((Subsystem_Fy_B.t662 +
    Subsystem_Fy_B.t663) - 1.0);
  Subsystem_Fy_B.t23 = Subsystem_Fy_B.t69 * (Subsystem_Fy_B.t780 +
    Subsystem_Fy_B.t783);
  Subsystem_Fy_B.t165 = Subsystem_Fy_B.t68 * (Subsystem_Fy_B.t781 +
    -Subsystem_Fy_B.t782);
  NominalState_predict[0] = ((((((0.0 - Subsystem_Fy_B.t17) +
    Subsystem_Fy_B.t169) - Subsystem_Fy_B.t18) * 0.0025 + x_nominal_prev[3] *
    6.0) + (((0.0 - Subsystem_Fy_B.t19) + Subsystem_Fy_B.t23) -
            Subsystem_Fy_B.t165) * 0.0025) + ((Subsystem_Fy_B.t168 -
    Subsystem_Fy_B.t122) - Subsystem_Fy_B.t22) * 0.0025) * 0.0025 / 6.0 +
    x_nominal_prev[0];
  Subsystem_Fy_B.t133 = Subsystem_Fy_B.t68 * ((Subsystem_Fy_B.t149 +
    Subsystem_Fy_B.t285) - 1.0);
  Subsystem_Fy_B.t137 = Subsystem_Fy_B.t67 * (Subsystem_Fy_B.t282 +
    Subsystem_Fy_B.t288);
  Subsystem_Fy_B.t154 = Subsystem_Fy_B.t69 * (Subsystem_Fy_B.t208 +
    -Subsystem_Fy_B.t198);
  Subsystem_Fy_B.t152 = Subsystem_Fy_B.t68 * ((Subsystem_Fy_B.t661 +
    Subsystem_Fy_B.t663) - 1.0);
  Subsystem_Fy_B.t160 = Subsystem_Fy_B.t67 * (Subsystem_Fy_B.t781 +
    Subsystem_Fy_B.t782);
  Subsystem_Fy_B.t164 = Subsystem_Fy_B.t69 * (Subsystem_Fy_B.t779 +
    -Subsystem_Fy_B.t784);
  NominalState_predict[1] = ((((((-9.82 - Subsystem_Fy_B.t133) +
    Subsystem_Fy_B.t137) - Subsystem_Fy_B.t154) * 0.0025 + x_nominal_prev[4] *
    6.0) + (((-9.82 - Subsystem_Fy_B.t152) + Subsystem_Fy_B.t160) -
            Subsystem_Fy_B.t164) * 0.0025) + (((Subsystem_Fy_B.t166 + -9.82) -
    Subsystem_Fy_B.t121) - Subsystem_Fy_B.t21) * 0.0025) * 0.0025 / 6.0 +
    x_nominal_prev[1];
  Subsystem_Fy_B.t24 = Subsystem_Fy_B.t69 * ((Subsystem_Fy_B.t149 +
    Subsystem_Fy_B.t150) - 1.0);
  Subsystem_Fy_B.t25 = Subsystem_Fy_B.t68 * (Subsystem_Fy_B.t208 +
    Subsystem_Fy_B.t198);
  Subsystem_Fy_B.t16 = Subsystem_Fy_B.t67 * (Subsystem_Fy_B.t206 +
    -Subsystem_Fy_B.t207);
  Subsystem_Fy_B.t15 = Subsystem_Fy_B.t69 * ((Subsystem_Fy_B.t661 +
    Subsystem_Fy_B.t662) - 1.0);
  Subsystem_Fy_B.t10 = Subsystem_Fy_B.t68 * (Subsystem_Fy_B.t779 +
    Subsystem_Fy_B.t784);
  Subsystem_Fy_B.t14 = Subsystem_Fy_B.t67 * (Subsystem_Fy_B.t780 +
    -Subsystem_Fy_B.t783);
  NominalState_predict[2] = ((((((0.0 - Subsystem_Fy_B.t24) + Subsystem_Fy_B.t25)
    - Subsystem_Fy_B.t16) * 0.0025 + x_nominal_prev[5] * 6.0) + (((0.0 -
    Subsystem_Fy_B.t15) + Subsystem_Fy_B.t10) - Subsystem_Fy_B.t14) * 0.0025) +
    ((Subsystem_Fy_B.t167 - Subsystem_Fy_B.t123) - Subsystem_Fy_B.t20) * 0.0025)
    * 0.0025 / 6.0 + x_nominal_prev[2];
  NominalState_predict[3] = x_nominal_prev[3] - ((((((((((((-0.0 -
    Subsystem_Fy_B.t168) + Subsystem_Fy_B.t122) + Subsystem_Fy_B.t22) +
    ((Subsystem_Fy_B.t839 + Subsystem_Fy_B.t840) - 1.0) * Subsystem_Fy_B.t67) -
    (Subsystem_Fy_B.t963 + Subsystem_Fy_B.t969) * Subsystem_Fy_B.t69) +
    Subsystem_Fy_B.t17 * 2.0) - Subsystem_Fy_B.t169 * 2.0) + Subsystem_Fy_B.t18 *
    2.0) + Subsystem_Fy_B.t19 * 2.0) - Subsystem_Fy_B.t23 * 2.0) +
    Subsystem_Fy_B.t165 * 2.0) - (Subsystem_Fy_B.t962 - Subsystem_Fy_B.t970) *
    Subsystem_Fy_B.t68) * 0.0025 / 6.0;
  NominalState_predict[4] = x_nominal_prev[4] - ((((((((((((58.92 -
    Subsystem_Fy_B.t166) + Subsystem_Fy_B.t121) + Subsystem_Fy_B.t21) +
    ((Subsystem_Fy_B.t838 + Subsystem_Fy_B.t840) - 1.0) * Subsystem_Fy_B.t68) -
    (Subsystem_Fy_B.t962 + Subsystem_Fy_B.t970) * Subsystem_Fy_B.t67) +
    Subsystem_Fy_B.t133 * 2.0) - Subsystem_Fy_B.t137 * 2.0) +
    Subsystem_Fy_B.t154 * 2.0) + Subsystem_Fy_B.t152 * 2.0) -
    Subsystem_Fy_B.t160 * 2.0) + Subsystem_Fy_B.t164 * 2.0) -
    (Subsystem_Fy_B.t964 - Subsystem_Fy_B.t968) * Subsystem_Fy_B.t69) * 0.0025 /
    6.0;
  NominalState_predict[5] = x_nominal_prev[5] - ((((((((((((-0.0 -
    Subsystem_Fy_B.t167) + Subsystem_Fy_B.t123) + Subsystem_Fy_B.t20) +
    ((Subsystem_Fy_B.t838 + Subsystem_Fy_B.t839) - 1.0) * Subsystem_Fy_B.t69) -
    (Subsystem_Fy_B.t964 + Subsystem_Fy_B.t968) * Subsystem_Fy_B.t68) +
    Subsystem_Fy_B.t24 * 2.0) - Subsystem_Fy_B.t25 * 2.0) + Subsystem_Fy_B.t16 *
    2.0) + Subsystem_Fy_B.t15 * 2.0) - Subsystem_Fy_B.t10 * 2.0) +
    Subsystem_Fy_B.t14 * 2.0) - (Subsystem_Fy_B.t963 - Subsystem_Fy_B.t969) *
    Subsystem_Fy_B.t67) * 0.0025 / 6.0;
  NominalState_predict[6] = x_nominal_prev[6] - (((((((((Subsystem_Fy_B.t474 +
    Subsystem_Fy_B.t45) + Subsystem_Fy_B.t358) + Subsystem_Fy_B.t447) +
    Subsystem_Fy_B.t652) + Subsystem_Fy_B.t653) + Subsystem_Fy_B.t654) +
    Subsystem_Fy_B.t70 * Subsystem_Fy_B.t886 / 2.0) + Subsystem_Fy_B.t71 *
    Subsystem_Fy_B.t844 / 2.0) + Subsystem_Fy_B.t72 * Subsystem_Fy_B.t865 / 2.0)
    * 0.0025 / 6.0;
  NominalState_predict[7] = (((((((((Subsystem_Fy_B.t446 + Subsystem_Fy_B.t50) -
    Subsystem_Fy_B.t276) + Subsystem_Fy_B.t280) + Subsystem_Fy_B.t880) -
    Subsystem_Fy_B.t842) + Subsystem_Fy_B.t885) - Subsystem_Fy_B.t71 *
    Subsystem_Fy_B.t865 / 2.0) + Subsystem_Fy_B.t72 * Subsystem_Fy_B.t844 / 2.0)
    + Subsystem_Fy_B.t70 * Subsystem_Fy_B.t828 / 2.0) * 0.0025 / 6.0 +
    x_nominal_prev[7];
  NominalState_predict[8] = (((((((((Subsystem_Fy_B.t229 + Subsystem_Fy_B.t51) +
    Subsystem_Fy_B.t363) - Subsystem_Fy_B.t364) + Subsystem_Fy_B.t845) +
    Subsystem_Fy_B.t867) - Subsystem_Fy_B.t869) + Subsystem_Fy_B.t70 *
    Subsystem_Fy_B.t865 / 2.0) - Subsystem_Fy_B.t72 * Subsystem_Fy_B.t886 / 2.0)
    + Subsystem_Fy_B.t71 * Subsystem_Fy_B.t828 / 2.0) * 0.0025 / 6.0 +
    x_nominal_prev[8];
  NominalState_predict[9] = (((((((((Subsystem_Fy_B.t43 - Subsystem_Fy_B.t377) +
    Subsystem_Fy_B.t376) + Subsystem_Fy_B.t375) - Subsystem_Fy_B.t864) +
    Subsystem_Fy_B.t868) + Subsystem_Fy_B.t651) - Subsystem_Fy_B.t70 *
    Subsystem_Fy_B.t844 / 2.0) + Subsystem_Fy_B.t71 * Subsystem_Fy_B.t886 / 2.0)
    + Subsystem_Fy_B.t72 * Subsystem_Fy_B.t828 / 2.0) * 0.0025 / 6.0 +
    x_nominal_prev[9];
  std::memcpy(&NominalState_predict[10], &x_nominal_prev[10], 9U * sizeof(real_T));
  Subsystem_Fy_B.t17 = Subsystem_Fy_B.t73 * ((Subsystem_Fy_B.t144 +
    Subsystem_Fy_B.t185) - 1.0);
  Subsystem_Fy_B.t169 = Subsystem_Fy_B.t75 * (Subsystem_Fy_B.t189 +
    Subsystem_Fy_B.t192);
  Subsystem_Fy_B.t18 = Subsystem_Fy_B.t74 * (Subsystem_Fy_B.t190 +
    -Subsystem_Fy_B.t191);
  Subsystem_Fy_B.t19 = Subsystem_Fy_B.t73 * ((Subsystem_Fy_B.t665 +
    Subsystem_Fy_B.t666) - 1.0);
  Subsystem_Fy_B.t23 = Subsystem_Fy_B.t75 * (Subsystem_Fy_B.t231 +
    Subsystem_Fy_B.t558);
  Subsystem_Fy_B.t165 = Subsystem_Fy_B.t74 * (Subsystem_Fy_B.t232 +
    -Subsystem_Fy_B.t450);
  NominalState_predict[19] = ((((((0.0 - Subsystem_Fy_B.t17) +
    Subsystem_Fy_B.t169) - Subsystem_Fy_B.t18) * 0.0025 + x_nominal_prev[22] *
    6.0) + (((0.0 - Subsystem_Fy_B.t19) + Subsystem_Fy_B.t23) -
            Subsystem_Fy_B.t165) * 0.0025) + ((Subsystem_Fy_B.t172 -
    Subsystem_Fy_B.t131) - Subsystem_Fy_B.t13) * 0.0025) * 0.0025 / 6.0 +
    x_nominal_prev[19];
  Subsystem_Fy_B.t133 = Subsystem_Fy_B.t74 * ((Subsystem_Fy_B.t186 +
    Subsystem_Fy_B.t185) - 1.0);
  Subsystem_Fy_B.t137 = Subsystem_Fy_B.t73 * (Subsystem_Fy_B.t190 +
    Subsystem_Fy_B.t191);
  Subsystem_Fy_B.t154 = Subsystem_Fy_B.t75 * (Subsystem_Fy_B.t188 +
    -Subsystem_Fy_B.t205);
  Subsystem_Fy_B.t152 = Subsystem_Fy_B.t74 * ((Subsystem_Fy_B.t664 +
    Subsystem_Fy_B.t666) - 1.0);
  Subsystem_Fy_B.t160 = Subsystem_Fy_B.t73 * (Subsystem_Fy_B.t232 +
    Subsystem_Fy_B.t450);
  Subsystem_Fy_B.t164 = Subsystem_Fy_B.t75 * (Subsystem_Fy_B.t230 +
    -Subsystem_Fy_B.t557);
  NominalState_predict[20] = ((((((-9.82 - Subsystem_Fy_B.t133) +
    Subsystem_Fy_B.t137) - Subsystem_Fy_B.t154) * 0.0025 + x_nominal_prev[23] *
    6.0) + (((-9.82 - Subsystem_Fy_B.t152) + Subsystem_Fy_B.t160) -
            Subsystem_Fy_B.t164) * 0.0025) + (((Subsystem_Fy_B.t170 + -9.82) -
    Subsystem_Fy_B.t130) - Subsystem_Fy_B.t12) * 0.0025) * 0.0025 / 6.0 +
    x_nominal_prev[20];
  Subsystem_Fy_B.t24 = Subsystem_Fy_B.t75 * ((Subsystem_Fy_B.t186 +
    Subsystem_Fy_B.t144) - 1.0);
  Subsystem_Fy_B.t25 = Subsystem_Fy_B.t74 * (Subsystem_Fy_B.t188 +
    Subsystem_Fy_B.t205);
  Subsystem_Fy_B.t16 = Subsystem_Fy_B.t73 * (Subsystem_Fy_B.t189 +
    -Subsystem_Fy_B.t192);
  Subsystem_Fy_B.t15 = Subsystem_Fy_B.t75 * ((Subsystem_Fy_B.t664 +
    Subsystem_Fy_B.t665) - 1.0);
  Subsystem_Fy_B.t10 = Subsystem_Fy_B.t74 * (Subsystem_Fy_B.t230 +
    Subsystem_Fy_B.t557);
  Subsystem_Fy_B.t14 = Subsystem_Fy_B.t73 * (Subsystem_Fy_B.t231 +
    -Subsystem_Fy_B.t558);
  NominalState_predict[21] = ((((((0.0 - Subsystem_Fy_B.t24) +
    Subsystem_Fy_B.t25) - Subsystem_Fy_B.t16) * 0.0025 + x_nominal_prev[24] *
    6.0) + (((0.0 - Subsystem_Fy_B.t15) + Subsystem_Fy_B.t10) -
            Subsystem_Fy_B.t14) * 0.0025) + ((Subsystem_Fy_B.t171 -
    Subsystem_Fy_B.t132) - Subsystem_Fy_B.t11) * 0.0025) * 0.0025 / 6.0 +
    x_nominal_prev[21];
  NominalState_predict[22] = x_nominal_prev[22] - ((((((((((((-0.0 -
    Subsystem_Fy_B.t172) + Subsystem_Fy_B.t131) + Subsystem_Fy_B.t13) +
    ((Subsystem_Fy_B.t836 + Subsystem_Fy_B.t837) - 1.0) * Subsystem_Fy_B.t73) -
    (Subsystem_Fy_B.t960 + Subsystem_Fy_B.t966) * Subsystem_Fy_B.t75) +
    Subsystem_Fy_B.t17 * 2.0) - Subsystem_Fy_B.t169 * 2.0) + Subsystem_Fy_B.t18 *
    2.0) + Subsystem_Fy_B.t19 * 2.0) - Subsystem_Fy_B.t23 * 2.0) +
    Subsystem_Fy_B.t165 * 2.0) - (Subsystem_Fy_B.t959 - Subsystem_Fy_B.t967) *
    Subsystem_Fy_B.t74) * 0.0025 / 6.0;
  NominalState_predict[23] = x_nominal_prev[23] - ((((((((((((58.92 -
    Subsystem_Fy_B.t170) + Subsystem_Fy_B.t130) + Subsystem_Fy_B.t12) +
    ((Subsystem_Fy_B.t835 + Subsystem_Fy_B.t837) - 1.0) * Subsystem_Fy_B.t74) -
    (Subsystem_Fy_B.t959 + Subsystem_Fy_B.t967) * Subsystem_Fy_B.t73) +
    Subsystem_Fy_B.t133 * 2.0) - Subsystem_Fy_B.t137 * 2.0) +
    Subsystem_Fy_B.t154 * 2.0) + Subsystem_Fy_B.t152 * 2.0) -
    Subsystem_Fy_B.t160 * 2.0) + Subsystem_Fy_B.t164 * 2.0) -
    (Subsystem_Fy_B.t961 - Subsystem_Fy_B.t965) * Subsystem_Fy_B.t75) * 0.0025 /
    6.0;
  NominalState_predict[24] = x_nominal_prev[24] - ((((((((((((-0.0 -
    Subsystem_Fy_B.t171) + Subsystem_Fy_B.t132) + Subsystem_Fy_B.t11) +
    ((Subsystem_Fy_B.t835 + Subsystem_Fy_B.t836) - 1.0) * Subsystem_Fy_B.t75) -
    (Subsystem_Fy_B.t961 + Subsystem_Fy_B.t965) * Subsystem_Fy_B.t74) +
    Subsystem_Fy_B.t24 * 2.0) - Subsystem_Fy_B.t25 * 2.0) + Subsystem_Fy_B.t16 *
    2.0) + Subsystem_Fy_B.t15 * 2.0) - Subsystem_Fy_B.t10 * 2.0) +
    Subsystem_Fy_B.t14 * 2.0) - (Subsystem_Fy_B.t960 - Subsystem_Fy_B.t966) *
    Subsystem_Fy_B.t73) * 0.0025 / 6.0;
  NominalState_predict[25] = x_nominal_prev[25] - (((((((((Subsystem_Fy_B.t228 +
    Subsystem_Fy_B.t44) + Subsystem_Fy_B.t227) + Subsystem_Fy_B.t226) +
    Subsystem_Fy_B.t648) + Subsystem_Fy_B.t649) + Subsystem_Fy_B.t650) +
    Subsystem_Fy_B.t76 * Subsystem_Fy_B.t843 / 2.0) + Subsystem_Fy_B.t77 *
    Subsystem_Fy_B.t875 / 2.0) + Subsystem_Fy_B.t78 * Subsystem_Fy_B.t876 / 2.0)
    * 0.0025 / 6.0;
  NominalState_predict[26] = (((((((((Subsystem_Fy_B.t221 + Subsystem_Fy_B.t5) -
    Subsystem_Fy_B.t220) + Subsystem_Fy_B.t219) + Subsystem_Fy_B.t218) -
    Subsystem_Fy_B.t254) + Subsystem_Fy_B.t255) - Subsystem_Fy_B.t77 *
    Subsystem_Fy_B.t876 / 2.0) + Subsystem_Fy_B.t78 * Subsystem_Fy_B.t875 / 2.0)
    + Subsystem_Fy_B.t76 * Subsystem_Fy_B.t890 / 2.0) * 0.0025 / 6.0 +
    x_nominal_prev[26];
  NominalState_predict[27] = (((((((((Subsystem_Fy_B.t472 + Subsystem_Fy_B.t6) +
    Subsystem_Fy_B.t241) - Subsystem_Fy_B.t471) + Subsystem_Fy_B.t879) +
    Subsystem_Fy_B.t881) - Subsystem_Fy_B.t841) + Subsystem_Fy_B.t76 *
    Subsystem_Fy_B.t876 / 2.0) - Subsystem_Fy_B.t78 * Subsystem_Fy_B.t843 / 2.0)
    + Subsystem_Fy_B.t77 * Subsystem_Fy_B.t890 / 2.0) * 0.0025 / 6.0 +
    x_nominal_prev[27];
  NominalState_predict[28] = (((((((((Subsystem_Fy_B.t49 - Subsystem_Fy_B.t356)
    + Subsystem_Fy_B.t357) + Subsystem_Fy_B.t362) - Subsystem_Fy_B.t877) +
    Subsystem_Fy_B.t884) + Subsystem_Fy_B.t647) - Subsystem_Fy_B.t76 *
    Subsystem_Fy_B.t875 / 2.0) + Subsystem_Fy_B.t77 * Subsystem_Fy_B.t843 / 2.0)
    + Subsystem_Fy_B.t78 * Subsystem_Fy_B.t890 / 2.0) * 0.0025 / 6.0 +
    x_nominal_prev[28];
  std::memcpy(&NominalState_predict[29], &x_nominal_prev[29], 13U * sizeof
              (real_T));
  for (i_0 = 0; i_0 < 40; i_0++) {
    std::memset(&Subsystem_Fy_B.b_m[i_0 * 40], 0, 40U * sizeof(real_T));
    for (i_1 = 0; i_1 < 40; i_1++) {
      for (i = 0; i < 40; i++) {
        b_tmp = 40 * i_0 + i;
        Subsystem_Fy_B.b_m[b_tmp] += Subsystem_Fy_B.b[40 * i_1 + i] *
          errorstate_prev[40 * i_0 + i_1];
      }
    }
  }

  for (i_0 = 0; i_0 < 30; i_0++) {
    std::memset(&Subsystem_Fy_B.c_b[i_0 * 40], 0, 40U * sizeof(real_T));
    for (i_1 = 0; i_1 < 30; i_1++) {
      for (i = 0; i < 40; i++) {
        b_tmp = 40 * i_0 + i;
        Subsystem_Fy_B.c_b[b_tmp] += Subsystem_Fy_B.c_c[40 * i_1 + i] *
          processNoise[30 * i_0 + i_1];
      }
    }
  }

  for (i_0 = 0; i_0 < 40; i_0++) {
    std::memset(&Subsystem_Fy_B.b_c[i_0 * 40], 0, 40U * sizeof(real_T));
  }

  for (i_1 = 0; i_1 < 40; i_1++) {
    for (i_0 = 0; i_0 < 40; i_0++) {
      for (i = 0; i < 40; i++) {
        b_tmp = 40 * i_0 + i;
        Subsystem_Fy_B.b_c[b_tmp] += Subsystem_Fy_B.b_m[40 * i_1 + i] *
          Subsystem_Fy_B.b[40 * i_1 + i_0];
      }

      Subsystem_Fy_B.c[i_0 + 40 * i_1] = 0.0;
    }
  }

  for (i_0 = 0; i_0 < 30; i_0++) {
    for (i = 0; i < 40; i++) {
      for (i_1 = 0; i_1 < 40; i_1++) {
        b_tmp = 40 * i + i_1;
        Subsystem_Fy_B.c[b_tmp] += Subsystem_Fy_B.c_b[40 * i_0 + i_1] *
          Subsystem_Fy_B.c_c[40 * i_0 + i];
      }
    }
  }

  for (i_0 = 0; i_0 < 1600; i_0++) {
    ErrorState_predict[i_0] = Subsystem_Fy_B.b_c[i_0] + Subsystem_Fy_B.c[i_0];
  }
}

// Function for MATLAB Function: '<S1>/MATLAB Function'
void Subsystem_FyModelClass::Subsystem_Fy_quatrotate_edit(const real_T q[4],
  const real_T r[3], real_T qout[3])
{
  real_T dcm_tmp;
  real_T dcm_tmp_0;
  real_T dcm_tmp_1;
  real_T dcm_tmp_2;
  real_T dcm_tmp_3;
  real_T dcm_tmp_4;
  real_T dcm_tmp_5;
  int32_T i;
  Subsystem_Fy_B.qin_idx_3 = std::sqrt(((q[0] * q[0] + q[1] * q[1]) + q[2] * q[2])
    + q[3] * q[3]);
  Subsystem_Fy_B.qin_idx_0 = q[0] / Subsystem_Fy_B.qin_idx_3;
  Subsystem_Fy_B.qin_idx_1 = q[1] / Subsystem_Fy_B.qin_idx_3;
  Subsystem_Fy_B.qin_idx_2 = q[2] / Subsystem_Fy_B.qin_idx_3;
  Subsystem_Fy_B.qin_idx_3 = q[3] / Subsystem_Fy_B.qin_idx_3;
  dcm_tmp_0 = Subsystem_Fy_B.qin_idx_0 * Subsystem_Fy_B.qin_idx_0;
  dcm_tmp_1 = Subsystem_Fy_B.qin_idx_1 * Subsystem_Fy_B.qin_idx_1;
  dcm_tmp_2 = Subsystem_Fy_B.qin_idx_2 * Subsystem_Fy_B.qin_idx_2;
  dcm_tmp_3 = Subsystem_Fy_B.qin_idx_3 * Subsystem_Fy_B.qin_idx_3;
  Subsystem_Fy_B.dcm_g[0] = ((dcm_tmp_0 + dcm_tmp_1) - dcm_tmp_2) - dcm_tmp_3;
  Subsystem_Fy_B.dcm_tmp_o = Subsystem_Fy_B.qin_idx_1 * Subsystem_Fy_B.qin_idx_2;
  dcm_tmp = Subsystem_Fy_B.qin_idx_0 * Subsystem_Fy_B.qin_idx_3;
  Subsystem_Fy_B.dcm_g[3] = (Subsystem_Fy_B.dcm_tmp_o + dcm_tmp) * 2.0;
  dcm_tmp_4 = Subsystem_Fy_B.qin_idx_1 * Subsystem_Fy_B.qin_idx_3;
  dcm_tmp_5 = Subsystem_Fy_B.qin_idx_0 * Subsystem_Fy_B.qin_idx_2;
  Subsystem_Fy_B.dcm_g[6] = (dcm_tmp_4 - dcm_tmp_5) * 2.0;
  Subsystem_Fy_B.dcm_g[1] = (Subsystem_Fy_B.dcm_tmp_o - dcm_tmp) * 2.0;
  dcm_tmp_0 -= dcm_tmp_1;
  Subsystem_Fy_B.dcm_g[4] = (dcm_tmp_0 + dcm_tmp_2) - dcm_tmp_3;
  dcm_tmp_1 = Subsystem_Fy_B.qin_idx_2 * Subsystem_Fy_B.qin_idx_3;
  Subsystem_Fy_B.dcm_tmp_o = Subsystem_Fy_B.qin_idx_0 * Subsystem_Fy_B.qin_idx_1;
  Subsystem_Fy_B.dcm_g[7] = (dcm_tmp_1 + Subsystem_Fy_B.dcm_tmp_o) * 2.0;
  Subsystem_Fy_B.dcm_g[2] = (dcm_tmp_4 + dcm_tmp_5) * 2.0;
  Subsystem_Fy_B.dcm_g[5] = (dcm_tmp_1 - Subsystem_Fy_B.dcm_tmp_o) * 2.0;
  Subsystem_Fy_B.dcm_g[8] = (dcm_tmp_0 - dcm_tmp_2) + dcm_tmp_3;
  for (i = 0; i < 3; i++) {
    qout[i] = Subsystem_Fy_B.dcm_g[i + 6] * r[2] + (Subsystem_Fy_B.dcm_g[i + 3] *
      r[1] + Subsystem_Fy_B.dcm_g[i] * r[0]);
  }
}

// Model step function
void Subsystem_FyModelClass::step()
{
  boolean_T b[3];
  boolean_T tmp2;
  static const real_T b_0[3] = { -0.0422477548569116, -0.0229749805022196,
    0.00224970158881409 };

  static const real_T c[3] = { 0.000481025186908457, 0.00104955038446195,
    -0.00186800748296271 };

  static const real_T d[3] = { -1.5, -48.6, -19.9 };

  static const real_T e[3] = { -0.00213087026849489, -0.101001412610616,
    0.000779715520280068 };

  static const real_T f[3] = { 0.00319134325340909, -0.00277695549666309,
    0.0234631224105483 };

  static const real_T g[3] = { -0.6023, 0.1729, -0.0705 };

  static const real_T h[169] = { 0.17606932693795824, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17606932693795824, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17606932693795824, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.19356880436665769, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.19356880436665769, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.19356880436665769, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0001 };

  static const real_T i[900] = { 0.0014383814760000001, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0014383814760000001, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0014383814760000001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.000125910841, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000125910841, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000125910841, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    5.58802321E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 5.58802321E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.58802321E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.6124240249999999E-7, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.6124240249999999E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.6124240249999999E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00015452678315783477, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.00015452678315783477, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.00015452678315783477, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0014383814760000001, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0014383814760000001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0014383814760000001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000125910841, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.000125910841, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.000125910841, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 5.58802321E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.58802321E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.58802321E-6, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.6124240249999999E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.6124240249999999E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.6124240249999999E-7, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.00016988515364193914, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.00016988515364193914, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00016988515364193914 };

  boolean_T exitg1;

  // Outputs for Atomic SubSystem: '<Root>/Subsystem_Fy'
  // Outport: '<Root>/aavel' incorporates:
  //   Inport: '<Root>/simu_w'
  //   MATLAB Function: '<S1>/MATLAB Function'

  Subsystem_Fy_Y.aavel[0] = Subsystem_Fy_U.simu_w[0];

  // Outport: '<Root>/meas' incorporates:
  //   Inport: '<Root>/fimu_a'
  //   Inport: '<Root>/fimu_m'
  //   Inport: '<Root>/fimu_w'
  //   Inport: '<Root>/simu_a'
  //   Inport: '<Root>/simu_m'
  //   Inport: '<Root>/simu_w'
  //   MATLAB Function: '<S1>/MATLAB Function'

  Subsystem_Fy_Y.meas[0] = Subsystem_Fy_U.simu_a[0];
  Subsystem_Fy_Y.meas[3] = Subsystem_Fy_U.simu_w[0];
  Subsystem_Fy_Y.meas[6] = Subsystem_Fy_U.simu_m[0];
  Subsystem_Fy_Y.meas[9] = Subsystem_Fy_U.fimu_a[0];
  Subsystem_Fy_Y.meas[12] = Subsystem_Fy_U.fimu_w[0];
  Subsystem_Fy_Y.meas[15] = Subsystem_Fy_U.fimu_m[0];

  // Outport: '<Root>/aavel' incorporates:
  //   Inport: '<Root>/simu_w'
  //   MATLAB Function: '<S1>/MATLAB Function'

  Subsystem_Fy_Y.aavel[1] = Subsystem_Fy_U.simu_w[1];

  // Outport: '<Root>/meas' incorporates:
  //   Inport: '<Root>/fimu_a'
  //   Inport: '<Root>/fimu_m'
  //   Inport: '<Root>/fimu_w'
  //   Inport: '<Root>/simu_a'
  //   Inport: '<Root>/simu_m'
  //   Inport: '<Root>/simu_w'
  //   MATLAB Function: '<S1>/MATLAB Function'

  Subsystem_Fy_Y.meas[1] = Subsystem_Fy_U.simu_a[1];
  Subsystem_Fy_Y.meas[4] = Subsystem_Fy_U.simu_w[1];
  Subsystem_Fy_Y.meas[7] = Subsystem_Fy_U.simu_m[1];
  Subsystem_Fy_Y.meas[10] = Subsystem_Fy_U.fimu_a[1];
  Subsystem_Fy_Y.meas[13] = Subsystem_Fy_U.fimu_w[1];
  Subsystem_Fy_Y.meas[16] = Subsystem_Fy_U.fimu_m[1];

  // Outport: '<Root>/aavel' incorporates:
  //   Inport: '<Root>/simu_w'
  //   MATLAB Function: '<S1>/MATLAB Function'

  Subsystem_Fy_Y.aavel[2] = Subsystem_Fy_U.simu_w[2];

  // Outport: '<Root>/meas' incorporates:
  //   Inport: '<Root>/Fy'
  //   Inport: '<Root>/fimu_a'
  //   Inport: '<Root>/fimu_m'
  //   Inport: '<Root>/fimu_w'
  //   Inport: '<Root>/simu_a'
  //   Inport: '<Root>/simu_m'
  //   Inport: '<Root>/simu_w'
  //   MATLAB Function: '<S1>/MATLAB Function'

  Subsystem_Fy_Y.meas[2] = Subsystem_Fy_U.simu_a[2];
  Subsystem_Fy_Y.meas[5] = Subsystem_Fy_U.simu_w[2];
  Subsystem_Fy_Y.meas[8] = Subsystem_Fy_U.simu_m[2];
  Subsystem_Fy_Y.meas[11] = Subsystem_Fy_U.fimu_a[2];
  Subsystem_Fy_Y.meas[14] = Subsystem_Fy_U.fimu_w[2];
  Subsystem_Fy_Y.meas[17] = Subsystem_Fy_U.fimu_m[2];
  Subsystem_Fy_Y.meas[18] = Subsystem_Fy_U.Fy;

  // MATLAB Function: '<S1>/MATLAB Function' incorporates:
  //   Inport: '<Root>/fimu_a'
  //   Inport: '<Root>/fimu_m'
  //   Inport: '<Root>/simu_a'
  //   Inport: '<Root>/simu_m'

  if (!Subsystem_Fy_DW.NominalState_not_empty) {
    project_state_to_constrain_fcn(b_0, c, d, e, f, g, Subsystem_Fy_U.fimu_a,
      Subsystem_Fy_U.fimu_m, Subsystem_Fy_U.simu_a, Subsystem_Fy_U.simu_m,
      Subsystem_Fy_DW.NominalState);
    Subsystem_Fy_DW.NominalState_not_empty = true;
  }

  // Outport: '<Root>/states' incorporates:
  //   MATLAB Function: '<S1>/MATLAB Function'

  std::memcpy(&Subsystem_Fy_Y.states[0], &Subsystem_Fy_DW.NominalState[0], 42U *
              sizeof(real_T));

  // MATLAB Function: '<S1>/MATLAB Function' incorporates:
  //   Inport: '<Root>/Fy'
  //   Inport: '<Root>/Fy_tol'
  //   Inport: '<Root>/fimu_a'
  //   Inport: '<Root>/fimu_m'
  //   Inport: '<Root>/fimu_w'
  //   Inport: '<Root>/simu_a'
  //   Inport: '<Root>/simu_m'
  //   Outport: '<Root>/aavel'

  std::memcpy(&Subsystem_Fy_B.State[0], &Subsystem_Fy_DW.NominalState[0], 42U *
              sizeof(real_T));
  std::memcpy(&Subsystem_Fy_B.dv[0], &Subsystem_Fy_DW.errorCov_[0], 1600U *
              sizeof(real_T));
  std::memcpy(&Subsystem_Fy_B.dv4[0], &Subsystem_Fy_DW.NominalState[0], 42U *
              sizeof(real_T));
  Subsystem_F_correct_measurement(Subsystem_Fy_B.dv, h, Subsystem_Fy_U.fimu_m,
    Subsystem_Fy_U.simu_m, Subsystem_Fy_U.Fy > Subsystem_Fy_U.Fy_tol,
    Subsystem_Fy_B.dv4, Subsystem_Fy_DW.NominalState, Subsystem_Fy_DW.errorCov_);
  Subsystem_Fy_B.dv5[0] = Subsystem_Fy_U.fimu_a[0];
  Subsystem_Fy_B.dv5[3] = Subsystem_Fy_U.fimu_w[0];
  Subsystem_Fy_B.dv5[6] = Subsystem_Fy_U.simu_a[0];
  Subsystem_Fy_B.dv5[9] = Subsystem_Fy_Y.aavel[0];
  Subsystem_Fy_B.dv5[1] = Subsystem_Fy_U.fimu_a[1];
  Subsystem_Fy_B.dv5[4] = Subsystem_Fy_U.fimu_w[1];
  Subsystem_Fy_B.dv5[7] = Subsystem_Fy_U.simu_a[1];
  Subsystem_Fy_B.dv5[10] = Subsystem_Fy_Y.aavel[1];
  Subsystem_Fy_B.dv5[2] = Subsystem_Fy_U.fimu_a[2];
  Subsystem_Fy_B.dv5[5] = Subsystem_Fy_U.fimu_w[2];
  Subsystem_Fy_B.dv5[8] = Subsystem_Fy_U.simu_a[2];
  Subsystem_Fy_B.dv5[11] = Subsystem_Fy_Y.aavel[2];
  std::memcpy(&Subsystem_Fy_B.dv4[0], &Subsystem_Fy_DW.NominalState[0], 42U *
              sizeof(real_T));
  std::memcpy(&Subsystem_Fy_B.dv[0], &Subsystem_Fy_DW.errorCov_[0], 1600U *
              sizeof(real_T));
  Subsystem_Fy_predict_eksf(i, Subsystem_Fy_B.dv5, Subsystem_Fy_B.dv4,
    Subsystem_Fy_B.dv, Subsystem_Fy_DW.NominalState, Subsystem_Fy_DW.errorCov_);
  Subsystem_Fy_B.q12_idx_1 = ((Subsystem_Fy_B.State[6] * Subsystem_Fy_B.State[6]
    + Subsystem_Fy_B.State[7] * Subsystem_Fy_B.State[7]) + Subsystem_Fy_B.State
    [8] * Subsystem_Fy_B.State[8]) + Subsystem_Fy_B.State[9] *
    Subsystem_Fy_B.State[9];
  Subsystem_Fy_B.State_g[0] = Subsystem_Fy_B.State[6] / Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_B.State_g[1] = -Subsystem_Fy_B.State[7] /
    Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_B.State_g[2] = -Subsystem_Fy_B.State[8] /
    Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_B.State_g[3] = -Subsystem_Fy_B.State[9] /
    Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_quatrotate_edit(Subsystem_Fy_B.State_g, Subsystem_Fy_U.fimu_w,
    Subsystem_Fy_B.dv7);
  Subsystem_Fy_quatrotate_edit(&Subsystem_Fy_B.State[25], Subsystem_Fy_B.dv7,
    Subsystem_Fy_B.dv6);

  // Outport: '<Root>/aavel' incorporates:
  //   MATLAB Function: '<S1>/MATLAB Function'

  Subsystem_Fy_Y.aavel[0] = -Subsystem_Fy_Y.aavel[0] + Subsystem_Fy_B.dv6[0];
  Subsystem_Fy_Y.aavel[1] = -Subsystem_Fy_Y.aavel[1] + Subsystem_Fy_B.dv6[1];
  Subsystem_Fy_Y.aavel[2] = -Subsystem_Fy_Y.aavel[2] + Subsystem_Fy_B.dv6[2];

  // Outport: '<Root>/aavel_IE_DP_rad' incorporates:
  //   MATLAB Function: '<S1>/MATLAB Function'
  //   Outport: '<Root>/aavel'

  Subsystem_Fy_Y.aavel_IE_DP_rad[0] = Subsystem_Fy_Y.aavel[0];
  Subsystem_Fy_Y.aavel_IE_DP_rad[1] = Subsystem_Fy_Y.aavel[2];

  // Outport: '<Root>/aavel_IE_DP_deg' incorporates:
  //   MATLAB Function: '<S1>/MATLAB Function'
  //   Outport: '<Root>/aavel_IE_DP_rad'

  Subsystem_Fy_Y.aavel_IE_DP_deg[0] = Subsystem_Fy_Y.aavel_IE_DP_rad[0] * 180.0 /
    3.1415926535897931;
  Subsystem_Fy_Y.aavel_IE_DP_deg[1] = Subsystem_Fy_Y.aavel_IE_DP_rad[1] * 180.0 /
    3.1415926535897931;

  // MATLAB Function: '<S1>/MATLAB Function1' incorporates:
  //   Outport: '<Root>/states'

  Subsystem_Fy_B.q12_idx_1 = ((Subsystem_Fy_Y.states[25] *
    Subsystem_Fy_Y.states[25] + Subsystem_Fy_Y.states[26] *
    Subsystem_Fy_Y.states[26]) + Subsystem_Fy_Y.states[27] *
    Subsystem_Fy_Y.states[27]) + Subsystem_Fy_Y.states[28] *
    Subsystem_Fy_Y.states[28];
  Subsystem_Fy_B.axang_idx_0 = Subsystem_Fy_Y.states[25] /
    Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_B.axang_idx_1 = -Subsystem_Fy_Y.states[26] /
    Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_B.axang_idx_2 = -Subsystem_Fy_Y.states[27] /
    Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_B.axang_idx_3 = -Subsystem_Fy_Y.states[28] /
    Subsystem_Fy_B.q12_idx_1;
  Subsystem_Fy_B.q12_idx_0 = ((Subsystem_Fy_B.axang_idx_0 *
    Subsystem_Fy_Y.states[6] - Subsystem_Fy_B.axang_idx_1 *
    Subsystem_Fy_Y.states[7]) - Subsystem_Fy_B.axang_idx_2 *
    Subsystem_Fy_Y.states[8]) - Subsystem_Fy_B.axang_idx_3 *
    Subsystem_Fy_Y.states[9];
  Subsystem_Fy_B.q12_idx_1 = (Subsystem_Fy_B.axang_idx_0 *
    Subsystem_Fy_Y.states[7] + Subsystem_Fy_Y.states[6] *
    Subsystem_Fy_B.axang_idx_1) + (Subsystem_Fy_B.axang_idx_2 *
    Subsystem_Fy_Y.states[9] - Subsystem_Fy_Y.states[8] *
    Subsystem_Fy_B.axang_idx_3);
  Subsystem_Fy_B.q12_idx_2 = (Subsystem_Fy_B.axang_idx_0 *
    Subsystem_Fy_Y.states[8] + Subsystem_Fy_Y.states[6] *
    Subsystem_Fy_B.axang_idx_2) + (Subsystem_Fy_Y.states[7] *
    Subsystem_Fy_B.axang_idx_3 - Subsystem_Fy_B.axang_idx_1 *
    Subsystem_Fy_Y.states[9]);
  Subsystem_Fy_B.axang_idx_0 = (Subsystem_Fy_B.axang_idx_0 *
    Subsystem_Fy_Y.states[9] + Subsystem_Fy_Y.states[6] *
    Subsystem_Fy_B.axang_idx_3) + (Subsystem_Fy_B.axang_idx_1 *
    Subsystem_Fy_Y.states[8] - Subsystem_Fy_Y.states[7] *
    Subsystem_Fy_B.axang_idx_2);
  Subsystem_Fy_B.axang_idx_1 = std::sqrt(((Subsystem_Fy_B.q12_idx_0 *
    Subsystem_Fy_B.q12_idx_0 + Subsystem_Fy_B.q12_idx_1 *
    Subsystem_Fy_B.q12_idx_1) + Subsystem_Fy_B.q12_idx_2 *
    Subsystem_Fy_B.q12_idx_2) + Subsystem_Fy_B.axang_idx_0 *
    Subsystem_Fy_B.axang_idx_0);
  Subsystem_Fy_B.axang_idx_2 = Subsystem_Fy_B.q12_idx_0 /
    Subsystem_Fy_B.axang_idx_1;
  Subsystem_Fy_B.q12_idx_0 = 2.0 * std::acos(Subsystem_Fy_B.axang_idx_2);
  if (Subsystem_Fy_B.q12_idx_0 > 3.1415926535897931) {
    Subsystem_Fy_B.q12_idx_0 -= 6.2831853071795862;
  }

  Subsystem_Fy_B.axang_idx_2 = std::sqrt(1.0 - Subsystem_Fy_B.axang_idx_2 *
    Subsystem_Fy_B.axang_idx_2);
  Subsystem_Fy_Y.aangleh[0] = Subsystem_Fy_B.q12_idx_1 /
    Subsystem_Fy_B.axang_idx_1 / Subsystem_Fy_B.axang_idx_2 *
    Subsystem_Fy_B.q12_idx_0;
  b[0] = rtIsNaN(Subsystem_Fy_Y.aangleh[0]);
  Subsystem_Fy_Y.aangleh[1] = Subsystem_Fy_B.q12_idx_2 /
    Subsystem_Fy_B.axang_idx_1 / Subsystem_Fy_B.axang_idx_2 *
    Subsystem_Fy_B.q12_idx_0;
  b[1] = rtIsNaN(Subsystem_Fy_Y.aangleh[1]);
  Subsystem_Fy_Y.aangleh[2] = Subsystem_Fy_B.axang_idx_0 /
    Subsystem_Fy_B.axang_idx_1 / Subsystem_Fy_B.axang_idx_2 *
    Subsystem_Fy_B.q12_idx_0;
  b[2] = rtIsNaN(Subsystem_Fy_Y.aangleh[2]);
  tmp2 = true;
  Subsystem_Fy_B.i = 0;
  exitg1 = false;
  while ((!exitg1) && (Subsystem_Fy_B.i < 3)) {
    if (!b[Subsystem_Fy_B.i]) {
      tmp2 = false;
      exitg1 = true;
    } else {
      Subsystem_Fy_B.i++;
    }
  }

  if (tmp2) {
    Subsystem_Fy_Y.aangleh[0] = 0.0;
    Subsystem_Fy_Y.aangleh[1] = 0.0;
    Subsystem_Fy_Y.aangleh[2] = 0.0;
  }

  // Outport: '<Root>/IE_DP_rad' incorporates:
  //   MATLAB Function: '<S1>/MATLAB Function1'

  Subsystem_Fy_Y.IE_DP_rad[0] = Subsystem_Fy_Y.aangleh[0];

  // MATLAB Function: '<S1>/MATLAB Function1'
  Subsystem_Fy_Y.aangleh[0] = Subsystem_Fy_Y.aangleh[0] * 180.0 /
    3.1415926535897931;
  Subsystem_Fy_Y.aangleh[1] = Subsystem_Fy_Y.aangleh[1] * 180.0 /
    3.1415926535897931;

  // Outport: '<Root>/IE_DP_rad' incorporates:
  //   MATLAB Function: '<S1>/MATLAB Function1'

  Subsystem_Fy_Y.IE_DP_rad[1] = Subsystem_Fy_Y.aangleh[2];

  // MATLAB Function: '<S1>/MATLAB Function1'
  Subsystem_Fy_Y.aangleh[2] = Subsystem_Fy_Y.aangleh[2] * 180.0 /
    3.1415926535897931;

  // Outport: '<Root>/IE_DP_deg' incorporates:
  //   MATLAB Function: '<S1>/MATLAB Function1'

  Subsystem_Fy_Y.IE_DP_deg[0] = Subsystem_Fy_Y.aangleh[0];
  Subsystem_Fy_Y.IE_DP_deg[1] = Subsystem_Fy_Y.aangleh[2];

  // End of Outputs for SubSystem: '<Root>/Subsystem_Fy'
}

// Model initialize function
void Subsystem_FyModelClass::initialize()
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  {
    int32_T i;
    static const real_T b[40] = { 0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0,
      5.0292208889999989E-5, 5.0292208889999989E-5, 5.0292208889999989E-5,
      1.4511816224999995E-6, 1.4511816224999995E-6, 1.4511816224999995E-6,
      0.0013907410484205129, 0.0013907410484205129, 0.0013907410484205129, 0.25,
      0.25, 0.25, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 5.0292208889999989E-5,
      5.0292208889999989E-5, 5.0292208889999989E-5, 1.4511816224999995E-6,
      1.4511816224999995E-6, 1.4511816224999995E-6, 0.0013907410484205129,
      0.0013907410484205129, 0.0013907410484205129, 0.0, 0.0, 0.0, 0.0 };

    // SystemInitialize for Atomic SubSystem: '<Root>/Subsystem_Fy'
    // SystemInitialize for MATLAB Function: '<S1>/MATLAB Function'
    std::memset(&Subsystem_Fy_DW.errorCov_[0], 0, 1600U * sizeof(real_T));
    for (i = 0; i < 40; i++) {
      Subsystem_Fy_DW.errorCov_[i + 40 * i] = b[i];
    }

    // End of SystemInitialize for MATLAB Function: '<S1>/MATLAB Function'
    // End of SystemInitialize for SubSystem: '<Root>/Subsystem_Fy'
  }
}

// Model terminate function
void Subsystem_FyModelClass::terminate()
{
  // (no terminate code required)
}

// Constructor
Subsystem_FyModelClass::Subsystem_FyModelClass() :
  Subsystem_Fy_B(),
  Subsystem_Fy_DW(),
  Subsystem_Fy_U(),
  Subsystem_Fy_Y(),
  Subsystem_Fy_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
Subsystem_FyModelClass::~Subsystem_FyModelClass()
{
  // Currently there is no destructor body generated.
}

// Real-Time Model get method
Subsystem_FyModelClass::RT_MODEL_Subsystem_Fy_T * Subsystem_FyModelClass::getRTM
  ()
{
  return (&Subsystem_Fy_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
