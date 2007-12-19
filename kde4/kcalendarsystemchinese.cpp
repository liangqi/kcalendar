/*
    Copyright (c) 2007 Liang Qi <cavendish.qi@gmail.com>
        Calendar conversion routines based on:
        1. ccal v2.4 by Zhuo Meng <zhuo@thunder.cwru.edu>
        2. NOVAS-C v2.0 (1 Nov 98) by
               U. S. Naval Observatory
               Astronomical Applications Dept.
               3450 Massachusetts Ave., NW
               Washington, DC  20392-5420
        3. Lunar Outreach Services by Christopher Osburn(1996)
 
    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.
 
    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.
 
    You should have received a copy of the GNU Library General Public License
    along with this library; see the file COPYING.LIB.  If not, write to
    the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
    Boston, MA 02110-1301, USA.
*/

// Derived chinese kde calendar class

#include "kcalendarsystemchinese.h"

#include "kdebug.h"
#include "klocale.h"

#include <QtCore/QDate>
#include <QtCore/QCharRef>
#include <QtCore/QHash>

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <vector>

using namespace std;

typedef std::vector<double> vdouble;

//===========================================================================
// novas.h
//
// NOVAS-C Version 2.0 (1 Nov 98)
// Header file for novas.c
//===========================================================================

/*
   Define "origin" constants.
*/

   #define BARYC  0
   #define HELIOC 1

/*
   Function prototypes
*/

static    void earthtilt (double tjd,
                   double *mobl, double *tobl, double *eqeq,
                   double *psi, double *eps);

static    short int aberration (double *pos, double *vel, double lighttime,
                         double *pos2);

static    void precession (double tjd1, double *pos, double tjd2,
                    double *pos2);

static    short int nutate (double tjd, short int fn1, double *pos,
                     double *pos2);

static    short int nutation_angles (double tdbtime,
                              double *longnutation,
                              double *obliqnutation);

static    void fund_args (double t,
                   double a[5]);

static    void tdb2tdt (double tdb,
                 double *tdtjd, double *secdiff);

static    void radec2vector (double ra, double dec, double dist,
                      double *vector);

static    short int solarsystem (double tjd, short int body, short int origin, 
                          double *pos, double *vel);

static    double julian_date (short int year, short int month, short int day,
                       double hour);

/*
static    void cal_date (double tjd,
                  short int *year, short int *month, short int *day,
                  double *hour);
*/

//===========================================================================
//  End of novas.h
//===========================================================================

//===========================================================================
// Copyright 1996, Christopher Osburn, Lunar Outreach Services,
//
// moon.h
//
// Headers and function defs for moon phase routines
//===========================================================================

static double moonphasebylunation(int lun, int phi);  /* moonphase.c */
static double moonphase(double k, int phi);  /* moonphase.c */
static double torad(double x);              /* misc.c      */

//===========================================================================
//  End of moon.h
//===========================================================================

//===========================================================================
// Copyright 1996, Christopher Osburn, Lunar Outreach Services,
//
// misc.cpp
// 1996/02/11
//
// Miscellaneous routines for moon phase programs
//===========================================================================

static double torad(double x)
/* convert x to radians */
{
  x = fmod(x, 360.0); /* normalize the angle */
  return ((x) * 0.01745329251994329576);
                    /* and return the result */
}

//===========================================================================
//  End of misc.cpp
//===========================================================================

//===========================================================================
// Copyright 1996, Christopher Osburn, Lunar Outreach Services,
//
// moonphase.cpp
// 1996/02/11
//
// calculate phase of the moon per Meeus Ch. 47
//
// Parameters:
//    int lun:  phase parameter.  This is the number of lunations
//              since the New Moon of 2000 January 6.
//
//    int phi:  another phase parameter, selecting the phase of the
//              moon.  0 = New, 1 = First Qtr, 2 = Full, 3 = Last Qtr
//
// Return:  Apparent JD of the needed phase
//===========================================================================

static double moonphase(double k, int phi)
{
  int i;                       /* iterator to be named later.  Every
                                  program needs an i */
  double T;                    /* time parameter, Julian Centuries since
                                  J2000 */
  double JDE;                  /* Julian Ephemeris Day of phase event */
  double E;                    /* Eccentricity anomaly */
  double M;                    /* Sun's mean anomaly */
  double M1;                   /* Moon's mean anomaly */
  double F;                    /* Moon's argument of latitude */
  double O;                    /* Moon's longitude of ascenfing node */
  double A[15];                /* planetary arguments */
  double W;                    /* added correction for quarter phases */

  T = k / 1236.85;                            /* (47.3) */

  /* this is the first approximation.  all else is for style points! */
  JDE = 2451550.09765 + (29.530588853 * k)    /* (47.1) */
        + T * T * (0.0001337 + T * (-0.000000150 + 0.00000000073 * T));

  /* these are correction parameters used below */
  E = 1.0                                     /* (45.6) */
      + T * (-0.002516 + -0.0000074 * T);
  M = 2.5534 + 29.10535669 * k                /* (47.4) */
      + T * T * (-0.0000218 + -0.00000011 * T);
  M1 = 201.5643 + 385.81693528 * k            /* (47.5) */
       + T * T * (0.0107438 + T * (0.00001239 + -0.000000058 * T));
  F = 160.7108 + 390.67050274 * k             /* (47.6) */
      + T * T * (-0.0016341 * T * (-0.00000227 + 0.000000011 * T));
  O = 124.7746 - 1.56375580 * k               /* (47.7) */
      + T * T * (0.0020691 + 0.00000215 * T);

  /* planetary arguments */
  A[0]  = 0; /* unused! */
  A[1]  = 299.77 +  0.107408 * k - 0.009173 * T * T;
  A[2]  = 251.88 +  0.016321 * k;
  A[3]  = 251.83 + 26.651886 * k;
  A[4]  = 349.42 + 36.412478 * k;
  A[5]  =  84.66 + 18.206239 * k;
  A[6]  = 141.74 + 53.303771 * k;
  A[7]  = 207.14 +  2.453732 * k;
  A[8]  = 154.84 +  7.306860 * k;
  A[9]  =  34.52 + 27.261239 * k;
  A[10] = 207.19 +  0.121824 * k;
  A[11] = 291.34 +  1.844379 * k;
  A[12] = 161.72 + 24.198154 * k;
  A[13] = 239.56 + 25.513099 * k;
  A[14] = 331.55 +  3.592518 * k;

  /* all of the above crap must be made into radians!!! */
  /* except for E... */

  M = torad(M);
  M1 = torad(M1);
  F = torad(F);
  O = torad(O);

  /* all those planetary arguments, too! */
  for (i=1; i<=14; i++)
    A[i] = torad(A[i]);

  /* ok, we have all the parameters, let's apply them to the JDE.
    (remember the JDE?  this is a program about the JDE...)        */
  
  switch(phi)
  {
    /* a special case for each different phase.  NOTE!, 
       I'm not treating these in a 0123 order!!!  Pay
       attention, there,  you!                         */

    case 0: /* New Moon */
      JDE = JDE
          - 0.40720         * sin (M1)
          + 0.17241 * E     * sin (M)
          + 0.01608         * sin (2.0 * M1)
          + 0.01039         * sin (2.0 * F)
          + 0.00739 * E     * sin (M1 - M)
          - 0.00514 * E     * sin (M1 + M)
          + 0.00208 * E * E * sin (2.0 * M)
          - 0.00111         * sin (M1 - 2.0 * F)
          - 0.00057         * sin (M1 + 2.0 * F)
          + 0.00056 * E     * sin (2.0 * M1 + M)
          - 0.00042         * sin (3.0 * M1)
          + 0.00042 * E     * sin (M + 2.0 * F)
          + 0.00038 * E     * sin (M - 2.0 * F)
          - 0.00024 * E     * sin (2.0 * M1 - M)
          - 0.00017         * sin (O)
          - 0.00007         * sin (M1 + 2.0 * M)
          + 0.00004         * sin (2.0 * M1 - 2.0 * F)
          + 0.00004         * sin (3.0 * M)
          + 0.00003         * sin (M1 + M - 2.0 * F)
          + 0.00003         * sin (2.0 * M1 + 2.0 * F)
          - 0.00003         * sin (M1 + M + 2.0 * F)
          + 0.00003         * sin (M1 - M + 2.0 * F)
          - 0.00002         * sin (M1 - M - 2.0 * F)
          - 0.00002         * sin (3.0 * M1 + M)
          + 0.00002         * sin (4.0 * M1);

          break;

    case 2: /* Full Moon */
      JDE = JDE
          - 0.40614         * sin (M1)
          + 0.17302 * E     * sin (M)
          + 0.01614         * sin (2.0 * M1)
          + 0.01043         * sin (2.0 * F)
          + 0.00734 * E     * sin (M1 - M)
          - 0.00515 * E     * sin (M1 + M)
          + 0.00209 * E * E * sin (2.0 * M)
          - 0.00111         * sin (M1 - 2.0 * F)
          - 0.00057         * sin (M1 + 2.0 * F)
          + 0.00056 * E     * sin (2.0 * M1 + M)
          - 0.00042         * sin (3.0 * M1)
          + 0.00042 * E     * sin (M + 2.0 * F)
          + 0.00038 * E     * sin (M - 2.0 * F)
          - 0.00024 * E     * sin (2.0 * M1 - M)
          - 0.00017         * sin (O)
          - 0.00007         * sin (M1 + 2.0 * M)
          + 0.00004         * sin (2.0 * M1 - 2.0 * F)
          + 0.00004         * sin (3.0 * M)
          + 0.00003         * sin (M1 + M - 2.0 * F)
          + 0.00003         * sin (2.0 * M1 + 2.0 * F)
          - 0.00003         * sin (M1 + M + 2.0 * F)
          + 0.00003         * sin (M1 - M + 2.0 * F)
          - 0.00002         * sin (M1 - M - 2.0 * F)
          - 0.00002         * sin (3.0 * M1 + M)
          + 0.00002         * sin (4.0 * M1);

          break;

    case 1: /* First Quarter */
    case 3: /* Last Quarter */
      JDE = JDE
          - 0.62801         * sin (M1)
          + 0.17172 * E     * sin (M)
          - 0.01183 * E     * sin (M1 + M)
          + 0.00862         * sin (2.0 * M1)
          + 0.00804         * sin (2.0 * F)
          + 0.00454 * E     * sin (M1 - M)
          + 0.00204 * E * E * sin (2.0 * M)
          - 0.00180         * sin (M1 - 2.0 * F)
          - 0.00070         * sin (M1 + 2.0 * F)
          - 0.00040         * sin (3.0 * M1)
          - 0.00034 * E     * sin (2.0 * M1 - M)
          + 0.00032 * E     * sin (M + 2.0 * F)
          + 0.00032 * E     * sin (M - 2.0 * F)
          - 0.00028 * E * E * sin (M1 + 2.0 * M)
          + 0.00027 * E     * sin (2.0 * M1 + M)
          - 0.00017         * sin (O)
          - 0.00005         * sin (M1 - M - 2.0 * F)
          + 0.00004         * sin (2.0 * M1 + 2.0 * F)
          - 0.00004         * sin (M1 + M + 2.0 * F)
          + 0.00004         * sin (M1 - 2.0 * M)
          + 0.00003         * sin (M1 + M - 2.0 * F)
          + 0.00003         * sin (3.0 * M)
          + 0.00002         * sin (2.0 * M1 - 2.0 * F)
          + 0.00002         * sin (M1 - M + 2.0 * F)
          - 0.00002         * sin (3.0 * M1 + M);
    
      W = 0.00306
        - 0.00038 * E * cos(M)
        + 0.00026 * cos(M1)
        - 0.00002 * cos(M1 - M)
        + 0.00002 * cos(M1 + M)
        + 0.00002 * cos(2.0 * F);
      if (phi == 3)
        W = -W;
      JDE += W;

      break;

    default: /* oops! */
      fprintf(stderr, "The Moon has exploded!\n");
      exit(1);
      break; /* unexecuted code */
  }
      /* now there are some final correction to everything */
  JDE = JDE
      + 0.000325 * sin(A[1])
      + 0.000165 * sin(A[2])
      + 0.000164 * sin(A[3])
      + 0.000126 * sin(A[4])
      + 0.000110 * sin(A[5])
      + 0.000062 * sin(A[6])
      + 0.000060 * sin(A[7])
      + 0.000056 * sin(A[8])
      + 0.000047 * sin(A[9])
      + 0.000042 * sin(A[10])
      + 0.000040 * sin(A[11])
      + 0.000037 * sin(A[12])
      + 0.000035 * sin(A[13])
      + 0.000023 * sin(A[14]);

  return JDE;
}

static double moonphasebylunation(int lun, int phi)
{
  double k;

  k = lun + phi / 4.0;
  return moonphase(k, phi);
}

//===========================================================================
//  End of moonphase.cpp
//===========================================================================

//===========================================================================
// Adapted from novascon.h from the NOVAS-C package
// The whole package can be obtained from
// http://aa.usno.navy.mil/AA/software/novas/novas_c/novasc_info.html
//
// Naval Observatory Vector Astrometry Subroutines
// C Version
//
// U. S. Naval Observatory
// Astronomical Applications Dept.
// 3450 Massachusetts Ave., NW
// Washington, DC  20392-5420
//===========================================================================

//===========================================================================
// novascon.h
//
// NOVAS-C Version 2.0 (1 Nov 98)
// Header file for novascon.c
//===========================================================================

//const short int FN0;

/*
   TDB Julian date of epoch J2000.0.
*/

//const double T0;

/*
   Speed of light in AU/Day.
*/

//const double C;

/*
   Value of 2.0 * pi in radians.
*/

//const double TWOPI;

/*
   Angle conversion constants.
*/

//const double RAD2SEC;
//const double DEG2RAD;
//const double RAD2DEG;

//===========================================================================
//  End of novascon.h
//===========================================================================

//===========================================================================
// novascon.cpp
//
// NOVAS-C Version 2.0 (1 Nov 98)
// Constants file
//===========================================================================

const short int FN0 = 0;

/*
   TDB Julian date of epoch J2000.0.
*/

const double T0 = 2451545.00000000;

/*
   Speed of light in AU/Day.
*/

const double C = 173.14463348;

/*
   Value of pi in radians.
*/

const double TWOPI = 6.28318530717958647692;

/*
   Angle conversion constants.
*/

const double RAD2SEC = 206264.806247096355;
const double DEG2RAD = 0.017453292519943296;
const double RAD2DEG = 57.295779513082321;

//===========================================================================
//  End of novascon.cpp
//===========================================================================

//===========================================================================
// solarsystem.h
//
// NOVAS-C Version 2.0 (1 Nov 98)
// Header file for all source files containing versions of
// NOVAS-C function 'solarsystem'
//===========================================================================

/*
   Function prototypes
*/

static short int solarsystem (double tjd, short int body, short int origin,
                          double *pos, double *vel);

//===========================================================================
//  End of solarsystem.h
//===========================================================================

//===========================================================================
// solsys3.cpp
//
// NOVAS-C Version 2.0 (1 Nov 98)
// Solar System function; version 3.
//===========================================================================

/*
   Additional function prototype.
*/

static void sun_eph (double jd,
              double *ra, double *dec, double *dis);

/********sun_eph */

static void sun_eph (double jd,

              double *ra, double *dec, double *dis)
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute equatorial spherical coordinates of Sun referred to
      the mean equator and equinox of date.

   REFERENCES:
      Bretagnon, P. and Simon, J.L. (1986).  Planetary Programs and
         Tables from -4000 to + 2800. (Richmond, VA: Willmann-Bell).

   INPUT
   ARGUMENTS:
      jd (double)
         Julian date on TDT or ET time scale.

   OUTPUT
   ARGUMENTS:
      ra (double)
         Right ascension referred to mean equator and equinox of date
         (hours).
      dec (double)
         Declination referred to mean equator and equinox of date 
         (degrees).
      dis (double)
         Geocentric distance (AU).

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0
      TWOPI
      RAD2DEG

   FUNCTIONS
   CALLED:
      sin           math.h
      cos           math.h
      asin          math.h
      atan2         math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-94/JAB (USNO/AA)
      V1.1/05-96/JAB (USNO/AA): Compute mean coordinates instead of
                                apparent.

   NOTES:
      1. Quoted accuracy is 2.0 + 0.03 * T^2 arcsec, where T is
      measured in units of 1000 years from J2000.0.  See reference.

------------------------------------------------------------------------
*/
{
   short int i;

   double sum_lon = 0.0;
   double sum_r = 0.0;
   const double factor = 1.0e-07;
   double u, arg, lon, lat, t, t2, emean, sin_lon;

   struct sun_con
   {
   double l;
   double r;
   double alpha;
   double nu;
   };

   static const struct sun_con con[50] =
      {{403406.0,      0.0, 4.721964,     1.621043},
       {195207.0, -97597.0, 5.937458, 62830.348067}, 
       {119433.0, -59715.0, 1.115589, 62830.821524}, 
       {112392.0, -56188.0, 5.781616, 62829.634302}, 
       {  3891.0,  -1556.0, 5.5474  , 125660.5691 }, 
       {  2819.0,  -1126.0, 1.5120  , 125660.9845 }, 
       {  1721.0,   -861.0, 4.1897  ,  62832.4766 }, 
       {     0.0,    941.0, 1.163   ,      0.813  }, 
       {   660.0,   -264.0, 5.415   , 125659.310  }, 
       {   350.0,   -163.0, 4.315   ,  57533.850  }, 
       {   334.0,      0.0, 4.553   ,    -33.931  }, 
       {   314.0,    309.0, 5.198   , 777137.715  }, 
       {   268.0,   -158.0, 5.989   ,  78604.191  }, 
       {   242.0,      0.0, 2.911   ,      5.412  }, 
       {   234.0,    -54.0, 1.423   ,  39302.098  }, 
       {   158.0,      0.0, 0.061   ,    -34.861  }, 
       {   132.0,    -93.0, 2.317   , 115067.698  }, 
       {   129.0,    -20.0, 3.193   ,  15774.337  }, 
       {   114.0,      0.0, 2.828   ,   5296.670  }, 
       {    99.0,    -47.0, 0.52    ,  58849.27   }, 
       {    93.0,      0.0, 4.65    ,   5296.11   }, 
       {    86.0,      0.0, 4.35    ,  -3980.70   }, 
       {    78.0,    -33.0, 2.75    ,  52237.69   }, 
       {    72.0,    -32.0, 4.50    ,  55076.47   }, 
       {    68.0,      0.0, 3.23    ,    261.08   }, 
       {    64.0,    -10.0, 1.22    ,  15773.85   }, 
       {    46.0,    -16.0, 0.14    ,  188491.03  }, 
       {    38.0,      0.0, 3.44    ,   -7756.55  }, 
       {    37.0,      0.0, 4.37    ,     264.89  }, 
       {    32.0,    -24.0, 1.14    ,  117906.27  }, 
       {    29.0,    -13.0, 2.84    ,   55075.75  }, 
       {    28.0,      0.0, 5.96    ,   -7961.39  }, 
       {    27.0,     -9.0, 5.09    ,  188489.81  }, 
       {    27.0,      0.0, 1.72    ,    2132.19  }, 
       {    25.0,    -17.0, 2.56    ,  109771.03  }, 
       {    24.0,    -11.0, 1.92    ,   54868.56  }, 
       {    21.0,      0.0, 0.09    ,   25443.93  }, 
       {    21.0,     31.0, 5.98    ,  -55731.43  }, 
       {    20.0,    -10.0, 4.03    ,   60697.74  }, 
       {    18.0,      0.0, 4.27    ,    2132.79  }, 
       {    17.0,    -12.0, 0.79    ,  109771.63  }, 
       {    14.0,      0.0, 4.24    ,   -7752.82  }, 
       {    13.0,     -5.0, 2.01    ,  188491.91  }, 
       {    13.0,      0.0, 2.65    ,     207.81  }, 
       {    13.0,      0.0, 4.98    ,   29424.63  }, 
       {    12.0,      0.0, 0.93    ,      -7.99  }, 
       {    10.0,      0.0, 2.21    ,   46941.14  }, 
       {    10.0,      0.0, 3.59    ,     -68.29  }, 
       {    10.0,      0.0, 1.50    ,   21463.25  }, 
       {    10.0,     -9.0, 2.55    ,  157208.40  }};

/*
   Define the time unit 'u', measured in units of 10000 Julian years
   from J2000.0.
*/

   u = (jd - T0) / 3652500.0;
   
/*
   Compute longitude and distance terms from the series.
*/

   for (i = 0; i < 50; i++)
   {
      arg = con[i].alpha + con[i].nu * u;
      sum_lon += con[i].l * sin (arg);
      sum_r += con[i].r * cos (arg);
   }

/*
   Compute longitude, latitude, and distance referred to mean equinox
   and ecliptic of date.
*/

   lon = 4.9353929 + 62833.1961680 * u + factor * sum_lon;

   lon = fmod (lon, TWOPI);
   if (lon < 0.0)
      lon += TWOPI;

   lat = 0.0;

   *dis = 1.0001026 + factor * sum_r;

/*
   Compute mean obliquity of the ecliptic.
*/

   t = u * 100.0;
   t2 = t * t;
   emean = (0.001813 * t2 * t - 0.00059 * t2 - 46.8150 * t +
      84381.448) / RAD2SEC;

/*
   Compute equatorial spherical coordinates referred to the mean equator 
   and equinox of date.
*/

   sin_lon = sin (lon);
   *ra = atan2 ((cos (emean) * sin_lon), cos (lon)) * RAD2DEG;
   *ra = fmod (*ra, 360.0);
   if (*ra < 0.0)
      *ra += 360.0;
   *ra = *ra / 15.0;

   *dec = asin (sin (emean) * sin_lon) * RAD2DEG;
   
   return;
}

/********solarsystem */

static short int solarsystem (double tjd, short int body, short int origin,
                       double *pos, double *vel)

/*
------------------------------------------------------------------------

   PURPOSE:    
      Provides the position and velocity of the Earth at epoch 'tjd'
      by evaluating a closed-form theory without reference to an 
      external file.  This function can also provide the position
      and velocity of the Sun.

   REFERENCES: 
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.
      Explanatory Supplement to The Astronomical Almanac (1992).

   INPUT
   ARGUMENTS:
      tjd (double)
         TDB Julian date.
      body (short int)
         Body identification number.
         Set 'body' = 0 or 'body' = 1 or 'body' = 10 for the Sun.
         Set 'body' = 2 or 'body' = 3 for the Earth.
      origin (short int)
         Origin code; solar system barycenter   = 0,
                      center of mass of the Sun = 1.

   OUTPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector of 'body' at 'tjd'; equatorial rectangular
         coordinates in AU referred to the mean equator and equinox
         of J2000.0.
      vel[3] (double)
         Velocity vector of 'body' at 'tjd'; equatorial rectangular
         system referred to the mean equator and equinox of J2000.0,
         in AU/Day.

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.
         1...Input Julian date ('tjd') out of range.
         2...Invalid value of 'body'.

   GLOBALS
   USED:
      T0, TWOPI.

   FUNCTIONS
   CALLED:
      sun_eph          solsys3.c
      radec2vector     novas.c
      precession       novas.c
      sin              math.h
      cos              math.h
      fabs             math.h
      fmod             math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/05-96/JAB (USNO/AA) Convert to C; substitute new theory of
                               Sun.
      V1.1/06-98/JAB (USNO/AA) Updated planetary masses & mean elements.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'solsys' version 3.

------------------------------------------------------------------------
*/
{
   short int ierr = 0;
   short int i;

/*
   The arrays below contain data for the four largest planets.  Masses
   are DE405 values; elements are from Explanatory Supplement, p. 316). 
   These data are used for barycenter computations only.
*/

   const double pm[4] = {1047.349, 3497.898, 22903.0, 19412.2};
   const double pa[4] = {5.203363, 9.537070, 19.191264, 30.068963};
   const double pl[4] = {0.600470, 0.871693, 5.466933, 5.321160};
   const double pn[4] = {1.450138e-3, 5.841727e-4, 2.047497e-4, 
                         1.043891e-4};

/*
   'obl' is the obliquity of ecliptic at epoch J2000.0 in degrees.
*/

   const double obl = 23.43929111;

   static double tlast = 0.0;
   static double sine, cose, tmass, pbary[3], vbary[3];

   double oblr, qjd, ras, decs, diss, pos1[3], p[3][3], dlon, sinl,
      cosl, x, y, z, xdot, ydot, zdot, f;

/*
   Initialize constants.
*/

   if (tlast == 0.0)
   {
      oblr = obl * TWOPI / 360.0;
      sine = sin (oblr);
      cose = cos (oblr);
      tmass = 1.0;
      for (i = 0; i < 4; i++)
         tmass += 1.0 / pm[i];
      tlast = 1.0;
   }

/*
   Check if input Julian date is within range.

   if ((tjd < 2340000.5) || (tjd > 2560000.5))
      return (ierr = 1);
*/

/*
   Form helicentric coordinates of the Sun or Earth, depending on
   'body'.
*/

   if ((body == 0) || (body == 1) || (body == 10))
      for (i = 0; i < 3; i++)
         pos[i] = vel[i] = 0.0;

    else if ((body == 2) || (body == 3))
    {
      for (i = 0; i < 3; i++)
      {
         qjd = tjd + (double) (i - 1) * 0.1;
         sun_eph (qjd, &ras,&decs,&diss);
         radec2vector (ras,decs,diss, pos1);
         precession (qjd,pos1,T0, pos);
         p[i][0] = -pos[0];
         p[i][1] = -pos[1];
         p[i][2] = -pos[2];
      }
      for (i = 0; i < 3; i++)
      {
         pos[i] = p[1][i];
         vel[i] = (p[2][i] - p[0][i]) / 0.2;
      }
    }

    else
      return (ierr = 2);
           
/*
   If 'origin' = 0, move origin to solar system barycenter.

   Solar system barycenter coordinates are computed from rough
   approximations of the coordinates of the four largest planets.
*/

   if (origin == 0)
   {
      if (fabs (tjd - tlast) >= 1.0e-06)
      {
         for (i = 0; i < 3; i++)
            pbary[i] = vbary[i] = 0.0;

/*
   The following loop cycles once for each of the four planets.

   'sinl' and 'cosl' are the sine and cosine of the planet's mean
   longitude.
*/

         for (i = 0; i < 4; i++)
         {
            dlon = pl[i] + pn[i] * (tjd - T0);
            dlon = fmod (dlon, TWOPI);
            sinl = sin (dlon);
            cosl = cos (dlon);

            x =  pa[i] * cosl;
            y =  pa[i] * sinl * cose;
            z =  pa[i] * sinl * sine;
            xdot = -pa[i] * pn[i] * sinl;
            ydot =  pa[i] * pn[i] * cosl * cose;
            zdot =  pa[i] * pn[i] * cosl * sine;

            f = 1.0 / (pm[i] * tmass);

            pbary[0] += x * f;
            pbary[1] += y * f;
            pbary[2] += z * f;
            vbary[0] += xdot * f;
            vbary[1] += ydot * f;
            vbary[2] += zdot * f;
         }

         tlast = tjd;
      }

      for (i = 0; i < 3; i++)
      {
         pos[i] -= pbary[i];
         vel[i] -= vbary[i];
      }
   }

   return (ierr);
}

//===========================================================================
//  End of solsys3.cpp
//===========================================================================

//===========================================================================
// novas.cpp
//
// NOVAS-C Version 2.0 (1 Nov 98)
//===========================================================================

/*
   Global variables.

   'PSI_COR' and 'EPS_COR' are celestial pole offsets for high-
   precision applications.  See function 'cel_pole' for more details.
*/

static double PSI_COR = 0.0;
static double EPS_COR = 0.0;


/********earthtilt */

static void earthtilt (double tjd, 

                double *mobl, double *tobl, double *eq, double *dpsi,
                double *deps)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Computes quantities related to the orientation of the Earth's
      rotation axis at Julian date 'tjd'.

   REFERENCES: 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.
      Transactions of the IAU (1994). Resolution C7; Vol. XXIIB, p. 59.
      McCarthy, D. D. (ed.) (1996). IERS Technical Note 21. IERS
         Central Bureau, Observatoire de Paris), pp. 21-22.

   INPUT
   ARGUMENTS:
      tjd (double)
         TDB Julian date of the desired time

   OUTPUT
   ARGUMENTS:
      *mobl (double)
         Mean obliquity of the ecliptic in degrees at 'tjd'.
      *tobl (double)
         True obliquity of the ecliptic in degrees at 'tjd'.
      *eq (double)
         Equation of the equinoxes in seconds of time at 'tjd'.
      *dpsi (double)
         Nutation in longitude in arcseconds at 'tjd'.
      *deps (double)
         Nutation in obliquity in arcseconds at 'tjd'.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      PSI_COR, EPS_COR, DEG2RAD 

   FUNCTIONS
   CALLED:
      nutation_angles  novas.c
      fund_args        novas.c
      fabs             math.h
      pow              math.h
      cos              math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/08-93/WTH (USNO/AA) Translate Fortran.
      V1.1/06-97/JAB (USNO/AA) Incorporate IAU (1994) and IERS (1996) 
                               adjustment to the "equation of the 
                               equinoxes".
      V1.2/10-97/JAB (USNO/AA) Implement function that computes 
                               arguments of the nutation series.
      V1.3/07-98/JAB (USNO/AA) Use global variables 'PSI_COR' and 
                               'EPS_COR' to apply celestial pole offsets
                               for high-precision applications.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'etilt'.
      2. Values of the celestial pole offsets 'PSI_COR' and 'EPS_COR'
      are set using function 'cel_pole', if desired.  See the prolog
      of 'cel_pole' for details.

------------------------------------------------------------------------
*/
{
   static double tjd_last = 0.0;
   static double t, dp, de;
   double d_psi, d_eps, mean_obliq, true_obliq, eq_eq, args[5];

/*
   Compute time in Julian centuries from epoch J2000.0.
*/

  t = (tjd - T0) / 36525.0;

/*
   Compute the nutation angles (arcseconds) from the standard nutation 
   model if the input Julian date is significantly different from the 
   last Julian date.
*/

  if (fabs (tjd - tjd_last) > 1.0e-6)
      nutation_angles (t, &dp,&de);

/*
   Apply observed celestial pole offsets.
*/

   d_psi = dp + PSI_COR;
   d_eps = de + EPS_COR;

/*
   Compute mean obliquity of the ecliptic in arcseconds.
*/

   mean_obliq = 84381.4480 - 46.8150 * t - 0.00059 * pow (t, 2.0)
      + 0.001813 * pow (t, 3.0);

/*
   Compute true obliquity of the ecliptic in arcseconds.
*/

   true_obliq = mean_obliq + d_eps;

/*
   Convert obliquity values to degrees.
*/

   mean_obliq /= 3600.0;
   true_obliq /= 3600.0;

/*
   Compute equation of the equinoxes in seconds of time.

   'args[4]' is "omega", the longitude of the ascending node of the 
   Moon's mean orbit on the ecliptic in radians.  This is also an 
   argument of the nutation series.
*/

   fund_args (t, args);

   eq_eq = d_psi * cos (mean_obliq * DEG2RAD) +
      (0.00264  * sin (args[4]) + 0.000063 * sin (2.0 * args[4]));

   eq_eq /= 15.0;
                           
/*
   Reset the value of the last Julian date and set the output values.
*/

   tjd_last = tjd;

   *dpsi = d_psi;
   *deps = d_eps;
   *eq = eq_eq;
   *mobl = mean_obliq;
   *tobl = true_obliq;

   return;
}

/********aberration */

short int aberration (double *pos, double *ve, double lighttime,

                      double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Corrects position vector for aberration of light.  Algorithm
      includes relativistic terms.

   REFERENCES: 
      Murray, C. A. (1981) Mon. Notices Royal Ast. Society 195, 639-648.
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      pos[3] (double)
         Position vector, referred to origin at center of mass of the
         Earth, components in AU.
      ve[3] (double)
         Velocity vector of center of mass of the Earth, referred to
         origin at solar system barycenter, components in AU/day.
      lighttime (double)
         Light time from body to Earth in days.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, referred to origin at center of mass of the
         Earth, corrected for aberration, components in AU

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.

   GLOBALS
   USED:
      C

   FUNCTIONS
   CALLED:
      sqrt      math.h
      pow       math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'aberat'.
      2. If 'lighttime' = 0 on input, this function will compute it.

------------------------------------------------------------------------
*/
{
   short int j;

   double p1mag, vemag, beta, dot,cosd, gammai, p, q, r;

   if (lighttime == 0.0)
   {
      p1mag = sqrt (pow (pos[0], 2.0) + pow (pos[1], 2.0)
                  + pow (pos[2], 2.0));
      lighttime = p1mag / C;
   }
    else
      p1mag = lighttime * C;

   vemag = sqrt (pow (ve[0], 2.0) + pow (ve[1], 2.0) 
               + pow (ve[2], 2.0));
   beta = vemag / C;
   dot = pos[0] * ve[0] + pos[1] * ve[1] + pos[2] * ve[2];

   cosd = dot / (p1mag * vemag);
   gammai = sqrt (1.0 - pow (beta, 2.0));
   p = beta * cosd;
   q = (1.0 + p / (1.0 + gammai)) * lighttime;
   r = 1.0 + p;

   for (j = 0; j < 3; j++)
      pos2[j] = (gammai * pos[j] + q * ve[j]) / r;

   return 0;
}

/********precession */

static void precession (double tjd1, double *pos, double tjd2,

                 double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Precesses equatorial rectangular coordinates from one epoch to
      another.  The coordinates are referred to the mean equator and
      equinox of the two respective epochs.

   REFERENCES:
      Explanatory Supplement to AE and AENA (1961); pp. 30-34.
      Lieske, J., et al. (1977). Astron. & Astrophys. 58, 1-16. 
      Lieske, J. (1979). Astron. & Astrophys. 73, 282-284. 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      tjd1 (double)
         TDB Julian date of first epoch.
      pos[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to mean equator and equinox of first epoch.
      tjd2 (double)
         TDB Julian date of second epoch.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to mean equator and equinox of second epoch.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      T0, RAD2SEC

   FUNCTIONS
   CALLED:
      sin    math.h
      cos    math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/03-98/JAB (USNO/AA) Change function type from 'short int' to
                               'void'.
      V1.3/12-99/JAB (USNO/AA) Precompute trig terms for greater
                               efficiency.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'preces'.

------------------------------------------------------------------------
*/
{
   double xx, yx, zx, xy, yy, zy, xz, yz, zz, t, t1, t02, t2, t3,
      zeta0, zee, theta, cz0, sz0, ct, st, cz, sz;

/*
   't' and 't1' below correspond to Lieske's "big T" and "little t".
*/

   t = (tjd1 - T0) / 36525.0;
   t1 = (tjd2 - tjd1) / 36525.0;
   t02 = t * t;
   t2 = t1 * t1;
   t3 = t2 * t1;

/*
   'zeta0', 'zee', 'theta' below correspond to Lieske's "zeta-sub-a",
   "z-sub-a", and "theta-sub-a".
*/

   zeta0 = (2306.2181 + 1.39656 * t - 0.000139 * t02) * t1
         + (0.30188 - 0.000344 * t) * t2 + 0.017998 * t3;

   zee = (2306.2181 + 1.39656 * t - 0.000139 * t02) * t1
       + (1.09468 + 0.000066 * t) * t2 + 0.018203 * t3;

   theta = (2004.3109 - 0.85330 * t - 0.000217 * t02) * t1
         + (-0.42665 - 0.000217 * t) * t2 - 0.041833 * t3;

   zeta0 /= RAD2SEC;
   zee /= RAD2SEC;
   theta /= RAD2SEC;

/*
   Precalculate trig terms.
*/

   cz0 = cos (zeta0);
   sz0 = sin (zeta0);
   ct = cos (theta);
   st = sin (theta);
   cz = cos (zee);
   sz = sin (zee);

/*
   Precession rotation matrix follows.
*/

   xx =  cz0 * ct * cz - sz0 * sz;
   yx = -sz0 * ct * cz - cz0 * sz;
   zx = -st * cz;
   xy = cz0 * ct * sz + sz0 * cz;
   yy = -sz0 * ct * sz + cz0 * cz;
   zy = -st * sz;
   xz = cz0 * st;
   yz = -sz0 * st;
   zz = ct;

/*
   Perform rotation.
*/

   pos2[0] = xx * pos[0] + yx * pos[1] + zx * pos[2];
   pos2[1] = xy * pos[0] + yy * pos[1] + zy * pos[2];
   pos2[2] = xz * pos[0] + yz * pos[1] + zz * pos[2];

   return;
}

/********nutate */

static short int nutate (double tjd, short int fn, double *pos, 

                  double *pos2)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Nutates equatorial rectangular coordinates from mean equator and
      equinox of epoch to true equator and equinox of epoch. Inverse
      transformation may be applied by setting flag 'fn'.

   REFERENCES: 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.

   INPUT
   ARGUMENTS:
      tdb (double)
         TDB julian date of epoch.
      fn (short int)
         Flag determining 'direction' of transformation;
            fn  = 0 transformation applied, mean to true.
            fn != 0 inverse transformation applied, true to mean.
      pos[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to mean equator and equinox of epoch.

   OUTPUT
   ARGUMENTS:
      pos2[3] (double)
         Position vector, geocentric equatorial rectangular coordinates,
         referred to true equator and equinox of epoch.

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.

   GLOBALS
   USED:
      DEG2RAD, RAD2SEC

   FUNCTIONS
   CALLED:
      earthtilt     novas.c
      cos           math.h
      sin           math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/01-93/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'nutate'.

------------------------------------------------------------------------
*/
{
   double cobm, sobm, cobt, sobt, cpsi, spsi, xx, yx, zx, xy, yy, zy,
      xz, yz, zz, oblm, oblt, eqeq, psi, eps;

   earthtilt (tjd, &oblm,&oblt,&eqeq,&psi,&eps);

   cobm = cos (oblm * DEG2RAD);
   sobm = sin (oblm * DEG2RAD);
   cobt = cos (oblt * DEG2RAD);
   sobt = sin (oblt * DEG2RAD);
   cpsi = cos (psi / RAD2SEC);
   spsi = sin (psi / RAD2SEC);

/*
   Nutation rotation matrix follows.
*/

   xx = cpsi;
   yx = -spsi * cobm;
   zx = -spsi * sobm;
   xy = spsi * cobt;
   yy = cpsi * cobm * cobt + sobm * sobt;
   zy = cpsi * sobm * cobt - cobm * sobt;
   xz = spsi * sobt;
   yz = cpsi * cobm * sobt - sobm * cobt;
   zz = cpsi * sobm * sobt + cobm * cobt;

   if (!fn)
   {

/*
   Perform rotation.
*/

      pos2[0] = xx * pos[0] + yx * pos[1] + zx * pos[2];
      pos2[1] = xy * pos[0] + yy * pos[1] + zy * pos[2];
      pos2[2] = xz * pos[0] + yz * pos[1] + zz * pos[2];
   }
    else
   {

/*
   Perform inverse rotation.
*/

      pos2[0] = xx * pos[0] + xy * pos[1] + xz * pos[2];
      pos2[1] = yx * pos[0] + yy * pos[1] + yz * pos[2];
      pos2[2] = zx * pos[0] + zy * pos[1] + zz * pos[2];
   }

   return 0;
}

/********nutation_angles */

static short int nutation_angles (double t,

                           double *longnutation, double *obliqnutation)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Provides fast evaluation of the nutation components according to
      the 1980 IAU Theory of Nutation.

   REFERENCES: 
      Kaplan, G. H. et. al. (1989). Astron. Journ. Vol. 97, 
         pp. 1197-1210, and references therein.
      Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
         Subroutines"; USNO internal document dated 20 Oct 1988;
         revised 15 Mar 1990.
      Miller, B. R. (1989). Proceedings of the ACM-SIGSAM International
         Symposium on Symbolic and Algebraic Computation; pp. 199-206.

   INPUT
   ARGUMENTS:
      t (double)
         TDB time in Julian centuries since J2000.0

   OUTPUT
   ARGUMENTS:
      *longnutation (double)
         Nutation in longitude in arcseconds.
      *obliqnutation (double)
         Nutation in obliquity in arcseconds.

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      fund_args         novas.c
      sin               math.h
      cos               math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/11-88/BRM (NIST)
      V1.1/08-93/WTH (USNO/AA): Translate Fortran.
      V1.2/10-97/JAB (USNO/AA): Add function to compute arguments.

   NOTES:
      1. This function is based on computer-generated Fortran code.
      Original Fortran code generated on 11/29/88 16:35:35 at the
      National Institutes of Standards and Technology (NIST), by 
      Bruce R. Miller.
      2. This function is the "C" version of Fortran NOVAS routine
      'nod', member 'vanut1f'.

------------------------------------------------------------------------
*/
{
   double clng[106] = {1.0,   1.0,  -1.0, -1.0,   1.0,  -1.0,  -1.0,
                      -1.0,  -1.0,  -1.0, -1.0,   1.0,  -1.0,   1.0,
                      -1.0,   1.0,   1.0, -1.0,  -1.0,   1.0,   1.0,
                      -1.0,   1.0,  -1.0,  1.0,  -1.0,  -1.0,  -1.0,
                       1.0,  -1.0,  -1.0,  1.0,  -1.0,   1.0,   2.0,
                       2.0,   2.0,   2.0,  2.0,  -2.0,   2.0,   2.0,
                       2.0,   3.0,  -3.0, -3.0,   3.0,  -3.0,   3.0,
                      -3.0,   3.0,   4.0,  4.0,  -4.0,  -4.0,   4.0,
                      -4.0,   5.0,   5.0,  5.0,  -5.0,   6.0,   6.0,
                       6.0,  -6.0,   6.0, -7.0,   7.0,   7.0,  -7.0,
                      -8.0,  10.0,  11.0, 12.0, -13.0, -15.0, -16.0,
                     -16.0,  17.0, -21.0,-22.0,  26.0,  29.0,  29.0,
                     -31.0, -38.0, -46.0, 48.0, -51.0,  58.0,  59.0,
                      63.0,  63.0,-123.0,129.0,-158.0,-217.0,-301.0,
                    -386.0,-517.0, 712.0,1426.0,2062.0,-2274.0,
                  -13187.0,-171996.0},
      clngx[14]={ 0.1,-0.1,0.1,0.1,0.1,0.1,0.2,-0.2,-0.4,0.5,1.2,
                 -1.6,-3.4,-174.2},
       cobl[64]={    1.0,    1.0,    1.0,   -1.0,   -1.0,   -1.0,
                     1.0,    1.0,    1.0,    1.0,    1.0,   -1.0,
                     1.0,   -1.0,    1.0,   -1.0,   -1.0,   -1.0,
                     1.0,   -1.0,    1.0,    1.0,   -1.0,   -2.0,
                    -2.0,   -2.0,    3.0,    3.0,   -3.0,    3.0,
                     3.0,   -3.0,    3.0,    3.0,   -3.0,    3.0,
                     3.0,    5.0,    6.0,    7.0,   -7.0,    7.0,
                    -8.0,    9.0,  -10.0,  -12.0,   13.0,   16.0,
                   -24.0,   26.0,   27.0,   32.0,  -33.0,  -53.0,
                    54.0,  -70.0,  -95.0,  129.0,  200.0,  224.0,
                  -895.0,  977.0, 5736.0,92025.0},
       coblx[8]={ -0.1, -0.1,  0.3,  0.5, -0.5, -0.6, -3.1,  8.9};

   short int i, ii, i1, i2, iop;
   short int nav1[10]={0,0,1,0,2,1,3,0,4,0},
       nav2[10]={ 0, 0, 0, 5, 1, 1, 3, 3, 4, 4},
       nav[183]={ 2, 0, 1, 1, 5, 2, 2, 0, 2, 1, 0, 3, 2, 5, 8, 1,17, 8,
                  1,18, 0, 2, 0, 8, 0, 1, 3, 2, 1, 8, 0,17, 1, 1,15, 1,
                  2,21, 1, 1, 2, 8, 2, 0,29, 1,21, 2, 2, 1,29, 2, 0, 9,
                  2, 5, 4, 2, 0, 4, 0, 1, 9, 2, 1, 4, 0, 2, 9, 2, 2, 4,
                  1,14,44, 2, 0,45, 2, 5,44, 2,50, 0, 1,36, 2, 2, 5,45,
                  1,37, 2, 2, 1,45, 2, 1,44, 2,53, 1, 2, 8, 4, 1,40, 3,
                  2,17, 4, 2, 0,64, 1,39, 8, 2,27, 4, 1,50,18, 1,21,47,
                  2,44, 3, 2,44, 8, 2,45, 8, 1,46, 8, 0,67, 2, 1, 5,74,
                  1, 0,74, 2,50, 8, 1, 5,78, 2,17,53, 2,53, 8, 2, 0,80,
                  2, 0,81, 0, 7,79, 1, 7,81, 2, 1,81, 2,24,44, 1, 1,79,
                  2,27,44},
      llng[106]={ 57, 25, 82, 34, 41, 66, 33, 36, 19, 88, 18,104, 93,
                  84, 47, 28, 83, 86, 69, 75, 89, 30, 58, 73, 46, 77,
                  23, 32, 59, 72, 31, 16, 74, 22, 98, 38, 62, 96, 37,
                  35,  6, 76, 85, 51, 26, 10, 13, 63,105, 52,102, 67,
                  99, 15, 24, 14,  3,100, 65, 11, 55, 68, 20, 87, 64,
                  95, 27, 60, 61, 80, 91, 94, 12, 43, 71, 42, 97, 70,
                   7, 49, 29,  2,  5, 92, 50, 78, 56, 17, 48, 40, 90,
                   8, 39, 54, 81, 21,103, 53, 45,101,  0,  1,  9, 44,
                  79,  4},
      llngx[14]={ 81, 7, 97, 0, 39, 40, 9, 44, 45,103,101, 79, 1, 4},
      lobl[64]={  51, 98, 17, 21,  5,  2, 63,105, 38, 52,102, 62, 96,
                  37, 35, 76, 36, 88, 85,104, 93, 84, 83, 67, 99,  8,
                  68,100, 60, 61, 91, 87, 64, 80, 95, 65, 55, 94, 43,
                  97,  0, 71, 70, 42, 49, 92, 50, 78, 56, 90, 48, 40,
                  39, 54,  1, 81,103, 53, 45,101,  9, 44, 79,  4},
      loblx[8] ={ 53,  1,103,  9, 44,101, 79,  4};

   double a[5], angle, cc, ss1, cs, sc, c[106], s[106], lng, lngx, obl,
      oblx;

/*
   Compute the arguments of the nutation series in radians.
*/

   fund_args (t, a);

/*
   Evaluate the series.
*/

   i = 0;
   for (ii = 0; ii < 10; ii += 2)
   {
      angle = a[nav1[ii]] * (double) (nav1[1+ii]+1);
      c[i] = cos (angle);
      s[i] = sin (angle);
      i += 1;
   }

   i = 5;
   for (ii = 0; ii < 10; ii += 2)
   {
      i1 = nav2[ii];
      i2 = nav2[1+ii];

      c[i] = c[i1] * c[i2] - s[i1] * s[i2];
      s[i] = s[i1] * c[i2] + c[i1] * s[i2];
      i += 1;
   }

   i = 10;
   for (ii = 0; ii < 183; ii += 3)
   {
      iop = nav[ii];
      i1 = nav[1+ii];
      i2 = nav[2+ii];
      switch (iop)
      {
         case 0:
            c[i] = c[i1] * c[i2] - s[i1] * s[i2];
            s[i] = s[i1] * c[i2] + c[i1] * s[i2];
            i += 1;
            break;
         case 1:
            c[i] = c[i1] * c[i2] + s[i1] * s[i2];
            s[i] = s[i1] * c[i2] - c[i1] * s[i2];
            i += 1;
            break;
         case 2:
            cc = c[i1] * c[i2];
            ss1 = s[i1] * s[i2];
            sc = s[i1] * c[i2];
            cs = c[i1] * s[i2];
            c[i] = cc - ss1;
            s[i] = sc + cs;
            i += 1;
            c[i] = cc + ss1;
            s[i] = sc - cs;
            i += 1;
            break;
      }
      if (iop == 3)
         break;
   }

   lng = 0.0;
   for (i = 0; i < 106; i++)
      lng += clng[i] * s[llng[i]];

   lngx = 0.0;
   for (i = 0; i < 14; i++)
      lngx += clngx[i] * s[llngx[i]];

   obl = 0.0;
   for (i = 0; i < 64; i++)
      obl += cobl[i] * c[lobl[i]];

   oblx = 0.0;
   for (i = 0; i < 8; i++)
      oblx += coblx[i] * c[loblx[i]];

   *longnutation = (lng + t * lngx) / 10000.0;
   *obliqnutation = (obl + t * oblx) / 10000.0;

   return 0;
}

/********fund_args */

static void fund_args (double t,

                double a[5])
/*
------------------------------------------------------------------------

   PURPOSE:
      To compute the fundamental arguments.

   REFERENCES:
      Seidelmann, P.K. (1982) Celestial Mechanics 27, 79-106 (1980 IAU 
         Theory of Nutation).

   INPUT
   ARGUMENTS:
      t (double)
         TDB time in Julian centuries since J2000.0

   OUTPUT
   ARGUMENTS:
      a[5] (double)
         Fundamental arguments, in radians:
          a[0] = l (mean anomaly of the Moon)
          a[1] = l' (mean anomaly of the Sun)
          a[2] = F (L - omega; L = mean longitude of the Moon)
          a[3] = D (mean elongation of the Moon from the Sun)
          a[4] = omega (mean longitude of the Moon's ascending node)

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      TWOPI

   FUNCTIONS
   CALLED:
      fmod     math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/10-97/JAB (USNO/AA)
      V1.1/07-98/JAB (USNO/AA): Place arguments in the range 0-TWOPI
                                radians.

   NOTES:
      1. The fundamental arguments are used in computing the nutation
      angles and in the expression for sidereal time.

------------------------------------------------------------------------
*/
{
   short int i;

   a[0] = 2.3555483935439407 + t * (8328.691422883896
                             + t * (1.517951635553957e-4
                             + 3.1028075591010306e-7 * t));
   a[1] = 6.240035939326023 + t * (628.3019560241842
                            + t * (-2.7973749400020225e-6
                            - 5.817764173314431e-8 * t));
   a[2] = 1.6279019339719611 + t * (8433.466158318453 
                             + t * (-6.427174970469119e-5
                             + 5.332950492204896e-8 * t));
   a[3] = 5.198469513579922 + t * (7771.377146170642
                            + t * (-3.340851076525812e-5
                            + 9.211459941081184e-8 * t));
   a[4] = 2.1824386243609943 + t * (-33.75704593375351
                             + t * (3.614285992671591e-5
                             + 3.878509448876288e-8 * t));

   for (i = 0; i < 5; i++)
   {
      a[i] = fmod (a[i],TWOPI);
      if (a[i] < 0.0)
         a[i] += TWOPI; 
   }

   return;
}

/********radec2vector */

static void radec2vector (double ra, double dec, double dist,
                   double *vector)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Converts equatorial spherical coordinates to a vector (equatorial
      rectangular coordinates).

   REFERENCES: 
      None.

   INPUT
   ARGUMENTS:
      ra (double)
         Right ascension (hours).
      dec (double)
         Declination (degrees).

   OUTPUT
   ARGUMENTS:
      vector[3] (double)
         Position vector, equatorial rectangular coordinates (AU).

   RETURNED
   VALUE:
      (short int)
         0...Everything OK.

   GLOBALS
   USED:
      DEG2RAD

   FUNCTIONS
   CALLED:
      cos     math.h
      sin     math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/05-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.

   NOTES:
      None.

------------------------------------------------------------------------
*/
{

   vector[0] = dist * cos (DEG2RAD * dec) * cos (DEG2RAD * 15.0 * ra);
   vector[1] = dist * cos (DEG2RAD * dec) * sin (DEG2RAD * 15.0 * ra);
   vector[2] = dist * sin (DEG2RAD * dec);

   return;
}

/********tdb2tdt */

static void tdb2tdt (double tdb,

              double *tdtjd, double *secdiff)
/*
------------------------------------------------------------------------

   PURPOSE:    
      Computes the terrestrial time (TT) or terrestrial dynamical time 
      (TDT) Julian date corresponding to a barycentric dynamical time 
      (TDB) Julian date.

   REFERENCES: 
      Explanatory Supplement to the Astronomical Almanac, pp. 42-44 and 
         p. 316.

   INPUT
   ARGUMENTS:
      tdb (double)
         TDB Julian date.

   OUTPUT
   ARGUMENTS:
      *tdtjd (double)
         TT (or TDT) Julian date.
      *secdiff (double)
         Difference tdbjd-tdtjd, in seconds.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      RAD2SEC, T0

   FUNCTIONS
   CALLED:
      sin   math.h
      fmod  math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/07-92/TKB (USNO/NRL Optical Interfer.) Translate Fortran.
      V1.1/08-93/WTH (USNO/AA) Update to C Standards.
      V1.2/06-98/JAB (USNO/AA) New algorithm (see reference).

   NOTES:
      1. Expressions used in this version are approximations resulting
      in accuracies of about 20 microseconds.
      2. This function is the "C" version of Fortran NOVAS routine
      'times'.
------------------------------------------------------------------------
*/
{

/*
   'ecc' = eccentricity of earth-moon barycenter orbit.
*/

   double ecc = 0.01671022;
   double rev = 1296000.0;
   double tdays, m, l, lj, e;

   tdays = tdb - T0;
   m = ( 357.51716 + 0.985599987 * tdays ) * 3600.0;
   l = ( 280.46435 + 0.985609100 * tdays ) * 3600.0;
   lj = ( 34.40438 + 0.083086762 * tdays ) * 3600.0;
   m = fmod (m,rev) / RAD2SEC;
   l = fmod (l,rev) / RAD2SEC;
   lj = fmod (lj,rev) / RAD2SEC;
   e = m + ecc * sin (m) + 0.5 * ecc * ecc * sin (2.0 * m);
   *secdiff = 1.658e-3 * sin (e) + 20.73e-6 * sin (l - lj);
   *tdtjd = tdb - *secdiff / 86400.0;

    return;
}

/********julian_date */

static double julian_date (short int year, short int month, short int day,
                    double hour)
/*
------------------------------------------------------------------------

   PURPOSE:
      This function will compute the Julian date for a given calendar
      date (year, month, day, hour).

   REFERENCES: 
      Fliegel & Van Flandern, Comm. of the ACM, Vol. 11, No. 10, October
      1968, p. 657.

   INPUT
   ARGUMENTS:
      year (short int)
         Year.
      month (short int)
         Month number.
      day (short int)
         Day-of-month.
      hour (double)
         Hour-of-day.

   OUTPUT
   ARGUMENTS:
      None.

   RETURNED
   VALUE:
      (double)
         Julian date.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      None.

   VER./DATE/
   PROGRAMMER:
      V1.0/06-98/JAB (USNO/AA)

   NOTES:
      1. This function is the "C" version of Fortran NOVAS routine
      'juldat'.
      2. This function makes no checks for a valid input calendar
      date.
------------------------------------------------------------------------
*/
{
   long int jd12h;

   double tjd;

   jd12h = (long) day - 32075L + 1461L * ((long) year + 4800L
      + ((long) month - 14L) / 12L) / 4L
      + 367L * ((long) month - 2L - ((long) month - 14L) / 12L * 12L)
      / 12L - 3L * (((long) year + 4900L + ((long) month - 14L) / 12L)
      / 100L) / 4L;
   tjd = (double) jd12h - 0.5 + hour / 24.0;

   return (tjd);
}

/********cal_date */

/*
static void cal_date (double tjd,
               short int *year, short int *month, short int *day,
               double *hour)
*/

/*
------------------------------------------------------------------------

   PURPOSE:    
      This function will compute a date on the Gregorian calendar given
      the Julian date.

   REFERENCES: 
      Fliegel & Van Flandern, Comm. of the ACM, Vol. 11, No. 10,
         October 1968, p. 657.

   INPUT
   ARGUMENTS:
      tjd (double)
         Julian date.

   OUTPUT
   ARGUMENTS:
      *year (short int)
         Year.
      *month (short int)
         Month number.
      *day (short int)
         Day-of-month.
      *hour (double)
         Hour-of-day.

   RETURNED
   VALUE:
      None.

   GLOBALS
   USED:
      None.

   FUNCTIONS
   CALLED:
      fmod     math.h

   VER./DATE/
   PROGRAMMER:
      V1.0/06-98/JAB (USNO/AA)

   NOTES:
      1. This routine valid for any 'jd' greater than zero.
      2. Input julian date can be based on any UT-like time scale
      (UTC, UT1, TT, etc.) - output time value will have same basis.
      3. This function is the "C" version of Fortran NOVAS routine
      'caldat'.


------------------------------------------------------------------------
*/
/*
{
   long int jd, k, m, n;

   double djd;

   djd = tjd + 0.5;
   jd = (long int) djd;

   *hour = fmod (djd,1.0) * 24.0;

   k     = jd + 68569L;
   n     = 4L * k / 146097L;

   k     = k - (146097L * n + 3L) / 4L;
   m     = 4000L * (k + 1L) / 1461001L;
   k     = k - 1461L * m / 4L + 31L;

   *month = (short int) (80L * k / 2447L);
   *day   = (short int) (k - 2447L * (long int) *month / 80L);
   k      = (long int) *month / 11L;

   *month = (short int) ((long int) *month + 2L - 12L * k);
   *year  = (short int) (100L * (n - 49L) + m + k);

   return;
}
*/

//===========================================================================
//  End of novas.cpp
//===========================================================================

//===========================================================================
// Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
//
// mphases.h
//
// Computes the julian dates (Beijing Time) of a given moon phase (0-3, for
// new moon to last quarter) between two given time marks
//===========================================================================

/* Inputs:
   tstart, start of the period in julian date (Beijing Time)
   tend, end of the period in julian date (Beijing Time)
   phase, the desired phase, 0-3 for new moon to last quarter
   Output:
   vjdphases, the julian dates (Beijing Time) of all occurings of the
              desired phase
*/
static void mphases(double tstart, double tend, int phase, vdouble& vjdphases);

//===========================================================================
//  End of mphases.h
//===========================================================================

//===========================================================================
// Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
//
// mphases.cpp
//===========================================================================

/* Inputs:
   tstart, start of the period in julian date (Beijing Time)
   tend, end of the period in julian date (Beijing Time)
   phase, the desired phase, 0-3 for new moon to last quarter
   Output:
   vjdphases, the julian dates (Beijing Time) of all occurings of the
              desired phase
*/
static void mphases(double tstart, double tend, int phase, vdouble& vjdphases)
{
    double offs;
    if (tstart < 2425245) /* Before 1928 */
        offs = (116.0 + 25.0 / 60.0) / 360.0;
    else
        offs = 120.0 / 360.0;
/* Find the lunation number of the first given phase after tstart */
    double period = 29.53058853;
    int lun = (int) ((tstart - offs - 2451550.09765 - phase * 7.375) / period);
    double jd1 = moonphasebylunation(lun, phase) + offs;
    while (jd1 - tstart > 29)
    {
        lun--;
        jd1 = moonphasebylunation(lun, phase) + offs;
    }
    while (tstart > jd1)
    {
        lun++;
        jd1 = moonphasebylunation(lun, phase) + offs;
    }
/* Compute subsequent phases until after tend */
    vjdphases.erase(vjdphases.begin(), vjdphases.end());
    vjdphases.push_back(jd1);
    while (jd1 < tend - 29)
    {
        lun++;
        jd1 = moonphasebylunation(lun, phase) + offs;
        if (jd1 < tend)
            vjdphases.push_back(jd1);
    }
}

//===========================================================================
//  End of mphases.cpp
//===========================================================================

//===========================================================================
// Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
//
// solarterm.h
//
// Computes Beijing time of solarterms for a given year
//===========================================================================

/* Given year, computes the julian date (Beijing Time) of the winter solstice
   of previous year (jdpws) and all solarterms of the given year (vjdterms)
*/
static void solarterm(short int year, double& jdpws, vdouble& vjdterms);

//===========================================================================
//  End of solarterm.h
//===========================================================================

//===========================================================================
// Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
//
// solarterm.cpp
//===========================================================================

#define PI 3.141592653589793

/* Computes the angle in degrees of given UTC time,
   with 0 being winter solstice
*/
static double timeangle(double t)
{
    double sposg[3], spos[3];
    double eposb[3], evelb[3], eposh[3], evelh[3];
    double y, angle;
    double dummy, secdiff, tdb, lighttime;

    tdb2tdt(t, &dummy, &secdiff);
    tdb = t + secdiff / 86400.0;
    solarsystem(tdb, 3, HELIOC, eposh, evelh);
    solarsystem(tdb, 3, BARYC, eposb, evelb);
    /* Convert to geocentric sun position */
    spos[0] = -eposh[0];
    spos[1] = -eposh[1];
    spos[2] = -eposh[2];
    lighttime = sqrt(spos[0]*spos[0] + spos[1]*spos[1] + spos[2]*spos[2]) / C;
    aberration(spos, evelb, lighttime, sposg);
    precession(T0, sposg, tdb, spos);
    nutate(tdb, FN0, spos, sposg);
    y = sqrt(sposg[1] * sposg[1] + sposg[2] * sposg[2]);
    y *= (sposg[1] > 0.0) ? 1.0 : -1.0;
    angle = atan2(y, sposg[0]) + PI / 2.0;   
    if (angle < 0.0)
        angle += 2.0 * PI;
    angle *= 180.0 / PI;
    return angle;
}

/* Computes the UTC time for given solar term:
   tstart:  search start time
   tend:    search end time
   termang: angle value of the given term
   returns the UTC time for the solar term
*/
static double termtime(double tstart, double tend, double termang)
{
    double errlimit;
    double tl, tu, tdif, f;
    double vs, ve, vl, vu;

    errlimit = 1.0 / 2880.0; /* Half a minute */
    f = (sqrt(5.0) - 1.0) / 2.0;
    vs = fabs(timeangle(tstart) - termang);
    ve = fabs(timeangle(tend) - termang);
    tdif = f * (tend - tstart);
    tl = tend - tdif;
    tu = tstart + tdif;
    vl = fabs(timeangle(tl) - termang);
    vu = fabs(timeangle(tu) - termang);
    while (tend - tstart > errlimit)
    {
        if (vl < vu && vl < vs && vl < ve ||
            vs < vu && vs < vl && vs < ve)
        {
            tend = tu;
            ve = vu;
            tu = tl;
            vu = vl;
            tdif = f * (tend - tstart);
            tl = tend - tdif;
            vl = fabs(timeangle(tl) - termang);
        }
        else
        {
            tstart = tl;
            vs = vl;
            tl = tu;
            vl = vu;
            tdif = f * (tend - tstart);
            tu = tstart + tdif;
            vu = fabs(timeangle(tu) - termang);
        }
    }
    return (tstart + tend) / 2.0;
}

/* Given year, computes the julian date (Beijing Time) of the winter solstice
   of previous year (jdpws) and all solarterms of the given year (vjdterms)
*/
static void solarterm(short int year, double& jdpws, vdouble& vjdterms)
{
    int dstart, dend;
    short int j, month, day;
    double angle, offs;

    if (year < 1928)
        offs = (116.0 + 25.0 / 60.0) / 360.0;
    else
        offs = 120.0 / 360.0;
    /* Determine the time of winter solstice of previous year */
    dstart = (int) julian_date(year - 1, 12, 18, 12.0);
    dend = (int) julian_date(year - 1, 12, 25, 12.0);
    jdpws = termtime(dstart, dend, 0.0) + offs;
    vjdterms.resize(24);
    for (j = 0; j < 24; j++)
    {
        month = j / 2 + 1;
        day = (j % 2) * 14;
        dstart = (int) julian_date(year, month, 1 + day, 12.0);
        dend = (int) julian_date(year, month, 10 + day, 12.0);
        angle = (j + 1) * 15.0;
        vjdterms[j] = termtime(dstart, dend, angle) + offs;
    }
}

//===========================================================================
//  End of solarterm.cpp
//===========================================================================

//===========================================================================
// Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
//
// lunaryear.h
//
// Determines the lunar month numbers for a given year
//===========================================================================

/* Inputs:
   year
   Outputs:
   vterms, vector of solarterm times for the given year
   lastnew, julian day of last new moon for previous year
   lastmon, month number for that started on last new moon day of previous year
   vmoons, vector of new moon times for the given year
   vmonth, vector of lunar month numbers for the given year with half
           advance for leap month, i.e. 7.5 for leap 7th month
   nextnew, julian day of first new moon for next year
   Returns:
   calendar month number in which a leap lunar month begins with 0.5 added if
   this leap lunar month runs into next calendar month
*/

static double lunaryear(short int year, vdouble& vterms, double& lastnew,
                 double& lastmon, vdouble& vmoons, vdouble& vmonth,
                 double& nextnew);

//===========================================================================
//  End of lunaryear.h
//===========================================================================

//===========================================================================
// Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
//
// lunaryear.cpp
//===========================================================================

/* Inputs:
   jdws, Time of winter solstice of given year
   vmoons, vector of new moon times since day after winter solstice of
           previous year
   Returns:
   true -- normal year, false -- leap year
*/
static bool IsNormalYear(double jdws, vdouble& vmoons)
{
    size_t i = 0, nmonth = 0;
    double jend = jdws + 0.5;
    jend = int(jend) + 0.5;
    while (vmoons[i++] < jend && i <= vmoons.size())
        nmonth++;
    //assert(nmonth == 12 || nmonth == 13);
    if (nmonth == 12) /* Normal year */
      return true;
    return false;
}

/* Inputs:
   jstart, time of new moon for the starting of the month
   jend, time of new moon for the starting of the next month
   vterms, vector of solarterms of the given year
   Returns:
   true -- a Zhongqi is in, false -- no Zhongqi
*/
static bool IsZhongQiInMonth(double jstart, double jend, vdouble& vterms)
{
    jstart -= 0.5;
    jstart = int(jstart) + 0.5;
    jend -= 0.5;
    jend = int(jend) + 0.5;
    size_t i;
    for (i = 1; i < vterms.size(); i += 2)
    {
        if (vterms[i] >= jstart && vterms[i] < jend)
            return true;
        if (vterms[i] >= jend)
            return false;
    }
    return false;
}

/* Inputs:
   vjds, vector of julian dates to have the hour portion trimmed
   Outputs:
   vjds, vector of julian dates with the hour portion trimmed
*/
static void TrimHour(vdouble& vjds)
{
    size_t i, t;
    for (i = 0; i < vjds.size(); i++)
    {
        t = size_t(vjds[i] + 0.5);
        vjds[i] = double(t);
    }
}

/* Inputs:
   year
   Outputs:
   vterms, vector of solarterm days for the given year
   lastnew, julian day for the last new moon of previous year
   lastmon, month number for that started on last new moon day of previous year
   vmoons, vector of new moon days for the given year
   vmonth, vector of lunar month numbers for the given year with half
           advance for leap month, i.e. 7.5 for leap 7th month
   nextnew, julian day for the first new moon of next year
   Returns:
   calendar month number in which a leap lunar month begins with 0.5 added if
   this leap lunar month runs into next calendar month
*/
static double lunaryear(short int year, vdouble& vterms, double& lastnew,
                 double& lastmon, vdouble& vmoons, vdouble& vmonth,
                 double& nextnew)
{
#ifdef USE_YEARCACHE
    /* Use cache if in range */
    if (YearCache::IsInRange(year)) /* Cached year */
    {
        YearCache::GetYear(year, vterms, lastnew, lastmon, vmoons, vmonth,
                           nextnew);
    }
    else
#endif
    {
        /* Determine solar terms */
        double jdpws;
        solarterm(year, jdpws, vterms);
        /* Determine new moons since day after previous winter solstice */
        double jstart = jdpws + 0.5;
        jstart = int(jstart) + 0.5;
        double jend = julian_date(year, 12, 31, 23.999);
        mphases(jstart, jend, 0, vmoons);
        /* Determine the month numbers */
        vmonth.resize(vmoons.size());
        size_t i;
        if (IsNormalYear(vterms[23], vmoons)) /* Normal year */
        {
            size_t n = 12;
            for (i = 0; i < vmoons.size(); i++)
            {
                vmonth[i] = n;
                if (n == 12)
                    n = 0;
                n++;
            }
            if (n == 1) /* Involves next lunar year */
            {
                double jdnwsp;
                vdouble vjdnterms;
                solarterm(year + 1, jdnwsp, vjdnterms);
                vdouble vnmoons;
                jstart = jdnwsp + 0.5;
                jstart = int(jstart) + 0.5;
                jend = julian_date(year + 1, 12, 31, 23.999);
                mphases(jstart, jend, 0, vnmoons);
                if (!IsNormalYear(vjdnterms[23], vnmoons))
                {
                    /* Check if a zhongqi falls inside the month */
                    if (!IsZhongQiInMonth(vnmoons[0], vnmoons[1], vjdnterms))
                        vmonth.back() = 11.5;
                }
            }
        }
        else /* Leap year */
        {
            bool bleaped = false;
            size_t n = 11;
            for (i = 0; i < vmoons.size() - 1; i++)
            {
                /* Check if a zhongqi falls inside the month */
                if (bleaped || IsZhongQiInMonth(vmoons[i], vmoons[i + 1], vterms))
                    vmonth[i] = ++n;
                else
                {
                    vmonth[i] = n + 0.5;
                    bleaped = true;
                }
                if (n == 12)
                    n = 0;
            }
            vmonth.back() = n + 1;
        }
        if (vmoons[0] < julian_date(year, 1, 1, 0.0))
        {
            lastnew = floor(vmoons.front() + 0.5);
            lastmon = vmonth.front();
            vmoons.erase(vmoons.begin());
            vmonth.erase(vmonth.begin());
        }
        else /* Need to find the last new moon for previous year */
        {
            vdouble vlastnew;
            mphases(vmoons[0] - 35.0, vmoons[0] - 25.0, 0, vlastnew);
            TrimHour(vlastnew);
            lastnew = vlastnew.back();
            lastmon = 11.0;
        }
        /* Need to find the first new moon for next year */
        vdouble vnextnew;
        mphases(vmoons.back() + 25.0, vmoons.back() + 35.0, 0, vnextnew);
        TrimHour(vnextnew);
        nextnew = vnextnew.back();
        /* Convert to whole day numbers */
        TrimHour(vmoons);
        TrimHour(vterms);
    }
    /* Scan for leap month and return the calendar month */
    if (int(lastmon + 0.9) != int(lastmon)) /* lastmon is a leap month */
    {
        if (julian_date(year, 1, 1, 12.0) < vmoons[0]) /* runs into new year */
            return 0.5;
    }
    short int monnum = 2;
  size_t i;
    for (i = 0; i < vmoons.size(); i++)
    {
        if (int(vmonth[i] + 0.9) != int(vmonth[i])) /* found leap month */
        {
            double jdfirst;
            while (monnum <= 12 && (jdfirst = julian_date(year, monnum, 1, 12.0)) < vmoons[i])
                monnum++;
            if (monnum == 13)
                return 12.0;
            /* See if leap month runs into next month */
            if (i != vmoons.size() - 1 && jdfirst < vmoons[i + 1])
                return (monnum - 0.5); /* Yes */
            else
                return (monnum - 1.0); /* No */
        }
    }
    return 0.0;
}

//===========================================================================
//  End of lunaryear.cpp
//===========================================================================

class CMonthInfo {
public:
  double month;
  int beginday; // julian day
  int endday; // julian day
  int total;
};

class CYearInfo {
public:
  int year;
  bool leap;
  double leapmonth;
  int beginday; // julian day
  int endday; // julian day
  int total;
/*
  int firstweekday;
  int lastweekday;
  int totalweeks;
*/
  vector<CMonthInfo> monthinfos;
};

#define YEAR_INTERVAL 2697

static class CYearInfo* getYearInfo( int year )
{
  static CYearInfo result;
  if ( year < 4342 || year > 9999 ) {
    result.year = -1;
    return (&result);
  }

  if ( year == result.year )
    return (&result);

  int gYear1 = year - YEAR_INTERVAL;
  result.year = year;

  vdouble vterms, vmoons, vmonth;
  double lastnew, lastmon, nextnew;
  double lmon = lunaryear(gYear1, vterms, lastnew, lastmon, vmoons, vmonth, nextnew);

  if ( lmon > 0.0 ) {
    result.leap = true;
  } else {
    result.leap = false;
  }

  result.monthinfos.clear();

  bool sign = false;
  for ( uint i = 0; i < vmoons.size(); ++i ) {
    if ( !sign && vmonth[i] > 10. ) {
      continue;
    } else {
      sign = true;
    }

    CMonthInfo mi;
    mi.month = vmonth[i];
    mi.beginday = vmoons[i];

    result.monthinfos.push_back(mi);
  }

  int idx = result.monthinfos.size();

  if ( result.monthinfos[idx-1].month == 12.0 ) {
    result.monthinfos[idx-1].endday = nextnew - 1;
  } else {
    int gYear2 = gYear1 + 1;
    vdouble vterms2, vmoons2, vmonth2;
    double lastnew2, lastmon2, nextnew2;
    lmon = lunaryear(gYear2, vterms2, lastnew2, lastmon2, vmoons2, vmonth2, nextnew2);

    CMonthInfo mi;
    mi.month = vmonth2[0];
    mi.beginday = vmoons2[0];
    mi.endday = vmoons2[1] - 1;

    result.monthinfos.push_back(mi);
  }

  idx = result.monthinfos.size();
  for ( int i = 0; i < idx - 1; ++i ) {
    result.monthinfos[i].endday = result.monthinfos[i+1].beginday - 1;
    result.monthinfos[i].total = result.monthinfos[i].endday - result.monthinfos[i].beginday + 1;
  }

  result.monthinfos[idx-1].total = result.monthinfos[idx-1].endday - result.monthinfos[idx-1].beginday + 1;

  result.beginday = result.monthinfos[0].beginday;
  result.endday = result.monthinfos[idx-1].endday;
  result.total = result.endday - result.beginday + 1;

  if ( result.leap ) {
    idx = result.monthinfos.size();
    for ( int i = 0; i < idx; ++i ) {
      double month = result.monthinfos[i].month;
      int mon = floor(month);
      if ( month == mon + 0.5 ) {
        result.leapmonth = month;
        break;
      }
    }
  }

  return (&result);
}

static bool gregorianToChinese( const QDate &date, int &year, int &month, int &day, double &rmonth )
{
  int jday = date.toJulianDay();

  year = -1;
  month = -1;
  day = -1;
  rmonth = 0.0;

  year = date.year() + YEAR_INTERVAL;

  class CYearInfo *ci = getYearInfo( year );

  if ( ci->year == -1 ) {
    year = year - 1;
    ci = getYearInfo( year );

    if ( ci->year == -1 ) {
      year = -1;
      return false;
    }

    if ( jday > ci->endday ) {
      year = -1;
      return false;
    }
  } else {
  if ( jday < ci->beginday ) {
    year = year - 1;
    ci = getYearInfo( year );

    if ( ci->year == -1 ) {
      year = -1;
      return false;
    }
  }
  }

  for ( uint i = 0; i < ci->monthinfos.size(); ++i ) {
    if ( jday >= ci->monthinfos[i].beginday && jday <= ci->monthinfos[i].endday ) {
      month = i + 1;
      day = jday - ci->monthinfos[i].beginday + 1;
      rmonth = ci->monthinfos[i].month;
      return true;
    }
  }
  return false;
}

static bool chineseToGregorian( QDate & date, int y, int m, int day, double rm = 0.0 )
{
  date = QDate();

  class CYearInfo *ci = getYearInfo( y );

  if ( ci->year == -1 )
    return false;

  int n = ci->monthinfos.size();

  int mon;
  if ( rm < 1.0 ) {
    mon = m;
  } else {
    mon = floor(rm);
    if ( ci->leap ) {
      if ( ci->leapmonth <= rm )
        mon = mon + 1;
    }
  }

  if ( mon < 1 || mon > n )
    return false;

  int total = ci->monthinfos[mon-1].total;

  if ( day < 1 || day > total )
    return false;

  int jday = ci->monthinfos[mon-1].beginday + day - 1;

  date = QDate::fromJulianDay( jday );
  return true;
}

KCalendarSystemChinese::KCalendarSystemChinese( const KLocale * locale )
        : KCalendarSystem( locale ), d( 0 )
{
}

KCalendarSystemChinese::~KCalendarSystemChinese()
{
}

QString KCalendarSystemChinese::calendarType() const
{
    return QLatin1String( "chinese" );
}

QDate KCalendarSystemChinese::epoch() const
{
    return QDate( 1645, 1, 28);
}

QDate KCalendarSystemChinese::earliestValidDate() const
{
    // The limitation is from ccal which is 1645-01-28/4342-01-01
    return QDate( 1645, 1, 28);
}

QDate KCalendarSystemChinese::latestValidDate() const
{
    // Set to last day of year 9999 until confirm date formats & widets support > 9999
    // Last day of Chinese year 9999 is 9999-12-30
    // Which in Gregorian is 7303-02-05
    // Which is jd xxxx FIXME Find out jd and use that instead
    // Can't call setDate( 9999, 12, 29 ) as it creates circular reference!
    return QDate( 7303, 2, 5 );
}

bool KCalendarSystemChinese::isValid( int y, int month, int day ) const
{
  QDate date;
  return chineseToGregorian( date, y, month, day );
}

bool KCalendarSystemChinese::isValid( const QDate &date ) const
{
    return KCalendarSystem::isValid( date );
}

bool KCalendarSystemChinese::setDate( QDate &date, int year, int month, int day ) const
{
  date = QDate();
  return KCalendarSystem::setDate( date, year, month, day );
}

bool KCalendarSystemChinese::setYMD( QDate & date, int y, int m, int day ) const
{
  return chineseToGregorian( date, y, m, day );
}

int KCalendarSystemChinese::year( const QDate &date ) const
{
  int year, month, day;
  double rmonth;

  if ( gregorianToChinese( date, year, month, day, rmonth ) )
    return year;

  return -1;
}

int KCalendarSystemChinese::month( const QDate &date ) const
{
  int year, month, day;
  double rmonth;

  if ( gregorianToChinese( date, year, month, day, rmonth ) )
    return month;

  return -1;
}

int KCalendarSystemChinese::day( const QDate &date ) const
{
  int year, month, day;
  double rmonth;

  if ( gregorianToChinese( date, year, month, day, rmonth ) )
    return day;

  return -1;
}

QDate KCalendarSystemChinese::addYears( const QDate &date, int nyears ) const
{
  QDate result;

  int year, month, day;
  double rmonth;

  if ( !(gregorianToChinese( date, year, month, day, rmonth )) ) {
    return QDate();
  }

  if ( !(chineseToGregorian( result, year+nyears, month, day, rmonth )) ) {
    return QDate();
  }

  return result;
}

QDate KCalendarSystemChinese::addMonths( const QDate &date, int nmonths ) const
{
  QDate result = date;

  while ( nmonths > 0 ) {
    result = addDays( result, daysInMonth( result ) );
    --nmonths;
  }

  while ( nmonths < 0 ) {
    // get the number of days in the previous month to be consistent with
    // addMonths where nmonths > 0
    int nDaysInMonth = daysInMonth( addDays( result, -day( result ) ) );
    result = addDays( result, -nDaysInMonth );
    ++nmonths;
  }

  return result;
}

QDate KCalendarSystemChinese::addDays( const QDate &date, int ndays ) const
{
  return date.addDays( ndays );
}

int KCalendarSystemChinese::monthsInYear( const QDate &date ) const
{
  int year, month, day;
  double rmonth;

  if ( !(gregorianToChinese( date, year, month, day, rmonth )) ) {
    return -1;
  }

  if ( isLeapYear( year ) ) {
    return 13;
  } else {
    return 12;
  }
}

int KCalendarSystemChinese::weeksInYear( const QDate &date ) const
{
  return KCalendarSystem::weeksInYear( year( date ) );
}

int KCalendarSystemChinese::weeksInYear( int year ) const
{
  class CYearInfo *ci = getYearInfo( year );

  if ( ci->year == -1 ) {
    year = -1;
    return -1;
  }

  // Makes assumption that Julian Day 0 was day 1 of week
  // This is true for Julian/Gregorian calendar with jd 0 being Monday
  // We add 1 for ISO compliant numbering for 7 day week
  // Assumes we've never skipped weekdays
  int daysinweek = daysInWeek( QDate() );

  int first = ( ci->beginday % daysinweek ) + 1;
  int last = ( ci->endday % daysinweek ) + 1;

  int firstweekday;
  int lastweekday;

  if ( first > 4 ) {
    firstweekday = ci->beginday - first + 8;
  } else {
    firstweekday = ci->beginday - first + 1;
  }

  if ( last < 4 ) {
    lastweekday = ci->endday - last;
  } else {
    lastweekday = ci->endday - last + 7;
  }

  return ( (lastweekday + 1 - firstweekday) / daysinweek );
}

int KCalendarSystemChinese::daysInYear( const QDate &date ) const
{
  int year = date.year() + YEAR_INTERVAL;

  class CYearInfo *ci = getYearInfo( year );

  if ( ci->year == -1 ) {
    year = -1;
    return -1;
  }

  int jday = date.toJulianDay();
  if ( jday < ci->beginday ) {
    year = year - 1;
    ci = getYearInfo( year );

    if ( ci->year == -1 ) {
      year = -1;
      return -1;
    }
  }

  return ci->total;
}

int KCalendarSystemChinese::daysInMonth( const QDate &date ) const
{
  int year = date.year() + YEAR_INTERVAL;

  class CYearInfo *ci = getYearInfo( year );

  if ( ci->year == -1 ) {
    year = -1;
    return -1;
  }

  int jday = date.toJulianDay();
  if ( jday < ci->beginday ) {
    year = year - 1;
    ci = getYearInfo( year );

    if ( ci->year == -1 ) {
      year = -1;
      return -1;
    }
  }

  for ( uint i = 0; i < ci->monthinfos.size(); ++i ) {
    if ( jday >= ci->monthinfos[i].beginday && jday <= ci->monthinfos[i].endday ) {
      return ci->monthinfos[i].total;
    }
  }

  return -1;
}

int KCalendarSystemChinese::daysInWeek( const QDate &date ) const
{
    return KCalendarSystem::daysInWeek( date );
}

int KCalendarSystemChinese::dayOfYear( const QDate &date ) const
{
  int year = date.year() + YEAR_INTERVAL;

  class CYearInfo *ci = getYearInfo( year );

  if ( ci->year == -1 ) {
    year = -1;
    return -1;
  }

  int jday = date.toJulianDay();
  if ( jday < ci->beginday ) {
    year = year - 1;
    ci = getYearInfo( year );

    if ( ci->year == -1 ) {
      year = -1;
      return -1;
    }
  }

  return (jday - ci->beginday + 1);
}

int KCalendarSystemChinese::dayOfWeek( const QDate &date ) const
{
  return KCalendarSystem::dayOfWeek( date );
}

int KCalendarSystemChinese::weekNumber( const QDate &date, int *yearNum ) const
{
  int y = year( date );

  class CYearInfo *ci = getYearInfo( y );

  if ( ci->year == -1 ) {
    if ( yearNum )
      *yearNum = -1;
    return -1;
  }

  // Makes assumption that Julian Day 0 was day 1 of week
  // This is true for Julian/Gregorian calendar with jd 0 being Monday
  // We add 1 for ISO compliant numbering for 7 day week
  // Assumes we've never skipped weekdays
  int daysinweek = daysInWeek( date );

  int first = ( ci->beginday % daysinweek ) + 1;
  int last = ( ci->endday % daysinweek ) + 1;

  int firstweekday;
  int lastweekday;

  if ( first > 4 ) {
    firstweekday = ci->beginday - first + 8;
  } else {
    firstweekday = ci->beginday - first + 1;
  }

  if ( last < 4 ) {
    lastweekday = ci->endday - last;
  } else {
    lastweekday = ci->endday - last + 7;
  }

  int jday = date.toJulianDay();

  if ( jday < firstweekday ) {
    if ( yearNum )
      * yearNum = y - 1;
    return weeksInYear( y - 1 );
  } else if ( jday > lastweekday ) {
    if ( yearNum )
      * yearNum = y + 1;
    return 1;
  } else { // firstweekday <= jday <= lastweekday
    if ( yearNum )
      * yearNum = y;
    return ((jday - firstweekday)/daysinweek) + 1;
  }

  return -1;
}

bool KCalendarSystemChinese::isLeapYear( int year ) const
{
  class CYearInfo *ci = getYearInfo( year );

  if ( ci->year == -1 )
    return false;

  return ci->leap;
}

bool KCalendarSystemChinese::isLeapYear( const QDate &date ) const
{
  return isLeapYear( year( date ) );
}

QString KCalendarSystemChinese::monthName( int month, int year, MonthNameFormat format ) const
{
  class CYearInfo *ci = getYearInfo( year );

  if ( month < 1 ) {
    return QString();
  }

  if ( (uint)month > ci->monthinfos.size() ) {
    return QString();
  }

  bool sign;
  int mon;
  if ( !ci->leap ) {
    sign = false;
    mon = month;
  } else {
    int lm = floor( ci->leapmonth );
    if ( month <= lm ) {
      sign = false;
      mon = month;
    } else if ( month == lm + 1 ) {
      sign = true;
      mon = month - 1;
    } else {
      sign = false;
      mon = month - 1;
    }
  }

  if ( !sign ) {
    if ( format == ShortNamePossessive ) {
        switch ( mon ) {
        case 1:
            return ki18nc( "of January",   "of Jan" ).toString( locale() );
        case 2:
            return ki18nc( "of February",  "of Feb" ).toString( locale() );
        case 3:
            return ki18nc( "of March",     "of Mar" ).toString( locale() );
        case 4:
            return ki18nc( "of April",     "of Apr" ).toString( locale() );
        case 5:
            return ki18nc( "of May short", "of May" ).toString( locale() );
        case 6:
            return ki18nc( "of June",      "of Jun" ).toString( locale() );
        case 7:
            return ki18nc( "of July",      "of Jul" ).toString( locale() );
        case 8:
            return ki18nc( "of August",    "of Aug" ).toString( locale() );
        case 9:
            return ki18nc( "of September", "of Sep" ).toString( locale() );
        case 10:
            return ki18nc( "of October",   "of Oct" ).toString( locale() );
        case 11:
            return ki18nc( "of November",  "of Nov" ).toString( locale() );
        case 12:
            return ki18nc( "of December",  "of Dec" ).toString( locale() );
        default:
            return QString();
        }
    }

    if ( format == LongNamePossessive ) {
        switch ( mon ) {
        case 1:
            return ki18n( "of January" ).toString( locale() );
        case 2:
            return ki18n( "of February" ).toString( locale() );
        case 3:
            return ki18n( "of March" ).toString( locale() );
        case 4:
            return ki18n( "of April" ).toString( locale() );
        case 5:
            return ki18nc( "of May long", "of May" ).toString( locale() );
        case 6:
            return ki18n( "of June" ).toString( locale() );
        case 7:
            return ki18n( "of July" ).toString( locale() );
        case 8:
            return ki18n( "of August" ).toString( locale() );
        case 9:
            return ki18n( "of September" ).toString( locale() );
        case 10:
            return ki18n( "of October" ).toString( locale() );
        case 11:
            return ki18n( "of November" ).toString( locale() );
        case 12:
            return ki18n( "of December" ).toString( locale() );
        default:
            return QString();
        }
    }

    if ( format == ShortName ) {
        switch ( mon ) {
        case 1:
            return ki18nc( "January", "Jan" ).toString( locale() );
        case 2:
            return ki18nc( "February", "Feb" ).toString( locale() );
        case 3:
            return ki18nc( "March", "Mar" ).toString( locale() );
        case 4:
            return ki18nc( "April", "Apr" ).toString( locale() );
        case 5:
            return ki18nc( "May short", "May" ).toString( locale() );
        case 6:
            return ki18nc( "June", "Jun" ).toString( locale() );
        case 7:
            return ki18nc( "July", "Jul" ).toString( locale() );
        case 8:
            return ki18nc( "August", "Aug" ).toString( locale() );
        case 9:
            return ki18nc( "September", "Sep" ).toString( locale() );
        case 10:
            return ki18nc( "October", "Oct" ).toString( locale() );
        case 11:
            return ki18nc( "November", "Nov" ).toString( locale() );
        case 12:
            return ki18nc( "December", "Dec" ).toString( locale() );
        default:
            return QString();
        }
    }

    // Default to LongName
    switch ( mon ) {
    case 1:
        return ki18n( "January" ).toString( locale() );
    case 2:
        return ki18n( "February" ).toString( locale() );
    case 3:
        return ki18nc( "March long", "March" ).toString( locale() );
    case 4:
        return ki18n( "April" ).toString( locale() );
    case 5:
        return ki18nc( "May long", "May" ).toString( locale() );
    case 6:
        return ki18n( "June" ).toString( locale() );
    case 7:
        return ki18n( "July" ).toString( locale() );
    case 8:
        return ki18nc( "August long", "August" ).toString( locale() );
    case 9:
        return ki18n( "September" ).toString( locale() );
    case 10:
        return ki18n( "October" ).toString( locale() );
    case 11:
        return ki18n( "November" ).toString( locale() );
    case 12:
        return ki18n( "December" ).toString( locale() );
    default:
        return QString();
    }
  } else {
    if ( format == ShortNamePossessive ) {
        switch ( mon ) {
        case 1:
            return ki18nc( "of Leap January",   "of Leap Jan" ).toString( locale() );
        case 2:
            return ki18nc( "of Leap February",  "of Leap Feb" ).toString( locale() );
        case 3:
            return ki18nc( "of Leap March",     "of Leap Mar" ).toString( locale() );
        case 4:
            return ki18nc( "of Leap April",     "of Leap Apr" ).toString( locale() );
        case 5:
            return ki18nc( "of Leap May short", "of Leap May" ).toString( locale() );
        case 6:
            return ki18nc( "of Leap June",      "of Leap Jun" ).toString( locale() );
        case 7:
            return ki18nc( "of Leap July",      "of Leap Jul" ).toString( locale() );
        case 8:
            return ki18nc( "of Leap August",    "of Leap Aug" ).toString( locale() );
        case 9:
            return ki18nc( "of Leap September", "of Leap Sep" ).toString( locale() );
        case 10:
            return ki18nc( "of Leap October",   "of Leap Oct" ).toString( locale() );
        case 11:
            return ki18nc( "of Leap November",  "of Leap Nov" ).toString( locale() );
        case 12:
            return ki18nc( "of Leap December",  "of Leap Dec" ).toString( locale() );
        default:
            return QString();
        }
    }

    if ( format == LongNamePossessive ) {
        switch ( mon ) {
        case 1:
            return ki18n( "of Leap January" ).toString( locale() );
        case 2:
            return ki18n( "of Leap February" ).toString( locale() );
        case 3:
            return ki18n( "of Leap March" ).toString( locale() );
        case 4:
            return ki18n( "of Leap April" ).toString( locale() );
        case 5:
            return ki18nc( "of Leap May long", "of Leap May" ).toString( locale() );
        case 6:
            return ki18n( "of Leap June" ).toString( locale() );
        case 7:
            return ki18n( "of Leap July" ).toString( locale() );
        case 8:
            return ki18n( "of Leap August" ).toString( locale() );
        case 9:
            return ki18n( "of Leap September" ).toString( locale() );
        case 10:
            return ki18n( "of Leap October" ).toString( locale() );
        case 11:
            return ki18n( "of Leap November" ).toString( locale() );
        case 12:
            return ki18n( "of Leap December" ).toString( locale() );
        default:
            return QString();
        }
    }

    if ( format == ShortName ) {
        switch ( mon ) {
        case 1:
            return ki18nc( "Leap January", "Leap Jan" ).toString( locale() );
        case 2:
            return ki18nc( "Leap February", "Leap Feb" ).toString( locale() );
        case 3:
            return ki18nc( "Leap March", "Leap Mar" ).toString( locale() );
        case 4:
            return ki18nc( "Leap April", "Leap Apr" ).toString( locale() );
        case 5:
            return ki18nc( "Leap May short", "Leap May" ).toString( locale() );
        case 6:
            return ki18nc( "Leap June", "Leap Jun" ).toString( locale() );
        case 7:
            return ki18nc( "Leap July", "Leap Jul" ).toString( locale() );
        case 8:
            return ki18nc( "Leap August", "Leap Aug" ).toString( locale() );
        case 9:
            return ki18nc( "Leap September", "Leap Sep" ).toString( locale() );
        case 10:
            return ki18nc( "Leap October", "Leap Oct" ).toString( locale() );
        case 11:
            return ki18nc( "Leap November", "Leap Nov" ).toString( locale() );
        case 12:
            return ki18nc( "Leap December", "Leap Dec" ).toString( locale() );
        default:
            return QString();
        }
    }

    // Default to LongName
    switch ( mon ) {
    case 1:
        return ki18n( "Leap January" ).toString( locale() );
    case 2:
        return ki18n( "Leap February" ).toString( locale() );
    case 3:
        return ki18nc( "Leap March long", "Leap March" ).toString( locale() );
    case 4:
        return ki18n( "Leap April" ).toString( locale() );
    case 5:
        return ki18nc( "Leap May long", "Leap May" ).toString( locale() );
    case 6:
        return ki18n( "Leap June" ).toString( locale() );
    case 7:
        return ki18n( "Leap July" ).toString( locale() );
    case 8:
        return ki18nc( "Leap August long", "Leap August" ).toString( locale() );
    case 9:
        return ki18n( "Leap September" ).toString( locale() );
    case 10:
        return ki18n( "Leap October" ).toString( locale() );
    case 11:
        return ki18n( "Leap November" ).toString( locale() );
    case 12:
        return ki18n( "Leap December" ).toString( locale() );
    default:
        return QString();
    }
  }
}

QString KCalendarSystemChinese::monthName( const QDate& date, MonthNameFormat format ) const
{
  return monthName( month( date ), year( date ), format );
}

QString KCalendarSystemChinese::weekDayName( int weekDay, WeekDayNameFormat format ) const
{
    if ( format == ShortDayName ) {
        switch ( weekDay ) {
        case 1:  return ki18nc( "Monday",    "Mon" ).toString( locale() );
        case 2:  return ki18nc( "Tuesday",   "Tue" ).toString( locale() );
        case 3:  return ki18nc( "Wednesday", "Wed" ).toString( locale() );
        case 4:  return ki18nc( "Thursday",  "Thu" ).toString( locale() );
        case 5:  return ki18nc( "Friday",    "Fri" ).toString( locale() );
        case 6:  return ki18nc( "Saturday",  "Sat" ).toString( locale() );
        case 7:  return ki18nc( "Sunday",    "Sun" ).toString( locale() );
        default: return QString();
        }
    }

    switch ( weekDay ) {
    case 1:  return ki18n( "Monday" ).toString( locale() );
    case 2:  return ki18n( "Tuesday" ).toString( locale() );
    case 3:  return ki18n( "Wednesday" ).toString( locale() );
    case 4:  return ki18n( "Thursday" ).toString( locale() );
    case 5:  return ki18n( "Friday" ).toString( locale() );
    case 6:  return ki18n( "Saturday" ).toString( locale() );
    case 7:  return ki18n( "Sunday" ).toString( locale() );
    default: return QString();
    }
}

QString KCalendarSystemChinese::weekDayName( const QDate &date, WeekDayNameFormat format ) const
{
  return weekDayName( dayOfWeek( date ), format );
}

QString KCalendarSystemChinese::yearString( const QDate &pDate, StringFormat format ) const
{
  return KCalendarSystem::yearString( pDate, format );
}

QString KCalendarSystemChinese::monthString( const QDate &pDate, StringFormat format ) const
{
    return KCalendarSystem::monthString( pDate, format );
}

QString KCalendarSystemChinese::dayString( const QDate &pDate, StringFormat format ) const
{
  return KCalendarSystem::dayString( pDate, format );
}

int KCalendarSystemChinese::yearStringToInteger( const QString &sNum, int &iLength ) const
{
  return KCalendarSystem::yearStringToInteger( sNum, iLength );
}

int KCalendarSystemChinese::monthStringToInteger( const QString &sNum, int &iLength ) const
{
  return KCalendarSystem::monthStringToInteger( sNum, iLength );
}

int KCalendarSystemChinese::dayStringToInteger( const QString &sNum, int &iLength ) const
{
  return KCalendarSystem::yearStringToInteger( sNum, iLength );
}

QString KCalendarSystemChinese::formatDate( const QDate &date, KLocale::DateFormat format ) const
{
  return KCalendarSystem::formatDate( date, format );
}

QDate KCalendarSystemChinese::readDate( const QString &str, bool *ok ) const
{
  return KCalendarSystem::readDate( str, ok );
}

QDate KCalendarSystemChinese::readDate( const QString &intstr, const QString &fmt, bool *ok ) const
{
  return KCalendarSystem::readDate( intstr, fmt, ok );
}

QDate KCalendarSystemChinese::readDate( const QString &str, KLocale::ReadDateFlags flags, bool *ok ) const
{
  return KCalendarSystem::readDate( str, flags, ok );
}

int KCalendarSystemChinese::weekDayOfPray() const
{
  return 7; // sunday
}

int KCalendarSystemChinese::weekStartDay() const
{
  return KCalendarSystem::weekStartDay();
}

bool KCalendarSystemChinese::isLunar() const
{
  return false;
}

bool KCalendarSystemChinese::isLunisolar() const
{
  return true;
}

bool KCalendarSystemChinese::isSolar() const
{
  return false;
}

bool KCalendarSystemChinese::isProleptic() const
{
  return false;
}

bool KCalendarSystemChinese::julianDayToDate( int jd, int &year, int &month, int &day ) const
{
  QDate date = QDate::fromJulianDay( jd );
  double rmonth;
  return gregorianToChinese( date, year, month, day, rmonth );
}

bool KCalendarSystemChinese::dateToJulianDay( int year, int month, int day, int &jd ) const
{
  QDate date;
  if ( !( chineseToGregorian( date, year, month, day ) ) ) {
    jd = -1;
    return false;
  }

  jd = date.toJulianDay();
  return true;
}
