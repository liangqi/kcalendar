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
// yearcache.h
//
// Cache results for some years for faster response 
//===========================================================================

#ifndef USE_YEARCACHE
#define USE_YEARCACHE
#endif

class YearCache
{
public:
    /* Input:
       year
       Returns:
       true -- year is in range of cache, false -- otherwise
    */
    static bool IsInRange(short int year);
    /* Inputs:
       year
       Outputs:
       vterms, vector of solarterm times for the given year
       lastnew, julian day of last new moon for previous year
       lastmon, month number for that started on last new moon day of 
                previous year
       vmoons, vector of new moon times for the given year
       vmonth, vector of lunar month numbers for the given year with half
               advance for leap month, i.e. 7.5 for leap 7th month
       nextnew, julian day of first new moon for next year
    */
    static void GetYear(short int year, vdouble& vterms, double& lastnew,
        double& lastmon, vdouble& vmoons, vdouble& vmonth, double& nextnew);
};

//===========================================================================
//  End of yearcache.h
//===========================================================================

//===========================================================================
// Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
//
// yearcache.cpp
//
// Cache results for some years for faster response 
//===========================================================================

/*
   Copyright (c) 2000-2006, by Zhuo Meng (zhuo@thunder.cwru.edu).
   All rights reserved.

   Distributed under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/* Cache results for some years for faster response */
static short int yearfrom = 1950, yearto = 2050;

static double allterms[] = {
2433288.0, 2433302.0, 2433317.0, 2433332.0, 2433347.0, 2433362.0,
2433377.0, 2433392.0, 2433408.0, 2433423.0, 2433439.0, 2433455.0,
2433471.0, 2433486.0, 2433502.0, 2433518.0, 2433533.0, 2433548.0,
2433564.0, 2433579.0, 2433594.0, 2433609.0, 2433624.0, 2433638.0,
2433653.0, 2433668.0, 2433682.0, 2433697.0, 2433712.0, 2433727.0,
2433742.0, 2433758.0, 2433773.0, 2433789.0, 2433804.0, 2433820.0,
2433836.0, 2433852.0, 2433867.0, 2433883.0, 2433898.0, 2433914.0,
2433929.0, 2433944.0, 2433959.0, 2433974.0, 2433989.0, 2434004.0,
2434018.0, 2434033.0, 2434048.0, 2434063.0, 2434077.0, 2434093.0,
2434108.0, 2434123.0, 2434138.0, 2434154.0, 2434170.0, 2434185.0,
2434201.0, 2434217.0, 2434232.0, 2434248.0, 2434264.0, 2434279.0,
2434294.0, 2434309.0, 2434324.0, 2434339.0, 2434354.0, 2434369.0,
2434383.0, 2434398.0, 2434413.0, 2434428.0, 2434443.0, 2434458.0,
2434473.0, 2434488.0, 2434504.0, 2434519.0, 2434535.0, 2434551.0,
2434566.0, 2434582.0, 2434598.0, 2434613.0, 2434629.0, 2434644.0,
2434659.0, 2434675.0, 2434690.0, 2434704.0, 2434719.0, 2434734.0,
2434749.0, 2434763.0, 2434778.0, 2434793.0, 2434808.0, 2434823.0,
2434838.0, 2434853.0, 2434869.0, 2434884.0, 2434900.0, 2434916.0,
2434932.0, 2434947.0, 2434963.0, 2434979.0, 2434994.0, 2435009.0,
2435025.0, 2435040.0, 2435055.0, 2435070.0, 2435084.0, 2435099.0,
2435114.0, 2435129.0, 2435143.0, 2435158.0, 2435173.0, 2435188.0,
2435203.0, 2435219.0, 2435234.0, 2435250.0, 2435265.0, 2435281.0,
2435297.0, 2435312.0, 2435328.0, 2435344.0, 2435359.0, 2435375.0,
2435390.0, 2435405.0, 2435420.0, 2435435.0, 2435450.0, 2435464.0,
2435479.0, 2435494.0, 2435509.0, 2435524.0, 2435538.0, 2435553.0,
2435569.0, 2435584.0, 2435599.0, 2435615.0, 2435631.0, 2435646.0,
2435662.0, 2435678.0, 2435693.0, 2435709.0, 2435725.0, 2435740.0,
2435755.0, 2435770.0, 2435785.0, 2435800.0, 2435815.0, 2435830.0,
2435844.0, 2435859.0, 2435874.0, 2435889.0, 2435904.0, 2435919.0,
2435934.0, 2435949.0, 2435965.0, 2435980.0, 2435996.0, 2436012.0,
2436027.0, 2436043.0, 2436059.0, 2436074.0, 2436090.0, 2436105.0,
2436120.0, 2436136.0, 2436151.0, 2436165.0, 2436180.0, 2436195.0,
2436210.0, 2436224.0, 2436239.0, 2436254.0, 2436269.0, 2436284.0,
2436299.0, 2436314.0, 2436330.0, 2436345.0, 2436361.0, 2436377.0,
2436392.0, 2436408.0, 2436424.0, 2436439.0, 2436455.0, 2436470.0,
2436486.0, 2436501.0, 2436516.0, 2436531.0, 2436545.0, 2436560.0,
2436575.0, 2436590.0, 2436604.0, 2436619.0, 2436634.0, 2436649.0,
2436664.0, 2436680.0, 2436695.0, 2436711.0, 2436726.0, 2436742.0,
2436758.0, 2436773.0, 2436789.0, 2436805.0, 2436820.0, 2436836.0,
2436851.0, 2436866.0, 2436881.0, 2436896.0, 2436911.0, 2436925.0,
2436940.0, 2436955.0, 2436970.0, 2436984.0, 2436999.0, 2437014.0,
2437030.0, 2437045.0, 2437060.0, 2437076.0, 2437092.0, 2437107.0,
2437123.0, 2437139.0, 2437154.0, 2437170.0, 2437185.0, 2437201.0,
2437216.0, 2437231.0, 2437246.0, 2437261.0, 2437276.0, 2437291.0,
2437305.0, 2437320.0, 2437335.0, 2437350.0, 2437365.0, 2437380.0,
2437395.0, 2437410.0, 2437426.0, 2437441.0, 2437457.0, 2437472.0,
2437488.0, 2437504.0, 2437520.0, 2437535.0, 2437551.0, 2437566.0,
2437581.0, 2437596.0, 2437611.0, 2437626.0, 2437641.0, 2437656.0,
2437671.0, 2437685.0, 2437700.0, 2437715.0, 2437730.0, 2437745.0,
2437760.0, 2437775.0, 2437791.0, 2437806.0, 2437822.0, 2437838.0,
2437853.0, 2437869.0, 2437885.0, 2437900.0, 2437916.0, 2437931.0,
2437947.0, 2437962.0, 2437977.0, 2437992.0, 2438006.0, 2438021.0,
2438036.0, 2438051.0, 2438065.0, 2438080.0, 2438095.0, 2438110.0,
2438125.0, 2438141.0, 2438156.0, 2438172.0, 2438187.0, 2438203.0,
2438219.0, 2438234.0, 2438250.0, 2438266.0, 2438281.0, 2438297.0,
2438312.0, 2438327.0, 2438342.0, 2438357.0, 2438372.0, 2438386.0,
2438401.0, 2438416.0, 2438431.0, 2438445.0, 2438460.0, 2438475.0,
2438491.0, 2438506.0, 2438521.0, 2438537.0, 2438553.0, 2438568.0,
2438584.0, 2438600.0, 2438615.0, 2438631.0, 2438646.0, 2438662.0,
2438677.0, 2438692.0, 2438707.0, 2438722.0, 2438737.0, 2438752.0,
2438766.0, 2438781.0, 2438796.0, 2438811.0, 2438826.0, 2438841.0,
2438856.0, 2438871.0, 2438887.0, 2438902.0, 2438918.0, 2438933.0,
2438949.0, 2438965.0, 2438981.0, 2438996.0, 2439012.0, 2439027.0,
2439042.0, 2439057.0, 2439072.0, 2439087.0, 2439102.0, 2439117.0,
2439132.0, 2439146.0, 2439161.0, 2439176.0, 2439191.0, 2439206.0,
2439221.0, 2439236.0, 2439252.0, 2439267.0, 2439283.0, 2439299.0,
2439314.0, 2439330.0, 2439346.0, 2439361.0, 2439377.0, 2439392.0,
2439408.0, 2439423.0, 2439438.0, 2439453.0, 2439467.0, 2439482.0,
2439497.0, 2439512.0, 2439526.0, 2439541.0, 2439556.0, 2439571.0,
2439586.0, 2439602.0, 2439617.0, 2439633.0, 2439648.0, 2439664.0,
2439680.0, 2439695.0, 2439711.0, 2439727.0, 2439742.0, 2439758.0,
2439773.0, 2439788.0, 2439803.0, 2439818.0, 2439833.0, 2439847.0,
2439862.0, 2439877.0, 2439892.0, 2439906.0, 2439921.0, 2439936.0,
2439952.0, 2439967.0, 2439982.0, 2439998.0, 2440013.0, 2440029.0,
2440045.0, 2440061.0, 2440076.0, 2440092.0, 2440107.0, 2440123.0,
2440138.0, 2440153.0, 2440168.0, 2440183.0, 2440198.0, 2440213.0,
2440227.0, 2440242.0, 2440257.0, 2440272.0, 2440287.0, 2440302.0,
2440317.0, 2440332.0, 2440348.0, 2440363.0, 2440379.0, 2440394.0,
2440410.0, 2440426.0, 2440442.0, 2440457.0, 2440473.0, 2440488.0,
2440503.0, 2440518.0, 2440533.0, 2440548.0, 2440563.0, 2440578.0,
2440593.0, 2440607.0, 2440622.0, 2440637.0, 2440652.0, 2440667.0,
2440682.0, 2440697.0, 2440713.0, 2440728.0, 2440744.0, 2440760.0,
2440775.0, 2440791.0, 2440807.0, 2440822.0, 2440838.0, 2440853.0,
2440869.0, 2440884.0, 2440899.0, 2440914.0, 2440928.0, 2440943.0,
2440958.0, 2440973.0, 2440987.0, 2441002.0, 2441017.0, 2441032.0,
2441047.0, 2441063.0, 2441078.0, 2441094.0, 2441109.0, 2441125.0,
2441141.0, 2441156.0, 2441172.0, 2441188.0, 2441203.0, 2441219.0,
2441234.0, 2441249.0, 2441264.0, 2441279.0, 2441294.0, 2441308.0,
2441323.0, 2441338.0, 2441353.0, 2441367.0, 2441382.0, 2441397.0,
2441413.0, 2441428.0, 2441443.0, 2441459.0, 2441474.0, 2441490.0,
2441506.0, 2441522.0, 2441537.0, 2441553.0, 2441568.0, 2441584.0,
2441599.0, 2441614.0, 2441629.0, 2441644.0, 2441659.0, 2441674.0,
2441688.0, 2441703.0, 2441718.0, 2441733.0, 2441748.0, 2441763.0,
2441778.0, 2441793.0, 2441808.0, 2441824.0, 2441840.0, 2441855.0,
2441871.0, 2441887.0, 2441903.0, 2441918.0, 2441934.0, 2441949.0,
2441964.0, 2441979.0, 2441994.0, 2442009.0, 2442024.0, 2442039.0,
2442054.0, 2442068.0, 2442083.0, 2442098.0, 2442113.0, 2442128.0,
2442143.0, 2442158.0, 2442174.0, 2442189.0, 2442205.0, 2442221.0,
2442236.0, 2442252.0, 2442268.0, 2442283.0, 2442299.0, 2442314.0,
2442330.0, 2442345.0, 2442360.0, 2442375.0, 2442389.0, 2442404.0,
2442419.0, 2442434.0, 2442448.0, 2442463.0, 2442478.0, 2442493.0,
2442508.0, 2442524.0, 2442539.0, 2442555.0, 2442570.0, 2442586.0,
2442602.0, 2442617.0, 2442633.0, 2442649.0, 2442664.0, 2442679.0,
2442695.0, 2442710.0, 2442725.0, 2442740.0, 2442755.0, 2442769.0,
2442784.0, 2442799.0, 2442814.0, 2442828.0, 2442843.0, 2442858.0,
2442873.0, 2442889.0, 2442904.0, 2442920.0, 2442935.0, 2442951.0,
2442967.0, 2442983.0, 2442998.0, 2443014.0, 2443029.0, 2443045.0,
2443060.0, 2443075.0, 2443090.0, 2443105.0, 2443120.0, 2443135.0,
2443149.0, 2443164.0, 2443179.0, 2443194.0, 2443209.0, 2443224.0,
2443239.0, 2443254.0, 2443269.0, 2443285.0, 2443301.0, 2443316.0,
2443332.0, 2443348.0, 2443363.0, 2443379.0, 2443395.0, 2443410.0,
2443425.0, 2443440.0, 2443455.0, 2443470.0, 2443485.0, 2443500.0,
2443515.0, 2443529.0, 2443544.0, 2443559.0, 2443574.0, 2443589.0,
2443604.0, 2443619.0, 2443635.0, 2443650.0, 2443666.0, 2443682.0,
2443697.0, 2443713.0, 2443729.0, 2443744.0, 2443760.0, 2443775.0,
2443790.0, 2443806.0, 2443821.0, 2443836.0, 2443850.0, 2443865.0,
2443880.0, 2443895.0, 2443909.0, 2443924.0, 2443939.0, 2443954.0,
2443969.0, 2443985.0, 2444000.0, 2444015.0, 2444031.0, 2444047.0,
2444063.0, 2444078.0, 2444094.0, 2444110.0, 2444125.0, 2444140.0,
2444156.0, 2444171.0, 2444186.0, 2444201.0, 2444216.0, 2444230.0,
2444245.0, 2444260.0, 2444275.0, 2444289.0, 2444304.0, 2444319.0,
2444334.0, 2444350.0, 2444365.0, 2444381.0, 2444396.0, 2444412.0,
2444428.0, 2444444.0, 2444459.0, 2444475.0, 2444490.0, 2444506.0,
2444521.0, 2444536.0, 2444551.0, 2444566.0, 2444581.0, 2444596.0,
2444610.0, 2444625.0, 2444640.0, 2444655.0, 2444670.0, 2444685.0,
2444700.0, 2444715.0, 2444730.0, 2444746.0, 2444762.0, 2444777.0,
2444793.0, 2444809.0, 2444824.0, 2444840.0, 2444856.0, 2444871.0,
2444886.0, 2444901.0, 2444916.0, 2444931.0, 2444946.0, 2444961.0,
2444976.0, 2444990.0, 2445005.0, 2445020.0, 2445035.0, 2445050.0,
2445065.0, 2445080.0, 2445096.0, 2445111.0, 2445127.0, 2445143.0,
2445158.0, 2445174.0, 2445190.0, 2445205.0, 2445221.0, 2445236.0,
2445251.0, 2445267.0, 2445282.0, 2445296.0, 2445311.0, 2445326.0,
2445341.0, 2445355.0, 2445370.0, 2445385.0, 2445400.0, 2445415.0,
2445430.0, 2445445.0, 2445461.0, 2445476.0, 2445492.0, 2445508.0,
2445524.0, 2445539.0, 2445555.0, 2445571.0, 2445586.0, 2445601.0,
2445617.0, 2445632.0, 2445647.0, 2445662.0, 2445677.0, 2445691.0,
2445706.0, 2445721.0, 2445735.0, 2445750.0, 2445765.0, 2445780.0,
2445795.0, 2445811.0, 2445826.0, 2445842.0, 2445857.0, 2445873.0,
2445889.0, 2445904.0, 2445920.0, 2445936.0, 2445951.0, 2445967.0,
2445982.0, 2445997.0, 2446012.0, 2446027.0, 2446042.0, 2446057.0,
2446071.0, 2446086.0, 2446101.0, 2446116.0, 2446130.0, 2446146.0,
2446161.0, 2446176.0, 2446191.0, 2446207.0, 2446223.0, 2446238.0,
2446254.0, 2446270.0, 2446285.0, 2446301.0, 2446317.0, 2446332.0,
2446347.0, 2446362.0, 2446377.0, 2446392.0, 2446407.0, 2446422.0,
2446436.0, 2446451.0, 2446466.0, 2446481.0, 2446496.0, 2446511.0,
2446526.0, 2446541.0, 2446557.0, 2446572.0, 2446588.0, 2446604.0,
2446619.0, 2446635.0, 2446651.0, 2446666.0, 2446682.0, 2446697.0,
2446712.0, 2446728.0, 2446743.0, 2446757.0, 2446772.0, 2446787.0,
2446802.0, 2446816.0, 2446831.0, 2446846.0, 2446861.0, 2446876.0,
2446891.0, 2446906.0, 2446922.0, 2446937.0, 2446953.0, 2446969.0,
2446984.0, 2447000.0, 2447016.0, 2447032.0, 2447047.0, 2447062.0,
2447078.0, 2447093.0, 2447108.0, 2447123.0, 2447137.0, 2447152.0,
2447167.0, 2447182.0, 2447196.0, 2447211.0, 2447226.0, 2447241.0,
2447256.0, 2447272.0, 2447287.0, 2447303.0, 2447318.0, 2447334.0,
2447350.0, 2447365.0, 2447381.0, 2447397.0, 2447412.0, 2447428.0,
2447443.0, 2447458.0, 2447473.0, 2447488.0, 2447503.0, 2447517.0,
2447532.0, 2447547.0, 2447562.0, 2447577.0, 2447591.0, 2447606.0,
2447622.0, 2447637.0, 2447652.0, 2447668.0, 2447684.0, 2447699.0,
2447715.0, 2447731.0, 2447746.0, 2447762.0, 2447777.0, 2447793.0,
2447808.0, 2447823.0, 2447838.0, 2447853.0, 2447868.0, 2447883.0,
2447897.0, 2447912.0, 2447927.0, 2447942.0, 2447957.0, 2447972.0,
2447987.0, 2448002.0, 2448018.0, 2448033.0, 2448049.0, 2448064.0,
2448080.0, 2448096.0, 2448112.0, 2448127.0, 2448143.0, 2448158.0,
2448173.0, 2448189.0, 2448204.0, 2448218.0, 2448233.0, 2448248.0,
2448263.0, 2448277.0, 2448292.0, 2448307.0, 2448322.0, 2448337.0,
2448352.0, 2448367.0, 2448383.0, 2448398.0, 2448414.0, 2448430.0,
2448445.0, 2448461.0, 2448477.0, 2448492.0, 2448508.0, 2448523.0,
2448539.0, 2448554.0, 2448569.0, 2448584.0, 2448598.0, 2448613.0,
2448628.0, 2448643.0, 2448657.0, 2448672.0, 2448687.0, 2448702.0,
2448717.0, 2448733.0, 2448748.0, 2448764.0, 2448779.0, 2448795.0,
2448811.0, 2448826.0, 2448842.0, 2448858.0, 2448873.0, 2448889.0,
2448904.0, 2448919.0, 2448934.0, 2448949.0, 2448964.0, 2448978.0,
2448993.0, 2449008.0, 2449023.0, 2449037.0, 2449052.0, 2449067.0,
2449083.0, 2449098.0, 2449113.0, 2449129.0, 2449145.0, 2449160.0,
2449176.0, 2449192.0, 2449207.0, 2449223.0, 2449238.0, 2449254.0,
2449269.0, 2449284.0, 2449299.0, 2449314.0, 2449329.0, 2449344.0,
2449358.0, 2449373.0, 2449388.0, 2449403.0, 2449418.0, 2449433.0,
2449448.0, 2449463.0, 2449479.0, 2449494.0, 2449510.0, 2449525.0,
2449541.0, 2449557.0, 2449573.0, 2449588.0, 2449604.0, 2449619.0,
2449634.0, 2449649.0, 2449664.0, 2449679.0, 2449694.0, 2449709.0,
2449724.0, 2449738.0, 2449753.0, 2449768.0, 2449783.0, 2449798.0,
2449813.0, 2449828.0, 2449844.0, 2449859.0, 2449875.0, 2449891.0,
2449906.0, 2449922.0, 2449938.0, 2449953.0, 2449969.0, 2449984.0,
2450000.0, 2450015.0, 2450030.0, 2450045.0, 2450059.0, 2450074.0,
2450089.0, 2450104.0, 2450118.0, 2450133.0, 2450148.0, 2450163.0,
2450178.0, 2450194.0, 2450209.0, 2450225.0, 2450240.0, 2450256.0,
2450272.0, 2450287.0, 2450303.0, 2450319.0, 2450334.0, 2450350.0,
2450365.0, 2450380.0, 2450395.0, 2450410.0, 2450425.0, 2450439.0,
2450454.0, 2450469.0, 2450484.0, 2450498.0, 2450513.0, 2450528.0,
2450544.0, 2450559.0, 2450574.0, 2450590.0, 2450605.0, 2450621.0,
2450637.0, 2450653.0, 2450668.0, 2450684.0, 2450699.0, 2450715.0,
2450730.0, 2450745.0, 2450760.0, 2450775.0, 2450790.0, 2450805.0,
2450819.0, 2450834.0, 2450849.0, 2450864.0, 2450879.0, 2450894.0,
2450909.0, 2450924.0, 2450940.0, 2450955.0, 2450971.0, 2450986.0,
2451002.0, 2451018.0, 2451034.0, 2451049.0, 2451065.0, 2451080.0,
2451095.0, 2451110.0, 2451125.0, 2451140.0, 2451155.0, 2451170.0,
2451185.0, 2451199.0, 2451214.0, 2451229.0, 2451244.0, 2451259.0,
2451274.0, 2451289.0, 2451305.0, 2451320.0, 2451336.0, 2451352.0,
2451367.0, 2451383.0, 2451399.0, 2451414.0, 2451430.0, 2451445.0,
2451461.0, 2451476.0, 2451491.0, 2451506.0, 2451520.0, 2451535.0,
2451550.0, 2451565.0, 2451579.0, 2451594.0, 2451609.0, 2451624.0,
2451639.0, 2451655.0, 2451670.0, 2451686.0, 2451701.0, 2451717.0,
2451733.0, 2451748.0, 2451764.0, 2451780.0, 2451795.0, 2451811.0,
2451826.0, 2451841.0, 2451856.0, 2451871.0, 2451886.0, 2451900.0,
2451915.0, 2451930.0, 2451945.0, 2451959.0, 2451974.0, 2451989.0,
2452005.0, 2452020.0, 2452035.0, 2452051.0, 2452066.0, 2452082.0,
2452098.0, 2452114.0, 2452129.0, 2452145.0, 2452160.0, 2452176.0,
2452191.0, 2452206.0, 2452221.0, 2452236.0, 2452251.0, 2452266.0,
2452280.0, 2452295.0, 2452310.0, 2452325.0, 2452340.0, 2452355.0,
2452370.0, 2452385.0, 2452401.0, 2452416.0, 2452432.0, 2452447.0,
2452463.0, 2452479.0, 2452495.0, 2452510.0, 2452526.0, 2452541.0,
2452556.0, 2452571.0, 2452586.0, 2452601.0, 2452616.0, 2452631.0,
2452646.0, 2452660.0, 2452675.0, 2452690.0, 2452705.0, 2452720.0,
2452735.0, 2452750.0, 2452766.0, 2452781.0, 2452797.0, 2452813.0,
2452828.0, 2452844.0, 2452860.0, 2452875.0, 2452891.0, 2452906.0,
2452922.0, 2452937.0, 2452952.0, 2452967.0, 2452981.0, 2452996.0,
2453011.0, 2453026.0, 2453040.0, 2453055.0, 2453070.0, 2453085.0,
2453100.0, 2453116.0, 2453131.0, 2453147.0, 2453162.0, 2453178.0,
2453194.0, 2453209.0, 2453225.0, 2453241.0, 2453256.0, 2453272.0,
2453287.0, 2453302.0, 2453317.0, 2453332.0, 2453347.0, 2453361.0,
2453376.0, 2453391.0, 2453406.0, 2453420.0, 2453435.0, 2453450.0,
2453466.0, 2453481.0, 2453496.0, 2453512.0, 2453527.0, 2453543.0,
2453559.0, 2453575.0, 2453590.0, 2453606.0, 2453621.0, 2453637.0,
2453652.0, 2453667.0, 2453682.0, 2453697.0, 2453712.0, 2453727.0,
2453741.0, 2453756.0, 2453771.0, 2453786.0, 2453801.0, 2453816.0,
2453831.0, 2453846.0, 2453861.0, 2453877.0, 2453893.0, 2453908.0,
2453924.0, 2453940.0, 2453955.0, 2453971.0, 2453987.0, 2454002.0,
2454017.0, 2454032.0, 2454047.0, 2454062.0, 2454077.0, 2454092.0,
2454107.0, 2454121.0, 2454136.0, 2454151.0, 2454166.0, 2454181.0,
2454196.0, 2454211.0, 2454227.0, 2454242.0, 2454258.0, 2454274.0,
2454289.0, 2454305.0, 2454321.0, 2454336.0, 2454352.0, 2454367.0,
2454383.0, 2454398.0, 2454413.0, 2454428.0, 2454442.0, 2454457.0,
2454472.0, 2454487.0, 2454501.0, 2454516.0, 2454531.0, 2454546.0,
2454561.0, 2454577.0, 2454592.0, 2454608.0, 2454623.0, 2454639.0,
2454655.0, 2454670.0, 2454686.0, 2454702.0, 2454717.0, 2454732.0,
2454748.0, 2454763.0, 2454778.0, 2454793.0, 2454808.0, 2454822.0,
2454837.0, 2454852.0, 2454867.0, 2454881.0, 2454896.0, 2454911.0,
2454926.0, 2454942.0, 2454957.0, 2454973.0, 2454988.0, 2455004.0,
2455020.0, 2455036.0, 2455051.0, 2455067.0, 2455082.0, 2455098.0,
2455113.0, 2455128.0, 2455143.0, 2455158.0, 2455173.0, 2455188.0,
2455202.0, 2455217.0, 2455232.0, 2455247.0, 2455262.0, 2455277.0,
2455292.0, 2455307.0, 2455322.0, 2455338.0, 2455354.0, 2455369.0,
2455385.0, 2455401.0, 2455416.0, 2455432.0, 2455448.0, 2455463.0,
2455478.0, 2455493.0, 2455508.0, 2455523.0, 2455538.0, 2455553.0,
2455568.0, 2455582.0, 2455597.0, 2455612.0, 2455627.0, 2455642.0,
2455657.0, 2455672.0, 2455688.0, 2455703.0, 2455719.0, 2455735.0,
2455750.0, 2455766.0, 2455782.0, 2455797.0, 2455813.0, 2455828.0,
2455843.0, 2455859.0, 2455874.0, 2455889.0, 2455903.0, 2455918.0,
2455933.0, 2455948.0, 2455962.0, 2455977.0, 2455992.0, 2456007.0,
2456022.0, 2456038.0, 2456053.0, 2456068.0, 2456084.0, 2456100.0,
2456116.0, 2456131.0, 2456147.0, 2456163.0, 2456178.0, 2456193.0,
2456209.0, 2456224.0, 2456239.0, 2456254.0, 2456269.0, 2456283.0,
2456298.0, 2456313.0, 2456328.0, 2456342.0, 2456357.0, 2456372.0,
2456387.0, 2456403.0, 2456418.0, 2456434.0, 2456449.0, 2456465.0,
2456481.0, 2456496.0, 2456512.0, 2456528.0, 2456543.0, 2456559.0,
2456574.0, 2456589.0, 2456604.0, 2456619.0, 2456634.0, 2456649.0,
2456663.0, 2456678.0, 2456693.0, 2456708.0, 2456723.0, 2456738.0,
2456753.0, 2456768.0, 2456783.0, 2456799.0, 2456815.0, 2456830.0,
2456846.0, 2456862.0, 2456877.0, 2456893.0, 2456909.0, 2456924.0,
2456939.0, 2456954.0, 2456969.0, 2456984.0, 2456999.0, 2457014.0,
2457029.0, 2457043.0, 2457058.0, 2457073.0, 2457088.0, 2457103.0,
2457118.0, 2457133.0, 2457149.0, 2457164.0, 2457180.0, 2457196.0,
2457211.0, 2457227.0, 2457243.0, 2457258.0, 2457274.0, 2457289.0,
2457304.0, 2457320.0, 2457335.0, 2457349.0, 2457364.0, 2457379.0,
2457394.0, 2457408.0, 2457423.0, 2457438.0, 2457453.0, 2457468.0,
2457483.0, 2457498.0, 2457514.0, 2457529.0, 2457545.0, 2457561.0,
2457577.0, 2457592.0, 2457608.0, 2457624.0, 2457639.0, 2457654.0,
2457670.0, 2457685.0, 2457700.0, 2457715.0, 2457730.0, 2457744.0,
2457759.0, 2457774.0, 2457788.0, 2457803.0, 2457818.0, 2457833.0,
2457848.0, 2457864.0, 2457879.0, 2457895.0, 2457910.0, 2457926.0,
2457942.0, 2457957.0, 2457973.0, 2457989.0, 2458004.0, 2458020.0,
2458035.0, 2458050.0, 2458065.0, 2458080.0, 2458095.0, 2458110.0,
2458124.0, 2458139.0, 2458154.0, 2458169.0, 2458183.0, 2458199.0,
2458214.0, 2458229.0, 2458244.0, 2458260.0, 2458276.0, 2458291.0,
2458307.0, 2458323.0, 2458338.0, 2458354.0, 2458370.0, 2458385.0,
2458400.0, 2458415.0, 2458430.0, 2458445.0, 2458460.0, 2458475.0,
2458489.0, 2458504.0, 2458519.0, 2458534.0, 2458549.0, 2458564.0,
2458579.0, 2458594.0, 2458610.0, 2458625.0, 2458641.0, 2458656.0,
2458672.0, 2458688.0, 2458704.0, 2458719.0, 2458735.0, 2458750.0,
2458765.0, 2458781.0, 2458796.0, 2458810.0, 2458825.0, 2458840.0,
2458855.0, 2458869.0, 2458884.0, 2458899.0, 2458914.0, 2458929.0,
2458944.0, 2458959.0, 2458975.0, 2458990.0, 2459006.0, 2459022.0,
2459037.0, 2459053.0, 2459069.0, 2459084.0, 2459100.0, 2459115.0,
2459131.0, 2459146.0, 2459161.0, 2459176.0, 2459191.0, 2459205.0,
2459220.0, 2459235.0, 2459249.0, 2459264.0, 2459279.0, 2459294.0,
2459309.0, 2459325.0, 2459340.0, 2459356.0, 2459371.0, 2459387.0,
2459403.0, 2459418.0, 2459434.0, 2459450.0, 2459465.0, 2459481.0,
2459496.0, 2459511.0, 2459526.0, 2459541.0, 2459556.0, 2459571.0,
2459585.0, 2459600.0, 2459615.0, 2459630.0, 2459644.0, 2459659.0,
2459675.0, 2459690.0, 2459705.0, 2459721.0, 2459737.0, 2459752.0,
2459768.0, 2459784.0, 2459799.0, 2459815.0, 2459830.0, 2459846.0,
2459861.0, 2459876.0, 2459891.0, 2459906.0, 2459921.0, 2459936.0,
2459950.0, 2459965.0, 2459980.0, 2459995.0, 2460010.0, 2460025.0,
2460040.0, 2460055.0, 2460071.0, 2460086.0, 2460102.0, 2460117.0,
2460133.0, 2460149.0, 2460165.0, 2460180.0, 2460196.0, 2460211.0,
2460226.0, 2460242.0, 2460257.0, 2460271.0, 2460286.0, 2460301.0,
2460316.0, 2460330.0, 2460345.0, 2460360.0, 2460375.0, 2460390.0,
2460405.0, 2460420.0, 2460436.0, 2460451.0, 2460467.0, 2460483.0,
2460498.0, 2460514.0, 2460530.0, 2460545.0, 2460561.0, 2460576.0,
2460592.0, 2460607.0, 2460622.0, 2460637.0, 2460651.0, 2460666.0,
2460681.0, 2460696.0, 2460710.0, 2460725.0, 2460740.0, 2460755.0,
2460770.0, 2460786.0, 2460801.0, 2460817.0, 2460832.0, 2460848.0,
2460864.0, 2460879.0, 2460895.0, 2460911.0, 2460926.0, 2460942.0,
2460957.0, 2460972.0, 2460987.0, 2461002.0, 2461017.0, 2461031.0,
2461046.0, 2461061.0, 2461076.0, 2461090.0, 2461105.0, 2461120.0,
2461136.0, 2461151.0, 2461166.0, 2461182.0, 2461197.0, 2461213.0,
2461229.0, 2461245.0, 2461260.0, 2461276.0, 2461291.0, 2461307.0,
2461322.0, 2461337.0, 2461352.0, 2461367.0, 2461382.0, 2461397.0,
2461411.0, 2461426.0, 2461441.0, 2461456.0, 2461471.0, 2461486.0,
2461501.0, 2461516.0, 2461532.0, 2461547.0, 2461563.0, 2461578.0,
2461594.0, 2461610.0, 2461626.0, 2461641.0, 2461657.0, 2461672.0,
2461687.0, 2461702.0, 2461717.0, 2461732.0, 2461747.0, 2461762.0,
2461777.0, 2461791.0, 2461806.0, 2461821.0, 2461836.0, 2461851.0,
2461866.0, 2461881.0, 2461897.0, 2461912.0, 2461928.0, 2461944.0,
2461959.0, 2461975.0, 2461991.0, 2462006.0, 2462022.0, 2462037.0,
2462053.0, 2462068.0, 2462083.0, 2462098.0, 2462112.0, 2462127.0,
2462142.0, 2462157.0, 2462171.0, 2462186.0, 2462201.0, 2462216.0,
2462231.0, 2462247.0, 2462262.0, 2462278.0, 2462293.0, 2462309.0,
2462325.0, 2462340.0, 2462356.0, 2462372.0, 2462387.0, 2462403.0,
2462418.0, 2462433.0, 2462448.0, 2462463.0, 2462478.0, 2462492.0,
2462507.0, 2462522.0, 2462537.0, 2462551.0, 2462566.0, 2462581.0,
2462597.0, 2462612.0, 2462627.0, 2462643.0, 2462658.0, 2462674.0,
2462690.0, 2462706.0, 2462721.0, 2462737.0, 2462752.0, 2462768.0,
2462783.0, 2462798.0, 2462813.0, 2462828.0, 2462843.0, 2462858.0,
2462872.0, 2462887.0, 2462902.0, 2462917.0, 2462932.0, 2462947.0,
2462962.0, 2462977.0, 2462993.0, 2463008.0, 2463024.0, 2463039.0,
2463055.0, 2463071.0, 2463087.0, 2463102.0, 2463118.0, 2463133.0,
2463148.0, 2463163.0, 2463178.0, 2463193.0, 2463208.0, 2463223.0,
2463238.0, 2463252.0, 2463267.0, 2463282.0, 2463297.0, 2463312.0,
2463327.0, 2463342.0, 2463358.0, 2463373.0, 2463389.0, 2463405.0,
2463420.0, 2463436.0, 2463452.0, 2463467.0, 2463483.0, 2463498.0,
2463514.0, 2463529.0, 2463544.0, 2463559.0, 2463573.0, 2463588.0,
2463603.0, 2463618.0, 2463632.0, 2463647.0, 2463662.0, 2463677.0,
2463692.0, 2463708.0, 2463723.0, 2463739.0, 2463754.0, 2463770.0,
2463786.0, 2463801.0, 2463817.0, 2463833.0, 2463848.0, 2463864.0,
2463879.0, 2463894.0, 2463909.0, 2463924.0, 2463939.0, 2463953.0,
2463968.0, 2463983.0, 2463998.0, 2464012.0, 2464027.0, 2464042.0,
2464058.0, 2464073.0, 2464088.0, 2464104.0, 2464119.0, 2464135.0,
2464151.0, 2464167.0, 2464182.0, 2464198.0, 2464213.0, 2464229.0,
2464244.0, 2464259.0, 2464274.0, 2464289.0, 2464304.0, 2464319.0,
2464333.0, 2464348.0, 2464363.0, 2464378.0, 2464393.0, 2464408.0,
2464423.0, 2464438.0, 2464453.0, 2464469.0, 2464485.0, 2464500.0,
2464516.0, 2464532.0, 2464547.0, 2464563.0, 2464579.0, 2464594.0,
2464609.0, 2464624.0, 2464639.0, 2464654.0, 2464669.0, 2464684.0,
2464699.0, 2464713.0, 2464728.0, 2464743.0, 2464758.0, 2464773.0,
2464788.0, 2464803.0, 2464819.0, 2464834.0, 2464850.0, 2464866.0,
2464881.0, 2464897.0, 2464913.0, 2464928.0, 2464944.0, 2464959.0,
2464975.0, 2464990.0, 2465005.0, 2465020.0, 2465034.0, 2465049.0,
2465064.0, 2465079.0, 2465093.0, 2465108.0, 2465123.0, 2465138.0,
2465153.0, 2465169.0, 2465184.0, 2465200.0, 2465215.0, 2465231.0,
2465247.0, 2465262.0, 2465278.0, 2465294.0, 2465309.0, 2465325.0,
2465340.0, 2465355.0, 2465370.0, 2465385.0, 2465400.0, 2465414.0,
2465429.0, 2465444.0, 2465459.0, 2465473.0, 2465488.0, 2465503.0,
2465519.0, 2465534.0, 2465549.0, 2465565.0, 2465580.0, 2465596.0,
2465612.0, 2465628.0, 2465643.0, 2465659.0, 2465674.0, 2465690.0,
2465705.0, 2465720.0, 2465735.0, 2465750.0, 2465765.0, 2465780.0,
2465794.0, 2465809.0, 2465824.0, 2465839.0, 2465854.0, 2465869.0,
2465884.0, 2465899.0, 2465914.0, 2465930.0, 2465946.0, 2465961.0,
2465977.0, 2465993.0, 2466008.0, 2466024.0, 2466040.0, 2466055.0,
2466070.0, 2466085.0, 2466100.0, 2466115.0, 2466130.0, 2466145.0,
2466160.0, 2466174.0, 2466189.0, 2466204.0, 2466219.0, 2466234.0,
2466249.0, 2466264.0, 2466280.0, 2466295.0, 2466311.0, 2466327.0,
2466342.0, 2466358.0, 2466374.0, 2466389.0, 2466405.0, 2466420.0,
2466436.0, 2466451.0, 2466466.0, 2466481.0, 2466495.0, 2466510.0,
2466525.0, 2466540.0, 2466554.0, 2466569.0, 2466584.0, 2466599.0,
2466614.0, 2466630.0, 2466645.0, 2466660.0, 2466676.0, 2466692.0,
2466708.0, 2466723.0, 2466739.0, 2466755.0, 2466770.0, 2466785.0,
2466801.0, 2466816.0, 2466831.0, 2466846.0, 2466861.0, 2466875.0,
2466890.0, 2466905.0, 2466920.0, 2466934.0, 2466949.0, 2466964.0,
2466979.0, 2466995.0, 2467010.0, 2467026.0, 2467041.0, 2467057.0,
2467073.0, 2467089.0, 2467104.0, 2467120.0, 2467135.0, 2467151.0,
2467166.0, 2467181.0, 2467196.0, 2467211.0, 2467226.0, 2467241.0,
2467255.0, 2467270.0, 2467285.0, 2467300.0, 2467315.0, 2467330.0,
2467345.0, 2467360.0, 2467375.0, 2467391.0, 2467407.0, 2467422.0,
2467438.0, 2467454.0, 2467469.0, 2467485.0, 2467501.0, 2467516.0,
2467531.0, 2467546.0, 2467561.0, 2467576.0, 2467591.0, 2467606.0,
2467621.0, 2467635.0, 2467650.0, 2467665.0, 2467680.0, 2467695.0,
2467710.0, 2467725.0, 2467741.0, 2467756.0, 2467772.0, 2467788.0,
2467803.0, 2467819.0, 2467835.0, 2467850.0, 2467866.0, 2467881.0,
2467896.0, 2467912.0, 2467927.0, 2467942.0, 2467956.0, 2467971.0,
2467986.0, 2468001.0, 2468015.0, 2468030.0, 2468045.0, 2468060.0,
2468075.0, 2468090.0, 2468106.0, 2468121.0, 2468137.0, 2468153.0,
2468169.0, 2468184.0, 2468200.0, 2468216.0, 2468231.0, 2468246.0,
2468262.0, 2468277.0, 2468292.0, 2468307.0, 2468322.0, 2468336.0,
2468351.0, 2468366.0, 2468381.0, 2468395.0, 2468410.0, 2468425.0,
2468440.0, 2468456.0, 2468471.0, 2468487.0, 2468502.0, 2468518.0,
2468534.0, 2468549.0, 2468565.0, 2468581.0, 2468596.0, 2468612.0,
2468627.0, 2468642.0, 2468657.0, 2468672.0, 2468687.0, 2468702.0,
2468716.0, 2468731.0, 2468746.0, 2468761.0, 2468776.0, 2468791.0,
2468806.0, 2468821.0, 2468836.0, 2468852.0, 2468868.0, 2468883.0,
2468899.0, 2468915.0, 2468930.0, 2468946.0, 2468962.0, 2468977.0,
2468992.0, 2469007.0, 2469022.0, 2469037.0, 2469052.0, 2469067.0,
2469082.0, 2469096.0, 2469111.0, 2469126.0, 2469141.0, 2469156.0,
2469171.0, 2469186.0, 2469202.0, 2469217.0, 2469233.0, 2469248.0,
2469264.0, 2469280.0, 2469296.0, 2469311.0, 2469327.0, 2469342.0,
2469357.0, 2469373.0, 2469388.0, 2469402.0, 2469417.0, 2469432.0,
2469447.0, 2469461.0, 2469476.0, 2469491.0, 2469506.0, 2469521.0,
2469536.0, 2469551.0, 2469567.0, 2469582.0, 2469598.0, 2469614.0,
2469629.0, 2469645.0, 2469661.0, 2469676.0, 2469692.0, 2469707.0,
2469723.0, 2469738.0, 2469753.0, 2469768.0, 2469783.0, 2469797.0,
2469812.0, 2469827.0, 2469841.0, 2469856.0, 2469871.0, 2469886.0,
2469901.0, 2469917.0, 2469932.0, 2469948.0, 2469963.0, 2469979.0,
2469995.0, 2470010.0, 2470026.0, 2470042.0, 2470057.0, 2470073.0,
2470088.0, 2470103.0, 2470118.0, 2470133.0, 2470148.0, 2470163.0
};

static double allmoons[] = {
2433271.0, 2433300.0, 2433330.0, 2433359.0, 2433389.0, 2433419.0,
2433448.0, 2433478.0, 2433508.0, 2433537.0, 2433566.0, 2433596.0,
2433625.0, 2433655.0, 2433684.0, 2433714.0, 2433743.0, 2433773.0,
2433803.0, 2433832.0, 2433862.0, 2433891.0, 2433921.0, 2433950.0,
2433980.0, 2434009.0, 2434039.0, 2434068.0, 2434098.0, 2434127.0,
2434157.0, 2434186.0, 2434216.0, 2434245.0, 2434275.0, 2434305.0,
2434334.0, 2434364.0, 2434393.0, 2434423.0, 2434452.0, 2434482.0,
2434511.0, 2434540.0, 2434570.0, 2434600.0, 2434629.0, 2434659.0,
2434689.0, 2434718.0, 2434748.0, 2434777.0, 2434807.0, 2434836.0,
2434866.0, 2434895.0, 2434924.0, 2434954.0, 2434983.0, 2435013.0,
2435043.0, 2435072.0, 2435102.0, 2435132.0, 2435161.0, 2435191.0,
2435220.0, 2435250.0, 2435279.0, 2435308.0, 2435338.0, 2435367.0,
2435397.0, 2435426.0, 2435456.0, 2435486.0, 2435516.0, 2435545.0,
2435575.0, 2435604.0, 2435634.0, 2435663.0, 2435692.0, 2435722.0,
2435751.0, 2435781.0, 2435810.0, 2435840.0, 2435870.0, 2435900.0,
2435929.0, 2435959.0, 2435988.0, 2436018.0, 2436047.0, 2436076.0,
2436106.0, 2436135.0, 2436165.0, 2436194.0, 2436224.0, 2436253.0,
2436283.0, 2436313.0, 2436343.0, 2436372.0, 2436402.0, 2436431.0,
2436460.0, 2436490.0, 2436519.0, 2436549.0, 2436578.0, 2436608.0,
2436637.0, 2436667.0, 2436697.0, 2436726.0, 2436756.0, 2436785.0,
2436815.0, 2436844.0, 2436874.0, 2436903.0, 2436933.0, 2436962.0,
2436992.0, 2437021.0, 2437051.0, 2437080.0, 2437110.0, 2437140.0,
2437169.0, 2437199.0, 2437228.0, 2437258.0, 2437287.0, 2437317.0,
2437346.0, 2437376.0, 2437405.0, 2437435.0, 2437464.0, 2437494.0,
2437523.0, 2437553.0, 2437583.0, 2437612.0, 2437642.0, 2437671.0,
2437701.0, 2437730.0, 2437760.0, 2437789.0, 2437818.0, 2437848.0,
2437877.0, 2437907.0, 2437937.0, 2437966.0, 2437996.0, 2438026.0,
2438055.0, 2438085.0, 2438114.0, 2438144.0, 2438173.0, 2438202.0,
2438232.0, 2438261.0, 2438291.0, 2438320.0, 2438350.0, 2438380.0,
2438410.0, 2438439.0, 2438469.0, 2438498.0, 2438528.0, 2438557.0,
2438586.0, 2438616.0, 2438645.0, 2438675.0, 2438704.0, 2438734.0,
2438764.0, 2438794.0, 2438823.0, 2438853.0, 2438882.0, 2438912.0,
2438941.0, 2438970.0, 2439000.0, 2439029.0, 2439058.0, 2439088.0,
2439118.0, 2439147.0, 2439177.0, 2439207.0, 2439237.0, 2439266.0,
2439296.0, 2439325.0, 2439354.0, 2439384.0, 2439413.0, 2439442.0,
2439472.0, 2439502.0, 2439531.0, 2439561.0, 2439591.0, 2439620.0,
2439650.0, 2439680.0, 2439709.0, 2439738.0, 2439768.0, 2439797.0,
2439827.0, 2439856.0, 2439886.0, 2439915.0, 2439945.0, 2439974.0,
2440004.0, 2440034.0, 2440063.0, 2440093.0, 2440122.0, 2440152.0,
2440181.0, 2440211.0, 2440240.0, 2440270.0, 2440299.0, 2440329.0,
2440358.0, 2440388.0, 2440417.0, 2440447.0, 2440477.0, 2440506.0,
2440536.0, 2440565.0, 2440595.0, 2440624.0, 2440654.0, 2440683.0,
2440712.0, 2440742.0, 2440771.0, 2440801.0, 2440831.0, 2440860.0,
2440890.0, 2440920.0, 2440949.0, 2440979.0, 2441008.0, 2441038.0,
2441067.0, 2441096.0, 2441126.0, 2441155.0, 2441185.0, 2441214.0,
2441244.0, 2441274.0, 2441304.0, 2441333.0, 2441363.0, 2441392.0,
2441422.0, 2441451.0, 2441480.0, 2441510.0, 2441539.0, 2441569.0,
2441598.0, 2441628.0, 2441658.0, 2441687.0, 2441717.0, 2441747.0,
2441776.0, 2441806.0, 2441835.0, 2441864.0, 2441894.0, 2441923.0,
2441952.0, 2441982.0, 2442012.0, 2442041.0, 2442071.0, 2442101.0,
2442131.0, 2442160.0, 2442190.0, 2442219.0, 2442248.0, 2442278.0,
2442307.0, 2442336.0, 2442366.0, 2442396.0, 2442425.0, 2442455.0,
2442485.0, 2442515.0, 2442544.0, 2442574.0, 2442603.0, 2442632.0,
2442662.0, 2442691.0, 2442720.0, 2442750.0, 2442779.0, 2442809.0,
2442839.0, 2442869.0, 2442898.0, 2442928.0, 2442957.0, 2442987.0,
2443016.0, 2443046.0, 2443075.0, 2443104.0, 2443134.0, 2443163.0,
2443193.0, 2443223.0, 2443252.0, 2443282.0, 2443312.0, 2443341.0,
2443371.0, 2443400.0, 2443430.0, 2443459.0, 2443489.0, 2443518.0,
2443547.0, 2443577.0, 2443606.0, 2443636.0, 2443666.0, 2443695.0,
2443725.0, 2443755.0, 2443784.0, 2443814.0, 2443843.0, 2443873.0,
2443902.0, 2443932.0, 2443961.0, 2443990.0, 2444020.0, 2444049.0,
2444079.0, 2444109.0, 2444138.0, 2444168.0, 2444198.0, 2444227.0,
2444257.0, 2444286.0, 2444316.0, 2444345.0, 2444374.0, 2444404.0,
2444433.0, 2444463.0, 2444492.0, 2444522.0, 2444552.0, 2444581.0,
2444611.0, 2444641.0, 2444670.0, 2444700.0, 2444729.0, 2444758.0,
2444788.0, 2444817.0, 2444846.0, 2444876.0, 2444906.0, 2444935.0,
2444965.0, 2444995.0, 2445025.0, 2445054.0, 2445084.0, 2445113.0,
2445142.0, 2445172.0, 2445201.0, 2445230.0, 2445260.0, 2445289.0,
2445319.0, 2445349.0, 2445379.0, 2445409.0, 2445438.0, 2445468.0,
2445497.0, 2445526.0, 2445556.0, 2445585.0, 2445614.0, 2445644.0,
2445673.0, 2445703.0, 2445733.0, 2445763.0, 2445792.0, 2445822.0,
2445852.0, 2445881.0, 2445910.0, 2445940.0, 2445969.0, 2445998.0,
2446028.0, 2446057.0, 2446087.0, 2446117.0, 2446146.0, 2446176.0,
2446206.0, 2446235.0, 2446265.0, 2446294.0, 2446324.0, 2446353.0,
2446382.0, 2446412.0, 2446441.0, 2446471.0, 2446500.0, 2446530.0,
2446560.0, 2446589.0, 2446619.0, 2446649.0, 2446678.0, 2446708.0,
2446737.0, 2446767.0, 2446796.0, 2446825.0, 2446855.0, 2446884.0,
2446914.0, 2446943.0, 2446973.0, 2447003.0, 2447032.0, 2447062.0,
2447092.0, 2447121.0, 2447151.0, 2447180.0, 2447209.0, 2447239.0,
2447268.0, 2447298.0, 2447327.0, 2447357.0, 2447386.0, 2447416.0,
2447446.0, 2447475.0, 2447505.0, 2447535.0, 2447564.0, 2447594.0,
2447623.0, 2447652.0, 2447682.0, 2447711.0, 2447741.0, 2447770.0,
2447800.0, 2447829.0, 2447859.0, 2447889.0, 2447919.0, 2447948.0,
2447978.0, 2448007.0, 2448036.0, 2448066.0, 2448095.0, 2448124.0,
2448154.0, 2448183.0, 2448213.0, 2448243.0, 2448273.0, 2448303.0,
2448332.0, 2448362.0, 2448391.0, 2448420.0, 2448450.0, 2448479.0,
2448508.0, 2448538.0, 2448567.0, 2448597.0, 2448627.0, 2448657.0,
2448686.0, 2448716.0, 2448746.0, 2448775.0, 2448804.0, 2448834.0,
2448863.0, 2448892.0, 2448922.0, 2448951.0, 2448981.0, 2449011.0,
2449040.0, 2449070.0, 2449100.0, 2449129.0, 2449159.0, 2449188.0,
2449218.0, 2449247.0, 2449276.0, 2449306.0, 2449335.0, 2449365.0,
2449394.0, 2449424.0, 2449454.0, 2449484.0, 2449513.0, 2449543.0,
2449572.0, 2449602.0, 2449631.0, 2449660.0, 2449690.0, 2449719.0,
2449749.0, 2449778.0, 2449808.0, 2449838.0, 2449867.0, 2449897.0,
2449926.0, 2449956.0, 2449986.0, 2450015.0, 2450044.0, 2450074.0,
2450103.0, 2450133.0, 2450162.0, 2450192.0, 2450221.0, 2450251.0,
2450281.0, 2450310.0, 2450340.0, 2450369.0, 2450399.0, 2450429.0,
2450458.0, 2450487.0, 2450517.0, 2450546.0, 2450576.0, 2450605.0,
2450635.0, 2450664.0, 2450694.0, 2450724.0, 2450753.0, 2450783.0,
2450813.0, 2450842.0, 2450872.0, 2450901.0, 2450930.0, 2450960.0,
2450989.0, 2451018.0, 2451048.0, 2451078.0, 2451107.0, 2451137.0,
2451167.0, 2451196.0, 2451226.0, 2451256.0, 2451285.0, 2451314.0,
2451344.0, 2451373.0, 2451402.0, 2451432.0, 2451461.0, 2451491.0,
2451521.0, 2451551.0, 2451580.0, 2451610.0, 2451640.0, 2451669.0,
2451698.0, 2451728.0, 2451757.0, 2451786.0, 2451816.0, 2451845.0,
2451875.0, 2451905.0, 2451934.0, 2451964.0, 2451994.0, 2452023.0,
2452053.0, 2452082.0, 2452112.0, 2452141.0, 2452170.0, 2452200.0,
2452229.0, 2452259.0, 2452288.0, 2452318.0, 2452348.0, 2452378.0,
2452407.0, 2452437.0, 2452466.0, 2452496.0, 2452525.0, 2452554.0,
2452584.0, 2452613.0, 2452643.0, 2452672.0, 2452702.0, 2452732.0,
2452761.0, 2452791.0, 2452821.0, 2452850.0, 2452880.0, 2452909.0,
2452938.0, 2452968.0, 2452997.0, 2453027.0, 2453056.0, 2453086.0,
2453115.0, 2453145.0, 2453175.0, 2453204.0, 2453234.0, 2453263.0,
2453293.0, 2453322.0, 2453352.0, 2453381.0, 2453411.0, 2453440.0,
2453470.0, 2453499.0, 2453529.0, 2453558.0, 2453588.0, 2453618.0,
2453647.0, 2453677.0, 2453706.0, 2453736.0, 2453765.0, 2453795.0,
2453824.0, 2453854.0, 2453883.0, 2453913.0, 2453942.0, 2453972.0,
2454001.0, 2454031.0, 2454061.0, 2454090.0, 2454120.0, 2454150.0,
2454179.0, 2454208.0, 2454238.0, 2454267.0, 2454296.0, 2454326.0,
2454355.0, 2454385.0, 2454415.0, 2454445.0, 2454474.0, 2454504.0,
2454534.0, 2454563.0, 2454592.0, 2454622.0, 2454651.0, 2454680.0,
2454710.0, 2454739.0, 2454769.0, 2454799.0, 2454828.0, 2454858.0,
2454888.0, 2454918.0, 2454947.0, 2454976.0, 2455006.0, 2455035.0,
2455064.0, 2455094.0, 2455123.0, 2455153.0, 2455182.0, 2455212.0,
2455242.0, 2455272.0, 2455301.0, 2455331.0, 2455360.0, 2455390.0,
2455419.0, 2455448.0, 2455478.0, 2455507.0, 2455537.0, 2455566.0,
2455596.0, 2455626.0, 2455655.0, 2455685.0, 2455715.0, 2455744.0,
2455774.0, 2455803.0, 2455832.0, 2455862.0, 2455891.0, 2455921.0,
2455950.0, 2455980.0, 2456009.0, 2456039.0, 2456069.0, 2456098.0,
2456128.0, 2456157.0, 2456187.0, 2456216.0, 2456246.0, 2456275.0,
2456305.0, 2456334.0, 2456364.0, 2456393.0, 2456423.0, 2456452.0,
2456482.0, 2456512.0, 2456541.0, 2456571.0, 2456600.0, 2456630.0,
2456659.0, 2456689.0, 2456718.0, 2456748.0, 2456777.0, 2456807.0,
2456836.0, 2456866.0, 2456895.0, 2456925.0, 2456955.0, 2456984.0,
2457014.0, 2457043.0, 2457073.0, 2457102.0, 2457132.0, 2457161.0,
2457190.0, 2457220.0, 2457249.0, 2457279.0, 2457309.0, 2457339.0,
2457368.0, 2457398.0, 2457427.0, 2457457.0, 2457486.0, 2457516.0,
2457545.0, 2457574.0, 2457604.0, 2457633.0, 2457663.0, 2457693.0,
2457722.0, 2457752.0, 2457782.0, 2457811.0, 2457841.0, 2457870.0,
2457900.0, 2457929.0, 2457958.0, 2457988.0, 2458017.0, 2458047.0,
2458076.0, 2458106.0, 2458136.0, 2458166.0, 2458195.0, 2458225.0,
2458254.0, 2458284.0, 2458313.0, 2458342.0, 2458372.0, 2458401.0,
2458431.0, 2458460.0, 2458490.0, 2458520.0, 2458550.0, 2458579.0,
2458609.0, 2458638.0, 2458668.0, 2458697.0, 2458726.0, 2458756.0,
2458785.0, 2458814.0, 2458844.0, 2458874.0, 2458903.0, 2458933.0,
2458963.0, 2458993.0, 2459022.0, 2459052.0, 2459081.0, 2459110.0,
2459140.0, 2459169.0, 2459199.0, 2459228.0, 2459258.0, 2459287.0,
2459317.0, 2459347.0, 2459376.0, 2459406.0, 2459435.0, 2459465.0,
2459494.0, 2459524.0, 2459553.0, 2459583.0, 2459612.0, 2459642.0,
2459671.0, 2459701.0, 2459730.0, 2459760.0, 2459790.0, 2459819.0,
2459849.0, 2459878.0, 2459908.0, 2459937.0, 2459967.0, 2459996.0,
2460026.0, 2460055.0, 2460084.0, 2460114.0, 2460144.0, 2460173.0,
2460203.0, 2460233.0, 2460262.0, 2460292.0, 2460321.0, 2460351.0,
2460380.0, 2460410.0, 2460439.0, 2460468.0, 2460498.0, 2460527.0,
2460557.0, 2460587.0, 2460616.0, 2460646.0, 2460676.0, 2460705.0,
2460735.0, 2460764.0, 2460794.0, 2460823.0, 2460852.0, 2460882.0,
2460911.0, 2460941.0, 2460970.0, 2461000.0, 2461030.0, 2461060.0,
2461089.0, 2461119.0, 2461148.0, 2461178.0, 2461207.0, 2461236.0,
2461266.0, 2461295.0, 2461324.0, 2461354.0, 2461384.0, 2461414.0,
2461443.0, 2461473.0, 2461503.0, 2461532.0, 2461562.0, 2461591.0,
2461620.0, 2461650.0, 2461679.0, 2461708.0, 2461738.0, 2461768.0,
2461797.0, 2461827.0, 2461857.0, 2461887.0, 2461916.0, 2461946.0,
2461975.0, 2462004.0, 2462034.0, 2462063.0, 2462092.0, 2462122.0,
2462152.0, 2462181.0, 2462211.0, 2462241.0, 2462270.0, 2462300.0,
2462329.0, 2462359.0, 2462388.0, 2462418.0, 2462447.0, 2462476.0,
2462506.0, 2462536.0, 2462565.0, 2462595.0, 2462624.0, 2462654.0,
2462684.0, 2462713.0, 2462743.0, 2462772.0, 2462802.0, 2462831.0,
2462861.0, 2462890.0, 2462919.0, 2462949.0, 2462979.0, 2463008.0,
2463038.0, 2463067.0, 2463097.0, 2463127.0, 2463156.0, 2463186.0,
2463215.0, 2463245.0, 2463274.0, 2463304.0, 2463333.0, 2463362.0,
2463392.0, 2463421.0, 2463451.0, 2463481.0, 2463510.0, 2463540.0,
2463570.0, 2463599.0, 2463629.0, 2463658.0, 2463688.0, 2463717.0,
2463746.0, 2463776.0, 2463805.0, 2463835.0, 2463864.0, 2463894.0,
2463924.0, 2463954.0, 2463983.0, 2464013.0, 2464042.0, 2464072.0,
2464101.0, 2464130.0, 2464160.0, 2464189.0, 2464219.0, 2464248.0,
2464278.0, 2464308.0, 2464337.0, 2464367.0, 2464397.0, 2464426.0,
2464456.0, 2464485.0, 2464514.0, 2464544.0, 2464573.0, 2464602.0,
2464632.0, 2464662.0, 2464691.0, 2464721.0, 2464751.0, 2464781.0,
2464810.0, 2464840.0, 2464869.0, 2464898.0, 2464928.0, 2464957.0,
2464986.0, 2465016.0, 2465045.0, 2465075.0, 2465105.0, 2465135.0,
2465165.0, 2465194.0, 2465224.0, 2465253.0, 2465282.0, 2465312.0,
2465341.0, 2465370.0, 2465400.0, 2465429.0, 2465459.0, 2465489.0,
2465519.0, 2465548.0, 2465578.0, 2465607.0, 2465637.0, 2465666.0,
2465696.0, 2465725.0, 2465754.0, 2465784.0, 2465813.0, 2465843.0,
2465873.0, 2465902.0, 2465932.0, 2465962.0, 2465991.0, 2466021.0,
2466050.0, 2466080.0, 2466109.0, 2466139.0, 2466168.0, 2466197.0,
2466227.0, 2466256.0, 2466286.0, 2466316.0, 2466345.0, 2466375.0,
2466404.0, 2466434.0, 2466464.0, 2466493.0, 2466523.0, 2466552.0,
2466581.0, 2466611.0, 2466640.0, 2466670.0, 2466699.0, 2466729.0,
2466759.0, 2466788.0, 2466818.0, 2466848.0, 2466877.0, 2466907.0,
2466936.0, 2466966.0, 2466995.0, 2467024.0, 2467054.0, 2467083.0,
2467113.0, 2467142.0, 2467172.0, 2467202.0, 2467231.0, 2467261.0,
2467291.0, 2467320.0, 2467350.0, 2467379.0, 2467408.0, 2467438.0,
2467467.0, 2467496.0, 2467526.0, 2467556.0, 2467585.0, 2467615.0,
2467645.0, 2467675.0, 2467704.0, 2467734.0, 2467763.0, 2467792.0,
2467822.0, 2467851.0, 2467880.0, 2467910.0, 2467939.0, 2467969.0,
2467999.0, 2468029.0, 2468059.0, 2468088.0, 2468118.0, 2468147.0,
2468176.0, 2468206.0, 2468235.0, 2468264.0, 2468294.0, 2468323.0,
2468353.0, 2468383.0, 2468413.0, 2468442.0, 2468472.0, 2468501.0,
2468531.0, 2468560.0, 2468590.0, 2468619.0, 2468648.0, 2468678.0,
2468707.0, 2468737.0, 2468767.0, 2468796.0, 2468826.0, 2468856.0,
2468885.0, 2468915.0, 2468944.0, 2468974.0, 2469003.0, 2469032.0,
2469062.0, 2469091.0, 2469121.0, 2469150.0, 2469180.0, 2469210.0,
2469239.0, 2469269.0, 2469299.0, 2469328.0, 2469358.0, 2469387.0,
2469416.0, 2469446.0, 2469475.0, 2469505.0, 2469534.0, 2469564.0,
2469593.0, 2469623.0, 2469653.0, 2469682.0, 2469712.0, 2469742.0,
2469771.0, 2469801.0, 2469830.0, 2469859.0, 2469889.0, 2469918.0,
2469948.0, 2469977.0, 2470007.0, 2470036.0, 2470066.0, 2470096.0,
2470125.0, 2470155.0, 2470185.0
};

static double allmonth[] = {
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
8.5, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 3.5, 4.0,
5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0,
5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0,
5.0, 6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
10.5, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 7.0,
8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
7.0, 8.0, 8.5, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0,
5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0,
5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0,
5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0,
4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 2.5,
3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.5, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0,
8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.5,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0,
6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 3.5, 4.0,
5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0,
5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0, 4.0,
5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 11.5, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0, 3.0,
4.0, 5.0, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0, 2.0,
3.0, 4.0, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 1.0,
2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
11.0, 12.0, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
10.0, 11.0, -1.0
};

bool
YearCache::IsInRange(short int year)
{
    return (year >= yearfrom && year <= yearto);
}

void
YearCache::GetYear(short int year, vdouble& vterms, double& lastnew,
    double& lastmon, vdouble& vmoons, vdouble& vmonth, double& nextnew)
{
    int yeardiff = year - yearfrom;
    int termidx = yeardiff * 24;
    vterms.erase(vterms.begin(), vterms.end());
    vterms.insert(vterms.end(), &allterms[termidx], &allterms[termidx + 24]);
    double tstart = julian_date(year, 1, 1, 0.0);
    double tend = julian_date(year + 1, 1, 1, 0.0);
    int moonidx = yeardiff * 12 + yeardiff * 7 / 19 + 1;
    while (allmoons[moonidx - 1] > tstart)
        moonidx--;
    while (allmoons[moonidx] < tstart)
        moonidx++;
    lastnew = allmoons[moonidx - 1];
    lastmon = allmonth[moonidx - 1];
    vmoons.erase(vmoons.begin(), vmoons.end());
    vmonth.erase(vmonth.begin(), vmonth.end());
    while (allmoons[moonidx] < tend)
    {
        vmoons.push_back(allmoons[moonidx]);
        vmonth.push_back(allmonth[moonidx]);
        moonidx++;
    }
    nextnew = allmoons[moonidx];
}

//===========================================================================
//  End of yearcache.cpp
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

  if ( nmonths == 0 )
    return result;

  int year, month, day;
  double rmonth;

  if ( !(gregorianToChinese( date, year, month, day, rmonth )) ) {
    return QDate();
  }

  nmonths = month + nmonths;
  int jday;
  class CYearInfo *ci;

  if ( nmonths <= 0 ) {
    while ( nmonths <= 0 ) {
      year = year - 1;

      ci = getYearInfo( year );

      if ( ci->year == -1 ) {
        return QDate();
      }

      int mtotal = ci->monthinfos.size();

      nmonths = nmonths + mtotal;

      if ( nmonths <= 0 ) {
        continue;
      } else {
        jday = ci->monthinfos[nmonths - 1].beginday + day - 1;
        result = QDate::fromJulianDay( jday );
        return result;
      }
    }
  }

  if ( nmonths > 0 ) {
    while ( nmonths > 0 ) {
      ci = getYearInfo( year );

      if ( ci->year == -1 ) {
        return QDate();
      }

      int mtotal = ci->monthinfos.size();

      if ( nmonths <= mtotal ) {
        jday = ci->monthinfos[nmonths - 1].beginday + day - 1;
        result = QDate::fromJulianDay( jday );
        return result;
      } else {
        year = year + 1;
        nmonths = nmonths - mtotal;
        continue;
      }
    }
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
