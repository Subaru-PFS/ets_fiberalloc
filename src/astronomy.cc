#include <cmath>
#include <string>
#include <memory>
#include <regex>
#include "lsconstants.h"
#include "math_utils.h"
#include "string_utils.h"
#include "vec3.h"
#include "pointing.h"
#include "rotmatrix.h"
#include "astronomy.h"

using namespace std;

namespace {

double greg2julian (int y, int m, int d)
  {
  if (m<=2) // Jan, Feb
    { --y; m+=12; }
  int a=y/100;
  int b=a/4;
  int c=2-a+b;
  int e=365.25*(y+4716);
  int f=30.6001*(m+1);
  return c+d+e+f-1524.5;
  }

void julian2greg (double jd, int *year, int *month, int *day)
  {
  double q=jd+0.5;
  int z=(int) q;
  int w=(z-1867216.25)/36524.25;
  int x=w/4;
  int a=z+1+w-x;
  int b=a+1524;
  int c=(b-122.1)/365.25;
  int d=365.25*c;
  int e=(b-d)/30.6001;
  int f=30.6001*e;
  *day=b-d-f+(q-z);
  *month=e-1; if (*month>12) month-=12;
  *year=c-4716; if (*month<=2) (*year)++;
  }
#if 0
double jd2gmst (double jd)
  {
  double jd0=(int)(jd+0.5)-0.5;
  double h=(jd-jd0)*24;
  double d=jd-2451545.0;
  double d0=jd0-2451545.0;
  double t=d/36525.;
  double res = 6.697374558 + 0.06570982441908*d0 + 1.00273790935*h + 0.000026*t*t;
  return fmodulo(res,24.);
  }
#else
double jd2gmst (double jd)
  {
  double jd0=jd-2451545.0;
  double t=jd0/36525.;
  double res=280.46061837 + 360.98564736629*jd0 + 3.87933e-4*t*t -t*t*t/3.871e7;
  res/=15;
  return fmodulo(res,24.);
  }
#endif
#if 0
double jd2gast (double jd)
  {
  double gmst=jd2gmst(jd);
  double d=jd-2451545.0;
  double omega=125.04-0.052954*d;
  double l = 280.47 + 0.98565*d;
  double eps = 23.4393 - 0.0000004*d;
  double dpsi=-0.000319*sin(omega*degr2rad) - 0.000024*sin(2*l*degr2rad);
  double res=gmst+dpsi*cos(eps*degr2rad);
  return fmodulo(res,24.);
  }
double jd2gast (double jd)
  {
  double gmst=jd2gmst(jd);
  double t=(jd-2451545.0)/36525.;
  double EPSILONm = 23.439291-0.0130111*t - 1.64E-07*t*t + 5.04E-07*t*t*t;
  double L = 280.4665 + 36000.7698*t;
  double dL = 218.3165 + 481267.8813*t;
  double OMEGA = 125.04452 - 1934.136261*t;

L*=degr2rad; dL*=degr2rad; OMEGA*=degr2rad;
  double dPSI = -17.20*sin(OMEGA) - 1.32*sin(2.*L) - .23*sin(2.*dL)
    + 21.*sin(2.*OMEGA);
  double dEPSILON = 9.20*cos(OMEGA) + .57*cos(2.*L) + .10*cos(2.*dL) -
    .09*cos(2.*OMEGA);

  dPSI = dPSI*(1./3600)*degr2rad;
  dEPSILON = dEPSILON*(1./3600)*degr2rad;
cout <<"BLAH" << dPSI << " " << dEPSILON << endl;
  return fmodulo(gmst + dPSI*cos(EPSILONm*degr2rad+dEPSILON)*rad2degr/15.,24.);
  }
#endif
double jd2gmst_approx (double jd)
  {
  double res = 18.697374558 + 24.06570982441908*(jd-2451545.0);
  return fmodulo(res,24.);
  }

double gmst2ha (double gmst, double lon, double ra) // time in h, angles in rad
  {
  return fmodulo(gmst*15*degr2rad+lon-ra, twopi);
  }
double iso8601toJD (const std::string &datetime)
  {
  regex reg_date(R"foo(^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})Z$)foo");
  std::smatch match;
  planck_assert(regex_search(datetime,match,reg_date),"unknown date format");
  planck_assert(match.size()==7,"unexpected number of matches");
  double jd0=greg2julian(stringToData<int>(match[1]),
                         stringToData<int>(match[2]),
                         stringToData<int>(match[3]));
  jd0 += stringToData<int>(match[4])/24. +
         stringToData<int>(match[5])/(24.*60.) +
         stringToData<double>(match[6])/(24.*60.*60.);
  return jd0;
  }

template<size_t n> double poly (double x, const array<double, n> &c)
  {
  double res=0., v=1.;
  for (size_t i=0; i<c.size(); ++i)
    {
    res+=c[i]*v;
    v*=x;
    }
  return res;
  }

void nutate (double jd, double &d_psi, double &d_eps)
  {
  //  form time in Julian centuries from 1900.0 Hmmm? Looks tather like 2000.0
  double t = (jd - 2451545.)/36525.;

  // Mean elongation of the Moon
  const array<double,4> coeff1 { 297.85036,  445267.111480, -0.0019142, 1./189474 };
  double d=fmodulo(poly(t,coeff1)*degr2rad,twopi);

  // Sun's mean anomaly
  const array<double,4> coeff2 { 357.52772, 35999.050340, -0.0001603, -1./3e5 };
  double m=fmodulo(poly(t,coeff2)*degr2rad,twopi);

  // Moon's mean anomaly
  const array<double,4> coeff3 { 134.96298, 477198.867398, 0.0086972, 1./5.625e4 };
  double mprime = fmodulo(poly(t,coeff3)*degr2rad,twopi);

  // Moon's argument of latitude
  const array<double,4> coeff4 { 93.27191, 483202.017538, -0.0036825, -1./3.27270e5 };
  double f = fmodulo(poly(t,coeff4)*degr2rad,twopi);

  // Longitude of the ascending node of the Moon's mean orbit on the ecliptic,
  // measured from the mean equinox of the date

  const array<double,4> coeff5 { 125.04452, -1934.136261, 0.0020708, 1./4.5e5 };
  double omega = fmodulo(poly(t,coeff5)*degr2rad,twopi);

  const array<double,63> d_lng { 0,-2,0,0,0,0,-2,0,0,-2,-2,-2,0,2,0,2,0,0,-2,0,2,0,0,-2,0,-2,0,0,2,
   -2,0,-2,0,0,2,2,0,-2,0,2,2,-2,-2,2,2,0,-2,-2,0,-2,-2,0,-1,-2,1,0,0,-1,0,0,
     2,0,2};
  const array<double,63> m_lng
 {0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,1,0,-1,0,0,0,1,1,-1,0,
  0,0,0,0,0,-1,-1,0,0,0,1,0,0,1,0,0,0,-1,1,-1,-1,0,-1};
  const array<double,63> mp_lng{0,0,0,0,0,1,0,0,1,0,1,0,-1,0,1,-1,-1,1,2,-2,0,2,2,1,0,0,-1,0,-1,
   0,0,1,0,2,-1,1,0,1,0,0,1,2,1,-2,0,1,0,0,2,2,0,1,1,0,0,1,-2,1,1,1,-1,3,0 };
  const array<double,63> f_lng {0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,2,2,2,2,0,0,2,0,0,
   0,-2,2,2,2,0,2,2,0,2,2,0,0,0,2,0,2,0,2,-2,0,0,0,2,2,0,0,2,2,2,2 };
  const array<double,63> om_lng {1,2,2,2,0,0,2,1,2,2,0,1,2,0,1,2,1,1,0,1,2,2,0,2,0,0,1,0,1,2,1,
   1,1,0,1,2,2,0,2,1,0,2,1,1,1,0,1,1,1,1,1,0,0,0,0,0,2,0,0,2,2,2,2 };
  const array<double,63> sin_lng {-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217,
    -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 26, -22,
     21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, -7,
     6,6,6,-6,-6,5,-5,-5,-5,4,4,4,-4,-4,-4,3,-3,-3,-3,-3,-3,-3,-3};
  const array<double,63> sdelt {-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1,
     0,0,0.1, 0,-0.1,0,0,0,0,0,0,0,0,0,0, -0.1, 0, 0.1, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  const array<double,63> cos_lng { 92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95,0,-70,-53,0,
    -33, 26, 32, 27, 0, -24, 16,13,0,-12,0,0,-10,0,-8,7,9,7,6,0,5,3,-3,0,3,3,
     0,-3,-3,3,3,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  const array<double,63> cdelt {8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  d_psi=d_eps=0.;
  // Sum the periodic terms
  for (size_t n=0; n<d_lng.size(); ++n)
    {
    double arg=d_lng[n]*d + m_lng[n]*m +mp_lng[n]*mprime + f_lng[n]*f +om_lng[n]*omega;
    double sarg = sin(arg);
    double carg = cos(arg);
    d_psi += 0.0001*(sdelt[n]*t + sin_lng[n])*sarg;
    d_eps += 0.0001*(cdelt[n]*t + cos_lng[n])*carg;
    }
  }

void co_nutate (double jd, double &ra, double &dec)
  {
  double d_psi, d_eps;
  nutate (jd,d_psi,d_eps);
  //  form time in Julian centuries from 1900.0 Hmmm? Looks tather like 2000.0
  double t = (jd - 2451545.)/36525.;
  double eps0 = 23.4392911*3600. - 46.8150*t - 0.00059*t*t + 0.001813*t*t*t;
  double eps = (eps0 + d_eps)/3600.*degr2rad; // true obliquity of the ecliptic in radians

  double ce = cos(eps);
  double se = sin(eps);

  // convert ra-dec to equatorial rectangular coordinates
  vec3 p1(pointing(halfpi-dec,ra));
const double d2as=pi/(180.*3600.);
  // apply corrections to each rectangular coordinate
  vec3 p2 (p1.x - (p1.y*ce + p1.z*se)*d_psi * d2as,
           p1.y + (p1.x*ce*d_psi - p1.z*d_eps) * d2as,
           p1.z + (p1.x*se*d_psi + p1.y*d_eps) * d2as);
  pointing pp2(p2);
  dec=halfpi-pp2.theta;
  ra=pp2.phi;
  }
double jd2gast (double jd)
  {
  double d_psi, d_eps;
  nutate (jd,d_psi,d_eps);
  //  form time in Julian centuries from 1900.0 Hmmm? Looks tather like 2000.0
  double t = (jd - 2451545.)/36525.;
  double eps0 = 23.4392911*3600. - 46.8150*t - 0.00059*t*t + 0.001813*t*t*t;
  double eps = (eps0 + d_eps)/3600.*degr2rad; // true obliquity of the ecliptic in radians
  double res=jd2gmst(jd) + d_psi/3600./15. *cos(eps);
  return fmodulo(res,24.);
  }

void precess (double &ra, double &dec, double equinox1, double equinox2)
  {
  const double sec2rad=degr2rad/3600.;
  vec3 x(pointing(halfpi-dec, ra));
  double t = 1e-3*(equinox2-equinox1);
  double st = 1e-3*(equinox1-2000.);
  double A=sec2rad*t*(23062.181 + st*(139.656 +0.0139*st)
    + t*(30.188 - 0.344*st+17.998*t));
  double B=sec2rad*t*t*(79.280 + 0.410*st + 0.205*t) + A;
  double C=sec2rad*t*(20043.109 - st*(85.33 + 0.217*st)
    + t*(-42.665 - 0.217*st -41.833*t));

  double sina = sin(A), sinb = sin(B), sinc = sin(C),
         cosa = cos(A), cosb = cos(B), cosc = cos(C);

  rotmatrix r(
    vec3( cosa*cosb*cosc-sina*sinb,sina*cosb+cosa*sinb*cosc, cosa*sinc),
    vec3(-cosa*sinb-sina*cosb*cosc,cosa*cosb-sina*sinb*cosc,-sina*sinc),
    vec3(-cosb*sinc, -sinb*sinc, cosc));

  vec3 x2 = r.Transform(x); //rotate to get output direction cosines

  pointing ptg(x2);
  ra = ptg.phi;
  ra += (ra<0.)*twopi;
  dec= halfpi-ptg.theta;
  }

double co_refract_forward (double a, double p, double t)
  {
  a*=rad2degr;
  // you have observed the altitude a, and would like to know what the "apparent"
  // altitude is (the one behind the atmosphere).
  bool w = a<15.;
  double r = 0.0166667/tan((a + 7.31/(a+4.4))*degr2rad);

  if (w)
    r = 3.569*(0.1594 + .0196*a + .00002*a*a)/(1.+.505*a+.0845*a*a);
  double tpcor = p/1010. * 283/t;
  r *= tpcor;
  return r*degr2rad;
  }

double co_refract (double alt_in, double alt_observatory)
  {
  double alpha = 0.0065; // temp lapse rate [deg C per meter]

  double temperature=211.5;
  if (alt_observatory<=11000)
    temperature = 283.0 - alpha*alt_observatory;

  // estimate Pressure based on altitude, using U.S. Standard Atmosphere formula.
  double pressure = 1010.*pow(1-6.5/288000*alt_observatory,5.255);
  double epsilon = 0.25; // accuracy of iteration for observed=1 case, in arcseconds

  // calculate initial refraction guess
  double dr = co_refract_forward(alt_in,pressure,temperature);
  double cur = alt_in + dr; //guess of observed location

  double last;
  do {
    last = cur;
    dr = co_refract_forward(cur,pressure,temperature);
    cur= alt_in + dr;
    } while (abs(last-cur)*rad2degr*3600. >= epsilon);
  return cur;
  }

void sunpos (double jd, double &ra, double &dec, double &longmed)
  {
  const double dtor=degr2rad;
  // form time in Julian centuries from 1900.0
  double t = (jd - 2415020.0)/36525.0;

  // form sun's mean longitude
  double l = (279.696678+fmodulo(36000.768925*t,360.0))*3600.0;

  // allow for ellipticity of the orbit (equation of centre)
  // using the Earth's mean anomaly ME

  double me = 358.475844 + (fmodulo(35999.049750*t,360.0));
  double ellcor  = (6910.1 - 17.2*t)*sin(me*degr2rad) + 72.3*sin(2.0*me*degr2rad);
  l+=ellcor;

  // allow for the Venus perturbations using the mean anomaly of Venus MV
  double mv = 212.603219 + fmodulo(58517.803875*t,360.0);
  double vencorr = 4.8 * cos((299.1017 + mv - me)*dtor) +
          5.5 * cos((148.3133 +  2.0 * mv  -  2.0 * me )*dtor) +
          2.5 * cos((315.9433 +  2.0 * mv  -  3.0 * me )*dtor) +
          1.6 * cos((345.2533 +  3.0 * mv  -  4.0 * me )*dtor) +
          1.0 * cos((318.15   +  3.0 * mv  -  5.0 * me )*dtor);
  l+=vencorr;

  // Allow for the Mars perturbations using the mean anomaly of Mars MM
  double mm = 319.529425  +  fmodulo(19139.858500*t,360.0);
  double marscorr = 2.0 * cos((343.8883 -  2.0 * mm  +  2.0 * me)*dtor ) +
            1.8 * cos((200.4017 -  2.0 * mm  + me) * dtor);
  l+=marscorr;

  // Allow for the Jupiter perturbations using the mean anomaly of
  // Jupiter MJ
  double mj = 225.328328  +  fmodulo(3034.6920239*t,360.0);
  double jupcorr = 7.2 * cos(( 179.5317 - mj + me )*dtor) +
          2.6 * cos((263.2167  -  mj ) *dtor) +
          2.7 * cos(( 87.1450  -  2.0 * mj  +  2.0 * me ) *dtor) +
          1.6 * cos((109.4933  -  2.0 * mj  +  me ) *dtor);
  l+=jupcorr;

  // Allow for the Moons perturbations using the mean elongation of
  // the Moon from the Sun D
  double  d = 350.7376814  + fmodulo(445267.11422*t,360.0);
  double mooncorr = 6.5 * sin(d*dtor);
  l+=mooncorr;

  // Allow for long period terms
  double longterm  = + 6.4 * sin(( 231.19  +  20.20 * t )*dtor);
  l+=longterm;
  l  = fmodulo(l+2592000.0,1296000.0);
  longmed = l/3600.0;

  // Allow for Aberration
  l-=20.5;

  // Allow for Nutation using the longitude of the Moons mean node OMEGA
  double omega = 259.183275 - fmodulo(1934.142008*t,360.0);
  l-=17.2 * sin(omega*dtor);

  // Form the True Obliquity
  double oblt  = 23.452294 - 0.0130125*t + (9.2*cos(omega*dtor))/3600.0;

  // Form Right Ascension and Declination
  l/=3600.0;
  ra  = atan2( sin(l*dtor) * cos(oblt*dtor) , cos(l*dtor) );
  if (ra<0) ra+=twopi;
  dec = asin(sin(l*dtor) * sin(oblt*dtor));
  oblt*=dtor;
  longmed*=dtor;
  }

void co_aberration (double &ra, double &dec, double jd)
  {
  const double d2r=degr2rad;
  double T = (jd -2451545.0)/36525.0; // julian centuries from J2000 of jd.
  double d_psi, d_epsilon;
  nutate (jd, d_psi, d_epsilon); // d_psi and d_epsilon in degrees?!
  double eps0 = (23+26/60.+21.448/3600.)*3600. - 46.8150*T - 0.00059*T*T +
               0.001813*T*T*T;
  double eps = (eps0 + d_epsilon)/3600.*degr2rad; // true obliquity of the ecliptic in radians

  double sunra, sundec, sunlon;
  sunpos (jd, sunra, sundec, sunlon);

  // Earth's orbital eccentricity
  double e = 0.016708634 - 0.000042037*T - 0.0000001267*T*T;
  // longitude of perihelion, in degrees
  double pi = 102.93735 + 1.71946*T + 0.00046*T*T;
  double k = 20.49552; //constant of aberration, in arcseconds

  //Useful Trig Functions
  double cd = cos(dec), sd = sin(dec);
  double ce = cos(eps), te = tan(eps);
  double cp = cos(pi*d2r), sp = sin(pi*d2r);
  double cs = cos(sunlon), ss = sin(sunlon);
  double ca = cos(ra), sa = sin(ra);

  double term1 = (ca*cs*ce+sa*ss)/cd,
         term2 = (ca*cp*ce+sa*sp)/cd,
         term3 = (cs*ce*(te*cd-sa*sd)+ca*sd*ss),
         term4 = (cp*ce*(te*cd-sa*sd)+ca*sd*sp);

  double d_ra = -k * term1 + e*k * term2;
  double d_dec = -k * term3 + e*k * term4;
  ra+=d_ra*degr2rad/3600;
  dec+=d_dec*degr2rad/3600;
  }


} // unnamed namespace

void eq2hor_subaru (double ra, double decl, const string &time,
  double &alt, double &az)
  {
  const double j2000= 2451545.0;
  double jd=iso8601toJD(time);
  cout << "jd="<< dataToString(jd) << endl;

  double lat=(19+49/60.+32/3600.)*degr2rad; //Subaru
  double lon=-(155+28/60.+34/3600.)*degr2rad; //Subaru
  double gast=jd2gast(jd);
  cout << dataToString(ra*rad2degr) << " " << dataToString(decl*rad2degr) << endl;
  precess(ra, decl, 2000., 2000. + (jd-j2000) / 365.25);
  cout << dataToString(ra*rad2degr) << " " << dataToString(decl*rad2degr) << endl;
  co_nutate(jd,ra,decl);
  cout << dataToString(ra*rad2degr) << " " << dataToString(decl*rad2degr) << endl;
  co_aberration(ra,decl,jd);
  cout << dataToString(ra*rad2degr) << " " << dataToString(decl*rad2degr) << endl;
  double ha=gmst2ha (gast,lon,ra);
  cout << "HA: " << dataToString(ha*rad2degr)<<endl;
  alt=asin(sin(decl)*sin(lat)+cos(decl)*cos(lat)*cos(ha));
  az=acos((sin(decl)-sin(alt)*sin(lat))/(cos(alt)*cos(lat)));
  if (sin(ha)>0) az=twopi-az;
  alt = co_refract(alt,4139.);
  cout << dataToString(rad2degr*alt) << " " << dataToString(rad2degr*az) << endl;
  }
