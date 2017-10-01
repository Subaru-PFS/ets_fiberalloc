#include <cmath>
#include <string>
#include <memory>
#include <regex>
#include <array>
#include "lsconstants.h"
#include "math_utils.h"
#include "string_utils.h"
#include "vec3.h"
#include "pointing.h"
#include "rotmatrix.h"

using namespace std;

namespace {

class eq2hor
  {
  private:
    const double j2000= 2451545.0;
    const double sec2rad=degr2rad/3600.;

    rotmatrix npmat;
    double lat, lon, gast, eps, jd, sunlon, altitude;

  public:
    eq2hor (double lat_obs, double lon_obs, double alt_obs, const std::string &time);
    pointing radec2altaz (const pointing &radec) const;
  };


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

double jd2gmst (double jd)
  {
  double jd0=jd-2451545.0;
  double t=jd0/36525.;
  double res=280.46061837 + 360.98564736629*jd0 + 3.87933e-4*t*t -t*t*t/3.871e7;
  res/=15;
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

void nutate (double jd, double &d_psi, double &d_eps) // result in arcsec
  {
  double t = (jd - 2451545.)/36525.;

  // Mean elongation of the Moon
  const array<double,4> coeff1 {{297.85036,445267.111480,-0.0019142,1./189474}};
  double d=fmodulo(poly(t,coeff1)*degr2rad,twopi);

  // Sun's mean anomaly
  const array<double,4> coeff2 {{357.52772, 35999.050340, -0.0001603, -1./3e5}};
  double m=fmodulo(poly(t,coeff2)*degr2rad,twopi);

  // Moon's mean anomaly
  const array<double,4> coeff3 {{134.96298,477198.867398,0.0086972,1./5.625e4}};
  double mprime = fmodulo(poly(t,coeff3)*degr2rad,twopi);

  // Moon's argument of latitude
  const array<double,4> coeff4
    {{93.27191,483202.017538,-0.0036825,-1./3.2727e5}};
  double f = fmodulo(poly(t,coeff4)*degr2rad,twopi);

  // Longitude of the ascending node of the Moon's mean orbit on the ecliptic,
  // measured from the mean equinox of the date
  const array<double,4> coeff5 {{125.04452, -1934.136261, 0.0020708, 1./4.5e5}};
  double omega = fmodulo(poly(t,coeff5)*degr2rad,twopi);

  const array<double,63> d_lng {{ 0,-2,0,0,0,0,-2,0,0,-2,-2,-2,0,2,0,2,0,0,-2,0,
    2,0,0,-2,0,-2,0,0,2,-2,0,-2,0,0,2,2,0,-2,0,2,2,-2,-2,2,2,0,-2,-2,0,-2,-2,0,
    -1,-2,1,0,0,-1,0,0, 2,0,2}};
  const array<double,63> m_lng {{ 0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,2,0,2,1,0,-1,0,0,0,1,1,-1,0,0,0,0,0,0,-1,-1,0,0,0,1,0,0,1,0,0,0,
    -1,1,-1,-1,0,-1}};
  const array<double,63> mp_lng {{0,0,0,0,0,1,0,0,1,0,1,0,-1,0,1,-1,-1,1,2,-2,0,
    2,2,1,0,0,-1,0,-1,0,0,1,0,2,-1,1,0,1,0,0,1,2,1,-2,0,1,0,0,2,2,0,1,1,0,0,1,
    -2,1,1,1,-1,3,0 }};
  const array<double,63> f_lng {{0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,
    2,2,2,2,0,0,2,0,0,0,-2,2,2,2,0,2,2,0,2,2,0,0,0,2,0,2,0,2,-2,0,0,0,2,2,0,0,2,
    2,2,2 }};
  const array<double,63> om_lng {{1,2,2,2,0,0,2,1,2,2,0,1,2,0,1,2,1,1,0,1,2,2,0,
    2,0,0,1,0,1,2,1,1,1,0,1,2,2,0,2,1,0,2,1,1,1,0,1,1,1,1,1,0,0,0,0,0,2,0,0,2,2,
    2,2 }};
  const array<double,63> sin_lng {{-171996, -13187, -2274, 2062, 1426, 712,
    -517, -386, -301, 217, -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38,
    -31, 29, 29, 26, -22, 21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7,
    -7, -7, 6,6,6,-6,-6,5,-5,-5,-5,4,4,4,-4,-4,-4,3,-3,-3,-3,-3,-3,-3,-3}};
  const array<double,63> sdelt {{-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4,
    0, -0.5, 0, 0.1, 0,0,0.1, 0,-0.1,0,0,0,0,0,0,0,0,0,0, -0.1, 0, 0.1,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }};
  const array<double,63> cos_lng {{ 92025, 5736, 977, -895, 54, -7, 224, 200,
    129, -95,0,-70,-53,0, -33, 26, 32, 27, 0, -24, 16,13,0,-12,0,0,-10,0,-8,7,9,
    7,6,0,5,3,-3,0,3,3,0,-3,-3,3,3,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }};
  const array<double,63> cdelt {{8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0,
    -0.1, 0.3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};

  d_psi=d_eps=0.;
  // Sum the periodic terms
  for (size_t n=0; n<d_lng.size(); ++n)
    {
    double arg=d_lng[n]*d + m_lng[n]*m + mp_lng[n]*mprime
               + f_lng[n]*f +om_lng[n]*omega;
    double sarg = sin(arg);
    double carg = cos(arg);
    d_psi += 0.0001*(sdelt[n]*t + sin_lng[n])*sarg;
    d_eps += 0.0001*(cdelt[n]*t + cos_lng[n])*carg;
    }
  }

double co_refract_forward (double a, double p, double t)
  {
  a*=rad2degr;
  // you have observed the altitude a, and would like to know what the
  // "apparent" altitude is (the one behind the atmosphere).
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
  double temperature = max(211.5,283.0 - 0.0065*alt_observatory);

  // estimate Pressure based on altitude, using U.S. Standard Atmosphere formula
  double pressure = 1010.*pow(1-6.5/288000*alt_observatory,5.255);
  double epsilon = 0.25; // accuracy of iteration for observed=1 case, in arcsec

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

double sunpos2 (double jd)
  {
  const double dtor=degr2rad;
  // form time in Julian centuries from 1900.0
  double t = (jd - 2415020.0)/36525.0;

  // form sun's mean longitude
  double l = (279.696678+fmodulo(36000.768925*t,360.0))*3600.0;

  // allow for ellipticity of the orbit (equation of centre)
  // using the Earth's mean anomaly ME
  double me = 358.475844 + (fmodulo(35999.049750*t,360.0));
  double ellcor =(6910.1 - 17.2*t)*sin(me*degr2rad) + 72.3*sin(2.0*me*degr2rad);
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
  return l/3600.0*dtor;
  }
void co_aberration2 (double &ra, double &dec, double jd, double sunlon, double eps)
  {
  const double d2r=degr2rad;
  double T = (jd -2451545.0)/36525.0; // julian centuries from J2000 of jd.

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


    static rotmatrix axis_rotation_matrix(const vec3 &axis, double angle)
      {
      rotmatrix res;
      res.Make_Axis_Rotation_Transform (axis, angle);
      return res;
      }

    static rotmatrix get_precession_matrix (double equinox1, double equinox2)
      {
      const double sec2rad=degr2rad/3600.;
      double t = 1e-3*(equinox2-equinox1);
      double st = 1e-3*(equinox1-2000.);
      double A=sec2rad*t*(23062.181 + st*(139.656 +0.0139*st)
        + t*(30.188 - 0.344*st+17.998*t));
      double B=sec2rad*t*t*(79.280 + 0.410*st + 0.205*t) + A;
      double C=sec2rad*t*(20043.109 - st*(85.33 + 0.217*st)
        + t*(-42.665 - 0.217*st -41.833*t));

      double sina = sin(A), sinb = sin(B), sinc = sin(C),
             cosa = cos(A), cosb = cos(B), cosc = cos(C);

      return rotmatrix(
        vec3( cosa*cosb*cosc-sina*sinb,sina*cosb+cosa*sinb*cosc, cosa*sinc),
        vec3(-cosa*sinb-sina*cosb*cosc,cosa*cosb-sina*sinb*cosc,-sina*sinc),
        vec3(-cosb*sinc, -sinb*sinc, cosc));
      }

    eq2hor::eq2hor (double lat_obs, double lon_obs, double alt_obs, const string &time)
      {
      lat=lat_obs; lon=lon_obs; altitude=alt_obs;
      jd=iso8601toJD(time);
      double d_psi, d_eps;
      nutate(jd, d_psi, d_eps);
      double t = (jd - 2451545.)/36525.;
      double eps0 = 23.4392911*3600. - 46.8150*t - 0.00059*t*t + 0.001813*t*t*t;
      // true obliquity of the ecliptic in radians
      eps = (eps0 + d_eps)/3600.*degr2rad;
      gast=fmodulo(jd2gmst(jd) + d_psi/3600./15. *cos(eps),24.);
      rotmatrix prec_mat=get_precession_matrix(2000., 2000. + (jd-j2000)/365.25);
      rotmatrix nut_mat =   axis_rotation_matrix(vec3(1,0,0),eps)
                           *(axis_rotation_matrix(vec3(0,0,1),d_psi*sec2rad)
                           *axis_rotation_matrix(vec3(1,0,0),-eps0*sec2rad));
      npmat=nut_mat*prec_mat;
      sunlon=sunpos2(jd);
      }
    pointing eq2hor::radec2altaz (const pointing &radec) const
      {
      pointing ptg(radec);
      ptg=pointing(npmat.Transform(vec3(ptg)));
      double ra=ptg.phi;
      ra+=(ra<0.)*twopi;
      double dec=halfpi-ptg.theta;
      co_aberration2(ra,dec,jd,sunlon,eps);
      double ha=gmst2ha (gast,lon,ra);
      double alt=asin(sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(ha));
      double az=acos((sin(dec)-sin(alt)*sin(lat))/(cos(alt)*cos(lat)));
      if (sin(ha)>0) az=twopi-az;
      alt=co_refract(alt,altitude);
      return pointing(halfpi-alt,az);
      }

const double obs_lat=(19+49/60.+32/3600.)*degr2rad,
                 obs_lon=-(155+28/60.+34/3600.)*degr2rad;
// Height of the observatory (meters above sea level)
const double obs_height=4139.;
vector<complex<double>> targetToPFI(const vector<pointing> &altaz, const pointing &los, double psi)
  {
  // altitude and azimuth of North celestial pole:
  pointing altaz_ncp (halfpi-obs_lat,0.);
  vec3 z{los}, skypole(altaz_ncp);
  vec3 x=(skypole-z*dotprod(z,skypole)).Norm();
  vec3 y=crossprod(z,x);
  xcomplex<double> cpsi(cos(psi),sin(psi));
  const double a0=0., a1=-3.2e2, a2=-1.37e1, a3=-7.45e0;
  vector<xcomplex<double>> res;
  for (auto&& t:altaz)
    {
    vec3 pos(t);
    vec3 xp=pos-y*dotprod(pos,y);
    vec3 yp=pos-x*dotprod(pos,x);
    xcomplex<double> pnew (atan2(dotprod(xp,x),dotprod(xp,z))*rad2degr,
               atan2(dotprod(yp,y),dotprod(yp,z))*rad2degr);
    pnew*=cpsi;
    double rsq=norm(pnew);
    res.emplace_back((a3*rsq*rsq+a2*rsq+a1)*pnew.real()+a0,
                 (-a3*rsq*rsq-a2*rsq-a1)*pnew.imag()+a0);
    }
  return res;
  }
/*! Converts RA/DEC in degrees to colatitude/longitude in radians. */
inline pointing radec2ptg (double ra, double dec)
  { return pointing((90-dec)*degr2rad,ra*degr2rad); }

} // unnamed namespace

vector<complex<double>> cconv (const vector<double> &ra, const vector<double> &dec,
  double tel_ra, double tel_dec, double psi, const string &time)
  {
  psi+=90.; // correction to match declination with PFI y-axis
  vector<pointing> t0;
  for (size_t i=0; i<ra.size(); ++i)
    t0.emplace_back(radec2ptg(ra[i],dec[i]));
  eq2hor eqtest(obs_lat, obs_lon, obs_height, time);
  for (auto &t:t0)
    t=eqtest.radec2altaz(t);
  vector<complex<double>> res;

  return targetToPFI(t0,eqtest.radec2altaz(radec2ptg(tel_ra,tel_dec)),psi*degr2rad);
  }
