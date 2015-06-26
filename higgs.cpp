#include <iostream>
#include <cmath>
#include <complex>
#include <ctime> //for seeding the random number with time
#include <cstdlib> 
#include <functional>
#include <vector>
#include "LHAPDF/LHAPDF.h"
#include "gsl/gsl_spline.h"

using namespace std;

complex<double> operator+(complex<double> a, double b){return {a.real()+b, a.imag()};}
complex<double> operator-(complex<double> a, double b){return {a.real()-b, a.imag()};}

class Basic{
private: 
  double Y; 
  complex<double> X;
  double Alpha{0.113671};
  double Gmu = 1.16637E-5;
public:
  Basic(double m0, double mhiggs): Y{m0*m0/mhiggs/mhiggs}, X{(sqrt(1-4*Y)-1.0)/sqrt(1-4*Y)+1.0} {}
  ~Basic()=default;

  complex<double> G0() {return 4.0*Y*pow((1.0+2.0*Y*log(X)),2.0)/2.0;};
  complex<double> G12() {-4.0*Y*(2.0-(1.0-4.0*Y))*pow(X,2.0)/2.0;};

  complex<double> Sigma(){return Gmu*Alpha*Alpha/128.0/sqrt(2.0)/M_PI*pow(abs(G12()),2);};

};

double pdf(double z, const gsl_spline* spline, gsl_interp_accel* acc)
{
  return gsl_spline_eval(spline, z, acc);
}

int main()
{
  double vlist[49],coor[49];

  LHAPDF::PDF *p = LHAPDF::mkPDF("MSTW2008lo68cl",0);

  for(int i=1; i<=10; i++)
    {
      int num1{9+i}; int num2{19+i}; int num3{29+i}; int num5{i-1};
      
      coor[num5] = 1.E-5*(i);
      coor[num1] = 1.E-4*(i+1);
      coor[num2] = 1.E-3*(i+1);
      coor[num3] = 1.E-2*(i+1);

      vlist[num5] = abs(p->xfxQ2(0, coor[num5],125.));
      vlist[num1] = abs(p->xfxQ2(0, coor[num1],125.));
      vlist[num2] = abs(p->xfxQ2(0, coor[num2],125.));
      vlist[num3] = abs(p->xfxQ2(0, coor[num3],125.));
      
    }
  
  for(int i=0; i<9; i++)
    {
      int num1 = i+40;
      coor[num1] = i*0.1+0.2;
      vlist[num1] = abs(p->xfxQ2(0, coor[num1],125.));
    }

   
  gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 49);

  gsl_spline_init (spline,coor,vlist,49);
  double yi;

  for(int i = 0; i < 49; i++)
    {
      yi = gsl_spline_eval (spline, coor[i], acc);
      cout << coor[i] << ","<< abs(yi) << endl;
    }
  double answer = pdf(0.034,spline,acc);
  cout << answer << endl;

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  delete p;
  return 0;
}
//compile with option `lhapdf-config --cflags --ldflags` -lgsl -lgslcblas -lm
