#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"

double SpinDensity(double *x_val, double *par)
{
  double x = x_val[0];
  double rho00 = par[0];
  double Norm = par[1]; // Norm = 0.75 <= sampling | Norm from fitting <= extract rho00

  double dNdCosThetaStar = Norm*((1.0-rho00)+(3.0*rho00-1)*x*x);

  return dNdCosThetaStar;
}

double PolyRes(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*(x-1.0/3.0) + 1.0/3.0;

  return y;
}

double flowSample(double *x_val, double *par) 
{
  double x, y, v2;
  x = x_val[0];
  v2 = par[0];
  y = 1.0 + 2.0*v2*cos(2.0*x);

  return y;
}

double flow(double *x_val, double *par) 
{
  double x, y, v2, Norm;
  x = x_val[0];
  v2 = par[0];
  Norm = par[1];
  y = Norm*(1.0 + 2.0*v2*cos(2.0*x));

  return y;
}

double EventPlaneGaus(double *var, double *par)
{
  double Psi2 = var[0];
  double res = par[0]; // res = <cos(2*(Psi2-Psi_RP))>

  double sigma = acos(res)/2.0;
  double sigmaSquare = sigma*sigma;
  double norm = 1.0/sqrt(2.0*sigmaSquare*TMath::Pi());
  double power = -1.0*Psi2*Psi2/(2.0*sigmaSquare);

  double y = norm*exp(power);

  return y;
}

double EventPlaneResolution(double *var, double *par)
{
  double chi = var[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  double y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

double EventPlaneDist(double *var, double *par)
{
  double DeltaPsi = var[0];
  double chi = par[0];
  double arg = chi/TMath::Sqrt(2.0);
  double arg2 = -0.5*chi*chi;
  double pi = TMath::Pi();
  double norm = par[1];

  double cos = TMath::Cos(2.0*DeltaPsi);
  double sin2 = TMath::Sin(2.0*DeltaPsi)*TMath::Sin(2.0*DeltaPsi);
  double y = norm*(TMath::Exp(arg2)+TMath::Sqrt(pi)*arg*cos*TMath::Exp(arg2*sin2)*(1.0+TMath::Erf(arg*cos)));

  return y;
}

double v2_pT_FitFunc(double *x_val, double *par)
{
  // Fit function for v2 vs. pT
  // From arXiv:nucl-th/0403030v5: Resonance decay effects on anisotrotpy parameters
  double v2, pT, a, b, c, d, n;
  pT = x_val[0];
  n  = par[0]; // number-of-constituent quarks
  a  = par[1];
  b  = par[2];
  c  = par[3];
  d  = par[4];

  if(c != 0.0)
  {
    v2 = a*n/(1.0 + exp(-(pT/n - b)/c)) - d*n;
  }
  else v2 = 0.0;

  return v2;
}

double Levy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = numer*power/denom;

  return y;
}

double pTLevy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = pT*numer*power/denom;

  return y;
}

double meanLevy(double *var, double *par)
{
  double const m0 = 1.01940; // phi-meson mass
  double pT   = var[0];
  double mT   = sqrt(pT*pT+m0*m0);
  double dNdy = par[0];
  double n    = par[1];
  double T    = par[2];

  double numer = dNdy*(n-1)*(n-2);
  double denom = n*T*(n*T+m0*(n-2));
  double power = pow(1+(mT-m0)/(n*T),-1.0*n);

  double y = pT*pT*numer*power/denom;

  return y;
}

//---------------Back Ground Subtraction-----------------------
double BreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  return BW;
}

double PolyBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double Poly = par[3] + par[4]*x;

  double y = BW + Poly;

  return y;
}

double Poly(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*x;

  return y;
}

double LegeBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double p1 = x;
  double p2 = 0.5*(3.0*x*x-1.0);
  double Para = par[3] + par[4]*p1 + par[5]*p2;
  // double Para = par[3]*(1.0 + par[4]*p1 + par[5]*p2);

  double y = BW + Para;

  return y;
}

double Lege(double *x_val, double *par)
{
  double x = x_val[0];
  double p1 = x;
  double p2 = 0.5*(3.0*x*x-1.0);
  double y = par[0] + par[1]*p1 + par[2]*p2;
  // double y = par[0]*(1.0 + par[1]*p1 + par[2]*p2);

  return y;
}

double ChebBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double t1 = x;
  double t2 = 2.0*x*x-1.0;
  double Cheb = par[3] + par[4]*t1 + par[5]*t2;
  // double Cheb = par[3]*(1.0 + par[4]*t1 + par[5]*t2);

  double y = BW + Cheb;

  return y;
}

double Cheb(double *x_val, double *par)
{
  double x = x_val[0];
  double t1 = x;
  double t2 = 2.0*x*x-1.0;
  double y = par[0] + par[1]*t1 + par[2]*t2;
  // double y = par[0]*(1.0 + par[1]*t1 + par[2]*t2);

  return y;
}

double LogaBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double log = par[3] + par[4]*TMath::Log(par[5]*x+par[6]);

  double y = BW + log;

  return y;
}

double Loga(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*TMath::Log(par[2]*x+par[3]);

  return y;
}

double SqRtBreitWigner(double *x_val, double *par)
{
  double x = x_val[0];
  double m0 = par[0];
  double Gamma = par[1];
  double Norm = par[2];

  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double BW = Norm*Gamma/denom;

  double sqrt = par[3] + par[4]*TMath::Sqrt(x+par[5]);

  double y = BW + sqrt;

  return y;
}

double SqRt(double *x_val, double *par)
{
  double x = x_val[0];
  double y = par[0] + par[1]*TMath::Sqrt(x+par[2]);

  return y;
}

double DELPHIBreitWigner(double *var, double *par)
{
  double x = var[0];
  double m0 = par[0];
  double Gamma = par[1];
  double nsig = par[2];
  double denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  double sig = nsig*Gamma/denom;

  double mMassKaon = 0.49368;
  double c0 = par[3];
  double c1 = par[4];
  double c2 = par[5];
  double c3 = par[6];
  double c4 = par[7];
  double nbkg = par[8];
  double x_var = x-mMassKaon*2.0;
  double bkg = nbkg*TMath::Power(x_var,c0)*TMath::Exp(c1*x+c2*x*x+c3*x*x*x+c4*x*x*x*x);

  return sig+bkg;
}

double DELPHI(double *var, double *par)
{
  double mMassKaon = 0.49368;
  double x = var[0];
  double c0 = par[0];
  double c1 = par[1];
  double c2 = par[2];
  double c3 = par[3];
  double c4 = par[4];
  double nbkg = par[5];
  double x_var = x-mMassKaon*2.0;
  double bkg = nbkg*TMath::Power(x_var,c0)*TMath::Exp(c1*x+c2*x*x+c3*x*x*x+c4*x*x*x*x);

  return bkg;
}

//---------------Back Ground Subtraction-----------------------

float ErrorAdd(float x, float y)
{
  return sqrt(x*x+y*y);
}

float ErrTimes(float x, float y, float dx, float dy)
{
  return x*y*ErrorAdd(dx/x,dy/y);
}

float ErrDiv(float x, float y, float dx, float dy)
{
  return x/y*ErrorAdd(dx/x,dy/y);
}


// full event plane resolution
double Resolution_Full(double *x_val, double *par)
{
  double y;
  double chi = x_val[0];
  double arg = chi*chi/4.0;
  double norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}
