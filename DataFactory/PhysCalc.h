#ifndef _PHYSCALC_H
#define _PHYSCALC_H

// The units here are second, centimeter, sr

double bogdanova(double energy, double theta){
  double es = 18 / (energy * TMath::Cos(theta) + 145) * pow((energy + 2.7 / TMath::Cos(theta)), -2.7) * (energy + 5) / (energy + 5 / TMath::Cos(theta));
  return es;
}


double bogdanova_1d(double *x, double *p){
  return bogdanova(x[0], p[0]);
}

double bogdanova_2d(double* x, double*p) { return bogdanova(x[0], x[1]); }

TF2* func2d_bogdanova(double elow, double eup, double theta_low, double theta_up){
  TF2 *f2 = new TF2("f2", bogdanova_2d, elow, eup, theta_low, theta_up, 0);
  return f2;
}

double reyna_bugaev(double energy, double theta){
  double xp = energy * TMath::Cos(theta);
  double power = 0.2455 + 1.288 * TMath::Log10(xp) - 0.2555 * TMath::Log10(xp) * TMath::Log10(xp) + 0.0209 * TMath::Log10(xp) * TMath::Log10(xp) * TMath::Log10(xp);
  double es = 0.00253 * pow(xp, -power) * pow(TMath::Cos(theta), 3);
  return es;
}

double reyna_bugaev_1d(double *x, double *p){
  return reyna_bugaev(x[0], p[0]);
}

double reyna_bugaev_2d(double* x, double*p) { return reyna_bugaev(x[0], x[1]); }

TF2* func2d_reyna_bugaev(double elow, double eup, double theta_low, double theta_up){
  TF2 *f2 = new TF2("f2", reyna_bugaev_2d, elow, eup, theta_low, theta_up, 0);
  return f2;
}

double reyna_hebbeker(double energy, double theta){
  double y = TMath::Log10(energy * TMath::Cos(theta));
  double power = 0.133 * (y * y * y / 2 - 5 * y * y / 2 + 3 * y) - 2.521 \
          * (-2 * y * y * y / 3 + 3 * y * y - 10 * y / 3 + 1) - 5.78 * (y * y * y / 6 - y * y / 2 + y / 3) - 2.11 \
          * (y * y * y / 3 - 2 * y * y + 11 * y / 3 - 2);
  double es = 0.86 * pow(10, power) * pow(TMath::Cos(theta), 3) / 10000; // /10000 means conversion from m-1 into cm-2
  return es;
}

double reyna_hebbeker_1d(double *x, double *p){
  return reyna_hebbeker(x[0], p[0]);
}

double reyna_hebbeker_2d(double* x, double*p) { return reyna_hebbeker(x[0], x[1]); }

TF2* func2d_reyna_hebbeker(double elow, double eup, double theta_low, double theta_up){
  TF2 *f2 = new TF2("f2", reyna_hebbeker_2d, elow, eup, theta_low, theta_up, 0);
  return f2;
}

double calc_ratio_with_emin(double *x, double *p)
{
  // x[0] = emin, p[0] = theta, p[1] = ecut, p[2] = flux_id
  double ecut = p[1];
  int flux_id = (int) (p[2]+0.1);

  if (x[0] < ecut)
    return 1.;

  TF1 *flux_model;
  if (flux_id == 1)
    flux_model = new TF1("flux_model", bogdanova_1d, ecut, 3000., 1);
  else if (flux_id == 2)
    flux_model = new TF1("flux_model", reyna_bugaev_1d, ecut, 3000., 1);
  else if (flux_id == 3)
    flux_model = new TF1("flux_model", reyna_hebbeker_1d, ecut, 3000., 1);

  flux_model->SetParameters(p);
  double intg = flux_model->Integral(x[0], 3000.0, 1e-12);
  double norm = flux_model->Integral(ecut, 3000.0, 1e-12);
  return intg / norm;
}

double ratio_from_length(double *x, double *p)
{
  // x[0] = length, p[0] = theta, p[1] = ecut, p[2] = flux_id
  // p[3] = _relation_id, p[3+] = parameters of relation

  double emin;
  if (p[3] == 1)
    emin = p[4] + x[0] * p[5];
  else if (p[3] == 2)
    emin = p[4] + x[0] * p[5] + x[0] * x[0] * p[6];
  else
    return -10000;

  x[0] = emin;
  double res = calc_ratio_with_emin(x, p);

  return res;
}


double calc_flux_with_function(double *x, double *p)
{
  // x[0] = emin, p[0] = theta, p[1] = ecut, p[2] = flux_id
  double ecut = p[1];
  int flux_id = (int) (p[2]+0.1);

  TF1 *flux_model;
  if (flux_id == 1)
    flux_model = new TF1("flux_model", bogdanova_1d, ecut, 3000., 1);
  else if (flux_id == 2)
    flux_model = new TF1("flux_model", reyna_bugaev_1d, ecut, 3000., 1);
  else if (flux_id == 3)
    flux_model = new TF1("flux_model", reyna_hebbeker_1d, ecut, 3000., 1);

  flux_model->SetParameters(p);

  if (x[0] < ecut)
    return flux_model->Integral(ecut, 3000.0, 1e-12);

  double intg = flux_model->Integral(x[0], 3000.0, 1e-12);
  return intg;
}


double get_emin_from_ratio(double theta, double ratio, double ecut, int flux_id)
{
  if (ratio > 1.)
    return ecut;

  TF1 *func_emin_ratio = new TF1("func_emin_ratio", calc_ratio_with_emin, ecut - 0.1, 3000., 3);

  double p[3];
  p[0] = theta;
  p[1] = ecut;
  p[2] = (double) flux_id;
  func_emin_ratio->SetParameters(p);

  double emin = -10;
  if (ratio > 1e-10)
    emin = func_emin_ratio->GetX(ratio, ecut - 0.1, 3000., 1e-6, 500);

  if (emin != emin)
  {
    std::cerr << "[get_emin_from_ratio] Cannot find root:  theta = "<< theta <<", ratio = " << ratio << std::endl;
    emin = -10;
   }
  return emin;
}


double leff_to_emin_linear(double *x, double *p){
  // x[0] = leff, p = parameters
  double emin = p[0] + p[1] * x[0];
  return emin;
}

double leff_to_emin_quadratic(double *x, double *p){
  double emin = p[0] + p[1] * x[0] + p[2] * x[0] * x[0];
  return emin;
}

double emin_to_leff_linear(double x, double *p){
  // x[0] = emin, p = parameters
  double emin = x;
  TF1 *func_leff_to_emin_linear = new TF1("func_leff_to_emin_linear", leff_to_emin_linear, 0, 5000, 2);
  func_leff_to_emin_linear->SetParameters(p);

  double leff = func_leff_to_emin_linear->GetX(emin, -20, 5000., 1e-6, 500);

  return leff;
}

double emin_to_leff_quadratic(double x, double *p){
  // x[0] = emin, p = parameters
  double emin = x;
  TF1 *func_leff_to_emin_quadratic = new TF1("func_leff_to_emin_quadratic", leff_to_emin_quadratic, 0, 5000, 3);
  func_leff_to_emin_quadratic->SetParameters(p);

  double leff = func_leff_to_emin_quadratic->GetX(emin, -20, 5000., 1e-6, 500);

  return leff;
}

double get_leff_from_emin(double emin, double ecut, int rel_id, double *p)
{
  if (emin < -9) return -10;
  if (emin < ecut) return 0;

  TF1 *func_leff_emin;
  if (rel_id == 1) func_leff_emin = new TF1("func_leff_emin", leff_to_emin_linear, 0, 3000., 2);
  if (rel_id == 2) func_leff_emin = new TF1("func_leff_emin", leff_to_emin_quadratic, 0, 3000., 3);

  func_leff_emin->SetParameters(p);

  double leff = func_leff_emin->GetX(emin, 0, 5000., 1e-6, 500);
  if (leff != leff){
    std::cerr << "[get_leff_from_emin] Cannot find root:  emin = "<< emin << std::endl;
    }

  return leff;
}

#endif