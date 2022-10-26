# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "log_normal_truncated_ab.hpp"

# include "pb/lib/common.h"
# include "pb/lib/normal.hpp"
# include "pb/lib/log_normal.hpp"

//****************************************************************************80

double log_normal_truncated_ab_cdf ( double x, double mu, double sigma, 
  double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_TRUNCATED_AB_CDF evaluates the Log Normal truncated AB CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 < X.
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Output, double LOG_NORMAL_TRUNCATED_AB_CDF, the value of the CDF.
//
{
  double cdf;
  bool check;
  double lncdf_a;
  double lncdf_b;
  double lncdf_x;

  check = log_normal_truncated_ab_check ( mu, sigma, a, b );

  if ( ! check )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_TRUNCATED_AB_CDF - Fatal error!\n";
    cerr << "  Parameters are not legal.\n";
    exit ( 1 );
  }

  if ( x <= a )
  {
    cdf = 0.0;
  }
  else if ( b <= x )
  {
    cdf = 1.0;
  }
  else
  {
    lncdf_a = log_normal_cdf ( a, mu, sigma );
    lncdf_b = log_normal_cdf ( b, mu, sigma );
    lncdf_x = log_normal_cdf ( x, mu, sigma );

    cdf = ( lncdf_x - lncdf_a ) / ( lncdf_b - lncdf_a );
  }

  return cdf;
}
//****************************************************************************80

double log_normal_truncated_ab_cdf_inv ( double cdf, double mu, double sigma, 
  double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_TRUNCATED_AB_CDF_INV inverts the Log Normal truncated AB CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Input, double LOG_NORMAL_TRUNCATED_AB_CDF_INV, the corresponding argument.
//
{
  bool check;
  double lncdf_a;
  double lncdf_b;
  double lncdf_x;
  double x;

  check = log_normal_truncated_ab_check ( mu, sigma, a, b );

  if ( ! check )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_TRUNCATED_AB_CDF_INV - Fatal error!\n";
    cerr << "  Parameters are not legal.\n";
    exit ( 1 );
  }

  if ( cdf <= 0.0 )
  {
    x = a;
  }
  else if ( 1.0 <= cdf )
  {
    x = b;
  }
  else
  {
    lncdf_a = log_normal_cdf ( a, mu, sigma );
    lncdf_b = log_normal_cdf ( b, mu, sigma );

    lncdf_x = lncdf_a + cdf * ( lncdf_b - lncdf_a );
    x = log_normal_cdf_inv ( lncdf_x, mu, sigma );
  }

  return x;
}
//****************************************************************************80

bool log_normal_truncated_ab_check ( double mu, double sigma, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_TRUNCATED_AB_CHECK checks the Log Normal truncated AB PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Output, bool LOG_NORMAL_TRUNCATED_AB_CHECK, is true if the parameters 
//    are legal.
//
{
  (void)(mu);
  (void)(a);
  (void)(b);
  bool check;

  check = true;

  if ( sigma <= 0.0 )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_TRUNCATED_AB_CHECK - Fatal error!\n";
    cerr << "  SIGMA <= 0.\n";
    check = false;
  }

  return check;
}
//****************************************************************************80

double log_normal_truncated_ab_mean ( double mu, double sigma, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_TRUNCATED_AB_MEAN: mean of the Log Normal truncated AB PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Output, double LOG_NORMAL_TRUNCATED_AB_MEAN, the mean of the PDF.
//
{
  double a0;
  double b0;
  double c1;
  double c2;
  double c3;
  double c4;
  bool check;
  double ln_mean;
  double mean;

  check = log_normal_truncated_ab_check ( mu, sigma, a, b );

  if ( ! check )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_TRUNCATED_AB_MEAN - Fatal error!\n";
    cerr << "  Parameters are not legal.\n";
    exit ( 1 );
  }

  a0 = ( log ( a ) - mu ) / sigma;
  b0 = ( log ( b ) - mu ) / sigma;

  c1 = normal_01_cdf ( sigma - a0 );
  c2 = normal_01_cdf ( sigma - b0 );
  c3 = normal_01_cdf ( + a0 );
  c4 = normal_01_cdf ( + b0 );

  ln_mean = exp ( mu + 0.5 * sigma * sigma );

  mean = ln_mean * ( c1 - c2 ) / ( c4 - c3 );

  return mean;
}
//****************************************************************************80

double log_normal_truncated_ab_pdf ( double x, double mu, double sigma, 
  double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_TRUNCATED_AB_PDF evaluates the Log Normal truncated AB PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 < X
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Output, double LOG_NORMAL_TRUNCATED_AB_PDF, the value of the PDF.
//
{
  bool check;
  double lncdf_a;
  double lncdf_b;
  double lnpdf_x;
  double pdf;

  check = log_normal_truncated_ab_check ( mu, sigma, a, b );

  if ( ! check )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_TRUNCATED_AB_PDF - Fatal error!\n";
    cerr << "  Parameters are not legal.\n";
    exit ( 1 );
  }

  if ( x <= a )
  {
    pdf = 0.0;
  }
  else if ( b <= x )
  {
    pdf = 0.0;
  }
  else
  {
    lncdf_a = log_normal_cdf ( a, mu, sigma );
    lncdf_b = log_normal_cdf ( b, mu, sigma );
    lnpdf_x = log_normal_pdf ( x, mu, sigma );

    pdf = lnpdf_x / ( lncdf_b - lncdf_a );
  }

  return pdf;
}
//****************************************************************************80

double log_normal_truncated_ab_sample ( double mu, double sigma, double a, 
  double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_TRUNCATED_AB_SAMPLE samples the Log Normal truncated AB PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double LOG_NORMAL_TRUNCATED_AB_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  bool check;
  double lncdf_a;
  double lncdf_b;
  double x;

  check = log_normal_truncated_ab_check ( mu, sigma, a, b );

  if ( ! check )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_TRUNCATED_AB_SAMPLE - Fatal error!\n";
    cerr << "  Parameters are not legal.\n";
    exit ( 1 );
  }

  lncdf_a = log_normal_cdf ( a, mu, sigma );
  lncdf_b = log_normal_cdf ( b, mu, sigma );

  cdf = r8_uniform_ab ( lncdf_a, lncdf_b, seed );

  x = log_normal_cdf_inv ( cdf, mu, sigma );

  return x;
}
//****************************************************************************80

double log_normal_truncated_ab_variance ( double mu, double sigma, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_TRUNCATED_AB_VARIANCE: variance of Log Normal truncated AB PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Output, double LOG_NORMAL_TRUNCATED_AB_VARIANCE, the variance of the PDF.
//
{
  double a0;
  double b0;
  double c1;
  double c2;
  double c3;
  double c4;
  bool check;
  double ln_xsquared;
  double lntab_xsquared;
  double mean;
  double variance;

  check = log_normal_truncated_ab_check ( mu, sigma, a, b );

  if ( ! check )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_TRUNCATED_AB_VARIANCE - Fatal error!\n";
    cerr << "  Parameters are not legal.\n";
    exit ( 1 );
  }

  mean = log_normal_truncated_ab_mean ( mu, sigma, a, b );

  a0 = ( log ( a ) - mu ) / sigma;
  b0 = ( log ( b ) - mu ) / sigma;

  c1 = normal_01_cdf ( 2.0 * sigma - a0 );
  c2 = normal_01_cdf ( 2.0 * sigma - b0 );
  c3 = normal_01_cdf ( + a0 );
  c4 = normal_01_cdf ( + b0 );

  ln_xsquared = exp ( 2.0 * mu + 2.0 * sigma * sigma );

  lntab_xsquared = ln_xsquared * ( c1 - c2 ) / ( c4 - c3 );

  variance = lntab_xsquared - mean * mean;

  return variance;
}
//****************************************************************************80
