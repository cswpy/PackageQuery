# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "log_normal.hpp"

# include "pb/lib/common.h"
# include "pb/lib/normal.hpp"

//****************************************************************************80

double log_normal_cdf ( double x, double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF evaluates the Log Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Output, double LOG_NORMAL_CDF, the value of the CDF.
//
{
  double cdf;
  double logx;

  if ( x <= 0.0 )
  {
    cdf = 0.0;
  }
  else
  {
    logx = log ( x );

    cdf = normal_cdf ( logx, mu, sigma );
  }

  return cdf;
}
//****************************************************************************80

double log_normal_cdf_inv ( double cdf, double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF_INV inverts the Log Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Input, double LOG_NORMAL_CDF_INV, the corresponding argument.
//
{
  double logx;
  double x;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cerr << "\n";
    cerr << "LOG_NORMAL_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  logx = normal_cdf_inv ( cdf, mu, sigma );

  x = exp ( logx );

  return x;
}
//****************************************************************************80

void log_normal_cdf_values ( int &n_data, double &mu, double &sigma,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = Log NormalDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &MU, the mean of the distribution.
//
//    Output, double &SIGMA, the shape parameter of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
     0.2275013194817921E-01,
     0.2697049307349095E+00,
     0.5781741008028732E+00,
     0.7801170895122241E+00,
     0.4390310097476894E+00,
     0.4592655190218048E+00,
     0.4694258497695908E+00,
     0.4755320473858733E+00,
     0.3261051056816658E+00,
     0.1708799040927608E+00,
     0.7343256357952060E-01,
     0.2554673736161761E-01 };

  static double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double sigma_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

bool log_normal_check ( double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CHECK checks the parameters of the Log Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Output, bool LOG_NORMAL_CHECK, is true if the parameters are legal.
//
{
  (void)(mu);
  bool check;

  check = true;

  if ( sigma <= 0.0 )
  {
    cout << "\n";
    cout << "LOG_NORMAL_CHECK - Fatal error!\n";
    cout << "  SIGMA <= 0.\n";
    check = false;
  }

  return check;
}
//****************************************************************************80

double log_normal_mean ( double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_MEAN returns the mean of the Log Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Output, double LOG_NORMAL_MEAN, the mean of the PDF.
//
{
  double mean;

  mean = exp ( mu + 0.5 * sigma * sigma );

  return mean;
}
//****************************************************************************80

double log_normal_pdf ( double x, double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_PDF evaluates the Log Normal PDF.
//
//  Discussion:
//
//    PDF(A,B;X)
//      = exp ( - 0.5 * ( ( log ( X ) - MU ) / SIGMA )^2 )
//        / ( SIGMA * X * sqrt ( 2 * PI ) )
//
//    The Log Normal PDF is also known as the Cobb-Douglas PDF,
//    and as the Antilog_normal PDF.
//
//    The Log Normal PDF describes a variable X whose logarithm
//    is normally distributed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Output, double LOG_NORMAL_PDF, the value of the PDF.
//
{
  double pdf;
  const double r8_pi = 3.14159265358979323;
  double y;

  if ( x <= 0.0 )
  {
    pdf = 0.0;
  }
  else
  {
    y = ( log ( x ) - mu ) / sigma;
    pdf = exp ( - 0.5 * y * y ) / ( sigma * x * sqrt ( 2.0 * r8_pi ) );
  }

  return pdf;
}
//****************************************************************************80

double log_normal_sample ( double mu, double sigma, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_SAMPLE samples the Log Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double LOG_NORMAL_SAMPLE, a sample of the PDF.
//
{
  double cdf;
  double x;

  cdf = r8_uniform_01 ( seed );

  x = log_normal_cdf_inv ( cdf, mu, sigma );

  return x;
}
//****************************************************************************80

double log_normal_variance ( double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_VARIANCE returns the variance of the Log Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Output, double LOG_NORMAL_VARIANCE, the variance of the PDF.
//
{
  double variance;

  variance = exp ( 2.0 * mu + sigma * sigma ) * ( exp ( sigma * sigma ) - 1.0 );

  return variance;
}
//****************************************************************************80

double normal_cdf ( double x, double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF evaluates the Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Output, double CDF, the value of the CDF.
//
{
  double cdf;
  double y;

  y = ( x - mu ) / sigma;

  cdf = normal_01_cdf ( y );

  return cdf;
}
//****************************************************************************80

double normal_cdf_inv ( double cdf, double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF_INV inverts the Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Algorithm AS 111,
//    Applied Statistics,
//    Volume 26, pages 118-121, 1977.
//
//  Parameters:
//
//    Input, double CDF, the value of the CDF.
//    0.0 <= CDF <= 1.0.
//
//    Input, double MU, SIGMA, the parameters of the PDF.
//    0.0 < SIGMA.
//
//    Output, double NORMAL_CDF_INV, the corresponding argument.
//
{
  double x;
  double x2;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cerr << "\n";
    cerr << "NORMAL_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  x2 = normal_01_cdf_inv ( cdf );

  x = mu + sigma * x2;

  return x;
}
//****************************************************************************80

bool normal_check ( double mu, double sigma )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CHECK checks the parameters of the Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2004
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
//    Output, bool NORMAL_CHECK, is true if the parameters are legal.
//
{
  (void)(mu);
  if ( sigma <= 0.0 )
  {
    cout << "\n";
    cout << "NORMAL_CHECK - Fatal error!\n";
    cout << "  SIGMA <= 0.\n";
    return false;
  }

  return true;
}
//****************************************************************************80