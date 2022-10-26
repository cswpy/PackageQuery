# include <cmath>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "truncated_normal.hpp"

# include "pb/lib/common.h"
# include "pb/lib/normal.hpp"

double truncated_normal_ab_cdf ( double x, double mu, double sigma, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_CDF evaluates the truncated Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, B, the lower and upper truncation limits.
//
//    Output, double TRUNCATED_NORMAL_AB_CDF, the value of the CDF.
//
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double cdf;
  double xi;
  double xi_cdf;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else if ( x <= b )
  {
    alpha = ( a - mu ) / sigma;
    beta = ( b - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    alpha_cdf = normal_01_cdf ( alpha );
    beta_cdf = normal_01_cdf ( beta );
    xi_cdf = normal_01_cdf ( xi );

    cdf = ( xi_cdf - alpha_cdf ) / ( beta_cdf - alpha_cdf );
  }
  else
  {
    cdf = 1.0;
  }
  
  return cdf;
}
//****************************************************************************80

void truncated_normal_ab_cdf_values ( int &n_data, double &mu, double &sigma, 
  double &a, double &b, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_CDF_VALUES: values of the Truncated Normal CDF.
//
//  Discussion:
//
//    The Normal distribution, with mean Mu and standard deviation Sigma,
//    is truncated to the interval [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Output, double &SIGMA, the standard deviation of the distribution.
//
//    Output, double &A, &B, the lower and upper truncation limits.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0 };

  static double b_vec[N_MAX] = {
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0 };

  static double fx_vec[N_MAX] = {
    0.3371694242213513, 
    0.3685009225506048, 
    0.4006444233448185, 
    0.4334107066903040, 
    0.4665988676496338, 
    0.5000000000000000, 
    0.5334011323503662, 
    0.5665892933096960, 
    0.5993555766551815, 
    0.6314990774493952, 
    0.6628305757786487 };

  static double mu_vec[N_MAX] = {
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0 };

  static double sigma_vec[N_MAX] = {
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0  };

  static double x_vec[N_MAX] = {
     90.0, 
     92.0, 
     94.0, 
     96.0, 
     98.0, 
    100.0, 
    102.0, 
    104.0, 
    106.0, 
    108.0, 
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_ab_cdf_inv ( double cdf, double mu, double sigma, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_CDF_INV inverts the truncated Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2013
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
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, B, the lower and upper truncation limits.
//
//    Output, double TRUNCATED_NORMAL_AB_CDF_INV, the corresponding argument.
//
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double x;
  double xi;
  double xi_cdf;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_AB_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );

  xi_cdf = ( beta_cdf - alpha_cdf ) * cdf + alpha_cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_ab_mean ( double mu, double sigma, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_MEAN returns the mean of the truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, B, the lower and upper truncation limits.
//
//    Output, double TRUNCATED_NORMAL_AB_MEAN, the mean of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double beta;
  double beta_cdf;
  double beta_pdf;
  double mean;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta  );

  alpha_pdf = normal_01_pdf ( alpha );
  beta_pdf = normal_01_pdf ( beta );

  mean = mu + sigma * ( alpha_pdf - beta_pdf ) / ( beta_cdf - alpha_cdf );

  return mean;
}
//****************************************************************************80

double truncated_normal_ab_moment ( int order, double mu, double sigma, 
  double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_MOMENT: moments of the truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Phoebus Dhrymes,
//    Moments of Truncated Normal Distributions,
//    May 2005.
//
//  Parameters:
//
//    Input, int ORDER, the order of the moment.
//    0 <= ORDER.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//    0.0 < S.
//
//    Input, double A, B, the lower and upper truncation limits.
//    A < B.
//
//    Output, double TRUNCATED_NORMAL_AB_MOMENT, the moment of the PDF.
//
{
  double a_cdf;
  double a_h;
  double a_pdf;
  double b_cdf;
  double b_h;
  double b_pdf;
  double ir;
  double irm1;
  double irm2;
  double moment;
  int r;

  if ( order < 0 )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n";
    cerr << "  ORDER < 0.\n";
    exit ( 1 );
  }

  if ( sigma <= 0.0 )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n";
    cerr << "  SIGMA <= 0.0.\n";
    exit ( 1 );
  }

  if ( b <= a )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n";
    cerr << "  B <= A.\n";
    exit ( 1 );
  }

  a_h = ( a - mu ) / sigma;
  a_pdf = normal_01_pdf ( a_h );
  a_cdf = normal_01_cdf ( a_h );

  if ( a_cdf == 0.0 )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n";
    cerr << "  PDF/CDF ratio fails, because A_CDF too small.\n";
    cerr << "  A_PDF = " << a_pdf << "\n";
    cerr << "  A_CDF = " << a_cdf << "\n";
    exit ( 1 );
  }

  b_h = ( b - mu ) / sigma;
  b_pdf = normal_01_pdf ( b_h );
  b_cdf = normal_01_cdf ( b_h );

  if ( b_cdf == 0.0 )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_AB_MOMENT - Fatal error!\n";
    cerr << "  PDF/CDF ratio fails, because B_CDF too small.\n";
    cerr << "  B_PDF = " << b_pdf << "\n";
    cerr << "  B_CDF = " << b_cdf << "\n";
    exit ( 1 );
  }

  moment = 0.0;
  irm2 = 0.0;
  irm1 = 0.0;

  for ( r = 0; r <= order; r++ )
  {
    if ( r == 0 )
    {
      ir = 1.0;
    }
    else if ( r == 1 )
    {
      ir = - ( b_pdf - a_pdf ) / ( b_cdf - a_cdf );
    }
    else
    {
      ir = ( double ) ( r - 1 ) * irm2 
        - ( pow ( b_h, r - 1 ) * b_pdf - pow ( a_h, r - 1 ) * a_pdf )
        / ( b_cdf - a_cdf );
    }

    moment = moment + r8_choose ( order, r ) * pow ( mu, order - r ) 
      * pow ( sigma, r ) * ir;

    irm2 = irm1;
    irm1 = ir;
  }

  return moment;
}
//****************************************************************************80

double truncated_normal_ab_pdf ( double x, double mu, double sigma, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_PDF evaluates the truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, B, the lower and upper truncation limits.
//
//    Output, double TRUNCATED_NORMAL_AB_PDF, the value of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double pdf;
  double xi;
  double xi_pdf;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else if ( x <= b )
  {
    alpha = ( a - mu ) / sigma;
    beta = ( b - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    alpha_cdf = normal_01_cdf ( alpha );
    beta_cdf = normal_01_cdf ( beta );
    xi_pdf = normal_01_pdf ( xi );

    pdf = xi_pdf / ( beta_cdf - alpha_cdf ) / sigma;
  }
  else
  {
    pdf = 0.0;
  }
  
  return pdf;
}
//****************************************************************************80

void truncated_normal_ab_pdf_values ( int &n_data, double &mu, double &sigma, 
  double &a, double &b, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_PDF_VALUES: values of the Truncated Normal PDF.
//
//  Discussion:
//
//    The Normal distribution, with mean Mu and standard deviation Sigma,
//    is truncated to the interval [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Output, double &SIGMA, the standard deviation of the distribution.
//
//    Output, double &A, &B, the lower and upper truncation limits.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0 };

  static double b_vec[N_MAX] = {
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0 };

  static double fx_vec[N_MAX] = {
    0.01543301171801836,
    0.01588394472270638,
    0.01624375997031919,
    0.01650575046469259,
    0.01666496869385951,
    0.01671838200940538,
    0.01666496869385951,
    0.01650575046469259,
    0.01624375997031919,
    0.01588394472270638,
    0.01543301171801836 };

  static double mu_vec[N_MAX] = {
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0 };

  static double sigma_vec[N_MAX] = {
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0  };

  static double x_vec[N_MAX] = {
     90.0, 
     92.0, 
     94.0, 
     96.0, 
     98.0, 
    100.0, 
    102.0, 
    104.0, 
    106.0, 
    108.0, 
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_ab_sample ( double mu, double sigma, double a, double b, 
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, B, the lower and upper truncation limits.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double TRUNCATED_NORMAL_AB_SAMPLE, a sample of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double beta;
  double beta_cdf;
  double u;
  double x;
  double xi;
  double xi_cdf;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );

  u = r8_uniform_01 ( seed );
  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf );
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_ab_variance ( double mu, double sigma, double a, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_VARIANCE returns the variance of the truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, B, the lower and upper truncation limits.
//
//    Output, double TRUNCATED_NORMAL_AB_VARIANCE, the variance of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double beta;
  double beta_cdf;
  double beta_pdf;
  double variance;

  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;

  alpha_pdf = normal_01_pdf ( alpha );
  beta_pdf = normal_01_pdf ( beta );

  alpha_cdf = normal_01_cdf ( alpha );
  beta_cdf = normal_01_cdf ( beta );

  variance = sigma * sigma * ( 1.0 
    + ( alpha * alpha_pdf - beta * beta_pdf ) / ( beta_cdf - alpha_cdf ) 
    - pow ( ( alpha_pdf - beta_pdf ) / ( beta_cdf - alpha_cdf ), 2 ) );

  return variance;
}
//****************************************************************************80

double truncated_normal_a_cdf ( double x, double mu, double sigma, double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_CDF evaluates the lower truncated Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_CDF, the value of the CDF.
//
{
  double alpha;
  double alpha_cdf;
  double cdf;
  double xi;
  double xi_cdf;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else
  {
    alpha = ( a - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    alpha_cdf = normal_01_cdf ( alpha );
    xi_cdf = normal_01_cdf ( xi );

    cdf = ( xi_cdf - alpha_cdf ) / ( 1.0 - alpha_cdf );
  }
  
  return cdf;
}
//****************************************************************************80

void truncated_normal_a_cdf_values ( int &n_data, double &mu, double &sigma, 
  double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_CDF_VALUES: values of the Lower Truncated Normal CDF.
//
//  Discussion:
//
//    The Normal distribution, with mean Mu and standard deviation Sigma,
//    is truncated to the interval [A,+oo).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Output, double &SIGMA, the standard deviation of the distribution.
//
//    Output, double &A, the lower truncation limit.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0 };

  static double fx_vec[N_MAX] = {
    0.3293202045481688, 
    0.3599223134505957, 
    0.3913175216041539, 
    0.4233210140873113, 
    0.4557365629792204, 
    0.4883601253415709, 
    0.5209836877039214, 
    0.5533992365958304, 
    0.5854027290789878, 
    0.6167979372325460, 
    0.6474000461349729 };

  static double mu_vec[N_MAX] = {
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0 };

  static double sigma_vec[N_MAX] = {
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0  };

  static double x_vec[N_MAX] = {
     90.0, 
     92.0, 
     94.0, 
     96.0, 
     98.0, 
    100.0, 
    102.0, 
    104.0, 
    106.0, 
    108.0, 
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_a_cdf_inv ( double cdf, double mu, double sigma, 
  double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_CDF_INV inverts the lower truncated Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2013
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
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_CDF_INV, the corresponding argument.
//
{
  double alpha;
  double alpha_cdf;
  double x;
  double xi;
  double xi_cdf;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_A_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  alpha = ( a - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );

  xi_cdf = ( 1.0 - alpha_cdf ) * cdf + alpha_cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_a_mean ( double mu, double sigma, double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_MEAN returns the mean of the lower truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_MEAN, the mean of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double mean;

  alpha = ( a - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );

  alpha_pdf = normal_01_pdf ( alpha );

  mean = mu + sigma * alpha_pdf / ( 1.0 - alpha_cdf );

  return mean;
}
//****************************************************************************80

double truncated_normal_a_moment ( int order, double mu, double sigma, 
  double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_MOMENT: moments of the lower truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Phoebus Dhrymes,
//    Moments of Truncated Normal Distributions,
//    May 2005.
//
//  Parameters:
//
//    Input, int ORDER, the order of the moment.
//    0 <= ORDER.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_MOMENT, the moment of the PDF.
//
{
  double moment;

  moment = r8_mop ( order )
    * truncated_normal_b_moment ( order, - mu, sigma, - a );

  return moment;
}
//****************************************************************************80

double truncated_normal_a_pdf ( double x, double mu, double sigma, double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_PDF evaluates the lower truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_PDF, the value of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double pdf;
  double xi;
  double xi_pdf;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else
  {
    alpha = ( a - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    alpha_cdf = normal_01_cdf ( alpha );
    xi_pdf = normal_01_pdf ( xi );

    pdf = xi_pdf / ( 1.0 - alpha_cdf ) / sigma;
  }
  
  return pdf;
}
//****************************************************************************80

void truncated_normal_a_pdf_values ( int &n_data, double &mu, double &sigma, 
  double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_PDF_VALUES: values of the Lower Truncated Normal PDF.
//
//  Discussion:
//
//    The Normal distribution, with mean Mu and standard deviation Sigma,
//    is truncated to the interval [A,+oo).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Output, double &SIGMA, the standard deviation of the distribution.
//
//    Output, double &A, the lower truncation limit.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0, 
     50.0 };

  static double fx_vec[N_MAX] = {
     0.01507373507401876,
     0.01551417047139894,
     0.01586560931024694,
     0.01612150073158793,
     0.01627701240029317,
     0.01632918226724295,
     0.01627701240029317,
     0.01612150073158793,
     0.01586560931024694,
     0.01551417047139894,
     0.01507373507401876 };

  static double mu_vec[N_MAX] = {
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0 };

  static double sigma_vec[N_MAX] = {
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0  };

  static double x_vec[N_MAX] = {
     90.0, 
     92.0, 
     94.0, 
     96.0, 
     98.0, 
    100.0, 
    102.0, 
    104.0, 
    106.0, 
    108.0, 
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_a_sample ( double mu, double sigma, double a, 
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_SAMPLE samples the lower truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double TRUNCATED_NORMAL_A_SAMPLE, a sample of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double u;
  double x;
  double xi;
  double xi_cdf;

  alpha = ( a - mu ) / sigma;

  alpha_cdf = normal_01_cdf ( alpha );

  u = r8_uniform_01 ( seed );
  xi_cdf = alpha_cdf + u * ( 1.0 - alpha_cdf );
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_a_variance ( double mu, double sigma, double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_VARIANCE: variance of the lower truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_VARIANCE, the variance of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double alpha_pdf;
  double variance;

  alpha = ( a - mu ) / sigma;

  alpha_pdf = normal_01_pdf ( alpha );

  alpha_cdf = normal_01_cdf ( alpha );

  variance = sigma * sigma * ( 1.0 
    + alpha * alpha_pdf / ( 1.0 - alpha_cdf ) 
    - pow ( alpha_pdf / ( 1.0 - alpha_cdf ), 2 ) );

  return variance;
}
//****************************************************************************80

double truncated_normal_b_cdf ( double x, double mu, double sigma, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_CDF evaluates the upper truncated Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double B, the upper truncation limit.
//
//    Output, double TRUNCATED_NORMAL_B_CDF, the value of the CDF.
//
{
  double beta;
  double beta_cdf;
  double cdf;
  double xi;
  double xi_cdf;

  if ( x <= b )
  {
    beta = ( b - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    beta_cdf = normal_01_cdf ( beta );
    xi_cdf = normal_01_cdf ( xi );

    cdf = xi_cdf / beta_cdf;
  }
  else
  {
    cdf = 1.0;
  }
  
  return cdf;
}
//****************************************************************************80

void truncated_normal_b_cdf_values ( int &n_data, double &mu, double &sigma, 
  double &b, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_CDF_VALUES: values of the upper Truncated Normal CDF.
//
//  Discussion:
//
//    The Normal distribution, with mean Mu and standard deviation Sigma,
//    is truncated to the interval (-oo,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Output, double &SIGMA, the standard deviation of the distribution.
//
//    Output, double &B, the upper truncation limit.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double b_vec[N_MAX] = {
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0 };

  static double fx_vec[N_MAX] = {
    0.3525999538650271, 
    0.3832020627674540, 
    0.4145972709210122, 
    0.4466007634041696, 
    0.4790163122960786, 
    0.5116398746584291, 
    0.5442634370207796, 
    0.5766789859126887, 
    0.6086824783958461, 
    0.6400776865494043, 
    0.6706797954518312 };

  static double mu_vec[N_MAX] = {
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0 };

  static double sigma_vec[N_MAX] = {
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0  };

  static double x_vec[N_MAX] = {
     90.0, 
     92.0, 
     94.0, 
     96.0, 
     98.0, 
    100.0, 
    102.0, 
    104.0, 
    106.0, 
    108.0, 
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_b_cdf_inv ( double cdf, double mu, double sigma, 
  double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_CDF_INV inverts the upper truncated Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2013
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
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double B, the upper truncation limit.
//
//    Output, double TRUNCATED_NORMAL_B_CDF_INV, the corresponding argument.
//
{
  double beta;
  double beta_cdf;
  double x;
  double xi;
  double xi_cdf;

  if ( cdf < 0.0 || 1.0 < cdf )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_B_CDF_INV - Fatal error!\n";
    cerr << "  CDF < 0 or 1 < CDF.\n";
    exit ( 1 );
  }

  beta = ( b - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );

  xi_cdf = beta_cdf * cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_b_mean ( double mu, double sigma, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_MEAN returns the mean of the upper truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double B, the upper truncation limit.
//
//    Output, double TRUNCATED_NORMAL_B_MEAN, the mean of the PDF.
//
{
  double beta;
  double beta_cdf;
  double beta_pdf;
  double mean;

  beta = ( b - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );

  beta_pdf = normal_01_pdf ( beta );

  mean = mu - sigma * beta_pdf / beta_cdf;

  return mean;
}
//****************************************************************************80

double truncated_normal_b_moment ( int order, double mu, double sigma, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_MOMENT: moments of the upper truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Phoebus Dhrymes,
//    Moments of Truncated Normal Distributions,
//    May 2005.
//
//  Parameters:
//
//    Input, int ORDER, the order of the moment.
//    0 <= ORDER.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double B, the upper truncation limit.
//
//    Output, double TRUNCATED_NORMAL_B_MOMENT, the moment of the PDF.
//
{
  double f;
  double h;
  double h_cdf;
  double h_pdf;
  double ir;
  double irm1;
  double irm2;
  double moment;
  int r;

  if ( order < 0 )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_B_MOMENT - Fatal error!\n";
    cerr << "  ORDER < 0.\n";
    exit ( 1 );
  }

  h = ( b - mu ) / sigma;
  h_pdf = normal_01_pdf ( h );
  h_cdf = normal_01_cdf ( h );

  if ( h_cdf == 0.0 )
  {
    cerr << "\n";
    cerr << "TRUNCATED_NORMAL_B_MOMENT - Fatal error!\n";
    cerr << "  CDF((B-MU)/SIGMA) = 0.\n";
    exit ( 1 );
  }

  f = h_pdf / h_cdf;

  moment = 0.0;
  irm2 = 0.0;
  irm1 = 0.0;

  for ( r = 0; r <= order; r++ )
  {
    if ( r == 0 )
    {
      ir = 1.0;
    }
    else if ( r == 1 )
    {
      ir = - f;
    }
    else
    {
      ir = - pow ( h, r - 1 ) * f + ( double ) ( r - 1 ) * irm2;
    }

    moment = moment + r8_choose ( order, r ) * pow ( mu, order - r ) 
      * pow ( sigma, r ) * ir;

    irm2 = irm1;
    irm1 = ir;
  }

  return moment;
}
//****************************************************************************80

double truncated_normal_b_pdf ( double x, double mu, double sigma, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_PDF evaluates the upper truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double B, the upper truncation limit.
//
//    Output, double TRUNCATED_NORMAL_B_PDF, the value of the PDF.
//
{
  double beta;
  double beta_cdf;
  double pdf;
  double xi;
  double xi_pdf;

  if ( x <= b )
  {
    beta = ( b - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    beta_cdf = normal_01_cdf ( beta );
    xi_pdf = normal_01_pdf ( xi );

    pdf = xi_pdf / beta_cdf / sigma;
  }
  else
  {
    pdf = 0.0;
  }
  
  return pdf;
}
//****************************************************************************80

void truncated_normal_b_pdf_values ( int &n_data, double &mu, double &sigma, 
  double &b, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_PDF_VALUES: values of the Upper Truncated Normal PDF.
//
//  Discussion:
//
//    The Normal distribution, with mean Mu and standard deviation Sigma,
//    is truncated to the interval (-oo,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Output, double &SIGMA, the standard deviation of the distribution.
//
//    Output, double &B, the upper truncation limit.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double b_vec[N_MAX] = {
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0, 
     150.0 };

  static double fx_vec[N_MAX] = {
    0.01507373507401876, 
    0.01551417047139894, 
    0.01586560931024694, 
    0.01612150073158793, 
    0.01627701240029317, 
    0.01632918226724295, 
    0.01627701240029317, 
    0.01612150073158793, 
    0.01586560931024694, 
    0.01551417047139894, 
    0.01507373507401876 };

  static double mu_vec[N_MAX] = {
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0, 
     100.0 };

  static double sigma_vec[N_MAX] = {
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0, 
    25.0  };

  static double x_vec[N_MAX] = {
     90.0, 
     92.0, 
     94.0, 
     96.0, 
     98.0, 
    100.0, 
    102.0, 
    104.0, 
    106.0, 
    108.0, 
    110.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    b = 0.0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    b = b_vec[n_data-1];
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double truncated_normal_b_sample ( double mu, double sigma, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_SAMPLE samples the upper truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double B, the upper truncation limit.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double TRUNCATED_NORMAL_B_SAMPLE, a sample of the PDF.
//
{
  double beta;
  double beta_cdf;
  double u;
  double x;
  double xi;
  double xi_cdf;

  beta = ( b - mu ) / sigma;

  beta_cdf = normal_01_cdf ( beta );

  u = r8_uniform_01 ( seed );
  xi_cdf = u * beta_cdf;
  xi = normal_01_cdf_inv ( xi_cdf );

  x = mu + sigma * xi;

  return x;
}
//****************************************************************************80

double truncated_normal_b_variance ( double mu, double sigma, double b )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_VARIANCE: variance of the upper truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double B, the upper truncation limit.
//
//    Output, double TRUNCATED_NORMAL_B_VARIANCE, the variance of the PDF.
//
{
  double beta;
  double beta_cdf;
  double beta_pdf;
  double variance;

  beta = ( b - mu ) / sigma;

  beta_pdf = normal_01_pdf ( beta );

  beta_cdf = normal_01_cdf ( beta );

  variance = sigma * sigma * ( 1.0 
    - beta * beta_pdf / beta_cdf 
    - pow ( beta_pdf / beta_cdf, 2 ) );

  return variance;
}
