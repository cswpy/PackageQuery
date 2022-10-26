#pragma once

double normal_01_cdf ( double x );
double normal_01_cdf_inv ( double cdf );
void normal_01_cdf_values ( int &n_data, double &x, double &fx );
double normal_01_mean ( );
double normal_01_moment ( int order );
double normal_01_pdf ( double x );
double normal_01_sample ( int &seed );
double normal_01_variance ( );

double normal_ms_cdf ( double x, double mu, double sigma );
double normal_ms_cdf_inv ( double cdf, double mu, double sigma );
double normal_ms_mean ( double mu, double sigma );
double normal_ms_moment ( int order, double mu, double sigma );
double normal_ms_moment_central ( int order, double mu, double sigma );
double normal_ms_moment_central_values ( int order, double mu, double sigma );
double normal_ms_moment_values ( int order, double mu, double sigma );
double normal_ms_pdf ( double x, double mu, double sigma );
double normal_ms_sample ( double mu, double sigma, int &seed );
double normal_ms_variance ( double mu, double sigma );