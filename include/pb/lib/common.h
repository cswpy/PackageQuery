#pragma once

int i4_uniform_ab ( int a, int b, int &seed );

double r8_abs ( double x );
double r8_choose ( int n, int k );
double r8_factorial2 ( int n );
void r8_factorial2_values ( int &n_data, int &n, double &f );
double r8_huge ( );
double r8_log_2 ( double x );
double r8_mop ( int i );
double r8_uniform_01 ( int &seed );

void r8poly_print ( int n, double a[], string title );
double r8poly_value_horner ( int n, double a[], double x );

double *r8vec_linspace_new ( int n, double a_first, double a_last );
double r8vec_max ( int n, double x[] );
double r8vec_mean ( int n, double x[] );
double r8vec_min ( int n, double x[] );
void r8vec_print ( int n, double a[], string title );
double r8vec_variance ( int n, double x[] );

void timestamp ( );

double r8_hypot ( double x, double y );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8mat_transpose_in_place ( int n, double a[] );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_norm ( int n, double a[] );
void svsort ( int n, double d[], double v[] ) ;

double r8_uniform_ab ( double a, double b, int &seed );