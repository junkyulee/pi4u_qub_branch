int i4_uniform_ab ( int a, int b, int *seed );

float normal_01_cdf ( float x );
float normal_01_cdf_inv ( float cdf );
void normal_01_cdf_values ( int *n_data, float *x, float *fx );
float normal_01_mean ( );
float normal_01_moment ( int order );
float normal_01_pdf ( float x );
float normal_01_sample ( int *seed );
float normal_01_variance ( );

float normal_ms_cdf ( float x, float mu, float sigma );
float normal_ms_cdf_inv ( float cdf, float mu, float sigma );
float normal_ms_mean ( float mu, float sigma );
float normal_ms_moment ( int order, float mu, float sigma );
float normal_ms_moment_central ( int order, float mu, float sigma );
float normal_ms_moment_central_values ( int order, float mu, float sigma );
float normal_ms_moment_values ( int order, float mu, float sigma );
float normal_ms_pdf ( float x, float mu, float sigma );
float normal_ms_sample ( float mu, float sigma, int *seed );
float normal_ms_variance ( float mu, float sigma );

float r8_abs ( float x );
float r8_choose ( int n, int k );
float r8_factorial2 ( int n );
void r8_factorial2_values ( int *n_data, int *n, float *f );
float r8_huge ( );
float r8_log_2 ( float x );
float r8_mop ( int i );
float r8_uniform_01 ( int *seed );

void r8poly_print ( int n, float a[], char *title );
float r8poly_value_horner ( int n, float a[], float x );

float *r8vec_linspace_new ( int n, float a, float b );
float r8vec_max ( int n, float x[] );
float r8vec_mean ( int n, float x[] );
float r8vec_min ( int n, float x[] );
void r8vec_print ( int n, float a[], char *title );
float r8vec_variance ( int n, float x[] );

void timestamp ( );

float truncated_normal_ab_cdf ( float x, float mu, float sigma, float a, 
  float b );
void truncated_normal_ab_cdf_values ( int *n_data, float *mu, float *sigma, 
  float *a, float *b, float *x, float *fx );
float truncated_normal_ab_cdf_inv ( float cdf, float mu, float sigma, float a, 
  float b );
float truncated_normal_ab_mean ( float mu, float sigma, float a, float b );
float truncated_normal_ab_moment ( int order, float mu, float sigma, float a, float b );
float truncated_normal_ab_pdf ( float x, float mu, float sigma, float a, 
  float b );
void truncated_normal_ab_pdf_values ( int *n_data, float *mu, float *sigma, 
  float *a, float *b, float *x, float *fx );
float truncated_normal_ab_sample ( float mu, float sigma, float a, float b, 
  int *seed );
float truncated_normal_ab_variance ( float mu, float sigma, float a, float b );

float truncated_normal_a_cdf ( float x, float mu, float sigma, float a );
void truncated_normal_a_cdf_values ( int *n_data, float *mu, float *sigma, 
  float *a, float *x, float *fx );
float truncated_normal_a_cdf_inv ( float cdf, float mu, float sigma, float a );
float truncated_normal_a_mean ( float mu, float sigma, float a );
float truncated_normal_a_moment ( int order, float mu, float sigma, float a );
float truncated_normal_a_pdf ( float x, float mu, float sigma, float a );
void truncated_normal_a_pdf_values ( int *n_data, float *mu, float *sigma, 
  float *a, float *x, float *fx );
float truncated_normal_a_sample ( float mu, float sigma, float a, int *seed );
float truncated_normal_a_variance ( float mu, float sigma, float a );

float truncated_normal_b_cdf ( float x, float mu, float sigma, float b );
void truncated_normal_b_cdf_values ( int *n_data, float *mu, float *sigma, 
  float *b, float *x, float *fx );
float truncated_normal_b_cdf_inv ( float cdf, float mu, float sigma, float b );
float truncated_normal_b_mean ( float mu, float sigma, float b );
float truncated_normal_b_moment ( int order, float mu, float sigma, float b );
float truncated_normal_b_pdf ( float x, float mu, float sigma, float b );
void truncated_normal_b_pdf_values ( int *n_data, float *mu, float *sigma, 
  float *b, float *x, float *fx );
float truncated_normal_b_sample ( float mu, float sigma, float b, int *seed );
float truncated_normal_b_variance ( float mu, float sigma, float b );


float log_normal_truncated_ab_pdf ( float x, float mu, float sigma, float a, float b );
