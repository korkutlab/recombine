double soft(double x, double d);

double sum(double *x, int n);

double max(double *x, int n);

void pmax(double *x, int n, double cutoff);

double max_fabs(double *x, int n);

double l2n(double *x, int n);

double l1n(double *x, int n);

double diff_l1n(double *x1, double *x2, int n);

void normalize_by_l2n(double *x, int n);

double l1n_after_normalized_by_l2n(double *x, int n);

double BinarySearch(double *argu, int n, double sumabs);

double rand_gen();

double normalRandom();

void matrix_multiply(double *A, double *B, int m, int n, int p, double *C);

double cross_prod(double *x1, double *x2, int n);

void v2w(double *v, int p, int G,
         int *group_ps, int *group_offsets, int *group_idxs,
         double *w);

void update_v(double *dtu, int G,
              int *group_ps, int *group_offsets, int *group_idxs,
              double lambda2,
              double *v);

double expectation_approx(double *beta, double a, double b, int p, int l);

double pstar(double x, double theta, double lambda1, double lambda0);

double lambdastar(double x, double theta, double lambda1, double lambda0);

double g(double x, double theta, double lambda1, double lambda0);

double threshold(double theta, double lambda1, double lambda0);

double crossprod_Dj_u(double *D, double *u, int n2, int j);

double SSL(double dtu, double w, double lambda0, double lambda1, double theta, double delta);

int checkConvergence(double *beta, double *beta_old, double eps, int p);

