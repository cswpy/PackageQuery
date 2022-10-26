#include "det_bound.h"

#include "pb/lib/normal.hpp"

void DetBound::computeRotationMatrix(){
  M.resize(m, m);
  double a = 1 / (m + sqrt(m));
  double b = 1 / sqrt(m);
  for (int i = 0; i < m; i ++){
    for (int j = 0; j < m; j ++){
      if (i == m - 1) M(j, i) = b;
      else if (j == m - 1) M(j, i) = -b;
      else if (i == j) M(j, i) = -a + 1;
      else M(j, i) = -a;
    }
  }
}

double DetBound::computeLowerBound(double E, double rho){
  double left = min(0.0, E);
  double right = max(0.0, E);
  double c = left;
  while (fabs(left - right) > eps){
    c = (right + left) / 2;
    double min_prob = 1.0;
    for (int i = 0; i < m; i ++){
      double ou = E * means[i];
      double ov = fabs(E) * vars[i];
      double u = c * means[i];
      double v = fabs(c) * vars[i];
      double sigma = sqrt(v);
      double osigma = sqrt(ov);
      switch (att_senses[i]){
        case LowerBounded:
          min_prob = min(min_prob, 1 - normal_ms_cdf(normal_ms_cdf_inv(1-rho, u, sigma), ou, osigma));
          break;
        case UpperBounded:
          min_prob = min(min_prob, normal_ms_cdf(normal_ms_cdf_inv(rho, u, sigma), ou, osigma));
          break;
        case Bounded:
          min_prob = min(min_prob, normal_ms_cdf(normal_ms_cdf_inv(0.5 + rho/2, u, sigma), ou, osigma) - normal_ms_cdf(normal_ms_cdf_inv(0.5 - rho/2, u, sigma), ou, osigma));
          break;
        default:
          throw invalid_argument("Invalid constraint sense");
      }
    }
    if (min_prob < pow(eps, 1.0/m)) left = c;
    else right = c;
  }
  return c;
}

DetBound::DetBound(){
}

DetBound::DetBound(vector<int> att_senses, vector<double> means, vector<double> vars, int seed, double eps): att_senses(att_senses), means(means), vars(vars), seed(seed), eps(eps){
  m = att_senses.size();
  setSeed(seed);
  computeRotationMatrix();
}

DetBound::DetBound(vector<int> att_senses, VectorXd means, VectorXd vars, int seed, double eps): att_senses(att_senses), eps(eps){
  m = att_senses.size();
  DetBound::means.resize(m);
  VectorXd::Map(&DetBound::means[0], m) = means;
  DetBound::vars.resize(m);
  VectorXd::Map(&DetBound::vars[0], m) = vars;
  setSeed(seed);
  computeRotationMatrix();
}

double DetBound::measureHardness(double E, const VectorXd& bl, const VectorXd& bu){
  double hardness = 0;
  for (int i = 0; i < m; i ++){
    double ou = E * means[i];
    double ov = fabs(E) * vars[i];
    double osigma = sqrt(ov);
    if (bl(i) != -DBL_MAX && bu(i) != DBL_MAX) hardness -= log10(max(normal_ms_cdf(bu(i), ou, osigma) - normal_ms_cdf(bl(i), ou, osigma), eps));
    else if (bl(i) != -DBL_MAX) hardness -= log10(max(1 - normal_ms_cdf(bl(i), ou, osigma), eps));
    else if (bu(i) != DBL_MAX) hardness -= log10(max(normal_ms_cdf(bu(i), ou, osigma), eps));
  }
  return hardness;
}

double DetBound::minHardness(double& minE, const VectorXd& bl, const VectorXd& bu){
  double left = DBL_MAX;
  double right = -DBL_MAX;
  for (int i = 0; i < m; i ++){
    if (bl(i) != -DBL_MAX) left = min(left, bl(i) / means[i]);
    if (bu(i) != DBL_MAX) right = max(right, bu(i) / means[i]);
  }
  while (fabs(left - right) > eps){
    double mid_left = (2*left + right) / 3;
    double mid_right = (left + 2*right) / 3;
    if (measureHardness(mid_left, bl, bu) < measureHardness(mid_right, bl, bu)) right = mid_right;
    else left = mid_left;
  }
  minE = (left + right) / 2;
  return measureHardness(minE, bl, bu);
}

VectorXd DetBound::sampleCenter(double E, double rho, double alpha){
  VectorXd center (m);
  if (m == 1){
    center(0) = E;
  } else{
    double c = computeLowerBound(E, rho);
    double V = alpha * (m - 1) * (E - c) * (E - c);
    normal_distribution<double> normal_dist(0.0, 1.0);
    VectorXd x (m); x.fill(0);
    double norm = 0;
    for (int i = 0; i < m-1; i ++){
      x(i) = normal_dist(gen);
      norm += x(i)*x(i);
    }
    x *= sqrt(m * V / norm);
    center = M * x;
    for (int i = 0; i < m; i ++) center(i) += E;
  }
  return center;
}

void DetBound::assignBound(double rho, const VectorXd& center, VectorXd& bl, VectorXd& bu){
  for (int i = 0; i < m; i ++){
    double u = center(i) * means[i];
    double v = fabs(center(i)) * vars[i];
    double sigma = sqrt(v);
    switch (att_senses[i]){
      case LowerBounded:
        bl(i) = normal_ms_cdf_inv(1-rho, u, sigma);
        bu(i) = DBL_MAX;
        break;
      case UpperBounded:
        bl(i) = -DBL_MAX;
        bu(i) = normal_ms_cdf_inv(rho, u, sigma);
        break;
      case Bounded:
        bl(i) = normal_ms_cdf_inv(0.5 - rho/2, u, sigma);
        bu(i) = normal_ms_cdf_inv(0.5 + rho/2, u, sigma);
        break;
      default:
        throw invalid_argument("Invalid constraint sense");
    }
  }
}

double DetBound::sample(double E, double rho, double alpha, VectorXd& bl, VectorXd& bu){
  assert(isLess(rho, 1));
  VectorXd center = sampleCenter(E, rho, alpha);
  assignBound(rho, center, bl, bu);
  return measureHardness(E, bl, bu);
}

double DetBound::sampleHardness(double E, double alpha, double hardness, VectorXd& bl, VectorXd& bu, double rho){
  VectorXd center = sampleCenter(E, rho, alpha);
  double left = 0.0;
  double right = 1.0;
  while (fabs(left - right) > eps){
    double mid = (left + right) / 2;
    assignBound(mid, center, bl, bu);
    if (measureHardness(E, bl, bu) < hardness) right = mid;
    else left = mid;
  }
  return (left + right) / 2;
}

void DetBound::setSeed(int seed){
  DetBound::seed = seed;
  if (seed < 0) {
    random_device rd;
    DetBound::seed = rd();
  }
  gen.seed(DetBound::seed);
}