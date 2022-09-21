#include "det_bound.h"

const double kRho = 0.5;

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
      switch (consSense[i]){
        case LowerBounded:
          min_prob = min(min_prob, 1 - normalCdf(ou, ov, normalQuantile(u, v, 1-rho)));
          break;
        case UpperBounded:
          min_prob = min(min_prob, normalCdf(ou, ov, normalQuantile(u, v, rho)));
          break;
        case Bounded:
          min_prob = min(min_prob, normalCdf(ou, ov, normalQuantile(u, v, 0.5 + rho/2)) - normalCdf(ou, ov, normalQuantile(u, v, 0.5 - rho/2)));
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

DetBound::DetBound(vector<int> consSense, vector<double> means, vector<double> vars, int seed, double eps): consSense(consSense), means(means), vars(vars), seed(seed), eps(eps){
  m = consSense.size();
  if (seed < 0) {
    random_device rd;
    DetBound::seed = rd();
  }
  gen.seed(DetBound::seed);
  computeRotationMatrix();
}

double DetBound::measureHardness(double E, const VectorXd& bl, const VectorXd& bu){
  double hardness = 0;
  for (int i = 0; i < m; i ++){
    double ou = E * means[i];
    double ov = fabs(E) * vars[i];
    if (bl(i) != -DBL_MAX && bu(i) != DBL_MAX) hardness -= log10(max(normalCdf(ou, ov, bu(i)) - normalCdf(ou, ov, bl(i)), eps));
    else if (bl(i) != -DBL_MAX) hardness -= log10(max(1 - normalCdf(ou, ov, bl(i)), eps));
    else if (bu(i) != DBL_MAX) hardness -= log10(max(normalCdf(ou, ov, bu(i)), eps));
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
    switch (consSense[i]){
      case LowerBounded:
        bl(i) = normalQuantile(u, v, 1-rho);
        bu(i) = DBL_MAX;
        break;
      case UpperBounded:
        bl(i) = -DBL_MAX;
        bu(i) = normalQuantile(u, v, rho);
        break;
      case Bounded:
        bl(i) = normalQuantile(u, v, 0.5 - rho/2);
        bu(i) = normalQuantile(u, v, 0.5 + rho/2);
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

double DetBound::sampleHardness(double E, double alpha, double hardness, VectorXd& bl, VectorXd& bu){
  VectorXd center = sampleCenter(E, kRho, alpha);
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