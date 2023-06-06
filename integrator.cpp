#include "integrator.h"

using namespace std;

Integrator::Integrator() {}

Integrator::Integrator(int n, int bisectionSteps, double bisectionEps)
    : n(n), bisectionSteps(bisectionSteps), bisectionEps(bisectionEps),
      nodes(initNodes()), weights(initWeights()) {}

double Integrator::legendre(int m, double x) {
  if (m == 0)
    return 1.0;
  if (m == 1)
    return x;

  std::vector<double> lCache(m + 1);

  lCache[0] = 1.0;
  lCache[1] = x;

  for (int i = 2; i <= m; i++)
    lCache[i] = ((2 * i - 1) * x * lCache[i - 1] - (i - 1) * lCache[i - 2]) / i;

  return lCache[m];
}

vector<double> Integrator::initNodes() {

  double pi = 4 * atan(1);

  vector<double> glNodes(n);

  int kStart, kEnd, kOffset;
  if (n % 2) {
    kStart = 1;
    kEnd = (n - 1) / 2 + 1;
    kOffset = 2;
    glNodes[kEnd - 1] = 0.0;
  } else {
    kStart = 0;
    kEnd = n / 2;
    kOffset = 1;
  }

  int root = 0;
  for (int k = kStart; k < kEnd; k++) {

    double theta = (ceil(0.5 * n) - 0.25 - k) * pi / (n + 0.5);
    double a = cos((ceil(0.5 * n) - k) * pi / (n + 1.0));
    double b = cos(theta);
    double c;

    // perform bisection

    int iter = 0;

    do {
      c = 0.5 * (a + b);
      if (legendre(n, c) == 0.0)
        break;
      if (legendre(n, a) * legendre(n, c) < 0)
        b = c;
      else
        a = c;
      iter++;
    } while (fabs(a - b) >= bisectionEps && iter <= bisectionSteps);

    glNodes[kEnd + root] = c;
    glNodes[kEnd - root - kOffset] = -c;
    root++;
  }

  return glNodes;
}

vector<double> Integrator::initWeights() {
  vector<double> glWeights(n);

  for (int i = 0; i < n; i++) {
    double x = nodes[i];
    double dlegendre =
        n * (x * legendre(n, x) - legendre(n - 1, x)) / (x * x - 1.0);

    glWeights[i] = 2.0 / ((1.0 - x * x) * dlegendre * dlegendre);
  }

  return glWeights;
}
