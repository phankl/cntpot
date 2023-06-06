#ifndef integrator_header
#define integrator_header

#include <cmath>
#include <vector>

using namespace std;

class Integrator {
private:
  int n, bisectionSteps;
  double bisectionEps;
  vector<double> nodes, weights;

  double legendre(int, double);
  vector<double> initNodes();
  vector<double> initWeights();

public:
  Integrator();
  Integrator(int, int, double);

  template <typename F> double integrate(F function, double a, double b) const {
    double mean = 0.5 * (a + b);
    double halfRange = 0.5 * (b - a);

    double sum = 0;
    for (int i = 0; i < n; i++)
      sum += weights[i] * function(mean + halfRange * nodes[i]);

    return halfRange * sum;
  }
};

#endif
