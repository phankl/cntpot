#include "potential.h"

CubicSpline uInfGeneration() {
  double hStart = 2 * R_CNT + DELTA0;
  double hEnd = 2 * R_CNT + R_C;
  double hSpacing = (hEnd - hStart) / (UINF_POINTS - 1);

  std::vector<double> uInfData(UINF_POINTS, 0);
  std::vector<std::vector<double>> fileData(UINF_POINTS,
                                            std::vector<double>(2, 0));

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < UINF_POINTS; i++) {
    double h = hStart + i * hSpacing;
    double uDensity = uExactInfDensity(h, 0, 0);
    uInfData[i] = uDensity;
    fileData[i][0] = h;
    fileData[i][1] = uDensity;
  }

  io::printData(fileData);

  CubicSpline spline(uInfData, hStart, hSpacing);
  return spline;
}

CubicSpline gammaOrthGeneration(const CubicSpline &uInfParaSpline) {
  double hStart = 0.0;
  double hEnd = 2 * R_CNT + R_C - DELTA0;
  double hSpacing = (hEnd - hStart) / (GAMMA_POINTS - 1);

  std::vector<double> gammaOrthData(GAMMA_POINTS, 0);
  std::vector<std::vector<double>> fileData(GAMMA_POINTS,
                                            std::vector<double>(3, 0));

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < GAMMA_POINTS; i++) {
    double h = hStart + i * hSpacing;

    // Define lambda functions for minimisation
    auto uInfLambda = [h](double xi) -> double {
      return uExactInfDensity(h, 0.5 * PI, xi);
    };
    auto uInfBarLambda = [h, uInfParaSpline](double xi) -> double {
      double hBar = sqrt(h * h + xi * xi);
      if (hBar >= 2 * R_CNT + R_C)
        return 0;
      else
        return uInfParaSpline(hBar);
    };

    int precision_bits = std::numeric_limits<double>::digits;
    double xiUpper = sqrt(pow(2 * R_CNT + R_C, 2) - h * h);

    std::pair<double, double> uInfMinimum =
        boost::math::tools::brent_find_minima(uInfLambda, -DELTA0, xiUpper,
                                              precision_bits);
    std::pair<double, double> uInfBarMinimum =
        boost::math::tools::brent_find_minima(uInfBarLambda, -DELTA0, xiUpper,
                                              precision_bits);

    double gamma = uInfMinimum.second / uInfBarMinimum.second;
    double omega;
    if (std::fabs(uInfMinimum.first) < 1.0e-3)
      omega = 1.0;
    else
      omega = std::fabs(uInfBarMinimum.first / uInfMinimum.first);
    if (omega < 1.0)
      omega = 1.0;

    gammaOrthData[i] = gamma;
    fileData[i][0] = h;
    fileData[i][1] = gamma;
    // fileData[i][2] = omega;
  }

  io::printData(fileData);

  CubicSpline spline(gammaOrthData, hStart, hSpacing);
  return spline;
}

BicubicSpline phiGeneration(const CubicSpline &uInfParaSpline) {
  double hStart = 0.0;
  double hEnd = 2 * R_CNT + R_C;
  double hSpacing = (hEnd - hStart) / (PHI_POINTS - 1);
  double psiSpacing = 1.0 / (PHI_POINTS - 1);

  std::vector<std::vector<double>> fileData(PHI_POINTS * PHI_POINTS,
                                            std::vector<double>(3, 0));
  std::vector<CubicSpline> hSplines(PHI_POINTS, CubicSpline());

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < PHI_POINTS; i++) {
    double h = hStart + i * hSpacing;

    auto zetaIntegrand = [h, uInfParaSpline](double zeta) -> double {
      double hBar = sqrt(h * h + zeta * zeta);
      if (hBar >= 2 * R_CNT + R_C)
        return 0;
      else
        return uInfParaSpline(hBar);
    };

    double hMax = 2 * R_CNT + R_C;
    double zetaLower = zetaMin(h);
    double zetaUpper = sqrt(hMax * hMax - h * h);
    double zetaSpacing = (zetaUpper - zetaLower) / (PHI_POINTS - 1);

    std::vector<double> constantHData(PHI_POINTS, 0);

    for (int j = 0; j < PHI_POINTS; j++) {
      double zeta = zetaLower + j * zetaSpacing;
      double psi = j * psiSpacing;
      double phi;
      if (zeta == zetaLower)
        phi = 0;
      else
        phi = LENGTH_INTEGRATOR.integrate(zetaIntegrand, zetaLower, zeta);

      constantHData[j] = phi;
      fileData[i * PHI_POINTS + j][0] = h;
      fileData[i * PHI_POINTS + j][1] = psi;
      fileData[i * PHI_POINTS + j][2] = phi;
    }

    hSplines[i] = CubicSpline(constantHData, 0.0, psiSpacing);
  }

  io::printData(fileData);

  BicubicSpline spline(hSplines, hStart, hSpacing);
  return spline;
}

BicubicSpline uSemiGeneration() {
  double hStart = 0;
  double hEnd = 2 * R_CNT + R_C;
  double hSpacing = (hEnd - hStart) / (USEMI_POINTS - 1);
  double xiStart = -R_C;
  double xiEnd = R_C;
  double xiSpacing = (xiEnd - xiStart) / (USEMI_POINTS - 1);

  // Switch indices for smooth extrapolation of data
  int hSwitch = ceil((2 * R_CNT + DELTA3) / hSpacing);
  int xiSwitch = ceil((R_C - DELTA3) / xiSpacing);

  std::vector<std::vector<double>> fileData(USEMI_POINTS * USEMI_POINTS,
                                            std::vector<double>(3, 0));
  std::vector<std::vector<double>> uSemiData(
      USEMI_POINTS, std::vector<double>(USEMI_POINTS, 0));
  std::vector<CubicSpline> hSplines(USEMI_POINTS, CubicSpline());

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < USEMI_POINTS; i++) {
    double h = hStart + i * hSpacing;
    for (int j = 0; j < USEMI_POINTS; j++) {
      double xi = xiStart + j * xiSpacing;
      double uSemi;
      if (i < hSwitch && j >= xiSwitch)
        uSemi = 0;
      else
        uSemi = uSemi = uExactSemiDensity(h, 0, xi, 0);

      uSemiData[i][j] = uSemi;
      fileData[i * USEMI_POINTS + j][0] = h;
      fileData[i * USEMI_POINTS + j][1] = xi;
      fileData[i * USEMI_POINTS + j][2] = uSemi;
    }
  }

// Smoothly extrapolate data in undefined domain
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < hSwitch; i++) {
    double uSemiXiSwitch = uSemiData[i][xiSwitch - 1];
    for (int j = xiSwitch; j < USEMI_POINTS; j++) {
      double uSemiHSwitch = uSemiData[hSwitch][j];
      double uSemi =
          ((hSwitch - i) * uSemiXiSwitch + (j - xiSwitch + 1) * uSemiHSwitch) /
          (hSwitch - i + j - xiSwitch + 1);
      uSemiData[i][j] = uSemi;
      fileData[i * USEMI_POINTS + j][2] = uSemi;
    }
  }

  io::printData(fileData);

// Create 1D splines at constant h
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < USEMI_POINTS; i++) {
    hSplines[i] = CubicSpline(uSemiData[i], xiStart, xiSpacing);
  }
  BicubicSpline spline(hSplines, hStart, hSpacing);
  return spline;
}

double uExactInfDensity(double h, double alpha, double xi) {
  double sinAlpha = sin(alpha);
  double axisDistance = sqrt(h * h + xi * xi * sinAlpha * sinAlpha);
  double deltaH = axisDistance - 2 * R_CNT;
  if (fabs(deltaH) >= R_C)
    return 0;

  // Compute eta integral range
  double etaMax;
  if (alpha == 0)
    etaMax = sqrt(R_C * R_C - deltaH * deltaH);
  // Conversative choice for general case
  else
    etaMax = 2 * R_CNT + R_C;

  // 3D integral with nested lambda functions
  auto phi1Integrand = [h, alpha, xi, etaMax](double phi1) -> double {
    auto etaIntegrand = [h, alpha, xi, phi1](double eta) -> double {
      auto phi2Integrand = [h, alpha, xi, eta, phi1](double phi2) -> double {
        double distance = surfaceElementDistance(h, alpha, xi, eta, phi1, phi2);
        return ljPair(distance);
      };
      return ANGLE_INTEGRATOR.integrate(phi2Integrand, 0, 2.0 * PI);
    };
    return LENGTH_INTEGRATOR.integrate(etaIntegrand, xi - etaMax, xi + etaMax);
  };
  return R_CNT * R_CNT * N_SIGMA * N_SIGMA *
         ANGLE_INTEGRATOR.integrate(phi1Integrand, 0, 2.0 * PI);
}

double uExactSemiDensity(double h, double alpha, double xi, double etaEnd) {
  double sinAlpha = sin(alpha);
  double axisDistance = sqrt(h * h + xi * xi * sinAlpha * sinAlpha);
  double deltaH = axisDistance - 2 * R_CNT;
  if (deltaH < 0)
    deltaH = 0;
  if (fabs(deltaH) >= R_C)
    return 0;

  // Compute integral range
  double etaMax;
  if (alpha == 0)
    etaMax = sqrt(R_C * R_C - deltaH * deltaH);
  // Conversative choice for general case
  else
    etaMax = 2 * R_CNT + R_C;

  if (etaEnd > xi - etaMax && etaEnd >= xi + etaMax)
    return 0;

  // 3D integral with nested lambda functions
  auto phi1Integrand = [h, alpha, xi, etaEnd, etaMax](double phi1) -> double {
    auto etaIntegrand = [h, alpha, xi, phi1](double eta) -> double {
      auto phi2Integrand = [h, alpha, xi, eta, phi1](double phi2) -> double {
        double distance = surfaceElementDistance(h, alpha, xi, eta, phi1, phi2);
        return ljPair(distance);
      };
      return ANGLE_INTEGRATOR.integrate(phi2Integrand, 0, 2.0 * PI);
    };
    if (etaEnd > xi - etaMax && etaEnd < xi + etaMax)
      return LENGTH_INTEGRATOR.integrate(etaIntegrand, etaEnd, xi + etaMax);
    else
      return LENGTH_INTEGRATOR.integrate(etaIntegrand, xi - etaMax,
                                         xi + etaMax);
  };
  return R_CNT * R_CNT * N_SIGMA * N_SIGMA *
         ANGLE_INTEGRATOR.integrate(phi1Integrand, 0, 2.0 * PI);
}

double uExactInf(double h, double alpha, double xi1, double xi2) {
  auto xiIntegrand = [h, alpha](double xi) -> double {
    return uExactInfDensity(h, alpha, xi);
  };
  return LENGTH_INTEGRATOR.integrate(xiIntegrand, xi1, xi2);
}

double uExactSemi(double h, double alpha, double xi1, double xi2,
                  double etaEnd) {
  auto xiIntegrand = [h, alpha, etaEnd](double xi) -> double {
    return uExactSemiDensity(h, alpha, xi, etaEnd);
  };
  return LENGTH_INTEGRATOR.integrate(xiIntegrand, xi1, xi2);
}

double uApproximateInf(double h, double alpha, double xi1, double xi2,
                       const CubicSpline &gammaOrthSpline,
                       const BicubicSpline &phiSpline) {
  double sinAlpha = sin(alpha);
  double gamma = gammaFunction(h, alpha, gammaOrthSpline);
  double omega = omegaFunction(alpha);
  double a = omega * sinAlpha;
  double hMax = 2 * R_CNT + R_C;

  double zeta1 = xi1 * a;
  double zeta2 = xi2 * a;
  double zetaLower = zetaMin(h);
  double zetaUpper = sqrt(hMax * hMax - h * h);
  double psi1 = (fabs(zeta1) - zetaLower) / (zetaUpper - zetaLower);
  double psi2 = (fabs(zeta2) - zetaLower) / (zetaUpper - zetaLower);
  double phi1;
  double phi2;
  if (zeta1 < 0)
    phi1 = -phiSpline(h, psi1);
  else
    phi1 = phiSpline(h, psi1);
  if (zeta2 < 0)
    phi2 = -phiSpline(h, psi2);
  else
    phi2 = phiSpline(h, psi2);

  return gamma / a * (phi2 - phi1);
}

double uApproximateSemi(double h, double alpha, double xi1, double xi2,
                        double etaEnd, const CubicSpline &gammaOrthSpline,
                        const BicubicSpline &uSemiParaSpline) {
  double sinAlpha = sin(alpha);
  double cosAlpha = cos(alpha);
  double gamma = gammaFunction(h, alpha, gammaOrthSpline);
  double omega = omegaFunction(alpha);
  double theta = thetaFunction(alpha);
  double a = omega * sinAlpha;
  double b = theta * etaEnd;
  double hSquared = h * h;

  double deltaXi = (xi2 - xi2) / (TRAPEZOIDAL_POINTS - 1);

  double uSemi = 0;

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < TRAPEZOIDAL_POINTS; i++) {
    double c = 1.0;
    if (i == 0 || i == TRAPEZOIDAL_POINTS - 1)
      c = 0.5;
    double xiBar = xi1 + i * deltaXi;
    double zetaBar = xiBar * a;
    double hBar = sqrt(hSquared + zetaBar * zetaBar);
    double etaBar = xiBar * cosAlpha - b;

    uSemi += c * uSemiParaSpline(hBar, etaBar);
  }

  uSemi *= deltaXi * gamma;
  return uSemi;
}
