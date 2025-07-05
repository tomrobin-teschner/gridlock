#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>
#include <map>
#include <memory>
#include <cmath>
#include <filesystem>
#include <chrono>
#include <sstream>

using FieldType = typename std::vector<std::vector<double>>;

#include "src/meshLooper.hpp"
#include "src/postProcessor.hpp"
#include "src/boundaryConditions.hpp"
#include "src/ui.hpp"

#include "Eigen/Eigen"
#include "nlohmann/json.hpp"

// helper functions
auto field = [](int numX, int numY) {
  return FieldType(numX, std::vector<double>(numY));
};

auto map2Dto1D = [](int i, int j, int numX, int numY) { return j * numX + i; };
auto map1Dto2D = [](int index, int numX, int numY) { return std::make_tuple(index % numX, index / numX); };

int main(int argv, char* argc[]) {

  // initialise UI
  UI ui;
  ui.draw(42, 2, "Gridlock 1.0");
  ui.draw(42, 3, "A 2D structured CFD solver for incompressible flows");
  ui.draw(42, 4, "Tom-Robin Teschner");

  // delete the output folder and then recreate it
  std::filesystem::remove_all("output/");
  std::filesystem::create_directories("output/");

  // extract second argument from command line
  if (argv != 2) {
    ui.draw(0, 7,"Error: I was expecting an input file to process, Correct usage:");
    ui.draw(0, 8, std::string(argc[0]) + " <filename>");
    ui.draw(0, 9, "<filename> is the JSON input file, located in the input/ folder");
    return -1;
  }
  std::string filename = argc[1];
  ui.draw(0, 7, "Case setup defined in: ");
  ui.draw(35, 7, filename);

  // reading JSON input file
  std::ifstream inputFile(filename);
  if (!inputFile.is_open()) {
    ui.draw(0, 8, "Error: Could not open input file: " + filename);
    ui.draw(0, 9, "Press any key to continue!");
    ui.block();
    return -1;
  }

  // in case file could be read, parse it from JSON
  nlohmann::json parameters;
  inputFile >> parameters;
  inputFile.close();

  // Mesh parameters
  int numX = parameters["mesh"]["numX"];
  int numY = parameters["mesh"]["numY"];
  double Lx = parameters["mesh"]["Lx"];
  double Ly = parameters["mesh"]["Ly"];
  double dx = Lx / (numX - 1); double dy = Ly / (numY - 1);
  int numGhostPoints = 1;
  int totalSizeX = numX + 2 * numGhostPoints;
  int totalSizeY = numY + 2 * numGhostPoints;
  
  // Time stepping parameters
  int timeSteps = parameters["time"]["timeSteps"];
  double CFL = parameters["time"]["CFL"];
  int outputFrequency = parameters["output"]["outputFrequency"];

  // residual and convergence checking
  double resU = 0.0; double resV = 0.0; double resP = 0.0;
  double resUNorm = 1.0; double resVNorm = 1.0; double resPNorm = 1.0;
  double epsU = parameters["convergence"]["u"];
  double epsV = parameters["convergence"]["v"];
  double epsP = parameters["convergence"]["p"];

  // input for the linear solver
  Eigen::VectorXd bp((numX) * (numY));
  Eigen::VectorXd xp((numX) * (numY));
  Eigen::SparseMatrix<double> Ap((numX) * (numY), (numX) * (numY));

  Eigen::VectorXd bu((numX) * (numY));
  Eigen::VectorXd xu((numX) * (numY));
  Eigen::SparseMatrix<double> Au((numX) * (numY), (numX) * (numY));

  Eigen::VectorXd bv((numX) * (numY));
  Eigen::VectorXd xv((numX) * (numY));
  Eigen::SparseMatrix<double> Av((numX) * (numY), (numX) * (numY));

  // physical properties
  double nu = parameters["fluid"]["nu"];

  // Mesh generation
  auto x = field(totalSizeX, totalSizeY);
  auto y = field(totalSizeX, totalSizeY);

  // mesh looper to facilitate looping over mesh
  MeshLooper looper(numX, numY, numGhostPoints);

  looper.loopAll([dx, dy, numGhostPoints, &x, &y](int i, int j) {
    x[i][j] = dx * i - dx * numGhostPoints;
    y[i][j] = dy * j - dy * numGhostPoints;
  });

  // solution vectors
  auto u = std::make_shared<std::vector<std::vector<double>>>(totalSizeX, std::vector<double>(totalSizeY));
  auto v = std::make_shared<std::vector<std::vector<double>>>(totalSizeX, std::vector<double>(totalSizeY));
  auto p = std::make_shared<std::vector<std::vector<double>>>(totalSizeX, std::vector<double>(totalSizeY));

  auto uOld = std::make_shared<std::vector<std::vector<double>>>(totalSizeX, std::vector<double>(totalSizeY));
  auto vOld = std::make_shared<std::vector<std::vector<double>>>(totalSizeX, std::vector<double>(totalSizeY));
  auto pOld = std::make_shared<std::vector<std::vector<double>>>(totalSizeX, std::vector<double>(totalSizeY));

  // Initial conditions
  looper.loopAll([u, v, p](int i, int j) {
    (*u)[i][j] = 0.0; (*v)[i][j] = 0.0; (*p)[i][j] = 0.0;
  });

  // Instantiate boundary condition class
  BoundaryConditions bc(u, v, p, x, y, looper, parameters);
  bc.applyBCs("u", u);
  bc.applyBCs("v", v);
  bc.applyBCs("p", p);

  // create post processing object
  PostProcessor output("solution", numX, numY, numGhostPoints, x, y);
  output.registerField("u", u);
  output.registerField("v", v);
  output.registerField("p", p);

  // create residual file for plotting later
  std::ofstream residualsFile;
  residualsFile.open("output/residual.txt");
  assert(residualsFile.is_open() && "Could not open residual file!");
  residualsFile << "time,u,v,p" << std::endl;

  // loop over time
  double totalTime = 0.0;
  int t = 0;
  auto timeStart = std::chrono::high_resolution_clock::now();
  for (t=0; t<timeSteps; ++t) {
    
    // create deep copy of old solution
    *uOld = *u; *vOld = *v; *pOld = *p;

    // determien stable timestep
    double dt = std::numeric_limits<double>::max();
    looper.loopWithBoundaries([CFL, nu, dx, dy, &dt, &x, &y, uOld, vOld](int i, int j) {
      auto minSpacing = std::min(dx, dy);
      auto velocityMag = std::sqrt((*uOld)[i][j] * (*uOld)[i][j] + (*vOld)[i][j] * (*vOld)[i][j]);
      auto inviscidTimeStep = CFL * minSpacing / velocityMag;
      auto viscousTimeStep = (0.5 * CFL) * minSpacing * minSpacing / nu;
      auto tempDt = std::min(inviscidTimeStep, viscousTimeStep);
      if (tempDt < dt) dt = tempDt;
    });

    // solve the u-momentum equations
    Au.setZero();
    looper.loopWithBoundaries([&](int i, int j) {
      // create zero-based inmdices for i and j direction
      auto idx = i - numGhostPoints;
      auto jdx = j - numGhostPoints;

      // map 2D indices to 1D indices. Subtract 2 from numX|Y to ignore boundary points
      auto ic  = map2Dto1D(idx, jdx, numX, numY);
      auto ip1 = map2Dto1D(idx + 1, jdx, numX, numY);
      auto im1 = map2Dto1D(idx - 1, jdx, numX, numY);
      auto jp1 = map2Dto1D(idx, jdx + 1, numX, numY);
      auto jm1 = map2Dto1D(idx, jdx - 1, numX, numY);

      // set up right-hand side
      bu(ic) = (*uOld)[i][j] / dt;
      
      // compute coefficients for matrix
      auto umax = std::max((*uOld)[i][j], 0.0) / dx;
      auto umin = std::min((*uOld)[i][j], 0.0) / dx;
      auto vmax = std::max((*vOld)[i][j], 0.0) / dy;
      auto vmin = std::min((*vOld)[i][j], 0.0) / dy;
      auto nudx2 = nu / std::pow(dx, 2);
      auto nudy2 = nu / std::pow(dy, 2);

      // Construct coefficient matrix
      Au.coeffRef(ic, ic) = (1.0 / dt) + umax - umin + vmax - vmin + 2.0 * nudx2 + 2.0 * nudy2;

      // west
      if (idx == 0) {
        auto westBCType = parameters["boundaries"]["west"]["u"][0];
        auto westBCValue = parameters["boundaries"]["west"]["u"][1];
        if (westBCType == "dirichlet") {
          bu(ic) -= 2.0 * westBCValue * (- umax - nudx2);
        } else if (westBCType == "neumann") {
          Au.coeffRef(ic, ip1) = 2.0 * (- umax - nudx2);
          bu(ic) += 2.0 * westBCValue / dx;
        }

      // east
      } else if (idx == numX - 1) {
        auto eastBCType = parameters["boundaries"]["east"]["u"][0];
        auto eastBCValue = parameters["boundaries"]["east"]["u"][1];
        if (eastBCType == "dirichlet") {
          bu(ic) -= 2.0 * eastBCValue * (umin - nudx2);
        } else if (eastBCType == "neumann") {
          Au.coeffRef(ic, im1) = 2.0 * (umin - nudx2);
          bu(ic) -= 2.0 * eastBCValue / dx;
        }
          
      } else {
        Au.coeffRef(ic, ip1) = umin - nudx2;
        Au.coeffRef(ic, im1) = - umax - nudx2;
      }

      // south
      if (jdx == 0) {
        auto southBCType = parameters["boundaries"]["south"]["u"][0];
        auto southBCValue = parameters["boundaries"]["south"]["u"][1];
        if (southBCType == "dirichlet") {
          bu(ic) -= 2.0 * southBCValue * (- vmax - nudy2);
        } else if (southBCType == "neumann") {
          Au.coeffRef(ic, jp1) = 2.0 * (- vmax - nudy2);
          bu(ic) += 2.0 * southBCValue / dy;
        }
      
      // north
      } else if (jdx == numY - 1) {
        auto northBCType = parameters["boundaries"]["north"]["u"][0];
        auto northBCValue = parameters["boundaries"]["north"]["u"][1];
        if (northBCType == "dirichlet") {
          bu(ic) -= 2.0 * northBCValue * (vmin - nudy2);
        } else if (northBCType == "neumann") {
          Au.coeffRef(ic, jm1) = 2.0 * (vmin - nudy2);
          bu(ic) -= 2.0 * northBCValue / dy;
        }  
        
      } else {
        Au.coeffRef(ic, jp1) = vmin - nudy2;
        Au.coeffRef(ic, jm1) = - vmax - nudy2;
      }
    });

    // initial residual
    auto initResU = bu - Au * xu;
    resU = initResU.lpNorm<2>();

    // normalise residuals
    if (t < 2) if (resU > 0.0) resUNorm = resU;
    resU /= resUNorm;

    // solve pressure poisson solver with the preconditioned Conjugate Gradient method
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstabU;
    bicgstabU.compute(Au);
    bicgstabU.setMaxIterations(parameters["linearSolver"]["u"]["maxIterations"]);
    bicgstabU.setTolerance(parameters["linearSolver"]["u"]["tolerance"]);

    xu = bicgstabU.solve(bu);

    // update pressure
    looper.loopWithBoundaries([&](int i, int j) {
      auto idx = i - numGhostPoints;
      auto jdx = j - numGhostPoints;
      (*u)[i][j] = xu(map2Dto1D(idx, jdx, numX, numY));
    });
    
    bc.applyBCs("u", u);

    // solve the v-momentum equations
    Av.setZero();
    looper.loopWithBoundaries([&](int i, int j) {
      // create zero-based inmdices for i and j direction
      auto idx = i - numGhostPoints;
      auto jdx = j - numGhostPoints;

      // map 2D indices to 1D indices. Subtract 2 from numX|Y to ignore boundary points
      auto ic  = map2Dto1D(idx, jdx, numX, numY);
      auto ip1 = map2Dto1D(idx + 1, jdx, numX, numY);
      auto im1 = map2Dto1D(idx - 1, jdx, numX, numY);
      auto jp1 = map2Dto1D(idx, jdx + 1, numX, numY);
      auto jm1 = map2Dto1D(idx, jdx - 1, numX, numY);

      // set up right-hand side
      bv(ic) = (*vOld)[i][j] / dt;
      
      // compute coefficients for matrix
      auto umax = std::max((*uOld)[i][j], 0.0) / dx;
      auto umin = std::min((*uOld)[i][j], 0.0) / dx;
      auto vmax = std::max((*vOld)[i][j], 0.0) / dy;
      auto vmin = std::min((*vOld)[i][j], 0.0) / dy;
      auto nudx2 = nu / std::pow(dx, 2);
      auto nudy2 = nu / std::pow(dy, 2);

      // Construct coefficient matrix
      Av.coeffRef(ic, ic) = (1.0 / dt) + umax - umin + vmax - vmin + 2.0 * nudx2 + 2.0 * nudy2;

      // west
      if (idx == 0) {
        auto westBCType = parameters["boundaries"]["west"]["v"][0];
        auto westBCValue = parameters["boundaries"]["west"]["v"][1];
        if (westBCType == "dirichlet") {
          bv(ic) -= 2.0 * westBCValue * (- umax - nudx2);
        } else if (westBCType == "neumann") {
          Av.coeffRef(ic, ip1) = 2.0 * (- umax - nudx2);
          bv(ic) += 2.0 * westBCValue / dx;
        }

      // east
      } else if (idx == numX - 1) {
        auto eastBCType = parameters["boundaries"]["east"]["v"][0];
        auto eastBCValue = parameters["boundaries"]["east"]["v"][1];
        if (eastBCType == "dirichlet") {
          bv(ic) -= 2.0 * eastBCValue * (umin - nudx2);
        } else if (eastBCType == "neumann") {
          Av.coeffRef(ic, im1) = 2.0 * (umin - nudx2);
          bv(ic) -= 2.0 * eastBCValue / dx;
        }
          
      } else {
        Av.coeffRef(ic, ip1) = umin - nudx2;
        Av.coeffRef(ic, im1) = - umax - nudx2;
      }

      // south
      if (jdx == 0) {
        auto southBCType = parameters["boundaries"]["south"]["v"][0];
        auto southBCValue = parameters["boundaries"]["south"]["v"][1];
        if (southBCType == "dirichlet") {
          bv(ic) -= 2.0 * southBCValue * (- vmax - nudy2);
        } else if (southBCType == "neumann") {
          Av.coeffRef(ic, jp1) = 2.0 * (- vmax - nudy2);
          bv(ic) += 2.0 * southBCValue / dy;
        }
      
      // north
      } else if (jdx == numY - 1) {
        auto northBCType = parameters["boundaries"]["north"]["v"][0];
        auto northBCValue = parameters["boundaries"]["north"]["v"][1];
        if (northBCType == "dirichlet") {
          bv(ic) -= 2.0 * northBCValue * (vmin - nudy2);
        } else if (northBCType == "neumann") {
          Av.coeffRef(ic, jm1) = 2.0 * (vmin - nudy2);
          bv(ic) -= 2.0 * northBCValue / dy;
        }  
        
      } else {
        Av.coeffRef(ic, jp1) = vmin - nudy2;
        Av.coeffRef(ic, jm1) = - vmax - nudy2;
      }
    });

    // initial residual
    auto initResV = bv - Av * xv;
    resV = initResV.lpNorm<2>();

    // normalise residuals
    if (t < 2) if (resV > 0.0) resVNorm = resV;
    resV /= resVNorm;

    // solve pressure poisson solver with the preconditioned Conjugate Gradient method
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstabV;
    bicgstabV.compute(Av);
    bicgstabV.setMaxIterations(parameters["linearSolver"]["v"]["maxIterations"]);
    bicgstabV.setTolerance(parameters["linearSolver"]["v"]["tolerance"]);

    xv = bicgstabV.solve(bv);

    // update pressure
    looper.loopWithBoundaries([&](int i, int j) {
      auto idx = i - numGhostPoints;
      auto jdx = j - numGhostPoints;
      (*v)[i][j] = xv(map2Dto1D(idx, jdx, numX, numY));
    });
    
    bc.applyBCs("v", v);

    // // solve momentum equations
    // looper.loopInterior([&](int i, int j) {
    //   // second-order central scheme for diffusive fluxes
    //   auto d2udx2 = ((*uOld)[i + 1][j] - 2.0 * (*uOld)[i][j] + (*uOld)[i - 1][j]) / dx / dx;
    //   auto d2udy2 = ((*uOld)[i][j + 1] - 2.0 * (*uOld)[i][j] + (*uOld)[i][j - 1]) / dy / dy;
    //   auto d2vdx2 = ((*vOld)[i + 1][j] - 2.0 * (*vOld)[i][j] + (*vOld)[i - 1][j]) / dx / dx;
    //   auto d2vdy2 = ((*vOld)[i][j + 1] - 2.0 * (*vOld)[i][j] + (*vOld)[i][j - 1]) / dy / dy;

    //   // first-order upwind for convective fluxes
    //   auto uplus  = std::max((*uOld)[i][j], 0.0);
    //   auto uminus = std::min((*uOld)[i][j], 0.0);
    //   auto vplus  = std::max((*vOld)[i][j], 0.0);
    //   auto vminus = std::min((*vOld)[i][j], 0.0);

    //   auto dudx_backward = ((*uOld)[i][j] - (*uOld)[i - 1][j]) / dx;
    //   auto dudx_forward  = ((*uOld)[i + 1][j] - (*uOld)[i][j]) / dx;
    //   auto dvdx_backward = ((*vOld)[i][j] - (*vOld)[i - 1][j]) / dx;
    //   auto dvdx_forward  = ((*vOld)[i + 1][j] - (*vOld)[i][j]) / dx;

    //   auto dudy_backward = ((*uOld)[i][j] - (*uOld)[i][j - 1]) / dy;
    //   auto dudy_forward  = ((*uOld)[i][j + 1] - (*uOld)[i][j]) / dy;
    //   auto dvdy_backward = ((*vOld)[i][j] - (*vOld)[i][j - 1]) / dy;
    //   auto dvdy_forward  = ((*vOld)[i][j + 1] - (*vOld)[i][j]) / dy;

    //   auto duudx = uplus * dudx_backward + uminus * dudx_forward;
    //   auto dvudy = vplus * dudy_backward + vminus * dudy_forward;
    //   auto duvdx = uplus * dvdx_backward + uminus * dvdx_forward;
    //   auto dvvdy = vplus * dvdy_backward + vminus * dvdy_forward;
      
    //   // solve momentum equations
    //   (*u)[i][j] = (*uOld)[i][j] + dt * (nu * (d2udx2 + d2udy2) - duudx - dvudy);
    //   (*v)[i][j] = (*vOld)[i][j] + dt * (nu * (d2vdx2 + d2vdy2) - duvdx - dvvdy);
    // });
   
    // // pressure projection step
    // bc.applyBCs("u", u);
    // bc.applyBCs("v", v);

    // set up right-hand side and coefficient matrix
    looper.loopWithBoundaries([&](int i, int j) {
      // create zero-based inmdices for i and j direction
      auto idx = i - numGhostPoints;
      auto jdx = j - numGhostPoints;

      // map 2D indices to 1D indices. Subtract 2 from numX|Y to ignore boundary points
      auto ic  = map2Dto1D(idx, jdx, numX, numY);
      auto ip1 = map2Dto1D(idx + 1, jdx, numX, numY);
      auto im1 = map2Dto1D(idx - 1, jdx, numX, numY);
      auto jp1 = map2Dto1D(idx, jdx + 1, numX, numY);
      auto jm1 = map2Dto1D(idx, jdx - 1, numX, numY);

      // set up right-hand side
      auto dudx = ((*u)[i + 1][j] - (*u)[i - 1][j]) / (2.0 * dx);
      auto dvdy = ((*v)[i][j + 1] - (*v)[i][j - 1]) / (2.0 * dy);
      auto div = (dudx + dvdy);
      bp(ic) = div / dt;

      // set up matrix coefficients
      auto dx2 = 1.0 / std::pow(dx, 2);
      auto dy2 = 1.0 / std::pow(dy, 2);

      // Construct coefficient matrix
      Ap.coeffRef(ic, ic) = - 2.0 * (dx2 + dy2);

      // west
      if (idx == 0) {
        auto westBCType = parameters["boundaries"]["west"]["p"][0];
        auto westBCValue = parameters["boundaries"]["west"]["p"][1];
        if (westBCType == "dirichlet") {
          bp(ic) -= 2.0 * westBCValue * dx2;
        } else if (westBCType == "neumann") {
          Ap.coeffRef(ic, ip1) = 2.0 * dx2;
          bp(ic) += 2.0 * westBCValue / dx;
        }

      // east
      } else if (idx == numX - 1) {
        auto eastBCType = parameters["boundaries"]["east"]["p"][0];
        auto eastBCValue = parameters["boundaries"]["east"]["p"][1];
        if (eastBCType == "dirichlet") {
          bp(ic) -= 2.0 * eastBCValue * dx2;
        } else if (eastBCType == "neumann") {
          Ap.coeffRef(ic, im1) = 2.0 * dx2;
          bp(ic) -= 2.0 * eastBCValue / dx;
        }
          
      } else {
        Ap.coeffRef(ic, ip1) = dx2;
        Ap.coeffRef(ic, im1) = dx2;
      }

      // south
      if (jdx == 0) {
        auto southBCType = parameters["boundaries"]["south"]["p"][0];
        auto southBCValue = parameters["boundaries"]["south"]["p"][1];
        if (southBCType == "dirichlet") {
          bp(ic) -= 2.0 * southBCValue * dy2;
        } else if (southBCType == "neumann") {
          Ap.coeffRef(ic, jp1) = 2.0 * dy2;
          bp(ic) += 2.0 * southBCValue / dy;
        }
      
      // north
      } else if (jdx == numY - 1) {
        auto northBCType = parameters["boundaries"]["north"]["p"][0];
        auto northBCValue = parameters["boundaries"]["north"]["p"][1];
        if (northBCType == "dirichlet") {
          bp(ic) -= 2.0 * northBCValue * dy2;
        } else if (northBCType == "neumann") {
          Ap.coeffRef(ic, jm1) = 2.0 * dy2;
          bp(ic) -= 2.0 * northBCValue / dy;
        }  
        
      } else {
        Ap.coeffRef(ic, jp1) = dy2;
        Ap.coeffRef(ic, jm1) = dy2;
      }
    });

    // initial residual
    auto initResP = bp - Ap * xp;
    resP = initResP.lpNorm<2>();

    // normalise residuals
    if (t < 2) if (resP > 0.0) resPNorm = resP;
    resP /= resPNorm;

    // solve pressure poisson solver with the preconditioned Conjugate Gradient method
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstabP;
    bicgstabP.compute(Ap);
    bicgstabP.setMaxIterations(parameters["linearSolver"]["p"]["maxIterations"]);
    bicgstabP.setTolerance(parameters["linearSolver"]["p"]["tolerance"]);

    xp = bicgstabP.solve(bp);

    // update pressure
    looper.loopWithBoundaries([&](int i, int j) {
      auto idx = i - numGhostPoints;
      auto jdx = j - numGhostPoints;
      (*p)[i][j] = xp(map2Dto1D(idx, jdx, numX, numY));
    });
    
    bc.applyBCs("p", p);

    // update velocity
    looper.loopInterior([&](int i, int j) {
      (*u)[i][j] = (*u)[i][j] - dt * ((*p)[i + 1][j] - (*p)[i - 1][j]) / (2.0 * dx);
      (*v)[i][j] = (*v)[i][j] - dt * ((*p)[i][j + 1] - (*p)[i][j - 1]) / (2.0 * dy);
    });

    bc.applyBCs("u", u);
    bc.applyBCs("v", v);

    if (outputFrequency != -1 && t % outputFrequency == 0) 
      output.write(t);

    totalTime = totalTime + dt;

    auto timeEndIteration = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(timeEndIteration - timeStart);
    auto averageDuration = 0.001 * duration.count() / (t + 1);
    auto timeRemaining = (timeSteps - (t + 1)) * averageDuration;

    int hh = std::floor(timeRemaining / 60.0 / 60.0);
    int mm = std::floor(timeRemaining / 60.0 - hh * 60.0);
    int ss = std::floor(timeRemaining - hh * 60.0 * 60.0 - mm * 60.0);

    std::stringstream pResSS; pResSS << std::scientific << std::setprecision(2) << resP;
    std::stringstream pIter; pIter << std::fixed << std::setw(4) << std::setfill(' ')  << bicgstabP.iterations();
    std::stringstream uResSS; uResSS << std::scientific << std::setprecision(2) << resU;
    std::stringstream uIter; uIter << std::fixed << std::setw(4) << std::setfill(' ')  << bicgstabU.iterations();
    std::stringstream vResSS; vResSS << std::scientific << std::setprecision(2) << resV;
    std::stringstream vIter; vIter << std::fixed << std::setw(4) << std::setfill(' ')  << bicgstabV.iterations();
    std::stringstream hhSS; hhSS << std::fixed << std::setw(2) << std::setfill('0')  << hh;
    std::stringstream mmSS; mmSS << std::fixed << std::setw(2) << std::setfill('0')  << mm;
    std::stringstream ssSS; ssSS << std::fixed << std::setw(2) << std::setfill('0')  << ss;
    std::stringstream averageDurationSS; averageDurationSS << std::scientific << std::setprecision(2) << averageDuration;
    std::stringstream totalTimeSS; totalTimeSS << std::scientific << std::setprecision(2) << totalTime;

    ui.draw(0, 8, "Time step:  " + std::to_string(t + 1) + " / " + std::to_string(timeSteps));
    ui.draw(35, 8, "Total time: " + totalTimeSS.str() + " [s]");
    ui.draw(0, 9, "Time / iteration: " + averageDurationSS.str() + " [s]");
    ui.draw(35, 9, "Remaining time: " + hhSS.str() + ":" + mmSS.str() + ":" + ssSS.str() + " [hh:mm:ss]");
    ui.draw(0, 10, "u residual: " + uResSS.str());
    ui.draw(35, 10, "u iterations: " + uIter.str());
    ui.draw(0, 11, "v residual: " + vResSS.str());
    ui.draw(35, 11, "v iterations: " + vIter.str());
    ui.draw(0, 12, "p residual: " + pResSS.str());
    ui.draw(35, 12, "p iterations: " + pIter.str());

    // write residuals to file
    residualsFile << std::fixed << t + 1 << ",";
    residualsFile << std::scientific << std::setprecision(5) << resU << "," << resV << "," << resP << std::endl;

    // check convergence
    if (resU < epsU && resV < epsV && resP < epsP) {
      ui.draw(0, 14, "Solution converged. Press any key to continue!");
      break;
    }
  }

  // close residual file
  residualsFile.close();

  // if (parameters["output"]["createRestartFile"] == true) {
  //   // open file in binary mode
  //   std::ofstream restartFile("output/restart.bin", std::ios::binary);

  //   // write mesh
  //   restartFile.write((char*)&numX, sizeof(int));
  //   restartFile.write((char*)&numY, sizeof(int));
  //   restartFile.write((char*)&numGhostPoints, sizeof(int));
  //   restartFile.write((char*)&t, sizeof(int));
  //   for (int i = 0; i < numX + 2 * numGhostPoints; ++i) {
  //     for (int j = 0; j < numY + 2 * numGhostPoints; ++j) {
  //       restartFile.write((char*)&x[i][j], sizeof(double));
  //       restartFile.write((char*)&y[i][j], sizeof(double));
  //       restartFile.write((char*)u[i][j], sizeof(double));
  //       restartFile.write((char*)v[i][j], sizeof(double));
  //       restartFile.write((char*)p[i][j], sizeof(double));
  //     }
  //   }
  // }

  // write output to file
  output.write();

  // draw end message if simulation has not converged
  if (epsU > resU || epsV > resV || epsP > resP)
    ui.draw(0, 14, "Simulation finished. Press any key to continue!");
  ui.block();

  return 0;
}