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

#include "src/meshLooper/meshLooper.hpp"
#include "src/boundaryConditions.hpp"
#include "src/timeStep/timeStep.hpp"
#include "src/infrastructure/parameterFile/parameterFile.hpp"
#include "src/infrastructure/ui/ui.hpp"
#include "src/infrastructure/stopWatch/stopWatch.hpp"
#include "src/infrastructure/utilities/data.hpp"
#include "src/postProcessing/postProcessing.hpp"
#include "src/fieldArray/fieldArray.hpp"
#include "src/residuals/residuals.hpp"

#include "Eigen/Eigen"
#include "nlohmann/json.hpp"

// helper functions
auto field = [](int numX, int numY) {
  return FieldType(numX, std::vector<double>(numY));
};

auto map2Dto1D = [](int i, int j, int numX, int numY) { return j * numX + i; };
auto map1Dto2D = [](int index, int numX, int numY) { return std::make_tuple(index % numX, index / numX); };

int main(int argv, char* argc[]) {

  // Read parameter files
  ParameterFile parameterFile("input/solver.json");
  auto solverParameters = parameterFile.getSolverParameters();
  auto meshParameters = parameterFile.getMeshParameters();
  auto bcParameters = parameterFile.getBoundaryConditionsParameters();

  // check if problem is a fully neumann boundary value problem for the pressure
  auto pWestBC = bcParameters["boundaries"]["west"]["p"][0];
  auto pEastBC = bcParameters["boundaries"]["east"]["p"][0];
  auto pSouthBC = bcParameters["boundaries"]["south"]["p"][0];
  auto pNorthBC = bcParameters["boundaries"]["north"]["p"][0];
  bool fullyNeumann = false;
  if (pWestBC == "neumann" && pEastBC == "neumann" && pSouthBC == "neumann" && pNorthBC == "neumann")
    fullyNeumann = true;

  // Mesh parameters
  int numX = meshParameters["mesh"]["numX"];
  int numY = meshParameters["mesh"]["numY"];
  double Lx = meshParameters["mesh"]["Lx"];
  double Ly = meshParameters["mesh"]["Ly"];

  double dx = Lx / (numX - 1); double dy = Ly / (numY - 1);
  int numGhostPoints = 1;
  int totalSizeX = numX + 2 * numGhostPoints;
  int totalSizeY = numY + 2 * numGhostPoints;
  
  // Time stepping parameters
  int timeSteps = solverParameters["time"]["timeSteps"];
  double CFL = solverParameters["time"]["CFL"];
  int outputFrequency = solverParameters["output"]["outputFrequency"];

  // residual and convergence checking
  int uIter = 0; int vIter = 0; int pIter = 0;
  double epsU = solverParameters["convergence"]["u"];
  double epsV = solverParameters["convergence"]["v"];
  double epsP = solverParameters["convergence"]["p"];

  // set up picard linearisation
  int maxPicardIterations = solverParameters["linearisation"]["maxPicardIterations"];
  double picardToleranceU = solverParameters["linearisation"]["tolerance"]["u"];
  double picardToleranceV = solverParameters["linearisation"]["tolerance"]["v"];
  double picardToleranceP = solverParameters["linearisation"]["tolerance"]["p"];

  // initialise UI
  UI ui(timeSteps, maxPicardIterations);
  ui.createSkeleton();

  // input for the linear solver
  Eigen::VectorXd bp((numX) * (numY));
  Eigen::VectorXd xp((numX) * (numY));
  Eigen::SparseMatrix<double> Ap((numX) * (numY), (numX) * (numY));
  double alphaP = solverParameters["linearSolver"]["p"]["underRelaxation"];
  
  Eigen::VectorXd bu((numX) * (numY));
  Eigen::VectorXd xu((numX) * (numY));
  Eigen::SparseMatrix<double> Au((numX) * (numY), (numX) * (numY));
  double alphaU = solverParameters["linearSolver"]["u"]["underRelaxation"];
  
  Eigen::VectorXd bv((numX) * (numY));
  Eigen::VectorXd xv((numX) * (numY));
  Eigen::SparseMatrix<double> Av((numX) * (numY), (numX) * (numY));
  double alphaV = solverParameters["linearSolver"]["v"]["underRelaxation"];

  // physical properties
  double nu = solverParameters["fluid"]["nu"];

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
  FieldArray u(totalSizeX, totalSizeY);
  FieldArray v(totalSizeX, totalSizeY);
  FieldArray p(totalSizeX, totalSizeY);

  FieldArray uOld(totalSizeX, totalSizeY);
  FieldArray vOld(totalSizeX, totalSizeY);
  FieldArray pOld(totalSizeX, totalSizeY);

  FieldArray uPicardOld(totalSizeX, totalSizeY);
  FieldArray vPicardOld(totalSizeX, totalSizeY);
  FieldArray pPicardOld(totalSizeX, totalSizeY);
  
  // create post processing object
  auto outputFileName = solverParameters["output"]["filename"];
  PostProcessing output(outputFileName, numX, numY, numGhostPoints, x, y);
  output.registerField("u", &u);
  output.registerField("v", &v);
  output.registerField("p", &p);

  Residuals outerResiduals(&u, &v, &p, looper, true);
  Residuals picardResiduals(&u, &v, &p, looper);

  // Instantiate boundary condition class
  BoundaryConditions bc(x, y, looper, bcParameters);
  bc.applyBCs("u", u);
  bc.applyBCs("v", v);
  bc.applyBCs("p", p);

  // create high-resolution stop watch
  StopWatch timer(timeSteps);
  
  // Create time step calculation object
  TimeStep timeStep(solverParameters, looper);


  // loop over time
  double totalTime = 0.0;
  int t = 0;
  for (t = 0; t<timeSteps; ++t) {
    
    // create deep copy of old solution
    uOld = u; vOld = v; pOld = p;

    // initialise residuals
    outerResiduals.init();

    // determine stable timestep
    auto dt = timeStep.getTimeStep(u, v, dx, dy);

    // get timings
    auto [elapsedHH, elapsedMM, elapsedSS] = timer.elapsed();
    auto [remainingHH, remainingMM, remainingSS] = timer.remaining(t);

    // update time information of what is already available
    ui.updateTimestatistics(t + 1, CFL, dt, totalTime, elapsedHH, elapsedMM, elapsedSS,
      remainingHH, remainingMM, remainingSS);    

    // outer (picard) loop
    for (int k = 0; k < maxPicardIterations; ++k) {

      // create deep copy of old solution
      uPicardOld = u; vPicardOld = v; pPicardOld = p;

      // auto resUPicardStart = 0.0;
      // looper.loopWithBoundaries([&resUPicardStart, &u](int i, int j) {
      //   resUPicardStart += std::fabs(u[i, j]);
      // });

      picardResiduals.init();

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
        
        // compute coefficients for matrix
        auto umax = std::max(uPicardOld[i, j], 0.0) / dx;
        auto umin = std::min(uPicardOld[i, j], 0.0) / dx;
        auto vmax = std::max(vPicardOld[i, j], 0.0) / dy;
        auto vmin = std::min(vPicardOld[i, j], 0.0) / dy;
        auto nudx2 = nu / std::pow(dx, 2);
        auto nudy2 = nu / std::pow(dy, 2);
        
        auto aP = (1.0 / dt) + umax - umin + vmax - vmin + 2.0 * nudx2 + 2.0 * nudy2;
        auto aE = umin - nudx2;
        auto aW = - umax - nudx2;
        auto aN = vmin - nudy2;
        auto aS = - vmax - nudy2; 

        // Construct coefficient matrix
        Au.coeffRef(ic, ic) = aP;

        // set up right-hand side
        bu(ic) = uOld[i, j] / dt;

        // west
        if (idx == 0) {
          auto westBCType = bcParameters["boundaries"]["west"]["u"][0];
          auto westBCValue = bcParameters["boundaries"]["west"]["u"][1];
          if (westBCType == "dirichlet") {
            Au.coeffRef(ic, ip1) = aE - aW;
            bu(ic) += 2.0 * westBCValue * umax + 2.0 * westBCValue * nudx2;
          } else if (westBCType == "neumann") {
            Au.coeffRef(ic, ip1) = aE + aW;
            bu(ic) += -2.0 * westBCValue * umax * dx - 2.0 * nu * westBCValue / dx;
          }

        // east
        } else if (idx == numX - 1) {
          auto eastBCType = bcParameters["boundaries"]["east"]["u"][0];
          auto eastBCValue = bcParameters["boundaries"]["east"]["u"][1];
          if (eastBCType == "dirichlet") {
            Au.coeffRef(ic, im1) = aW - aE;
            bu(ic) += -2.0 * eastBCValue * umin + 2.0 * eastBCValue * nudx2;
          } else if (eastBCType == "neumann") {
            Au.coeffRef(ic, im1) = aW + aE;
            bu(ic) += - 2.0 * eastBCValue * umin * dx+ 2.0 * nu * eastBCValue / dx;
          }
            
        } else {
          Au.coeffRef(ic, ip1) = aE;
          Au.coeffRef(ic, im1) = aW;
        }

        // south
        if (jdx == 0) {
          auto southBCType = bcParameters["boundaries"]["south"]["u"][0];
          auto southBCValue = bcParameters["boundaries"]["south"]["u"][1];
          if (southBCType == "dirichlet") {
            Au.coeffRef(ic, jp1) = aN - aS;
            bu(ic) += 2.0 * southBCValue * vmax + 2.0 * southBCValue * nudy2;
          } else if (southBCType == "neumann") {
            Au.coeffRef(ic, jp1) = aN + aS;
            bu(ic) += -2.0 * southBCValue * vmax * dy- 2.0 * nu * southBCValue / dy;
          }
        
        // north
        } else if (jdx == numY - 1) {
          auto northBCType = bcParameters["boundaries"]["north"]["u"][0];
          auto northBCValue = bcParameters["boundaries"]["north"]["u"][1];
          if (northBCType == "dirichlet") {
            Au.coeffRef(ic, jm1) = aS - aN;
            bu(ic) += -2.0 * northBCValue * vmin + 2.0 * northBCValue * nudy2;
          } else if (northBCType == "neumann") {
            Au.coeffRef(ic, jm1) = aS + aN;
            bu(ic) += - 2.0 * northBCValue * vmin * dy + 2.0 * nu * northBCValue / dy;
          }
          
        } else {
          Au.coeffRef(ic, jp1) = aN;
          Au.coeffRef(ic, jm1) = aS;
        }
      });

      // solve pressure poisson solver with the preconditioned Conjugate Gradient method
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstabU;
      bicgstabU.compute(Au);
      bicgstabU.setMaxIterations(solverParameters["linearSolver"]["u"]["maxIterations"]);
      bicgstabU.setTolerance(solverParameters["linearSolver"]["u"]["tolerance"]);

      xu = bicgstabU.solve(bu);
      uIter = bicgstabU.iterations();
      
      // update u-velocity with under-relaxation applied
      looper.loopWithBoundaries([&](int i, int j) {
        auto idx = i - numGhostPoints;
        auto jdx = j - numGhostPoints;
        u[i, j] = alphaU * xu(map2Dto1D(idx, jdx, numX, numY)) + (1.0 - alphaU) * uPicardOld[i, j];
      });
      
      // ensure ghost cells receive most up to date values
      bc.applyBCs("u", u);

      auto resVPicardStart = 0.0;
      looper.loopWithBoundaries([&](int i, int j) {
        resVPicardStart += std::fabs(v[i, j]);
      });

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

        // compute coefficients for matrix
        auto umax = std::max(uPicardOld[i, j], 0.0) / dx;
        auto umin = std::min(uPicardOld[i, j], 0.0) / dx;
        auto vmax = std::max(vPicardOld[i, j], 0.0) / dy;
        auto vmin = std::min(vPicardOld[i, j], 0.0) / dy;
        auto nudx2 = nu / std::pow(dx, 2);
        auto nudy2 = nu / std::pow(dy, 2);

        auto aP = (1.0 / dt) + umax - umin + vmax - vmin + 2.0 * nudx2 + 2.0 * nudy2;
        auto aE = umin - nudx2;
        auto aW = - umax - nudx2;
        auto aN = vmin - nudy2;
        auto aS = - vmax - nudy2; 

        // Construct coefficient matrix
        Av.coeffRef(ic, ic) = aP;

        // set up right-hand side
        bv(ic) = vOld[i, j] / dt;

        // west
        if (idx == 0) {
          auto westBCType = bcParameters["boundaries"]["west"]["v"][0];
          auto westBCValue = bcParameters["boundaries"]["west"]["v"][1];
          if (westBCType == "dirichlet") {
            Av.coeffRef(ic, ip1) = aE - aW;
            bv(ic) += 2.0 * westBCValue * umax + 2.0 * westBCValue * nudx2;
          } else if (westBCType == "neumann") {
            Av.coeffRef(ic, ip1) = aE + aW;
            bv(ic) += -2.0 * westBCValue * umax * dx - 2.0 * nu * westBCValue / dx;
          }

        // east
        } else if (idx == numX - 1) {
          auto eastBCType = bcParameters["boundaries"]["east"]["v"][0];
          auto eastBCValue = bcParameters["boundaries"]["east"]["v"][1];
          if (eastBCType == "dirichlet") {
            Av.coeffRef(ic, im1) = aW - aE;
            bv(ic) += -2.0 * eastBCValue * umin + 2.0 * eastBCValue * nudx2;
          } else if (eastBCType == "neumann") {
            Av.coeffRef(ic, im1) = aW + aE;
            bv(ic) += - 2.0 * eastBCValue * umin * dx + 2.0 * nu * eastBCValue / dx;
          }
            
        } else {
          Av.coeffRef(ic, ip1) = aE;
          Av.coeffRef(ic, im1) = aW;
        }

        // south
        if (jdx == 0) {
          auto southBCType = bcParameters["boundaries"]["south"]["v"][0];
          auto southBCValue = bcParameters["boundaries"]["south"]["v"][1];
          if (southBCType == "dirichlet") {
            Av.coeffRef(ic, jp1) = aN - aS;
            bv(ic) += 2.0 * southBCValue * vmax + 2.0 * southBCValue * nudy2;
          } else if (southBCType == "neumann") {
            Av.coeffRef(ic, jp1) = aN + aS;
            bv(ic) += -2.0 * southBCValue * vmax * dy - 2.0 * nu * southBCValue / dy;
          }
        
        // north
        } else if (jdx == numY - 1) {
          auto northBCType = bcParameters["boundaries"]["north"]["v"][0];
          auto northBCValue = bcParameters["boundaries"]["north"]["v"][1];
          if (northBCType == "dirichlet") {
            Av.coeffRef(ic, jm1) = aS - aN;
            bv(ic) += -2.0 * northBCValue * vmin + 2.0 * northBCValue * nudy2;
          } else if (northBCType == "neumann") {
            Av.coeffRef(ic, jm1) = aS + aN;
            bv(ic) += - 2.0 * northBCValue * vmin * dy + 2.0 * nu * northBCValue / dy;
          }
          
        } else {
          Av.coeffRef(ic, jp1) = aN;
          Av.coeffRef(ic, jm1) = aS;
        }
      });

      // solve pressure poisson solver with the preconditioned Conjugate Gradient method
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstabV;
      bicgstabV.compute(Av);
      bicgstabV.setMaxIterations(solverParameters["linearSolver"]["v"]["maxIterations"]);
      bicgstabV.setTolerance(solverParameters["linearSolver"]["v"]["tolerance"]);

      xv = bicgstabV.solve(bv);
      vIter = bicgstabV.iterations();

      // update v-velocity with under-relaxation applied
      looper.loopWithBoundaries([&](int i, int j) {
        auto idx = i - numGhostPoints;
        auto jdx = j - numGhostPoints;
        v[i, j] = alphaV * xv(map2Dto1D(idx, jdx, numX, numY)) + (1.0 - alphaV) * vPicardOld[i, j];
      });
      
      // ensure ghost cells receive most up to date values
      bc.applyBCs("v", v);

      auto resPPicardStart = 0.0;
      looper.loopWithBoundaries([&](int i, int j) {
        resPPicardStart += std::fabs(p[i, j]);
      });

      // set up right-hand side and coefficient matrix
      Ap.setZero();
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
        auto dudx = (u[i + 1, j] - u[i - 1, j]) / (2.0 * dx);
        auto dvdy = (v[i, j + 1] - v[i, j - 1]) / (2.0 * dy);
        auto div = (dudx + dvdy);

        // set up matrix coefficients
        auto dx2 = 1.0 / std::pow(dx, 2);
        auto dy2 = 1.0 / std::pow(dy, 2);

        auto aP = - 2.0 * dx2 - 2.0 * dy2;
        auto aE = dx2;
        auto aW = dx2;
        auto aN = dy2;
        auto aS = dy2; 

        // Construct coefficient matrix
        Ap.coeffRef(ic, ic) = aP;

        // set up right-hand side
        bp(ic) = div / dt;

        // west
        if (idx == 0) {
          auto westBCType = bcParameters["boundaries"]["west"]["p"][0];
          auto westBCValue = bcParameters["boundaries"]["west"]["p"][1];
          if (westBCType == "dirichlet") {
            Ap.coeffRef(ic, ip1) = aE - aW;
            bp(ic) -= 2.0 * westBCValue * dx2;
          } else if (westBCType == "neumann") {
            Ap.coeffRef(ic, ip1) = aE + aW;
            bp(ic) += 2.0 * westBCValue / dx;
          }

        // east
        } else if (idx == numX - 1) {
          auto eastBCType = bcParameters["boundaries"]["east"]["p"][0];
          auto eastBCValue = bcParameters["boundaries"]["east"]["p"][1];
          if (eastBCType == "dirichlet") {
            Ap.coeffRef(ic, im1) = aW - aE;
            bp(ic) -= 2.0 * eastBCValue * dx2;
          } else if (eastBCType == "neumann") {
            Ap.coeffRef(ic, im1) = aW + aE;
            bp(ic) -= 2.0 * eastBCValue / dx;
          }
            
        } else {
          Ap.coeffRef(ic, ip1) = aE;
          Ap.coeffRef(ic, im1) = aW;
        }

        // south
        if (jdx == 0) {
          auto southBCType = bcParameters["boundaries"]["south"]["p"][0];
          auto southBCValue = bcParameters["boundaries"]["south"]["p"][1];
          if (southBCType == "dirichlet") {
            Ap.coeffRef(ic, jp1) = aN - aS;
            bp(ic) -= 2.0 * southBCValue * dy2;
          } else if (southBCType == "neumann") {
            Ap.coeffRef(ic, jp1) = aN + aS;
            bp(ic) += 2.0 * southBCValue / dy;
          }
        
        // north
        } else if (jdx == numY - 1) {
          auto northBCType = bcParameters["boundaries"]["north"]["p"][0];
          auto northBCValue = bcParameters["boundaries"]["north"]["p"][1];
          if (northBCType == "dirichlet") {
            Ap.coeffRef(ic, jm1) = aS - aN;
            bp(ic) -= 2.0 * northBCValue * dy2;
          } else if (northBCType == "neumann") {
            Ap.coeffRef(ic, jm1) = aS + aN;
            bp(ic) -= 2.0 * northBCValue / dy;
          }
          
        } else {
          Ap.coeffRef(ic, jp1) = aN;
          Ap.coeffRef(ic, jm1) = aS;
        }

        // check if pressure has fully neumann boundary condition, if so, compute average and subtract
        if (idx == 0 && jdx == 0) {
          if (fullyNeumann) {
            Ap.coeffRef(ic, ic) = 1.0;
            bp(ic) = 0.0;            
          }
        }
      });

      // solve pressure poisson solver with the preconditioned Conjugate Gradient method
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> bicgstabP;
      bicgstabP.compute(Ap);
      bicgstabP.setMaxIterations(solverParameters["linearSolver"]["p"]["maxIterations"]);
      bicgstabP.setTolerance(solverParameters["linearSolver"]["p"]["tolerance"]);

      xp = bicgstabP.solve(bp);
      pIter = bicgstabP.iterations();

      // update pressure with under-relaxation applied
      looper.loopWithBoundaries([&](int i, int j) {
        auto idx = i - numGhostPoints;
        auto jdx = j - numGhostPoints;
        p[i, j] = alphaP * xp(map2Dto1D(idx, jdx, numX, numY)) + (1.0 - alphaP) * pPicardOld[i, j];
      });

      if (fullyNeumann) p[numGhostPoints, numGhostPoints] = 0.0;
      
      // ensure ghost cells receive most up to date values
      bc.applyBCs("p", p);
      
      // update velocity
      looper.loopWithBoundaries([&](int i, int j) {
        u[i, j] = u[i, j] - dt * (p[i + 1, j] - p[i - 1, j]) / (2.0 * dx);
        v[i, j] = v[i, j] - dt * (p[i, j + 1] - p[i, j - 1]) / (2.0 * dy);
      });
      
      bc.applyBCs("u", u);
      bc.applyBCs("v", v);
      
      // update residual
      auto [uPicardResValue, vPicardResValue, pPicardResValue] = picardResiduals.getResidual(k);
      
      // output picard iteration statistics
      ui.updatePicardIteration(k + 1, uPicardResValue, vPicardResValue, pPicardResValue, uIter, vIter, pIter);

      // break out of picard iterations if linearisation loop has converged
      if (uPicardResValue < picardToleranceU && vPicardResValue < picardToleranceV) break;
    } // end outer picard loop

    auto [uResidualValue, vResidualValue, pResidualValue] = outerResiduals.getResidual(t);

    if (outputFrequency != -1 && (t + 1) % outputFrequency == 0) 
      output.write(t + 1);

    // determine total time at the end of time step
    totalTime = totalTime + dt;
    
    // check divergence of velocity field
    double divU = 0.0;
    looper.loopInterior([&](int i, int j) { 
      auto dudx = (u[i + 1, j] - u[i - 1, j]) / (2.0 * dx);
      auto dvdy = (v[i, j + 1] - v[i, j - 1]) / (2.0 * dy);
      divU = dudx + dvdy;
    });

    // update residuals information
    ui.updateIteration(uResidualValue, vResidualValue, pResidualValue, divU);
    
    // check convergence
    if (uResidualValue < epsU && vResidualValue < epsV && pResidualValue < epsP) {
      ui.draw(17, 0, "Solution converged. Press any key to continue!");
      break;
    }
  }

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
  auto [uResidualValue, vResidualValue, pResidualValue] = picardResiduals.getResidual(t);
  if (uResidualValue > epsU || vResidualValue > epsV || pResidualValue > epsP)
    ui.draw(17, 0, "Simulation finished but did not converge. Press any key to continue!");
  ui.block();

  return 0;
}