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
#include "src/boundaryConditions/boundaryConditions.hpp"
#include "src/timeStep/timeStep.hpp"
#include "src/infrastructure/parameterFile/parameterFile.hpp"
#include "src/infrastructure/ui/ui.hpp"
#include "src/infrastructure/stopWatch/stopWatch.hpp"
#include "src/infrastructure/utilities/data.hpp"
#include "src/postProcessing/postProcessing.hpp"
#include "src/fieldArray/fieldArray.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"
#include "src/residuals/residuals.hpp"
#include "src/linearSolver/linearSolverEigen.hpp"

#include "src/governingEquations/pressureProjection/momentumU.hpp"

#include "Eigen/Eigen"
#include "nlohmann/json.hpp"
#include "toml++/toml.hpp"

int main(int argv, char* argc[]) {

  // Read parameter files
  ParameterFile parameterFile("input/solver.toml");
  auto solverParameters = parameterFile.getSolverParameters();
  auto meshParameters = parameterFile.getMeshParameters();
  auto bcParameters = parameterFile.getBoundaryConditionsParameters();

  // check if problem is a fully neumann boundary value problem for the pressure
  std::string pWestBC = bcParameters["boundaries"]["west"]["p"][0].value_or("");
  std::string pEastBC = bcParameters["boundaries"]["east"]["p"][0].value_or("");
  std::string pSouthBC = bcParameters["boundaries"]["south"]["p"][0].value_or("");
  std::string pNorthBC = bcParameters["boundaries"]["north"]["p"][0].value_or("");

  bool fullyNeumann = false;
  if (pWestBC == "neumann" && pEastBC == "neumann" && pSouthBC == "neumann" && pNorthBC == "neumann")
    fullyNeumann = true;

    std::cout << solverParameters["config_files"]["meshFile"] << std::endl;
  std::cout << pWestBC << " " << pEastBC << " " << pSouthBC << " " << pNorthBC << " " << fullyNeumann << std::endl;

  int numGhostPoints = 1;
  int totalSizeX = meshParameters["mesh"]["numX"].value_or(0) + 2 * numGhostPoints;
  int totalSizeY = meshParameters["mesh"]["numY"].value_or(0) + 2 * numGhostPoints;
  
  // Time stepping parameters
  int timeSteps = solverParameters["time"]["timeSteps"].value_or(0);
  double CFL = solverParameters["time"]["CFL"].value_or(0.0);
  int outputFrequency = solverParameters["output"]["outputFrequency"].value_or(0);

  // residual and convergence checking
  int uIter = 0; int vIter = 0; int pIter = 0;
  double epsU = solverParameters["convergence"]["u"].value_or(1e-3);
  double epsV = solverParameters["convergence"]["v"].value_or(1e-3);
  double epsP = solverParameters["convergence"]["p"].value_or(1e-3);
  
  // set up picard linearisation
  int maxPicardIterations = solverParameters["linearisation"]["maxPicardIterations"].value_or(0);
  double picardToleranceU = solverParameters["linearisation"]["tolerance"]["u"].value_or(1e-3);
  double picardToleranceV = solverParameters["linearisation"]["tolerance"]["v"].value_or(1e-3);
  double picardToleranceP = solverParameters["linearisation"]["tolerance"]["p"].value_or(1e-3);

  // initialise UI
  UI ui(timeSteps, maxPicardIterations);
  ui.createSkeleton();

  // Mesh generation
  Mesh mesh(meshParameters, numGhostPoints);

  // input for the linear solver
  LinearSolverEigen<Eigen::SparseMatrix<double>, Eigen::VectorXd,
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>> uSolver(mesh.numX(), mesh.numY());
  
  LinearSolverEigen<Eigen::SparseMatrix<double>, Eigen::VectorXd,
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>> vSolver(mesh.numX(), mesh.numY());
  
  LinearSolverEigen<Eigen::SparseMatrix<double>, Eigen::VectorXd,
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>> pSolver(mesh.numX(), mesh.numY());
  
  double alphaU = solverParameters["linearSolver"]["u"]["underRelaxation"].value_or(1.0);
  double alphaV = solverParameters["linearSolver"]["v"]["underRelaxation"].value_or(1.0);
  double alphaP = solverParameters["linearSolver"]["p"]["underRelaxation"].value_or(1.0);
  
  // physical properties
  double nu = solverParameters["fluid"]["nu"].value_or(1.0);
  
  // solution vector
  FieldArrayManager pv(totalSizeX, totalSizeY);

  // Set up the governing equations
  MomentumU momentumU(pv, solverParameters, mesh);

  // create post processing object
  std::string outputFileName = solverParameters["output"]["filename"].value_or("");
  PostProcessing output(pv, outputFileName, mesh);
  output.registerFields({PV::U, PV::V, PV::P});

  Residuals outerResiduals(pv, mesh, {PV::U, PV::V, PV::P}, true);
  Residuals picardResiduals(pv, mesh, {PV::U, PV::V, PV::P});
  
  // Instantiate boundary condition class
  BoundaryConditions bc(mesh, pv, bcParameters);
  bc.updateGhostPoints({PV::U, PV::V, PV::P});

  // create high-resolution stop watch
  StopWatch timer(timeSteps);

  // Create time step calculation object
  TimeStep timeStep(solverParameters, mesh, pv);

  // loop over time
  double totalTime = 0.0;
  int t = 0;
  for (t = 0; t<timeSteps; ++t) {
    
    // create deep copy of old solution
    pv.storeOldFields();
    
    // initialise residuals
    outerResiduals.init();
    
    // determine stable timestep
    auto dt = timeStep.getTimeStep();
    
    // get timings
    auto [elapsedHH, elapsedMM, elapsedSS] = timer.elapsed();
    auto [remainingHH, remainingMM, remainingSS] = timer.remaining(t);

    // update time information of what is already available
    ui.updateTimestatistics(t + 1, CFL, dt, totalTime, elapsedHH, elapsedMM, elapsedSS,
      remainingHH, remainingMM, remainingSS);    

    // outer (picard) loop
    for (int k = 0; k < maxPicardIterations; ++k) {

      // create deep copy of old solution
      pv.storePicardOldFields();

      picardResiduals.init();
    
      // solve the u-momentum equations
      uSolver.setZero();
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        // create indices
        auto [idx, jdx] = mesh.loop().zeroBasedIndices(i, j);
        auto [ic, ip1, im1, jp1, jm1] = mesh.loop().getMatrixIndices(idx, jdx);

        // initialise matrix coefficients
        double aP = 0.0, aE = 0.0, aW = 0.0, aN = 0.0, aS = 0.0, b = 0.0;

        // contributions due to du/dt
        aP += (1.0 / dt);
        b  += pv(PV::U_OLD)[i, j] / dt;

        // contributions due to u*div(u)
        auto umax = std::max(pv(PV::U_PICARD_OLD)[i, j], 0.0) / mesh.dx();
        auto umin = std::min(pv(PV::U_PICARD_OLD)[i, j], 0.0) / mesh.dx();
        auto vmax = std::max(pv(PV::V_PICARD_OLD)[i, j], 0.0) / mesh.dy();
        auto vmin = std::min(pv(PV::V_PICARD_OLD)[i, j], 0.0) / mesh.dy();

        aP += umax - umin + vmax - vmin;
        aE += umin; aW += - umax; aN += vmin; aS += - vmax;

        // contributions due to nu*div(grad(u))
        auto nudx2 = nu / std::pow(mesh.dx(), 2);
        auto nudy2 = nu / std::pow(mesh.dy(), 2);
        
        aP += 2.0 * nudx2 + 2.0 * nudy2;
        aE -= nudx2; aW -= nudx2; aN -= nudy2; aS -= nudy2; 

        // Construct coefficient matrix
        uSolver.setMatrixAt(ic, ic, aP);

        // set up right-hand side
        uSolver.setRHSAt(ic, b);

        // west
        if (idx == 0) {
          std::string westBCType = bcParameters["boundaries"]["west"]["u"][0].value_or("");
          double westBCValue = bcParameters["boundaries"]["west"]["u"][1].value_or(0.0);
          if (westBCType == "dirichlet") {
            uSolver.setMatrixAt(ic, ip1, aE - aW);
            uSolver.addRHSAt(ic, 2.0 * westBCValue * umax + 2.0 * westBCValue * nudx2);
          } else if (westBCType == "neumann") {
            uSolver.setMatrixAt(ic, ip1, aE + aW);
            uSolver.addRHSAt(ic, -2.0 * westBCValue * umax * mesh.dx() - 2.0 * nu * westBCValue / mesh.dx());
          }

        // east
        } else if (idx == mesh.numX() - 1) {
          std::string eastBCType = bcParameters["boundaries"]["east"]["u"][0].value_or("");
          double eastBCValue = bcParameters["boundaries"]["east"]["u"][1].value_or(0.0);
          if (eastBCType == "dirichlet") {
            uSolver.setMatrixAt(ic, im1, aW - aE);
            uSolver.addRHSAt(ic, -2.0 * eastBCValue * umin + 2.0 * eastBCValue * nudx2);
          } else if (eastBCType == "neumann") {
            uSolver.setMatrixAt(ic, im1, aW + aE);
            uSolver.addRHSAt(ic, - 2.0 * eastBCValue * umin * mesh.dx() + 2.0 * nu * eastBCValue / mesh.dx());
          }
            
        } else {
          uSolver.setMatrixAt(ic, ip1, aE);
          uSolver.setMatrixAt(ic, im1, aW);
        }

        // south
        if (jdx == 0) {
          std::string southBCType = bcParameters["boundaries"]["south"]["u"][0].value_or("");
          double southBCValue = bcParameters["boundaries"]["south"]["u"][1].value_or(0.0);
          if (southBCType == "dirichlet") {
            uSolver.setMatrixAt(ic, jp1, aN - aS);
            uSolver.addRHSAt(ic, 2.0 * southBCValue * vmax + 2.0 * southBCValue * nudy2);
          } else if (southBCType == "neumann") {
            uSolver.setMatrixAt(ic, jp1, aN + aS);
            uSolver.addRHSAt(ic, -2.0 * southBCValue * vmax * mesh.dy() - 2.0 * nu * southBCValue / mesh.dy());
          }
        
        // north
        } else if (jdx == mesh.numY() - 1) {
          std::string northBCType = bcParameters["boundaries"]["north"]["u"][0].value_or("");
          double northBCValue = bcParameters["boundaries"]["north"]["u"][1].value_or(0.0);
          if (northBCType == "dirichlet") {
            uSolver.setMatrixAt(ic, jm1, aS - aN);
            uSolver.addRHSAt(ic, -2.0 * northBCValue * vmin + 2.0 * northBCValue * nudy2);
          } else if (northBCType == "neumann") {
            uSolver.setMatrixAt(ic, jm1, aS + aN);
            uSolver.addRHSAt(ic, - 2.0 * northBCValue * vmin * mesh.dy() + 2.0 * nu * northBCValue / mesh.dy());
          }
          
        } else {
          uSolver.setMatrixAt(ic, jp1, aN);
          uSolver.setMatrixAt(ic, jm1, aS);
        }
      });

      // solve pressure poisson solver with the preconditioned Conjugate Gradient method
      auto [uIter, xu] = uSolver.solve(solverParameters["linearSolver"]["u"]["maxIterations"].value_or(0), 
        solverParameters["linearSolver"]["u"]["tolerance"].value_or(1e-3));

      // update u-velocity with under-relaxation applied
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        auto idx = i - numGhostPoints;
        auto jdx = j - numGhostPoints;
        pv(PV::U)[i, j] = alphaU * xu(mesh.loop().map2Dto1D(idx, jdx)) + (1.0 - alphaU) * pv(PV::U_PICARD_OLD)[i, j];
      });
      
      // ensure ghost cells receive most up to date values
      bc.updateGhostPoints(PV::U);

      // solve the v-momentum equations
      vSolver.setZero();
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        // create indices
        auto [idx, jdx] = mesh.loop().zeroBasedIndices(i, j);
        auto [ic, ip1, im1, jp1, jm1] = mesh.loop().getMatrixIndices(idx, jdx);

        // compute coefficients for matrix
        auto umax = std::max(pv(PV::U_PICARD_OLD)[i, j], 0.0) / mesh.dx();
        auto umin = std::min(pv(PV::U_PICARD_OLD)[i, j], 0.0) / mesh.dx();
        auto vmax = std::max(pv(PV::V_PICARD_OLD)[i, j], 0.0) / mesh.dy();
        auto vmin = std::min(pv(PV::V_PICARD_OLD)[i, j], 0.0) / mesh.dy();
        auto nudx2 = nu / std::pow(mesh.dx(), 2);
        auto nudy2 = nu / std::pow(mesh.dy(), 2);

        auto aP = (1.0 / dt) + umax - umin + vmax - vmin + 2.0 * nudx2 + 2.0 * nudy2;
        auto aE = umin - nudx2;
        auto aW = - umax - nudx2;
        auto aN = vmin - nudy2;
        auto aS = - vmax - nudy2; 

        // Construct coefficient matrix
        vSolver.setMatrixAt(ic, ic, aP);

        // set up right-hand side
        vSolver.setRHSAt(ic, pv(PV::V_OLD)[i, j] / dt);

        // west
        if (idx == 0) {
          std::string westBCType = bcParameters["boundaries"]["west"]["v"][0].value_or("");
          double westBCValue = bcParameters["boundaries"]["west"]["v"][1].value_or(0.0);
          if (westBCType == "dirichlet") {
            vSolver.setMatrixAt(ic, ip1, aE - aW);
            vSolver.addRHSAt(ic, 2.0 * westBCValue * umax + 2.0 * westBCValue * nudx2);
          } else if (westBCType == "neumann") {
            vSolver.setMatrixAt(ic, ip1, aE + aW);
            vSolver.addRHSAt(ic, -2.0 * westBCValue * umax * mesh.dx() - 2.0 * nu * westBCValue / mesh.dx());
          }

        // east
        } else if (idx == mesh.numX() - 1) {
          std::string eastBCType = bcParameters["boundaries"]["east"]["v"][0].value_or("");
          double eastBCValue = bcParameters["boundaries"]["east"]["v"][1].value_or(0.0);
          if (eastBCType == "dirichlet") {
            vSolver.setMatrixAt(ic, im1, aW - aE);
            vSolver.addRHSAt(ic, -2.0 * eastBCValue * umin + 2.0 * eastBCValue * nudx2);
          } else if (eastBCType == "neumann") {
            vSolver.setMatrixAt(ic, im1, aW + aE);
            vSolver.addRHSAt(ic, - 2.0 * eastBCValue * umin * mesh.dx() + 2.0 * nu * eastBCValue / mesh.dx());
          }
            
        } else {
          vSolver.setMatrixAt(ic, ip1, aE);
          vSolver.setMatrixAt(ic, im1, aW);
        }

        // south
        if (jdx == 0) {
          std::string southBCType = bcParameters["boundaries"]["south"]["v"][0].value_or("");
          double southBCValue = bcParameters["boundaries"]["south"]["v"][1].value_or(0.0);
          if (southBCType == "dirichlet") {
            vSolver.setMatrixAt(ic, jp1, aN - aS);
            vSolver.addRHSAt(ic, 2.0 * southBCValue * vmax + 2.0 * southBCValue * nudy2);
          } else if (southBCType == "neumann") {
            vSolver.setMatrixAt(ic, jp1, aN + aS);
            vSolver.addRHSAt(ic, -2.0 * southBCValue * vmax * mesh.dy() - 2.0 * nu * southBCValue / mesh.dy());
          }
        
        // north
        } else if (jdx == mesh.numY() - 1) {
          std::string northBCType = bcParameters["boundaries"]["north"]["v"][0].value_or("");
          double northBCValue = bcParameters["boundaries"]["north"]["v"][1].value_or(0.0);
          if (northBCType == "dirichlet") {
            vSolver.setMatrixAt(ic, jm1, aS - aN);
            vSolver.addRHSAt(ic, -2.0 * northBCValue * vmin + 2.0 * northBCValue * nudy2);
          } else if (northBCType == "neumann") {
            vSolver.setMatrixAt(ic, jm1, aS + aN);
            vSolver.addRHSAt(ic, - 2.0 * northBCValue * vmin * mesh.dy() + 2.0 * nu * northBCValue / mesh.dy());
          }
          
        } else {
          vSolver.setMatrixAt(ic, jp1, aN);
          vSolver.setMatrixAt(ic, jm1, aS);
        }
      });

      // solve pressure poisson solver with the preconditioned Conjugate Gradient method
      auto [vIter, xv] = vSolver.solve(solverParameters["linearSolver"]["v"]["maxIterations"].value_or(0),
        solverParameters["linearSolver"]["v"]["tolerance"].value_or(1e-3));

      // update v-velocity with under-relaxation applied
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        auto idx = i - numGhostPoints;
        auto jdx = j - numGhostPoints;
        pv(PV::V)[i, j] = alphaV * xv(mesh.loop().map2Dto1D(idx, jdx)) + (1.0 - alphaV) * pv(PV::V_PICARD_OLD)[i, j];
      });
      
      // ensure ghost cells receive most up to date values
      bc.updateGhostPoints(PV::V);

      // set up right-hand side and coefficient matrix
      pSolver.setZero();
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        // create indices
        auto [idx, jdx] = mesh.loop().zeroBasedIndices(i, j);
        auto [ic, ip1, im1, jp1, jm1] = mesh.loop().getMatrixIndices(idx, jdx);

        // set up right-hand side
        auto dudx = (pv(PV::U)[i + 1, j] - pv(PV::U)[i - 1, j]) / (2.0 * mesh.dx());
        auto dvdy = (pv(PV::V)[i, j + 1] - pv(PV::V)[i, j - 1]) / (2.0 * mesh.dy());
        auto div = (dudx + dvdy);

        // set up matrix coefficients
        auto dx2 = 1.0 / std::pow(mesh.dx(), 2);
        auto dy2 = 1.0 / std::pow(mesh.dy(), 2);

        auto aP = - 2.0 * dx2 - 2.0 * dy2;
        auto aE = dx2;
        auto aW = dx2;
        auto aN = dy2;
        auto aS = dy2; 

        // Construct coefficient matrix
        pSolver.setMatrixAt(ic, ic, aP);

        // set up right-hand side
        pSolver.setRHSAt(ic, div / dt);

        // west
        if (idx == 0) {
          std::string westBCType = bcParameters["boundaries"]["west"]["p"][0].value_or("");
          double westBCValue = bcParameters["boundaries"]["west"]["p"][1].value_or(0.0);
          if (westBCType == "dirichlet") {
            pSolver.setMatrixAt(ic, ip1, aE - aW);
            pSolver.addRHSAt(ic, -2.0 * westBCValue * dx2);
          } else if (westBCType == "neumann") {
            pSolver.setMatrixAt(ic, ip1, aE + aW);
            pSolver.addRHSAt(ic, 2.0 * westBCValue / mesh.dx());
          }

        // east
        } else if (idx == mesh.numX() - 1) {
          std::string eastBCType = bcParameters["boundaries"]["east"]["p"][0].value_or("");
          double eastBCValue = bcParameters["boundaries"]["east"]["p"][1].value_or(0.0);
          if (eastBCType == "dirichlet") {
            pSolver.setMatrixAt(ic, im1, aW - aE);
            pSolver.addRHSAt(ic, -2.0 * eastBCValue * dx2);
          } else if (eastBCType == "neumann") {
            pSolver.setMatrixAt(ic, im1, aW + aE);
            pSolver.addRHSAt(ic, -2.0 * eastBCValue / mesh.dx());
          }
            
        } else {
          pSolver.setMatrixAt(ic, ip1, aE);
          pSolver.setMatrixAt(ic, im1, aW);
        }

        // south
        if (jdx == 0) {
          std::string southBCType = bcParameters["boundaries"]["south"]["p"][0].value_or("");
          double southBCValue = bcParameters["boundaries"]["south"]["p"][1].value_or(0.0);
          if (southBCType == "dirichlet") {
            pSolver.setMatrixAt(ic, jp1, aN - aS);
            pSolver.addRHSAt(ic, -2.0 * southBCValue * dy2);
          } else if (southBCType == "neumann") {
            pSolver.setMatrixAt(ic, jp1, aN + aS);
            pSolver.addRHSAt(ic, 2.0 * southBCValue / mesh.dy());
          }
        
        // north
        } else if (jdx == mesh.numY() - 1) {
          std::string northBCType = bcParameters["boundaries"]["north"]["p"][0].value_or("");
          double northBCValue = bcParameters["boundaries"]["north"]["p"][1].value_or(0.0);
          if (northBCType == "dirichlet") {
            pSolver.setMatrixAt(ic, jm1, aS - aN);
            pSolver.addRHSAt(ic, -2.0 * northBCValue * dy2);
          } else if (northBCType == "neumann") {
            pSolver.setMatrixAt(ic, jm1, aS + aN);
            pSolver.addRHSAt(ic, -2.0 * northBCValue / mesh.dy());
          }
          
        } else {
          pSolver.setMatrixAt(ic, jp1, aN);
          pSolver.setMatrixAt(ic, jm1, aS);
        }

        // check if pressure has fully neumann boundary condition, if so, compute average and subtract
        if (idx == 0 && jdx == 0) {
          if (fullyNeumann) {
            pSolver.setMatrixAt(ic, ic, 1.0);
            pSolver.setRHSAt(ic, 0.0);            
          }
        }
      });

      // solve pressure poisson solver with the preconditioned Conjugate Gradient method
      auto [pIter, xp] = pSolver.solve(solverParameters["linearSolver"]["p"]["maxIterations"].value_or(0),
        solverParameters["linearSolver"]["p"]["tolerance"].value_or(1e-3));

      // update pressure with under-relaxation applied
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        auto idx = i - numGhostPoints;
        auto jdx = j - numGhostPoints;
        pv(PV::P)[i, j] = alphaP * xp(mesh.loop().map2Dto1D(idx, jdx)) + (1.0 - alphaP) * pv(PV::P_PICARD_OLD)[i, j];
      });

      if (fullyNeumann) pv(PV::P)[numGhostPoints, numGhostPoints] = 0.0;
      
      // ensure ghost cells receive most up to date values
      bc.updateGhostPoints(PV::P);
      
      // update velocity
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        pv(PV::U)[i, j] = pv(PV::U)[i, j] - dt * (pv(PV::P)[i + 1, j] - pv(PV::P)[i - 1, j]) / (2.0 * mesh.dx());
        pv(PV::V)[i, j] = pv(PV::V)[i, j] - dt * (pv(PV::P)[i, j + 1] - pv(PV::P)[i, j - 1]) / (2.0 * mesh.dy());
      });
      
      bc.updateGhostPoints({PV::U, PV::V});
      
      // update residual
      auto resPicard = picardResiduals.getResidual(k);
      
      // output picard iteration statistics
      ui.updatePicardIteration(k + 1, resPicard[PV::U], resPicard[PV::V], resPicard[PV::P], uIter, vIter, pIter);

      // break out of picard iterations if linearisation loop has converged
      if (resPicard[PV::U] < picardToleranceU && resPicard[PV::V] < picardToleranceV) break;
    } // end outer picard loop

    auto res = outerResiduals.getResidual(t);

    if (outputFrequency != -1 && (t + 1) % outputFrequency == 0) 
      output.write(t + 1);

    // determine total time at the end of time step
    totalTime = totalTime + dt;
    
    // check divergence of velocity field
    double divU = 0.0;
    mesh.loop().loopInterior([&](int i, int j) { 
      auto dudx = (pv(PV::U)[i + 1, j] - pv(PV::U)[i - 1, j]) / (2.0 * mesh.dx());
      auto dvdy = (pv(PV::V)[i, j + 1] - pv(PV::V)[i, j - 1]) / (2.0 * mesh.dy());
      divU = dudx + dvdy;
    });

    // update residuals information
    ui.updateIteration(res[PV::U], res[PV::V], res[PV::P], divU);
    
    // check convergence
    if (res[PV::U] < epsU && res[PV::V] < epsV && res[PV::P] < epsP) {
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
  auto res = outerResiduals.getResidual(t);
  if (res[PV::U] > epsU || res[PV::V] > epsV || res[PV::P] > epsP)
    ui.draw(17, 0, "Simulation finished but did not converge. Press any key to continue!");
  ui.block();

  return 0;
}