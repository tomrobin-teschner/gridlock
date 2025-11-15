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
#include "src/mesh/mesh.hpp"
#include "src/boundaryConditions/boundaryConditions.hpp"
#include "src/timeInfo/timeInfo.hpp"
#include "src/infrastructure/parameters/parameters.hpp"
#include "src/infrastructure/ui/ui.hpp"
#include "src/infrastructure/stopWatch/stopWatch.hpp"
#include "src/infrastructure/utilities/data.hpp"
#include "src/postProcessing/postProcessing.hpp"
#include "src/fieldArray/fieldArray.hpp"
#include "src/fieldArray/fieldArrayManager.hpp"
#include "src/residuals/residuals.hpp"
#include "src/linearSolver/linearSolverEigen.hpp"

// #include "src/governingEquations/pressureProjection/momentumU.hpp"

#include "Eigen/Eigen"

int main(int argv, char* argc[]) {

  // Read parameter files
  Parameters params;

  // initialise UI
  UI ui(params);
  ui.createSkeleton();

  // Mesh generation
  Mesh mesh(params);

  // input for the linear solver
  LinearSolverEigen<Eigen::SparseMatrix<double>, Eigen::VectorXd,
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>> uSolver(params, mesh, "u");
  
  LinearSolverEigen<Eigen::SparseMatrix<double>, Eigen::VectorXd,
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>> vSolver(params, mesh, "v");
  
  LinearSolverEigen<Eigen::SparseMatrix<double>, Eigen::VectorXd,
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>>> pSolver(params, mesh, "p");
  
  // physical properties
  auto nu = params.solver<double>("fluid", "nu");
  
  // solution vector
  FieldArrayManager pv(mesh);

  // // Set up the governing equations
  // MomentumU momentumU(pv, solverParameters, mesh);

  // create post processing object
  PostProcessing output(pv, params, mesh);
  output.registerFields({PV::U, PV::V, PV::P});

  Residuals outerResiduals(params, pv, mesh, {PV::U, PV::V, PV::P}, "linearSolver", true);
  Residuals picardResiduals(params, pv, mesh, {PV::U, PV::V, PV::P}, "linearisation");
  
  // Instantiate boundary condition class
  BoundaryConditions bc(params, mesh, pv);
  bc.updateGhostPoints({PV::U, PV::V, PV::P});

  // create high-resolution stop watch
  StopWatch timer(params);

  // Create time step calculation object
  TimeInfo timeInfo(params, mesh, pv);

  // loop over time
  double totalTime = 0.0;
  int t = 0;
  for (t = 0; t < timeInfo.getNumTimeSteps(); ++t) {
    
    // create deep copy of old solution
    pv.storeOldFields();
    
    // initialise residuals
    outerResiduals.init();
    
    // determine stable timestep
    auto dt = timeInfo.getTimeStep();
    
    // get timings
    auto [elapsedHH, elapsedMM, elapsedSS] = timer.elapsed();
    auto [remainingHH, remainingMM, remainingSS] = timer.remaining(t);

    // update time information of what is already available
    ui.updateTimestatistics(t + 1, timeInfo.getCFL(), dt, totalTime, elapsedHH, elapsedMM, elapsedSS,
      remainingHH, remainingMM, remainingSS);    

    // outer (picard) loop
    for (int k = 0; k < params.solver<int>("linearisation", "maxPicardIterations"); ++k) {

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
          auto westBCType = params.bcs<std::string>("boundaries", "west", "u", 0);
          auto westBCValue = params.bcs<double>("boundaries", "west", "u", 1);
          if (westBCType == "dirichlet") {
            uSolver.setMatrixAt(ic, ip1, aE - aW);
            uSolver.addRHSAt(ic, 2.0 * westBCValue * umax + 2.0 * westBCValue * nudx2);
          } else if (westBCType == "neumann") {
            uSolver.setMatrixAt(ic, ip1, aE + aW);
            uSolver.addRHSAt(ic, -2.0 * westBCValue * umax * mesh.dx() - 2.0 * nu * westBCValue / mesh.dx());
          }

        // east
        } else if (idx == mesh.numX() - 1) {
          auto eastBCType = params.bcs<std::string>("boundaries", "east", "u", 0);
          auto eastBCValue = params.bcs<double>("boundaries", "east", "u", 1);
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
          auto southBCType = params.bcs<std::string>("boundaries", "south", "u", 0);
          auto southBCValue = params.bcs<double>("boundaries", "south", "u", 1);
          if (southBCType == "dirichlet") {
            uSolver.setMatrixAt(ic, jp1, aN - aS);
            uSolver.addRHSAt(ic, 2.0 * southBCValue * vmax + 2.0 * southBCValue * nudy2);
          } else if (southBCType == "neumann") {
            uSolver.setMatrixAt(ic, jp1, aN + aS);
            uSolver.addRHSAt(ic, -2.0 * southBCValue * vmax * mesh.dy() - 2.0 * nu * southBCValue / mesh.dy());
          }
        
        // north
        } else if (jdx == mesh.numY() - 1) {
          auto northBCType = params.bcs<std::string>("boundaries", "north", "u", 0);
          auto northBCValue = params.bcs<double>("boundaries", "north", "u", 1);
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
      auto [uIter, xu] = uSolver.solve();

      // update u-velocity with under-relaxation applied
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        auto idx = i - mesh.numGhostPoints();
        auto jdx = j - mesh.numGhostPoints();
        auto alphaU = uSolver.getUnderRelaxation();
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
          auto westBCType = params.bcs<std::string>("boundaries", "west", "v", 0);
          auto westBCValue = params.bcs<double>("boundaries", "west", "v", 1);
          if (westBCType == "dirichlet") {
            vSolver.setMatrixAt(ic, ip1, aE - aW);
            vSolver.addRHSAt(ic, 2.0 * westBCValue * umax + 2.0 * westBCValue * nudx2);
          } else if (westBCType == "neumann") {
            vSolver.setMatrixAt(ic, ip1, aE + aW);
            vSolver.addRHSAt(ic, -2.0 * westBCValue * umax * mesh.dx() - 2.0 * nu * westBCValue / mesh.dx());
          }

        // east
        } else if (idx == mesh.numX() - 1) {
          auto eastBCType = params.bcs<std::string>("boundaries", "east", "v", 0);
          auto eastBCValue = params.bcs<double>("boundaries", "east", "v", 1);
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
          auto southBCType = params.bcs<std::string>("boundaries", "south", "v", 0);
          auto southBCValue = params.bcs<double>("boundaries", "south", "v", 1);
          if (southBCType == "dirichlet") {
            vSolver.setMatrixAt(ic, jp1, aN - aS);
            vSolver.addRHSAt(ic, 2.0 * southBCValue * vmax + 2.0 * southBCValue * nudy2);
          } else if (southBCType == "neumann") {
            vSolver.setMatrixAt(ic, jp1, aN + aS);
            vSolver.addRHSAt(ic, -2.0 * southBCValue * vmax * mesh.dy() - 2.0 * nu * southBCValue / mesh.dy());
          }
        
        // north
        } else if (jdx == mesh.numY() - 1) {
          auto northBCType = params.bcs<std::string>("boundaries", "north", "v", 0);
          auto northBCValue = params.bcs<double>("boundaries", "north", "v", 1);
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
      auto [vIter, xv] = vSolver.solve();

      // update v-velocity with under-relaxation applied
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        auto idx = i - mesh.numGhostPoints();
        auto jdx = j - mesh.numGhostPoints();
        auto alphaV = vSolver.getUnderRelaxation();
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
          auto westBCType = params.bcs<std::string>("boundaries", "west", "p", 0);
          auto westBCValue = params.bcs<double>("boundaries", "west", "p", 1);
          if (westBCType == "dirichlet") {
            pSolver.setMatrixAt(ic, ip1, aE - aW);
            pSolver.addRHSAt(ic, -2.0 * westBCValue * dx2);
          } else if (westBCType == "neumann") {
            pSolver.setMatrixAt(ic, ip1, aE + aW);
            pSolver.addRHSAt(ic, 2.0 * westBCValue / mesh.dx());
          }

        // east
        } else if (idx == mesh.numX() - 1) {
          auto eastBCType = params.bcs<std::string>("boundaries", "east", "p", 0);
          auto eastBCValue = params.bcs<double>("boundaries", "east", "p", 1);
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
          auto southBCType = params.bcs<std::string>("boundaries", "south", "p", 0);
          auto southBCValue = params.bcs<double>("boundaries", "south", "p", 1);
          if (southBCType == "dirichlet") {
            pSolver.setMatrixAt(ic, jp1, aN - aS);
            pSolver.addRHSAt(ic, -2.0 * southBCValue * dy2);
          } else if (southBCType == "neumann") {
            pSolver.setMatrixAt(ic, jp1, aN + aS);
            pSolver.addRHSAt(ic, 2.0 * southBCValue / mesh.dy());
          }
        
        // north
        } else if (jdx == mesh.numY() - 1) {
          auto northBCType = params.bcs<std::string>("boundaries", "north", "p", 0);
          auto northBCValue = params.bcs<double>("boundaries", "north", "p", 1);
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
          if (bc.fullyNeumann()) {
            pSolver.setMatrixAt(ic, ic, 1.0);
            pSolver.setRHSAt(ic, 0.0);            
          }
        }
      });

      // solve pressure poisson solver with the preconditioned Conjugate Gradient method
      auto [pIter, xp] = pSolver.solve();

      // update pressure with under-relaxation applied
      mesh.loop().loopWithBoundaries([&](int i, int j) {
        auto idx = i - mesh.numGhostPoints();
        auto jdx = j - mesh.numGhostPoints();
        auto alphaP = pSolver.getUnderRelaxation();
        pv(PV::P)[i, j] = alphaP * xp(mesh.loop().map2Dto1D(idx, jdx)) + (1.0 - alphaP) * pv(PV::P_PICARD_OLD)[i, j];
      });

      if (bc.fullyNeumann()) pv(PV::P)[mesh.numGhostPoints(), mesh.numGhostPoints()] = 0.0;
      
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
      if (resPicard[PV::U] < picardResiduals.getTolerance(PV::U) && resPicard[PV::V] < picardResiduals.getTolerance(PV::V)) break;
    } // end outer picard loop

    auto res = outerResiduals.getResidual(t);

    if (timeInfo.getOutputFrequency() != -1 && (t + 1) % timeInfo.getOutputFrequency() == 0) 
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
    if (res[PV::U] < outerResiduals.getTolerance(PV::U) && res[PV::V] < outerResiduals.getTolerance(PV::V) && res[PV::P] < outerResiduals.getTolerance(PV::P)) {
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
  if (res[PV::U] > outerResiduals.getTolerance(PV::U) || res[PV::V] > outerResiduals.getTolerance(PV::V) ||
    res[PV::P] > outerResiduals.getTolerance(PV::P))
    ui.draw(17, 0, "Simulation finished but did not converge. Press any key to continue!");
  ui.block();

  return 0;
}