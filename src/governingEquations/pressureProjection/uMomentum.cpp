#include "src/governingEquations/pressureProjection/uMomentum.hpp"

uMomentum::uMomentum(nlohmann::json parameters, MeshLooper &meshLooper)
  : GoverningEquationsBase(parameters, meshLooper) {}


void uMomentum::solve() {
  // get mesh parameters
  auto numGhostpoints = _meshLooper.getNumGhostPoints();
  auto numX = _meshLooper.getNumX();
  auto numY = _meshLooper.getNumY();

  Eigen::VectorXd bu((numX) * (numY));
  Eigen::VectorXd xu((numX) * (numY));
  Eigen::SparseMatrix<double> Au((numX) * (numY), (numX) * (numY));
  double alphaU = _parameters["linearSolver"]["u"]["underRelaxation"];

  Au.setZero();
  _meshLooper.loopWithBoundaries([&](int i, int j) {
    // create zero-based inmdices for i and j direction
    auto idx = i - numGhostpoints;
    auto jdx = j - numGhostpoints;

    // map 2D indices to 1D indices. Subtract 2 from numX|Y to ignore boundary points
    auto ic  = _meshLooper.map2Dto1D(idx, jdx);
    auto ip1 = _meshLooper.map2Dto1D(idx + 1, jdx);
    auto im1 = _meshLooper.map2Dto1D(idx - 1, jdx);
    auto jp1 = _meshLooper.map2Dto1D(idx, jdx + 1);
    auto jm1 = _meshLooper.map2Dto1D(idx, jdx - 1);
    
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
  bicgstabU.setMaxIterations(_parameters["linearSolver"]["u"]["maxIterations"]);
  bicgstabU.setTolerance(_parameters["linearSolver"]["u"]["tolerance"]);

  xu = bicgstabU.solve(bu);
  uIter = bicgstabU.iterations();
  
  // update u-velocity with under-relaxation applied
  _meshLooper.loopWithBoundaries([&](int i, int j) {
    auto idx = i - numGhostPoints;
    auto jdx = j - numGhostPoints;
    u[i, j] = alphaU * xu(_meshLooper.map2Dto1D(idx, jdx)) + (1.0 - alphaU) * uPicardOld[i, j];
  });
  
  // ensure ghost cells receive most up to date values
  bc.updateGhostPoints("u", u);
}