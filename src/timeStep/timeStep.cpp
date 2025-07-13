#include "src/timeStep/timeStep.hpp"

TimeStep::TimeStep(nlohmann::json parameters, MeshLooper looper) : _looper(looper), _CFL(parameters["time"]["CFL"]),
  _maxCFL(parameters["time"]["maxCFL"]), _useSER(parameters["time"]["SER"]), _alphaSER(parameters["time"]["alphaSER"])
  { }

TimeStep::TimeStepType TimeStep::getTimeStep(std::shared_ptr<FieldType> u, std::shared_ptr<FieldType> v, std::shared_ptr<FieldType> p,
  double dx, double dy, double resU, double resV, double resP) {
  double CFL = _CFL;
  if (_useSER) {
    // https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=87f810d91a9094924942a60f072b6ffb5e94918d
    double res = std::max(resU, std::max(resV, resP));
    CFL = std::min(_CFL * std::pow((1.0 / res), _alphaSER), _maxCFL); 
  }

  double dt = std::numeric_limits<double>::max();
  _looper.loopWithBoundaries([CFL, dx, dy, &dt, u, v](int i, int j) {
    auto minSpacing = std::min(dx, dy);
    auto velocityMag = std::sqrt((*u)[i][j] * (*u)[i][j] + (*v)[i][j] * (*v)[i][j]);
    auto inviscidTimeStep = CFL * minSpacing / velocityMag;
    if (inviscidTimeStep < dt) dt = inviscidTimeStep;
  });
  return TimeStep::TimeStepType(dt, CFL);
}