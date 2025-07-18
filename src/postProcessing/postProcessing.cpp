#include "src/postProcessing/postProcessing.hpp"

PostProcessing::PostProcessing(std::string filename, int numX, int numY, int numGhostPoints, FieldType &x, FieldType &y)
  : _filename(filename), _numX(numX), _numY(numY), _numGhostPoints(numGhostPoints), _x(x), _y(y) {
  // delete the output folder and then recreate it
  std::filesystem::remove_all("output/");
  std::filesystem::create_directories("output/");
}
  
void PostProcessing::registerField(std::string name, FieldArray *field) {
  _fields[name] = field;
}

void PostProcessing::write(int iteration) {
  
  // construct filename
  std::stringstream file;
  file << "output/";
  if (iteration != -1) {
    file << _filename << "_" << std::setw(10) << std::setfill('0') << iteration << ".dat";
  } else {
    file << _filename << ".dat";
  }

  std::ofstream solution(file.str());
  solution << "TITLE = \"Solution\"\n";
  solution << "VARIABLES = \"x\", \"y\"";
  for (const auto &fields : _fields)
    solution << ", \"" << fields.first << "\"";
  solution << "\n";

  solution << "ZONE T=\"Solution\", I=" << _numX << ", J=" << _numY << ", F=POINT\n";

  // write coordinates and data
  solution << std::scientific << std::setprecision(6);
  for (int j = _numGhostPoints; j < _numY + _numGhostPoints; ++j)
    for (int i = _numGhostPoints; i < _numX + _numGhostPoints; ++i) {
      solution << _x[i][j] << " " << _y[i][j];
      for (const auto &fields : _fields)
        solution << " " << (*fields.second)[i, j];
      solution << "\n";
    }
}