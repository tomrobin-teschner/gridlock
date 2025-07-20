#include "src/postProcessing/postProcessing.hpp"

PostProcessing::PostProcessing(std::string filename, const Mesh& mesh)
  : _filename(filename), _mesh(mesh) {
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

  solution << "ZONE T=\"Solution\", I=" << _mesh.numX() << ", J=" << _mesh.numY() << ", F=POINT\n";

  // write coordinates and data
  _mesh.loop().loopWithBoundariesReversed([this, &solution](int i, int j) {
    solution << std::scientific << std::setprecision(6);
    solution << _mesh.x(i, j) << " " << _mesh.y(i, j);
    for (const auto &fields : _fields)
      solution << " " << (*fields.second)[i, j];
    solution << "\n";
  });
}