#include <gmsh.h>
#include <cmath>
#include <vector>

int main() {
  gmsh::initialize();
  gmsh::model::add("turtle");

  gmsh::merge("turtle.stl");

  gmsh::model::mesh::classifySurfaces(80 * M_PI / 180, true, true, M_PI);
  gmsh::model::mesh::createGeometry();

  std::vector<std::pair<int, int>> surfaces;
  gmsh::model::getEntities(surfaces, 2);

  std::vector<int> tags;
  for (auto s : surfaces) {
    tags.push_back(s.second);
    gmsh::model::setColor({s}, 30, 200, 0);
  }

  gmsh::model::geo::synchronize();

  int surfaceLoop = gmsh::model::geo::addSurfaceLoop(tags);
  gmsh::model::geo::addVolume({surfaceLoop});

  gmsh::model::geo::synchronize();

  gmsh::option::setNumber("Mesh.MeshSizeMin", 2.0);
  gmsh::option::setNumber("Mesh.MeshSizeMax", 3.0);
  gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 2.0);
  gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 3.0);

  gmsh::model::mesh::generate(3);
  gmsh::write("turtle.msh");
  gmsh::fltk::run();
  gmsh::finalize();
}
