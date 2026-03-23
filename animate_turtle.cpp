#include <gmsh.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <cmath>
#include <vector>
#include <string>

int main() {
    gmsh::initialize();
    gmsh::model::add("turtle_anim");

    gmsh::merge("turtle.stl");
    gmsh::model::mesh::classifySurfaces(80 * M_PI / 180, true, true, M_PI);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);

    std::vector<int> surfaceTags;
    for (auto s : surfaces) surfaceTags.push_back(s.second);

    gmsh::model::geo::synchronize();
    int loop = gmsh::model::geo::addSurfaceLoop(surfaceTags);
    gmsh::model::geo::addVolume({loop});
    gmsh::model::geo::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMin", 1.4);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 2.4);
    gmsh::model::mesh::generate(3);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, tmp;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, tmp);

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 3);

    int tetraIndex = -1;
    for (int i = 0; i < elementTypes.size(); i++)
        if (elementTypes[i] == 4) { tetraIndex = i; break; }

    int pointCount = nodeTags.size();

    std::size_t maxTag = 0;
    for (int i = 0; i < pointCount; i++)
        if (nodeTags[i] > maxTag) maxTag = nodeTags[i];
    std::vector<int> tagToIndex(maxTag + 1, -1);
    for (int i = 0; i < pointCount; i++) tagToIndex[nodeTags[i]] = i;

    double minX = 1e9, maxX = -1e9;
    double minY = 1e9, maxY = -1e9;
    double minZ = 1e9, maxZ = -1e9;
    for (int i = 0; i < pointCount; i++) {
        double x = nodeCoords[3 * i];
        double y = nodeCoords[3 * i + 1];
        double z = nodeCoords[3 * i + 2];
        if (x < minX) minX = x; if (x > maxX) maxX = x;
        if (y < minY) minY = y; if (y > maxY) maxY = y;
        if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
    }
    double centerX = (minX + maxX) / 2;
    double centerY = (minY + maxY) / 2;
    double centerZ = (minZ + maxZ) / 2;
    double sizeX = maxX - minX;
    double sizeY = maxY - minY;
    double sizeZ = maxZ - minZ;

    int frames = 240;
    double totalTime = 3.0;
    double omega = 2 * M_PI * 1.5;

    system("mkdir -p output");

    for (int f = 0; f < frames; f++) {
        double t = totalTime * f / (frames - 1);

        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(pointCount);

        auto velocity = vtkSmartPointer<vtkDoubleArray>::New();
        velocity->SetName("Velocity");
        velocity->SetNumberOfComponents(3);
        velocity->SetNumberOfTuples(pointCount);

        auto wave = vtkSmartPointer<vtkDoubleArray>::New();
        wave->SetName("Wave");
        wave->SetNumberOfTuples(pointCount);

        auto color = vtkSmartPointer<vtkUnsignedCharArray>::New();
        color->SetName("Color");
        color->SetNumberOfComponents(3);
        color->SetNumberOfTuples(pointCount);

        for (int i = 0; i < pointCount; i++) {
            double x = nodeCoords[3 * i];
            double y = nodeCoords[3 * i + 1];
            double z = nodeCoords[3 * i + 2];

            double normX = (x - centerX) / (sizeX / 2);
            double normY = (y - centerY) / (sizeY / 2);
            double normZ = (z - centerZ) / (sizeZ / 2);

            double legFactor = 0.0;
            if (normX > 0.3 && normY < -0.3) legFactor = 1.0;
            if (normX < -0.3 && normY < -0.3) legFactor = 1.0;

            double amplitude = 0.15 * sizeY * legFactor;
            double phase = omega * t;
            if (normX > 0) phase += M_PI;

            double shiftY = amplitude * sin(phase);
            double shiftX = 0.05 * amplitude * cos(phase);

            double vx = 0.05 * amplitude * (-omega) * sin(phase);
            double vy = amplitude * omega * cos(phase);
            double vz = 0.0;

            if (legFactor == 0) {
                shiftX = shiftY = 0;
                vx = vy = vz = 0;
            }

            points->SetPoint(i, x + shiftX, y + shiftY, z);
            velocity->SetTuple3(i, vx, vy, vz);

            double radius = sqrt((x - centerX)*(x - centerX) + (y - centerY)*(y - centerY) + (z - centerZ)*(z - centerZ));
            double scalar = sin(6 * radius / (sizeX + sizeY) - 2 * omega * t);
            wave->SetValue(i, scalar);

            double blend = 0.5 + 0.5 * scalar;
            if (blend < 0) blend = 0;
            if (blend > 1) blend = 1;
            unsigned char r = (1 - blend) * 60 + blend * 150;
            unsigned char g = (1 - blend) * 120 + blend * 100;
            unsigned char b = (1 - blend) * 40 + blend * 50;
            color->SetTuple3(i, r, g, b);
        }

        auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        grid->SetPoints(points);

        std::vector<std::size_t> &tetraNodes = elementNodeTags[tetraIndex];
        int tetraCount = tetraNodes.size() / 4;
        for (int i = 0; i < tetraCount; i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            for (int j = 0; j < 4; j++) {
                std::size_t gmshTag = tetraNodes[4 * i + j];
                tetra->GetPointIds()->SetId(j, tagToIndex[gmshTag]);
            }
            grid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        grid->GetPointData()->AddArray(velocity);
        grid->GetPointData()->SetActiveVectors("Velocity");
        grid->GetPointData()->AddArray(wave);
        grid->GetPointData()->SetActiveScalars("Wave");
        grid->GetPointData()->AddArray(color);

        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(("output/frame_" + std::to_string(f) + ".vtu").c_str());
        writer->SetInputData(grid);
        writer->SetDataModeToAscii();
        writer->SetCompressor(nullptr);
        writer->Write();
    }

    gmsh::finalize();
    return 0;
}
