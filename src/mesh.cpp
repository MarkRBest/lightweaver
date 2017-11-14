#include "mesh.hpp"

lw::mesh::mesh(){}

void lw::mesh::addTriangle(trianglePtr& tri) {
    tris.push_back(tri);
}
