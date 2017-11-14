#pragma once

#include "vec3f.hpp"
#include "triangle.hpp"

#include <math.h>
#include <limits>
#include <iostream>

namespace lw {


class BSPTree {

public:

    BSPTree();

    inline int64_t numNodes() const { return num_nodes; }
    inline int64_t numLeaves() const { return num_leaves; }

    void build(std::vector<trianglePtr>& _tris);
    float raytrace(vec3f& orig, vec3f& dir, int64_t plane, int64_t address, trianglePtr& hittri );

private:

    std::vector<float> splits;
    std::vector<std::vector<trianglePtr> > tris;  // half this space is wasted

    std::vector<float> workspace;  // used for spliting triangles into left and right

    int64_t num_nodes, num_leaves;

    void pbuild(int64_t plane, int64_t depth, int64_t address, std::vector<trianglePtr>& _tris);

};


}// namespace lw
