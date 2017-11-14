#pragma once

#include <vector>
#include <boost/shared_ptr.hpp>
#include "triangle.hpp"

namespace lw {

class mesh {
private:
    std::vector<trianglePtr> tris;
public:
    mesh();
    void addTriangle(trianglePtr& tri);
    inline std::vector<trianglePtr>& getTris(){return tris;}
    inline int triCount(){ return tris.size(); }

    friend std::ostream& operator << (std::ostream& stream, mesh& v)
    {
        stream << "mesh {tris: " << v.tris.size() << "}";
        return stream;
    }

};

typedef boost::shared_ptr<mesh> meshPtr;

} // namespace lw
