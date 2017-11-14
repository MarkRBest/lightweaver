#pragma once

#include <vector>
#include "camera.hpp"
#include "mesh.hpp"
#include "light.hpp"

namespace lw{

class sceen{
private:
    cameraPtr c;
    std::vector<lightPtr> l;
    std::vector<meshPtr> m;
public:
    sceen(){}
    void addCamera(cameraPtr _c);
    cameraPtr getCam();

    void addLight(lightPtr _l);
    std::vector<lightPtr>& getLights(){return l;}

    void addMesh(meshPtr _m);
    std::vector<meshPtr>& getMeshs(){return m;}

    friend std::ostream& operator << (std::ostream& stream, sceen& v)
    {
        stream << "sceen" << std::endl;
        stream << *v.c << std::endl;

        for(std::vector<lightPtr>::iterator it = v.l.begin(); it != v.l.end(); ++it) {
            stream << *(*it) << std::endl;
        }

        for(std::vector<meshPtr>::iterator it = v.m.begin(); it != v.m.end(); ++it) {
            stream << *(*it) << std::endl;
        }

        return stream;
    }
};

} // namespace lw