#pragma once

#include <boost/shared_ptr.hpp>
#include "vec3f.hpp"

namespace lw {

class light{

private:
    vec3f point;

public:
    light(float x, float y, float z);
    vec3f& pnt(){return point;}

    friend std::ostream& operator << (std::ostream& stream, light& v)
    {
        stream << "light {x: " << v.point[0];
        stream << " y: " << v.point[1];
        stream << " z: " << v.point[2] <<"}";
        return stream;
    }

};

typedef boost::shared_ptr<light> lightPtr;

} // namespace lw
