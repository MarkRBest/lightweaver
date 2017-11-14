#pragma once

#include <boost/shared_ptr.hpp>
#include "vec3f.hpp"

namespace lw{

class camera{
private:
    float m_ax, m_ay, m_az;
    vec3f point;
public:
    camera(float ax, float ay, float az, float x, float y, float z);
    float ax(){return m_ax;}
    float ay(){return m_ay;}
    float az(){return m_az;}
    vec3f&  pnt(){return point;}
    friend std::ostream& operator << (std::ostream& stream, camera& v)
    {
        stream << "camera {angle_x: " << v.m_ax;
        stream << " angle_y: " << v.m_ay;
        stream << " angle_z: " << v.m_az;
        stream << " x: " << v.point[0];
        stream << " y: " << v.point[1];
        stream << " z: " << v.point[2] <<"}";
        return stream;
    }
};

typedef boost::shared_ptr<camera> cameraPtr;

} // namespace lw
