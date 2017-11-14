#pragma once

#include "vec3f.hpp"

namespace lw{

class ray {

private:
    vec3f o, d;

public:
    ray(float o0, float o1, float o2, float d0, float d1, float d2);
    ray(float* _o, float* _d);

    void setOrig(float x, float y, float z);
    void setOrig(vec3f& orig);
    void setDir(float x, float y, float z);
    void setDir(vec3f& dir);

    inline vec3f& orig() {return o;}
    inline vec3f& dir() {return d;}

    friend std::ostream& operator << (std::ostream& stream, ray& v)
    {
        stream << "ray {";
        stream << "orig: " << v.o;
        stream << " dir: " << v.d <<"}";
        return stream;
    }
};

} //namespace lw
