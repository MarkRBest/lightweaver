#pragma once

#include <math.h>
#include <iostream>
#include <boost/shared_ptr.hpp>

namespace lw {

class vec2f {

private:
    float data[2];

public:
    vec2f(){}
    vec2f(float x, float y) {
        data[0] = x;
        data[1] = y;
    }

    inline float& operator [](int i) {
        return data[i];
    }


    friend std::ostream& operator <<(std::ostream& stream, const vec2f& v) {
        stream << "vec2f {x: " << v.data[0];
        stream << " y: " << v.data[1];
        stream << "}";
        return stream;
    }
};

typedef boost::shared_ptr<vec2f> vec2fPtr;

} // namespace lw

