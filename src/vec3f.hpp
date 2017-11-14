#pragma once

#include <math.h>
#include <iostream>
#include <boost/shared_ptr.hpp>

namespace lw {

class vec3f {

private:
    float data[3];

public:
    vec3f(){}
    vec3f(float x, float y, float z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }

    inline float& operator [](int i) {
        return data[i];
    }

    inline const float operator [](int i) const {
        return data[i];
    }

    inline vec3f operator -() const {
        return vec3f(-data[0], -data[1], -data[2]);
    }

    inline float length() {
        float len = data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
        return pow(len, 0.5);
    }

    inline void norm() {
        float len = length();
        for (int i = 0; i < 3; ++i)
            data[i] /= len;
    }

    inline void operator =(const vec3f& v) {
        data[0] = v.data[0];
        data[1] = v.data[1];
        data[2] = v.data[2];
    }

    inline void operator +=(const vec3f& v) {
        data[0] += v.data[0];
        data[1] += v.data[1];
        data[2] += v.data[2];
    }

    inline void operator -=(const vec3f& v) {
        data[0] -= v.data[0];
        data[1] -= v.data[1];
        data[2] -= v.data[2];
    }

    inline void operator *=(const vec3f& v) {
        data[0] *= v.data[0];
        data[1] *= v.data[1];
        data[2] *= v.data[2];
    }

    inline vec3f operator +(const vec3f& v) {
        return vec3f(data[0] + v.data[0], data[1] + v.data[1],
                data[2] + v.data[2]);
    }

    inline vec3f operator -(const vec3f& v) {
        return vec3f(data[0] - v.data[0], data[1] - v.data[1],
                data[2] - v.data[2]);
    }

    inline vec3f operator /(float v) {
        return vec3f(data[0] / v, data[1] / v, data[2] / v);
    }

    inline vec3f operator *(float v) {
        return vec3f(data[0] * v, data[1] * v, data[2] * v);
    }

    inline vec3f operator /(int v) {
        return vec3f(data[0] / v, data[1] / v, data[2] / v);
    }

    inline vec3f operator *(int v) {
        return vec3f(data[0] * v, data[1] * v, data[2] * v);
    }

    inline vec3f operator -(vec3f& v) {
        return vec3f(data[0] - v[0], data[1] - v[1], data[2] - v[2]);
    }

    inline vec3f operator +(vec3f& v) {
        return vec3f(data[0] + v[0], data[1] + v[1], data[2] + v[2]);
    }

    inline void operator *=(const float v) {
        data[0] *= v;
        data[1] *= v;
        data[2] *= v;
    }

    inline void operator /=(const float v) {
        data[0] /= v;
        data[1] /= v;
        data[2] /= v;
    }

    inline void operator +=(const float v) {
        data[0] += v;
        data[1] += v;
        data[2] += v;
    }

    inline void operator -=(const float v) {
        data[0] -= v;
        data[1] -= v;
        data[2] -= v;
    }

    friend bool operator <(const vec3f& t1, const vec3f& t2) {
        for (int i = 0; i < 3; i++) {
            if (t1.data[i] < t2.data[i])
                return true;
            else if (t1.data[i] > t2.data[i]) {
                return false;
            }
        }
        return false;
    }

    friend std::ostream& operator <<(std::ostream& stream, const vec3f& v) {
        stream << "vec3f {x: " << v.data[0];
        stream << " y: " << v.data[1];
        stream << " z: " << v.data[2] << "}";
        return stream;
    }

};

typedef boost::shared_ptr<vec3f> vec3fPtr;

} // namespace lw
