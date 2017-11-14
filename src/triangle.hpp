#pragma once

#include <boost/shared_ptr.hpp>
#include "vec3f.hpp"
#include "vec2f.hpp"

namespace lw {

class triangle {
private:
    vec3f vert[3];
    vec3f norm[3];
    vec2f uv[3];

    float areaT;
    float area(vec3f& va, vec3f& vb, vec3f& vc);

public:
    triangle();
    triangle(vec3fPtr v1, vec3fPtr v2, vec3fPtr v3);
    triangle(vec3fPtr v1, vec3fPtr v2, vec3fPtr v3, vec3fPtr n1, vec3fPtr n2,
            vec3fPtr n3);
    triangle(vec3fPtr v1, vec3fPtr v2, vec3fPtr v3,
             vec3fPtr n1, vec3fPtr n2, vec3fPtr n3,
             vec2fPtr uv1, vec2fPtr uv2, vec2fPtr uv3);

    inline vec3f& getVert(int i) {
        return vert[i];
    }

    inline vec3f& getNorm(int i) {
        return norm[i];
    }

    inline void setNorm(int i, const vec3f& _norm) {
            return norm[i] = _norm;
    }

    inline vec2f getUv(int i){
        return uv[i];
    }

    void calcWeights(vec3f& pnt, vec3f& weights);
    void getInterpVert(vec3f& weights, vec3f& vert);
    void getInterpNorm(vec3f& weights, vec3f& norm);
    void getInterpUv(vec3f& weights, vec2f& uv);

    float intersect(vec3f& p, vec3f& d);

    friend std::ostream& operator <<(std::ostream& stream, triangle& v) {
        stream << "tri {";
        stream << "v1:" << v.vert[0];
        stream << " v2: " << v.vert[1];
        stream << " v3: " << v.vert[2] << "}";
        return stream;
    }
};

typedef boost::shared_ptr<triangle> trianglePtr;

} // namespace lw
