#include "triangle.hpp"

#include <iostream>

using lw::vec3f;

lw::triangle::triangle() : areaT(0) {}

lw::triangle::triangle(vec3fPtr v1, vec3fPtr v2, vec3fPtr v3) {
    for (int i = 0; i < 3; i++) {
        vert[0][i] = (*v1)[i];
        vert[1][i] = (*v2)[i];
        vert[2][i] = (*v3)[i];
    }
    areaT = area(vert[0],vert[1],vert[2]);
}

lw::triangle::triangle(vec3fPtr v1, vec3fPtr v2, vec3fPtr v3, vec3fPtr n1,
        vec3fPtr n2, vec3fPtr n3) {
    for (int i = 0; i < 3; i++) {
        vert[0][i] = (*v1)[i];
        vert[1][i] = (*v2)[i];
        vert[2][i] = (*v3)[i];

        norm[0][i] = (*n1)[i];
        norm[1][i] = (*n2)[i];
        norm[2][i] = (*n3)[i];
    }
    areaT = area(vert[0],vert[1],vert[2]);
}

lw::triangle::triangle(
        vec3fPtr v1, vec3fPtr v2, vec3fPtr v3,
        vec3fPtr n1, vec3fPtr n2, vec3fPtr n3,
        vec2fPtr uv1, vec2fPtr uv2, vec2fPtr uv3) {
    for (int i = 0; i < 3; i++) {
        vert[0][i] = (*v1)[i];
        vert[1][i] = (*v2)[i];
        vert[2][i] = (*v3)[i];

        norm[0][i] = (*n1)[i];
        norm[1][i] = (*n2)[i];
        norm[2][i] = (*n3)[i];
    }
    // copy uv coordinates
    for (int i = 0; i < 2; i++) {
        uv[0][i] = (*uv1)[i];
        uv[1][i] = (*uv2)[i];
        uv[2][i] = (*uv3)[i];
    }
    areaT = area(vert[0],vert[1],vert[2]);
}

#define innerProduct(v,q) \
        ((v)[0] * (q)[0] + \
        (v)[1] * (q)[1] + \
        (v)[2] * (q)[2])

#define crossProduct(a,b,c) \
    (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
    (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
    (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

/* a = b - c */
#define vector(a,b,c) \
    (a)[0] = (b)[0] - (c)[0];	\
    (a)[1] = (b)[1] - (c)[1];	\
    (a)[2] = (b)[2] - (c)[2];


float lw::triangle::intersect(vec3f& p, vec3f& d) {
    //return 0.0;
    float h[3], q[3];
    float a, f, u, v;

    vec3f e1 = vert[1] - vert[0];
    vec3f e2 = vert[2] - vert[0];

    crossProduct (h, d, e2);
    a = innerProduct(e1,h);

    if (a > -0.00001 && a < 0.00001)
        return (false);

    f = 1 / a;
    vec3f s = p - vert[0];
    u = f * (innerProduct(s,h));

    if (u < 0.0 || u > 1.0)
        return (false);

    crossProduct(q, s, e1);
    v = f * innerProduct(d,q);

    if (v < 0.0 || u + v > 1.0)
        return (false);

    // at this stage we can compute t to find out where
    // the intersection point is on the line
    return f * innerProduct(e2,q);
}

float lw::triangle::area(vec3f& va, vec3f& vb, vec3f& vc) {
    vec3f ea = vb - va;
    vec3f eb = vc - va;
    vec3f ec = vc - vb;
    float a = innerProduct(ea, ea);
    float b = innerProduct(eb, eb);
    float c = innerProduct(ec, ec);
    float area = (2 * a * b + 2 * b * c + 2 * c * a - (a * a) - (b * b)
            - (c * c)) / 16.0;
    return pow(area, 0.5);
}

void lw::triangle::calcWeights(vec3f& pnt, vec3f& weights) {
    float areaA, areaB, areaC;
    areaB = area(vert[0], pnt, vert[2]);
    if (std::isnan(areaB))
        areaB = 0.0;
    areaC = area(vert[0], pnt, vert[1]);
    if (std::isnan(areaC))
        areaC = 0.0;
    areaA = areaT - areaB - areaC;

    weights[0] = areaA / areaT;
    weights[1] = areaB / areaT;
    weights[2] = areaC / areaT;

//    std::cout << areaT << " "<< areaA << " " << areaB << " " << areaC << std::endl;
}

void lw::triangle::getInterpVert(vec3f& weights, vec3f& _vert) {
    for (int i = 0; i < 3; i++){
        _vert[i] = 0.0;
        for (int j = 0; j < 3; j++){
            _vert[i] += weights[j] * vert[j][i];
        }
    }
}


void lw::triangle::getInterpNorm(vec3f& weights, vec3f& _norm) {

    for (int i = 0; i < 3; i++){
        _norm[i] = 0.0;
        for (int j = 0; j < 3; j++){
            _norm[i] += weights[j] * norm[j][i];
        }
    }
}

void lw::triangle::getInterpUv(vec3f& weights, vec2f& _uv) {
    for (int i = 0; i < 2; i++) {
        _uv[i] = 0.0;
        for (int j = 0; j < 3; j++) {
            _uv[i] += weights[j] * uv[j][i];
        }
    }
}
