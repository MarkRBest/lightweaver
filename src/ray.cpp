#include "ray.hpp"

lw::ray::ray(float o0,float o1, float o2, float d0, float d1, float d2)
:o(o0,o1,o2),d(d0,d1,d2) {}

lw::ray::ray(float* _o, float* _d)
:o(_o[0],_o[1],_o[2]),d(_d[0],_d[1],_d[2]){}

void lw::ray::setOrig(float x, float y, float z)
{
    o[0]=x; o[1]=y; o[2]=z;
}

void lw::ray::setOrig(vec3f& _orig)
{
    for(int i=0;i<3;i++)
        o[i] = _orig[i];
}

void lw::ray::setDir(float x, float y, float z)
{
    d[0]=x; d[1]=y; d[2]=z;
}

void lw::ray::setDir(vec3f& _dir)
{
    for(int i=0;i<3;i++)
        d[i] = _dir[i];
}
