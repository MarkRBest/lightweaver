#include "sceen.hpp"

void lw::sceen::addCamera(cameraPtr _c)
{
 c=_c;
}

lw::cameraPtr lw::sceen::getCam()
{
    return c;
}

void lw::sceen::addLight(lightPtr _l)
{
 l.push_back(_l);
}

void lw::sceen::addMesh(meshPtr _m)
{
  m.push_back(_m);
}