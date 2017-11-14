#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <list>
#include <map>
#include <string>
#include <algorithm>

#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>

#include <boost/make_shared.hpp>
#include <boost/limits.hpp>

#include "vec3f.hpp"
#include "ray.hpp"
#include "light.hpp"
#include "mesh.hpp"
#include "sceen.hpp"
#include "camera.hpp"
#include "triangle.hpp"
#include "texture.hpp"
#include "Timer.hpp"
#include "Octree.hpp"
#include "PhotonMap.hpp"
#include "Utils.hpp"
#include "bsptree.hpp"

using namespace lw;

sceen s;
Octnode root;
PhotonMap pMap;
Texture* tex1;
BSPTree tree;
OctTree octtree;

#define HARD_SHADOWS 1
#define SOFT_SHADOWS 2
#define SHADOWS HARD_SHADOWS

void loadSceen(sceen& s, const char* f) {
    std::ifstream myfile;
    std::string line;
    myfile.open(f);
    int64_t ntris =0;
    if (myfile.is_open()) {
        meshPtr m;
        int index = 0;
        vec3fPtr verts[3];
        vec3fPtr norms[3];
        vec2fPtr uv[3];

        while (myfile.good()) {
            std::getline(myfile, line);
            std::stringstream ss(line);
            std::string tok;
            ss >> tok;
            if (!strcmp(tok.c_str(), "Camera")) {
                std::string tmp;
                float ax, ay, az=0.0, x, y, z;
                ss >> tmp >> ax >> tmp >> ay;
                ss >> tmp >> x >> tmp >> y >> tmp >> z;
                //std::cout << camera << std::endl;
                s.addCamera(boost::make_shared<camera>(ax, ay, az, x, y, z));
            } else if (!strcmp(tok.c_str(), "Light")) {
                std::string tmp;
                float x, y, z;
                ss >> tmp >> x >> tmp >> y >> tmp >> z;
                s.addLight(boost::make_shared<light>(x, y, z));
            } else if (!strcmp(tok.c_str(), "Mesh")) {
                m = boost::make_shared<mesh>();
            } else if (!strcmp(tok.c_str(), "Tri:")) {
                index = 0;
                ntris ++;
            } else if (!strcmp(tok.c_str(), "MeshEnd")) {
                //std::cout<< *m <<std::endl;
                s.addMesh(m);
            } else if (!strcmp(tok.c_str(), "vert:")) {
                std::string tmp;
                float vx, vy, vz;
                ss >> tmp >> vx >> tmp >> vy >> tmp>> vz;
                verts[index] = boost::make_shared<vec3f>(vx, vy, vz);
                float nx, ny, nz;
                ss >> tmp >> nx >> tmp >> ny >> tmp>> nz;
                norms[index] = boost::make_shared<vec3f>(nx, ny, nz);
                float tu, tv;
                ss >> tmp >> tu >> tmp >> tv;
                uv[index] = boost::make_shared<vec2f>(tu, tv);
                index++;
            }

            if (index > 2) {
                trianglePtr tri = boost::make_shared<triangle>(
                        verts[0], verts[1], verts[2],
                        norms[0], norms[1], norms[2],
                        uv[0], uv[1], uv[2]);
                m->addTriangle(tri);
            }
        }
    }
    std::cout << "Loaded: " << f << " " << s.getMeshs().size() << " Meshes";
    std::cout << " " << s.getLights().size() << " Lights";
    std::cout << " " << ntris << " Triangles" << std::endl;

    myfile.close();
}

typedef std::map<vec3f, std::vector<vec3f> > NormsMap;
NormsMap normsMap;

void smoothNormals(sceen& sceen) {
    std::vector<lw::meshPtr>& meshList = sceen.getMeshs();
    std::vector<lw::meshPtr>::iterator mit;
    for (mit = meshList.begin(); mit != meshList.end(); ++mit) {
        std::vector<trianglePtr>& tris = (*mit)->getTris();
        std::vector<lw::trianglePtr>::iterator tIter;
        for (tIter = tris.begin(); tIter != tris.end(); ++tIter) {
            trianglePtr& tri = *tIter;
            for (int i = 0; i < 3; i++) {
                vec3f& vert = tri->getVert(i);
                NormsMap::iterator nmIter = normsMap.find(vert);

                // if there is no current norm list then add empty one
                if (nmIter == normsMap.end()) {
                    nmIter = normsMap.insert(nmIter,
                            std::pair<vec3f, std::vector<vec3f> >(vert,
                                    std::vector<vec3f>()));
                }

                nmIter->second.push_back(tri->getNorm(i));
            }
        }
    }

    // average normals
    NormsMap::iterator nmIter;
    for(nmIter = normsMap.begin(); nmIter != normsMap.end(); ++nmIter){
        std::vector<vec3f>& normList =  nmIter->second;

        // compute the average normal at the vertex
        vec3f avg(0,0,0);
        std::vector<vec3f>::iterator nIter;
        for(nIter= normList.begin(); nIter != normList.end(); ++nIter){
            avg+=*nIter;
        }
        avg /= normList.size();
        normList.clear();

        // replace list with the average of the normals
        normList.push_back(avg);
    }

    // apply averaged normals
    for (mit = meshList.begin(); mit != meshList.end(); ++mit) {
        std::vector<trianglePtr>& tris = (*mit)->getTris();
        std::vector<lw::trianglePtr>::iterator tIter;
        for (tIter = tris.begin(); tIter != tris.end(); ++tIter) {
            trianglePtr& tri = *tIter;
            for (int i = 0; i < 3; i++) {
                vec3f& vert = tri->getVert(i);
                NormsMap::iterator nmIter = normsMap.find(vert);
                tri->setNorm(i, nmIter->second[0]);
            }
        }
    }
}

triangle* shadowCache = NULL;

bool notInShadow(vec3f& orig, vec3f& dir) {
    float t = MAX_FLOAT;
    if(shadowCache){
        float tmp = shadowCache->intersect(orig, dir);
        if (0.0 < tmp && tmp < t) {
            return false;
        }
    }

//    triangle* hit = NULL;
//    t = octtree.raytrace(orig, dir, camNode, hit);
//
//    return t == MAX_FLOAT ? false : true;

//    // check all meshs
//    std::vector<meshPtr>& meshs = s.getMeshs();
//    std::vector<meshPtr>::iterator it;
//    for (it = meshs.begin(); it != meshs.end(); ++it) {
//        std::vector<trianglePtr>& tris = (*it)->getTris();
//        // check each triangle in the mesh
//        for (std::vector<trianglePtr>::iterator it2 = tris.begin();
//                it2 != tris.end(); ++it2) {
//            float tmp = (*it2)->intersect(orig, dir);
//            if (0.0 < tmp && tmp < t) {
//                // there was an intersection
//                t = tmp;
//                shadowCache = &(*(*it2));
//                return false;
//            }
//        }
//    }

    return true;
}

inline float frand(){
    return ((float) rand()/(float) RAND_MAX) - 0.5;
}

void raytrace(vec3f& c, ray& r, Octnode* camNode) {
    float t = MAX_FLOAT;
//    triangle* hit = NULL;

//    trianglePtr hit = boost::make_shared<triangle>(triangle());
//    t = tree.raytrace(r.orig(), r.dir(), 0, 0, hit);
//    if (t > 0) { std::cout << "bsptree hit " << time << " " << *hit << std::endl; }


    triangle* hit = NULL;
    t = octtree.raytrace(r.orig(), r.dir(), camNode, hit);
//    if (t > MAX_FLOAT) { std::cout << "octtree hit " << t << " " << *hit << std::endl; }

    // clear pixel
    c[0] = c[1] = c[2] = 0.0;

    // shade pixel if hit
    if (t != MAX_FLOAT) {
        std::vector<lightPtr>& lights = s.getLights();
        std::vector<lightPtr>::iterator light_iter;
        for(light_iter = lights.begin(); light_iter != lights.end(); ++light_iter)
        {

            lightPtr& light = *light_iter;

            vec3f hitPnt = r.dir();
            hitPnt *= t-0.0001;
            hitPnt += r.orig();

            vec3f lDir = light->pnt() - hitPnt;
            lDir.norm();

            vec3f w, triNorm;
            vec2f uv;
            hit->calcWeights(hitPnt, w);
            //w.norm();
            hit->getInterpNorm(w, triNorm);
//            hit->getInterpUv(w, uv);
//            triNorm = hit->getNorm(0);
            // get colour from uv
    //        int h_idx = (tex1->height()-1) * uv[0];
    //        int w_idx = (tex1->width()-1) * uv[1];
    //        unsigned char* tex = tex1-> pixel(h_idx,w_idx);

            float ilum = 0.05; // ambient
            vec3f slDir;

            float sCount = 0;

    #if SHADOWS == SOFT_SHADOWS
            const int NumTestRays = 25, NumFullRays = 250;

            int countRays = 0;
            for (int i = 0; i < NumTestRays; i++) {
                slDir = light->pnt();
                slDir[0] += frand()/2;
                slDir[1] += frand()/2;
                slDir[2] += frand()/2;
                slDir -= hitPnt;
                if (notInShadow(hitPnt, slDir)) {
                    sCount++;
                }
                countRays++;
            }

            // if part in shadow then fully test with a lot of rays
            if (sCount > 0 && sCount < NumTestRays) {
                for (int i = NumTestRays; i < NumFullRays; i++) {
                    slDir = light->pnt();
                    slDir[0] += frand()/2;
                    slDir[1] += frand()/2;
                    slDir[2] += frand()/2;
                    slDir -= hitPnt;
                    if (notInShadow(hitPnt, slDir)) {
                        sCount++;
                    }
                    countRays++;
                }
            }
            sCount /= countRays;
    #else
            sCount = 1.0;
//            if (notInShadow(hitPnt, lDir)) {
//                sCount++;
//            }
    #endif
            // shade the image
            //float pVal = 1.0;//pMap.gatherPhotons(hitPnt, 1.0);

            float gourad = innerProduct(lDir, triNorm);
            if (gourad > 0.0) {
                ilum += 0.8 / s.getLights().size() * gourad * sCount;
            }

    //        std::cout << ":" << pVal << " " << gourad <<std::endl;

            c[0] += ilum;//* ((float)tex[0]/256);
            c[1] += ilum;//* ((float)tex[1]/256);
            c[2] += ilum;//* ((float)tex[2]/256);
        }

        c[0] = std::min(c[0], (float)0.9999999);
        c[1] = std::min(c[1], (float)0.9999999);
        c[2] = std::min(c[2], (float)0.9999999);
    }
}

int scale = 5;
int window_width = 128*scale, window_height = 128*scale;

struct RGB {
    uint8_t r, g, b;
};

void draw() {
    lw::Timer timer;
    //std::cout << "Casting photons"<<std::endl;
    //pMap.build(s);
    //std::cout << "Photon map built in "<< timer.elapsedTime() / 1000000 <<" ms" <<std::endl;


    cameraPtr cam = s.getCam();
    std::cout << "Camera " << *cam << std::endl;
    Octnode* camNode = octtree.getLeafNode(cam->pnt());
    //std::cout << cam->pnt() << std::endl;
    std::string camadr;
    printAddress(camNode->depth(), camNode->address(), &camadr);
    std::cout << "Cam Node " << camadr << " " << *camNode << std::endl;

    vec3f& o = cam->pnt();
    float ax = cam->ax(), ay = cam->ay();

    vec3f c;
    ray r(o[0], o[1], o[2], 0, 0, 0);

//    float pi_6 = 3.14159265359 / 4.0;
//    float x_step = pi_6 / (window_width / 2.00);
//    float y_step = pi_6 / (window_height / 2.00);
    timer.reset();

    GLubyte pixels[window_width][window_height][3];

    vec3f up (0, 1, 0);
//    vec3f camDir (-(cos(ax) + cos(ay)) / 2.0, sin(ax), -sin(ay));
    vec3f camDir (-cam->pnt()[0], -cam->pnt()[1], -cam->pnt()[2]);
//    vec3f camDir (cos(cam->ax()), cos(cam->ay()), cos(cam->az()));
    camDir.norm();
    std::cout << "cam pos " << cam->pnt() << std::endl;
    std::cout << "cam dir " << camDir << std::endl;
    vec3f side;
    crossProduct(side , camDir, up);
    crossProduct(up , camDir, side);
    vec3f pnt = cam->pnt() + camDir;
    vec3f stepx = (side/(int)(window_width/1.0));
    vec3f stepy = (up/(int)(window_height/1.0));
    for(int i=0;i<(window_width/2);i++)
        pnt -= stepx;
    for(int i=0;i<(window_height/2);i++)
        pnt -= stepy;

    for (int i = 0; i < window_width; ++i) {
        std::cout << "raster " << i << " of " << window_width << std::endl;
        for (int j = 0; j < window_height; ++j) {
            vec3f dir = pnt + stepx*i + stepy*j - cam->pnt();
            dir.norm();
            r.setDir(dir[0], dir[1], dir[2]);
            c[0]=0.0; c[1]=0.0; c[2]=0.0;

            raytrace(c, r, camNode);

            int j_idx = window_height-j-1;
            int i_idx = window_width-i-1;
            pixels[i_idx][j_idx][0] = GLubyte(c[0] * 256);
            pixels[i_idx][j_idx][1] = GLubyte(c[1] * 256);
            pixels[i_idx][j_idx][2] = GLubyte(c[2] * 256);
        }
    }

    std::cout <<"Time: " << timer.elapsedTime() / 1000000 <<std::endl;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawPixels(window_width, window_height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glutSwapBuffers();
}

int main(int argc, const char* argv[]) {
    glutInit(&argc, const_cast<char**>(argv));

    lw::Timer timer;
    std::string f;

    if (argc != 2){
        std::cout << "usage LightWeaver <filename>" << std::endl;
        exit(0);
    }
    f = std::string(argv[1]);

    loadSceen(s, f.c_str());
    octtree.build(s, 0, 10);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("Lightweaver Raytracer");

    glutDisplayFunc(draw);

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0, 0.0, 0.0, 1.0);

    glutMainLoop();

    return 0;
}
