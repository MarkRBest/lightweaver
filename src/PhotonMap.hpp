#pragma once

#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include <list>
#include <forward_list>

#include "sceen.hpp"
#include "vec3f.hpp"
#include "Utils.hpp"
#include "kdNode.hpp"

namespace lw {

#define NUM_PHOTONS (1<<20)

class Photon {
private:
public:
    vec3f m_pnt;     // impact point
    vec3f m_dir;     // inverse incident direction
    vec3f m_norm;    // surface normal
    float m_energy;  // photon energy

    Photon(vec3f& pnt, vec3f& dir, vec3f& norm, float energy) :
            m_pnt(pnt), m_dir(dir), m_energy(energy) {
    }
    inline vec3f& pnt() {
        return m_pnt;
    }
    inline vec3f& dir() {
        return m_dir;
    }
    inline vec3f& norm() {
        return m_norm;
    }
    inline float energy() {
        return m_energy;
    }

    friend std::ostream& operator <<(std::ostream& stream, Photon& v) {
        stream << "Photon {";
        stream << "pnt: " << v.m_pnt;
        stream << ", dir: " << v.m_dir;
        stream << ", norm: " << v.m_norm;
        stream << ", energy: " << v.m_energy;
        stream << "}";
        return stream;
    }
};

struct cmpPhoton {
    int plane;
    cmpPhoton(int plane) :
            plane(plane) {
    }
    ;

    inline bool operator()(const Photon& x, const Photon& y) {
        return x.m_pnt[plane] < y.m_pnt[plane];
    }
};

typedef std::vector<Photon> PhotonVec;

class PhotonMap {

private:
    PhotonVec photons;
    std::vector<trianglePtr> triStore;
    kdNode<Photon>* root;

    void unpack(sceen& s) {
        // check all meshs
        std::vector<meshPtr>::iterator it;
        std::vector<meshPtr>& meshs = s.getMeshs();
        for (it = meshs.begin(); it != meshs.end(); ++it) {
            // check each triangle in the mesh
            std::vector<trianglePtr>& tris = (*it)->getTris();
            std::vector<trianglePtr>::iterator it2;
            for (it2 = tris.begin(); it2 != tris.end(); ++it2) {
                triStore.push_back(*it2);
            }
        }
    }

    void castPhoton(vec3f& o, vec3f d, vec3f& hitPnt, vec3f& norm) {
        float t = MAX_FLOAT;
        triangle* hit = NULL;
        std::vector<trianglePtr>::iterator tIt;
        for (tIt = triStore.begin(); tIt != triStore.end(); ++tIt) {
            float tmp = (*tIt)->intersect(o, d);
            if (0.0 < tmp && tmp < t) {
                // there was an intersection
                t = tmp;
                hit = &(*(*tIt));
            }
        }
        if (hit != NULL) {
            hitPnt = o + d * t;
            vec3f w;
            hit->calcWeights(hitPnt, w);
            hit->getInterpNorm(w, norm);
        }
    }


    void generatePhotons(sceen& s) {
        // cast photons
        int32_t n = triStore.size();
        std::vector<lightPtr>& lights = s.getLights();

        std::vector<lightPtr>::iterator lIt;
        for (lIt = lights.begin(); lIt != lights.end(); ++lIt) {
            for (int32_t i = 0; i < NUM_PHOTONS; ++i) {
                int rIdx = rand() % n;
                trianglePtr tri = triStore[rIdx];
                vec3f w(randf(), randf(),randf());
                float sum = w[0] + w[1] + w[2];
                w /= sum;
                vec3f dir;
                vec3f& lightPnt = (*lIt)->pnt();
                tri->getInterpVert(w, dir);
                dir = dir-lightPnt;
                dir.norm();

                // cast the photon
                vec3f hitPnt, norm;
                // cast ray at the triangle and see what is hit
                castPhoton(lightPnt, dir, hitPnt, norm);
                // store the photon
                dir = -dir;
                photons.push_back(Photon(hitPnt, dir, norm, 1.0));
            }
        }
    }
    void findMedian(PhotonVec photons, kdNode<Photon>* node, int depth) {
//        std::cout << "D: "<<depth <<" "<< photons.size() << std::endl;
        if (photons.size() == 1) {
            Photon* p = new Photon(photons[0]);
            node->setData(p);
            node->setLeft(NULL);
            node->setRight(NULL);
        } else {
            node->setData(NULL);

            // split photons
            std::sort(photons.begin(), photons.end(), cmpPhoton(node->plane()));
            PhotonVec left, right;

            uint32_t mIdx = photons.size() / 2;
            for (uint32_t i = 0; i < mIdx; ++i) {
                left.push_back(photons[i]);
            }

            for (uint32_t i = mIdx; i < photons.size(); ++i) {
                right.push_back(photons[i]);
            }
            photons.clear();

            int32_t nPlane =(node->plane() + 1) % 3;
            node->setLeft(new kdNode<Photon>(nPlane));
            node->setRight(new kdNode<Photon>(nPlane));

            // recurse
            findMedian(left, node->left(), depth+1);
            findMedian(right, node->right(),depth+1);
        }
    }



public:

    void build(sceen& s) {
        unpack(s);
        generatePhotons(s); // cast photons
            root = new kdNode<Photon>(0);
            if(photons.size()){
                findMedian(photons, root, 0);
            }
        }

    float gatherPhotons(vec3f& hitPnt, float radius) {
        int total = 0;
        float v = 0.0;
        std::vector<kdNode<Photon>*> nodes;
        nodes.push_back(root);

        while (!nodes.empty()) {
            kdNode<Photon>* n = nodes.back();
            nodes.pop_back();

            if (n->isLeaf()) {
                Photon* p = n->data();
                vec3f dv = (p->m_pnt - hitPnt);
                float dist = dv.length();
                if (dist < radius) {
                    dist = 1.0;//(radius - dist)/radius;
                    v += dist * innerProduct(p->dir(), p->norm());
                    total++;
                }
            } else {
                if (hitPnt[n->plane()] < n->split()) {
                    nodes.push_back(n->left());
                    if (n->split() - hitPnt[n->plane()] < radius)
                        nodes.push_back(n->right());
                } else {
                    nodes.push_back(n->right());
                    if (hitPnt[n->plane()] - n->split() < radius)
                        nodes.push_back(n->left());
                }
            }
        }
        return v;
    }
};

}
