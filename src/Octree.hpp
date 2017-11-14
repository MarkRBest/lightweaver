#pragma once

#include <stdint.h>
#include <map>

#include "vec3f.hpp"
#include "triangle.hpp"
#include "ray.hpp"
#include "sceen.hpp"
#include "Utils.hpp"
#include "Timer.hpp"

namespace lw {

void printAddress(int64_t depth, int64_t address, std::string* s); // pretty print octtree address

class Octnode{
private:
    uint32_t m_depth;
    uint64_t m_address;
    vec3f m_min, m_max, m_mid; // boundaries of the octree
    uint64_t m_data; // pointer for either children or triangles
    Octnode** m_neighbours; // nieghbour octrees for 3d-dda

public:
    Octnode() : m_depth(0), m_address(0), m_min(), m_max(), m_mid(),
    m_data(0), m_neighbours(NULL) {}
    Octnode(uint32_t depth, uint64_t address, const vec3f& min, const vec3f& max) :
            m_depth(depth), m_address(address), m_min(min), m_max(max), m_mid(),
            m_data(0), m_neighbours(NULL)
    {
        m_mid = min;
        m_mid += max;
        m_mid /= 2.0;
    }

    inline uint32_t depth() const { return m_depth; }
    inline uint64_t address() const { return m_address; }
    inline void setDepth(uint64_t depth) { m_depth = depth; }
    inline void setAddress(uint64_t address) { m_address = address; }
    inline vec3f& min(){return m_min;}
    inline vec3f& max(){return m_max;}
    inline vec3f& mid(){return m_mid;}
    inline void setMinMax(const vec3f& min,const vec3f& max){m_min =min; m_max=max; m_mid=(m_min+m_max)/(float)2.0;}

    inline Octnode** children() {return (Octnode**)(m_data);}
    inline std::vector<triangle>* tris() const { return (std::vector<triangle>*)(m_data); }
    inline void setData(Octnode** data){m_data= (uint64_t)data;}
    inline void setData(std::vector<triangle>* data){m_data= (uint64_t)data;}

    inline bool isLeaf() const {return m_neighbours != NULL;}
    inline Octnode** neighbours(){return m_neighbours;}
    inline void setNeighbours(Octnode** neighbour){m_neighbours = neighbour;}

    float intersectRay(ray& ray);
    bool intersectTri(triangle& tri);
    int planeBoxOverlap(float normal[3], float vert[3], vec3f& maxbox);

    friend std::ostream& operator <<(std::ostream& stream, const Octnode& v);
};

class OctTree {

public:

    static int64_t masks[3];// = {1,2,4};

    OctTree (uint64_t _max_depth=0, uint64_t _max_tris=0) : max_depth(_max_depth), max_tris(_max_tris) { };
    void build(sceen& sceen, uint64_t _max_depth, uint64_t _max_tris);
    Octnode* getNode(int64_t depth, int64_t address);
    Octnode* getLeafNode(vec3f& pnt);
    int64_t shiftAddr(int64_t address, int mask, bool fwd);

    float raytrace(vec3f& orig, vec3f& dir, Octnode* camnode, triangle*& hittri);

private:
    uint64_t max_depth;
    uint64_t max_tris;

    float dists[16][3];
    float ray_times[16][3];
    Octnode root;


    void buildRootnode(Octnode& root, sceen& sceen);
    void buildOctree(Octnode* node, uint64_t depth, uint64_t address, std::vector<triangle>* tris);
    void setNeighbours(Octnode* node, int64_t address);
};

} // namespace lw
