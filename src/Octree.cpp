#include "Octree.hpp"

namespace lw {

int64_t OctTree::masks[3] = {1,2,4};

void printAddress(int64_t depth, int64_t address, std::string* s) {

    if (address == -1)
    {
        s->append("-1");
    }
    else
    {
        while (depth > 0)
        {
            for (int i = 0; i < 3; ++i)
            {
                int mask = 1 << i;
                s->append(((mask & address) ? "1" : "0"));
            }
            s->append("|");
            address >>= 3;
            depth--;
        }
    }
}

#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
        dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
        dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
        dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
        dest[0]=v1[0]-v2[0]; \
        dest[1]=v1[1]-v2[1]; \
        dest[2]=v1[2]-v2[2];

#define FINDMINMAX(x0,x1,x2,min,max) \
        min = max = x0;   \
        if(x1<min) min=x1;\
        if(x1>max) max=x1;\
        if(x2<min) min=x2;\
        if(x2>max) max=x2;


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)             \
        p0 = a*v0[Y] - b*v0[Z];                    \
        p2 = a*v2[Y] - b*v2[Z];                    \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
        rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)              \
        p0 = a*v0[Y] - b*v0[Z];                    \
        p1 = a*v1[Y] - b*v1[Z];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
        rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)             \
        p0 = -a*v0[X] + b*v0[Z];                   \
        p2 = -a*v2[X] + b*v2[Z];                       \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)              \
        p0 = -a*v0[X] + b*v0[Z];                   \
        p1 = -a*v1[X] + b*v1[Z];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
        if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/
#define AXISTEST_Z12(a, b, fa, fb)             \
        p1 = a*v1[X] - b*v1[Y];                    \
        p2 = a*v2[X] - b*v2[Y];                    \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
        if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)              \
        p0 = a*v0[X] - b*v0[Y];                \
        p1 = a*v1[X] - b*v1[Y];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
        rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
        if(min>rad || max<-rad) return 0;


typedef uint32_t udword;
//! Integer representation of a floating-point value.
#define IR(x)   ((udword&)x)
#define RAYAABB_EPSILON 0.0000001

void swap(float& tmin, float& tmax){
    float tmp = tmin;
    tmin = tmax;
    tmax = tmp;
}

float Octnode::intersectRay(ray& ray) {
    /*bool Inside = true;
    float MaxT[3], coord[3];

    vec3f& origin = ray.orig(), dir = ray.dir();

    MaxT[0] = MaxT[1] = MaxT[2] = -1.0f;

    // Find candidate planes.
    for (int i = 0; i < 3; i++) {
        if (origin[i] < m_min[i]) {
            coord[i] = m_min[i];
            Inside = false;

            // Calculate T distances to candidate planes
            if (dir[i] > 0)
                MaxT[i] = (m_min[i] - origin[i]) / dir[i];
        } else if (origin[i] > m_max[i]) {
            coord[i] = m_max[i];
            Inside = false;
            // Calculate T distances to candidate planes
            if (dir[i] < 0)
                MaxT[i] = (origin[i] - m_max[i]) / -dir[i];
        }
    }

    // Ray origin inside bounding box
    if (Inside) {
        return 0.0;
    }

    // Get largest of the maxT's for final choice of intersection
    udword WhichPlane = 0;
    if (MaxT[1] > MaxT[WhichPlane])
        WhichPlane = 1;
    if (MaxT[2] > MaxT[WhichPlane])
        WhichPlane = 2;

    // Check final candidate actually inside box
    if (MaxT[WhichPlane] < 0)
        return -1.0f;

    for (uint32_t i = 0; i < 3; ++i) {
        if (i != WhichPlane) {
            coord[i] = origin[i] + MaxT[WhichPlane] * dir[i];
#ifdef RAYAABB_EPSILON
            if(coord[i] < m_min[i] - RAYAABB_EPSILON || coord[i] > m_max[i] + RAYAABB_EPSILON) return -1.0f;
#else
            if (coord[i] < m_min[i] || coord[i] > m_max[i])
                return -1.0f;
#endif
        }
    }
    return MaxT[WhichPlane];    // ray hits box
     */
    vec3f&  orig = ray.orig();
    vec3f& dir = ray.dir();

    float tmin = (m_min[0] - orig[0]) / dir[0];
    float tmax = (m_max[0] - orig[0]) / dir[0];
    if (tmin > tmax)
        swap(tmin, tmax);
    float tymin = (m_min[1] - orig[1]) / dir[1];
    float tymax = (m_max[1] - orig[1]) / dir[1];
    if (tymin > tymax)
        swap(tymin, tymax);
    if ((tmin > tymax) || (tymin > tmax))
        return -1.0;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;
    float tzmin = (m_min[2] - orig[2]) / dir[2];
    float tzmax = (m_max[2] - orig[2]) / dir[2];
    if (tzmin > tzmax)
        swap(tzmin, tzmax);
    if ((tmin > tzmax) || (tzmin > tmax))
        return -1.0;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    //if ((tmin > r.tmax) || (tmax < r.tmin)) return false;
    //if (r.tmin < tmin) r.tmin = tmin;
    //if (r.tmax > tmax) r.tmax = tmax;
    return tmin;
}



std::ostream& operator <<(std::ostream& stream, const Octnode& v)
{
    std::string s;
    printAddress(v.depth(), v.address(), &s);
    stream << "Octnode{";
    stream << "depth: " << v.depth() << " ";
    stream << "adr: " << s << " ";
    stream <<" min: (" << v.m_min[0] << ", "<< v.m_min[1] << ", "<< v.m_min[2] << ")";
    stream <<" max: (" << v.m_max[0] << ", "<< v.m_max[1] << ", "<< v.m_max[2] << ")";
    stream <<" mid: (" << v.m_mid[0] << ", "<< v.m_mid[1] << ", "<< v.m_mid[2] << ")";
    if (v.isLeaf())
        stream << " leaf: 1, tris: " << v.tris()->size();
    else
        stream << " leaf: 0, tris: 0";
    stream << "}";
    return stream;
}


bool Octnode::intersectTri(triangle& tri)
{
    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */

    vec3f boxcenter((m_min[0]+m_max[0])/2,(m_min[1]+m_max[1])/2,(m_min[2]+m_max[2])/2);
    vec3f boxhalfsize ((m_max[0]-m_min[0])/2,(m_max[1]-m_min[1])/2,(m_max[2]-m_min[2])/2);
    /*   float triverts[3][3] = {
       {tri.getVert(0)[0],tri->vertex[0].pos.y,tri->vertex[0].pos.z},
       {tri->vertex[1].pos.x,tri->vertex[1].pos.y,tri->vertex[1].pos.z},
       {tri->vertex[2].pos.x,tri->vertex[2].pos.y,tri->vertex[2].pos.z}
   };
     */
    float v0[3],v1[3],v2[3];
    float min,max,p0,p1,p2,rad,fex,fey,fez;      // -NJMP- "d" local variable removed
    float normal[3],e0[3],e1[3],e2[3];

    /* This is the fastest branch on Sun */
    /* move everything so that the boxcenter is in (0,0,0) */
    SUB(v0,tri.getVert(0),boxcenter);
    SUB(v1,tri.getVert(1),boxcenter);
    SUB(v2,tri.getVert(2),boxcenter);

    /* compute triangle edges */
    SUB(e0,v1,v0);      /* tri edge 0 */
    SUB(e1,v2,v1);      /* tri edge 1 */
    SUB(e2,v0,v2);      /* tri edge 2 */

    /* Bullet 3:  */
    /*  test the 9 tests first (this was faster) */
    fex = fabsf(e0[X]);
    fey = fabsf(e0[Y]);
    fez = fabsf(e0[Z]);

    AXISTEST_X01(e0[Z], e0[Y], fez, fey);
    AXISTEST_Y02(e0[Z], e0[X], fez, fex);
    AXISTEST_Z12(e0[Y], e0[X], fey, fex);

    fex = fabsf(e1[X]);
    fey = fabsf(e1[Y]);
    fez = fabsf(e1[Z]);

    AXISTEST_X01(e1[Z], e1[Y], fez, fey);
    AXISTEST_Y02(e1[Z], e1[X], fez, fex);
    AXISTEST_Z0(e1[Y], e1[X], fey, fex);

    fex = fabsf(e2[X]);
    fey = fabsf(e2[Y]);
    fez = fabsf(e2[Z]);

    AXISTEST_X2(e2[Z], e2[Y], fez, fey);
    AXISTEST_Y1(e2[Z], e2[X], fez, fex);
    AXISTEST_Z12(e2[Y], e2[X], fey, fex);

    /* Bullet 1: */
    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */
    /* test in X-direction */
    FINDMINMAX(v0[X],v1[X],v2[X],min,max);
    if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;

    /* test in Y-direction */
    FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
    if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;

    /* test in Z-direction */
    FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
    if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;

    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
    CROSS(normal,e0,e1);
    // -NJMP- (line removed here)
    if(!planeBoxOverlap(normal,v0, boxhalfsize)) return 0;    // -NJMP-
    return 1;   /* box and triangle overlaps */
}

int Octnode::planeBoxOverlap(float normal[3], float vert[3],
        vec3f& maxbox) {  // -NJMP-
    int q;
    float vmin[3], vmax[3], v;
    for (q = X; q <= Z; q++) {
        v = vert[q];                  // -NJMP-
        if (normal[q] > 0.0f) {
            vmin[q] = -maxbox[q] - v;   // -NJMP-
            vmax[q] = maxbox[q] - v;   // -NJMP-
        } else {
            vmin[q] = maxbox[q] - v;   // -NJMP-
            vmax[q] = -maxbox[q] - v;   // -NJMP-
        }
    }
    if (DOT(normal,vmin) > 0.0f)
        return 0;   // -NJMP-
    if (DOT(normal,vmax) >= 0.0f)
        return 1;  // -NJMP-
    return 0;
}



void OctTree::build(sceen& sceen, uint64_t _max_depth, uint64_t _max_tris)
{
    Timer timer;

    max_depth = _max_depth;
    max_tris =_max_tris;

    buildRootnode(root, sceen);

    vec3f dist = root.max() - root.min();
    for(uint64_t i=0; i<=max_depth; ++i)
    {
        std::cout << "dist " << i << " " << dist <<std::endl;
        dists[i][0] = dist[0];
        dists[i][1] = dist[1];
        dists[i][2] = dist[2];
        dist /= 2.0;
    }

    buildOctree(&root, 0, 0, root.tris());
    setNeighbours(&root, 0);

    std::cout << "Octree Built in "<< timer.elapsedTime()/1000000 << " ms"<< std::endl;
}


// get the address of adjacent node
int64_t OctTree::shiftAddr(int64_t address, int mask, bool fwd)
{
    int addrCpy = address;
    do
    {
        if (fwd)
        {
            if (address & mask)
            {
                address -= mask;
            }
            else
            {
                address += mask;
                return address;
            }
            mask <<= 3;
            addrCpy >>= 3;
        }
        else
        {
            if (address & mask)
            {
                address -= mask;
                return address;
            }
            else
            {
                address += mask;
            }
            mask <<= 3;
            addrCpy >>= 3;
        }
    }
    while (addrCpy != 0);
    return -1;
}

Octnode* OctTree::getLeafNode(vec3f& pnt) {
    Octnode* node = &root;

    while (!node->isLeaf()) {
        int index = 0;

        //Octnode** children = node->children();
        vec3f& mid = node->children()[0]->max();
        for (int i = 0; i < 3; ++i)
            if (pnt[i] > mid[i])
                index += (1 << i);

        node = node->children()[index];
    }
    return node;
}

// get a node from a root node given an address
Octnode* OctTree::getNode(int64_t depth, int64_t address)
{
    Octnode* node = &root;

    if(address == -1){
        return NULL;
    }

    // todo - can be sped up with  proper masking
    std::vector<int64_t> tmp;
    for(int i=0; i<depth; ++i) {
        tmp.push_back((address & 7));
        address >>= 3;
    }

    while(!node->isLeaf() && depth > 0){
        node = node->children()[tmp[depth-1]];
        depth--;
    }
    return node;
}


void OctTree::buildRootnode(Octnode& root, sceen& sceen)
{
    std::vector<triangle>* rootTris = new std::vector<triangle>();
    vec3f min(MAX_FLOAT, MAX_FLOAT, MAX_FLOAT);
    vec3f max(MIN_FLOAT, MIN_FLOAT, MIN_FLOAT);

    std::vector<meshPtr>& meshList = sceen.getMeshs();
    std::vector<meshPtr>::iterator mit;
    for (mit = meshList.begin(); mit != meshList.end(); ++mit) {
        std::vector<trianglePtr>& tris = (*mit)->getTris();
        std::vector<trianglePtr>::iterator tIter;
        for (tIter = tris.begin(); tIter != tris.end(); ++tIter) {
            trianglePtr& tri = *tIter;

            // record triangle
            rootTris->push_back(*tri);

            // find min/max dimensions
            for (int i = 0; i < 3; i++) {
                vec3f& vert = tri->getVert(i);
                for (int j = 0; j < 3; j++) {
                    if (min[j] > vert[j])
                        min[j] = vert[j];
                }
                for (int j = 0; j < 3; j++) {
                    if (max[j] < vert[j])
                        max[j] = vert[j];
                }
            }
        }
    }

    // check boundaries for the camera
    vec3f& campnt = sceen.getCam()->pnt();
    for(uint64_t i=0; i<3; ++i){
        if(min[i] > campnt[i])
            min[i] = campnt[i];
        if(max[i] < campnt[i])
            max[i] = campnt[i];
    }

    for(uint64_t i=0; i<3; ++i)
    {
        min[i]-= 1.0;
        max[i]+= 1.0;
    }

    root.setMinMax(min,max);
    root.setData(rootTris);
    root.setNeighbours(NULL);
}

void OctTree::buildOctree(Octnode* node, uint64_t depth, uint64_t address, std::vector<triangle>* tris)
{
    node->setDepth(depth);
    node->setAddress(address);

    if (depth >= max_depth) // || tris->size() <= max_tris)
    {
        node->setData(tris);

        // mark node as a leaf
        Octnode** n = new Octnode*[6];
        node->setNeighbours(n);
    }
    else
    {
        // not a leaf
        node->setNeighbours(NULL);

        Octnode** children = new Octnode*[8];
        node->setData(children);

        for (int i = 0; i < 8; ++i) {
            children[i] = new Octnode();
            vec3f min, max;

            for (int j = 0; j < 3; ++j)
            {
                int64_t mask = 1 << j;
                if (i & mask)
                {
                    min[j] = node->mid()[j];
                    max[j] = node->max()[j];
                }
                else
                {
                    min[j] = node->min()[j];
                    max[j] = node->mid()[j];
                }
            }

            children[i]->setMinMax(min, max);

            // split the triangles
            std::vector<triangle>* subset = new std::vector<triangle>();
            std::vector<triangle>::iterator tIter;
            for (tIter = tris->begin(); tIter != tris->end(); ++tIter)
            {
                if (children[i]->intersectTri(*tIter))
                {
                    subset->push_back(*tIter);
                }
            }

            // recurse
            buildOctree(children[i], depth + 1, (address<<3)+i, subset);
        }

        // full list of triangles is not needed anymore so free it up
        delete (tris);
    }

    std::cout << "Depth: " << depth << " ";
    std::string adr;
    printAddress(depth, address, &adr);
    std::cout << "Adr: "<< adr << " ";
    std::cout << tris->size() << " " << *node << std::endl;
}

void OctTree::setNeighbours(Octnode* node, int64_t address)
{
    if (node->isLeaf())
    {
        Octnode** n = node->neighbours();

        n[0] = getNode(node->depth(), shiftAddr(address, 1, 0));
        n[1] = getNode(node->depth(), shiftAddr(address, 1, 1));
        n[2] = getNode(node->depth(), shiftAddr(address, 2, 0));
        n[3] = getNode(node->depth(), shiftAddr(address, 2, 1));
        n[4] = getNode(node->depth(), shiftAddr(address, 4, 0));
        n[5] = getNode(node->depth(), shiftAddr(address, 4, 1));
    }
    else
    {
        for (int i=0; i<8; i++)
        {
            setNeighbours(node->children()[i], (address<<3)+i);
        }
    }
}


float OctTree::raytrace(vec3f& orig, vec3f& dir, Octnode* camnode, triangle*& hittri)
{
    float hit_time = MAX_FLOAT;


    // calculate the time it takes to cross an oct node at each level
    for (uint64_t i=0; i<=max_depth; ++i)
    {
        for (int64_t j=0; j<3; ++j)
        {
            ray_times[i][j] = fabs(dists[i][j]/dir[j]);
        }
//        std::cout << "ray_times "<< i << " " << ray_times[i][0] << " " << ray_times[i][1] << " " << ray_times[i][2] << std::endl;
    }

    // set up initial times to next node
    int64_t move_idx[3];
    float node_times[3];
    for (int64_t i=0; i<3; ++i)
    {
        if(dir[i] > 0)
        {
            move_idx[i] = (i*2)+1;
            node_times[i] = (camnode->max()[i] - orig[i]) / dir[i];
        }
        else
        {
            move_idx[i] = (i*2);
            node_times[i] = (camnode->min()[i] - orig[i]) / dir[i];
        }
    }

//    std::cout << "============================================================ " << ::std::endl;
//    std::cout << "orig: " << orig[0] << " " << orig[1] << " " << orig[2] << std::endl;
//    std::cout << "dir: " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;

    Octnode* node = camnode;
    float total_time = 0;
    while(node != NULL)
    {
//        std::cout << "Node time: " << node_times[0] << " " << node_times[1] << " " << node_times[2] << std::endl;

        // move to next node
        int64_t min_index;
        if(node_times[0] < node_times[1]) {
            if(node_times[0] < node_times[2]) {
                min_index = 0;
            } else {
                min_index = 2;
            }
        } else {
            if(node_times[1] < node_times[2]) {
                min_index = 1;
            } else {
                min_index = 2;
            }
        }

        // record how long the ray has been travel through the tree
        total_time += node_times[min_index];

        // check triangles in camera node
        std::vector<triangle>* tris = node->tris();
        std::vector<triangle>::iterator iter;
        for(iter = tris->begin(); iter != tris->end(); ++iter)
        {
            float time = iter->intersect(orig, dir);
            if (time > 0 && time < hit_time)// && time < total_time)
            {
                hit_time = time;
                hittri = &(*iter);
            }
        }

        if (hit_time < MAX_FLOAT) {
            return hit_time;
        }

        float min_time = node_times[min_index];
        for (int64_t i=0; i<3; ++i)
        {
            node_times[i] -= min_time;
        }


        Octnode* newnode = node->neighbours()[move_idx[min_index]];

//        std::string adr;
//        uint64_t shifted = shiftAddr(node->address(), pow(2,min_index), dir[min_index]>0);
//        Octnode* tmpNode = getNode(node->depth(), shifted);
//        if (tmpNode != NULL)
//        {
//            printAddress(tmpNode->depth(), tmpNode->address(), &adr);
//            std::cout << "shifted " << adr << " " << *tmpNode <<std::endl;
//        }
//        adr.clear();
//        printAddress(node->depth(), shifted, &adr);
//        std::cout << "total_time " << total_time << " " << min_index << " " << move_idx[min_index] << " " << adr << std::endl;
//        adr.clear();
//        printAddress (node->depth(), node->address(), &adr);
//        std::cout << "node " << adr << " " << *node << std::endl;

        if (newnode != NULL)
        {
//            std::string adr;
//            printAddress (newnode->depth(), newnode->address(),&adr);
//            std::cout << "next " << adr << " " <<  *newnode << std::endl;
            if (newnode->depth() < node->depth())
            {
                // going up levels
                // we went up a level in the tree so need to know which node we left so we can adjust the time
                uint64_t tmp_addr = node->address();
                for(int64_t i=newnode->depth(); i<node->depth(); ++i)
                {
                    for(int64_t j=0; j<3; ++j)
                    {
                        if(dir[j] > 0)
                        {
                            if((tmp_addr & masks[j]) == 0)
                            {
                                node_times[i] += ray_times[i+1][j];
                            }
                        }
                        else
                        {
                            if((tmp_addr & masks[j]) != 0)
                            {
                                node_times[i] += ray_times[i+1][j];
                            }
                        }
                    }
                    tmp_addr >>= 3;
                }
            }
            else
            {
                if (!newnode->isLeaf())
                {
                    // going down in levels
                    // for(int64_t i=newnode->depth(); i<node->depth(); ++i)
                    while(!newnode->isLeaf())
                    {
                        // we went down a level in the tree
                        int64_t oct = 0;
                        for (int64_t i=0; i<3; ++i)
                        {
                            if(dir[i] > 0)
                            {
                                if(node_times[i] > ray_times[node->depth()+1][i])
                                {
                                    node_times[i] -= ray_times[node->depth()+1][i];
                                    oct += pow(2, i);
                                }
                                else { }
                            }
                            else
                            {
                                if(node_times[i] > ray_times[node->depth()+1][i])
                                {
                                    node_times[i] -= ray_times[node->depth()+1][i];
                                }
                                else
                                {
                                    oct += pow(2, i);
                                }
                            }
                        }
                        newnode = newnode->children()[oct];
                    }
                }
                else
                {
                    // the nodes are the same level so we are done
                    node_times[min_index] += ray_times[newnode->depth()][min_index];
                }


            }
        }
        node = newnode;
    }

//    abort();

    return MAX_FLOAT;
}


} // namespace lw
