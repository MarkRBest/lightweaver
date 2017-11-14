#pragma once

#include <iostream>
#include <stdint.h>
#include <vector>

namespace lw {

template<class X>
class kdNode;

template<class X>
class kdNode {
private:
    int32_t m_plane;
    float m_split;
    X* m_data;
    kdNode *m_l, *m_r;

public:

    kdNode(int32_t plane) :
            m_plane(plane), m_split(0.0), m_data(NULL), m_l(NULL), m_r(NULL) {
    }

    inline bool isLeaf(){
        return m_data != NULL;
    }

    inline int32_t plane() {
        return m_plane;
    }

    inline void setPlane(int32_t plane) {
        m_plane = plane;
    }

    inline void setSplit(float v) {
        m_split = v;
    }

    inline float split() {
        return m_split;
    }

    inline void setLeft(kdNode* left) {
        m_l = left;
    }

    inline kdNode* left() {
        return m_l;
    }

    inline void setRight(kdNode* right) {
        m_r = right;
    }

    inline kdNode* right() {
        return m_r;
    }

    inline void setData(X* x) {
        m_data = x;
    }

    inline X* data() {
        return m_data;
    }
};

} // namespace lw
