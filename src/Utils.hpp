#pragma once

#include <limits>

namespace lw {

#define MAX_FLOAT std::numeric_limits<float>::max()
#define MIN_FLOAT std::numeric_limits<float>::min()

#define innerProduct(v,q) \
        ((v)[0] * (q)[0] + \
        (v)[1] * (q)[1] + \
        (v)[2] * (q)[2])

#define crossProduct(a,b,c) \
    (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
    (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
    (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

// generate random float within [0,1]
#define randf() ((float) rand() / RAND_MAX)

}
