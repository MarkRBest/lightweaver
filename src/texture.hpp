#pragma once

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include <jconfig.h>
#include <jpeglib.h>
#include <setjmp.h>

#include <vector>

namespace lw {

class Texture {
private:
    unsigned char*** data;
    int row_count;
    int w, h;

    int read_JPEG_file(const char * filename);
    void store_data(JSAMPROW pixelrow, int row_width);

public:
    Texture(const char* filename) :
        row_count(0), w(0), h(0) {
        read_JPEG_file(filename);
    }
    inline int height() { return h; }
    inline int width() { return w; }
    inline unsigned char* pixel(int h, int w){
        return data[h][w];
    }
};

} // namespace lw
