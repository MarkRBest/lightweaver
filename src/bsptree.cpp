#include "bsptree.hpp"


lw::BSPTree::BSPTree() : num_nodes(0), num_leaves(0) {}

void lw::BSPTree::build(std::vector<trianglePtr>& _tris)
{
    for(uint64_t i=0; i<2*_tris.size(); ++i)
    {
        tris.push_back(std::vector<trianglePtr>());
        splits.push_back(0.0);
    }

    // recursively build the trees
    pbuild(0, 0, 0, _tris);
}


float lw::BSPTree::raytrace(vec3f& orig, vec3f& dir, int64_t plane, int64_t address, trianglePtr& hittri)
{
    float split = splits[address];

    // std::cout << "plane " << plane << " " << address << " " << split << std::endl;

    if(split == -666)
    {
        float min_hit_time = std::numeric_limits<float>::max();
        std::vector<trianglePtr>& tmp = tris[address];
        for(auto tri : tmp)
        {
            float hit_time = tri->intersect(orig, dir);
            if(hit_time > 0 && hit_time < min_hit_time)
            {
                hittri = tri;
                min_hit_time = hit_time;
            }

        }

        return min_hit_time;
    }

    int64_t order[2];
    if(orig[plane] < split)
    {
        order[0] = 1;
        order[1] = 2;
    } else {
        order[0] = 2;
        order[1] = 1;
    }

    float plane_hit_time = (split - orig[plane]) / dir[plane];
    float hit_time = raytrace(orig, dir, (plane+1)%3, (address*2)+order[0], hittri);
    if (hit_time != std::numeric_limits<float>::max() && hit_time < plane_hit_time)
    {
        return hit_time;
    }

    // check if the ray goes through the
    if (plane_hit_time > 0)
    {
        hit_time = raytrace(orig, dir, (plane+1)%3, (address*2)+order[1], hittri);
    }

    return hit_time;
}


void lw::BSPTree::pbuild(int64_t plane, int64_t depth, int64_t address, std::vector<trianglePtr>& _tris)
    {
        num_nodes++;
//        std::cout << "plane " << plane << " depth " << depth << " addr " << address << " tris " << _tris.size() << " 2x " << 2*_tris.size() << std::endl;
        if(_tris.size() == 1)
        {
            splits[address] = -666;
            std::vector<trianglePtr>& tmp = tris[address];
            tmp.insert(tmp.begin(), _tris.begin(), _tris.end());
        }
        else
        {
            workspace.clear();
            for(auto & tri : _tris)
            {
                workspace.push_back(tri->getVert(0)[plane]);
                workspace.push_back(tri->getVert(1)[plane]);
                workspace.push_back(tri->getVert(2)[plane]);
            }
            std::sort(workspace.begin(), workspace.end());

            int64_t n = workspace.size()/2;

            std::vector<trianglePtr> best_left, best_right;
            float best_split;
            int32_t best_delta = std::numeric_limits<int32_t>::max();
            std::vector<trianglePtr> left, right;
            int jitter = 0, sign = +1;
            float proposal;
            int cnt=0;
            do {
                cnt++;

                proposal = workspace[n + jitter*sign];
                if (sign == -1) {
                    jitter++;
                }
                sign *= -1;

                left.clear();
                right.clear();

                for (auto triptr: _tris)
                {
                    float min_plane = std::min(std::min(triptr->getVert(0)[plane], triptr->getVert(1)[plane]),triptr->getVert(2)[plane]);
                    float max_plane = std::max(std::max(triptr->getVert(0)[plane], triptr->getVert(1)[plane]),triptr->getVert(2)[plane]);

                    if (max_plane <= proposal)
                    {
                        left.push_back(triptr);
                    }
                    else if (min_plane >= proposal)
                    {
                        right.push_back(triptr);
                    }
                    else
                    {
                        left.push_back(triptr);
                        right.push_back(triptr);
                    }
                }

                int delta = abs(left.size() - right.size());
                if (best_delta > delta)
                {
                    best_delta = delta;
                    best_split = proposal;

                    best_left.clear();
                    best_left.insert(best_left.begin(), left.begin(), left.end());

                    best_right.clear();
                    best_right.insert(best_right.begin(), right.begin(), right.end());
                }

            } while(jitter < _tris.size() * 0.2 && abs(left.size() - right.size()) > 1);


            if (
                _tris.size()*1.5 < (best_left.size()+best_right.size()) ||
                std::min(best_left.size(), best_right.size()) == 0  ||
                std::max(best_left.size(), best_right.size()) <= 1000000 )
            {
                // if we couldn't sub divided the lists then we are at a node
                splits[address] = -666;
                std::vector<trianglePtr>& tmp = tris[address];
                tmp.insert(tmp.begin(), _tris.begin(), _tris.end());

//                std::cout << "Leaf " << depth << " " << address <<" " << best_split << " " << _tris.size() << " " << best_left.size() << " " << best_right.size() << " " << best_left.size() + best_right.size() << std::endl;
                num_leaves++;
            }
            else
            {
//                std::cout << "Split " << depth << " " << address <<" " << best_split << " " << _tris.size() << " " << best_left.size() << " " << best_right.size() << " " << best_left.size() + best_right.size() << std::endl;
                splits[address] = best_split;

                // split the set of triangles
                pbuild((plane+1)%3, depth+1, address*2+1, best_left);
                pbuild((plane+1)%3, depth+1, address*2+2, best_right);
            }
        }
    }
