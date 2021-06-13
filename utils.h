#ifndef UTILS_H
#define UTILS_H

#include <iostream>

#include "physical_const.h"

namespace smoke_simulation{
//グリッドを初期化する関数
inline int get_voxel_face_index_x(int ix, int iy){
    if(ix<0){
        std::cout<<"範囲外参照x0"<<std::endl;
    }
    if(ix>smoke_simulation::physical_const::kGrid_num_x){
        std::cout<<"範囲外参照x1"<<std::endl;
    }
    if(iy<0){
        std::cout<<"範囲外参照x2"<<std::endl;
    }
    if(iy>smoke_simulation::physical_const::kGrid_num_y-1){
        std::cout<<"範囲外参照x3"<<std::endl;
    }
    return smoke_simulation::physical_const::kGrid_num_y*ix+iy;
}
inline int get_voxel_face_index_y(int ix, int iy){
    if(ix<0){
        std::cout<<"範囲外参照y0"<<std::endl;
    }
    if(ix>smoke_simulation::physical_const::kGrid_num_x-1){
        std::cout<<"範囲外参照y1"<<std::endl;
    }
    if(iy<0){
        std::cout<<"範囲外参照y2"<<std::endl;
    }
    if(iy>smoke_simulation::physical_const::kGrid_num_y){
        std::cout<<"範囲外参照y3"<<std::endl;
    }
    return (smoke_simulation::physical_const::kGrid_num_y+1)*ix+iy;
}
inline int get_voxel_center_index(int ix, int iy){
    if(ix<0){
        std::cout<<"範囲外参照c0"<<std::endl;
    }
    if(ix>smoke_simulation::physical_const::kGrid_num_x-1){
        std::cout<<"範囲外参照c1"<<std::endl;
    }
    if(iy<0){
        std::cout<<"範囲外参照c2"<<std::endl;
    }
    if(iy>smoke_simulation::physical_const::kGrid_num_y-1){
        std::cout<<"範囲外参照c3"<<std::endl;
    }
    return smoke_simulation::physical_const::kGrid_num_y*ix+iy;
}

}//namespace smoke_simulation

#endif //INITIALIZE_GRID_H
