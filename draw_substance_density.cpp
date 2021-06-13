#include "draw_substance_density.h"

#include "physical_const.h"
#include "utils.h"
#include <iostream>

namespace smoke_simulation{
//物質密度を描画する関数
void draw_substance_density(Grid& all_grid, int scale, cv::VideoWriter& writer){
    cv::Mat src(scale*(smoke_simulation::physical_const::kGrid_num_y), scale*(smoke_simulation::physical_const::kGrid_num_x), CV_8UC1);
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            for(int j=0;j<scale;j++){
                for(int k=0;k<scale;k++){
                    if(all_grid.substance_density[get_voxel_center_index(ix,iy)]>0.0){
                        if(all_grid.substance_density[get_voxel_center_index(ix,iy)]>1.0){
                            src.at<unsigned char>(scale*(smoke_simulation::physical_const::kGrid_num_y-1-iy)+j, scale*(ix)+k)=(unsigned char)255;
                        }
                        else{
                            src.at<unsigned char>(scale*(smoke_simulation::physical_const::kGrid_num_y-1-iy)+j, scale*(ix)+k)=(unsigned char)(all_grid.substance_density[get_voxel_center_index(ix,iy)]*255.9);
                        }
                    }
                    else{
                        src.at<unsigned char>(scale*(smoke_simulation::physical_const::kGrid_num_y-1-iy)+j, scale*(ix)+k)=(unsigned char)0;
                    }
                }
            }
        }
    }
    cv::imshow(" ",   src);
    writer << src;
    cv::waitKey(1);
}
}// namespace smoke_simulation
