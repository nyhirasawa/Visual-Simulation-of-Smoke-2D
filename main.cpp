//clang++ main.cpp grid.cpp initialize_grid.cpp utils.cpp draw_substance_density.cpp update_fluid_velocity.cpp linear_solver.cpp sparse_matrix.cpp move_substances.cpp `pkg-config --cflags --libs opencv4` -std=c++11 -O3
#include <opencv2/opencv.hpp>

#include <iostream>
#include <vector>

#include "grid.h"
#include "initialize_grid.h"
#include "physical_const.h"
#include "draw_substance_density.h"
#include "update_fluid_velocity.h"
#include "move_substances.h"
#include "utils.h"

int main(int argc, char* argv[]){
    //系の情報が乗ったグリッド
    smoke_simulation::Grid all_grid(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y);

    //グリッドを初期化
    smoke_simulation::initialize_grid(all_grid);

    //出力する動画のフォーマットを設定
    //フレームレート
    double fps=1.0/smoke_simulation::physical_const::kDt;
    //出力形式('M', 'P', '4', 'V' ならmp4形式で出力する)
    int fourcc=cv::VideoWriter::fourcc('M', 'P', '4', 'V');
    const int scale=4;
    int width=smoke_simulation::physical_const::kGrid_num_x*scale;
    int height=smoke_simulation::physical_const::kGrid_num_y*scale;
    cv::VideoWriter writer;
    writer.open("result.mp4", fourcc, fps, cv::Size(width, height), false);


    //メインの計算部分
    const int num_frame=300;
    for(int i=0;i<num_frame;i++){
        std::cout<<i<<"/"<<num_frame<<std::endl;
        //substance densityの描画
        smoke_simulation::draw_substance_density(all_grid, scale, writer);
        //流体の速度場の更新
        smoke_simulation::update_fluid_velocity(all_grid);
        //substance densityの更新
        smoke_simulation::move_substances(all_grid);
    }

    return 0;
}
