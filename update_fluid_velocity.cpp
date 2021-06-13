#include "update_fluid_velocity.h"

#include <math.h>
#include <iostream>
#include <chrono>//時間計測用
#include <fstream>//ファイル書き出し用

#include "linear_solver.h"
#include "physical_const.h"
#include "sparse_matrix.h"
#include "utils.h"

namespace smoke_simulation{

void add_source(Grid& all_grid){
/*
    for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;++iy){
        if(iy>(7*smoke_simulation::physical_const::kGrid_num_y/15)
         &&iy<(8*smoke_simulation::physical_const::kGrid_num_y/15)){
            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0, iy)]=10.0;
            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(0, iy)]=0.0;
            all_grid.substance_density[get_voxel_center_index(0, iy)]=5.0;
//            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(1, iy)]=10.0;
//            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(1, iy)]=0.0;
//            all_grid.substance_density[get_voxel_center_index(1, iy)]=5.0;
        }
    }
*/
}

//壁における速度場をセットする関数
void set_boundary_velocity(Grid& all_grid){
    //y軸に垂直な壁の速度をセット
    for(int i=1;i<smoke_simulation::physical_const::kGrid_num_x-1;i++){
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(i,0)]=-all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(i,1)];
        all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(i,smoke_simulation::physical_const::kGrid_num_y)]
            =-all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(i,smoke_simulation::physical_const::kGrid_num_y-1)];
    }
    //x軸に垂直な壁の速度をセット
    for(int i=1;i<smoke_simulation::physical_const::kGrid_num_y-1;i++){
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0,i)]=-all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(1,i)];
        all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x,i)]
            =-all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x-1,i)];
    }
    //四隅の速度場は周辺から線形補間
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0,0)]=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(1,0)]+all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0,1)])/2.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x,0)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x-1,0)]
        +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x,1)])/2.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0,smoke_simulation::physical_const::kGrid_num_y-1)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(0,smoke_simulation::physical_const::kGrid_num_y-2)]
        +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(1, smoke_simulation::physical_const::kGrid_num_y-1)])/2.0;
    all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-1)]
        =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1)]
        +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(smoke_simulation::physical_const::kGrid_num_x, smoke_simulation::physical_const::kGrid_num_y-2)])/2.0;

    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(0,0)]=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(1,0)]+all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(0,1)])/2.0;
    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(smoke_simulation::physical_const::kGrid_num_x-1,0)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(smoke_simulation::physical_const::kGrid_num_x-2,0)]
        +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(smoke_simulation::physical_const::kGrid_num_x-1,1)])/2.0;
    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(0,smoke_simulation::physical_const::kGrid_num_y)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(0,smoke_simulation::physical_const::kGrid_num_y-1)]
        +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(1,smoke_simulation::physical_const::kGrid_num_y)])/2.0;
    all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(smoke_simulation::physical_const::kGrid_num_x-1,smoke_simulation::physical_const::kGrid_num_y)]
        =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(smoke_simulation::physical_const::kGrid_num_x-2,smoke_simulation::physical_const::kGrid_num_y)]
        +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(smoke_simulation::physical_const::kGrid_num_x-1,smoke_simulation::physical_const::kGrid_num_y-1)])/2.0;
}

//壁における圧力場をセットする関数
void set_boundary_pressure(Grid& all_grid){
    //y軸に垂直な壁の圧力をセット
    for(int i=1;i<smoke_simulation::physical_const::kGrid_num_x-1;i++){
        all_grid.pressure[get_voxel_center_index(i,0)]
            =all_grid.pressure[get_voxel_center_index(i,1)];
        all_grid.pressure[get_voxel_center_index(i,smoke_simulation::physical_const::kGrid_num_y-1)]
            =all_grid.pressure[get_voxel_center_index(i,smoke_simulation::physical_const::kGrid_num_y-2)];
    }
    //x軸に垂直な壁の圧力をセット
    for(int i=1;i<smoke_simulation::physical_const::kGrid_num_y-1;i++){
        all_grid.pressure[get_voxel_center_index(0,i)]
            =all_grid.pressure[get_voxel_center_index(1,i)];
        all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-1,i)]
            =all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-2,i)];
    }
    //四隅の圧力場は周辺から線形補間
    all_grid.pressure[get_voxel_center_index(0,0)]
        =(all_grid.pressure[get_voxel_center_index(1,0)]
        +all_grid.pressure[get_voxel_center_index(0,1)])/2.0;
    all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-1,0)]
        =(all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-2,0)]
        +all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-1,1)])/2.0;
    all_grid.pressure[get_voxel_center_index(0,smoke_simulation::physical_const::kGrid_num_y-1)]
        =(all_grid.pressure[get_voxel_center_index(0,smoke_simulation::physical_const::kGrid_num_y-2)]
        +all_grid.pressure[get_voxel_center_index(1, smoke_simulation::physical_const::kGrid_num_y-1)])/2.0;
    all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-1)]
        =(all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-2, smoke_simulation::physical_const::kGrid_num_y-1)]
        +all_grid.pressure[get_voxel_center_index(smoke_simulation::physical_const::kGrid_num_x-1, smoke_simulation::physical_const::kGrid_num_y-2)])/2.0;
}


//(vorticityの計算に使う)cell centered velocityの計算
void calc_cell_centered_velocity(Grid& all_grid){
    //cell_centered_velocityの計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy)][0]
                =(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]+all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix+1,iy)])/2.0;
            all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy)][1]
                =(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]+all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy+1)])/2.0;
            all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy)][2]=0.0;
        }
    }
}

//vorticity confinement termの計算
void calc_vorticity_confinement(Grid& all_grid){
    //vorticityの計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            //2次元バージョン
            all_grid.vorticity[get_voxel_center_index(ix,iy)][0]=0.0;
            all_grid.vorticity[get_voxel_center_index(ix,iy)][1]=0.0;
            //if~else if はgrid cellの範囲外を参照しないようにする処方
            if(ix+1>=smoke_simulation::physical_const::kGrid_num_x && iy+1>=smoke_simulation::physical_const::kGrid_num_y){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =(-all_grid.velocity_in_cell_center[get_voxel_center_index(ix-1,iy)][1]
                      +all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy-1)][0])
                      /smoke_simulation::physical_const::kCell_length;
            }
            else if(ix+1>=smoke_simulation::physical_const::kGrid_num_x && iy-1<0){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =(-all_grid.velocity_in_cell_center[get_voxel_center_index(ix-1,iy)][1]
                      -all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy+1)][0])
                      /smoke_simulation::physical_const::kCell_length;
            }
            else if(ix+1>=smoke_simulation::physical_const::kGrid_num_x){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =(-all_grid.velocity_in_cell_center[get_voxel_center_index(ix-1,iy)][1]
                      /smoke_simulation::physical_const::kCell_length)
                    +((-all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy+1)][0]
                       +all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy-1)][0])
                       /(2.0*smoke_simulation::physical_const::kCell_length));
            }
            else if(ix-1<0 && iy+1>=smoke_simulation::physical_const::kGrid_num_y){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =(all_grid.velocity_in_cell_center[get_voxel_center_index(ix+1,iy)][1]
                     +all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy-1)][0])
                     /smoke_simulation::physical_const::kCell_length;
            }
            else if(ix-1<0 && iy-1<0){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =(all_grid.velocity_in_cell_center[get_voxel_center_index(ix+1,iy)][1]
                     -all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy+1)][0])
                     /smoke_simulation::physical_const::kCell_length;
            }
            else if(ix-1<0){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =(all_grid.velocity_in_cell_center[get_voxel_center_index(ix+1,iy)][1]
                     /smoke_simulation::physical_const::kCell_length)
                    +(-all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy+1)][0]
                      +all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy-1)][0])
                     /(2.0*smoke_simulation::physical_const::kCell_length);
            }
            else if(iy+1>=smoke_simulation::physical_const::kGrid_num_y){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =((all_grid.velocity_in_cell_center[get_voxel_center_index(ix+1,iy)][1]
                      -all_grid.velocity_in_cell_center[get_voxel_center_index(ix-1,iy)][1])
                      /(2.0*smoke_simulation::physical_const::kCell_length))
                     +(all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy-1)][0]
                      /smoke_simulation::physical_const::kCell_length);
            }
            else if(iy-1<0){
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =((all_grid.velocity_in_cell_center[get_voxel_center_index(ix+1,iy)][1]
                      -all_grid.velocity_in_cell_center[get_voxel_center_index(ix-1,iy)][1])
                      /(2.0*smoke_simulation::physical_const::kCell_length))
                     +(-all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy+1)][0]
                      /smoke_simulation::physical_const::kCell_length);
            }
            else{
                all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                    =(all_grid.velocity_in_cell_center[get_voxel_center_index(ix+1,iy)][1]
                     -all_grid.velocity_in_cell_center[get_voxel_center_index(ix-1,iy)][1]
                     -all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy+1)][0]
                     +all_grid.velocity_in_cell_center[get_voxel_center_index(ix,iy-1)][0])
                     /(2.0*smoke_simulation::physical_const::kCell_length);
            }
        }
    }
    //ノルムを1に正規化したvorticity amplitude gradient(論文のN)の計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            //if~else if はgrid cellの範囲外を参照しないようにする処方
            if(ix+1>=smoke_simulation::physical_const::kGrid_num_x && iy+1>=smoke_simulation::physical_const::kGrid_num_y){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else if(ix+1>=smoke_simulation::physical_const::kGrid_num_x && iy-1<0){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else if(ix+1>=smoke_simulation::physical_const::kGrid_num_x){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/(2.0*smoke_simulation::physical_const::kCell_length);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else if(ix-1<0 && iy+1>=smoke_simulation::physical_const::kGrid_num_y){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else if(ix-1<0 && iy-1<0){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else if(ix-1<0){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/(2.0*smoke_simulation::physical_const::kCell_length);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else if(iy+1>=smoke_simulation::physical_const::kGrid_num_y){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/(2.0*smoke_simulation::physical_const::kCell_length);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else if(iy-1<0){
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/(2.0*smoke_simulation::physical_const::kCell_length);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/smoke_simulation::physical_const::kCell_length;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            else{
                double vorticity_amplitude_x0=sqrt(all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix-1,iy)][2]);
                double vorticity_amplitude_x1=sqrt(all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]*all_grid.vorticity[get_voxel_center_index(ix+1,iy)][2]);
                double vorticity_amplitude_y0=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy-1)][2]);
                double vorticity_amplitude_y1=sqrt(all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][0]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][1]
                                                  +all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]*all_grid.vorticity[get_voxel_center_index(ix,iy+1)][2]);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]=(vorticity_amplitude_x1-vorticity_amplitude_x0)/(2.0*smoke_simulation::physical_const::kCell_length);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]=(vorticity_amplitude_y1-vorticity_amplitude_y0)/(2.0*smoke_simulation::physical_const::kCell_length);
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]=0.0;
            }
            double normalize_factor=sqrt(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]*all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]
                                        +all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]*all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]
                                        +all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]*all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]);
            //0除算を回避するための処理
            if(normalize_factor>0.0001){
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][0]/=normalize_factor;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][1]/=normalize_factor;
                all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix, iy)][2]/=normalize_factor;
            }
        }
    }
    //vorticity confinment term の計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][0]
                +=smoke_simulation::physical_const::kConfinement_amplitude
                 *smoke_simulation::physical_const::kCell_length
                 *(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix,iy)][1]
                  *all_grid.vorticity[get_voxel_center_index(ix,iy)][2]
                  -all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix,iy)][2]
                  *all_grid.vorticity[get_voxel_center_index(ix,iy)][1]);
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][1]
                +=smoke_simulation::physical_const::kConfinement_amplitude
                 *smoke_simulation::physical_const::kCell_length
                 *(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix,iy)][2]
                  *all_grid.vorticity[get_voxel_center_index(ix,iy)][0]
                  -all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix,iy)][0]
                  *all_grid.vorticity[get_voxel_center_index(ix,iy)][2]);
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][2]
                +=smoke_simulation::physical_const::kConfinement_amplitude
                 *smoke_simulation::physical_const::kCell_length
                 *(all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix,iy)][0]
                  *all_grid.vorticity[get_voxel_center_index(ix,iy)][1]
                  -all_grid.normalized_vorticity_amplitude_gradient[get_voxel_center_index(ix,iy)][1]
                  *all_grid.vorticity[get_voxel_center_index(ix,iy)][0]);
        }
    }
}

//外力場による速度場の更新
void update_fluid_velocity_by_external_force(Grid& all_grid){
    //計算した外力場により速度場を更新
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x+1;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            if(ix-1<0){
                all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]
                    +=all_grid.external_force_field[get_voxel_center_index(ix,iy)][0]
                     *smoke_simulation::physical_const::kDt;
            }
            else if(ix>=smoke_simulation::physical_const::kGrid_num_x){
                all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]
                    +=all_grid.external_force_field[get_voxel_center_index(ix-1,iy)][0]
                     *smoke_simulation::physical_const::kDt;
            }
            else{
                all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]
                    +=((all_grid.external_force_field[get_voxel_center_index(ix-1,iy)][0]
                     +all_grid.external_force_field[get_voxel_center_index(ix,iy)][0])/2.0)
                     *smoke_simulation::physical_const::kDt;
            }
        }
    }
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y+1;iy++){
            if(iy-1<0){
                all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]
                    +=all_grid.external_force_field[get_voxel_center_index(ix,iy)][1]
                     *smoke_simulation::physical_const::kDt;
            }
            else if(iy>=smoke_simulation::physical_const::kGrid_num_y){
                all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]
                    +=all_grid.external_force_field[get_voxel_center_index(ix,iy-1)][1]
                     *smoke_simulation::physical_const::kDt;
            }
            else{
                all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]
                    +=((all_grid.external_force_field[get_voxel_center_index(ix,iy-1)][1]
                       +all_grid.external_force_field[get_voxel_center_index(ix,iy)][1])/2.0)
                      *smoke_simulation::physical_const::kDt;
            }
        }
    }
}


//外力項の計算
void add_force_fluid(Grid& all_grid){
    //外力場をリセット
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][0]=0.0;
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][1]=0.0;
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][2]=0.0;
        }
    }
    calc_cell_centered_velocity(all_grid);
/*
    //重力の効果
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.external_force_field[get_voxel_center_index(ix,iy)][1]-=9.8*all_grid.substance_density[get_voxel_center_index(ix, iy)]*smoke_simulation::physical_const::kDt;
        }
    }
*/
    calc_vorticity_confinement(all_grid);
    update_fluid_velocity_by_external_force(all_grid);
}

//x成分のvelocityのlinear interpolation
double linear_interpolation_x(double advected_x, double advected_y, Grid& all_grid){
    //バックトレース先の座標のindex
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(double)(advected_index_y);
    }
    //バックトレース先の速度を線形補間する
    double a0, a1;
    double b0, b1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
//    return all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x,advected_index_y)];

//    if(advected_index_y+1>smoke_simulation::physical_const::kGrid_num_y-1
//     &&advected_index_x+1>smoke_simulation::physical_const::kGrid_num_x){
//         return all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x,advected_index_y)];
//     }
//    else if(advected_index_y+1>smoke_simulation::physical_const::kGrid_num_y-1){
//        return a1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x,advected_index_y)]
//              +a0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x+1,advected_index_y)];
//    }
//    else if(advected_index_x+1>smoke_simulation::physical_const::kGrid_num_x){
//        return b1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x,advected_index_y)]
//              +b0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x,advected_index_y+1)];
//    }
//    else{
        return a1*(b1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x,advected_index_y)]
              +b0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x,advected_index_y+1)])
              +a0*(b1*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x+1,advected_index_y)]
              +b0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x+1,advected_index_y+1)]);
//    }
}

//y成分のvelocityのlinear interpolation
double linear_interpolation_y(double advected_x, double advected_y, Grid& all_grid){
    //バックトレース先の座標のindex
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(double)(advected_index_x);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(double)(advected_index_x);
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(double)advected_index_y;
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-1;
        advected_y=(double)advected_index_y;
    }
    //バックトレース先の速度を線形補間する
    double a0, a1;
    double b0, b1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
//    return all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x,advected_index_y)];

//    if(advected_index_y+1>smoke_simulation::physical_const::kGrid_num_y
//     &&advected_index_x+1>smoke_simulation::physical_const::kGrid_num_x-1){
//         return all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x,advected_index_y)];
//     }
//    else if(advected_index_y+1>smoke_simulation::physical_const::kGrid_num_y){
//        return a1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x,advected_index_y)]
//              +a0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x+1,advected_index_y)];
//    }
//    else if(advected_index_x+1>smoke_simulation::physical_const::kGrid_num_x-1){
//        return b1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x,advected_index_y)]
//              +b0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x,advected_index_y+1)];
//    }
//    else{
        return a1*(b1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x,advected_index_y)]
              +b0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x,advected_index_y+1)])
              +a0*(b1*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x+1,advected_index_y)]
              +b0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x+1,advected_index_y+1)]);
//    }
}

//x成分のvelocityのmonotonic cubic interpolation
double monotonic_cubic_interpolation_x(double advected_x, double advected_y, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<0){
        advected_index_x=0;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x+1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_y<0){
        advected_index_y=0;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-1;
        advected_y=(double)(advected_index_y);
    }
    //範囲外を参照しないようにする処置
    if(advected_index_x<=1||advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1
     ||advected_index_y<=1||advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-2){
         return linear_interpolation_x(advected_x, advected_y, all_grid);
    }
    //monotonic cubic interpolation の処理
    else{
        double dk_0;
        double dk_1;
        double delta_k;
        double interpolated_velocity_x[4];
        //y軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            dk_0  =  (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y+1)]
                     -all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y-1)])/2.0;
            dk_1  =  (all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y+2)]
                     -all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)])/2.0;
            delta_k = all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y+1)]
                     -all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)];
            double b0 = advected_y-advected_index_y;
            //3つのスロープの符号が異なる場合
            if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                interpolated_velocity_x[i]=(1.0-b0)*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i,advected_index_y)]
                                          +b0*all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i,advected_index_y+1)];
//                interpolated_velocity_x[i]=all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)];
            }

            else{
                //論文の式は間違っている。正しくはこっち
                interpolated_velocity_x[i]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                          +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)];
            }

        }
        //x軸方向に関するinterpolation
        dk_0=(interpolated_velocity_x[2]-interpolated_velocity_x[0])/2.0;
        dk_1=(interpolated_velocity_x[3]-interpolated_velocity_x[1])/2.0;
        delta_k=(interpolated_velocity_x[2]-interpolated_velocity_x[1]);
        double a0=advected_x-advected_index_x;
        //3つのスロープの符号が異なる場合
        if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
            return (1.0-a0)*interpolated_velocity_x[1]+a0*interpolated_velocity_x[2];
        }
        else{
            //論文の式は間違っている。正しくはこっち
            return (dk_0+dk_1-2*delta_k)*(a0*a0*a0)+(3*delta_k-2*dk_0-dk_1)*(a0*a0)+dk_0*a0
                    +interpolated_velocity_x[1];
        }
    }
}

//y成分のvelocityのmonotonic cubic interpolation
double monotonic_cubic_interpolation_y(double advected_x, double advected_y, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    //バックトレース先の座標のindexが系の外に出てしまった場合の処理
    if(advected_index_x<0){
        advected_index_x=0;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-1;
        advected_x=(double)advected_index_x;
    }
    if(advected_index_y<0){
        advected_index_y=0;
        advected_y=(double)(advected_index_y);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y+1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y;
        advected_y=(double)(advected_index_y);

    }
    //範囲外を参照しないようにする処置
    if(advected_index_x<=1||advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-2
     ||advected_index_y<=1||advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
         return linear_interpolation_y(advected_x, advected_y, all_grid);
    }
    //monotonic cubic interpolation の処理
    else{
        double dk_0;
        double dk_1;
        double delta_k;
        double interpolated_velocity_y[4];
        //y軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            dk_0  =  (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i, advected_index_y+1)]
                     -all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i, advected_index_y-1)])/2.0;
            dk_1  =  (all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i, advected_index_y+2)]
                     -all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i, advected_index_y)])/2.0;
            delta_k = all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i, advected_index_y+1)]
                     -all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i, advected_index_y)];
            double b0 = advected_y-advected_index_y;
            //3つのスロープの符号が異なる場合
            if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                interpolated_velocity_y[i]=(1.0-b0)*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i,advected_index_y)]
                                          +b0*all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i,advected_index_y+1)];
//                interpolated_velocity_x[i]=all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(advected_index_x-1+i, advected_index_y)];
            }
            else{
                //論文の式は間違っている。正しくはこっち
                interpolated_velocity_y[i]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                          +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(advected_index_x-1+i, advected_index_y)];
            }
        }
        //x軸方向に関するinterpolation
        dk_0=(interpolated_velocity_y[2]-interpolated_velocity_y[0])/2.0;
        dk_1=(interpolated_velocity_y[3]-interpolated_velocity_y[1])/2.0;
        delta_k=(interpolated_velocity_y[2]-interpolated_velocity_y[1]);
        double a0=advected_x-advected_index_x;
        //3つのスロープの符号が異なる場合
        if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
            return (1.0-a0)*interpolated_velocity_y[1]+a0*interpolated_velocity_y[2];
        }
        else{
            //論文の式は間違っている。正しくはこっち
            return (dk_0+dk_1-2*delta_k)*(a0*a0*a0)+(3*delta_k-2*dk_0-dk_1)*(a0*a0)+dk_0*a0
                    +interpolated_velocity_y[1];
        }
    }
}


//advect項の計算
//速度場を時間 -dt だけバックトレースしてadvect項を計算する
void advect_fluid(Grid& all_grid){
    double velocity_after_advect[smoke_simulation::physical_const::kGrid_num_x+1][smoke_simulation::physical_const::kGrid_num_y+1][2];
    //velocityのx成分を計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x+1;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            double velocity_y;
            if(ix<=0){
                velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]+all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy+1)])/2.0;
            }
            else if(ix>=smoke_simulation::physical_const::kGrid_num_x){
                velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix-1,iy)]
                           +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix-1,iy+1)])/2.0;
            }
            else{
                velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix-1,iy)]
                           +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix-1,iy+1)]
                           +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]
                           +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy+1)])/4.0;
            }
            //バックトレース先の位置
            double advected_x=(double)ix-((all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
            double advected_y=(double)iy-((velocity_y*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
            //バックトレース先の速度を補間(linear interpolationの場合)
//            velocity_after_advect[ix][iy][0]=linear_interpolation_x(advected_x, advected_y, all_grid);
            //バックトレース先の速度を補間(monotonic cubic interpolationの場合)
            velocity_after_advect[ix][iy][0]=monotonic_cubic_interpolation_x(advected_x, advected_y, all_grid);
        }
    }

    //velocityのy成分を計算
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y+1;iy++){
            double velocity_x;
            if(iy<=0){
                velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]+all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix+1,iy)])/2.0;
            }
            else if(iy>=smoke_simulation::physical_const::kGrid_num_y){
                velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy-1)]
                           +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix+1,iy-1)])/2.0;
            }
            else{
                velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy-1)]
                           +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]
                           +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix+1,iy-1)]
                           +all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix+1,iy)])/4.0;
            }
            //バックトレース先の位置
            double advected_x=(double)ix-((velocity_x*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
            double advected_y=(double)iy-((all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
            //バックトレース先の速度を補間(linear interpolationの場合)
//            velocity_after_advect[ix][iy][1]=linear_interpolation_y(advected_x, advected_y, all_grid);
            //バックトレース先の速度を補間(monotonic cubic interpolationの場合)
            velocity_after_advect[ix][iy][1]=monotonic_cubic_interpolation_y(advected_x, advected_y, all_grid);
        }
    }
    //計算結果をコピー
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x+1;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]=velocity_after_advect[ix][iy][0];
        }
    }
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y+1;iy++){
            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]=velocity_after_advect[ix][iy][1];
        }
    }
}

//Poisson eq. を解くことによる圧力の計算
void calc_pressure(Grid& all_grid){
    //グリッドの総数
    const int N=smoke_simulation::physical_const::kGrid_num_x
               *smoke_simulation::physical_const::kGrid_num_y;
    //係数行列の計算
    linear_algebra::sparse_matrix A(N,N);
//    linear_algebra::sparse_matrix_with_diagonal_element A(N,N);
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            if(ix-1>=0){
                A.input_element(get_voxel_center_index(ix, iy),get_voxel_center_index(ix-1, iy),1);
            }
            if(iy-1>=0){
                A.input_element(get_voxel_center_index(ix, iy),get_voxel_center_index(ix, iy-1),1);
            }
            A.input_element(get_voxel_center_index(ix, iy),get_voxel_center_index(ix, iy),-4);
            if(iy+1<=smoke_simulation::physical_const::kGrid_num_y-1){
                A.input_element(get_voxel_center_index(ix, iy),get_voxel_center_index(ix, iy+1),1);
            }
            if(ix+1<=smoke_simulation::physical_const::kGrid_num_x-1){
                A.input_element(get_voxel_center_index(ix, iy),get_voxel_center_index(ix+1, iy),1);
            }
        }
    }
    //連立方程式の右辺のベクトルの計算
    std::vector<double> b(N);
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
                 b[get_voxel_center_index(ix, iy)]=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix+1,iy)]-all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]
                                            +all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy+1)]-all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)])
                                            *smoke_simulation::physical_const::kCell_length
                                            /(smoke_simulation::physical_const::kDt);
        }
    }
    //時間計測用
/*
    std::chrono::system_clock::time_point  start, end;
    start = std::chrono::system_clock::now(); // 時間計測開始
*/
    //CG法により圧力場を得る
    linear_algebra::conjugate_gradient(A, b, all_grid.pressure, N, 10000, 0.0001);
//    linear_algebra::incomplete_cholesky_conjugate_gradient(A, b, all_grid.pressure, N, 10000, 0.0001);
/*
    end = std::chrono::system_clock::now();  // 時間計測終了
    std::ofstream writing_file;
    writing_file.open("length_of_time_CG_64_2d.dat", std::ios::app);
    writing_file << (std::chrono::duration_cast<std::chrono::microseconds>(end-start).count()) << std::endl;
    writing_file.close();
*/
    //gauss seidel法を使う場合(Aはsparse_matrix_with_diagonal_elementにする)
//    gauss_seidel(A,b,all_grid.pressure,N,200);
}


//pressure gradient termの計算
void calc_pressure_gradient_term(Grid& all_grid){
    //圧力の計算
    calc_pressure(all_grid);
    //圧力場に境界条件をセット
    set_boundary_pressure(all_grid);
    //pressure gradeint term によって速度場を更新
    for(int ix=1;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=1;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]-=(all_grid.pressure[get_voxel_center_index(ix, iy)]-all_grid.pressure[get_voxel_center_index(ix-1, iy)])
            *smoke_simulation::physical_const::kDt
            /(smoke_simulation::physical_const::kCell_length);
            all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]-=(all_grid.pressure[get_voxel_center_index(ix, iy)]-all_grid.pressure[get_voxel_center_index(ix, iy-1)])
            *smoke_simulation::physical_const::kDt
            /(smoke_simulation::physical_const::kCell_length);
        }
    }
}

//流体の 1 time step
void update_fluid_velocity(Grid& all_grid){
//固定境界条件の場合
//    add_source(all_grid);
    add_force_fluid(all_grid);
    set_boundary_velocity(all_grid);
    advect_fluid(all_grid);
    set_boundary_velocity(all_grid);
    calc_pressure_gradient_term(all_grid);
    set_boundary_velocity(all_grid);
}
}
