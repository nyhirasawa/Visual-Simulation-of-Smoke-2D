#include "move_substances.h"

#include <iostream>

#include "linear_solver.h"
#include "physical_const.h"
#include "utils.h"

namespace smoke_simulation{
//外力項の計算
void add_force_substances(Grid& all_grid){

}

//substance densityのlinear interpolation
double linear_interpolation_substances(double advected_x, double advected_y, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    if(advected_index_x<1){
        advected_index_x=1;
        advected_x=(advected_index_x);
//        advected_x=(advected_index_x+0.5);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-1){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-2;
        advected_x=(advected_index_x);
//        advected_x=(advected_index_x+0.5);
    }
    if(advected_index_y<1){
        advected_index_y=1;
        advected_y=(advected_index_y);
//        advected_y=(advected_index_y+0.5);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-1){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-2;
        advected_y=(advected_index_y);
//        advected_y=(advected_index_y+0.5);
    }
    double a0, a1, b0, b1;
    a0=advected_x-advected_index_x;
    a1=1.0-a0;
    b0=advected_y-advected_index_y;
    b1=1.0-b0;
//    if(advected_index_x+1>=smoke_simulation::physical_const::kGrid_num_x
//     &&advected_index_y+1>=smoke_simulation::physical_const::kGrid_num_y){
//        return all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y)];
//    }
//    else if(advected_index_x+1>=smoke_simulation::physical_const::kGrid_num_x){
//        return b1*all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y)]
//              +b0*all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y+1)];
//    }
//    else if(advected_index_y+1>=smoke_simulation::physical_const::kGrid_num_y){
//        return a1*all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y)]
//              +a0*all_grid.substance_density[get_voxel_center_index(advected_index_x+1,advected_index_y)];
//    }
//    else{
        return a1*(b1*all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y)]
              +b0*all_grid.substance_density[get_voxel_center_index(advected_index_x,advected_index_y+1)])
              +a0*(b1*all_grid.substance_density[get_voxel_center_index(advected_index_x+1,advected_index_y)]
              +b0*all_grid.substance_density[get_voxel_center_index(advected_index_x+1,advected_index_y+1)]);
//    }
}

//substance densityのmonotonic cubic interpolation
double monotonic_cubic_substances(double advected_x, double advected_y, Grid& all_grid){
    int advected_index_x=(int)(advected_x);
    int advected_index_y=(int)(advected_y);
    if(advected_index_x<0){
        advected_index_x=0;
        advected_x=(double)(advected_index_x+0.5);
    }
    if(advected_index_x>=smoke_simulation::physical_const::kGrid_num_x){
        advected_index_x=smoke_simulation::physical_const::kGrid_num_x-1;
        advected_x=(double)(advected_index_x+0.5);
    }
    if(advected_index_y<0){
        advected_index_y=0;
        advected_y=(double)(advected_index_y+0.5);
    }
    if(advected_index_y>=smoke_simulation::physical_const::kGrid_num_y){
        advected_index_y=smoke_simulation::physical_const::kGrid_num_y-1;
        advected_y=(double)(advected_index_y+0.5);
    }
    //範囲外を参照しないようにする処置
    if(advected_index_x<=1||advected_index_x>=smoke_simulation::physical_const::kGrid_num_x-2
     ||advected_index_y<=1||advected_index_y>=smoke_simulation::physical_const::kGrid_num_y-2){
         return linear_interpolation_substances(advected_x, advected_y, all_grid);
    }
    //monotonic cubic interpolation の処理
    else{
        double dk_0;
        double dk_1;
        double delta_k;
        double interpolated_substance_density[4];
        //y軸方向に関するinterpolation
        for(int i=0;i<4;i++){
            dk_0  =  (all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i, advected_index_y+1)]
                     -all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i, advected_index_y-1)])/2.0;
            dk_1  =  (all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i, advected_index_y+2)]
                     -all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i, advected_index_y)])/2.0;
            delta_k = all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i, advected_index_y+1)]
                     -all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i, advected_index_y)];
            double b0 = advected_y-advected_index_y;
            if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
                interpolated_substance_density[i]=(1.0-b0)*all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i,advected_index_y)]
                                                 +b0      *all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i,advected_index_y+1)];
            }

            else{
                //論文の式は間違っている。正しくはこっち
                interpolated_substance_density[i]=(dk_0+dk_1-2*delta_k)*(b0*b0*b0)+(3*delta_k-2*dk_0-dk_1)*(b0*b0)+dk_0*b0
                                                 +all_grid.substance_density[get_voxel_center_index(advected_index_x-1+i, advected_index_y)];
            }

        }
        //x軸方向に関するinterpolation
        dk_0=(interpolated_substance_density[2]-interpolated_substance_density[0])/2.0;
        dk_1=(interpolated_substance_density[3]-interpolated_substance_density[1])/2.0;
        delta_k=(interpolated_substance_density[2]-interpolated_substance_density[1]);
        double a0=advected_x-advected_index_x;
        //3つのスロープの符号が異なる場合
        if(dk_0*dk_1<0.0||dk_0*delta_k<0.0||dk_1*delta_k<0.0){
            return (1.0-a0)*interpolated_substance_density[1]
                   +a0     *interpolated_substance_density[2];
        }
        else{
            //論文の式は間違っている。正しくはこっち
            return (dk_0+dk_1-2*delta_k)*(a0*a0*a0)+(3*delta_k-2*dk_0-dk_1)*(a0*a0)+dk_0*a0
                    +interpolated_substance_density[1];
        }
    }
}

//流体の速度場によってsubstanceが運ばれる項(流体のadvect項に相当)
//速度場を時間 -dt だけバックトレースして計算する
void transport_substances(Grid& all_grid){
    double substance_density_after_advect[smoke_simulation::physical_const::kGrid_num_x][smoke_simulation::physical_const::kGrid_num_y];
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            double velocity_x=(all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix,iy)]+all_grid.velocity_in_voxel_face_x[get_voxel_face_index_x(ix+1,iy)])/2.0;
            double velocity_y=(all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy)]+all_grid.velocity_in_voxel_face_y[get_voxel_face_index_y(ix,iy+1)])/2.0;
            //バックトレース先の座標
            double advected_x=ix-((velocity_x*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
            double advected_y=iy-((velocity_y*smoke_simulation::physical_const::kDt)/smoke_simulation::physical_const::kCell_length);
            //バックトレース先の座標のindex
            //バックトレース先の速度を補間する
//            substance_density_after_advect[ix][iy]=linear_interpolation_substances(advected_x, advected_y, all_grid);
            substance_density_after_advect[ix][iy]=monotonic_cubic_substances(advected_x, advected_y, all_grid);
        }
    }
    //計算結果をコピー
    for(int ix=0;ix<smoke_simulation::physical_const::kGrid_num_x;ix++){
        for(int iy=0;iy<smoke_simulation::physical_const::kGrid_num_y;iy++){
            all_grid.substance_density[get_voxel_center_index(ix,iy)]=substance_density_after_advect[ix][iy];
        }
    }
}

//上の4ステップをまとめただけの関数(substance densityの1時間ステップ分の更新に相当)
void move_substances(Grid& all_grid){
//    add_force_substances(all_grid);
//    diffuse_substances(all_grid);
    transport_substances(all_grid);
//    dissipate_substances(all_grid);
}
}
