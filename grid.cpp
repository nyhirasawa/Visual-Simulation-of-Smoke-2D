#include "grid.h"

namespace smoke_simulation{

Grid::Grid(int nx, int ny): Grid_num_x(nx), Grid_num_y(ny){
    //substance_densityのメモリ確保
    substance_density = new double[nx*ny];

    //BoundaryConditionのメモリ確保
    boundary_condition = new BoundaryCondition[nx*ny];

    //velocityのメモリ確保
    velocity_in_voxel_face_x = new double[(nx+1)*ny];
    velocity_in_voxel_face_y = new double[nx*(ny+1)];
    velocity_in_cell_center  = new double*[nx*ny];
    for(int i=0; i<nx*ny;i++){
        velocity_in_cell_center[i]=new double[3];
    }

    //vorticityのメモリ確保
    vorticity = new double*[nx*ny];
    for(int i=0; i<nx*ny;i++){
        vorticity[i] = new double[3];
    }

    //normalized_vorticity_amplitude_gradientのメモリ確保
    normalized_vorticity_amplitude_gradient = new double*[nx*ny];
    for(int i=0; i<nx*ny;i++){
        normalized_vorticity_amplitude_gradient[i] = new double[3];
    }

    //温度場temperatureのメモリ確保
    temperature = new double[nx*ny];

    //pressureのメモリ確保
    pressure.resize(nx*ny);


    //cell center上のexternal force fieldのメモリ確保
//    external_force_field = new double[(2*nx+1)*(2*ny+1)];
    external_force_field = new double*[nx*ny];
    for(int i=0; i<nx*ny;i++){
        external_force_field[i] = new double[3];
    }
}

Grid::~Grid(){
    //substance_densityのメモリ解放
    delete[] substance_density;

    //BoundaryConditionのメモリ解放
    delete[] boundary_condition;

    //velocityのメモリ解放
    delete[] velocity_in_voxel_face_x;
    delete[] velocity_in_voxel_face_y;
    for(int i=0; i<Grid_num_x*Grid_num_y;i++){
        delete[] velocity_in_cell_center[i];
    }
    delete[] velocity_in_cell_center;

    //vorticityのメモリ解放
    for(int i=0; i<Grid_num_x*Grid_num_y;i++){
        delete[] vorticity[i];
    }
    delete[] vorticity;

    //normalized_vorticity_amplitude_gradientのメモリ解放
    for(int i=0; i<Grid_num_x*Grid_num_y;i++){
        delete[] normalized_vorticity_amplitude_gradient[i];
    }
    delete[] normalized_vorticity_amplitude_gradient;

    //温度場temperatureのメモリ確保
    delete[] temperature;

//    delete[] pressure;

    //external force fieldのメモリ確保
    for(int i=0; i<Grid_num_x*Grid_num_y;i++){
        delete[] external_force_field[i];
    }
    delete[] external_force_field;
}

}
