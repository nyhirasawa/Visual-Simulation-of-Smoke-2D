#ifndef SMOKE_SIMULATION_GRID_H
#define SMOKE_SIMULATION_GRID_H

#include <vector>

namespace smoke_simulation{
//セルの境界条件
enum BoundaryCondition{
    FLUID=0,
    WALL=1,
};

//系の物理量が乗るグリッドの定義
class Grid{
public:
    const int Grid_num_x, Grid_num_y;
    double *velocity_in_voxel_face_x, *velocity_in_voxel_face_y;
    double **velocity_in_cell_center;
    double **vorticity;
    double **normalized_vorticity_amplitude_gradient;
    std::vector<double> pressure;
    double *temperature;
    double *substance_density;
    double **external_force_field;
    BoundaryCondition *boundary_condition;
    Grid(int nx, int ny);
    ~Grid();
};
}
#endif//SMOKE_SIMULATION_GRID_H
