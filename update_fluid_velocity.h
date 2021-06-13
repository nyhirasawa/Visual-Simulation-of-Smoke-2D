#ifndef UPDATE_FLUID_VELOCITY_H
#define UPDATE_FLUID_VELOCITY_H
#include <vector>
#include "fftw3.h"

#include "grid.h"
#include "physical_const.h"

namespace smoke_simulation{
void add_source(Grid& all_grid);

//壁における速度場をセットする関数
void set_boundary_velocity(Grid& all_grid);

//壁における圧力場をセットする関数
void set_boundary_pressure(Grid& all_grid);

//(vorticityの計算に使う)cell centered velocityの計算
void calc_cell_centered_velocity(Grid& all_grid);
//vorticity confinement termの計算
void calc_vorticity_confinement(Grid& all_grid);
//外力場による速度場の更新
void update_fluid_velocity_by_external_force(Grid& all_grid);
//外力項の計算(vorticity confinement termの計算と外力場による速度場の更新をまとめただけの関数)
void add_force_fluid(Grid& all_grid);

//速度場のlinear interpolation
double linear_interpolation_x(double advected_x, double advected_y, Grid& all_grid);
double linear_interpolation_y(double advected_x, double advected_y, Grid& all_grid);
//速度場のmonotonic cubic interpolation
double monotonic_cubic_interpolation_x(double advected_x, double advected_y, Grid& all_grid);
double monotonic_cubic_interpolation_y(double advected_x, double advected_y, Grid& all_grid);

//advect項の計算
void advect_fluid(Grid& all_grid);

//Poisson eq. を解くことによる圧力の計算
void calc_pressure(Grid& all_grid);
//pressure gradient termの計算
void calc_pressure_gradient_term(Grid& all_grid);
//上記4ステップをまとめただけの関数
void update_fluid_velocity(Grid& all_grid);
}

#endif //UPDATE_FLUID_VELOCITY_H
