#ifndef MOVING_SUBSTANCES_H
#define MOVING_SUBSTANCES_H
#include <vector>
#include "grid.h"
#include "physical_const.h"

namespace smoke_simulation{
//外力項の計算
void add_force_substances(Grid& all_grid);
//substance densityのlinear interpolation
double linear_interpolation_substances(double advected_x, double advected_y, Grid& all_grid);
//substance densityのmonotonic cubic interpolation
double monotonic_cubic_substances(double advected_x, double advected_y, Grid& all_grid);
//流体の速度場によってsubstanceが運ばれる項(流体のadvect項に相当)
void transport_substances(Grid& all_grid);
//diffuse項の計算
void diffuse_substances(Grid& all_grid);
//dissipate項の計算
void dissipate_substances(Grid& all_grid);
//上の4ステップをまとめただけの関数(substance densityの1時間ステップ分の更新に相当)
void move_substances(Grid& all_grid);
}//namespace smoke_simulation

#endif //MOVING_SUBSTANCES_H
