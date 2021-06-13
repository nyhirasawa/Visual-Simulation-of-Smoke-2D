smoke_simulation: main.cpp grid.cpp initialize_grid.cpp utils.cpp draw_substance_density.cpp update_fluid_velocity.cpp linear_solver.cpp sparse_matrix.cpp move_substances.cpp
	clang++ main.cpp grid.cpp initialize_grid.cpp utils.cpp draw_substance_density.cpp update_fluid_velocity.cpp linear_solver.cpp sparse_matrix.cpp move_substances.cpp `pkg-config --cflags --libs opencv4` -std=c++11 -O3
