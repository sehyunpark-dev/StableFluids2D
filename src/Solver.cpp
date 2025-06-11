#include "Solver.h"

Solver::Solver(MACGrid2D *grid, int grid_res, float dx, float dt) : 
    grid_(grid), grid_res_(grid->getRes()), grid_center_(grid->getGridCenter()),
    dx_(grid->getCellSize()), dt_(0.03), density_(1.225f), viscosity_(0.01f),
    solid_wall_max_(glm::vec2(grid_res_ / 2.0f, grid_res_ / 2.0f) + grid_center_),
    solid_wall_min_(glm::vec2(-grid_res_ / 2.0f, -grid_res_ / 2.0f) + grid_center_)
{
    u_.resize(grid_res_ * (grid_res_ + 1), 0.0f);
    v_.resize((grid_res_ + 1) * grid_res_, 0.0f);
    
    pressure_.resize(grid_res_ * grid_res_, 0.0f);
    velocity_.resize(grid_res_ * grid_res_, glm::vec2(0.0f, 0.0f));
    body_forces_.resize(grid_res_ * grid_res_, glm::vec2(0.0f, 0.0f));
    smoke_density_.resize(grid_res_ * grid_res_, 0.0f);

    initScene();
}

void Solver::initScene()
{
    // Add initial smoke to the first frame of the scene
    glm::vec2 source_center = glm::vec2(-grid_res_ / 4.0f * dx_, 0.0f) + grid_center_;
    int source_center_idx = grid_->getCellIndex(source_center);
    int source_center_idx_x = source_center_idx % grid_res_;
    int source_center_idx_y = source_center_idx / grid_res_;

    for (int i = -2; i <= 2; i++)
    {
        for (int j = -2; j <= 2; j++)
        {
            int x = source_center_idx_x + i;
            int y = source_center_idx_y + j;
            if (0 <= x && x < grid_res_ && 0 <= y && y < grid_res_)
            {
                // Set smoke density to 1.0 in the source area
                smoke_density_[y * grid_res_ + x] = 1.0f;
                std::cout << "Smoke density set at (" << x << ", " << y << ") to " << smoke_density_[y * grid_res_ + x] << " / real index: " << y * grid_res_ + x << std::endl;
            }
        }
    }
}