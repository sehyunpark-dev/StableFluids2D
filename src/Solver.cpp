#include "Solver.h"
#include <cmath>

Solver::Solver(MACGrid2D *grid, int grid_res, float dx, float dt) : 
    grid_(grid), grid_res_(grid->getRes()), grid_center_(grid->getGridCenter()),
    dx_(grid->getCellSize()), dt_(0.03), isSimulating_(false),
    density_(1.225f), viscosity_(1e-6f),
    solid_wall_max_(glm::vec2(grid_res_ / 2.0f, grid_res_ / 2.0f) + grid_center_),
    solid_wall_min_(glm::vec2(-grid_res_ / 2.0f, -grid_res_ / 2.0f) + grid_center_)
{
    u_.resize(grid_res_ * (grid_res_ + 1), 0.0f);
    v_.resize((grid_res_ + 1) * grid_res_, 0.0f);
    
    pressure_.resize(grid_res_ * grid_res_, 0.0f);
    velocity_.resize(grid_res_ * grid_res_, glm::vec2(0.0f, 0.0f));
    smoke_density_.resize(grid_res_ * grid_res_, 0.0f);

    source_idx_.clear();

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
                source_idx_.push_back(y * grid_res_ + x);
            }
        }
    }
}

void Solver::step()
{
    if (isSimulating_)
    {
        addBodyForce();
        setBoundaryCondition();

        advect();


        diffuse();
        project();
    }
}

void Solver::reset()
{

}

void Solver::addBodyForce()
{
    glm::vec2 source_velocity = glm::vec2(10.0f, 0.0f);

    for (int i : source_idx_)
    {
        // Keep smoke density in the source 1.0
        smoke_density_[i] = 1.0f;

        int x = i % grid_res_;
        int y = i / grid_res_;

        u_[y * (grid_res_ + 1) + x] += source_velocity.x * dt_;
        u_[y * (grid_res_ + 1) + (x + 1)] += source_velocity.x * dt_;
        v_[y * grid_res_ + x] += source_velocity.y * dt_;
        v_[(y + 1) * grid_res_ + x] += source_velocity.y * dt_;
    }
}

void Solver::advect()
{
    std::vector<float> u_new = u_;
    std::vector<float> v_new = v_;
    std::vector<float> smoke_density_new = smoke_density_;
    
    // Advect smoke density
    for (int y = 0; y < grid_res_; y++)
    {
        for (int x = 0; x < grid_res_; x++)
        {
            glm::vec2 cell_coord = grid_->getCellCoord(x, y);
            glm::vec2 cur_velocity = getVelocity(cell_coord);
        }
    }
}

void Solver::diffuse()
{

}

void Solver::project()
{

}

void Solver::setBoundaryCondition()
{
    // No-stick condition : make velocity ⋅ solid wall normal = 0
    // That is, we need to make u or v component zero!
    
    for (int x = 0; x < grid_res_; x++)
    {
        v_[0 * grid_res_ + x] = 0.0f;           // Top wall
        v_[grid_res_ * grid_res_ + x] = 0.0f;   // Bottom wall
    }

    for (int y = 0; y < grid_res_; y++)
    {
        u_[y * (grid_res_ + 1) + 0] = 0.0f;         // Left wall
        u_[y * (grid_res_ + 1) + grid_res_] = 0.0f; // Right wall
    }
}

std::vector<float> &Solver::getSmokeDensityVector()
{
    return smoke_density_;
}

glm::vec2 Solver::getVelocity(glm::vec2 &coord)
{
    
}

float Solver::getSmokeBilerpValue(const glm::vec2 &pos, const std::vector<float> &smoke_vector)
{
    int idx = grid_->getCellIndex(pos);
    int x = idx % grid_res_;
    int y = idx % grid_res_;

    glm::vec2 cell_center = grid_->getCellCoord(x, y);
    glm::vec2 diff = pos - cell_center;

    
    float neighbors[4];
    
    // For boundary condition...
    int x_plus  = x == grid_res_ ? x : x+1;
    int x_minus = x == 0         ? 0 : x-1;
    int y_plus  = y == grid_res_ ? y : y+1;
    int y_minus = y == 0         ? 0 : y-1;
    
    if (diff.x >= 0 && diff.y >= 0)     // the pos is located in the right lower
    {
        neighbors[0] = smoke_vector[y * grid_res_ + x];
        neighbors[1] = smoke_vector[y * grid_res_ + x+1];
        neighbors[2] = smoke_vector[(y+1) * grid_res_ + x];
        neighbors[3] = smoke_vector[(y+1) * grid_res_ + x+1];
    }
    else if (diff.x >= 0 && diff.y < 0) // the pos is located in the right upper
    {
        neighbors[0] = smoke_vector[y * grid_res_ + x];
        neighbors[1] = smoke_vector[y * grid_res_ + x+1];
        neighbors[2] = smoke_vector[(y-1) * grid_res_ + x];
        neighbors[3] = smoke_vector[(y-1) * grid_res_ + x+1];
    }
    else if (diff.x < 0 && diff.y >= 0) // the pos is located in the left lower
    {
        neighbors[0] = smoke_vector[y * grid_res_ + x];
        neighbors[1] = smoke_vector[y * grid_res_ + x-1];
        neighbors[2] = smoke_vector[(y+1) * grid_res_ + x];
        neighbors[3] = smoke_vector[(y+1) * grid_res_ + x-1];
    }
    else                                // the pos is located in the left upper
    {
        neighbors[0] = smoke_vector[y * grid_res_ + x];
        neighbors[1] = smoke_vector[y * grid_res_ + x-1];
        neighbors[2] = smoke_vector[(y-1) * grid_res_ + x];
        neighbors[3] = smoke_vector[(y-1) * grid_res_ + x-1];
    }

    float x_frac = diff.x / dx_;
    float y_frac = diff.y / dx_;

    float sx_0 = (1.0f - x_frac) * neighbors[0] + x_frac * neighbors[1];
    float sx_1 = (1.0f - x_frac) * neighbors[2] + x_frac * neighbors[3];
    float val = (1.0f - y_frac) * sx_0 + y_frac * sx_1;

    return val;
}