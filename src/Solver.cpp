#include "Solver.h"
#include <cmath>

Solver::Solver(MACGrid2D *grid, int grid_res, float dx, float dt) : 
    grid_(grid), grid_res_(grid->getRes()), grid_center_(grid->getGridCenter()),
    dx_(grid->getCellSize()), dt_(0.03), isSimulating_(false),
    density_(1.225f), viscosity_(1e-6f)
{
    u_vec_.resize(grid_res_ * (grid_res_ + 1), 0.0f);
    v_vec_.resize((grid_res_ + 1) * grid_res_, 0.0f);
    
    pressure_vec_.resize(grid_res_ * grid_res_, 0.0f);
    divergence_vec_.resize(grid_res_ * grid_res_, 0.0f);

    smoke_vec_.resize(grid_res_ * grid_res_, 0.0f);

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

    for (int i = -5; i <= 5; i++)
    {
        for (int j = -5; j <= 5; j++)
        {
            int x = source_center_idx_x + i;
            int y = source_center_idx_y + j;
            if (0 <= x && x < grid_res_ && 0 <= y && y < grid_res_)
            {
                // Set smoke density to 1.0 in the source area
                smoke_vec_[y * grid_res_ + x] = 1.0f;
                source_idx_.push_back(y * grid_res_ + x);
            }
        }
    }
}

void Solver::processSimulator()
{
    isSimulating_ = true;
}

void Solver::step()
{
    if (isSimulating_)
    {
        addBodyForce();
        setVelocityBoundaryCondition();

        advect();
        setVelocityBoundaryCondition();

        diffuse();
        setVelocityBoundaryCondition();
        
        project();
        setVelocityBoundaryCondition();

        frame_count_++;
    }
}

void Solver::reset()
{

}

void Solver::addBodyForce()
{
    glm::vec2 source_velocity = glm::vec2(3.0f, 0.0f);

    for (int i : source_idx_)
    {
        // Keep smoke density in the source 1.0
        smoke_vec_[i] = 1.0f;

        int x = i % grid_res_;
        int y = i / grid_res_;

        u_vec_[y * (grid_res_ + 1) + x] += source_velocity.x * dt_;
        u_vec_[y * (grid_res_ + 1) + (x+1)] += source_velocity.x * dt_;
        v_vec_[y * grid_res_ + x] += source_velocity.y * dt_;
        v_vec_[(y+1) * grid_res_ + x] += source_velocity.y * dt_;
    }
}

void Solver::advect()
{
    std::vector<float> u_new = u_vec_;
    std::vector<float> v_new = v_vec_;
    std::vector<float> smoke_new = smoke_vec_;
    
    // Advect smoke
    for (int y = 0; y < grid_res_; y++)
    {
        for (int x = 0; x < grid_res_; x++)
        {
            glm::vec2 cell_coord = grid_->getCellCoord(x, y);
            glm::vec2 cur_velocity = getVelocity(cell_coord, u_vec_, v_vec_);
            
            glm::vec2 prev_coord = cell_coord - cur_velocity * dt_;
            float prev_smoke = getSmokeBilerpValue(prev_coord, smoke_vec_);

            smoke_new[y * grid_res_ + x] = prev_smoke;
        }
    }

    // Advect U
    for (int y = 0; y < grid_res_; y++)
    {
        for (int x = 0; x < grid_res_ + 1; x++)
        {
            glm::vec2 u_coord = grid_->getUCoord(x, y);
            glm::vec2 cur_velocity = getVelocity(u_coord, u_vec_, v_vec_);

            glm::vec2 prev_coord = u_coord - cur_velocity * dt_;
            float prev_u = getUBilerpValue(prev_coord, u_vec_);

            u_new[y * (grid_res_ + 1) + x] = prev_u;
        }
    }
    
    // Advect V
    for (int y = 0; y < grid_res_ + 1; y++)
    {
        for (int x = 0; x < grid_res_; x++)
        {
            glm::vec2 v_coord = grid_->getVCoord(x, y);
            glm::vec2 cur_velocity = getVelocity(v_coord, u_vec_, v_vec_);

            glm::vec2 prev_coord = v_coord - cur_velocity * dt_;
            float prev_v = getVBilerpValue(prev_coord, v_vec_);

            v_new[y * grid_res_ + x] = prev_v;
        }
    }

    smoke_vec_ = smoke_new;
    u_vec_ = u_new;
    v_vec_ = v_new;
}

void Solver::diffuse()
{
    std::vector<float> smoke_new = smoke_vec_;
    std::vector<float> u_new = u_vec_;
    std::vector<float> v_new = v_vec_;

    int iter = 20;
    float constant = dt_ * viscosity_ / (dx_ * dx_);

    // Gauss-Seidel method
    for (int k = 0; k < iter; k++)
    {
        // Smoke diffusion
        for (int y = 1; y < grid_res_ - 1; y++) 
        {
            for (int x = 1; x < grid_res_ - 1; x++)
            {
                float neighbor_sum = 
                    smoke_new[y * grid_res_ + (x - 1)] + // left
                    smoke_new[y * grid_res_ + (x + 1)] + // right
                    smoke_new[(y - 1) * grid_res_ + x] + // top
                    smoke_new[(y + 1) * grid_res_ + x];  // bottom
                
                int idx = y * grid_res_ + x;
                smoke_new[idx] = (smoke_vec_[idx] + constant * neighbor_sum) / (1 + 4 * constant);
            }
        }

        // U diffusion
        for (int y = 1; y < grid_res_ - 1; y++) 
        {
            for (int x = 1; x < grid_res_; x++)
            {
                float neighbor_sum = 
                    u_new[y * (grid_res_ + 1) + (x - 1)] + // left
                    u_new[y * (grid_res_ + 1) + (x + 1)] + // right
                    u_new[(y - 1) * (grid_res_ + 1) + x] + // top
                    u_new[(y + 1) * (grid_res_ + 1) + x];  // bottom
                
                int idx = y * (grid_res_ + 1) + x;
                u_new[idx] = (u_vec_[idx] + constant * neighbor_sum) / (1 + 4 * constant);
            }
        }

        // V diffusion
        for (int y = 1; y < grid_res_; y++)
        {
            for (int x = 1; x < grid_res_ - 1; x++)
            {
                float neighbor_sum = 
                    v_new[y * grid_res_ + (x - 1)] + // left
                    v_new[y * grid_res_ + (x + 1)] + // right
                    v_new[(y - 1) * grid_res_ + x] + // top
                    v_new[(y + 1) * grid_res_ + x];  // bottom

                int idx = y * grid_res_ + x;
                v_new[idx] = (v_vec_[idx] + constant * neighbor_sum) / (1 + 4 * constant);
            }
        }
    }

    smoke_vec_ = smoke_new;
    u_vec_ = u_new;
    v_vec_ = v_new;
}

void Solver::project()
{
    pressure_vec_.resize(grid_res_ * grid_res_, 0.0f);
    divergence_vec_.resize(grid_res_ * grid_res_, 0.0f);

    // 1. Calculate divergence field : ∇⋅u = ∂u/dx + ∂v/dx
    for (int y = 0; y < grid_res_; y++) 
    {
        for (int x = 0; x < grid_res_; x++)
        {
            float du = u_vec_[y * (grid_res_ + 1) + (x + 1)] - u_vec_[y * (grid_res_ + 1) + x];
            float dv = v_vec_[y * grid_res_ + x] - v_vec_[(y + 1) * grid_res_ + x];

            divergence_vec_[y * grid_res_ + x] = (du + dv) / dx_;
        }
    }

    // 2. Solve Poisson equation with Gauss-Seidel method : (∇^2)p = ∇⋅u
    int iter = 20;
    pressure_vec_.resize(grid_res_ * grid_res_, 0.0f);

    for (int i = 0; i < iter; i++)
    {
        for (int y = 1; y < grid_res_ - 1; y++)
        {
            for (int x = 1; x < grid_res_ - 1; x++)
            {
                float p_right  = pressure_vec_[y * grid_res_ + (x + 1)];
                float p_left   = pressure_vec_[y * grid_res_ + (x - 1)];
                float p_top    = pressure_vec_[(y - 1) * grid_res_ + x];
                float p_bottom = pressure_vec_[(y + 1) * grid_res_ + x];

                float neighbor_sum = p_right + p_left + p_top + p_bottom;
                
                pressure_vec_[y * grid_res_ + x] = 
                    (neighbor_sum - (divergence_vec_[y * grid_res_ + x] * dx_ * dx_)) / 4.0f;
            }
        }

        // Apply Neumann boundary condition : ∂p/∂n = 0, 
        // which means that pressure of boundary cells are equal to their "one neighbor"!
        for (int y = 0; y < grid_res_; y++) {
            pressure_vec_[y * grid_res_ + 0] = pressure_vec_[y * grid_res_ + 1]; // Left
            pressure_vec_[y * grid_res_ + (grid_res_ - 1)] = pressure_vec_[y * grid_res_ + (grid_res_ - 2)]; // Right
        }

        for (int x = 0; x < grid_res_; x++) {
            pressure_vec_[0 * grid_res_ + x] = pressure_vec_[1 * grid_res_ + x]; // Top
            pressure_vec_[(grid_res_ - 1) * grid_res_ + x] = pressure_vec_[(grid_res_ - 2) * grid_res_ + x]; // Bottom
        }
    }

    // 3. Correct the velocity component : vel -= (dt/ρ) * ∇p
    // Let boundary u, v values free - they will be zero in setVelocityBoundaryCondition()
    for (int y = 0; y < grid_res_; y++) 
    {
        for (int x = 1; x < grid_res_; x++)
        {
            float p_left   = pressure_vec_[y * grid_res_ + (x - 1)];
            float p_right  = pressure_vec_[y * grid_res_ + x];
            float grad_p_x = (p_right - p_left) / dx_;
            u_vec_[y * (grid_res_ + 1) + x] -= (dt_ / density_) * grad_p_x;
        }
    }

    for (int y = 1; y < grid_res_; y++)
    {
        for (int x = 0; x < grid_res_; x++)
        {
            float p_top    = pressure_vec_[(y - 1) * grid_res_ + x];
            float p_bottom = pressure_vec_[y * grid_res_ + x];
            float grad_p_y = (p_top - p_bottom) / dx_;
            v_vec_[y * grid_res_ + x] -= (dt_ / density_) * grad_p_y;
        }
    }
}

void Solver::setVelocityBoundaryCondition()
{
    // No-stick condition : make velocity ⋅ solid wall normal = 0
    // That is, we need to make u or v component zero!
    
    for (int x = 0; x < grid_res_; x++)
    {
        v_vec_[0 * grid_res_ + x] = 0.0f;           // Top wall
        v_vec_[grid_res_ * grid_res_ + x] = 0.0f;   // Bottom wall
    }

    for (int y = 0; y < grid_res_; y++)
    {
        u_vec_[y * (grid_res_ + 1) + 0] = 0.0f;         // Left wall
        u_vec_[y * (grid_res_ + 1) + grid_res_] = 0.0f; // Right wall
    }
}

glm::vec2 Solver::getVelocity(const glm::vec2 &pos, 
    const std::vector<float> &u_vec, const std::vector<float> &v_vec)
{
    float u = getUBilerpValue(pos, u_vec);
    float v = getVBilerpValue(pos, v_vec);
    return glm::vec2(u, v);
}

float Solver::getSmokeBilerpValue(const glm::vec2 &pos, const std::vector<float> &smoke_vector)
{
    int idx = grid_->getCellIndex(pos);
    int x = idx % grid_res_;
    int y = idx / grid_res_;

    // pick the base cell which is located on the "left upper"
    int base_x = x;
    int base_y = y;

    glm::vec2 cell_center = grid_->getCellCoord(x, y);
    glm::vec2 diff = pos - cell_center;

    if (diff.x < 0) { base_x = x - 1; }
    if (diff.y > 0) { base_y = y - 1; }

    // Deal with the boundary condition
    base_x = glm::clamp(base_x, 0, grid_res_ - 1);
    base_y = glm::clamp(base_y, 0, grid_res_ - 1);
    
    int x_l = base_x;                              // Left x
    int x_r = glm::min(base_x + 1, grid_res_ - 1); // Right x
    int y_t = base_y;                              // Top y
    int y_b = glm::min(base_y + 1, grid_res_ - 1); // Bottom y

    float smoke_tl = smoke_vector[y_t * grid_res_ + x_l]; // Top-Left smoke value
    float smoke_tr = smoke_vector[y_t * grid_res_ + x_r]; // Top-Right smoke value
    float smoke_bl = smoke_vector[y_b * grid_res_ + x_l]; // Bottom-Left smoke value
    float smoke_br = smoke_vector[y_b * grid_res_ + x_r]; // Bottom-Right smoke value

    glm::vec2 tl_pos = grid_->getCellCoord(x_l, y_t);
    float tx = (pos.x - tl_pos.x) / dx_;
    float ty = (tl_pos.y - pos.y) / dx_;
    tx = glm::clamp(tx, 0.0f, 1.0f);
    ty = glm::clamp(ty, 0.0f, 1.0f);

    float val_top    = (1.0f - tx) * smoke_tl + tx * smoke_tr;
    float val_bottom = (1.0f - tx) * smoke_bl + tx * smoke_br;
    float val = (1.0f - ty) * val_top + ty * val_bottom;

    return val;
}

float Solver::getUBilerpValue(const glm::vec2 &pos, const std::vector<float> &u_vec)
{
    int idx = grid_->getUIndex(pos);
    int x = idx % (grid_res_ + 1);
    int y = idx / (grid_res_ + 1);

    // pick the base u index which is located on the "left upper"
    int base_x = x;
    int base_y = y;

    glm::vec2 u_pos = grid_->getUCoord(x, y);
    glm::vec2 diff = pos - u_pos;

    if (diff.x < 0) { base_x = x - 1; }
    if (diff.y > 0) { base_y = y - 1; }

    // Deal with the boundary condition
    base_x = glm::clamp(base_x, 0, grid_res_);
    base_y = glm::clamp(base_y, 0, grid_res_ - 1);

    int x_l = base_x;
    int x_r = glm::min(base_x + 1, grid_res_);
    int y_t = base_y;
    int y_b = glm::min(base_y + 1, grid_res_ - 1);

    float u_tl = u_vec[y_t * (grid_res_ + 1) + x_l]; // Top-Left u value
    float u_tr = u_vec[y_t * (grid_res_ + 1) + x_r]; // Top-Right u value
    float u_bl = u_vec[y_b * (grid_res_ + 1) + x_l]; // Bottom-Left u value
    float u_br = u_vec[y_b * (grid_res_ + 1) + x_r]; // Bottom-Right u value

    glm::vec2 tl_pos = grid_->getUCoord(x_l, y_t);
    float tx = (pos.x - tl_pos.x) / dx_;
    float ty = (tl_pos.y - pos.y) / dx_;
    tx = glm::clamp(tx, 0.0f, 1.0f);
    ty = glm::clamp(ty, 0.0f, 1.0f);

    float val_top    = (1.0f - tx) * u_tl + tx * u_tr;
    float val_bottom = (1.0f - tx) * u_bl + tx * u_br;
    float val = (1.0f - ty) * val_top + ty * val_bottom;

    return val;
}

float Solver::getVBilerpValue(const glm::vec2 &pos, const std::vector<float> &v_vec)
{
    int idx = grid_->getVIndex(pos);
    int x = idx % grid_res_;
    int y = idx / grid_res_;

    // pick the base v index which is located on the "left upper"
    int base_x = x;
    int base_y = y;

    glm::vec2 v_pos = grid_->getVCoord(x, y);
    glm::vec2 diff = pos - v_pos;

    if (diff.x < 0) { base_x = x - 1; }
    if (diff.y > 0) { base_y = y - 1; }

    // Deal with the boundary condition
    base_x = glm::clamp(base_x, 0, grid_res_ - 1);
    base_y = glm::clamp(base_y, 0, grid_res_);

    int x_l = base_x;
    int x_r = glm::min(base_x + 1, grid_res_ - 1);
    int y_t = base_y;
    int y_b = glm::min(base_y + 1, grid_res_);

    float v_tl = v_vec[y_t * grid_res_ + x_l]; // Top-Left v value
    float v_tr = v_vec[y_t * grid_res_ + x_r]; // Top-Right v value
    float v_bl = v_vec[y_b * grid_res_ + x_l]; // Bottom-Left v value
    float v_br = v_vec[y_b * grid_res_ + x_r]; // Bottom-Right v value

    glm::vec2 tl_pos = grid_->getVCoord(x_l, y_t);
    float tx = (pos.x - tl_pos.x) / dx_;
    float ty = (tl_pos.y - pos.y) / dx_;
    tx = glm::clamp(tx, 0.0f, 1.0f);
    ty = glm::clamp(ty, 0.0f, 1.0f);

    float val_top    = (1.0f - tx) * v_tl + tx * v_tr;
    float val_bottom = (1.0f - tx) * v_bl + tx * v_br;
    float val = (1.0f - ty) * val_top + ty * val_bottom;

    return val;
}