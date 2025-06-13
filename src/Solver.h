#ifndef SOLVER_H
#define SOLVER_H

#include "MACGrid2D.h"
#include <vector>
#include "glm/glm.hpp"

class Solver
{
private:
    // Parameters of the grid and simulation
    const MACGrid2D *grid_;
    const int       grid_res_;      // Resolution of the grid
    const glm::vec2 grid_center_;   // Center of the grid in world coordinates
    const float     dx_;
    const float     dt_ = 0.03f;    // Time step for the simulation
    bool            isSimulating_ = false;

    // Physical properties of the fluid (Constants)
    const float density_   = 1.225f;    // Density(rho) at the cell center
    const float viscosity_ = 1e-6f;     // Viscosity of the fluid

    unsigned int frame_count_ = 0; // Frame count for the simulation

    // Velocity of cell faces (not cell centers!)
    std::vector<float> u_vec_; // the x-component of velocity at the vertical faces (size : height * (width + 1)))
    std::vector<float> v_vec_; // the y-component of velocity at the horizontal faces (size : (height + 1) * width))
    
    std::vector<float> pressure_vec_;   // Pressure(p) of each cell
    std::vector<float> divergence_vec_; // divergence(∇⋅u = ∂u/dx + ∂v/dx) of each cell
    
    // Smoke density of each sell (size : height * width, Range [0, 1])
    // This is used to visualize the smoke density in the scene
    // It is not a physical quantity, but a visualization aid.
    std::vector<float> smoke_vec_;

    // the vector that contains indices of source(smoke) cells
    std::vector<int> source_idx_;

    void addBodyForce();
    void advect();
    void diffuse();
    void project();
    void setVelocityBoundaryCondition();

    glm::vec2 getVelocity(const glm::vec2 &pos, 
        const std::vector<float> &u_vec, const std::vector<float> &v_vec);
    float getSmokeBilerpValue(const glm::vec2 &pos, const std::vector<float> &smoke_vector);
    float getUBilerpValue(const glm::vec2 &pos, const std::vector<float> &u);
    float getVBilerpValue(const glm::vec2 &pos, const std::vector<float> &v);

public:
    Solver(MACGrid2D *grid, int grid_res, float dx, float dt);
    ~Solver() = default;

    void initScene();
    void processSimulator();
    void step();
    void reset();
    
    // Getters & Setters

    inline std::vector<float> &getSmokeVector()
    {
        return smoke_vec_;
    }

    inline float getU(int x, int y) const
    {
        if (0 <= y && y < grid_res_ && 0 <= x && x < (grid_res_ + 1))
        {
            return u_vec_[y * (grid_res_ + 1) + x];
        }
        else
        {
            std::cerr << "[Error] Invalid U index" << std::endl;
            return -1;
        }
    }

    inline float getV(int x, int y) const
    {
        if (0 <= y && y < (grid_res_ + 1) && 0 <= x && x < grid_res_)
        {
            return v_vec_[y * grid_res_ + x];
        }
        else
        {
            std::cerr << "[Error] Invalid V index" << std::endl;
            return -1;
        }
    }

    inline float getP(int x, int y) const
    {
        if (0 <= y && y < grid_res_ && 0 <= x && x < grid_res_)
        {
            return pressure_vec_[y * grid_res_ + x];
        }
        else
        {
            std::cerr << "[Error] Invalid Pressure index" << std::endl;
            return -1;
        }
    }

    inline void setU(int x, int y, float val)
    {
        if (0 <= y && y < grid_res_ && 0 <= x && x < (grid_res_ + 1))
        {
            u_vec_[y * (grid_res_ + 1) + x] = val;
        }
        else
        {
            std::cerr << "[Error] Invalid U index" << std::endl;
        }
    }

    inline void setV(int x, int y, float val)
    {
        if (0 <= y && y < (grid_res_ + 1) && 0 <= x && x < grid_res_)
        {
            v_vec_[y * grid_res_ + x] = val;
        }
        else
        {
            std::cerr << "[Error] Invalid V index" << std::endl;
        }
    }

    inline void setP(int x, int y, float val)
    {
        if (0 <= y && y < grid_res_ && 0 <= x && x < grid_res_)
        {
            pressure_vec_[y * grid_res_ + x] = val;
        }
        else
        {
            std::cerr << "[Error] Invalid Pressure index" << std::endl;
        }
    }

    inline unsigned int getFrameCount()
    {
        return frame_count_;
    }
};

#endif