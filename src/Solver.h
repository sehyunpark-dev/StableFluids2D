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
    const float     dt_ = 0.03f;     // Time step for the simulation

    // Physical properties of the fluid (Constants)
    const float density_   = 1.225f;    // Density(rho) at the cell center
    const float viscosity_ = 0.01f;     // Viscosity of the fluid

    // For Boundary conditions
    const glm::vec2 solid_wall_max_;    // Maximum coordinates of the solid wall (max_x, max_y)
    const glm::vec2 solid_wall_min_;    // Minimum coordinates of the solid wall (min_x, min_y)

    unsigned int frame_count_ = 0; // Frame count for the simulation

    // Velocity of cell faces (not cell centers!)
    std::vector<float> u_; // the x-component of velocity at the vertical faces (size : height * (width + 1)))
    std::vector<float> v_; // the y-component of velocity at the horizontal faces (size : (height + 1) * width))
    
    // Physical quantities of each cell (size : height * width)
    std::vector<float>     pressure_;   // Pressure(p) at the cell center
    std::vector<glm::vec2> velocity_;   // Velocity(u) at the cell center
    std::vector<glm::vec2> body_forces_; // Body forces(f) at the cell center (e.g. gravity)
    
    // Smoke density of each sell (size : height * width, Range [0, 1])
    // This is used to visualize the smoke density in the scene
    // It is not a physical quantity, but a visualization aid.
    std::vector<float> smoke_density_;

public:
    Solver(MACGrid2D *grid, int grid_res, float dx, float dt);
    ~Solver() = default;

    void initScene();

    void addBodyForce();
    void advect();
    void diffuse();
    void project();
    
    // Getters & Setters

    inline float getU(int x, int y) const
    {
        if (0 <= y && y < grid_res_ && 0 <= x && x < (grid_res_ + 1))
        {
            return u_[y * (grid_res_ + 1) + x];
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
            return v_[y * grid_res_ + x];
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
            return pressure_[y * grid_res_ + x];
        }
        else
        {
            std::cerr << "[Error] Invalid Pressure index" << std::endl;
            return -1;
        }
    }

    std::vector<float>& getSmokeDensityVector()
    {
        return smoke_density_;
    }

    inline void setU(int x, int y, float val)
    {
        if (0 <= y && y < grid_res_ && 0 <= x && x < (grid_res_ + 1))
        {
            u_[y * (grid_res_ + 1) + x] = val;
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
            v_[y * grid_res_ + x] = val;
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
            pressure_[y * grid_res_ + x] = val;
        }
        else
        {
            std::cerr << "[Error] Invalid Pressure index" << std::endl;
        }
    }
};

#endif