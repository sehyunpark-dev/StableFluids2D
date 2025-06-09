#ifndef MAC_GRID_2D_H
#define MAC_GRID_2D_H

#include <vector>
#include <glm/glm.hpp>

class MACGrid2D
{
private:
    int res_x_, res_y_; // Resolution of the grid in x and y directions (e.g. 128x128)
    float cell_size_; // Size of each cell in the grid (e.g. 0.01))

    std::vector<double> u_; // x-component of velocity (size : (width + 1) * height))
    std::vector<double> v_; // y-component of velocity (size : width * (height + 1)))
    std::vector<double> p_; // pressure of each cell (size : width * height)

public:
    MACGrid2D(int res_x, int res_y, double cell_size);

    double getU(int i, int j) const;
    double getV(int i, int j) const;
    double getP(int i, int j) const;

    void setU(int i, int j, double val);
    void setV(int i, int j, double val);
    void setP(int i, int j, double val);
};

#endif