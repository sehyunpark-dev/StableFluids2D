#include <iostream>

#include <MACGrid2D.h>

MACGrid2D::MACGrid2D(
    int res_x, int res_y, double cell_size) : 
        res_x_(res_x), res_y_(res_y), cell_size_(cell_size)
{
    // Initialize the grid parameters
    u_.resize((res_x_ + 1) * res_y_, 0.0);
    v_.resize(res_x_ * (res_y_ + 1), 0.0);
    p_.resize(res_x_ * res_y_, 0.0);
}

double MACGrid2D::getU(int i, int j) const
{
    if (0 <= i && i < (res_x_ + 1) && 0 <= j && j < res_y_)
    {
        return u_[i + j * (res_x_ + 1)];
    }
    else
    {
        std::cerr << "[Error] Invalid U index" << std::endl;
        return -1;
    }
}

double MACGrid2D::getV(int i, int j) const
{
    if (0 <= i && i < res_x_ && 0 <= j && j < (res_y_ + 1))
    {
        return v_[i + j * res_x_];
    }
    else
    {
        std::cerr << "[Error] Invalid V index" << std::endl;
        return -1;
    }
}

double MACGrid2D::getP(int i, int j) const
{
    if (0 <= i && i < res_x_ && 0 <= j && j < res_y_)
    {
        return p_[i + j * res_x_];
    }
    else
    {
        std::cerr << "[Error] Invalid Pressure index" << std::endl;
        return -1;
    }
}

void MACGrid2D::setU(int i, int j, double val)
{
    if (0 <= i && i < (res_x_ + 1) && 0 <= j && j < res_y_)
    {
        u_[i + j * (res_x_ + 1)] = val;
    }
    else
    {
        std::cerr << "[Error] Invalid U index" << std::endl;
    }
}

void MACGrid2D::setV(int i, int j, double val)
{
    if (0 <= i && i < res_x_ && 0 <= j && j < (res_y_ + 1))
    {
        v_[i + j * res_x_] = val;
    }
    else
    {
        std::cerr << "[Error] Invalid V index" << std::endl;
    }
}

void MACGrid2D::setP(int i, int j, double val)
{
    if (0 <= i && i < res_x_ && 0 <= j && j < res_y_)
    {
        p_[i + j * res_x_] = val;
    }
    else
    {
        std::cerr << "[Error] Invalid Pressure index" << std::endl;
    }
}