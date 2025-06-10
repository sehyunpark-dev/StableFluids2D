#include "MACGrid2D.h"

MACGrid2D::MACGrid2D(int res, float cell_size) : res_(res), cell_size_(cell_size)
{
    // Initialize the grid parameters
    u_.resize((res_ + 1) * res_, 0.0);
    v_.resize(res_ * (res_ + 1), 0.0);
    p_.resize(res_ * res_, 0.0);
    cell_coord_.resize(res_ * res_, glm::vec2(0.0, 0.0));
    
    float center_idx = (res_ - 1) / 2.0f;
    for (int i = 0; i < res_; i++)
    {
        for (int j = 0; j < res_; j++)
        {
            float y = (center_idx - i) * cell_size;
            float x = (j - center_idx) * cell_size;
            cell_coord_[i * res_ + j] = glm::vec2(x, y);
        }
    }
}

int MACGrid2D::getRes() const
{
    return res_;
}

float MACGrid2D::getCellSize() const
{
    return cell_size_;
}

glm::vec2 MACGrid2D::getGridCenter() const
{
    return grid_center_;
}

std::vector<glm::vec2> MACGrid2D::getCellCoord() const
{
    return cell_coord_;
}

glm::vec2 MACGrid2D::getCellCoord(int i, int j) const
{
    if (0 <= i && i < res_ && 0 <= j && j < res_)
    {
        return cell_coord_[i * res_ + j];
    }
    else
    {
        std::cerr << "[Error] Invalid cell coordinate index" << std::endl;
        return glm::vec2(-1.0f, -1.0f);
    }
}


double MACGrid2D::getU(int i, int j) const
{
    if (0 <= i && i < res_ && 0 <= j && j < (res_ + 1))
    {
        return u_[i * (res_ + 1) + j];
    }
    else
    {
        std::cerr << "[Error] Invalid U index" << std::endl;
        return -1;
    }
}

double MACGrid2D::getV(int i, int j) const
{
    if (0 <= i && i < (res_ + 1) && 0 <= j && j < res_)
    {
        return v_[i * res_ + j];
    }
    else
    {
        std::cerr << "[Error] Invalid V index" << std::endl;
        return -1;
    }
}

double MACGrid2D::getP(int i, int j) const
{
    if (0 <= i && i < res_ && 0 <= j && j < res_)
    {
        return p_[i * res_ + j];
    }
    else
    {
        std::cerr << "[Error] Invalid Pressure index" << std::endl;
        return -1;
    }
}

void MACGrid2D::setU(int i, int j, double val)
{
    if (0 <= i && i < res_ && 0 <= j && j < (res_ + 1))
    {
        u_[i * (res_ + 1) + j] = val;
    }
    else
    {
        std::cerr << "[Error] Invalid U index" << std::endl;
    }
}

void MACGrid2D::setV(int i, int j, double val)
{
    if (0 <= i && i < (res_ + 1) && 0 <= j && j < res_)
    {
        v_[i * res_ + j] = val;
    }
    else
    {
        std::cerr << "[Error] Invalid V index" << std::endl;
    }
}

void MACGrid2D::setP(int i, int j, double val)
{
    if (0 <= i && i < res_ && 0 <= j && j < res_)
    {
        p_[i * res_ + j] = val;
    }
    else
    {
        std::cerr << "[Error] Invalid Pressure index" << std::endl;
    }
}