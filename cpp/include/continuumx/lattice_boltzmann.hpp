#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <algorithm>

namespace continuumx {

// D2Q9 lattice directions and weights
constexpr int Q = 9;
constexpr std::array<int, Q> cx = {0, 1, 0, -1, 0, 1, -1, -1, 1};
constexpr std::array<int, Q> cy = {0, 0, 1, 0, -1, 1, 1, -1, -1};
constexpr std::array<double, Q> w = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
                                     1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
constexpr double cs2 = 1.0 / 3.0;

class LatticeBoltzmannD2Q9 {
public:
    LatticeBoltzmannD2Q9(int nx, int ny, double tau = 0.6);
    
    void initialize(double rho0 = 1.0, double ux = 0.0, double uy = 0.0);
    void collide();
    void stream();
    void step();
    
    void apply_bounce_back(int x, int y);
    void apply_inlet(int x, int y, double ux, double uy);
    void apply_outlet(int x, int y);
    
    double density(int x, int y) const;
    std::array<double, 2> velocity(int x, int y) const;
    double vorticity(int x, int y) const;
    
    void add_force(int x, int y, double fx, double fy);
    
    // Accessors
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    double tau() const { return tau_; }
    double omega() const { return omega_; }
    
    const std::vector<double>& distributions() const { return f_; }
    
private:
    int nx_, ny_;
    double tau_, omega_;
    std::vector<double> f_;      // current distributions [nx*ny*Q]
    std::vector<double> f_tmp_;  // temporary for streaming
    std::vector<double> force_x_, force_y_;
    
    int idx(int x, int y, int q) const { return (y * nx_ + x) * Q + q; }
    double equilibrium(double rho, double ux, double uy, int q) const;
};

} // namespace continuumx
