#include "continuumx/lattice_boltzmann.hpp"
#include <cstring>

namespace continuumx {

LatticeBoltzmannD2Q9::LatticeBoltzmannD2Q9(int nx, int ny, double tau)
    : nx_(nx), ny_(ny), tau_(tau), omega_(1.0 / tau)
{
    if (nx <= 0 || ny <= 0) throw std::invalid_argument("Invalid grid size");
    if (tau <= 0.5) throw std::invalid_argument("tau must be > 0.5 for stability");
    f_.resize(nx * ny * Q, 0.0);
    f_tmp_.resize(nx * ny * Q, 0.0);
    force_x_.resize(nx * ny, 0.0);
    force_y_.resize(nx * ny, 0.0);
}

void LatticeBoltzmannD2Q9::initialize(double rho0, double ux, double uy) {
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            for (int q = 0; q < Q; ++q) {
                f_[idx(x, y, q)] = equilibrium(rho0, ux, uy, q);
            }
        }
    }
}

double LatticeBoltzmannD2Q9::equilibrium(double rho, double ux, double uy, int q) const {
    double cu = cx[q] * ux + cy[q] * uy;
    double uu = ux * ux + uy * uy;
    return w[q] * rho * (1.0 + cu / cs2 + (cu * cu) / (2.0 * cs2 * cs2) - uu / (2.0 * cs2));
}

void LatticeBoltzmannD2Q9::collide() {
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            double rho = density(x, y);
            auto [ux, uy] = velocity(x, y);
            // Simple Guo forcing
            ux += 0.5 * force_x_[y * nx_ + x] / rho;
            uy += 0.5 * force_y_[y * nx_ + x] / rho;
            for (int q = 0; q < Q; ++q) {
                double feq = equilibrium(rho, ux, uy, q);
                f_[idx(x, y, q)] = f_[idx(x, y, q)] - omega_ * (f_[idx(x, y, q)] - feq);
            }
        }
    }
    // Reset forces after use
    std::fill(force_x_.begin(), force_x_.end(), 0.0);
    std::fill(force_y_.begin(), force_y_.end(), 0.0);
}

void LatticeBoltzmannD2Q9::stream() {
    std::fill(f_tmp_.begin(), f_tmp_.end(), 0.0);
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            for (int q = 0; q < Q; ++q) {
                int x_new = x + cx[q];
                int y_new = y + cy[q];
                // Periodic for now (basic); walls handled separately
                if (x_new >= 0 && x_new < nx_ && y_new >= 0 && y_new < ny_) {
                    f_tmp_[idx(x_new, y_new, q)] = f_[idx(x, y, q)];
                }
            }
        }
    }
    f_.swap(f_tmp_);
}

void LatticeBoltzmannD2Q9::step() {
    collide();
    stream();
}

void LatticeBoltzmannD2Q9::apply_bounce_back(int x, int y) {
    // Simple bounce-back: reverse directions
    static const std::array<int, Q> opp = {0, 3, 4, 1, 2, 7, 8, 5, 6};
    for (int q = 1; q < Q; ++q) {
        f_[idx(x, y, opp[q])] = f_[idx(x, y, q)];
    }
}

void LatticeBoltzmannD2Q9::apply_inlet(int x, int y, double ux, double uy) {
    double rho = density(x, y);
    for (int q = 0; q < Q; ++q) {
        f_[idx(x, y, q)] = equilibrium(rho, ux, uy, q);
    }
}

void LatticeBoltzmannD2Q9::apply_outlet(int x, int y) {
    // Zero-gradient approximation
    if (x > 0) {
        for (int q = 0; q < Q; ++q) {
            f_[idx(x, y, q)] = f_[idx(x - 1, y, q)];
        }
    }
}

double LatticeBoltzmannD2Q9::density(int x, int y) const {
    double rho = 0.0;
    for (int q = 0; q < Q; ++q) rho += f_[idx(x, y, q)];
    return rho;
}

std::array<double, 2> LatticeBoltzmannD2Q9::velocity(int x, int y) const {
    double rho = density(x, y);
    if (rho < 1e-10) return {0.0, 0.0};
    double ux = 0.0, uy = 0.0;
    for (int q = 0; q < Q; ++q) {
        ux += f_[idx(x, y, q)] * cx[q];
        uy += f_[idx(x, y, q)] * cy[q];
    }
    return {ux / rho, uy / rho};
}

double LatticeBoltzmannD2Q9::vorticity(int x, int y) const {
    if (x <= 0 || x >= nx_-1 || y <= 0 || y >= ny_-1) return 0.0;
    auto [ux_yp, uy_yp] = velocity(x, y+1);
    auto [ux_ym, uy_ym] = velocity(x, y-1);
    auto [ux_xp, uy_xp] = velocity(x+1, y);
    auto [ux_xm, uy_xm] = velocity(x-1, y);
    return 0.5 * ((uy_xp - uy_xm) - (ux_yp - ux_ym));
}

void LatticeBoltzmannD2Q9::add_force(int x, int y, double fx, double fy) {
    if (x >= 0 && x < nx_ && y >= 0 && y < ny_) {
        force_x_[y * nx_ + x] += fx;
        force_y_[y * nx_ + x] += fy;
    }
}

} // namespace continuumx
