#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

struct OrbitState {
    double t;
    double x;
    double y;
    double vx;
    double vy;
};

// forward declarations

OrbitState rhs(const OrbitState& state);

void write_history(const std::vector<OrbitState>& history);

OrbitState update_state(const OrbitState& state,
                        const double dt,
                        const OrbitState& state_derivs);

std::vector<OrbitState> integrate(const double a,
                                  const double tmax,
                                  const double dt_in);

std::vector<OrbitState> integrate_rk2(const double a,
                                      const double tmax,
                                      const double dt_in);

double error(const std::vector<OrbitState>& history);


const double GM = 4.0 * M_PI * M_PI;   // G * Mass in AU, year, solar mass units

OrbitState rhs(const OrbitState& state) {

    OrbitState dodt{};

    // dx/dt = vx; dy/dt = vy

    dodt.x = state.vx;
    dodt.y = state.vy;

    // d(vx)/dt = - GMx/r**3; d(vy)/dt = - GMy/r**3

    double r = std::sqrt(state.x * state.x + state.y * state.y);

    dodt.vx = - GM * state.x / std::pow(r, 3);
    dodt.vy = - GM * state.y / std::pow(r, 3);

    return dodt;

}

void write_history(const std::vector<OrbitState>& history) {

    std::ofstream outfile("orbit.dat");

    for (auto o : history) {
        outfile << std::setw(12) << o.t
                << std::setw(12) << o.x
                << std::setw(12) << o.y
                << std::setw(12) << o.vx
                << std::setw(12) << o.vy << std::endl;

    }

}

OrbitState update_state(const OrbitState& state, const double dt,
                        const OrbitState& state_derivs) {

    OrbitState state_new{};

    state_new.t = state.t + dt;
    state_new.x = state.x + dt * state_derivs.x;
    state_new.y = state.y + dt * state_derivs.y;
    state_new.vx = state.vx + dt * state_derivs.vx;
    state_new.vy = state.vy + dt * state_derivs.vy;

    return state_new;
}

std::vector<OrbitState> integrate(const double a,
                                  const double tmax, const double dt_in) {

    // how the history of the orbit

    std::vector<OrbitState> orbit_history{};

    // set initial conditions
    OrbitState state{};

    // assume circular orbit on the x-axis, counter-clockwise orbit

    state.t = 0.0;
    state.x = a;
    state.y = 0.0;
    state.vx = 0.0;
    state.vy = std::sqrt(GM / a);

    orbit_history.push_back(state);

    double dt = dt_in;

    // integration loop
    while (state.t < tmax) {

        if (state.t + dt > tmax) {
            dt = tmax - state.t;
        }

        // get the derivatives
        auto state_derivs = rhs(state);

        // update the state
        state = update_state(state, dt, state_derivs);

        orbit_history.push_back(state);
    }

    return orbit_history;

}

std::vector<OrbitState> integrate_rk2(const double a,
                                      const double tmax, const double dt_in) {

    // how the history of the orbit

    std::vector<OrbitState> orbit_history{};

    // set initial conditions
    OrbitState state{};

    // assume circular orbit on the x-axis, counter-clockwise orbit

    state.t = 0.0;
    state.x = a;
    state.y = 0.0;
    state.vx = 0.0;
    state.vy = std::sqrt(GM / a);

    orbit_history.push_back(state);

    double dt = dt_in;

    // integration loop
    while (state.t < tmax) {

        if (state.t + dt > tmax) {
            dt = tmax - state.t;
        }

        // get the derivatives
        auto state_derivs = rhs(state);

        // get the derivatives at the midpoint in time
        auto state_star = update_state(state, 0.5 * dt, state_derivs);
        state_derivs = rhs(state_star);

        // update the state
        state = update_state(state, dt, state_derivs);

        orbit_history.push_back(state);
    }

    return orbit_history;

}

double error(const std::vector<OrbitState>& history) {

    const OrbitState& first = history[0];
    const OrbitState& last = history.back();

    double r_init = std::sqrt(first.x * first.x +
                              first.y * first.y);
    double r_final = std::sqrt(last.x * last.x +
                               last.y * last.y);

    return std::abs(r_final - r_init);
}

int main() {

    double tmax = 1.0;
    double a = 1.0;      // 1 AU

    // convergence comparing Euler and RK2

    for (double dt : {0.1, 0.05, 0.025, 0.0125}) {
        auto orbit_history = integrate(a, tmax, dt);
        auto orbit_history2 = integrate_rk2(a, tmax, dt);
        std::cout << std::setw(10) << dt << " "
                  << std::setw(10) << error(orbit_history) << " "
                  << std::setw(10) << error(orbit_history2) << std::endl;
    };


    // demonstration of writing
    double dt = 0.01;

    auto orbit_history = integrate_rk2(a, tmax, dt);
    write_history(orbit_history);

}
