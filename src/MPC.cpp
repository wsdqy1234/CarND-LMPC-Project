#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "./Eigen-3.3/Eigen/Core"
// TODO: Set the timestep length and duration
using CppAD::AD;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
/*const double Lf = 2.67;
double ref_v = 70;
size_t N = 10;
double dt = 0.1;*/

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  double Lf;
  double ref_v;
  size_t N;
  double dt;
  FG_eval(Eigen::VectorXd coeffs,double Lf, double ref_v, size_t N, double dt) {
      this->coeffs = coeffs;
      this->Lf = Lf;
      this->ref_v = ref_v;
      this->N = N;
      this->dt = dt;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    fg[0] = 0;
    int cte_penalty = 3000;
    int epsi_penalty = 3000;
    int v_penalty = 1;
    int delta_penalty = 5;
    int a_penalty = 5;
    int Ddelta_penalty = 200;
    int Da_penalty = 10;
      size_t x_start = 0;
      size_t y_start = x_start + N;
      size_t psi_start = y_start + N;
      size_t v_start = psi_start + N;
      size_t cte_start = v_start + N;
      size_t epsi_start = cte_start + N;
      size_t delta_start = epsi_start + N;
      size_t a_start = delta_start + N - 1;

    for (size_t i = 0; i < N; i++) {
      fg[0] += cte_penalty * CppAD::pow(vars[cte_start + i], 2);
      fg[0] += epsi_penalty * CppAD::pow(vars[epsi_start + i], 2);
      fg[0] += v_penalty * CppAD::pow(vars[v_start + i] - ref_v, 2);
    }
    for (size_t i = 0; i < N - 1; i++) {
      fg[0] += delta_penalty * CppAD::pow(vars[delta_start + i], 2);
      fg[0] += a_penalty * CppAD::pow(vars[a_start + i], 2);
      // try adding penalty for speed + steer
      fg[0] += 700*CppAD::pow(vars[delta_start + i] * vars[v_start+i], 2);
    }
    for (size_t i = 0; i < N - 2; i++) {
      fg[0] += Ddelta_penalty * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += Da_penalty * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (size_t t = 1; t < N; t++) {
      AD<double> x1 = vars[x_start + t];
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y1 = vars[y_start + t];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v1 = vars[v_start + t];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi1 = vars[epsi_start + t];
      AD<double> epsi0 = vars[epsi_start + t - 1];
      AD<double> a = vars[a_start + t - 1];
      AD<double> delta = vars[delta_start + t - 1];
      if (t > 1) {   // use previous actuations (to account for latency)
        a = vars[a_start + t - 2];
        delta = vars[delta_start + t - 2];
      }

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));

      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 - v0/Lf * delta * dt);
      fg[1 + v_start + t] = v1 - (v0 + a * dt);
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) - v0/Lf * delta * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::MPC(size_t N_now) {
    N = N_now;
}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;
    size_t x_start = 0;
    size_t y_start = x_start + N;
    size_t psi_start = y_start + N;
    size_t v_start = psi_start + N;
    size_t cte_start = v_start + N;
    size_t epsi_start = cte_start + N;
    size_t delta_start = epsi_start + N;
    size_t a_start = delta_start + N - 1;

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  size_t n_vars = N * 6 + (N - 1) * 2;
  size_t n_constraints = N * 6;

  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  for (size_t i = 0; i < cte_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  for (size_t i = cte_start; i < epsi_start; i++) {
    vars_lowerbound[i] = -2;
    vars_upperbound[i] = 2;
  }
  for (size_t i = epsi_start; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  for (size_t i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;//-25<delta<25
    vars_upperbound[i] = 0.436332;
  }


  for (size_t i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;//-1<a<1
    vars_upperbound[i] = 1.0;
  }


  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  FG_eval fg_eval(coeffs,Lf,ref_v,N,dt);


  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  options += "Numeric max_cpu_time          0.5\n";

  CppAD::ipopt::solve_result<Dvector> solution;

  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  vector<double> result;

  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  for (size_t i = 0; i < N-1; i++) {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }


  return result;
}
