#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
private:
    double Lf = 2.67;
    double ref_v = 70;
    double dt = 0.1;
public:
    size_t N = 10;
    MPC();
    MPC(size_t N_now);
    virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
    vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
