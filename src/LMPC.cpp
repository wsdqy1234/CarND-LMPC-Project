//
// Created by lin on 2020/7/14.
//

#include "LMPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "./Eigen-3.3/Eigen/Core"


LMPC::LMPC() {}

LMPC::~LMPC() {}

using CppAD::AD;
Config config;

class FG_eval1 {
public:
    Eigen::VectorXd coeffs;
    vector<State> local_SS;
    vector<size_t> local_q_function;
    FG_eval1(Eigen::VectorXd coeffs, vector<State> local_SS, vector<size_t> local_q_function) {
        this->coeffs = coeffs;
        this->local_SS = local_SS;
        this->local_q_function = local_q_function;
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    void operator()(ADvector& fg, const ADvector& vars) {
        // LMPC uses IPOPT to achieve solution

        // Set the constraint on the variable
        // The constraint at the initial moment
        fg[1 + config.x_start] = vars[config.x_start];
        fg[1 + config.y_start] = vars[config.y_start];
        fg[1 + config.psi_start] = vars[config.psi_start];
        fg[1 + config.v_start] = vars[config.v_start];
        fg[1 + config.cte_start] = vars[config.cte_start];
        fg[1 + config.epsi_start] = vars[config.epsi_start];
        // The constraints in the rolling optimization process
        for (size_t t = 1; t < config.N; t++) { // N is a constant value representing the roll optimazation step.
            AD<double> x1 = vars[config.x_start + t];
            AD<double> x0 = vars[config.x_start + t - 1];
            AD<double> y1 = vars[config.y_start + t];
            AD<double> y0 = vars[config.y_start + t - 1];
            AD<double> psi1 = vars[config.psi_start + t];
            AD<double> psi0 = vars[config.psi_start + t - 1];
            AD<double> v1 = vars[config.v_start + t];
            AD<double> v0 = vars[config.v_start + t - 1];
            AD<double> cte1 = vars[config.cte_start + t];
            AD<double> cte0 = vars[config.cte_start + t - 1];
            AD<double> epsi1 = vars[config.epsi_start + t];
            AD<double> epsi0 = vars[config.epsi_start + t - 1];
            AD<double> a = vars[config.a_start + t - 1];
            AD<double> delta = vars[config.delta_start + t - 1];
            if (t > 1) {   // use previous actuations (to account for latency)
                a = vars[config.a_start + t - 2];
                delta = vars[config.delta_start + t - 2];
            }

            //求t-1时刻的变量
            //f0为t-1时刻的真值y，psides0为t-1时刻的真值psi
            AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
            AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));

            // Dynamic constraints on the model
            fg[1 + config.x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * config.dt);
            fg[1 + config.y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * config.dt);
            fg[1 + config.psi_start + t] = psi1 - (psi0 - v0/config.Lf * delta * config.dt);
            fg[1 + config.v_start + t] = v1 - (v0 + a * config.dt);
            fg[1 + config.cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * config.dt));
            fg[1 + config.epsi_start + t] = epsi1 - ((psi0 - psides0) - v0/config.Lf * delta * config.dt);
        }
        AD<double> xN = vars[config.x_start + config.N - 1];
        AD<double> yN = vars[config.y_start + config.N - 1];
        AD<double> psiN = vars[config.psi_start + config.N - 1];
        AD<double> vN = vars[config.v_start + config.N - 1];
        AD<double> xN_bar = 0;
        AD<double> yN_bar = 0;
        AD<double> psiN_bar = 0;
        AD<double> vN_bar = 0;
        AD<double> lambda_sum = 0;

        // Definition and calculation of Cost
        fg[0] = 0;
        size_t SS_len = local_SS.size();
        size_t n_vars = vars.size();
        size_t slack_start = config.lambda_start + SS_len;
        for (size_t i = 0; i < SS_len; i++) {
            lambda_sum = lambda_sum + vars[config.lambda_start+i];
            xN_bar = xN_bar + vars[config.lambda_start+i]*local_SS.at(i).x;
            yN_bar = yN_bar + vars[config.lambda_start+i]*local_SS.at(i).y;
            psiN_bar = psiN_bar + vars[config.lambda_start+i]*local_SS.at(i).psi;
            vN_bar = vN_bar + vars[config.lambda_start+i]*local_SS.at(i).v;
            fg[0] += vars[config.lambda_start+i]*local_q_function.at(i);
        }
        for (size_t i = slack_start; i < n_vars; i++) {
            fg[0] += 1000 * vars[i] * vars[i];
        }
        fg[1 + config.epsi_start + config.N] = lambda_sum - 1;
        fg[1 + config.epsi_start + config.N + 1] = xN - xN_bar + vars[slack_start];
        fg[1 + config.epsi_start + config.N + 2] = yN - yN_bar + vars[slack_start + 1];
        std::cout << "fg.size() = " << fg.size() << std::endl;
    }
};

vector<double> LMPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, double s, double origin_x, double origin_y, double origin_psi) {
    bool ok = true;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];

    UpdateCoeffs(coeffs);
    UpdateLocal(s, origin_x, origin_y, origin_psi);// Update the local collection
    

    // Set the number of model variables (including 6 observation variables and 2 execution variables)
    size_t n_lambda = local_SS_.size();
    size_t n_slack = config_.extend_constraints - 1;// Add slack variables
    size_t n_vars = config.N * 6 + (config.N - 1) * 2 + n_lambda + n_slack;

    size_t n_constraints = config.N * 6 + config_.extend_constraints;// Set the number of constraints for observation variables
    std::cout << "local_SS_.size() = " << local_SS_.size() << ", n_vars = " << n_vars << ", n_constraints = " << n_constraints << std::endl;
    // Initialization of model variables (all time set to 0)
    Dvector vars(n_vars);
    for (size_t i = 0; i < n_vars; i++) {
        vars[i] = 0;
    }

    // The upper and lower bounds of each variable
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);


    // States setting of model variables at time 0 (load input state)
    vars[config.x_start] = x;
    vars[config.y_start] = y;
    vars[config.psi_start] = psi;
    vars[config.v_start] = v;
    vars[config.cte_start] = cte;
    vars[config.epsi_start] = epsi;

    // Set the upper and lower bounds of all observed variables to the maximum range
    for (size_t i = 0; i < config.cte_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }
    for (size_t i = config.cte_start; i < config.epsi_start; i++) {
        vars_lowerbound[i] = -config_.road_half_width;
        vars_upperbound[i] = config_.road_half_width;
    }
    for (size_t i = config.epsi_start; i < config.delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    //Set the upper and lower bounds of all execution variables
    for (size_t i = config.delta_start; i < config.a_start; i++) {
        vars_lowerbound[i] = -0.436332;//-25<delta<25
        vars_upperbound[i] = 0.436332;
    }

    // Acceleration/decceleration upper and lower limits.
    for (size_t i = config.a_start; i < config.lambda_start; i++) {
        vars_lowerbound[i] = -1.0;//-1<a<1
        vars_upperbound[i] = 1.0;
    }
    for (size_t i = config_.lambda_start; i < n_vars; i++) {
        vars_lowerbound[i] = 0;
        vars_upperbound[i] = 1;
    }

    // The upper and lower limits of the observed variable constraints
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    //The upper and lower bounds of the initial state are the ground truth of the observed variables
    constraints_lowerbound[config.x_start] = x;
    constraints_lowerbound[config.y_start] = y;
    constraints_lowerbound[config.psi_start] = psi;
    constraints_lowerbound[config.v_start] = v;
    constraints_lowerbound[config.cte_start] = cte;
    constraints_lowerbound[config.epsi_start] = epsi;

    constraints_upperbound[config.x_start] = x;
    constraints_upperbound[config.y_start] = y;
    constraints_upperbound[config.psi_start] = psi;
    constraints_upperbound[config.v_start] = v;
    constraints_upperbound[config.cte_start] = cte;
    constraints_upperbound[config.epsi_start] = epsi;

    FG_eval1 fg_eval(coeffs, GetLocalSS(), GetLocalQFunc());// Generate IPOPT object


    // IPOPT is used to solve LMPC nonlinear optimization problem
    std::string options;
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval1>(
            options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
            constraints_upperbound, fg_eval, solution);

    // Check some of the solution values
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    // Returns the calculated execution variable at time 1, and the predicted (x,y)
    vector<double> result;
    result.push_back(solution.x[config.delta_start]);
    result.push_back(solution.x[config.a_start]);

    for (size_t i = 0; i < config.N-1; i++) {
        result.push_back(solution.x[config.x_start + i + 1]);
        result.push_back(solution.x[config.y_start + i + 1]);
    }
    std::cout << "LMPC:\n";
    std::cout << "status == success: " << (solution.status == CppAD::ipopt::solve_result<Dvector>::success) << std::endl;
    std::cout << "status == local_infeasibility: " << (solution.status == CppAD::ipopt::solve_result<Dvector>::local_infeasibility) << std::endl;
    std::cout << "result.size() = " << result.size() << std::endl;
    if (!result.empty()) {
        std::cout << "cost = " << solution.obj_value << std::endl;
        std::cout << "slack[0] = " << solution.x[solution.x.size()-2] << ", slack[1] = " << solution.x[solution.x.size()-1] << std::endl;
    }
    return result;
}

void LMPC::CollectX(State x) {
    x_collect_.push_back(x);
}

void LMPC::CollectU(Input u) {
    u_collect_.push_back(u);
}

void LMPC::AddTraj() {
    if (!x_collect_.empty()) {
        ite_index_.push_back(x_safe_set_.size());
        x_safe_set_.insert(x_safe_set_.end(), x_collect_.begin(), x_collect_.end());
        u_safe_set_.insert(u_safe_set_.end(), u_collect_.begin(), u_collect_.end());
        size_t n = x_collect_.size();
        for (size_t i = 0; i < n; i++) {
            //Here Qfun is defined as the thought of recursive strategy
            //the number of states in the last lap minus the number of states that have passed through the current iteration
            // i.e., the number of states in two adjacent laps is assumed to be almost constant
            q_function_.push_back(n-i-1);
        }
        ite_++;
    }
}

void LMPC::EraseCollect() {
    vector<State>().swap(x_collect_);
    vector<Input>().swap(u_collect_);
}

void LMPC::UpdateLocal(double s, double px, double py, double psi) {
    //Use Local SS to optimize, and here's an update to Local SS
    vector<State>().swap(local_SS_);
    vector<size_t>().swap(local_q_function_);
    size_t len = x_safe_set_.size();
    size_t count = local_SS_.size();
    int min_ite = max(ite_+1-config.num_ite_build_local_SS,0);
    for(size_t i = ite_index_.at(min_ite); i < len; i++) {
        if (count < config_.num_local_SS_point) {
            if (x_safe_set_.at(i).s > s) {
                local_SS_.push_back(x_safe_set_.at(i));
                local_q_function_.push_back(q_function_.at(i));
                double dx = local_SS_.at(local_SS_.size()-1).x - px;
                double dy = local_SS_.at(local_SS_.size()-1).y - py;
                local_SS_.at(local_SS_.size()-1).x = dx * cos(-psi) - dy * sin(-psi);
                local_SS_.at(local_SS_.size()-1).y = dx * sin(-psi) + dy * cos(-psi);
                count++;
            }
        } else break;
    }

    ofstream SS_out;
    SS_out.open("SS.csv",ios::app);
    SS_out<<setiosflags(ios::fixed)<<setprecision(5);
    SS_out<<local_SS_.at(0).s<<", "<<local_q_function_.at(0)<<", "<<local_SS_.at(config_.num_local_SS_point-1).s<<", "<<local_q_function_.at(config_.num_local_SS_point-1)<<std::endl;
    SS_out.close();
}

void LMPC::UpdateCoeffs(Eigen::VectorXd coeffs) {
    coeffs_ = coeffs;
}
