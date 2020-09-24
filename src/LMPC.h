//
// Created by lin on 2020/7/14.
//

#ifndef LMPC_H
#define LMPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "FTOCP.h"

using namespace std;

class LMPC {
private:
    vector<State> x_safe_set_;//All collected State vectors, States
    vector<Input> u_safe_set_;//All the collected Execute vectors, Input
    vector<size_t> ite_index_;//Record the ending index of (j-1)-th LMPC, which is the location of Safe Set
    vector<State> x_collect_;//Collect the state vectors in the current lap
    vector<Input> u_collect_;//Collect the execute vectors in the current lap
    vector<size_t> q_function_;//Qfun of LMPC
    vector<State> local_SS_;//Local variable collection used to solve LMPC problem, local Safe Set
    vector<size_t> local_q_function_;// Local Qfun collection used to solve LMPC problem, local Qfun.
    Eigen::VectorXd coeffs_;//Polynomial coefficient of the reference path
    int ite_ = 0;//iteration index, equal to lap index
    FTOCP ftocp_;
    Config config_;
    void UpdateLocal(double s, double px, double py, double psi);//Update local variables
    void UpdateCoeffs(Eigen::VectorXd coeffs);
public:
    LMPC();
    virtual ~LMPC();
    void CollectX(State x);
    void CollectU(Input u);
    void AddTraj();//After one lap, x_collect_ and u_collect_ are added to safe_set
    void EraseCollect();//Empty this circle of collected data
    vector<State> GetLocalSS() {
        vector<State> local_SS(local_SS_);
        return local_SS;
    }
    vector<size_t> GetLocalQFunc() {
        vector<size_t> local_q_function(local_q_function_);
        return local_q_function;
    }
    Eigen::VectorXd GetCoeffs() {
        Eigen::VectorXd coeffs(coeffs_);
        return coeffs;
    }
    inline size_t GetCollectLength() { return x_collect_.size(); }
    inline size_t GetSafeSetLength() { return x_safe_set_.size(); }
    inline size_t GetUSafeSetLength() { return u_safe_set_.size(); }
    inline size_t GetQFunctionLength() { return q_function_.size(); }
    vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs, double s, double origin_x, double origin_y, double origin_psi);
};

#endif //LMPC_H
