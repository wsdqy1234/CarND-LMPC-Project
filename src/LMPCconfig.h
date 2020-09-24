//
// Created by lin on 2020/7/15.
//

#ifndef LMPCCONFIG_H
#define LMPCCONFIG_H

#include <cstdio>

using namespace std;

class Config {
public:
    const size_t N = 10;
    const double dt = 0.1;
    const size_t num_local_SS_point = 30;
    const int num_ite_build_local_SS = 3;//用
    const size_t extend_constraints = 3;
    const double road_half_width = 2.5;
    const size_t x_start = 0;
    const size_t y_start = x_start + N;
    const size_t psi_start = y_start + N;
    const size_t v_start = psi_start + N;
    const size_t cte_start = v_start + N;
    const size_t epsi_start = cte_start + N;
    const size_t delta_start = epsi_start + N;
    const size_t a_start = delta_start + N - 1;
    const size_t lambda_start = a_start + N - 1;
    const double Lf = 2.67;
    double ref_v = 70;
};

#endif //LMPCCONFIG_H
