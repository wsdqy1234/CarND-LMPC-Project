//
// Created by lin on 2020/7/14.
//

#ifndef FTOCP_H
#define FTOCP_H

#include "LMPCconfig.h"

class State {
public:
    double s;
    double x;
    double y;
    double psi;
    double v;
    double cte;
    double epsi;
    State(double S, double X, double Y, double Psi, double V, double Cte, double Epsi) : s(S), x(X), y(Y), psi(Psi), v(V), cte(Cte), epsi(Epsi) {}
    ~State() {}
};

class Input {
public:
    double delta;
    double acceleration;
    Input(double Delta, double Acceleration) : delta(Delta), acceleration(Acceleration) {}
    ~Input() {}
};

class FTOCP {
private:
    Config config_;
};

#endif //FTOCP_H
