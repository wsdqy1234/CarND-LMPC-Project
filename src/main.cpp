#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "LMPC.h"
#include "json.hpp"
#define MAX_LAP_NUM 8

using json = nlohmann::json;
// The Angle and radian conversion function
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

//Check the Socket Event for JSON-type data
string hasData(string s) {
	auto found_null = s.find("null");
	auto b1 = s.find_first_of("[");
	auto b2 = s.rfind("}]");
	if (found_null != string::npos) {
		return "";
	} else if (b1 != string::npos && b2 != string::npos) {
		return s.substr(b1, b2 - b1 + 2);
	}
	return "";
}

// Evaluate a polynomial.
// result = coeffs[0]*x^0+...+coeffs[size-1]*x^(size-1)
double polyeval(Eigen::VectorXd coeffs, double x) {
	double result = 0.0;
	for (int i = 0; i < coeffs.size(); i++) {
		result += coeffs[i] * pow(x, i);
	}
	return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
												int order) {
	assert(xvals.size() == yvals.size());
	assert(order >= 1 && order <= xvals.size() - 1);
	Eigen::MatrixXd A(xvals.size(), order + 1);

	for (int i = 0; i < xvals.size(); i++) {
		A(i, 0) = 1.0;
	}

	for (int j = 0; j < xvals.size(); j++) {
		for (int i = 0; i < order; i++) {
			A(j, i + 1) = A(j, i) * xvals(j);
		}
	}

	auto Q = A.householderQr();
	auto result = Q.solve(yvals);
	return result;
}

bool checkIsGoal(double px_old, double py_old, double px, double py){
	double goalx = -40.62;
	double goaly = 108.73;
	bool flag;
    flag = (px_old > goalx && px < goalx) || (py_old > goaly && py < goaly);
	return flag;
}



int main() {
	uWS::Hub h;

	// Initialize the LMPC controller
	LMPC lmpc;
	h.onMessage([&lmpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
										 uWS::OpCode opCode) {
        string sdata = string(data).substr(0, length);
		static MPC mpc;
		static int lap = 0;
		static int feasible_lap = 3;
		static size_t N_initial = 7; //Rolling optimization step N during initial iteration
		static double s_total = 0;
		static double px_old = -40.62;
		static double py_old = 108.73;
//		static double v_old = 0;
		static double s_old = 0;
		static chrono::steady_clock::time_point lap_time_start = chrono::steady_clock::now();
		static chrono::steady_clock::time_point lap_time_end = chrono::steady_clock::now();
		static chrono::steady_clock::time_point pre_message_time = chrono::steady_clock::now();
		static chrono::steady_clock::time_point current_message_time = chrono::steady_clock::now();

        // A "42" at the beginning of a message represents a new message from a Websocket
		if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2' && lap <= MAX_LAP_NUM) {//If the message is valid
		    string s = hasData(sdata);
			if (s != "") {
			    auto j = json::parse(s);
				string event = j[0].get<string>();
				if (event == "telemetry") {// j[1] is the data JSON object
					vector<double> ptsx = j[1]["ptsx"];
					vector<double> ptsy = j[1]["ptsy"];
					double px = j[1]["x"];
					double py = j[1]["y"];
					double psi = j[1]["psi"];
					double v = j[1]["speed"];
					double steer_angle = j[1]["steering_angle"];	// steering angle is in the opposite direction
					double acceleration = j[1]["throttle"];
                    current_message_time = chrono::steady_clock::now();
                    double delta_t = chrono::duration_cast<chrono::duration<double>>(current_message_time-pre_message_time).count();
                    pre_message_time = current_message_time;

					// Output lap time
					if (checkIsGoal(px_old, py_old, px, py)) { //If the endpoint is reached, Lap Time is calculated and the iteration information is output
						lap_time_end = chrono::steady_clock::now();
						chrono::duration<double> lap_time_used = chrono::duration_cast<chrono::duration<double>>(lap_time_end-lap_time_start);
						if (lap_time_used.count() > 1){
							lap++;
							s_total = s_old;
							s_old = s_old - s_total;
							lmpc.AddTraj();//Add all the state vectors and input vectors in this lap to the Safe Set
                            std::cout<<"lap = "<<lap<<", "<<"time use:"<<lap_time_used.count()<<"s"<<", "<<"s_total:"<<s_total<<std::endl;
							std::cout<<"x_collect_length: "<<lmpc.GetCollectLength()<<std::endl;
							std::cout<<"safe_set_length: "<<lmpc.GetSafeSetLength()<<std::endl;
							std::cout<<"u_safe_set_length: "<<lmpc.GetUSafeSetLength()<<std::endl;
							std::cout<<"q_function_length: "<<lmpc.GetQFunctionLength()<<std::endl;
							lmpc.EraseCollect();//Collect again
                            lap_time_start = chrono::steady_clock::now();
						}
					}
					double delta_s = sqrt(pow((px-px_old),2)+pow((py-py_old),2));
                    s_old += delta_s; // the same as approximating s with a straight line, as long as delta T goes to 0
                    px_old = px;
                    py_old = py;

                    // transform waypoints to be from car's perspective
                    // this means we can consider px = 0, py = 0, and psi = 0
                    // greatly simplifying future calculations
					vector<double> waypoints_x;
					vector<double> waypoints_y;
					for (unsigned int i = 0; i < ptsx.size(); i++) {
						double dx = ptsx[i] - px;
						double dy = ptsy[i] - py;
						waypoints_x.push_back(dx * cos(-psi) - dy * sin(-psi));
						waypoints_y.push_back(dx * sin(-psi) + dy * cos(-psi));
					}
					double* ptrx = &waypoints_x[0];
					double* ptry = &waypoints_y[0];
					Eigen::Map<Eigen::VectorXd> waypoints_x_eig(ptrx, 6);
					Eigen::Map<Eigen::VectorXd> waypoints_y_eig(ptry, 6);

					auto coeffs = polyfit(waypoints_x_eig, waypoints_y_eig, 3);
					double cte = polyeval(coeffs, 0);	// px = 0, py = 0
					double epsi = -atan(coeffs[1]);	// p

                    //Feasible_traj is feasible_TRAj for MPC and FEASJ for LMPC. Write the pre-iteration state variable to CSV file
                    if (lap <= feasible_lap){  //Feasible solution for accumulation feasible LAP per iteration of MPC
                        ofstream traj_out;
                        traj_out.open("feasible_traj.csv",ios::app);
                        traj_out<<setiosflags(ios::fixed)<<setprecision(5);
                        traj_out<<px<<", "<<py<<", "<<psi<<", "<<v<<", "<<cte<<", "<<epsi<<", "<<steer_angle<<", "<<acceleration<<"\n";
                        traj_out.close();
                    }
                    else{ //Start accumulating LMPC solutions
                        ofstream traj_out;
                        traj_out.open("traj.csv",ios::app);
                        traj_out<<setiosflags(ios::fixed)<<setprecision(5);
                        traj_out<<delta_t<<"; "<<s_old<<", "<<px<<", "<<py<<", "<<psi<<", "<<v<<", "<<cte<<", "<<epsi<<"; "<<steer_angle<<", "<<acceleration<<"\n";
                        traj_out.close();
                    }

					// Writes the error to a local file
					ofstream in;
					in.open("error.txt",ios::app);
					in<<cte<<", "<<epsi<<"\n";
					in.close();

					// Use MPC and LMPC respectively
					double steer_value = j[1]["steering_angle"];
					double throttle_value = j[1]["throttle"];

                    vector<double> mpc_x_vals;//MPC预测的轨迹
                    vector<double> mpc_y_vals;
                    vector<double> next_x_vals;//参考轨迹线
                    vector<double> next_y_vals;

                    Eigen::VectorXd state(6);
					state << 0, 0, 0, v, cte, epsi;
					if (lap <= feasible_lap) {
                        if(lap == 0){
                            mpc.N = N_initial;
                            auto vars = mpc.Solve(state, coeffs);
                            //Change the MPC's rolling optimization step size
                            steer_value = vars[0];
                            throttle_value = vars[1];
                            //Add the MPC predicted points to the List, which represents the predicted trajectory of the system
                            //That corresponds to the green line on the simulator
                            for (unsigned int i = 2; i < vars.size(); i ++) {
                                if (i%2 == 0) {
                                    mpc_x_vals.push_back(vars[i]);
                                }
                                else {
                                    mpc_y_vals.push_back(vars[i]);
                                }
                            }

                        }else{//In the first lap, MPC was used for solution optimization and feasible solutions were accumulated
                            mpc.N = N_initial+lap-1;
                            auto vars = mpc.Solve(state, coeffs);
                            steer_value = vars[0];
                            throttle_value = vars[1];
                            lmpc.CollectX(State(s_old, px, py, psi, v, cte, epsi));
                            lmpc.CollectU(Input(steer_value, throttle_value));
                            for (unsigned int i = 2; i < vars.size(); i ++) {
                                if (i%2 == 0) {
                                    mpc_x_vals.push_back(vars[i]);
                                }
                                else {
                                    mpc_y_vals.push_back(vars[i]);
                                }
                            }
                        }

					}
					else {
					    auto vars = lmpc.Solve(state, coeffs, s_old, px, py, psi);
					    steer_value = vars[0];
					    throttle_value = vars[1];
                        lmpc.CollectX(State(s_old, px, py, psi, v, cte, epsi));
                        lmpc.CollectU(Input(steer_value, throttle_value));
                        for (unsigned int i = 2; i < vars.size(); i ++) {
                            if (i%2 == 0) {
                                mpc_x_vals.push_back(vars[i]);
                            }
                            else {
                                mpc_y_vals.push_back(vars[i]);
                            }
                        }
					}

                    // Add the reference trajectory points to the List, which represents the system's reference trajectory
                    // Corresponds to the yellow line on the simulator
                    for (double i = 0; i < 100; i += 3){
                        next_x_vals.push_back(i);
                        next_y_vals.push_back(polyeval(coeffs, i));
                    }
					//Send steering Angle and throttle information to the simulator
					json msgJson;
					// NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
					// Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
					msgJson["steering_angle"] = steer_value/(deg2rad(25));
					msgJson["throttle"] = throttle_value;
                    msgJson["mpc_x"] = mpc_x_vals;
                    msgJson["mpc_y"] = mpc_y_vals;
					msgJson["next_x"] = next_x_vals;
					msgJson["next_y"] = next_y_vals;


					auto msg = "42[\"steer\"," + msgJson.dump() + "]";

                    //Delay module
                    // The aim is to simulate a real driving environment in which the car can indeed give immediate instructions
					this_thread::sleep_for(chrono::milliseconds(100));
					ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
				}
			} else {
				//Manual
				std::string msg = "42[\"manual\",{}]";
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
			}
		}
	});


	h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
										 size_t, size_t) {
		const std::string s = "<h1>Hello world!</h1>";
		if (req.getUrl().valueLength == 1) {
			res->end(s.data(), s.length());
		} else {
			// i guess this should be done more gracefully?
			res->end(nullptr, 0);
		}
	});

	h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
		std::cout << "Connected!!!" << std::endl;
	});

	h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
												 char *message, size_t length) {
		ws.close();
		std::cout << "Disconnected" << std::endl;
	});

	int port = 4567;
	if (h.listen(port)) {
		std::cout << "Listening to port " << port << std::endl;
	} else {
		std::cerr << "Failed to listen to port" << std::endl;
		return -1;
	}

	h.run();
}
