#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen/Core"
#include "Eigen/QR"
#include "mpc.h"
#include "json.hpp"
#include "base_type.h"

using std::cout;
using std::endl;
using std::string;

using json = nlohmann::json;

/**
 * Check if the SocketIO event has JSON data.
 * If there is data the JSON object in string format will be returned,
 * else the empty string "" will be returned.
 */
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

/** Evaluate a polynomial. */
T polyeval(const std::vector<T> &coeffs, T x) {
    T result = 0.0;
    for (size_t i = 0; i < coeffs.size(); i++)
        result += coeffs[i] * pow(x, i);
    return result;
}

/**
 * Fit a polynomial.
 * Adapted from
 * https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
 */
template<size_t order>
Eigen::Matrix<T, order + 1, 1> polyfit(std::vector<T> &xvals, std::vector<T> &yvals) {

    /* Make sure that we have the same number of x and y values */
    assert(xvals.size() == yvals.size());

    /* Make sure that the user specified a positive order, and that the order is less than the number of data points */
    assert(order >= 1 && order <= xvals.size() - 1);

    /* Now, let's fit */
    Eigen::Matrix<T, Eigen::Dynamic, order + 1> A(xvals.size(), order + 1);
    for (size_t i = 0; i < xvals.size(); i++)
        A(i, 0) = 1.0;
    for (size_t j = 0; j < xvals.size(); j++)
        for (size_t i = 0; i < order; i++)
            A(j, i + 1) = A(j, i) * xvals[j];

    auto Q = A.householderQr();
    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> yvals_(yvals.data(), yvals.size());
    auto result = Q.solve(yvals_);
    return result;
}

/** Convert a string to the specified datatype */
template<typename T>
T stringTo(const std::string &str) {
    std::istringstream ss(str);
    T num;
    ss >> num;
    return num;
}

int main(int argc, char *argv[]) {

    /* See if the user passed default values */
    const size_t N = argc > 1 ? stringTo<size_t>(argv[1]) : 10;
    const T dt = argc > 2 ? stringTo<T>(argv[2]) : .1;
    const T reference_velocity = argc > 3 ? stringTo<T>(argv[3]) : 40;

    uWS::Hub h;

    /* Create the model predictive controller */
    const T initial_time_delay_estimate = 0.1;
    MPC mpc(N, dt, reference_velocity, initial_time_delay_estimate);

    /* Handle a new message */
    h.onMessage([&mpc](uWS::WebSocket <uWS::SERVER> ws,
                       char *data,
                       size_t length,
                       uWS::OpCode opCode) {

        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        string sdata = string(data).substr(0, length);
        if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
            string s = hasData(sdata);
            if (s != "") {
                auto j = json::parse(s);
                string event = j[0].get<string>();
                if (event == "telemetry") {

                    // j[1] is the data JSON object
                    std::vector<double> ptsx = j[1]["ptsx"];
                    std::vector<double> ptsy = j[1]["ptsy"];
                    T px = j[1]["x"];
                    T py = j[1]["y"];
                    T psi = j[1]["psi"];
                    const T v = j[1]["speed"];
                    const T previous_angle = j[1]["steering_angle"];
                    const T previous_throttle = j[1]["throttle"];
                    const T global_psi = psi;

                    /* Translate the global coordinate system to the origin of the car's reference frame. */
                    for (size_t i = 0; i < ptsx.size(); ++i) {
                        ptsx[i] -= px;
                        ptsy[i] -= py;
                    }
                    px = 0;
                    py = 0;

                    /* Now rotate the coordinates */
                    for (size_t i = 0; i < ptsx.size(); ++i) {
                        const T x = ptsx[i];
                        const T y = ptsy[i];
                        ptsx[i] = +cos(psi) * x + sin(psi) * y;
                        ptsy[i] = -sin(psi) * x + cos(psi) * y;
                    }
                    psi = 0;

                    /* Fit a polynomial through the waypoints. This requires us to translate the coordinates from
                     * the global coordinate system to the local one. */
                    const auto coeffs = polyfit<3>(ptsx, ptsy);
                    const std::vector<T> reference_polynomial(coeffs.data(), coeffs.data() + coeffs.size());

                    /* Calculate steering angle and throttle using MPC.
                     * The simulator wants these values in the range [-1,1].
                     * After computing the solution, the MPC saves these normalized
                     * values in the fields normalized_xxx */
                    const std::vector<T> state = {px, py, psi, v};
                    const std::vector<T> previous_controls = {-previous_angle, previous_throttle};
                    mpc.solve(state, reference_polynomial, previous_controls, global_psi);

                    /* Now put all of the data in a message that we will pass back to the simulator. */
                    json msgJson;
                    msgJson["steering_angle"] = -mpc.normalized_steering_angle;
                    msgJson["throttle"] = mpc.normalized_throttle;

                    /* Display the MPC predicted trajectory as a green line */
                    msgJson["mpc_x"] = mpc.x_vals;
                    msgJson["mpc_y"] = mpc.y_vals;

                    /* Display the reference line as a yellow line */
                    std::vector<T> ptsy_fit(ptsx.size());
                    for (size_t i = 0; i < ptsx.size(); ++i)
                        ptsy_fit[i] = polyeval(reference_polynomial, ptsx[i]);
                    msgJson["next_x"] = ptsx;
                    msgJson["next_y"] = ptsy_fit;

                    auto msg = "42[\"steer\"," + msgJson.dump() + "]";
//                    cout << endl << endl << "Message : " << endl << msg << endl << endl;

                    // Latency
                    // The purpose is to mimic real driving conditions where
                    // the car does actuate the commands instantly.
                    //
                    // Feel free to play around with this value but should be to drive
                    // around the track with 100ms latency.
                    //
                    // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
                    // SUBMITTING.
                    // TODO
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
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

    h.onConnection([&h](uWS::WebSocket <uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket <uWS::SERVER> ws, int code,
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