# Model Predictive Control for Autonomous Vehicles
In this project, we will use [*model predictive control* (MPC)](https://en.wikipedia.org/wiki/Model_predictive_control) to steer a car around a track in `Project 5: MPC Controller` of [this simulator](https://github.com/udacity/self-driving-car-sim/releases). Specifically, we will develop a model predictive controller in C++ which uses [IPOPT](https://projects.coin-or.org/Ipopt) to solve the MPC optimization problem, and then uses [micro WebSockets (µWS)](https://github.com/uNetworking/uWebSockets) to communicate with the simulator. 

## Problem Setup

At each iteration, the simulator sends us a message containing the following data expressed in the simulator's inertial reference frame (not the car's reference frame):

1. A vector of `(x,y)` coordinates for some upcoming waypoints that we should try to pass through.
2. The current `(x,y)` coordinates of the car. 
3. The orientation `psi` of the car. 
4. The velocity `v` of the car. 

After sending us a message, the simulator waits for us to return a steering angle and throttle to use during the next iteration. 

One of the difficulties in this project is that we will be artificially introducing a time delay for the control signal. Specifically, the main function of our program looks like this:

```
int main(int argc, char *argv[]) {

    /* Set up the websockets */
    ...   
    
    /* Create the model predictive controller */
    MPC mpc(...);
    
    /* Handle a new message from the simulator */
    h.onMessage([&mpc](uWS::WebSocket <uWS::SERVER> ws,
                       char *data,
                       size_t length,
                       uWS::OpCode opCode) {

		/* Get the waypoints, position, orientation, and velocity */
		waypoints = ...

	    /* Use the model predictive controller to compute the control values */
	    mpc.solve(...);
        
        /* Put this thread to sleep to simulate a time delay */
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        
        /* Now send our control commands */
        send( mpc.steering_angle, mpc.throttle );
    });
    
    return 0;
}
```

Stably controlling the car despite this time delay is the primary challenge of this project.

## The Model  

As the name suggests, a model predictive controller needs a model of the system, which it will use to predict how the system will respond to hypothetical inputs given the current state of the system. The model predictive control then defines a cost function over some finite horizon of future inputs, and chooses the control which it believes will minimize the cost function over that horizon.  

In our case, we will use a very simple model of the car's dynamics, where the state vector contains the following values in the car's coordinate system:

1. The current `(x,y)` coordinates of the car. 
2. The orientation `psi` of the car. 
3. The velocity `v` of the car. 

Specifically,
$$
\dfrac{d}{dt}\left[ \begin{array}{ccc} x(t) \\ y(t) \\ \psi(t) \\ v(t) \end{array} \right] = 
\left[ \begin{array}{ccc} v(t) \cos(\psi(t)) \\ v(t) \sin(\psi(t)) \\ 0 \\ 0 \end{array} \right] +
\left[ \begin{array}{ccc} 0 & 0 \\ 0 & 0  \\ v(t) / L_f & 0 \\ 0 & 1 \end{array} \right] 
\left[ \begin{array}{ccc} \texttt{steering_angle}(t)  \\ \texttt{throttle}(t) \end{array} \right]
$$
where $L_f$ is an empirically-determined constant related to the turning radius of the car. For instance, in our simulations, we will use a value of $L_f = 2.67$.

Actually, since we are implementing a digital controller, we need the discrete-time form of these equations. Hence we use the following approximation
$$
\left[ \begin{array}{ccc} x(t + h) \\ y(t+ h) \\ \psi(t+ h) \\ v(t+ h) \end{array} \right] = 
\left[ \begin{array}{ccc} x(t) + h v(t) \cos(\psi(t)) \\ y(t) + h v(t) \sin(\psi(t)) \\ \psi(t) \\ v(t) \end{array} \right] +
\left[ \begin{array}{ccc} 0 & 0 \\ 0 & 0  \\ h v(t) / L_f & 0 \\ 0 & h \end{array} \right] 
\left[ \begin{array}{ccc} \texttt{steering_angle}(t)  \\ \texttt{throttle}(t) \end{array} \right]
$$
where $h$ denotes the sampling period.

## MPC Cost Function

Cost function design is more of an art than a science. In our case, we want a cost function that balances the following objectives:

1. Stay as close as possible to the waypoints.
2. Go as fast as possible around the track. 
3. Minimize the change in steering angle from timestep to timestep. 

The first two points should be obvious. While you might understand why the third point is important, you probably are underestimating its significance... 

Due to the fact that we might have a large time delay in our problem, it is vitally import that the controller does not make any sharp movements. This is important because if there is a time delay, but the control signal is smooth (in a qualitative, not mathematical sense), then we won't really feel the effect of that time delay because the control that we want to use at the current time step is very similar to the one we used at the previous time step. This is the core principle we will use to overcome time-delay errors. 

As for how we will achieve objective 1 - well, in that regard, we will simply fit a $3^{rd}$ order polynomial $f(x)$ to the waypoints, and then try to minimize the predicted crosstrack error, that is, the difference between $f(x)$ and predicted $y$-positions of the car over some finite horizon.

So let's get to to it... To compute the cost function, we will simulate the system $N$ time steps into the future with sample period $h$. Specifically, the model predictive controller internally has an estimate $d$ of the time delay. When the simulator gives the MPC an initial state $\texttt{state}(t-d)$, the MPC propagates that state to time $t$, after which it predicts the values $\texttt{state}(t+h), \ldots, \texttt{state}(t+Nh)$ using the above-mentioned model. Then we define the cost function $J$ to be 
$$
J =  \left( \sum_{k=1}^N \Big[  f(x(t_k)) - y(t_k) \Big]^2 + \dfrac{\alpha}{\max( v(t_k), 0.01 )^2}  \right) + \left( \sum_{k=1}^{N-1} \Big[ u(t_{k+1}) - u(t_k) \Big]^2 v^2(t_k) \right)
$$
where $t_k = t + kh$. The first term ensures that the car will stay close to the waypoints, the second term penalizes the car for going slow, and the third term ensures that the result is smooth. 

There are two things I want to draw your attention to here:

1. The coefficient $\alpha$ is a design parameter which controls how much we will penalize the car for going slow. In the simulator, the default value is 1, which means "cautious driver". I have tried values up to 1000, which means "race car driver". This almost always works out well, and is really fun to watch :).
2. The third term, which smoothes out the control is weighted by the square of the velocity. This means that when going fast, we really want to make sure that the control is smooth. If we are going slow, then we don't really care. This multiplier is significant because we will pay more dearly for errors at high speeds, than at low speeds. You know that...

> We solve this optimization problem with IPOPT at each time step, and implement only the subsequent control value, despite the fact that we have predicted $N$ time steps into the future. 

## The Code

To compile the code, you should create a `build` folder in your repo. Then 

```
cd build
cmake .. 
make
```

Once compiled, you can run our controller for the command line with no arguments:

```
./mpc
```

However, if you want to tune this controller for a specific scenario, you can specify any number of the following arguments:

1. The `velocity_scale`, $\alpha$, was defined in the previous section. The default value is 1 ("cautious driver"), but you should try running this with 1000 to see how fast we can go:

   ```
   ./mpc 1000
   ```

2. The time delay assumed by the model predictive controller. Internally, the model predictive has some hard-coded estimate of the time delay. To account for time delay, we simply take the initial state provided by the simulator and propagate that state using our dynamic model into the future by the specified time delay. We have to do this because our control will only take effect after this delay, so the value of the state after this propagation represents the actual initial state for the system when the first control command arrives. To change this value, you simply specify a second argument. For instance, the default MPC code uses the same time delay 0.1s, as the amount of time that we put the code to sleep. If you want to see what happens when you assume that there is no time delay and a velocity scale of 30, then you could run our program with 

   ```
   ./mpc 30 0
   ```

3. The number of time steps $N$ to perform the optimization over. The default value is 10. The danger in modifying this is that when $N$ is too large, then computing the optimal solution may take too long. In this case, you are introducing your own new time delay into the problem due to the computation time. Even increasing this value to 15 seems to have a significant effect on the quality of the solution. If you want to specify $\alpha=55$, time delay of 0.02, and $N = 15$, you would use:

   ```
   ./mpc 55 0.02 15
   ```

4. The sample time step $h$ used in computing the prediction horizon. To compute the cost function, we predict what the states will be at times $t+h,\ldots,t+Nh$. You can change this value of $h$. The default value is 0.66. If you want to try 1-second sampling intervals, then you could run 

   ```
   ./mpc 55 0.02 15 1
   ```

With that being said, I would recommend that you only mess with the velocity scale and the time delay assumed by the MPC algorithm. The other values are highly tuned, and changing them does not appear to improve performance. 

## Results

All of the following results are shown on the "Fantastic" video setting. You should use that setting to replicate our results. 

### Cautious Driver, Correct Time Delay

First, let's see how the MPC does when we use all of the default values ($\alpha = 1$, and time delay = 0.1s):

<video controls="controls" width=100%>
  <source type="video/mp4" src="videos/mpc_1.mp4"></source>
</video>

As you can see, this controller slows down significantly at the corners, but stays very close to the center of the lane throughout the simulation.

### Race Car Driver, Correct Time Delay 

If you want to go faster, then you can increase the `velocity_scale`, $\alpha$. In the following video, we increase $\alpha$ to 1000, and use the default time delay (0.1s).

<video controls="controls" width=100%>
  <source type="video/mp4" src="videos/mpc_1000.mp4"></source>
</video>

Now you can see the car really flying around the track, and staying pretty close to the center line while doing so.

### Cautious Driver, Assuming No Time Delay

The code is set up so that there will be a time delay of at least 0.1s. Internally, the MPC has it's own estimate of this value. In the following video, we use the "cautious driver" setting ($\alpha = 1$), and we tell the MPC that there is no time delay (even though there definitely is). 

<video controls="controls" width=100%>
  <source type="video/mp4" src="videos/mpc_1_0.mp4"></source>
</video>

In this video, we don't see that much performance degradation. That is probably because the car is going too slow to really notice it. 

### Race Car Driver, Assuming No Time Delay

As a final test case, we will use the race car driver setting ($\alpha=1000$), and again tell the MPC that there is no time delay (even though there definitely is one). In the following video, we see that this introduces some oscillatory behavior into the controller, but doesn't completely destabilize it. 

<video controls="controls" width=100%>
  <source type="video/mp4" src="videos/mpc_1000_0.mp4"></source>
</video>



