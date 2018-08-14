# Diffusion process simulator

Simulates a variety of 1D and 2D diffusion processes for user chosen parameter values. See below for a description of the diffusion processes.  

To run this app through R, input the following commands into your R console

library(shiny)

runGitHub("diffusion_process_simulator", "jensenmelanie")


## Diffusion Process
While there is a variety of phenoma that can be modeled the following diffusion processes, the following processes are introduced in terms of modeling the position of a microparticle. We assume that the process is observed for time <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\large&space;t&space;\in&space;[0,T]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\large&space;t&space;\in&space;[0,T]" title="\large t \in [0,T]" /></a>

### Brownian Motion/Pure Diffusion 
A pure diffusion process satisfies the following Stochastic Differential Equation (SDE):

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;d&space;X_t&space;&&space;=&space;\sqrt{2D}dW_t&space;\\&space;X_0&space;&&space;=&space;0&space;\end{align*}" title="\begin{align*} d X_t & = \sqrt{2D}dW_t \\ X_0 & = 0 \end{align*}" />
</p>
  
where D > 0 is the diffusion constant and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\large&space;W_t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\large&space;W_t" title="\large W_t" /></a> is standard Brownian Motion. 

### Ornstein-Uhlenbeck (OU) process with time trend
An Ornstein-Uhlenbeck process with a time trend satisfies the following SDE:
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;d&space;X_t&space;&&space;=&space;-&space;\kappa&space;(X_t&space;-&space;\nu&space;t&space;)&space;dt&space;&plus;&space;\sqrt{2D}dW_t&space;\\&space;X_0&space;&&space;=&space;0&space;\end{align*}" title="\begin{align*} d X_t & = - \kappa (X_t - \nu t ) dt + \sqrt{2D}dW_t \\ X_0 & = 0 \end{align*}" />
</p>
where the model parameters have the following interpretation:

- <a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\kappa" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;\kappa" title="\large \kappa" /></a>: how "strongly" the particle reacts to the random perturbations 
- <a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;\nu" title="\large \nu" /></a>: the speed of the time trend
- D: the size of the noise due to thermal flucation 

By setting <a href="https://www.codecogs.com/eqnedit.php?latex=\large&space;\nu" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\large&space;\nu" title="\large \nu" /></a> = 0, then a centered (around 0) Ornstein-Uhlenbeck process can be simulated.

### Switching Ornstein Uhlenbeck with a time trend
Assume there are k switch points,<a html="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\large&space;\tau_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\large&space;\tau_i" title="\large \tau_i" /></a>, for i = 0, 1, ..., k+1, where <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\large&space;\tau_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\large&space;\tau_0" title="\large \tau_0" /></a>
= 0 and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\large&space;\tau_{k&plus;1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\large&space;\tau_{k&plus;1}" title="\large \tau_{k+1}" /></a> = T. 
We define a switching Ornstein Uhlenbeck with a time trend to satisfy the following SDE:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?\begin{align*}&space;d&space;X_t&space;=&space;-&space;\kappa_i\big(X_t&space;-&space;\nu_i&space;(t&space;-&space;\tau_{i-1})&space;&plus;&space;X_{\tau_{i-1}}\big)&space;&&space;dt&space;&plus;&space;\sqrt{2D_i}dW_t&space;\quad&space;\text{for&space;}&space;\tau_{i-1}&space;<&space;t&space;\leq&space;\tau_i\\&space;X_0&space;&&space;=&space;0&space;\end{align*}" title="\begin{align*} d X_t = - \kappa_i\big(X_t - \nu_i (t - \tau_{i-1}) + X_{\tau_{i-1}}\big) & dt + \sqrt{2D_i}dW_t \quad \text{for } \tau_{i-1} < t \leq \tau_i\\ X_0 & = 0 \end{align*}" />
</p>
where the parameters can take different values during each state.


## Restriction on Parameters

- D > 0 
- <img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{120}&space;\large&space;\kappa" title="\large \kappa" /> > 0
- <img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{120}&space;\large&space;\nu&space;\in&space;\mathbb{R}" title="\large \nu \in \mathbb{R}" />


### Switching OU process
For a given number of switch points k, the user can either specify 1 value or k + 1 values for each of the three parameters, D, <img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{120}&space;\large&space;\kappa" title="\large \kappa" /> , and <img src="https://latex.codecogs.com/gif.latex?\inline&space;\dpi{120}&space;\large&space;\nu" title="\large \nu" /> .
In the case when only 1 value is given, the process will be simulated assuming the given parameter value remains the constant for each state.


### 2D OU process and switchin OU process
An additional parameter is required for the simulations. The user needs to specify the angle of motion in radians. 


