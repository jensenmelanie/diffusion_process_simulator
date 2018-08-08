library(RColorBrewer)
mellow.color.pal = colorRampPalette(rev(brewer.pal(11,'Spectral')))
##Multiple Velocities
##=================================================
##1D processes 
##=================================================

##-------------------------------------------------
##BM
##-------------------------------------------------

BM_1D = function(time_seq, true_D){
	particle_data = list()
	
	Xn = c(0)
	dt_seq = c(0,diff(time_seq))

	for (n in 2:length(time_seq)){
		Xn[n] = Xn[n-1] + sqrt(2*dt_seq[n]*true_D)*rnorm(n = 1, 0, 1)		
	}
	
	particle_data$trajectory = cbind(time_seq, Xn)
	particle_data$parameters = matrix(c(true_D), nrow = 1, ncol = 1 )
	colnames(particle_data$parameters) = c("D")
	return(particle_data)

}

##-------------------------------------------------
##OU with drift
##-------------------------------------------------
OUdrift_1D = function(time_seq, true_D, true_kg, true_nu){

particle_data = list()

dt_seq = c(0,diff(time_seq))

OU_variance = (true_D/true_kg)*(1 - exp(-2*dt_seq*true_kg))

Xn = c(0); 


path_length = length(time_seq)



for(n in 2:(path_length) ){
		
	rho = exp(-true_kg*dt_seq); 
	An = (time_seq[n] - time_seq[n-1]*rho[n]) - (1/true_kg)*(1- rho[n])
	
	Xn[n] =  rho[n]*Xn[n-1] + true_nu*An +  sqrt(OU_variance[n])*rnorm(1, mean = 0, sd = 1)
	
	}

#pos_limits = range(Xn)
#plot(time_seq, Xn , type = "l", xlim = range(time_seq), ylim = pos_limits, xlab = "Time (seconds)", ylab = expression(paste("Position (", mu, "m)")))

particle_data$trajectory = cbind(time_seq, Xn)
particle_data$parameters = matrix(c(true_D, true_kg, true_nu), nrow = 1, ncol = 3 )
colnames(particle_data$parameters) = c("D", "kg", "nu" )
return(particle_data)

}


##-------------------------------------------------
##Switching OU with drift
##-------------------------------------------------

##Either the switch points can be randomized or chosen
## Can handle different parameter values for d, kg and nu


OUdrift_switching_1D = function(time_seq, switch_points = list("random",2), true_D, true_kg_vec, true_nu_vec){

particle_data = list()

num_steps = length(time_seq) -1
path_length = num_steps+1


if(class(switch_points) == "list"){
	num_switch = switch_points[[2]]
	possible_switch_pts = ceiling((path_length*0.1)):floor(path_length*0.9)
	switch_pts_actual = sort((sample(possible_switch_pts, num_switch)))
}else if (class(switch_points) == "numeric" || class(switch_points) == "integer"){
	num_switch = length(switch_points)
	switch_pts_actual = switch_points
}else{
	stop("Either a list of form: list('random', #switch points), or actual switch points")
}
	
switch_pts = c(0, switch_pts_actual, path_length)
	
my_colors = mellow.color.pal(num_switch+2)
if(sum(as.numeric(my_colors == "#FFFFBF")) >0 ){
	my_colors = my_colors[- which(my_colors == "#FFFFBF")]
}else{
	my_colors = my_colors[1:(num_switch+1)]
	}

## Setting up the state parameters

if(length(true_D) != (num_switch +1)){
	if( length(true_D) == 1){
		true_D = rep(true_D, times = (num_switch +1))
	}else{
		stop(paste("Dimension of initial D values should be one or ", num_switch +1))}
}


if(length(true_kg_vec) != (num_switch +1)){
	if( length(true_kg_vec) == 1){
		true_kg_vec = rep(true_kg_vec, times = (num_switch +1))
	}else{
		stop(paste("Dimension of initial kg values should be one or ", num_switch +1))}
}

if( length(true_nu_vec) != (num_switch +1)){
	if( length(true_nu_vec) == 1){
		true_nu_vec = rep(true_nu_vec, num_switch +1)
	}else{
		stop(paste("Dimension of initial velocities should be one or ", num_switch +1))
	}
}

state_vec = c(); switch_pt_vec = c(); 
k  = 0
for(tt in 1:(num_switch+1)){
k = k+1
state_index = (switch_pts[tt]+1):switch_pts[tt+1]
state_vec[state_index] =  k
switch_pt_vec[state_index] = switch_pts[k]
}

##Starting as OU with drift and OU -- in this case both processes have the same  variance but they have different centered terms

dt_vec = c(0,diff(time_seq))


Xn = c(0); 

for(n in 2:(path_length) ){

	switch_pt_index = max(switch_pt_vec[n],1)
	state_index = state_vec[n]
	
	OU_variance = (true_D[state_index]/true_kg_vec[state_index])*(1 - exp(-2*dt_vec[n]*true_kg_vec[state_index]))
	
	rho = exp(-true_kg_vec[state_index]*dt_vec[n]); 
	An = (time_seq[n] - time_seq[n-1]*rho) - (1/true_kg_vec[state_index])*(1- rho)
	OUdrift_mean_X = rho*Xn[n-1] + true_nu_vec[state_index]*(An - time_seq[switch_pt_index]*(1- rho) ) + Xn[switch_pt_index]*(1- rho)

	Xn[n] = OUdrift_mean_X + sqrt(OU_variance)*rnorm(1, mean = 0, sd = 1)
			
	}

particle_data$trajectory = cbind(time_seq, Xn)
particle_data$parameters = matrix(nrow = 4, ncol = length(true_kg_vec))
particle_data$parameters[1,] = true_D
particle_data$parameters[2,] = true_kg_vec
particle_data$parameters[3,] = true_nu_vec
particle_data$parameters[4,] = c(0,switch_pts_actual)
rownames(particle_data$parameters) = c("D","kg","nu", "tau" )
return(particle_data)
}






##=================================================
##2D processes 
##=================================================


##-------------------------------------------------
##BM
##-------------------------------------------------

BM_2D = function(time_seq, true_D){
	particle_data = list()
	
	Xn = c(0); Yn = c(0)
	dt_seq = c(0,diff(time_seq))

	for (n in 2:length(time_seq)){
		Xn[n] = Xn[n-1] + sqrt(2*dt_seq[n]*true_D)*rnorm(n = 1, 0, 1)		
		Yn[n] = Yn[n-1] + sqrt(2*dt_seq[n]*true_D)*rnorm(n = 1, 0, 1)		
		
	}
	
	particle_data$trajectory = cbind(time_seq,Xn, Yn)
	particle_data$parameters = matrix(c(true_D), nrow = 1, ncol = 1 )
	colnames(particle_data$parameters) = c("D")
	return(particle_data)

}
##-------------------------------------------------
##Switching OU with drift
##-------------------------------------------------


OUdrift_2D = function(time_seq, true_D, true_kg, true_nu, true_r){

particle_data = list()

path_length = length(time_seq)
dt_seq = c(0, diff(time_seq))

OU_variance = (true_D/true_kg)*(1 - exp(-2*dt_seq*true_kg))

Xn = c(0); 
Yn = c(0)

for(n in 2:(path_length) ){


	
	rho = exp(-true_kg*dt_seq); 
	An = (time_seq[n] - time_seq[n-1]*rho[n]) - (1/true_kg)*(1- rho[n])
	
	
	OUdrift_mean_X = rho[n]*Xn[n-1] + true_nu*cos(true_r)*An
	OUdrift_mean_Y = rho[n]*Yn[n-1] + true_nu*sin(true_r)*An

	

	Xn[n] = OUdrift_mean_X + sqrt(OU_variance[n])*rnorm(1, mean = 0, sd = 1)
			
	Yn[n] = OUdrift_mean_Y + sqrt(OU_variance[n])*rnorm(1, mean = 0, sd = 1)
	}

particle_data$trajectory = cbind(time_seq, Xn, Yn)
particle_data$parameters = matrix(c(true_D, true_kg, true_nu, true_r), nrow = 1, ncol = 4)
colnames(particle_data$parameters) = c("D","kg","nu", "r" )
return(particle_data)

}



##-------------------------------------------------
##Switching OU with drift
##-------------------------------------------------
##Either the switch points can be randomized or chosen
## Can handle different parameter values for d, kg and nu


OUdrift_switching_2D = function(time_seq, switch_points = list("random",2), true_D, true_kg_vec, true_nu_vec, true_r){

particle_data = list()

num_steps = length(time_seq) -1
path_length = num_steps+1


if(class(switch_points) == "list"){
	num_switch = switch_points[[2]]
	possible_switch_pts = ceiling((path_length*0.1)):floor(path_length*0.9)
	switch_pts_actual = sort((sample(possible_switch_pts, num_switch)))
}else if (class(switch_points) == "numeric"  || class(switch_points) == "integer"){
	num_switch = length(switch_points)
	switch_pts_actual = switch_points
}else{
	stop("Either a list of form: list('random', #switch points), or actual switch points")
}
	
switch_pts = c(0, switch_pts_actual, path_length)
	
my_colors = mellow.color.pal(num_switch+2)
if(sum(as.numeric(my_colors == "#FFFFBF")) >0 ){
	my_colors = my_colors[- which(my_colors == "#FFFFBF")]
}else{
	my_colors = my_colors[1:(num_switch+1)]
	}



## Setting up the state parameters

if(length(true_D) != (num_switch +1)){
	if( length(true_D) == 1){
		true_D = rep(true_D, times = (num_switch +1))
	}else{
		stop(paste("Dimension of initial D values should be one or ", num_switch +1))}
}


if(length(true_kg_vec) != (num_switch +1)){
	if( length(true_kg_vec) == 1){
		true_kg_vec = rep(true_kg_vec, times = (num_switch +1))
	}else{
		stop(paste("Dimension of initial kg values should be one or ", num_switch +1))}
}

if( length(true_nu_vec) != (num_switch +1)){
	if( length(true_nu_vec) == 1){
		true_nu_vec = rep(true_nu_vec, num_switch +1)
	}else{
		stop(paste("Dimension of initial velocities should be one or ", num_switch +1))
	}
}




state_vec = c(); switch_pt_vec = c(); 
k  = 0
for(tt in 1:(num_switch+1)){
k = k+1
state_index = (switch_pts[tt]+1):switch_pts[tt+1]
state_vec[state_index] =  k
switch_pt_vec[state_index] = switch_pts[k]
}

##Starting as OU with drift and OU -- in this case both processes have the same  variance but they have different centered terms

dt_vec = c(0,diff(time_seq))


Xn = c(0); 
Yn = c(0)

for(n in 2:(path_length) ){

	switch_pt_index = max(switch_pt_vec[n],1)
	state_index = state_vec[n]
	
	OU_variance = (true_D[state_index]/true_kg_vec[state_index])*(1 - exp(-2*dt_vec[n]*true_kg_vec[state_index]))
	
	rho = exp(-true_kg_vec[state_index]*dt_vec[n]); 
	An = (time_seq[n] - time_seq[n-1]*rho) - (1/true_kg_vec[state_index])*(1- rho)
	OUdrift_mean_X = rho*Xn[n-1] + true_nu_vec[state_index]*cos(true_r)*(An - time_seq[switch_pt_index]*(1- rho) ) + Xn[switch_pt_index]*(1- rho)
	OUdrift_mean_Y = rho*Yn[n-1] + true_nu_vec[state_index]*sin(true_r)*(An - time_seq[switch_pt_index]*(1- rho) ) + Yn[switch_pt_index]*(1- rho)

	

	Xn[n] = OUdrift_mean_X + sqrt(OU_variance)*rnorm(1, mean = 0, sd = 1)
			
	Yn[n] = OUdrift_mean_Y + sqrt(OU_variance)*rnorm(1, mean = 0, sd = 1)
	}

particle_data$trajectory = cbind(time_seq, Xn, Yn)
particle_data$parameters = matrix(nrow = 4, ncol = length(true_kg_vec))
particle_data$parameters[1,] = true_D
particle_data$parameters[2,] = true_kg_vec
particle_data$parameters[3,] = true_nu_vec
particle_data$parameters[4,] = c(0,switch_pts_actual)
rownames(particle_data$parameters) = c("D","kg","nu", "tau" )
return(particle_data)
}




##=================================================
##PLotting functions
##=================================================
mywidth = 2
##-------------------------------------------------
##1D plot
##-------------------------------------------------


plot_sde_1D = function(time_seq, Xn,pos_limits= NA){
	if(is.na(pos_limits)[1] == TRUE){
		pos_limits = c(-1,1)*max(abs(Xn))
		}
	plot(time_seq, Xn , type = "l", xlim = range(time_seq), ylim = pos_limits, xlab = "Time (seconds)", ylab = expression(paste("Position (", mu, "m)")), lwd = mywidth )
	abline(h = 0, col = "firebrick3", lty = 2)
}


##-------------------------------------------------
##2D plot
##-------------------------------------------------

plot_sde_2D = function(time_seq, Xn, Yn, pos_limits = NA){
	if( is.na(pos_limits)[1] == TRUE){
		pos_limits = range(c(Xn, Yn))
	}
par(mfrow = c(1,3))

plot(Xn, Yn , type = "l", xlim = range(pos_limits), ylim = pos_limits, xlab = expression(paste("X Position (", mu, "m)")), ylab = expression(paste("Y Position (", mu, "m)")),lwd = mywidth)
abline(h = 0, col = "firebrick3", lty = 2)
abline(v = 0, col = "firebrick3", lty = 2)
points(Xn[1], Yn[1], pch = 16, col = "springgreen3")



plot(time_seq, Xn , type = "l", xlim = range(time_seq), ylim = pos_limits, xlab = "Time (seconds)", ylab = expression(paste("X Position (", mu, "m)")),lwd = mywidth)
abline(h = 0, col = "firebrick3", lty = 2)

plot(time_seq, Yn , type = "l", xlim = range(time_seq), ylim = pos_limits, xlab = "Time (seconds)", ylab = expression(paste("Y Position (", mu, "m)")),lwd = mywidth)
abline(h = 0, col = "firebrick3", lty = 2)

}


plot_sde_switch= function(time_seq, switch_points,Xn, Yn = NA, pos_limits = NA){

	switch_pts = c(switch_points, length(time_seq))
	num_switch = length(switch_pts)-2
	my_colors = mellow.color.pal(num_switch +2)

	if(sum(as.numeric(my_colors == "#FFFFBF")) >0 ){
		my_colors = my_colors[- which(my_colors == "#FFFFBF")]
	}else{
		my_colors = my_colors[1:(num_switch+1)]
	}


	if (is.na(Yn[1]) == TRUE){
		if(is.na(pos_limits)[1] == TRUE){pos_limits = c(-1,1)*max(abs(Xn))}
		plot(time_seq[(switch_pts[1]+1):switch_pts[2]], Xn[(switch_pts[1]+1):switch_pts[2]] , type = "l", xlim = range(time_seq), ylim = pos_limits, xlab = "Time (seconds)", ylab = expression(paste("X Position (", mu, "m)")), col = my_colors[1],lwd = mywidth)
	k  = 0
	for(tt in 2:(num_switch+1)){
		state_index = (switch_pts[tt]):switch_pts[tt+1]
		lines(time_seq[state_index], Xn[state_index], col = my_colors[tt],lwd = mywidth)
	}
abline(v = time_seq[switch_pts[-c(1, length(switch_pts))]], lty = 4,  col = my_colors[1:(length(switch_pts)-1)])
abline(h = 0, col = "firebrick3", lty = 2)


#-------------------------------------------------		
## FOR 2D SDEs		
#-------------------------------------------------		

	}else{
		
	par(mfrow = c(1,3))
	if(is.na(pos_limits)[1] == TRUE){pos_limits = range(c(Xn, Yn))}
	plot(Xn[(switch_pts[1]+1):switch_pts[2]], Yn[(switch_pts[1]+1):switch_pts[2]] , type = "l", xlim = range(pos_limits), ylim = pos_limits, xlab = expression(paste("X Position (", mu, "m)")), ylab = expression(paste("Position (", mu, "m)")), col = my_colors[1],lwd = mywidth)
	for(tt in 2:(num_switch+1)){
		state_index = (switch_pts[tt]):switch_pts[tt+1]
		lines(Xn[state_index], Yn[state_index], col = my_colors[tt],lwd = mywidth)
	}
abline(h = 0, col = "firebrick3", lty = 2)
abline(v = 0, col = "firebrick3", lty = 2)


	plot(time_seq[(switch_pts[1]+1):switch_pts[2]], Xn[(switch_pts[1]+1):switch_pts[2]] , type = "l", xlim = range(time_seq), ylim = pos_limits, xlab = "Time (seconds)", ylab = expression(paste("X Position (", mu, "m)")), col = my_colors[1],lwd = mywidth)
	for(tt in 2:(num_switch+1)){
		state_index = (switch_pts[tt]):switch_pts[tt+1]
		lines(time_seq[state_index], Xn[state_index], col = my_colors[tt],lwd = mywidth)
	}
abline(v = time_seq[switch_pts[-c(1, length(switch_pts))]], lty = 4,  col = my_colors[1:(length(switch_pts)-1)])
abline(h = 0, col = "firebrick3", lty = 2)


	plot(time_seq[(switch_pts[1]+1):switch_pts[2]], Yn[(switch_pts[1]+1):switch_pts[2]] , type = "l", xlim = range(time_seq), ylim = pos_limits, xlab = "Time (seconds)", ylab = expression(paste("Y Position (", mu, "m)")), col = my_colors[1],lwd = mywidth)
	for(tt in 2:(num_switch+1)){
		state_index = (switch_pts[tt]):switch_pts[tt+1]
		lines(time_seq[state_index], Yn[state_index], col = my_colors[tt], lwd = mywidth )
}
abline(v = time_seq[switch_pts[-c(1, length(switch_pts))]], lty = 4,  col = my_colors[1:(length(switch_pts)-1)])		
abline(h = 0, col = "firebrick3", lty = 2)
	
		
		
	}
	
	
}



find_time = function(time_seq, specific_time){
	return(which(round(time_seq - specific_time,5) == 0))
}

