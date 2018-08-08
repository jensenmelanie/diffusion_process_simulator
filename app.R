library(shiny)
source("helper_functions.R")

##=============================================================================
## user interface
##=============================================================================

ui <- fluidPage(

	# App Title-----
	titlePanel("Simulating Diffusion Processes"),
	
	
	
	#output
	plotOutput(outputId = "distPlot"),
	##Sidebar Layout with input and output definitions
	#sidebarLayout(
	hr(),
	
	fluidRow(style = "padding-bottom:200px;",
	
		column(3,
		#Sidebar panel for inputs
		#sidebarPanel(
			
			#Input:Select box for type
			selectInput("SDEclass", label = h3("Type of Diffusion Process"), choices = list("Brownian Motion" = "BM", "Ornstein Uhlenbeck" = "OU", "Switching Ornstein Uhlenbeck" = "switch"), selected = "BM"),
			
			#Input:Select box for SDE Dimension
			selectInput("SDEdim", label = h4("Dimension"), choices = list("1D" = 1, 			"2D " = 2), selected = 1),
			
				#Conditional Input: angle if 2D
				conditionalPanel(
				condition = "input.SDEdim == 2 && input.SDEclass != 'BM' ", sliderInput("SDEangle", label = h5("Angle of motion (radians)"), min = 0, max = pi, value = pi/4))
			
		),
		
		column(3, offset = 1,
			
				#Conditional Input: BM Parameters
				conditionalPanel(condition = "input.SDEclass == 'BM'", 	numericInput("BM_dPar", label = h5("Diffusivity Coefficient"), min = 0, value = 1)),
			
				#Conditional Input: OU parameters
				conditionalPanel(condition = "input.SDEclass == 'OU'", 	numericInput("dPar", label = h5("Diffusivity Coefficient"), min = 0, value = 1, step = 0.5)),
				conditionalPanel(condition = "input.SDEclass == 'OU'", 	numericInput("kgPar", label = h5("Spring Constant-ish"), min = 0, value = 10)),
				conditionalPanel(condition = "input.SDEclass == 'OU'", 	numericInput("nuPar", label = h5("Motor Velocity"), value = 0, step = 0.5)),
				
				#Conditional Input: Switch parameters
				
				conditionalPanel(condition = "input.SDEclass == 'switch'", 	selectInput("switchType", label = h5("Switch Points"),choices = list("Random" = 'random', "User selected" = "given"),selected = "random")),
					conditionalPanel(condition ="input.SDEclass == 'switch' && input.switchType== 'random'", numericInput("numSwitch", label = h6("Number of Switch Points"), min = 1, value = 2 )),
					conditionalPanel(condition ="input.SDEclass == 'switch' &&  input.switchType== 'given'", textInput("timeSwitch", label = h6("Time of Switch Points"), value = "2.5, 7.5" )),			
				conditionalPanel(condition = "input.SDEclass == 'switch'", 	textInput("switch_dPar", label = h5("Diffusivity Coefficients for each state"), value = "1,1,1")),
				conditionalPanel(condition = "input.SDEclass == 'switch'", 	textInput("switch_kgPar", label = h5("Spring Constant-ish for each state"),  value = "10,10,10")),
				conditionalPanel(condition = "input.SDEclass == 'switch'", 	textInput("switch_nuPar", label = h5("Motor Velocity for each state"), value = "1,0, -1"))
		),
		
		column(3,offset = 1,		

		##Time of the process
		numericInput("dt", label = h5("Time Step"), min = 0, value = 0.1, step = 0.1),
		sliderInput("SDEtime", label = h5("Time (seconds)"), min = 0, max = 100, value = c(0, 10)),
		sliderInput("posLimits", label = h4("Plot Position Limits"), min = -50, max = 50, value = c(-10,10)),
		##Action Button
		actionButton("simSDE", "Simulate my diffusion process", style = 'background-color:#84CF41;color:white;text-align:center;padding:15px 32px; font-size:16px; margin-top:10px;margin-left:10px')
		) #end of the last column

		
			) #fluid row
	
	#Main panel for displaying the outputs---
	#mainPanel(
	


		
	
	
	#) #main panel
	
	
	
	#) #sidebar Layout
	

) # fluid Page




##=============================================================================
## Server
##=============================================================================

server <- function(input,output){

	
	output$distPlot <-renderPlot({
	
	input$simSDE
	
	isolate({
		
	time_seq = seq(input$SDEtime[1], input$SDEtime[2], by = input$dt)
	path_length = length(time_seq)
	
		if(input$SDEdim == 1){
			if(input$SDEclass == "BM"){
				mysde = BM_1D(time_seq,input$BM_dPar)
				plot_sde_1D(mysde$trajectory[,1], mysde$trajectory[,2],pos_limits = input$posLimits)
				
			}else if (input$SDEclass == "OU"){
				mysde = OUdrift_1D(time_seq,input$dPar, input$kgPar, input$nuPar)
				plot_sde_1D(mysde$trajectory[,1], mysde$trajectory[,2],pos_limits = input$posLimits)
				
			}else if(input$SDEclass == "switch"){
				if(input$switchType == "random"){
					switch_pts = c(sort(sample(5:(path_length -5), input$numSwitch)))
				}else{
					given_tau = as.numeric(unlist(strsplit(as.character(input$timeSwitch), ",")))
					switch_pts = unlist(lapply(as.list(given_tau), find_time, 	time_seq = time_seq))
				}
				
				D_vec = as.numeric(unlist(strsplit(input$switch_dPar, ",")))
				
				kg_vec = as.numeric(unlist(strsplit(input$switch_kgPar, ",")))
				
				nu_vec = as.numeric(unlist(strsplit(input$switch_nuPar, ",")))

				mysde = OUdrift_switching_1D(time_seq, switch_pts, D_vec, kg_vec, nu_vec)
				
				
				plot_sde_switch(mysde$trajectory[,1],mysde$parameters[4,],mysde$trajectory[,2], pos_limits = input$posLimits)
			}
		
			
		}else if(input$SDEdim ==2){
		
		if(input$SDEclass == "BM"){
				mysde = BM_2D(time_seq,input$dPar)
				plot_sde_2D(mysde$trajectory[,1], mysde$trajectory[,2],mysde$trajectory[,3], pos_limits = input$posLimits)
			}else if (input$SDEclass == "OU"){
				mysde = OUdrift_2D(time_seq,input$dPar, input$kgPar, input$nuPar, input$SDEangle)
				plot_sde_2D(mysde$trajectory[,1], mysde$trajectory[,2],mysde$trajectory[,3], pos_limits = input$posLimits)
	
			}else if(input$SDEclass == "switch"){
							if(input$switchType == "random"){
					switch_pts = sort(sample(5:(path_length -5), input$numSwitch))
				}else{
					given_tau = as.numeric(unlist(strsplit(input$timeSwitch, ",")))
					switch_pts = unlist(lapply(as.list(given_tau), find_time, 	time_seq = time_seq))
				}
				
				D_vec = as.numeric(unlist(strsplit(input$switch_dPar, ",")))
				
				kg_vec = as.numeric(unlist(strsplit(input$switch_kgPar, ",")))

				nu_vec = as.numeric(unlist(strsplit(input$switch_nuPar, ",")))
				
				mysde = OUdrift_switching_2D(time_seq, switch_pts, D_vec, kg_vec, nu_vec, input$SDEangle)
				
				plot_sde_switch(mysde$trajectory[,1],mysde$parameters[4,],mysde$trajectory[,2],mysde$trajectory[,3], pos_limits = input$posLimits)
		
				
			} #End of the class type for 2d
		
		
		} #end of the dimension loop
		
	}) #end of isolate
	
	}) #end of render plot
	
	
	
} # end of the server

shinyApp(ui = ui, server = server)
