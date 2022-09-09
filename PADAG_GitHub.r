install.packages("tidyr")
install.packages("dplyr")
install.packages("fpow")
install.packages("shiny")
install.packages("plotly")
install.packages("DT")

library(fpow)
library(dplyr) # for mutate()
library(shiny)
library(plotly)
library(DT)

#library(utils)
library(tidyr)
#library(rgl)
#library(car)

ui <- fluidPage(
# Application title	
titlePanel(tags$h4(p(strong("PADAG (Power Analysis and Data Augmentation with Patients Missing Genotypic Data)")))),
br(),
# Sidebar 
sidebarLayout(
        sidebarPanel(
		textInput("n", "Number of Families**:", placeholder = "Enter one value or values separated by a comma...",width = "500px"),
textInput("r", "Coefficient of Relationship**:", placeholder = "Enter one value or values separated by a comma...",width = "500px"),						 						 
radioButtons("outcome",
                         "Outcome**:",
                         c("Continuous" = 1,
                           "Dichotomous" = 0),inline = T),
textInput("heritability_t","Heritability**:", placeholder = "Enter one value or values separated by a comma...",width = "500px"),
radioButtons("power_or_ncp",
                         "Is baseline power or noncentrality parameter (NCP) available?**",
                         c("Baseline Power" = 0,
                           "Baseline NCP" = 1
						   ),width="500px", inline = T),
textInput("power_t",
                        "Baseline Power (e.g.,0.75)**:", placeholder = "Enter one decimal or decimals separated by a comma...",width = "500px"),

textInput("ncp1_t",
                        "Baseline NCP:", placeholder = "Enter one value or values separated by a comma...",width = "500px"),
textInput("rw",
                        "Weighted Coefficient of Relationship:", placeholder = "Enter one value or values separated by a comma...",width = "500px"),
radioButtons("side",                      
						 "Alternative Hypothesis**:",
                         c("Two-Sided" = 1,
                           "One-Sided" = 0),inline = T),
numericInput("alpha",
                         "Type I error rate**:",
                         value="0.05", max=1,
                         min=0), 
						 
						 actionButton("go" ,"Calculate Power with Data Augmentation")# class = "btn btn-primary")	
		),
		
#actionButton("go" ,"Calculate Power with Data Augmentation", 
#                         icon("arrow-circle-right")),# class = "btn btn-primary")	
#		),					 
  # Show a plot of the augumented power 
        #mainPanel(
        # rglwidgetOutput("plot",  width = 800, height = 600)
        #)
		
		mainPanel(

          # Output: Tabset w/ plot, summary, and table ----
          tabsetPanel(type = "tabs",
                      tabPanel("Data Output", 
                               # Output: HTML table with requested number of observations ----
                               h3("Input values (**) are used to obtain statistical power with data augmentation:"),
                               #tableOutput("table")
							   DT::dataTableOutput("table")
                               #h3("II. Second title:"),
                               #tableOutput("VUL")   
                      ),
                  tabPanel("Plot", plotlyOutput("plot"))   
         )
        )
		
		
		)

)

server <- function(input, output, session) {

observe({
    if (input$power_or_ncp=="1") {#Baseline NCP
      isolate(updateTextInput(session,"ncp1_t",label="Baseline NCP**:", placeholder = "Enter one value or values separated by a comma..."))
      isolate(updateTextInput(session,"power_t",label="Baseline Power (e.g.,0.75):", placeholder = "Enter one decimal or decimals separated by a comma..."))
    } else if (input$power_or_ncp=="0") {
      isolate(updateTextInput(session,"ncp1_t",label="Baseline NCP:", placeholder = "Enter one value or values separated by a comma..."))
      isolate(updateTextInput(session,"power_t",label="Baseline Power (e.g.,0.75)**:", placeholder = "Enter one decimal or decimals separated by a comma..."))
    } 
 })
 
observe({
	    if (input$outcome=="1") {					
      isolate(updateTextInput(session,"rw",label="Weighted Coefficient of Relationship:", placeholder = "Enter one value or values separated by a comma..."))
    } else if (input$outcome=="0") {
      isolate(updateTextInput(session,"rw",label="Weighted Coefficient of Relationship**:", placeholder = "Enter one value or values separated by a comma..."))
    } 
 })
  
Output <-  reactiveValues(datasetInput = NULL)
storage <- reactiveValues()
	  

	   
   re<- observeEvent(input$go,{
 if (input$outcome=="1") {	#"Continuous" = 1
	if (input$power_or_ncp=="0") { #Baseline power
	storage$para <- data.frame(n=input$n,power1=input$power_t,r=input$r,heritability=input$heritability_t, alpha_F=input$alpha)
	   storage$para2<-storage$para %>% separate_rows(n) %>% separate_rows(power1) %>% separate_rows(r) %>% separate_rows(heritability) 
	   storage$para3<-storage$para2 %>%  mutate(n =as.numeric(n)) %>% mutate(power1=as.numeric(power1)) %>% mutate(r=as.numeric(r)) %>% mutate(heritability=as.numeric(heritability))  
	   storage$para3<-storage$para3  %>%  mutate(rho =r*heritability) %>% mutate(df2=n-2) %>% mutate(type2=1-power1) 
	   storage$para3<-storage$para3  %>% mutate(type2=as.numeric(type2)) %>% mutate(df2=as.numeric(df2)) %>% mutate(rho=as.numeric(rho))
	   
	   if (input$side=="1"){ #"Two-Sided" = 1	
	   for(i in 1:dim(storage$para3)[1]){
		storage$para3[i,"ncp1"] <-  ncparamF(storage$para3[i,"alpha_F"],storage$para3[i,"type2"],1,storage$para3[i,"df2"])
		}
	   storage$para3<-storage$para3  %>% mutate(ncp2=ncp1*(1-2*rho*r+r^2)/(1-rho^2)) 
	   storage$para3<-storage$para3 %>% mutate(ncp2_ncp1=ncp2/ncp1) 
	   storage$para3<-storage$para3  %>% mutate(power2=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F)) 
	   } else if (input$side=="0"){ #"Two-Sided" = 0 
	   storage$para3<-storage$para3 %>% mutate(alpha_F=alpha_F*2)
	   for(i in 1:dim(storage$para3)[1]){
		storage$para3[i,"ncp1"] <-  ncparamF(storage$para3[i,"alpha_F"],storage$para3[i,"type2"],1,storage$para3[i,"df2"])
		}
	   storage$para3<-storage$para3  %>% mutate(ncp2=ncp1*(1-2*rho*r+r^2)/(1-rho^2)) 
	   storage$para3<-storage$para3 %>% mutate(ncp2_ncp1=ncp2/ncp1) 
	   storage$para3<-storage$para3  %>% mutate(power2=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F)) 
	   }
	   } 
else if (input$power_or_ncp=="1"){ #Baseline NCP
       storage$para <- data.frame(n=input$n,ncp1=input$ncp1_t,r=input$r,heritability=input$heritability_t) 
	   storage$para2<-storage$para %>% separate_rows(n) %>% separate_rows(r) %>% separate_rows(heritability) %>% separate_rows(ncp1)
	   storage$para3<-storage$para2 %>%  mutate(n =as.numeric(n)) %>% mutate(heritability=as.numeric(heritability)) %>% mutate(ncp1=as.numeric(ncp1)) %>% mutate(r=as.numeric(r)) %>%
	   mutate(rho =r*heritability) %>% mutate(ncp2=ncp1*(1-2*rho*r+r^2)/(1-rho^2)) %>% mutate(ncp2_ncp1=ncp2/ncp1) 

		if (input$side=="1"){ #"Two-Sided" = 1					   
	   storage$para3<-storage$para3 %>% mutate(df2=n-2)  %>% mutate(power1=pf(qf(input$alpha,1,df2,lower.tail=F), 1, df2, ncp1, lower.tail = F))  %>% 
	   mutate(power2=pf(qf(input$alpha,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F)) 
	   } else if (input$side=="0"){ #"Two-Sided" = 0
	   storage$para3<-storage$para3 %>% mutate(alpha_F=input$alpha*2)
	   storage$para3<-storage$para3 %>% mutate(df2=n-2)  %>% mutate(power1=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp1, lower.tail = F))  %>% 
	   mutate(power2=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F)) 
	   }
	   }	   
	   }
	   
else if (input$outcome=="0") {#Dichotomous
	   if (input$power_or_ncp=="0") {
	   req(input$n,input$power_t,input$r,input$heritability_t,input$alpha)
	   storage$para <- data.frame(n=input$n,power1=input$power_t,r=input$r,heritability=input$heritability_t,rw=input$rw,alpha_F=input$alpha)
	   #storage$para2<-separate_rows(storage$para,n) 
	   storage$para2<-storage$para %>% separate_rows(n) %>% separate_rows(power1) %>% separate_rows(r) %>% separate_rows(heritability) %>% separate_rows(rw)
	   storage$para3<-storage$para2 %>%  mutate(n=as.numeric(n)) %>% mutate(power1=as.numeric(power1)) %>% mutate(r=as.numeric(r))  %>% mutate(rw=as.numeric(rw)) %>% mutate(heritability=as.numeric(heritability))  
	   storage$para3<-storage$para3  %>%  mutate(rho =r*heritability) %>% mutate(df2=n-2) %>% mutate(type2=1-power1) 
	   storage$para3<-storage$para3  %>%  mutate(type2=as.numeric(type2)) %>% mutate(df2=as.numeric(df2)) %>% mutate(rho=as.numeric(rho))
		if (input$side=="1"){ #"Two-Sided" = 1	
		for(i in 1:dim(storage$para3)[1]){
		storage$para3[i,"ncp1"] <-  ncparamF(storage$para3[i,"alpha_F"],storage$para3[i,"type2"],1,storage$para3[i,"df2"])
		} 
	   storage$para3<-storage$para3  %>% mutate(ncp2=ncp1*(1-2*r*rho*rw+r^2)/(1-(rho*rw)^2)) 
	   storage$para3<-storage$para3 %>%  mutate(ncp2_ncp1=ncp2/ncp1) 
	   storage$para3<-storage$para3  %>% mutate(power2=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F)) 
		} else if (input$side=="0"){ #"Two-Sided" = 0 
		storage$para3<-storage$para3 %>% mutate(alpha_F=alpha_F*2)
			for(i in 1:dim(storage$para3)[1]){
			storage$para3[i,"ncp1"] <-  ncparamF(storage$para3[i,"alpha_F"],storage$para3[i,"type2"],1,storage$para3[i,"df2"])
			}
		storage$para3<-storage$para3  %>% mutate(ncp2=ncp1*(1-2*r*rho*rw+r^2)/(1-(rho*rw)^2)) 
		storage$para3<-storage$para3 %>%  mutate(ncp2_ncp1=ncp2/ncp1) 
		storage$para3<-storage$para3  %>% mutate(power2=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F)) 
		}
		}
	    else if (input$power_or_ncp=="1"){#Baseline NCP
	   storage$para <- data.frame(n=input$n,ncp1=input$ncp1_t,r=input$r,heritability=input$heritability_t,rw=input$rw) 
	   storage$para2<-storage$para %>% separate_rows(n) %>% separate_rows(r) %>% separate_rows(heritability) %>% separate_rows(ncp1) %>% separate_rows(rw)
	   storage$para3<-storage$para2 %>%  mutate(n =as.numeric(n)) %>% mutate(heritability=as.numeric(heritability)) %>% mutate(ncp1=as.numeric(ncp1)) %>% mutate(rw=as.numeric(rw)) %>% mutate(r=as.numeric(r)) %>%
	   mutate(rho =r*heritability) %>% mutate(ncp2=ncp1*(1-2*r*rho*rw+r^2)/(1-(rho*rw)^2)) %>% mutate(ncp2_ncp1=ncp2/ncp1) 
	   if (input$side=="1"){ #"Two-Sided" = 1		
	   storage$para3<-storage$para3 %>% mutate(df2=n-2)  %>% mutate(power1=pf(qf(input$alpha,1,df2,lower.tail=F), 1, df2, ncp1, lower.tail = F))  %>% 
	   mutate(power2=pf(qf(input$alpha,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F))
	   } else if (input$side=="0"){ #"Two-Sided" = 0
	   storage$para3<-storage$para3 %>% mutate(alpha_F=input$alpha*2)
	   storage$para3<-storage$para3 %>% mutate(df2=n-2)  %>% mutate(power1=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp1, lower.tail = F))  %>% 
	   mutate(power2=pf(qf(alpha_F,1,df2,lower.tail=F), 1, df2, ncp2, lower.tail = F)) 
	   }
	   } 
	   }
		
Output$datasetInput<-storage$para3

output$table=DT::renderDataTable({
       datatable(Output$datasetInput) %>% 
           formatPercentage(c("power1", "power2"), 2) %>%
		   formatRound(columns=c('ncp2','ncp1', 'ncp2_ncp1'), digits=3) 
		   })

output$plot<-renderPlotly({
plot_ly(storage$para3, x = ~n,y = ~power1, z = ~power2) %>%
add_markers(color=~heritability) 
})
})
}

shinyApp(ui, server)
