#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
        
        # Application title
        titlePanel("Analysis of Resonator ImmunoAssays (ARIA)"),
        
        # Sidebar with a slider input for number of bins 
        sidebarLayout(
                sidebarPanel(
                        textInput("plotName", label = "Plot Name", 
                                  placeholder = "e.g., XPP-01a"),
                        textInput("cntl", label = "Control Rings", 
                                  placeholder = "e.g., thermal, raw, blank"),
                        fileInput("chipLayout", "Select Chip Layout",
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        fileInput("filesList", multiple = TRUE, "Select Rings",
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        radioButtons("ch", label = "Experiment Type:", 
                                     choices = c("U-Channel" = "u.ch",
                                                 "2 Channel" = "two.ch")),
                        radioButtons("avg", label = "Average Clusters:",
                                     choices = c("No", "Yes")),
                        textInput("time2", label = "Start Time (min)",
                                  placeholder = "e.g., 22"),
                        textInput("time1", label = "End Time (min)",
                                  placeholder = "e.g., 42"),
                        actionButton("keepData", "Download Data"),
                        actionButton("goButton", "Go!")
                ),
                
                # Show a plot of the generated distribution
                mainPanel(
                        
                        plotOutput("run"), plotOutput("netShifts")
                )
        )
)

# Define server logic
server <- function(input, output) {
        
        output$run <- renderPlot({
                
                req(input$goButton)
                req(input$plotName)
                req(input$cntl)
                library(tidyverse)
                library(ggthemes)
                
                # define recipe as global variable for use in other functions
                recipe <<- read_csv(input$chipLayout[,4])
                colnames(recipe)[1] <<- "Target" # rename col; remove byte mark
                targets <- recipe$Target
                
                # generate list of rings to analyze (gets all *.csv files)
                rings <- input$filesList
                
                # create empty data frame to store data
                df <- data.frame()
                
                # add data to data frame corresponding for each ring in rings
                for (i in 1:nrow(rings)) {
                        ring <- as.vector(rings[i,1])
                        dat <- read_csv(rings[i,4], col_names = FALSE)
                        time_shift <- dat[ ,1]
                        shift <- dat[ ,2]
                        ringStr <- strsplit(ring, "\\.")[[1]]
                        ringNum <- as.numeric(ringStr[1])
                        recipe.col <- which(recipe$Ring == ringNum)
                        groupNum <- recipe$Group[recipe.col]
                        ring <- rep(ringNum, nrow(dat))
                        group <- rep(groupNum, nrow(dat))
                        groupName <- as.character(recipe$Target[[recipe.col]])
                        groupName <- rep(groupName, nrow(dat))
                        channel <- recipe$Channel[[recipe.col]]
                        channel <- rep(channel, nrow(dat))
                        run <- rep(input$plotName, nrow(dat))
                        time_point <- seq(1:nrow(dat))
                        tmp <- data.frame(ring, group, time_shift, shift, 
                                          groupName, channel, run, time_point)
                        df <- rbind(df, tmp)
                }
                
                # renames columns in df
                names(df) <- c("Ring", "Group", "Time", "Shift", "Target", 
                               "Channel", "Experiment", "Time Point")
                
                if (input$cntl != "raw"){
                        # get thermal control averages
                        controls <- filter(df, Target == input$cntl)
                        ringList <- unique(controls$Ring)
                        
                        # gets times from first thermal control
                        times <- filter(controls, 
                                        Ring == ringList[1]) %>% select(Time)
                        df.controls <- data.frame(times)
                        
                        # create dataframe with all controls
                        for (i in ringList){
                                ringShift <- filter(controls, Ring == i) %>% 
                                        select(Shift)
                                names(ringShift) <- paste('Ring', i, sep='')
                                df.controls <- cbind(df.controls, ringShift)
                        }
                        
                        # averages thermal controls
                        cols <- ncol(df.controls)
                        if (length(unique(controls$Ring)) != 1) {
                                df.controls$avgControls <- 
                                        rowMeans(df.controls[,c(2:cols)])
                        } else {
                                df.controls$avgControls <- 
                                        df.controls[,c(2:cols)]
                        }
                        avgControls <- as.vector(df.controls$avgControls)
                        
                        #subtracts thermal controls from each ring
                        ringNames <- unique(df$Ring)
                        for(i in ringNames){
                                ringDat <- filter(df, Ring == i) %>% 
                                        select(Shift)
                                ringTC <- ringDat - avgControls
                                df[df$Ring == i, 4] <- ringTC
                        }
                }
                
                if (input$avg == "Yes"){
                        if (input$ch == "two.ch"){
                                dat.avg <- df %>%
                                        group_by(Target, `Time Point`, 
                                                 Channel) %>% 
                                        summarise_each(funs(mean, sd), 
                                                       c(Time, Shift))
                        } else {
                                dat.avg <- df %>%
                                        group_by(Target, `Time Point`) %>% 
                                        summarise_each(funs(mean, sd), 
                                                       c(Time, Shift))
                        }
                        
                        
                        plots <- ggplot(dat.avg, aes(Time_mean, Shift_mean, 
                                                     color = Target)) + 
                                geom_line() + 
                                geom_ribbon(aes(ymin = Shift_mean - Shift_sd, 
                                                ymax = Shift_mean + Shift_sd, 
                                                linetype = NA), 
                                            fill = "slategrey", alpha = 1/8) +
                                labs(x = "Time (min)",
                                     y = expression(paste("Relative Shift (",
                                                          Delta,"pm)")),
                                     color = "Target") +
                                ggtitle(input$plotName)
                        theme_few(base_size = 22)
                        
                } else {
                        plots <- ggplot(df, aes(Time, Shift, 
                                                group = factor(Ring),
                                                color = factor(Target))) + 
                                labs(x = "Time (min)", 
                                     y = expression(paste("Relative Shift (",
                                                          Delta,"pm)")),
                                     color = "Target") +
                                geom_line() + 
                                ggtitle(input$plotName) +
                                theme_few(base_size = 22)
                }
                
                
                if (input$ch == "two.ch") {
                        plots <- plots + facet_grid(.~Channel)
                }
                
                #ggsave(plots, filename = paste0(input$plotName, "_",
                #                                input$ch, "_", 
                #                                input$cntl, ".png"),
                #       device = "png", width = 8, height = 6)
                
                #write_csv(dat, paste0(input$plotName, "_", input$ch, "_",
                #                      input$cntl, ".csv"))
                
                plots
                
        })
        output$downloadPlot <- downloadHandler({
                filename = function() { paste(input$plotName, '.png', sep='') }
                content = function(file) {
                        ggsave(plots, plot = plotInput(), device = "png")
                }
        }
        )
}

# Run the application 
shinyApp(ui = ui, server = server)

