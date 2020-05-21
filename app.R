library(shiny)
library(readxl)
library(ggplot2)
library(kableExtra)
library(tidyr)
library(data.table)
library(tidyverse)

Ginty_2019_Neuron_Receptor_RNASeq <- read_excel("Ginty_2019_Neuron_Receptor_RNASeq.xlsx")
Ginty_2019_Neuron_Receptor_RNASeq_Mean <- read_excel("Ginty_2019_Neuron_Receptor_RNASeq.xlsx", 
                                                     sheet = "Mean_Values")
Ginty_2019_Neuron_Receptor_RNASeq_SD <- read_excel("Ginty_2019_Neuron_Receptor_RNASeq.xlsx",    
                                                   sheet = "SD_Values")

ui <- fluidPage(
    textInput("text", label = h3("Type exact gene name (check attached spreadsheet first)"), value = "Cacna1h"),
    hr(),
    fluidRow(column(3, verbatimTextOutput("value"))),
    plotOutput("RNAPlot")
    
)


# Define server logic required to draw a histogram
server <- function(input, output) {

    
    output$RNAPlot <- renderPlot({Chosen_Genes <- input$text
    Restricted_Mean <- Ginty_2019_Neuron_Receptor_RNASeq_Mean %>%
        filter(Gene %in% Chosen_Genes)
    Restricted_SD <-  Ginty_2019_Neuron_Receptor_RNASeq_SD %>%
        filter(Gene %in% Chosen_Genes)
    Restricted_SDTable <- pivot_longer(
        data = Restricted_SD, 
        cols = c("NonPeptidergic_Noc","Peptidergic_Noc","C","Adelta", "ABetaRA", "ABetaSA", "Abetafield"),
        names_to = "cell_type",
        values_to = "transcript_rpkm_sd" 
    )
    n <- c(3, 3, 3, 3, 5, 3, 3)
    SD <- Restricted_SDTable$transcript_rpkm_sd
    SEM <- SD/sqrt(n)
    Restricted_MeanTable <- pivot_longer(
        data = Restricted_Mean, 
        cols = c("NonPeptidergic_Noc","Peptidergic_Noc","C","Adelta", "ABetaRA", "ABetaSA", "Abetafield"),
        names_to = "cell_type",
        values_to = "transcript_rpkm"
    ) %>%
        add_column(SD) %>% 
        add_column(SEM) %>% 
        add_column(Cell_Type = rep(c("Non Peptidergic Nociceptor", "Peptidergic Nociceptor", "C-LTMR", "A-Delta LTMR", "A-Beta Rapidly Adapting", "A-Beta Slowly Adapting", "A-Beta Field")))
        ggplot(data = Restricted_MeanTable, aes(Cell_Type, transcript_rpkm, fill = Cell_Type)) + 
            facet_wrap(ncol = 3,
                       ~fct_rev(Gene)) +
            geom_col(width = 0.2) + 
            ylab("Transcript RPKM value") + 
            xlab("Cell Type") +
            geom_errorbar(data = Restricted_MeanTable, aes(ymin = transcript_rpkm-SEM, ymax = transcript_rpkm+SEM), width = 0.1) +
            theme(axis.text.x = element_text(angle = 60, hjust = 1))+ theme(panel.border = element_blank(), panel.grid.major =   element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
            geom_text(aes(label = round(transcript_rpkm, 1)), nudge_y = 50, angle = 60)
    })
}
   
# Run the application 
shinyApp(ui = ui, server = server)



