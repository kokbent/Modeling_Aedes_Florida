library(data.table)
library(dplyr)
library(shiny)
library(leaflet)
library(shinydashboard)
tmp <- fread("county_lookup.csv")
county_lookup <- tmp$county
names(county_lookup) <- tmp$name

header <- dashboardHeader(
    title = "Occurrence of Aedes in Florida",
    titleWidth = 450
)

sidebar <- dashboardSidebar(
    width = 400,
    selectInput(
        "county", 
        "County",
        choices = county_lookup
    ),
    
    numericInput(
        "hum_pop_dens", 
        "(Local) Human population density (units?)",
        min = 0,
        max = 5225,
        step = 0.5,
        value = 480.5
    ),
    
    sliderInput(
        "wind", 
        "Wind speed (m/s)",
        min = 2.5,
        max = 17.5,
        step = 0.1,
        value = 5.5
    ),
    
    sliderInput(
        "tmin", 
        "Minimum temperature (celcius)",
        min = 3.5,
        max = 28.5,
        step = 0.1,
        value = 22.5
    ),
    
    sliderInput(
        "delta_tmax", 
        "Residual of maximum temperature (celcius)",
        min = -4,
        max = 5,
        step = 0.1,
        value = 0
    ),
    
    sliderInput(
        "rh", 
        "Relative humidity (%)",
        min = 48,
        max = 93,
        step = 0.5,
        value = 77
    )
)

body <- dashboardBody(
    
    fluidRow(
        
        # output map,
        box(
            title = "Geographic location", 
            status = "primary", 
            solidHeader = TRUE,
            width = 6,
            leafletOutput("cntmap")
        ),
        
        # output occurrence
        box(
            title = "Occurence", 
            status = "primary",
            solidHeader = T,
            width = 6,
            plotOutput("occurrence_plot")
        )
    ),
    
    fluidRow(
        # output occurrence
        box(
            title = "Abundance if occur", 
            status = "primary",
            solidHeader = T,
            width = 6,
            plotOutput("cond_plot")
        ),
        
        # output map
        box(
            title = "Overall abundance", 
            status = "primary",
            solidHeader = T,
            width = 6,
            plotOutput("overall_plot")
        )
        
    ),
    
)

ui <- dashboardPage(
    header,
    sidebar,
    body
)