library(shiny)
library(dplyr)
library(readr)
library(ggplot2)
library(mvtnorm)
library(rgdal)
library(raster)

### Read model and shapefile
ae_alb <- read_rds("data/ae_albopictus_mod.rds")
ae_aeg <- read_rds("data/ae_aegypti_mod.rds")
tmp <- fread("county_lookup.csv")
county_lookup <- tmp$county
names(county_lookup) <- tmp$name
ct_shp <- shapefile("shp/fl_cnt.shp")
availability <- ct_shp$County %in% names(county_lookup)
availability <- ifelse(availability, "yes", "no")

### Define colour for plotting
pal <- function (x) {
    case_when(
        x == "yes" ~ "#6f68ff",
        x == "selected" ~ "#ff5e92",
        T ~ "#777777"
    )
}

### Define server logic
shinyServer(function(input, output) {
    ## Colour for county selected vs available
    available <- reactive(
        {
            av <- availability
            cty_selected <- names(which(input$county == county_lookup))
            ind <- which(ct_shp$County == cty_selected)
            av[ind] <- "selected"
            av
        }
    )
    
    ## Calculate predicted values and make into a dataframe
    out_df <- reactive(
        {
            # Manipulate input into a matrix (of covariates)
            form <- ~ hum_pop_dens + # human population density
                wind + tmin + delta_tmax + rh + # climatic variables
                trap
            in_l <- list(hum_pop_dens = input$hum_pop_dens,
                         wind = input$wind,
                         tmin = input$tmin,
                         delta_tmax = input$delta_tmax,
                         rh = input$rh,
                         county = input$county)
            in_df <- as.data.frame(in_l, stringAsFactors = F)
            in_df$county <- as.character(in_df$county)
            in_df <- in_df[rep(1, 3),]
            in_df$trap <- factor(c("baseline", "light", "other"), 
                                 levels = c("baseline", "light", "other"))
            in_mat <- model.matrix(form, in_df)
            
            ## Aedes albopictus
            # Extract parameters from model and Monte Carlo
            betas <- t(rmvnorm(1000, c(ae_alb$fe$cond, ae_alb$fe$zi), ae_alb$fvcov[1:16, 1:16]))
            re_cond <- rnorm(1000, ae_alb$re$cond$county[in_df$county,], ae_alb$rsd[1,3])
            zi_cond <- rnorm(1000, ae_alb$re$zi$county[in_df$county,], ae_alb$rsd[2,3])
            occur <- (1 - arm::invlogit(t(in_mat %*% (betas[9:16,])) + zi_cond))
            cond <- exp(t(in_mat %*% (betas[1:8,])) + re_cond)
            pred_overall <- occur * cond
            
            # Calculate mean and CI
            df_ae_alb <- data.frame(
                species = "Aedes albopictus",
                trap = in_df$trap,
                occ_mean = colMeans(occur),
                occ_lci = apply(occur, 2, quantile, 0.025),
                occ_uci = apply(occur, 2, quantile, 0.975),
                cond_mean = colMeans(cond),
                cond_lci = apply(cond, 2, quantile, 0.025),
                cond_uci = apply(cond, 2, quantile, 0.975),
                overall_mean = colMeans(pred_overall),
                overall_lci = apply(pred_overall, 2, quantile, 0.025),
                overall_uci = apply(pred_overall, 2, quantile, 0.975),
                stringsAsFactors = F
            )
            
            ## Aedes aegypti
            betas <- t(rmvnorm(1000, c(ae_aeg$fe$cond, ae_aeg$fe$zi), ae_aeg$fvcov[1:16, 1:16]))
            re_cond <- rnorm(1000, ae_aeg$re$cond$county[in_df$county,], ae_aeg$rsd[1,3])
            zi_cond <- rnorm(1000, ae_aeg$re$zi$county[in_df$county,], ae_aeg$rsd[2,3])
            occur <- (1 - arm::invlogit(t(in_mat %*% (betas[9:16,])) + zi_cond))
            cond <- exp(t(in_mat %*% (betas[1:8,])) + re_cond)
            pred_overall <- occur * cond
            
            df_ae_aeg <- data.frame(
                species = "Aedes aegypti",
                trap = in_df$trap,
                occ_mean = colMeans(occur),
                occ_lci = apply(occur, 2, quantile, 0.025),
                occ_uci = apply(occur, 2, quantile, 0.975),
                cond_mean = colMeans(cond),
                cond_lci = apply(cond, 2, quantile, 0.025),
                cond_uci = apply(cond, 2, quantile, 0.975),
                overall_mean = colMeans(pred_overall),
                overall_lci = apply(pred_overall, 2, quantile, 0.025),
                overall_uci = apply(pred_overall, 2, quantile, 0.975),
                stringsAsFactors = F
            )
            
            bind_rows(df_ae_alb, df_ae_aeg)
        }
        
    )
    
    ## Make plots
    output$occurrence_plot <- renderPlot({
        ggplot(out_df(), aes(x = species, y = occ_mean, colour = trap, fill = trap)) +
            geom_errorbar(aes(ymin = occ_mean, ymax = occ_mean),
                          width = 0.5, size = 2,
                          position = position_dodge(width = 0.6)) +
            geom_crossbar(aes(ymin = occ_lci, ymax = occ_uci),
                          width = 0.5,
                          alpha = 0.4,
                          linetype = 0,
                          position = position_dodge(width = 0.6)) +
            labs(x = "Species", y = "Occurrence Probability") +
            theme_bw()

        })
    
    output$cond_plot <- renderPlot({
        ggplot(out_df(), aes(x = species, y = cond_mean, colour = trap, fill = trap)) +
            geom_errorbar(aes(ymin = cond_mean, ymax = cond_mean),
                          width = 0.5, size = 2,
                          position = position_dodge(width = 0.6)) +
            geom_crossbar(aes(ymin = cond_lci, ymax = cond_uci),
                          width = 0.5,
                          alpha = 0.4,
                          linetype = 0,
                          position = position_dodge(width = 0.6)) +
            labs(x = "Species", y = "Abundance if occur") +
            theme_bw()
        
    })
    
    output$overall_plot <- renderPlot({
        ggplot(out_df(), aes(x = species, y = overall_mean, colour = trap, fill = trap)) +
            geom_errorbar(aes(ymin = overall_mean, ymax = overall_mean),
                          width = 0.5, size = 2,
                          position = position_dodge(width = 0.6)) +
            geom_crossbar(aes(ymin = overall_lci, ymax = overall_uci),
                          width = 0.5,
                          alpha = 0.4,
                          linetype = 0,
                          position = position_dodge(width = 0.6)) +
            labs(x = "Species", y = "Overall abundance") +
            theme_bw()
        
    })
    
    ## Make map
    output$cntmap <- renderLeaflet(
        leaflet(data = ct_shp) %>%
            addPolygons(
                fillColor = ~pal(available()),
                color = "#000", weight = 1,
                opacity = 0.8,
                fillOpacity = 0.8,
                label = ~htmltools::htmlEscape(paste0(County))
            )
    )

})
