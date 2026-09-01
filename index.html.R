library(shiny)
library(shinydashboard)
library(tidyLPA)
library(dplyr)
library(ggplot2)
library(nnet)
library(DT)
library(tidyr)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "SPPLab - LPA Studio", titleWidth = 280),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("1. Data Import", tabName = "data_import", icon = icon("file-upload")),
      menuItem("2. LPA Analysis", tabName = "lpa_analysis", icon = icon("chart-line")),
      menuItem("3. Group Comparisons", tabName = "demographics", icon = icon("users")),
      menuItem("4. Logistic Regression", tabName = "regression", icon = icon("sliders-h")),
      menuItem("5. Longitudinal LPTA", tabName = "lpta_module", icon = icon("exchange-alt"))
    ),
    hr(),
    div(
      style = "padding: 15px; color: #b8c7ce; font-size: 12px; text-align: center;",
      p(strong("Sport and Performance Psychology Lab")),
      p("Open Source Analytics Tools")
    )
  ),
  dashboardBody(
    tabItems(
      # TAB 1: Data Import
      tabItem(tabName = "data_import",
        fluidRow(
          box(title = "Select Data File (.csv)", width = 6, status = "primary", solidHeader = TRUE,
            fileInput("file", "Upload Dataset", accept = c(".csv")),
            helpText("Shinylive ortamı için CSV formatı önerilir.")
          ),
          box(title = "Dataset Preview", width = 12, status = "info",
            DTOutput("raw_data_table")
          )
        )
      ),
      
      # TAB 2: LPA Analysis
      tabItem(tabName = "lpa_analysis",
        fluidRow(
          box(title = "LPA Parameters", width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("lpa_var_select"),
            numericInput("max_k", "Maximum Number of Profiles (K):", value = 4, min = 2, max = 6),
            checkboxInput("scale_vars", "Standardize Variables (Z-Score)", value = TRUE),
            actionButton("run_lpa", "Run LPA Analysis", class = "btn-success", icon = icon("play")),
            hr(),
            p(tags$small(em("Methodological Specification: Model 1 (equal variances, covariances fixed to 0) via tidyLPA.")))
          ),
          box(title = "Profile Plot (Means)", width = 8, status = "info",
            plotOutput("lpa_plot")
          )
        ),
        fluidRow(
          box(title = "Model Fit Indices & Comparative Statistics", width = 12, status = "warning",
            DTOutput("lpa_compare_table")
          )
        ),
        fluidRow(
          box(title = "Class Counts & Proportions", width = 6, status = "info",
            DTOutput("lpa_counts_table")
          ),
          box(title = "Average Latent Class Probabilities", width = 6, status = "info",
            DTOutput("lpa_prob_matrix_table")
          )
        ),
        fluidRow(
          box(title = "Suggested Methodological Citation", width = 12, status = "success", solidHeader = TRUE,
            verbatimTextOutput("citation_text")
          )
        )
      ),
      
      # TAB 3: Group Comparisons
      tabItem(tabName = "demographics",
        fluidRow(
          box(title = "Comparison Variable & Test Selection", width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("demo_var_select"),
            selectInput("demo_type", "Variable Type / Test:", 
                        choices = c("Continuous (Independent t-Test / ANOVA)" = "continuous", 
                                    "Categorical (Chi-Square Test)" = "categorical")),
            actionButton("run_demo", "Compare Groups", class = "btn-primary")
          ),
          box(title = "Statistical Test Output", width = 8, status = "info",
            verbatimTextOutput("demo_results")
          )
        )
      ),
      
      # TAB 4: Logistic Regression
      tabItem(tabName = "regression",
        fluidRow(
          box(title = "Predict Profile Membership", width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("reg_pred_select"),
            actionButton("run_reg", "Run Multinomial Regression", class = "btn-success")
          ),
          box(title = "Regression Coefficients & Odds Ratios", width = 8, status = "info",
            verbatimTextOutput("reg_results")
          )
        )
      ),

      # TAB 5: Longitudinal LPTA
      tabItem(tabName = "lpta_module",
        fluidRow(
          box(title = "LPTA Parameters", width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("lpta_t1_select"),
            uiOutput("lpta_t2_select"),
            numericInput("lpta_k", "Number of Profiles (K):", value = 3, min = 2, max = 5),
            actionButton("run_lpta_trans", "Run Transition Analysis", class = "btn-success", icon = icon("random"))
          ),
          box(title = "Profile Transition Matrix", width = 8, status = "info",
            DTOutput("lpta_transition_table")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  raw_data <- reactive({
    req(input$file)
    df <- read.csv(input$file$datapath, stringsAsFactors = FALSE)
    return(df)
  })
  
  output$raw_data_table <- renderDT({
    req(raw_data())
    datatable(raw_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$lpa_var_select <- renderUI({
    req(raw_data())
    selectInput("lpa_vars", "Select LPA Indicator Variables (at least 2):", 
                choices = names(raw_data()), multiple = TRUE)
  })
  
  lpa_res <- eventReactive(input$run_lpa, {
    req(input$lpa_vars)
    validate(need(length(input$lpa_vars) >= 2, "Please select at least 2 indicator variables."))
    
    df_working <- raw_data()
    for(col in input$lpa_vars) {
      df_working[[col]] <- as.numeric(df_working[[col]])
    }
    
    df_clean <- df_working[complete.cases(df_working[, input$lpa_vars, drop = FALSE]), , drop = FALSE]
    lpa_input <- df_clean[, input$lpa_vars, drop = FALSE]
    
    models_enum <- lpa_input %>%
      single_imputation() %>%
      estimate_profiles(1:input$max_k)
    
    fit_table <- get_fit(models_enum)
    selected_model <- models_enum[[input$max_k]]
    class_preds <- get_data(selected_model)$Class
    df_clean$LPA_Class <- as.factor(class_preds)
    
    list(models = models_enum, selected_model = selected_model, data_with_class = df_clean, fit_table = fit_table)
  })
  
  output$lpa_plot <- renderPlot({
    req(lpa_res())
    df <- lpa_res()$data_with_class
    vars <- input$lpa_vars
    
    plot_df <- df
    if(input$scale_vars) {
      plot_df[vars] <- lapply(plot_df[vars], function(x) as.numeric(scale(x)))
    }
    
    profile_means <- plot_df %>%
      group_by(LPA_Class) %>%
      summarise(across(all_of(vars), mean, .names = "{.col}")) %>%
      pivot_longer(cols = -LPA_Class, names_to = "Variable", values_to = "Mean")
    
    ggplot(profile_means, aes(x = Variable, y = Mean, group = LPA_Class, color = LPA_Class)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5) +
      theme_minimal(base_size = 14) +
      labs(title = paste(input$max_k, "Profile Model Plot"), x = "Variables", y = "Mean") +
      theme(legend.position = "bottom")
  })
  
  output$lpa_compare_table <- renderDT({
    req(lpa_res())
    fit_df <- lpa_res()$fit_table
    fit_df$Model <- paste(fit_df$classes_number, "Profile(s)")
    datatable(fit_df[, c("Model", "AIC", "BIC", "SABIC", "Entropy")], options = list(dom = 't'), rownames = FALSE) %>%
      formatRound(columns = c("AIC", "BIC", "SABIC", "Entropy"), digits = 3)
  })
  
  output$lpa_counts_table <- renderDT({
    req(lpa_res())
    df <- lpa_res()$data_with_class
    counts_df <- df %>% group_by(LPA_Class) %>% summarise(Count = n(), Proportion = paste0(round(n()/nrow(df)*100, 2), " %"))
    datatable(counts_df, options = list(dom = 't'), rownames = FALSE)
  })
  
  output$lpa_prob_matrix_table <- renderDT({
    req(lpa_res())
    data_preds <- get_data(lpa_res()$selected_model)
    prob_cols <- grep("^Cprob", names(data_preds), value = TRUE)
    if(length(prob_cols) > 0) {
      prob_matrix <- data_preds %>% group_by(Class) %>% summarise(across(all_of(prob_cols), mean))
      datatable(prob_matrix, options = list(dom = 't'), rownames = FALSE) %>% formatRound(columns = -1, digits = 3)
    } else {
      datatable(data.frame(Message = "Not available."))
    }
  })
  
  output$citation_text <- renderText({
    "Senel, E. (2026). LPA Studio: Open source web application for Latent Profile Analysis."
  })
  
  output$demo_var_select <- renderUI({
    req(lpa_res())
    selectInput("demo_var", "Select Comparison Variable:", choices = setdiff(names(lpa_res()$data_with_class), "LPA_Class"))
  })
  
  output$reg_pred_select <- renderUI({
    req(lpa_res())
    selectInput("reg_preds", "Select Predictor Variables:", choices = setdiff(names(lpa_res()$data_with_class), "LPA_Class"), multiple = TRUE)
  })
  
  observeEvent(input$run_demo, {
    output$demo_results <- renderPrint({
      req(lpa_res(), input$demo_var)
      df <- lpa_res()$data_with_class
      formula <- as.formula(paste(input$demo_var, "~ LPA_Class"))
      if(input$demo_type == "continuous") {
        print(summary(aov(formula, data = df)))
      } else {
        print(chisq.test(table(df[[input$demo_var]], df$LPA_Class)))
      }
    })
  })
  
  observeEvent(input$run_reg, {
    output$reg_results <- renderPrint({
      req(lpa_res(), input$reg_preds)
      df <- lpa_res()$data_with_class
      model <- multinom(as.formula(paste("LPA_Class ~", paste(input$reg_preds, collapse = " + "))), data = df, trace = FALSE)
      print(summary(model))
    })
  })

  output$lpta_t1_select <- renderUI({ req(raw_data()); selectInput("lpta_t1_vars", "Time 1 Variables:", choices = names(raw_data()), multiple = TRUE) })
  output$lpta_t2_select <- renderUI({ req(raw_data()); selectInput("lpta_t2_vars", "Time 2 Variables:", choices = names(raw_data()), multiple = TRUE) })
  
  lpta_res <- eventReactive(input$run_lpta_trans, {
    req(input$lpta_t1_vars, input$lpta_t2_vars)
    df <- raw_data()
    m1 <- estimate_profiles(df[, input$lpta_t1_vars], input$lpta_k)
    m2 <- estimate_profiles(df[, input$lpta_t2_vars], input$lpta_k)
    list(t1 = get_data(m1)$Class, t2 = get_data(m2)$Class)
  })

  output$lpta_transition_table <- renderDT({
    req(lpta_res())
    res <- lpta_res()
    datatable(as.data.frame(table(Time1 = res$t1, Time2 = res$t2)), options = list(dom = 't'), rownames = FALSE)
  })
}

shinyApp(ui = ui, server = server)
