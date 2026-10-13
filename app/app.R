# ============================================================================
# ctDNA Time-Resolved Absolute Recurrence Risk Calculator (Shiny)
# Engine: data/continuous_curves.csv (M1c one-stage Royston-Parmar continuous
# absolute-risk curves). Only curative-intent scenarios that pass the M4
# dual-evidence reliability gate, with >=3 independent studies and grade-A
# curves, are included.
# For research and educational use only; not individual medical advice.
# PROSPERO CRD420261516792.
# ============================================================================
library(shiny)
library(readr)
library(dplyr)
library(ggplot2)

curves <- read_csv("data/continuous_curves.csv", show_col_types = FALSE)

SET <- list(
  "Non-metastatic CRC"             = "MRD",
  "Locally-advanced rectal cancer" = "LARC",
  "Resectable colorectal liver metastases" = "CRLM"
)
WIN <- list(
  "Pre-operative"                   = "preop",
  "Post-operative MRD (4-12 weeks)" = "postopMRD",
  "After adjuvant/systemic therapy" = "postACT"
)
WIN_FOR_SET <- list(
  MRD  = c("postopMRD", "postACT"),
  LARC = c("postopMRD"),
  CRLM = c("preop", "postopMRD", "postACT")
)
RED <- "#C0392B"; BLUE <- "#2C3E50"

cell <- function(scen, win) curves %>% filter(scen_code == scen, win == win, endp_code == "DFS")
at_t <- function(d, t) {
  d <- d %>% arrange(time_month)
  if (t <= max(d$time_month)) {
    z <- d[which.min(abs(d$time_month - t)), ]
    return(c(CIR = z$CIR, lo = z$CIRlo, hi = z$CIRhi))
  }
  c(NA, NA, NA)
}
fmt <- function(x) ifelse(is.na(x), "--", sprintf("%.0f%%", 100 * x))
fmtci <- function(l, h) ifelse(is.na(l), "--", sprintf("%.0f%%-%.0f%%", 100 * l, 100 * h))

ui <- fluidPage(
  titlePanel(
    div(
      strong("ctDNA Time-Resolved Absolute Recurrence Risk Calculator"),
      div(style = "font-size:13px;color:#555;font-weight:400;",
          "Cumulative incidence of recurrence by ctDNA status, disease setting and assessment window")
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      selectInput("setting", "Disease setting", choices = names(SET)),
      selectInput("window", "Assessment window", choices = NULL),
      radioButtons("status", "ctDNA status",
                   choices = c("ctDNA-positive" = 1, "ctDNA-negative" = 0), inline = TRUE),
      sliderInput("horizon", "Time horizon (months)", min = 1, max = 24, value = 12, step = 1),
      hr(),
      div(style = "font-size:12px;color:#666;",
          "Pooled absolute risks from reconstructed patient-level data, fitted with one-stage
          Royston-Parmar splines. Bands are pointwise (individual-level) confidence intervals.
          For research and educational use only; not for individual clinical decisions.")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel(
          "Risk estimate",
          fluidRow(
            column(
              4,
              wellPanel(
                style = "padding:12px;",
                div(style = "font-size:13px;color:#555;", textOutput("lab_sel")),
                div(style = "font-size:34px;line-height:1.1;", textOutput("cir_sel")),
                div(style = "font-size:12px;color:#555;", textOutput("ci_sel"))
              )
            ),
            column(
              8,
              wellPanel(
                style = "padding:12px;",
                tableOutput("contrast_tbl")
              )
            )
          ),
          plotOutput("curve_plot", height = "430px"),
          h5("Landmark cumulative incidence"),
          tableOutput("landmark_tbl")
        ),
        tabPanel(
          "Methods & reliability",
          h5("How estimates were built"),
          p("Kaplan-Meier curves from 24 published studies were digitised and reconstructed into
            patient-level data using the Guyot iterative algorithm, with interval-conservation
            checks (reconstructed numbers at risk exactly match the published tables; 192/192 curves,
            11,481 reconstructed patients). Absolute risks were pooled with two-stage random-effects
            models and refitted as continuous curves using one-stage Royston-Parmar splines."),
          h5("Reliability gate (why only these scenarios are shown)"),
          p("A parameter-simulation grid (320 cells x 500 replicates) and a closed-loop plasmode
            validation on three published appendix individual-level datasets calibrated the
            reconstruction. Acceptance: mean absolute error <= 0.05, relative bias <= 0.10, coverage
            >= 0.90. Only cells with >=3 independent studies, grade-A curves and complete risk tables
            are exposed. Advanced disease and anti-EGFR rechallenge (<3 studies) are excluded."),
          h5("External calibration"),
          p("Against the published Li 2022 stage III post-chemotherapy cohort, predicted 12-month
            risk was 4.7% vs observed 6.2% (ctDNA-negative) and 62% vs observed 45.5% (ctDNA-positive,
            conservative). The Ogaard advanced-disease cohort was an external anchor, not pooled."),
          h5("Interpretation"),
          p("Positive ctDNA indicates a markedly elevated, front-loaded risk (~80% of recurrences
            accrue within 12 months and ~90% within 24 months); published reports place ctDNA
            detection a median of ~7-10 months ahead of imaging. Negative ctDNA confers a low but
            non-zero risk. Pre-operative ctDNA in resectable liver metastases showed little
            discrimination and should not be used alone."),
          h5("Limitations"),
          p("Confidence bands are individual-level and narrower than random-effects CIs; use the
            pooled random-effects CIs in the manuscript for primary claims. Assays and surveillance
            schedules differ across studies; this is an absolute-risk synthesis, not a covariate
            predictive model.")
        ),
        tabPanel(
          "Data & sources",
          h5("Data statement"),
          p("Reconstruction used published Kaplan-Meier curves and risk tables. The published appendix
            individual-level datasets of Li 2022 and Ogaard 2022 were used for validation only and
            are NOT redistributed; obtain them from the original publications. A reproducible data-raw
            script and the study list with DOIs are provided in the repository."),
          h5("Reuse & registration"),
          p("Protocol: PROSPERO CRD420261516792. Please cite the accompanying manuscript and the
            original contributing studies. This application is for reproducibility and academic use;
            it is not a medical device and does not replace clinical judgement.")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$setting, {
    scen <- SET[[input$setting]]
    wins <- WIN_FOR_SET[[scen]]
    updateSelectInput(session, "window",
                      choices = setNames(wins, names(WIN)[match(wins, WIN)]))
  })

  dat <- reactive({
    scen <- SET[[input$setting]]
    req(input$window)
    d <- cell(scen, input$window)
    req(nrow(d) > 0)
    d
  })
  pos <- reactive(dat() %>% filter(ctdna == 1))
  neg <- reactive(dat() %>% filter(ctdna == 0))

  observeEvent(dat(), {
    tmax <- floor(min(dat() %>% group_by(ctdna) %>%
                        summarise(tm = max(time_month)) %>% dplyr::pull(tm)))
    updateSliderInput(session, "horizon", min = 1, max = max(12, tmax),
                      value = min(input$horizon, tmax))
  })

  output$lab_sel <- renderText({
    sprintf("%s | %s | %d months", input$setting,
            names(WIN)[match(input$window, WIN)], input$horizon)
  })
  output$cir_sel <- renderText({
    d <- if (input$status == 1) pos() else neg()
    fmt(at_t(d, input$horizon)["CIR"])
  })
  output$ci_sel <- renderText({
    d <- if (input$status == 1) pos() else neg()
    a <- at_t(d, input$horizon)
    paste("95% CI", fmtci(a["lo"], a["hi"]))
  })

  output$contrast_tbl <- renderTable({
    h <- input$horizon
    ap <- at_t(pos(), h); an <- at_t(neg(), h)
    rd <- ap["CIR"] - an["CIR"]
    rr <- if (!is.na(an["CIR"]) && an["CIR"] > 0.001) ap["CIR"] / an["CIR"] else NA
    data.frame(
      `At selected horizon` = c("ctDNA-positive", "ctDNA-negative", "Risk difference", "Risk ratio"),
      Estimate = c(fmt(ap["CIR"]), fmt(an["CIR"]),
                   sprintf("%+.0f pp", 100 * rd),
                   ifelse(is.na(rr), "--", sprintf("%.1f-fold", rr))),
      `95.CI` = c(fmtci(ap["lo"], ap["hi"]), fmtci(an["lo"], an["hi"]), "", ""),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)

  output$landmark_tbl <- renderTable({
    rows <- lapply(c(1, 0), function(s) {
      d <- if (s == 1) pos() else neg()
      lab <- if (s == 1) "ctDNA-positive" else "ctDNA-negative"
      data.frame(
        Group = lab,
        `12 mo` = fmt(at_t(d, 12)["CIR"]),
        `24 mo` = fmt(at_t(d, 24)["CIR"]),
        `36 mo` = fmt(at_t(d, 36)["CIR"]),
        check.names = FALSE
      )
    })
    bind_rows(rows)
  }, striped = TRUE, hover = TRUE)

  output$curve_plot <- renderPlot({
    d <- dat() %>% mutate(grp = ifelse(ctdna == 1, "ctDNA-positive", "ctDNA-negative"))
    h <- input$horizon
    ds <- if (input$status == 1) pos() else neg()
    a <- at_t(ds, h)
    ggplot(d, aes(x = time_month, y = CIR * 100, color = grp, fill = grp)) +
      geom_ribbon(aes(ymin = CIRlo * 100, ymax = CIRhi * 100), alpha = 0.12, linewidth = 0,
                  show.legend = FALSE) +
      geom_line(linewidth = 1.2) +
      geom_vline(xintercept = h, linetype = "dashed", color = "grey50") +
      annotate("point", x = h, y = 100 * a["CIR"], color = "black", size = 3) +
      annotate("text", x = h, y = 100 * a["CIR"] + 4, label = fmt(a["CIR"]),
               size = 4.5, fontface = "bold") +
      scale_color_manual(values = c("ctDNA-positive" = RED, "ctDNA-negative" = BLUE), name = NULL) +
      scale_fill_manual(values = c("ctDNA-positive" = RED, "ctDNA-negative" = BLUE), name = NULL) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(x = "Months since landmark", y = "Cumulative incidence of recurrence, CIR (%)",
           title = paste(input$setting, "-", names(WIN)[match(input$window, WIN)])) +
      theme_bw(base_size = 14) +
      theme(plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"),
            legend.position = "bottom", text = element_text(family = "sans"))
  })
}

shinyApp(ui, server)
