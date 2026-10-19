suppressMessages({library(shiny); library(readr); library(dplyr); library(ggplot2)})
app <- source("app.R")$value
cases <- list(c("Non-metastatic CRC","postACT"),
              c("Non-metastatic CRC","postopMRD"),
              c("Locally-advanced rectal cancer","postopMRD"),
              c("Locally-advanced rectal cancer","postNAC"),
              c("Resectable colorectal liver metastases","postopMRD"),
              c("Resectable colorectal liver metastases","preop"),
              c("Resectable colorectal liver metastases","postACT"))
for (sw in cases) {
  testServer(app, {
    session$setInputs(setting = sw[1], window = sw[2], status = 1, horizon = 12)
    h <- output$kflag
    txt <- if (is.null(h)) "NULL" else paste(as.character(h), collapse = " ")
    txt <- gsub("<[^>]+>", " ", txt)
    txt <- gsub("\\s+", " ", txt)
    cat(sw[1], "|", sw[2], "->", substr(trimws(txt), 1, 78), "\n")
  })
}
cat("TEST DONE\n")
