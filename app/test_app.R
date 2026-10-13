cat("shiny installed:", requireNamespace("shiny", quietly=TRUE), "\n")
p <- try(parse(file="app.R"), silent=TRUE)
cat("app.R parse OK:", inherits(p, "expression"), "\n")
suppressPackageStartupMessages({library(readr);library(dplyr)})
curves <- read_csv("data/continuous_curves.csv", show_col_types=FALSE)
at_t <- function(d,t){ d <- arrange(d, time_month)
  if (t <= max(d$time_month)){ z <- d[which.min(abs(d$time_month-t)),]; c(CIR=z$CIR,lo=z$CIRlo,hi=z$CIRhi)
  } else c(NA,NA,NA) }
cells <- distinct(curves, scen_code, win)
for (i in seq_len(nrow(cells))){
  ce <- cells[i,]
  d <- filter(curves, scen_code==ce$scen_code, win==ce$win, endp_code=="DFS")
  pp <- filter(d, ctdna==1); nn <- filter(d, ctdna==0)
  cat(ce$scen_code, ce$win,
      "| pos12", sprintf("%.0f%%",100*at_t(pp,12)["CIR"]),
      "neg12", sprintf("%.0f%%",100*at_t(nn,12)["CIR"]),
      "| pos24", sprintf("%.0f%%",100*at_t(pp,24)["CIR"]),
      "neg24", sprintf("%.0f%%",100*at_t(nn,24)["CIR"]), "\n")
}
