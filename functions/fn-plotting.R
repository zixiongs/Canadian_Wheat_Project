# Covert P-values to significance levels

sig_level <- function(p.value) {
  case_when(
    p.value > 0.05 ~ "ns",
    p.value <= 0.001 ~ "***",
    p.value <= 0.01 ~ "**",
    p.value <= 0.05 ~ "*"
  )
}