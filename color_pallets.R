source("git_modules/scRNA_tools/scRNA_functions.R")

#color palette for the morphotypes (from L. major atlas paper - Carolina M. C. Catta-Preta et al 2024)
morphotype_colors <-  c(
  amastigotes = "#669b4d",
  procyclics = "#8b5b95",
  `procyclics-to-nectomonads` = "#47a3c1",
  nectomonads = "#a94989",
  leptomonads = "#bd6831",
  haptomonads = "#4c80b3",
  `early-metacyclics` = "#42a06e",
  `late-metacyclics` = "#e9e23a"
)

#replace it by more diverse, color-blind friendly colors
morphotype_colors <- create_colors(c(names(morphotype_colors), "gamete-like"), palette = "colorblind friendly")

#color palette for the cell cycle phases (from T. brucei scRNA paper)
cell_cycle_colors <- c(
  early_g1 = "#dc746d",
  late_g1 = "#6bbdc1",
  s_phase = "#b47bfb",
  g2m = "#8bad28",
  non_cycling = "#a9a9a8"
)
