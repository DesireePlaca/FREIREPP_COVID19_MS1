# the main packages are DFA.CANCOR, psych, glmnet, whitening and ggplot2 (and its extensions).

# PACKAGES ---------------------------------------------------------------------

load_lib <- c(
  # data processing:
  "tools",
  "purrr",
  "data.table",
  "plyr",
  "dplyr",
  "Formula",
  "rlist",
  "car",
  "MASS",
  # graphs:
  "ggplot2",
  "ggrepel",
  "gridExtra",
  "reshape2",
  "scales",
  "qgraph",
  "corrplot",
  "RColorBrewer",
  "GGally",
  "sna",
  "network",
  "ggnetwork",
  "stringr",
  # descriptive statistics:
  "EnvStats",
  "CCA",
  "whitening",
  "DFA.CANCOR",
  "pls",
  "psych",
  # inferencial statistics:
  "npmv",
  "CCP",
  "glmnet",
  "boot",
  "ppcor",
  "splines",
  "LaplacesDemon",
  "sandwich",
  "heplots",
  # math:
  "expm"
)

install_lib <- load_lib[!load_lib %in% installed.packages()]
for (lib in install_lib) install.packages(lib, dependencies = TRUE)
sapply(load_lib, require, character = TRUE)