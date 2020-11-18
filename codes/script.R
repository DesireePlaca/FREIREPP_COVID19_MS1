##################################################################################
#                                                                                #
# WARNING: set the folder "R_statistical_analysis" as your working directory or  #
# uncomment and run the code on the next line for having access to all objets:   #
#load("script.RData")                                                            #
#                                                                                #
##################################################################################

# SOURCES ----------------------------------------------------------------------
# calling packages and functions from external sources

source(paste(getwd(), "codes", "packages.R", sep = "/"))
source(paste(getwd(), "codes", "functions.R", sep = "/"))


# DATA -------------------------------------------------------------------------
# loading data from spreadsheets

wd_data <- paste(getwd(), "data", sep = "/")

data_files <- list.files(wd_data, pattern = "\\.csv$")
data_paths <- paste(wd_data, data_files, sep = "/")

data <- purrr::map(data_paths[-1], ~ read.table(.x,
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
))
names(data) <- gsub("(\\..+)$", "", data_files[-1])

data_all <- read.table(data_paths[1],
  header = FALSE, sep = "\t", stringsAsFactors = FALSE)


# DATA PROCESSING --------------------------------------------------------------
# changing data sets for analysis

df_sample_details <- data$sample_condition

df_selected_genes <- data$DEGs_Neutrophil
df_selected_genes[, 1] <- toupper(gsub(" +", "", df_selected_genes[, 1]))
df_selected_genes[, 2] <- toupper(gsub(" +", "", df_selected_genes[, 2]))

v_selected_genes <- unname(unlist(data$DEGs_Neutrophil))
v_empty <- which(v_selected_genes == "")
v_selected_genes <- v_selected_genes[-v_empty]

target_genes <- unique(unlist(data$target_genes))
target_neutrophil <- !is.na(match(target_genes, df_selected_genes[, 1]))
target_cytokine <- !is.na(match(target_genes, df_selected_genes[, 2]))

df_expr <- data$expr0
col_names <- df_expr$Genes
row_names <- colnames(df_expr)[-1]
df_expr <- data.table::transpose(df_expr[, -1])
rownames(df_expr) <- row_names
colnames(df_expr) <- col_names

conditions <- unique(df_sample_details$Class)
df_conditions <- lapply(conditions, function(i) {
  id_condition_i <- subset(df_sample_details, Class == i)$SampleName
  have_condition_i <- which(rownames(df_expr) %in% id_condition_i)
  df_condition_i <- df_expr[have_condition_i, ]
  data.frame(df_condition_i, group = i)
})
names(df_conditions) <- conditions

df_all <- data_all
col_names_all <- df_all[, 1]
row_names_all <- df_all[7, -1]
df_all <- data.table::transpose(df_all[, -1])
rownames(df_all) <- row_names_all
colnames(df_all) <- col_names_all

df_r <- df_all[, c("n1_ct", "n1_ct_Class", "Group", "Gender", "Age", "Sample")]
df_r$n1_ct <- as.numeric(df_r$n1_ct)

df_R <- join(
  data.frame(df_expr, Sample = rownames(df_expr)),
  df_r,
  by = "Sample"
)

df_R$Age <- as.numeric(df_R$Age)

df_R <- df_R[, c(
  "n1_ct", "n1_ct_Class", "Group", "Gender", "Age", v_selected_genes
)]

df_R_reg <- df_R[-which(is.na(df_R), arr.ind = TRUE)[, 1], ]

bolean_neutrophil  <- v_selected_genes %in% df_selected_genes[, 1]
bolean_cytokine  <- v_selected_genes %in% df_selected_genes[, 2]

df_CCA_neutrophil <- df_R[, c("Group", v_selected_genes[bolean_neutrophil])]
df_CCA_cytokine <- df_R[, c("Group", v_selected_genes[bolean_cytokine])]

neutrophil_x <- as.matrix(subset(df_CCA_neutrophil, Group == "COVID")[, -1])
cytokine_y <- as.matrix(subset(df_CCA_cytokine, Group == "COVID")[, -1])

X <- scale(log2(neutrophil_x + 1))
Y <- scale(log2(cytokine_y + 1))


# DATA ANALYSIS ----------------------------------------------------------------

#Parameters---------------------------------------------------------------------

set.seed(2020)

n_boot <- 1000

col_neutrophil <- "#5f69f7dc"
col_cytokine <- "#e26e51fb"

n_genes <- 10

facet_label <- c(
  COVID_F = "Female - SC2",
  CONTROL_F = "Female - Neg. SC2",
  COVID_M = "Male - SC2",
  CONTROL_M = "Male - Neg. SC2",
  COVID_young = "Young - SC2",
  CONTROL_young = "Young - Neg. SC2",
  COVID_elderly = "Elderly - SC2",
  CONTROL_elderly = "Elderly - Neg. SC2",
  COVID_low = "Low Load - SC2",
  COVID_high = "High Load - SC2",
  COVID_medium = "Medium Load - SC2",
  CONTROL_high = "Neg. SC2",
  Neutrophil = "Neutrophil mediated immunity",
  Cytokine = "Cytokine-mediated signaling pathway"
)

sel_cond <- seq_len(length(show_conditions()))[-5]

# Neutrophils ------------------------------------------------------------------

# segment plot:
stat_neutrophil <- summary_groups(sel_cond, "neutrophil")

neutrophil_median <- plot_segments(genes = "neutrophil", var = "median",
  stat_groups = stat_neutrophil, save = TRUE, facet_label = facet_label,
  x_title = "Median gene expression", xoff = -35, yoff = 200)

neutrophil_iqr <- plot_segments(genes = "neutrophil", var = "iqr",
  stat_groups = stat_neutrophil, save = TRUE, facet_label = facet_label,
  x_title = "Interquartile range of gene expression", xoff = -30, yoff = 800)

# non-parametric manova:
neutrophil_g1 <- merge_conditions(show_conditions()[1:4], genes = "neutrophil")

neutrophil_g1 <- stat_supergroups(neutrophil_g1,
  supergroups = list(
    Covid = c("COVID_F", "COVID_M"),
    Control = c("CONTROL_F", "CONTROL_M")
  )
)
n_neutrophil <- nrow(neutrophil_g1)

list_boot_neutrophil <- list()
for (i in 0:n_boot) {
  if (i == 0) {
    df_sample <- 1:n_neutrophil
  } else {
    df_sample <- sample(1:n_neutrophil, n_neutrophil, replace = TRUE)
  }

manova_neutrophil <- nonpartest(
  ABCA13 | ACAA1 | ACTR1B | AGA | AGPAT2 | ALDH3B1 | ALDOC | AMPD3 |
    ANPEP | AP2A2 | APEH | APRT | ATP6V0C | ATP8B4 | CANT1 | CAT | CD14 |
    CD59 | CEP290 | CKAP4 | COMMD3 | CST3 | CTSD | CTSH | CXCL8 | CYBA | DBNL |
    DDOST | DGAT1 | DNAJC5 | DNASE1L1 | DPP7 | EEF2 | ENPP4 | FOS | FUCA1 |
    GAA | GALNS | GSDMD | GSTP1 | GUSB | HMOX2 | HSPA1A | HSPA1B | IMPDH1 |
    IMPDH2 | JUN | JUNB | LAMP1 | LAMTOR1 | LRG1 | LTF | MAN2B1 | METTL7A |
    MIF | MLEC | MVP | NAPRT | NBEAL2 | NDUFC2 | NIT2 | NME2 | ORMDL3 | PDAP1 |
    PDXK | PFKL | PGM1 | PGRMC1 | PKP1 | PPIE | PRKCD | PSMA2 | PSMB1 | PSMC3 |
    PTGES2 | PTPN6 | PTPRN2 | PYCARD | PYGB | QSOX1 | RAB24 | RAB4B | S100A9 |
    S100P | SERPINB6 | SLC27A2 | SLPI | STK11IP | STXBP2 | SVIP | SYNGR1 |
    TCIRG1 | TMC6 | TMEM30A | TOLLIP | TOM1 | TUBB4B | TXNDC5 |
    VAT1 ~ supergroup,
  data = neutrophil_g1[df_sample, ], permtest = FALSE, plots = FALSE,
  permreps = 1000, tests = c(1, 0, 0, 0)
)

  if (i == 0) obs_neutrophil <- manova_neutrophil

  list_boot_neutrophil[[i + 1]] <- manova_neutrophil

  cat("ANOVA/ boot/ neutrophil/ i = ", i, "\n")
}

releffect_neutrophil <- plot_releffects(list_boot_neutrophil,
  genes = "neutrophil", save = TRUE)

# Cytokines --------------------------------------------------------------------

# segment plot:
stat_cytokine <- summary_groups(sel_cond, "cytokine")

cytokine_median <- plot_segments(genes = "cytokine", var = "median",
  stat_groups = stat_cytokine, save = TRUE, facet_label = facet_label,
  x_title = "Median gene expression", xoff = -25, yoff = 65)

cytokine_iqr <- plot_segments(genes = "cytokine", var = "iqr",
  stat_groups = stat_cytokine, save = TRUE, facet_label = facet_label,
  x_title = "Interquartile range of gene expression", xoff = -35, yoff = 75)

# non-parametric manova:
cytokine_g1 <- merge_conditions(show_conditions()[1:4], genes = "cytokine")

cytokine_g1 <- stat_supergroups(cytokine_g1,
  supergroups = list(
    Covid = c("COVID_F", "COVID_M"),
    Control = c("CONTROL_F", "CONTROL_M")
  )
)
n_cytokine <- nrow(cytokine_g1)

list_boot_cytokine <- list()
for (i in 0:n_boot) {
  if (i == 0) {
    df_sample <- 1:n_cytokine
  } else {
    df_sample <- sample(1:n_cytokine, n_cytokine, replace = TRUE)
  }

  manova_cytokine <- nonpartest(
    ADAR | AIM2 | BCL2 | BIRC3 | BST2 | BTK | CASP1 | CCL2 | CCL22 | CCL3 |
      CCL4 | CCL5 | CCR1 | CCR5 | CD300LF | CD36 | CD4 | CD40LG | CD80 | CD86 |
      CEACAM1 | CIITA | CNTF | CSF1 | CSF1R | CSF2RA | CSF2RB | CXCL10 |
      CXCL11 | CXCL13 | CXCL6 | CXCL9 | CXCR1 | CXCR2 | DUOX2 | F13A1 |
      FASLG | FCGR1A | FCGR1B | FLRT2 | FYN | GBP1 | GBP2 | IFI30 | IFI35 |
      IFI6 | IFIT1 | IFIT2 | IFIT3 | IFITM1 | IFITM2 | IFITM3 | IL10RA |
      IL12RB1 | IL15RA | IL18BP | IL18RAP | IL1B | IL1RAP | IL1RN |
      IL20RB | IL21R | IL2RB | IL2RG | IL7R | INPP5D | IRAK2 | IRAK3 |
      IRF4 | IRF7 | IRF8 | ISG15 | ISG20 | ITGAM | ITGAX | JAK3 |
      LCP1 | LEPR | LMNB1 | LRRK2 | MAP2K6 | MSN | MTAP | MX1 | MX2 |
      NANOG | NR4A3 | OAS1 | OAS2 | OAS3 | OASL | OPRD1 | OPRM1 | PML |
      PSMB9 | RNASEL | RSAD2 | SOCS2 | ST18 | STAT2 | STAT4 | TEC |
      TNFRSF9 | TNFSF13B | TNFSF14 | TNFSF15 | TNFSF8 | TRIM21 |
      TRIM31 | TRIM34 | TRIM5 | XAF1 | ZEB1 | ZNF675 ~ supergroup,
    data = cytokine_g1[df_sample, ], permtest = FALSE, plots = FALSE,
    permreps = 1000, tests = c(1, 0, 0, 0)
  )

  if (i == 0) obs_cytokine <- manova_cytokine

  list_boot_cytokine[[i + 1]] <- manova_cytokine

  cat("ANOVA/ boot/ cytokine/ i = ", i, "\n")
}

releffect_cytokine <- plot_releffects(list_boot_cytokine, genes = "cytokine",
  save = TRUE, width = 14)

# Neutrophils x Cytokines -----------------------------------------------------

# pca aiming regression covariates:

var_exclude <- which(colnames(df_R_reg) %in% c(
  target_genes, "Age", "Group", "n1_ct", "n1_ct_Class", "Gender"
))

df_pca_reg <- log2(
  subset(df_R_reg, Gender != "" & Group == "COVID")[, -var_exclude] + 1
)

pca_covar <- principal(df_pca_reg,
  nfactors = 1, rotate = "varimax",
  scores = TRUE
)

genes_covar <- pca_covar$scores
colnames(genes_covar) <- gsub("R", "P", colnames(genes_covar))

load_genes <- loadings(pca_covar)[, 1]
numeric_load <- load_genes
load_genes[load_genes < 0.7] <- 0
load_genes[load_genes >= 0.7] <- 1

type.genes <- ifelse(
  colnames(df_pca_reg) %in% df_selected_genes[, 1],
  "Neutrophil mediated immunity",
  "Cytokine-mediated signaling pathway"
)
name.genes <- colnames(df_pca_reg)
name.genes[load_genes == 0] <- ""

k <- 0.7
cor_genes <- cor(df_pca_reg)
cor_genes_k <- cor_genes
cor_genes_k[abs(cor_genes_k) < k] <- 0
cor_genes_k[abs(cor_genes_k) >= k] <- 1

net <- network(cor_genes_k, directed = FALSE)
net %v% "type" <- type.genes
net %v% "load" <- numeric_load
net %v% "name" <- name.genes

set.seed(100)
p1_pca <- ggplot(data = ggnetwork(net, layout = "kamadakawai"),
  aes(x, y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.8, color = "grey") +
  geom_nodes(aes(color = type, size = load^3), alpha = 0.5,
    show.legend = FALSE) +
  geom_nodetext(aes(label = name, size = load^7), fontface = "bold",
    show.legend = FALSE) +
  scale_colour_manual("", values = c(
    "Neutrophil mediated immunity" = col_neutrophil,
    "Cytokine-mediated signaling pathway" = col_cytokine)) +
  theme_blank()

  save_plot(TRUE,
      "PCA_LogData",
      "Frist-Principal-Component-PC1_For-Interpretation-Of-PC1", "",
      p1_pca,
      width = 12, height = 7
  )

df_mlm <- cbind(
  subset(
    data.frame(
      df_R_reg[, c("Group", "Gender")],
      scale(df_R_reg[, c("Age", "n1_ct")]),
      log2(1 + df_R_reg[, target_genes])
    ), Gender != "" & Group == "COVID"
  )[, -1],
  genes_covar
)
df_mlm$Gender <- ifelse(df_mlm$Gender == "F", 1, 0)

# cca (DFA.CANCOR):

signed_cor <- whitening::cca(X, Y)$lambda

cca3 <- new_CANCOR(data.frame(X, Y), colnames(X), colnames(Y), FALSE, 1,
  "structure", FALSE)

cca_xy_12 <- plot_cca(cancor = cca3, n = 1, save = TRUE, width = 13,
  plot_type = 1, k_min = 0.7)

k <- 0.7
seed <- 1
net_cca_2 <- network_cca(X, Y, cca3, signed_cor, 2, k, seed, 2.2, save = TRUE)

# ridge regression -----------------------------------------------------------

set.seed(2020)

var_exclude <- which(colnames(df_R_reg) %in% c(
  "n1_ct", "Age", "Gender", "Group", "n1_ct_Class")
)

df_ridge_reg <- cbind(
  subset(df_R_reg, n1_ct > 0 & Gender != "")["n1_ct"],
  scale(subset(df_R_reg, n1_ct > 0 & Gender != "")["Age"]),
  subset(df_R_reg, n1_ct > 0 & Gender != "")["Gender"],
  log2(subset(df_R_reg, n1_ct > 0 & Gender != "")[, -var_exclude] + 1)
)
df_ridge_reg$Gender <- ifelse(df_ridge_reg$Gender == "F", 1, 0)

n_obs <- nrow(df_ridge_reg)
index_sample <- sample(seq_len(n_obs), floor((2 / 3) * n_obs), replace = FALSE)
df_train <- df_ridge_reg[index_sample, ]
df_test <- df_ridge_reg[-index_sample, ]

fit_ridge <- cv.glmnet(as.matrix(df_train[, -1]), df_train[, "n1_ct"],
  type.measure = "mse", alpha = 0, family = "gaussian")

predict_ridge <- predict(fit_ridge, s = fit_ridge$lambda.min,
  newx = as.matrix(df_test[, -1]))

mse_ridge <- mean((df_test[, "n1_ct"] - c(predict_ridge))^2)

genes_nonzero_coef <- coef(fit_ridge)[-(1:3), ]
genes_nonzero_name <- names(genes_nonzero_coef)
genes_nonzero_name <- factor(genes_nonzero_name, ordered = TRUE,
  levels = genes_nonzero_name[order(genes_nonzero_coef)])

df_plot_ridge <- data.frame(
  names = genes_nonzero_name,
  coefs = genes_nonzero_coef,
  genes = ifelse(genes_nonzero_name %in% df_selected_genes[, 1],
    "Neutrophil", "Cytokine")
)

p_ridge_1 <- ggplot(data = df_plot_ridge, aes(xend = names, x = names,
  yend = coefs, y = 0, colour = genes)) +
  geom_hline(yintercept = 0, alpha = 0.4) +
  geom_segment(show.legend = FALSE) +
  geom_point(aes(y = coefs), show.legend = FALSE) +
  theme_classic() +
  scale_colour_manual("", values = c(Neutrophil = col_neutrophil,
    Cytokine = col_cytokine)) +
  facet_wrap(~ genes, ncol = 1, scales = "free_x",
    labeller = as_labeller(facet_label)) +
  theme(
    legend.text = element_text(size = 11),
    axis.text.x = element_text(angle = 90, size = 7),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 13),
    panel.grid.major = element_line(),
    panel.grid.minor.y = element_line(),
    strip.text = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(-0.07, 0.07), breaks = seq(-0.1, 0.1, 0.02),
    labels = scales::number_format(accuracy = 0.001, decimal.mark = '.')
  ) +
  labs(y = "Regression coefficient estimate", x = "Genes")

p_ridge_2 <- ggplot(data = data.frame(y_test = df_test[, "n1_ct"],
  pred = as.vector(predict_ridge)), aes(x = pred, y = y_test)) +
  geom_abline(slope = 1, intercept = 0, col = "gold", alpha = 0.6, size = 1.2) +
  geom_point(alpha = 0.6, size = 3) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.line.y = element_line(),
    axis.title = element_text(size = 13)
  ) +
  scale_y_continuous(limits = c(15, 30)) +
  scale_x_continuous(limits = c(15, 30)) +
  labs(y = "Observed response from test set", x = "Model prediction")

save_plot(
    save = TRUE,
    "RidgeReg_Log-Data",
    "Regression-Coefficients-Estimates", "",
    p_ridge_1, 10, 6
)

save_plot(
    save = TRUE,
    "RidgeReg_Log-Data",
    "Measure-Of-Prediction-Accuracy", "",
    p_ridge_2, 10, 6
)

### Multivariate Regression-------------------------------------------------------------------

formula <- cbind(
  CD14, CXCR1, CXCR2, IL1B, NBEAL2, S100A9
) ~ Gender + Age + n1_ct + PC1 + I(Gender * Age) +
    I(Gender * n1_ct) + I(Gender * PC1)

data <- df_mlm

##### Fit:

fit <- lm.t(formula, data, dDist = "dmvn", k = log2(nrow(data)), PLOT = FALSE)

Est <- lapply(fit$Est$Summ, round, digits = 3)

over_MD <- sum((fit$Diagnostics$MD > fit$Diagnostics$cut.MD)) / nrow(data)

##### Hypothesis Testing:

n_var <- ncol(fit$Matrices$B)
n_covar <- nrow(fit$Matrices$B)

D <- diag(n_var)
D <- sapply(seq_len(n_var - 1), function(i) D[, i] - D[, i + 1])

set.seed(123)

###### B:
fit.fast <- lm.t.fast(formula, data, dDist = "dmvn")
B_coef <- boot_B(fit = fit.fast, PLOT = FALSE)

###### VCOV-HC:
fit_mlm <- lm(formula, data)
sig <- sqrt(diag(vcovHC(fit_mlm)))

###### Prediction:

n_boot <- 1
df_plot_reg <- df_mlm[, c("Age", "n1_ct", "PC1")]
pred_Age <- prediction_data("Age", "Age (scaled)", df_plot_reg, x_pos = -2.8, y_pos = 13, save = TRUE)
pred_n1_ct <- prediction_data("n1_ct", "Viral load (scaled)", df_plot_reg, x_pos = -.95, y_pos = 13, save = TRUE)
pred_PC1 <- prediction_data("PC1", "(PC1) Principal component score based on the remaining genes", df_plot_reg, x_pos = -2.2, y_pos = 13, save = TRUE)

###### Omnibus test:

n_boot <- 1000

Omnibus <- boot_H0(line = 0, PLOT = FALSE, fit = fit.fast)
Omnibus$model$fit$Matrices$B
Omnibus$model$fit.0$Matrices$B
Omnibus$model$fit$Test
Omnibus$result$p.boot
Omnibus$result$p.obs

###### Intercept:
Intercept <- boot_H0(line = 1, PLOT = FALSE, fit = fit.fast)
Intercept$model$fit$Matrices$B
Intercept$model$fit.0$Matrices$B
Intercept$model$fit$Test
Intercept$result$p.boot
Intercept$result$p.obs

###### Gender:
Gender <- boot_H0(line = 2, PLOT = FALSE, fit = fit.fast)
Gender$model$fit$Matrices$B
Gender$model$fit.0$Matrices$B
Gender$model$fit$Test
Gender$result$p.boot
Gender$result$p.obs

###### Age:
Age <- boot_H0(line = 3, PLOT = FALSE, fit = fit.fast)
Age$model$fit$Matrices$B
Age$model$fit.0$Matrices$B
Age$model$fit$Test
Age$result$p.boot
Age$result$p.obs

###### n1_ct:
n1_ct <- boot_H0(line = 4, PLOT = FALSE, fit = fit.fast)
n1_ct$model$fit$Matrices$B
n1_ct$model$fit.0$Matrices$B
n1_ct$model$fit$Test
n1_ct$result$p.boot
n1_ct$result$p.obs

###### PC1:
PC1 <- boot_H0(line = 5, PLOT = FALSE, fit = fit.fast)
PC1$model$fit$Matrices$B
PC1$model$fit.0$Matrices$B
PC1$model$fit$Test
PC1$result$p.boot
PC1$result$p.obs

###### Gender * Age:
GenderAge <- boot_H0(line = 6, PLOT = FALSE, fit = fit.fast)
GenderAge$model$fit$Matrices$B
GenderAge$model$fit.0$Matrices$B
GenderAge$model$fit$Test
GenderAge$result$p.boot
GenderAge$result$p.obs

###### Gender * n1_ct:
Gendern1_ct <- boot_H0(line = 7, PLOT = FALSE, fit = fit.fast)
Gendern1_ct$model$fit$Matrices$B
Gendern1_ct$model$fit.0$Matrices$B
Gendern1_ct$model$fit$Test
Gendern1_ct$result$p.boot
Gendern1_ct$result$p.obs

###### Gender * PC1:
GenderPC1 <- boot_H0(line = 8, PLOT = FALSE, fit = fit.fast)
GenderPC1$model$fit$Matrices$B
GenderPC1$model$fit.0$Matrices$B
GenderPC1$model$fit$Test
GenderPC1$result$p.boot
GenderPC1$result$p.obs

##### Diagnostics:

n <- nrow(data)
MD <- data.frame(x = seq_len(n), y = fit$Diagnostics$MD, type = "MD",
  cut1 = c(fit$Diagnostics$cut.MD, rep(NA, n -1)),
  cut2 = c(fit$Diagnostics$cut.MD, rep(NA, n -1)))
SCR <- data.frame(x = seq_len(n), y = fit$Diagnostics$SCRi, type = "SCR",
  cut1 = NA, cut2 = NA)
StudRes <- data.frame(x = rep(seq_len(n), n_var), 
  y = c(fit$Diagnostics$stud.res), type = "StudRes",
  cut1 = c(-2, rep(NA, n - 1)), cut2 = c(2, rep(NA, n - 1)))
StudFit <- data.frame(x = c(fit$Diagnostics$XB), 
  y = c(fit$Diagnostics$stud.res), type = "StudFit",
  cut1 = c(-2, rep(NA, n - 1)), cut2 = c(2, rep(NA, n - 1)))
AgeRes <- data.frame(x = rep(data[, "Age"], n_var),
  y = c(fit$Diagnostics$stud.res), type = "AgeRes",
  cut1 = c(-2, rep(NA, n - 1)), cut2 = c(2, rep(NA, n - 1)))
CtRes <- data.frame(x = rep(data[, "n1_ct"], n_var),
  y = c(fit$Diagnostics$stud.res), type = "CtRes",
  cut1 = c(-2, rep(NA, n - 1)), cut2 = c(2, rep(NA, n - 1)))
PC1Res <- data.frame(x = rep(data[, "PC1"], n_var),
  y = c(fit$Diagnostics$stud.res), type = "PC1Res",
  cut1 = c(-2, rep(NA, n - 1)), cut2 = c(2, rep(NA, n - 1)))
GenderRes <- data.frame(x = rep(data[, "Gender"], n_var),
  y = c(fit$Diagnostics$stud.res), type = "GenderRes",
  cut1 = c(-2, rep(NA, n - 1)), cut2 = c(2, rep(NA, n - 1)))

df_diagnostics <- rbind(StudRes, StudFit, MD, SCR, AgeRes, CtRes, PC1Res,
  GenderRes)
df_diagnostics$type <- factor(df_diagnostics$type, ordered = TRUE,
  levels = unique(df_diagnostics$type))

facet_label <- c(
  StudRes = "Studentized Resdiuals Index of Observations",
  StudFit = "Studantized Resdiuals x Fitted Values",
  MD = "Mahalanobis Distance x Index of Observations",
  SCR = "Scale Ratio x Index of observations",
  AgeRes = "Studentized Resdiuals x Age (scaled)",
  CtRes = "Studentized Resdiuals x Viral load (scaled)",
  PC1Res = "Studentized Resdiuals x PC1",
  GenderRes = "Studentized Resdiuals x Gender"
)

p_diag <- ggplot(data = subset(df_diagnostics, !(type %in% c("MD", "SCR"))),
  aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_segment(inherit.aes = FALSE,
    data = subset(df_diagnostics, type %in% c("MD", "SCR")),
    aes(x = x, y = 0, xend = x, yend = y), alpha = 0.7) +
  geom_hline(aes(yintercept = cut1), linetype = "dashed", show.legend = FALSE) +
  geom_hline(aes(yintercept = cut2), linetype = "dashed", show.legend = FALSE) +
  facet_wrap(~ type, ncol = 4, scales = "free",
    labeller = as_labeller(facet_label)) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 11),
    panel.border = element_rect(fill = "transparent"),
    axis.line = element_blank()
  ) + xlab("") + ylab("")

save_plot(TRUE, "MLM_LogData", "Model-Diagnostics", "", p_diag, width = 12, height = 7)


## Tables------------------------------------------------------------------------

#### Estimates:

write.table(round(fit$Est$Cor, 3), file = paste(getwd(), "tables", "MRV_Matrix-Of-Correlation-Estimate-For-Response-Variables.txt", sep = "/"), sep = "\t", quote = FALSE)

D_Est <- do.call(rbind, lapply(Est, round, digits = 3))
D_Est <- cbind(D_Est, rob_sd = round(sig, 3), p_boot = round(c(B_coef$p.boot), 3))
write.table(D_Est, file = paste(getwd(), "tables", "MRV_Coefficient-Estimates.txt", sep = "/"), sep = "\t", quote = FALSE)

#### Hypothesis Testing:

D_Test1 <- do.call(rbind, lapply(list(
   Omnibus = c(Omnibus$model$fit$Test[3], list(p.boot = Omnibus$result$p.boot), Omnibus$model$fit$Test[6]),
   Intercept = c(Intercept$model$fit$Test[3], list(p.boot = Intercept$result$p.boot), Intercept$model$fit$Test[6]),
   Gender = c(Gender$model$fit$Test[3], list(p.boot = Gender$result$p.boot), Gender$model$fit$Test[6]),
   Age = c(Age$model$fit$Test[3], list(p.boot = Age$result$p.boot), Age$model$fit$Test[6]),
   n1_ct = c(n1_ct$model$fit$Test[3], list(p.boot = n1_ct$result$p.boot), n1_ct$model$fit$Test[6]),
   PC1 = c(PC1$model$fit$Test[3], list(p.boot = PC1$result$p.boot), PC1$model$fit$Test[6]),
   GenderAge = c(GenderAge$model$fit$Test[3], list(p.boot = GenderAge$result$p.boot), GenderAge$model$fit$Test[6]),
   Gendern1_ct = c(Gendern1_ct$model$fit$Test[3], list(p.boot = Gendern1_ct$result$p.boot), Gendern1_ct$model$fit$Test[6]),
   GenderPC1 = c(GenderPC1$model$fit$Test[3], list(p.boot = GenderPC1$result$p.boot), GenderPC1$model$fit$Test[6])
  ), data.frame
 )
)
D_Test1 <- data.frame(round(D_Test1[, 1:2], 3), test = D_Test1[, 3])
write.table(D_Test1, file = paste(getwd(), "tables", "MRV_Likelihood-Ratio-Test.txt", sep = "/"), sep = "\t", quote = FALSE)