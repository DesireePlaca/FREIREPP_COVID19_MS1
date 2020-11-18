# these are not pure functions, it depends on global variables.

# FUNCTIONS --------------------------------------------------------------------

show_conditions <- function() unique(df_sample_details$Class)

merge_conditions <- function(conditions, genes = "all") {
  if (genes %in% c("all", "neutrophil", "cytokine", "other")) {
    df <- do.call(rbind, df_conditions[conditions])
    k <- ncol(df)

    if (genes == "neutrophil") {
      genes <- which(colnames(df[, -k]) %in% df_selected_genes[, 1])
      df <- df[, c(genes, k)]
    } else if (genes == "cytokine") {
      genes <- which(colnames(df[, -k]) %in% df_selected_genes[, 2])
      df <- df[, c(genes, k)]
    } else if (genes == "other") {
      genes <- which(!(colnames(df[, -k]) %in% unlist(df_selected_genes)))
      df <- df[, c(genes, k)]
    }

    return(df)
  } else {
    warnings("choose genes = c('all', 'neutrophil', 'cytokine', 'other')")
  }
}

labeller <- function(df, var_value, var_label, size = 5) {
  variable <- as.character(df[, var_label])
  value <- df[, var_value]
  sort_value <- rev(order(value))
  positions <- which(value %in% value[sort_value[1:size]])
  to_label <- variable[positions]
  labels <- rep(NA, length(variable))
  labels[positions] <- to_label

  return(labels)
}

set_col_label <- function(df, var, n_genes) {
  df$lab_var <- unlist(
    lapply(unique(df$group), function(i) {
      data_group <- subset(df, group == i)
      labeller(data_group, var, "variable", size = n_genes)
    })
  )
  return(df)
}

summary_groups <- function(conditions, genes) {
  df <- merge_conditions(conditions, genes = genes)

  melt_df <- reshape2::melt(df, id.vars = "group")

  if (ncol(df) < 500) {
    descr_df <- ddply(melt_df, ~ group + variable, summarise,
      mean = mean(value), median = median(value), sd = sd(value),
      iqr = iqr(value))

    output <- list(df = df, descr = descr_df)
    return(output)
  } else {
    stop("number of genes exceeds 500. STOP!")
  }
}

save_plot <- function(save, genes, var_name, append, p, width, height) {
  if (save) {
    pdf_name <- paste0(genes, "_", var_name, append, ".pdf")
    png_name <- paste0(genes, "_", var_name, append, ".png")

    pdf(paste(getwd(), "figures", pdf_name, sep = "/"),
      width = width, height = height
    )
      print(p)
    dev.off()

    png(paste(getwd(), "figures", png_name, sep = "/"),
      width = width, height = height, units = "in", res = 300
    )
      print(p)
    dev.off()

    print("plot saved!")
  } else {
    print("plot NOT saved.")
  }
}

plot_segments <- function(genes, var, stat_groups, save = FALSE, width = 12,
  height = 7, append = "", n_genes = 10, facet_label = NULL, col = NULL,
  x_title = "statistics", xoff = -20, yoff = 100) {
  df <- stat_groups$df
  descr_df <- stat_groups$descr

  descr_df <- set_col_label(descr_df, var, n_genes)

  descr_output <- descr_df
  var_name <- var
  colnames(descr_df)[which(colnames(descr_df) == var)] <- "var"

  if (is.null(col)) {
    if (genes == "neutrophil") {
      col  <- col_neutrophil
    } else if (genes == "cytokine") {
      col <- col_cytokine
    } else if (genes == "other") {
      col <- "pink"
    } else {
      col <- "gold"
    }
  }

  descr_df <- arrange(descr_df, group, var) 

  descr_df$variable <- factor(descr_df$variable, ordered = TRUE,
    levels = unique(subset(descr_df, group == "COVID_F")$variable))

  descr_df$group <- factor(descr_df$group, ordered = TRUE,
    levels = names(facet_label))

  descr_df <- subset(descr_df, group != "COVID_medium") 

  xx_title <- ifelse(tolower(genes) == "neutrophil",
    "Neutrophil mediated immunity genes",
    "Cytokine-mediated signaling pathway genes")

  p <- ggplot(data = descr_df, aes(x = variable, y = 0, xend = variable,
  yend = var, group = group)) +
  geom_segment(colour = col) +
  geom_point(aes(y = ifelse(!is.na(lab_var), var, NA)), size = 1) +
  geom_text_repel(aes(y = var, label = lab_var), size = 2.9, force = 20,
    segment.alpha = 0.4, segment.size = 0.3, min.segment.length = 0,
    nudge_x = xoff, nudge_y = yoff, fontface = "bold") +
  theme_classic() +
  theme(
    legend.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 10)
  ) +
  labs(y = x_title, x = xx_title)

  if (is.null(facet_label)) {
    p <- p +  facet_wrap(~ group)
  } else {
    p <- p +  facet_wrap(~ group, labeller = as_labeller(facet_label))
  }

  save_plot(
      save,
      "Desc_Raw-Data",
      paste0(str_to_title(genes), "-"), gsub(" +", "-", x_title),
      p, width, height
  )

  list(p = p, df = df, descr = descr_output)
}

stat_supergroups <- function(stat_group, supergroups) {
  list_super <- lapply(seq_len(length(supergroups)), function(i) {
    join_groups <- subset(stat_group, group %in% supergroups[[i]])
    stat_super <- data.frame(join_groups, supergroup = names(supergroups[i]))
  })

  do.call(rbind, list_super)
}

plot_releffects <- function(manova, genes, text_x = 8, save = FALSE,
  width = 12, height = 7, append = "", col = NULL) {
  obs_manova <- manova[[1]]$twogroupreleffects
  row_manova <- rownames(obs_manova)

  releffects <- lapply(manova[-1], "[[", 2)

  df_perc_boot <- lapply(seq_len(nrow(obs_manova)), function(i) {
    q <- sapply(seq_len(ncol(obs_manova)), function(j) {
      v <- sapply(releffects, function(r) {
        row <- rownames(r)
        if(row[1] != row_manova[1]) r <- r[2:1, ]
        r[i, j]
      })
      quantile(v, probs = c(0.025, 0.975))
    })
    colnames(q) <- colnames(obs_manova)
    df <- data.frame(Var2 = colnames(q), t(q),
      group = ifelse(i == 1, row_manova[1], row_manova[2]))
    df$Var2 <- factor(df$Var2, levels = sort(unique(df$Var2)), ordered = TRUE)
    return(df)
  })
  names(df_perc_boot) <- row_manova

  df <- melt(
    cbind(obs_manova, group = row_manova),
    id.vars = "group"
  )
  df$Var1 <- "point_est"
  df <- df[, c(4, 2, 3, 1)]
  colnames(df)[2] <- "Var2"

  df <- dplyr::arrange(df, group, desc(value), Var2)
  df$Var2 <- as.character(df$Var2)
  df$Var2 <- factor(df$Var2, levels = unique(df$Var2), ordered = TRUE)

  for (k in seq_len(2)) {
    df_perc_boot[[k]]$Var2 <- factor(df_perc_boot[[k]]$Var2, levels = levels(df$Var2), ordered = TRUE)
  }

  if (is.null(col)) {
    if (genes == "neutrophil") {
      col  <- col_neutrophil
    } else if (genes == "cytokine") {
      col <- col_cytokine
    } else if (genes == "other") {
      col <- "pink"
    } else {
      col <- "gold"
    }
  }

  if (names(df_perc_boot)[1] == "Control") {
    col1 <- "#555252"
    col2 <- col
    k <- 2
  } else {
    col1 <- col
    col2 <- "#555252"
    k <- 1
  }

  x_title <- ifelse(tolower(genes) == "neutrophil",
    "Neutrophil mediated immunity genes",
    "Cytokine-mediated signaling pathway genes")

  p1 <- ggplot(data = subset(df, group == "Covid"), aes(x = Var2, y = value,
    group = group, color = group)) +
    geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[k]],
      aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = col,
      alpha = 0.3) +
    geom_line(show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    geom_hline(yintercept = 0.5, alpha = 0.6, linetype = "dashed", size = 1) +
    scale_colour_manual("", values = c(Covid = col,
      Control = "#555252")) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 90, size = text_x),
      axis.line.y = element_blank(),
      axis.title = element_text(size = 15),
      panel.grid.major = element_line()
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(y = "Relative effect", x = x_title)

  p2 <- ggplot(data = df, aes(x = Var2, y = value, group = group,
    color = group)) +
    geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[1]],
      aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = col1,
      alpha = 0.3) +
    geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[2]],
      aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group), fill = col2,
      alpha = 0.3) +
    geom_line() +
    geom_point() +
    scale_colour_manual("", values = c(Covid = col,
      Control = "#555252"),
      labels = c(Covid = "SARS-CoV-2", Control = "negative SARS-CoV-2")) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 12),
      axis.text.x = element_text(angle = 90, size = text_x),
      axis.line.y = element_blank(),
      axis.title = element_text(size = 15),
      panel.grid.major = element_line(),
      legend.position = "top",
      legend.box.background = element_rect()
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(y = "Relative effect", x = x_title)

  save_plot(
      save,
      "MANOVA_Rank-Data", str_to_title(genes),
      "_Relative-Effects_Comparison-Between-Control-And-Covid",
      p2, width, height
  )

  return(list(p1 = p1, p2 = p2))
}

circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100){
    r = diameter
    tt <- seq(0, 2 * pi, length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

plot_cca <- function(cancor, n, save = FALSE, width = 12, height = 7,
  plot_type = 1, k_min = 0.6) {
cca3 <- cancor

cca_x <- cca3$CoefStruct11[, n:(n + 1)]
df_cca_x <- data.frame(cca_x, genes = rownames(cca_x), group = "Neutrophil",
  col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_x[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 1))
  }), check.names = FALSE
)
colnames(df_cca_x)[1:2] <- paste0("CV", 1:2)

cca_y <- cca3$CoefStruct22[, n:(n + 1)]
df_cca_y <- data.frame(cca_y, genes = rownames(cca_y), group = "Cytokine",
  col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_y[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 2))
  })
)
colnames(df_cca_y)[1:2] <- paste0("CV", 1:2)

cca_xy <- cca3$CoefStruct12[, n:(n + 1)]
df_cca_xy <- data.frame(cca_xy, genes = rownames(cca_xy), group = "Cytokine",
  col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_xy[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 2))
  })
)
colnames(df_cca_xy)[1:2] <- paste0("CV", 1:2)

cca_yx <- cca3$CoefStruct21[, n:(n + 1)]
df_cca_yx <- data.frame(cca_yx, genes = rownames(cca_yx), group = "Cytokine",
  col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_yx[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 2))
  })
)
colnames(df_cca_yx)[1:2] <- paste0("CV", 1:2)

if(plot_type == 2) {
  df_cca_x$CV2 <- df_cca_xy$CV1
  df_cca_y$CV2 <- df_cca_yx$CV1
}

df_cca <- rbind(df_cca_x, df_cca_y)

facet_label <- c(
  Neutrophil = "Neutrophil mediated immunity",
  Cytokine = "Cytokine-mediated signaling pathway"
)

df_cca$label_x <- df_cca$label_y <- NA
df_cca$label_x[grep("Neutrophil", df_cca$group)[1]] <- "Ne-CV1"
df_cca$label_y[grep("Neutrophil", df_cca$group)[1]] <- "Ne-CV2"
df_cca$label_x[grep("Cytokine", df_cca$group)[1]] <- "Cy-CV1"
df_cca$label_y[grep("Cytokine", df_cca$group)[1]] <- "Cy-CV2"

p <- ggplot(data = df_cca, aes(x = CV1, y = CV2, colour = col, group = group)) +
 geom_point(alpha = 1, size = 2, show.legend = FALSE) +
 geom_path(data = circleFun(), aes(x, y), inherit.aes = FALSE, alpha = 0.6) +
 geom_hline(yintercept = 0, alpha = 0.6) +
 geom_vline(xintercept = 0, alpha = 0.6) +
 geom_text(aes(label = label_x), x = 0.2, y = 0.05, color = "gold",
  size = 4.2, fontface = "bold") +
 geom_text(aes(label = label_y), x = -0.05, y = -0.2, color = "gold", 
  angle = 90, size = 4.2, fontface = "bold") +
 facet_wrap(~ group, labeller = as_labeller(facet_label)) +
 coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
 scale_colour_manual(values = c("0" = "grey", "1" = col_neutrophil,
  "2" = col_cytokine)) +
 geom_text_repel(aes(label = ifelse(col != "0", genes, NA)), size = 3.3,
  segment.alpha = 0.4, segment.size = 0.3, min.segment.length = 0,
  fontface = "bold", show.legend = FALSE) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 11),
    axis.line = element_blank(),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 13)
  ) +
  labs(y = paste("Correlation with canonical variate", n + 1),
    x = paste("Correlation with canonical variate", n))

  save_plot(
      save,
      "CCA_Log-Data",
      "Correlation-With-Two-Pairs-Of-Canonical-Variates", "",
      p, width, height
  )

  output <- list(p = p, df = df_cca)

  return(output)
}

network_cca <- function(X, Y, cca3, signed_cor, n = 2, k = 0.86, seed = 100,
  label_size = 2.3, save = FALSE) {
  cor_genes_cca <- cor(cbind(X, Y))

  n <- min(n, length(signed_cor))
  n_cor <- seq_len(n)

  cor_cv <- signed_cor[n_cor]

  u_cv <- do.call(cbind, lapply(n_cor, function(i) {
    cor_ux_cv <- cca3$CoefStruct11[, i]
    cor_uy_cv <- cca3$CoefStruct21[, i]
    u_cv <- c(cor_ux_cv, cor_uy_cv)
  }))

  v_cv <- do.call(cbind, lapply(n_cor, function(i) {
    cor_vx_cv <- cca3$CoefStruct12[, i]
    cor_vy_cv <- cca3$CoefStruct22[, i]
    v_cv <- c(cor_vx_cv, cor_vy_cv)
  }))

  diag_cor <- diag(cor_cv)
  diag_1 <- diag(2 * n)
  diag_1[1:n, (n + 1):ncol(diag_1)] <- diag_cor
  diag_1[(n + 1):nrow(diag_1), 1:n] <- diag_cor

  cv_cols <- cbind(u_cv, v_cv)
  cv_rows <- cbind(t(cv_cols), diag_1)
  bind_cols <- cbind(cor_genes_cca, cv_cols)
  bind_rows <- rbind(bind_cols, cv_rows)

  colnames(bind_rows) <- rownames(bind_rows) <- c(
    colnames(cor_genes_cca), paste0("Ne-CV", n_cor), paste0("Cy-CV", n_cor)
  )

  diag(bind_rows) <- 0

  cat("is symmetric? ", is.symmetric.matrix(bind_rows), "\n")

  cor_genes_cv_k <- bind_rows
  cor_genes_cv_k[abs(cor_genes_cv_k) < k] <- 0
  cor_genes_cv_k[abs(cor_genes_cv_k) >= k] <- 1

  type.genes <- ifelse(
    colnames(bind_rows) %in% df_selected_genes[, 1],
    "Neutrophil mediated immunity",
    ifelse(colnames(bind_rows) %in% df_selected_genes[, 2],
      "Cytokine-mediated signaling pathway",
      "Canonical variate"
    )
  )
  name.genes <- colnames(bind_rows)
  name.genes[colSums(cor_genes_cv_k) == 0] <- ""

  net <- network(cor_genes_cv_k, directed = FALSE)
  net %v% "type" <- type.genes
  net %v% "name" <- name.genes

  set.seed(seed)
  df_network <- ggnetwork(net, layout = "kamadakawai")
  p1 <- ggplot(data = df_network, aes(x, y, xend = xend, yend = yend)) +
    geom_edges(alpha = 0.8, color = "grey") +
    geom_nodes(aes(color = type), alpha = 0.5, size = 5,
      show.legend = FALSE) +
    geom_nodetext(aes(label = name), fontface = "bold", size = label_size,
      alpha = 0.7) +
    scale_colour_manual("", values = c(
      "Neutrophil mediated immunity" = col_neutrophil,
      "Cytokine-mediated signaling pathway" = col_cytokine,
      "Canonical variate" = "gold")) +
    theme_blank()

  save_plot(save,
      "CCA_Log-Data",
      "Correlation-Between-Genes-And-Two-Pairs-Of-Canonical-Variates", "",
      p1,
      width = 12, height = 7
  )

  return(p1)
}

A_fun_1 <- function(pos, n_col) {
  matrix(c(rep(0, pos - 1), 1, rep(0, n_col - pos)), nrow = 1)
}

A_fun_2 <- function(pos, n_col, n_coef) {
  cbind(
    matrix(0, n_coef, pos - 1),
    diag(n_coef),
    matrix(0, n_coef, n_col - pos - n_coef + 1)
  )
}

H0_fun_A <- function(A = NULL, B_ind_zero = NULL) {
  H0 <- lm.t.fast(formula, data,
    dDist = "dmvn", A = A,
    Chisq = FALSE, PLOT = FALSE,
    B_ind_zero = B_ind_zero
  )
}

H0_fun_A.boot <- function(A = NULL, B_null, f_u) {
  H0.boot <- lm.t.boot(formula, data,
    dDist = "dmvn", A = A,
    Chisq = FALSE, PLOT = FALSE,
    B_null = B_null, f_u = f_u
  )
}

boot_H0 <- function(line = 0, PLOT = FALSE, fit, fit_null, j) {
  if (line == 0) {
    A.H0 <- NULL
    fit <- H0_fun_A(A = A.H0)
    B_ind_zero <- seq_len(length(fit$Matrices$B))
    fit.0 <- H0_fun_A(A = A.H0, B_ind_zero = B_ind_zero)
  } else if (line > 0){
    A.H0 <- A_fun_1(line, n_covar)
    fit <- H0_fun_A(A = A.H0)
    B_ind_zero <- matrix(c(rep(line, n_var), 1:n_var), ncol = 2)
    fit.0 <- H0_fun_A(A = A.H0, B_ind_zero = B_ind_zero)
  } else {
    A.H0 <- NULL
    fit <- fit
    fit.0 <- fit_null
  }

  fit.boot <- list()
  for (i in 1:n_boot) {
    fit.boot[[i]] <- H0_fun_A.boot(
      A = A.H0,
      B_null = fit.0$to.boot$B_null, f_u = fit.0$to.boo$w2
    )

    cat(
      "fit.boot/", "i = ", i,
      " / est.obs = ", fit.boot[[i]]$Test$est.obs, "\n"
    )
  }

  if(line >= 0) {
    boot.est <- sapply(fit.boot, function(i) i$Test$est.obs)
    p.boot <- sum(boot.est > fit$Test$est.obs) / n_boot
    p.obs <- fit$Test$p.value
  } else {
    boot.est <- sapply(fit.boot, function(i) i$Estimates$Summary[j])
    t.obs <- fit$Estimates$Summary[j]
    upr <- sum(boot.est > t.obs) / n_boot
    lwr <- sum(boot.est <= t.obs) / n_boot
    p.boot <- 2*min(lwr, upr)
    p.obs <- fit$Estimates$t.pvalue[j]
  }

  if (PLOT) {
    dev.new()
    hist(boot.est)
    abline(v = fit$Test$est.obs, col = 2, lwd = 2)
  }

  output <- list(
    model = list(fit = fit, fit.0 = fit.0, fit.boot = fit.boot),
    result = list(boot.est = boot.est, p.boot = p.boot, p.obs = p.obs)
  )

  return(output)
}

boot_B <- function(fit, PLOT = FALSE) {
  B <- fit$Matrices$B
  ind_coef <- seq_len(length(B))

  boot.simu <- lapply(ind_coef, function(j) {
    cat("\nind_coef = ", j, "-----------------------------------------", "\n")
    fit.0 <- H0_fun_A(A = NULL, B_ind_zero = j)
    b <- boot_H0(line = -1, PLOT = PLOT, fit = fit, fit_null = fit.0, j = j)
  })

  p.boot <- matrix(sapply(boot.simu, function(k) k$result$p.boot), nrow(B))
  p.obs <- matrix(sapply(boot.simu, function(k) k$result$p.obs), nrow(B))
  colnames(p.obs) <- colnames(p.boot) <- colnames(B)
  rownames(p.obs) <- rownames(p.boot) <- rownames(B)

  output <- list(B = B, p.boot = p.boot, p.obs = p.obs)
  return(output)
}

prediction_data <- function(v_name, x_title, df, x_pos, y_pos, save = TRUE) {
  x1 <- data.frame(Gender = rep(0:1, each = 100), summarise_all(df, mean))
  x2 <- do.call(data.frame, lapply(x1[, -1], function(i) x1$Gender * i))
  newdata <- data.frame(x1, x2)
  k_df <- colnames(df)
  colnames(newdata) <- c("Gender", k_df, paste0("Gender_", k_df))

  v_range <- range(df_mlm[, v_name])
  newdata[, v_name] <- rep(seq(v_range[1], v_range[2], length.out = 100), 2)
  newdata[, paste0("Gender_", v_name)] <- newdata[, "Gender"] * newdata[, v_name]

  B_boot <- lapply(seq_len(n_boot), function(i) {
    cat("pred.boot / i = ", i, "\n")
    est <- lm.t.boot(formula, df_mlm,
      B_null = fit$to.boot$B_null, f_u = fit$to.boot$w2
    )
    est$Boot.par$B_boot
  })
  pred_boot <- pred.mlm.boot(fit.fast, B_boot, newdata, alpha = 0.05)

  pred <- pred_boot
  df_pred <- data.frame(
    var_x = newdata[, v_name],
    group = factor(newdata[, "Gender"], levels = 0:1, labels = c("Male", "Female")),
    fit = pred[, , 1]#,
    #lwr = pred[, , 2],
    #upr = pred[, , 3]
  )

  melt_pred <- melt(df_pred, id.vars = c("group", "var_x"))
  melt_pred$var_y <- substring(melt_pred$variable, 5)
  melt_pred$p1 <- melt_pred$p2 <- melt_pred$p3 <- NA
  melt_pred$p1[melt_pred$var_y == "CD14"][1] <- format(B_coef$p.boot[v_name, "CD14"], nsmall = 2)
  melt_pred$p2[melt_pred$var_y == "CD14"][1] <- format(B_coef$p.boot[paste0("I(Gender * ", v_name, ")"), "CD14"], nsmall = 2)
  melt_pred$p3[melt_pred$var_y == "CD14"][1] <- format(B_coef$p.boot["Gender", "CD14"], nsmall = 2)
  melt_pred$p1[melt_pred$var_y == "CXCR1"][1] <- format(B_coef$p.boot[v_name, "CXCR1"], nsmall = 2)
  melt_pred$p2[melt_pred$var_y == "CXCR1"][1] <- format(B_coef$p.boot[paste0("I(Gender * ", v_name, ")"), "CXCR1"], nsmall = 2)
  melt_pred$p3[melt_pred$var_y == "CXCR1"][1] <- format(B_coef$p.boot["Gender", "CXCR1"], nsmall = 2)
  melt_pred$p1[melt_pred$var_y == "CXCR2"][1] <- format(B_coef$p.boot[v_name, "CXCR2"], nsmall = 2)
  melt_pred$p2[melt_pred$var_y == "CXCR2"][1] <- format(B_coef$p.boot[paste0("I(Gender * ", v_name, ")"), "CXCR2"], nsmall = 2)
  melt_pred$p3[melt_pred$var_y == "CXCR2"][1] <- format(B_coef$p.boot["Gender", "CXCR2"], nsmall = 2)
  melt_pred$p1[melt_pred$var_y == "IL1B"][1] <- format(B_coef$p.boot[v_name, "IL1B"], nsmall = 2)
  melt_pred$p2[melt_pred$var_y == "IL1B"][1] <- format(B_coef$p.boot[paste0("I(Gender * ", v_name, ")"), "IL1B"], nsmall = 2)
  melt_pred$p3[melt_pred$var_y == "IL1B"][1] <- format(B_coef$p.boot["Gender", "IL1B"], nsmall = 2)
  melt_pred$p1[melt_pred$var_y == "NBEAL2"][1] <- format(B_coef$p.boot[v_name, "NBEAL2"], nsmall = 2)
  melt_pred$p2[melt_pred$var_y == "NBEAL2"][1] <- format(B_coef$p.boot[paste0("I(Gender * ", v_name, ")"), "NBEAL2"], nsmall = 2)
  melt_pred$p3[melt_pred$var_y == "NBEAL2"][1] <- format(B_coef$p.boot["Gender", "NBEAL2"], nsmall = 2)
  melt_pred$p1[melt_pred$var_y == "S100A9"][1] <- format(B_coef$p.boot[v_name, "S100A9"], nsmall = 2)
  melt_pred$p2[melt_pred$var_y == "S100A9"][1] <- format(B_coef$p.boot[paste0("I(Gender * ", v_name, ")"), "S100A9"], nsmall = 2)
  melt_pred$p3[melt_pred$var_y == "S100A9"][1] <- format(B_coef$p.boot["Gender", "S100A9"], nsmall = 2)

  df_p <- melt(df_mlm[, c(v_name, "Gender", target_genes)],
    id.vars = c(v_name, "Gender"))
  colnames(df_p)[3] <- "var_y"
  colnames(df_p)[colnames(df_p) == v_name] <- "v_point"
  df_p[, "Gender"] <- factor(df_p[, "Gender"], levels = 0:1, labels = c("Male", "Female"))

  p1 <- (ggplot(data = melt_pred,
    aes(x = var_x, y = value, colour = group, group = group)) +
    geom_point(
      inherit.aes = FALSE, data = df_p, aes(x = v_point, y = value,
        colour = factor(Gender)), alpha = 0.3, size = 1) +
    scale_colour_manual("", values = c(Male = "#24e43e",
      Female = "#b520f0"), labels = c(Male = "Male - SC2", Female = "Female - SC2")) +
    geom_line(size = 1) +
    geom_text(colour = "black", size = ifelse(eval(v_name) == "n1_ct", 2.8, 3), x = x_pos, y = y_pos,
      vjust = "inward", hjust = "inward",
      aes(label = ifelse(!is.na(p1), paste0(
        ifelse(eval(v_name) == "n1_ct", "Viral load", eval(v_name)), ": p-value = ", p1,
        "\nGender: p-value = ", p3,
        "\nGender x ", ifelse(eval(v_name) == "n1_ct", "Viral load", eval(v_name)), ": p-value = ", p2),
        NA))) +
    facet_wrap(~ var_y, nrow = 1) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 10),
      axis.line.y = element_blank(),
      axis.title = element_text(size = 13),
      panel.grid.major = element_line(),
      legend.position = "top",
      legend.box.background = element_rect()
    ) +
    labs(y = "Transformed gene expression", x = x_title)
  )

  save_plot(save, "MLM_LogData", "Covariate-", v_name, p1, width = 12, height = 5)

  return(p1)
}

new_CANCOR <- function (data, set1, set2, plot = TRUE, plotCV = 1, plotcoefs = "structure", 
    verbose = TRUE) 
{
    data <- as.data.frame(data[, c(set1, set2)])
    if (anyNA(data) == TRUE) {
        data <- na.omit(data)
        NAflag = TRUE
    }
    else {
        NAflag = FALSE
    }
    set1data <- as.data.frame(data[, set1])
    set2data <- as.data.frame(data[, set2])
    Ncases <- nrow(set1data)
    NVset1 <- ncol(set1data)
    NVset2 <- ncol(set2data)
    CorrelSet1 <- stats::cor(set1data)
    CorrelSet2 <- stats::cor(set2data)
    CorrelSet1n2 <- stats::cor(set2data, set1data)
    output <- DFA.CANCOR:::canonical.cor(set1data, set2data)
    cancorrels <- output$cancorrels
    mv_Wilk <- DFA.CANCOR:::Wilk(rho = cancorrels[, 2], Ncases = Ncases, p = NVset1, 
        q = NVset2)
    colnames(mv_Wilk) <- c("            Wilk's Lambda", "   F-approx. ", 
        "  df1", "    df2", "         p")
    rownames(mv_Wilk) <- paste(1:nrow(mv_Wilk), paste("through ", 
        nrow(mv_Wilk), sep = ""))
    mv_Pillai <- DFA.CANCOR:::Pillai(rho = cancorrels[, 2], Ncases = Ncases, 
        p = NVset1, q = NVset2)
    colnames(mv_Pillai) <- c("  Pillai-Bartlett Trace", "   F-approx. ", 
        "  df1", "    df2", "         p")
    rownames(mv_Pillai) <- paste(1:nrow(mv_Pillai), paste("through ", 
        nrow(mv_Pillai), sep = ""))
    mv_Hotelling <- DFA.CANCOR:::Hotelling(rho = cancorrels[, 2], Ncases = Ncases, 
        p = NVset1, q = NVset2)
    colnames(mv_Hotelling) <- c("   Hotelling-Lawley Trace", 
        "   F-approx. ", "  df1", "    df2", "         p")
    rownames(mv_Hotelling) <- paste(1:nrow(mv_Hotelling), paste("through ", 
        nrow(mv_Hotelling), sep = ""))
    mv_Roy <- DFA.CANCOR:::RoyRoot(rho = cancorrels[, 2], Ncases = Ncases, 
        p = NVset1, q = NVset2)
    colnames(mv_Roy) <- c("      Roy's Largest Root", "   F-approx. ", 
        "  df1", "    df2", "         p")
    rownames(mv_Roy) <- paste(1:nrow(mv_Roy), paste("through ", 
        nrow(mv_Roy), sep = ""))
    mv_BartlettV <- DFA.CANCOR:::BartlettV(rho = cancorrels[, 2], Ncases = Ncases, 
        p = NVset1, q = NVset2)
    colnames(mv_BartlettV) <- c("            Wilk's Lambda", 
        "   F-approx. ", "  df", "         p")
    rownames(mv_BartlettV) <- paste(1:nrow(mv_BartlettV), paste("through ", 
        nrow(mv_BartlettV), sep = ""))
    mv_Rao <- DFA.CANCOR:::Rao(rho = cancorrels[, 2], Ncases = Ncases, p = NVset1, 
        q = NVset2)
    colnames(mv_Rao) <- c("   Wilk's Lambda", "   F-approx. ", 
        "    df1", "       df2", "         p")
    rownames(mv_Rao) <- paste(1:nrow(mv_Rao), paste("through ", 
        nrow(mv_Rao), sep = ""))
    colnames(output$cancorrels) <- c("  Eigenvalue", "  Canonical r", 
        "   Canonical r sq.", "         t", "     df", "   p value")
    rownames(output$cancorrels) <- paste(" Canonical function ", 
        1:nrow(output$cancorrels), sep = "")
    colnames(output$raw1) <- paste("     CV", 1:ncol(output$raw1), 
        sep = "")
    colnames(output$raw2) <- paste("     CV", 1:ncol(output$raw2), 
        sep = "")
    colnames(output$struct11) <- paste("     CV", 1:ncol(output$struct11), 
        sep = "")
    colnames(output$struct21) <- paste("     CV", 1:ncol(output$struct21), 
        sep = "")
    colnames(output$struct12) <- paste("     CV", 1:ncol(output$struct12), 
        sep = "")
    colnames(output$struct22) <- paste("     CV", 1:ncol(output$struct22), 
        sep = "")
    colnames(output$stand1) <- paste("     CV", 1:ncol(output$stand1), 
        sep = "")
    colnames(output$stand2) <- paste("     CV", 1:ncol(output$stand2), 
        sep = "")
    if (plot == TRUE | is.null(plot)) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        if (is.null(plot)) 
            plotCV = 1
        if (plotcoefs == "structure" | is.null(plotcoefs)) {
            layout(matrix(c(1, 2, 1, 2), nrow = 2, byrow = T))
            barplot(output$struct11[, plotCV], ylim = c(-1, 1), 
                ylab = "CV1", main = "Set 1 Structure Coefficients", 
                col = "blue", las = 2)
            box()
            barplot(output$struct22[, plotCV], ylim = c(-1, 1), 
                ylab = "CV1", main = "Set 2 Structure Coefficients", 
                col = "blue", las = 2)
            box()
        }
        if (plotcoefs == "standardized") {
            layout(matrix(c(1, 2, 1, 2), nrow = 2, byrow = T))
            barplot(output$stand1[, plotCV], ylim = c(-1.3, 1.3), 
                ylab = "CV1", main = "Set 1 Standardized Coefficients", 
                col = "blue", las = 2)
            box()
            barplot(output$stand2[, plotCV], ylim = c(-1.3, 1.3), 
                ylab = "CV1", main = "Set 2 Standardized Coefficients", 
                col = "blue", las = 2)
            box()
        }
    }
    if (verbose == TRUE) {
        cat("\n\n\nCanonical Correlation Analysis\n")
        if (NAflag == TRUE) 
            cat("\n\nCases with missing values were found and removed from the data matrix.\n\n")
        cat("\n\nPearson correlations for Set 1:\n\n")
        print(round(CorrelSet1, 2))
        cat("\n\nPearson correlations for Set 2:\n\n")
        print(round(CorrelSet2, 2))
        cat("\n\nPearson correlations between Set 1 & Set 2:\n\n")
        print(round(CorrelSet1n2, 2))
        cat("\n\n\nMultivariate peel-down significance tests:\n\n")
        print(round(mv_Wilk, 4))
        cat("\n\n")
        print(round(mv_Pillai, 4))
        cat("\n\n")
        print(round(mv_Hotelling, 4))
        cat("\n\n")
        print(round(mv_Roy, 4))
        cat("\n\n")
        cat("\n\n\nBartlett's V test:\n")
        print(round(mv_BartlettV, 4))
        cat("\n\n")
        cat("\n\n\nRao's V test:\n")
        print(round(mv_Rao, 4))
        cat("\n\n")
        cat("\n\n\nCanonical correlations:\n\n")
        print(round(output$cancorrels, 3))
        cat("\nThe above t tests are for the Pearson correlations between the canonical variate scores")
        cat("\nfor each function, i.e., they are not the multivariate significance tests.\n\n")
        cat("\n\n\n\nRaw canonical coefficients for Set 1:\n\n")
        print(round(output$raw1, 2))
        cat("\n\nRaw canonical coefficients for Set 2:\n\n")
        print(round(output$raw2, 2))
        cat("\n\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 1 variates:\n\n")
        print(round(output$struct11, 2))
        cat("\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 1 variates:\n\n")
        print(round(output$struct21, 2))
        cat("\n\nStructure coefficients (Pearson correlations) for Set 1 variables with the Set 2 variates:\n\n")
        print(round(output$struct12, 2))
        cat("\n\nStructure coefficients (Pearson correlations) for Set 2 variables with the Set 2 variates:\n\n")
        print(round(output$struct22, 2))
        cat("\n\n\nStandardized coefficients for Set 1 variables:\n\n")
        print(round(output$stand1, 2))
        cat("\n\nStandardized coefficients for Set 2 variables:\n\n")
        print(round(output$stand2, 2))
        cat("\n\n\n\n")
    }
    CANCORoutput <- list(cancorrels = cancorrels, CoefRawSet1 = output$raw1, 
        CoefRawSet2 = output$raw2, CoefStruct11 = output$struct11, 
        CoefStruct21 = output$struct21, CoefStruct12 = output$struct12, 
        CoefStruct22 = output$struct22, CoefStandSet1 = output$stand1, 
        CoefStandSet2 = output$stand2, mv_Wilk = mv_Wilk, mv_Pillai = mv_Pillai, 
        mv_Hotelling = mv_Hotelling, mv_Roy = mv_Roy, mv_BartlettV = mv_BartlettV, 
        mv_Rao = mv_Rao, CorrelSet1 = CorrelSet1, CorrelSet2 = CorrelSet2, 
        CorrelSet1n2 = CorrelSet1n2)
    return(CANCORoutput)
}

lm.t <- function(formula, data, dDist = "dmvn", v.par = NULL, k = 2, 
                 A = NULL, D = NULL, C = NULL, alpha = 0.05, 
                 PLOT = TRUE, Chisq = FALSE, Intercept = TRUE, 
                 B_ind_zero = NULL){
 ## Warnings:
 if((is.null(v.par) | length(v.par) != 1 | !is.numeric(v.par)) & dDist != "dmvn"){
  dDist <- "dmvn"
  cat("dDist automaticaly changed to dmvn.\n")
 }

 ## Iterations:
 n <- nrow(data)
 Form <- Formula(formula)

 R <- lapply(c(0, 1:n), function(i){
   if(i) data <- data[-i, ]
   dt <- model.frame(Form, data = data)

   Y <- as.matrix(model.part(Form, data = dt, lhs = 1))
   if(any(grep("\\.", colnames(Y)))) colnames(Y) <- sapply(strsplit(colnames(Y), "\\."), "[[", 2)
   X <- as.matrix(data.frame(Intercept = 1, model.part(Form, data = dt, rhs = 1), check.names = FALSE))
   if(!Intercept) X <- X[, -which(colnames(X) == "Intercept")]
   n <- nrow(Y)
   p <- ncol(Y)

   B <- solve(t(X)%*%X)%*%t(X)%*%Y
    if(!is.null(B_ind_zero)) B[B_ind_zero] <- 0
   QB <- t(Y - X%*%B)%*%(Y - X%*%B)

    if(dDist == "dmvt"){
     name.par <- "df"
     u0 <- 1/n
     scl <- v.par/(v.par - 2)
     Sig <- u0*QB
     Cov <- Sig
    } else if(dDist == "dmvpe"){
       name.par <- "kappa"
       u0 <- ((v.par*p^(v.par - 1))/n)^(1/v.par)
       scl <- (2^(1/v.par)*gamma((p + 2)/(2*v.par)))/(p*gamma(p/(2*v.par)))
       Sig <- u0*QB
       Cov <- u0*scl*QB
      } else if(dDist == "dmvn"){
         name.par <- "empty"
         u0 <- 1/n #(n - nrow(B))
         scl <- 1
         Sig <- u0*QB
         Cov <- Sig
        } else{
           stop("Choose a suitable distribution!")
          }

   list.args <- list(x = Y, mu = X%*%B, S = Cov, name.par = v.par, log = TRUE)
   if(!(name.par == "empty")) names(list.args)[4] <- name.par
   if(name.par == "empty") list.args <- list.args[-4]
   logL <- sum(do.call(dDist, list.args))

   r <- list(Y = Y, X = X, B = B, QB = QB, u0 = u0, Sig = Sig, Cov = Cov,
             scl = scl, logL = logL, n = n, p = p, name.par = name.par)
  }
 )

 ## Estimates:
 for(i in 1:length(R[[1]])) assign(names(R[[1]])[i], R[[1]][[i]])
 P <- X%*%solve(t(X)%*%X)%*%t(X)
 E <- (diag(n) - P)%*%Y
 Yfit <- X%*%B

 ## GAIC:
 n.pars <- prod(dim(B)) #+ prod(dim(Sig))
 gaic <- k*n.pars - 2*logL

 ## Hypothesis Testing:
 q <- nrow(B)
 if(is.null(A)) A <- diag(q)
 m <- nrow(A)
 if(is.null(D)) D <- diag(p)
 u <- ncol(D)
 if(is.null(C)) C <- matrix(0, nrow = m, ncol = u)

 ### General Linear Hypothesis:
 T <- t(A%*%B%*%D - C)%*%solve(A%*%solve(t(X)%*%X)%*%t(A))%*%(A%*%B%*%D - C)
 lambda <- det(t(D)%*%QB%*%D)/det(t(D)%*%QB%*%D + T)
 kk <- min(dim(C))
  if(p == 1 & !Chisq){
   df1 <- kk
   df2 <- n
   df <- c(df1 = df1, df2 = df2)
   name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
   est.obs <- ((1 - lambda)/lambda)*df2/df1
   est.q <- qf(1 - alpha, df1, df2)
   est.pvalue <- pf(est.obs, df1, df2, lower.tail = FALSE)
  } else if(kk == 1 & !Chisq){
     df1 <- p
     df2 <- n + 1 - p
     df <- c(df1 = df1, df2 = df2)
     name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
     est.obs <- ((1 - lambda)/lambda)*df2/df1
     est.q <- qf(1 - alpha, df1, df2)
     est.pvalue <- pf(est.obs, df1, df2, lower.tail = FALSE)
    } else if(p == 2 & !Chisq){
       df1 <- kk
       df2 <- n - 1
       df <- c(df1 = df1, df2 = df2)
       name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
       est.obs <- ((1 - sqrt(lambda))/sqrt(lambda))*df2/df1
       est.q <- qf(1 - alpha, 2*df1, 2*df2)
       est.pvalue <- pf(est.obs, 2*df1, 2*df2, lower.tail = FALSE)
      } else if(kk == 2 & !Chisq){
         df1 <- p
         df2 <- n + 1 - p
         df <- c(df1 = df1, df2 = df2)
         name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
         est.obs <- ((1 - sqrt(lambda))/sqrt(lambda))*df2/df1
         est.q <- qf(1 - alpha, 2*df1, 2*df2)
         est.pvalue <- pf(est.obs, 2*df1, 2*df2, lower.tail = FALSE)
        } else{
           df1 <- p*kk
           df2 <- NULL
           df <- c(df1 = df1, df2 = df2)
           name.test <- paste("Chisq", df1, paste("(", alpha, ")", collapse = ""), sep = "_")
           est.obs <- -(n - 0.5*(p - kk + 1))*log(lambda)
           est.q <- qchisq(1 - alpha, df1)
           est.pvalue <- pchisq(est.obs, df1, lower.tail = FALSE)
          }

 ### Multivariate Student-t Test:
 if(dDist == "dmvt"){
  Bt.Ct <- matrix(c(A%*%B%*%D - C), ncol = 1)
  Dt <- scl*t(Bt.Ct)%*%solve((t(D)%*%Sig%*%D)%x%(A%*%solve(t(X)%*%X)%*%t(A)))%*%(Bt.Ct)
  kp <- nrow(Bt.Ct)
  est.t.obs <- Dt/kp
  est.t.q <- qf(1 - alpha, kp, v.par)
  est.t.pvalue <- pf(est.t.obs, kp, v.par, lower.tail = FALSE)
  name.t.test <- paste("F", paste(c(kp, v.par), collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
 } else{
    est.t.obs <- "empty"
    est.t.q <- "empty"
    est.t.pvalue <- "empty"
    name.t.test <- "empty"
   }

 ### T-test:
 sdB <- matrix(sqrt(diag(scl*Sig%x%solve(t(X)%*%X))), nrow = q, byrow = FALSE)
 dimnames(sdB) <- dimnames(B)
 t.obs <- B/sdB
 t.pvalue <- 2*pt(abs(t.obs), n - q, lower.tail = FALSE)

 ## Diagnostics:
 pii <- diag(P)
 ei <- lapply(1:n, function(i) matrix(E[i, ], ncol = 1))
 sig2i <- do.call(rbind, lapply(2:(n+1), function(i) diag(R[[i]]$scl*R[[i]]$Sig)))
 w2 <- E/sqrt(1 - replicate(ncol(E), pii))

 ### Residuals:
 Sinv <- solve(QB/(n - q))
 ri2 <- sapply(1:n, function(i) 1/(1 - pii[i])*t(ei[[i]])%*%Sinv%*%ei[[i]])
 bi <- ri2/(n - q)
 cut.bi <- qbeta(1 - alpha, p/2, (n - p - q)/2)

 Sinv.i <- lapply(2:(n + 1), function(i) solve(R[[i]]$Sig/(R[[i]]$u0*(n - q - 1))))
 ti2 <- sapply(1:n, function(i) 1/(1 - pii[i])*t(ei[[i]])%*%Sinv.i[[i]]%*%ei[[i]])
 fi <- ((n - p - q)*ti2)/(p*(n - q - 1))
 cut.fi <- qf(1 - alpha, p, n - p - q)

 stud.res <- sapply(1:ncol(Y), function(i) E[, i]/sqrt(sig2i[, i]*(1 - pii)))
 colnames(stud.res) <- colnames(Y)

 ### Cook's Distance:
 Di <- (1/q)*(pii/(1-pii))*ri2
 Di.star <- (1/q)*sapply(1:n, function(i) sum(diag(Sinv.i[[i]]%*%t(B - R[[i + 1]]$B)%*%t(X)%*%X%*%(B - R[[i + 1]]$B))))

 ### Scale Ratio:
 SCRi <- ((R[[2]]$u0/u0)^p*(1 - bi))^q*(1/(1 - pii)^p)

 ### Likelihood Displacement:
 LDi <- 2*(logL - sapply(R[-1], function(i) i$logL))

 ### Mahalanobis Distance:
 MD <- sapply(1:n, function(i) matrix(Y[i, ] - (X%*%B)[i, ], nrow = 1)%*%solve(Sig)%*%matrix(Y[i, ] - (X%*%B)[i, ], ncol = 1))
 if(dDist == "dmvn"){
  cut.MD <- qchisq(1 - alpha, p)
 } else if(dDist == "dmvt"){
    MD <- MD/p
    cut.MD <- qf(1 - alpha, p, v.par)
   } else if(dDist == "dmvpe"){
      MD <- MD^v.par
      cut.MD <- qgamma(1 - alpha, shape = p/(2*v.par), rate = 0.5)
     } else{
        stop("Choose a suitable distribution!")
       }

 ### Plots:
 if(PLOT){
  sub <- paste(dDist, ifelse(name.par == "empty", "NULL", v.par), sep = ", ")
 
  dev.new()
  par(mfrow = c(1, 3))

  plot(bi, type = "h", sub = paste("(a)", sub, sep = ": "), main = expression(frac(r[i]^2,n - q)%~%beta))
  abline(h = cut.bi, lty = 2)
  if(any(bi > cut.bi)) text(y = bi[bi > cut.bi], x = (1:n)[bi > cut.bi], labels = which(bi > cut.bi))

  plot(fi, type = "h", sub = paste("(b)", sub, sep = ": "), main = expression(frac(n - p - q, p(n - q - 1))~t[i]^2%~%F))
  abline(h = cut.fi, lty = 2)
  if(any(fi > cut.fi)) text(y = fi[fi > cut.fi], x = (1:n)[fi > cut.fi], labels = which(fi > cut.fi))
  
  plot(MD, type = "h", sub = paste("(c)", sub, sep = ": "), main = "Mahalanobis Distance")
  abline(h = cut.MD, lty = 2)
  if(any(MD > cut.MD)) text(y = MD[MD > cut.MD], x = (1:n)[MD > cut.MD], labels = which(MD > cut.MD))

  dev.new()
  par(mfrow = c(2, 2))

  plot(Di, type = "h", sub = paste("(d)", sub, sep = ": "), main = "Cook's Distance")
  plot(Di.star, type = "h", sub = paste("(e)", sub, sep = ": "), main = "Modified Cook's Distance")

  plot(SCRi, type = "h", sub = paste("(f)", sub, sep = ": "), main = "Scale Ratio")
  plot(LDi, type = "h", sub = paste("(g)", sub, sep = ": "), main = "Likelihood Displacement")

  dev.new()
  par(mfrow = c(1, 2))
  plot(unlist(stud.res) ~ unlist(Yfit), type = "p", pch = 20, xlab = "Fitted Values", ylab = "Studentized Residuals", sub = paste("(h)", sub, sep = ": "), ylim = c(min(stud.res, -2, -max(stud.res)), max(stud.res, 2, -min(stud.res))))
  abline(h = c(-2, 0, 2), lty = 2)
  plot(NA, type = "n", pch = 20, xlab = "Index", ylab = "Studentized Residuals", sub = paste("(i)", sub, sep = ": "), xlim = c(1, nrow(stud.res)), ylim = c(min(stud.res, -2, -max(stud.res)), max(stud.res, 2, -min(stud.res))))
  sapply(1:nrow(stud.res), function(i){
    a <- stud.res[i, ]
    points(a ~ rep(i, ncol(stud.res)), pch = 20)
   }
  )
  abline(h = c(-2, 0, 2), lty = 2)
 }

 ## Outputs:
 L.r <- lapply(1:p, function(i){
   a <- cbind(B[, i], sdB[, i], t.obs[, i], t.pvalue[, i])
   colnames(a) <- c("Estimate", "sd", "t", "pvalue")
   a
  }
 )
 names(L.r) <- colnames(B)

 r <- list(Estimates = list(B = B, Summary = L.r, Cov = scl*Sig, Cor = cov2cor(scl*Sig)), 
           GAIC = list(GAIC = gaic, k = k, logL = logL, n.pars = n.pars), 
           Test = list("U.p.kk.n (kk = min(m, u))" = c(p = p, kk = kk, n = n), lambda= lambda, est.obs = est.obs, est.q = est.q, p.value = est.pvalue, name.test = name.test),
           StudentTest = list("df = m*u"= c(m = m, u = u), est.obs = est.t.obs, est.q = est.t.q, p.value = est.t.pvalue, name.test = name.t.test),
           Matrices = list(A = A, D = D, C = C, B = B, AB = A%*%B, BD = B%*%D, ABD = A%*%B%*%D),
           Diagnostics = list(XB = Yfit, E = E, pii = unname(pii), sig2i = sig2i, stud.res = stud.res, ri2 = ri2, ti2 = ti2, bi = bi, cut.bi = cut.bi, fi = fi, cut.fi = cut.fi, MD = MD, cut.MD = cut.MD, Di = unname(Di), Di.star = Di.star, SCRi = unname(SCRi), LDi = LDi),
           v.par = list(par = ifelse(name.par == "empty", NA, v.par), dDist = dDist),
           to.boot = list(B_null = B, w2 = w2)
 )

 return(r)
}

lm.t.boot <- function(formula, data, dDist = "dmvn", v.par = NULL, k = 2,
                      A = NULL, D = NULL, C = NULL, alpha = 0.05,
                      PLOT = TRUE, Chisq = FALSE, Intercept = TRUE,
                      B_null, f_u) {
  ## Warnings:
  if ((is.null(v.par) | length(v.par) != 1 | !is.numeric(v.par)) & dDist != "dmvn") {
    dDist <- "dmvn"
    cat("dDist automaticaly changed to dmvn.\n")
  }

  ## Iterations:
  n <- nrow(data)
  Form <- Formula(formula)

  dt <- model.frame(Form, data = data)

  X <- as.matrix(data.frame(Intercept = 1, model.part(Form, data = dt, rhs = 1), check.names = FALSE))
  if (!Intercept) X <- X[, -which(colnames(X) == "Intercept")]

  v_boot <- matrix(ifelse(runif(length(f_u)) > 0.5, -1, 1), ncol = ncol(f_u))
  Y <- X %*% B_null + f_u * v_boot

  n <- nrow(Y)
  p <- ncol(Y)

  B <- solve(t(X) %*% X) %*% t(X) %*% Y
  QB <- t(Y - X %*% B) %*% (Y - X %*% B)

  name.par <- "empty"
  u0 <- 1 / n # (n - nrow(B))
  scl <- 1
  Sig <- u0 * QB
  Cov <- Sig

  list.args <- list(x = Y, mu = X %*% B, S = Cov, name.par = v.par, log = TRUE)
  if (!(name.par == "empty")) names(list.args)[4] <- name.par
  if (name.par == "empty") list.args <- list.args[-4]
  logL <- 1 # sum(do.call(dDist, list.args))

  ## Estimates:
  P <- X %*% solve(t(X) %*% X) %*% t(X)
  E <- (diag(n) - P) %*% Y
  Yfit <- X %*% B
  pii <- diag(P)

  ## GAIC:
  n.pars <- prod(dim(B)) #+ prod(dim(Sig))
  gaic <- k * n.pars - 2 * logL

  ## Hypothesis Testing:
  q <- nrow(B)
  if (is.null(A)) A <- diag(q)
  m <- nrow(A)
  if (is.null(D)) D <- diag(p)
  u <- ncol(D)
  if (is.null(C)) C <- matrix(0, nrow = m, ncol = u)

  ### General Linear Hypothesis:
  T <- t(A %*% B %*% D - C) %*% solve(A %*% solve(t(X) %*% X) %*% t(A)) %*% (A %*% B %*% D - C)
  lambda <- det(t(D) %*% QB %*% D) / det(t(D) %*% QB %*% D + T)
  kk <- min(dim(C))
  if (p == 1 & !Chisq) {
    df1 <- kk
    df2 <- n
    df <- c(df1 = df1, df2 = df2)
    name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
    est.obs <- ((1 - lambda) / lambda) * df2 / df1
  } else if (kk == 1 & !Chisq) {
    df1 <- p
    df2 <- n + 1 - p
    df <- c(df1 = df1, df2 = df2)
    name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
    est.obs <- ((1 - lambda) / lambda) * df2 / df1
  } else if (p == 2 & !Chisq) {
    df1 <- kk
    df2 <- n - 1
    df <- c(df1 = df1, df2 = df2)
    name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
    est.obs <- ((1 - sqrt(lambda)) / sqrt(lambda)) * df2 / df1
  } else if (kk == 2 & !Chisq) {
    df1 <- p
    df2 <- n + 1 - p
    df <- c(df1 = df1, df2 = df2)
    name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
    est.obs <- ((1 - sqrt(lambda)) / sqrt(lambda)) * df2 / df1
  } else {
    df1 <- p * kk
    df2 <- NULL
    df <- c(df1 = df1, df2 = df2)
    name.test <- paste("Chisq", df1, paste("(", alpha, ")", collapse = ""), sep = "_")
    est.obs <- -(n - 0.5 * (p - kk + 1)) * log(lambda)
  }

  ### T-test:
  sdB <- matrix(sqrt(diag(scl * Sig %x% solve(t(X) %*% X))), nrow = q, byrow = FALSE)
  dimnames(sdB) <- dimnames(B)
  t.obs <- B / sdB

  r <- list(
    Estimates = list(Summary = t.obs),
    Test = list("U.p.kk.n (kk = min(m, u))" = c(p = p, kk = kk, n = n), lambda = lambda, est.obs = est.obs, name.test = name.test),
    Diagnostics = list(E = E, pii = unname(pii)),
    Boot.par = list(Y_boot = Y, v_boot = v_boot, B_boot = B)
  )

  return(r)
}

lm.t.fast <- function(formula, data, dDist = "dmvn", v.par = NULL, k = 2, 
                 A = NULL, D = NULL, C = NULL, alpha = 0.05, 
                 PLOT = TRUE, Chisq = FALSE, Intercept = TRUE,
                 B_ind_zero = NULL){
 ## Warnings:
 if((is.null(v.par) | length(v.par) != 1 | !is.numeric(v.par)) & dDist != "dmvn"){
  dDist <- "dmvn"
  cat("dDist automaticaly changed to dmvn.\n")
 }

 ## Iterations:
 n <- nrow(data)
 Form <- Formula(formula)

 dt <- model.frame(Form, data = data)

 Y <- as.matrix(model.part(Form, data = dt, lhs = 1))
 if (any(grep("\\.", colnames(Y)))) colnames(Y) <- sapply(strsplit(colnames(Y), "\\."), "[[", 2)
 X <- as.matrix(data.frame(Intercept = 1, model.part(Form, data = dt, rhs = 1), check.names = FALSE))
 if (!Intercept) X <- X[, -which(colnames(X) == "Intercept")]
 n <- nrow(Y)
 p <- ncol(Y)

 B <- solve(t(X) %*% X) %*% t(X) %*% Y
 if (!is.null(B_ind_zero)) B[B_ind_zero] <- 0
 QB <- t(Y - X %*% B) %*% (Y - X %*% B)

 name.par <- "empty"
 u0 <- 1 / n # (n - nrow(B))
 scl <- 1
 Sig <- u0 * QB
 Cov <- Sig

 list.args <- list(x = Y, mu = X %*% B, S = Cov, name.par = v.par, log = TRUE)
 if (!(name.par == "empty")) names(list.args)[4] <- name.par
 if (name.par == "empty") list.args <- list.args[-4]
 logL <- 1 # sum(do.call(dDist, list.args))

 ## Estimates:
 P <- X%*%solve(t(X)%*%X)%*%t(X)
 E <- (diag(n) - P)%*%Y
 Yfit <- X%*%B
 pii <- diag(P)
 w2 <- E/sqrt(1 - replicate(ncol(E), pii))

 ## GAIC:
 n.pars <- prod(dim(B)) #+ prod(dim(Sig))
 gaic <- k*n.pars - 2*logL

 ## Hypothesis Testing:
 q <- nrow(B)
 if(is.null(A)) A <- diag(q)
 m <- nrow(A)
 if(is.null(D)) D <- diag(p)
 u <- ncol(D)
 if(is.null(C)) C <- matrix(0, nrow = m, ncol = u)

 ### General Linear Hypothesis:
 T <- t(A%*%B%*%D - C)%*%solve(A%*%solve(t(X)%*%X)%*%t(A))%*%(A%*%B%*%D - C)
 lambda <- det(t(D)%*%QB%*%D)/det(t(D)%*%QB%*%D + T)
 kk <- min(dim(C))
  if(p == 1 & !Chisq){
   df1 <- kk
   df2 <- n
   df <- c(df1 = df1, df2 = df2)
   name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
   est.obs <- ((1 - lambda)/lambda)*df2/df1
   est.q <- qf(1 - alpha, df1, df2)
   est.pvalue <- pf(est.obs, df1, df2, lower.tail = FALSE)
  } else if(kk == 1 & !Chisq){
     df1 <- p
     df2 <- n + 1 - p
     df <- c(df1 = df1, df2 = df2)
     name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
     est.obs <- ((1 - lambda)/lambda)*df2/df1
     est.q <- qf(1 - alpha, df1, df2)
     est.pvalue <- pf(est.obs, df1, df2, lower.tail = FALSE)
    } else if(p == 2 & !Chisq){
       df1 <- kk
       df2 <- n - 1
       df <- c(df1 = df1, df2 = df2)
       name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
       est.obs <- ((1 - sqrt(lambda))/sqrt(lambda))*df2/df1
       est.q <- qf(1 - alpha, 2*df1, 2*df2)
       est.pvalue <- pf(est.obs, 2*df1, 2*df2, lower.tail = FALSE)
      } else if(kk == 2 & !Chisq){
         df1 <- p
         df2 <- n + 1 - p
         df <- c(df1 = df1, df2 = df2)
         name.test <- paste("F", paste(df, collapse = ","), paste("(", alpha, ")", collapse = ""), sep = "_")
         est.obs <- ((1 - sqrt(lambda))/sqrt(lambda))*df2/df1
         est.q <- qf(1 - alpha, 2*df1, 2*df2)
         est.pvalue <- pf(est.obs, 2*df1, 2*df2, lower.tail = FALSE)
        } else{
           df1 <- p*kk
           df2 <- NULL
           df <- c(df1 = df1, df2 = df2)
           name.test <- paste("Chisq", df1, paste("(", alpha, ")", collapse = ""), sep = "_")
           est.obs <- -(n - 0.5*(p - kk + 1))*log(lambda)
           est.q <- qchisq(1 - alpha, df1)
           est.pvalue <- pchisq(est.obs, df1, lower.tail = FALSE)
          }

 ### T-test:
 sdB <- matrix(sqrt(diag(scl*Sig%x%solve(t(X)%*%X))), nrow = q, byrow = FALSE)
 dimnames(sdB) <- dimnames(B)
 t.obs <- B/sdB
 t.pvalue <- 2*pt(abs(t.obs), n - q, lower.tail = FALSE)

 r <- list(Estimates = list(Summary = t.obs,  t.pvalue =  t.pvalue),
           Test = list("U.p.kk.n (kk = min(m, u))" = c(p = p, kk = kk, n = n), lambda= lambda, est.obs = est.obs, est.q = est.q, p.value = est.pvalue, name.test = name.test),
           Matrices = list(A = A, D = D, C = C, B = B, AB = A%*%B, BD = B%*%D, ABD = A%*%B%*%D),
           Diagnostics = list(E = E, pii = unname(pii)),
           to.boot = list(B_null = B, w2 = w2)
 )

 return(r)
}

pred.mlm.boot <- function(fit, B_boot, newdata, alpha = 0.05) {
  newdata <- as.matrix(cbind(Intercept = 1, newdata))
  n <- nrow(newdata)
  k <- ncol(B_boot[[1]])
  fit <- newdata %*% fit$Matrices$B
  pred_list <- lapply(B_boot, function(i) newdata %*% i)
  limits <- sapply(seq_len(n * k), function(i) {
    v <- unlist(Map(function(x) x[i], pred_list))
    q <- quantile(v, probs = c(alpha / 2, 1 - alpha / 2))
  })
  lwr <- matrix(limits[1, ], n, k)
  upr <- matrix(limits[2, ], n, k)
  if (nrow(newdata) == 1L) {
    ci <- rbind(fit, lwr, upr)
    rownames(ci) <- c("fit", "lwr", "upr")
  } else {
    ci <- array(0, dim = c(n, k, 3))
    dimnames(ci) <- list(1:n, colnames(B_boot[[1]]), c("fit", "lwr", "upr"))
    ci[, , 1] <- fit
    ci[, , 2] <- lwr
    ci[, , 3] <- upr
  }

  return(ci)
}