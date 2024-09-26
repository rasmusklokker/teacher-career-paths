## 4. Examples with life course data
## 4.1. Sequence data

library("seqHMM")

data("biofam", package = "TraMineR")
biofam_seq <- seqdef(biofam[, 10:25], start = 15, labels = c("parent",
  "left", "married", "left+marr", "child", "left+child", "left+marr+ch",
  "divorced"))

data("biofam3c", package = "seqHMM")
marr_seq <- seqdef(biofam3c$married, start = 15, alphabet = c("single",
  "married", "divorced"))
child_seq <- seqdef(biofam3c$children, start = 15,
  alphabet = c("childless", "children"))
left_seq <- seqdef(biofam3c$left, start = 15, alphabet = c("with parents",
  "left home"))

attr(marr_seq, "cpal") <- c("violetred2", "darkgoldenrod2", "darkmagenta")
attr(child_seq, "cpal") <- c("darkseagreen1", "coral3")
attr(left_seq, "cpal") <- c("lightblue", "red3")


## Figure 2
ssplot(list("Marriage" = marr_seq, "Parenthood" = child_seq,
  "Residence" = left_seq))

## Figure 1
seq_data <- list(biofam_seq[1:10, ], marr_seq[1:10, ], child_seq[1:10, ],
  left_seq[1:10, ])
ssplot(seq_data, type = "I", sortv = "from.start", sort.channel = 1,
  ylab = c("Original", "Marriage", "Parenthood", "Residence"),
  xtlab = 15:30, xlab = "Age", ylab.pos = c(1, 1.5), title.n = FALSE,
  title = "Ten first sequences", legend.prop = 0.63,
  ncol.legend = c(3, 1, 1, 1))

## Figure 3
ssp_f <- ssp(list(marr_seq[biofam3c$covariates$sex == "woman", ],
    child_seq[biofam3c$covariates$sex == "woman", ],
    left_seq[biofam3c$covariates$sex == "woman", ]),
  type = "I", sortv = "mds.obs", with.legend = FALSE, title = "Women",
  ylab.pos = c(1, 2, 1), xtlab = 15:30, ylab = c("Married", "Children",
    "Residence"))

ssp_m <- update(ssp_f, title = "Men",
  x = list(marr_seq[biofam3c$covariates$sex == "man", ],
    child_seq[biofam3c$covariates$sex == "man", ],
    left_seq[biofam3c$covariates$sex == "man", ]))

gridplot(list(ssp_f, ssp_m), ncol = 2, nrow = 2, byrow = TRUE,
  legend.pos = "bottom", legend.pos2 = "top", row.prop = c(0.65, 0.35))


## 4.2. Hidden Markov models

## Single-channel model with automatic starting values
sc_initmod_random <- build_hmm(observations = biofam_seq, n_states = 5)

## Single-channel model with user-defined starting values
sc_init <- c(0.9, 0.06, 0.02, 0.01, 0.01)

sc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0.02, 0.80, 0.10,
  0.05, 0.03, 0.02, 0.03, 0.80, 0.10, 0.05, 0.02, 0.03, 0.05, 0.80, 0.10,
  0.02, 0.03, 0.05, 0.05, 0.85), nrow = 5, ncol = 5, byrow = TRUE)

sc_emiss <- matrix(NA, nrow = 5, ncol = 8)
sc_emiss[1, ] <- seqstatf(biofam_seq[, 1:4])[, 2] + 0.1
sc_emiss[2, ] <- seqstatf(biofam_seq[, 5:7])[, 2] + 0.1
sc_emiss[3, ] <- seqstatf(biofam_seq[, 8:10])[, 2] + 0.1
sc_emiss[4, ] <- seqstatf(biofam_seq[, 11:13])[, 2] + 0.1
sc_emiss[5, ] <- seqstatf(biofam_seq[, 14:16])[, 2] + 0.1
sc_emiss <- sc_emiss / rowSums(sc_emiss)

rownames(sc_trans) <- colnames(sc_trans) <- rownames(sc_emiss) <-
  paste("State", 1:5)

colnames(sc_emiss) <- attr(biofam_seq, "labels")

sc_trans
round(sc_emiss, 3)

sc_initmod <- build_hmm(observations = biofam_seq, initial_probs = sc_init,
  transition_probs = sc_trans, emission_probs = sc_emiss)

sc_fit <- fit_model(sc_initmod)

sc_fit$logLik

sc_fit$model

## Multi-channel model
mc_init <- c(0.9, 0.05, 0.02, 0.02, 0.01)

mc_trans <- matrix(c(0.80, 0.10, 0.05, 0.03, 0.02, 0, 0.90, 0.05, 0.03,
  0.02, 0, 0, 0.90, 0.07, 0.03, 0, 0, 0, 0.90, 0.10, 0, 0, 0, 0, 1),
  nrow = 5, ncol = 5, byrow = TRUE)

mc_emiss_marr <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05, 0.90,
  0.05, 0.05, 0.90, 0.05, 0.30, 0.30, 0.40), nrow = 5, ncol = 3,
  byrow = TRUE)

mc_emiss_child <- matrix(c(0.9, 0.1, 0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.5,
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_emiss_left <- matrix(c(0.9, 0.1, 0.1, 0.9, 0.1, 0.9, 0.1, 0.9, 0.5,
  0.5), nrow = 5, ncol = 2, byrow = TRUE)

mc_obs <- list(marr_seq, child_seq, left_seq)

mc_emiss <- list(mc_emiss_marr, mc_emiss_child, mc_emiss_left)

mc_initmod <- build_hmm(observations = mc_obs, initial_probs = mc_init,
  transition_probs = mc_trans, emission_probs = mc_emiss,
  channel_names = c("Marriage", "Parenthood", "Residence"))

mc_initmod

mc_fit <- fit_model(mc_initmod, em_step = FALSE, local_step = TRUE,
  threads = 4)

hmm_biofam <- mc_fit$model
BIC(hmm_biofam)

## 4.3. Clustering and mixture hidden Markov models
mc_init2 <- c(0.9, 0.05, 0.03, 0.02)

mc_trans2 <- matrix(c(0.85, 0.05, 0.05, 0.05, 0, 0.90, 0.05, 0.05, 0, 0,
  0.95, 0.05, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = TRUE)

mc_emiss_marr2 <- matrix(c(0.90, 0.05, 0.05, 0.90, 0.05, 0.05, 0.05,
  0.85, 0.10, 0.05, 0.80, 0.15), nrow = 4, ncol = 3, byrow = TRUE)

mc_emiss_child2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mc_emiss_left2 <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  nrow = 4, ncol = 2, byrow = TRUE)

mhmm_init <- list(mc_init, mc_init2)

mhmm_trans <- list(mc_trans, mc_trans2)

mhmm_emiss <- list(list(mc_emiss_marr, mc_emiss_child, mc_emiss_left),
  list(mc_emiss_marr2, mc_emiss_child2, mc_emiss_left2))

biofam3c$covariates$cohort <- cut(biofam3c$covariates$birthyr,
  c(1908, 1935, 1945, 1957))
biofam3c$covariates$cohort <- factor(biofam3c$covariates$cohort,
  labels=c("1909-1935", "1936-1945", "1946-1957"))

init_mhmm <- build_mhmm(observations = mc_obs, initial_probs = mhmm_init,
  transition_probs = mhmm_trans, emission_probs = mhmm_emiss,
  formula = ~sex + cohort, data = biofam3c$covariates,
  channel_names = c("Marriage", "Parenthood", "Residence"),
  cluster_names = c("Cluster 1", "Cluster 2"))

set.seed(1011)
mhmm_fit <- fit_model(init_mhmm, local_step = TRUE, threads = 4,
  control_em = list(restart = list(times = 100)))
mhmm <- mhmm_fit$model

summary(mhmm, conditional_se = FALSE)

## 4.4. Visualizing hidden Markov models
plot(hmm_biofam)

## Figure 4
plot(hmm_biofam, vertex.size = 50, vertex.label.dist = 4.5,
     edge.curved = c(0, 0.6, -0.8, 0.6, 0, 0.6, 0), legend.prop = 0.3,
     combined.slice.label = "States with prob. < 0.05")

## Figure 5
vertex_layout <- matrix(c(1, 2, 2, 3, 1, 0, 0.5, -0.5, 0, -1),
                        ncol = 2)
plot(hmm_biofam, layout = vertex_layout, xlim = c(0.5, 3.5),
     ylim = c(-1.5, 1), rescale = FALSE, vertex.size = 50,
     vertex.label.pos = c("left", "top", "bottom", "right", "left"),
     edge.curved = FALSE, edge.width = 1, edge.arrow.size = 1,
     with.legend = "left", legend.prop = 0.4, label.signif = 1,
     combine.slices = 0, cpal = colorpalette[[30]][14:5])

## Figure 6
ssplot(hmm_biofam, plots = "both", type = "I", sortv = "mds.hidden",
       title = "Observed and hidden state sequences", xtlab = 15:30,
       xlab = "Age")

## 4.5. Visualizing mixture hidden Markov models

## MHMM plot in a grid
plot(mhmm, interactive = FALSE, nrow = 2, legend.prop = 0.45,
  vertex.size = 50, vertex.label.cex = 1.3, cex.legend = 1.3,
  edge.curved = 0.65, edge.label.cex = 1.3, cex.edge.width = 0.8,
  edge.arrow.size = 0.8)

## Interactive MHMM plot
if (interactive()) {
    mssplot(mhmm, ask = TRUE)
}

