
#LCA model construction----
library(poLCA)
FEATURE_COLS<-c("Obesity","Central_Obesity","Visceral_adiposity",
                "low_HDL","high_TG","high_LDL",
                "high_FPG","high_HbA1c","IR",
                "high_SBP","high_DBP",
                "Liver_Steatosis","high_hsCRP","Hyperuricemia")
df_lca <- df_hcpc[, FEATURE_COLS]


df_lca[FEATURE_COLS] <- lapply(df_lca[FEATURE_COLS], function(x){
  factor(x, levels = c(0,1))
})


f_lca <- as.formula(
  paste0("cbind(", paste(FEATURE_COLS, collapse = ","), ") ~ 1")
)

set.seed(2024)

maxK <- 6
lca_fits <- vector("list", maxK)

for(k in 2:maxK){
  lca_fits[[k]] <- poLCA(
    f_lca,
    data    = df_lca,
    nclass  = k,
    nrep    = 30,     
    maxiter = 5000,
    verbose = FALSE
  )
}


fit_stats <- data.frame(
  Classes = 2:maxK,
  logLik  = sapply(2:maxK, function(k) lca_fits[[k]]$llik),
  AIC     = sapply(2:maxK, function(k) lca_fits[[k]]$aic),
  BIC     = sapply(2:maxK, function(k) lca_fits[[k]]$bic),
  Gsq     = sapply(2:maxK, function(k) lca_fits[[k]]$Gsq),
  Chisq  = sapply(2:maxK, function(k) lca_fits[[k]]$Chisq)
)

fit_stats


bestK <- fit_stats$Classes[which.min(fit_stats$BIC)]
best_model <- lca_fits[[5]]

#numbers of risk factors across LCAs----
#number of risk factors----
df_hcpc$risk_factors_numbers <- rowSums(
  sapply(df_hcpc[, FEATURE_COLS], function(x) as.numeric(as.character(x))),
  na.rm = TRUE
)

summary(df_hcpc$risk_factors_numbers)

#


summary(df_hcpc$LCA_class)

plot_df <- df_hcpc %>%
  dplyr::select(LCA_class, risk_factors_numbers ) %>%
  pivot_longer(
    cols = risk_factors_numbers,
    names_to = "Variable",
    values_to = "Value"
  )

library(ggplot2)

plot_df$LCA_class <- factor(plot_df$LCA_class,
                            levels = c("1","2","3","4","5"),
                            labels = c("Class_1","Class_2","Class_3","Class_4","Class_5"))
class_cols5 <- c(
  "Class_1" = "#1f77b4",
  "Class_2" = "#e377c2",
  "Class_3" = "#ff7f0e",
  "Class_4" = "#17becf",
  "Class_5" = "#2ca02c" 
)

violin5<-ggplot(plot_df, aes(x = LCA_class, y = Value, fill = LCA_class)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  #facet_wrap(~ Variable, scales = "free_y", nrow = 1) +
  stat_compare_means(
    comparisons = list(
      c("Class_1","Class_2"),
      c("Class_1","Class_3"),
      c("Class_1","Class_4"),
      c("Class_1","Class_5")
    ),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  scale_fill_manual(values = class_cols5)+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 12,14))+
  theme_minimal(base_size = 14) +
  ylab("Numbers of risk factors")+
  # labs(
  #   title = "Comparison of Numbers of risk factors"
  # ) +
  theme(legend.title = element_blank(),
        plot.title = element_text(
          size = 14,
          face = "bold",
          hjust = 0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12,angle = 30, hjust = 1, color = "black"),
        strip.text = element_text(face = "bold")
  )

violin5

#Heatmap and radar plot for LCA----

library(dplyr)
prob_mat <- prob_df %>%
  dplyr::select(Feature, Class, Prob) %>%          #  3 cols
  distinct() %>%                             
  pivot_wider(
    id_cols = Feature,                       
    names_from = Class,
    values_from = Prob
  ) %>%
  column_to_rownames("Feature") %>%         
  as.matrix()

prob_mat<-prob_mat[,c(1,2,4,5,3)]


number_mat <- apply(prob_mat,2, function(x){sprintf("%.0f%%", x * 100)})

number_mat
set.seed(2024)



library(ComplexHeatmap)
library(circlize)
library(grid)

colnames(prob_mat)<-paste0(colnames(prob_mat),"\n",pct_label)

m <- prob_mat  # Feature × Class probability matrix  (0-1)

# 1) color
col_fun <- colorRamp2(
  c(0, 0.5, 1),
  c("#F7FBFF", "#2171B5", "#08306B")
)

# 2) row cluster：correlation distance + ward.D2
d_rows  <- as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
hc_rows <- hclust(d_rows, method = "ward.D2")

# 3) cut tree for features
row_k <- 5
row_split <- cutree(hc_rows, k = row_k)

# 4) plot

# pct_label <- c("23.6%", "23.6%", "13.5%", "20.4%", "18.9%")


heatmap_5c<-Heatmap(
  m,
  name = "Prob",
  col = col_fun,
  
  cluster_rows = hc_rows,
  cluster_columns = FALSE,
  
  # #  cutree_rows=5：split
  row_split = 5,
  row_gap = unit(2, "mm"),
  row_title = NULL, 
  
  # row names
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12, just = "centre"),
  column_names_rot = 45, column_names_centered = TRUE,
  column_names_max_height = unit(50, "mm"),
  #bottom_annotation = ha_bottom,
  # 
  rect_gp = gpar(col = "white", lwd = 1),
  
  # title
  column_title = "A.Latent Class Feature Profiles (Conditional probability, %)",
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold",
    just = "left"   # 👈 核心
  ),
  
  # annotation percentage
  cell_fun = function(j, i, x, y, w, h, fill) {
    v <- m[i, j]
    txt <- number_mat[i, j]  # 
    # text color 
    lab_col <- ifelse(v >= 0.60, "white", "black")  # 0.5~0.7 
    grid.text(txt, x, y, gp = gpar(col = lab_col, fontsize = 12))
  }
)
# 1) ComplexHeatmap----grob
library(ComplexHeatmap)
library(grid)
library(ggpubr)
ht_grob <- grid.grabExpr(
  draw(heatmap_5c, padding = unit(c(5, 5, 5, 5), "mm"))  #
)

# 2) grab----ggplot
ht_plot <- as_ggplot(ht_grob) + theme_void() + theme(plot.margin = margin(5,5,5,5))

dev.off()
#radar-----
# 1) 
FEATURE_COLS

label_map <- c(
  Obesity                 = "BMI↑",
  Central_Obesity         = "WC↑",
  Visceral_adiposity      = "VAT↑",
  low_HDL                 = "HDL↓",
  high_TG                 = "TG↑",
  high_LDL                = "LDL↑",
  high_FPG                = "FPG↑",
  high_HbA1c              = "HbA1c↑",
  IR                      = "HOMA-IR↑",
  high_SBP                = "SBP↑",
  high_DBP                = "DBP↑",
  Liver_Steatosis         = "CAP↑",
  high_hsCRP              = "hsCRP↑",
  Hyperuricemia           = "UA↑"
)

FEATURES <- rownames(prob_mat)


freq_df <- df_hcpc[,c(FEATURE_COLS,"LCA_class")] %>%
  tidyr::pivot_longer(
    cols = -LCA_class,
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  dplyr::group_by(LCA_class, Feature) %>%
  dplyr::summarise(
    Freq = mean(Value == 1, na.rm = TRUE),  # 👈 关键
    .groups = "drop"
  )

rad <- freq_df %>%
  dplyr::mutate(
    Feature = factor(Feature, levels = FEATURES),
    Feature_lab = factor(label_map[as.character(Feature)],
                         levels = label_map[FEATURES]),
    Class = factor(LCA_class, levels = 1:5, labels = paste0("Class_", 1:5)),
    Prob_pct = Freq * 100
  ) %>%
  dplyr::arrange(Class, Feature_lab)



# # 3) 
# rad <- prob_df %>%
#   dplyr::mutate(
#     Feature = factor(Feature, levels = FEATURES),
#     Feature_lab = factor(Feature,
#                          levels = FEATURES),
#     Class = factor(Class, levels = paste0("Class_", 1:5)),
#     Prob_pct = Prob * 100
#   ) %>%
#   arrange(Class, Feature_lab)

# 4) 

rad_closed <- rad %>%
  dplyr::group_by(Class) %>%
  dplyr::arrange(Feature_lab, .by_group = TRUE) %>%
  dplyr::reframe(
    Feature_lab = c(as.character(Feature_lab), as.character(Feature_lab)[1]),
    Prob_pct    = c(Prob_pct, Prob_pct[1])
  ) %>%
  dplyr::group_by(Class) %>%              # 👈 关键：重新 group
  dplyr::mutate(idx = dplyr::row_number()) %>%  # 👈 连续编号
  dplyr::ungroup() %>%
  dplyr::mutate(
    Feature_lab = factor(Feature_lab, levels = levels(rad$Feature_lab))
  )

#write.csv(rad,"/Users/jiayufeng/Desktop/NUSWORK/RESET Manuscript/updated_data/rad_reset.csv")
#label
class_labeller <- c(
  Class_1 = "Class 1\nMetabolic preserved",
  Class_2 = "Class 2\nIsolated Hypertension",
  Class_3 = "Class 3\nLean IR/Hyperglycemia",
  Class_4 = "Class 4\nObese IS/Normoglycemia",
  Class_5 = "Class 5\nObese IR/Hyperglycemia"
)



# 5) 
class_cols <- c(
  "Class_1" = "#1f77b4",
  "Class_2" = "#e377c2",
  "Class_3" = "#ff7f0e",
  "Class_4" = "#17becf",
  "Class_5" = "#2ca02c" 
)

# 6) 
K <- length(levels(rad$Feature_lab))   

#install.packages("systemfonts")  
library(systemfonts)

systemfonts::system_fonts() |> head()

#install.packages("showtext")
library(showtext)
showtext_auto() 

p <- ggplot(rad_closed, aes(x = idx, y = Prob_pct, group = Class, color = Class)) +
  geom_polygon(aes(fill = Class), alpha = 0.08, linewidth = 0) +
  geom_line(linewidth = 1.2, lineend = "round") +
  geom_point(size = 1.7) +
  coord_polar(start = 0) +
  #new name
  facet_wrap(
    ~ Class,
    nrow = 3,
    ncol = 2,
    labeller = as_labeller(class_labeller)
  )+
  scale_color_manual(values = class_cols) +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  # 
  scale_x_continuous(
    breaks = 1:K,
    labels = levels(rad$Feature_lab),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.6),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.4),
    axis.title = element_blank(),
    axis.text.y = element_text(color = "grey40", size = 9),
    axis.text.x = element_text(color = "#213a66", size = 8),
    strip.text = element_text(color = "black", size = 12),
    legend.position = "none",
    
    # 
    panel.spacing.x = unit(4, "pt"),
    panel.spacing.y = unit(6, "pt"),
    plot.margin = margin(5, 2, 2, 2),
    
    plot.title = element_text(size = 14, face = "bold",hjust = 0.5)
  )+
  ggtitle("B.Latent Class Feature Profiles (Prevalence, %)")
p



#decision tree model for replicating the LCA----
library(rpart)

df_hcpc$LCA_classC<-as.factor(paste0("Class_",df_hcpc$LCA_class))
summary(df_hcpc$LCA_classC)

tree_fit <- rpart(
  LCA_classC ~ .,
  data = df_hcpc[, c("LCA_classC", FEATURE_COLS)],
  method = "class",
  control = rpart.control(
    maxdepth = 3,      
    minbucket = 30,
    cp = 0.01
  )
)


con_vars <- c("VFat_cm2_1", "Systolic_mmHg_avg", "HOMA_IR_blood", "HbA1c_pct")
summary(df_hcpc[,con_vars])

bin_vars <- c("Visceral_adiposity", "high_SBP", "IR", "high_HbA1c")  # 你树里出现的这些
summary(df_hcpc[,bin_vars])

df_hcpc$pred_class <- predict(tree_fit, type = "class")
summary(df_hcpc$pred_class)

df_tree <- df_hcpc[, c("LCA_classC", FEATURE_COLS)]
df_tree[bin_vars] <- lapply(df_tree[bin_vars], function(x) factor(x, levels = c(0,1), labels = c("No","Yes")))

tree_fit2 <- rpart(
  LCA_classC ~ .,
  data = df_tree[, c("LCA_classC", FEATURE_COLS)],
  method = "class",
  control = rpart.control(maxdepth = 3, minbucket = 30, cp = 0.01)
)

library(rpart.plot)
class_cols <- list(
  "#1f77b4",  # Class 1
  "#e377c2",  # Class 2
  "#ff7f0e",  # Class 3
  "#17becf",  # Class 4
  "#2ca02c"   # Class 5
)


rpart.plot::rpart.plot(
  tree_fit2,
  type = 3,
  extra = 104,
  fallen.leaves = TRUE,
  cex = 1,
  tweak = 0.9,
  box.palette = class_cols,
  shadow.col = "gray",
  main = "Decision tree for latent class classification"
  #nn = TRUE
)

#Validation in PICMAN using Decision tree identified classes----
library(dplyr)

df_hcpc<-df_forest %>% drop_na(all_of(FEATURE_COLS))

df <- df_hcpc %>%
  dplyr::mutate(
    Class5 = case_when(
      # Visceral_adiposity = No
      Visceral_adiposity == 0  & high_SBP == 0  & IR == 0  ~ "1",
      Visceral_adiposity == 0  & high_SBP == 0  & IR == 1 ~ "3",
      Visceral_adiposity == 0  & high_SBP == 1 & IR == 0  ~ "2",
      Visceral_adiposity == 0  & high_SBP == 1 & IR == 1 ~ "3",
      
      # Visceral_adiposity = Yes
      Visceral_adiposity == 1 & IR == 0 ~ "4",
      Visceral_adiposity == 1 & IR == 1 & high_HbA1c == 0  ~ "4",
      Visceral_adiposity == 1 & IR == 1 & high_HbA1c == 1 ~ "5"
    ),
    Class5 = factor(Class5, levels = 1:5)
  )

summary(df$Class5)



#Limma test for identifying proteomic signature associated with VAT, HOMA-IR, SBP or HbA1c----
limma_linear <- function(data, x) {
  design <- model.matrix(~x)
  fit <- lmFit(data, design)
  fit2 <- eBayes(fit)
  res <- topTable(fit2, sort.by = "logFC", number = Inf)
  sigres <- subset(res, P.Value < 0.05)
  sigres <- sigres[order(-abs(sigres$logFC),sigres$P.Value), ]
  return(sigres)
}



summary(multi_omic[,c("HOMA_IR","vat","sbp","hb_a1c_base_measurements")])


results <- list()
cont_vars_IR_obese<-c("HOMA_IR","vat","sbp","hb_a1c_base_measurements")

for (data_name in names(omics_list)) {
  data <- omics_list[[data_name]]
  
  for (var in cont_vars_IR_obese) {
    x <- as.vector(multi_omic[[var]])  
    key <- paste(data_name, var, sep = "_")
    results[[key]] <- limma_linear(data, x)
  }
}

#Limma test for 5 classes proteomic comparison, using class 1 as reference----


limma_pairwise_cov <- function(data, pheno) {

  design <- model.matrix(
    as.formula(paste("~ 0 + group")),
    data = pheno
  )
  colnames(design)[1:length(levels(pheno$group))] <- levels(pheno$group)
  

  contrast_matrix <- makeContrasts(
    Isolated_HTN_vs_Metabolic_Preserved = Class_2 - Class_1,
    Lean_IR_vs_Metabolic_Preserved = Class_3 - Class_1,
    Obese_IS_vs_Metabolic_Preserved = Class_4 - Class_1,
    Lean_IR_vs_Obese_IS = Class_3 - Class_4,
    Obese_IR_vs_Metabolic_Preserved = Class_5 - Class_1,
    Obese_IR_vs_Lean_IR = Class_5 - Class_3,
    Obese_IR_vs_Obese_IS = Class_5 - Class_4,
    levels = design
  )
  

  fit <- lmFit(data, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  

  res_Isolated_HTN_vs_Metabolic_Preserved<- topTable(fit2, coef = "Isolated_HTN_vs_Metabolic_Preserved", number = Inf)
  res_Lean_IR_vs_Metabolic_Preserved  <- topTable(fit2, coef = "Lean_IR_vs_Metabolic_Preserved", number = Inf)
  res_Obese_IS_vs_Metabolic_Preserved <- topTable(fit2, coef = "Obese_IS_vs_Metabolic_Preserved", number = Inf)
  res_Lean_IR_vs_Obese_IS <- topTable(fit2, coef = "Lean_IR_vs_Obese_IS", number = Inf)
  res_Obese_IR_vs_Metabolic_Preserved <- topTable(fit2, coef = "Obese_IR_vs_Metabolic_Preserved", number = Inf)
  res_Obese_IR_vs_Lean_IR <- topTable(fit2, coef = "Obese_IR_vs_Lean_IR", number = Inf)
  res_Obese_IR_vs_Obese_IS <- topTable(fit2, coef = "Obese_IR_vs_Obese_IS", number = Inf)
  
  return(list(
    Isolated_HTN_vs_Metabolic_Preserved = res_Isolated_HTN_vs_Metabolic_Preserved[order(-abs(res_Isolated_HTN_vs_Metabolic_Preserved$logFC), res_Isolated_HTN_vs_Metabolic_Preserved$adj.P.Val), ],
    Lean_IR_vs_Metabolic_Preserved = res_Lean_IR_vs_Metabolic_Preserved[order(-abs(res_Lean_IR_vs_Metabolic_Preserved$logFC), res_Lean_IR_vs_Metabolic_Preserved$adj.P.Val), ],
    Obese_IS_vs_Metabolic_Preserved = res_Obese_IS_vs_Metabolic_Preserved[order(-abs(res_Obese_IS_vs_Metabolic_Preserved$logFC), res_Obese_IS_vs_Metabolic_Preserved$adj.P.Val), ],
    Lean_IR_vs_Obese_IS = res_Lean_IR_vs_Obese_IS[order(-abs(res_Lean_IR_vs_Obese_IS$logFC), res_Lean_IR_vs_Obese_IS$adj.P.Val), ],
    Obese_IR_vs_Metabolic_Preserved = res_Obese_IR_vs_Metabolic_Preserved[order(-abs(res_Obese_IR_vs_Metabolic_Preserved$logFC), res_Obese_IR_vs_Metabolic_Preserved$adj.P.Val), ],
    Obese_IR_vs_Lean_IR = res_Obese_IR_vs_Lean_IR[order(-abs(res_Obese_IR_vs_Lean_IR$logFC), res_Obese_IR_vs_Lean_IR$adj.P.Val), ],
    Obese_IR_vs_Obese_IS = res_Obese_IR_vs_Obese_IS[order(-abs(res_Obese_IR_vs_Obese_IS$logFC), res_Obese_IR_vs_Obese_IS$adj.P.Val), ]
  ))
}



omics_list_linear <- Map(
  function(data_matrix, linear_var) {
    data_matrix[row.names(data_matrix) %in% linear_var, , drop = FALSE]
  },
  omics_list,
  all_list_linear
)



#all, 
results_list_pair <- lapply(omics_list_linear, function(data_matrix) {
  limma_pairwise_cov(data_matrix, pheno)
})

results_list_pair_all2<-unlist(results_list_pair, recursive = FALSE)

#sig

results_list_pair_sig2<-lapply(results_list_pair_all2, function(df){
  df <- subset(df, P.Value < 0.05)
  df <- df[order(-abs(df$logFC),df$P.Value), ]
})

# class 3 and class 4, discordant phenotypes proteomic signature (validated in UKB)----


ukb_leanir <- results_list_ukb_sig2$res_B_Lean_IR_vs_Metabolic_Preserved %>%
  select(Variable = Protein, logFC_UKB_LeanIR = logFC, adj.P.Val_UKB_LeanIR = adj.P.Val)

ukb_obeseis <- results_list_ukb_sig2$res_C_Obese_IS_vs_Metabolic_Preserved %>%
  select(Variable = Protein, logFC_UKB_ObeseIS = logFC, adj.P.Val_UKB_ObeseIS = adj.P.Val)

# ---- 1. join UKB results onto your discovery discordance table ----
discordance_ukb <- discordance %>%
  left_join(ukb_leanir, by = "Variable") %>%
  left_join(ukb_obeseis, by = "Variable")

# ---- 2. category-aware validation flag ----
discordance_ukb <- discordance_ukb %>%
  mutate(
    sig_ukb_leanir  = !is.na(adj.P.Val_UKB_LeanIR)  & adj.P.Val_UKB_LeanIR  < 0.05,
    sig_ukb_obeseis = !is.na(adj.P.Val_UKB_ObeseIS) & adj.P.Val_UKB_ObeseIS < 0.05,

    dir_concord_leanir  = sign(logFC_LeanIR)  == sign(logFC_UKB_LeanIR),
    dir_concord_obeseis = sign(logFC_ObeseIS) == sign(logFC_UKB_ObeseIS),
    validated_ukb = case_when(
      category == "Unique_LeanIR"  ~ sig_ukb_leanir  & dir_concord_leanir,
      category == "Unique_ObeseIS" ~ sig_ukb_obeseis & dir_concord_obeseis,
      startsWith(category, "Common") ~ sig_ukb_leanir & sig_ukb_obeseis &
        dir_concord_leanir & dir_concord_obeseis,
      TRUE ~ FALSE
    )
  )

# ---- 3. final validated discordance protein set ----
validated_discordance <- discordance_ukb %>%
  filter(validated_ukb) %>%
  select(Variable, category, logFC_LeanIR, logFC_ObeseIS,
         P.Value_LeanIR, P.Value_ObeseIS, validated_ukb)

validated_discordance



# ----  curated biological theme lookup (edit/expand as needed) ----
module_lookup <- tibble::tribble(
  ~Variable,   ~module,
  # --- Hepatic metabolism & injury---
  "ADH4",      "Hepatic_Metabolism_Injury",
  "ADH1B",     "Hepatic_Metabolism_Injury",
  "GSTA1",     "Hepatic_Metabolism_Injury",
  "CA5A",      "Hepatic_Metabolism_Injury",
  "FTCD",      "Hepatic_Metabolism_Injury",
  "DCXR",      "Hepatic_Metabolism_Injury",
  "HAO1",      "Hepatic_Metabolism_Injury",
  "KRT8",      "Hepatic_Metabolism_Injury",   # 
  "KRT18",     "Hepatic_Metabolism_Injury",   # 
  # --- Inflammation & immune(Lean-IR )---
  "PLCB2",     "Inflammatory_signaling",
  "IL6",       "Inflammatory_signaling",
  "SERPINE1",  "Inflammatory_signaling",
  # --- Insulin & endocrine(SHARED)---
  "LEP",       "Insulin_Endocrine_Signaling",
  "IGFBP1",    "Insulin_Endocrine_Signaling",
  "FGF21",     "Insulin_Endocrine_Signaling",
  
  # --- Mitochondrial & energy metabolism ---
  "RTN4IP1",   "Other_Cellular",
  "ECHS1",     "Other_Cellular",         
  "ECHDC3",    "Other_Cellular",
  "AIFM1",     "Other_Cellular",
  
  # --- Nucleic-acid & protein regulation---
  "GTPBP2",    "Other_Cellular",
  "PARP1",     "Other_Cellular",
  "EIF4E",     "Other_Cellular",
  "ELOA",      "Other_Cellular",
  "ITPA",      "Other_Cellular",   
  # --- Cytoskeletal & structural ---
  "IGSF9",     "Other_Cellular",
  "ESYT2",     "Other_Cellular",
  "GAS2",      "Other_Cellular",
  # --- Other / unclassified---
  "OXT",       "Other_Cellular",
  "DBH",       "Other_Cellular",          
  "DKKL1",     "Other_Cellular",          
  "TREH",      "Other_Cellular"           
)

module_labels <- c(
  "Hepatic_Metabolism_Injury"  = "Hepatic metabolism & injury",
  "Inflammatory_signaling"            = "Inflammatory signaling",
  "Insulin_Endocrine_Signaling"    = "Insulin & endocrine signaling",
  "Other_Cellular"            = "Other cellular processes"
)
module_colors <- c(
  "Hepatic_Metabolism_Injury"  = "#16a085",  # orange — pulled away from red
  "Inflammatory_signaling"            = "#e63946",  # crimson red — kept distinct from orange
  "Insulin_Endocrine_Signaling"    = "#8e44ad",  # purple
  "Other_Cellular"   = "#7f8c8d"
)

# ---- 2. build long panel data ----
make_panel_df <- function(data, logfc_col, pval_col) {
  data %>%
    filter(!is.na(.data[[logfc_col]]), !is.na(.data[[pval_col]])) %>%
    transmute(
      Variable,
      logFC     = .data[[logfc_col]],
      neglogp   = -log10(.data[[pval_col]]),
      is_common = startsWith(category, "Common")
    ) %>%
    left_join(module_lookup, by = "Variable")
}

lean_df  <- make_panel_df(df31, "logFC_LeanIR",  "P.Value_LeanIR")
obese_df <- make_panel_df(df31, "logFC_ObeseIS", "P.Value_ObeseIS")



core_modules <- c("Hepatic_Metabolism_Injury",
                  "Inflammatory_signaling",
                  "Insulin_Endocrine_Signaling")

core_df <- df31 %>%
  dplyr::left_join(module_lookup, by = "Variable") %>%
  dplyr::filter(module %in% core_modules) %>%
  dplyr::mutate(
    module_label = factor(module_labels[module], levels = module_labels[core_modules]),
    x = dplyr::coalesce(logFC_LeanIR,  0),          # NA(非显著)→ 0,贴 x 轴
    y = dplyr::coalesce(logFC_ObeseIS, 0),          # NA(非显著)→ 0,贴 y 轴
    sig_grp = dplyr::case_when(
      category == "Common_Same"     ~ "Both",
      category == "Common_Opposite" ~ "Both",
      category == "Unique_LeanIR"   ~ "Lean-IR only",
      category == "Unique_ObeseIS"  ~ "Obese-IS only",
      TRUE ~ "NS")
  )

rng <- max(abs(c(core_df$x, core_df$y))) * 1.05

px<-ggplot(core_df, aes(x, y)) +
  geom_hline(yintercept = 0, color = "grey75") +
  geom_vline(xintercept = 0, color = "grey75") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +  # 共享同向
  geom_point(aes(color = module_label, shape = sig_grp), size = 2, stroke = 0.7) +
  geom_text_repel(aes(label = Variable, color = module_label),
                  size = 4.5, max.overlaps = Inf, box.padding = 0.4,
                  segment.size = 0.3, min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual(values = setNames(module_colors[core_modules], module_labels[core_modules]),
                     name = "Pathway") +
  scale_shape_manual(values = c(Both = 16, `Lean-IR only` = 17, `Obese-IS only` = 15),
                     name = "Significant in") +
  coord_equal(xlim = c(-rng, rng), ylim = c(-rng, rng)) +
  labs(x = expression(log[2]*"FC  Lean-IR vs Metabolically Preserved"),
       y = expression(log[2]*"FC  Obese-IS vs Metabolically Preserved"),
       #title = "Core-pathway proteins across both discordant phenotypes"
  )+
  scale_x_continuous(breaks = seq(-1.5, 1.5, by = 0.5))+
  scale_y_continuous(breaks = seq(-1.5, 1.5, by = 0.5))+
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.location = "panel",
        legend.justification = "center",
        legend.spacing.y = unit(0.2, "cm"),           
        legend.margin = margin(t = 0.1, b = 0.1),       
        legend.box.spacing = unit(0.2, "cm"),     
        plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"))+ guides(color = guide_legend(nrow = 2, byrow = TRUE))

px


# Heatmap core-Proteomic signature for Class_1 / Class_3 / Class_4 ----


core_keep <- df_plot$module %in% core_modules
core_prot <- df_plot$protein[core_keep]


sel_cols <- c("Class_1","Class_3","Class_4")

mat_sub  <- mat_scaled[core_prot, sel_cols, drop = FALSE]


row_split_sub <- droplevels(df_plot$module_label[core_keep])


sig_star <- function(p) ifelse(is.na(p), "",
                               ifelse(p < 0.001, "***",
                                      ifelse(p < 0.01 , "**",
                                             ifelse(p < 0.05 , "*", ""))))
prot <- rownames(mat_sub)
star_mat <- matrix("", nrow = nrow(mat_sub), ncol = 3,
                   dimnames = list(prot, sel_cols))
star_mat[, "Class_3"] <- sig_star(df31$P.Value_LeanIR [match(prot, df31$Variable)])  # Lean-IR
star_mat[, "Class_4"] <- sig_star(df31$P.Value_ObeseIS[match(prot, df31$Variable)])  # Obese-IS



labels_wrapped <- gsub(" & ", " &\n", levels(row_split_sub))   
ha_row <- rowAnnotation(
  ModuleColor = row_split_sub,
  ModuleLabel = anno_block(
    labels = labels_wrapped,
    gp = gpar(fill = NA, col = NA),
    labels_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_rot = 0
  ),
  col = list(ModuleColor = module_col),  
  show_annotation_name = FALSE,
  show_legend = FALSE,
  annotation_width = unit.c(unit(4, "mm"), unit(45, "mm"))
)

ht <- Heatmap(
  mat_sub,
  name = "z-score",
  col  = colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b40426")),

  right_annotation = ha_row,
  row_split = row_split_sub,
  row_title = NULL,
  show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 12),
  

  show_column_names = TRUE,
  column_labels = c("Metabolically\nPreserved", "Lean-IR", "Obese-IS"),  
  column_names_rot = 15,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  
  cluster_rows = TRUE, cluster_row_slices = FALSE, cluster_columns = FALSE,
  row_gap = unit(1, "mm"), border = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(star_mat[i, j], x, y, gp = gpar(fontsize = 12, fontface = "bold"))
  }
  #column_title = "Proteomic signatures in discordant phenotypes",
  #column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
ht

#Partial Linear correlation for core proteins and VAT+HOMA-IR----

df_linear <- as.data.frame(expression_matrix_prot) %>%
  rownames_to_column("protein") %>%
  dplyr::left_join(module_lookup[module_lookup$module %in% core_modules,], by = c("protein" = "Variable")) %>%
  dplyr::filter(!is.na(module)) %>%
  dplyr::mutate(module_label = factor(module_labels[module],
                                      levels = display_levels)) %>%
  dplyr::arrange(module_label, protein)

colnames(df_forest)

df_outcome<-df_forest[,c("subject_id","vat","HOMA_IR")]


# 1. Extract expression data from df_linear
# Identify sample columns (exclude protein and module columns)
meta_cols <- c("protein", "module", "module_label")
sample_cols <- setdiff(colnames(df_linear), meta_cols)

# Extract just the expression data
expr_data <- df_linear[, sample_cols]
rownames(expr_data) <- df_linear$protein

# 2. Transpose and merge with phenotype
expr_transposed <- t(expr_data) %>% as.data.frame() %>%
  rownames_to_column("subject_id") %>%
  dplyr::left_join(df_outcome, by = "subject_id")

# 3. Calculate correlation matrix
protein_names <- df_linear$protein
cor_matrix <- cor(expr_transposed[, protein_names],
                  expr_transposed[, c("vat", "HOMA_IR")],
                  use = "complete.obs",
                  method = "pearson")

# Calculate p-values
library(Hmisc)
p_values <- rcorr(as.matrix(expr_transposed[, protein_names]),
                  as.matrix(expr_transposed[, c("vat", "HOMA_IR")]))$P[protein_names, c("vat", "HOMA_IR")]


# 4. Partial correlations (protein vs each trait, adjusting for the other)
library(ppcor)
vars2 <- c("vat", "HOMA_IR")

cor_matrix <- matrix(NA_real_, nrow = length(protein_names), ncol = 2,
                     dimnames = list(protein_names, vars2))
p_values   <- matrix(NA_real_, nrow = length(protein_names), ncol = 2,
                     dimnames = list(protein_names, vars2))

for (p in protein_names) {
  sub <- expr_transposed[, c(p, "vat", "HOMA_IR")]
  sub <- sub[complete.cases(sub), , drop = FALSE]
  if (nrow(sub) < 10) next
  xp <- as.numeric(sub[, p]); xv <- as.numeric(sub[, "vat"]); xi <- as.numeric(sub[, "HOMA_IR"])
  pc_vat <- ppcor::pcor.test(xp, xv, xi)   # protein vs VAT | HOMA_IR
  pc_ir  <- ppcor::pcor.test(xp, xi, xv)   # protein vs HOMA_IR | VAT
  cor_matrix[p, "vat"]     <- pc_vat$estimate
  cor_matrix[p, "HOMA_IR"] <- pc_ir$estimate
  p_values[p, "vat"]       <- pc_vat$p.value
  p_values[p, "HOMA_IR"]   <- pc_ir$p.value
}


