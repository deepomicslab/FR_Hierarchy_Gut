# TOPOSCORE paper - Analysis script
# v0.2 - JW - Jun 2023

# load required libraries and helper functions
source('tp_helper.R')
library(pROC)
library(randomForest)

## 1. Discovery analysis set ----
log_msg('####### Discovery analysis set #########')
clin_disc <- load_clin(cohort = 'Disc', dir = '../data/immu')
### 1.5 Toposcoring ----
#pred_disc <- read.table('../result/immu/SE/pred_disc.tsv', header = TRUE, sep = '\t')
pred_disc <- read.table('../result/immu/SE/pred_disc.tsv', header = TRUE, sep = '\t')
pred_disc = mutate(pred_disc, OS12=factor(OS12, levels = c('NR', 'R')))
roc <- calc_roc(pred_disc$OS12, pred_disc$TOPOB01, verbose = TRUE)
log_msg('ROC AUC = %.2f [%.2f - %.2f]', roc$AUC[1], roc$AUC[2], roc$AUC[3]) # ROC AUC = 0.54 [0.62 - 0.69]
youden <- roc$ROC_DF %>% mutate(SENS = TPR, SPEC = 1 - FPR) %>% mutate(J = SENS + SPEC - 1)
ggplot(youden, aes(x = THRESHOLD, y = J)) + geom_point()
ycut_nr <- youden[which.max(youden$J), ] # 0.5123106 0.498
ycut_r <- youden[which(youden$THRESHOLD > 0.6 & youden$J > 0.17), ] # 0.5871212 0.6291866
log_msg('Cut-off thresholds = %.4f and %.4f', ycut_nr$THRESHOLD, ycut_r$THRESHOLD) #Cut-off thresholds = 0.5123 and 0.5871  Cut-off thresholds = 0.4976 and 0.6292
plt_roc <- plot_roc(roc$ROC_DF) + 
  geom_point(data = ycut_nr, color = 'red') +
  geom_point(data = ycut_r, color = 'green')
ggsave(plt_roc, filename = '../result/immu/SE/fig_roc.pdf', width = 10, height = 10, units = "cm")
scores_disc = pred_disc %>% select(-OS12)
plt_kde_disc <- plot_toposcoreb01_density(scores_disc, clin_disc, 
                                lims = c(ycut_r$THRESHOLD, ycut_nr$THRESHOLD))
ggsave(plt_kde_disc, filename = '../result/immu/SE/fig_kde_disc.pdf', width = 12, height = 10, units = "cm")


### 1.6 Prediction in discovery cohort (full signature) ----
pred_disc <- read.table('../result/immu/SE/pred_valid.tsv', header = TRUE, sep = '\t')
pred_disc <- assign_prediction(pred_disc, ycut_r$THRESHOLD, ycut_nr$THRESHOLD) # Prediction: NR [47] < 0.50 < gray_zone [93] < 0.63 < R [90]
hr_disc <- get_hr(pred_disc, type = 'OS', by = 'PRED')
#rediction: NR [41] < 0.51 < gray_zone [69] < 0.59 < R [120]
log_msg('Prediction discovery: HR = %.2f [%.2f-%.2f], p = %.1e', hr_disc[1], hr_disc[2], hr_disc[3], hr_disc[4])
# Prediction discovery: HR = 0.47 [0.31-0.72], p = 4.2e-04
#fig_km_disc <- print_plot(plot_mykm(pred_disc, type = 'OS', by = 'PRED'))
#ggsave(fig_km_disc, filename = '../result/immu/SE/fig_km_disc.pdf', width = 12, height = 10, units = "cm")


pred_disc = mutate(pred_disc, OS12=factor(OS12, levels = c('NR', 'R')))
rf_model <- randomForest(OS12 ~ PRED, data = pred_disc)
predictions <- predict(rf_model, pred_disc, type = "response")

roc_obj <- roc(as.numeric(pred_disc$OS12), as.numeric(predictions))
auc_score <- auc(roc_obj)
plot(roc_obj, main = paste0("AUC = ", round(auc_score, 4)))
roc_obj <- pROC::roc(as.numeric(pred_disc$OS12), as.numeric(predictions),
auc = TRUE, ci = TRUE)
roc_df <- data.frame(
THRESHOLD = rev(roc_obj$thresholds),
TPR = rev(roc_obj$sensitivities),
FPR = rev(1 - roc_obj$specificities)
)
list(AUC = as.numeric(roc_obj$ci), ROC_DF = roc_df)
p = ggplot(roc_df, aes(x = FPR, y = TPR)) +
geom_step() +
geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'blue') +
xlim(0, 1) + ylim(0, 1) +
xlab('1 - Specificity') + ylab('Sensitivity') +
GRAPHICS$THEME
ggsave(p, filename = '../result/immu/SE/result.pdf', width = 12, height = 10, units = "cm")