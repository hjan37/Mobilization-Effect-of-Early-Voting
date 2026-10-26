# =============================================================================
# 표 8 재현: 사전투표자의 동원 결정요인 — 이항 로지스틱 GLMM
# =============================================================================
# 논문: 사전투표제도의 투표율 제고 효과
#       제22대 국회의원선거에 대한 가설적 설문 문항과
#       위계적 일반화 선형모형(GLMM)을 활용한 실증 분석
# 투고 저널: 한국정치연구
# =============================================================================

# --- 0. 패키지 ---
if (!require("lme4")) install.packages("lme4", repos = "https://cloud.r-project.org")
library(lme4)

# --- 1. 데이터 읽기 ---
d <- read.csv("data/early_voters_n319.csv", fileEncoding = "UTF-8")

cat("N =", nrow(d), "\n")
cat("동원:", sum(d$mobilized == 1), "/ 대체:", sum(d$mobilized == 0), "\n")
cat("동원 비율:", round(mean(d$mobilized), 3), "\n")
cat("시·도 수:", length(unique(d$region)), "\n\n")

# 그룹 변수를 factor로
d$region_f <- factor(d$region)

# =============================================================================
# 모형 1: 정당일체감 유무 (Q19 기준)
# =============================================================================
cat(strrep("=", 60), "\n")
cat("모형 1: 정당일체감 유무\n")
cat(strrep("=", 60), "\n\n")

m1 <- glmer(
  mobilized ~ female + age10 + pol_interest + elec_interest +
    education + income + ideology + party_ident_q19 +
    (1 | region_f),
  data   = d,
  family = binomial(link = "logit"),
  nAGQ   = 1  # Laplace approximation
)

print(summary(m1), correlation = FALSE)

cat("\n오즈비 (OR):\n")
or1 <- exp(fixef(m1))
ci1 <- exp(confint(m1, parm = "beta_", method = "Wald"))
print(round(cbind(OR = or1, ci1), 3))

cat("\n적합도 지표:\n")
cat("  Log-Likelihood:", round(as.numeric(logLik(m1)), 1), "\n")
cat("  AIC:           ", round(AIC(m1), 1), "\n")
cat("  BIC:           ", round(BIC(m1), 1), "\n")

vc1 <- as.data.frame(VarCorr(m1))
cat("  τ₀₀ (region):  ", round(vc1$vcov, 6), "\n")
icc1 <- vc1$vcov / (vc1$vcov + pi^2 / 3)
cat("  ICC:           ", round(icc1, 4), "\n")

# =============================================================================
# 모형 2: 정당별 일체감 (Q19B 기준)
# =============================================================================
cat("\n")
cat(strrep("=", 60), "\n")
cat("모형 2: 정당별 일체감\n")
cat(strrep("=", 60), "\n\n")

m2 <- glmer(
  mobilized ~ female + age10 + pol_interest + elec_interest +
    education + income + ideology + pid_dem + pid_ppp + pid_other +
    (1 | region_f),
  data   = d,
  family = binomial(link = "logit"),
  nAGQ   = 1
)

print(summary(m2), correlation = FALSE)

cat("\n오즈비 (OR):\n")
or2 <- exp(fixef(m2))
ci2 <- exp(confint(m2, parm = "beta_", method = "Wald"))
print(round(cbind(OR = or2, ci2), 3))

cat("\n적합도 지표:\n")
cat("  Log-Likelihood:", round(as.numeric(logLik(m2)), 1), "\n")
cat("  AIC:           ", round(AIC(m2), 1), "\n")
cat("  BIC:           ", round(BIC(m2), 1), "\n")

vc2 <- as.data.frame(VarCorr(m2))
cat("  τ₀₀ (region):  ", round(vc2$vcov, 6), "\n")
icc2 <- vc2$vcov / (vc2$vcov + pi^2 / 3)
cat("  ICC:           ", round(icc2, 4), "\n")

# =============================================================================
# 모형 비교 및 논문 수치 대조
# =============================================================================
cat("\n")
cat(strrep("=", 60), "\n")
cat("논문 표 8 수치 대조\n")
cat(strrep("=", 60), "\n\n")

cat("적합도 비교:\n")
cat(sprintf("  모형 1  LL = %6.1f (논문: -200.9)  AIC = %5.1f (논문: 421.8)  BIC = %5.1f (논문: 459.5)\n",
            as.numeric(logLik(m1)), AIC(m1), BIC(m1)))
cat(sprintf("  모형 2  LL = %6.1f (논문: -201.5)  AIC = %5.1f (논문: 427.0)  BIC = %5.1f (논문: 472.2)\n",
            as.numeric(logLik(m2)), AIC(m2), BIC(m2)))

# 핵심 계수 비교
cat("\n주요 계수 비교 (모형 1):\n")
cat(sprintf("  정당일체감 있음  b = %6.3f (논문: -0.483)  OR = %4.2f (논문: 0.62)\n",
            fixef(m1)["party_ident_q19"], exp(fixef(m1)["party_ident_q19"])))

cat("\n주요 계수 비교 (모형 2):\n")
cat(sprintf("  민주당           b = %6.3f (논문: -0.106)  OR = %4.2f (논문: 0.90)\n",
            fixef(m2)["pid_dem"], exp(fixef(m2)["pid_dem"])))
cat(sprintf("  국민의힘         b = %6.3f (논문: -0.298)  OR = %4.2f (논문: 0.74)\n",
            fixef(m2)["pid_ppp"], exp(fixef(m2)["pid_ppp"])))
cat(sprintf("  기타정당         b = %6.3f (논문:  0.522)  OR = %4.2f (논문: 1.69)\n",
            fixef(m2)["pid_other"], exp(fixef(m2)["pid_other"])))

# 비내포 구조 설명
cat("\n")
cat(strrep("-", 60), "\n")
cat("비내포(non-nested) 구조 확인\n")
cat(strrep("-", 60), "\n")
cat("  모형 1 '정당일체감 유무' (Q19 기준):     ", sum(d$party_ident_q19), "명\n")
cat("  모형 2  정당별 더미     (Q19B 기준):     ", sum(d$pid_dem + d$pid_ppp + d$pid_other), "명\n")
diff_n <- sum(d$pid_dem + d$pid_ppp + d$pid_other) - sum(d$party_ident_q19)
cat("  Q19='없다' & Q19A='있다' 불일치 인원:   ", diff_n, "명\n")
cat("  → 두 모형은 엄밀한 내포관계가 아님\n")
cat("  → LL(모형 1) > LL(모형 2)는 통계적으로 문제없음\n")
cat("  → AIC/BIC 기준 비교가 적절\n")

cat("\n=== 재현 완료 ===\n")

# =============================================================================
# 결과 저장
# =============================================================================
if (!dir.exists("output")) dir.create("output")

# 고정효과 테이블
fe1 <- data.frame(
  model    = "Model 1",
  variable = names(fixef(m1)),
  b        = round(fixef(m1), 3),
  se       = round(sqrt(diag(vcov(m1))), 3),
  OR       = round(exp(fixef(m1)), 2),
  p        = round(summary(m1)$coefficients[, "Pr(>|z|)"], 3)
)
fe2 <- data.frame(
  model    = "Model 2",
  variable = names(fixef(m2)),
  b        = round(fixef(m2), 3),
  se       = round(sqrt(diag(vcov(m2))), 3),
  OR       = round(exp(fixef(m2)), 2),
  p        = round(summary(m2)$coefficients[, "Pr(>|z|)"], 3)
)
results <- rbind(fe1, fe2)
row.names(results) <- NULL
write.csv(results, "output/table8_coefficients.csv", row.names = FALSE)

# 적합도 지표
fit <- data.frame(
  model = c("Model 1", "Model 2"),
  LL    = round(c(as.numeric(logLik(m1)), as.numeric(logLik(m2))), 1),
  AIC   = round(c(AIC(m1), AIC(m2)), 1),
  BIC   = round(c(BIC(m1), BIC(m2)), 1),
  tau00 = round(c(vc1$vcov, vc2$vcov), 6),
  ICC   = round(c(icc1, icc2), 4)
)
write.csv(fit, "output/table8_fit_indices.csv", row.names = FALSE)

cat("\n결과 저장 완료:\n")
cat("  output/table8_coefficients.csv\n")
cat("  output/table8_fit_indices.csv\n")

# =============================================================================
#coded by Ahn 2026.10.26
