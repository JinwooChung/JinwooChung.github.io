# =============================================================================
# 04_forecast_validate.R
# Bass 확산 모형 학습 — 4단계: 예측 검증(holdout)과 논문 모형 재현
#
# Part A. 표본 밖 예측 검증 (순인원)
#   "표본 내 적합도(R2)가 좋다" ≠ "예측을 잘한다". 일부 연도를 숨겨 두고 맞히는지 본다.
#     H1: 2002~2015 로 추정 → 2016~2019 예측   (코로나 없는 깨끗한 검증)
#     H2: 2002~2019 로 추정 → 2022~2025 예측   (코로나 이후 회복기; 2020~21 은 오차 계산에서 제외)
#   비교 기준(벤치마크): Bass 가 "아무 모형도 안 쓴 것"보다 나은가?
#     Naive : 마지막 관측값이 그대로 유지된다고 예측
#     Drift : 마지막 5년의 평균 증가분만큼 매년 늘어난다고 예측
#   예측 방식
#     OLS      : 마지막 실제 누적 N_T 에서 출발해 n = a + bN + cN^2 를 한 해씩 반복 적용
#     NLS 2종  : 추정된 곡선 n(t) = m[F(t) - F(t-1)] 을 그대로 연장
#
# Part B. 논문 모형 재현 (연인원 2002~2023)
#   자료: 연인원 2002~2021 (20년사) + 2022·2023 (논문 표1 누적값의 차분: 43.0, 54.6 만)
#   목표: 논문의 m=907.5만, p=.01, q=.18 과 표1·표2 예측치를 재현 → 논문이 쓴 추정법을 역추적
#
# 필요 파일: bass_functions.R, templestay_2002to2025.csv
# 출력: out/04_*.png, out/04_holdout.csv, out/04_paper_replication.csv
# =============================================================================

source("bass_functions.R")

# OLS 예측: 출발 누적 N_start 에서 h 년 앞까지 이산형 Bass 식을 반복 적용
ols_forecast <- function(fit, N_start, h) {
  cf <- coef(fit)[1:3]; out <- numeric(h); N <- N_start
  for (i in 1:h) { n <- cf[1] + cf[2] * N + cf[3] * N^2; out[i] <- n; N <- N + n }
  out
}
mape <- function(actual, pred) mean(abs(pred - actual) / actual) * 100

# =============================================================================
# Part A. 표본 밖 예측 검증
# =============================================================================
d <- load_templestay(start_year = 2002)            # 순인원, t=1 ↔ 2002

run_holdout <- function(train_end, test_years) {
  tr <- d[d$year <= train_end, ]
  fyears <- (train_end + 1):max(test_years)          # 예측 구간 (코로나 연도 포함해 연속으로)
  h <- length(fyears); tf <- fyears - 2001
  preds <- list()
  f_ols <- fit_ols(tr)
  preds$OLS     <- ols_forecast(f_ols$fit, tail(tr$N, 1), h)
  f_cum <- fit_nls_cum(tr); preds$NLS_cum <- bass_n(tf, f_cum$p, f_cum$q, f_cum$m)
  f_inc <- fit_nls_inc(tr); preds$NLS_inc <- bass_n(tf, f_inc$p, f_inc$q, f_inc$m)
  last  <- tail(tr$n, 1)
  preds$Naive <- rep(last, h)
  preds$Drift <- last + (1:h) * mean(diff(tail(tr$n, 6)))
  pars <- rbind(OLS = c(f_ols$p, f_ols$q, f_ols$m), NLS_cum = c(f_cum$p, f_cum$q, f_cum$m),
                NLS_inc = c(f_inc$p, f_inc$q, f_inc$m))
  colnames(pars) <- c("p", "q", "m")
  idx <- fyears %in% test_years
  act <- d$n[match(fyears[idx], d$year)]
  tab <- data.frame(holdout = sprintf("~%d → %d-%d", train_end, min(test_years), max(test_years)),
                    model = names(preds),
                    MAPE = sapply(preds, function(p) mape(act, p[idx])),
                    pred_last = sapply(preds, function(p) p[length(p)]),
                    actual_last = tail(act, 1), row.names = NULL)
  list(tab = tab, preds = preds, fyears = fyears, pars = pars, train_end = train_end, test = test_years)
}

H1 <- run_holdout(2015, 2016:2019)
H2 <- run_holdout(2019, 2022:2025)
hold <- rbind(H1$tab, H2$tab)
write.csv(hold, file.path("out", "04_holdout.csv"), row.names = FALSE)

cat("\n[A-1] 학습 구간별 추정 모수 (순인원, 만 명)\n")
cat("  H1 (2002~2015)\n"); print(round(H1$pars, 4))
cat("  H2 (2002~2019)\n"); print(round(H2$pars, 4))
cat("\n[A-2] 예측 오차 (MAPE, 연간 순인원 기준; 작을수록 좋음)\n")
print(transform(hold, MAPE = sprintf("%.1f%%", MAPE), pred_last = round(pred_last, 1),
                actual_last = round(actual_last, 1)), row.names = FALSE)

MODEL_COL <- c(METHOD_COL, Naive = "black", Drift = "#7b3294")
MODEL_LTY <- c(OLS = 1, NLS_cum = 1, NLS_inc = 1, Naive = 3, Drift = 2)
open_png("04_fig1_holdout.png", 1900, 800)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
for (H in list(H1, H2)) {
  plot(d$year, d$n, type = "n", xlim = c(2002, 2025), ylim = c(0, 40), xlab = "", ylab = "만 명",
       main = sprintf("%d년까지로 추정 → %d~%d 예측", H$train_end, min(H$test), max(H$test)))
  u <- par("usr"); rect(min(H$test) - .5, u[3], max(H$test) + .5, u[4], col = "grey93", border = NA)
  if (H$train_end == 2019) rect(2019.5, u[3], 2021.5, u[4], col = adjustcolor("red", .08), border = NA)
  for (k in names(H$preds)) lines(H$fyears, H$preds[[k]], col = MODEL_COL[k], lty = MODEL_LTY[k], lwd = 2)
  tr <- d$year <= H$train_end
  points(d$year[tr], d$n[tr], pch = 16)
  points(d$year[!tr], d$n[!tr], pch = 1, cex = 1.1)
  legend("topleft", bty = "n", cex = .75, ncol = 2,
         legend = c("학습 자료", "검증 자료", names(H$preds)),
         col = c("black", "black", MODEL_COL[names(H$preds)]), pch = c(16, 1, rep(NA, 5)),
         lty = c(NA, NA, MODEL_LTY[names(H$preds)]), lwd = 2)
}
dev.off()

# =============================================================================
# Part B. 논문 모형 재현 (연인원)
# =============================================================================
raw <- read.csv("templestay_2002to2025.csv", fileEncoding = "UTF-8-BOM")
pd  <- data.frame(year = 2002:2023,
                  n = c(raw$pday_total[raw$year <= 2021] / 1e4, 43.0, 54.6))   # 2022·23: 논문 표1 차분
pd$t <- pd$year - 2001; pd$N <- cumsum(pd$n)

# 논문 수치 (표1: 2002~2023 자료 모형의 예측 누적, 표2: 향후 예측 누적)
paper_t1 <- data.frame(year = 2011:2023,
  actual = c(204.9, 240.5, 278.9, 316.2, 357.3, 398.8, 447.5, 499.1, 552.0, 575.8, 601.6, 644.6, 699.2),
  pred   = c(216.2, 254.0, 294.5, 337.5, 382.3, 428.3, 474.5, 520.3, 564.7, 607.0, 646.7, 683.1, 716.1))
paper_t2 <- data.frame(year = c(2024, 2029, 2034, 2039, 2044, 2049),
                       pred = c(745.4, 842.9, 883.9, 899.2, 904.6, 906.5))

# 세 가지 방법으로 추정해 논문 모수와 비교
cat("\n[B-1] 연인원 2002~2023 추정 — 논문: m=907.5, p=.01, q=.18\n")
fits_pd <- lapply(BASS_METHODS, function(f) f(pd))
print(round(t(sapply(fits_pd, function(r) c(p = r$p, q = r$q, m = r$m))), 4))
cat("   (코로나 이전 2002~2019) — 논문: m=954.3, p=.01, q=.18\n")
print(round(t(sapply(BASS_METHODS, function(f) { r <- f(pd[pd$year <= 2019, ]); c(p = r$p, q = r$q, m = r$m) })), 4))

# OLS 로 확인된 경우: 논문의 "예측 누적"은 N=0 에서 출발해 이산형 식을 반복 적용한 경로와 일치하는가?
ols_pd <- fits_pd$OLS
path   <- ols_forecast(ols_pd$fit, 0, length(2002:2049))        # 2002~2049 연간 예측 경로
cum_path <- setNames(cumsum(path), 2002:2049)
rep_t1 <- transform(paper_t1, replicated = round(cum_path[as.character(year)], 1))
rep_t1$diff_pct_paper <- round(100 * (paper_t1$actual - paper_t1$pred) / paper_t1$actual, 2)
rep_t2 <- transform(paper_t2, replicated = round(cum_path[as.character(year)], 1))
cat("\n[B-2] 표1 재현: 논문 예측 누적 vs 재현 (OLS 계수, N=0 에서 반복)\n"); print(rep_t1, row.names = FALSE)
cat("\n[B-3] 표2 재현: 향후 예측 누적\n"); print(rep_t2, row.names = FALSE)
write.csv(rbind(cbind(table = "표1", rep_t1[, c("year", "pred", "replicated")]),
                cbind(table = "표2", rep_t2[, c("year", "pred", "replicated")])),
          file.path("out", "04_paper_replication.csv"), row.names = FALSE)

# 논문 모형이 말한 2024·2025 연간 수요 vs 실제 순인원 추세
n_path <- setNames(path, 2002:2049)
cat(sprintf("\n[B-4] 논문 모형의 연간 연인원 예측: 2023 %.1f만 → 2024 %.1f만 → 2025 %.1f만 (감소)\n",
            n_path["2023"], n_path["2024"], n_path["2025"]))
cat(sprintf("      실제 연인원 2023 = %.1f만. 실제 순인원은 2024 %+.1f%%, 2025 %+.1f%% (역대 최다 경신)\n",
            54.6, 100 * (d$n[d$year == 2024] / d$n[d$year == 2023] - 1),
            100 * (d$n[d$year == 2025] / d$n[d$year == 2024] - 1)))

open_png("04_fig2_paper_replication.png", 1900, 800)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
yy <- 2002:2049
plot(pd$year, pd$N, pch = 16, xlim = range(yy), ylim = c(0, 1000), xlab = "", ylab = "만 명(연인원)",
     main = "누적 연인원: 논문 표1·표2 vs 재현")
lines(yy, cum_path, col = METHOD_COL["OLS"], lwd = 2)
points(paper_t1$year, paper_t1$pred, pch = 4, col = "red", cex = 1.2)
points(paper_t2$year, paper_t2$pred, pch = 4, col = "red", cex = 1.2)
abline(h = ols_pd$m, lty = 2, col = METHOD_COL["OLS"])
legend("bottomright", bty = "n", cex = .8,
       legend = c("실제 누적 연인원", "재현: OLS 이산형 반복 경로", "논문 표1·표2 수치", "재현 m"),
       col = c("black", METHOD_COL["OLS"], "red", METHOD_COL["OLS"]),
       pch = c(16, NA, 4, NA), lty = c(NA, 1, NA, 2), lwd = 2)

plot(pd$year, pd$n, pch = 16, xlim = c(2002, 2040), ylim = c(0, 60), xlab = "", ylab = "만 명(연인원)",
     main = "연간 연인원: 논문 모형의 예측 경로")
lines(yy, n_path, col = METHOD_COL["OLS"], lwd = 2)
abline(v = 2023.5, lty = 3)
text(2024, 54, "논문 모형: 2024년 이후 감소\n실제 순인원:\n2024 +15.9%, 2025 +5.1%", pos = 4, cex = .75)
dev.off()

cat("\n그래프와 결과표를 'out' 폴더에 저장했습니다.\n")
