# =============================================================================
# 02_bass_fit.R
# Bass 확산 모형 학습 — 2단계: 기본 추정 (세 가지 추정법 비교)
#
# 설정 (의도적으로 "손대지 않은" 기본형 — 3단계에서 하나씩 바꿔 본다)
#   - 자료: 순인원 합계 2002~2025 (24개 연도), 단위 만 명
#   - 출발점: 2002년 = t=1
#   - 코로나 2020~21: 그대로 포함
#
# 보는 것
#   (1) 추정법(OLS / 누적 NLS / 증분 NLS)에 따라 p, q, m 이 얼마나 다른가
#   (2) 표준오차: 각 모수가 얼마나 "확실하게" 정해졌는가
#   (3) 모형이 말하는 정점 시점·규모, 2025년까지 m 의 몇 %에 도달했는가
#   (4) 적합 곡선과 잔차: 어디서 체계적으로 어긋나는가
#
# 필요 파일: bass_functions.R, templestay_2002to2025.csv
# 출력: out/02_*.png, out/02_estimates.csv
# =============================================================================

source("bass_functions.R")

d  <- load_templestay(start_year = 2002)
T0 <- 2001                                   # t=1 ↔ 2002년

# -----------------------------------------------------------------------------
# 1. 세 가지 방법으로 추정
# -----------------------------------------------------------------------------
fits <- lapply(BASS_METHODS, function(f) f(d))
est  <- do.call(rbind, lapply(fits, summarise_fit, d = d, t0 = T0))

cat("\n[1] 추정 결과 (순인원 2002~2025, 코로나 포함, 단위: 만 명)\n")
show <- with(est, data.frame(
  방법 = method,
  p = sprintf("%.4f", p), `se(p)` = ifelse(is.na(se_p), "-", sprintf("%.4f", se_p)),
  q = sprintf("%.3f", q), `se(q)` = ifelse(is.na(se_q), "-", sprintf("%.3f", se_q)),
  m = sprintf("%.0f", m), `se(m)` = ifelse(is.na(se_m), "-", sprintf("%.0f", se_m)),
  `q/p` = sprintf("%.1f", q_over_p),
  정점연도 = sprintf("%.1f", peak_year), `정점 연간(만)` = sprintf("%.1f", peak_n),
  `2025 도달률` = sprintf("%.0f%%", 100 * F_last),
  RMSE = sprintf("%.2f", rmse_n), R2 = sprintf("%.2f", r2_n), check.names = FALSE))
print(show, row.names = FALSE)
write.csv(est, file.path("out", "02_estimates.csv"), row.names = FALSE)

# OLS 회귀계수도 따로 보기 (m, p, q 로 역산되기 전의 원래 모습)
cat("\n[2] OLS 회귀계수 (n_t = a + b N_{t-1} + c N_{t-1}^2)\n")
print(round(summary(fits$OLS$fit)$coefficients, 5))

# -----------------------------------------------------------------------------
# 2. 그림 1 — 적합 곡선: 연간(좌) · 누적(우), 2060년까지 연장
# -----------------------------------------------------------------------------
yrs_f <- 2002:2060; tt <- yrs_f - T0
covid <- d$year %in% 2020:2021

open_png("02_fig1_fit.png", 1900, 800)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

ymax <- max(d$n, sapply(fits, function(r) if (is.na(r$m)) 0 else max(bass_n(tt, r$p, r$q, r$m))))
plot(d$year, d$n, pch = 16, col = ifelse(covid, "red", "black"), xlim = range(yrs_f),
     ylim = c(0, ymax * 1.05), xlab = "", ylab = "만 명", main = "연간 순인원: 실제 vs Bass 적합")
for (k in names(fits)) { r <- fits[[k]]; if (!is.na(r$m))
  lines(yrs_f, bass_n(tt, r$p, r$q, r$m), col = METHOD_COL[k], lwd = 2) }
abline(v = 2025.5, lty = 3)
legend("topright", bty = "n", cex = .8, legend = c("실제", "실제(코로나)", names(fits)),
       col = c("black", "red", METHOD_COL), pch = c(16, 16, NA, NA, NA), lty = c(NA, NA, 1, 1, 1), lwd = 2)

mmax <- max(d$N, est$m, na.rm = TRUE)
plot(d$year, d$N, pch = 16, xlim = range(yrs_f), ylim = c(0, mmax * 1.05),
     xlab = "", ylab = "만 명", main = "누적 순인원: 실제 vs Bass 적합 (점선 = 추정 m)")
for (k in names(fits)) { r <- fits[[k]]; if (!is.na(r$m)) {
  lines(yrs_f, bass_N(tt, r$p, r$q, r$m), col = METHOD_COL[k], lwd = 2)
  abline(h = r$m, col = METHOD_COL[k], lty = 2) } }
abline(v = 2025.5, lty = 3)
dev.off()

# -----------------------------------------------------------------------------
# 3. 그림 2 — 연간값 잔차 (실제 - 적합): 체계적인 패턴이 있는가
# -----------------------------------------------------------------------------
open_png("02_fig2_residuals.png", 1500, 700)
par(mar = c(4, 4.5, 3, 1))
res <- sapply(fits, function(r) if (is.na(r$m)) rep(NA, nrow(d)) else d$n - bass_n(d$t, r$p, r$q, r$m))
matplot(d$year, res, type = "b", pch = 16, lty = 1, col = METHOD_COL[colnames(res)],
        xlab = "", ylab = "실제 - 적합 (만 명)", main = "연간 순인원 잔차")
abline(h = 0, col = "grey50")
u <- par("usr"); rect(2019.5, u[3], 2021.5, u[4], col = adjustcolor("red", .08), border = NA)
legend("bottomleft", bty = "n", cex = .8, legend = colnames(res), col = METHOD_COL[colnames(res)], lty = 1, pch = 16)
dev.off()

cat("\n그래프와 추정표를 'out' 폴더에 저장했습니다.\n")
