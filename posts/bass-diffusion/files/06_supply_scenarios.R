# =============================================================================
# 06_supply_scenarios.R
# Bass 확산 모형 학습 — 6단계: 공급(사찰 수)을 반영한 확장 모형과 정책 시나리오
#
# 정책 질문: "사찰 수(공급)를 그대로 두면 수요는 어디서 멈추고,
#             공급을 늘리면 그 천장이 얼마나 올라가는가?"
#
# 모형 (모두 순인원, 코로나 2020~21 은 추정에서 제외)
#   [M0] 고정 m Bass (이산형 NLS)
#        n_t = [p + q N_{t-1}/m] (m - N_{t-1})
#   [M1] 동적 m Bass: 공급이 "천장"을 올린다 (Mahajan & Peterson 1978 계열)
#        n_t = [p + q N_{t-1}/m_t] (m_t - N_{t-1}),   m_t = k × 사찰수_t
#        k = 사찰 1곳이 떠받치는 잠재 참가자 규모(만 명)
#   [M2] Generalized Bass: 공급이 "속도"를 올린다 (Bass, Krishnan & Jain 1994)
#        n_t = m [F(τ_t) - F(τ_{t-1})],  τ_t = Σ x_s,  x_s = 1 + β Δln(사찰수_s)
#        비교 기준은 같은 곡선형인 [M2-0] 증분 NLS (β = 0 인 경우)
#   * M0·M1 은 실제 전년 누적 N_{t-1} 을 쓰는 "한 해씩" 형태, M2 는 곡선형이라
#     적합도(RSS)는 같은 형태끼리 비교한다: M0 vs M1, M2-0 vs M2
#
# 시나리오 (2026~2035 사찰 수)
#   동결 : 2025년 158곳 유지
#   추세 : 2015~2025 연평균 증가 폭만큼 매년 증가
#   확대 : 매년 +6곳 (정책 지원 가정)
#
# 민감도 분석 (5절): 공급 확대의 "질"이 떨어지면?
#   α     : 신규 사찰의 효율 (기존 사찰 대비 잠재시장 기여 비율)  m_t = k × [158 + α × (사찰_t - 158)]
#   q 배수 : 구전 효과 변화 (신규 사찰 만족도 하락 → 부정적 입소문이 전체 q 를 낮춤)
#            신규 사찰이 늘어나는 만큼 점진적으로 하락해 2035년에 해당 배수에 도달한다고 가정
#   비교 기준은 '동결'(α, q 변화 없음)
#
# 필요 파일: bass_functions.R, templestay_2002to2025.csv
# 출력: out/06_*.png, out/06_models.csv, out/06_scenarios.csv, out/06_sensitivity.csv
# =============================================================================

source("bass_functions.R")
set.seed(2026)

d <- load_templestay(start_year = 2002)
d$S     <- d$temples
d$N_lag <- c(0, head(d$N, -1))
d$covid <- d$year %in% 2020:2021
d$dlS   <- c(0, diff(log(d$S)))
e <- d[!d$covid, ]                                   # 추정용 자료

# -----------------------------------------------------------------------------
# 1. 모형 추정
# -----------------------------------------------------------------------------
multi_nls <- function(formula, data, starts, lower, upper) {
  best <- NULL
  for (st in starts) {
    f <- tryCatch(nls(formula, data, start = st, algorithm = "port", lower = lower, upper = upper,
                      control = nls.control(maxiter = 1000)), error = function(e) NULL)
    if (!is.null(f) && (is.null(best) || deviance(f) < deviance(best))) best <- f
  }
  best
}
fit_M0 <- function(data) multi_nls(n ~ (p + q * N_lag / m) * (m - N_lag), data,
  starts = lapply(c(1.5, 2, 3), function(k) list(p = .01, q = .15, m = k * max(data$N))),
  lower = c(p = 1e-6, q = 1e-5, m = max(data$N_lag) + 1), upper = c(p = 1, q = 3, m = 1e5))
fit_M1 <- function(data) multi_nls(n ~ (p + q * N_lag / (k * S)) * (k * S - N_lag), data,
  starts = lapply(c(1.2, 2, 4), function(a) list(p = .01, q = .15, k = a * max(data$N_lag / data$S) + .1)),
  lower = c(p = 1e-6, q = 1e-5, k = max(data$N_lag / data$S) + 1e-3), upper = c(p = 1, q = 3, k = 1e3))

M0 <- fit_M0(e); M1 <- fit_M1(e)

# M2: 누적 시계 τ 에 공급 증가율을 반영 (nls 대신 optim — τ 가 모든 연도에 걸친 누적이므로)
gbm_path <- function(par, dlS) {                       # par = (p, q, m, beta)
  x <- pmax(1 + par[4] * dlS, 0.01); tau <- cumsum(x)
  Fv <- bass_F(tau, par[1], par[2]); par[3] * (Fv - c(0, head(Fv, -1)))
}
fit_curve <- function(with_beta) {
  obj <- function(par) { if (!with_beta) par[4] <- 0
    sum((d$n[!d$covid] - gbm_path(par, d$dlS)[!d$covid])^2) }
  best <- NULL
  for (b0 in if (with_beta) c(0, .5, 1, 2) else 0) for (m0 in c(500, 900, 1500)) {
    o <- optim(c(.007, .15, m0, b0), obj, method = "L-BFGS-B",
               lower = c(1e-5, 1e-4, max(d$N), -2), upper = c(.5, 2, 1e5, 10))
    if (is.null(best) || o$value < best$value) best <- o }
  par <- best$par; if (!with_beta) par[4] <- 0
  list(par = setNames(par, c("p", "q", "m", "beta")), rss = best$value)
}
M2_0 <- fit_curve(FALSE); M2 <- fit_curve(TRUE)

n_obs <- nrow(e)
aic <- function(rss, k) n_obs * log(rss / n_obs) + 2 * k
models <- data.frame(
  model = c("M0 고정m", "M1 동적m(k×사찰)", "M2-0 증분NLS", "M2 GBM(공급→속도)"),
  p = c(coef(M0)["p"], coef(M1)["p"], M2_0$par["p"], M2$par["p"]),
  q = c(coef(M0)["q"], coef(M1)["q"], M2_0$par["q"], M2$par["q"]),
  m_or_k = c(coef(M0)["m"], coef(M1)["k"], M2_0$par["m"], M2$par["m"]),
  beta = c(NA, NA, 0, M2$par["beta"]),
  RSS = c(deviance(M0), deviance(M1), M2_0$rss, M2$rss),
  AIC = c(aic(deviance(M0), 3), aic(deviance(M1), 3), aic(M2_0$rss, 3), aic(M2$rss, 4)), row.names = NULL)
write.csv(models, file.path("out", "06_models.csv"), row.names = FALSE)

cat("\n[1] 모형 비교 (순인원, 코로나 제외; M1 의 m_or_k 는 k = 사찰 1곳당 잠재 참가자, 만 명)\n")
print(transform(models, p = round(p, 4), q = round(q, 3), m_or_k = round(m_or_k, 2),
                beta = round(beta, 3), RSS = round(RSS, 1), AIC = round(AIC, 1)), row.names = FALSE)
cat("\n   M1 표준오차:\n"); print(round(summary(M1)$coefficients[, 1:2], 4))
k_hat <- coef(M1)["k"]
cat(sprintf("\n   M1 해석: 2025년 사찰 158곳 → 잠재시장 m = %.2f × 158 = %.0f만 명 (2025 누적 %.0f만 = %.0f%%)\n",
            k_hat, k_hat * 158, tail(d$N, 1), 100 * tail(d$N, 1) / (k_hat * 158)))

# -----------------------------------------------------------------------------
# 2. 검증: 1년 앞 rolling 예측 (M1 은 다음 해 실제 사찰 수를 안다고 가정한 "조건부 예측")
# -----------------------------------------------------------------------------
step_M0 <- function(cf, N) (cf["p"] + cf["q"] * N / cf["m"]) * (cf["m"] - N)
step_M1 <- function(cf, N, S) (cf["p"] + cf["q"] * N / (cf["k"] * S)) * (cf["k"] * S - N)
roll <- list()
for (T in 2010:2024) {
  y <- T + 1; if (y %in% 2020:2021) next
  tr <- e[e$year <= T, ]; a0 <- fit_M0(tr); a1 <- fit_M1(tr)
  N_T <- d$N[d$year == T]; act <- d$n[d$year == y]
  roll[[length(roll) + 1]] <- data.frame(origin = T, target = y, actual = act,
    M0 = if (is.null(a0)) NA else step_M0(coef(a0), N_T),
    M1 = if (is.null(a1)) NA else step_M1(coef(a1), N_T, d$S[d$year == y]),
    Naive = d$n[d$year == T])
}
roll <- do.call(rbind, roll)
ape <- function(p) round(mean(abs(p - roll$actual) / roll$actual, na.rm = TRUE) * 100, 1)
cat(sprintf("\n[2] 1년 앞 rolling 예측 MAPE (%d회): M0 고정m %.1f%% | M1 동적m %.1f%% | Naive %.1f%%\n",
            nrow(roll), ape(roll$M0), ape(roll$M1), ape(roll$Naive)))

# -----------------------------------------------------------------------------
# 3. 시나리오 예측 (M1) + 부트스트랩 불확실성
# -----------------------------------------------------------------------------
fy <- 2026:2035; h <- length(fy)
trend <- (d$S[d$year == 2025] - d$S[d$year == 2015]) / 10
scen <- list(
  "동결"   = rep(158, h),
  "추세"   = 158 + trend * (1:h),
  "확대"   = 158 + 6 * (1:h))
project_M1 <- function(cf, Spath) {
  N <- tail(d$N, 1); out <- numeric(h)
  for (i in 1:h) { n <- max(step_M1(cf, N, Spath[i]), 0); out[i] <- n; N <- N + n }
  out
}
project_M0 <- function(cf) {
  N <- tail(d$N, 1); out <- numeric(h)
  for (i in 1:h) { n <- max(step_M0(cf, N), 0); out[i] <- n; N <- N + n }
  out
}
point <- sapply(scen, function(S) project_M1(coef(M1), S))
ref_M0 <- project_M0(coef(M0))

# 잔차 부트스트랩: 적합값 + 재표집 잔차로 자료를 다시 만들어 M1 재추정 → 시나리오 재예측
B <- 300; fit_e <- fitted(M1); res_e <- resid(M1)
boot <- array(NA, c(B, h, length(scen)))
for (b in 1:B) {
  eb <- e; eb$n <- fit_e + sample(res_e, replace = TRUE)
  fb <- fit_M1(eb); if (is.null(fb)) next
  for (s in seq_along(scen)) boot[b, , s] <- project_M1(coef(fb), scen[[s]])
}
lo <- apply(boot, c(2, 3), quantile, .1, na.rm = TRUE)
hi <- apply(boot, c(2, 3), quantile, .9, na.rm = TRUE)

tab <- do.call(rbind, lapply(seq_along(scen), function(s) data.frame(
  scenario = names(scen)[s], year = fy, temples = round(scen[[s]], 1),
  n = point[, s], n_lo80 = lo[, s], n_hi80 = hi[, s],
  N = tail(d$N, 1) + cumsum(point[, s]))))
write.csv(tab, file.path("out", "06_scenarios.csv"), row.names = FALSE)

cat(sprintf("\n[3] 시나리오 (M1 동적m). 추세 = 2015~2025 평균 +%.1f곳/년. 괄호 = 80%% 부트스트랩 구간\n", trend))
for (s in names(scen)) {
  x <- tab[tab$scenario == s, ]
  cat(sprintf("   %-4s 2035 사찰 %5.1f곳 | 연간 순인원 2030 %.1f만 (%.1f~%.1f), 2035 %.1f만 (%.1f~%.1f) | 누적 2035 %.0f만\n",
              s, x$temples[h], x$n[5], x$n_lo80[5], x$n_hi80[5], x$n[h], x$n_lo80[h], x$n_hi80[h], x$N[h]))
}
cat(sprintf("   참고: M0 고정m 이면 연간 순인원 2030 %.1f만, 2035 %.1f만 (공급과 무관하게 감소)\n", ref_M0[5], ref_M0[h]))

cat("\n[4] 공급 확대의 한계 효과 (동결 대비)\n")
for (s in c("추세", "확대")) {
  dS  <- scen[[s]][h] - 158
  dN  <- sum(point[, s]) - sum(point[, "동결"])
  dn  <- point[h, s] - point[h, "동결"]
  cat(sprintf("   %s: 2035년까지 사찰 +%.0f곳 → 10년 누적 +%.1f만 명, 2035년 연간 +%.1f만 명 (사찰 1곳당 연 %.0f명)\n",
              s, dS, dN, dn, 1e4 * dn / dS))
}
cat(sprintf("   비교: 2025년 실제 사찰당 순인원 = %.0f명\n", d$uniq_total[d$year == 2025] / 158))

# -----------------------------------------------------------------------------
# 4. 그림
# -----------------------------------------------------------------------------
SC_COL <- c("동결" = "#8a8a8a", "추세" = "#2f6db5", "확대" = "#d9822b")

# 그림 1: 적합 비교 (M0 vs M1, 한 해씩 적합값)
open_png("06_fig1_fit.png", 1900, 800)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 4))
f0 <- step_M0(coef(M0), d$N_lag); f1 <- step_M1(coef(M1), d$N_lag, d$S)
plot(d$year, d$n, pch = 16, col = ifelse(d$covid, "red", "black"), ylim = c(0, 40),
     xlab = "", ylab = "만 명", main = "연간 순인원: 고정 m(M0) vs 동적 m(M1) 적합")
lines(d$year, f0, col = "#8a8a8a", lwd = 2); lines(d$year, f1, col = "#2f6db5", lwd = 2)
legend("topleft", bty = "n", cex = .8, legend = c("실제", "실제(코로나, 추정 제외)",
       sprintf("M0 고정 m (RSS %.0f)", deviance(M0)), sprintf("M1 동적 m (RSS %.0f)", deviance(M1))),
       col = c("black", "red", "#8a8a8a", "#2f6db5"), pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1), lwd = 2)

plot(d$year, d$N, pch = 16, ylim = c(0, max(k_hat * d$S, coef(M0)["m"]) * 1.05), xlab = "",
     ylab = "만 명", main = "누적 순인원과 잠재시장: 고정 m vs 동적 m = k × 사찰수")
lines(d$year, k_hat * d$S, col = "#2f6db5", lwd = 2, type = "s")
abline(h = coef(M0)["m"], col = "#8a8a8a", lwd = 2, lty = 2)
legend("topleft", bty = "n", cex = .8, legend = c("실제 누적", "M1 동적 m (k×사찰수)", "M0 고정 m"),
       col = c("black", "#2f6db5", "#8a8a8a"), pch = c(16, NA, NA), lty = c(NA, 1, 2), lwd = 2)
dev.off()

# 그림 2: 시나리오
open_png("06_fig2_scenarios.png", 1900, 800)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
plot(d$year, d$n, pch = 16, col = ifelse(d$covid, "red", "black"), xlim = c(2002, 2035),
     ylim = c(0, max(hi, na.rm = TRUE) * 1.05), xlab = "", ylab = "만 명",
     main = "연간 순인원: 사찰 수 시나리오별 예측 (음영 = 80% 구간)")
for (s in names(scen)) {
  i <- match(s, names(scen))
  polygon(c(fy, rev(fy)), c(lo[, i], rev(hi[, i])), col = adjustcolor(SC_COL[s], .18), border = NA)
  lines(c(2025, fy), c(tail(d$n, 1), point[, s]), col = SC_COL[s], lwd = 2)
}
lines(c(2025, fy), c(tail(d$n, 1), ref_M0), col = "black", lty = 3, lwd = 2)
abline(v = 2025.5, lty = 3, col = "grey60")
legend("topleft", bty = "n", cex = .8, legend = c(sprintf("동결 (158곳)"), sprintf("추세 (+%.1f곳/년)", trend),
       "확대 (+6곳/년)", "참고: 고정 m 모형"), col = c(SC_COL, "black"), lty = c(1, 1, 1, 3), lwd = 2)

plot(d$year, d$S, type = "b", pch = 16, xlim = c(2002, 2035), ylim = c(0, max(unlist(scen)) * 1.05),
     xlab = "", ylab = "사찰 수", main = "시나리오별 사찰 수 경로")
for (s in names(scen)) lines(c(2025, fy), c(158, scen[[s]]), col = SC_COL[s], lwd = 2)
dev.off()

# -----------------------------------------------------------------------------
# 5. 민감도 분석: 신규 사찰 효율(α) × 구전 효과 변화(q 배수)
# -----------------------------------------------------------------------------
# q 하락은 신규 사찰이 쌓이는 만큼 점진적으로 나타난다고 가정:
#   q 배수_t = 1 - (1 - qmult) × (신규 사찰_t / 2035년 신규 사찰)  → 2035년에 qmult 도달
project_M1_sens <- function(cf, Spath, alpha = 1, qmult = 1) {
  N <- tail(d$N, 1); out <- numeric(h)
  new_share <- if (Spath[h] > 158) (Spath - 158) / (Spath[h] - 158) else rep(0, h)
  for (i in 1:h) {
    m  <- cf["k"] * (158 + alpha * (Spath[i] - 158))
    qm <- 1 - (1 - qmult) * new_share[i]
    n <- max((cf["p"] + cf["q"] * qm * N / m) * (m - N), 0)
    out[i] <- n; N <- N + n
  }
  out
}
alphas <- c(1, .75, .5, .25); qmults <- c(1, .9, .8)
base_freeze <- project_M1_sens(coef(M1), scen[["동결"]])
sens <- do.call(rbind, lapply(c("추세", "확대"), function(s) do.call(rbind, lapply(qmults, function(qm)
  do.call(rbind, lapply(alphas, function(a) {
    x <- project_M1_sens(coef(M1), scen[[s]], a, qm)
    dS <- scen[[s]][h] - 158
    data.frame(scenario = s, alpha = a, q_mult = qm, n_2035 = x[h],
               gain_cum10 = sum(x) - sum(base_freeze),
               gain_per_temple_2035 = 1e4 * (x[h] - base_freeze[h]) / dS)
  }))))))
write.csv(sens, file.path("out", "06_sensitivity.csv"), row.names = FALSE)

cat(sprintf("\n[5] 민감도: 확대 시나리오 (동결 대비, 동결 2035 연간 = %.1f만)\n", base_freeze[h]))
sx <- sens[sens$scenario == "확대", ]
cat("   10년 누적 효과(만 명) — 행: q 배수, 열: 신규 사찰 효율 α\n")
g <- xtabs(gain_cum10 ~ q_mult + alpha, sx); g <- g[order(-as.numeric(rownames(g))), order(-as.numeric(colnames(g)))]
print(round(g, 1))
cat("   2035년 연간 순인원(만 명)\n")
g2 <- xtabs(n_2035 ~ q_mult + alpha, sx); g2 <- g2[order(-as.numeric(rownames(g2))), order(-as.numeric(colnames(g2)))]
print(round(g2, 1))
neg <- sx[sx$gain_cum10 < 0, ]
if (nrow(neg)) cat(sprintf("   ※ %d개 조합에서 확대가 동결보다 오히려 손해 (구전 효과 하락이 공급 효과를 상쇄)\n", nrow(neg)))

# 그림 3: 민감도
open_png("06_fig3_sensitivity.png", 1900, 800)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
QCOL <- c("1" = "#2f6db5", "0.9" = "#d9822b", "0.8" = "#c0392b")
plot(NA, xlim = c(1, .25), ylim = range(c(sx$gain_cum10, 0)), xlab = "신규 사찰 효율 α (기존 = 1)",
     ylab = "10년 누적 효과 (만 명, 동결 대비)", main = "확대 시나리오의 효과: α × 구전 효과")
abline(h = 0, lty = 3)
for (qm in qmults) with(sx[sx$q_mult == qm, ], lines(alpha, gain_cum10, type = "b", pch = 16, lwd = 2, col = QCOL[as.character(qm)]))
legend("topright", bty = "n", cex = .8, title = "q 배수", legend = c("1.0 (변화 없음)", "0.9", "0.8"), col = QCOL, lwd = 2, pch = 16)

plot(d$year, d$n, pch = 16, col = ifelse(d$covid, "red", "black"), xlim = c(2002, 2035), ylim = c(0, 45),
     xlab = "", ylab = "만 명", main = "확대(+6곳/년) 경로: 질 관리에 따라")
combos <- list(c(1, 1), c(.75, 1), c(.5, 1), c(.75, .8))
ccol <- c("#d9822b", "#2f6db5", "#7b3294", "#c0392b")
for (i in seq_along(combos)) lines(c(2025, fy), c(tail(d$n, 1), project_M1_sens(coef(M1), scen[["확대"]], combos[[i]][1], combos[[i]][2])), col = ccol[i], lwd = 2)
lines(c(2025, fy), c(tail(d$n, 1), base_freeze), col = "grey50", lwd = 2, lty = 2)
legend("topleft", bty = "n", cex = .8, col = c(ccol, "grey50"), lwd = 2, lty = c(1, 1, 1, 1, 2),
       legend = c(sprintf("α=%.2f, q×%.1f", sapply(combos, `[`, 1), sapply(combos, `[`, 2)), "동결"))
dev.off()

cat("\n그래프와 결과표를 'out' 폴더에 저장했습니다.\n")
