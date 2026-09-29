# =============================================================================
# 05_rolling_origin.R
# Bass 확산 모형 학습 — 5단계: 관측 시점을 늘려 가며 다시 추정 (rolling origin)
#
# 질문
#   (1) 자료가 한 해씩 쌓일 때 m 과 "정점 연도" 추정은 어떻게 움직이는가?
#       → m 이 수렴하는가, 아니면 누적값을 뒤쫓아 계속 커지는가?
#       → 모형이 말하는 정점이 관측 끝 연도보다 앞(= 이미 지났다)에 놓이는가?
#   참고: 코로나 연도를 추정에서 빼므로 T=2020, 2021 의 추정은 T=2019 와 같다
#   (2) 여러 예측 시작 시점에서 평가하면 Bass 와 Naive/Drift 중 누가 나은가?
#       (4단계의 검증 2회 → 여기서는 최대 16회. 운의 영향을 줄인다)
#
# 설정
#   관측 끝 연도(origin) T = 2008, 2009, ..., 2025  (각 T 마다 2002~T 자료로 추정)
#   코로나 2020~21 은 추정에서 제외 (3단계 '제외' 설정; 누적값에는 포함)
#   예측 평가: 1년 앞, 3년 앞 연간 순인원. 대상 연도가 2020~21 이면 평가에서 제외
#
# 필요 파일: bass_functions.R, templestay_2002to2025.csv
# 출력: out/05_*.png, out/05_rolling_estimates.csv, out/05_rolling_errors.csv
# =============================================================================

source("bass_functions.R")

d <- load_templestay(start_year = 2002)
d$covid <- as.integer(d$year %in% 2020:2021)
d$N_lag <- c(0, head(d$N, -1))

ols_forecast <- function(fit, N_start, h) {
  cf <- coef(fit)[1:3]; out <- numeric(h); N <- N_start
  for (i in 1:h) { n <- cf[1] + cf[2] * N + cf[3] * N^2; out[i] <- n; N <- N + n }
  out
}

origins <- 2008:2025
H <- 3                                                   # 최대 예측 시계
est_rows <- list(); err_rows <- list()

for (T in origins) {
  tr  <- d[d$year <= T & d$covid == 0, ]
  N_T <- d$N[d$year == T]
  fits <- lapply(BASS_METHODS, function(f) f(tr))

  # (1) 추정치 기록
  for (k in names(fits)) {
    r <- fits[[k]]
    pk <- if (is.na(r$m)) c(t_star = NA) else bass_peak(r$p, r$q, r$m)
    est_rows[[length(est_rows) + 1]] <- data.frame(
      origin = T, method = k, p = r$p, q = r$q, m = r$m,
      peak_year = 2001 + pk["t_star"], N_T = N_T, m_over_N = r$m / N_T, row.names = NULL)
  }

  # (2) 1~3년 앞 예측
  fy <- (T + 1):(T + H); tf <- fy - 2001
  preds <- list(
    OLS     = if (is.na(fits$OLS$m)) rep(NA, H) else ols_forecast(fits$OLS$fit, N_T, H),
    NLS_cum = with(fits$NLS_cum, if (is.na(m)) rep(NA, H) else bass_n(tf, p, q, m)),
    NLS_inc = with(fits$NLS_inc, if (is.na(m)) rep(NA, H) else bass_n(tf, p, q, m)),
    Naive   = rep(d$n[d$year == T], H),
    Drift   = d$n[d$year == T] + (1:H) * mean(diff(tail(d$n[d$year <= T], 6))))
  for (h in c(1, 3)) {
    y <- T + h
    if (y > max(d$year) || y %in% 2020:2021) next
    act <- d$n[d$year == y]
    for (k in names(preds))
      err_rows[[length(err_rows) + 1]] <- data.frame(
        origin = T, horizon = h, target = y, model = k,
        pred = preds[[k]][h], actual = act, ape = abs(preds[[k]][h] - act) / act * 100)
  }
}
est <- do.call(rbind, est_rows); err <- do.call(rbind, err_rows)
write.csv(est, file.path("out", "05_rolling_estimates.csv"), row.names = FALSE)
write.csv(err, file.path("out", "05_rolling_errors.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 결과 출력
# -----------------------------------------------------------------------------
cat("\n[1] 관측 끝 연도별 m 추정 (만 명) — 괄호: m / 그 시점 실제 누적\n")
wide <- reshape(est[, c("origin", "method", "m", "m_over_N")], idvar = "origin",
                timevar = "method", direction = "wide")
show <- data.frame(origin = wide$origin, 누적 = round(est$N_T[est$method == "OLS"], 1))
for (k in names(BASS_METHODS))
  show[[k]] <- sprintf("%6.0f (%.1f배)", wide[[paste0("m.", k)]], wide[[paste0("m_over_N.", k)]])
print(show, row.names = FALSE)

cat("\n[2] 관측 끝 연도별 '모형이 말하는 정점 연도'\n")
pk <- reshape(est[, c("origin", "method", "peak_year")], idvar = "origin", timevar = "method", direction = "wide")
names(pk) <- sub("peak_year.", "", names(pk)); pk[-1] <- round(pk[-1], 1)
print(pk, row.names = FALSE)

cat("\n[3] Rolling origin 예측 오차 (MAPE %, 괄호: 평가 횟수)\n")
agg <- aggregate(ape ~ model + horizon, err, function(v) c(mape = mean(v), n = length(v)))
agg <- data.frame(model = agg$model, horizon = agg$horizon,
                  MAPE = round(agg$ape[, "mape"], 1), n = agg$ape[, "n"])
tab <- reshape(agg, idvar = "model", timevar = "horizon", direction = "wide")
tab <- tab[match(c("OLS", "NLS_cum", "NLS_inc", "Naive", "Drift"), tab$model), ]
names(tab) <- c("model", "1년앞 MAPE", "1년앞 n", "3년앞 MAPE", "3년앞 n")
print(tab, row.names = FALSE)

# 각 평가에서 1등을 몇 번 했나 (1년 앞)
win <- do.call(rbind, lapply(split(err[err$horizon == 1, ], err$origin[err$horizon == 1]),
                             function(x) x[which.min(x$ape), c("origin", "model")]))
cat("\n[4] 1년 앞 예측에서 가장 정확했던 모형 (횟수)\n"); print(table(win$model))

# -----------------------------------------------------------------------------
# 그림 1 — m 추정치의 궤적 vs 실제 누적
# -----------------------------------------------------------------------------
open_png("05_fig1_m_path.png", 1900, 800)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
plot(NA, xlim = range(origins), ylim = c(0, max(est$m, na.rm = TRUE) * 1.05),
     xlab = "관측 끝 연도 T", ylab = "만 명", main = "추정된 m 의 궤적 (자료를 T 년까지 사용)")
u <- par("usr"); rect(2019.5, u[3], 2021.5, u[4], col = adjustcolor("red", .08), border = NA)
for (k in names(BASS_METHODS)) with(est[est$method == k, ], lines(origin, m, type = "b", pch = 16, col = METHOD_COL[k], lwd = 2))
with(est[est$method == "OLS", ], lines(origin, N_T, lwd = 2, lty = 2))
legend("topleft", bty = "n", cex = .8, legend = c(names(BASS_METHODS), "그 시점 실제 누적 N(T)"),
       col = c(METHOD_COL, "black"), lty = c(1, 1, 1, 2), pch = c(16, 16, 16, NA), lwd = 2)

plot(NA, xlim = range(origins), ylim = c(1, max(est$m_over_N, na.rm = TRUE) * 1.05), log = "y",
     xlab = "관측 끝 연도 T", ylab = "m / N(T)  (로그 눈금)", main = "m 은 누적의 몇 배? (1 = 이미 포화)")
u <- par("usr"); rect(2019.5, u[3], 2021.5, u[4], col = adjustcolor("red", .08), border = NA)
for (k in names(BASS_METHODS)) with(est[est$method == k, ], lines(origin, m_over_N, type = "b", pch = 16, col = METHOD_COL[k], lwd = 2))
abline(h = 1, lty = 3)
dev.off()

# -----------------------------------------------------------------------------
# 그림 2 — 정점 연도 추정 vs 관측 끝 연도 (대각선 = "지금이 정점")
# -----------------------------------------------------------------------------
open_png("05_fig2_peak_path.png", 1100, 900)
par(mar = c(4.5, 4.5, 3, 1))
yl <- range(c(est$peak_year, origins), na.rm = TRUE)
plot(NA, xlim = range(origins), ylim = yl, xlab = "관측 끝 연도 T", ylab = "모형이 말하는 정점 연도",
     main = "모형은 거의 늘 '정점은 이미 지났다'고 말한다")
abline(0, 1, lty = 2, col = "grey40")
text(2021, 2021.6, "정점 = 관측 끝 연도\n(선 아래 = 이미 지났다고 판단)", pos = 2, cex = .75, col = "grey40")
for (k in names(BASS_METHODS)) with(est[est$method == k, ], lines(origin, peak_year, type = "b", pch = 16, col = METHOD_COL[k], lwd = 2))
legend("topleft", bty = "n", cex = .8, legend = names(BASS_METHODS), col = METHOD_COL, lwd = 2, pch = 16)
dev.off()

# -----------------------------------------------------------------------------
# 그림 3 — 대표 시점(2010·2013·2016·2019·2025)의 증분 NLS 곡선
# -----------------------------------------------------------------------------
sel <- c(2010, 2013, 2016, 2019, 2025)
scol <- c("#fdae61", "#f46d43", "#d73027", "#74add1", "#313695")
yy <- 2002:2050
open_png("05_fig3_curves.png", 1500, 800)
par(mar = c(4, 4.5, 3, 1))
plot(d$year, d$n, pch = 16, col = ifelse(d$covid == 1, "red", "black"), xlim = range(yy), ylim = c(0, 40),
     xlab = "", ylab = "만 명", main = "증분 NLS: 관측 끝 연도별 적합 곡선 (세로 점선 = 관측 끝)")
for (i in seq_along(sel)) {
  r <- est[est$method == "NLS_inc" & est$origin == sel[i], ]
  lines(yy, bass_n(yy - 2001, r$p, r$q, r$m), col = scol[i], lwd = 2)
  abline(v = sel[i] + .5, col = scol[i], lty = 3)
}
legend("topright", bty = "n", cex = .8, legend = paste0("~", sel, "  m=", round(est$m[est$method == "NLS_inc" & est$origin %in% sel])),
       col = scol, lwd = 2)
dev.off()

cat("\n그래프와 결과표를 'out' 폴더에 저장했습니다.\n")
