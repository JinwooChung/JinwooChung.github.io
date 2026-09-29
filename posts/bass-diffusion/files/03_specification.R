# =============================================================================
# 03_specification.R
# Bass 확산 모형 학습 — 3단계: 분석 설정(specification) 비교
#
# 질문: "분석자가 내리는 선택"이 p, q, m 과 정점 판단을 얼마나 바꾸는가?
#
#   결정 1. 출발점:  2002년 = t1  vs  2003년 = t1  (2002는 월드컵 특수로 운영된 해)
#   결정 2. 코로나 2020~21 처리
#       (a) 포함   : 아무것도 하지 않음 (2단계와 같음)
#       (b) 제외   : 두 해를 추정에서 뺌. 단, 누적값 N 에는 그 해 참가자가 그대로 들어 있음
#       (c) 더미   : 코로나 연도 효과를 모수로 따로 추정
#            - OLS    : n_t = a + b N_{t-1} + c N_{t-1}^2 + δ·covid     (δ = 줄어든 인원, 만 명)
#            - 증분NLS: n_t = m[F(t)-F(t-1)] × (1 - δ·covid)            (δ = 줄어든 비율)
#            - 누적NLS: 누적식에는 더미를 자연스럽게 넣기 어려워 생략
#
#   → 2 × 3 × 3방법 = 최대 18개 추정
#
# 필요 파일: bass_functions.R, templestay_2002to2025.csv
# 출력: out/03_*.png, out/03_specs.csv
# =============================================================================

source("bass_functions.R")

run_spec <- function(start_year, covid_mode) {
  d  <- load_templestay(start_year = start_year)
  t0 <- start_year - 1
  d$covid <- as.integer(d$year %in% 2020:2021)
  d$N_lag <- c(0, head(d$N, -1))              # 실제 전년 누적 (행을 빼도 유지되도록 미리 계산)

  d_est <- if (covid_mode == "exclude") d[d$covid == 0, ] else d

  fits <- switch(covid_mode,
    include = list(OLS = fit_ols(d_est), NLS_cum = fit_nls_cum(d_est), NLS_inc = fit_nls_inc(d_est)),
    exclude = list(OLS = fit_ols(d_est), NLS_cum = fit_nls_cum(d_est), NLS_inc = fit_nls_inc(d_est)),
    dummy   = list(OLS = fit_ols(d_est, covid_dummy = TRUE),
                   NLS_cum = na_result("NLS_cum"),
                   NLS_inc = fit_nls_inc_dummy(d_est)))

  out <- do.call(rbind, lapply(fits, function(r) {
    s <- summarise_fit(r, d, t0)             # 적합도·도달률은 전체 자료 기준으로 계산
    s$delta <- NA
    if (covid_mode == "dummy" && !is.na(r$m)) {
      s$delta <- if (r$method == "OLS") coef(r$fit)["covid"] else r$extra["delta"]
    }
    s
  }))
  cbind(start = start_year, covid = covid_mode, out)
}

specs <- expand.grid(start = c(2002, 2003), covid = c("include", "exclude", "dummy"),
                     stringsAsFactors = FALSE)
res <- do.call(rbind, Map(run_spec, specs$start, specs$covid))
res <- res[!(res$covid == "dummy" & res$method == "NLS_cum"), ]
rownames(res) <- NULL
write.csv(res, file.path("out", "03_specs.csv"), row.names = FALSE)

cat("\n[1] 설정별 추정 결과 (단위: 만 명; δ: OLS=줄어든 인원(만), NLS_inc=줄어든 비율)\n")
show <- with(res, data.frame(
  시작 = start, 코로나 = covid, 방법 = method,
  p = sprintf("%.4f", p), q = sprintf("%.3f", q),
  m = sprintf("%.0f", m), `se(m)` = ifelse(is.na(se_m), "-", sprintf("%.0f", se_m)),
  정점 = sprintf("%.1f", peak_year), `2025도달` = sprintf("%.0f%%", 100 * F_last),
  δ = ifelse(is.na(delta), "", sprintf("%.2f", delta)), R2 = sprintf("%.2f", r2_n),
  check.names = FALSE))
print(show, row.names = FALSE)

cat("\n[2] 방법별 m 범위 (모든 설정을 통틀어)\n")
print(aggregate(m ~ method, res, function(v) round(c(min = min(v), max = max(v)))))

# -----------------------------------------------------------------------------
# 그림 1 — 설정별 m 추정치 (점) ± 1.96 se (선), 방법별 색
# -----------------------------------------------------------------------------
res$label <- sprintf("%d·%s", res$start, c(include = "포함", exclude = "제외", dummy = "더미")[res$covid])
labs <- unique(res$label)
open_png("03_fig1_m_by_spec.png", 1500, 850)
par(mfrow = c(1, 2), mar = c(5, 7, 3, 1))
for (what in c("m", "peak_year")) {
  xr <- if (what == "m") range(c(res$m - 1.96 * ifelse(is.na(res$se_m), 0, res$se_m),
                                 res$m + 1.96 * ifelse(is.na(res$se_m), 0, res$se_m)), na.rm = TRUE)
        else range(res$peak_year, 2025, na.rm = TRUE)
  xr <- if (what == "m") c(max(0, xr[1]), min(xr[2], 3000)) else xr
  plot(NA, xlim = xr, ylim = c(.5, length(labs) + .5), yaxt = "n", ylab = "",
       xlab = if (what == "m") "m (만 명)" else "정점 연도",
       main = if (what == "m") "잠재시장 m (± 1.96 se)" else "모형이 말하는 정점 연도")
  axis(2, at = seq_along(labs), labels = labs, las = 1, cex.axis = .8)
  off <- c(OLS = -.2, NLS_cum = 0, NLS_inc = .2)
  for (i in seq_len(nrow(res))) {
    y <- match(res$label[i], labs) + off[res$method[i]]
    x <- res[[what]][i]; col <- METHOD_COL[res$method[i]]
    if (what == "m" && !is.na(res$se_m[i]))
      segments(max(0, x - 1.96 * res$se_m[i]), y, x + 1.96 * res$se_m[i], y, col = col, lwd = 2)
    points(x, y, pch = 16, col = col, cex = 1.2)
  }
  if (what == "m") {                                   # 2025년 말 실제 누적 = m 의 하한
    abline(v = max(res$m * res$F_last, na.rm = TRUE), lty = 3)
    text(max(res$m * res$F_last, na.rm = TRUE), length(labs) + .45, "2025 누적", pos = 4, cex = .7)
  }
  if (what == "peak_year") abline(v = 2025, lty = 3)
  legend("bottomright", bty = "n", cex = .75, legend = names(METHOD_COL), col = METHOD_COL, pch = 16)
}
dev.off()

# -----------------------------------------------------------------------------
# 그림 2 — 증분 NLS 의 연간 적합 곡선을 설정별로 겹쳐 보기
# -----------------------------------------------------------------------------
d_all <- load_templestay(start_year = 2002)
yrs <- 2002:2060
open_png("03_fig2_curves_nls_inc.png", 1400, 800)
par(mar = c(4, 4.5, 3, 1))
plot(d_all$year, d_all$n, pch = 16, col = ifelse(d_all$year %in% 2020:2021, "red", "black"),
     xlim = range(yrs), ylim = c(0, 40), xlab = "", ylab = "만 명",
     main = "증분 NLS 적합 곡선: 설정별 비교 (더미 모형은 코로나가 없었을 경우의 곡선)")
sub <- res[res$method == "NLS_inc", ]
lcol <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
for (i in seq_len(nrow(sub))) {
  t0 <- sub$start[i] - 1
  yy <- yrs[yrs >= sub$start[i]]
  lines(yy, bass_n(yy - t0, sub$p[i], sub$q[i], sub$m[i]), col = lcol[i], lwd = 2)
}
abline(v = 2025.5, lty = 3)
legend("topright", bty = "n", cex = .8,
       legend = ifelse(sub$covid == "dummy", paste0(sub$label, " (제외와 사실상 겹침)"), sub$label),
       col = lcol[seq_len(nrow(sub))], lwd = 2)
dev.off()

cat("\n그래프와 결과표를 'out' 폴더에 저장했습니다.\n")
