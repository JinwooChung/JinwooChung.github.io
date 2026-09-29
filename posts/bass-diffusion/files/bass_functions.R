# =============================================================================
# bass_functions.R — Bass 확산 모형 공용 함수 (2단계 이후 모든 스크립트에서 source)
#
#   모형:  누적 채택 비율  F(t) = (1 - e^{-(p+q)t}) / (1 + (q/p) e^{-(p+q)t})
#          누적 채택자     N(t) = m F(t)
#          연간 신규       n(t) = N(t) - N(t-1)      (연 단위 집계 자료용)
#   위험함수: f(t) / (1 - F(t)) = p + q F(t)
#
#   추정 방법 (입력 d 는 t, n, N 열을 가진 data.frame)
#     fit_ols     : Bass(1969) 이산형 OLS  n_t = a + b N_{t-1} + c N_{t-1}^2
#     fit_nls_cum : 누적 NLS               N_t = m F(t)
#     fit_nls_inc : 증분 NLS (Srinivasan & Mason 1986)  n_t = m [F(t) - F(t-1)]
#   각 함수는 list(method, p, q, m, se, fit) 를 돌려준다. (se: 표준오차, OLS 는 NA)
# =============================================================================

bass_F <- function(t, p, q) (1 - exp(-(p + q) * t)) / (1 + (q / p) * exp(-(p + q) * t))
bass_N <- function(t, p, q, m) m * bass_F(t, p, q)
bass_n <- function(t, p, q, m) bass_N(t, p, q, m) - bass_N(t - 1, p, q, m)

# 정점: 연속형 근사 t* = ln(q/p)/(p+q),  n* = m (p+q)^2 / 4q,  정점에서의 누적 비율 F* = (q-p)/2q
bass_peak <- function(p, q, m) {
  if (q <= p) return(c(t_star = 0, n_star = m * p, F_star = 0))  # q<=p 이면 처음부터 감소
  c(t_star = log(q / p) / (p + q), n_star = m * (p + q)^2 / (4 * q), F_star = (q - p) / (2 * q))
}

na_result <- function(method) list(method = method, p = NA, q = NA, m = NA,
                                   se = c(p = NA, q = NA, m = NA), fit = NULL)

# ---- (A) OLS ---------------------------------------------------------------
# n_t = p m + (q - p) N_{t-1} - (q/m) N_{t-1}^2  → 계수 (a, b, c) 로부터
#   m = [-b - sqrt(b^2 - 4ac)] / 2c,  p = a/m,  q = -c m
fit_ols <- function(d, N0 = 0) {
  d$N_lag <- c(N0, head(d$N, -1))
  fit <- lm(n ~ N_lag + I(N_lag^2), data = d)
  cf <- coef(fit); a <- cf[1]; b <- cf[2]; c <- cf[3]
  disc <- b^2 - 4 * a * c
  if (is.na(c) || c >= 0 || disc < 0) {
    r <- na_result("OLS"); r$fit <- fit; r$note <- "c>=0 또는 판별식<0: m 역산 불가"; return(r)
  }
  m <- unname((-b - sqrt(disc)) / (2 * c))
  list(method = "OLS", p = unname(a / m), q = unname(-c * m), m = m,
       se = c(p = NA, q = NA, m = NA), fit = fit)
}

# ---- NLS 공통: 다중 시작점, 경계 제약(port), RSS 최소 해 채택 ---------------
nls_multistart <- function(formula, d, method, N0 = 0) {
  starts <- list()
  o <- fit_ols(d, N0)
  if (!is.na(o$m) && o$p > 0 && o$q > 0 && o$m > max(d$N))
    starts[[1]] <- list(p = o$p, q = o$q, m = o$m)
  for (pq in list(c(.03, .38), c(.01, .15), c(.005, .5), c(.001, .3)))
    for (k in c(1.5, 3, 6))
      starts[[length(starts) + 1]] <- list(p = pq[1], q = pq[2], m = k * max(d$N))
  best <- NULL
  for (st in starts) {
    f <- tryCatch(nls(formula, data = d, start = st, algorithm = "port",
                      lower = c(p = 1e-6, q = 1e-5, m = max(d$N)),
                      upper = c(p = 1, q = 3, m = 1000 * max(d$N)),
                      control = nls.control(maxiter = 1000)),
                  error = function(e) NULL)
    if (!is.null(f) && (is.null(best) || deviance(f) < deviance(best))) best <- f
  }
  if (is.null(best)) return(na_result(method))
  cf <- coef(best)
  se <- tryCatch(summary(best)$coefficients[c("p", "q", "m"), "Std. Error"],
                 error = function(e) c(p = NA, q = NA, m = NA))
  list(method = method, p = unname(cf["p"]), q = unname(cf["q"]), m = unname(cf["m"]),
       se = se, fit = best)
}

# ---- (B) 누적 NLS ----------------------------------------------------------
fit_nls_cum <- function(d, N0 = 0)
  nls_multistart(N ~ m * bass_F(t, p, q), d, "NLS_cum", N0)

# ---- (C) 증분 NLS ----------------------------------------------------------
fit_nls_inc <- function(d, N0 = 0)
  nls_multistart(n ~ m * (bass_F(t, p, q) - bass_F(t - 1, p, q)), d, "NLS_inc", N0)

BASS_METHODS <- list(OLS = fit_ols, NLS_cum = fit_nls_cum, NLS_inc = fit_nls_inc)

# ---- 결과 요약 한 줄 -------------------------------------------------------
#  t0: t=1 에 해당하는 연도 - 1 (예: 2002년이 t=1 이면 t0 = 2001)
summarise_fit <- function(r, d, t0) {
  if (is.na(r$m)) return(data.frame(method = r$method, p = NA, q = NA, m = NA,
                                    se_p = NA, se_q = NA, se_m = NA, q_over_p = NA,
                                    peak_year = NA, peak_n = NA, F_last = NA,
                                    rmse_n = NA, r2_n = NA))
  pk <- bass_peak(r$p, r$q, r$m)
  n_hat <- bass_n(d$t, r$p, r$q, r$m)
  data.frame(method = r$method, p = r$p, q = r$q, m = r$m,
             se_p = r$se["p"], se_q = r$se["q"], se_m = r$se["m"],
             q_over_p = r$q / r$p,
             peak_year = t0 + pk["t_star"], peak_n = pk["n_star"],
             F_last = tail(d$N, 1) / r$m,                     # 마지막 관측 연도까지 도달한 비율
             rmse_n = sqrt(mean((d$n - n_hat)^2)),             # 연간값 기준 적합 오차
             r2_n = 1 - sum((d$n - n_hat)^2) / sum((d$n - mean(d$n))^2),
             row.names = NULL)
}

# ---- 한글 글꼴·PNG 헬퍼 ----------------------------------------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(KR = windowsFont("Malgun Gothic")); FONT <- "KR"
} else FONT <- "NanumGothic"
open_png <- function(file, w, h, dir = "out") {
  dir.create(dir, showWarnings = FALSE)
  png(file.path(dir, file), width = w, height = h, res = 150)
  par(family = FONT)
}
METHOD_COL <- c(OLS = "#8a8a8a", NLS_cum = "#2f6db5", NLS_inc = "#d9822b")

# ---- 자료 읽기 (순인원, 만 명 단위) ----------------------------------------
load_templestay <- function(file = "templestay_2002to2025.csv", start_year = 2002) {
  raw <- read.csv(file, fileEncoding = "UTF-8-BOM")
  d <- raw[raw$year >= start_year, ]
  d$t <- d$year - (start_year - 1)
  d$n <- d$uniq_total / 1e4
  d$N <- cumsum(d$n)
  d
}
