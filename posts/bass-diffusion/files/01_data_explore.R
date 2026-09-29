# =============================================================================
# 01_data_explore.R
# Bass 확산 모형 학습 — 1단계: 데이터 탐색
#
# 목적: 모형을 돌리기 전에 자료의 "모양"을 먼저 본다.
#   - 연간 순인원이 S자 확산의 종 모양(증가→정점→감소)을 보이는가?
#   - 누적 순인원이 S자 곡선의 어디쯤에 와 있는가? (정점을 지났는가?)
#   - 이상 연도(2002 월드컵 특수, 2020~21 코로나)는 어디인가?
#   - 공급(사찰 수) 확대가 참가자 증가와 얼마나 겹치는가?
#   - Bass 위상도(n_t vs N_{t-1}): 2단계 OLS 추정이 실제로 하는 일을 미리 본다.
#
# 입력: templestay_2002to2025.csv  (순인원 2002~2025, 연인원 2002~2021, 사찰 수)
# 출력: out/01_*.png, out/01_summary.csv
# 필요 패키지: 없음 (base R)
# 실행: 작업 폴더를 BassModel 로 두고 source("01_data_explore.R")
# =============================================================================

out_dir <- "out"
dir.create(out_dir, showWarnings = FALSE)

# 한글 글꼴: Windows 는 맑은 고딕, 그 외(Linux/Mac)는 나눔고딕
if (.Platform$OS.type == "windows") {
  windowsFonts(KR = windowsFont("Malgun Gothic")); FONT <- "KR"
} else FONT <- "NanumGothic"
open_png <- function(file, w, h) {
  png(file.path(out_dir, file), width = w, height = h, res = 150)
  par(family = FONT)
}

# -----------------------------------------------------------------------------
# 1. 자료 읽기와 파생 변수
# -----------------------------------------------------------------------------
d <- read.csv("templestay_2002to2025.csv", fileEncoding = "UTF-8-BOM")

d$t        <- d$year - 2001                        # 도입 후 연차 (2002 = 1)
d$n        <- d$uniq_total / 1e4                   # 연간 순인원 (만 명) — Bass 의 n(t)
d$N        <- cumsum(d$n)                          # 누적 순인원 (만 명) — Bass 의 N(t)
d$N_lag    <- c(0, head(d$N, -1))                  # 전년도 말 누적 N(t-1)
d$growth   <- c(NA, diff(d$n) / head(d$n, -1))     # 전년 대비 증감률
d$per_temple <- d$uniq_total / d$temples           # 사찰당 순인원 (명)
d$for_share  <- d$uniq_for / d$uniq_total          # 외국인 비중
d$pday_ratio <- d$pday_total / d$uniq_total        # 연인원 / 순인원 (평균 체류 일수의 대리 지표)

covid <- d$year %in% 2020:2021
wc02  <- d$year == 2002

# 콘솔 요약표
show <- data.frame(
  연도 = d$year, 사찰 = d$temples,
  `순인원(만)` = round(d$n, 1), `누적(만)` = round(d$N, 1),
  증감률 = ifelse(is.na(d$growth), "", sprintf("%+.1f%%", 100 * d$growth)),
  사찰당 = round(d$per_temple), 외국인비중 = sprintf("%.0f%%", 100 * d$for_share),
  `연/순` = ifelse(is.na(d$pday_ratio), "", sprintf("%.2f", d$pday_ratio)),
  check.names = FALSE)
cat("\n[1] 연도별 요약\n"); print(show, row.names = FALSE)
write.csv(d, file.path(out_dir, "01_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 몇 가지 핵심 수치
i_max <- which.max(d$n)
cat(sprintf("\n[2] 연간 순인원 최대: %d년 %.1f만 명 (마지막 관측 연도 = %d)\n",
            d$year[i_max], d$n[i_max], max(d$year)))
cat(sprintf("    누적 순인원 2025년 말: %.1f만 명\n", tail(d$N, 1)))
pre <- d$year %in% 2003:2019
cat(sprintf("    2003→2019 연평균 증가율: 순인원 %.1f%%, 사찰 수 %.1f%%\n",
            100 * ((d$n[d$year == 2019] / d$n[d$year == 2003])^(1/16) - 1),
            100 * ((d$temples[d$year == 2019] / d$temples[d$year == 2003])^(1/16) - 1)))
cat(sprintf("    연도별 순인원과 사찰 수의 상관 (2003~2019): r = %.2f\n",
            cor(d$n[pre], d$temples[pre])))

# -----------------------------------------------------------------------------
# 2. 그림 1 — 네 가지 기본 모양 (2x2)
# -----------------------------------------------------------------------------
C_DOM <- "#2f6db5"; C_FOR <- "#d9822b"; C_GREY <- "grey45"
shade_covid <- function() {
  u <- par("usr"); rect(2019.5, u[3], 2021.5, u[4], col = adjustcolor("red", .08), border = NA)
}

open_png("01_fig1_overview.png", 1800, 1300)
par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 4))

# (a) 연간 순인원: 내국인 + 외국인 누적 막대
m <- rbind(d$uniq_dom, d$uniq_for) / 1e4
bp <- barplot(m, names.arg = d$year, col = c(C_DOM, C_FOR), border = NA, las = 2,
              cex.names = .7, ylab = "만 명", main = "(a) 연간 순인원 n(t)")
legend("topleft", c("내국인", "외국인"), fill = c(C_DOM, C_FOR), bty = "n", border = NA)
u <- par("usr"); rect(bp[19] - .6, u[3], bp[20] + .6, u[4], col = adjustcolor("red", .08), border = NA)
text(mean(bp[19:20]), max(colSums(m)) * .97, "코로나", cex = .8, col = "red")

# (b) 누적 순인원
plot(d$year, d$N, type = "b", pch = 16, col = C_DOM, xlab = "", ylab = "만 명",
     main = "(b) 누적 순인원 N(t)")
shade_covid(); grid(col = "grey90")

# (c) 공급: 사찰 수와 사찰당 참가자
plot(d$year, d$temples, type = "b", pch = 16, col = C_GREY, xlab = "", ylab = "사찰 수",
     main = "(c) 공급: 사찰 수(회색) · 사찰당 순인원(주황)")
shade_covid()
par(new = TRUE)
plot(d$year, d$per_temple, type = "b", pch = 17, col = C_FOR, axes = FALSE, xlab = "", ylab = "")
axis(4, col.axis = C_FOR); mtext("사찰당 순인원(명)", side = 4, line = 2.5, col = C_FOR, cex = .8)

# (d) 외국인 비중과 연/순 비율
plot(d$year, 100 * d$for_share, type = "b", pch = 16, col = C_FOR, xlab = "", ylab = "외국인 비중(%)",
     main = "(d) 외국인 비중(주황) · 연인원/순인원(회색)")
shade_covid()
par(new = TRUE)
plot(d$year, d$pday_ratio, type = "b", pch = 17, col = C_GREY, axes = FALSE, xlab = "", ylab = "",
     xlim = range(d$year))
axis(4, col.axis = C_GREY); mtext("연인원/순인원", side = 4, line = 2.5, col = C_GREY, cex = .8)
dev.off()

# -----------------------------------------------------------------------------
# 3. 그림 2 — Bass 위상도: n(t) 대 N(t-1)
#    Bass 이산형: n_t = p*m + (q-p)*N_{t-1} - (q/m)*N_{t-1}^2
#    → 모형이 맞다면 점들이 "뒤집힌 U자(포물선)" 위에 놓인다.
#      - 포물선의 꼭짓점 = 정점, 오른쪽 끝이 가로축과 만나는 곳 = m
#      - 점들이 아직 오르막에만 있으면 꼭짓점·m 을 자료가 직접 보여주지 못한다.
# -----------------------------------------------------------------------------
open_png("01_fig2_phase.png", 1300, 1000)
par(mar = c(4.5, 4.5, 3, 1))
col_pt <- ifelse(covid, "red", ifelse(wc02, C_FOR, C_DOM))
plot(d$N_lag, d$n, pch = 16, col = col_pt, cex = 1.2,
     xlim = c(0, max(d$N_lag) * 1.05), ylim = c(0, max(d$n) * 1.1),
     xlab = "전년도 말 누적 순인원 N(t-1) (만 명)", ylab = "당해 연도 순인원 n(t) (만 명)",
     main = "Bass 위상도: n(t) vs N(t-1)")
lines(d$N_lag, d$n, col = "grey80")
text(d$N_lag, d$n, sprintf("%02d", d$year %% 100), cex = .65, col = col_pt,
     pos = ifelse(d$year == 2003, 2, 3))   # 2003·2004 라벨 겹침 방지
# 참고용: 코로나 2년을 뺀 점들에 2차식을 그대로 맞춘 곡선 (= 2단계 OLS 의 예고편)
ok <- !covid
fit_q <- lm(n ~ N_lag + I(N_lag^2), data = d[ok, ])
xx <- seq(0, max(d$N_lag) * 1.05, length.out = 200)
lines(xx, predict(fit_q, data.frame(N_lag = xx)), lty = 2, col = C_GREY)
legend("topleft", bty = "n", cex = .8,
       legend = c("일반 연도", "2002 (월드컵 특수)", "2020-21 (코로나)", "2차식 적합(코로나 제외, 참고)"),
       col = c(C_DOM, C_FOR, "red", C_GREY), pch = c(16, 16, 16, NA), lty = c(NA, NA, NA, 2))
dev.off()
cat("\n[3] 위상도 2차식(코로나 제외) 계수: ",
    paste(names(coef(fit_q)), signif(coef(fit_q), 3), sep = "=", collapse = ", "), "\n")

cat(sprintf("\n그래프와 요약표를 '%s' 폴더에 저장했습니다.\n", out_dir))
