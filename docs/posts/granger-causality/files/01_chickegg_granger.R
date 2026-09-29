# =============================================================================
# 그레인저 인과 검정 실습 1단계: 닭과 달걀 (ChickEgg)          [v3, 2026-09-29]
# -----------------------------------------------------------------------------
# 자료 : lmtest 패키지의 ChickEgg
#        미국 연간 닭 개체 수(chicken)와 달걀 생산량(egg), 1930~1983, T = 54
# 출처 : Thurman & Fisher (1988), "Chickens, Eggs, and Causality, or Which
#        Came First?", American Journal of Agricultural Economics 70(2).
#
# 이 사례의 성격 (먼저 읽어 주세요)
#   이 스크립트는 "표준 절차를 따르면 결론이 강건해 보이지만, 잔차 진단을 하면
#   경고가 나오는 사례"입니다. 1~9절의 표준 절차는 "달걀 → 닭"이 모든 방법에서
#   유의하다고 말하지만, 10절의 진단에서 이분산과 관계의 불안정성이 드러나고,
#   이를 고려하면 결론이 크게 약해집니다(11절). 유명한 결과도 진단을 거치면
#   흔들릴 수 있다는 점, 그리고 잔차 진단·하위표본 분석은 "추가분석"이 아니라
#   본 분석의 일부라는 점을 보여 주는 것이 이 사례의 목적입니다.
#
# 이 스크립트에서 익힐 것
#   1) 정상성 점검 (ADF / KPSS)
#   2) 시차 선택 (정보기준) → 선택된 시차를 이후 모든 절에서 자동으로 사용
#   3) grangertest()로 양방향 검정
#   4) model0(제약) vs model1(비제약)을 직접 만들어 같은 결과가 나오는지 확인
#   5) 시차를 바꾸면 결론이 흔들리는지 (민감도)
#   6) 비정상 수준값을 올바르게 검정하는 Toda–Yamamoto 방식
#   7) VAR 틀의 검정과 부트스트랩 p값 (분포 가정 없이)
#   8) 표준 절차의 요약표 (잠정 결론)
#   9) 잔차 진단 · 이분산 강건 검정 · 하위표본 · 이동창 분석
#  10) 최종 결론과 보고 문장
#
# 실행 방법
#   - 0절(패키지·자료·도우미 함수)과 3절(시차 선택)을 먼저 실행하세요.
#     이후 절들은 3절에서 정한 p_lvl, p_dif를 사용합니다.
#   - 섹션 단위로 한 블록씩 실행하면서 출력과 주석을 비교해 보세요.
#
# v2 변경 사항
#   - 시차를 VARselect() 결과에서 자동으로 가져옴 (기존: 3으로 고정)
#   - 5절을 시차가 바뀌어도 동작하도록 일반화
#   - 8절 부트스트랩 10,000회, 차분 자료 부트스트랩 비교 추가
#   - 9절 요약표 추가
# v3 변경 사항
#   - 10절 추가: 잔차 진단(이분산·자기상관), HC3 강건 검정, 하위표본, 이동창
#   - 11절 추가: 진단을 반영한 최종 결론과 보고 문장
#   - 8절 해석 수정: 부트스트랩 괴리의 원인은 비정상성이 아니라 이분산으로 판단
#   - 9절을 "표준 절차의 잠정 결론"으로 재정의 (기존 보고 문장은 11절로 이동·수정)
# =============================================================================


# ---- 0. 패키지 · 자료 · 도우미 함수 ------------------------------------------
pkgs <- c("lmtest", "tseries", "urca", "vars", "sandwich")
for (pk in pkgs) {
  if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
}
library(lmtest)
library(tseries)
library(sandwich)   # 이분산 강건 표준오차 (10절)
has_vars <- requireNamespace("vars", quietly = TRUE)
if (has_vars) library(vars)

data(ChickEgg)                 # ts 행렬: 열 이름 "chicken", "egg"
ChickEgg_d <- diff(ChickEgg)   # 1차 차분 자료

# 시차 변수를 만드는 도우미: embed()로 [t, t-1, ..., t-p] 열을 만든다
make_lags <- function(x, p, name) {
  m <- embed(as.numeric(x), p + 1)
  colnames(m) <- c(name, paste0(name, "_L", 1:p))
  m
}

# 여러 계열의 시차 변수를 한 데이터프레임으로 (앞의 p개 관측치는 자동으로 빠짐)
lagged_df <- function(tsmat, p) {
  parts <- lapply(colnames(tsmat), function(v) make_lags(tsmat[, v], p, v))
  as.data.frame(do.call(cbind, parts))
}

# Toda–Yamamoto 검정
#   1) 적정 시차 p에 최대 적분차수 dmax를 더한 p + dmax 시차로 추정
#   2) 원래의 p개 X 시차만 0인지 Wald(χ²) 검정 (추가 시차는 검정에서 제외)
toda_yamamoto <- function(tsmat, y, x, p, dmax = 1) {
  k   <- p + dmax
  dd  <- lagged_df(tsmat, k)
  rhs <- c(paste0(y, "_L", 1:k), paste0(x, "_L", 1:k))
  fit <- lm(reformulate(rhs, response = y), data = dd)
  w   <- waldtest(fit, paste0(x, "_L", 1:p), test = "Chisq")
  data.frame(direction = paste(x, "->", y), p = p, dmax = dmax,
             Chisq = round(w$Chisq[2], 3), df = abs(w$Df[2]),
             p_value = round(w$`Pr(>Chisq)`[2], 4))
}

# grangertest()의 p값만 뽑는 도우미
gp <- function(f, k, dat) grangertest(f, order = k, data = dat)$`Pr(>F)`[2]


# ---- 1. 자료 살펴보기 ---------------------------------------------------------
head(ChickEgg)
summary(ChickEgg)

# 두 계열을 각각 표준화해서 한 그림에 그리기 (단위가 크게 달라서)
ts.plot(scale(ChickEgg), col = c("firebrick", "steelblue"), lwd = 2,
        main = "ChickEgg (standardized)", ylab = "z-score")
legend("topright", c("chicken", "egg"),
       col = c("firebrick", "steelblue"), lwd = 2, bty = "n")

# [관찰 포인트]
# - 1940년대 전반(2차대전기)에 달걀 생산이 크게 늘었다가, 이후 닭 개체 수가
#   뒤따라 움직이는 구간이 보입니다. "달걀이 먼저 움직인다"는 느낌이 그래프에서
#   드는지 먼저 눈으로 확인해 두세요. 검정은 그 느낌을 "닭 자신의 과거를
#   넘어서는 정보인가"로 엄밀하게 따지는 일입니다.


# ---- 2. 정상성 점검 -----------------------------------------------------------
# ADF  : 귀무가설 = 단위근 있음(비정상).  p가 작으면 정상.
# KPSS : 귀무가설 = 정상.                p가 작으면 비정상.
# 두 검정을 함께 보는 이유: 소표본에서는 ADF 검정력이 낮아서 한쪽만 보면 헷갈림.

stationarity <- function(x, name) {
  data.frame(
    series = name,
    ADF_p  = round(suppressWarnings(adf.test(x)$p.value), 3),
    KPSS_p = round(suppressWarnings(kpss.test(x)$p.value), 3)
  )
}
rbind(
  stationarity(ChickEgg[, "chicken"],   "chicken (level)"),
  stationarity(ChickEgg[, "egg"],       "egg (level)"),
  stationarity(ChickEgg_d[, "chicken"], "chicken (diff1)"),
  stationarity(ChickEgg_d[, "egg"],     "egg (diff1)")
)

# [해석]
# - 수준값: 두 계열 모두 ADF로는 단위근을 기각 못 하고, KPSS로는 정상성이 기각됨
#   → 수준값은 비정상(I(1))으로 보는 게 안전합니다.
# - 1차 차분: KPSS로는 정상. egg는 ADF도 정상, chicken은 ADF가 애매(T=53의
#   낮은 검정력). 대체로 "1차 차분하면 정상", 즉 최대 적분차수 d_max = 1.
# - 그래서 수준값으로 그대로 검정하면 표준 F분포 가정이 흔들립니다.
#   아래에서는 수준값 F검정, 차분값 F검정, Toda–Yamamoto, 부트스트랩을 모두
#   돌려 결론이 같은지 확인합니다.


# ---- 3. 시차 선택 ------------------------------------------------------------
# 연간 자료이므로 최대 시차 6 안에서 정보기준으로 고릅니다.
# 여기서 정한 p_lvl(수준값용), p_dif(차분값용)을 이후 모든 절에서 사용합니다.
crit <- "SC(n)"   # 기본 기준: SC(=BIC). 짧은 시차를 선호해 소표본에 안전한 편.
                  # "AIC(n)", "HQ(n)", "FPE(n)"로 바꿔 볼 수 있음

if (has_vars) {
  sel_lvl <- VARselect(ChickEgg,   lag.max = 6, type = "const")$selection
  sel_dif <- VARselect(ChickEgg_d, lag.max = 6, type = "const")$selection
  print(sel_lvl)
  print(sel_dif)
  p_lvl <- as.integer(sel_lvl[crit])
  p_dif <- as.integer(sel_dif[crit])
} else {
  p_lvl <- 2L; p_dif <- 1L   # vars가 없을 때: 이 자료의 선택 결과를 직접 입력
  message("vars 패키지가 없어 시차를 직접 지정했습니다: p_lvl = 2, p_dif = 1")
}
cat("선택된 시차  수준값 p_lvl =", p_lvl, " / 차분값 p_dif =", p_dif, "\n")

# [해석]
# - 이 자료에서는 네 기준(AIC, HQ, SC, FPE)이 모두 수준값 2, 차분값 1을 고릅니다.
# - 차분하면 시차가 하나 줄어드는 것은 자연스러운 일입니다. 수준값 VAR(2)는
#   "차분값 시차 1 + 수준 간 장기관계 항"(오차수정모형)으로 다시 쓸 수 있습니다.
# - 일반적으로는 AIC가 긴 시차, SC가 짧은 시차를 고르는 경향이 있어 기준끼리
#   엇갈리는 일이 흔합니다. 그럴 때는 6절처럼 여러 시차에서 결론이 유지되는지
#   보는 게 더 중요합니다.
# - 선택된 시차는 "모형이 과거를 몇 기까지 볼지"에 대한 설정이지,
#   "달걀이 몇 년 뒤에 닭에 영향을 준다"는 뜻이 아닙니다.


# ---- 4. grangertest()로 양방향 검정 (수준값) ---------------------------------
# 공식 Y ~ X 는 "X가 Y를 그레인저 인과하는가"를 검정합니다. (방향 주의!)

# (a) 달걀 → 닭 : 닭 자신의 과거에 달걀의 과거를 더하면 닭 예측이 나아지는가
grangertest(chicken ~ egg, order = p_lvl, data = ChickEgg)

# (b) 닭 → 달걀
grangertest(egg ~ chicken, order = p_lvl, data = ChickEgg)

# [해석] (p_lvl = 2 기준)
# - (a) F(2, 47) = 8.82, p ≈ 0.0006 → 기각. "달걀은 닭을 그레인저 인과한다."
# - (b) F(2, 47) = 0.88, p ≈ 0.42   → 기각 못 함.
# - 결론: 한 방향만 성립 → "달걀이 먼저"라는 논문의 유명한 결론.
# - 출력 읽는 법
#     Model 1 = 비제약모형(X 시차 포함), Model 2 = 제약모형
#     Res.Df  = 잔차 자유도 (관측치 52 − 계수 수)
#     Df      = 두 모형의 계수 차이 = 검정하는 X 시차 수
# - (참고) 교과서 예시의 시차 3에서는 (a) p ≈ 0.003, (b) p ≈ 0.62. 6절 참조.
# - 이 F검정은 비정상 수준값에 F분포를 적용한 것이라 엄밀하지 않습니다.
#   타당성은 7절(Toda–Yamamoto)과 8절(부트스트랩)에서 보완합니다.


# ---- 5. model0 / model1을 직접 만들어 같은 검정 재현 -------------------------
# grangertest() 안에서 일어나는 일을 그대로 풀어 쓴 것입니다. (학습용)
d <- lagged_df(ChickEgg, p_lvl)
nrow(d)   # 54 − p_lvl : 두 모형이 똑같은 관측치를 씀

own_lags   <- paste0("chicken_L", 1:p_lvl)   # 닭 자신의 과거
cause_lags <- paste0("egg_L",     1:p_lvl)   # 달걀의 과거

m0 <- lm(reformulate(own_lags,                response = "chicken"), data = d)  # 제약
m1 <- lm(reformulate(c(own_lags, cause_lags), response = "chicken"), data = d)  # 비제약
formula(m0); formula(m1)   # 어떤 식이 만들어졌는지 확인

anova(m0, m1)   # 잔차제곱합(RSS) 비교 F검정 → 4절 (a)와 F, p가 똑같아야 함

c(RSS_model0 = deviance(m0), RSS_model1 = deviance(m1),
  RSS_reduction = 1 - deviance(m1) / deviance(m0))

summary(m1)$coefficients

# [해석] (p_lvl = 2 기준)
# - anova()의 F(8.82)와 p(0.0006)가 grangertest()와 정확히 일치합니다.
#   달걀 시차를 넣자 잔차제곱합이 약 27% 줄었습니다.
# - 그레인저 검정은 egg_L1, egg_L2를 "한꺼번에" 0으로 놓는 결합검정입니다.
# - 계수를 보면 egg_L1 ≈ +89, egg_L2 ≈ −94로 크기가 비슷하고 부호가 반대입니다.
#   즉 닭 개체 수를 예측하는 것은 달걀 "수준"보다 달걀의 "전년 대비 증감"
#   (egg_{t-1} − egg_{t-2})에 가깝다는 뜻입니다. 개별 계수는 서로 상관이 커서
#   불안정하므로 참고로만 보세요.


# ---- 6. 시차 민감도 -----------------------------------------------------------
# 시차를 1~6으로 바꿔 가며, 수준값과 차분값 양쪽에서 p값을 비교합니다.
sens <- do.call(rbind, lapply(1:6, function(k) {
  data.frame(
    lag = k,
    lvl_egg_to_chicken = gp(chicken ~ egg, k, ChickEgg),
    lvl_chicken_to_egg = gp(egg ~ chicken, k, ChickEgg),
    dif_egg_to_chicken = gp(chicken ~ egg, k, ChickEgg_d),
    dif_chicken_to_egg = gp(egg ~ chicken, k, ChickEgg_d),
    selected = trimws(paste(ifelse(k == p_lvl, "lvl", ""),
                            ifelse(k == p_dif, "dif", "")))
  )
}))
sens[, 2:5] <- round(sens[, 2:5], 4)
print(sens)

# [해석]
# - 수준값 lag = 1 에서는 달걀 → 닭이 유의하지 않고(p ≈ 0.28), lag ≥ 2부터 유의.
#   → 시차를 하나만 보고 결론 내리면 위험하다는 전형적인 예.
# - 차분값에서는 lag 1부터 유의. 시차가 길어지면 p가 커지는 구간도 있는데,
#   모수가 늘면서 자유도가 깎여 검정력이 떨어지는 효과입니다.
# - 닭 → 달걀은 어떤 설정에서도 유의하지 않음. 결론의 방향은 안정적입니다.
# - selected 열은 3절에서 선택된 시차를 표시합니다.


# ---- 7. Toda–Yamamoto: 비정상 수준값을 올바르게 검정하기 ------------------------
# 단위근·공적분 여부와 관계없이 Wald 통계량이 점근적으로 χ²(p)를 따르도록
# 시차를 d_max만큼 더 넣고, 원래 p개 시차만 검정합니다.
ty <- rbind(
  toda_yamamoto(ChickEgg, y = "chicken", x = "egg",     p = p_lvl),  # 본 결과
  toda_yamamoto(ChickEgg, y = "egg",     x = "chicken", p = p_lvl),
  toda_yamamoto(ChickEgg, y = "chicken", x = "egg",     p = 3),      # 민감도
  toda_yamamoto(ChickEgg, y = "egg",     x = "chicken", p = 3)
)
print(ty)

# [해석]
# - p = 2: 달걀 → 닭 χ²(2) = 10.80, p ≈ 0.0045 / 닭 → 달걀 p ≈ 0.63
# - p = 3: 달걀 → 닭 χ²(3) = 11.87, p ≈ 0.008  / 닭 → 달걀 p ≈ 0.76
# - 수준값 F검정보다 p가 조금 큰 것은 추가 시차 때문에 관측치·자유도를 더 쓰기
#   때문입니다. 약간의 검정력을 내주고 검정의 타당성을 얻는 셈입니다.
# - 비정상성을 교정해도 결론이 같으므로, 결론은 검정 방법에 좌우되지 않습니다.


# ---- 8. VAR 틀의 검정과 부트스트랩 p값 ----------------------------------------
# (1) VAR는 두 식을 한 시스템으로 추정합니다. 검정은 여전히 방향별로 합니다.
#     2변수에서는 4절과 같은 F값이 나오고 자유도 계산만 달라집니다.
# (2) 부트스트랩은 "인과가 없는 가상 자료"를 반복 생성해 F의 분포를 직접
#     만들고, 실제 F 이상이 나오는 비율로 p값을 구합니다. 분포 가정이 없습니다.
# 주의: 변수가 3개 이상이면 causality(cause = "X")는 "X → 나머지 변수 전체"의
#       결합검정이 됩니다. 특정 Y 하나에 대한 조건부 검정은 5절 방식으로 하세요.
#
# 부트스트랩 반복 수: p값의 모의실험 오차 ≈ sqrt(p(1−p)/B).
#   p ≈ 0.04, B = 10,000이면 ±0.002. p가 유의수준에 가까울수록 B를 늘리세요.
#   시드는 재현성을 위해 고정하고, 결과를 좋게 하려고 바꾸지 않습니다.
#   10,000회 × 2번이라 PC에 따라 몇 분 걸릴 수 있습니다.
B <- 10000

if (has_vars) {
  # (a) 수준값 VAR (시차 p_lvl)
  v_lvl <- VAR(ChickEgg, p = p_lvl, type = "const")
  print(causality(v_lvl, cause = "egg")$Granger)
  print(causality(v_lvl, cause = "chicken")$Granger)
  set.seed(1)
  boot_lvl <- causality(v_lvl, cause = "egg", boot = TRUE, boot.runs = B)$Granger
  print(boot_lvl)

  # (b) 차분값 VAR (시차 p_dif): 정상 자료에서도 두 p값이 벌어지는지 비교
  v_dif <- VAR(ChickEgg_d, p = p_dif, type = "const")
  print(causality(v_dif, cause = "egg")$Granger)
  set.seed(1)
  boot_dif <- causality(v_dif, cause = "egg", boot = TRUE, boot.runs = B)$Granger
  print(boot_dif)
}

# [해석]
# - 실제 결과 (B = 10,000, seed 1)
#     수준값 VAR(2): 점근 p ≈ 0.0003 → 부트스트랩 p = 0.030
#     차분값 VAR(1): 점근 p ≈ 0.0017 → 부트스트랩 p = 0.055  (5% 수준에서 비유의)
#   (참고: 교과서 시차 3의 수준값 VAR에서는 0.0018 → 0.036)
# - 점근 p값의 기준값 (Python statsmodels로 교차 확인)
#     수준값 VAR(2): 달걀 → 닭 F(2, 94) = 8.82 / 닭 → 달걀 F(2, 94) = 0.88, p ≈ 0.42
#     차분값 VAR(1): 달걀 → 닭 F(1, 98) = 10.37
# - 처음에는 "괴리의 원인은 수준값의 비정상성이고, 정상인 차분 자료에서는 두 p값이
#   가까울 것"이라고 예상했습니다. 그런데 차분 자료에서도 20배 이상 벌어졌으므로
#   이 가설은 틀렸습니다.
# - vars의 부트스트랩은 이분산을 허용하는 wild bootstrap 방식으로 알려져 있습니다
#   (?causality 도움말에서 확인해 보세요).
#   점근 F검정은 잔차 분산이 일정하다고 가정하므로, 이 차이는 잔차의 이분산을
#   의심하게 합니다. → 10절에서 직접 확인합니다.


# ---- 9. 표준 절차의 요약표 (잠정 결론) ----------------------------------------
summary_tbl <- data.frame(
  method = c("F test, level", "F test, diff", "Toda-Yamamoto (dmax = 1)"),
  lag    = c(p_lvl, p_dif, p_lvl),
  egg_to_chicken = round(c(gp(chicken ~ egg, p_lvl, ChickEgg),
                           gp(chicken ~ egg, p_dif, ChickEgg_d),
                           ty$p_value[1]), 4),
  chicken_to_egg = round(c(gp(egg ~ chicken, p_lvl, ChickEgg),
                           gp(egg ~ chicken, p_dif, ChickEgg_d),
                           ty$p_value[2]), 4)
)
if (has_vars) {
  summary_tbl <- rbind(summary_tbl, data.frame(
    method = paste0(c("VAR bootstrap, level", "VAR bootstrap, diff"), " (B = ", B, ")"),
    lag    = c(p_lvl, p_dif),
    egg_to_chicken = round(c(boot_lvl$p.value, boot_dif$p.value), 4),
    chicken_to_egg = NA))
}
print(summary_tbl, row.names = FALSE)

# [잠정 결론]
# - 여기까지만 보면 "달걀 → 닭은 모든 방법에서 유의(부트스트랩 차분만 경계선),
#   닭 → 달걀은 모두 비유의"이므로 결론이 강건해 보입니다.
# - 그러나 4~8절의 검정은 부트스트랩을 제외하면 모두 "잔차 분산이 일정하다"는
#   같은 가정 위에 있습니다. 같은 가정을 공유하는 방법끼리 결과가 일치하는 것은
#   그 가정이 깨졌을 때의 강건성을 보장하지 않습니다. → 10절에서 가정을 점검합니다.


# ---- 10. 잔차 진단과 강건성 재점검 --------------------------------------------
# (1) 잔차 진단 : 비제약모형(model1)의 잔차가 회귀 가정을 만족하는가
#     - Breusch–Pagan  : 귀무가설 = 등분산.     p가 작으면 이분산.
#     - Breusch–Godfrey: 귀무가설 = 자기상관 없음. p가 작으면 시차가 부족하다는 신호.
fit_pair <- function(tsmat, y, x, p) {
  d   <- lagged_df(tsmat, p)
  own <- paste0(y, "_L", 1:p); cz <- paste0(x, "_L", 1:p)
  list(m0 = lm(reformulate(own, response = y), data = d),
       m1 = lm(reformulate(c(own, cz), response = y), data = d),
       year = as.numeric(time(tsmat))[-(1:p)])
}

diag_row <- function(label, tsmat, y, x, p) {
  f <- fit_pair(tsmat, y, x, p)
  data.frame(
    model      = label, lag = p,
    BP_p       = round(bptest(f$m1)$p.value, 4),               # 이분산
    BG_p       = round(bgtest(f$m1, order = 2)$p.value, 4),    # 잔차 자기상관
    classic_p  = round(waldtest(f$m0, f$m1)$`Pr(>F)`[2], 4),   # 일반 F검정
    HC3_p      = round(waldtest(f$m0, f$m1,                    # 이분산 강건 F검정
                          vcov = vcovHC(f$m1, type = "HC3"))$`Pr(>F)`[2], 4)
  )
}

diag_tbl <- rbind(
  diag_row("level: egg -> chicken", ChickEgg,   "chicken", "egg",     p_lvl),
  diag_row("level: chicken -> egg", ChickEgg,   "egg",     "chicken", p_lvl),
  diag_row("diff : egg -> chicken", ChickEgg_d, "chicken", "egg",     p_dif),
  diag_row("diff : chicken -> egg", ChickEgg_d, "egg",     "chicken", p_dif)
)
print(diag_tbl, row.names = FALSE)

# 잔차를 연도별로 그려 보기 (수준값, 달걀 → 닭 식)
f_lvl <- fit_pair(ChickEgg, "chicken", "egg", p_lvl)
plot(f_lvl$year, resid(f_lvl$m1), type = "h", lwd = 2, col = "grey40",
     xlab = "year", ylab = "residual", main = "Residuals: chicken ~ lags (level)")
abline(h = 0); abline(v = 1949.5, lty = 2, col = "firebrick")

tapply(resid(f_lvl$m1), f_lvl$year < 1950, sd)   # TRUE = 1930~49, FALSE = 1950~83

# [해석]
# - BP_p: 달걀 → 닭 식에서 p < 0.001 (수준·차분 모두) → 강한 이분산.
#   잔차 표준편차가 1930~49년 약 30,000, 1950~83년 약 15,000으로 2배 차이.
#   대공황·2차대전기의 큰 충격 때문입니다. 그래프에서 점선 왼쪽이 크게 튑니다.
# - BG_p: 모두 0.3 이상 → 잔차 자기상관 없음. 시차 선택(2, 1)은 적절했습니다.
# - HC3_p: 이분산을 반영하면 달걀 → 닭의 p가
#     수준값 0.0006 → 0.040,  차분값 0.0023 → 0.040
#   으로 커집니다. 8절 부트스트랩 결과(0.030, 0.055)와 비슷한 수준입니다.
#   → 부트스트랩 괴리의 주원인은 이분산이었다는 해석이 뒷받침됩니다.


# (2) 이분산 강건 Toda–Yamamoto : 7절과 같은 검정에 HC3 공분산만 바꿔 끼움
toda_yamamoto_hc <- function(tsmat, y, x, p, dmax = 1) {
  k   <- p + dmax
  dd  <- lagged_df(tsmat, k)
  fit <- lm(reformulate(c(paste0(y, "_L", 1:k), paste0(x, "_L", 1:k)), response = y),
            data = dd)
  w   <- waldtest(fit, paste0(x, "_L", 1:p), test = "Chisq",
                  vcov = vcovHC(fit, type = "HC3"))
  data.frame(direction = paste(x, "->", y), p = p, dmax = dmax,
             Chisq = round(w$Chisq[2], 3), df = abs(w$Df[2]),
             p_value = round(w$`Pr(>Chisq)`[2], 4))
}
ty_hc <- rbind(
  toda_yamamoto_hc(ChickEgg, "chicken", "egg",     p_lvl),
  toda_yamamoto_hc(ChickEgg, "egg",     "chicken", p_lvl)
)
print(ty_hc)

# [해석]
# - 달걀 → 닭: χ²(2) = 10.80, p = 0.0045 (7절)  →  χ²(2) = 4.53, p ≈ 0.10 (HC3)
#   이분산을 고려하면 본 결과로 삼았던 검정이 5% 수준에서 유의하지 않게 됩니다.
# - 닭 → 달걀: p ≈ 0.76으로 여전히 비유의.


# (3) 하위표본 : 관계가 기간 전체에 걸쳐 있는가, 특정 시기에만 있는가
sub_row <- function(from, to) {
  ce <- window(ChickEgg, start = from, end = to)
  data.frame(period = paste0(from, "-", to), T = nrow(ce),
             lvl_egg_to_chicken = round(gp(chicken ~ egg, p_lvl, ce), 4),
             lvl_chicken_to_egg = round(gp(egg ~ chicken, p_lvl, ce), 4),
             dif_egg_to_chicken = round(gp(chicken ~ egg, p_dif, diff(ce)), 4))
}
print(rbind(sub_row(1930, 1983), sub_row(1930, 1949), sub_row(1950, 1983)),
      row.names = FALSE)

# [해석]
# - 1930~49 (T = 20): 달걀 → 닭 p ≈ 0.01 (수준·차분). 관측치가 적은데도 유의.
# - 1950~83 (T = 34): 달걀 → 닭 p ≈ 0.93 (수준), 0.80 (차분). 흔적이 없음.
#   T = 34면 검정력이 낮긴 하지만, p가 0.9 근처라 "놓쳤다"고 보기도 어렵습니다.
# - 전체 기간의 유의한 결과는 1950년 이전 구간이 만들어 낸 것으로 보입니다.


# (4) 이동창(rolling window) : 창을 한 해씩 옮기며 p값이 어떻게 변하는지
W <- 25   # 창 길이(년)
starts <- start(ChickEgg)[1]:(end(ChickEgg)[1] - W + 1)
roll <- data.frame(
  start = starts, end = starts + W - 1,
  p = sapply(starts, function(s) {
    ce <- window(ChickEgg, start = s, end = s + W - 1)
    gp(chicken ~ egg, p_dif, diff(ce))           # 차분값, 달걀 → 닭
  })
)
plot(roll$end, roll$p, type = "b", pch = 19, log = "y",
     xlab = "window end year", ylab = "p-value (log scale)",
     main = paste0("Rolling ", W, "-year Granger test: egg -> chicken (diff)"))
abline(h = 0.05, lty = 2, col = "firebrick")
print(transform(roll, p = round(p, 3)))

# [해석]
# - 1930~1954 창부터 1941~1965 창까지는 p ≤ 0.005로 매우 유의합니다.
# - 창의 시작이 1943년 이후로 넘어가 1940년대 전반이 빠지는 순간 p가 0.05 위로
#   올라가고(1944·1945 시작 창은 0.05~0.09로 경계선), 이후로는 0.2~0.9 사이입니다.
# - 즉 "달걀 → 닭" 신호는 2차대전기 달걀 증산과 그 뒤의 닭 개체 수 증가라는
#   특정 에피소드에 집중되어 있습니다.


# ---- 11. 최종 결론과 보고 문장 -------------------------------------------------
# [최종 결론]
# - 표준 절차(4~7절)에서는 "달걀 → 닭"이 강하게 유의했지만(p < 0.01),
#   잔차에 강한 이분산이 있었고, 이를 반영하면 증거는 경계선 수준으로 약해집니다
#   (HC3 F검정 p ≈ 0.04, 부트스트랩 p ≈ 0.03~0.055, HC3 Toda–Yamamoto p ≈ 0.10).
# - 관계는 1950년 이전, 특히 1940년대 전반에 집중되어 있고 이후에는 나타나지
#   않습니다. 기간 전체에 걸친 안정적인 관계로 보기 어렵습니다.
# - 닭 → 달걀은 어떤 방법·기간에서도 유의하지 않았습니다(1930~49 수준값 p ≈ 0.05
#   하나를 제외). 이것이 "닭이 달걀에 영향이 없다"는 뜻은 아닙니다. 연 단위
#   자료로는 같은 해 안의 동시점 효과를 잡을 수 없습니다.
#
# [보고 문장 예시]
#   두 계열이 모두 I(1)이어서 Toda–Yamamoto 방식으로 검정했다(시차 2, d_max = 1).
#   일반 Wald 검정에서는 달걀 생산량이 닭 개체 수를 그레인저 인과하는 것으로
#   나타났으나(χ²(2) = 10.80, p = 0.005), 잔차에 강한 이분산이 확인되어
#   (Breusch–Pagan p < 0.001) 이분산 강건 표준오차로 다시 검정하자 유의성이
#   약해졌다(χ²(2) = 4.53, p = 0.10). 하위표본과 이동창 분석에서도 이 관계는
#   1950년 이전, 특히 1940년대 전반에만 나타났다. 따라서 달걀 → 닭의 그레인저
#   인과는 기간 전체에 걸친 안정적 관계라기보다 2차대전기 에피소드에 의존한
#   결과로 해석하는 것이 타당하다. 반대 방향의 인과는 확인되지 않았다.
#
# [다른 프로젝트에 적용할 때 체크리스트]
#   □ 정상성 점검 → 적분차수 d_max 결정
#   □ 정보기준으로 시차 선택 + 주변 시차 민감도
#   □ 양방향 검정 (Toda–Yamamoto 또는 차분 검정)
#   □ 잔차 진단: 이분산(BP), 자기상관(BG)
#   □ 이분산이 있으면 HC 강건 검정 또는 wild bootstrap
#   □ 하위표본·이동창으로 관계의 시간적 안정성 확인
#     (코로나, 금융위기, 제도 변경 등 알려진 구조적 단절 전후 비교)
#   □ 결론에 "통제변수", "기간", "정보적 선행" 단서 붙이기


# =============================================================================
# 스스로 해 볼 과제
#   1) 3절의 crit를 "AIC(n)"로 바꾸고 전체를 다시 실행해 보기 (이 자료에선 동일)
#   2) 10절 (4)의 창 길이 W를 20, 30으로 바꿔 보기. 결론이 창 길이에 민감한가?
#   3) 동시점 관계 보기: lm(chicken ~ egg, data = ChickEgg)의 R²는 높은데
#      그레인저 결과와 어떻게 다른 이야기를 하는지 비교
#   4) "달걀 → 닭"이 성립해도 닭이 달걀을 낳지 않는다는 뜻이 아닙니다.
#      닭 → 달걀 효과는 같은 해 안에 일어나는 동시점 효과라 연 단위 시차로는
#      잡히지 않을 수 있습니다. 측정된 것은 '정보적 선행'일 뿐입니다.
# =============================================================================
