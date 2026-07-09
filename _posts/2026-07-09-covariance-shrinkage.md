---
title: "공분산 행렬 Shrinkage 추정 완벽 정리: Linear에서 Nonlinear까지"
date: 2026-07-09 20:00:00 +0900
categories: [Statistics, Finance]
tags: [shrinkage, covariance-matrix, random-matrix-theory, portfolio-optimization, ledoit-wolf]
math: true
---

포트폴리오 최적화, PCA, 판별분석 등 다변량 통계의 거의 모든 곳에서 공분산 행렬 추정이 필요합니다. 그런데 변수의 개수(차원)가 표본 크기에 비해 크지 않더라도, 표본공분산행렬은 생각보다 훨씬 부정확한 추정량입니다. 이 글에서는 이 문제를 해결하는 **shrinkage(축소) 추정** 기법을, Ledoit과 Wolf의 두 논문을 바탕으로 정리합니다.

- Ledoit & Wolf, *"Honey, I Shrunk the Sample Covariance Matrix"* (2003) — **Linear shrinkage**의 원리와 실전 포트폴리오 적용
- Ledoit & Wolf, *"Nonlinear Shrinkage Estimation of Large-Dimensional Covariance Matrices"*, Annals of Statistics (2012) — **Nonlinear shrinkage**로의 확장

수학적으로 유도가 필요한 두 가지 핵심 결과(최적 축소 강도 공식, 오라클 공식)는 본문에서는 결과와 의미만 설명하고, 전체 유도 과정은 **부록**에 따로 정리했습니다.

**목차**
1. 서론: 왜 공분산 행렬 추정이 어려운가
2. 표본공분산행렬의 문제
3. Shrinkage의 기본 아이디어
4. 실전 검증: 포트폴리오 최적화에서의 Linear Shrinkage
5. 왜 Nonlinear Shrinkage가 필요한가
6. 회전등변 추정량(Rotation-Equivariant Estimator)과 오라클
7. 큰 차원 점근이론: Marchenko-Pastur 방정식
8. Oracle Nonlinear Shrinkage 공식
9. Bona Fide 추정량: 실제로 계산 가능하게 만들기
10. 몬테카를로 시뮬레이션으로 본 성능
11. Linear vs Nonlinear 요약 비교
12. 마무리
- 부록 A. 최적 선형 축소 강도(δ*) 유도
- 부록 B. 회전등변 오라클 공식 유도
- 참고문헌

## 1. 서론: 왜 공분산 행렬 추정이 어려운가

$p$개의 변수(예: 주식 종목)에 대해 $n$번의 관측치가 있다고 합시다. 가장 자연스러운 공분산 행렬 추정량은 표본공분산행렬

$$
S = \frac{1}{n}\sum_{t=1}^{n}(y_t-\bar{y})(y_t-\bar{y})'
$$

입니다. $S$는 불편추정량(unbiased estimator)이라는 좋은 성질을 갖고 있습니다. 문제는 **차원 $p$가 표본 크기 $n$에 비해 작지 않을 때** 발생합니다. 금융에서는 이런 상황이 오히려 일반적입니다 — 수백 개 종목의 리스크를 추정하려는데, 안정적인 시계열은 몇 년 치(수십~수백 개월)밖에 없는 경우가 흔하니까요.

## 2. 표본공분산행렬의 문제

$p$가 $n$에 비해 작지 않으면 $S$의 개별 원소들은 심한 추정오차를 갖게 됩니다. 특히 중요한 것은 **고유값(eigenvalue)의 왜곡**입니다: $S$의 가장 큰 고유값들은 실제(모집단) 값보다 체계적으로 더 크게, 가장 작은 고유값들은 더 작게 추정되는 경향이 있습니다. 즉 표본공분산행렬의 고유값 스펙트럼은 실제보다 더 넓게 퍼져 보입니다.

이게 왜 문제냐면, 포트폴리오 최적화 같은 절차는 공분산 행렬의 정보를 그대로 신뢰해서 의사결정을 내립니다. 만약 특정 종목 조합의 분산이 (추정오차 때문에) 우연히 아주 작게 나왔다면, 평균-분산 최적화 알고리즘은 바로 그 조합에 가장 큰 베팅을 겁니다. 즉 **가장 신뢰할 수 없는 추정치에 가장 크게 의존하는** 역설적인 상황이 발생합니다. 이를 Michaud(1989)는 "error-maximization"이라 불렀습니다. 결과적으로 매니저의 실현 성과는 실제 종목 선정 능력보다 체계적으로 낮게 나타납니다.

## 3. Shrinkage의 기본 아이디어

Shrinkage의 핵심 발상은 간단합니다: 극단적인 값으로 치우친 추정치를 좀 더 중심 쪽으로 당겨오자는 것입니다. 통계학적으로 이는 전형적인 **편향-분산 트레이드오프(bias-variance tradeoff)** 문제입니다.

비유하자면, 다트를 던지는 두 사람이 있다고 생각해봅시다. 한 명은 평균적으로는 정중앙을 맞히지만(불편) 매번 크게 흩어져 던지는 사람($S$에 해당), 다른 한 명은 항상 비슷한 자리에 몰아 던지지만 정중앙에서 살짝 벗어난 곳을 겨냥하는 사람(구조화된 추정량 $F$에 해당)입니다. 이 둘의 조준을 적절히 평균 내면, 둘 중 어느 한쪽만 쓰는 것보다 더 정확한 결과를 얻을 수 있습니다. Shrinkage는 정확히 이 원리를 공분산 행렬 추정에 적용한 것입니다.

### 3.1 Shrinkage 추정량의 세 가지 구성요소

어떤 shrinkage 추정량이든 세 가지 요소로 구성됩니다.

1. **구조가 없는 추정량**: 불편이지만 분산(추정오차)이 큰 추정량. 여기서는 표본공분산행렬 $S$.
2. **구조화된 추정량(shrinkage target)**: 자유모수가 적어 추정오차는 작지만, 모델이 틀렸을 경우 편향이 큰 추정량. $F$로 표기.
3. **축소 강도(shrinkage intensity)** $\delta \in [0,1]$: 두 추정량을 얼마나 섞을지 결정하는 가중치.

이 세 요소를 결합한 최종 추정량은 다음과 같은 볼록결합(convex combination)입니다.

$$
\hat{\Sigma}_{\text{Shrink}} = \delta F + (1-\delta) S
$$

### 3.2 Shrinkage Target: Constant Correlation Model

좋은 target이 되려면 (1) 자유모수가 적어야 하고 (2) 동시에 데이터의 중요한 특징은 반영해야 합니다. Ledoit & Wolf(2003)가 제안한 target은 **등상관모형(constant correlation model)**입니다: 모든 종목 쌍의 상관계수가 동일하다고 가정하는 것입니다.

구체적으로, 표본상관계수 $r_{ij} = s_{ij}/\sqrt{s_{ii}s_{jj}}$들의 평균을 $\bar r$이라 하면, target 행렬 $F$는

$$
f_{ii} = s_{ii}, \qquad f_{ij} = \bar r \sqrt{s_{ii}s_{jj}} \quad (i \neq j)
$$

로 정의됩니다. 즉 **분산은 표본값을 그대로 쓰고, 상관계수만 하나의 평균값으로 통일**하는 것입니다. 이렇게 하면 추정해야 할 상관계수의 개수가 $p(p-1)/2$개에서 단 1개($\bar r$)로 줄어들어, 추정오차가 극적으로 감소합니다. 물론 실제 상관구조가 종목마다 다르다면 이 모형은 어느 정도 "틀린" 모형이 되지만, 그 대가로 얻는 분산 감소가 훨씬 크다는 것이 shrinkage의 핵심 논리입니다.

### 3.3 최적 Shrinkage 강도 $\delta^*$

$\delta$가 클수록 target 쪽으로 더 많이 끌려가고, 작을수록 표본공분산행렬에 가까워집니다. "최적의" $\delta$는 참값 $\Sigma$와의 거리(Frobenius norm 기준 제곱오차)의 기댓값을 최소화하는 값으로 정의합니다.

$$
\delta^* = \underset{\delta}{\arg\min}\; \mathbb{E}\left[\|\delta F + (1-\delta)S - \Sigma\|^2\right]
$$

이 최적화 문제를 풀면 (표본 크기 $n\to\infty$, 차원 $p$ 고정인 고전적 점근이론 하에서) 다음과 같은 아름다운 결과가 나옵니다.

$$
\delta^* \approx \frac{\kappa}{n}, \qquad \kappa = \frac{\pi - \rho}{\gamma}
$$

여기서 세 양은 각각 다음을 의미합니다.

- $\pi$: $S$의 각 원소가 갖는 **표본변동성(추정오차의 크기)**의 총합. $S$가 얼마나 "시끄러운(noisy)" 추정량인지를 나타냅니다.
- $\gamma$: target $F$가 수렴하는 모집단 값 $\Phi$가 진짜 $\Sigma$와 얼마나 다른지, 즉 **target의 모형 오설정(misspecification) 정도**.
- $\rho$: $F$와 $S$가 공유하는 표본변동성. $F$ 자체도 같은 데이터로 만들어지기 때문에 $S$의 노이즈와 어느 정도 얽혀 있는데, 이를 보정해주는 항입니다.

직관적으로 해석하면: **$S$가 시끄러울수록($\pi$↑) 더 많이 축소하고, target이 잘못됐을수록($\gamma$↑) 덜 축소한다**는 것입니다. 실제 데이터에서는 $\pi,\rho,\gamma$를 각각 일관추정량(consistent estimator)으로 추정한 뒤

$$
\hat\delta^* = \max\left\{0,\ \min\left\{\frac{\hat\kappa}{n},\ 1\right\}\right\}
$$

를 사용합니다 (0과 1 사이를 벗어나지 않도록 자르는 안전장치입니다). 이 공식의 전체 유도 과정은 **부록 A**에서 단계별로 다룹니다.

## 4. 실전 검증: 포트폴리오 최적화에서의 Linear Shrinkage

Ledoit & Wolf(2003)는 이 방법을 실제 미국 주식 데이터(1983~2002년, 월간)에 적용해 검증했습니다. 매달 시가총액 상위 $p$개 종목(=30, 50, 100, 225, 500)으로 벤치마크를 구성하고, 직전 $n=60$개월 데이터로 공분산 행렬을 추정한 뒤, 이를 이용해 평균-분산 최적화로 포트폴리오를 구성합니다. 성과는 연율화된 **정보비율(Information Ratio, IR)** — 벤치마크 대비 초과수익을 그 변동성으로 나눈 값 — 로 측정합니다.

결과를 요약하면 다음과 같습니다 (연율화 300bp 목표수익 기준, 50회 반복 평균).

| 종목 수 $p$ | 표본공분산행렬 IR | Shrinkage(등상관) IR |
|---|---|---|
| 30 | 0.97 | **1.24** |
| 50 | 0.79 | **1.14** |
| 100 | 0.59 | **0.91** |
| 225 | 0.37 | **0.54** |
| 500 | 0.20 | **0.30** |

모든 $p$에서 shrinkage 추정량이 표본공분산행렬을 뚜렷하게 앞섭니다. 또한 shrinkage 추정량은 매달 포트폴리오를 재조정할 때의 회전율(turnover)도 표본공분산행렬보다 낮게 나타나 거래비용 측면에서도 유리했습니다. 흥미로운 점은 $p$가 커질수록(즉 종목 수가 늘어나 $p/n$ 비율이 커질수록) 모든 방법의 절대적인 IR은 낮아지지만, shrinkage의 **상대적** 우위는 오히려 더 뚜렷해진다는 것입니다.

## 5. 왜 Nonlinear Shrinkage가 필요한가

Linear shrinkage는 하나의 $\delta$를 **모든 고유값에 똑같은 비율로** 적용합니다. 예를 들어 $\delta=0.5$라면, 모든 표본 고유값이 전체 평균 쪽으로 정확히 절반씩 이동합니다. 그런데 이것이 정말 최선일까요?

이론적으로, Marchenko-Pastur 방정식(뒤에서 다룹니다)이 시사하는 진짜 최적의 변환은 **비선형(nonlinear)** 입니다. Linear shrinkage는 이 비선형 문제에 대한 1차 근사에 불과합니다. 이 근사가 얼마나 좋은지는 상황에 따라 크게 갈립니다.

- $p/n$ 비율이 크거나(표본 크기가 차원에 비해 상대적으로 작음), 모집단 고유값들이 서로 가까이 모여 있는 경우 → linear shrinkage가 잠재적 개선분의 대부분을 이미 포착합니다.
- 반대로 $p/n$ 비율이 작거나(표본 크기가 차원에 비해 충분히 큼), 모집단 고유값들이 넓게 퍼져 있는 경우 → linear shrinkage는 거의 개선을 이루지 못합니다.

두 번째 경우가 직관적으로 이상하게 느껴질 수 있습니다 — 데이터가 많은데 왜 linear 방법이 잘 못할까요? 이유는, 데이터가 충분히 많으면 모집단 고유값 분포의 **정교한(비선형적인) 형태**까지 정밀하게 복원할 수 있는 정보가 이미 데이터 안에 있는데, linear shrinkage는 그 정보를 하나의 숫자 $\delta$로 뭉개버리기 때문입니다. 반대로 데이터가 부족하면 애초에 모든 것이 흐릿해서, 정교한 비선형 보정을 해봤자 얻는 이득 자체가 크지 않습니다.

이 문제를 해결하려면, 고유값의 **위치에 따라 개별적으로 다른 강도**의 축소를 적용해야 합니다. 이것이 nonlinear shrinkage의 핵심 아이디어입니다.

## 6. 회전등변 추정량(Rotation-Equivariant Estimator)과 오라클

특정 좌표축을 우대할 이유가 없다면, 데이터를 회전시켰을 때 추정량도 같이 회전해야 한다는 요구는 매우 자연스럽습니다. 임의의 직교행렬(회전행렬) $W$에 대해

$$
\hat\Sigma(Y_nW) = W'\hat\Sigma(Y_n)W
$$

를 만족하는 추정량을 **회전등변(rotation-equivariant)** 추정량이라 부릅니다. 이 조건을 만족하는 추정량은 반드시 표본공분산행렬 $S$와 **같은 고유벡터**를 가져야 하고, 오직 고유값들($D_n = \text{diag}(d_1,\dots,d_p)$)만 다르게 바꿀 수 있습니다. 즉

$$
\hat\Sigma = U_n D_n U_n'
$$

의 형태여야 합니다 ($U_n$은 $S$의 고유벡터 행렬). Linear shrinkage도, 앞으로 다룰 nonlinear shrinkage도 모두 이 틀 안에 있는 방법입니다 — 고유벡터는 그대로 두고 고유값만 손을 보는 것이죠.

그렇다면 질문: 만약 진짜 $\Sigma$를 알고 있다면(현실에서는 불가능하지만), 이 클래스 안에서 **가장 좋은** $D_n$은 무엇일까요? 답은 놀랍도록 단순합니다.

$$
d_i^{*} = u_i' \Sigma_n u_i \qquad (i=1,\dots,p)
$$

즉 각 고유값을 그 고유벡터 방향으로 $\Sigma$를 "투영"한 값으로 바꾸는 것이 최적입니다. 이를 **오라클(oracle) 추정량**이라 부릅니다 — 실제로는 계산할 수 없지만(모집단 $\Sigma$를 알아야 하므로), 회전등변 추정량 중 이론적으로 도달 가능한 최선의 성능을 알려주는 기준점 역할을 합니다. 이 결과의 유도는 **부록 B**에서 다룹니다.

## 7. 큰 차원 점근이론: Marchenko-Pastur 방정식

오라클 공식 $d_i^* = u_i'\Sigma u_i$는 $\Sigma$를 모르니 그대로 쓸 수 없습니다. Ledoit & Péché(2011)는 **큰 차원 점근이론(large-dimensional asymptotics)** — 차원 $p$와 표본 크기 $n$이 함께 무한대로 가되 비율 $c = p/n$이 일정한 상수로 수렴하는 점근 체계 — 를 이용해, $d_i^*$를 **관측 가능한 양만으로 근사**하는 방법을 제시했습니다.

이 체계에서는 표본 고유값들의 분포 $F_n$이 어떤 (확정적인) 극한분포 $F$로 수렴한다는 사실이 알려져 있습니다. 그리고 이 극한분포 $F$는 모집단 고유값 분포 $H$와 **Marchenko-Pastur 방정식**이라는 관계식으로 연결되어 있습니다 (Marchenko & Pastur, 1967).

이 방정식을 다루는 핵심 도구가 **Stieltjes 변환**입니다. 어떤 분포 $G$에 대해

$$
m_G(z) \equiv \int \frac{1}{\lambda - z}\, dG(\lambda), \qquad z \in \mathbb{C}^+
$$

로 정의되며, 확률분포를 복소평면 위의 해석함수로 바꿔주는 변환입니다 (푸리에 변환과 비슷한 역할을 한다고 생각하면 됩니다). Marchenko-Pastur 방정식은 $F$의 Stieltjes 변환과 $H$의 관계를 적분방정식 형태로 표현합니다.

> 이 방정식 자체의 완전한 증명은 확률론적 랜덤행렬이론(Random Matrix Theory)의 깊은 결과이며, Marchenko & Pastur(1967), Silverstein & Choi(1995) 등 여러 논문에 걸쳐 확립된 것입니다. 이 글에서는 결과를 어떻게 활용하는지에 집중하고, 이론 자체의 증명은 다루지 않습니다.

## 8. Oracle Nonlinear Shrinkage 공식

Marchenko-Pastur 이론을 이용하면, $d_i^*$를 다음과 같이 근사할 수 있습니다 (Ledoit & Péché, 2011).

$$
d_i^{or} \equiv \frac{\lambda_i}{\left|1 - c - c\lambda_i\, \breve{m}_F(\lambda_i)\right|^2}
$$

여기서 $\lambda_i$는 표본 고유값, $c = p/n$은 concentration ratio, $\breve{m}_F(\lambda)$는 표본 고유값의 극한분포 $F$의 Stieltjes 변환을 실수축 경계로 확장한 값입니다.

이 공식이 중요한 이유는, **더 이상 $\Sigma$(모집단 공분산 행렬) 자체를 직접 요구하지 않는다**는 데 있습니다. 대신 표본 고유값들이 (극한적으로) 어떤 분포를 이루는지에만 의존합니다. 이는 원리적으로 데이터로부터 추정 가능한 양입니다 — 다음 절의 주제입니다.

같은 논리를 정밀도행렬(precision matrix, 공분산행렬의 역행렬) $\Sigma_n^{-1}$에도 적용할 수 있으며, 그 오라클 공식은

$$
a_i^{or} \equiv \lambda_i^{-1}\left(1 - c - 2c\lambda_i\, \text{Re}[\breve{m}_F(\lambda_i)]\right)
$$

입니다. 흥미롭게도 $(S_n^{or})^{-1} \neq P_n^{or}$입니다 — 즉 공분산 행렬을 잘 추정한 뒤 역행렬을 취하는 것과, 정밀도행렬을 처음부터 직접 추정하는 것은 다른 결과를 줍니다. (뒤의 시뮬레이션에서 직접 추정이 더 낫다는 것을 확인합니다.)

## 9. Bona Fide 추정량: 실제로 계산 가능하게 만들기

오라클 공식은 여전히 $\breve{m}_F(\lambda_i)$ — 극한분포 $F$의 정보 — 를 요구합니다. $F$는 모집단 고유값 분포 $H$에 의해 결정되는데, $H$ 역시 우리가 모르는 대상입니다. Ledoit & Wolf(2012) 논문의 핵심 기여는 바로 이 **$H$를 데이터로부터 일관되게(consistently) 추정하는 방법**을 제시한 것입니다.

방법의 개요는 다음과 같습니다.

1. 후보 분포 $\tilde H$들을 (그리드 위의 유한개 기저함수들의 가중합으로) 표현합니다.
2. 각 후보 $\tilde H$를 Marchenko-Pastur 방정식에 대입하면, 이론적으로 예측되는 표본 고유값 분포 $F_{\tilde H}$를 계산할 수 있습니다.
3. 이 이론적 분포 $F_{\tilde H}$가 **실제 관측된** 표본 고유값들의 분포 $F_n$과 가장 가까워지도록 $\tilde H$(정확히는 기저함수 가중치들)를 최적화합니다.

이 최적화 문제는 비볼록(non-convex)이지만, Sequential Linear Programming(SLP)이라는 수치기법으로 빠르고 안정적으로 풀 수 있음이 실험적으로 확인되었습니다. 최적화가 끝나면 추정된 $\breve{m}_F^{*}(\lambda_i)$를 얻고, 이를 오라클 공식에 대입해 최종 **bona fide(실제 사용 가능한)** 추정량을 얻습니다.

$$
\hat{d}_i = \frac{\lambda_i}{\left|1 - \frac{p}{n} - \frac{p}{n}\lambda_i\, \breve{m}_F^{*}(\lambda_i)\right|^2}
$$

이렇게 얻은 $\hat S_n = U_n \hat D_n U_n'$은, 이론적으로 표본 크기가 커질수록 오라클 추정량 $S_n^{or}$에 수렴한다는 것이 증명되어 있습니다(Corollary 5.2). 즉 "모르는 것을 아는 것처럼" 흉내 낼 수 있게 된 것입니다.

## 10. 몬테카를로 시뮬레이션으로 본 성능

Ledoit & Wolf(2012)는 성능을 비교하기 위해 **PRIAL(Percentage Relative Improvement in Average Loss)** 지표를 사용합니다.

$$
\text{PRIAL}(\hat\Sigma_n) = 100 \times \left\{1 - \frac{\mathbb{E}[\|\hat\Sigma_n - S_n^{*}\|^2]}{\mathbb{E}[\|S_n - S_n^{*}\|^2]}\right\}\%
$$

즉 표본공분산행렬 $S_n$ 대비 얼마나 유한표본 오라클 $S_n^*$에 가까워졌는지를 %로 나타낸 것입니다 ($S_n$은 0%, $S_n^*$는 100%). 주요 결과를 요약하면 다음과 같습니다 (모집단 고유값이 20%는 1, 40%는 3, 40%는 10인 대표적 시나리오 기준).

- **차원에 따른 수렴**: $c=p/n=1/3$ 고정 시, $p=30$인 작은 행렬에서도 nonlinear shrinkage는 이미 잠재적 개선분의 **88%**를 달성하며, $p$가 커질수록 오라클에 빠르게 수렴합니다. Linear shrinkage는 65~70% 선에서 정체됩니다.
- **Concentration ratio($c=p/n$) 효과**: nonlinear shrinkage는 $c$가 0.1~0.9 전 구간에서 항상 **90% 이상**을 유지합니다. Linear shrinkage는 $c$가 작을 때(데이터가 상대적으로 풍부할 때) 오히려 30%대까지 떨어지고, $c$가 클 때는 93%까지 회복됩니다 — 5절에서 설명한 "역설"이 그대로 확인됩니다.
- **고유값 분산(dispersion) 효과**: 모집단 고유값들이 거의 동일한 극단적인 경우(항등행렬, $d=0$)에는 오히려 linear shrinkage(99.9%)가 nonlinear(99.4%)를 근소하게 앞섭니다 — 이 경우 선형 근사가 이미 정확하기 때문입니다. 하지만 고유값이 퍼질수록 linear shrinkage의 성능은 급격히 나빠져 55%까지 떨어지는 반면, nonlinear shrinkage는 96% 이상을 유지합니다.
- **정밀도행렬(precision matrix) 추정**: 공분산 행렬과 유사하게 큰 개선을 보이며, 정밀도행렬을 직접 nonlinear shrinkage로 추정하는 것이 공분산 행렬을 추정한 뒤 역행렬을 취하는 것보다 뚜렷하게 우수합니다.
- **비정규분포(fat tail)에 대한 강건성**: 정규분포 대신 자유도 3인 t-분포로 데이터를 생성해도 성능에 유의미한 차이가 없었습니다.
- **다양한 스펙트럼 형태**: 균등·U자형·종형·좌우 비대칭 등 다양한 모집단 고유값 분포 형태(베타분포로 생성)에서도 nonlinear shrinkage는 항상 88% 이상, 대부분 95% 이상을 기록했습니다.
- **고전적 점근이론(고정 차원, $n\to\infty$)**: $p=100$을 고정하고 $n$을 10,000까지 늘려도, linear shrinkage의 이득은 이론상 0으로 수렴하는 반면(최적 $\delta^*$ 자체가 0으로 수렴하므로), nonlinear shrinkage와 오라클은 여전히 **약 60%의 개선**을 유지합니다. 이는 $n$이 커질수록 표본 고유값들이 각 모집단 고유값 주변으로 뭉쳐 군집(cluster)을 이루는 "spectral separation" 현상 때문인데, 군집 내부에서는 여전히 정교한 비선형 보정이 필요하기 때문입니다.
- **다른 추정량과의 비교**: Stein(1975), Haff(1980), Won et al.(2009, 조건수 제약 MLE) 등 기존 문헌의 추정량들과 비교했을 때도 nonlinear shrinkage가 가장 우수했으며, 이 방법의 교차검증(cross-validation) 버전이 근소한 차이로 2위를 차지했습니다.

## 11. Linear vs Nonlinear 요약 비교

| 구분 | Linear Shrinkage (2003) | Nonlinear Shrinkage (2012) |
|---|---|---|
| 축소 방식 | 모든 고유값에 동일한 비율 $\delta$ 적용 | 고유값 위치별로 다른 비율 적용 |
| 계산 방법 | 닫힌 형태(closed-form) 공식 | 비볼록 최적화(SLP) 필요 |
| 구현 난이도 | 쉬움 | 상대적으로 복잡, 전용 솔버 필요 |
| 고유값이 밀집된 경우 | 거의 최적에 가까움 | 거의 동일한 성능 |
| 고유값이 넓게 퍼진 경우 | 성능 급격히 저하 | 안정적으로 높은 성능 유지 |
| $p/n$이 작은 경우(데이터 풍부) | 상대적으로 성능 저조 | 안정적 |
| 정밀도행렬 직접 추정 | 별도 이론 없음 | 지원 (더 우수) |

## 12. 마무리

두 논문이 공통적으로 전하는 메시지는 명확합니다: **차원이 크지 않은 극히 예외적인 상황이 아니라면, 표본공분산행렬을 있는 그대로 쓰지 말라**는 것입니다. Linear shrinkage는 구현이 간단하고 닫힌 형태 공식이 있어 빠르게 적용할 수 있으며, 모집단 고유값이 서로 비슷하거나 표본 크기가 차원에 비해 크지 않은 상황에서는 이미 개선분의 대부분을 가져다줍니다. 반면 nonlinear shrinkage는 계산이 더 복잡하지만, 고유값이 넓게 퍼져 있거나 표본이 상대적으로 풍부한 경우까지 포함해 전 구간에서 안정적으로 우수한 성능을 보장합니다. 두 방법 모두 표본공분산행렬을 압도적으로 능가한다는 점에서, "표본공분산행렬을 그대로 쓰는 것"은 더 이상 정당화되기 어렵습니다.

---

## 부록 A. 최적 선형 축소 강도($\delta^*$) 유도

목표는 손실함수

$$
L(\delta) = \|\delta F + (1-\delta)S - \Sigma\|^2 = \sum_{i,j}\left(\delta f_{ij} + (1-\delta)s_{ij} - \sigma_{ij}\right)^2
$$

의 기댓값(위험, risk) $R(\delta) = \mathbb{E}[L(\delta)]$을 최소화하는 $\delta$를 찾는 것입니다.

**1단계. 항 재정리.** 각 원소 $(i,j)$에 대해

$$
\delta f_{ij} + (1-\delta)s_{ij} - \sigma_{ij} = (s_{ij}-\sigma_{ij}) + \delta(f_{ij}-s_{ij})
$$

로 쓸 수 있습니다. $a_{ij} := s_{ij}-\sigma_{ij}$ (평균 0인 $S$의 추정오차), $b_{ij} := f_{ij}-s_{ij}$로 놓으면, $(i,j)$ 원소의 손실은 $(a_{ij}+\delta b_{ij})^2$이고, 기댓값을 취하면

$$
\mathbb{E}\left[(a_{ij}+\delta b_{ij})^2\right] = \text{Var}(s_{ij}) + 2\delta\, \text{Cov}(s_{ij}-\sigma_{ij},\, f_{ij}-s_{ij}) + \delta^2\, \mathbb{E}\left[(f_{ij}-s_{ij})^2\right]
$$

**2단계. 점근적 크기 분석 ($n\to\infty$, $p$ 고정).** $S$는 $\sqrt n$-일관추정량이므로 $\text{Var}(s_{ij}) = O(1/n)$입니다. $F$ 역시 모집단 등상관행렬 $\Phi$의 $\sqrt n$-일관추정량이므로, $f_{ij}-s_{ij}$를 다음과 같이 분해할 수 있습니다.

$$
f_{ij}-s_{ij} = \underbrace{(f_{ij}-\phi_{ij})}_{O_p(1/\sqrt n)} - \underbrace{(s_{ij}-\sigma_{ij})}_{O_p(1/\sqrt n)} + \underbrace{(\phi_{ij}-\sigma_{ij})}_{O(1),\ \text{고정된 상수}}
$$

앞의 두 항은 $n\to\infty$일 때 0으로 수렴하지만, 마지막 항 $(\phi_{ij}-\sigma_{ij})$는 **target의 모형 오설정**을 나타내는 고정된 값으로, $n$과 무관합니다. 따라서

$$
\mathbb{E}\left[(f_{ij}-s_{ij})^2\right] \;\approx\; (\phi_{ij}-\sigma_{ij})^2 + O(1/n) \;=:\; \gamma_{ij} + O(1/n)
$$

한편 $S$는 불편($\mathbb{E}[s_{ij}]=\sigma_{ij}$)이므로

$$
\text{Cov}(s_{ij}-\sigma_{ij},\,f_{ij}-s_{ij}) = \text{Cov}(s_{ij},f_{ij}) - \text{Var}(s_{ij})
$$

이고, 두 항 모두 $O(1/n)$ 크기입니다. 다음과 같이 정의합니다.

$$
\pi_{ij} := \lim_{n\to\infty} n\cdot\text{Var}(s_{ij}), \qquad \rho_{ij} := \lim_{n\to\infty} n\cdot\text{Cov}(f_{ij},s_{ij})
$$

**3단계. $\delta$에 대한 이차식으로 정리.** 위 결과들을 종합하면, $(i,j)$ 원소의 위험은

$$
R_{ij}(\delta) \approx \gamma_{ij}\,\delta^2 + \frac{2\delta(\rho_{ij}-\pi_{ij}) + \pi_{ij}}{n}
$$

모든 $(i,j)$에 대해 합산하고 $\gamma=\sum\gamma_{ij}$, $\rho=\sum\rho_{ij}$, $\pi=\sum\pi_{ij}$로 정의하면,

$$
R(\delta) \approx \gamma\delta^2 + \frac{2(\rho-\pi)}{n}\delta + \frac{\pi}{n}
$$

$\gamma = \sum_{ij}(\phi_{ij}-\sigma_{ij})^2 \geq 0$이므로 이는 아래로 볼록한(convex) 이차함수입니다.

**4단계. 최소화.** $\dfrac{dR}{d\delta} = 2\gamma\delta + \dfrac{2(\rho-\pi)}{n} = 0$을 풀면

$$
\delta^* = \frac{\pi-\rho}{n\gamma} = \frac{\kappa}{n}, \qquad \kappa := \frac{\pi-\rho}{\gamma}
$$

가 유도됩니다. 이것이 본문에서 소개한 공식입니다.

**실전 추정.** $\pi,\rho,\gamma$는 표본으로부터 다음과 같은 일관추정량으로 계산합니다 (Ledoit & Wolf, 2003).

$$
\hat\pi = \sum_{i,j}\hat\pi_{ij}, \qquad \hat\pi_{ij} = \frac{1}{n}\sum_{t=1}^{n}\left\{(y_{it}-\bar y_i)(y_{jt}-\bar y_j) - s_{ij}\right\}^2
$$

$$
\hat\gamma = \sum_{i,j}(f_{ij}-s_{ij})^2
$$

$\hat\rho$는 대각/비대각 항을 나누어 델타 방법(delta method)으로 유도되는 다소 복잡한 형태를 가지며, 대각 성분은 $\hat\pi_{ii}$와 같고 비대각 성분은 $\bar r$과 $s_{ii}, s_{jj}, s_{ij}$의 표본 변동 간 공분산 추정치들의 조합으로 계산됩니다. 이 세 추정량을 결합하면 $\hat\kappa=(\hat\pi-\hat\rho)/\hat\gamma$를 얻고, 최종적으로

$$
\hat\delta^* = \max\left\{0,\ \min\left\{\frac{\hat\kappa}{n},\ 1\right\}\right\}
$$

를 사용합니다 (이론상 $\hat\kappa/n$이 $[0,1]$을 벗어나는 것은 매우 드물지만, 유한표본에서 발생할 경우를 대비한 안전장치입니다).

## 부록 B. 회전등변 오라클 공식 유도

목표는, 고정된 직교행렬 $U$ (예: $S$의 고유벡터 행렬, $U'U=UU'=I_p$)에 대해, 대각행렬 $D=\text{diag}(d_1,\dots,d_p)$를 움직여

$$
f(D) = \|UDU' - \Sigma\|^2_F = \text{trace}\left[(UDU'-\Sigma)^2\right]
$$

를 최소화하는 것입니다 ($UDU'-\Sigma$는 대칭행렬이므로 전치해도 자기 자신).

**1단계. 전개.**

$$
f(D) = \text{trace}[UDU'UDU'] - 2\,\text{trace}[UDU'\Sigma] + \text{trace}[\Sigma^2]
$$

$U'U=I_p$이므로 $UDU'UDU' = UD^2U'$이고, trace의 순환성(cyclic property)에 의해

$$
\text{trace}[UD^2U'] = \text{trace}[D^2U'U] = \text{trace}[D^2] = \sum_i d_i^2
$$

또한 순환성을 이용해

$$
\text{trace}[UDU'\Sigma] = \text{trace}[DU'\Sigma U]
$$

**2단계. $\tilde\Sigma := U'\Sigma U$ 도입.** $\tilde\Sigma$는 $p\times p$ 대칭행렬이며(일반적으로 대각행렬은 아닙니다 — $U$가 $\Sigma$의 고유벡터가 아니라 $S$의 고유벡터이기 때문입니다), 그 대각원소를 $\tilde\sigma_i := \tilde\Sigma_{ii} = u_i'\Sigma u_i$라 합시다. $D$가 대각행렬이므로

$$
\text{trace}[D\tilde\Sigma] = \sum_i d_i \tilde\sigma_i
$$

**3단계. 완전제곱식.** 지금까지의 결과를 모으면

$$
f(D) = \sum_i d_i^2 - 2\sum_i d_i\tilde\sigma_i + \text{trace}[\Sigma^2] = \sum_i (d_i-\tilde\sigma_i)^2 - \sum_i \tilde\sigma_i^2 + \text{trace}[\Sigma^2]
$$

**4단계. 최소화.** $d_i$에 의존하는 항은 오직 $\sum_i(d_i-\tilde\sigma_i)^2$ 뿐이며, 이는 각 $i$에 대해 독립적으로 $d_i = \tilde\sigma_i$일 때 (즉 0일 때) 최소가 됩니다. 따라서

$$
d_i^{*} = \tilde\sigma_i = u_i'\Sigma u_i \qquad (i=1,\dots,p)
$$

가 $f(D)$를 최소화하는 유일해입니다. 즉 회전등변 추정량 클래스 안에서 프로베니우스 노름 기준 최적의 대각원소는 정확히 $u_i'\Sigma u_i$입니다.

같은 논리를 $\Sigma$ 대신 $\Sigma^{-1}$에 적용하면(전개 과정이 완전히 동일합니다), 정밀도행렬에 대한 오라클도

$$
a_i^{*} = u_i'\Sigma^{-1}u_i
$$

로 얻어짐을 바로 알 수 있습니다.

---

## 참고문헌

- Ledoit, O. and Wolf, M. (2003). *Honey, I Shrunk the Sample Covariance Matrix.*
- Ledoit, O. and Wolf, M. (2004). A well-conditioned estimator for large-dimensional covariance matrices. *Journal of Multivariate Analysis*, 88(2), 365–411.
- Ledoit, O. and Wolf, M. (2012). Nonlinear shrinkage estimation of large-dimensional covariance matrices. *The Annals of Statistics*, 40(2), 1024–1060.
- Ledoit, O. and Péché, S. (2011). Eigenvectors of some large sample covariance matrix ensembles. *Probability Theory and Related Fields*, 151, 233–264.
- Marchenko, V. A. and Pastur, L. A. (1967). Distribution of eigenvalues for some sets of random matrices. *Sbornik: Mathematics*, 1, 457–483.
- Silverstein, J. W. and Choi, S.-I. (1995). Analysis of the limiting spectral distribution of large-dimensional random matrices. *Journal of Multivariate Analysis*, 54, 295–309.
- Stein, C. (1956). Inadmissibility of the usual estimator for the mean of a multivariate normal distribution. *Proc. 3rd Berkeley Symposium*, 197–206.
- Michaud, R. (1989). The Markowitz optimization enigma: Is optimized optimal? *Financial Analysts Journal*, 45, 31–42.
