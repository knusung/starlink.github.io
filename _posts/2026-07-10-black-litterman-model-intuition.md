---
title: "Black-Litterman 모델의 직관: 균형에서 시작해 뷰(View)로 기울이기"
date: 2026-07-10 09:00:00 +0900
categories: [Quant Finance]
tags: [Black-Litterman, Bayesian Statistics, Portfolio Optimization, Asset Allocation, CAPM]
math: true
---

> *원문: He, Guangliang and Litterman, Robert (1999). "The Intuition Behind Black-Litterman Model Portfolios." Goldman Sachs Asset Management. SSRN: [https://ssrn.com/abstract=334304](https://ssrn.com/abstract=334304)*

Black-Litterman 모델은 1990년 Goldman Sachs의 Fischer Black과 Robert Litterman이 발표한 이후 자산배분 실무에서 가장 널리 쓰이는 모델 중 하나가 되었다. 하지만 원 논문들은 결과가 "왜" 그렇게 나오는지에 대한 직관을 충분히 설명하지 않았고, 많은 실무자들이 이 모델을 일종의 블랙박스로 취급해왔다. 이 글에서 정리하는 He and Litterman(1999)의 논문은 정확히 그 틈을 메운다. 이들은 제약이 없는 투자자의 최적 포트폴리오가 아주 단순한 형태 — **스케일된 시장균형 포트폴리오 + 뷰 포트폴리오들의 가중합** — 로 정확히 분해된다는 것을 보인다.

이 글은 통계/ML을 공부하는 대학원생 수준의 독자를 대상으로, Black-Litterman 모델을 베이지안 통계의 언어로 풀어 설명한다. 수식 전개가 필요한 부분은 결과와 직관 위주로 본문에서 다루고, 세부적인 유도·증명 과정은 모두 부록으로 분리했다.

**핵심 요약**
- 무제약 최적 포트폴리오는 $w^{*} = \dfrac{1}{1+\tau}\left(w_{eq} + P'\Lambda\right)$ 로, **"스케일된 시장균형 포트폴리오" + "투자자 뷰 포트폴리오들의 가중합"** 으로 정확히 분해된다.
- 각 뷰의 가중치 $\lambda_k$ 는 **뷰의 강도(기대수익률·확신도)에 비례**하고, **시장균형 및 다른 뷰와의 공분산에는 반비례(페널티)** 한다.
- 새로운 뷰가 시장에 이미 반영된 정보와 정확히 일치한다면(균형이 내재하는 기대수익률과 같다면), 그 뷰의 가중치는 정확히 0이 된다.

## 목차
1. [Black-Litterman 모델이 필요한 이유](#1-black-litterman-모델이-필요한-이유)
2. [모델 설정: 베이지안 프레임워크](#2-모델-설정-베이지안-프레임워크)
3. [핵심 결과: 무제약 최적 포트폴리오](#3-핵심-결과-무제약-최적-포트폴리오)
4. [뷰 가중치의 성질](#4-뷰-가중치의-성질)
5. [수치 예제로 확인하기](#5-수치-예제로-확인하기)
6. [정리하며](#6-정리하며)
7. [참고문헌](#7-참고문헌)
8. [부록: 수식 유도와 증명](#부록-수식-유도와-증명)

---

## 1. Black-Litterman 모델이 필요한 이유

표준적인 평균-분산(mean-variance) 최적화는 사용자가 전체 $N$개 자산에 대한 기대수익률 벡터 $\mu$를 직접 입력하고, 옵티마이저가 최적 비중을 계산해주는 방식이다. 문제는 기대수익률과 포트폴리오 비중 사이의 매핑이 매우 복잡하고 비선형적이어서, 기대수익률에 대한 아주 작은 가정 변화가 포트폴리오 비중에는 매우 극단적이고 비직관적인 결과를 낳는다는 점이다(이른바 "오차 극대화(error maximization)" 문제). 게다가 기대수익률 값을 어디서부터 시작해야 할지 자연스러운 출발점도 없다.

Black-Litterman 모델은 이 문제를 베이지안 방식으로 해결한다. 사용자가 전체 기대수익률을 직접 입력하는 대신, **시장균형(CAPM equilibrium)** 을 사전분포(prior)로 삼고, 투자자가 확신하는 **일부 뷰(view)** 만 입력하면, 모델이 이 둘을 결합해 전체 자산의 기대수익률과 최적 포트폴리오 비중을 동시에 산출해준다.

| 구분 | 전통적 평균-분산 최적화 | Black-Litterman 모델 |
|---|---|---|
| 필요한 입력 | 전체 $N$개 자산의 기대수익률 벡터 $\mu$ | 시장균형($\Pi$, 자동 계산) + 투자자가 확신하는 일부 뷰 $(P,Q,\Omega)$ |
| 출발점(중립 비중) | 명시적인 중립 출발점이 없음 | 시가총액 가중 시장균형 포트폴리오 $w_{eq}$ |
| "오차 극대화" 문제 | 기대수익률의 작은 오차가 비중에 극단적으로 증폭됨 | 뷰가 균형 대비 상대적으로 표현되어 안정적 |
| 결과 해석 가능성 | 왜 이런 비중이 나왔는지 역추적이 어려움 | $w^{*}=$ 균형 + 뷰의 가중합으로 분해되어 해석이 쉬움 |
| 불확실성 처리 | 기대수익률을 확정값(포인트 추정치)으로 취급 | 기대수익률 자체를 확률변수로 취급(베이지안) |
| 뷰가 하나도 없을 때 | 정의되지 않음(애초에 입력이 없으므로) | 정확히 (스케일된) 시장균형 포트폴리오로 수렴 |

예를 들어 "독일 주식이 나머지 유럽 시장보다 좋을 것 같다"는 견해를 표준 옵티마이저에 반영하려면, 독일의 기대수익률만 올리고 나머지는 그대로 두는 식으로 접근하기 쉽다. 그런데 뒤에서(5장) 직접 확인하겠지만, 그 결과로 나오는 포트폴리오는 독일 비중이 튀는 것은 물론, 아무 견해도 표명하지 않은 호주·캐나다·일본·미국의 비중까지 크게 흔들어 놓는다. 왜 그런 일이 벌어지는지 표준 옵티마이저의 산출물만 봐서는 파악하기 어렵다. Black-Litterman은 바로 이 "왜"에 대한 답을 준다.

## 2. 모델 설정: 베이지안 프레임워크

### 표기법

| 기호 | 의미 |
|---|---|
| $N$ | 전체 자산 개수 |
| $K$ | 투자자가 가진 뷰의 개수 ($K \le N$) |
| $r$ | $N\times1$ 자산 수익률 벡터, $r \sim N(\mu,\Sigma)$ |
| $\mu$ | (미지의, 확률변수로 취급되는) $N\times1$ 기대수익률 벡터 |
| $\Sigma$ | $N\times N$ 자산 수익률 공분산행렬 (기지의 상수로 가정) |
| $\delta$ | 시장 평균 위험회피계수 |
| $w_{eq}$ | 시가총액 기준 시장균형 포트폴리오 |
| $\Pi$ | 균형 초과수익률(리스크 프리미엄) 벡터 |
| $\tau$ | 균형 추정치 $\Pi$에 대한 불확실성 크기(스칼라) |
| $P$ | $K\times N$ 뷰 포트폴리오 행렬 (각 행이 하나의 뷰) |
| $Q$ | $K\times1$ 뷰의 기대수익률 벡터 |
| $\Omega$ | $K\times K$ 뷰의 불확실성(오차) 공분산행렬, 대각행렬 |
| $\bar\mu$ | 사후 기대수익률 벡터 |
| $\bar M^{-1}$ | 사후 공분산행렬 |
| $\bar\Sigma$ | 수익률 자체의 총분산, $\Sigma+\bar M^{-1}$ |
| $w^{*}$ | 무제약 투자자의 최적 포트폴리오 |
| $\Lambda$ | 각 뷰 포트폴리오에 대한 가중치 벡터 ($K\times1$) |
| $A$ | $K\times K$ 행렬, $A=\Omega/\tau+P\Sigma P'/(1+\tau)$ |

### 사전분포: 시장균형

자산 수익률은 다음과 같이 정규분포를 따른다고 가정한다.

$$ r \sim N(\mu, \Sigma) $$

시장 전체가 시가총액 가중 포트폴리오 $w_{eq}$를 보유하는 균형 상태에서, CAPM 논리에 따라 균형 리스크 프리미엄은 다음과 같이 역산된다(이를 **역최적화, reverse optimization** 라고 부른다).

$$ \Pi = \delta\, \Sigma\, w_{eq} $$

Black-Litterman의 베이지안 사전분포는 진짜 기대수익률 $\mu$가 이 균형값 $\Pi$를 중심으로 분포한다고 본다.

$$ \mu = \Pi + \epsilon^{(e)}, \qquad \epsilon^{(e)} \sim N(0,\ \tau\Sigma) $$

여기서 $\tau$는 "균형 추정치를 얼마나 믿는가"를 나타내는 스칼라다. $\tau$가 작을수록 사전분포가 $\Pi$ 주변에 촘촘히 몰려 있다는 뜻이고, 이는 곧 균형에 대한 확신이 크다는 의미다.

### 뷰: 새로운 정보

투자자는 여기에 더해 $K$개의 주관적인 뷰를 갖고 있다. 하나의 뷰는 "어떤 포트폴리오 $p$의 기대수익률이 평균 $q$, 표준편차 $\omega$인 정규분포를 따른다"는 진술이다. $K$개의 뷰를 모아 행렬로 쓰면,

$$ P\mu = Q + \epsilon^{(v)}, \qquad \epsilon^{(v)} \sim N(0,\ \Omega) $$

이며 $\Omega$는 (일반성을 잃지 않고) 대각행렬로 둘 수 있다. 사전분포의 오차 $\epsilon^{(e)}$와 뷰의 오차 $\epsilon^{(v)}$는 서로 독립이라고 가정한다.

$$ \begin{pmatrix}\epsilon^{(e)}\\ \epsilon^{(v)}\end{pmatrix} \sim N\!\left(0,\ \begin{pmatrix}\tau\Sigma & 0\\ 0 & \Omega\end{pmatrix}\right) $$

### 사후분포: 베이지안 결합

사전분포와 뷰를 결합하면, 기대수익률의 사후분포는 $N(\bar\mu, \bar M^{-1})$ 이며,

$$ \bar\mu = \Big[(\tau\Sigma)^{-1} + P'\Omega^{-1}P\Big]^{-1}\Big[(\tau\Sigma)^{-1}\Pi + P'\Omega^{-1}Q\Big] $$

$$ \bar M^{-1} = \Big[(\tau\Sigma)^{-1} + P'\Omega^{-1}P\Big]^{-1} $$

기대수익률 자체가 불확실한 확률변수이므로, 수익률의 총분산은 원래의 $\Sigma$에 이 사후 불확실성 $\bar M^{-1}$이 더해진다.

$$ r \sim N(\bar\mu, \bar\Sigma), \qquad \bar\Sigma = \Sigma + \bar M^{-1} $$

> **베이즈 통계 관점에서 보면.** 위 사후분포 식은 낯설지 않다. 서로 독립인 두 개의 가우시안 정보원 — ① 사전분포 $N(\Pi,\tau\Sigma)$, ② 선형관측모형 $P\mu=Q+\epsilon^{(v)}$ — 를 결합하는, 정밀도(precision, 공분산의 역행렬) 가중 평균 공식 그 자체다. 정확히 Kalman 필터의 측정 갱신(measurement update)이나, 가우시안 사전분포를 갖는 베이지안 선형회귀(리지 회귀의 MAP 추정과 동일한 구조)에서 등장하는 식과 같다. Black-Litterman의 "새로움"은 새로운 통계 기법이 아니라, **자산배분 문제에 이 표준적인 베이지안 결합을 적용하되, "관측"을 개별 자산이 아니라 임의의 포트폴리오(뷰)에 대해 정의**했다는 데 있다.

| 구분 | 사전분포 (Prior) | 뷰 (Views) |
|---|---|---|
| 정보의 출처 | CAPM 시장균형(전체 투자자의 집단적 판단) | 개별 투자자의 주관적 견해 |
| 확률적 표현 | $\mu \sim N(\Pi,\tau\Sigma)$ | $P\mu = Q+\epsilon^{(v)},\ \epsilon^{(v)}\sim N(0,\Omega)$ |
| 불확실성 파라미터 | $\tau$ (스칼라, 균형 추정 전체의 불확실성) | $\Omega$ ($K\times K$ 대각행렬, 뷰별 개별 불확실성) |
| 다루는 대상 | $N$개 자산 전체의 기대수익률 | $K$개 뷰 포트폴리오에 대한 진술 ($K\le N$) |
| 없을 경우 | 모델 성립 불가(반드시 필요) | 없어도 무방 → 사후분포=사전분포, $w^{*}=w_{eq}/(1+\tau)$ |

이 전체 흐름을 그림으로 요약하면 다음과 같다.

![Black-Litterman 베이지안 프레임워크](/assets/img/posts/black-litterman-model-intuition/bl-framework-diagram.png)
*그림 1. 사전분포(시장균형)와 뷰가 베이지안 결합을 통해 사후분포로, 다시 최적 포트폴리오로 이어지는 흐름.*

## 3. 핵심 결과: 무제약 최적 포트폴리오

### 최적화 문제

사후분포 $N(\bar\mu,\bar\Sigma)$가 주어졌을 때, 위험회피계수 $\delta$를 가진 무제약 투자자는 표준적인 평균-분산 최적화를 푼다.

$$ \max_{w}\ w'\bar\mu - \frac{\delta}{2}w'\bar\Sigma w $$

1계 조건(FOC)으로부터 익숙한 형태의 해를 얻는다.

$$ w^{*} = \frac{1}{\delta}\bar\Sigma^{-1}\bar\mu $$

이 식 자체는 표준 평균-분산 최적화와 다를 게 없다. He and Litterman(1999)의 기여는 여기서 멈추지 않고, $\bar\mu$와 $\bar\Sigma$의 구체적 형태(식 8~10)를 대입해 이 해를 훨씬 더 해석하기 쉬운 형태로 정리했다는 데 있다. (대입·정리 과정은 다소 손이 많이 가는 행렬대수라, 결과만 아래에 제시하고 전체 유도는 [부록 A](#부록-a-무제약-최적-포트폴리오-공식-유도)로 미룬다.)

### 결과: 균형 + 뷰의 가중합

$$ w^{*} = \frac{1}{1+\tau}\Big(w_{eq} + P'\Lambda\Big) $$

여기서 $\Lambda$는 각 뷰 포트폴리오에 대한 가중치를 담은 $K\times1$ 벡터로,

$$ \Lambda = \frac{\tau\Omega^{-1}Q}{\delta} \ -\ A^{-1}\frac{P\Sigma}{1+\tau}w_{eq} \ -\ A^{-1}\frac{P\Sigma}{1+\tau}P'\frac{\tau\Omega^{-1}Q}{\delta}, \qquad A = \frac{\Omega}{\tau} + \frac{P\Sigma P'}{1+\tau} $$

이제 경제적 직관이 아주 분명해진다. **투자자는 먼저 스케일된 시장균형 포트폴리오 $w_{eq}/(1+\tau)$를 보유하고, 그 위에 자신이 가진 각 뷰를 표현하는 포트폴리오를 $\Lambda$만큼씩 추가로 편입한다.** $\Lambda$의 세 항은 각각 명확한 의미를 갖는다.

1. **뷰의 강도** ($\tau\Omega^{-1}Q/\delta$): 뷰의 기대수익률 $q_k$가 높을수록, 또는 확신도(정밀도) $\omega_k^{-1}/\tau$가 높을수록 가중치가 커진다.
2. **균형과의 공분산에 대한 페널티** ($-A^{-1}P\Sigma w_{eq}/(1+\tau)$): 뷰 포트폴리오가 시장균형 포트폴리오와 이미 비슷한 방향이라면, 그 뷰는 "새로운 정보"를 덜 담고 있는 셈이므로 가중치가 깎인다.
3. **다른 뷰들과의 공분산에 대한 페널티** (마지막 항): 여러 뷰가 서로 상관되어 있다면 정보가 중복 계산되는 셈이므로, 그만큼 가중치가 깎인다.

즉 Black-Litterman은 "얼마나 강하게 확신하는가"와 "그 정보가 얼마나 새로운가(중복되지 않는가)"를 동시에 반영해 각 뷰에 자금을 배분해주는 셈이다.

## 4. 뷰 가중치의 성질

$\Lambda$가 포트폴리오 구성의 핵심이므로, 저자들은 이 가중치가 언제 양(+)이 되고 어떻게 변하는지를 두 가지 명제로 정리한다. 두 명제 모두 결과와 직관만 아래에 제시하고, 증명 전체는 [부록 B](#부록-b-property-31의-증명)와 [부록 C](#부록-c-property-32의-증명)에서 다룬다.

### Property 3.1 — 새로운 뷰가 추가될 때

기존 $K$개 뷰 $(P,Q,\Omega)$에 새로운 뷰 $(p,q,\omega)$ 하나를 추가한다고 하자. 새 뷰의 가중치 $\hat\lambda_{K+1}$는 다음 부호 관계를 만족한다.

$$ \operatorname{sign}(\hat\lambda_{K+1}) = \operatorname{sign}\big(q - p'\tilde\mu\big), \qquad \tilde\mu \equiv \Sigma\bar\Sigma^{-1}\bar\mu = \delta\,\Sigma w^{*} $$

여기서 $\tilde\mu$는 "기존 $K$개 뷰만 있었을 때의 최적 포트폴리오가 내재하고 있는 기대수익률"이다. 직관은 명확하다: **새 뷰가 기존 모델이 이미 암시하고 있던 값보다 더 강세(약세)라면 양(음)의 가중치를 받고, 정확히 일치한다면($q=p'\tilde\mu$) 가중치는 0이 된다** — 이미 시장에 반영된 정보를 다시 말한 것에 불과하기 때문이다. 이 경우 기존 뷰들의 가중치도 전혀 바뀌지 않는다.

### Property 3.2 — 뷰의 강도와 확신도가 변할 때

특정 뷰 $k$에 대해,

- $\lambda_k$ 는 그 뷰의 기대수익률 $q_k$에 대해 **증가함수**다. (더 강세로 볼수록 더 많이 투자한다.)
- $|\lambda_k|$ 는 그 뷰의 확신도(정밀도) $\omega_k^{-1}$에 대해 **증가함수**다. (더 확신할수록 그 뷰에 더 크게 베팅한다.)

| 구분 | Property 3.1 | Property 3.2 |
|---|---|---|
| 질문 | 새 뷰를 추가하면 가중치는 어떻게 정해지나? | 기존 뷰의 파라미터가 바뀌면 가중치는 어떻게 변하나? |
| 핵심 변수 | 새 뷰의 기대수익률 $q$ vs 내재 기대수익률 $p'\tilde\mu$ | 뷰 $k$의 기대수익률 $q_k$, 확신도 $\omega_k^{-1}$ |
| 결과 | $\operatorname{sign}(\lambda_{K+1})=\operatorname{sign}(q-p'\tilde\mu)$ | $\lambda_k$는 $q_k$의 증가함수, $\lvert\lambda_k\rvert$는 $\omega_k^{-1}$의 증가함수 |
| 특수 사례 | $q=p'\tilde\mu \Rightarrow \lambda_{K+1}=0$, 다른 뷰 가중치 불변 | — |
| 5장의 대응 예시 | Table 8 (캐나다>일본 뷰가 이미 내재된 경우) | Table 6 (강세↑), Table 7 (확신도↓) |

한 가지 덧붙이면, $A^{-1}$은 일반적으로 완전한(대각이 아닌) 행렬이므로, 뷰 $k$의 파라미터를 바꾸면 이론적으로는 $\lambda_k$뿐 아니라 다른 모든 $\lambda_i$도 함께 미세하게 움직일 수 있다 ($\partial\lambda_i/\partial q_k = (A^{-1})_{ik}/\delta$). 이 교차 효과는 뷰 포트폴리오들끼리 서로 상관되어 있을 때(자산이 겹치지 않아도 공분산 $\Sigma$를 통해 상관될 수 있다) 나타나며, 5장의 수치 예제에서 실제로 관찰된다.

## 5. 수치 예제로 확인하기

이제 위 결과를 원 논문의 수치 예제로 직접 확인해본다. 시장은 7개 주요 산업국의 주가지수로 구성되며, 위험회피계수는 $\delta=2.5$, 사전 불확실성은 $\tau=0.05$(20년치 데이터로 균형을 추정했을 때의 신뢰수준에 대응)로 둔다.

### 시장 데이터

**표 1. 7개국 주가지수 수익률 상관계수**

| | Australia | Canada | France | Germany | Japan | UK | USA |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Australia** | 1.000 | 0.488 | 0.478 | 0.515 | 0.439 | 0.512 | 0.491 |
| **Canada** | 0.488 | 1.000 | 0.664 | 0.655 | 0.310 | 0.608 | 0.779 |
| **France** | 0.478 | 0.664 | 1.000 | 0.861 | 0.355 | 0.783 | 0.668 |
| **Germany** | 0.515 | 0.655 | 0.861 | 1.000 | 0.354 | 0.777 | 0.653 |
| **Japan** | 0.439 | 0.310 | 0.355 | 0.354 | 1.000 | 0.405 | 0.306 |
| **UK** | 0.512 | 0.608 | 0.783 | 0.777 | 0.405 | 1.000 | 0.652 |
| **USA** | 0.491 | 0.779 | 0.668 | 0.653 | 0.306 | 0.652 | 1.000 |

![7개국 상관계수 히트맵](/assets/img/posts/black-litterman-model-intuition/correlation-heatmap.png)
*그림 2. 유럽 3국(France·Germany·UK)의 상관관계가 특히 높고, 일본은 상대적으로 독립적으로 움직인다.*

**표 2. 변동성, 시가총액 비중, 균형 리스크 프리미엄**

| 국가 | 변동성 $\sigma$ (%) | 시가총액 비중 $w_{eq}$ (%) | 균형 리스크 프리미엄 $\Pi$ (%) |
|---|---:|---:|---:|
| Australia | 16.0 | 1.6 | 3.9 |
| Canada | 20.3 | 2.2 | 6.9 |
| France | 24.8 | 5.2 | 8.4 |
| Germany | 27.1 | 5.5 | 9.0 |
| Japan | 21.0 | 11.6 | 4.3 |
| UK | 20.0 | 12.4 | 6.8 |
| USA | 18.7 | 61.5 | 7.6 |

### 왜 전통적 접근은 문제가 되는가

투자자가 "독일 주식이 나머지 유럽 시장 대비 연 5% 초과 수익을 낼 것"이라는 견해를 가졌다고 하자. Black-Litterman을 모르는 투자자가 표준 옵티마이저에 이 견해를 최대한 신중하게 반영하려 해도(독일의 기대수익률만 살짝 올리고 나머지는 균형값 유지), 결과는 다음과 같이 나온다.

**표 3. 전통적 평균-분산 접근으로 독일 뷰를 반영한 결과**

| 국가 | 기대수익률 $\mu$ (%) | 최적비중 $w_{opt}$ (%) | $\mu-\Pi$ | $w_{opt}-w_{eq}$ |
|---|---:|---:|---:|---:|
| Australia | 3.9 | −5.1 | 0.0 | −6.7 |
| Canada | 6.9 | −2.3 | 0.0 | −4.5 |
| France | 7.6 | −50.1 | −0.8 | −55.3 |
| Germany | 11.5 | 83.6 | 2.4 | 78.1 |
| Japan | 4.3 | 14.9 | 0.0 | 3.3 |
| UK | 6.0 | −22.8 | −0.8 | −35.2 |
| USA | 7.6 | 66.6 | 0.0 | 5.1 |

기대수익률은 아주 조금만 바뀌었는데($\mu-\Pi$ 열이 대부분 0이거나 매우 작다), 포트폴리오 비중은 극단적으로 흔들린다. 게다가 아무 견해도 표명하지 않은 호주·캐나다·일본·미국의 비중까지 크게 달라진다. 이것이 바로 "오차 극대화" 문제다.

### Black-Litterman으로 같은 뷰를 표현하면

같은 견해를 "독일 주식 롱 + 나머지 유럽 시장 시가총액 가중 숏"으로 구성된 포트폴리오의 기대수익률이 5%라는 뷰로 표현하면 결과는 다음과 같다.

**표 4. 뷰 1개(독일 vs 유럽)를 반영한 Black-Litterman 결과**

| 국가 | 뷰 포트폴리오 $p$ (%) | 사후 기대수익률 $\bar\mu$ (%) | 최적비중 $w^{*}$ (%) | 편차 $w^{*}-\frac{w_{eq}}{1+\tau}$ |
|---|---:|---:|---:|---:|
| Australia | 0.0 | 4.3 | 1.5 | 0.0 |
| Canada | 0.0 | 7.6 | 2.1 | 0.0 |
| France | −29.5 | 9.3 | −4.0 | −8.9 |
| Germany | 100.0 | 11.0 | 35.4 | 30.2 |
| Japan | 0.0 | 4.5 | 11.0 | 0.0 |
| UK | −70.5 | 7.0 | −9.5 | −21.3 |
| USA | 0.0 | 8.1 | 58.6 | 0.0 |

($q=5.00,\ \omega/\tau=0.021,\ \lambda\approx0.302$)

훨씬 직관적이다. 뷰가 없는 국가(호주·캐나다·일본·미국)의 편차는 정확히 0이고, 편차가 있는 국가들은 정확히 뷰 포트폴리오 $p$에 $\lambda$를 곱한 값이다 (예: 독일 $100.0\times0.302\approx30.2$, 프랑스 $-29.5\times0.302\approx-8.9$). 3장에서 말한 "균형 + 뷰 포트폴리오 $\times\ \lambda$"라는 구조가 숫자로 정확히 확인된다.

### 뷰가 여러 개일 때

이제 두 번째 뷰 — "캐나다 주식이 미국 대비 연 3% 초과 수익" — 를 추가한다.

**표 5. 뷰 2개(독일 vs 유럽, 캐나다 vs 미국)를 반영한 결과**

| 국가 | 뷰1: 독일-유럽 | 뷰2: 캐나다-미국 | $\bar\mu$ (%) | $w^{*}$ (%) | 편차 |
|---|---:|---:|---:|---:|---:|
| Australia | 0.0 | 0.0 | 4.4 | 1.5 | 0.0 |
| Canada | 0.0 | 100.0 | 8.7 | 41.9 | 39.8 |
| France | −29.5 | 0.0 | 9.5 | −3.4 | −8.4 |
| Germany | 100.0 | 0.0 | 11.2 | 33.6 | 28.4 |
| Japan | 0.0 | 0.0 | 4.6 | 11.0 | 0.0 |
| UK | −70.5 | 0.0 | 7.0 | −8.2 | −20.0 |
| USA | 0.0 | −100.0 | 7.5 | 18.8 | −39.8 |
| **$q$** | 5.00 | 3.00 | | | |
| **$\omega/\tau$** | 0.021 | 0.017 | | | |
| **$\lambda$** | 0.298 | 0.418 | | | |

두 개의 독립적인 뷰를 추가했을 뿐인데, 캐나다·미국·프랑스·독일·영국의 비중이 모두 크게 바뀐다 — 이는 문제가 아니라 정확히 의도된 결과다. 각 뷰 포트폴리오가 서로 다른 자산에 걸쳐 있고, 자산 간 상관관계(표 1)를 통해 영향이 전파되기 때문이다. 흥미로운 점은 뷰 1의 가중치가 표 4의 0.302에서 0.298로 살짝 낮아졌다는 것 — 캐나다/미국 포트폴리오가 프랑스·독일·영국과 상관되어 있어(표 1), 정보가 부분적으로 겹치기 때문이다.

![균형 비중 vs Black-Litterman 최적 비중](/assets/img/posts/black-litterman-model-intuition/equilibrium-vs-bl-weights.png)
*그림 3. 단 두 개의 뷰만으로 균형 대비 비중이 얼마나 크게 기울어지는지(tilt) 확인할 수 있다.*

### 뷰가 바뀌면 가중치는 어떻게 달라지는가

Property 3.1, 3.2를 아래 세 가지 변형을 통해 직접 확인해보자.

- **표 6**: 캐나다/미국 뷰를 더 강세로 바꾼다 (3% → 4%).
- **표 7**: 독일/유럽 뷰의 확신도를 낮춘다 ($\omega/\tau$를 0.021 → 0.043으로, 절반의 확신도).
- **표 8**: "캐나다가 일본 대비 연 4.12% 초과 수익"이라는 세 번째 뷰를 추가한다. 그런데 4.12%는 표 7의 모델이 이미 암시하고 있던 캐나다-일본 스프레드와 정확히 일치하는 값이다.

**표 9. 세 시나리오의 가중치 변화 요약**

| 시나리오 | 뷰1: 독일>유럽 $(q_1,\ \omega/\tau_1)$ | 뷰2: 캐나다>미국 $(q_2,\ \omega/\tau_2)$ | 뷰3: 캐나다>일본 $(q_3,\ \omega/\tau_3)$ | $\lambda_1$ | $\lambda_2$ | $\lambda_3$ | 비고 |
|---|---|---|---|---:|---:|---:|---|
| 표 5 (기본, 뷰 2개) | 5.00%, 0.021 | 3.00%, 0.017 | — | 0.298 | 0.418 | — | 기준 케이스 |
| 표 6 (뷰2 강세 ↑) | 5.00%, 0.021 | 4.00%, 0.017 | — | 0.292 | 0.538 | — | $q_2$: 3%→4%, $\lambda_2$ 증가 (Property 3.2) |
| 표 7 (뷰1 확신도 ↓) | 5.00%, 0.043 | 4.00%, 0.017 | — | 0.193 | 0.544 | — | $\omega/\tau_1$: 0.021→0.043, $\lvert\lambda_1\rvert$ 감소 (Property 3.2) |
| 표 8 (뷰3 추가, 이미 내재) | 5.00%, 0.043 | 4.00%, 0.017 | 4.12%, 0.059 | 0.193 | 0.544 | **0.000** | $q_3$가 균형+기존 뷰의 내재 기대수익률과 정확히 일치 → $\lambda_3=0$ (Property 3.1) |

표 8이 특히 인상적이다. 세 번째 뷰를 추가했는데도 포트폴리오는 단 1bp도 바뀌지 않는다 — 이미 시장(균형 + 기존 두 뷰)이 그 정보를 완전히 반영하고 있었기 때문이다. Property 3.1이 예측한 그대로다.

## 6. 정리하며

Black-Litterman 모델을 "기대수익률을 만들어내는 신비한 블랙박스"로 보는 대신, He and Litterman(1999)은 다음과 같은 아주 단순한 그림을 제시한다.

> 투자자는 먼저 (자신의 불확실성을 반영해 스케일된) 시장균형 포트폴리오를 보유한다. 그리고 자신이 가진 각각의 뷰를 표현하는 포트폴리오를, 그 뷰의 강도와 새로움(다른 정보와 얼마나 겹치지 않는가)에 비례하는 크기로 추가 편입한다.

제약이 있거나 위험 허용도가 세계 평균과 다른 투자자라도, 사후 기대수익률 $\bar\mu$와 공분산 $\bar\Sigma$만 표준 포트폴리오 최적화 패키지에 입력하면 동일한 논리를 적용할 수 있다. 표준 평균-분산 최적화와 달리, 제대로 구현된 Black-Litterman 모델은 언제나 "왜 이런 비중이 나왔는가"를 비교적 쉽게 설명할 수 있는 포트폴리오를 만들어준다.

## 7. 참고문헌

- He, Guangliang and Litterman, Robert (1999). *The Intuition Behind Black-Litterman Model Portfolios.* Goldman Sachs Asset Management. SSRN: [https://ssrn.com/abstract=334304](https://ssrn.com/abstract=334304)
- Black, Fischer and Litterman, Robert (1990). *Asset Allocation: Combining Investor Views with Market Equilibrium.* Fixed Income Research, Goldman Sachs.
- Black, Fischer and Litterman, Robert (1992). *Global Portfolio Optimization.* Financial Analysts Journal, 48(5), 28–43.
- Black, Fischer (1989). *Universal Hedging: Optimizing Currency Risk and Reward in International Equity Portfolios.* Financial Analysts Journal, 45(4), 16–22.
- Satchell, S. and Scowcroft, A. (1997). *A Demystification of the Black-Litterman Model: Managing Quantitative and Traditional Portfolio Construction.* Research Papers in Management Studies, Judge Institute of Management Studies, University of Cambridge.

---

## 부록: 수식 유도와 증명

### 부록 A. 무제약 최적 포트폴리오 공식 유도

FOC로부터 얻은 $w^{*}=\frac{1}{\delta}\bar\Sigma^{-1}\bar\mu$ 에, 사후분포 식(8),(9)를 대입하면

$$ w^{*} = \frac{1}{\delta}\bar\Sigma^{-1}\bar M^{-1}\Big[(\tau\Sigma)^{-1}\Pi + P'\Omega^{-1}Q\Big] $$

여기서 다음 행렬 항등식을 이용한다 — 임의의 가역행렬 $X,Y$(단, $X+Y$도 가역)에 대해 성립하는 관계다.

$$ (X+Y)^{-1} = Y^{-1} - Y^{-1}(X^{-1}+Y^{-1})^{-1}Y^{-1} $$

(양변에 $(X+Y)$를 곱해 직접 전개하면 항등식이 성립함을 바로 확인할 수 있다. Woodbury 항등식 계열의 한 특수한 형태로, Kalman 필터·가우시안 프로세스 회귀 등에서 흔히 쓰인다.) $X=\Sigma,\ Y=\bar M^{-1}$로 두면,

$$ \bar\Sigma^{-1} = (\Sigma+\bar M^{-1})^{-1} = \bar M - \bar M(\bar M+\Sigma^{-1})^{-1}\bar M $$

이를 이용해 $\bar\Sigma^{-1}\bar M^{-1}$ 항을 정리하면 (다소 지루하지만 기계적인 행렬대수를 거치면)

$$ \bar\Sigma^{-1}\bar M^{-1} = \frac{\tau}{1+\tau}\left[I - A^{-1}P\frac{\Sigma}{1+\tau}\right],\qquad A=\frac{\Omega}{\tau}+\frac{P\Sigma P'}{1+\tau} $$

가 되고, 이를 $w^{*}$ 식에 대입해 $\Pi=\delta\Sigma w_{eq}$ 를 이용해 정리하면 본문의 핵심 결과에 도달한다.

$$ w^{*} = \frac{1}{1+\tau}\big(w_{eq} + P'\Lambda\big) $$

$$ \Lambda = \frac{\tau\Omega^{-1}Q}{\delta} - A^{-1}\frac{P\Sigma}{1+\tau}w_{eq} - A^{-1}\frac{P\Sigma}{1+\tau}P'\frac{\tau\Omega^{-1}Q}{\delta} $$

이 식은 $A$의 정의를 이용해 대수적으로 다음과 같이 더 간단한 동치 형태로도 쓸 수 있다(부록 C의 Property 3.2 증명에서 이 형태를 사용한다).

$$ \Lambda = A^{-1}\left[\frac{Q}{\delta} - \frac{P\Sigma}{1+\tau}w_{eq}\right] $$

두 형태가 같다는 것은 $A^{-1}A=I$, 즉 $A^{-1}\big[\Omega/\tau + P\Sigma P'/(1+\tau)\big]=I$ 를 이용해 $A^{-1}\Omega/\tau = I - A^{-1}P\Sigma P'/(1+\tau)$ 로 바꿔 대입하면 확인할 수 있다.

또한 위험회피계수가 $\hat\delta$인 다른 투자자는 $\hat w^{*}=(\delta/\hat\delta)w^{*}$로, 목표 변동성이 $\sigma$인 투자자는 $\tilde w^{*}=\big(\sigma\delta/\sqrt{\bar\mu'\bar\Sigma^{-1}\bar\mu}\big)\,w^{*}$로 단순히 스케일링해서 얻을 수 있다. 제약이 있는 경우에는 $\bar\mu,\bar\Sigma$를 표준 포트폴리오 최적화 패키지에 입력하면 된다.

### 부록 B. Property 3.1의 증명

기존 $K$개 뷰에 새로운 뷰 $(p,q,\omega)$를 추가하면, 확장된 행렬들은 다음과 같다.

$$ \hat P=\begin{pmatrix}P\\p'\end{pmatrix},\ \ \hat Q=\begin{pmatrix}Q\\q\end{pmatrix},\ \ \hat\Omega=\begin{pmatrix}\Omega&0\\0&\omega\end{pmatrix},\ \ \hat A=\begin{pmatrix}A&b\\b'&c\end{pmatrix} $$

여기서 $b=\dfrac{P\Sigma}{1+\tau}p,\quad c=\dfrac{\omega}{\tau}+p'\dfrac{\Sigma}{1+\tau}p$ 이다. 블록행렬의 역행렬 공식을 쓰면

$$ \hat A^{-1} = \begin{pmatrix} A^{-1}+A^{-1}bb'A^{-1}/d & -A^{-1}b/d\\ -b'A^{-1}/d & 1/d \end{pmatrix},\qquad d=c-b'A^{-1}b $$

$\Lambda$의 확장 버전 식에 $\hat A^{-1}$의 마지막 행을 적용하고, 행렬 곱을 전개해 정리하면(논문에서 "몇 단계의 대수 전개"라고만 언급하는, 다소 지루하지만 기계적인 계산이다)

$$ \hat\lambda_{K+1} = \frac{q - p'\Sigma\bar\Sigma^{-1}\bar\mu}{(c-b'A^{-1}b)\,\delta} $$

를 얻고, $\hat A^{-1}$의 앞 $K$개 행을 적용하면 나머지 뷰들의 가중치가

$$ \hat\Lambda_{1:K} = \Lambda - \hat\lambda_{K+1}A^{-1}b $$

로 정리된다. 즉 $\hat\Lambda = \begin{pmatrix}\Lambda-\hat\lambda_{K+1}A^{-1}b\\ \hat\lambda_{K+1}\end{pmatrix}$.

이제 부호를 따져보자. $\Omega$는 양정치(positive definite), $\hat P\Sigma\hat P'$는 양반정치(positive semi-definite)이므로 이 둘의 합인 $\hat A$는 양정치이고, 따라서 $\hat A^{-1}$도 양정치다. $d=c-b'A^{-1}b$는 $\hat A^{-1}$의 한 대각원소이므로 반드시 양수다. 따라서

$$ \operatorname{sign}(\hat\lambda_{K+1}) = \operatorname{sign}\big(q-p'\Sigma\bar\Sigma^{-1}\bar\mu\big) = \operatorname{sign}\big(q-p'\tilde\mu\big),\qquad \tilde\mu\equiv\Sigma\bar\Sigma^{-1}\bar\mu=\delta\Sigma w^{*} $$

가 성립한다. $\tau,\Omega\to 0$이면서 $\Omega/\tau$가 유한하게 유지되는 극한에서는 $\tilde\mu\to\bar\mu$이므로, 이 조건은 $q-p'\bar\mu$의 부호와 같아진다. $\blacksquare$

### 부록 C. Property 3.2의 증명

부록 A에서 얻은 간단한 형태 $\Lambda = A^{-1}\big[Q/\delta - P\Sigma w_{eq}/(1+\tau)\big]$ 에서 출발한다.

**$q_k$에 대한 증가성.** $Q$의 $k$번째 원소에 대해 미분하면

$$ \frac{\partial\lambda_i}{\partial q_k} = \frac{1}{\delta}\big(A^{-1}\big)_{ik}, \qquad \text{특히}\quad \frac{\partial\lambda_k}{\partial q_k} = \frac{1}{\delta}\big(A^{-1}\big)_{kk} $$

$A^{-1}$은 양정치 행렬이므로 그 대각원소 $(A^{-1})_{kk}$는 항상 양수다. 따라서 $\lambda_k$는 $q_k$의 증가함수다. (동시에 $i\ne k$인 경우의 교차편미분 $(A^{-1})_{ik}/\delta$가 0이 아닐 수 있다는 점도 이 식에서 바로 드러난다 — 본문 4장에서 언급한 "다른 뷰에 대한 파급 효과"의 근거다.)

**확신도 $\omega_k^{-1}$에 대한 증가성.** $A=\Omega/\tau+P\Sigma P'/(1+\tau)$ 이므로 $\partial A/\partial\omega_k^{-1} = -\omega_k^{2}\,\iota_{kk}/\tau$ (단, $\iota_{kk}$는 $(k,k)$ 원소만 1이고 나머지는 0인 행렬). 역행렬의 미분 공식 $\partial A^{-1}/\partial x = -A^{-1}(\partial A/\partial x)A^{-1}$을 쓰면

$$ \frac{\partial\Lambda}{\partial\omega_k^{-1}} = \big(\omega_k^{2}/\tau\big)\,A^{-1}\iota_{kk}\,\Lambda $$

이 벡터의 $k$번째 원소는

$$ \frac{\partial\lambda_k}{\partial\omega_k^{-1}} = \big(\omega_k^{2}/\tau\big)\big(A^{-1}\big)_{kk}\,\lambda_k $$

$(\omega_k^2/\tau)(A^{-1})_{kk} > 0$ 이므로, 이 미분값은 $\lambda_k$와 항상 같은 부호를 갖는다. 즉 $\lambda_k>0$이면 확신도가 높아질수록 $\lambda_k$는 더 커지고, $\lambda_k<0$이면 더 작아진다 — 어느 경우든 $|\lambda_k|$는 확신도 $\omega_k^{-1}$의 증가함수다. $\blacksquare$
