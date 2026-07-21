---
layout: post
title: "[스펙트럴 그래프 이론 #2] 고유값은 어떻게 최적화 문제의 답이 되는가 — Courant-Fischer 정리"
date: 2026-07-21
categories: [math, graph-theory]
tags: [spectral-graph-theory, linear-algebra, courant-fischer]
math: true
---

> Daniel Spielman의 *Spectral and Algebraic Graph Theory* 강의노트를 공부하고 제 방식대로 재구성한 글입니다. 정확한 원문은 저자의 공개 강의노트를 참고해주세요.

지난 글에서 라플라시안의 고유값이 왜 의미를 갖는지 감을 잡았다면, 이번 글에서는 그 "의미"를 엄밀하게 만들어주는 정리를 다룹니다. 바로 **Courant-Fischer 정리**입니다. 한 줄로 요약하면: *대칭행렬의 고유값들은 특정한 최적화 문제의 답으로 정확히 표현된다.*

## 레일리 지수 (Rayleigh Quotient)

행렬 $M$과 벡터 $x$에 대해 **레일리 지수**를 다음과 같이 정의합니다.

$$\frac{x^T M x}{x^T x}$$

만약 $x=\psi$가 고유값 $\mu$의 고유벡터라면

$$\frac{\psi^T M \psi}{\psi^T \psi} = \frac{\psi^T \mu \psi}{\psi^T \psi} = \mu$$

즉 고유벡터의 레일리 지수는 정확히 그 고유값입니다. Courant-Fischer 정리는 여기서 한 걸음 더 나아가, **레일리 지수를 최대화(혹은 최소화)하는 벡터가 바로 고유벡터**임을 말해줍니다.

아래 그림은 $2\times 2$ 대칭행렬에 대해, 단위원 위의 점 $x=(\cos\theta,\sin\theta)$마다 $x^TMx$ 값이 어떻게 변하는지 보여줍니다. 값이 최대가 되는 방향과 최소가 되는 방향이 정확히 두 고유벡터의 방향과 일치하는 것을 볼 수 있습니다.

![레일리 지수의 기하학적 의미](assets/img/posts/라플라시안-행렬-ch2-courant-fischer/ch2_rayleigh_quotient.png)

## Courant-Fischer 정리

**정리.** 대칭행렬 $M$의 고유값을 $\mu_1 \ge \mu_2 \ge \cdots \ge \mu_n$이라 하면,

$$\mu_1 = \max_{x\ne 0} \frac{x^TMx}{x^Tx}, \qquad \mu_n = \min_{x\ne 0} \frac{x^TMx}{x^Tx}$$

이고, 일반적으로 $k$번째 고유값은

$$\mu_k = \max_{\substack{S\subseteq\mathbb{R}^n \\ \dim(S)=k}} \ \min_{\substack{x\in S \\ x\ne 0}} \frac{x^TMx}{x^Tx} = \min_{\substack{T\subseteq\mathbb{R}^n \\ \dim(T)=n-k+1}} \ \max_{\substack{x\in T \\ x\ne 0}} \frac{x^TMx}{x^Tx}$$

로 표현됩니다.

말로 풀면 이렇습니다. 가장 큰 고유값은 "전체 공간에서 레일리 지수를 최대화한 값"이고, 가장 작은 고유값은 "전체 공간에서 최소화한 값"입니다. 중간 고유값 $\mu_k$는 좀 더 미묘한데, "$k$차원 부분공간을 잘 고르면 그 안에서의 최솟값을 최대한 크게 만들 수 있고, 그 최댓값이 바로 $\mu_k$"라는 뜻입니다 (혹은 대칭적으로, $(n-k+1)$차원 부분공간에서 최댓값을 최소화).

이 정리가 왜 유용할까요? 대표적으로, 라플라시안의 $\lambda_2$(Fiedler 값)에 상한을 걸고 싶을 때 굳이 최적의 부분공간을 찾을 필요 없이, **아무 벡터나 하나 골라서** ($1$에 직교하도록) 그 레일리 지수를 계산하기만 하면 $\lambda_2$의 상한을 즉시 얻습니다. 이런 벡터를 **테스트 벡터(test vector)**라 부르며, 5장에서 본격적으로 활용됩니다.

## 증명의 핵심 아이디어

증명은 어떤 벡터 $x$든 고유벡터들의 정규직교기저로 전개할 수 있다는 사실에서 출발합니다.

$$x = \sum_i c_i \psi_i, \qquad c_i = \psi_i^T x$$

이 전개를 이차형식에 대입하면 다음과 같은 깔끔한 보조정리를 얻습니다.

**보조정리.** $x = \sum_i c_i \psi_i$이면 $x^TMx = \sum_i c_i^2 \mu_i$이다.

즉 이차형식은 각 고유공간 방향의 성분 제곱에 그 방향의 고유값을 곱해 더한 것과 같습니다. 이 사실 하나로 정리의 대부분이 자연스럽게 풀립니다 — $\mu_k$ 이상인 성분에 가중치를 몰아주면 값이 커지고, 반대로 하면 작아지기 때문입니다. 완전한 증명은 부록 A에 정리했습니다.

## Courant-Fischer로 스펙트럴 정리를 거꾸로 증명하기

흥미로운 점은, 앞선 글에서 사용했던 스펙트럴 정리(대칭행렬은 항상 직교하는 고유벡터들을 갖는다)를 **거꾸로 Courant-Fischer의 특수한 경우(가장 큰/작은 고유값에 대한 부분)만으로부터 유도**할 수 있다는 것입니다. 순서를 다음과 같이 뒤집는 셈입니다.

1. 레일리 지수를 최대화하는 단위벡터 $x$가 (콤팩트집합 위 연속함수이므로) 반드시 존재한다는 해석학적 사실만 받아들인다.
2. 그 최댓값 지점에서 그래디언트가 0이라는 미적분 조건으로부터, $x$가 실제로 $M$의 고유벡터임을 보인다.
3. $M$에서 이 고유벡터의 기여분을 "빼낸" 행렬 $\hat M = M - \mu\psi\psi^T$을 만들고, 그 랭크가 정확히 1 작아짐을 보인 뒤 귀납법을 적용한다.

이렇게 하면 대칭행렬이라면 무조건 정규직교 고유벡터 집합을 가진다는 사실이 순수하게 최적화 논증만으로 증명됩니다. 자세한 단계는 부록 B에 있습니다.

## 비대칭 행렬의 경우: 특이값 분해

여기까지 다룬 레일리 지수 논증은 **대칭행렬**에서만 성립합니다. 정사각형이 아니거나 대칭이 아닌 행렬 $A$에 대응하는 개념은 **특이값(singular value)**입니다.

$$A = U\Sigma V^T$$

로 분해할 때, $\Sigma$의 대각원소들이 특이값이고 $U, V$의 열들이 각각 왼쪽/오른쪽 특이벡터입니다. 특이값은 $AA^T$ 혹은 $A^TA$의 고유값의 제곱근으로 정의되며, Courant-Fischer와 똑같은 형태의 최적화적 특성화가 성립합니다:

$$\sigma_1 = \max_{\|u\|=1,\|v\|=1} u^TAv$$

이 개념은 4장에서 이분그래프(bipartite graph)의 인접행렬을 다룰 때 다시 등장합니다.

---

## 부록: 증명 모음

### 부록 A. Courant-Fischer 정리의 증명

**보조정리 2.1.1.** $M$이 대칭행렬이고 고유값 $\mu_1,\dots,\mu_n$에 대응하는 정규직교 고유벡터 $\psi_1,\dots,\psi_n$이 있다고 하자. $x=\sum_i c_i\psi_i$이면

$$x^TMx = \sum_i c_i^2\mu_i$$

*증명.* $M\psi_j = \mu_j\psi_j$이고 $\psi_i^T\psi_j$는 $i=j$일 때 1, 아니면 0이므로

$$x^TMx = \Big(\sum_i c_i\psi_i\Big)^T M\Big(\sum_j c_j\psi_j\Big) = \sum_{i,j} c_ic_j\mu_j\,\psi_i^T\psi_j = \sum_i c_i^2\mu_i \qquad\blacksquare$$

**정리의 증명 (일반적인 $\mu_k$에 대해).**

*(도달 가능성)* $S = \text{span}(\psi_1,\dots,\psi_k)$로 두면, $S$ 위의 임의의 $x=\sum_{i\le k}c_i\psi_i$에 대해

$$\frac{x^TMx}{x^Tx} = \frac{\sum_{i\le k}\mu_i c_i^2}{\sum_{i\le k}c_i^2} \ge \mu_k \cdot \frac{\sum_{i\le k}c_i^2}{\sum_{i\le k}c_i^2} = \mu_k$$

(각 $\mu_i \ge \mu_k$이기 때문) 이므로 $\min_{x\in S} \frac{x^TMx}{x^Tx} \ge \mu_k$가 성립하고, $x=\psi_k$일 때 등호가 성립하므로 이 최솟값은 정확히 $\mu_k$입니다.

*(최적성)* 임의의 $k$차원 부분공간 $S$에 대해 $\min_{x\in S}\frac{x^TMx}{x^Tx} \le \mu_k$임을 보이면 됩니다. $T=\text{span}(\psi_k,\dots,\psi_n)$은 차원이 $n-k+1$이므로, 임의의 $k$차원 $S$와는 반드시 최소 1차원의 교집합을 가집니다 ($k + (n-k+1) > n$이기 때문). 이 교집합에서 벡터 하나를 고르면

$$\min_{x\in S}\frac{x^TMx}{x^Tx} \le \min_{x\in S\cap T}\frac{x^TMx}{x^Tx} \le \max_{x\in T}\frac{x^TMx}{x^Tx}$$

이제 $x=\sum_{i\ge k} c_i\psi_i \in T$에 대해 (같은 보조정리로)

$$\frac{x^TMx}{x^Tx} = \frac{\sum_{i\ge k}\mu_ic_i^2}{\sum_{i\ge k}c_i^2} \le \mu_k$$

이므로 $\max_{x\in T}\frac{x^TMx}{x^Tx} \le \mu_k$. 두 부등식을 결합하면 원하는 결과를 얻습니다. $\blacksquare$

### 부록 B. Courant-Fischer의 특수한 경우로부터 스펙트럴 정리 증명하기

**1단계 (극값의 존재와 고유벡터 성질).** $x^TMx$를 단위벡터 위에서 최대화하는 벡터 $x$가 극값이라면, 함수 $f(x) = \frac{x^TMx}{x^Tx}$의 그래디언트가 0이어야 합니다. $\nabla(x^Tx) = 2x$, $\nabla(x^TMx)=2Mx$이므로 몫의 미분법으로

$$\nabla f(x) = \frac{2(x^Tx)Mx - 2(x^TMx)x}{(x^Tx)^2} = 0 \ \Longrightarrow\ Mx = \frac{x^TMx}{x^Tx}\,x$$

즉 극값을 이루는 $x$는 자기 자신의 레일리 지수를 고유값으로 갖는 고유벡터입니다. 최댓값에서 이 값은 가장 큰 고유값 $\mu_1$이 됩니다.

**2단계 (귀납적으로 나머지 고유벡터 뽑아내기).** $\psi$가 고유값 $\mu\ne 0$인 단위고유벡터일 때 $\hat M = M-\mu\psi\psi^T$를 정의하면:

- $\hat M \psi = M\psi - \mu\psi(\psi^T\psi) = \mu\psi-\mu\psi = 0$ (즉 $\psi$는 $\hat M$의 영공간에 속함)
- $M$의 영공간에 속하는 임의의 $x$는 $\psi$와 직교하므로 ($M$의 상(span)과 영공간은 서로 직교) $\hat Mx = Mx - \mu\psi(\psi^Tx) = 0$
- 대칭행렬의 상은 영공간과 직교하므로, $M$의 상은 $\hat M$의 상과 $\psi$가 张(span)하는 1차원 공간의 합이 됨 → $\text{rank}(\hat M) = \text{rank}(M) - 1$

이 사실로 랭크에 대한 귀납법을 돌리면: $\hat M$에 대해 귀납가설을 적용해 정규직교 고유벡터 $\psi_1,\dots,\psi_r$과 0이 아닌 고유값 $\mu_1,\dots,\mu_r$을 얻고, 여기에 $\psi_{r+1}=\psi$, $\mu_{r+1}=\mu$를 추가하면

$$M = \sum_{i=1}^{r+1}\mu_i\psi_i\psi_i^T$$

를 얻습니다. $\psi_{r+1}$이 나머지와 직교함은, $\psi_i$ ($i\le r$)가 $\hat M$의 상에 있고 $\psi_{r+1}$은 $\hat M$의 영공간에 있다는 사실에서 바로 나옵니다.

**남은 조각: 0이 아닌 고유벡터가 항상 존재하는가?** $M\ne 0$인 대칭행렬이라면, 어떤 벡터 $x$에 대해 $x^TMx\ne 0$임을 보일 수 있습니다 (대각원소가 0이 아니면 표준기저벡터를, 모두 0이면 0이 아닌 비대각원소 위치의 두 기저벡터의 합을 사용). 이 $x$에서 시작해 레일리 지수를 최대화(혹은 부호가 반대면 $-M$에 대해 최대화)하면 0이 아닌 고유값의 고유벡터를 얻습니다. 이렇게 랭크에 대한 귀납의 기초 단계가 성립합니다. $\blacksquare$
