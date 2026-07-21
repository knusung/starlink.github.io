---
layout: post
title: "[스펙트럴 그래프 이론 #5] 기본 그래프들의 고유값을 직접 계산해보기"
date: 2026-07-21
categories: [math, graph-theory]
tags: [spectral-graph-theory, laplacian, eigenvalues]
math: true
---

> Daniel Spielman의 *Spectral and Algebraic Graph Theory* 강의노트를 공부하며 제 나름대로 재구성한 글입니다.

이론을 아무리 익혀도, 실제로 손에 잡히는 예제 몇 개를 직접 계산해보지 않으면 감이 잘 안 잡힙니다. 이번 글에서는 완전그래프, 별그래프, 그래프의 곱(격자·하이퍼큐브), 링그래프, 경로그래프의 라플라시안 고유값을 하나씩 구해봅니다. 이 그래프들은 모두 연결그래프이므로 고유값 0은 항상 중복도 1입니다.

![기본 그래프들](assets/img/posts/라플라시안-행렬-ch5-fundamental-graphs/ch5_fundamental_graphs.png)

## 등주비율(isoperimetric ratio) 복습

계산에 들어가기 전에, 뒤에서 계속 등장할 개념 하나를 미리 정의해둡니다. 정점 집합 $S$의 **경계** $\partial(S)$는 정확히 한쪽 끝만 $S$에 속하는 간선들의 집합이고, **등주비율**은

$$\theta(S) \stackrel{\text{def}}{=} \frac{|\partial(S)|}{|S|}$$

입니다. 20장에서 증명되는 부등식 $\theta(S)\ge\lambda_2(1-s)$ (단 $s=|S|/n$)가 이번 글 곳곳에서 기준점 역할을 합니다.

## 완전그래프 $K_n$

모든 정점 쌍이 간선으로 이어진 그래프입니다.

**고유값: $0$(중복도 1), $n$(중복도 $n-1$).**

가장 쉬운 방법은 $L_{K_n} = nI - \mathbf1\mathbf1^T$라는 관찰입니다. 상수벡터에 직교하는 임의의 벡터 $x$에 대해 $\mathbf1^Tx=0$이므로 $L_{K_n}x = nx - \mathbf1(\mathbf1^Tx) = nx$가 바로 나옵니다.

완전그래프는 등주부등식이 **꽉 차는(tight)** 경우이기도 합니다. $S\subset[n]$이면 $\theta(S) = n-|S| = \lambda_2(1-s)$가 정확히 성립하는데, 이는 $S$ 안의 모든 정점이 $S$ 밖의 모든 정점과 연결되어 있어 경계가 최대한 "빽빽"하기 때문입니다.

## 별그래프 $S_n$

정점 1이 중심이고, 나머지 $2,\dots,n$이 모두 정점 1에만 연결된 그래프입니다.

**고유값: $0$(중복도 1), $1$(중복도 $n-2$), $n$(중복도 1).**

핵심 관찰은 다음 보조정리입니다.

> 차수가 1인 두 정점 $a,b$가 같은 이웃을 공유하면, $\delta_a-\delta_b$는 항상 고유값 1의 고유벡터다.

별그래프의 잎(leaf) 정점들은 모두 중심을 공유 이웃으로 가지므로, 임의의 두 잎 $i,i+1$에 대해 $\delta_i-\delta_{i+1}$이 고유값 1의 고유벡터가 되고, 이런 식으로 서로 독립인 $n-2$개의 고유벡터를 얻습니다. 나머지 고유값은 트레이스(대각합 = 고유값의 합)를 이용해서 구합니다: $\text{tr}(L_{S_n}) = 2n-2$이고 이미 찾은 $n-1$개 고유값의 합이 $n-2$이므로, 남은 고유값은 $n$이어야 합니다.

## 그래프의 곱: 격자와 하이퍼큐브

두 그래프 $G=(V,E,v), H=(W,F,w)$의 **곱** $G\times H$는 정점집합이 $V\times W$이고, 한쪽 좌표가 고정된 채 다른 쪽 좌표만 $G$ 또는 $H$의 간선을 따라 바뀌는 그래프입니다. 경로 두 개를 곱하면 격자가, 한 개의 간선을 반복해서 곱하면 하이퍼큐브가 나옵니다.

**정리.** $G,H$의 라플라시안 고유값/고유벡터가 각각 $(\lambda_i,\alpha_i)$, $(\mu_j,\beta_j)$라면, $G\times H$는 고유값 $\lambda_i+\mu_j$를 갖는 고유벡터 $\gamma_{i,j}(a,b)=\alpha_i(a)\beta_j(b)$를 갖습니다.

즉 **곱 그래프의 고유값은 두 원래 그래프의 고유값의 모든 쌍별 합**이고, 고유벡터는 두 원래 고유벡터의 곱(외적)입니다. 증명은 라플라시안을 연산자로 직접 적용해서 항을 정리하면 나옵니다 (부록 A). 크로네커 곱을 써서 $L_{G\times H} = L_G\otimes I_W + I_V\otimes L_H$로 표현할 수도 있습니다.

### 하이퍼큐브 $H_d$

한 개의 간선 그래프 $H_1$을 자기 자신과 $d-1$번 곱하면 $d$차원 하이퍼큐브 $H_d$가 나옵니다. $H_1$의 고유값은 $0, 2$이므로, 곱 정리를 반복 적용하면:

- $H_d$는 고유값 $2i$ ($0\le i\le d$)를 중복도 $\binom{d}{i}$로 가집니다.
- 각 고유벡터는 $y\in\{0,1\}^d$로 인덱싱되며 $\psi_y(x) = (-1)^{y^Tx}$의 형태를 갖습니다.

$\lambda_2(H_d)=2$라는 사실과 20장의 등주부등식을 결합하면 다음 따름정리가 즉시 나옵니다.

**따름정리.** 하이퍼큐브에서 $\theta_{H_d}\ge 1$, 즉 전체 정점의 절반 이하인 임의의 정점 집합에 대해, 경계 간선의 수는 그 집합의 크기보다 작지 않다.

이 부등식은 정점 전체의 절반, 즉 "한쪽 면"(예: 라벨이 0으로 시작하는 모든 정점)을 잡으면 등호가 성립할 정도로 꽉 찬 결과입니다.

## 테스트 벡터로 $\lambda_2$의 상한 구하기

Courant-Fischer 정리에 의해, $1$에 직교하는 **아무 벡터나** 하나 넣으면 $\lambda_2$의 상한을 얻습니다 ($\lambda_2 \le v^TLv/v^Tv$). 이런 용도로 쓰는 벡터를 **테스트 벡터**라 부릅니다.

경로그래프 $P_n$에서, "정점 번호 자체"를 테스트 벡터로 쓰고 싶지만 그건 $1$에 직교하지 않으므로, 대신 $x(a) = (n+1)-2a$를 사용합니다. 이 벡터는 정점 번호에 선형이면서 $1^Tx=0$을 만족합니다. 계산하면

$$\lambda_2(P_n) \le \frac{\sum_{a<n}(x(a)-x(a+1))^2}{\sum_a x(a)^2} = \frac{4(n-1)}{(n+1)n(n-1)/3} = \frac{12}{n(n+1)}$$

아래 그림은 이 상한과 실제 $\lambda_2$ 값을 $n$을 늘려가며 비교한 것입니다. 상한이 실제 값과 같은 $O(1/n^2)$ 오더를 정확히 잡아내는 것을 볼 수 있습니다.

![경로그래프의 λ2: 실제값과 테스트 벡터 상한](assets/img/posts/라플라시안-행렬-ch5-fundamental-graphs/ch5_path_lambda2_bound.png)

흥미로운 점 하나: 20장의 등주부등식 (5.1)을 경로그래프에 적용하면 훨씬 느슨한 결과가 나옵니다. $S=\{1,\dots,n/2\}$에서 $\theta(S)=2/n$인데, 부등식이 주는 하한은 $\lambda_2(1-s)$ 형태로 $O(1/n^2)$ 오더라서, 실제 등주비율 $2/n$과는 자릿수 하나 정도 차이가 납니다. 21장에서 다루는 **Cheeger 부등식**이 바로 이 간극이 "제곱보다 나쁘지는 않다"는 것을 보장해줍니다.

한편 Courant-Fischer의 다른 절반(최댓값-최솟값 형태)은 $\lambda_2$의 **하한**을 증명할 때는 별 도움이 안 됩니다. $n-1$차원 부분공간 위에서의 최솟값의 최댓값이라는 형태는 다루기가 훨씬 어렵기 때문입니다. 이런 하한은 대개 다른 기법(비교 논증 등, 6장에서 다룸)으로 증명합니다.

## 링그래프 $R_n$

정점을 $0,\dots,n-1$로 보고 $\pmod n$으로 인접한 링(사이클) 그래프입니다.

**정리.** $R_n$의 라플라시안은 다음 고유벡터를 갖는다.

$$x_k(a)=\cos(2\pi ka/n),\ 0\le k\le n/2, \qquad y_k(a)=\sin(2\pi ka/n),\ 1\le k<n/2$$

각각 고유값 $2-2\cos(2\pi k/n)$을 갖습니다.

이 결과를 눈으로 확인하는 좋은 방법은, 각 정점을 단위원 위의 점 $(\cos(2\pi a/n),\sin(2\pi a/n))$에 그려보는 것입니다. 한 정점의 두 이웃의 평균은 그 정점 자신과 같은 방향을 가리키는 벡터가 되는데, 이는 곧 원 위의 $x,y$좌표 자체가 (같은 고유값의) 고유벡터임을 시각적으로 보여줍니다. 대수적으로는 삼각함수의 덧셈정리를 이용해 직접 검증할 수 있습니다 (부록 B).

## 경로그래프 $P_n$을 링그래프로부터 유도하기

마지막으로, 경로그래프의 정확한 고유값·고유벡터를 (테스트 벡터의 근사가 아니라) 완전히 계산해봅니다. 트릭은 $P_n$을 $R_{2n}$의 "몫(quotient)"으로 보는 것입니다.

링 $R_{2n}$의 정점을 그림처럼 재배열하면(정점 $a$가 정점 $n+a$의 바로 위에 오도록), $P_n$의 각 정점 $a$를 $R_{2n}$의 정점 $a$와 $a+n$을 "접어서 합친" 것으로 볼 수 있습니다. 이때 $\psi(a)=\psi(a+n)$을 만족하는 $R_{2n}$의 고유벡터들이 $P_n$의 고유벡터로 정확히 대응됩니다.

**정리.** $P_n$의 라플라시안은 고유값 $2(1-\cos(\pi k/n))$, $0\le k<n$을 가지며, 대응 고유벡터는

$$v_k(a) = \cos\!\left(\frac{\pi ka}{n}-\frac{\pi k}{2n}\right)$$

이 유도는 $R_{2n}$의 고유벡터 쌍 $x_k,y_k$를 각도 $\theta=-\pi k/(2n)$만큼 회전시켜서, "접었을 때" 값이 일치하는 벡터를 만드는 방식으로 진행됩니다. 자세한 계산은 부록 C에 정리했습니다. 이런 식으로 큰 그래프의 대칭성을 이용해 작은 그래프의 스펙트럼을 유도하는 기법을 **등가분할(equitable partition)**이라 부르며, Godsil의 저서에 훨씬 자세한 이론이 있다고 원저자는 언급합니다.

---

## 부록: 증명 모음

### 부록 A. 그래프 곱의 고유값 정리 증명

$\alpha$가 $L_G$의 고유값 $\lambda$의 고유벡터, $\beta$가 $L_H$의 고유값 $\mu$의 고유벡터일 때, $\gamma(a,b)=\alpha(a)\beta(b)$가 $L_{G\times H}$의 고유값 $\lambda+\mu$의 고유벡터임을 보입니다.

$$
\begin{aligned}
(L_{G\times H}\gamma)(a,b) &= \sum_{(a,\hat a)\in E} v_{a,\hat a}(\gamma(a,b)-\gamma(\hat a,b)) + \sum_{(b,\hat b)\in F} w_{b,\hat b}(\gamma(a,b)-\gamma(a,\hat b))\\
&= \sum_{(a,\hat a)\in E} v_{a,\hat a}\beta(b)(\alpha(a)-\alpha(\hat a)) + \sum_{(b,\hat b)\in F} w_{b,\hat b}\alpha(a)(\beta(b)-\beta(\hat b))\\
&= \beta(b)\cdot\lambda\alpha(a) + \alpha(a)\cdot\mu\beta(b)\\
&= (\lambda+\mu)\,\alpha(a)\beta(b)
\end{aligned}
$$

세 번째 줄에서 $\sum_{(a,\hat a)\in E}v_{a,\hat a}(\alpha(a)-\alpha(\hat a)) = (L_G\alpha)(a) = \lambda\alpha(a)$라는 라플라시안 연산자의 정의를 그대로 사용했습니다. $\blacksquare$

### 부록 B. 링그래프 고유벡터의 대수적 검증

삼각함수 덧셈정리 $\cos(a+b)=\cos a\cos b - \sin a\sin b$를 이용해 직접 계산합니다.

$$
\begin{aligned}
(L_{R_n}x_k)(a) &= 2x_k(a) - x_k(a+1) - x_k(a-1)\\
&= 2\cos(2\pi ka/n) - \cos(2\pi ka/n)\cos(2\pi k/n) + \sin(2\pi ka/n)\sin(2\pi k/n)\\
&\quad -\cos(2\pi ka/n)\cos(2\pi k/n) - \sin(2\pi ka/n)\sin(2\pi k/n)\\
&= 2\cos(2\pi ka/n) - 2\cos(2\pi ka/n)\cos(2\pi k/n)\\
&= (2-2\cos(2\pi k/n))\,x_k(a)
\end{aligned}
$$

$y_k$에 대한 계산도 정확히 같은 방식으로 진행됩니다. 두 고유벡터 쌍을 각도 $\theta$만큼 "회전"시켜도 여전히 같은 고유값의 고유벡터가 된다는 사실은

$$\cos(2\pi ka/n+\theta) = (\cos\theta)\,x_k(a) - (\sin\theta)\,y_k(a)$$

로 확인됩니다 — 우변이 $x_k,y_k$가 张(span)하는 2차원 고유공간 안에 있기 때문입니다. $\blacksquare$

### 부록 C. 경로그래프 고유값을 링그래프로부터 유도하기

$I_n$을 $n$차원 항등행렬이라 하면, 정점을 접는 연산이 다음 행렬 등식으로 표현됩니다.

$$\begin{pmatrix}I_n & I_n\end{pmatrix} L_{R_{2n}} \begin{pmatrix}I_n\\I_n\end{pmatrix} = 2L_{P_n}$$

$\psi(a)=\psi(a+n)$을 만족하는 $R_{2n}$의 고유값 $\lambda$짜리 고유벡터 $\psi$가 있으면, $\varphi(a)=\psi(a)$ ($1\le a\le n$)로 정의할 때 $\binom{I_n}{I_n}\varphi = \psi$이므로

$$\begin{pmatrix}I_n&I_n\end{pmatrix}L_{R_{2n}}\begin{pmatrix}I_n\\I_n\end{pmatrix}\varphi = \begin{pmatrix}I_n&I_n\end{pmatrix}L_{R_{2n}}\psi = \begin{pmatrix}I_n&I_n\end{pmatrix}\lambda\psi = 2\lambda\varphi$$

즉 $\varphi$는 $L_{P_n}$의 고유값 $\lambda$짜리 고유벡터입니다.

이제 이런 $\psi$를 실제로 찾습니다. 부록 B의 $x_k,y_k$를 각도 $\theta=-\pi k/(2n)$만큼 회전시킨 벡터

$$\psi = (\cos\theta)x_k - (\sin\theta)y_k$$

를 새 인덱싱 하에서 쓰면

$$\psi(a) = \begin{cases}\cos(2\pi ka/(2n) - \pi k/2n) & 1\le a\le n\\ \cos(2\pi k(2n+1-a)/(2n) - \pi k/2n) & n<a\le 2n\end{cases}$$

$\cos(-x)=\cos x$를 이용해 $\psi(a+n)=\psi(a)$를 직접 확인할 수 있습니다:

$$
\begin{aligned}
\cos\!\left(\frac{2\pi k(2n+1-a)}{2n}-\frac{\pi k}{2n}\right) &= \cos\!\left(2\pi - \frac{2\pi ka}{2n} + \frac{\pi k}{n} - \frac{\pi k}{2n}\right)\\
&= \cos\!\left(-\frac{2\pi ka}{2n}+\frac{\pi k}{2n}\right) = \cos\!\left(\frac{2\pi ka}{2n}-\frac{\pi k}{2n}\right)
\end{aligned}
$$

따라서 $\psi(a)=\psi(a+n)$이 확인되고, 대응하는 $P_n$의 고유값은 $R_{2n}$에서의 고유값 $2-2\cos(2\pi k/(2n)) = 2(1-\cos(\pi k/n))$과 정확히 같습니다. $1\le k<n$에 대해 이렇게 $n-1$개의 서로 다른 고유값을 얻고, 마지막 남은 고유값 0은 상수벡터에서 나옵니다. $\blacksquare$
