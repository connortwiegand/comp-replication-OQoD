---
title : Optimum Debt replication
author : Evan Hart
date: June 2021
---
```julia
# using Weave; weave("/Users/modelt/Documents/Julia/Julia Spring 2021/Computation/optimum debt/outline2.jmd",doctype="md2html")
```
# Summary
Using an augmented Aigari model calibrated to the US economy the authors calculate the optimum quantity of debt.  This optimum of debt/GDP ratio is equal to 2/3. 

# Baby Model 
Version of the 1994a Aigari model, augmented to permit growth and to include government debt, lump sum taxes and government consumption. The consumer chooses stochastic processes for consumption and asset holdings in order to maximize the expected discounted sum of utilities of consumption subject to the sequence of budget constraints and nonnegativity constraints, i.e.

$$max_{c_t,a_{t+1}} E[\sum_{t=0}^\infty \beta^t c_t^{1-v}/(1-v) |a_0,e_0]$$\
subject to
$$c_t +a_{t+1} \leq (1+r)a_t +w_te_t-T_t$$
$$c_t\geq0, a_t\geq0,t\geq0.$$




# 3.1 Technical Appendix: Benchmark Model 
## Maximization problem
The consumer chooses sequences of consumption, asset holdings, and leisure to maximize expected utility, presented in stationary form.
$$max_{\tilde{c_t},l_t,\tilde{a}_{t+1}} E[ \sum_{t=0}^\infty {[\beta(1+g)^{\eta(1-\mu)}]}^t {(\tilde{c_t}^\eta l_t^{1-\eta})}^{1-\mu}/(1-\mu) |\tilde{a_0},e_0]$$
subject to\
$$\tilde{c_t}+ (1+g)\tilde{a}_{t+1} \leq (1+\bar{r})\tilde{a}_t +\bar{w}_te_t(1-l_t)+\chi$$
$$\tilde{a}_t\geq0,$$
$$l_t \leq1$$


## FOC
Because the decision rules are stationary, we can compute the functions $c(x, i), \alpha(x, i)$, and $l(x, i)$ that satisfy the following first-order conditions:

$$\eta(1+g)c(x,i)^{\eta(1-\mu)-1}l(x,i)^{(1-\eta)(1-\mu)}= \beta(1+g)^{\eta(1-\mu)}\{\sum_j\pi_{i,j}\eta(1+\bar{r})c(x',j)^{\eta(1-\mu)-1}l(x',j)^{(1-\eta)(1-\mu)}+\zeta min(\alpha(x,i),0)^2\}$$ (1)

$$(1-\eta)c(x,i)^{\eta(1-\mu)}l(x,i)^{(1-\eta)(1-\mu)-1} =\eta\bar{w}e(i)c(x,i)^{\eta(1-\mu)-1}l(x,i)^{(1-\eta)(1-\mu)}+\zeta min(1-l(x,i),0)^2$$ (2)

$$c(x,i)+(1+g)\alpha(x,i)= (1+\bar{r})x+\bar{w}e(i)(1-l(x,i))+\chi$$ (3)

$$x'=\alpha(x,i)$$ (4)


Notice that $c(x, i)$ is a function of $x, i, e(i)$, the parameters, and the function $\alpha(·)$. Therefore, with a candidate solution for $\alpha,$ we can back out $c$ using (3). This is not the case for leisure, however, since (2) is an implicit function of $l(x, i)$ (if we assume that $x, i, e(i)$, the parameters, $\alpha(·)$, and $c(·)$ are known)

## Numerical Solution for $l(x,i)$
We can construct a numerical solution for $l(x,i)$ by applying a Newton method with given values for $\eta, \mu, \zeta , \bar{r}, \bar{w}, \chi, g, x, i, e(i)$, and $\alpha(x, i)$. We simply iterate on

$$l^{k+1}=l^k-f(l^k)/j(l^k),\hspace{10mm} k=0,1,...$$ (5)
where
$$f(l)=(1-\eta)c(l)^{\eta(1-\mu)}l^{(1-\eta)(1-\mu)-1}-\zeta min(1-l(x,i),0)^2-\bar{w}e(i)\eta c(l)^{\eta(1-\mu)-1}l^{(1-\eta)(1-\mu)}$$ 

$$c(l)=(1+\bar{r})x+\bar{w}e(i)(1-l)+\chi-(1+g)\alpha(x,i)$$ 


$$j(l) = -2\bar{w}e(i)(1-\eta)\eta(1-\mu)c(l)^{\eta(1-\mu)-1}l^{(1-\eta)(1-\mu)-1}$$
$$+(1-\eta)((1-\eta)(1-\mu)-1)c(l)^{\eta(1-\mu)}l^{(1-\eta)(1-\mu)-2}$$
$$+\bar{w}^2e(i)^2\eta(\eta(1-\mu)-1)c(l)^{\eta(1-\mu)-2}l^{(1-\eta)(1-\mu)}$$
$$+2\zeta min(1-l,0)$$

The function f is the Euler equation in (2), the function c is consumption derived from the budget constraint, and the function J is the derivative of f with respect to l. We start the iterations in (5) with an initial guess $l_0$, and we stop when $|l_{k+1} − l_k|$ is less than some tolerance parameter.





## Residual
Let $l^*(x,i,\alpha)$ be the solution to (5). Then we can write the first-order conditions in (1)-(3) in terms of one residual;
$$R(x,i;\alpha)= \eta(1+g)c(x,i)^{\eta(1-\mu)-1}l^*(x,i)^{(1-\eta)(1-\mu)}$$
$$\hspace{30mm}-\beta(1+g)^{\eta(1-\mu)}\{\sum_j\pi_{i,j}\eta(1+\bar{r})c(x',j)^{\eta(1-\mu)-1}$$
$$\hspace{20mm}\cdot l^*(x',j)^{(1-\eta)(1-\mu)}+\zeta min(\alpha(x,i),0)^2\}$$


Then we can write the first-order conditions in (1)-(3) in terms of one residual;
$$R(x,i;\alpha)= \eta(1+g)c(x,i)^{\eta(1-\mu)-1}l^*(x,i)^{(1-\eta)(1-\mu)}$$
$$\hspace{30mm}-\beta(1+g)^{\eta(1-\mu)}\{\sum_j\pi_{i,j}\eta(1+\bar{r})c(x',j)^{\eta(1-\mu)-1}$$
$$\hspace{20mm}\cdot l^*(x',j)^{(1-\eta)(1-\mu)}+\zeta min(\alpha(x,i),0)^2\}$$


If we ignore the constraint on leisure, the residual is given by
$$R(x,i;\alpha)= \eta(1+g)\{\frac{1-\eta}{\eta\bar{w}e(i)}\}^{(1-\eta)(1-\mu)}$$
$$\hspace{30mm}-\beta(1+g)^{\eta(1-\mu)}\{\sum_j\pi_{i,j}\eta(1+\bar{r})c(x',j)^{\eta(1-\mu)-1}$$
$$\hspace{20mm}\cdot l^*(x',j)^{(1-\eta)(1-\mu)}+\zeta min(\alpha(x,i),0)^2\}$$

```julia
function R(BenchmarkParameters, e, r)
@unpack g,χ,γ,b,θ = BenchmarkParameters
Ψ_e = 0.2
Ψ_e2 = 0.1
x = 2
xvec = [1 , 2]

t_y = (γ+ χ+(r_bar-g)*b)/(1-δ*k) # pg 458 eq (15) check t against equation (14) in appendix
wr = 1-θ# pg 461
kr = θ/(r+δ)# pg 461
A = (1-η)/[ η* (1-ty)*(1- θ) ] # pg 461
B = 1 - γ- θ*(g+ δ)/(r+ δ) # pg 461
N = 1/(1+A*B) # pg 461
r_bar = (1-t_y)*r # eq (12) from main paper
w_bar = (1-t_y)*wr/N# eq (16) from main paper
# α_bar = kr+b
sum = sum(πi .*  (1+r_bar)*(1+r_bar) *alpha_approximation(Ψ_e,Ψ_e2,x,xvec) * w_bar*e 
.+  χ - (1+g) *alpha_approximation(Ψ_e,Ψ_e2,alpha_approximation(Ψ_e,Ψ_e2,x,xvec),xvec))

R = (1+g) * [ (1+r_bar)*x + w_bar* e[i] + χ - (1+g)* alpha_approximation(Ψ_e,Ψ_e2,x,xvec) ]^(-υ)   
- β*(1+g)^(1-υ)* (sum^(-υ)  + ζ* minimum( alpha_approximation(Ψ_e,Ψ_e2,x,xvec),0))^2  
end


```

$\int R(x,i;\alpha^h) *N_e(x) dx = 0$