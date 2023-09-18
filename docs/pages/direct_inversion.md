## Direct amplitude inversion

Mean square voltage: $$\langle V^{2} \rangle = \left(\frac{K_{\rm s} K_{\rm t}}{r\psi_{(r)}}\right) C e^{-4[r\alpha_{\rm w} + \alpha_{\rm s}]}\text{, with }\alpha_{\rm s} = \int_{0}^{r}\xi C$$

### 1. Explicit solution

Lets define : $$f(r) = \langle V^{2} \rangle \left(\frac{r\psi_{(r)}}{K_{\rm s} K_{\rm t}}\right) e^{4r\alpha_{\rm w}} = Ce^{-4\alpha_{\rm s}}.$$

Then after logarithmic differentiation: $$\frac{f'}{f}  = \frac{C'}{C} - 4\xi C,$$

Under the form $$C' - \frac{f'}{f}C = 4\xi C,$$
we define $u = 1/C$, $u' = -C'/C^{2}$, which leads to:
$$u' + \frac{f'}{f}u = -4\xi.$$

A solution can be written as:

$$u(r) = \frac{1}{f}\left(-4\xi\int^{r} f + B \right),$$ where $B$ is a constant to be determined, which leads to

$$C(r) = \frac{f}{B - 4\xi\displaystyle\int^{r}f} $$

### 2. Constant determination

#### A. with $C(r_{0}) = C_{0}$ known

It leads to: $$B = \frac{f(r_{0})}{C_{0}} + 4\xi\int^{r_{0}}f.$$

#### B. with $\displaystyle\frac{1}{2\delta r}\int_{r_{0} - \delta r}^{r_{0} + \delta r}C = C_{0}$ known

We use $$\frac{1}{2\delta r}\int_{r_{0} - \delta r}^{r_{0} + \delta r} \frac{f}{B - 4\xi\int^{r}f} = C_{0}.$$

By using $u = \int^{r}f$, and $u' = f$, we integrate the previous equation, leading to:

$$\frac{-1}{2\delta r}\frac{1}{4\xi}\ln\left[\frac{B - 4\xi u_{r_{0} + \delta r}}{B - 4\xi u_{r_{0} - \delta r}}\right] = C_{0}$$

Then:

$$B = 4\xi\frac{u_{r_{0} + \delta r} - e^{-8\delta r\xi C_{0}}u_{r_{0} - \delta r}}{1 - e^{-8\delta r\xi C_{0}}}.$$

Note that when $\delta r \to 0$, we recover the solution for a point-imposed concentration $C(r_{0}) = C_{0}$.