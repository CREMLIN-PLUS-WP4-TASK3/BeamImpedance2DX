<!-- pandoc -f markdown -t html -s  --mathjax=https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js --metadata title='Derivation' -o derivation.html derivation.md -->

Based on paper [Space charge and resistive wall impedance computation in the frequency domain using
the finite element method by Uwe Niedermayer, Oliver Boine-Frankenheim, and Herbert De Gersem](https://doi.org/10.1103/PhysRevSTAB.18.032001).
Equations in this paper will be referred to as `Uwe:1`, `Uwe:2` etc.

---

Our goal is to obtain longitudinal and transverse impedance (Uwe:4-6)

$$\underline{Z}_\parallel(\omega) = -\frac{1}{q^2}\int_{beam}{\vec{\underline{E}}\vec{\underline{J}}^{*}_\parallel\;dV}$$
$$\underline{Z}_\perp(\omega) = -\frac{\beta c}{(q d_r)^2\omega}\int_{beam}{\vec{\underline{E}}\vec{\underline{J}}^{*}_\perp\;dV}$$

*Underlined variables denote complex values.*

For this we need to solve Maxwell's equations in 2D.

$$
\begin{cases}
\nabla\cdot\vec{\underline{D}} = \rho \\
\nabla\cdot\vec{\underline{B}} = 0 \\
\nabla\times\vec{\underline{E}} = -j\omega\vec{\underline{B}} \\
\nabla\times\vec{\underline{H}} = J_s+j\omega\vec{\underline{D}} \\
\vec{\underline{D}} = [\underline{\varepsilon}]\vec{\underline{E}} \\
\vec{\underline{B}} = [\underline{\mu}]\vec{\underline{H}}
\end{cases}
$$

Simplify by setting $\rho$ to zero and reducing tensor values $[\underline{\varepsilon}]$ and $[\underline{\mu}]$ to scalars
$\underline{\varepsilon}$ and $\underline{\mu}$.

$$
\begin{cases}
\nabla\cdot\underline{\varepsilon}\vec{\underline{E}} = \rho \\
\nabla\cdot\underline{\mu}\vec{\underline{H}} = 0 \\
\nabla\times\vec{\underline{E}} = -j\omega\underline{\mu}\vec{\underline{H}} \\
\nabla\times\vec{\underline{H}} = J_s+j\omega\underline{\varepsilon}\vec{\underline{E}}
\end{cases}
$$

$$\vec{\underline{H}} = j\frac{1}{\omega\underline{\mu}}\nabla\times\vec{\underline{E}}$$
$$j\frac{1}{\omega}\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}} = J_s+j\omega\underline{\varepsilon}\vec{\underline{E}}$$
$$\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}} = -j\omega J_s+\omega^2\underline{\varepsilon}\vec{\underline{E}}$$

The governing wave equation is (Uwe:7)

$$\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}} - \omega^2\underline{\varepsilon}\vec{\underline{E}} = -j\omega J_s$$

$$\underline{\nu} = \frac{1}{\underline{\mu}} =
\frac{\mu^{'}+j\mu^{''}}{(\mu^{'}-j\mu^{''})(\mu^{'}+j\mu^{''})} =
\frac{\mu^{'}+j\mu^{''}}{|\underline{\mu}|^2}$$

$$\underline{\varepsilon} = \varepsilon_0\varepsilon_r - j\frac{\sigma}{\omega}$$

Using Helmholtz decomposition for splitting electric field into solenoidal and irrotational parts.

$$\vec{\underline{E}} = \vec{\underline{E}}_{curl} + \vec{\underline{E}}_{div}$$

$$\nabla\times\underline{\nu}\nabla\times(\vec{\underline{E}}_{curl}+\vec{\underline{E}}_{div}) = -j\omega J_s+\omega^2\underline{\varepsilon}(\vec{\underline{E}}_{curl}+\vec{\underline{E}}_{div})$$

$$\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}_{curl} = -j\omega
J_s+\omega^2\underline{\varepsilon}\vec{\underline{E}}_{curl}+
\omega^2\underline{\varepsilon}\vec{\underline{E}}_{div}$$

Since problem is 2D, also splitting electric field (Uwe:8)

$$\vec{\underline{E}} =
\begin{pmatrix}
\vec{E}_\perp^\Re \\
E_z^\Re
\end{pmatrix}
+j
\begin{pmatrix}
\vec{E}_\perp^\Im \\
E_z^\Im
\end{pmatrix}$$

Fields in $z$ direction propagate without loss

$$F(z)=F_0 \exp{\left(\frac{-j\omega z}{\beta c_0}\right)}$$

Therefore, partial derivative in $\vec{e}_z$ direction gives

$$\frac{\partial F(z)}{\partial z}=-j \frac{\omega}{\beta c_0} F(z)$$

## Solving for irrotational field

Using Gauss's law

$$\nabla\cdot\underline{\varepsilon}\vec{\underline{E}} = \rho$$

And solving for the potential.

$$-\nabla\Phi=\vec{E}_{div}$$

This way solenoidal part will be zeroed out because of the vector identity
$\nabla\times\left(-\nabla\Phi\right)=0$

$$-\nabla\cdot\underline{\varepsilon}\nabla\underline{\Phi} = \frac{1}{\beta c_0}J_s$$
$$-\nabla\cdot\left(\varepsilon_0\varepsilon_r -
j\frac{\sigma}{\omega}\right)\nabla\Phi = \frac{1}{\beta c_0}J_s$$

$$-\nabla\cdot\varepsilon_0\varepsilon_r\nabla\Phi
+j\nabla\cdot\frac{\sigma}{\omega}\nabla\Phi=\frac{1}{\beta c_0}J_s$$

$\Phi$ is a scalar field.

For 2D problem gradient may be decomposed as follows.

$$\nabla\Phi=
\frac{\partial\Phi}{\partial x}\vec{x}+
\frac{\partial\Phi}{\partial y}\vec{y}+
\frac{\partial\Phi}{\partial z}\vec{z}=
\begin{pmatrix}
\frac{\partial}{\partial x} \\
\frac{\partial}{\partial y} \\
\end{pmatrix}\Phi\vec{n}+
\frac{\partial\Phi}{\partial z}\Phi\vec{z}=
\nabla_\perp\Phi\vec{n}-j\frac{\omega}{\beta c_0}\Phi\vec{e}_{z}
$$

Material properties are constant across elements.
So, in the calculations they are handled as if they were constant.

$$-\varepsilon_0\varepsilon_r\nabla\cdot
\left(\nabla_\perp\Phi\vec{n} - j \frac{\omega}{\beta c_0}\Phi\vec{e}_{z}\right)
+j\frac{\sigma}{\omega}\nabla\cdot
\left(\nabla_\perp\Phi\vec{n} - j \frac{\omega}{\beta c_0}\Phi\vec{e}_{z}\right)=\frac{1}{\beta c_0}J_s$$

$$-\varepsilon_0\varepsilon_r\nabla\cdot\nabla_\perp\Phi\vec{n}
+j\frac{\omega\varepsilon}{\beta c_0}\nabla\Phi\vec{e}_z
+j\frac{\sigma}{\omega}\nabla\cdot\nabla_\perp\Phi\vec{n}
+\frac{\sigma}{\beta c_0}\nabla\Phi\vec{e}_{z}
=\frac{1}{\beta c_0}J_s$$

$$
\begin{align}
&-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi\vec{n}
-\varepsilon_0\varepsilon_r\frac{\partial}{\partial z}\nabla_\perp\Phi\vec{n}\vec{e}_{z}
+j\frac{\omega\varepsilon_0\varepsilon_r}{\beta c_0}\nabla_\perp\Phi\vec{e}_{z}\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi\vec{e}_{z}
+\\&+j\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi\vec{n}
+j\frac{\sigma}{\omega}\frac{\partial}{\partial z}\nabla_\perp\Phi\vec{n}\vec{e}_{z}
+\frac{\sigma}{\beta c_0}\nabla_\perp\Phi\vec{e}_{z}\vec{n}
-j\frac{\omega\sigma}{\beta^2 c_0^2}\Phi\vec{e}_{z}
=\frac{1}{\beta c_0}J_s
\end{align}
$$

$\vec{n}$ and $\vec{e}_{z}$ are orthogonal, so their product $\vec{n} \vec{e}_{z} = 0$.

$$-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi\vec{e}_{z}
+j\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi\vec{n}
-j\frac{\omega\sigma}{\beta^2 c_0^2}\Phi\vec{e}_{z}
=\frac{1}{\beta c_0}J_s$$

Split $\underline{\Phi}$ into real and imaginary parts.

$$
\begin{align}
&-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\left(\Phi^\Re+j\Phi^\Im\right)\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\left(\Phi^\Re+j\Phi^\Im\right)\vec{e}_{z}
+\\&+j\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\left(\Phi^\Re+j\Phi^\Im\right)\vec{n}
-j\frac{\omega\sigma}{\beta^2 c_0^2}\left(\Phi^\Re+j\Phi^\Im\right)\vec{e}_{z}
=\frac{1}{\beta c_0}J_s
\end{align}
$$

$$
\begin{align}
&-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
-j\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}
+j\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}
+\\&+j\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
-\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
-j\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}
+\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}
=\frac{1}{\beta c_0}J_s
\end{align}
$$

$$\begin{cases}
\Re\left\{
-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
-j\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}
+j\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}
+\right.
\\
\left.
+j\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
-\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
-j\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}
+\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}
\right\}=
\Re\left\{
\frac{1}{\beta c_0}J_s
\right\} \\
\Im\left\{
-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
-j\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}
+j\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}
+\right.
\\
\left.
+j\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
-\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
-j\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}
+\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}
\right\}=
\Im\left\{
\frac{1}{\beta c_0}J_s
\right\}
\end{cases}$$

Arriving to the equation (Uwe:19)

$$\begin{cases}
-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}
-\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
+\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}=
\frac{1}{\beta c_0}J_s \\
-\varepsilon_0\varepsilon_r\nabla_\perp\cdot\nabla_\perp\Phi^\Im\vec{n}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\Phi^\Im\vec{e}_{z}
+\frac{\sigma}{\omega}\nabla_\perp\cdot\nabla_\perp\Phi^\Re\vec{n}
-\frac{\omega\sigma}{\beta^2 c_0^2}\Phi^\Re\vec{e}_{z}=0
\end{cases}$$

Using Galerkin method lowering equation order by multiplying equation by the test functions $v^\Re$ for real equation
part and $v^\Im$ for imaginary equation part **from the left** and
[integrating by parts](https://en.wikipedia.org/wiki/Integration_by_parts#Higher_dimensions).
Writing down separately for each term.

$$
-\varepsilon_0\varepsilon_r\int_\Omega{v^\Re\nabla_\perp^2\Phi^\Re \;d\Omega}
=\varepsilon_0\varepsilon_r\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Re \;d\Omega}
-\varepsilon_0\varepsilon_r\int_{\partial\Omega}{v^\Re\nabla_\perp\Phi^\Re\vec{n}\;dS}
$$

$$
\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\int_\Omega{v^\Re\Phi^\Re \;d\Omega}
$$

$$
-\frac{\sigma}{\omega}\int_\Omega{v^\Re\nabla_\perp^2\Phi^\Im \;d\Omega}
=\frac{\sigma}{\omega}\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Im \;d\Omega}
-\frac{\sigma}{\omega}\int_{\partial\Omega}{v^\Re\nabla_\perp\Phi^\Im\vec{n}\;dS}
$$

$$
\frac{\omega\sigma}{\beta^2 c_0^2}\int_\Omega{v^\Re\Phi^\Im \;d\Omega}
$$

$$
\frac{1}{\beta c_0}\int_\Omega{v^\Re J_s \;d\Omega}
$$

$$
-\varepsilon_0\varepsilon_r\int_\Omega{v^\Im\nabla_\perp^2\Phi^\Im \;d\Omega}
=\varepsilon_0\varepsilon_r\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Im \;d\Omega}
-\varepsilon_0\varepsilon_r\int_{\partial\Omega}{v^\Im\nabla_\perp\Phi^\Im\vec{n}\;dS}
$$

$$
\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}
$$

$$
\frac{\sigma}{\omega}\int_\Omega{v^\Im\nabla_\perp^2\Phi^\Re \;d\Omega}
=-\frac{\sigma}{\omega}\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Re \;d\Omega}
+\frac{\sigma}{\omega}\int_{\partial\Omega}{v^\Im\nabla_\perp\Phi^\Re\vec{n}\;dS}
$$

$$
-\frac{\omega\sigma}{\beta^2 c_0^2}\int_\Omega{v^\Im\Phi^\Re \;d\Omega}
$$

Using PEC boundary condition thus setting boundary field values to zero.
For trial functions $\Phi^\Re$, $\Phi^\Im$ and test functions $v^\Re$, $v^\Im$
equation looks like this

$$
\begin{align}
&\varepsilon_0\varepsilon_r\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Re \;d\Omega}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\int_\Omega{v^\Re\Phi^\Re \;d\Omega}+
\frac{\sigma}{\omega}\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Im \;d\Omega}+\\
&+\frac{\omega\sigma}{\beta^2 c_0^2}\int_\Omega{v^\Re\Phi^\Im \;d\Omega}+
\varepsilon_0\varepsilon_r\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Im \;d\Omega}+
\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta^2 c_0^2}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}-\\
&-\frac{\sigma}{\omega}\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Re \;d\Omega}
-\frac{\omega\sigma}{\beta^2 c_0^2}\int_\Omega{v^\Im\Phi^\Re \;d\Omega}
=\frac{1}{\beta c_0}\int_\Omega{v^\Re J_s \;d\Omega}
\end{align}
$$

After multiplying both sides of the equation by $\beta c_0$ we get

$$
\begin{align}
&\varepsilon_0\varepsilon_r\beta c_0\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Re \;d\Omega}
+\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta c_0}\int_\Omega{v^\Re\Phi^\Re \;d\Omega}+
\frac{\sigma\beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Re \cdot \nabla_\perp\Phi^\Im \;d\Omega}+\\
&+\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Re\Phi^\Im \;d\Omega}+
\varepsilon_0\varepsilon_r\beta c_0\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Im \;d\Omega}+
\frac{\omega^2\varepsilon_0\varepsilon_r}{\beta c_0}\int_\Omega{v^\Im\Phi^\Im \;d\Omega}-\\
&-\frac{\sigma\beta c_0}{\omega}\int_\Omega{\nabla_\perp v^\Im \cdot \nabla_\perp\Phi^\Re \;d\Omega}
-\frac{\omega\sigma}{\beta c_0}\int_\Omega{v^\Im\Phi^\Re \;d\Omega}
=\int_\Omega{v^\Re J_s \;d\Omega}
\end{align}
$$


After solving for potential we may get the irrotational part of electric field from the equations

$$E_{div\,\perp}^\Re = -\nabla\Phi^\Re$$
$$E_{div\,\perp}^\Im = -\nabla\Phi^\Im$$
$$E_{div\,z}^\Re = -\frac{\omega}{\beta c_0}\Phi^\Im$$
$$E_{div\,z}^\Im = \frac{\omega}{\beta c_0}\Phi^\Re$$

## Solving for solenoidal field

Orienting coordinate system so that section slice is placed on XY plane.
Z axis is oriented according to the right hand rule.
Enclosing surface normal is opposite to Z axis ($\vec{e} = -\vec{n}$).
Section boundary normal points inward, boundary tangential vector is oriented clockwise relative to
the enclosing surface normal direction.

$$\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{E}=-j\omega\vec{J}_S$$

$$\vec{\underline{E}} = \vec{\underline{E}}_{curl} + \vec{\underline{E}}_{div}$$

Since $\nabla\times\vec{\underline{E}}_{div}=0$

$$\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}_{curl}
-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{curl}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{div}
=-j\omega\vec{J}_S$$

$$\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}_{curl}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{curl}
=\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{div}-j\omega\vec{J}_S$$

$$\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}_{curl}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{curl}=\vec{\underline{R}}$$

$$\vec{\underline{R}}=\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{div}-j\omega\vec{\underline{J}}_{s}$$

Using decomposed curl operator (Uwe:24).

$$\nabla\times\vec{F}=
\left(
\begin{array}{cc|c}
0 & -\frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\
\frac{\partial}{\partial z} & 0 & -\frac{\partial}{\partial x} \\ \hline
-\frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0 \\
\end{array}
\right)
\begin{pmatrix}
F_x \\ F_y \\ F_z
\end{pmatrix}
=\left(
\begin{array}{cc|c}
j\hat{\operatorname{Z}} &  & \hat{\operatorname{A}} \\
&  &  \\ \hline
\hat{\operatorname{B}} &  & 0 \\
\end{array}
\right)
\begin{pmatrix}
\vec{F}_\perp \\ F_z
\end{pmatrix}
$$

#### Operator $\hat{\operatorname{Z}}$
operates on $\vec{e}_x$ and $\vec{e}_y$ field components and produces vector with no $\vec{e}_z$ component

$$j\hat{\operatorname{Z}}\vec{F}=
\frac{\partial\left(\vec{e}_y\cdot\vec{F}\right)}{\partial z}\vec{e}_x
-\frac{\partial\left(\vec{e}_x\cdot\vec{F}\right)}{\partial z}\vec{e}_y
=\frac{\partial}{\partial z}\vec{e}_z\times\vec{F}$$

Due to the uniformity of the structure in $\vec{e}_z$ direction
$$
j\hat{\operatorname{Z}}\vec{F}
=\frac{\partial}{\partial z}\vec{e}_z\times\vec{F}
=-j\frac{\omega}{\beta c_0}\vec{e}_z\times\vec{F}
$$

#### Operator $\hat{\operatorname{A}}$
operates on scalar and produces vector with no $\vec{e}_z$ component

$$
\hat{\operatorname{A}}\varphi=
\frac{\partial\varphi}{\partial y}\vec{e}_x
-\frac{\partial\varphi}{\partial x}\vec{e}_y
=
-\vec{e}_z\times\left(
\frac{\partial\varphi}{\partial x}\vec{e}_x
+\frac{\partial\varphi}{\partial y}\vec{e}_y
+\frac{\partial\varphi}{\partial z}\vec{e}_z
\right)
=-\vec{e}_z\times\nabla\varphi
$$

Therefore

$$\hat{\operatorname{A}}\hat{\operatorname{A}}\varphi=0$$
$$\hat{\operatorname{A}}\hat{\operatorname{Z}}\vec{F}=0$$

#### Operator $\hat{\operatorname{B}}$
operates on $\vec{e}_x$ and $\vec{e}_y$ field components and produces vector with only $\vec{e}_z$ component

$$\hat{\operatorname{B}}\vec{F}=
\left(
-\frac{\partial\left(\vec{e}_x\cdot\vec{F}\right)}{\partial y}
+\frac{\partial\left(\vec{e}_y\cdot\vec{F}\right)}{\partial x}
\right)\vec{e}_z
=\vec{e}_z\nabla\times\vec{F}
$$

Therefore

$$\hat{\operatorname{B}}\hat{\operatorname{B}}\vec{F}=0$$
$$\hat{\operatorname{Z}}\hat{\operatorname{B}}\vec{F}=0$$

Decomposing curl operators

$$
\begin{align}
\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}
=&
\nabla\times j\underline{\nu}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
\nabla\times\underline{\nu}\hat{\operatorname{A}}E_z +
\nabla\times\underline{\nu}\hat{\operatorname{B}}\vec{\underline{E}}_\perp \\
=&
-\hat{\operatorname{Z}}\underline{\nu}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
j\hat{\operatorname{A}}\underline{\nu}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
j\hat{\operatorname{B}}\underline{\nu}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp + \\
&j\hat{\operatorname{Z}}\underline{\nu}\hat{\operatorname{A}}E_z +
\hat{\operatorname{A}}\underline{\nu}\hat{\operatorname{A}}E_z +
\hat{\operatorname{B}}\underline{\nu}\hat{\operatorname{A}}E_z + \\
&j\hat{\operatorname{Z}}\underline{\nu}\hat{\operatorname{B}}\vec{\underline{E}}_\perp +
\hat{\operatorname{A}}\underline{\nu}\hat{\operatorname{B}}\vec{\underline{E}}_\perp +
\hat{\operatorname{B}}\underline{\nu}\hat{\operatorname{B}}\vec{\underline{E}}_\perp
\end{align}
$$

Using mentioned above equations

$$\hat{\operatorname{A}}\hat{\operatorname{A}}\varphi=0$$
$$\hat{\operatorname{A}}\hat{\operatorname{Z}}\vec{F}=0$$
$$\hat{\operatorname{B}}\hat{\operatorname{B}}\vec{F}=0$$
$$\hat{\operatorname{Z}}\hat{\operatorname{B}}\vec{F}=0$$

removing zero terms

$$
\begin{align}
\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}
=&
-\hat{\operatorname{Z}}\underline{\nu}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
j\hat{\operatorname{B}}\underline{\nu}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
j\hat{\operatorname{Z}}\underline{\nu}\hat{\operatorname{A}}E_z +
\hat{\operatorname{B}}\underline{\nu}\hat{\operatorname{A}}E_z +
\hat{\operatorname{A}}\underline{\nu}\hat{\operatorname{B}}\vec{\underline{E}}_\perp
\end{align}
$$

The next step is to split this equation into real and imaginary parts.
Firstly, split $\nu = \nu^\Re + j\nu^\Im$.
$$
\begin{align}
\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}
=&
-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}}\vec{\underline{E}}_\perp
-j\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}}\vec{\underline{E}}_\perp
+j\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}}\vec{\underline{E}}_\perp
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}}\vec{\underline{E}}_\perp
+j\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\underline{E}_z \\
&-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\underline{E}_z
+\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}}\underline{E}_z
+j\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}}\underline{E}_z
+\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}\vec{\underline{E}}_\perp
+j\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}\vec{\underline{E}}_\perp
\end{align}
$$

Now, splitting $\vec{\underline{E}}_\perp=\vec{E}_\perp^\Re+j\vec{E}_\perp^\Im$ and $\underline{E}_z=E_z^\Re+j E_z^\Im$

$$
\begin{align}
\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}
=&
-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}}\vec{E}_\perp^\Re
-j\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}}\vec{E}_\perp^\Im
-j\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}}\vec{E}_\perp^\Re
+\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}}\vec{E}_\perp^\Im
+j\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}}\vec{E}_\perp^\Re
-\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}}\vec{E}_\perp^\Im \\
&-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}}\vec{E}_\perp^\Re
-j\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}}\vec{E}_\perp^\Im
+j\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}E_z^\Re
-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}E_z^\Im
-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}E_z^\Re
-j\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}E_z^\Im \\
&+\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}}E_z^\Re
+j\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}}E_z^\Im
+j\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}}E_z^\Re
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}}E_z^\Im
+\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re
+j\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im \\
&+j\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re
-\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im
\end{align}
$$

Grouping terms (Uwe:27)

$$
\begin{align}
\nabla\times\underline{\nu}\nabla\times\vec{\underline{E}}
=
&\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}} & -\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}} \\
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Re \\
E_z^\Re \\
\end{pmatrix}
-\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}} \\
\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Im \\
E_z^\Im \\
\end{pmatrix} \\
+j
&\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}} & -\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}} \\
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Im \\
E_z^\Im \\
\end{pmatrix}
+j
\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}} \\
\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Re \\
E_z^\Re \\
\end{pmatrix}
\end{align}
$$



Matrix equations

$$
\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}} & -\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}} \\
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Re \\
E_z^\Re \\
\end{pmatrix}
$$

and

$$
-\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}} \\
\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Im \\
E_z^\Im \\
\end{pmatrix}
$$

produce real valued fields.
Therefore they will live in "real" function space.

Matrix equations

$$
j\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}} & -\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}} \\
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Im \\
E_z^\Im \\
\end{pmatrix}
$$

and

$$
j
\begin{pmatrix}
\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}} & \hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}} \\
\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}} & \hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}} \\
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp^\Re \\
E_z^\Re \\
\end{pmatrix}
$$

produce imaginary valued fields.
Therefore they will live in "imaginary" function space.

Let's rewrite the first equation term in matrix form

$$
\nabla\times\underline{\nu}\nabla\times\vec{E}
=[H]
\begin{pmatrix}
e^\Re_\perp \\
e^\Im_\perp \\
e^\Re_z \\
e^\Im_z \\
\end{pmatrix}
$$

$$
[H]=
\begin{pmatrix}
H^{\Re\Re}_{\perp\perp} &
H^{\Im\Re}_{\perp\perp} &
H^{\Re\Re}_{z\perp} &
H^{\Im\Re}_{z\perp} \\
H^{\Re\Im}_{\perp\perp} &
H^{\Im\Im}_{\perp\perp} &
H^{\Re\Im}_{z\perp} &
H^{\Im\Im}_{z\perp} \\
H^{\Re\Re}_{\perp z} &
H^{\Im\Re}_{\perp z} &
H^{\Re\Re}_{z z} &
H^{\Im\Re}_{z z} \\
H^{\Re\Im}_{\perp z} &
H^{\Im\Im}_{\perp z} &
H^{\Re\Im}_{z z} &
H^{\Im\Im}_{z z} \\
\end{pmatrix}
$$

Elements of matrix $[H]$ are:

<!-- 1 -->
$$
H^{\Re\Re}_{\perp\perp}
\vec{E}^\Re_\perp=
\left(\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}}\right)
\vec{E}^\Re_\perp
$$

<!-- 2 -->
$$
H^{\Im\Re}_{\perp\perp}
\vec{E}^\Im_\perp=
-\left(\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}}\right)
\vec{E}^\Im_\perp
$$

<!-- 3 -->
$$
H^{\Re\Re}_{z\perp}
\vec{E}^\Re_z=
-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}
\vec{E}^\Re_z
$$

<!-- 4 -->
$$
H^{\Im\Re}_{z\perp}
\vec{E}^\Im_z=
-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}
\vec{E}^\Im_z
$$

<!-- 5 -->
$$
H^{\Re\Im}_{\perp\perp}
\vec{E}^\Re_\perp=
j\left(\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}}\right)
\vec{E}^\Re_\perp
$$

<!-- 6 -->
$$
H^{\Im\Im}_{\perp\perp}
\vec{E}^\Im_\perp=
j\left(\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}-\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}}\right)
\vec{E}^\Im_\perp
$$

<!-- 7 -->
$$
H^{\Re\Im}_{z\perp}
\vec{E}^\Re_z=
j\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}
\vec{E}^\Re_z
$$

<!-- 8 -->
$$
H^{\Im\Im}_{z\perp}
\vec{E}^\Im_z=
-j\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}
\vec{E}^\Im_z
$$

<!-- 9 -->
$$
H^{\Re\Re}_{\perp z}
\vec{E}^\Re_\perp=
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}}
\vec{E}^\Re_\perp
$$

<!-- 10 -->
$$
H^{\Im\Re}_{\perp z}
\vec{E}^\Im_\perp=
-\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}}
\vec{E}^\Im_\perp
$$

<!-- 11 -->
$$
H^{\Re\Re}_{z z}
\vec{E}^\Re_z=
\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}}
\vec{E}^\Re_z
$$

<!-- 12 -->
$$
H^{\Im\Re}_{z z}
\vec{E}^\Im_z=
-\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}}
\vec{E}^\Im_z
$$

<!-- 13 -->
$$
H^{\Re\Im}_{\perp z}
\vec{E}^\Re_\perp=
j\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}}
\vec{E}^\Re_\perp
$$

<!-- 14 -->
$$
H^{\Im\Im}_{\perp z}
\vec{E}^\Im_\perp=
-j\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}}
\vec{E}^\Im_\perp
$$

<!-- 15 -->
$$
H^{\Re\Im}_{z z}
\vec{E}^\Re_z=
j\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}}
\vec{E}^\Re_z
$$

<!-- 16 -->
$$
H^{\Im\Im}_{z z}
\vec{E}^\Im_z=
j\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}}
\vec{E}^\Im_z
$$

Now we need to derive a weak form of the $[H]$ matrix.

#### Decomposed curl operator identities

Before deriving weak form of the $[H]$ matrix we need to derive some useful identities.

Using vector calculus identity
$$
\nabla\times\left(\varphi\vec{F}\right)
=\varphi\left(\nabla\times\vec{F}\right)
+\left(\nabla\varphi\right)\times\vec{F}
=\varphi\left(\nabla\times\vec{F}\right)
-\vec{F}\times\left(\nabla\varphi\right)
$$

$$
\varphi\left(\nabla\times\vec{F}\right)
=\nabla\times\left(\varphi\vec{F}\right)
+\vec{F}\times\left(\nabla\varphi\right)
$$

and Stokes' theorem

$$\int_\Omega{\nabla\times\vec{F}\vec{n}\;d\Omega}=\int_{\partial\Omega}{\vec{F}\vec{\tau}\;dS}$$

we may write (Uwe:28)

$$
\begin{align}
\int_\Omega{\varphi\left(\hat{\operatorname{B}}\vec{F}\right)\;d\Omega}
=&\int_\Omega{\varphi\left(\vec{e}_z\nabla\times\vec{F}\right)\;d\Omega}
=-\int_\Omega{\vec{n}\varphi\nabla\times\vec{F}\;d\Omega} \\
=&-\int_\Omega{\vec{n}\nabla\times\left(\varphi\vec{F}\right)\;d\Omega}
-\int_\Omega{\vec{n}\vec{F}\times\left(\nabla\varphi\right)\;d\Omega} \\
=&\int_\Omega{\vec{F}\left(-\vec{e}_z\times\left(\nabla\varphi\right)\right)\;d\Omega}
-\int_\Omega{\vec{n}\nabla\times\left(\varphi\vec{F}\right)\;d\Omega} \\
=&\int_\Omega{\vec{F}\hat{\operatorname{A}}\varphi \;d\Omega}
-\int_{\partial\Omega}{\varphi\vec{F}\vec{\tau}\;dS}
\end{align}
$$

Using vector identities

$$
\nabla\times\left(\varphi\vec{A}\right)
=\varphi\left(\nabla\times\vec{A}\right)
+\left(\nabla\varphi\right)\times\vec{A}
$$

$$
\vec{a}\left(\vec{b}\times\vec{c}\right)
=\vec{b}\left(\vec{c}\times\vec{a}\right)
=\vec{c}\left(\vec{a}\times\vec{b}\right)
$$

Stokes' theorem

$$
\int_{\partial\Omega}{\left(\varphi\vec{A}\right)\vec{\tau}\;dS}
=\int_\Omega{\nabla\times\left(\varphi\vec{A}\right)\vec{n}\;d\Omega}
=\int_{\Omega}{\varphi\left(\nabla\times\vec{A}\right)\vec{n}\;d\Omega}
+\int_{\Omega}{\left(\nabla\varphi\right)\times\vec{A}\vec{n}\;d\Omega}
$$

$$
\int_{\Omega}{\varphi\left(\nabla\times\vec{A}\right)\vec{n}\;d\Omega}
=-\int_{\Omega}{\left(\nabla\varphi\right)\times\vec{A}\vec{n}\;d\Omega}
+\int_{\partial\Omega}{\left(\varphi\vec{A}\right)\vec{\tau}\;dS}
$$

we may write

$$
\begin{align}
\int_{\Omega}{\varphi\hat{\operatorname{B}}\vec{F}\;d\Omega}
=&\int_{\Omega}{\varphi\left(\vec{e}_z\nabla\times\vec{F}\right)\;d\Omega} \\
=&-\int_{\Omega}{\vec{n}\varphi\left(\nabla\times\vec{F}\right)\;d\Omega} \\
=&\int_{\Omega}{\vec{n}\left(\nabla\varphi\right)\times\vec{F}\;d\Omega}
-\int_{\partial\Omega}{\left(\varphi\vec{F}\right)\vec{\tau}\;dS} \\
=&\int_{\Omega}{\vec{F}\left(-\vec{e}_z\times\nabla\varphi\right)\;d\Omega}
-\int_{\partial\Omega}{\left(\varphi\vec{F}\right)\vec{\tau}\;dS} \\
=&\int_{\Omega}{\vec{F}\left(\hat{\operatorname{A}}\varphi\right)\;d\Omega}
-\int_{\partial\Omega}{\left(\varphi\vec{F}\right)\vec{\tau}\;dS}
\end{align}
$$

Using Ostrogradsky-Gauss' theorem

$$\int_\Omega{\nabla\cdot\vec{F}\;d\Omega}=\int_{\partial\Omega}{\vec{F}\vec{n}\;dS}$$

and vector identities

$$
\nabla\left(\varphi\vec{A}\right)=\left(\nabla\varphi\right)^T\vec{A}+\varphi\nabla\vec{A}
$$

$$
\int_{\partial\Omega}{\varphi\vec{A}\vec{n}\;dS}
=\int_\Omega{\nabla\varphi\vec{A}\;d\Omega}
=\int_{\Omega}{\left(\nabla\varphi\right)^T\vec{A}\;d\Omega}
+\int_{\Omega}{\varphi\nabla\vec{A}\;d\Omega}
$$

$$
\int_{\Omega}{\varphi\nabla\vec{A}\;d\Omega}
=-\int_{\Omega}{\left(\nabla\varphi\right)^T\vec{A}\;d\Omega}
+\int_{\partial\Omega}{\varphi\vec{A}\vec{n}\;dS}
$$

#### Converting term $\hat{\operatorname{A}}\nu\hat{\operatorname{B}}\vec{E}_\perp$
into [weak
form](https://en.wikipedia.org/wiki/Weak_formulation) by multiplying it by the test function
$\vec{w}$ and
integrating over the calculation domain $\Omega$.

$$
\begin{align}
\int_\Omega{
\vec{w}
\left(
\hat{\operatorname{A}}\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
=&\int_\Omega{
\vec{w}
\left(
-\vec{e}_z\times\nabla\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
} \\
=&-\vec{e}_z\times
\int_\Omega{
\vec{w}
\left(
\nabla\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
} \\
=&\vec{e}_z\times
\int_\Omega{
\left(
\nabla
\vec{w}
\right)^T
\left(
\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
-\vec{e}_z\times
\int_{\partial\Omega}{
\vec{w}
\left(
\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\vec{n}
\;dS
} \\
=&\int_\Omega{
\left(
\vec{e}_z\times\nabla\vec{w}
\right)^T
\left(
\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
-\int_{\partial\Omega}{
\vec{w}
\left(
\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\left(
\vec{e}_z\times\vec{n}
\right)
\;dS
} \\
=&\int_\Omega{
\left(
\hat{\operatorname{B}}\vec{w}
\right)
\left(
\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
+\int_{\partial\Omega}{
\vec{w}
\left(
\nu\hat{\operatorname{B}}\vec{E}_\perp
\right)
\vec{\tau}
\;dS
}
\end{align}
$$

#### Converting term $\hat{\operatorname{B}}\nu\hat{\operatorname{A}}E_z$

into [weak
form](https://en.wikipedia.org/wiki/Weak_formulation) by multiplying it by the test function
$v$ and
integrating over the calculation domain $\Omega$.

$$
\begin{align}
\int_\Omega{
v
\left(
\hat{\operatorname{B}}\nu\hat{\operatorname{A}}E_z
\right)
\;d\Omega
}
=&\int_\Omega{
v
\left(
\vec{e}_z\nabla\times\nu\hat{\operatorname{A}}E_z
\right)
\;d\Omega
} \\
=&-\int_\Omega{
v
\left(
\nabla\times\nu\hat{\operatorname{A}}E_z
\right)
\vec{n}
\;d\Omega
} \\
=&\int_\Omega{
\left(
\nabla v
\right)
\times
\left(
\nu\hat{\operatorname{A}}E_z
\right)
\vec{n}
\;d\Omega
}
-\int_{\partial\Omega}{
v
\left(
\nu\hat{\operatorname{A}}E_z
\right)
\vec{\tau}
\;dS
} \\
=&\int_\Omega{
\left(
-\vec{e}_z\times\nabla v
\right)
\left(
\nu\hat{\operatorname{A}}E_z
\right)
\;d\Omega
}
-\int_{\partial\Omega}{
v
\left(
\nu\hat{\operatorname{A}}E_z
\right)
\vec{\tau}
\;dS
} \\
=&\int_\Omega{
\left(
\hat{\operatorname{A}}v
\right)
\left(
\nu\hat{\operatorname{A}}E_z
\right)
\;d\Omega
}
-\int_{\partial\Omega}{
v
\left(
\nu\hat{\operatorname{A}}E_z
\right)
\vec{\tau}
\;dS
}
\end{align}
$$

We may now use decomposed curl operator identities to convert $[H]$ matrix into weak form.
Using Galerkin method lowering equations order by multiplying them by the test functions and
[integrating by parts](https://en.wikipedia.org/wiki/Integration_by_parts#Higher_dimensions) (where necessary).
Test functions for real and imaginary transverse field, real and imaginary longitudinal field are
$\vec{w}^\Re$, $\vec{w}^\Im$, $v^\Re$ and $v^\Im$ respectively (Uwe:38).

$$
[\hat{H}]=
\begin{pmatrix}
\hat{H}^{\Re\Re}_{\perp\perp} &
\hat{H}^{\Im\Re}_{\perp\perp} &
\hat{H}^{\Re\Re}_{z\perp} &
\hat{H}^{\Im\Re}_{z\perp} \\
\hat{H}^{\Re\Im}_{\perp\perp} &
\hat{H}^{\Im\Im}_{\perp\perp} &
\hat{H}^{\Re\Im}_{z\perp} &
\hat{H}^{\Im\Im}_{z\perp} \\
\hat{H}^{\Re\Re}_{\perp z} &
\hat{H}^{\Im\Re}_{\perp z} &
\hat{H}^{\Re\Re}_{z z} &
\hat{H}^{\Im\Re}_{z z} \\
\hat{H}^{\Re\Im}_{\perp z} &
\hat{H}^{\Im\Im}_{\perp z} &
\hat{H}^{\Re\Im}_{z z} &
\hat{H}^{\Im\Im}_{z z} \\
\end{pmatrix}
$$

<!-- 1 -->
$$
\begin{align}
\hat{H}^{\Re\Re}_{\perp\perp}
&=\int_\Omega{\vec{w}^\Re H^{\Re\Re}_{\perp\perp}\vec{E}^\Re_\perp\;d\Omega} \\
&=\int_\Omega{\vec{w}^\Re\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}\vec{E}^\Re_\perp\;d\Omega}
-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
+\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)
\vec{\tau}\;dS}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Re\vec{E}^\Re_\perp\;d\Omega}
\end{align}
$$

<!-- 2 -->
$$
\begin{align}
\hat{H}^{\Im\Re}_{\perp\perp}
&=\int_\Omega{\vec{w}^\Re H^{\Im\Re}_{\perp\perp}
\vec{E}^\Im_\perp\;d\Omega} \\
&=-\int_\Omega{\vec{w}^\Re\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}\vec{E}^\Im_\perp\;d\Omega}
+\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\;d\Omega} \\
&=-\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
-\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
-\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Im\vec{E}^\Im_\perp\;d\Omega}
\end{align}
$$

<!-- 3 -->
$$
\hat{H}^{\Re\Re}_{z\perp}
=\int_\Omega{\vec{w}^\Re H^{\Re\Re}_{z\perp}
\vec{E}^\Re_z\;d\Omega}
=-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
$$

<!-- 4 -->
$$
\hat{H}^{\Im\Re}_{z\perp}
=\int_\Omega{\vec{w}^\Re H^{\Im\Re}_{z\perp}
\vec{E}^\Im_z\;d\Omega}
=-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
$$

<!-- 5 -->
$$
\begin{align}
\hat{H}^{\Re\Im}_{\perp\perp}
&=\int_\Omega{\vec{w}^\Im H^{\Re\Im}_{\perp\perp}
\vec{E}^\Re_\perp\;d\Omega} \\
&=\int_\Omega{\vec{w}^\Im\hat{\operatorname{A}}\nu^\Im\hat{\operatorname{B}}\vec{E}^\Re_\perp\;d\Omega}
-\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
+\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\vec{\tau}\;dS}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Im\vec{E}^\Re_\perp\;d\Omega}
\end{align}
$$

<!-- 6 -->
$$
\begin{align}
\hat{H}^{\Im\Im}_{\perp\perp}
&=\int_\Omega{\vec{w}^\Im H^{\Im\Im}_{\perp\perp}
\vec{E}^\Im_\perp\;d\Omega} \\
&=\int_\Omega{\vec{w}^\Im\hat{\operatorname{A}}\nu^\Re\hat{\operatorname{B}}\vec{E}^\Im_\perp\;d\Omega}
-\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
+\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Re\vec{E}^\Im_\perp\;d\Omega}
\end{align}
$$

<!-- 7 -->
$$
\hat{H}^{\Re\Im}_{z\perp}
=\int_\Omega{\vec{w}^\Im H^{\Re\Im}_{z\perp}\vec{E}^\Re_z\;d\Omega}
=\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
$$

<!-- 8 -->
$$
\hat{H}^{\Im\Im}_{z\perp}
=\int_\Omega{\vec{w}^\Im H^{\Im\Im}_{z\perp}\vec{E}^\Im_z\;d\Omega}
=-\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
$$

<!-- 9 -->
$$
\begin{align}
\hat{H}^{\Re\Re}_{\perp z}
&=\int_\Omega{v^\Re H^{\Re\Re}_{\perp z}\vec{E}^\Re_\perp\;d\Omega} \\
&=-\int_\Omega{v^\Re\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\;d\Omega} \\
&=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
\left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega} +
\int_{\partial\Omega}{v^\Re\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\vec{\tau}\;dS}
\end{align}
$$

<!-- 10 -->
$$
\begin{align}
\hat{H}^{\Im\Re}_{\perp z}
&=\int_\Omega{v^\Re H^{\Im\Re}_{\perp z} \vec{E}^\Im_\perp\;d\Omega} \\
&=-\int_\Omega{v^\Re\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\;d\Omega} \\
&=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
\left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega} +
\int_{\partial\Omega}{v^\Re\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\vec{\tau}\;dS}
\end{align}
$$

<!-- 11 -->
$$
\begin{align}
\hat{H}^{\Re\Re}_{z z}
&=\int_\Omega{v^\Re H^{\Re\Re}_{z z}\vec{E}^\Re_z\;d\Omega} \\
&=\int_\Omega{v^\Re\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{A}} v^\Re\right)
\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
-\int_{\partial\Omega}{v^\Re\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\vec{\tau}\;dS}
\end{align}
$$

<!-- 12 -->
$$
\begin{align}
\hat{H}^{\Im\Re}_{z z}
&=\int_\Omega{v^\Re H^{\Im\Re}_{z z}\vec{E}^\Im_z\;d\Omega} \\
&=-\int_\Omega{v^\Re\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega} \\
&=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
+\int_{\partial\Omega}{v^\Re\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\vec{\tau}\;dS}
\end{align}
$$

<!-- 13 -->
$$
\begin{align}
\hat{H}^{\Re\Im}_{\perp z}
&=\int_\Omega{v^\Im H^{\Re\Im}_{\perp z}\vec{E}^\Re_\perp\;d\Omega} \\
&=\int_\Omega{v^\Im\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
\left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega} -
\int_{\partial\Omega}{v^\Im\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\vec{\tau}\;dS}
\end{align}
$$

<!-- 14 -->
$$
\begin{align}
\hat{H}^{\Im\Im}_{\perp z}
&=\int_\Omega{v^\Im H^{\Im\Im}_{\perp z}\vec{E}^\Im_\perp\;d\Omega} \\
&=-\int_\Omega{v^\Im\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\;d\Omega} \\
&=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
\left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega} +
\int_{\partial\Omega}{v^\Im\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\vec{\tau}\;dS}
\end{align}
$$

<!-- 15 -->
$$
\begin{align}
\hat{H}^{\Re\Im}_{z z}
&=\int_\Omega{v^\Im H^{\Re\Im}_{z z}\vec{E}^\Re_z\;d\Omega} \\
&=\int_\Omega{v^\Im\hat{\operatorname{B}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
-\int_{\partial\Omega}{v^\Im\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\vec{\tau}\;dS}
\end{align}
$$

<!-- 16 -->
$$
\begin{align}
\hat{H}^{\Im\Im}_{z z}
&=\int_\Omega{v^\Im H^{\Im\Im}_{z z}\vec{E}^\Im_z\;d\Omega} \\
&=\int_\Omega{v^\Im\hat{\operatorname{B}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
-\int_{\partial\Omega}{v^\Im\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\vec{\tau}\;dS}
\end{align}
$$

Writing whole equation in matrix form (Uwe:36)
$$
[S_{curlcurl}+M_{\varepsilon}+M_{SIBC}] E_{curl}=
[N_{\varepsilon}] + [J_s]
$$

$S_{curlcurl}$ term is equal to the $\hat{H}$ minus the boundary integrals.
$M_{SIBC}$ contains all the boundary integrals from the $\hat{H}$.

$M_{\varepsilon}$ and $N_{\varepsilon}$ are terms $\omega^2\underline{\varepsilon}\vec{\underline{E}}_{curl}$
and $\omega^2\underline{\varepsilon}\vec{\underline{E}}_{div}$
written in weak from.

Let's split these terms in real/imaginary and longitudinal and transverse parts

$$
\begin{align}
\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}
=&\omega^2\varepsilon_0\varepsilon_r\vec{\underline{E}}
-j\sigma\omega\vec{\underline{E}} \\
=&\omega^2\varepsilon_0\varepsilon_r\vec{\underline{E}}_\perp
+\omega^2\varepsilon_0\varepsilon_r\underline{E}_z
-j\sigma\omega\vec{\underline{E}}_\perp
-j\sigma\omega\underline{E}_z \\
=&\omega^2\varepsilon_0\varepsilon_r\vec{E}^\Re_\perp
+j\omega^2\varepsilon_0\varepsilon_r\vec{E}^\Im_\perp
+\omega^2\varepsilon_0\varepsilon_r E_z^\Re
+j\omega^2\varepsilon_0\varepsilon_r E_z^\Im \\
&-j\sigma\omega\vec{E}^\Re_\perp
+\sigma\omega\vec{E}^\Im_\perp
-j\sigma\omega E_z^\Re
+\sigma\omega E_z^\Im
\end{align}
$$

$J_s$ is a source term


The whole matrix is (Uwe:37-40)
$$
\left[
\begin{pmatrix}
S^{\Re\Re}_{\perp\perp} &
S^{\Im\Re}_{\perp\perp} &
S^{\Re\Re}_{z\perp} &
S^{\Im\Re}_{z\perp} \\
S^{\Re\Im}_{\perp\perp} &
S^{\Im\Im}_{\perp\perp} &
S^{\Re\Im}_{z\perp} &
S^{\Im\Im}_{z\perp} \\
S^{\Re\Re}_{\perp z} &
S^{\Im\Re}_{\perp z} &
S^{\Re\Re}_{z z} &
S^{\Im\Re}_{z z} \\
S^{\Re\Im}_{\perp z} &
S^{\Im\Im}_{\perp z} &
S^{\Re\Im}_{z z} &
S^{\Im\Im}_{z z} \\
\end{pmatrix}
+
\begin{pmatrix}
M^{\Re\Re}_{\varepsilon\perp} &
M^{\Im\Re}_{\sigma\perp} & 0 & 0 \\
M^{\Re\Im}_{\sigma\perp} &
M^{\Im\Im}_{\varepsilon\perp} & 0 & 0 \\
0 & 0 &
M^{\Re\Re}_{\varepsilon z} &
M^{\Im\Re}_{\sigma z} \\
0 & 0 &
M^{\Re\Im}_{\sigma z} &
M^{\Im\Im}_{\varepsilon z} \\
\end{pmatrix}
+
\begin{pmatrix}
B^{\Re\Re}_{\perp\perp} &
B^{\Im\Re}_{\perp\perp} & 0 & 0 \\
B^{\Re\Im}_{\perp\perp} &
B^{\Im\Im}_{\perp\perp} & 0 & 0 \\
B^{\Re\Re}_{\perp z} &
B^{\Im\Re}_{\perp z} &
B^{\Re\Re}_{zz} &
B^{\Im\Re}_{zz} \\
B^{\Re\Re}_{\perp z} &
B^{\Im\Re}_{\perp z} &
B^{\Re\Im}_{zz} &
B^{\Im\Im}_{zz} \\
\end{pmatrix}
\right]
E_{curl}=R
$$

$E_{curl}$ terms

<!-- 1 -->
$$
S^{\Re\Re}_{\perp\perp}
=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Re\vec{E}^\Re_\perp\;d\Omega}
$$

<!-- 2 -->
$$
S^{\Im\Re}_{\perp\perp}
=-\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Re\right)
\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
-\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Re\nu^\Im\vec{E}^\Im_\perp\;d\Omega}
$$

<!-- 3 -->
$$
S^{\Re\Re}_{z\perp}
=-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
$$

<!-- 4 -->
$$
S^{\Im\Re}_{z\perp}
=-\int_\Omega{\vec{w}^\Re\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
$$

<!-- 5 -->
$$
S^{\Re\Im}_{\perp\perp}
=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\;d\Omega}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Im\vec{E}^\Re_\perp\;d\Omega}
$$

<!-- 6 -->
$$
S^{\Im\Im}_{\perp\perp}
=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}^\Im\right)
\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\;d\Omega}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}^\Im\nu^\Re\vec{E}^\Im_\perp\;d\Omega}
$$

<!-- 7 -->
$$
S^{\Re\Im}_{z\perp}
=\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Re\hat{\operatorname{A}}\vec{E}^\Re_z\;d\Omega}
$$

<!-- 8 -->
$$
S^{\Im\Im}_{z\perp}
=-\int_\Omega{\vec{w}^\Im\hat{\operatorname{Z}}\nu^\Im\hat{\operatorname{A}}\vec{E}^\Im_z\;d\Omega}
$$

<!-- 9 -->
$$
S^{\Re\Re}_{\perp z}
=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
\left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega}
$$

<!-- 10 -->
$$
S^{\Im\Re}_{\perp z}
=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
\left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega}
$$

<!-- 11 -->
$$
S^{\Re\Re}_{z z}
=\int_\Omega{\left(\hat{\operatorname{A}} v^\Re\right)
\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
$$

<!-- 12 -->
$$
S^{\Im\Re}_{z z}
=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Re\right)
\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
$$

<!-- 13 -->
$$
S^{\Re\Im}_{\perp z}
=\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
\left(\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\right)\;d\Omega}
$$

<!-- 14 -->
$$
S^{\Im\Im}_{\perp z}
=-\int_\Omega{\left(\hat{\operatorname{A}}v^\Im\right)
\left(\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\right)\;d\Omega}
$$

<!-- 15 -->
$$
S^{\Re\Im}_{z z}
=\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;d\Omega}
$$

<!-- 16 -->
$$
S^{\Im\Im}_{z z}
=\int_\Omega{\left(\hat{\operatorname{A}} v^\Im\right)
\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;d\Omega}
$$

<!-- 1 -->
$$
M^{\Re\Re}_{\varepsilon\perp}=
-\omega^2\int_{\Omega}{\vec{w}^\Re\varepsilon_0\varepsilon_r\vec{E}^\Re_\perp\;d\Omega}
$$

<!-- 2 -->
$$
M^{\Im\Re}_{\sigma\perp}=
-\omega\int_{\Omega}{\vec{w}^\Re\sigma\vec{E}^\Im_\perp\;d\Omega}
$$

<!-- 5 -->
$$
M^{\Re\Im}_{\sigma\perp}=
\omega\int_{\Omega}{\vec{w}^\Im\sigma\vec{E}^\Re_\perp\;d\Omega}
$$

<!-- 6 -->
$$
M^{\Im\Im}_{\varepsilon\perp}=
-\omega^2\int_{\Omega}{\vec{w}^\Im\varepsilon_0\varepsilon_r\vec{E}^\Im_\perp\;d\Omega}
$$

<!-- 11 -->
$$
M^{\Re\Re}_{\varepsilon z}=
-\omega^2\int_{\Omega}{v^\Re\varepsilon_0\varepsilon_r E^\Re_z\;d\Omega}
$$

<!-- 12 -->
$$
M^{\Im\Re}_{\sigma z}=
-\omega\int_{\Omega}{v^\Re\sigma E^\Im_z\;d\Omega}
$$

<!-- 15 -->
$$
M^{\Re\Im}_{\sigma z}=
\omega\int_{\Omega}{v^\Im\sigma E^\Re_z\;d\Omega}
$$

<!-- 16 -->
$$
M^{\Im\Im}_{\varepsilon z}=
-\omega^2\int_{\Omega}{v^\Im\varepsilon_0\varepsilon_r E^\Im_z\;d\Omega}
$$

<!-- 1 -->
$$
B^{\Re\Re}_{\perp\perp}
=\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)
\vec{\tau}\;dS}
$$

<!-- 2 -->
$$
B^{\Im\Re}_{\perp\perp}
=-\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
$$

<!-- 5 -->
$$
B^{\Re\Im}_{\perp\perp}
=\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\vec{\tau}\;dS}
$$

<!-- 6 -->
$$
B^{\Im\Im}_{\perp\perp}
=\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
$$

<!-- 9 -->
$$
B^{\Re\Re}_{\perp z}
=\int_{\partial\Omega}{v^\Re\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\;dS}
$$

<!-- 10 -->
$$
B^{\Im\Re}_{\perp z}
=\int_{\partial\Omega}{v^\Re\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\;dS}
$$

<!-- 11 -->
$$
B^{\Re\Re}_{z z}
=-\int_{\partial\Omega}{v^\Re\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;dS}
$$

<!-- 12 -->
$$
B^{\Im\Re}_{z z}
=\int_{\partial\Omega}{v^\Re\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;dS}
$$

<!-- 13 -->
$$
B^{\Re\Im}_{\perp z}
=-\int_{\partial\Omega}{v^\Im\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\;dS}
$$

<!-- 14 -->
$$
B^{\Im\Im}_{\perp z}
=\int_{\partial\Omega}{v^\Im\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\;dS}
$$

<!-- 15 -->
$$
B^{\Re\Im}_{z z}
=-\int_{\partial\Omega}{v^\Im\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;dS}
$$

<!-- 16 -->
$$
B^{\Im\Im}_{z z}
=-\int_{\partial\Omega}{v^\Im\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;dS}
$$

$E_{div}$ terms

<!-- 1 -->
$$
N^{\Re\Re}_{\varepsilon\perp}=
\omega^2\int_{\Omega}{\vec{w}^\Re\varepsilon_0\varepsilon_r\vec{E}^\Re_\perp\;d\Omega}
$$

<!-- 2 -->
$$
N^{\Im\Re}_{\sigma\perp}=
\omega\int_{\Omega}{\vec{w}^\Re\sigma\vec{E}^\Im_\perp\;d\Omega}
$$

<!-- 5 -->
$$
N^{\Re\Im}_{\sigma\perp}=
-\omega\int_{\Omega}{\vec{w}^\Im\sigma\vec{E}^\Re_\perp\;d\Omega}
$$

<!-- 6 -->
$$
N^{\Im\Im}_{\varepsilon\perp}=
\omega^2\int_{\Omega}{\vec{w}^\Im\varepsilon_0\varepsilon_r\vec{E}^\Im_\perp\;d\Omega}
$$

<!-- 11 -->
$$
N^{\Re\Re}_{\varepsilon z}=
\omega^2\int_{\Omega}{v^\Re\varepsilon_0\varepsilon_r E^\Re_z\;d\Omega}
$$

<!-- 12 -->
$$
N^{\Im\Re}_{\sigma z}=
\omega\int_{\Omega}{v^\Re\sigma E^\Im_z\;d\Omega}
$$

<!-- 15 -->
$$
N^{\Re\Im}_{\sigma z}=
-\omega\int_{\Omega}{v^\Im\sigma E^\Re_z\;d\Omega}
$$

<!-- 16 -->
$$
N^{\Im\Im}_{\varepsilon z}=
\omega^2\int_{\Omega}{v^\Im\varepsilon_0\varepsilon_r E^\Im_z\;d\Omega}
$$

$J_{s}$ term

$$
J_s^\Im =
-\omega\int_{\Omega}{v^\Im J_s\;d\Omega}
$$

Now we have all the equations to obtain $\vec{\underline{E}}=\vec{\underline{E}}_{curl}+\vec{\underline{E}}_{div}$ for
the chosen source function $J_{S}$.
These can be inserted into the equations for impedance (Uwe:4-6)
$$\begin{align}\underline{Z}_\parallel(\omega)
&=-\frac{1}{q^2}\int_{beam}{\vec{\underline{E}}\vec{\underline{J}}^{*}_\parallel\;dV} \\
&=-\frac{1}{q^2}\int_{beam}{\left(\vec{E}^\Re + j\vec{E}^\Im\right)
\left(\vec{J}_\parallel^\Re-j\vec{J}_\parallel^\Im\right)\;dV} \\
&=-\frac{1}{q^2}\left(
\int_{beam}{\vec{E}^\Re \vec{J}_\parallel^\Re\;dV}
-j\int_{beam}{\vec{E}^\Re \vec{J}_\parallel^\Im\;dV}
+j\int_{beam}{\vec{E}^\Im \vec{J}_\parallel^\Re\;dV}
+\int_{beam}{\vec{E}^\Im \vec{J}_\parallel^\Im\;dV}
\right)
\end{align}$$

For the purely real source function $\vec{J}$ and a step function
$$
\delta=
\begin{cases}
1 \quad (x,y)\in\Omega_{beam} \\
0 \quad (x,y)\notin\Omega_{beam}
\end{cases}
$$
$$\underline{Z}_\parallel(\omega)
=-\frac{1}{q^2}\left(
\int_{\Omega}{\vec{E}^\Re \vec{J}_\parallel^\Re \delta \;d\Omega}
+j\int_{\Omega}{\vec{E}^\Im \vec{J}_\parallel^\Re \delta \;d\Omega}
\right)$$

Similarly
$$\begin{align}
\underline{Z}_\perp(\omega)
&=-\frac{\beta c}{(q d_r)^2\omega}\int_{beam}{\vec{\underline{E}}\vec{\underline{J}}^{*}_\perp\;dV} \\
&=-\frac{\beta c}{(q d_r)^2\omega}\left(
\int_{\Omega}{\vec{E}^\Re \vec{J}_\perp^\Re \delta \;d\Omega}
+j\int_{\Omega}{\vec{E}^\Im \vec{J}_\perp^\Re \delta \;d\Omega}
\right)
\end{align}$$

## SIBC

When calculation domain is enclosed in perfectly conducting material (PEC), tangential electric field
$\vec{\underline{E}}_\tau=0$.
When using [Nédélec (first kind)](https://defelement.com/elements/nedelec1.html) finite elements, this condition may be
enforces by setting zero Dirichlet boundary conditions.
This also eliminates all the boundary integrals.

Skin depth of a metal is given by
$$\delta=\sqrt{\frac{2}{\mu_0\omega\sigma}}$$
Impedance of the metal surface is
$$Z=\frac{E_x}{H_y}=(1+j)\sqrt{\frac{\mu\omega}{2\sigma}}$$
Writing down Faraday's law and applying curl decomposition
$$\nabla\times\vec{\underline{E}}=-j\omega\mu\vec{\underline{H}}$$
$$j\hat{Z}\vec{\underline{E}}_\perp+\hat{A}\vec{\underline{E}}_\perp+\hat{B}\underline{E}_z
=-j\omega\mu\vec{\underline{H}}$$
Extracting tangential and $\vec{e}_z$ field parts
$$j\hat{Z}\vec{\underline{E}}_\tau+\hat{A}\underline{E}_z
=-j\omega\mu\vec{\underline{H}}_\tau$$
$$\hat{B}\vec{\underline{E}}_\tau=-j\omega\mu\underline{H}_z$$
Noticing that term $j\hat{Z}\vec{E}_\tau$ is normal to the boundary and therefore will not contribute to the boundary integral.
$$\begin{align}
\hat{A}\underline{E}_z=-j\omega\mu\vec{\underline{H}}_\tau
&=-\frac{j\omega\mu}{(1+j)\sqrt{\frac{\mu\omega}{2\sigma}}}\vec{\underline{E}}_z
=-\frac{j\omega\mu\sqrt{2\sigma}}{(1+j)\sqrt{\omega\mu}}\vec{\underline{E}}_z
=-\frac{2j\sqrt{\omega\mu\sigma}}{(1+j)\sqrt{2}}\vec{\underline{E}}_z
=-\frac{2j}{(1+j)\delta}\vec{\underline{E}}_z \\
&=-\frac{2j(1-j)}{(1+j)(1-j)\delta}\vec{\underline{E}}_z
=-\frac{2j+2}{(1^2-j^2)\delta}\vec{\underline{E}}_z
=-\frac{j+1}{\delta}\vec{\underline{E}}_z
\end{align}$$
$$\hat{B}\vec{\underline{E}}_\perp=-j\omega\mu\underline{H}_z
=-\frac{j+1}{\delta}\underline{E}_\perp$$
where $\delta=\sqrt{\frac{2}{\omega\mu\sigma}}$ is a skin depth.
Now, applying these equations to $[M_{SIBC}]$ matrix terms.

<!-- 1 -->
$$
\int_{\partial\Omega}{\left(\vec{w}^\Re\vec{\tau}\right)\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Re\vec{\tau}\right)\;dS}
=-\int_{\partial\Omega}{\left(\vec{w}^\Re\vec{\tau}\right)\left(\frac{\nu^\Re}{\delta}\vec{E}_\perp^\Re\vec{\tau}\right)\;dS}
-\int_{\partial\Omega}{\left(\vec{w}^\Im\vec{\tau}\right)\left(\frac{\nu^\Re}{\delta}\vec{E}_\perp^\Re\vec{\tau}\right)\;dS}
$$

<!-- 2 -->
$$
-\int_{\partial\Omega}{\vec{w}^\Re\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
=\int_{\partial\Omega}{\left(\vec{w}^\Re\vec{\tau}\right)\left(\frac{\nu^\Im}{\delta}\vec{E}_\perp^\Im\vec{\tau}\right)\;dS}
+\int_{\partial\Omega}{\left(\vec{w}^\Im\vec{\tau}\right)\left(\frac{\nu^\Im}{\delta}\vec{E}_\perp^\Im\vec{\tau}\right)\;dS}
$$

<!-- 5 -->
$$
\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Im\hat{\operatorname{B}}\vec{E}_\perp^\Re\right)\vec{\tau}\;dS}
=-\int_{\partial\Omega}{\left(\vec{w}^\Im\vec{\tau}\right)\left(\frac{\nu^\Im}{\delta}\vec{E}_\perp^\Re\vec{\tau}\right)\;dS}
+\int_{\partial\Omega}{\left(\vec{w}^\Re\vec{\tau}\right)\left(\frac{\nu^\Im}{\delta}\vec{E}_\perp^\Re\vec{\tau}\right)\;dS}
$$

<!-- 6 -->
$$
\int_{\partial\Omega}{\vec{w}^\Im\left(\nu^\Re\hat{\operatorname{B}}\vec{E}_\perp^\Im\right)\vec{\tau}\;dS}
=-\int_{\partial\Omega}{\left(\vec{w}^\Im\vec{\tau}\right)\left(\frac{\nu^\Re}{\delta}\vec{E}_\perp^\Im\vec{\tau}\right)\;dS}
+\int_{\partial\Omega}{\left(\vec{w}^\Re\vec{\tau}\right)\left(\frac{\nu^\Re}{\delta}\vec{E}_\perp^\Im\vec{\tau}\right)\;dS}
$$

<!-- 9 -->
$$
\int_{\partial\Omega}{v^\Re\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Re_\perp\vec{\tau}\;dS}=0
$$

<!-- 10 -->
$$
\int_{\partial\Omega}{v^\Re\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Im_\perp\vec{\tau}\;dS}=0
$$

<!-- 11 -->
$$
-\int_{\partial\Omega}{v^\Re\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;dS}
=\int_{\partial\Omega}{v^\Re\left(\frac{\nu^\Re}{\delta}\vec{E}_z^\Re\right)\;dS}
+\int_{\partial\Omega}{v^\Im\left(\frac{\nu^\Re}{\delta}\vec{E}_z^\Re\right)\;dS}
$$

<!-- 12 -->
$$
\int_{\partial\Omega}{v^\Re\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;dS}
=-\int_{\partial\Omega}{v^\Re\left(\frac{\nu^\Im}{\delta}\vec{E}_z^\Im\right)\;dS}
-\int_{\partial\Omega}{v^\Im\left(\frac{\nu^\Im}{\delta}\vec{E}_z^\Im\right)\;dS}
$$

<!-- 13 -->
$$
-\int_{\partial\Omega}{v^\Im\nu^\Re\hat{\operatorname{Z}}\vec{E}^\Re_\perp\;dS}
=0
$$

<!-- 14 -->
$$
\int_{\partial\Omega}{v^\Im\nu^\Im\hat{\operatorname{Z}}\vec{E}^\Im_\perp\;dS}
=0
$$

<!-- 15 -->
$$
-\int_{\partial\Omega}{v^\Im\left(\nu^\Im\hat{\operatorname{A}}\vec{E}_z^\Re\right)\;dS}
=\int_{\partial\Omega}{v^\Im\left(\frac{\nu^\Im}{\delta}\vec{E}_z^\Re\right)\;dS}
-\int_{\partial\Omega}{v^\Re\left(\frac{\nu^\Im}{\delta}\vec{E}_z^\Re\right)\;dS}
$$

<!-- 16! -->
$$
-\int_{\partial\Omega}{v^\Im\left(\nu^\Re\hat{\operatorname{A}}\vec{E}_z^\Im\right)\;dS}
=\int_{\partial\Omega}{v^\Im\left(\frac{\nu^\Re}{\delta}\vec{E}_z^\Im\right)\;dS}
-\int_{\partial\Omega}{v^\Re\left(\frac{\nu^\Re}{\delta}\vec{E}_z^\Im\right)\;dS}
$$

These boundary conditions are needed to simulate metal walls, so terms with imaginary permeability may be dropped.
Real permeability $1/\mu^\Re$ may be substituted with constant $1/\mu_0$.
