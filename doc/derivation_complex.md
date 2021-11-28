<!-- pandoc -f markdown -t html -s  --mathjax=https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js --metadata title='Derivation' -o derivation_complex.html derivation_complex.md -->

Based on paper [Space charge and resistive wall impedance computation in the frequency domain using
the finite element method by Uwe Niedermayer, Oliver Boine-Frankenheim, and Herbert De Gersem](https://doi.org/10.1103/PhysRevSTAB.18.032001).
Equations in this paper will be referred to as `Uwe:1`, `Uwe:2` etc.

---

Our goal is to obtain longitudinal and transverse impedance (Uwe:4-6)

$$\underline{Z}_\parallel(\omega) = -\frac{1}{q^2}\int_{beam}{\vec{\underline{E}}\underline{J}^{*}_\parallel\;dV}$$
$$\underline{Z}_\perp(\omega) = -\frac{\beta c}{(q d_r)^2\omega}\int_{beam}{\vec{\underline{E}}\underline{J}^{*}_\perp\;dV}$$

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

$$\nabla\times\frac{1}{\mu}\nabla\times(\vec{\underline{E}}_{curl}+\vec{\underline{E}}_{div}) = -j\omega J_s+\omega^2\underline{\varepsilon}(\vec{\underline{E}}_{curl}+\vec{\underline{E}}_{div})$$

$$\nabla\times\frac{1}{\mu}\nabla\times\vec{\underline{E}}_{curl} = -j\omega
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
$\underline{\Phi}$ is a scalar field.

For 2D problem gradient may be decomposed as follows.

$$\nabla\underline{\Phi}=
\frac{\partial\underline{\Phi}}{\partial x}\vec{x}+
\frac{\partial\underline{\Phi}}{\partial y}\vec{y}+
\frac{\partial\underline{\Phi}}{\partial z}\vec{z}=
\begin{pmatrix}
\frac{\partial}{\partial x} \\
\frac{\partial}{\partial y} \\
\end{pmatrix}\underline{\Phi}\vec{n}+
\frac{\partial}{\partial z}\underline{\Phi}\vec{z}=
\nabla_\perp\underline{\Phi}\vec{n}-j\frac{\omega}{\beta c_0}\underline{\Phi}\vec{e}_{z}
$$

Material properties are constant across elements.
So, in the calculations they are handled as if they were constant.

$$-\underline{\varepsilon}\nabla\cdot
\left(\nabla_\perp\underline{\Phi}\vec{n} - j \frac{\omega}{\beta c_0}\underline{\Phi}\vec{e}_{z}\right)
=\frac{1}{\beta c_0}J_s$$

$$-\underline{\varepsilon}\nabla\cdot\nabla_\perp\underline{\Phi}\vec{n}
+j\frac{\omega\underline{\varepsilon}}{\beta c_0}\nabla\underline{\Phi}\vec{e}_z
=\frac{1}{\beta c_0}J_s$$

$$
-\underline{\varepsilon}\nabla_\perp\cdot\nabla_\perp\underline{\Phi}\vec{n}
-\underline{\varepsilon}\frac{\partial}{\partial z}\nabla_\perp\underline{\Phi}\vec{n}\vec{e}_{z}
+j\frac{\omega\underline{\varepsilon}}{\beta c_0}\nabla_\perp\underline{\Phi}\vec{e}_{z}\vec{n}
+\frac{\omega^2\underline{\varepsilon}}{\beta^2 c_0^2}\underline{\Phi}\vec{e}_{z}
=\frac{1}{\beta c_0}J_s
$$

$\vec{n}$ and $\vec{e}_{z}$ are orthogonal, so their product $\vec{n} \vec{e}_{z} = 0$.

$$-\underline{\varepsilon}\nabla_\perp\cdot\nabla_\perp\underline{\Phi}\vec{n}
+\frac{\omega^2\underline{\varepsilon}}{\beta^2 c_0^2}\underline{\Phi}\vec{e}_{z}
=\frac{1}{\beta c_0}J_s$$

Using Galerkin method lowering equation order by multiplying equation by the test function $v$ and
[integrating by parts](https://en.wikipedia.org/wiki/Integration_by_parts#Higher_dimensions).
Writing down separately for each term.

$$
-\underline{\varepsilon}\int_\Omega{v\nabla_\perp^2\Phi \;d\Omega}
=\underline{\varepsilon}\int_\Omega{\nabla_\perp v \cdot \nabla_\perp\Phi \;d\Omega}
-\underline{\varepsilon}\int_{\partial\Omega}{v\nabla_\perp\Phi\vec{n}\;dS}
$$

$$
\frac{\omega^2\underline{\varepsilon}}{\beta^2 c_0^2}\int_\Omega{v\Phi \;d\Omega}
$$

$$
\frac{1}{\beta c_0}\int_\Omega{v J_s \;d\Omega}
$$

Using PEC boundary condition thus setting boundary field values to zero.
For trial functions $\Phi$ and test functions $v$
equation looks like this

$$
\underline{\varepsilon}\int_\Omega{\nabla_\perp v \cdot \nabla_\perp\Phi \;d\Omega}
+\frac{\omega^2\underline{\varepsilon}}{\beta^2 c_0^2}\int_\Omega{v\Phi \;d\Omega}
=\frac{1}{\beta c_0}\int_\Omega{v J_s \;d\Omega}
$$

After multiplying both sides of the equation by $\beta c_0$ we get

$$
\underline{\varepsilon}\beta c_0\int_\Omega{\nabla_\perp v \cdot \nabla_\perp\Phi \;d\Omega}
+\frac{\omega^2\underline{\varepsilon}}{\beta c_0}\int_\Omega{v\Phi \;d\Omega}
=\int_\Omega{v J_s \;d\Omega}
$$

After solving for potential we may get the irrotational part of electric field from the equations

$$E_{div\,\perp} = -\nabla\Phi$$
$$E_{div\,z} = j\frac{\omega}{\beta c_0}\Phi$$

## Solving for solenoidal field

Orienting coordinate system so that section slice is placed on XY plane.
Z axis is oriented according to the right hand rule.
Enclosing surface normal is opposite to Z axis ($\vec{e} = -\vec{n}$).
Section boundary normal points inward, boundary tangential vector is oriented clockwise relative to
the enclosing surface normal direction.

$$\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{E}=-j\omega{J}_S$$

$$\vec{\underline{E}} = \vec{\underline{E}}_{curl} + \vec{\underline{E}}_{div}$$

Since $\nabla\times\vec{\underline{E}}_{div}=0$

$$\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}}_{curl}
-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{curl}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{div}
=-j\omega{J}_S$$

$$\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}}_{curl}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{curl}
=\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{div}-j\omega{J}_S$$

$$\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}}_{curl}-\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{curl}=\vec{\underline{R}}$$

$$\vec{\underline{R}}=\omega^2\varepsilon_0\underline{\varepsilon}_r\vec{\underline{E}}_{div}-j\omega\underline{J}_{s}$$

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
\hat{\operatorname{Z}} &  & \hat{\operatorname{A}} \\
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

$$\hat{\operatorname{Z}}\vec{F}=
\frac{\partial\left(\vec{e}_y\cdot\vec{F}\right)}{\partial z}\vec{e}_x
-\frac{\partial\left(\vec{e}_x\cdot\vec{F}\right)}{\partial z}\vec{e}_y
=\frac{\partial}{\partial z}\vec{e}_z\times\vec{F}$$

Due to the uniformity of the structure in $\vec{e}_z$ direction
$$
\hat{\operatorname{Z}}\vec{F}
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
\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}}
=&
\nabla\times\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
\nabla\times\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z +
\nabla\times\frac{1}{\underline{\mu}}\hat{\operatorname{B}}\vec{\underline{E}}_\perp \\
=&
\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
\hat{\operatorname{A}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp + \\
&\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z +
\hat{\operatorname{A}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z +
\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z + \\
&\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{B}}\vec{\underline{E}}_\perp +
\hat{\operatorname{A}}\frac{1}{\underline{\mu}}\hat{\operatorname{B}}\vec{\underline{E}}_\perp +
\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{B}}\vec{\underline{E}}_\perp
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
\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}}
=&
\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{\underline{E}}_\perp +
\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z +
\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z +
\hat{\operatorname{A}}\frac{1}{\underline{\mu}}\hat{\operatorname{B}}\vec{\underline{E}}_\perp
\end{align}
$$

Grouping terms

$$
\begin{align}
\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{\underline{E}}
=
&\begin{pmatrix}
\hat{\operatorname{A}}\frac{1}{\underline{\mu}}\hat{\operatorname{B}}
+\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}} &
\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}} \\
\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}} &
\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}
\end{pmatrix}
\begin{pmatrix}
\vec{E}_\perp \\
E_z \\
\end{pmatrix}
\end{align}
$$

Let's rewrite the first equation term in matrix form

$$
\nabla\times\frac{1}{\underline{\mu}}\nabla\times\vec{E}
=[H]
\begin{pmatrix}
e_\perp \\
e_z \\
\end{pmatrix}
$$

$$
[H]=
\begin{pmatrix}
H_{\perp\perp} &
H_{\perp z} \\
H_{z\perp} &
H_{zz}
\end{pmatrix}
$$

Elements of matrix $[H]$ are:

<!-- 1 -->
$$H_{\perp\perp}\vec{E}_\perp=
\hat{\operatorname{A}}\frac{1}{\underline{\mu}}\hat{\operatorname{B}}\vec{E}_\perp
+\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{E}_\perp$$

<!-- 2 -->
$$H_{\perp z}E_z=\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z$$

<!-- 3 -->
$$H_{z\perp}\vec{E}_\perp=\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\vec{E}_\perp$$

<!-- 4 -->
$$H_{zz}{E}_z=\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}{E}_z$$

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

#### Converting term $\hat{\operatorname{A}}\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp$
into [weak
form](https://en.wikipedia.org/wiki/Weak_formulation) by multiplying it by the test function
$\vec{w}$ and
integrating over the calculation domain $\Omega$.

$$
\begin{align}
\int_\Omega{
\vec{w}
\left(
\hat{\operatorname{A}}\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
=&\int_\Omega{
\vec{w}
\left(
-\vec{e}_z\times\nabla\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
} \\
=&-\vec{e}_z\times
\int_\Omega{
\vec{w}
\left(
\nabla\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
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
\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
-\vec{e}_z\times
\int_{\partial\Omega}{
\vec{w}
\left(
\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
\right)
\vec{n}
\;dS
} \\
=&\int_\Omega{
\left(
\vec{e}_z\times\nabla\vec{w}
\right)^T
\left(
\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
-\int_{\partial\Omega}{
\vec{w}
\left(
\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
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
\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
\right)
\;d\Omega
}
+\int_{\partial\Omega}{
\vec{w}
\left(
\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp
\right)
\vec{\tau}
\;dS
}
\end{align}
$$

#### Converting term $\hat{\operatorname{B}}\frac{1}{\mu}\hat{\operatorname{A}}E_z$

into [weak
form](https://en.wikipedia.org/wiki/Weak_formulation) by multiplying it by the test function
$v$ and
integrating over the calculation domain $\Omega$.

$$
\begin{align}
\int_\Omega{
v
\left(
\hat{\operatorname{B}}\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\;d\Omega
}
=&\int_\Omega{
v
\left(
\vec{e}_z\nabla\times\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\;d\Omega
} \\
=&-\int_\Omega{
v
\left(
\nabla\times\frac{1}{\mu}\hat{\operatorname{A}}E_z
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
\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\vec{n}
\;d\Omega
}
-\int_{\partial\Omega}{
v
\left(
\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\vec{\tau}
\;dS
} \\
=&\int_\Omega{
\left(
-\vec{e}_z\times\nabla v
\right)
\left(
\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\;d\Omega
}
-\int_{\partial\Omega}{
v
\left(
\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\vec{\tau}
\;dS
} \\
=&\int_\Omega{
\left(
\hat{\operatorname{A}}v
\right)
\left(
\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\;d\Omega
}
-\int_{\partial\Omega}{
v
\left(
\frac{1}{\mu}\hat{\operatorname{A}}E_z
\right)
\vec{\tau}
\;dS
}
\end{align}
$$

We may now use decomposed curl operator identities to convert $[H]$ matrix into weak form.
Using Galerkin method lowering equations order by multiplying them by the test functions and
[integrating by parts](https://en.wikipedia.org/wiki/Integration_by_parts#Higher_dimensions) (where necessary).
Test functions for transverse and longitudinal fields $\vec{w}$ and $v$ respectively.

$$
[\hat{H}]=
\begin{pmatrix}
\hat{H}_{\perp\perp} &
\hat{H}_{\perp z} \\
\hat{H}_{z\perp} &
\hat{H}_{zz}
\end{pmatrix}
$$

<!-- 1 -->
$$\begin{align}
\hat{H}_{\perp\perp}\underline{\vec{E}}_\perp
&=\int_{\Omega}{\vec{w}H_{\perp\perp}\underline{\vec{E}}_\perp\;d\Omega} \\
&=\int_{\Omega}{\vec{w}\hat{\operatorname{A}}\frac{1}{\underline{\mu}}\hat{\operatorname{B}}\underline{\vec{E}}_\perp\;d\Omega}
+\int_{\Omega}{\vec{w}\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{Z}}\underline{\vec{E}}_\perp\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}\right)\left(\frac{1}{\mu}\hat{\operatorname{B}}\underline{\vec{E}}_\perp\right)\;d\Omega}
+\int_{\partial\Omega}{\vec{w}\left(\frac{1}{\mu}\hat{\operatorname{B}}\underline{\vec{E}}_\perp\right)\vec{\tau}\;dS}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}\frac{1}{\mu}\underline{\vec{E}}_\perp\;d\Omega}
\end{align}$$

<!-- 2 -->
$$\hat{H}_{\perp z}E_z
=\int_{\Omega}{\vec{w}H_{\perp z}E_z\;d\Omega}
=\int_{\Omega}{\vec{w}\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}E_z\;d\Omega}$$

<!-- 3 -->
$$\begin{align}
\hat{H}_{z\perp}\underline{\vec{E}}_\perp
&=\int_{\Omega}{vH_{z\perp}\underline{\vec{E}}_\perp\;d\Omega}
=\int_\Omega{v\hat{\operatorname{B}}\frac{1}{\mu}\hat{\operatorname{Z}}\underline{\vec{E}}_\perp\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{A}}v\right)
\left(\frac{1}{\mu}\hat{\operatorname{Z}}\underline{\vec{E}}_\perp\right)\;d\Omega}
-\int_{\partial\Omega}{v\frac{1}{\mu}\hat{\operatorname{Z}}\underline{\vec{E}}_\perp\vec{\tau}\;dS}
\end{align}$$

<!-- 4 -->
$$\begin{align}
\hat{H}_{zz}\underline{E}_z
&=\int_{\Omega}{v\hat{\operatorname{B}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}\underline{E}_z\;d\Omega} \\
&=\int_\Omega{\left(\hat{\operatorname{A}} v\right)
\left(\frac{1}{\mu}\hat{\operatorname{A}}\underline{E}_z\right)\;d\Omega}
-\int_{\partial\Omega}{v\left(\frac{1}{\mu}\hat{\operatorname{A}}\underline{E}_z\right)\vec{\tau}\;dS}
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

Let's split these terms in longitudinal and transverse parts

$$\omega^2\underline{\varepsilon}\vec{\underline{E}}
=\omega^2\underline{\varepsilon}\vec{\underline{E}}_\perp
+\omega^2\underline{\varepsilon}\underline{E}_z$$

$J_s$ is a source term


The whole matrix is (Uwe:37-40)
$$
\left[
\begin{pmatrix}
S_{\perp\perp} &
S_{\perp z} \\
S_{z \perp} &
S_{zz}
\end{pmatrix}
+
\begin{pmatrix}
M_{\perp} & 0 \\
0 & M_{z}
\end{pmatrix}
+
\begin{pmatrix}
B_{\perp\perp} & 0 \\
B_{z\perp} & B_{zz}
\end{pmatrix}
\right]
E_{curl}=R
$$

$E_{curl}$ terms

<!-- 1 -->
$$\hat{S}_{\perp\perp}
=\int_\Omega{\left(\hat{\operatorname{B}}\vec{w}\right)\left(\frac{1}{\mu}\hat{\operatorname{B}}\underline{\vec{E}}_\perp\right)\;d\Omega}
+\frac{\omega^2}{\beta^2 c_0^2}\int_\Omega{\vec{w}\frac{1}{\mu}\underline{\vec{E}}_\perp\;d\Omega}$$

<!-- 2 -->
$$\hat{S}_{\perp z}
=\int_{\Omega}{\vec{w}\hat{\operatorname{Z}}\frac{1}{\underline{\mu}}\hat{\operatorname{A}}\underline{E}_z\;d\Omega}$$

<!-- 3 -->
$$\hat{S}_{z\perp}
=\int_\Omega{\left(\hat{\operatorname{A}}v\right)
\left(\frac{1}{\mu}\hat{\operatorname{Z}}\underline{\vec{E}}_\perp\right)\;d\Omega}$$

<!-- 4 -->
$$\hat{S}_{zz}\underline{E}_z
=\int_\Omega{\left(\hat{\operatorname{A}} v\right)
\left(\frac{1}{\mu}\hat{\operatorname{A}}\underline{E}_z\right)\;d\Omega}$$

<!-- 1 -->
$$
M_{\perp}=
-\omega^2\int_{\Omega}{\vec{w}\underline{\varepsilon}\underline{\vec{E}}_\perp\;d\Omega}
$$

<!-- 4 -->
$$
M_{z}=
-\omega^2\int_{\Omega}{v\underline{\varepsilon}\underline{E}_z\;d\Omega}
$$

<!-- 1 -->
$$B_{\perp\perp}
=\int_{\partial\Omega}{\vec{w}\left(\frac{1}{\mu}\hat{\operatorname{B}}\underline{\vec{E}}_\perp\right)\vec{\tau}\;dS}$$

<!-- 3 -->
$$
B_{z\perp}
=-\int_{\partial\Omega}{\vec{w}\frac{1}{\mu}\hat{\operatorname{Z}}\underline{\vec{E}}_\perp\vec{\tau}\;dS}
$$

<!-- 4 -->
$$
B_{zz}
=-\int_{\partial\Omega}{v\left(\frac{1}{\mu}\hat{\operatorname{A}}\underline{E}_z\right)\vec{\tau}\;dS}
$$

$E_{div}$ terms

<!-- 1 -->
$$
N_{\perp}=
\omega^2\int_{\Omega}{\vec{w}\underline{\varepsilon}\underline{\vec{E}}_\perp\;d\Omega}
$$

<!-- 4 -->
$$
N_{z}=
\omega^2\int_{\Omega}{v\underline{\varepsilon}\underline{E}_z\;d\Omega}
$$

$J_{s}$ term

$$
J_s =
-j\omega\int_{\Omega}{v \underline{J}_s\;d\Omega}
$$

Now we have all the equations to obtain $\vec{\underline{E}}=\vec{\underline{E}}_{curl}+\vec{\underline{E}}_{div}$ for
the chosen source function $J_{S}$.
These can be inserted into the equations for impedance (Uwe:4-6)
$$\underline{Z}_\parallel(\omega)
=-\frac{1}{q^2}\int_{beam}{\vec{\underline{E}}\underline{J}^{*}_\parallel\;dV}$$

For the purely real source function ${J}$ and a step function
$$
\delta=
\begin{cases}
1 \quad (x,y)\in\Omega_{beam} \\
0 \quad (x,y)\notin\Omega_{beam}
\end{cases}
$$
$$\underline{Z}_\parallel(\omega)
=-j\frac{1}{q^2}\int_{\Omega}{\underline{\vec{E}} {J}_\parallel \delta \;d\Omega}$$

Similarly
$$\underline{Z}_\perp(\omega)
=-\frac{\beta c}{(q d_r)^2\omega}\int_{beam}{\vec{\underline{E}}\underline{J}^{*}_\perp\;dV}
=-j\frac{\beta c}{(q d_r)^2\omega}\int_{beam}{\vec{\underline{E}}J_\perp\;dV}
$$

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
$$\hat{Z}\vec{\underline{E}}_\perp+\hat{A}\vec{\underline{E}}_\perp+\hat{B}\underline{E}_z
=-j\omega\mu\vec{\underline{H}}$$
Extracting tangential and $\vec{e}_z$ field parts
$$\hat{Z}\vec{\underline{E}}_\tau+\hat{A}\underline{E}_z
=-j\omega\mu\vec{\underline{H}}_\tau$$
$$\hat{B}\vec{\underline{E}}_\tau=-j\omega\mu\underline{H}_z$$
Noticing that term $\hat{Z}\vec{E}_\tau$ is normal to the boundary and therefore will not contribute to the boundary integral.
$$\begin{align}
\hat{A}\underline{E}_z=-j\omega\mu\vec{\underline{H}}_\tau
&=-\frac{j\omega\mu}{(1+j)\sqrt{\frac{\mu\omega}{2\sigma}}}\underline{E}_z
=-\frac{j\omega\mu\sqrt{2\sigma}}{(1+j)\sqrt{\omega\mu}}\underline{E}_z
=-\frac{2j\sqrt{\omega\mu\sigma}}{(1+j)\sqrt{2}}\underline{E}_z
=-\frac{2j}{(1+j)\delta}\underline{E}_z \\
&=-\frac{2j(1-j)}{(1+j)(1-j)\delta}\underline{E}_z
=-\frac{2j+2}{(1^2-j^2)\delta}\underline{E}_z
=-\frac{1+j}{\delta}\underline{E}_z
\end{align}$$
$$\hat{B}\vec{\underline{E}}_\perp=-j\omega\mu\underline{H}_z
=-\frac{1+j}{\delta}\vec{\underline{E}}_\perp$$
where $\delta=\sqrt{\frac{2}{\omega\mu\sigma}}$ is a skin depth.
Now, applying these equations to $[M_{SIBC}]$ matrix terms.

<!-- 1 -->
$$B_{\perp\perp}
=\int_{\partial\Omega}{\vec{w}\left(\frac{1}{\mu}\hat{\operatorname{B}}\vec{E}_\perp\right)\vec{\tau}\;dS}
=-\int_{\partial\Omega}{\left(\vec{w}\vec{\tau}\right)\left(\frac{1+j}{\delta\mu}\vec{E}_\perp\vec{\tau}\right)\;dS}$$

<!-- 3 -->
$$
B_{z\perp}
=-\int_{\partial\Omega}{\vec{w}\frac{1}{\mu}\hat{\operatorname{Z}}\vec{E}_\perp\vec{\tau}\;dS}=0
$$

<!-- 4 -->
$$
B_{zz}
=-\int_{\partial\Omega}{v\left(\frac{1}{\mu}\hat{\operatorname{A}}{E}_z\right)\vec{\tau}\;dS}
=\int_{\partial\Omega}{v\frac{1+j}{\delta\mu}\underline{E}_z\;dS}
$$

These boundary conditions are needed to simulate metal walls, so permeability $1/\mu$ may be
substituted with constant $1/\mu_0$.
