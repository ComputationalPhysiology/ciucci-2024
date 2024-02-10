---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
title: Simulating mechanical forces
description: Simple example of a contracting cylinder in FEniCS
keywords: ["cardiac mechancis", "elasticity", "stress"]
date: "2024-01-09"
license: CC-BY-4.0
authors:
  - name: Henrik Finsberg
    affiliations:
      - Simula Research Laboratory
  - name: Samuel Wall
    affiliations:
      - Simula Research Laboratory
github: ComputationalPhysiology/ciucci-2024
exports:
  - format: pdf
    template: curvenote
    output: report.pdf
  - format: tex
    template: curvenote
    output: report.tex
---
## Meshes and microstructure


### Cylindrical domain

We consider a cylindrical domain $\Omega(a, r) = [-a, a] \times \omega(r) $ with

```{math}
\omega(r) = \{ (y, z) \in \mathbb{R}^2 : y^2 + z^2 < r \}
```

We let $a = 2000 \mu m$ and $r = 300 \mu m$ and will refer to $\Omega_0$ as $\Omega(300, 2000)$


```{figure} ./figures_static/cylinder/mesh.png
:name: cyl_domain
:alt: Cylindrical domain
:align: center

Cylindrical mesh with $a = 2000 \mu m$ and $r = 300 \mu m$
```

For the cylinder we also define a local coordinate system in the longitudinal (aka fiber) direction

In our case we orient the fibers along the height the cylinder, i.e

```{math}
\mathbf{f}_0 = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}
```


```{figure} ./figures_static/cylinder/long.png
:name: long
:alt: long
:align: center

Showing the orientation of the cardiac fibers for the cylinder along the longitudinal direction
```

In a cylinder there is a natural decomposition into radial and circumferential components. We have the radial basis function given by

```{math}
\mathbf{e}_r = \begin{pmatrix}
0 \\ y / (y^2 + z^2) \\ z / (y^2 + z^2)
\end{pmatrix}
```

```{figure} ./figures_static/cylinder/rad.png
:name: radial
:alt: Radial
:align: center

Showing the orientation of radial components
```

and the circumferential basis function given by

```{math}
\mathbf{e}_c = \begin{pmatrix}
0 \\ z / (y^2 + z^2) \\ -y / (y^2 + z^2)
\end{pmatrix}
```

```{figure} ./figures_static/cylinder/circ.png
:name: circ
:alt: Circ
:align: center

Showing the orientation of circular components
```



### LV domain

We will also test out an idealized LV domain generated using [`cardiac-geometries`](https://computationalphysiology.github.io/cardiac_geometries)

```{figure} ./figures_static/lv/mesh.png
:name: lv_mesh
:alt: LV mesh
:align: center

Idealized LV mesh with width of 1.5 mm, a short axis radius of 4 mm and a long axis radius of 8 mm.
```

For the LV we also have the rule based fiber orientations

```{figure} ./figures_static/lv/fiber.png
:name: lv_fiber
:alt: LV fiber
:align: center

Fiber orientations in the LV with and endo / epi fiber angle of 60/-60.
```


```{figure} ./figures_static/lv/sheet.png
:name: lv_sheet
:alt: LV sheet
:align: center

Sheet orientations in the LV.
```

```{figure} ./figures_static/lv/sheet-normal.png
:name: lv_sheet_norma
:alt: LV sheet normal
:align: center

Sheet-normal orientations in the LV.
```




## Material model

We use the the transversely isotropic version of the [Holzapfel Ogden model](https://doi.org/10.1098/rsta.2009.0091), i.e

```{math}
  \Psi(\mathbf{F})
  = \frac{a}{2 b} \left( e^{ b (I_1 - 3)}  -1 \right)
  + \frac{a_f}{2 b_f} \mathcal{H}(I_{4\mathbf{f}_0} - 1)
  \left( e^{ b_f (I_{4\mathbf{f}_0} - 1)_+^2} -1 \right)
  + \kappa (J \mathrm{ln}(J) - J + 1)
```

with

```{math}
(x)_+ = \max\{x,0\}
```

and

```{math}
\mathcal{H}(x) = \begin{cases}
    1, & \text{if $x > 0$} \\
    0, & \text{if $x \leq 0$}
\end{cases}
```

is the Heaviside function. Here

```{math}
I_1 = \mathrm{tr}(\mathbf{F}^T\mathbf{F})
```

and

```{math}
I_{4\mathbf{f}_0} = (\mathbf{F} \mathbf{f}_0)^T \cdot \mathbf{F} \mathbf{f}_0
```

and
```{math}
J = \mathrm{det}(\mathbf{F})
```

with $\mathbf{f}_0$ being the direction the muscle fibers.

### Material parameters

The material parameter are

| Parameter | Value       |
|-----------|-------------|
| $a$       | 2280 Pa     |
| $b$       | 9.726       |
| $a_f$     | 1685 Pa     |
| $b_f$     | 15.779      |

### Modeling of compressibility




## Modeling of active contraction

Similar to [Finsberg et. al](doi.org/10.1002/cnm.2982) we use an active strain formulation and decompose the deformation gradient into an active and an elastic part

```{math}
\mathbf{F} = \mathbf{F}_e \mathbf{F}_a
```

with

```{math}
\mathbf{F}_a = (1 - \gamma) \mathbf{f} \otimes \mathbf{f}  + \frac{1}{\sqrt{1 - \gamma}} (\mathbf{I} - \mathbf{f} \otimes \mathbf{f}).
```

In these experiments we use $\gamma = 0.2$ to represent end systole.

## Variational formulation



We model the myocardium as incompressible using a two field variational approach and $\mathbb{P}_2-\mathbb{P}_1$ finite elements for the displacement $\mathbf{u}$ and hydrostatic pressure $p$.
m

The Euler‐Lagrange equations in the Lagrangian form reads: find $(\mathbf{u},p) \in H^1(\Omega_0) \times L^2(\Omega_0)$ such that for all $(\delta \mathbf{u},\delta \mathbf{u}) \in H^1(\Omega_0) \times L^2(\Omega_0)$ we have

```{math}
\delta \Pi(\mathbf{u}, p) = 0.
```

For the cylinder we have

```{math}
\delta \Pi(\mathbf{u}, p) = \int_{\Omega_0}\left[ \mathbf{P} : \nabla \delta \mathbf{u} - \delta p (J - 1) - pJ \mathbf{F}^{-T}: \delta \mathbf{u} \right] \mathrm{d}V + \int_{\partial \Omega_0^{a} \cup \Omega_0^{-a}} k \mathbf{u} \cdot \delta \mathbf{u}  \mathrm{d}S
```
where $ \Omega_0^{a} = \{ a \} \times \omega(r) $ represents the boundaries at each end of the cylinder. Here we enforce a Robin type boundary condition at both ends of the cylinder with a spring $k$. This parameter will be used to represent the loading conditions and we choose $k=$ 1 Pa / $\mu$m to represent standard loading conditions.

For the left ventricle we have


```{math}
\delta \Pi(\mathbf{u}, p) = \int_{\Omega_0}\left[ \mathbf{P} : \nabla \delta \mathbf{u} - \delta p (J - 1) - pJ \mathbf{F}^{-T}: \delta \mathbf{u} \right] \mathrm{d}V  + \int_{\partial \Omega_0^{\mathrm{epi}}} k \mathbf{u} \cdot \delta \mathbf{u}  \mathrm{d}S + \int_{\partial \Omega_0^{\mathrm{endo}}} p_{\mathrm{lv}} J \mathbf{F}^T \mathbf{N} \cdot \delta \mathbf{u}  \mathrm{d}S.
```


Here we enforce a Robin type boundary condition at the epicardium, $ \Omega_0^{\mathrm{epi}}$, with a spring $k = 0.5$ kPa / mm, and an endocardial pressure $p_{\mathrm{lv}}$ on the endocardium $ \Omega_0^{\mathrm{endo}}$. In addition we fix the basal plane in the longitudinal direction.



## Cauchy stress

The Cauchy stress tensor is given by

```{math}
\sigma = J^{-1} \frac{\partial \Psi }{\partial \mathbf{F}}  \mathbf{F}^T
```

We can extract different components of the Cauchy stress tensor. For example the fiber component can be extracted  using the following formula

```{math}
\sigma_{ff} = \mathbf{f}^T \cdot \sigma \mathbf{f}
```
where $\mathbf{f}$ is the vector field displayed in {ref}`lv_fiber` in the current configuration.


## Numerical experiments

### Cylinder

For the experiments with the cylinder we simulated two states; one relaxed state and one contracted state. For the contracted state we used $\gamma = 0.2$

### Cylinder varying spring

We tested different values of the spring constant $k$, with $k = 3$ Pa / $\mu$m, representing standard load, $k = 0.1$ Pa / $\mu$m representing unloaded and $k = 10$ kPa / $\mu$m representing increased load.

```{figure} figures/cylinder/diameter.svg
:name: spring_diam
:alt: spring Diameter
:align: center

Resulting diameter at the center for the cylinder in relaxed and contracted state for spring constants
```

```{figure} figures/cylinder/length.svg
:name: spring_length
:alt: spring Length
:align: center

Resulting length of the cylinder in relaxed and contracted state for spring constants
```


```{figure} figures/cylinder/frac.svg
:name: spring_length
:alt: spring Length
:align: center

Fractional shortening
```


```{figure} figures/cylinder/strain.svg
:name: spring_strain
:alt: spring Strain
:align: center

Resulting strain the cylinder in the contracted state for different spring constants in the longitudinal $E_{xx}$, circumferential $E_{c}$ and radial $E_r$ direction
```

```{figure} figures/cylinder/stress.svg
:name: spring_stress
:alt: spring Stress
:align: center

Resulting stress the cylinder in the contracted state for different spring constants in the longitudinal $\sigma_{xx}$, circumferential $\sigma_{c}$ and radial $\sigma_r$ direction
```


```{figure} figures/cylinder/stress_dev.svg
:name: spring_stress
:alt: spring Stress
:align: center

Resulting deviatoric and hydrostatic stress the cylinder in the contracted state for different spring constants in the longitudinal $\sigma_{xx}$, circumferential $\sigma_{c}$ and radial $\sigma_r$ direction
```

## Twitch

To simulate twitch we create a syntetic curve for $\gamma$,

```{math}
\gamma(t) = \begin{cases}
\gamma_{\mathrm{min}} & \text{ if } t <= t_0, \\
\frac{\gamma_{\mathrm{max}} - \gamma_{\mathrm{min}}}{\beta} \left( e^{-\frac{t - t_0}{\tau_1}} - e^{-\frac{t - t_0}{\tau_2}} \right) + \gamma_{\mathrm{min}}  & \text{ otherwise} \\
\end{cases}
```

and

```{math}
\beta = \left(\frac{\tau_1}{\tau_2}\right)^{- \frac{1}{\frac{\tau_1}{\tau_2} - 1}} - \left(\frac{\tau_1}{\tau_2}\right)^{- \frac{1}{1 - \frac{\tau_2}{\tau_1}}}.
```



```{figure} figures/cylinder-twitch/diameter.svg
:name: cylinder_twitch_diam
:alt: cylinder_twitch Diameter
:align: center

Resulting diameter at the center for the cylinder in relaxed and contracted state for spring constants
```

```{figure} figures/cylinder-twitch/length.svg
:name: cylinder_twitch_length
:alt: cylinder_twitch Length
:align: center

Resulting length of the cylinder in relaxed and contracted state for spring constants
```


```{figure} figures/cylinder-twitch/frac.svg
:name: cylinder_twitch_length
:alt: cylinder_twitch Length
:align: center

Fractional shortening
```


```{figure} figures/cylinder-twitch/stress.svg
:name: cylinder_twitch_stress
:alt: cylinder_twitch Stress
:align: center

Resulting stress the cylinder in the contracted state for different spring constants in the longitudinal $\sigma_{xx}$, circumferential $\sigma_{c}$ and radial $\sigma_r$ direction
```

```{figure} figures/cylinder-twitch/stress_polar_label.svg
:name: cylinder_twitch_stress_label
:alt: cylinder_twitch Stress label
:align: center

Label corresponding to figure {ref}`cylinder_twitch_stress` showing which regions we average the stress
```





## LV

For the LV we have two different states; one where we have no endocardial pressure (referred to as the unloaded state) and one where we apply an endocardial pressure of 15 kPa (referred to as the loaded state). In both cases we set $\gamma = 0.2$.

```{figure} figures/lv/strain.svg
:name: lv_strain
:alt: LV strain
:align: center

Strain in the fiber, sheet and sheet normal direction for the unloaded and loaded state
```

```{figure} figures/lv/stress.svg
:name: lv_stress
:alt: LV stress
:align: center

Stress in the fiber, sheet and sheet normal direction for the unloaded and loaded state
```


```{figure} figures/lv/stress_dev.svg
:name: lv_stress
:alt: LV stress
:align: center

Stress in the fiber, sheet and sheet normal direction for the unloaded and loaded state
```
