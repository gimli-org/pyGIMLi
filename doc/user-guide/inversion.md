---
file_format: mystnb
kernelspec:
  name: python3
---

# Inversion

## Theory

Inversion frameworks are generalized, abstract approaches to solve a specific inversion problem without specifying the appropriate geophysical methods.
This can be a specific regularization strategy, an alternative formulation of the inverse problem or algorithms of routine inversion.
It is initialized by specific forward operators or managers that provide them.

### Minimization

For all inversion frameworks the minimization of an objective function is required.
The most common approaches are based on least-squares minimization of a data misfit term between the individual data $d_i$ and the corresponding forwarc response $f_i(\mathbf{m})$ of the model $\mathbf{m})$.
In order to account for measuring errors, we weight the residuals with the inverse of the error $\epsilon_i$ so that the overall data objective function $\Phi_\text{d}$ reads:

$$ \Phi_\text{d} = \sum_{i=1}^{N} \left( \frac{f_i(\mathbf{m}) - d_i}{\epsilon_i} \right)^2 = \| \mathbf{W}_\text{d} (\mathbf{\mathcal{F}}(\mathbf{m})-\mathbf{d}) \|^2_2 $$

where $\mathbf{W}_\text{d}$ is the data weighting matrix containing the inverse data errors and $\mathbf{\mathcal{F}}(\mathbf{m})$ is the forward response vector for the model $\mathbf{m}$.
As the forward operator is in general non-linear, the minimization of the data misfit term requires iterative linearization approaches, starting from an initial model $\mathbf{m}^0$ and updating the model in each iteration $k$ by $\Delta \mathbf{m}^k$:

$$ \mathbf{m}^{k+1} = \mathbf{m}^k + \tau^k\Delta \mathbf{m}^k $$

The step length $\tau^k$ can be determined by line search strategies to ensure sufficient decrease of the objective function in each iteration.

For inversion of a limited number of parameters (e.g., 1D layer models, spectral parameters) the minimization is performed directly ($\Phi_d\rightarrow\min$) by some stabilizing damping terms, as with the Levenberg-Marquardt {cite}`Levenberg1944, Marquardt1963` method.
However, for large-scale problems (e.g., 2D/3D mesh-based inversion) the direct minimization is often not feasible.

### Regularization

To stabilize the inversion of ill-posed geophysical problems, additional constraints on the model parameters are required.
A common approach is to add a regularization term $\Phi_\text{m}$ to the objective function that penalizes undesirable model features.
A very common choice for 2D/3D problems is the roughness of the model distribution, but there is a wide range of different regularization methods (different kinds of smoothness and damping, mixed operators, anisotropic smoothing). We express this term by the matrix $\mathbf{W}_\text{m}$ acting on the model parameters $\mathbf{m}$ and possibly a reference model $\mathbf{m_0}$:

$$ \Phi=\Phi_\text{d}+\lambda\Phi_\text{m} = \mathbf{W}_\text{d} (\mathbf{\mathcal{F}}(\mathbf{m})-\mathbf{d}) \|^2_2 + \lambda \| \mathbf{W}_\text{m} (\mathbf{m}-\mathbf{m_0}) \|^2_2 \rightarrow\min $$ (eq:min)

The dimensionless factor $\lambda$ scales the influence of the regularization term $\Phi_\text{m}$ (model objective function).

For minimizing the total objective function $\Phi$, there are different available strategies:

- Steepest-descent methods
- Conjugate-gradient methods
- Gauss-Newton methods
- Quasi-Newton methods
- Stochastic methods (e.g., genetic algorithms, Markov-Chain Monte Carlo)

### Gauss-Newton inversion

The default inversion framework is based on the generalized Gauss-Newton minimization scheme leading to the model update $\Delta\mathbf{m}^k$ in the $k^\text{th}$ iteration {cite}`ParkVan1991`:

$$ ({\mathbf{J}}^{\text{T}}\mathbf{W_d}^{\text{T}}\mathbf{W_d}\mathbf{J} + \lambda {\mathbf{W_m}}^{\text{T}}\mathbf{W_m}
)\Delta\mathbf{m}^k =
{\mathbf{J}}^{\text{T}} {\mathbf{W_d}}^{\text{T}} \mathbf{W_d}(\Delta\mathbf{d}^k)
− \lambda{\mathbf{W_m}}^{\text{T}}\mathbf{W_m}(\mathbf{m}^k − \mathbf{m}^0 )
$$

with $\Delta \mathbf{d}^k = \mathbf{d} − \mathcal{F}(\mathbf{m}^k)$
and $\Delta \mathbf{m}^k = \mathbf{m}^k - \mathbf{m}^{k-1}
$

which is solved using a conjugate-gradient least-squares solver {cite}`Guenther2006`.
The inversion process including the region-specific regularization is sketched in Fig.~\ref{fig:InversionBase}.

\begin{figure}
\centering\includegraphics[width=1\columnwidth]{gimli-fig-2.pdf}
\caption{Generalized inversion scheme. Already implemented (\autoref{tab:methods}) or custom forward operators can be used that provide the problem specific response function and its Jacobian. Various strategies are available to regularize the inverse problem. \label{fig:InversionBase}}
\end{figure}

All matrices of the inversion formulation can be directly accessed from Python and thereby offer opportunities for uncertainty and resolution analysis as well as experimental design {cite}`{e.g., }Wagner2015`.
Beyond different inversion approaches there are so-called frameworks for typical inversion (mostly regularization) tasks.
Examples that are already implemented in pyGIMLi are for example:

- **Marquardt scheme** inversion of few independent parameters, e.g., fitting of spectra {cite}`Loewer2016`
- **Soil-physical model reduction** incorporating soil-physical functions {cite}`Igel2016, Costabel2014`
- **Classical joint inversion** of two data sets for the same parameter like DC and EM {cite}`Guenther2013NSG`
- **Block joint inversion** of several 1D data using common layers, e.g., MRS+VES {cite}`Guenther2012`
- **Sequential (constrained) inversion** successive independent inversion of data sets, e.g., classic time-lapse inversion {cite}`{e.g., }Bechtold2012`
- **Simultaneous constrained inversion** of data sets of data neighbored in space LCI, e.g., {cite}`Costabel2016`, time (full time-lapse) or frequency {cite}`Guenther2016`
- **Structurally coupled cooperative inversion** of disparate data based on structural similarity (e.g., {cite}`Ronczka2017`
- **Structure-based inversion** using layered 2D models {cite}`Attwa2014`

- lambda
- chi^2

+++

## Input data

### Data weights / errors

### Data transforms

+++

## Model parametrization

### Mesh-free inversion (0-D)

### Mesh inversion

#### 1-D

#### 2-D

#### 3-D

+++

## Regularization - Including prior information

### Starting model

### Reference model

### Parameter limits

### Damping

### Smoothing

### Advanced regularization

+++

## Region concept

+++

## Model appraisal

### Data misfit

### Cumulative sensitivity

### Resolution
