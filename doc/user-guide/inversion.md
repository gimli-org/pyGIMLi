---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
kernelspec:
  language: python
  name: python3
---

# Inversion

+++

## Theory
Inversion frameworks are generalized, abstract approaches to solve a specific inversion problem without specifying the appropriate geophysical methods.
This can be a specific regularization strategy, an alternative formulation of the inverse problem or algorithms of routine inversion.
It is initialized by specific forward operators or managers that provide them.

### Gauss-Newton inversion

The default inversion framework is based on the generalized Gauss-Newton method and is compatible with any given forward operator and thus applicable to various physical problems.
We state the inversion problem as minimization of an objective function consisting of data misfit and model constraints:

$$ \mathbf{W}_\text{d} (\mathbf{\mathcal{F}}(\mathbf{m})-\mathbf{d}) \|^2_2 + \lambda \| \mathbf{W}_\text{m} (\mathbf{m}-\mathbf{m_0}) \|^2_2 \rightarrow\min $$ (eq:min)

Note that we do not include inequality constraints in the minimization but use transformations to restrict parameters to reasonable ranges {cite}`{e.g., }KimKim2011`.
$\mathbf{W}_\text{d}$ is the data weighting matrix containing the inverse data errors, $\mathbf{W}_\text{m}$ is the model constraint matrix (e.g., a first-order roughness operator), and $\mathbf{m}_0$ is a reference model.
The dimensionless factor $\lambda$ scales the influence of the regularization term.
There is a wide range of different regularization methods (different kinds of smoothness and damping, mixed operators, anisotropic smoothing)
The existing ones can be used flexibly to constrain different model parameters or subsurface parts (regions), but also be extended by own functions.
The application of the Gauss-Newton scheme on minimizing {eq}`eq:min` yields the model update $\Delta\mathbf{m}^k$ in the $k^\text{th}$ iteration {cite}`ParkVan1991`:

$$ ({\mathbf{J}}^{\text{T}}\mathbf{W_d}^{\text{T}}\mathbf{W_d}\mathbf{J} + \lambda {\mathbf{W_m}}^{\text{T}}\mathbf{W_m}
)\Delta\mathbf{m}^k =\: &
{\mathbf{J}}^{\text{T}} {\mathbf{W_d}}^{\text{T}} \mathbf{W_d}(\Delta\mathbf{d}^k) \notag \\
−& \lambda{\mathbf{W_m}}^{\text{T}}\mathbf{W_m}(\mathbf{m}^k − \mathbf{m}^0 ) \\
\quad\text{with}\quad \Delta \mathbf{d}^k = \mathbf{d} − \mathcal{F}(\mathbf{m}^k)
\quad\text{and}&\quad\Delta \mathbf{m}^k = \mathbf{m}^k - \mathbf{m}^{k-1}\notag
$$

which is solved using a conjugate-gradient least-squares solver {cite}`GuentherRueSpi2006`.
The inversion process including the region-specific regularization is sketched in Fig.~\ref{fig:InversionBase}.

\begin{figure}
\centering\includegraphics[width=1\columnwidth]{gimli-fig-2.pdf}
\caption{Generalized inversion scheme. Already implemented (\autoref{tab:methods}) or custom forward operators can be used that provide the problem specific response function and its Jacobian. Various strategies are available to regularize the inverse problem. \label{fig:InversionBase}}
\end{figure}

All matrices of the inversion formulation can be directly accessed from Python and thereby offer opportunities for uncertainty and resolution analysis as well as experimental design {cite}`{e.g., }WagnerGueSchMau2015GEO`.
Beyond different inversion approaches there are so-called frameworks for typical inversion (mostly regularization) tasks.
Examples that are already implemented in pyGIMLi are for example:

 * **Marquardt scheme** inversion of few independent parameters, e.g., fitting of spectra {cite}`LoewerIgeWag2016`
 * **Soil-physical model reduction** incorporating soil-physical functions {cite}`IgelStaGue2016SAGEEP, CostabelGue2014VZJ`
 * **Classical joint inversion** of two data sets for the same parameter like DC and EM {cite}`Guenther2013NSG`
 * **Block joint inversion** of several 1D data using common layers, e.g., MRS+VES {cite}`GuentherMue2012HESS`
 * **Sequential (constrained) inversion** successive independent inversion of data sets, e.g., classic time-lapse inversion {cite}`{e.g., }BechtoldVanWei2012`
 * **Simultaneous constrained inversion** of data sets of data neighbored in space LCI, e.g., {cite}`CostabelGueDluMue2016`, time (full time-lapse) or frequency {cite}`guenther2016JoAG`
 * **Structurally coupled cooperative inversion** of disparate data based on structural similarity (e.g., {cite}`Ronczka_2016`
 * **Structure-based inversion** using layered 2D models {cite}`AttwaAkcBas2014`










- Gauss-Newton formulation as in the paper
- lambda
- chi^2
- Terms should be explained in here

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
