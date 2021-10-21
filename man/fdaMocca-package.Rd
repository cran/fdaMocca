\name{fdaMocca-package}
%\Rdversion{0.1-0}
\alias{fdaMocca-package}
\docType{package}
\title{Model-based clustering for functional data with covariates
}
\description{
\code{fdaMocca} provides functions for model-based functional cluster analysis for functional data with optional covariates. The aim is to cluster a set of independent functional subjects (often called functional objects) into homogenous groups by using basis function representation of the functional data and allowing scalar covariates. A functional subject is defined as a curve and covariates. The spline coefficients and the (potential) covariates are modelled as a multivariate Gaussian mixture model, where the number of mixtures corresponds to the number of (predefined) clusters. The model allows for different cluster covariance structures for the basis coefficients and for the covariates. The parameters of the model are estimated by maximizing the observed mixture likelihood using an EM-type algorithm (Arnqvist and Sjöstedt de Luna, 2019).
}
\details{
\tabular{ll}{
Package: \tab fdaMocca\cr
Type: \tab Package\cr
License: \tab GPL (>= 2)\cr
}}
\author{
Per Arnqvist, Sara \enc{Sjöstedt}{Sjostedt} de Luna, Natalya Pya Arnqvist

Maintainer: Natalya Pya Arnqvist<nat.pya@gmail.com>

}
%\note{
% } 
\references{
Arnqvist, P., Bigler, C., Renberg, I., \enc{Sjöstedt}{Sjostedt} de Luna, S. (2016). Functional clustering of varved lake sediment to reconstruct past seasonal climate. \emph{Environmental and Ecological Statistics}, \bold{23}(4), 513--529.

Abramowicz, K., Arnqvist, P., Secchi, P., \enc{Sjöstedt}{Sjostedt} de Luna, S., Vantini, S., Vitelli, V. (2017). Clustering misaligned dependent curves applied to varved lake sediment for climate reconstruction. \emph{Stochastic Environmental Research and Risk Assessment}. Volume \bold{31.1}, 71--85.

Arnqvist, P., and \enc{Sjöstedt}{Sjostedt} de Luna, S. (2019). Model based functional clustering of varved lake sediments. \emph{arXiv preprint arXiv:1904.10265}.
}

