\name{varve}
\alias{varve}
\alias{varve_full}
\docType{data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Varved sediment data from lake Kassjön}

\description{Annually varved lake sediment data from the lake Kassjön in Northern Sweden. The Kassjön data are used to illustrate the ideas of the model-based functional clustering with covariates. The varved sediment of lake Kassjön covers approximately 6400 years and is believed to be historical records of weather and climate.  The varves (years) are clustered into similar groups based on their seasonal patterns (functional data) and additional covariates, all potentially carrying information on past climate/weather. 

The Kassjön data has been analyzed in several papers. In Petterson et al. (2010, 1999, 1993) the sediment data was captured with image analysis and after preprocessing, the data was recorded as gray scale values with yearly deliminators, thus giving 6388 years (-4386 -- 1901), or varves with 4--36 gray scale values per year. In Arnqvist et al. (2016) the shape/form of the
yearly grey scale observations was modeled as curve functions and analyzed in a non-parametric functional data analysis approach. In Abramowicz et al. (2016) a Bagging Voronoi K-Medoid Alignment (BVKMA) method was proposed to group the varves into different "climates". The suggested procedure simultaneously clusters and aligns spatially dependent curves and is a nonparametric statistical method that does not
rely on distributional or dependency structure assumptions.

}
%\usage{
%  load("covariates.rda")
%}
%- maybe also `usage' for other objects documented here.

\format{ \code{varve} data is a list containing six objects named as \code{x}, \code{time}, \code{timeindex}, \code{curve}, \code{grid}, \code{covariates}. See \code{\link{mocca}} for explanation of these objects. 

\code{varve_full} has \eqn{N=6326} observed subjects (years/varve), where for each varve we observed one function (seasonal pattern) and four covariates. \code{varve} is simply a reduced data set with only \eqn{N=1493} subjects.

} 


\details{The varve patterns have the following origin. During
spring, in connection to snow melt and spring runoff, minerogenic material is
transported from the catchment area into the lake through four small streams, which
gives rise to a bright colored layer, giving high gray-scale values (Petterson et al., 2010). During summer, autochthonous organic matter, sinks to the bottom and creates a darker layer (lower gray-scale values). During the winter, when the lake is
ice-covered, fine organic material is deposited, resulting in a thin blackish winter layer
(lowest gray-scale values). There is substantial within- and between year
variation, reflecting the balance between minerogenic and organic material. The
properties of each varve reflect, to a large extent, weather conditions and internal
biological processes in the lake the year that the varve was deposited. The minerogenic
input reflects the intensity of the spring run-off, which is dependent on the amount of
snow accumulated during the winter, and hence the variability in past winter climate.

The data consists of \eqn{N = 6326} (subjects) years and the \eqn{n_i} observations per year
ranges from 4 to 37. A few years were missing, see Arnqvist et al. (2016) for details. For
each year \eqn{i} we observe the centered seasonal pattern in terms of grey scale values (\eqn{y_i})'s
at \eqn{n_i} time points (pixels). We also register (the four covariates) the mean grey scale
within each year, the varve width \eqn{n_i} and the minerogenic accumulation rate (\eqn{mg/cm^2})
corresponding to the total amount of minerogenic material per \eqn{cm^2} in the varve (year) \eqn{i},
see Petterson et al. (2010) for details. In order to make the seasonal patterns comparable
we first put them on the same time scale \eqn{[0,1]}, such that pixel position \eqn{j} at year \eqn{i}
corresponds to position \eqn{\tilde{t}_{ij} = (j-1)/(n_i-1),~ j = 1, ..., n _i, ~i = 1, ..., N.} To make
the patterns more comparable (with respect to climate) they were further aligned by
landmark registration, synchronizing the first spring peaks, that are directly related to
the spring flood that occurs approximately the same time each year.

As in previous analysis (Arnqvist et al., 2016) the first peak of each year is aligned
towards a common spring peak with an affine warping, that is, if we denote the common
spring peak as \eqn{M_L} and the yearly spring peak as \eqn{L_i,~ i = 1, ..., N} and let \eqn{b = M_L/L_i},
\eqn{d=(1-M_L)/(1-L_i)} then we will have the warped time points according to \eqn{w(t_{ij}) = t_{ij}b}
for \eqn{t_{ij} < L_i} and \eqn{w(t_{ij}) = 1 + d(t_{ij}-1)} for \eqn{t_{ij}\geq L_i}. The common spring peak and the
yearly spring peaks are given in Arnqvist et al. (2016).

Focusing on the functional forms of the seasonal patterns we finally centered them
within years and worked with (the centered values) \eqn{y_i(t_{ij}) -\bar{y}_i,~ j = 1, ..., n_i,~ i = 1,...,N}
where \eqn{\bar{y}_i=\sum_{j=1}^{n_i}y_i(t_{ij})/n_i} is the mean grey scale value of varve (year) \eqn{i}. In addition
to the seasonal patterns we also include four covariates: i) \eqn{x_{1i}=\bar{y}_i}, the mean grey scale; ii) \eqn{x_{2i} = n_i}, the varve width (proportional to \eqn{n_i}); iii) \eqn{x_{3i}}, the minerogenic accumulation rate
corresponding to the accumulated amount of minerogenic material per \eqn{cm^2} in varve \eqn{i}; and iv) \eqn{x_{4i}}, the landmark which is the distance from the start of the year to the first peak, interpreted as the start of the spring, for details see (Petterson et al., 2010, and and Arnqvist et al. 2016).

\code{varve_full} is a full data set with \eqn{N=6326} years/curves spanning the time period 4486 B.C. until 1901 A.D..

\code{varve} is a reduced data set with \eqn{N=1493} years/curves covering the time period 750 BC to 750 AD.

}

\references{
Arnqvist, P., Bigler, C., Renberg, I., \enc{Sjöstedt}{Sjostedt} de Luna, S. (2016). Functional clustering of varved lake sediment to reconstruct past seasonal climate. \emph{Environmental and Ecological Statistics}, \bold{23}(4), 513--529.

Abramowicz, K., Arnqvist, P., Secchi, P., \enc{Sjöstedt}{Sjostedt} de Luna, S., Vantini, S., Vitelli, V. (2017). Clustering misaligned dependent curves applied to varved lake sediment for climate reconstruction. \emph{Stochastic Environmental Research and Risk Assessment}. Volume \bold{31.1}, 71--85.

Arnqvist, P., and \enc{Sjöstedt}{Sjostedt} de Luna, S. (2019). Model based functional clustering of varved lake sediments. \emph{arXiv preprint arXiv:1904.10265}.

Petterson, G., Renberg, I., \enc{Sjöstedt}{Sjostedt} de Luna, S., Arnqvist, P., and Anderson, N. J. (2010). Climatic influence on the inter-annual variability of late-Holocene minerogenic sediment supply in a boreal forest catchment. \emph{Earth surface processes and landforms}. \bold{35}(4), 390--398.

Petterson, G., B. Odgaard, and I. Renberg (1999). Image analysis as a method to quantify sediment components. \emph{Journal of Paleolimnology 22}. (4), 443--455.

Petterson, G., I. Renberg, P. Geladi, A. Lindberg, and F. Lindgren (1993). Spatial uniformity of sediment accumulation in varved lake sediments in Northern Sweden. \emph{Journal of Paleolimnology}. \bold{9}(3), 195--208.

}

%\source{
%The data are adapted from ...., where the original source is given as: 
%}

%\examples{
%}








