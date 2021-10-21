\name{plot.mocca}
%\Rdversion{1.0}
\alias{plot.mocca}
\title{mocca plotting
}
\description{The function
   takes a \code{mocca} object produced by \code{mocca()} and creates cluster means plots or covariance structure within each cluster.}

\usage{
\method{plot}{mocca}(x,type=1, select =NULL,transform=FALSE,covariance=TRUE,
    covariates =FALSE,lwd=2,ylab="",xlab="",main="",ylim=NULL,
    ncolors=NULL,probs=FALSE,pts=FALSE,size=50,
    years=NULL, years.names=NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a \code{mocca} object as produced by \code{mocca()}.}

   \item{type}{determines what type of plots to print. For \code{type=1} (default) cluster mean curves are shown in one plot on one page together with the overall mean curve; \code{type=2} produces the trend of the frequencies of the different clusters, together with mean probabilites (if \code{probs=TRUE}), the mean value of the included covariates (if present) within each cluster (not the model estimated covariate values) are also shown, if \code{pts=TRUE} points of the frequency trend are plotted, cluster means are shown on separate plots; \code{type=3} illustrates the covariance (or correlation) structure within each cluster. \code{type=2} is used with annual data.}
  
   \item{select}{allows the plot for a single cluster mean to be selected for printing with \code{type=1} or \code{type=2}. it can also be the order of the cluster means to be printed. If \code{NULL} (default), the cluster mean curves are in \code{{1,2,...,K}} order, where \code{K} is the number of clusters.  If you just want the plot for the cluster mean of the second cluster set \code{select=2}.}

  \item{transform}{logical, informs whether svd back-transformation of the spline model matrix should be applied (see Arnqvist and Sjöstedt de Luna, 2019).}
  \item{covariance}{logical, informs whether covariance (\code{TRUE}) or correlation (\code{FALSE}) matrices should be plotted}
  \item{covariates}{logical, informs whether covariates should be added when printing the covariance structure of the spline coefficients}
  \item{lwd}{defines the line width.}

  \item{ylab}{If supplied then this will be used as the \code{y} label for all plots.}
  \item{xlab}{If supplied then this will be used as the \code{x} label for all plots.}
  \item{main}{Used as title for plots if supplied.}
  \item{ylim}{If supplied then this pair of numbers are used as the \code{y} limits for each plot. Default \code{ylim=c(-45, 55)}.}
  \item{ncolors}{defines the number of colors (\eqn{\geq 1}) to be in the palette, used with the \code{rainbow()} function. 
   If \code{NULL} (default), \code{ncolors} equals the number of clusters \code{K}.}
   \item{probs}{logical, used with \code{type=2}, informs whether the mean probabilites should be printed.}
   \item{pts}{logical, used with \code{type=2}, if \code{TRUE} (default) points of the frequency trend are shown.}
   \item{size}{the bin size used with \code{type=2} (default: 50 years), the bin size of how many of those years belong to a specific cluster.}
  % \item{select}{}
  % \item{stand.cov}{logical, used with \code{type=2}, tells whether standardized covariates were used in the modelling.}
   \item{years}{a vector of years used with annual data and needed for \code{type=2} plot to calculate frequencies in the bins of size provided by the \code{size} argument.}
   \item{years.names}{a character vector that gives names of the years needed for \code{type=2} plot. This can be also supplied with \code{data}. With \code{varve} data \code{years.names} are supplied as rownames of the matrix of covariates. if \code{years.names=NULL} (default) then \code{years} are converted to the character vector and used as \code{years.names}.}
  \item{...}{ other graphics parameters to pass on to plotting commands.}

}

% \details{ to be added...
%%  ~~ If necessary, more details than the description above ~~
% }
\value{The function generates plots.

}

\author{
Per Arnqvist, Sara \enc{Sjöstedt}{Sjostedt} de Luna,  Natalya Pya Arnqvist
}
\references{
 Arnqvist, P., and \enc{Sjöstedt}{Sjostedt} de Luna, S. (2019). Model based functional clustering of varved lake sediments. \emph{arXiv preprint arXiv:1904.10265}.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{  \code{\link{mocca}}}

\examples{
## see ?mocca help files
}


