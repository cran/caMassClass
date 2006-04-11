\name{caMassClass-package}
\alias{caMassClass-package}
\alias{caMassClass}
\docType{package}
\title{Processing and Classification of Protein Mass Spectra Data}
\description{
Functions for processing and classification of protein
mass spectra data. Includes support for I/O in mzXML and CSV formats.
}
\details{
\tabular{ll}{
Package: \tab caMassClass \cr
Version: \tab 1.5 \cr
Date: \tab 2006-04-11 \cr
Depends: \tab R (>= 2.0.0), PROcess, e1071, nnet, rpart, caTools, XML, digest \cr
License: \tab The caMassClass Software License, Version 1.0 (See COPYING file \cr
 \tab or \url{http://ncicb.nci.nih.gov/download/camassclasslicense.jsp})\cr
URL: \tab \url{http://ncicb.nci.nih.gov/download/index.jsp} \cr
}

Index of Preprocessing Functions:
\tabular{ll}{
 \code{\link{msc.project.read}} \tab Read and Manage a Batch of Protein Mass Spectra \cr
 \code{\link{msc.project.run}} \tab Read and Preprocess Protein Mass Spectra \cr
 \code{\link{msc.preprocess.run}} \tab Preprocessing Pipeline of Protein Mass Spectra \cr
 \code{\link{msc.baseline.subtract}} \tab Baseline Subtraction for Mass Spectra Data \cr
 \code{\link{msc.mass.cut}} \tab Remove Low Mass Portion of the Mass Spectra Data. \cr
 \code{\link{msc.mass.adjust}} \tab Perform Normalization and Mass Drift Adjustment \cr
 \code{\link{msc.peaks.find}} \tab Find Peaks of Mass Spectra \cr
 \code{\link{msc.peaks.clust}} \tab Clusters Peaks of Mass Spectra (low level)\cr
 \code{\link{msc.peaks.align}} \tab Align Peaks of Mass Spectra into Biomarker Matrix \cr
 \code{\link{msc.biomarkers.fill}} \tab Fill Empty Spaces in Biomarker Matrix \cr
 \code{\link{msc.copies.merge}} \tab Merge Multiple Copies of Mass Spectra Samples \cr
 \code{\link{msc.sample.correlation}} \tab Sample Correlation (low level) \cr
}

Index of Classification Functions:
\tabular{ll}{
 \code{\link{msc.classifier.test}} \tab Test a Classifier through Cross-validation \cr
 \code{\link{msc.features.remove}} \tab Remove Highly Correlated Features \cr
 \code{\link{msc.features.scale}} \tab Scale Classification Data \cr
 \code{\link{msc.features.select}} \tab Reduce Number of Features Prior to Classification \cr
 \code{\link{msc.classifier.run}} \tab Train and Test Chosen Classifier. \cr
}

Index of I/O Functions:
\tabular{ll}{
 \code{\link{msc.rawMS.read.csv}} \tab Read Protein Mass Spectra from CSV files \cr
 \code{\link{msc.rawMS.read.mzXML}} \tab Read raw Mass Spectra from mzXML Files \cr
 \code{\link{msc.rawMS.write.mzXML}} \tab Write raw Mass Spectra to mzXML Files \cr
 \code{\link{msc.peaks.read.csv}} \tab Read Mass Spectra Peaks in CSV Format \cr
 \code{\link{msc.peaks.write.csv}} \tab Write Mass Spectra Peaks in CSV Format \cr
 \code{\link{msc.peaks.read.mzXML}} \tab Read peaks heights and positions from mzXML File \cr
 \code{\link{msc.peaks.write.mzXML}} \tab Write peaks heights and positions to mzXML File \cr
 \code{\link{msc.biomarkers.read.csv}} \tab Read biomarker matrix in CSV format \cr
 \code{\link{msc.biomarkers.write.csv}} \tab Write biomarker matrix in CSV format \cr
 \code{\link{msc.project.read}} \tab Read and Manage a Batch of Protein Mass Spectra Data \cr
 \code{\link{read.mzXML}} \tab Read mzXML Files into generic data structure \cr
 \code{\link{write.mzXML}} \tab Write mzXML Files from generic data structure \cr
}
}
\author{Jarek Tuszynski <jaroslaw.w.tuszynski@saic.com>}

\keyword{ package }
\seealso{
 See additional documentation in "R/library/caMassClass/doc" directory.
 
 Se also Packages: PROcess 
 (\url{http://www.bioconductor.org/packages/bioc/1.8/html/PROcess.html}), xcms
 (\url{http://www.bioconductor.org/packages/bioc/1.8/html/xcms.html}), 
 ppc (CRAN).
}

