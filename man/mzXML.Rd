\name{read.mzXML}
\alias{read.mzXML}
\alias{write.mzXML}
\alias{new.mzXML}
\title{Read and Write mzXML Files}
\description{
 Read and write mzXML files.
 }
\usage{
  mzXML = new.mzXML()
  mzXML = read.mzXML(fileName)
  write.mzXML(mzXML, fileName, precision=c(32, 64)) 
}

\arguments{
  \item{mzXML}{class storing partially parsed mzXML data}
  \item{fileName}{character string with name of the file (connection)}
  \item{precision}{precision to be used in saving \code{scan} data. Save double 
    (floating point) array using 32 or 64 bits?}
}

\details{
  The main task of \code{read.mzXML} and \code{write.mzXML} functions is to 
  extract and save scan data of mzXML files. In addition attempt is made to keep
  all other sections of mzXML file as unparsed XML code, so the data can be 
  extracted latter or saved into new mzXML files. Those unparsed sections are 
  stored as XML text 
}

\value{
  Function \code{read.mzXML} returns object of type 
  \code{mzXML}, containing: 
  \item{scan}{List of Mass Spectra scans. Each element of the list contain the 
    following elements: 
    \itemize{
      \item \code{peaks} - intensities or peaks of the scan 
      \item \code{mass}  - masses (m/z) coresponding to \code{peaks}. Vectors 
      \code{mass} and \code{peaks} have the same length.
      \item \code{num} - scan number
      \item \code{parentNum} - scan number of parent scan in case of recursively
      stored scans (\code{msLevel>1})
      \item \code{msLevel} - 1- means MS scan, 2- means MS/MS scan,
      etc.
      \item \code{header} - xml code of <scan> header might contain
      other useful attributes
      \item \code{maldi} - optional - acquisition dependent properties of a 
      MALDI experiment
      \item \code{scanOrigin} - optional - name of parent file(?)
      \item \code{precursorMz} - optional - information about the precursor ion
      \item \code{nameValue} - optional - properties of the scan not included 
      elsewhere 
    }
    All optional elements contain unparsed XML code, if corresponding 
    sections are present, or \code{NULL}. See mzXML schema and documentation for 
    more details
   }
  \item{header}{Stores header of <mzXML> section containing information about
    namespace and schema file location.}
  \item{msInstrument}{General information about the MS instrument. Stored as XML.}
  \item{parentFile}{Path to all the ancestor files. Stored as XML.}
  \item{dataProcessing}{Description of any data manipulation. Stored as XML.}
  \item{separation}{Information about the separation technique. Stored as XML.}
  \item{spotting}{Acquisition independent properties of a MALDI experiment. Stored as XML.}
  \item{indexOffset}{Offset of the index element. Either 0 or a vector.}
  
  Function \code{new.mzXML} returns the same object as \code{read.mzXML} but 
  with all fields equal to \code{NULL}.
  Function \code{write.mzXML} does not return anything. 
}

\references{ Definition of \code{mzXML} format:
  \url{http://tools.proteomecenter.org/mzXMLschema.php} }
\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{\code{\link[XML]{xmlTreeParse}} from \pkg{XML} }
\examples{
  directory = system.file("Test", package = "caMassClass")
  FileName = file.path(directory,"test1.xml")
  xml  = read.mzXML(FileName)
  xml
  
  # test reading/writing
  write.mzXML(xml, "temp.xml")
  xml2 = read.mzXML("temp.xml")
  file.remove("temp.xml")
  stopifnot(all(xml$scan[[1]]$peaks == xml2$scan[[1]]$peaks))
  stopifnot(xml$msInstrument == xml2$msInstrument)
  
  # extracting scan data from the output
  FileName = file.path(directory,"test2.xml")
  xml  = read.mzXML(FileName)
  plot(xml$scan[[1]]$mass, xml$scan[[1]]$peaks, type="l")
  
  # extracting data from unparsed sections
  tree = xmlTreeParse(xml$msInstrument, asText=TRUE, asTree=TRUE)
  x = xmlRoot(tree)
  xmlName(x)
  xmlAttrs(x[["msManufacturer"]]) ["value"]
  xmlAttrs(x[["software"]])
}
\keyword{file}

