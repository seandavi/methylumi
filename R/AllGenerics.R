
## Generic methods
if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object) standardGeneric("getHistory"))
if (is.null(getGeneric("betas"))) setGeneric("betas", function(object) standardGeneric("betas"))
if (is.null(getGeneric("betas<-"))) setGeneric("betas<-", function(object, value) standardGeneric("betas<-"))
if (is.null(getGeneric("mvals"))) setGeneric("mvals", function(object) standardGeneric("mvals"))
if (is.null(getGeneric("mvals<-"))) setGeneric("mvals<-", function(object, value) standardGeneric("mvals<-"))
if (is.null(getGeneric("pvals"))) setGeneric("pvals", function(object) standardGeneric("pvals"))
if (is.null(getGeneric("pvals<-"))) setGeneric("pvals<-", function(object, value) standardGeneric("pvals<-"))
if (is.null(getGeneric("unmethylated"))) setGeneric("unmethylated", function(object) standardGeneric("unmethylated"))
if (is.null(getGeneric("unmethylated<-"))) setGeneric("unmethylated<-", function(object, value) standardGeneric("unmethylated<-"))
if (is.null(getGeneric("unmethylated.OOB"))) setGeneric("unmethylated.OOB", function(object) standardGeneric("unmethylated.OOB"))
if (is.null(getGeneric("unmethylated.OOB<-"))) setGeneric("unmethylated.OOB<-", function(object,value) standardGeneric("unmethylated.OOB<-"))
if (is.null(getGeneric("methylated"))) setGeneric("methylated", function(object) standardGeneric("methylated"))
if (is.null(getGeneric("methylated<-"))) setGeneric("methylated<-", function(object, value) standardGeneric("methylated<-"))
if (is.null(getGeneric("methylated.OOB"))) setGeneric("methylated.OOB", function(object) standardGeneric("methylated.OOB"))
if (is.null(getGeneric("methylated.OOB<-"))) setGeneric("methylated.OOB<-", function(object,value) standardGeneric("methylated.OOB<-"))
if (is.null(getGeneric("hist"))) setGeneric("hist", function(x,...) standardGeneric("hist"))
if (is.null(getGeneric("corplot"))) setGeneric("corplot", function(x,...) standardGeneric("corplot"))
if (is.null(getGeneric("plotDensity"))) setGeneric("plotDensity", function(x,...) standardGeneric("plotDensity"))
if (is.null(getGeneric("QCdata"))) setGeneric("QCdata", function(object) standardGeneric("QCdata"))
if (is.null(getGeneric("QCdata<-"))) setGeneric("QCdata<-", function(object, value) standardGeneric("QCdata<-"))

# generic methods specific for MethyLumiM class
if (is.null(getGeneric("methylated.N"))) setGeneric("methylated.N", function(object) standardGeneric("methylated.N"))
if (is.null(getGeneric("methylated.N<-"))) setGeneric("methylated.N<-", function(object, value) standardGeneric("methylated.N<-"))
if (is.null(getGeneric("unmethylated.N"))) setGeneric("unmethylated.N", function(object) standardGeneric("unmethylated.N"))
if (is.null(getGeneric("unmethylated.N<-"))) setGeneric("unmethylated.N<-", function(object, value) standardGeneric("unmethylated.N<-"))
if (is.null(getGeneric("controlData"))) setGeneric("controlData", function(object) standardGeneric("controlData"))
if (is.null(getGeneric("controlData<-"))) setGeneric("controlData<-", function(object, value) standardGeneric("controlData<-"))
if (is.null(getGeneric("detection"))) setGeneric("detection", function(object) standardGeneric("detection"))
if (is.null(getGeneric("detection<-"))) setGeneric("detection<-", function(object, value) standardGeneric("detection<-"))
if (is.null(getGeneric("dataType"))) setGeneric("dataType", function(object) standardGeneric("dataType"))
if (is.null(getGeneric("dataType<-"))) setGeneric("dataType<-", function(object, value) standardGeneric("dataType<-"))


if (is.null(getGeneric("featureFilter"))) # {{{
    setGeneric("featureFilter", 
               function(eset,
                        require.entrez=FALSE,
                        require.GOBP=FALSE,
                        require.GOCC=FALSE,
                        require.GOMF=FALSE,
                        exclude.ChrX=FALSE,
                        require.closeToTSS=FALSE,
                        range.DistToTSS=c(-500, 300),
                        require.CpGisland=FALSE, ...)
               standardGeneric("featureFilter")) # }}}

if (is.null(getGeneric("varFilter"))) # {{{
    setGeneric("varFilter", 
               function(eset,
                        var.func=IQR, var.cutoff=0.5,
                        filterByQuantile=TRUE, ...)
               standardGeneric("varFilter")) # }}}

setGeneric('total.intensity', # {{{
             function(object) standardGeneric('total.intensity')) # }}}

setGeneric('unmethylated.N', # {{{
              function(object) standardGeneric('unmethylated.N')) # }}}

setGeneric('unmethylated.N<-', # {{{
              function(object,value) standardGeneric('unmethylated.N<-')) # }}}
setGeneric('methylated.N', # {{{
              function(object) standardGeneric('methylated.N')) # }}}
setGeneric('methylated.N<-', # {{{
              function(object, value) standardGeneric('methylated.N<-')) # }}}

## FIXME: use COLOR_CHANNEL for this, and add 'both'/'d2', or else retire it
setGeneric('getProbesByChannel', # {{{
  function(object, ...) standardGeneric('getProbesByChannel')) # }}}

setGeneric('intensitiesByChannel', # {{{
  function(object, ...) standardGeneric('intensitiesByChannel')) # }}}

setGeneric('Cy3', # {{{
             function(object) standardGeneric('Cy3')) # }}}
setGeneric("Cy3<-", # {{{
              function(object, value) standardGeneric('Cy3<-')) # }}}

setGeneric('Cy3.SD', # {{{
           function(object) standardGeneric('Cy3.SD')) # }}}
setGeneric('Cy3.N', # {{{
           function(object) standardGeneric('Cy3.N')) # }}}

setGeneric('Cy5', # {{{
             function(object) standardGeneric('Cy5')) # }}}
setGeneric("Cy5<-", # {{{
              function(object, value) standardGeneric('Cy5<-')) # }}}


setGeneric('Cy5.SD', # {{{
          function(object) standardGeneric('Cy5.SD')) # }}}
setGeneric('Cy5.N', # {{{
           function(object) standardGeneric('Cy5.N')) # }}}

setGeneric('negctls', # {{{
           function(object, channel) standardGeneric('negctls')) # }}}

setGeneric('negctls.SD', # {{{
  function(object, channel) standardGeneric('negctls.SD')) # }}}

setGeneric('negctls.stderr', # {{{
  function(object, channel) standardGeneric('negctls.stderr')) # }}}

setGeneric('negnorm', # {{{
           function(object, channel) standardGeneric('negnorm')) # }}}

setGeneric('normctls', # {{{
  function(object, ...) standardGeneric('normctls')) # }}}
