setClass('MethyLumi',
         contains='eSet')

setClass('MethyLumiQC',
         contains='MethyLumi')

setClassUnion("QCDataOrNULL",c('NULL',"MethyLumiQC"))

## Generic methods
if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object) standardGeneric("getHistory"))
if (is.null(getGeneric("betas"))) setGeneric("betas", function(object) standardGeneric("betas"))
if (is.null(getGeneric("betas<-"))) setGeneric("betas<-", function(object, value) standardGeneric("betas<-"))
if (is.null(getGeneric("pvals"))) setGeneric("pvals", function(object) standardGeneric("pvals"))
if (is.null(getGeneric("pvals<-"))) setGeneric("pvals<-", function(object, value) standardGeneric("pvals<-"))
if (is.null(getGeneric("exprs"))) setGeneric("exprs", function(object) standardGeneric("exprs"))
if (is.null(getGeneric("exprs<-"))) setGeneric("exprs<-", function(object, value) standardGeneric("exprs<-"))
if (is.null(getGeneric("unmethylated"))) setGeneric("unmethylated", function(object) standardGeneric("unmethylated"))
if (is.null(getGeneric("unmethylated<-"))) setGeneric("unmethylated<-", function(object, value) standardGeneric("unmethylated<-"))
if (is.null(getGeneric("methylated"))) setGeneric("methylated", function(object) standardGeneric("methylated"))
if (is.null(getGeneric("methylated<-"))) setGeneric("methylated<-", function(object, value) standardGeneric("methylated<-"))
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


