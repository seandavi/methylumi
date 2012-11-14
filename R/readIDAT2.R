readIDAT2 <- function(idatFile){
  fileSize <- file.info(idatFile)$size

  tempCon <- file(idatFile,"rb")
  prefixCheck <- readChar(tempCon,4)
  if(prefixCheck != "IDAT"){

  }

  versionNumber <- readBin(tempCon, "integer", n=1, size=8, endian="little")

  if(versionNumber<3)
	  stop("Older style IDAT files not supported:  consider updating your scanner settings")

  nFields <- readBin(tempCon, "integer", n=1, size=4, endian="little")

  fields <- matrix(0,nFields,3);
  colnames(fields) <- c("Field Code", "Byte Offset", "Bytes")
  for(i1 in 1:nFields){
    fields[i1,"Field Code"] <-
      readBin(tempCon, "integer", n=1, size=2, endian="little", signed=FALSE)
    fields[i1,"Byte Offset"] <-
      readBin(tempCon, "integer", n=1, size=8, endian="little")
  }

  knownCodes <-
    c("nSNPsRead"  = 1000,
      "IlluminaID" =  102,
      "SD"         =  103,
      "Mean"       =  104,
      "NBeads"     =  107,
      "MidBlock"   =  200,
      "RunInfo"    =  300,
      "RedGreen"   =  400,
      "MostlyNull" =  401,
      "Barcode"    =  402,
      "ChipType"   =  403,
      "MostlyA"    =  404,
      "Unknown.1"  =  405,
      "Unknown.2"  =  406,
      "Unknown.3"  =  407,
      "Unknown.4"  =  408,
      "Unknown.5"  =  409,
      "Unknown.6"  =  410,
      "Unknown.7"  =  510
      )

  nNewFields <- 1
  rownames(fields) <- paste("Null", 1:nFields)
  for(i1 in 1:nFields){
    temp <- match(fields[i1,"Field Code"], knownCodes)
    if(!is.na(temp)){
      rownames(fields)[i1] <- names(knownCodes)[temp]
    }else{
      rownames(fields)[i1] <- paste("newField", nNewFields, sep=".")
      nNewFields <- nNewFields + 1
    }
  }

  fields <- fields[order(fields[, "Byte Offset"]),]

  seek(tempCon, fields["nSNPsRead", "Byte Offset"])
  nSNPsRead <- readBin(tempCon, "integer", n=1, size=4, endian="little")

  readBlock <- function(nam) {
      switch(nam,
             "IlluminaID" = {
                 seek(tempCon, fields["IlluminaID", "Byte Offset"])
                 IlluminaID <- readBin(tempCon, "integer", n=nSNPsRead, size=4, endian="little")
                 IlluminaID
             },
             "SD" = {
                 seek(tempCon, fields["SD", "Byte Offset"])
                 SD <- readBin(tempCon, "integer", n=nSNPsRead, size=2, endian="little", signed=FALSE)
                 SD
             },
             "Mean" = {
                 seek(tempCon, fields["Mean", "Byte Offset"])
                 Mean <- readBin(tempCon, "integer", n=nSNPsRead, size=2, endian="little", signed=FALSE)
                 Mean
             },
             "NBeads" = {
                 seek(tempCon, fields["NBeads", "Byte Offset"])
                 NBeads <- readBin(tempCon, "integer", n=nSNPsRead, size=1, signed=FALSE)
                 NBeads
             },
             "MidBlock" = {
                 seek(tempCon, fields["MidBlock", "Byte Offset"])
                 nMidBlockEntries <- readBin(tempCon, "integer", n=1, size=4, endian="little")
                 MidBlock <- readBin(tempCon, "integer", n=nMidBlockEntries, size=4,
                                     endian="little")
                 MidBlock
             },
             "RunInfo" = {
                 seek(tempCon, fields["RunInfo", "Byte Offset"])
                 nRunInfoBlocks <- readBin(tempCon, "integer", n=1, size=4, endian="little")
                 RunInfo <- matrix(NA, nRunInfoBlocks, 5)
                 colnames(RunInfo) <- c("RunTime", "BlockType", "BlockPars",
                                        "BlockCode", "CodeVersion")
                 for(i1 in 1:2) { #nRunInfoBlocks){  ## MR edit
                     for(i2 in 1:5){
                         nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                         RunInfo[i1,i2] <- readChar(tempCon, nChars)
                     }
                 }
                 RunInfo
             },
             "RedGreen" = {
                 seek(tempCon, fields["RedGreen", "Byte Offset"])
                 RedGreen <- readBin(tempCon, "numeric", n=1, size=4,
                                     endian="little")
                 RedGreen
             },
             "MostlyNull" = {
                 seek(tempCon, fields["MostlyNull", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 MostlyNull <- readChar(tempCon, nChars)
                 MostlyNull
             },
             "Barcode" = {                 
                 seek(tempCon, fields["Barcode", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Barcode <- readChar(tempCon, nChars)
                 Barcode
             },
             "ChipType" = {
                 seek(tempCon, fields["ChipType", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 ChipType <- readChar(tempCon, nChars)
                 ChipType
             },
             "MostlyA" = {
                 seek(tempCon, fields["MostlyA", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 MostlyA <- readChar(tempCon, nChars)
             },
             "Unknown.1" = {
                 seek(tempCon, fields["Unknown.1", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Unknown.1 <- readChar(tempCon, nChars)
                 Unknown.1
             },
             "Unknown.2" = {
                 seek(tempCon, fields["Unknown.2", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Unknown.2 <- readChar(tempCon, nChars)
                 Unknown.2
             },
             "Unknown.3" = {
                 seek(tempCon, fields["Unknown.3", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Unknown.3 <- readChar(tempCon, nChars)
                 Unknown.3
             },
             "Unknown.4" = {
                 seek(tempCon, fields["Unknown.4", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Unknown.4 <- readChar(tempCon, nChars)
                 Unknown.4
             },
             "Unknown.5" = {
                 seek(tempCon, fields["Unknown.5", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Unknown.5 <- readChar(tempCon, nChars)
                 Unknown.5
             },
             "Unknown.6" = {
                 seek(tempCon, fields["Unknown.6", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Unknown.6 <- readChar(tempCon, nChars)
                 Unknown.6
             },
             "Unknown.7" = {
                 seek(tempCon, fields["Unknown.7", "Byte Offset"])
                 nChars <- readBin(tempCon, "integer", n=1, size=1, signed=FALSE)
                 Unknown.7 <- readChar(tempCon, nChars)
                 Unknown.7
             })
  }

  readFields <- setdiff(rownames(fields), "nSNPsRead")
  names(readFields) <- readFields
  
  allFields <- lapply(readFields, readBlock)

  close(tempCon)

  UnknownNames <- c("MostlyNull", "MostlyA", "Unknown.1",
                    "Unknown.2", "Unknown.3", "Unknown.4",
                    "Unknown.5", "Unknown.6", "Unknown.7")
  Unknowns <- allFields[intersect(names(allFields), UnknownNames)]

  Quants <- cbind(allFields$Mean, allFields$SD, allFields$NBeads)
  colnames(Quants) <- c("Mean", "SD", "NBeads")
  rownames(Quants) <- as.character(allFields$IlluminaID)

  InfoNames <- c("MidBlock", "RunInfo", "RedGreen", "Barcode", "ChipType")
  Info <- allFields[intersect(names(allFields), InfoNames)]
  
  idatValues <-
    list(fileSize=fileSize,
         versionNumber=versionNumber,
         nFields=nFields,
         fields=fields,
         nSNPsRead=nSNPsRead,
         Quants=Quants)
  idatValues <- c(idatValues, Info, list(Unknowns = Unknowns))
  idatValues
}
