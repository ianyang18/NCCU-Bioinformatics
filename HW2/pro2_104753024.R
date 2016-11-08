#################################################################
# Goal: Caculating the sequence alignment
# Input-1 : Testing file, see the format in "test.fasta" file
# Input-2 : Scoring table
# Input-3 : Strategies, global/local alignment
# Input-4 : gap_open
# Output  : Stroing the sequence alignment result
##################################################################
# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript pro2_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- tolower(args[i+1])
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<-strtoi(args[i+1])
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-strtoi(args[i+1])
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

#print("PARAMETERS")
#print(paste("input file         :", i_f))
#print(paste("output file        :", o_f))
#print(paste("score file         :", s_f))
#print(paste("gap open penalty   :", g_o))
#print(paste("gap extend penalty :", g_e))
#print(paste("alignment          :", aln_mode))

######################################
# main
######################################

########################################################
#
# Initialize the substitution/deletion/insertion matrice, which is used to record the optimal score
# Parameters: sequence_A, sequence_B, alignment_strategy, gap scoring scheme
# Output: list of matrices
#
########################################################
initializeScore <- function(sa, sb, method, gapOE) {
    la <- nchar(sa)
    mIndex <- la + 1
    lb <- nchar(sb)
    nIndex <- lb + 1
    # Initialize the table to record the score, first sequence(sa) in row, second sequence(sb) in column
    # table index starts from 1
    if(method == "global") {
        if(gapOE$open & gapOE$extend){
            dMatrix <- matrix(rep(0, mIndex*nIndex), nrow=mIndex, ncol=nIndex, byrow=TRUE)
            if(la != 0) {
                for(i in (1:mIndex)) {
                    dMatrix[i, 1] <- (i-1) * gap + gapOE$open
                }
            }
            if(lb != 0) {
                dMatrix[1, 2:nIndex] <- (-Inf)
            }

            iMatrix <<- matrix(rep(0, mIndex*nIndex), nrow=mIndex, ncol=nIndex, byrow=TRUE)
            if(la != 0) {
                iMatrix[2:mIndex, 1] <- (-Inf)
            }
            if(lb != 0) {
                for(j in (1:nIndex)) {
                    iMatrix[1, j] <- (j-1) * gap + gapOE$open
                }
            }

            sMatrix <- matrix(rep(0, mIndex*nIndex), nrow=mIndex, ncol=nIndex, byrow=TRUE)
            if(la != 0) {
                for(i in (1:mIndex)) {
                    sMatrix[i, 1] <- dMatrix[i, 1]
                }
            }
            if(lb != 0) {
                for(j in (1:nIndex)) {
                    sMatrix[1, j] <- iMatrix[1, j]
                }
            }
            sMatrix[1, 1] <- initialValue
        }
        else if(!gapOE$open & !gapOE$extend) {
            sMatrix <- matrix(rep(0, mIndex*nIndex), nrow=mIndex, ncol=nIndex, byrow=TRUE)
            if(la != 0) {
                for(i in (1:mIndex)) {
                    sMatrix[i, 1] <- (i-1) * gap + initialValue
                }
            }
            if(lb != 0) {
                for(j in (1:nIndex)) {
                    sMatrix[1, j] <- (j-1) * gap + initialValue
                }
            }
            #sMatrix[1, 1] <- initialValue
        }
    }
    else if(method == "local") {
        dMatrix <- matrix(rep(0, mIndex*nIndex), nrow=mIndex, ncol=nIndex, byrow=TRUE)
        if(la != 0) {
            for(i in (1:mIndex)) {
                dMatrix[i, 1] <- (i-1) * gap + gapOE$open
            }
        }
        if(lb != 0) {
            dMatrix[1, 2:nIndex] <- (-Inf)
        }

        iMatrix <<- matrix(rep(0, mIndex*nIndex), nrow=mIndex, ncol=nIndex, byrow=TRUE)
        if(la != 0) {
            iMatrix[2:mIndex, 1] <- (-Inf)
        }
        if(lb != 0) {
            for(j in (1:nIndex)) {
                iMatrix[1, j] <- (j-1) * gap + gapOE$open
            }
        }

        sMatrix <- matrix(rep(initialValue, mIndex*nIndex), nrow=mIndex, ncol=nIndex, byrow=TRUE)
    }
    if((method == "global") & !gapOE$open & !gapOE$extend) {
        list("sMatrix"=sMatrix)
    }
    else {
        list("sMatrix"=sMatrix, "dMatrix"=dMatrix, "iMatrix"=iMatrix)
    }
}

##########################################################
#
# Calculating the sequence alignment score according to the substitution/deletion/insertion table
# Parameters: sequence_A, sequence_B, alignment_strategy, gap scoring scheme, score table
# Output: Score
#
##########################################################
calScore <- function(sa, sb, method, gapOE, matrices) {
    la <- nchar(sa)
    mIndex <- la + 1
    lb <- nchar(sb)
    nIndex <- lb + 1
    if(method == "local") {
        bestScore <- initialValue
        iEndPos <- 0
        jEndPos <- 0
    }
    if(la != 0 & lb != 0) {
        for(i in (2:mIndex)) {
            for(j in (2:nIndex)) {
                charA <- substr(sa, i-1, i-1)
                charB <- substr(sb, j-1, j-1)
                # gap_open, gap_extend scheme
                if(gapOE$open & gapOE$extend) {
                    dExtendScore <- matrices$dMatrix[i-1, j] + gap
                    dOpenScore <- matrices$sMatrix[i-1, j] + gapOE$open + gap
                    matrices$dMatrix[i, j] = max(dExtendScore, dOpenScore)

                    iExtendScore <- matrices$iMatrix[i, j-1] + gap
                    iOpenScore <- matrices$sMatrix[i, j-1] + gapOE$open + gap
                    matrices$iMatrix[i, j] = max(iExtendScore, iOpenScore)
                    
                    if(exists("s_f")) {
                        sScore <- matrices$sMatrix[i-1, j-1] + scheme[charA, charB]
                    }
                    else {
                        sScore <- matrices$sMatrix[i-1, j-1] + if(charA == charB) match else miss
                    }

                    if(method == "global") {
                        matrices$sMatrix[i, j] <- max(sScore, matrices$iMatrix[i, j], matrices$dMatrix[i, j])
                    }
                    else if(method == "local") {
                        matrices$sMatrix[i, j] <- max(initialValue, sScore, matrices$iMatrix[i, j], matrices$dMatrix[i, j])
                        if(matrices$sMatrix[i, j] > bestScore) {
                            bestScore <- matrices$sMatrix[i, j]
                            iEndPos <- i
                            jEndPos <- j
                        }
                    }
                }
                # normal scheme
                else if(!gapOE$open & !gapOE$extend) {
                    if(exists("s_f")) {
                        sScore <- matrices$sMatrix[i-1, j-1] + scheme[charA, charB]
                    }
                    else {
                        sScore <- matrices$sMatrix[i-1, j-1] + if(charA == charB) match else miss
                    }
                    iScore <- matrices$sMatrix[i, j-1] + gap
                    dScore <- matrices$sMatrix[i-1, j] + gap
                    if(method == "global") {
                        matrices$sMatrix[i, j] <- max(sScore, iScore, dScore)
                    }
                    else if(method == "local") {
                        matrices$sMatrix[i, j] <- max(initialValue, sScore, iScore, dScore)
                        if(matrices$sMatrix[i, j] > bestScore) {
                            bestScore <- matrices$sMatrix[i, j]
                            iEndPos <- i
                            jEndPos <- j
                        }
                    }
                }
            }
        }
    }
    if(method == "global") {
        list("score"=matrices$sMatrix[mIndex, nIndex], "matrix"=matrices, "position"=c(mIndex, nIndex))
    }
    else if(method == "local") {
        list("score"=bestScore, "matrix"=matrices, "position"=c(iEndPos, jEndPos))
    }
}

####################################################################
#
# Given the index and the score to backtrace.
# Input:  Sequence A, Sequence B, score table, high score's index position, strategy
# Output: Alignment sequences
#
####################################################################
output <- function(sa, sb, matrices, mIndex, nIndex, method) {
    sMatrix <- matrices$sMatrix
    # i, j are the index for the string
    i <- mIndex - 1
    j <- nIndex -1
    charA <- substr(sa, i, i)
    charB <- substr(sb, j, j)
    #i, j plus one because of matrix index from 1
    if(method == "global" & (mIndex == 1 | nIndex == 1)) {
        if(mIndex > 1) {
            stringA <<- substr(sa, 1, i)
            stringB <<- paste(stringB, rep("-", i), collapse="", sep="")
        }
        if(nIndex > 1) {
            stringB <<- substr(sb, 1, j)
            stringA <<- paste(stringA, rep("-", j), collapse="", sep="")
        }
    }
    else if(method == "local" & sMatrix[mIndex, nIndex] == initialValue) {
    }
    else {
        if(exists("s_f")) {
            if((sMatrix[mIndex, nIndex] == sMatrix[mIndex-1, nIndex-1] + scheme[charA, charB])) {
                output(sa, sb, matrices, mIndex-1, nIndex-1, method)
                stringA <<- paste(stringA, substr(sa, i, i), sep="")
                stringB <<- paste(stringB, substr(sb, j, j), sep="")
            }
            else if(sMatrix[mIndex, nIndex] == sMatrix[mIndex-1, nIndex] + gap) {
                output(sa, sb, matrices, mIndex-1, nIndex, method)
                stringA <<- paste(stringA, substr(sa, i, i), sep="")
                stringB <<- paste(stringB, "-", sep="")
            } 
            else {
                output(sa, sb, matrices, mIndex, nIndex-1, method)
                stringA <<- paste(stringA, "-", sep="")
                stringB <<- paste(stringB, substr(sb, j, j), sep="")
            }
        }           
        else {
            if ((sMatrix[mIndex, nIndex] == sMatrix[mIndex-1, nIndex-1] + if (charA == charB) match else miss)) {
                output(sa, sb, matrices, mIndex-1, nIndex-1, method)
                stringA <<- paste(stringA, substr(sa, i, i), sep="")
                stringB <<- paste(stringB, substr(sb, j, j), sep="")
            } 
            else if(sMatrix[mIndex, nIndex] == sMatrix[mIndex-1, nIndex] + gap) {
                output(sa, sb, matrices, mIndex-1, nIndex, method)
                stringA <<- paste(stringA, substr(sa, i, i), sep="")
                stringB <<- paste(stringB, "-", sep="")
            } 
            else {
                output(sa, sb, matrices, mIndex, nIndex-1, method)
                stringA <<- paste(stringA, "-", sep="")
                stringB <<- paste(stringB, substr(sb, j, j), sep="")
            }
        }
    }
}


data <- read.table(i_f)
name <- names(data)
if(nrow(data) != 4 & !grepl("^[>]", data[[name]][1]) & !grepl("^[>]", data[[name]][3]) & !grepl("^[a-zA-Z]*$", data[[name]][2]) & !grepl("^[a-zA-Z]*$", data[[name]][4])) {
    print("Check your input file's format!")
} else {
    seqA <- toupper(lapply(data[[name]][2], as.character))
    seqB <- toupper(lapply(data[[name]][4], as.character))
    
    # load the score table
    if(exists("s_f")) {
        scheme <- read.table(s_f)
        gap <- scheme[nrow(scheme), 1]    # in PAM205.txt, assuming char([A-Z]) alignment with char(*) can be seen as the gap penalty
    }
    else {
        match <- 8    # same character
        miss <- -5    # different character
        gap <- -3     # gap penalty
    }

    gapOE <- list("open"=0, "extend"=0)
    if(exists("g_o") & exists("g_e")) {
        gapOE <- list("open"=g_o, "extend"=g_e)
        if(gapOE$open & gapOE$extend) {
            gap <- gapOE$extend
        }
    }
    
    # initial value of the table
    #initialValue <- scheme[nrow(scheme), ncol(scheme)]
    initialValue <- 0

    stringA <<- ""
    stringB <<- ""

    matrices <- initializeScore(seqA, seqB, aln_mode, gapOE)
    score <- calScore(seqA, seqB, aln_mode, gapOE, matrices)

    if(aln_mode == "global") {
        output(seqA, seqB, score$matrix, score$position[1], score$position[2], aln_mode)
    } else if(aln_mode == "local") {
        output(seqA, seqB, score$matrix, score$position[1], score$position[2], aln_mode)
    }

    # for debug use
    #print(score$matrix["sMatrix"])
    #print(score$score)
    #print(stringA)
    #print(stringB)
    write.table(c(lapply(data[[name]][1], as.character), "\n", stringA, "\n", lapply(data[[name]][3], as.character), "\n", stringB), o_f, quote=F, col.names=FALSE, row.names=FALSE)
}
