score <- function(sa, sb) {
    gap <- -3
    match <- 8
    miss <- -5
    la <- nchar(sa)
    lb <- nchar(sb)
    # Initialize the table to record the score, first string(sa) in row, second string(sb) in column
    # table index starts from 1
    sMatrix <- matrix(rep(0, (la+1)*(lb+1)), nrow=la+1, ncol=lb+1, byrow=TRUE)
    if(la != 0) {
        for(i in (1:la+1)) {
            sMatrix[i, 1] <- (i-1) * gap
        }
    }
    if(lb != 0) {
        for(j in (1:lb+1)) {
            sMatrix[1, j] <- (j-1) * gap
        }
    }
    #print(sMatrix)
    if(la != 0 & lb != 0) {
        for(i in (1:la+1)) {
            for(j in (1:lb+1)) {
                sScore <- sMatrix[i-1, j-1] + if (substr(sa, i-1, i-1) == substr(sb, j-1, j-1)) match else miss # string index starts from 0
                iScore <- sMatrix[i, j-1] + gap
                dScore <- sMatrix[i-1, j] + gap
                sMatrix[i, j] <- max(sScore, iScore, dScore)
            }
        }
    }
    #print(sMatrix)
    sMatrix
    #sMatrix[la+1, lb+1]
}

output <- function(sa, sb, sMatrix, i, j) {
    mIndex <- i+1
    nIndex <- j+1
    match <- 8
    miss <- -5
    gap <- -3
    #i, j plus one because of matrix index from 1
    if(mIndex == 1 | nIndex == 1) {
        if(mIndex > 1) {
            stringA <<- substr(sa, 0, i)
            stringB <<- paste(stringB, rep("-", i), collapse="", sep="")
        }
        if(nIndex > 1) {
            stringB <<- substr(sb, 0, j)
            stringA <<- paste(stringA, rep("-", j), collapse="", sep="")
        }
    } 
    else {
        if(sMatrix[mIndex, nIndex] == sMatrix[mIndex-1, nIndex-1] + if (substr(sa, i, i) == substr(sb, j, j)) match else miss) {
            output(sa, sb, sMatrix, i-1, j-1)
            stringA <<- paste(stringA, substr(sa, i, i), sep="")
            stringB <<- paste(stringB, substr(sb, j, j), sep="")
        } 
        else if(sMatrix[mIndex, nIndex] == sMatrix[mIndex-1, nIndex] + gap) {
            output(sa, sb, sMatrix, i-1, j)
            stringA <<- paste(stringA, substr(sa, i, i), sep="")
            stringB <<- paste(stringB, "-", sep="")
        } 
        else {
            output(sa, sb, sMatrix, i, j-1)
            stringA <<- paste(stringA, "-", sep="")
            stringB <<- paste(stringB, substr(sb, j, j), sep="")
        }
    }
}

stringA <<- ""
stringB <<- ""

str1="ATACATGTCT"
str2="GTACGTCGG"
#score("ATACATGTCT", "GTACGTCGG")
print(score(str1, str2))
output(str1, str2, score(str1, str2), nchar(str1), nchar(str2))

print(stringA)
print(stringB)
