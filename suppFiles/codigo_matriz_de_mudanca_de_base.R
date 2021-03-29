rm(list=ls())
library(splines)

chgbasismat <- function(knot, ord)
{
  dimmat <- length(knot) - ord
  answer <- matrix(0, nrow = dimmat, ncol = dimmat)
  for (j in 0:(ord-1))
  {
    brow <- splineDesign(knot, knot[1], ord, j)
    brow <- as.vector(brow/factorial(j))
    answer[j + 1, ] <- brow
  }
  nknot <- sort(-1*knot)
  for (j in 1:(dimmat - ord))
  {
    brow <- splineDesign(knot, knot[ord + j], ord, ord - 1)
    brow2 <- splineDesign(nknot, nknot[length(knot) - ord - (j - 1)],
                          ord, ord - 1)
    brow2 <- brow2[dimmat:1]
    brow <- brow + (-1)^ord * brow2
    brow <- as.vector(brow/factorial(ord - 1))
    answer[ord + j, ] <- brow
  }
  return(answer)
}

# Operacao de mudanca de base ilustrada no trabalho
solve(chgbasismat(c(0,0,0,0,1,2,3,4,4,4,4),4))%*%c(1,-1,1,-0.1,-0.8,2.2,-1.4)