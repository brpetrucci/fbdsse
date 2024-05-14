modelComb <- c(1, 2, 3)
parComb <- c(1, 2, 3, 4, 5, 6)
traitComb <- c(1, 2, 3, 4)
psiComb <- c(1, 2, 3)

count <- 1
for (m in modelComb) {
  if (m == 1) {
    for (p in parComb) {
      for (r in psiComb) {
        lines <- c(paste0("modelComb <- ", m),
                   paste0("parComb <- ", p),
                   paste0("traitComb <- 1"),
                   paste0("psiComb <- ", r))
        
        fileConn <- file(paste0("refs_", count, ".Rev"))
        writeLines(lines, fileConn)
        close(fileConn)
        
        count <- count + 1
      }
    }
  }
  
  else if (m == 2) {
    for (p in parComb) {
      for (t in traitComb) {
        lines <- c(paste0("modelComb <- ", m),
                   paste0("parComb <- ", p),
                   paste0("traitComb <- ", t),
                   paste0("psiComb <- 1"))
        
        fileConn <- file(paste0("refs_", count, ".Rev"))
        writeLines(lines, fileConn)
        close(fileConn)
        
        count <- count + 1
      }
    }
  }
  
  else {
    for (p in parComb) {
      for (t in traitComb) {
        for (r in psiComb) {
          lines <- c(paste0("modelComb <- ", m),
                     paste0("parComb <- ", p),
                     paste0("traitComb <- ", t),
                     paste0("psiComb <- ", r))
          fileConn <- file(paste0("refs_", count, ".Rev"))
          writeLines(lines, fileConn)
          close(fileConn)
          
          count <- count + 1
        }
      }
    }
  }
}
