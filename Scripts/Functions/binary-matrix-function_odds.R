bin.fun <- function(r1, r2, t, A, B, C){
  
  # Time matrix
  A.bin = 1*(A <= t)
  # Distance matrix
  B.bin = 1*(B < r2 & B > r1)
  
  
  ## Estimation
  a = sum(A.bin * B.bin * C, na.rm = TRUE) 
  b = sum(A.bin * B.bin * (1 - C), na.rm = TRUE)
  c = sum(A.bin * C, na.rm = TRUE)
  d = sum(A.bin * (1 - C), na.rm = TRUE)
  
  out <- data.frame(a = a,
                    b = b, 
                    c = c,
                    d = d,
                    tau = (a/b)/(c/d),
                    r_lower = r1,
                    r_upper = r2)
  return(out)
}
