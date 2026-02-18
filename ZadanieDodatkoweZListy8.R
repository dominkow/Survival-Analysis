czas1 <- c(28,89,175,195,309,377,393,421,447,462,709,744,770,1106,1206)
delta1 <- c(1,1,1,1,1,0,0,0,0,1,0,0,0,0,0)

czas2 <- c(34,88,137,199,280,291,299,300,309,351,358,369,369,370,375,382,392,429,451,1119)
delta2 <- c(1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,0)

dane <- data.frame(
  czas   = c(czas1, czas2),
  status = c(delta1, delta2),
  grupa  = factor(c(rep("Niski",15), rep("Wysoki",20)))
)

t_i <- sort(unique(dane$czas[dane$status == 1]))
D <- length(t_i)                

k <- nlevels(dane$grupa)        
grupy <- levels(dane$grupa)

d_ij <- r_ij <- matrix(0, D, k)

for(i in 1:D){
  ti <- t_i[i]
  for(j in seq_len(k)){
    g <- grupy[j]
    r_ij[i, j] <- sum(dane$czas[dane$grupa == g] >= ti)
    d_ij[i, j] <- sum(dane$czas == ti & dane$status == 1 & dane$grupa == g)
  }
}

r_i <- rowSums(r_ij)      
d_i <- rowSums(d_ij)  

W_logrank <- rep(1, D)                     
W_gehan <- r_i                              
W_tarone <- sqrt(r_i)                        
W_peto <- cumprod(1 - d_i/(r_i + 1))       

Z_funkcja <- function(W){
  colSums( W * (d_ij - r_ij * d_i / r_i) )  
}

# Budujemy macierz E 

macierz_E <- function(W){
  E <- matrix(0, k, k)
  for(i in 1:D){
    w2   <- W[i]^2
    wspolny_czynnik <- d_i[i] * (r_i[i] - d_i[i]) / (r_i[i] - 1)  
    for(j in 1:k){
      E[j,j] <- E[j,j] + w2 * (r_ij[i,j] / r_i[i]) * (1 - (r_ij[i,j] / r_i[i])) * wspolny_czynnik #przekątna
      for(g in 1:k){
        if(g == j) next
        E[j,g] <- E[j,g] - w2 * r_ij[i,j] * r_ij[i,g] / (r_i[i]^2 )* wspolny_czynnik
      }
    }
  }
  E
}

test_dla_danej_wagi <- function(W){
  Z  <- Z_funkcja(W)
  E<- macierz_E(W)
  Macierz <- E[1, 1]  # Macierz tylko (k-1)x(k-1) = 1x1
  Z_k_1 <- Z[1]     # wektor Z (k−1 elementów)
  Z2 <- as.numeric( Z_k_1 %*% solve(Macierz) %*% t(Z_k_1) )
  
  p <- 1 - pchisq(Z2, df = k-1)
  
  print("Statystyka Z2")
  print(Z2)
  print("p-value")
  print(p)
  
}

print("Log-rank")
test_dla_danej_wagi(W_logrank)
print("Gehan-Breslow")
test_dla_danej_wagi(W_gehan)
print("Tarone-Ware")
test_dla_danej_wagi(W_tarone)
print("Peto-Peto")
test_dla_danej_wagi(W_peto)