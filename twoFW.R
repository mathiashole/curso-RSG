#################################################################
###-| MODELADO DE REDES ECOLÓGICAS--CLASE 2: REDES TRÓFICAS |-###
#################################################################

library(deSolve)

#---------------------------------------------------------------------------------#
#---| (1) Modelo Rosenzweig-MacArthur de consumidor-recurso (RM_CR) (pág. 22) |---#
#---------------------------------------------------------------------------------#

#----------------------#
#---| (1.1) MODELO |---#
#----------------------#

# R: biomasa/abundancia del recurso R
# C: biomasa/abundancia del consumidor C
# r: tasa de crecimiento per cápita de R
# K: capacidad de carga de R
# a: tasa máxima de ataque (consumo) de C sobre R (#R/tiempo)
# e: eficiencia de conversion de R a C ("qué tanta biomasa consumida de R se conveirte en biomasa de C")
# R0: biomasa de saturación media ("valor de R al cual la tasa máxima de ataque adquiere la mitad de su valor: a/2")
# m: tasa de muerte de C (#C/tiempo)

dR.dt <- r*R*(1-R/K) - (a*C*R)/(R + R0)    # Tasa de crecimiento del recurso R (asume respuesta funcional Tipo II)

dC.dt <- e*(a*C*R)/(R + R0) - m*C  # Tasa de crecimiento del consumidor C (asume respuesta funcional Tipo II)


#-----------------------#
#--| (1.2) DINÁMICAS |--#
#-----------------------#

RM_CR <- function(tiempo, y, parms){
  
  # Biomasas iniciales
  R <- y[1]
  C <- y[2]
  
  # Parámetros
  r <- parms[1]
  K <- parms[2]
  a <- parms[3]
  e <- parms[4]
  R0 <- parms[5]
  m <- parms[6]
  
  # Modelos
  dR.dt <- r*R*(1-R/K) - (a*C*R)/(R + R0)
  dC.dt <- e*(a*C*R)/(R + R0) - m*C
  
  #Salida del modelo (biomasas a lo largo del tiempo)
  return(list(c("R" = dR.dt, "C" = dC.dt)))
  
}

#--------------------------------------#
#--| (1.3) ESCENARIOS DE SIMULACIÓN |--#
#--------------------------------------#

par(mfrow = c(2, 3)) # Prepara el área de gráficos para mostrar 6 plots distribuidos en 3 filas y dos colmunas
par(mar = c(4,4,4,4))

Tiempo <- 1:100

# Esc.1: ¿qué pasa con la dinámica de R cuando no hay consumidor (C.t0 = 0) #
r <- 1
K <- 2
a <- 1.3
e <- 0.5
R0 <- 0.5
m <- 0.5

R.t0 <- 0.75
C.t0 <- 0

Esc.1 <- ode(y = c(R.t0, C.t0), times = Tiempo, func = RM_CR, parms = c(r, K, a, e, R0, m))
Esc.1
matplot(x = Esc.1[,"time"], y = Esc.1[, c("1", "2")], type = "l", ylab = "C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red"), lty = 1, ylim = c(0, 2.5), main = "Esc.1: R sin C (C.t0 = 0)")
legend(x = "topright", legend = c(expression("R"), expression("C")), 
       lty = c(1, 1), bty = "n", col = c("black", "red"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.2: ¿qué pasa con la dinámica de R cuando el consumidor SÍ está presente (C.t0 > 0) #

Tiempo <- 1:300
C.t0 <- 0.05

Esc.2 <- ode(y = c(R.t0, C.t0), times = Tiempo, func = RM_CR, parms = c(r, K, a, e, R0, m))
Esc.2
matplot(x = Esc.2[,"time"], y = Esc.2[, c("1", "2")], type = "l", ylab = "C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red"), lty = 1, ylim = c(0, 2.5), main = "Esc.2: R + C bajo (C.t0 > 0)")
legend(x = "topright", legend = c(expression("R"), expression("C")), 
       lty = c(1, 1), bty = "n", col = c("black", "red"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.3: ¿qué pasa con la dinámica de R cuando el consumidor está presente en mayor biomasa (C >> 0) #

C.t0 <- 0.5

Esc.3 <- ode(y = c(R.t0, C.t0), times = Tiempo, func = RM_CR, parms = c(r, K, a, e, R0, m))
Esc.3
matplot(x = Esc.3[,"time"], y = Esc.3[, c("1", "2")], type = "l", ylab = "C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red"), lty = 1, ylim = c(0, 2.5), main = "Esc.3: R + C alto (C.t0 >> 0)")
legend(x = "topright", legend = c(expression("R"), expression("C")), 
       lty = c(1, 1), bty = "n", col = c("black", "red"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.4: ¿qué pasa con la dinámica de R y C si aumenta la fuerza de interacción (FI) entre R y C 
# (mayor eficiencia de conversión: e > 0.5)? # McCann 2011 (CAP. 5)

Tiempo <- 1:300
e <- 0.6 # antes, e = 0.5

Esc.4 <- ode(y = c(R.t0, C.t0), times = Tiempo, func = RM_CR, parms = c(r, K, a, e, R0, m))
Esc.4
matplot(x = Esc.4[,"time"], y = Esc.4[, c("1", "2")], type = "l", ylab = "C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red"), lty = 1, ylim = c(0, 2.5), main = "Esc.4: FI media (e = 0.6)")
legend(x = "topright", legend = c(expression("R"), expression("C")), 
       lty = c(1, 1), bty = "n", col = c("black", "red"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.5: ¿y si FI entre R y C es aún mayor (eficiencia de conversión alta: e >> 0.5)? #

Tiempo <- 1:300
e <- 0.65 # antes, e = 0.5 

Esc.5 <- ode(y = c(R.t0, C.t0), times = Tiempo, func = RM_CR, parms = c(r, K, a, e, R0, m))
Esc.5
matplot(x = Esc.5[,"time"], y = Esc.5[, c("1", "2")], type = "l", ylab = "C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red"), lty = 1, ylim = c(0, 2.5), main = "Esc.5: FI moderada (e = 0.65)")
legend(x = "topright", legend = c(expression("R"), expression("C")), 
       lty = c(1, 1), bty = "n", col = c("black", "red"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.6: ¿y si FI entre R y C es demasiado grande (eficiencia de conversión muy alta: e >>> 0.5)? #

Tiempo <- 1:300
e <- 0.7 # antes, e = 0.5 

Esc.6 <- ode(y = c(R.t0, C.t0), times = Tiempo, func = RM_CR, parms = c(r, K, a, e, R0, m))
Esc.6
matplot(x = Esc.6[,"time"], y = Esc.6[, c("1", "2")], type = "l", ylab = "C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red"), lty = 1, ylim = c(0, 2.5), main = "Esc.6: FI alta (e = 0.7)")
legend(x = "topright", legend = c(expression("R"), expression("C")), 
       lty = c(1,1,1), bty = "n", col = c("black", "red"), lwd = 3)
abline(h = c(0, 2), lty = 2)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

#--------------------------------------------------------------------------#
#--| Módulos complejos: cadena trófica & cadena trófica con competencia |--#
#--------------------------------------------------------------------------# McCann 2011 (CAP. 7)


#---------------------------------------------------------------------------------------------#
#---| (2) Modelo Rosenzweig-MacArthur de depredador-consumidor-recurso (RM_CT) (pág. 108) |---#
#---------------------------------------------------------------------------------------------#

#----------------------#
#---| (2.1) MODELO |---#
#----------------------#

# R: biomasa/abundancia del recurso R
# C: biomasa/abundancia del consumidor C
# P: biomasa/abundancia del depredador P
# r: tasa de crecimiento per cápita de R
# K: capacidad de carga de R
# aC: tasa máxima de ataque (consumo) de C sobre R (#R/tiempo)
# aP: tasa máxima de ataque (consumo) de  P sobre C (#C/tiempo)
# eC: eficiencia de conversion de R a C ("qué tanta biomasa consumida de R se conveirte en biomasa de C")
# eP: eficiencia de conversion de C a P ("qué tanta biomasa consumida de C se conveirte en biomasa de P")
# R0: biomasa de saturación media para C ("valor de R al cual la tasa máxima de ataque de C adquiere la mitad de su valor: aC/2")
# C0: biomasa de saturación media ("valor de C al cual la tasa máxima de ataque de P adquiere la mitad de su valor: aP/2")
# mC: tasa de muerte de C (#C/tiempo)
# mP: tasa de muerte de P (#P/tiempo)

dR.dt <- r*R*(1-R/K) - (aC*C*R)/(R + R0)    # Tasa de crecimiento del recurso R (asume respuesta funcional Tipo II)

dC.dt <- eC*(aC*C*R)/(R + R0) - mC*C - (aP*C*P)/(C + C0) # Tasa de crecimiento del consumidor C (asume respuesta funcional Tipo II)

dP.dt <- eP*(aP*P*C)/(C + C0) - mP*P  # Tasa de crecimiento del depredador P (asume respuesta funcional Tipo II)

#-----------------------#
#--| (2.2) DINÁMICAS |--#
#-----------------------#

RM_CT <- function(tiempo, y, parms){

  # Abundancias iniciales
  R <- y[1]
  C <- y[2]
  P <- y[3]
  
  # Parámetros
  r <- parms[1]
  K <- parms[2]
  aC <- parms[3]
  aP <- parms[4]
  eC <- parms[5]
  eP <- parms[6]
  R0 <- parms[7]
  C0 <- parms[8]
  mC <- parms[9]
  mP <- parms[10]
  
  # Modelos
  dR.dt <- r*R*(1-R/K) - (aC*C*R)/(R + R0)
  dC.dt <- eC*(aC*C*R)/(R + R0) - mC*C - (aP*C*P)/(C + C0)
  dP.dt <- eP*(aP*P*C)/(C + C0) - mP*P
  
  #Salida del modelo (abundancias a lo largo del tiempo)
  return(list(c("R" = dR.dt, "C" = dC.dt, "P" = dP.dt)))
  
}

#--------------------------------------#
#--| (2.3) ESCENARIOS DE SIMULACIÓN |--#
#--------------------------------------#

par(mfrow = c(1, 2)) # Prepara el área de gráficos para mostrar 6 plots distribuidos en 3 filas y dos colmunas
par(mar = c(4, 4, 4, 4))

# Esc.7: recordemos qué pasa con la dinámica C-R cuando la FI entre ellos es alta (eC = 0.7) y no hay depredador P (P.t0 = 0) # 
# (NOTA: estamos recreamos el Esc. 6 pero ahora con un modelo de cadena trófica)

Tiempo <- 1:300

r <- 1
K <- 2
aC <- 1.3
aP <- 1.3
eC <- 0.7
eP <- 0.8
R0 <- 0.5
C0 <- 0.5
mC <- 0.5
mP <- 0.5

R.t0 <- 0.75
C.t0 <- 0.5
P.t0 <- 0

matplot(x = Esc.6[,"time"], y = Esc.6[, c("1", "2")], type = "l", ylab = "C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red"), lty = 1, ylim = c(0, 2.5), main = "Esc.6: FI alta (e = 0.7)")
legend(x = "topright", legend = c(expression("R"), expression("C")), 
       lty = c(1,1,1), bty = "n", col = c("black", "red"), lwd = 3)
abline(h = c(0, 2), lty = 2)

Esc.7 <- ode(y = c(R.t0, C.t0, P.t0), times = Tiempo, func = RM_CT, parms = c(r, K, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.7
matplot(x = Esc.7[,"time"], y = Esc.7[, c("1", "2", "3")], type = "l", ylab = "P-C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.7: C & R sin P (P.t0 = 0)")
legend(x = "topright", legend = c(expression("R"), expression("C"), expression("P")), 
       lty = c(1, 1, 1), bty = "n", col = c("black", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.8: ¿qué pasa con la dinámica C-R (inestable) en presencia del depredador P (P.t0 > 0) # 

par(mfrow = c(1, 3)) 
par(mar = c(4, 4, 4, 4))

P.t0 <- 0.25

matplot(x = Esc.7[,"time"], y = Esc.7[, c("1", "2", "3")], type = "l", ylab = "P-C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.7: C & R sin P (P.t0 = 0)")
legend(x = "topright", legend = c(expression("R"), expression("C"), expression("P")), 
       lty = c(1, 1, 1), bty = "n", col = c("black", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

Esc.8 <- ode(y = c(R.t0, C.t0, P.t0), times = Tiempo, func = RM_CT, parms = c(r, K, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.8
matplot(x = Esc.8[,"time"], y = Esc.8[, c("1", "2", "3")], type = "l", ylab = "P-C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.8: C & R con P (P.t0 > 0)")
legend(x = "topright", legend = c(expression("R"), expression("C"), expression("P")), 
       lty = c(1, 1, 1), bty = "n", col = c("black", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.9: ¿qué pasa con la dinámica P-C-R  al aumentar la FI entre P-C (eP > 0.8) # 

eP <- 1 # antes, ep = 0.8

Esc.9 <- ode(y = c(R.t0, C.t0, P.t0), times = Tiempo, func = RM_CT, parms = c(r, K, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.9
matplot(x = Esc.9[,"time"], y = Esc.9[, c("1", "2", "3")], type = "l", ylab = "P-C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.9: FI entre P-C alta (eP = 1)")
legend(x = "topright", legend = c(expression("R"), expression("C"), expression("P")), 
       lty = c(1, 1, 1), bty = "n", col = c("black", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

#-------------------------------------------------------------------------------------------------------------------------------#
#---| (3) Modelo Rosenzweig-MacArthur de depredador-consumidor-recurso con competencia entre recursos (R1 y R2) (RM_CT.Comp)|---#
#-------------------------------------------------------------------------------------------------------------------------------#

#----------------------#
#---| (3.1) MODELO |---#
#----------------------#

# R1: abundancia/biomasa del recurso R1
# R2: abundancia/biomasa del recurso R2
# C: abundancia/biomasa del consumidor C
# P: abundancia/biomasa del depredador P
# r: tasa de crecimiento per cápita de R1 y R2
# K: capacidad de carga de R1 y R2
# alpha21: efecto de R2 sobre la tasa de crecimiento per cápita (r) de R1 ("qué tan importante es la competencia (per cápita) interespecífica en relación a la intraespecífica")
# alpha12: efecto de R1 sobre la tasa de crecimiento per cápita (r) de R2 ("qué tan importante es la competencia (per cápita) interespecífica en relación a la intraespecífica")
# aC: tasa máxima de ataque (consumo) de C sobre R1 y R2 (#R/tiempo)
# aP: tasa máxima de ataque (consumo) de  P sobre C (#C/tiempo)
# eC: eficiencia de conversion de R1 y R2 a C ("qué tanta biomasa consumida de R1 y R2 se conveirte en biomasa de C")
# eP: eficiencia de conversion de C a P ("qué tanta biomasa consumida de C se conveirte en biomasa de P")
# R0: biomasa de saturación media para C ("valor de R1 y R2 al cual la tasa máxima de ataque de C adquiere la mitad de su valor: aC/2")
# C0: biomasa de saturación media ("valor de C al cual la tasa máxima de ataque de P adquiere la mitad de su valor: aP/2")
# mC: tasa de muerte de C (#C/tiempo)
# mP: tasa de muerte de P (#P/tiempo)

dR1.dt <- r1*R1*(K1-R1-alpha21*R2)/K1 - (aC*C*R1)/(R1 + R0)  # Tasa de crecimiento del recurso R1

dR2.dt <- r2*R2*(K2-R2-alpha12*R1)/K2                        # Tasa de crecimiento del recurso R2

dC.dt <- eC*(aC*C*R1)/(R1 + R0) - mC*C - (aP*C*P)/(C + C0)   # Tasa de crecimiento del consumidor C

dP.dt <- eP*(aP*P*C)/(C + C0) - mP*P                         # Tasa de crecimiento del deprdador P

#-----------------------#
#--| (3.2) DINÁMICAS |--#
#-----------------------#

RM_CT.Comp <- function(tiempo, y, parms){
  
  # Abundancias iniciales
  R1 <- y[1]
  R2 <- y[2]
  C <- y[3]
  P <- y[4]
  
  # Parámetros
  r1 <- parms[1]
  r2 <- parms[2]
  K1 <- parms[3]
  K2 <- parms[4]
  alpha21 <- parms[5]
  alpha12 <- parms[6]
  aC <- parms[7]
  aP <- parms[8]
  eC <- parms[9]
  eP <- parms[10]
  R0 <- parms[11]
  C0 <- parms[12]
  mC <- parms[13]
  mP <- parms[14]
  
  # Modelos
  dR1.dt <- r1*R1*(K1-R1-alpha21*R2)/K1 - (aC*C*R1)/(R1 + R0) 
  dR2.dt <- r2*R2*(K2-R2-alpha12*R1)/K2   
  dC.dt <- eC*(aC*C*R1)/(R1 + R0) - mC*C - (aP*C*P)/(C + C0) 
  dP.dt <- eP*(aP*P*C)/(C + C0) - mP*P
  
  #Salida del modelo (abundancias a lo largo del tiempo)
  return(list(c("R1" = dR1.dt, "R2" = dR2.dt, "C" = dC.dt, "P" = dP.dt)))
  
}

#--------------------------------------#
#--| (3.3) ESCENARIOS DE SIMULACIÓN |--#
#--------------------------------------#

par(mfrow = c(1, 3)) 
par(mar = c(4, 4, 4, 4))

# Esc.10: recreemos el Esc. 8 (cadena trófica) pero con el modelo MR_CT.Comp en ausencia de R2 (R2.t0 = 0)) # 

Tiempo <- 1:300

r1 <- 1
r2 <- 1
K1 <- 2
K2 <- 2
alpha21 <- 0 
alpha12 <- 0 
aC <- 1.3
aP <- 1.3
eC <- 0.7
eP <- 0.8
R0 <- 0.5
C0 <- 0.5
mC <- 0.5
mP <- 0.5

R1.t0 <- 0.75
R2.t0 <- 0 # <---
C.t0 <- 0.5
P.t0 <- 0.25

matplot(x = Esc.8[,"time"], y = Esc.8[, c("1", "2", "3")], type = "l", ylab = "P-C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.8: C & R con P (P.t0 > 0)")
legend(x = "topright", legend = c(expression("R"), expression("C"), expression("P")), 
       lty = c(1, 1, 1), bty = "n", col = c("black", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

Esc.10 <- ode(y = c(R1.t0, R2.t0, C.t0, P.t0), times = Tiempo, func = RM_CT.Comp, parms = c(r1, r2, K1, K2, alpha21, alpha12, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.10
matplot(x = Esc.10[,"time"], y = Esc.10[, c("1", "2", "3", "4")], type = "l", ylab = "P-C-R", xlab = "Tiempo", 
        lwd = 3, col = c("black", "orange", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.10: P, C, R1 sin R2")   
legend(x = "topright", legend = c(expression("R1"), expression("R2"), expression("C"), expression("P")), 
       lty = c(1, 1, 1, 1), bty = "n", col = c("black", "orange", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.11: ¿qué pasa con la dinámica P-C-R1 en presencia del recurso R2 (R2.t0 > 0) pero sin competencia entre R1 y R2 (alpha21 = alpha12 = 0) # 

R2.t0 <- 0.1 # antes era 0

alpha21 <- 0  
alpha12 <- 0    

Esc.11 <- ode(y = c(R1.t0, R2.t0, C.t0, P.t0), times = Tiempo, func = RM_CT.Comp, parms = c(r1, r2, K1, K2, alpha21, alpha12, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.11
matplot(x = Esc.11[,"time"], y = Esc.11[, c("1", "2", "3", "4")], type = "l", ylab = "P-C-R1 & R2", xlab = "Tiempo", 
        lwd = 3, col = c("black", "orange", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.11: P, C, R1 & R2 sin competencia")   
legend(x = "topright", legend = c(expression("R1"), expression("R2"), expression("C"), expression("P")), 
       lty = c(1, 1, 1, 1), bty = "n", col = c("black", "orange", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.12: ¿qué pasa con la dinámica P-C-R1 en presencia de R2 (R2.t0 > 0) con competencia entre R1 y R2 (alpha21 > 0; alpha12 > 0) # 
# a medida que la competencia interespecífica aumenta                                                                              #

par(mfrow = c(1, 4))
par(mar = c(4, 4, 4, 4))

matplot(x = Esc.11[,"time"], y = Esc.11[, c("1", "2", "3", "4")], type = "l", ylab = "P-C-R1 & R2", xlab = "Tiempo", 
        lwd = 3, col = c("black", "orange", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.11: P, C, R1 & R2 sin competencia")   
legend(x = "topright", legend = c(expression("R1"), expression("R2"), expression("C"), expression("P")), 
       lty = c(1, 1, 1, 1), bty = "n", col = c("black", "orange", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.12.1: competencia interespecífica débil

alpha21 <- 1
alpha12 <- 2

Esc.12.1 <- ode(y = c(R1.t0, R2.t0, C.t0, P.t0), times = Tiempo, func = RM_CT.Comp, parms = c(r1, r2, K1, K2, alpha21, alpha12, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.12.1
matplot(x = Esc.12.1[,"time"], y = Esc.12.1[, c("1", "2", "3", "4")], type = "l", ylab = "P-C-R1 & R2", xlab = "Tiempo", 
        lwd = 3, col = c("black", "orange", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.12.1: P, C, R1 & R2 (alpha12 > alpha21)")    
legend(x = "topright", legend = c(expression("R1"), expression("R2"), expression("C"), expression("P")), 
       lty = c(1, 1, 1, 1), bty = "n", col = c("black", "orange", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.12.2: competencia interespecífica media

alpha21 <- 3
alpha12 <- 2

Esc.12.2 <- ode(y = c(R1.t0, R2.t0, C.t0, P.t0), times = Tiempo, func = RM_CT.Comp, parms = c(r1, r2, K1, K2, alpha21, alpha12, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.12.2
matplot(x = Esc.12.2[,"time"], y = Esc.12.2[, c("1", "2", "3", "4")], type = "l", ylab = "P-C-R1 & R2", xlab = "Tiempo", 
        lwd = 3, col = c("black", "orange", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.12.2: P, C, R1 & R2 (alpha12 < alpha21)")   
legend(x = "topright", legend = c(expression("R1"), expression("R2"), expression("C"), expression("P")), 
       lty = c(1, 1, 1, 1), bty = "n", col = c("black", "orange", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

# Esc.12.3: competencia interespecífica alta

alpha21 <- 3.5
alpha12 <- 2

Esc.12.3 <- ode(y = c(R1.t0, R2.t0, C.t0, P.t0), times = Tiempo, func = RM_CT.Comp, parms = c(r1, r2, K1, K2, alpha21, alpha12, aC, aP, eC, eP, R0, C0, mC, mP))
Esc.12.3
matplot(x = Esc.12.3[,"time"], y = Esc.12.3[, c("1", "2", "3", "4")], type = "l", ylab = "P-C-R1 & R2", xlab = "Tiempo", 
        lwd = 3, col = c("black", "orange", "red", "darkgreen"), lty = 1, ylim = c(0, 2.5), main = "Esc.12.3: P, C, R1 & R2 (alpha12 << alpha21)")   
legend(x = "topright", legend = c(expression("R1"), expression("R2"), expression("C"), expression("P")), 
       lty = c(1, 1, 1, 1), bty = "n", col = c("black", "orange", "red", "darkgreen"), lwd = 3)
abline(h = c(0, 2), lty = 2)

####################################
### FIN DEL SCRIPT: ¡GRACIAS! :) ###
####################################




