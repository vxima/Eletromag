#includes
import numpy as np
import matplotlib.pyplot as plt
#constantes (tudo ja esta em SI)
e0 = 8.85e-12
#VARIAVEIS 
N1 = 6
N2 = 5
N = N1 + N2
#permissividade = Er*E0
e1 = 2
e2 = 4
e1 = e1*e0
e2 = e2*e0
#comprimento dieletrico
d1 = 1e-3
d2 = 1e-3
d = d1 + d2
l1 = d1/N1
l2 = d2/N2
L = 2e-2
#potencial nas bordas
V0 = 0
Vn = 1
#K*V = D
def makeK_e1(N1 , N , e1 , l):
    #Inicializar matriz com 0s
    K1 = np.zeros((N+1, N+1))
    #Adicionar equações lineares, inserindo o bloco de 4 equaçoes para cada iteraçao
    for i in range(0 , N1):
        #diagonal principal
        K1[i][i] = K1[i][i]  + (-e1/l)
        K1[i+1][i+1] = K1[i+1][i+1]  + (-e1/l)
        #diagonal oposta
        K1[i+1][i] = K1[i+1][i] + (e1/l)
        K1[i][i+1] = K1[i][i+1] + (e1/l)
    return K1

def makeK_e2(N2 , N , e2 , l):
    #Inicializar matriz com 0s
    K2 = np.zeros((N+1, N+1))
    #Adicionar equações lineares, inserindo o bloco de 4 equaçoes para cada iteraçao
    for i in range((N-N2) , N): #continuando de N1
        #diagonal principal
        K2[i][i] = K2[i][i]  + (-e2/l)
        K2[i+1][i+1] = K2[i+1][i+1]  + (-e2/l)
        #diagonal oposta
        K2[i+1][i] = K2[i+1][i] + (e2/l)
        K2[i][i+1] = K2[i][i+1] + (e2/l)
    return K2
K1 = makeK_e1(N1 , N , e1 , l1 ) #calcula a matriz de cada dieletrico
K2 = makeK_e2(N2 , N , e2 , l2 )
K = K1+K2 #soma as matrizes dos dieletrico , tendo a matriz global

subK = K[1:N , 1:N] #pegar submatriz cortando a primeira e ultima linha da matriz
def makeD(N  , V0 , Vn , K):
    Dn = np.zeros(N+1)
    Dn[1] = Dn[0] - K[1][0]*V0 #substituo o valor de V0 e altero no D[1]
    Dn[N-1] = Dn[N] - K[N-1][N]*Vn #substituo o valor de Vn e altero no D[n-1]
    Dn = Dn[1:N] #pega o vetor D com os cortes 
    
    return Dn

D = makeD(N  , V0 , Vn , K)
V = np.linalg.solve(subK , D) #resolver sistema linar K*V = D

##capacitancia teorica
C1 = e1*(L*L)/d1
C2 = e2*(L*L)/d2
#capacitancia em serie
Ceq = (C1*C2)/(C2+C1)
##capacitancia numerica
Qsuperficie = L*L*(e2*(Vn - V[len(V) - 1 ]))/l2
Vtotal = Vn - V0
Cnumerica = Qsuperficie/Vtotal
print("Valor de Capacitancia numerica: " , Cnumerica)
print("Valor de Capacitancia teórica: " , Ceq)
#plotar grafico
xlab = []
ylab = []
x = 0
xlab.append(x)
for i in range(0 , N1): #seleciono os nós para o eixo x
    x = x + l1
    xlab.append(x)
for i in range(0 , N2):
    x = x + l2
    xlab.append(x)
xlab = [x*1e3 for x in xlab] # converto os valores do eixo x para milimetro
ylab.append(V0)
ylab.extend(V.tolist())
ylab.append(Vn)
print("Valores de Vi: " , ylab)

plt.xlabel('z(mm)')
plt.ylabel('V(z)')
plt.title('Potencial eletrico no capacitor')
plt.plot(xlab , ylab , color = 'red' , label="V(z)")
plt.legend()
plt.show()
