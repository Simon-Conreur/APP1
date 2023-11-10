from math import sin, cos, atan, sqrt
import matplotlib.pyplot as plt
import numpy as np


"""
TODO :
Calculer les valeurs de hl et hr
G_0
"""


# Rotation point
def rotation(p, c, angle):
    """
    :param p: Point qui tourne (y,z)
    :param c: Centre de rotation (y, z)
    :return: rotation(p)
    """
    np.array(p)
    np.array(c)
    p -= c
    # matrice de rotation
    p[0] = cos(angle)*p[0] - sin(angle)*p[1]
    p[1] = sin(angle)*p[0] + cos(angle)*p[1]
    return [p[0], p[1]]


# constantes en SI ---------------------------------------------------------------->
#  constantes physiques
g = 9.81
rho = 1000
D = 0  # constante d'ammortissemnt
I = 1  # Inertie de la structure

# constantes dimensions
L = 0.60  # Largeur de la barge
h1 = 0  # hauteur de la barge
h2 = 0  # hauteur du mas
d1 = 0  # distance centre de la barge/ centre du mas

#  masses
m1 = 0  # masse barge
m2 = 0  # masse mat + 1er bras
m_charge = 0  # Masse le la charge portée
m_bras_grap = 0  # Masse bras 2 et grappin
m3 = m_charge + m_bras_grap  # masse grappin + 2 eme bras + charge
m_tot = m1 + m2 + m3 # masse totale de la structure
# masses relatives
m_c = m2  # masse charge totale
m_s = m1  # masse de la stucture

# constantes calculées ---------------------------------------------------------------->

# hauteur flottaison moyenne
hc = m_tot/(L**2*rho)

# centres de gravité en Y initiaux
G1_0 = []  # centre de gravité de la barge
G2_0 = []  # centre de gravité du mas
G3_0 = []  # centre de gravité du bras 2, de la charge utile et du grappin
# Centres relatifs
G_0 = [0, (m1*(-hc+h1/2)+m2*(h1-hc+h2/2))/(m1+m2)]  # centre de gravité /!\ Change enfontction de ce qu'on veut simuler
G_c_0 = []  # centre de gravité de la charge (ce qui va appliquer un couple sans être calculé dans le cnetre de gravité globale)
P_0 = [G_0[0], -hc/2]  # change en fonction de simulation

# valeur initiales des variables principales
angl_0 = 0  # angle initiale
v_angl_0 = 0  # vitesse angulaire initilale
hl_0 = hc  # hauteur d'eau à gauche ! dépend de l'angle et de la position de G par rapport aux bords de la barge/!\ à modifier
hr_0 = hc  # hauteur d'eau à droite

# constantes simulation
step = 0.001  # pas de calcul
end = 20  # temps de calcul

# Variables ---------------------------------------------------------------->

#  variables de simulations
t = np.arange(0, end, step)  # temps
angl = np.empty_like(t)  # angle d'inclinaison en fonction du temps
v_angl = np.empty_like(t)  # Vitesse angulaire
a_angl = np.empty_like(t)  # Acceleration angulaire
hl = np.empty_like(t)  # hauteur d'eau sur la gauche
hr = np.empty_like(t)  # hauteur d'eau sur la droite

# constantes de simulation
f_p = rho*g*L**2*hc  # force pousée
f_G = m_s*g  # gravité appliquée sur centre de gravité
f_charge = m_c*g  # gravité appliquée sur G3

c_r = np.empty_like(t)  # couple redressement gravité poussé charge

# centre de gravité
G = np.empty_like(t)
for i in range(len(G)):
    G[i] = [0,0]  # Remplir la liste de [0,0] : coordonées nulles
print("G initial vaut :", G)

#  Centre de poussée
P = np.empty_like(t)
for i in range(len(P)):
    P[i] = [0,0]  # Remplir la liste de [0,0] : coordonées nulles

# centre de gravité de la charge
G_c = np.empty_like(t)
for i in range(len(G_c)):
    G_c[i] = [0,0]  # Remplir la liste de [0,0] : coordonées nulles


c_a = m_c*g*G_c[0][0]  # couple destabilisteur (constant)


# simulation ---------------------------------------------------------------->
def simulation():
    angl[0] = angl_0
    v_angl[0] = v_angl_0
    a_angl[0] = 0
    dt = step
    P[0] = P_0
    for i in range(len(t)-1):
        G[i] = rotation(G_0, [0, 0], angl[i])
        #il faut encore calculer P
        c_r[i] = f_p*(P[i][0]-G[i][0])
        # F = Ia
        #  a_angl[i+1] = (c_a[i]+c_r[i])/I
        v_angl[i+1] = v_angl[i]*(1-D*dt/I) + (c_r[i] - c_a)*dt/I
        angl[i+1] = angl[i]+v_angl[i+1]*dt


# Representation graphique
def graphiques():
    plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(t, angl, label="theta")
    plt.legend()
    plt.subplot(2, 1, 2)
    plt.plot(t, v_angl, label="omega")
    plt.legend()
    plt.show()










