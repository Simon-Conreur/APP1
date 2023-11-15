"""
Groupe 11.14 LEPL 1501-APP1
Programme de simulation de la barge flottante
TODO :
Vérifier que les données sont cohérentes avant de lancer la simulation
Fonction Inertie
"""

from math import sin, cos, atan, sqrt, tan
import matplotlib.pyplot as plt
import numpy as np


# Rotation point
def rotation(p, c, angle):
    """
    :param p: Point qui tourne (y,z)
    :param c: Centre de rotation (y, z)
    :param angle: Angle de rotation (radiants)
    :return: les coordonnées du points p apres rotation de centre c d'angle angle
    """
    p = np.array(p)
    c = np.array(c)
    p -= c
    y_p = p[0]
    z_p = p[1]
    # matrice de rotation
    p[0] = cos(angle)*y_p - sin(angle)*z_p
    p[1] = sin(angle)*y_p + cos(angle)*z_p

    p += c
    return [p[0], p[1]]


# r1 to r2
def translation(point, vect):
    """
    Déplace un point
    :param point: point à déplacer
    :param vect: vecteur de déplacement
    :return: les coordonnées de point après translation de vecteur vect
    """
    point[0] += vect[0]
    point[1] += vect[1]


def r1_2_0(point):
    """
    Repère 1 vers repère 2 à l'angle initial
    :param point: point à modifier
    :return: Le point dans r2
    """
    point[0], point[1] = point[0] + v1_2[0], point[1] + v1_2[1]  # translation vers R1
    rotation(point, [0,0], angl_0)  # rotation selon l'angle initial
    return point


def g_trapeze(theta, l, hc, v):
    """
    nb: C est le centre de poussée voir notes pour plus d'infos
    :param theta: Angle de rotation
    :param l: largeur de la barge
    :param hc: hauteur de l'eau au niveau du centre de gravité (hc)
    :param v: composante y du vecteur R1->R2 (= v1_2)
    :return: le centre de poussée P
    """
    hl = hc + tan(theta)*(l + (-v))
    hr = hc - tan(theta)*(l - (-v))
    lc_x = l*(hl + 2*hr)/(3*(hl + hr))
    lc_y = (hl**2 + hl*hr + hr**2)/(3*(hl + hr))
    Yc = -l/2 + lc_x + v
    Zc = -hc + lc_y
    C = rotation([Yc, Zc], [0,0], theta)
    return C


# /!\ les points sont donnés dans le repère 1
# constantes en SI ---------------------------------------------------------------->
#  constantes physiques
g = 9.81  # accélération due à la gravité de la Terre
rho = 1000  # masse volumique de l'eau
D = 1  # constante d'ammortissemnt

# constantes dimensions
L = 0.60  # Largeur de la barge
L2 = 0.05  # largeur du mat
h1 = 0.07  # hauteur de la barge
h2 = 0.50  # hauteur du mat
d1 = 0.10  # distance centre de la barge/ centre du mat

#  masses
m1 = 0.700  # masse barge
m2 = 0.250  # masse mat + 1er bras
m_charge = 0.200  # Masse le la charge portée
m_bras_grap = 0.100  # Masse bras 2 et grappin
m3 = m_charge + m_bras_grap  # masse grappin + 2 eme bras + charge
# masses relatives
m_tot = m1 + m2 + m3  # masse totale de la structure
m_c = m2  # masse charge totale
m_s = m1  # masse de la stucture

# constantes calculées ---------------------------------------------------------------->

# hauteur flottaison moyenne
hc = m_tot/(L**2*rho)

# centres de gravité en Y initiaux
G1_0 = [0, -hc+h1/2]  # centre de gravité de la barge
G2_0 = [d1, -hc +h1 + h2/2]  # centre de gravité du mas
G3_0 = []  # centre de gravité du bras 2, de la charge utile et du grappin

# Vecteur de translation pour passer du repère 1 au repère 2
# Les points sont calculés pour correspondre au repère 2 apd ici ---------------------------->
v1_2 = [-d1*m2/m_s, 0] # vecteur 1 vers 2
v2_1 = [d1*m2/m_s, 0]  # vecteur 2 vers 1

# Inertie de la structure
I = (m1/12)*(h1**2 + L**2 + (hc - h1/2)**2) + m2*(h2**2 + L2**2 + d1**2 + (hc - h1 - h2/2)**2)

# Centres relatifs
G_0 = [0, (m1*(-hc+h1/2)+m2*(h1-hc+h2/2))/(m1+m2)]  # centre de gravité /!\ Change enfontction de ce qu'on veut simuler
G_c_0 = []  # centre de gravité de la charge (ce qui va appliquer un couple sans être calculé dans le cnetre de gravité globale)
P_0 = [G_0[0], -hc/2]  # change en fonction de simulation

# valeur initiales des variables principales
angl_0 = 0  # angle initial
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

c_r = np.empty_like(t)  # couple redressement gravité poussée

c_a = np.empty_like(t)  # Couple appliqué par la charge

c_d = np.empty_like(t)  # Couple d'amortidssenment(v_angl)

# centre de gravité
tab = []
for i in range(len(t)):
    tab.append([0, 0])  # Remplir la liste de [0,0] : coordonées nulles
G = np.array(tab)

# centre de poussée rempli de coordonnées nulles
P = np.array(tab)
# centre de gravité de la charge rempli de coordonnées nulles
G_c = np.array(tab)
print("G initial vaut :", G)

# Déplacement des points initiaux vers le repère
G_0 = r1_2_0(G_0)
G_c_0 = r1_2_0(G_c_0)


# simulation ---------------------------------------------------------------->
def simulation():
    # G_c[0] = G_c_0 inutile car G_c[0] est défini dans la premiere itération de la boucle
    angl[0] = angl_0
    v_angl[0] = v_angl_0
    a_angl[0] = 0
    dt = step
    P[0] = P_0
    for i in range(len(t)-1):
        # Points au temps i
        G[i] = rotation(G_0, [0, 0], angl[i])  # centre de gravité
        G_c[i] = rotation(G_c_0, [0, 0], angl[i])  # centre de gravité de la charge
        P[i] = g_trapeze(angl[i], L, hc, v1_2[0])  # a vérifier car la manip est complquée

        # Couple au temps i
        c_a[i] = m_c * g * G_c[i][0]  # couple destabilisteur en fonction de G_c_Y
        c_r[i] = f_p*(P[i][0]-G[i][0])

        # vitesse et accélération au temps i+1  avec la méthode d'Euler depuis F = Ia
        v_angl[i+1] = v_angl[i] + (-D*v_angl[i] + c_r[i] - c_a[i])*dt/I
        angl[i+1] = angl[i]+v_angl[i]*dt


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










