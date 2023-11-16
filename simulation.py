"""
Groupe 11.14 LEPL 1501-APP1
Programme de simulation de la barge flottante
TODO :
Calculer mouvement charge
Calculer Submersion Decollement
Calculer les énergies
"""

from math import sin, cos, atan, sqrt, tan, pi
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


def r1_2_0(point):
    """
    Repère 1 vers repère 2 à l'angle initial plus de détails sur mes notes
    :param point: point à déplacer dans R1
    :return: Le point dans R2 (Translation de vecteur (X_G, 0) PUIS rotation de centre (0,0) et d'angle angl_0)
    """
    point[0], point[1] = point[0] + v1_2[0], point[1] + v1_2[1]  # translation vers R1
    rotation(point, [0,0], angl_0)  # rotation selon l'angle initial
    return point


def g_trapeze(theta, l, hc, v):
    """
    nb: C est le centre de poussée
    Calculs détaillés dans les slides du CM de projet S8
    :param theta: Angle de rotation
    :param l: largeur de la barge
    :param hc: hauteur de l'eau au niveau du centre de gravité (hc)
    :param v: composante y du vecteur R1->R2 (= v1_2)
    :return: le centre de poussée P du trapèze de bases hl, hr et de hauteur l
    """
    hl = hc + tan(theta)*(l/2 + (-v))
    hr = hc - tan(theta)*(l/2 - (-v))
    lc_x = l*(hl + 2*hr)/(3*(hl + hr))
    lc_y = (hl**2 + hl*hr + hr**2)/(3*(hl + hr))
    Yc = -l/2 + lc_x + v
    Zc = -hc + lc_y
    C = rotation([Yc, Zc], [0,0], theta)
    return C


# constantes en SI ---------------------------------------------------------------->
#  constantes physiques
g = 9.81  # accélération due à la gravité de la Terre
rho = 1000  # masse volumique de l'eau
D = 2  # constante d'ammortissemnt

# constantes dimensions
L = 0.60  # Largeur de la barge
L2 = 0.05  # largeur du mat
h1 = 0.10  # hauteur de la barge
h2 = 0.50  # hauteur du mat
d1 = 0.25  # distance centre de la barge/ centre du mat

#  masses
m1 = 0.700  # masse barge
m2 = 0.2  # masse mat + 1er bras
m_charge = 0.100  # Masse le la charge portée
m_bras_grap = 0.200  # Masse bras 2 et grappin
m3 = m_charge + m_bras_grap  # masse grappin + 2 eme bras + charge

# masses relatives
# (ces données sont à privilégier lors des calculs
# car si on change la configuration de la grue il suffit de changer les variables ici)
m_tot = m1 + m2 + m3  # masse totale de la structure
m_c = 0.3  # masse charge totale
m_s = m1 + m2  # masse de la stucture

# constantes calculées ---------------------------------------------------------------->

# hauteur flottaison moyenne au niveau du centre de gravité G
hc = m_tot/(L**2*rho)

# Inertie de la structure
I = (m1/12)*(h1**2 + L**2 + (hc - h1/2)**2) + m2*(h2**2 + L2**2 + d1**2 + (hc - h1 - h2/2)**2)

# Constantes valeurs initiales

# centres de gravité initiaux
G1_0 = [0, -hc+h1/2]  # centre de gravité de la barge
G2_0 = [d1, -hc + h1 + h2/2]  # centre de gravité du mas
G3_0 = [0.70, 0.50]  # centre de gravité du bras 2, de la charge utile et du grappin

# Variables cinématiques initiales
angl_0 = 0  # angle initial
v_angl_0 = 0  # vitesse angulaire initiale

# Centres relatifs initiaux
G_0 = [(G1_0[0]*m1 + G2_0[0]*m2)/m_s, (G1_0[1]*m1 + (G2_0[1])*m2)/m_s]  # centre de gravité de la structure initiale
G_c_0 = G3_0  # centre de gravité de la charge (ce qui applique c_a mais pas calculé dans le centre de gravité global)
"""P_0 = [0, hc/2]  # change en fonction de simulation cette ligne est inutile"""

# constantes simulation --------------------------------------------------->
step = 0.01  # pas de calcul
end = 10  # temps de calcul

# forces (pour calculer les énergies)
f_p = rho*g*L**2*hc  # force pousée
f_G = m_s*g  # gravité appliquée sur centre de gravité
f_charge = m_c*g  # gravité appliquée sur G3

# Variables de simulation ---------------------------------------------------------------->

# Variables une composante

# Cinématique
t = np.arange(0, end, step)  # temps
angl = np.empty_like(t)  # angle d'inclinaison en fonction du temps
v_angl = np.empty_like(t)  # Vitesse angulaire
a_angl = np.empty_like(t)  # Acceleration angulaire
hl = np.empty_like(t)  # hauteur d'eau sur la gauche
hr = np.empty_like(t)  # hauteur d'eau sur la droite

# couples
c_r = np.empty_like(t)  # couple redressement gravité poussée
c_a = np.empty_like(t)  # Couple appliqué par la charge
c_d = np.empty_like(t)  # Couple d'amortidssenment(v_angl)

# centres de force
G = np.zeros((len(t), 2))  # Remplir centre de gravité de coordonnées nulles
P = np.zeros_like(G)  # Remplir centre de poussée de coordonnées nulles
G_c = np.zeros_like(G)  # remplir centre de gravité de la charge de coordonnées nulles


# correspondance au repère 2 apd ici ---------------------------------------->

# Vecteur de translation pour passer du repère 1 au repère 2 (voir shchéma notes)
v1_2 = [-G_0[0], 0]  # vecteur 1 vers 2
v2_1 = [G_0[0], 0]  # vecteur 2 vers 1

# Déplacement des points initiaux vers le repère 2
G_0 = r1_2_0(G_0)  # G_0 exprimé dans R1 -> G_0 exprimé dans R2
G_c_0 = r1_2_0(G_c_0)  # G_c_0 exprimé dans R1 -> G_c_0 exprimé dans R2
"""P_0 = g_trapeze(angl_0, L, hc, v1_2[0])"""


# simulation ---------------------------------------------------------------->
def simulation():
    dt = step
    # initialisation des variables cinématiques
    angl[0] = angl_0
    v_angl[0] = v_angl_0
    # Boucle de simulation (in range(0, end, step))
    for i in range(len(t)-1):
        # Points au temps i
        G[i] = rotation(G_0, [0, 0], angl[i])  # centre de gravité
        G_c[i] = rotation(G_c_0, [0, 0], angl[i])  # centre de gravité de la charge
        P[i] = g_trapeze(angl[i], L, hc, v1_2[0])  # a vérifier car la manip est complquée

        # Couples au temps i
        c_a[i] = m_c * g * G_c[i][0]  # couple destabilisteur en fonction de G_c_Y
        c_r[i] = f_p*(P[i][0]-G[i][0])

        # vitesse au temps i+1 avec la méthode d'Euler depuis F = Ia
        v_angl[i+1] = v_angl[i] + (-D*v_angl[i] + c_r[i] + c_a[i])*dt/I
        # vItesse au temps i+1 en fonction de la vitesse au temps i+1
        angl[i+1] = angl[i]+v_angl[i]*dt


# Representation graphique
def graphiques():
    plt.figure(1)
    plt.subplot(4, 1, 1)
    plt.plot(t, angl*180/pi, label="theta")
    plt.legend()
    plt.subplot(4, 1, 2)
    plt.plot(t, v_angl*180/pi, label="omega")
    plt.legend()
    plt.subplot(4, 1, 3)
    plt.plot(t, G, label="c_r")
    plt.legend()
    plt.subplot(4, 1, 4)
    plt.plot(t, P, label="c_a")
    plt.legend()
    plt.show()

simulation()
graphiques()








