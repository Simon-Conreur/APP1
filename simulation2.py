"""
Groupe 11.14 LEPL 1501-APP1
Programme de simulation de la barge flottante simplifiée
TODO
Recalculer l'inertie
Rendre le déplacement de g_c plus intuitif
Diagramme de phase
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
    p -= c      # les points tournent autour de (0;0) donc on peut négliger cette ligne
    y_p = p[0]  # y_p = composante y du point qui tourne
    z_p = p[1]  # z_p = composante z du point qui tourne
    # matrice de rotation
    p[0] = cos(angle)*y_p - sin(angle)*z_p  # matrice de rotation ( y*cos(angle) - z*sin(angle) )
    p[1] = sin(angle)*y_p + cos(angle)*z_p  # ------------------- ( y*sin(angle) + z*cos(angle) )

    p += c
    return [p[0], p[1]]  # retourne le centre de rotation


def g_trapeze(theta, l, hc):
    """
    Calcul du centre de gravité d'un trapèze de bases hl et hr, le plus grand côté non parallèle vaut l.
    Theta est l'inclinaison qui engendre le trapèze
    Calculs détaillés dans les slides du CM de projet S8
    :param theta: Angle de rotation
    :param l: largeur de la barge
    :param hc: hauteur de l'eau au niveau du centre de gravité (hc)
    :param v: composante y du vecteur R1->R2 (= v1_2)
    :return: le centre de poussée P du trapèze de bases hl, hr et de hauteur l
    """
    hl = hc + tan(theta)*l/2
    hr = hc - tan(theta)*l/2
    lc_x = l * (hl + 2 * hr) / (3 * (hl + hr))
    lc_y = (hl ** 2 + hl * hr + hr ** 2) / (3 * (hl + hr))
    Yc = -l/2 + lc_x
    Zc = -hc + lc_y
    C = rotation([Yc, Zc], [0, 0], -theta)
    return C  # le centre de poussée


# constantes en SI ---------------------------------------------------------------->
#  constantes physiques
g = -9.81   # accélération due à la gravité de la Terre
rho = 1000  # masse volumique de l'eau
D = 1       # constante d'ammortissemnt

# constantes dimensions
L = 0.60   # Largeur de la barge
h1 = 0.03  # hauteur de la barge

#  masses
m1 = 3                       # masse barge
m_charge = 0.2               # Masse le la charge portée
m_bras_grap = 0              # Masse bras 2 et grappin
m3 = m_charge + m_bras_grap  # masse grappin + 2 eme bras + charge

# masses relatives
# (ces données sont à privilégier lors des calculs
# car si on change la configuration de la grue il suffit de changer les variables ici)
m_tot = m1 + m3  # masse totale de la structure
m_c = m3         # masse charge totale
m_s = m1         # masse de la stucture

# constantes calculées ---------------------------------------------------------------->

# hauteur flottaison moyenne au niveau du centre de gravité G
hc = m_tot/(rho*L**2)

# Inertie de la structure
I = (m1/12)*(h1**2 + L**2) + m1*(hc - h1/2)**2
I = 1.18
# Constantes valeurs initiales

# centres de gravité initiaux
G1_0 = np.array([0, -hc+h1/2])  # centre de gravité de la barge
G3_0 = np.array([0, 0.50])      # centre de gravité du bras 2, de la charge utile et du grappin

# Variables cinématiques initiales
angl_0 = 0*pi/180  # angle initial
v_angl_0 = 0       # vitesse angulaire initiale

# Centres relatifs initiaux
G_0 = [(G1_0[0]*m1)/m_s, (G1_0[1]*m1)/m_s]  # centre de gravité de la structure initiale
G_c_0 = G3_0                                # centre de gravité de la charge

# constantes simulation --------------------------------------------------->
step = 0.001  # pas de calcul
end = 10  # temps de calcul

# forces (pour calculer les énergies)
f_p = rho*(-g)*L**2*hc  # force pousée 0.25
f_G = m_s*g             # gravité appliquée sur centre de gravité
f_charge = m_c*g        # gravité appliquée sur G3

# Variables de simulation ---------------------------------------------------------------->

# Variables une composante

# Variables cinématiques
t = np.arange(0, end, step)  # temps
angl = np.empty_like(t)      # angle d'inclinaison en fonction du temps
v_angl = np.empty_like(t)    # Vitesse angulaire
a_angl = np.empty_like(t)    # Acceleration angulaire

# couples
c_r = np.empty_like(t)  # couple redressement gravité poussée
c_a = np.empty_like(t)  # Couple appliqué par la charge
c_d = np.empty_like(t)  # Couple d'amortidssenment(v_angl)


# Energies
Ek = np.empty_like(t)   # Energie cinétique
E_a = np.empty_like(t)  # Energie potentielle de la charge
E_P = np.empty_like(t)  # Energie potentielle due à la poussée d'archimède
E_G = np.empty_like(t)  # Energie potentielle due à la gravité
Em = np.empty_like(t)   # Energie mecanique

# centres de force
G = np.zeros((len(t), 2))  # Remplir centre de gravité de coordonnées nulles
P = np.zeros_like(G)       # Remplir centre de poussée de coordonnées nulles
G_c = np.zeros_like(G)     # remplir centre de gravité de la charge de coordonnées nulles


# déplacement de la charge
def G_c_deplacement(t, drop=False):
    time_drop = 4000  # nombre d'itération avant le laché du poids
    """
    Fonction qui gère les déplacements de la charge et son laché
    :param t: donne le nombre ed'itérations déjà réalisées
    :param drop: lache le poids si on le paramètre sur True
    :return: le vecteur de translation de la charge
    """
    global f_charge             # Permet de modifier la variable f_charge depuis la fonction
    if drop or t == time_drop:  # t = nombre d'itérations avant le lacher du poids
        f_charge = 0            # annule la force de la charge = lache la charge
    return np.array([t*0.0002, 0])     # une  fonction de x(t) y(t)


# simulation ---------------------------------------------------------------->
def simulation():
    dt = step  # delta t = le pas

    # initialisation des variables cinématiques
    angl[0] = angl_0      # angle initial
    v_angl[0] = v_angl_0  # vitesse angulaire initiale

    # Boucle de simulation (in range(0, end, step))
    for i in range(len(t)-1):
        # Points au temps i
        G[i] = rotation(G_0, [0, 0], angl[i])                           # centre de gravité de la barge
        G_c[i] = rotation(G_c_0 + G_c_deplacement(i), [0, 0], angl[i])  # centre de gravité de la charge
        P[i] = g_trapeze(angl[i], L, hc)                                   # Centre de poussée = centre de masse d'un trapèze


        # Couples au temps i
        c_a[i] = f_charge * G_c[i][0]   # couple destabilisteur appliqué par la charge
        c_r[i] = f_p*(P[i][0]-G[i][0])  # couple de redressement appliqué par la gravité et la poussée d'archimède

        # Vitesse au temps i+1 : Methode d'Euler qui provient de F = Ia
        v_angl[i+1] = v_angl[i] + (-D*v_angl[i] + c_r[i] + c_a[i])*dt/I

        # Angle au temps i+1 en fonction de la vitesse au temps i+1 : Methode d'Euler
        angl[i+1] = angl[i]+v_angl[i]*dt

        # Energies au temps i
        Ek[i] = I * m_tot * v_angl[i] ** 2 / 2      # Energie cinétique
        E_a[i] = -c_a[i]*angl[i]                    # Energie potentielle de la charge
        E_G[i] = -m_tot*g*(G[i][1] - G[0][1])       # Energie potentielle de gravité
        E_P[i] = -m_tot*g*(P[i][1] - P[0][1])       # Energie potentialle de flottason
        Em[i] = Ek[i] + E_a[i] + E_G[i] + E_P[i]    # Energie mecanique


# Representation graphique
def graphiques():
    # Angle en fonction du temps
    angl_submersion = atan((h1 - hc)*2/L)*180/pi  # angle maximal avant submersion
    angl_soulevement = atan(-hc*2/L)*180/pi       # angle maximal avant soulevement

    plt.figure(1)               # Fenetre 1
    plt.subplot(5, 1, 1)  # Premier graphique
    plt.plot(t, angl*180/pi, label="theta")
    plt.plot([0, end], [angl_submersion, angl_submersion], '--r', label="Submersion")     # Submersion : borne +
    plt.plot([0, end], [-angl_submersion, -angl_submersion], '--r')                       # Submersion : borne -
    plt.plot([0, end], [angl_soulevement, angl_soulevement], '--g', label="Soulèvement")  # Soulevement : borne +
    plt.plot([0, end], [-angl_soulevement, -angl_soulevement], '--g')                     # Soulevement : borne -
    plt.xlabel("[ s ]")  # légende unité axe x
    plt.ylabel("[ ° ]")  # légende unité axe y
    plt.legend()         # Graduation des axes

    # Vitesse angulaire en fonction du temps
    plt.subplot(5, 1, 2)                       # deuxième graphique
    plt.plot(t, v_angl*180/pi, label="omega")  # Vitesse angulaire
    plt.xlabel("[ s ]")                              # légende unité axe x
    plt.ylabel("[ °/s ]")                            # légende unité axe y
    plt.legend()                                     # Graduation des axes
    plt.subplot(5, 1, 3)
    plt.plot(t, P, label="c_r")
    plt.legend()

    # autre
    plt.subplot(5, 1, 4)
    plt.plot(t, c_r, label="c_a")
    plt.xlabel("[ s ]")
    plt.legend()
    plt.subplot(5, 1, 5)
    plt.plot(t, c_a, label="c_r")
    plt.xlabel("[ s ]")
    plt.legend()
    plt.show()


def graphiques_energie():
    plt.figure(2)  # Deuxieme fenetre
    plt.subplot(1, 1, 1)  # Premier graphique

    plt.plot(t, Ek, label="Energie Cinétique")                # Energie cinétique
    plt.plot(t, E_G, label="Energie potentielle de gravité")  # Energie potentille de gravité
    plt.plot(t, E_P, label="Energie potentielle de poussée")  # Energie potentielle de poussée
    plt.plot(t, E_a, label="Energie potentielle de charge")   # Energie potentielle de gravité de la charge
    plt.plot(t, Em, label="Energie mécanique")                # Energie mécanique

    plt.ylabel("[ J ]")  # légende unité axe y
    plt.xlabel("[ s ]")  # légende unité axe x
    plt.legend()         # Graduation des axes
    plt.show()           # Affichage des graphiques

def diagramme_de_phase():
    plt.figure(3)
    plt.title("Diagramme de phase")
    plt.plot(angl*180/pi, v_angl*180/pi)
    plt.xlabel("theta (o)")
    plt.ylabel("omega (o/s)")
    plt.legend()
    plt.show()
    
    
    
simulation()
graphiques()          # Affichage des graphiques cinématiques
graphiques_energie()  # Affichages des graphiques énergétiques
diagramme_de_phase()  # Affichage du diagramme de phase
