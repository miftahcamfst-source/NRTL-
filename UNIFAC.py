import numpy as np
import matplotlib.pyplot as plt


class MoleculeUNIFAC:
    """
    Classe representant une molecule selon la methode UNIFAC.
    Chaque molecule est decomposee en groupes fonctionnels
    auxquels sont associes des parametres R (volume) et Q (surface).
    """

    def _init_(self, name, groups):
        """
        Initialisation de la molecule.

        Parametres
        ----------
        name : str
            Nom de la molecule
        groups : dict
            Dictionnaire {nom_groupe: nombre d'occurrences}
        """
        self.name = name
        self.groups = groups

        # Parametres de volume R des groupes UNIFAC
        self.R = {
            'CH3': 0.9011,
            'CH2': 0.6744,
            'CH': 0.4469,
            'OH': 1.0000,
            'H2O': 0.9200,
            'CH3CO': 1.6724
        }

        # Parametres de surface Q des groupes UNIFAC
        self.Q = {
            'CH3': 0.848,
            'CH2': 0.540,
            'CH': 0.228,
            'OH': 1.200,
            'H2O': 1.400,
            'CH3CO': 1.488
        }

        # Calcul du parametre r (volume moleculaire total)
        self.r = sum(n * self.R[grp] for grp, n in groups.items())

        # Calcul du parametre q (surface moleculaire totale)
        self.q = sum(n * self.Q[grp] for grp, n in groups.items())

    def _repr_(self):
        """
        Representation textuelle de la molecule
        (utile pour l'affichage des parametres calcules)
        """
        return (f"Molecule({self.name}): r={self.r:.4f}, q={self.q:.4f}, "
                f"groupes={self.groups}")


def unifac_combinatorial(x1, mol1, mol2):
    """
    Calcul de la contribution combinatorielle du modele UNIFAC.

    Cette partie prend en compte :
    - la taille des molecules
    - leur forme
    - les effets entropiques de melange

    Parametres
    ----------
    x1 : float ou array
        Fraction molaire du constituant 1
    mol1, mol2 : MoleculeUNIFAC
        Objets molecules

    Retour
    ------
    ln_gamma1_comb, ln_gamma2_comb : float ou array
        Logarithmes des coefficients d'activite (partie combinatorielle)
    """

    # Fraction molaire du constituant 2
    x2 = 1 - x1

    # Nombre de coordination (constant en UNIFAC)
    z = 10

    # Parametres moleculaires
    r1, r2 = mol1.r, mol2.r
    q1, q2 = mol1.q, mol2.q

    # Fractions de volume (occupation de l'espace)
    phi1 = x1 * r1 / (x1 * r1 + x2 * r2)
    phi2 = x2 * r2 / (x1 * r1 + x2 * r2)

    # Fractions de surface (interactions de surface)
    theta1 = x1 * q1 / (x1 * q1 + x2 * q2)
    theta2 = x2 * q2 / (x1 * q1 + x2 * q2)

    # Parametres l (correction de forme)
    l1 = (z / 2) * (r1 - q1) - (r1 - 1)
    l2 = (z / 2) * (r2 - q2) - (r2 - 1)

    # Expression UNIFAC de la partie combinatorielle
    ln_gamma1_comb = (
        np.log(phi1 / x1)
        + (z / 2) * q1 * np.log(theta1 / phi1)
        + l1
        - (phi1 / x1) * (x1 * l1 + x2 * l2)
    )

    ln_gamma2_comb = (
        np.log(phi2 / x2)
        + (z / 2) * q2 * np.log(theta2 / phi2)
        + l2
        - (phi2 / x2) * (x1 * l1 + x2 * l2)
    )

    return ln_gamma1_comb, ln_gamma2_comb


# =======================
# Definition des molecules
# =======================

# Ethanol : CH3–CH2–OH
ethanol = MoleculeUNIFAC('Ethanol', {'CH3': 1, 'CH2': 1, 'OH': 1})

# Eau : H2O
eau = MoleculeUNIFAC('Eau', {'H2O': 1})

# Acetone : CH3–CO–CH3
acetone = MoleculeUNIFAC('Acetone', {'CH3': 2, 'CH3CO': 1})


# Affichage des parametres moleculaires calcules
print("Parametres moleculaires UNIFAC:")
print("=" * 50)
print(ethanol)
print(eau)
print(acetone)
print()


# ==========================================
# Calcul UNIFAC (partie combinatorielle seule)
# pour le melange Ethanol – Eau
# ==========================================

# Domaine de composition
x1_range = np.linspace(0.01, 0.99, 99)

# Calcul des logarithmes des coefficients d'activite
ln_gamma1_comb, ln_gamma2_comb = unifac_combinatorial(
    x1_range, ethanol, eau
)

# Passage aux coefficients d'activite
gamma1_comb = np.exp(ln_gamma1_comb)
gamma2_comb = np.exp(ln_gamma2_comb)


# =======================
# Affichage numerique
# =======================

print("Contribution combinatorielle UNIFAC – Ethanol(1) / Eau(2):")
print("=" * 60)
print(f"{'x_ethanol':<12} {'gamma1_comb':<15} {'gamma2_comb':<15}")
print("-" * 60)

for i in [0, 24, 49, 74, 98]:
    print(f"{x1_range[i]:<12.3f} {gamma1_comb[i]:<15.4f} {gamma2_comb[i]:<15.4f}")


# =======================
# Representation graphique
# =======================

plt.figure(figsize=(10, 6))

# Coefficients d'activite combinatoriels
plt.plot(x1_range, gamma1_comb, 'b-', linewidth=2,
         label='Ethanol (combinatoriel)')
plt.plot(x1_range, gamma2_comb, 'r-', linewidth=2,
         label='Eau (combinatoriel)')

# Reference ideale (gamma = 1)
plt.axhline(y=1, color='k', linestyle='--', alpha=0.3)

plt.xlabel('Fraction molaire ethanol, x', fontsize=12)
plt.ylabel("Coefficient d'activité (partie combinatorielle)", fontsize=12)
plt.title('UNIFAC – Contribution combinatorielle Ethanol–Eau',
          fontsize=13, fontweight='bold')

plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.xlim(0, 1)
plt.tight_layout()
plt.savefig('unifac_combinatorial.png', dpi=300, bbox_inches='tight')
plt.show()


# =======================
# Remarque importante
# =======================

print("\nNote : Le calcul complet UNIFAC necessite la partie residuelle,")
print("qui depend des parametres d'interaction entre groupes")
print("(tables UNIFAC / DECHEMA).")