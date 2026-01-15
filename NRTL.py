# =========================
# Import des bibliothèques
# =========================
import numpy as np              # Calcul numérique (tableaux, fonctions math)
import matplotlib.pyplot as plt # Tracés graphiques


# ======================================================
# Fonction NRTL pour un système binaire
# ======================================================
def nrtl_binary(x1, T, a12, a21, alpha):
    """
    Calcul des coefficients d'activité γ1 et γ2
    avec le modèle NRTL (Non-Random Two Liquids)
    pour un système binaire.

    Paramètres
    ----------
    x1 : float ou array
        Fraction molaire du constituant 1 (éthanol)
    T : float
        Température en Kelvin
    a12, a21 : float
        Paramètres d'interaction énergétique NRTL
    alpha : float
        Paramètre de non-randomness (0.2–0.47 en général)

    Retour
    ------
    gamma1, gamma2 : float ou array
        Coefficients d’activité du constituant 1 et 2
    """

    # Fraction molaire du constituant 2 (eau)
    x2 = 1 - x1

    # =========================
    # Calcul des paramètres τ
    # =========================
    # τij = aij / T
    # Ils représentent l'énergie d'interaction relative
    tau12 = a12 / T
    tau21 = a21 / T

    # =========================
    # Calcul des paramètres G
    # =========================
    # Gij = exp(-alpha * tauij)
    # alpha décrit le degré de non-randomness du mélange
    G12 = np.exp(-alpha * tau12)
    G21 = np.exp(-alpha * tau21)

    # ======================================================
    # Calcul de ln(gamma1) – coefficient d'activité éthanol
    # ======================================================
    term1_gamma1 = tau21 * (G21 / (x1 + x2 * G21))**2
    term2_gamma1 = tau12 * G12 / (x2 + x1 * G12)**2

    # Expression complète de ln(γ1)
    ln_gamma1 = x2**2 * (term1_gamma1 + term2_gamma1)

    # ======================================================
    # Calcul de ln(gamma2) – coefficient d'activité eau
    # ======================================================
    term1_gamma2 = tau12 * (G12 / (x2 + x1 * G12))**2
    term2_gamma2 = tau21 * G21 / (x1 + x2 * G21)**2

    # Expression complète de ln(γ2)
    ln_gamma2 = x1**2 * (term1_gamma2 + term2_gamma2)

    # Passage de ln(γ) à γ
    gamma1 = np.exp(ln_gamma1)
    gamma2 = np.exp(ln_gamma2)

    return gamma1, gamma2


# ======================================================
# Paramètres du système Éthanol (1) – Eau (2)
# ======================================================
T = 298.15        # Température (K) → 25 °C
a12 = -0.8009     # Paramètre NRTL éthanol → eau
a21 = 1239.5      # Paramètre NRTL eau → éthanol
alpha = 0.3       # Paramètre de non-randomness


# ======================================================
# Domaine de composition
# ======================================================
# Fraction molaire d’éthanol entre 0 et 1
x1_range = np.linspace(0.001, 0.999, 100)


# ======================================================
# Calcul des coefficients d’activité
# ======================================================
gamma1_vals, gamma2_vals = nrtl_binary(x1_range, T, a12, a21, alpha)


# ======================================================
# Affichage de quelques valeurs numériques
# ======================================================
print("Résultats NRTL pour Éthanol (1) - Eau (2) à T = 25 °C")
print("=" * 60)
print(f"{'x_ethanol':<12} {'gamma_ethanol':<18} {'gamma_eau':<18}")
print("-" * 60)

# Affichage de points représentatifs
for i in [0, 24, 49, 74, 99]:
    print(f"{x1_range[i]:<12.3f} {gamma1_vals[i]:<18.4f} {gamma2_vals[i]:<18.4f}")


# ======================================================
# Coefficients d’activité à dilution infinie
# ======================================================
# γ1∞ : éthanol très dilué dans l'eau
gamma1_inf = nrtl_binary(0.001, T, a12, a21, alpha)[0]

# γ2∞ : eau très diluée dans l’éthanol
gamma2_inf = nrtl_binary(0.999, T, a12, a21, alpha)[1]

print("\nCoefficients à dilution infinie :")
print(f"gamma1^inf (éthanol dans eau) = {gamma1_inf:.4f}")
print(f"gamma2^inf (eau dans éthanol) = {gamma2_inf:.4f}")


# ======================================================
# Tracé des coefficients d’activité
# ======================================================
plt.figure(figsize=(10, 6))

# Courbe γ éthanol
plt.plot(x1_range, gamma1_vals, 'b-', linewidth=2, label='γ éthanol')

# Courbe γ eau
plt.plot(x1_range, gamma2_vals, 'r-', linewidth=2, label='γ eau')

# Ligne γ = 1 (comportement idéal)
plt.axhline(y=1, linestyle='--', alpha=0.4)

# Mise en forme du graphique
plt.xlabel("Fraction molaire d'éthanol, x₁", fontsize=12)
plt.ylabel("Coefficient d'activité γ", fontsize=12)
plt.title("Coefficients d'activité – Système Éthanol / Eau (NRTL, 25 °C)",
          fontsize=13, fontweight='bold')

plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.xlim(0, 1)
plt.ylim(0, max(gamma1_inf, gamma2_inf) * 1.1)

# Sauvegarde et affichage
plt.tight_layout()
plt.savefig("nrtl_ethanol_eau.png", dpi=300, bbox_inches='tight')
plt.show()