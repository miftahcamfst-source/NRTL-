import numpy as np
import matplotlib.pyplot as plt

def uniquac_binary(x1, T, r1, r2, q1, q2, a12, a21):
    """
    Calcul des coefficients d’activité γ1 et γ2
    avec le modèle UNIQUAC pour un système binaire.

    UNIQUAC = Universal Quasi-Chemical
    Le modèle décompose ln(γ) en :
      - une partie combinatorielle (effets de taille et forme)
      - une partie résiduelle (interactions énergétiques)
    """

    # Constante des gaz parfaits
    R = 8.314  # J/(mol·K)

    # Nombre de coordination (hypothèse UNIQUAC)
    z = 10

    # Fraction molaire du composant 2
    x2 = 1 - x1

    # ======================================================
    # PARTIE COMBINATORIELLE
    # (effets entropiques : taille et forme des molécules)
    # ======================================================

    # Fractions de volume (occupation volumique)
    phi1 = x1 * r1 / (x1 * r1 + x2 * r2)
    phi2 = x2 * r2 / (x1 * r1 + x2 * r2)

    # Fractions de surface (contacts de surface)
    theta1 = x1 * q1 / (x1 * q1 + x2 * q2)
    theta2 = x2 * q2 / (x1 * q1 + x2 * q2)

    # Paramètres l (correction de coordination)
    l1 = (z / 2) * (r1 - q1) - (r1 - 1)
    l2 = (z / 2) * (r2 - q2) - (r2 - 1)

    # Contribution combinatorielle à ln(gamma1)
    ln_gamma1_comb = (
        np.log(phi1 / x1)                           # effet volume
        + (z / 2) * q1 * np.log(theta1 / phi1)      # effet surface
        + l1
        - (phi1 / x1) * (x1 * l1 + x2 * l2)
    )

    # Contribution combinatorielle à ln(gamma2)
    ln_gamma2_comb = (
        np.log(phi2 / x2)
        + (z / 2) * q2 * np.log(theta2 / phi2)
        + l2
        - (phi2 / x2) * (x1 * l1 + x2 * l2)
    )

    # ======================================================
    # PARTIE RÉSIDUELLE
    # (effets enthalpiques : interactions moléculaires)
    # ======================================================

    # Paramètres d’interaction réduits
    # tau_ij = exp(-a_ij / T)
    tau12 = np.exp(-a12 / T)
    tau21 = np.exp(-a21 / T)

    # Contribution résiduelle à ln(gamma1)
    ln_gamma1_res = q1 * (
        1
        - np.log(theta1 + theta2 * tau21)
        - theta1 / (theta1 + theta2 * tau21)
        - theta2 * tau12 / (theta2 + theta1 * tau12)
    )

    # Contribution résiduelle à ln(gamma2)
    ln_gamma2_res = q2 * (
        1
        - np.log(theta2 + theta1 * tau12)
        - theta2 / (theta2 + theta1 * tau12)
        - theta1 * tau21 / (theta1 + theta2 * tau21)
    )

    # ======================================================
    # COEFFICIENTS D’ACTIVITÉ TOTAUX
    # ======================================================

    ln_gamma1 = ln_gamma1_comb + ln_gamma1_res
    ln_gamma2 = ln_gamma2_comb + ln_gamma2_res

    gamma1 = np.exp(ln_gamma1)
    gamma2 = np.exp(ln_gamma2)

    return gamma1, gamma2


# ======================================================
# APPLICATION : Acétone (1) – Chloroforme (2)
# ======================================================

# Température (50 °C)
T = 323.15  # K

# Paramètres structurels UNIQUAC
# r : volume relatif
# q : surface relative
r1, q1 = 2.574, 2.336   # Acétone
r2, q2 = 2.870, 2.410   # Chloroforme

# Paramètres d’interaction énergétique
a12 = -52.39
a21 = 340.35

# Domaine de composition
x1_range = np.linspace(0.001, 0.999, 100)

# Calcul des coefficients d’activité
gamma1_vals, gamma2_vals = uniquac_binary(
    x1_range, T, r1, r2, q1, q2, a12, a21
)

# ======================================================
# AFFICHAGE NUMÉRIQUE
# ======================================================
print("Résultats UNIQUAC pour Acétone (1) – Chloroforme (2) à 50 °C")
print("=" * 60)
print(f"{'x_acetone':<12} {'gamma_acetone':<15} {'gamma_chloroforme':<15}")
print("-" * 60)

# Affichage de quelques points représentatifs
for i in [0, 24, 49, 74, 99]:
    print(f"{x1_range[i]:<12.3f} {gamma1_vals[i]:<15.4f} {gamma2_vals[i]:<15.4f}")

# ======================================================
# VÉRIFICATION DE GIBBS–DUHEM
# (condition thermodynamique de cohérence)
# ======================================================
from scipy.integrate import simpson

area1 = simpson(np.log(gamma1_vals) * x1_range, x1_range)
area2 = simpson(np.log(gamma2_vals) * (1 - x1_range), x1_range)

print("\nVérification Gibbs–Duhem :")
print(f"Aire ln(gamma1) : {area1:.6f}")
print(f"Aire ln(gamma2) : {area2:.6f}")
print(f"Différence : {abs(area1 - area2):.6f} (≈ 0 attendu)")

# ======================================================
# REPRÉSENTATIONS GRAPHIQUES
# ======================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Coefficients d’activité
ax1.plot(x1_range, gamma1_vals, 'b-', linewidth=2, label='γ acétone')
ax1.plot(x1_range, gamma2_vals, 'r-', linewidth=2, label='γ chloroforme')
ax1.axhline(y=1, linestyle='--', alpha=0.3)
ax1.set_xlabel("Fraction molaire acétone, x₁")
ax1.set_ylabel("Coefficient d’activité γ")
ax1.set_title("UNIQUAC – Coefficients d’activité", fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Énergie de Gibbs excédentaire (sans unité)
GE_RT = x1_range * (1 - x1_range) * (
    np.log(gamma1_vals) + np.log(gamma2_vals)
)

ax2.plot(x1_range, GE_RT, 'g-', linewidth=2)
ax2.axhline(y=0, linewidth=0.5)
ax2.set_xlabel("Fraction molaire acétone, x₁")
ax2.set_ylabel(r"$G^E / RT$")
ax2.set_title("Énergie de Gibbs excédentaire", fontweight='bold')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("uniquac_acetone_chloroforme.png", dpi=300, bbox_inches='tight')
plt.show()