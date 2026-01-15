import numpy as np
 import matplotlib . pyplot as plt
 def nrtl_binary (x1 , T, a12 , a21 , alpha ):
 """
 Calcul des coefficients d’activite avec NRTL binaire7
 Parametres :
 -----------
 x1 : float ou array
 Fraction molaire du constituant 1
   T : float
Temperature (K)
 a12 , a21 : float
 Parametres d’energie (K)
 alpha : float
 Parametre de non - randomness

 Retourne :
 ---------
 gamma1 , gamma2 : float ou array
 Coefficients d’activite
 """
 R = 8.314 # J/( mol .K)
 x2 = 1 - x1

 # Calcul des tau
 tau12 = a12 / T
tau21 = a21 / T

 # Calcul des G
 G12 = np.exp(- alpha * tau12 )
 G21 = np.exp(- alpha * tau21 )

 # Coefficients d’activite
term1_gamma1 = tau21 * (G21 / (x1 + x2* G21)) **2 term2_gamma1 = tau12 * G12 / (x2 + x1*G12 )**2 ln_gamma1 = x2 **2 * ( term1_gamma1 + term2_gamma1 )
term1_gamma2 = tau12 * (G12 / (x2 + x1* G12)) **2
term2_gamma2 = tau21 * G21 / (x1 + x2*G21 )**2
ln_gamma2 = x1 **2 * ( term1_gamma2 + term2_gamma2 )

 gamma1 = np. exp( ln_gamma1 )
 gamma2 = np. exp( ln_gamma2 )

 return gamma1 , gamma2

 # Parametres pour ethanol (1) -eau (2)
 #  T = 298.15 # K (25 C )
 a12 = -0.8009
 a21 = 1239.5
 alpha = 0.3

 # Calcul sur une plage de compositions
 x1_range = np. linspace (0.001 , 0.999 , 100)
 gamma1_vals , gamma2_vals = nrtl_binary ( x1_range , T, a12 , a21 ,
alpha )

 # Affichage des resultats
 print (" Resultats NRTL pour Ethanol (1) -Eau (2) a T = 25 C ")
 print ("=" *60)
 print (f"{’ x_ethanol ’: <12} {’ gamma_ethanol ’: <15} {’ gamma_eau
’: <15}")
 print ("-" *60)
 for i in [0, 24, 49, 74, 99]:
 print (f"{ x1_range [i ]: <12.3 f} { gamma1_vals [i ]: <15.4 f} {
gamma2_vals [i ]: <15.4 f}")

 # Coefficients a dilution infinie
 gamma1_inf = nrtl_binary (0.001 , T, a12 , a21 , alpha ) [0]
 gamma2_inf = nrtl_binary (0.999 , T, a12 , a21 , alpha ) [1]
 print ("\ nCoefficients a dilution infinie :")
 print (f" gamma_1 ^inf ( ethanol dans eau) = { gamma1_inf :.4 f}")
 print (f" gamma_2 ^inf (eau dans ethanol ) = { gamma2_inf :.4 f}")

 # Graphique
 plt. figure ( figsize =(10 , 6))
 plt. plot ( x1_range , gamma1_vals , ’b-’, linewidth =2, label =’
( thanol )’)
 plt. plot ( x1_range , gamma2_vals , ’r-’, linewidth =2, label =’
(eau )’)
78 plt. axhline (y=1, color =’k’, linestyle =’--’, alpha =0.3)
79 plt. xlabel (’Fraction molaire thanol , x ’, fontsize =12)
80 plt. ylabel (’Coefficient d\’ a c t i v i t , ’, fontsize =12)
81 plt. title (’ Coefficients d\’ a c t i v i t - S y s t m e thanol -Eau (
NRTL , 25 C )’,
82 fontsize =13 , fontweight =’bold ’)
83 plt. legend ( fontsize =11)
84 plt. grid (True , alpha =0.3)
85 plt. xlim (0, 1)
86 plt. ylim (0, max( gamma1_inf , gamma2_inf ) * 1.1)
8
plt. tight_layout ()
88 plt. savefig (’ nrtl_ethanol_eau .png ’, dpi =300 , bbox_inches =’tight ’
)
89 plt. show ()