Cette page contient les codes Matlab utilisé pour obtenir les résultats d'optimisation pour ce Challenge.

### Programmes principaux
1. `Pareto.m` permet de tracer le front de Pareto en fonction du bruit et de la puissance consommée.
2. `Optimisation_corde_twist.m` permet de obtenir un seul résultat d'optimisation.

### Fonctions diverses
1. `contraint_aerody_corde_twist.m` calcule la poussée de l'hélice. Cette fonction est appelée comme contrainte de l'optimisation.
2. `contraint_noice_corde_twist.m` calcule le bruit de l'hélice. Cette fonction est appelée comme contrainte de l'optimisation.
3. `contraint_noice_corde_twist_figure.m` permet de tracer la répartition du bruit dans différentes directions.
4. `contraint_noice_direction.m` calcule le bruit dans différentes directions, donc elle est appelée dans `contraint_noice_corde_twist_figure.m.
5. `fun_aerody_corde_twist.m` calcule la puissance consommée par l'hélice. Elle est appelée comme la fonction à minimiser de l'optimisation.
6. `noice_ele_corde_twist.m` calcule le bruit généré par un segment de pale.
7. `Viterna_Corrigan.m` sert à faire l'extrapolation de $C_l$ et $C_d$.
