### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 21305de4-1efd-11ec-071a-fd5d4d81dd65
md"""
!!! warning "Remarque importante"

	Les questions de ce *notebook* doivent être traitées de manière séquentielle : Q1, Q2... jusqu'à Q9.

"""

# ╔═╡ 50aedf4c-904a-4ef3-a113-6f21c6fc2d33
md"""
# Programmation

Quelques points vont être abordés en début de séance :

1. Deux nouveaux *packages* : `Plots` et `NLsolve` ;
1. Fonctions et types d'arguments ;
1. La notation `.` (*broadcast) ;
1. Premier (`first`) et dernier (`last`) éléments d'un tableau.

"""

# ╔═╡ 4c2ed0bd-ea08-495c-94c5-2c7e92ff29d6
md"""
# Recherche de la racine d'une fonction

Nous avons vu lors de la deuxième séance que les éléments
```math
\begin{aligned}
y_1 & \simeq y \left ( t_1 \right ), \\
y_2 & \simeq y \left ( t_2 \right ), \\
& \ldots
\end{aligned}
```
de la solution numérique du problème de Cauchy
```math
\left \{ \begin{aligned}
\dot{y} \left ( t \right ) & = f \left [ t, y \left ( t \right ) \right ], \\
y \left ( 0 \right ) & = y_0
\end{aligned} \right .
```
sont définis implicitement, c'est à dire comme racines de fonctions.

L'objectif de cette première partie est de se familiariser avec le *package* `NLsolve.jl` que nous utiliserons afin de résoudre des équations non-linéaires.

Les cellules suivantes décrivent comment obtenir la racine de la function
```math
\left ( \begin{matrix}
x_1 \\
x_2 \end{matrix} \right ) \mapsto \left ( \begin{matrix}
\left ( x_1 + 3 \right ) \left ( x_2 ^ 3 - 7 \right ) + 18 \\
\sin \left [ x_2 \exp \left ( x_1 \right ) - 1 \right ]
\end{matrix} \right )
```
à partir de la donnée initiale
```math
\left ( x_1, x_2 \right ) = \left ( 0.1, 1.2 \right ).
```

"""

# ╔═╡ 688b2980-43d3-4df1-825a-2da885beb08f
md"""
1. **Cas scalaire** -- Modifier l'exemple précédent afin de résoudre l'équation de Kepler
```math
10 - x + e \sin \left ( x \right ) = 0.
```
On définira dans un premier temps la fonction `kepler!` définie ci-dessous.

"""

# ╔═╡ 586ecaf0-6a36-468f-998c-d22cbb869990
md"""
# Modèle et solution exacte

On se concentre pour l'instant sur le modèle linéaire homogène pour lequel le second membre de l'EDO s'écrit
```math
f \colon \left ( t, y \right ) \mapsto \lambda y.
```

L'équation différentielle à résoudre s'écrit alors,
```math
\left \{ \begin{aligned}
\dot{y} \left ( t \right ) & = \lambda y \left ( t \right ), \\
y \left ( 0 \right ) & = y_0
\end{aligned} \right .
```
et la solution exacte est donnée sous la forme :
```math
y \colon t \mapsto \exp \left ( \lambda t \right ) y_0.
```

3. Implémenter la fonction ``f``, appelée ci-dessous `linear`, dans le cas ``\lambda = -1``.

"""

# ╔═╡ 344724c5-f849-4234-aa21-2c6e48a042c0
md"""
4. Implémenter la fonction `solution` qui correspond à la solution analytique dans le cas ``\lambda = -1`` et ``y_0 = 1``.

"""

# ╔═╡ f0559c7f-83ed-4797-8fe1-d87ed27419fa
md"""
# Schéma numérique

On rappelle que lors du cours précédent, tois schémas numériques ont été présentés, à savoir :
```math
y_{n + 1} - y_n - \tau f \left ( t_n, y_n \right ) = 0 \quad \text{(Euler explicite)},
```
```math
y_{n + 1} - y_n - \tau f \left ( t_{n + 1}, y_{n + 1} \right ) = 0 \quad \text{(Euler implicite)},
```
et
```math
y_{n + 1} - y_n - \tau f \left ( \frac{t_n + t_{n + 1}}{2}, \frac{y_n + y_{n + 1}}{2} \right ) = 0 \quad \text{(point milieu)}.
```

5. En vous inspirant de l'implémentation `explicit!` du schéma explicite d'Euler présentée ci-dessous, implémenter les fonctions `implicit!` (a) et `midpoint!` (b) dont la racine est ``y_{n + 1}``. On préservera le nombre et l'ordre des paramètres, au nombre de 6, à savoir

* `res` -- la valeur de la fonction implicite ;
* `x` -- la solution mise à jour (``y _ {n + 1}``) ;
* `y` -- la solution précédente (``y _ n``) ;
* `τ` -- le pas de temps (``\tau``) ;
* `f` -- le modèle (``f``) ;
* `t` -- l'instant précédent, (``t _ n``).

Pour le schéma d'Euler explicite, la solution mise à jour ``y_{n + 1}`` est donc définie comme la racine de la fonction implicite suivante
```math
F \left ( x, y, \tau, f, t \right ) = x - y - \tau f \left ( t, y \right ),
```
implémentée à l'aide de la fonction `explicit!` ci-dessous.

"""

# ╔═╡ 293407bc-2958-407f-b82b-d28f94023b60
md"""
!!! note "De l'usage de `!` en Julia"

	Par convention, `!` (*bang* en anglais) est ajouté à la fin du nom d'une fonction lorsque celle-ci modifie son premier argument (ici, `res`).

"""

# ╔═╡ c00f201c-b2b9-4226-ae7f-82269a399199
md"""
Modifier la fonction `implicit!` ci-dessous pour qu'elle corresponde au schéma d'Euler implicite
```math
F \left ( x, y, \tau, f, t \right ) = x - y - \tau f \left ( t + \tau, x \right ).
```

"""

# ╔═╡ 8825a1ac-5c1c-47fb-b16d-aade27463ad3
md"""
De même, modifier la fonction `midpoint!` ci-dessous pour lui faire correspondre le schéma du point milieu, qui s'écrira
```math
F \left ( x, y, \tau, f, t \right ) = x - y - \tau f \left ( t + \frac{\tau}{2}, \frac{x + y}{2} \right ).
```

"""

# ╔═╡ 690dbae8-a720-4d51-87a2-7834d84ee60b
md"""
# Intégration temporelle

Il reste à présent à assembler à implémenter la bouche d'intégration temporelle. Étant donnés

* Un modèle `f` ;
* Un pas de temps `τ` ;
* Et un instant `s`

la fonction `cauchy` implémentée ci-dessous retourne deux vecteurs, le premier contenant les instants
```math
t_0 \quad t_1 \quad \cdots \quad t_N = s
```
et le second la solution numérique, à savoir
```math
y_0 \quad y_1 \quad \cdots \quad y_N.
```

"""

# ╔═╡ 99b7cb68-508a-40ad-9007-4fbe9eb676ce
md"""
La solution numérique peut être obtenue et visualisée comme suit.

"""

# ╔═╡ a89a5fee-a28d-43df-be08-a146c476a5db
md"""
6. Utiliser les fonctions `linear` et `solution` définie précédemment dans l'implémentation de la fonction `error` ci-dessous, qui calcule l'erreur
```math
y_N - y \left ( t_N \right )
```
en fonction du schéma (`scheme!`) et du pas en temps (`τ`).

"""

# ╔═╡ b796ff78-0741-42b3-8285-b01a042feabc
md"""
7. Calculer (en utilisant la fonction `error`) et reporter les erreurs à l'instant `s = 1.0` dans le tableau ci-dessous. Commenter.

|             | `explicit!` | `implicit!` | `midpoint!` |
|:-----------:|:-----------:|:-----------:|:-----------:|
| `τ = 0.125` |             |             |             |
| `τ = 0.25`  |             |             |             |
| `τ = 0.5`   |             |             |             |
| `τ = 1.0`   |             |             |             |

8. On se place maintenant sur un horizon temporel plus long (`s = 10.0`). Augmenter la taille du pas de temps et commenter.

"""

# ╔═╡ a8242689-1552-4261-b1b8-767fa80b74ea
md"""

Tout l'intérêt de l'utilisation du package `NLsolve.jl` est que notre implémentation fonctionne pour les problèmes scalaires **non-linéaires**, ainsi que pour les cas **vectoriels**.

# Au delà du cas linéaire

9. Utiliser ou imaginer un modèle scalaire non-linéaire en modifier la fonction `nonlinear` ci-dessous, et visualiser votre solution numérique pour chacun des trois schémas sur le même graphique. Vous pourrez par exemple utiliser la question 3 de l'exercice vu en TD :
```math
f \colon \left ( t, y \right ) \mapsto 2t - y ^ 2.
```

!!! note "De l'usage du point"

	En Julia, le point (`.`) permet d'appliquer une fonction à chaque élément d'un tableau. Par exemple, la commande suivante élève chaque élément du tableau `y` au carré :
	```julia
	y .^ 2
	```

	Dans le doute, on peut aussi utiliser la *macro* `@.` comme suit :
	```julia
	@. y ^ 2
	```
    En un sens, elle "saupoudre" l'expression qui la suit de points.

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-beta4"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─21305de4-1efd-11ec-071a-fd5d4d81dd65
# ╟─50aedf4c-904a-4ef3-a113-6f21c6fc2d33
# ╟─4c2ed0bd-ea08-495c-94c5-2c7e92ff29d6
# ╟─688b2980-43d3-4df1-825a-2da885beb08f
# ╟─586ecaf0-6a36-468f-998c-d22cbb869990
# ╟─344724c5-f849-4234-aa21-2c6e48a042c0
# ╟─f0559c7f-83ed-4797-8fe1-d87ed27419fa
# ╟─293407bc-2958-407f-b82b-d28f94023b60
# ╟─c00f201c-b2b9-4226-ae7f-82269a399199
# ╟─8825a1ac-5c1c-47fb-b16d-aade27463ad3
# ╟─690dbae8-a720-4d51-87a2-7834d84ee60b
# ╟─99b7cb68-508a-40ad-9007-4fbe9eb676ce
# ╟─a89a5fee-a28d-43df-be08-a146c476a5db
# ╟─b796ff78-0741-42b3-8285-b01a042feabc
# ╟─a8242689-1552-4261-b1b8-767fa80b74ea
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
