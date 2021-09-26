### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 9e164dc1-55f8-4f75-b050-3f566f3cf328
using LinearAlgebra, NLsolve

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

# ╔═╡ d3588870-17f8-467b-928b-1b6781218296
function example!(res, x)
    res[1] = (x[1] + 3) * (x[2] ^ 3 - 7) + 18
    res[2] = sin(x[2] * exp(x[1]) - 1)
	nothing
end

# ╔═╡ 3fc5c2e0-d05f-400c-a250-6fb7207ff250
getproperty(nlsolve(example!, [0.1; 1.2]), :zero)

# ╔═╡ 688b2980-43d3-4df1-825a-2da885beb08f
md"""
1. **Cas scalaire** -- Modifier l'exemple précédent afin de résoudre l'équation de Kepler
```math
10 - x + e \sin \left ( x \right ) = 0.
```
On définira dans un premier temps la fonction `kepler!` définie ci-dessous.

"""

# ╔═╡ bd298666-6284-4b63-a217-1a0d46f2d7af
# Q1 -- À MODIFIER
function kepler!(res, x)
	res[1] = x[1]
	nothing
end

# ╔═╡ 8de162b4-2a13-48c1-b9e9-133035d842ce
if norm(nlsolve(kepler!, [0.0]).zero - [5.08912]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `kepler!` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `kepler!`.
	"""
end

# ╔═╡ e08458c6-78f3-4b34-b8ba-79dd83441f68
md"""
2. **Cas vectoriel** -- Modifier la fonction `system!` ci-dessous afin de résoudre le système d'équations
```math
\left \{ \begin{aligned}
x_1 + x_2 + x_3 ^ 2 & = 12, \\
x_1 ^ 2 - x_2 + x_3 & = 2, \\
2x_1 - x_2 ^ 2 + x_3 & = 1.
\end{aligned} \right .
```

"""

# ╔═╡ 74375a4a-8aea-43e3-b293-562ae8d5c114
# Q2 -- À MODIFIER
function system!(res, x)
	res[1] = x[1] + 1
	res[2] = x[2] + 2
	res[3] = x[3] + 3
	nothing
end

# ╔═╡ 2b611f9e-17c7-45d3-ba04-5563e9bb7ed3
if norm(nlsolve(system!, [0.0; 0.0; 0.0]).zero - [1.0; 2.0; 3.0]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `system!` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `system!`.
	"""
end

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

# ╔═╡ b147630e-e17e-48ca-8733-a0588f6e6d3b
# Q3 -- À MODIFIER
linear(t, y) = zero(y)

# ╔═╡ 293dea21-2128-41d3-b11c-47225d27d168
if norm(linear(nothing, [π; 2 // 3]) + [π; 2 // 3]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `linear` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `linear`.
	"""
end

# ╔═╡ 344724c5-f849-4234-aa21-2c6e48a042c0
md"""
4. Implémenter la fonction `solution` qui correspond à la solution analytique dans le cas ``\lambda = -1`` et ``y_0 = 1``.

"""

# ╔═╡ 4dafb136-cadb-4689-8511-01b72b308515
# Q4 -- À MODIFIER
solution(t, y = ones(1)) = zero(y)

# ╔═╡ 957c5d5f-df77-463e-8cc8-3affec13a636
if norm(solution(1.0) - [0.36788]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `solution` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `solution`.
	"""
end

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

# ╔═╡ 02ef58f1-235a-4108-917d-5c2e7e0527a1
function explicit!(res, x, y, τ, f, t)
	res .= x - y - τ * f(t, y)
end

# ╔═╡ c00f201c-b2b9-4226-ae7f-82269a399199
md"""
Modifier la fonction `implicit!` ci-dessous pour qu'elle corresponde au schéma d'Euler implicite
```math
F \left ( x, y, \tau, f, t \right ) = x - y - \tau f \left ( t + \tau, x \right ).
```

"""

# ╔═╡ fb414e6a-52ee-419c-b2dc-c320c27f17e2
# Q5a -- À MODIFIER
function implicit!(res, x, y, τ, f, t)
	res .= x - y - τ * f(t, y)
end

# ╔═╡ c0a78db4-794e-48f8-a34d-d7fcbddb279b
if norm(implicit!(zeros(1), ones(1), ones(1), √2, (t, y) -> t .+ y, π) - [-7.8571]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `implicit!` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `implicit!`.
	"""
end

# ╔═╡ 8825a1ac-5c1c-47fb-b16d-aade27463ad3
md"""
De même, modifier la fonction `midpoint!` ci-dessous pour lui faire correspondre le schéma du point milieu, qui s'écrira
```math
F \left ( x, y, \tau, f, t \right ) = x - y - \tau f \left ( t + \frac{\tau}{2}, \frac{x + y}{2} \right ).
```

"""

# ╔═╡ 1e4def72-64a4-4b2e-bf0f-86a4ee57bc10
# Q5b -- À MODIFIER
function midpoint!(res, x, y, τ, f, t)
	res .= x - y - τ * f(t, y)
end

# ╔═╡ a37aa68d-f4ff-4fb2-b8d2-ff80bd67b6a1
if norm(midpoint!(zeros(1), ones(1), ones(1), √2, (t, y) -> t .+ y, π) - [-6.8571]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `midpoint!` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `midpoint!`.
	"""
end

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

# ╔═╡ 14baffbc-2cf0-4b6c-8b35-bfaf0605575c
function cauchy(scheme!, f, τ, s, y₀, t₀ = zero(τ))
	t, y = t₀, y₀
    T, Y = [t], [y]

	while t < (1 - √eps(t)) * s
		y = getproperty(
			nlsolve(y) do res, x
				scheme!(res, x, y, τ, f, t)
			end,
			:zero
		)
        t += τ

        push!(Y, y)
        push!(T, t)
	end

	T, Y
end

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

# ╔═╡ c1c0ea70-fc6b-4928-adff-e58edd22d7e0
# Q6 -- À MODIFIER
function error(scheme!, τ, s)
	T, num = cauchy(scheme!, linear, τ, s, ones(1))
	exact = solution.(T)
	norm(last(num))
end

# ╔═╡ f5c664fc-9877-4800-b261-d9233d4c006a
if norm(error(explicit!, 0.2, 1.0) - 0.0401994) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `error` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `error`.
	"""
end

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

# ╔═╡ 1c422e7a-3e59-4b4a-af58-cfb1e5818c53
# Q9 -- À MODIFIER
nonlinear(t, y) = @. y

# ╔═╡ 778c181e-7253-4666-9075-14aa9f0315fa
if norm(nonlinear(π, [√2]) - [4.28319]) ≤ 1e-4
	md"""
	!!! tip "😃 Bonne réponse"

		Votre implémentation de `nonlinear` est correcte.
	"""
else
	md"""
	!!! danger "😡 Mauvaise réponse"

		Vérifier votre implémentation de `nonlinear`.
	"""
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"

[compat]
NLsolve = "~4.5.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-beta4"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "b8d49c34c3da35f220e7295659cd0bab8e739fed"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.33"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "bd4afa1fdeec0c8b89dad3c6e92bc6e3b0fec9ce"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.6.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "7220bc21c33e990c14f4a9a319b1d242ebc5b269"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "9f46deb4d4ee4494ffb5a40a27a2aced67bdd838"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.4"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "b5e930ac60b613ef3406da6d4f42c35d8dc51419"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.19"

[[deps.IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "34dc30f868e368f8a17b728a1238f3fcda43931a"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.3"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "144bab5b1443545bc4e791536c9f1eacb4eed06a"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.1"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ad42c30a6204c74d264692e633133dcea0e8b14e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.2"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "a8f30abc7c64a39d389680b74e749cf33f872a70"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─21305de4-1efd-11ec-071a-fd5d4d81dd65
# ╠═9e164dc1-55f8-4f75-b050-3f566f3cf328
# ╟─50aedf4c-904a-4ef3-a113-6f21c6fc2d33
# ╟─4c2ed0bd-ea08-495c-94c5-2c7e92ff29d6
# ╠═d3588870-17f8-467b-928b-1b6781218296
# ╠═3fc5c2e0-d05f-400c-a250-6fb7207ff250
# ╟─688b2980-43d3-4df1-825a-2da885beb08f
# ╠═bd298666-6284-4b63-a217-1a0d46f2d7af
# ╟─8de162b4-2a13-48c1-b9e9-133035d842ce
# ╟─e08458c6-78f3-4b34-b8ba-79dd83441f68
# ╠═74375a4a-8aea-43e3-b293-562ae8d5c114
# ╟─2b611f9e-17c7-45d3-ba04-5563e9bb7ed3
# ╟─586ecaf0-6a36-468f-998c-d22cbb869990
# ╠═b147630e-e17e-48ca-8733-a0588f6e6d3b
# ╟─293dea21-2128-41d3-b11c-47225d27d168
# ╟─344724c5-f849-4234-aa21-2c6e48a042c0
# ╠═4dafb136-cadb-4689-8511-01b72b308515
# ╟─957c5d5f-df77-463e-8cc8-3affec13a636
# ╟─f0559c7f-83ed-4797-8fe1-d87ed27419fa
# ╟─293407bc-2958-407f-b82b-d28f94023b60
# ╠═02ef58f1-235a-4108-917d-5c2e7e0527a1
# ╟─c00f201c-b2b9-4226-ae7f-82269a399199
# ╠═fb414e6a-52ee-419c-b2dc-c320c27f17e2
# ╟─c0a78db4-794e-48f8-a34d-d7fcbddb279b
# ╟─8825a1ac-5c1c-47fb-b16d-aade27463ad3
# ╠═1e4def72-64a4-4b2e-bf0f-86a4ee57bc10
# ╟─a37aa68d-f4ff-4fb2-b8d2-ff80bd67b6a1
# ╟─690dbae8-a720-4d51-87a2-7834d84ee60b
# ╠═14baffbc-2cf0-4b6c-8b35-bfaf0605575c
# ╟─99b7cb68-508a-40ad-9007-4fbe9eb676ce
# ╟─a89a5fee-a28d-43df-be08-a146c476a5db
# ╠═c1c0ea70-fc6b-4928-adff-e58edd22d7e0
# ╟─f5c664fc-9877-4800-b261-d9233d4c006a
# ╟─b796ff78-0741-42b3-8285-b01a042feabc
# ╟─a8242689-1552-4261-b1b8-767fa80b74ea
# ╠═1c422e7a-3e59-4b4a-af58-cfb1e5818c53
# ╟─778c181e-7253-4666-9075-14aa9f0315fa
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
