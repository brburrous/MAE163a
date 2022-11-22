### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 625a9d5c-62c0-11ed-3363-bf1b1417eda1
using SymPy, Latexify, PlutoUI, LaTeXStrings

# ╔═╡ 34827925-8a91-4ae9-ac10-2b9f412ac97c
@vars θ1 θ2 θ3 θ4 L1 L2 L3 L4 d2 d3

# ╔═╡ 3aa6b114-e38f-4f35-80ab-490d0fbe94de


# ╔═╡ 12bb4f91-d0cd-41ba-9d2f-284ab51eb458
begin 
	function Transform(α, a, d, θ)
		[
			cos(θ) 			-sin(θ) 		0 		a
			sin(θ)*cos(α) 	cos(θ)*cos(α) 	-sin(α)	-sin(α)*d
			sin(θ)*sin(α) 	cos(θ)*sin(α) 	cos(α) 	cos(α)*d
			0 				0 				0 		1
		]
	end
	Transform(V::Vector) = Transform(V...)
	
	function transforms(dh)
		[Transform(dh[i, :]) for i = 1:size(dh)[1]]
	end
	
	function inverse(hT)
		r = hT[1:3, 1:3]
		r_trans = [r[j,i] for i = 1:3, j = 1:3]
		p = hT[1:3, 4]
		p_trans = [-r_trans*p; 1]
		hcat(vcat(r_trans, [0 0 0]), p_trans)
	end
	
	function split(Ts::Vector, i, target)
		if i != 0
			invT = [inverse(Ts[j]) for j = i:-1:1]
			lh = reduce(*, invT)*target
			lh = simplify.(lh)
			rh = simplify.(reduce(*, Ts[i+1:end]))
		else
			rh = reduce(*, Ts)
			rh = simplify.(rh)
			lh = target
		end
		return [lh, rh]
	end
	
end

# ╔═╡ f3a37ec8-347b-4221-866f-688bc6781f16


# ╔═╡ 30d171e2-f450-4aff-bfa3-52b372236dc5
dhParams = [
	0 	0 	0 	θ1
	0 	L1 	d2 	θ2
	PI/2 	0 	d3 	θ3
	0 	L3 	0 	θ4
	0 	L4 	0 	0
]

# ╔═╡ 16e67336-e47b-4e56-bf6e-526c50c7854b
Ts = transforms(dhParams)

# ╔═╡ c2568846-c99e-4e5c-89e3-34ab422b764f
@vars r₁₁ r₁₂ r₁₃ r₂₁ r₂₂ r₂₃ r₃₁ r₃₂ r₃₃ p_x p_y p_z

# ╔═╡ 75f669cb-4313-4e35-a5b9-6c7fb28455bf
target = [
	r₁₁ r₁₂ r₁₃ p_x
	r₂₁ r₂₂ r₂₃ p_y
	r₃₁ r₃₂ r₃₃ p_z
	0 0 0 1
]

# ╔═╡ d938636d-f678-4210-993a-d5a16f0ff2fa
@bind i Slider(0:length(Ts)-1, show_value=true)

# ╔═╡ 8726790c-f378-49df-897c-10dd3f07d1d4
splitEq = split(Ts, i, target)

# ╔═╡ d3cddad4-c4c2-4e36-bcef-1f2b261e51dd
testEq = splitEq[1][1]

# ╔═╡ 27ab9d60-c4f1-4560-a316-088df657f55b
string(testEq)

# ╔═╡ 22430021-1347-4d11-a5e7-159915eb9b70
function getEq(Ts, i, target, n)
	splitEq = split(Ts, i, target)
	eqs = collect(zip(splitEq[1], splitEq[2]))
	eq = simplify.(eqs[n])
end

# ╔═╡ 4b1d01d1-546e-4121-8eab-3c7162054fe6
T01 = [
	 cos(θ1) -sin(θ1) 0 0
	sin(θ1) cos(θ1) 0 0
 	 0 		 0  	 1 0
	 0 		 0 		 0 1
]

# ╔═╡ 7a1f1531-f99c-4906-93b6-8b52701eb5fe
T12 = [
	cos(θ2) -sin(θ2) 0 0
	0 0 1 0 
	-sin(θ2) -cos(θ2) 0 0
	 0 		 0 		 0 1
]

# ╔═╡ 7f8a3a00-77c4-4761-88f4-933c51e2d6b1
@vars a2 a3 d4 θ5 θ6

# ╔═╡ 208e0373-199c-4240-8da0-79eb588f9f7b
T23 = [
	 cos(θ3) -sin(θ3) 0 a2
	sin(θ3) cos(θ3) 0 0
 	 0 		 0  	 1 d3
	 0 		 0 		 0 1
]


# ╔═╡ aba4e8d1-9c4f-4dbc-938a-9697ce19a6dc
T34 = [
	cos(θ4) -sin(θ4) 0 a3
	0 0 1 d4 
	-sin(θ4) -cos(θ4) 0 0
	 0 		 0 		 0 1
]

# ╔═╡ 97b3e997-f3e9-4fea-9863-c892557a338d
T45 = [
	cos(θ5) -sin(θ5) 0 0
	0 0 -1 0
	sin(θ5) cos(θ5) 0 0
	 0 		 0 		 0 1
]

# ╔═╡ 3ac89c31-f821-4217-822b-64f77729967b
T56 = [
	cos(θ6) -sin(θ6) 0 0
	0 0 1 0
	-sin(θ6) -cos(θ6) 0 0
	 0 		 0 		 0 1
]

# ╔═╡ e3d67624-4784-46a8-a495-10631bcc3e4a
Tlist = [T01, T12, T23, T34, T45, T56]

# ╔═╡ b343011d-8f6e-4b3f-8da2-9be085b45804
sympy.simplify(T01^-1)

# ╔═╡ bdcdcb96-d24d-4455-9e16-e3737d865dba
T12^-1

# ╔═╡ cea3d0ce-d717-49a1-942f-6a97d2b94ae4
simplify.(reduce(*, Tlist))

# ╔═╡ e5077ec3-6b48-4865-945d-088b8eadfcef
begin 
	function shorten_exp(exp)
		# regex for extracting args from sin and cos
		trig_arg = r"(?<=cos|sin)\(.+?\)"

		
		exp1 = replace(exp, trig_arg=>thetaReplace)
		exp2 = replace(exp1, "cos" => "c", "sin" => "s")
	end

	# returns comma seperated numbers corresponding to theta vals
	function thetaReplace(str)
		str_cpy = str
		str_cpy = Base.split(str_cpy, !isnumeric, keepempty=false)
		strfin = reduce((x,y)->x*","*y, string.(str_cpy))
		"_{"*strfin*"}"
	end
	shortEq(exp) = exp |> string |> shorten_exp |> latexstring
end

# ╔═╡ 08d75333-d848-41f9-b49a-9239b7654658
function eqsList(Ts, i, target)
	splitEq = split(Ts, i, target)
	eqs = collect(zip(splitEq[1], splitEq[2]))
	eqs2 = [string(x[1])*" = "*string(x[2]) for x in eqs]
	shortEq.(eqs2)[:]
end

# ╔═╡ 1c62a30e-2e44-4c55-b6f7-c483746847f2
eqsList(Ts, 3, target)

# ╔═╡ 2b87ed37-9daf-455e-9f17-9036cb4f5be5
eqsList(Tlist, 1, target)

# ╔═╡ 3285d327-f799-4abd-9e70-d4535f6bac88
typeof(shortEq(simplify(cos(θ1)cos(θ2)-sin(θ2)sin(θ1))))

# ╔═╡ 9a7b3b03-83d3-46fd-b9c1-5d6c3344e7d8


# ╔═╡ 83cf5f43-4691-4bf0-ab09-286f7600f547


# ╔═╡ 352dfa4b-859c-45b2-9133-1e28504857e0


# ╔═╡ e8ab5a3b-a114-452f-962f-5f35a3064ec9


# ╔═╡ 061f31b8-4778-4df0-9178-7837601664ac


# ╔═╡ a1435da6-fd8b-49b4-aaa6-0750bd305ba1
function invKine(x, y, z)
	θ0 = atan(y, x)
	xp = sqrt(x^2 + y^2)
	yp = z

	c2 = (xp^2 + yp^2 -l1^2 - l2^2)
	s2_1 = sqrt(1- c2^2)
	s2_2 = sqrt(1- c2^2)

	θ2_1 = atan(c2, s2_1)
	θ2_2 = atan(c2, s2_2)

	
	
	k1 = l1 + l2*c2
	k2_1 = l2*s2_1
	k2_2 = l2*s2_2

	k1 = r*cos(γ)
	θ1_1 = atan(y, x) - atan(k2_1, k1)
	θ1_2 = atan(y, x) - atan(k2_2, k1)

	θ3_1 = ϕ - θ1_1 - θ2_1
	θ3_1 = ϕ - θ1_2 - θ2_2

	return ([θ0, θ1_1, θ2_1, θ3_1],[θ0, θ1_2, θ2_2, θ3_2])
end

# ╔═╡ 8691754d-8e1a-4e09-b6d2-d46b2542871d
@vars d5

# ╔═╡ e6e7084f-44c6-4fb9-844b-4f727c7e124a
projDHs = [
	0 		0 		0 		θ1
	PI/2 	0	 	L1 		θ2
	0 		L2 		0 		θ3
	0 		L3 		0 		θ4
	0 		L4 		d5 		0
]

# ╔═╡ ae721655-fa6e-4b3b-98f5-baa0853968be
T_proj = transforms(projDHs)

# ╔═╡ b1218e42-ab66-4697-9b9c-4e8847a296bd
Tfinal = simplify.(reduce(*, T_proj))

# ╔═╡ 5139a3e1-da60-465e-b9af-140a3edc484e
eqsList(T_proj, 0, target)

# ╔═╡ 5c289151-54cf-48e6-9be4-e7955790f391
# Solution for theta 3
let
	eq1 = (Tfinal[1, 4], target[1, 4]).-(L1*sin(θ1) + L4*cos(θ1)*cos(θ2 + θ3 + θ4) + d5*sin(θ1))
	eq2 = (Tfinal[2, 4], target[2, 4]).-(-L1*cos(θ1) + L4*sin(θ1)*cos(θ2 + θ3 + θ4) - d5*cos(θ1))
	eq3 = (Tfinal[3, 4], target[3, 4]).-L4*sin(θ2 + θ3 + θ4)
	shortEqs = [shortEq.(x) for x in [eq1, eq2, eq3]]

	lh = simplify.(eq1.^2 .+ eq2.^2 .+ eq3.^2)[1]
	rh = simplify.(eq1.^2 .+ eq2.^2 .+ eq3.^2)[2]
	c3 = simplify((rh - (L2^2 + L3^2))/(2*L2*L3))
	s3 = sqrt(1-c3^2)
	theta3 = atan(c3, s3)
end

# ╔═╡ 6be2891e-a17a-442d-85cd-16ae138850fb
let
	s234 = getEq(T_proj, 0, target, 3)[1]
	c234 = getEq(T_proj, 0, target, 7)[1]
	theta234 = atan(c234, s234)
	theta4 = theta234 - θ3 - θ2
end

# ╔═╡ a928afbc-0c50-4d1a-8351-44cdf5d5fd8d
# Solution for theta 2
let
	eq1 = (Tfinal[1, 4], target[1, 4]).-(L1*sin(θ1) + L4*cos(θ1)*cos(θ2 + θ3 + θ4) + d5*sin(θ1))

	eq3 = (Tfinal[3, 4], target[3, 4]).-L4*sin(θ2 + θ3 + θ4)
	# shortEqs = [shortEq.(x) for x in [eq1, eq2, eq3]]

	a = L2*cos(θ1) + L3*cos(θ1)*cos(θ3)
	c = -L3*cos(θ1)*sin(θ3)
	e = L3*sin(θ3)
	f = L2 + L3*cos(θ3)

	d = eq1[2]
	g = eq3[2]
	# simplify(a*f-c*e)

	if (a*f - c*e) > 0
		theta2 = atan(a*g-d*e, d*f-c*g)
	else
		theta2 = atan(d*e-a*g, c*d-d*f)
	end

	theta2
end

# ╔═╡ 91e7cc11-ead7-4357-bcf6-450aa8482036


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
LaTeXStrings = "~1.3.0"
Latexify = "~0.15.17"
PlutoUI = "~0.7.48"
SymPy = "~1.1.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "05707eabb8148640957732d131ada7cdb7369c00"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "aaabba4ce1b7f8a9b34c015053d3b1edf60fa49c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.4.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "53b8b07b721b77144a0fbbbc2675222ebf40a02d"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.94.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "d12e612bba40d189cead6ff857ddb67bd2e6a387"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "de83b8c89b2744fee5279326fe8e3f4a9b94d1e1"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.1.7"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═625a9d5c-62c0-11ed-3363-bf1b1417eda1
# ╠═34827925-8a91-4ae9-ac10-2b9f412ac97c
# ╠═3aa6b114-e38f-4f35-80ab-490d0fbe94de
# ╠═12bb4f91-d0cd-41ba-9d2f-284ab51eb458
# ╠═f3a37ec8-347b-4221-866f-688bc6781f16
# ╠═30d171e2-f450-4aff-bfa3-52b372236dc5
# ╠═16e67336-e47b-4e56-bf6e-526c50c7854b
# ╠═c2568846-c99e-4e5c-89e3-34ab422b764f
# ╠═75f669cb-4313-4e35-a5b9-6c7fb28455bf
# ╠═d938636d-f678-4210-993a-d5a16f0ff2fa
# ╠═8726790c-f378-49df-897c-10dd3f07d1d4
# ╠═d3cddad4-c4c2-4e36-bcef-1f2b261e51dd
# ╠═27ab9d60-c4f1-4560-a316-088df657f55b
# ╠═08d75333-d848-41f9-b49a-9239b7654658
# ╠═22430021-1347-4d11-a5e7-159915eb9b70
# ╠═1c62a30e-2e44-4c55-b6f7-c483746847f2
# ╠═4b1d01d1-546e-4121-8eab-3c7162054fe6
# ╠═7a1f1531-f99c-4906-93b6-8b52701eb5fe
# ╠═7f8a3a00-77c4-4761-88f4-933c51e2d6b1
# ╠═208e0373-199c-4240-8da0-79eb588f9f7b
# ╠═aba4e8d1-9c4f-4dbc-938a-9697ce19a6dc
# ╠═97b3e997-f3e9-4fea-9863-c892557a338d
# ╠═3ac89c31-f821-4217-822b-64f77729967b
# ╠═e3d67624-4784-46a8-a495-10631bcc3e4a
# ╠═b343011d-8f6e-4b3f-8da2-9be085b45804
# ╠═bdcdcb96-d24d-4455-9e16-e3737d865dba
# ╠═2b87ed37-9daf-455e-9f17-9036cb4f5be5
# ╠═cea3d0ce-d717-49a1-942f-6a97d2b94ae4
# ╠═e5077ec3-6b48-4865-945d-088b8eadfcef
# ╠═3285d327-f799-4abd-9e70-d4535f6bac88
# ╠═9a7b3b03-83d3-46fd-b9c1-5d6c3344e7d8
# ╠═83cf5f43-4691-4bf0-ab09-286f7600f547
# ╠═352dfa4b-859c-45b2-9133-1e28504857e0
# ╠═e8ab5a3b-a114-452f-962f-5f35a3064ec9
# ╠═061f31b8-4778-4df0-9178-7837601664ac
# ╠═a1435da6-fd8b-49b4-aaa6-0750bd305ba1
# ╠═8691754d-8e1a-4e09-b6d2-d46b2542871d
# ╠═e6e7084f-44c6-4fb9-844b-4f727c7e124a
# ╠═ae721655-fa6e-4b3b-98f5-baa0853968be
# ╠═b1218e42-ab66-4697-9b9c-4e8847a296bd
# ╠═5139a3e1-da60-465e-b9af-140a3edc484e
# ╠═5c289151-54cf-48e6-9be4-e7955790f391
# ╠═6be2891e-a17a-442d-85cd-16ae138850fb
# ╠═a928afbc-0c50-4d1a-8351-44cdf5d5fd8d
# ╠═91e7cc11-ead7-4357-bcf6-450aa8482036
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
