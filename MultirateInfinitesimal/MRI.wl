(* ::Package:: *)

BeginPackage["MRI`"];
MRI::usage = "Package for analyzing MRI-GARK methods";

MRI::usage =
	"MRI[\[CapitalGamma]0, \[CapitalGamma]1, \[Ellipsis]] constructs an MRI-GARK method from coefficient matrices \[CapitalGamma]0, \[CapitalGamma]1, \[Ellipsis].\n" <>
	"MRI[] returns a Dataset containing a catalog of MRI-GARK methods.\n" <>
	"MRI[\"name\"] loads a method called name from the catalog.";
MRIGamma::usage =
	"MRIGamma[mri] gets the list of coefficients {\[CapitalGamma]0, \[CapitalGamma]1, \[Ellipsis]} of mri.\n" <>
	"MRIGamma[mri, k] gets coefficient matrix \[CapitalGamma]k. Note k is zero-indexed.";
MRIDegree::usage = "MRIDegree[mri] get the degree of the forcing polynomial of mri.";
MRIGammaBar::usage = "MRIGammaBar[mri] computes the integral of the forcing polynomial of mri.";
MRIStages::usage = "MRIStages[mri] returns the number of stages of mri.";
MRIA::usage = "MRIA[mri] returns the A coefficient matrix of the underlying Runge\[Dash]Kutta method of mri.";
MRIB::usage = "MRIB[mri] returns the b coefficients of the underlying Runge\[Dash]Kutta method of mri.";
MRIC::usage = "MRIC[mri] returns the c coefficients of the underlying Runge\[Dash]Kutta method of mri.";
MRIDeltaC::usage = "MRIDeltaC[mri] returns the distance between c's of mri.";
MRIInternalConsistency::usage = "MRIInternalConsistency[mri] computes the residual of the internal consistency condition for mri.";
MRIDecoupled::usage = "MRIDecoupled[mri] computes the residual of the decoupled stages condition for mri.";
MRIOrderConditions::usage =
	"MRIOrderConditions[mri, p] generate order conditions up to order p for mri\n" <>
	"MRIOrderConditions[mri, {p}] generates only the conditions of order p.";
MRIStability::usage =
	"MRIStability[mri, zf, zs] computes the linear stability function for mri applied to y'=\[Lambda]f*y+\[Lambda]s*y, where zf=h*\[Lambda]f and zs=h*\[Lambda]s.\n" <>
	"MRIStability[mri, zf, zs, s] computes the linear stability function for stage s.";


Begin["`Private`"];
(* Depends on Integreat 0.4.0 *)
Needs["Integreat`RK`"];
Needs["Integreat`BTrees`"];

Afs[m_MRI, -1, e_] := {MRIA[m] . e, \[FormalY]};
Afs[m_MRI, k_Integer, e_] := {MRIGamma[m, k] . e, \[FormalF][\[FormalY]*\[FormalF]^k]};
Afs[m_MRI, e_] := Table[Afs[m, i, e], {i, -1, MRIDegree[m]}];

LC[m_] := LowerTriangularize[ConstantArray[MRIDeltaC[m], MRIStages[m]], -1];

(* Root *)
SetAttributes[\[CapitalPhi], Listable];
\[CapitalPhi][m_MRI, BTreeN[Subscript[\[FormalF], 1], 2]] := 1;
\[CapitalPhi][m_MRI, BTreeN[Subscript[\[FormalF], 1][t_], 2]] := Total[Map[
	MRIDeltaC[m] . First[#] / BTreeGamma[BTree[\[FormalF][Last[#]]]] &,
	\[CapitalPhi][m, t, 1]
]];
\[CapitalPhi][m_MRI, BTreeN[Subscript[\[FormalF], 2], 2]] := Total[MRIB[m]];
\[CapitalPhi][m_MRI, BTreeN[Subscript[\[FormalF], 2][t_], 2]] := MRIB[m] . Total[\[CapitalPhi][m, t, 2][[All, 1]]];

(* Leaves *)
\[CapitalPhi][m_MRI, Subscript[\[FormalF], 1], 1] := {{MRIC[m], \[FormalY]}, {MRIDeltaC[m], \[FormalF]}};
\[CapitalPhi][m_MRI, Subscript[\[FormalF], 2], 1] := Afs[m, ConstantArray[1, MRIStages[m]]];
\[CapitalPhi][m_MRI, Subscript[\[FormalF], _], 2] := {{MRIC[m], \[FormalY]}};

(* Interior *)
\[CapitalPhi][m_MRI, Subscript[\[FormalF], 1][t_], 1] := Flatten[Map[
	{
		{LC[m] . First[#] / BTreeGamma[BTree[\[FormalF][Last[#]]]], \[FormalY]},
		{MRIDeltaC[m] * First[#], \[FormalF][Last[#]]}
	} &,
	\[CapitalPhi][m, t, 1]
], 1];
\[CapitalPhi][m_MRI, Subscript[\[FormalF], 1][t_], 2] := Map[
	{LC[m] . First[#] / BTreeGamma[BTree[\[FormalF][Last[#]]]], \[FormalY]} &,
	\[CapitalPhi][m, t, 1]
];
\[CapitalPhi][m_MRI, Subscript[\[FormalF], 2][t_], 1] := Flatten[
	Map[Afs[m, First[#]] &, \[CapitalPhi][m, t, 2]],
	1
];
\[CapitalPhi][m_MRI, Subscript[\[FormalF], 2][t_], 2] := Map[{MRIA[m] . First[#], \[FormalY]} &, \[CapitalPhi][m, t, 2]];

(* Operations *)
\[CapitalPhi][m_MRI, t_^n_, p_] := Times @@@ Tuples[ConstantArray[\[CapitalPhi][m, t, p], n]];
\[CapitalPhi][m_MRI, t_Times, p_] := Times @@@ Tuples[Map[\[CapitalPhi][m, #, p] &, List @@ t]];

SetAttributes[\[Phi], Listable];
\[Phi][z_, 0] := Exp[z];
\[Phi][z_, k_Integer] := Exp[z] * (Gamma[k] - Gamma[k, z]) / z^k;

\[Mu][m_MRI, z_] := Sum[MRIGamma[m, k] * \[Phi][z, k + 1], {k, 0, MRIDegree[m]}];


MRI[c__?SquareMatrixQ] := MRI[{c}];

Integreat`Internal`Catalog`AddCatalog[MRI,
	{
		"Backward Euler",
		MRI[{{1, 0, 0}, {-1, 0, 1}, {0, 0, 0}}]
	},
	{
		"MRI-GARK-IRK21a",
		MRI[{{1, 0, 0}, {-1/2, 0, 1/2}, {0, 0, 0}}]
	},
	{
		"MRI-GARK-ERK22a",
		MRI[{{1/2, 0}, {-1/2, 1}}]
	},
	{
		"MRI-GARK-ERK22b",
		MRI[{{1, 0}, {-1/2, 1/2}}]
	},
	{
		"MRI-GARK-ERK33a",
		MRI[
			{{1/3, 0, 0}, {-(1/3), 2/3, 0}, {0, -(2/3), 1}},
			{{0, 0, 0}, {0, 0, 0}, {1/2, 0, -(1/2)}}
		]
	},
	{
		"MRI-GARK-ERK45a",
		MRI[
			{
				{1/5, 0, 0, 0, 0},
				{-(53/16), 281/80, 0, 0, 0},
				{-(36562993/71394880), 34903117/17848720, -(88770499/71394880), 0, 0},
				{-(7631593/71394880), -(166232021/35697440), 6068517/1519040, 8644289/8924360, 0},
				{277061/303808, -(209323/1139280), -(1360217/1139280), -(148789/56964), 147889/45120}
			},
			{
				{0, 0, 0, 0, 0},
				{503/80, -(503/80), 0, 0, 0},
				{-(1365537/35697440), 4963773/7139488, -(1465833/2231090), 0, 0},
				{66974357/35697440, 21445367/7139488, -3, -(8388609/4462180), 0},
				{-(18227/7520), 2, 1, 5, -(41933/7520)}
			}
		]
	}
];

MRIGamma[m_MRI] := m[[1]];
MRIGamma[m_MRI, k_Integer?NonNegative] := MRIGamma[m][[k + 1]];
MRIDegree[m_MRI] := Length[MRIGamma[m]] - 1;
MRIGammaBar[m_MRI] := (1 / Range[MRIDegree[m] + 1]) . MRIGamma[m];
MRIStages[m_MRI] := Length[First[MRIGamma[m]]];
MRIA[m_MRI] := Table[Boole[i > j], {i, MRIStages[m]}, {j, MRIStages[m]}] . MRIGammaBar[m];
MRIB[m_MRI] := Total[MRIGammaBar[m]];
MRIC[m_MRI] := Total[MRIA[m], {2}];
MRIDeltaC[m_MRI] := Differences[Append[MRIC[m], 1]];

MRIInternalConsistency[m_MRI] := {
	Total[MRIGamma[m, 0], {2}] - MRIDeltaC[m],
	Total[Rest[MRIGamma[m]], {3}]
};

MRIDecoupled[m_MRI] := Most[MRIDeltaC[m]] * Diagonal[MRIGammaBar[m], 1];

MRIOrderConditions[m_MRI, p:_Integer?Positive | {_Integer?NonNegative}] := With[{
		t = BTreeN[p]
	},
	(\[CapitalPhi][m, t] - 1 / BTreeGamma[t]) / BTreeSigma[t]
];

MRIStability[m_MRI, zf_, zs_, s_] := With[{
		X = DiagonalMatrix[\[Phi][MRIDeltaC[m] * zf, 0]] + zs * \[Mu][m, MRIDeltaC[m] * zf],
		L = Table[Boole[i == j + 1], {i, MRIStages[m]}, {j, MRIStages[m]}]
	},
	Part[Inverse[IdentityMatrix[MRIStages[m]] - X . L] . X[[All, 1]], s]
];
MRIStability[m_MRI, zf_, zs_] := MRIStability[m, zf, zs, -1];

MRI /: RK[m_MRI] := RK[MRIA[m], MRIB[m], MRIC[m]];

MRI /: Variables[m_MRI] := Variables[MRIGamma[m]];

MRI /: MakeBoxes[m:MRI[__?SquareMatrixQ], format_] := With[{
		es = Append[ConstantArray[False, MRIStages[m] - 1], True]
	},
	GridBox[
		Map[MakeBoxes[#, format] &, ArrayFlatten[{MRIGamma[m]}], {2}],
		ColumnLines -> Join @@ ConstantArray[es, MRIDegree[m] + 1]
	]
];


End[];
EndPackage[];
