BeginPackage["INVARIANT`Structure`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
CifWyckoffImages    ::usage = "CifWyckoffImages[file, pos]"
GetLatticeVectors      ::usage = "GetLatticeVectors[dat]"
ImportIsodistortCIF ::usage = "ImportIsodistortCIF[\*StyleBox[\"filename\",\"TI\"]] opens .cif file and imports relevant data."
GetStrainTensor     ::usage = "GetStrainTensor[parent, sub]"
PbcDiff             ::usage = "PbcDiff[diff]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)
Options[ImportIsodistortCIF]    = {Fractional->False, CorrectLabels->True, Tolerance->10^-6}

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
CifWyckoffImages[file_, pos_] := Module[{CifData, CifFlags, xyzName, xyzStrData, xyzExpData, WyckoffImages},
     If[!FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[] ];
     CifData = Import[file, "string"] <> "\n"; (*--- fix for some strange behaviour ---*)
     CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
     CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData,RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] -> ""]];
 
     CifData = ImportString[CifData, "CIF"];
     CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
 
     xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
     xyzStrData = Flatten[xyzName /. CifData];
     xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
 
     WyckoffImages = Flatten[Table[Sort[Union[{#, pos[[i, 2]]} & /@ (xyzExpData /. Thread[ToExpression[{"x", "y", "z"}] -> pos[[i, 1]]])], Total[#1[[1]]] > Total[#2[[1]]] &], {i, Length[pos]}], 1];
     Return[WyckoffImages];
]

GetLatticeVectors[dat_] := Module[{a1, b1, c1, \[Alpha]1, \[Beta]1, \[Gamma]1,bt},
     a1 = dat[[1]]; b1 = dat[[2]]; c1 = dat[[3]]; \[Alpha]1 = dat[[4]]; \[Beta]1 = dat[[5]]; \[Gamma]1 = dat[[6]];
     bt = {{a1, 0, 0},
           {b1 Cos[\[Gamma]1], b1 Sin[\[Gamma]1], 0},
           {c1 Cos[\[Beta]1], c1 (Cos[\[Alpha]1] - Cos[\[Beta]1] Cos[\[Gamma]1])/Sin[\[Gamma]1],
            c1 Sqrt[1 - (Cos[\[Alpha]1]^2 + Cos[\[Beta]1]^2 + Cos[\[Gamma]1]^2) + 2 Cos[\[Alpha]1] Cos[\[Beta]1] Cos[\[Gamma]1]]/Sin[\[Gamma]1]}};
     Return[bt]
]

ImportIsodistortCIF[file_,OptionsPattern[]] := Module[{CifData, CifFlags, xyzName, ElemSym, ElemName, SpgName,
  xyzStrData, xyzExpData, xyzTranslation, Atoms, EveryAtomPos, LatticeVectors, LatticeScales, sol, i, j, k, failed,
  x,y,z, SpgNumber, aa,bb,cc,alpha,beta,gamma,FracOutText,Wyckoff,IsoDim,IsoDispModeRow,IsoDispModeCol,IsoDispModeVal,IsoDispModeMatrix,
  IsoDispModesID,IsoDispModesNam,IsoDispModesAmp,IsoDispModes,IsoDispModeNormVal,IsoDeltaCoordID,IsoDeltaCoordLabel,IsoDeltaCoordVal,IsoDeltaCoord,
  IsoCoordLabel,IsoCoordFormula,IsoCoord},

    If[!FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[] ];
	CifData = Import[file, "string"] <> "\n"; (*--- fix for some strange behaviour ---*)
    CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
    CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData,RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] -> ""]];
    IsoDispModes = StringCases[CifData, RegularExpression["([[:upper:]]\\W|[[:upper:]])_\\d+"]]; 
    CifData = StringReplace[CifData, Thread[IsoDispModes -> StringReplace[IsoDispModes, {"*" -> "conj", "_" -> ""}]]];

	CifData = ImportString[CifData, "CIF"];
	CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
	ElemSym = If[Position[CifFlags, "_chemical_formula_sum"] != {},"_chemical_formula_sum" /. CifData,""];
	ElemName = If[Position[CifFlags, "_chemical_name_mineral"] != {},"_chemical_name_mineral" /. CifData,""];

    IsoDim = If[Position[CifFlags, "_iso_displacivemode_number"] != {}, "_iso_displacivemode_number" /. CifData, 0];
    IsoDispModeRow = If[Position[CifFlags, "_iso_displacivemodematrix_row"] != {}, "_iso_displacivemodematrix_row" /. CifData, ""];
    IsoDispModeCol = If[Position[CifFlags, "_iso_displacivemodematrix_col"] != {}, "_iso_displacivemodematrix_col" /. CifData, ""];
    IsoDispModeVal = If[Position[CifFlags, "_iso_displacivemodematrix_value"] != {}, Rationalize["_iso_displacivemodematrix_value" /. CifData], ""];
    (*IsoDispModeMatrix = Normal[SparseArray[Thread[Transpose[{IsoDispModeRow, IsoDispModeCol}] -> IsoDispModeVal], {IsoDim, IsoDim}]];*)

    IsoDispModesID = If[Position[CifFlags, "_iso_displacivemode_ID"] != {}, "_iso_displacivemode_ID" /. CifData, ""];
    IsoDispModesNam = If[Position[CifFlags, "_iso_displacivemode_label"] != {}, "_iso_displacivemode_label" /. CifData, ""];
    IsoDispModesAmp = If[Position[CifFlags, "_iso_displacivemode_value"] != {}, "_iso_displacivemode_value" /. CifData, ""];
    IsoDispModeNormVal = If[Position[CifFlags, "_iso_displacivemodenorm_value"] != {},"_iso_displacivemodenorm_value" /. CifData, ""];
    IsoDispModes = Transpose[{IsoDispModesID, IsoDispModesNam, IsoDispModeNormVal, IsoDispModesAmp}];

    IsoDeltaCoordID = If[Position[CifFlags, "_iso_deltacoordinate_ID"] != {}, "_iso_deltacoordinate_ID" /. CifData, ""];
    IsoDeltaCoordLabel = If[Position[CifFlags, "_iso_deltacoordinate_label"] != {}, "_iso_deltacoordinate_label" /. CifData, ""];
    IsoDeltaCoordVal = If[Position[CifFlags, "_iso_deltacoordinate_value"] != {}, "_iso_deltacoordinate_value" /. CifData, ""];
    IsoDeltaCoord = Transpose[{IsoDeltaCoordID, IsoDeltaCoordLabel, IsoDeltaCoordVal}];

    IsoCoordLabel = If[Position[CifFlags, "_iso_coordinate_label"] != {}, "_iso_coordinate_label" /. CifData, ""];
    IsoCoordFormula = If[Position[CifFlags, "_iso_coordinate_formula"] != {}, "_iso_coordinate_formula" /. CifData, ""];
    IsoCoord = Transpose[{IsoCoordLabel, IsoCoordFormula}];

	(*--- get lattice ---*)
	aa = If[Position[CifFlags, "_cell_length_a"] != {},Rationalize["_cell_length_a" /. CifData], Print["lattice error"]; Abort[]];
	bb = If[Position[CifFlags, "_cell_length_b"] != {},Rationalize["_cell_length_b" /. CifData], Print["lattice error"]; Abort[]];
	cc = If[Position[CifFlags, "_cell_length_c"] != {},Rationalize["_cell_length_c" /. CifData], Print["lattice error"]; Abort[]];
	alpha = If[Position[CifFlags, "_cell_angle_alpha"] != {},Rationalize["_cell_angle_alpha" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	beta = If[Position[CifFlags, "_cell_angle_beta"] != {},Rationalize["_cell_angle_beta" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	gamma = If[Position[CifFlags, "_cell_angle_gamma"] != {},Rationalize["_cell_angle_gamma" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	LatticeVectors = GetLatticeVectors[{aa,bb,cc,alpha,beta,gamma}];
	LatticeScales = Table[ToExpression[FromCharacterCode[96+i]] -> Simplify@Norm[LatticeVectors[[i]]] ,{i,Length[LatticeVectors]}];
	LatticeVectors = Table[ToExpression[FromCharacterCode[96+i]]* Simplify@Normalize[LatticeVectors[[i]]] ,{i,Length[LatticeVectors]}];

    If[Position[CifFlags,"_atom_site_fract_x" | "_atom_site_fract_y" |"_atom_site_fract_z"] != {},
       Wyckoff = Transpose[{Transpose[({"_atom_site_fract_x","_atom_site_fract_y", "_atom_site_fract_z"} /. CifData)],"_atom_site_label" /. CifData}];,
       If[Position[CifFlags,"_atom_site_Cartn_x" | "_atom_site_Cartn_y" |"_atom_site_Cartn_z"] != {},
          Print["cartesian coords found, not implemented"];Abort[], 
          Print["no atom_side_fract"];Abort[]
       ];
    ];

	(*--- get multiplicity and transform to cartesian ---*)
    xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
    xyzStrData = Flatten[xyzName /. CifData];
    xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];

	(*--- apply multiplicity ---*)
    EveryAtomPos = Flatten[Table[Sort[Union[{Mod[#, 1], Wyckoff[[i, 2]]} & /@ (xyzExpData /. Thread[ToExpression[{"x", "y", "z"}] -> Wyckoff[[i, 1]]])], Total[#1[[1]]] > Total[#2[[1]]] &], {i, Length[("_atom_site_label" /. CifData)]}], 1];
    Atoms = {#[[1]], ElementLabel = Characters[#[[2]]]; ToUpperCase[ElementLabel[[1]]] <> If[Length[ElementLabel] > 1 && LetterQ[ElementLabel[[2]]], ElementLabel[[2]], ""]} & /@ EveryAtomPos;

	(*--- filter equal length ---*)
	{x,y,z} = {a,b,c} /. LatticeScales;
	If[ PossibleZeroQ[y-z],
		LatticeScales = Delete[LatticeScales,3];
		Atoms = Atoms /.{ c -> b};
		LatticeVectors = LatticeVectors /.{c->b},
	  If[ PossibleZeroQ[x-z],
		LatticeScales = Delete[LatticeScales,3];
		Atoms = Atoms /.{c->a};
		LatticeVectors = LatticeVectors /.{c->a}]
	];
	If[ PossibleZeroQ[x-y],
		LatticeScales = Delete[LatticeScales,2];
		Atoms = Atoms /.{b->a};
		LatticeVectors = LatticeVectors /.{b->a}];
	Clear[x,y,z];

	(*--- space group names ---*)
	SpgNumber = If[Position[CifFlags, "_space_group_IT_number"] != {},"_space_group_IT_number" /. CifData,0];
	If[SpgNumber!=0 && IntegerQ[SpgNumber],
		SpgName = GTSpaceGroups[SpgNumber, GOVerbose -> False][[4]]
	,
		If[Position[CifFlags, "_symmetry_space_group_name_H-M"] != {},
			SpgName = "_symmetry_space_group_name_H-M" /. CifData;
			SpgName = StringJoin[DeleteCases[Characters[SpgName], " "]];
			SpgName = StringReplace[SpgName, {"m3m" -> "m-3m"}];
			SpgName = BracketingBar[ToExpression[SpgName]],
			Print["space group name error"];
			SpgName = ""
		]
	];

    FracOutText="fractional";

	If[ValueQ[failed], Print[failed]];
	{{ElemSym, ElemName}, "","", SpgName, SpgNumber, LatticeVectors, Atoms, LatticeScales, FracOutText, Wyckoff, IsoDim, IsoDispModes, {IsoDispModeRow, IsoDispModeCol, IsoDispModeVal}, IsoDeltaCoord, IsoCoord}
]

GetStrainTensor[parent_, sub_] := Module[{e, strain},
  e = N[sub[[6]] /. sub[[8]]].Inverse[N[parent[[6]] /. parent[[8]]]] - IdentityMatrix[3];
  strain = 1/2 (e + Transpose[e] + Transpose[e].e)
]

PbcDiff[diff_] := Module[{}, 
  Which[-0.5 <= # < 0.5, #, # >= 0.5, # - 1, # < -0.5, # + 1] & /@ diff
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
