BeginPackage["INVARIANT`Structure`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
exportXYZ                    ::usage = "exportXYZ[filepath_, {comment_, vertices_, atomcoordinates_}]"
CifWyckoffImages             ::usage = "CifWyckoffImages[file, pos]"
cif2mcif                     ::usage = "cif2mcif[id, IsoMatrix, pos0]"
xyz2Matrix                   ::usage = "xyz2Matrix[expr]"
SymmetryOpOrigin             ::usage = "SymmetryOpOrigin[file, set]"
SymmetryOpVectorField        ::usage = "SymmetryOpVectorField[file, set]"
GetLatticeVectors            ::usage = "GetLatticeVectors[dat]"
ImportIsodistortCIF          ::usage = "ImportIsodistortCIF[\*StyleBox[\"filename\",\"TI\"]] opens .cif file and imports relevant data."
PosMatchTo                   ::usage = "PosMatchTo[spos, dpos, tol]"
GetStrainTensor              ::usage = "GetStrainTensor[parent, sub]"
GridPbc                      ::usage = "GridPbc[ixyz, Lxyz]"
GridNeighbors                ::usage = "GridNeighbors[r0, MeshDim]"
PbcDiff                      ::usage = "PbcDiff[diff]"
SiteJij                      ::usage = "SiteJij[file, PosVec]"
GetUpTo3rdNeighbors          ::usage = "GetUpTo3rdNeighbors[]"
MakeSuperCell                ::usage = "MakeSuperCell[Crys, {Nx, Ny, Nz}]"

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
exportXYZ[filepath_, {comment_, vertices_, atomcoordinates_}] := Module[{outputfilestream, outputstring},
  If[Or[Not@StringQ[filepath], Not@StringQ[comment]], Message[exportXYZ::usage]; Abort[]];
  If[Length@vertices != First@Dimensions@atomcoordinates, Message[exportXYZ::mismatch]; Abort[]];
  outputstring = ToString@Length[vertices] <> "\n" <> comment <> "\n" <> MapThread[(#1 <> "\t" <> ExportString[{#2}, "TSV"] <> "" &), {vertices, atomcoordinates}];
  outputfilestream = OpenWrite[filepath];
  WriteString[outputfilestream, outputstring];
  Close[outputfilestream];
Return[outputstring]]

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

GetIRvector[id_, IsoMatrix_, pos_] := Module[{IRvector},
  IRvector = {If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], # & /@ (pos\[Transpose][[2]])}\[Transpose];
  Return[IRvector]
]

cif2mcif[file_, id_, IsoMatrix_, pos0_] := Module[{CifTagVec, mcif, cif, vectors},
  CifTagVec = {{}, {"loop_"}, {"_atom_site_moment.label"}, {"_atom_site_moment.crystalaxis_x"}, {"_atom_site_moment.crystalaxis_y"}, {"_atom_site_moment.crystalaxis_z"}, {"_atom_site_moment.symmform"}};
  cif = Import[file, "Table"];
  vectors = {First@StringCases[#2, RegularExpression["[[:upper:]][[:lower:]]*"]] <> If[dd = StringCases[#2, RegularExpression["\\d+"]]; dd === {}, Unevaluated[Sequence[]], dd], #1[[1]], #1[[2]], #1[[3]], "mx,my,mz"} & @@@ GetIRvector[id, IsoMatrix, pos0];
  mcif = Join[cif, CifTagVec, vectors];
  Return[mcif]
  ]

SymmetryOpOrigin[file_, set_] := Module[{CifData, CifFlags, xyzName, xyzStrData, xyzExpData, newsets, op, posmap, diff, i, j},
  If[! FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[]];
  CifData = Import[file, "string"] <> "\n";
  (*---fix for some strange behaviour---*)
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] -> ""]];
  CifData = ImportString[CifData, "CIF"];
  CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
  xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
  xyzStrData = Flatten[xyzName /. CifData];
  xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  newsets = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #1] & @@@set, {op, xyzExpData}];
  posmap = Table[Flatten[Table[diff = Mod[newsets[[op]][[i]] - set[[j]][[1]], 1]; If[diff == {0., 0., 0.}, {j, i, set\[Transpose][[2]][[j]]},Unevaluated[Sequence[]]], {i, Length@set}, {j, Length@set}], 1], {op, Length@xyzExpData}];
  (*posmap = Table[SortBy[op, First], {op, posmap}];*)
  newsets = Table[{newsets[[op]], Transpose[posmap[[op]]][[3]]}\[Transpose], {op, Length@xyzExpData}];
  Return[newsets]
]

DistMatrixOld = Compile[{{pos1, _Real, 2}, {pos2, _Real, 2}}, Module[{i, j},
   Return[Table[Norm[Mod[Round[pos1[[i]] - pos2[[j]], 10^-8], 1]], {i, Length@pos1}, {j, Length@pos2}]]
]]

DistMatrix = Compile[{{pos1, _Real, 2}, {pos2, _Real, 2}}, Module[{i, j},
   Return[Table[Norm[Round[Which[# > 0.5, # - 1, # <= -0.5, # + 1, True, #] & /@ (pos1[[i]] - pos2[[j]]), 10^-8]^2], {i, Length@pos1}, {j, Length@pos2}]]
]]

PosMatchTo[spos_, dpos_, tol_, OptionsPattern["shift"->False]] := Module[{difftable, posmap, i, newpos},
  difftable = DistMatrix[spos, dpos];
  posmap = Position[difftable, x_ /; TrueQ[x <= tol]];
  newpos = If[OptionValue["shift"], PbcDiff[#]&/@(Table[dpos[[posmap[[i]][[2]]]], {i, Length@spos}] - spos) + spos,Table[dpos[[posmap[[i]][[2]]]], {i, Length@spos}]];
  Return[{posmap, newpos}]
]

SymmetryOpVectorField[file_, pos_, vec_, OptionsPattern["spin" -> False]] := Block[{originshift, CifData, CifFlags, xyzName, xyzStrData, xyzRotTranData, xyzTranslation, xyzRotData, field, newpos, newvec, difftable, diff, posmap, i, j, axial}, 
  If[! FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[]];
  CifData = Import[file, "string"] <> "\n";
  originshift = ToExpression[StringReplace[StringCases[CifData,RegularExpression["origin=\\W\\d\\.*\\d*,\\d\\.*\\d*,\\d\\.*\\d*\\W"]][[1]], {"origin=" -> "", "(" -> "{", ")" -> "}"}]];
  (*---fix for some strange behaviour---*)
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] ->""]];
  IsoDispModes = StringCases[CifData, RegularExpression["([[:upper:]][[:lower:]]*\\W*|[[:upper:]]*\\d*[[:lower:]]*)_\\d+"]];
  CifData = StringReplace[CifData, Thread[IsoDispModes -> StringReplace[IsoDispModes, {"*" -> "conj", "_" -> ""}]]];
  CifData = ImportString[CifData, "CIF"];
  CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];

  xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
  xyzStrData = Flatten[xyzName /. CifData];
  xyzRotTranData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  xyzTranslation = xyzRotTranData  /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0} ;
  xyzRotData = xyzRotTranData - xyzTranslation;
  
  axial = If[OptionValue["spin"], 1, 2];
  
  newvec = Table[{Det[xyz2Matrix[op]]^axial*N[op /. Thread[ToExpression[{"x", "y", "z"}] -> #1]], #2} & @@@ vec, {op, xyzRotData}];
  newpos = Table[{N[op /. Thread[ToExpression[{"x", "y", "z"}] -> (pos[[i]][[1]]+originshift)]], i}, {op, xyzRotTranData}, {i, Length@pos}];

  (*difftable = ParallelTable[DistMatrix[#+originshift&/@(pos\[Transpose][[1]]), newpos[[op]]\[Transpose][[1]]], {op, Length@xyzRotTranData}, DistributedContexts -> {"INVARIANT`Structure`Private`"}];*)
  difftable = Table[DistMatrix[#+originshift&/@(pos\[Transpose][[1]]), newpos[[op]]\[Transpose][[1]]], {op, Length@xyzRotTranData}];
  posmap = Table[Position[difftable[[op]], x_ /; TrueQ[x == 0]], {op, Length@xyzRotTranData}];
  field = Table[Table[{newvec[[i]][[posmap[[i]][[j]][[2]]]][[1]], pos\[Transpose][[2]][[j]]}, {j, Length@pos}], {i, Length@posmap}];
  Return[field]
]

xyz2Matrix[expr_] := Module[{m},
  m = Coefficient[#, {ToExpression["x"], ToExpression["y"], ToExpression["z"]}] & /@ expr;
  If[Abs[Det[Normalize[#]&/@m]] != 1, Print["xyz2Matrix gives wrong determinant"]; Abort[], Unevaluate[Sequence[]]];
  Return[m]
]


GetUpTo3rdNeighbors[] := Module[{FirstNeighborList, SecondNeighborList, ThirdNeighborList},
  FirstNeighborList = Flatten[{-1, 1}\[TensorProduct]IdentityMatrix[3], 1];
  SecondNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {2}], {0, 0, 0}];
  ThirdNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {3}], a_ /; MemberQ[FirstNeighborList, a]];
  Return[{FirstNeighborList, SecondNeighborList, ThirdNeighborList}]
]

SiteJij[file_, PosVec_] := Block[{CifData, CifFlags, xyzName, xyzStrData, xyzRotTranData, xyzTranslation, xyzRotData, field, newpos, newvec, difftable, diff, posmap, i, j}, 
  If[! FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[]];
  CifData = Import[file, "string"] <> "\n";
  (*---fix for some strange behaviour---*)
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] -> ""]];
  IsoDispModes = StringCases[CifData, RegularExpression["([[:upper:]][[:lower:]]*\\W*|[[:upper:]]*\\d*[[:lower:]]*)_\\d+\
"]];
  CifData = StringReplace[CifData, Thread[IsoDispModes -> StringReplace[IsoDispModes, {"*" -> "conj", "_" -> ""}]]];
  CifData = ImportString[CifData, "CIF"];
  CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
  xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
  xyzStrData = Flatten[xyzName /. CifData];
  xyzRotTranData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  xyzTranslation = xyzRotTranData /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0};
  xyzRotData = xyzRotTranData - xyzTranslation;
  
  newvec = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[2]]), {op, xyzRotData}];
  newpos = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[1]]), {op, xyzRotData}];
  
  Return[{#[[1]], #[[2]]}\[Transpose] &/@ ({newpos, newvec}\[Transpose])]
]

SymmetryOpVectorFieldOld[file_, pos_, vec_] := Block[{originshift, CifData, CifFlags, xyzName, xyzStrData, xyzRotTranData, xyzTranslation, xyzRotData, field, newpos, newvec, diff, posmap, i, j},
  If[! FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[]];
  CifData = Import[file, "string"] <> "\n";
  originshift = ToExpression[StringReplace[StringCases[CifData,RegularExpression["origin=\\W\\d\\.*\\d*,\\d\\.*\\d*,\\d\\.*\\d*\\W"]][[1]], {"origin=" -> "", "(" ->   "{", ")" -> "}"}]];
  (*---fix for some strange behaviour---*)
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] ->""]];
  IsoDispModes = StringCases[CifData, RegularExpression["([[:upper:]][[:lower:]]*\\W*|[[:upper:]]*\\d*[[:lower:]]*)_\\d+"]];
  CifData = StringReplace[CifData, Thread[IsoDispModes -> StringReplace[IsoDispModes, {"*" -> "conj", "_" -> ""}]]];
  CifData = ImportString[CifData, "CIF"];
  CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];

  xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
  xyzStrData = Flatten[xyzName /. CifData];
  xyzRotTranData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  xyzTranslation = xyzRotTranData  /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0} ;
  xyzRotData = xyzRotTranData - xyzTranslation;

  newvec = Table[{N[op /. Thread[ToExpression[{"x", "y", "z"}] -> #1]], #2} & @@@ vec, {op, xyzRotData}];
  newpos = Table[{N[op /. Thread[ToExpression[{"x", "y", "z"}] -> (#1+originshift)]], First@StringCases[#2, RegularExpression["[[:upper:]][[:lower:]]*"]]} & @@@ pos,  {op, xyzRotTranData}];
  posmap = ParallelTable[Flatten[Table[diff = Norm[Mod[Round[newpos[[op]][[i]][[1]] - N[(pos[[j]][[1]]+originshift)], 10^-10], 1]]; If[diff < 10^-10, {j, i, pos\[Transpose][[2]][[j]]}, Unevaluated[Sequence[]]], {i, Length@pos}, {j, Length@pos}], 1], {op, Length@xyzRotTranData}, DistributedContexts -> {"INVARIANT`Structure`Private`"}];
  posmap = Table[SortBy[op, First], {op, posmap}];
  field = Table[Table[{newvec[[i]][[posmap[[i]][[j]][[2]]]][[1]], posmap[[i]][[j]][[3]]}, {j, Length@pos}], {i, Length@posmap}];
  Return[field]
]

GetLatticeVectors[dat_] := Module[{a1, b1, c1, \[Alpha]1, \[Beta]1, \[Gamma]1,bt},
     a1 = dat[[1]]; b1 = dat[[2]]; c1 = dat[[3]]; \[Alpha]1 = dat[[4]]; \[Beta]1 = dat[[5]]; \[Gamma]1 = dat[[6]];
     bt = {{a1, 0, 0},
           {b1 Cos[\[Gamma]1], b1 Sin[\[Gamma]1], 0},
           {c1 Cos[\[Beta]1], c1 (Cos[\[Alpha]1] - Cos[\[Beta]1] Cos[\[Gamma]1])/Sin[\[Gamma]1],
            c1 Sqrt[1 - (Cos[\[Alpha]1]^2 + Cos[\[Beta]1]^2 + Cos[\[Gamma]1]^2) + 2 Cos[\[Alpha]1] Cos[\[Beta]1] Cos[\[Gamma]1]]/Sin[\[Gamma]1]}};
     Return[bt]
]

ImportIsodistortCIF[file_,OptionsPattern[]] := Module[{CifData, CifFlags, xyzName, ElemSym, ElemName, SpgName, xyzStrData, xyzExpData, xyzTranslation, Atoms, EveryAtomPos, LatticeVectors, LatticeScales, sol, i, j, k, failed, x,y,z, SpgNumber, aa,bb,cc,alpha,beta,gamma,FracOutText,Wyckoff,IsoDim,IsoDispModeRow,IsoDispModeCol,IsoDispModeVal,IsoDispModeMatrix, IsoDispModesID,IsoDispModesNam,IsoDispModesAmp,IsoDispModes,IsoDispModeNormVal,IsoDeltaCoordID,IsoDeltaCoordLabel,IsoDeltaCoordVal,IsoDeltaCoord, IsoCoordLabel,IsoCoordFormula,IsoCoord, IsoStrainModesID,IsoStrainModesNam,IsoStrainModesAmp,IsoStrainModes, IsoStrainModeRow, IsoStrainModeCol, IsoStrainModeVal},

    If[!FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[] ];
	CifData = Import[file, "string"] <> "\n"; (*--- fix for some strange behaviour ---*)
    CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
    CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData,RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] -> ""]];
    IsoDispModes = StringCases[CifData, RegularExpression["([[:upper:]][[:lower:]]*\\W*|[[:upper:]]*\\d*[[:lower:]]*)_\\d+"]]; 
    CifData = StringReplace[CifData, Thread[IsoDispModes -> StringReplace[IsoDispModes, {"*" -> "conj", "'" -> "Prime", "_" -> ""}]]];

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

    IsoStrainModeRow = If[Position[CifFlags, "_iso_strainmodematrix_row"] != {}, "_iso_strainmodematrix_row" /. CifData,   ""];
    IsoStrainModeCol = If[Position[CifFlags, "_iso_strainmodematrix_col"] != {}, "_iso_strainmodematrix_col" /. CifData,   ""];
    IsoStrainModeVal = If[Position[CifFlags, "_iso_strainmodematrix_value"] != {},                                             Rationalize["_iso_strainmodematrix_value" /. CifData], ""];

    IsoStrainModesID = If[Position[CifFlags, "_iso_strainmode_ID"] != {}, "_iso_strainmode_ID" /. CifData, ""];
    IsoStrainModesNam = If[Position[CifFlags, "_iso_strainmode_label"] != {}, "_iso_strainmode_label" /. CifData, ""];
    IsoStrainModesAmp = If[Position[CifFlags, "_iso_strainmode_value"] != {}, "_iso_strainmode_value" /. CifData, ""];
    IsoStrainModes = Transpose[{IsoStrainModesID, IsoStrainModesNam, IsoStrainModesAmp}];

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
	{{ElemSym, ElemName}, "","", SpgName, SpgNumber, LatticeVectors, Atoms, LatticeScales, FracOutText, Wyckoff, IsoDim, IsoDispModes, {IsoDispModeRow, IsoDispModeCol, IsoDispModeVal}, IsoStrainModes, {IsoStrainModeRow, IsoStrainModeCol, IsoStrainModeVal}, IsoDeltaCoord, IsoCoord}
]

GetStrainTensor[parent_, sub_, OptionsPattern["iso"->True]] := Module[{e, strain},
  strain= If[OptionValue["iso"],
             e = N[sub[[6]] /. sub[[8]]].Inverse[N[parent[[6]] /. parent[[8]]]] - IdentityMatrix[3];
             1/2 (e + Transpose[e] + Transpose[e].e),
             e = sub.Inverse[parent] - IdentityMatrix[3];
             1/2 (e + Transpose[e] + Transpose[e].e)]
]

GridPbc[ixyz_, Lxyz_] := Module[{}, Mod[ixyz - 1, Lxyz] + 1]

GridNeighbors[r0_, MeshDim_] := Module[{r00, list, FirstNeighborList, SecondNeighborList, ThirdNeighborList, Neighbours, a},
  FirstNeighborList = Flatten[{-1, 1}\[TensorProduct]IdentityMatrix[3], 1];
  SecondNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {2}], {0, 0, 0}];
  ThirdNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {3}], a_ /; MemberQ[FirstNeighborList, a]];
  r00 = GridPbc[r0, MeshDim];
  Neighbours = Table[GridPbc[r00 + #, MeshDim] & /@ list, {list, {FirstNeighborList, SecondNeighborList, ThirdNeighborList}}];
  Return[Neighbours]
]

PbcDiff[diff_] := Module[{}, 
  Which[-0.5 <= # < 0.5, #, # >= 0.5, # - 1, # < -0.5, # + 1] & /@ diff
]

MakeSuperCell[Crys_, mode_, k_, {Nx_, Ny_, Nz_}] := Module[{R, pos, ix, iy, iz, spos, sR},
  {R, pos} = Crys;
  spos = Flatten[Table[{Chop[(#1 + {ix, iy, iz})/{Nx, Ny, Nz} + #3 Exp[I 2 Pi {ix, iy, iz}.k]], #2}, {ix, 0, Nx - 1}, {iy, 0, Ny - 1}, {iz, 0, Nz - 1}] & @@@ (Append[pos\[Transpose], mode]\[Transpose]), 3];
  sR = DiagonalMatrix[{Nx, Ny, Nz}].R;
  Return[{sR, spos}]
]


(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
