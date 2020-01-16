BeginPackage["INVARIANT`ISODISTORT`",{"INVARIANT`Structure`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ShowIsoModes                 ::usage = "ShowIsoModes[PosVec]"
ISODISTORT                   ::usage = "ISODISTORT[R0, pos0, pos, IsoMatrix, label]"
ImposeIsoMode                ::usage = "ImposeIsoMode[Wyckoff0, IsoMatrix, modeset, s]"
GetIRvector                  ::usage = "GetIRvector[id, pos]"
GetOpMatrix                  ::usage = "GetOpMatrix[SymOpFile, pos, IsoMatrix, Modes]"
GetIsoVars                   ::usage = "GetIsoVars[IsoDispModes]"
GetIsoTransformRules         ::usage = "GetIsoDispTransformRules[OpMatrix, IsoDispModes, TranType]"
GetIsoStrainTransformRules   ::usage = "GetIsoStrainTransformRules[GridSymFile]"
Epsilon2Field                ::usage = "Epsilon2Field[strain]"
Field2Epsilon                ::usage = "Field2Epsilon[field]"
GetEpsilonijRule             ::usage = "GetEpsilonijRule[symfile]"
GetInvariants                ::usage = "GetInvariants[seeds, order, AllModes, OpMatrix, GridSymFile]"
ImposeDW                     ::usage = "ImposeDW[Wyckoff0, IsoMatrix, modeset, {Nx, Ny, Nz}]"
ImposeIsoStrainVariedDspMode ::usage = "ImposeIsoStrainVariedDspMode[Wyckoff0, IsoMatrix, modeset, LV]"
SimplifyElementSymbol        ::usage = "SimplifyElementSymbol[ele]"
ShowInvariantTable           ::usage = "ShowInvariantTable[TableData]"
Jij                          ::usage = "Jij[r0, MeshDim]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
mIso
Iso
Epsilon
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)
(*Options[ImportIsodistortCIF]    = {Fractional->False, CorrectLabels->True, Tolerance->10^-6}*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ShowIsoModes[PosVec_] := Module[{StructData},
  StructData = Table[{ElementData[i[[3]], "IconColor"], 
                      Sphere[i[[1]], QuantityMagnitude[ElementData[i[[3]], "AtomicRadius"], "Angstroms"]/2], 
                      Black, 
                      Text[i[[3]]<>ToString@Position[PosVec, i][[1, 1]], 
                      i[[1]]]}, {i, PosVec}];

  ArrowData = Table[{Green, 
                     Arrowheads[0.03], 
                     Arrow@Tube[{i[[1]], 
                     i[[1]] + i[[2]]}]}, {i, PosVec}];

  Print[Graphics3D[{StructData, ArrowData}, 
                   ImageSize -> 500, 
                   Axes -> True, 
                   AxesLabel -> (Style[#, Bold, 64]&/@{"a", "b", "c"}), 
                   ViewPoint -> {0, 0, \[Infinity]}]];
]

ISODISTORT[R0_, pos0_, pos_, IsoMatrix_, label_] := Module[{imode, Amp, NN},
  Amp = Table[Flatten[R0.PbcDiff[#] & /@ (Transpose[pos][[1]] - Transpose[pos0][[1]])].Normalize@Normal[IsoMatrix[[;; , imode]]], {imode, Length@label}];
  NN = Table[1/Norm@Flatten[R0.# & /@ Partition[Normal[IsoMatrix][[;; , imode]], 3]], {imode, Length@label}];
  Return[{Range[Length@label], label, NN, Chop[Amp, 2 10^-4]}\[Transpose]]
]

ImposeIsoMode[Wyckoff0_, IsoMatrix_, modeset_, s_] := Module[{mode, id, Amp, pos},
  pos = Wyckoff0;
  Do[id = mode[[1]]; Amp = s mode[[3]] mode[[4]];
     pos = Transpose[{#[[1]] & /@ pos + Amp If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], First@StringCases[#,RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}], {mode, modeset}];
  Return[pos]
]

Jij[symfile_] := Module[{tij, i}, 
  tij = DeleteDuplicates[DeleteCases[Flatten[Table[Total[Sum[Normalize[#[[1]][[2]]][[i]] Subscript[Iso, i, 0, 0, 0], {i, 3}]*Sum[Normalize[#[[2]][[2]]][[i]] Subscript[Iso, i, #[[2]][[1]][[1]], #[[2]][[1]][[2]], #[[2]][[1]][[3]]], {i, 3}] & /@ SiteJij[symfile, {{{0, 0, 0}, {0, 0, 1}}, {neighbor, vec}}]], {neighbor, Flatten[GetUpTo3rdNeighbors[], 1]}, {vec, IdentityMatrix[3]}]], 0, 2], #1 === -#2 || #1 === #2 &];
  Return[Table[tij[[i]]/Total[DeleteDuplicates[Abs[Flatten[CoefficientList[tij[[i]], Variables[tij[[i]]]]]]]] // Expand, {i, Length@tij}]]]

Epsilon2Field[strain_] := Module[{}, 
  {{{0, 0, 0}, IdentityMatrix[3][[strain[[2]]]]}, {IdentityMatrix[3][[strain[[3]]]], IdentityMatrix[3][[strain[[2]]]]}, {{0, 0, 0}, IdentityMatrix[3][[strain[[3]]]]}, {IdentityMatrix[3][[strain[[2]]]], IdentityMatrix[3][[strain[[3]]]]}}]

Field2Epsilon[field_] := Module[{}, 
  1/2 (Sum[Normalize[field[[2]][[1]]][[j]] Normalize[field[[1]][[2]]][[i]] Subscript[Epsilon, i, j], {j, 3}, {i, 3}] + Sum[Normalize[field[[4]][[1]]][[j]] Normalize[field[[3]][[2]]][[i]] Subscript[Epsilon, i, j], {j, 3}, {i, 3}])/.{Subscript[Epsilon, 2, 1] -> Subscript[Epsilon, 1, 2], Subscript[Epsilon, 3, 1] -> Subscript[Epsilon, 1, 3], Subscript[Epsilon, 3, 2] -> Subscript[Epsilon, 2, 3]}]

(*Epsilon2Field[strain_] := Module[{},
  {{{0, 0, 0}, IdentityMatrix[3][[strain[[2]]]]}, {IdentityMatrix[3][[strain[[3]]]], IdentityMatrix[3][[strain[[2]]]]}}
]

Field2Epsilon[field_] := Module[{},
  Sum[Normalize[field[[2]][[1]]][[j]] Normalize[field[[1]][[2]]][[i]] Subscript[Epsilon, i, j], {j, 3}, {i, 3}] /. {Subscript[Epsilon, 2, 1] -> Subscript[Epsilon, 1, 2], Subscript[Epsilon, 3, 1] -> Subscript[Epsilon, 1, 3], Subscript[Epsilon, 3, 2] -> Subscript[Epsilon, 2, 3]}
]*)

GetEpsilonijRule[symfile_] := Module[{tij, i},
  strains = DeleteDuplicates@Flatten[SparseArray[{{i_, j_} /; i == j -> Subscript[Epsilon, i, j], {i_, j_} /; i < j -> Subscript[Epsilon, i, j], {i_, j_} /; i > j -> Subscript[Epsilon, j, i]}, {3, 3}] // Normal];
  Thread[strains -> #] & /@ (Table[Field2Epsilon[#] & /@ SiteJij[symfile, Epsilon2Field[ep]], {ep, strains}]\[Transpose])
]

GetIsoVars[IsoDispModes_] := Module[{VarString, Var},
  VarString = {#1 & @@@ IsoDispModes, StringReplace[First@StringCases[#2, RegularExpression["[[:upper:]]*\\d[+-]"]], {"+" -> "Plus", "-" -> "Minus"}] & @@@ IsoDispModes, StringPart[#2, -2] & @@@ IsoDispModes}\[Transpose];
  Var = {#1, Subscript[ToExpression[#2], ToExpression[#3]]} & @@@ VarString;
  Return[Var]
]

GetIRvector[id_, IsoMatrix_, pos_] := Module[{IRvector},
  IRvector = {If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], # & /@ (pos\[Transpose][[2]])}\[Transpose];
  Return[IRvector]
]

GetOpMatrix[SymOpFile_, pos_, IsoMatrix_, Modes_, OptionsPattern["spin" -> False]] := Module[{op2, ir1, ir2, proj, mat, OpMatrix, AllIRvectors, AllTransformedIRvectors, IsoDim, t, CifData, CifFlags, xyzName, xyzStrData},
  CifData = Import[SymOpFile, "string"] <> "\n";
  originshift = ToExpression[StringReplace[StringCases[CifData,RegularExpression["origin=\\W\\d\\.*\\d*,\\d\\.*\\d*,\\d\\.*\\d*\\W"]][[1]], {"origin=" -> "", "(" -> "{", ")" -> "}"}]];
  (*---fix for some strange behaviour---*)
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] ->""]];
  CifData = ImportString[CifData, "CIF"];
  CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
 
  xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
  xyzStrData = Flatten[xyzName /. CifData];

  IsoDim = Length@Modes;
  OpMatrix = If[IsoDim != 0,
  AllIRvectors = Flatten[GetIRvector[#1, IsoMatrix, pos]\[Transpose][[1]]] & @@@ Modes;
  AllTransformedIRvectors = Transpose[ParallelTable[SymmetryOpVectorField[SymOpFile, pos, GetIRvector[id, IsoMatrix, pos], "spin" -> OptionValue["spin"]], {id, IsoDim}, DistributedContexts -> {"INVARIANT`ISODISTORT`Private`"}]];
  (*AllTransformedIRvectors = Transpose[Table[SymmetryOpVectorField[SymOpFile, pos, GetIRvector[id, IsoMatrix, pos]], {id, IsoDim}]];*)
  ParallelTable[mat=Table[Flatten[op2[[ir1]]\[Transpose][[1]]].AllIRvectors[[ir2]], {ir1, Range@IsoDim}, {ir2, Range@IsoDim}]; SparseArray[Normalize[#] & /@ mat], {op2, AllTransformedIRvectors}, DistributedContexts -> {"INVARIANT`ISODISTORT`Private`"}],
  Table[{}, {Length@xyzStrData}]
  ];
  Return[OpMatrix]
]

GetIsoTransformRules[OpMatrix_, IsoModes_, TranType_] := Module[{IsoDim, IsoVars, SpinDispRules, rules, i, var},
  IsoDim = Length@IsoModes;
  rules = If[IsoDim != 0,
  IsoVars = Which[TranType == "disp", Subscript[Iso, ToExpression[#1]] &@@@ IsoModes, TranType == "spin", Subscript[mIso, ToExpression[#1]] &@@@ IsoModes];
  SpinDispRules = IsoVars[[#1[[1]]]] -> #2 IsoVars[[#1[[2]]]] & @@@ Drop[ArrayRules[OpMatrix], -1];
  Table[First@DeleteDuplicates[Keys[SpinDispRules[[#]]] & /@ i] -> Total[Values[SpinDispRules[[#]]] & /@ i], {i, Table[Flatten[Position[Keys@SpinDispRules, var]], {var, Which[TranType == "disp", Subscript[Iso, #] & /@ Range[IsoDim], TranType == "spin", Subscript[mIso, #] & /@ Range[IsoDim]]}]}], {}];
  Return[rules]
]

GetIsoStrainTransformRules[GridSymFile_] := Module[{StrainRules},
  StrainRules = GetEpsilonijRule[GridSymFile];
  Return[StrainRules]
]

GetInvariants[seeds_, order_, DispModes_, OpDispMatrix_, SpinModes_, OpSpinMatrix_, GridSymFile_, OptionsPattern[{"strain"->False, "counting"->False}]] := Module[{monomials, invariant, TransformRules, i, ss, strains},
  strains = DeleteDuplicates@Flatten[SparseArray[{{i_, j_} /; i == j -> Subscript[Epsilon, i, j], {i_, j_} /; i < j -> Subscript[Epsilon, i, j], {i_, j_} /; i > j -> Subscript[Epsilon, j, i]}, {3, 3}] // Normal];
  CoeffNorm = Thread[seeds -> ConstantArray[1, Length@seeds]];
  monomials = MonomialList[Total[seeds]^order];
  monomials = If[OptionValue["strain"], Flatten[Table[# ss & /@ monomials, {ss, strains}]], monomials];
  TransformRules = Join[#1, #2, #3] &@@@ ({GetIsoTransformRules[#, DispModes, "disp"] & /@ OpDispMatrix, GetIsoTransformRules[#, SpinModes, "spin"] & /@ OpSpinMatrix, GetIsoStrainTransformRules[GridSymFile]}\[Transpose]);

  invariant = DeleteDuplicates[DeleteCases[Union[Total[(monomials /. # & /@ TransformRules)]], i_/;i==0], (#1 -#2 == 0 || #1 + #2 == 0) &];
  Return[Table[Rationalize[Simplify[invariant[[i]]/If[OptionValue["counting"], Total[Abs[Flatten[CoefficientList[invariant[[i]], Variables[invariant[[i]]]]]]], Total[DeleteDuplicates[Abs[Flatten[CoefficientList[invariant[[i]], Variables[invariant[[i]]]]]]]]]]], {i, Length@invariant}]]
]

ImposeDW[Wyckoff0_, IsoMatrix_, modeset_, {Nx_, Ny_, Nz_}] := Module[{mode, id, Amp, pos, s, ix, iy, iz, Superpos},
  Superpos = Table[{#1 + {ix, iy, iz}, #2} & @@@ Wyckoff0, {ix, 0, Nx - 1}, {iy, 0, Ny - 1}, {iz, 0, Nz - 1}];
  Do[pos = Superpos[[ix]][[iy]][[iz]];
     s = Cos[2 Pi {1/Nx, 1/Ny, 1/Nz}.{ix, iy, iz}];
     Do[id = mode[[1]]; Amp = s mode[[3]] mode[[4]];
        pos = Transpose[{#[[1]] & /@ pos + Amp If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], First@StringCases[#, RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}], {mode, modeset}
        ];
     Superpos[[ix]][[iy]][[iz]] = pos, {ix, Nx}, {iy, Ny}, {iz, Nz}
     ];
  Return[Superpos]
]

ImposeIsoStrainVariedDspMode[Wyckoff0_, IsoMatrix_, modeset_, LV_] := Module[{mode, modesetnew, id, Amp, pos, NN},
  pos = Wyckoff0;
  Do[id = mode[[1]];
     NN = 1/Norm@Flatten[LV.# & /@ Partition[Normal[IsoMatrix][[;; , id]], 3]];
     Amp = NN mode[[4]];
     pos = Transpose[{#[[1]] & /@ pos + Amp If[IntegerQ[id] , # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]],
                      First@StringCases[#, RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}],
   {mode, modeset}];
  Return[pos]
]

SimplifyElementSymbol[ele_] := Module[{SimplifyList},
  SimplifyList = First@StringCases[#, RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ ele;
  Return[SimplifyList]
]

ShowInvariantTable[TableData_, param_, OptionsPattern["FontSize" -> 12]] := Module[{m, n},
  Print[Rotate[Grid[Table[Style[Rotate[# // Expand, -270 Degree], Black, Bold, OptionValue["FontSize"]] & /@ (Flatten[Table[{Flatten[{"param", param[[n]]}], Prepend[TableData[[n]], n]}, {n, Length@TableData}], 1][[m]]), {m, 2 Length@TableData}], Alignment -> Left, Frame -> All], 270 Degree]]]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
