BeginPackage["INVARIANT`ISODISTORT`",{"INVARIANT`Structure`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ShowIsoModes                 ::usage = "ShowIsoModes[PosVec]"
ISODISTORT                   ::usage = "ISODISTORT[R0, pos0, pos, IsoMatrix, label]"
ImposeIsoMode                ::usage = "ImposeIsoMode[Wyckoff0, IsoMatrix, modeset, s]"
GetIRvector                  ::usage = "GetIRvector[id, pos]"
GetOpMatrix                  ::usage = "GetOpMatrix[SymOpFile, pos, IsoMatrix, Modes]"
GetIsoVars                   ::usage = "GetIsoVars[IsoDispModes]"
GetIsoTransformRules         ::usage = "GetIsoTransformRules[OpMatrix, IsoDispModes]"
GetInvariants                ::usage = "GetInvariants[seeds, order, AllModes, OpMatrix]"
ImposeDW                     ::usage = "ImposeDW[Wyckoff0, IsoMatrix, modeset, {Nx, Ny, Nz}]"
ImposeIsoStrainVariedDspMode ::usage = "ImposeIsoStrainVariedDspMode[Wyckoff0, IsoMatrix, modeset, LV]"
SimplifyElementSymbol        ::usage = "SimplifyElementSymbol[ele]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
Iso

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
     pos = Transpose[{#[[1]] & /@ pos + Amp If[IntegerQ[id] , # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], 
                      First@StringCases[#,RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}], 
     {mode, modeset}];
  Return[pos]
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

GetOpMatrix[SymOpFile_, pos_, IsoMatrix_, Modes_] := Module[{op2, ir1, ir2, proj, mat, OpMatrix, AllIRvectors, AllTransformedIRvectors, IsoDim, t},
  IsoDim = Length@Modes;
  AllIRvectors = Flatten[GetIRvector[#1, IsoMatrix, pos]\[Transpose][[1]]] & @@@ Modes;
  AllTransformedIRvectors = Transpose[ParallelTable[SymmetryOpVectorField[SymOpFile, pos, GetIRvector[id, IsoMatrix, pos]], {id, IsoDim}, DistributedContexts -> {"INVARIANT`ISODISTORT`Private`"}]];
  (*AllTransformedIRvectors = Transpose[Table[SymmetryOpVectorField[SymOpFile, pos, GetIRvector[id, IsoMatrix, pos]], {id, IsoDim}]];*)
  OpMatrix = ParallelTable[mat=Table[Flatten[op2[[ir1]]\[Transpose][[1]]].AllIRvectors[[ir2]], {ir1, Range@IsoDim}, {ir2, Range@IsoDim}];
  SparseArray[Normalize[#] & /@ mat], {op2, AllTransformedIRvectors}, DistributedContexts -> {"INVARIANT`ISODISTORT`Private`"}];
  Return[OpMatrix]
]

GetIsoTransformRules[OpMatrix_, IsoDispModes_] := Module[{IsoDim, IsoVars, rules, i, var},
  IsoDim = Length@IsoDispModes;
  IsoVars = Subscript[Iso, ToExpression[#1]] &@@@ IsoDispModes;
  rules = IsoVars[[#1[[1]]]] -> #2 IsoVars[[#1[[2]]]] & @@@ Drop[ArrayRules[OpMatrix], -1];
  rules = Table[First@DeleteDuplicates[Keys[rules[[#]]] & /@ i] -> Total[Values[rules[[#]]] & /@ i], {i, Table[Flatten[Position[Keys@rules, var]], {var, Subscript[Iso, #] & /@ Range[IsoDim]}]}];
  Return[rules]
]

GetInvariants[seeds_, order_, AllModes_, OpMatrix_, OptionsPattern[{"ctol" -> 10^-6}]] := Module[{monomials, invariant, i},
  CoeffNorm = Thread[seeds -> ConstantArray[1, Length@seeds]];
  monomials = MonomialList[Total[seeds]^order];
  invariant = DeleteDuplicates[DeleteCases[Union[Total[(monomials /. GetIsoTransformRules[#, AllModes]) & /@ OpMatrix]], 0.], (#1 -#2 == 0 || #1 + #2 == 0) &];
  Return[Table[Simplify[Rationalize[invariant[[i]]/Total[Abs[Flatten[CoefficientList[invariant[[i]], Variables[invariant[[i]]]]]]], OptionValue["ctol"]]], {i, Length@invariant}]]
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
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
