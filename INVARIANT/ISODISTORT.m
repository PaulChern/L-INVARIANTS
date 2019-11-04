BeginPackage["INVARIANT`ISODISTORT`",{"INVARIANT`Structure`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ISODISTORT                   ::usage = "ISODISTORT[R0, pos0, pos, IsoMatrix, label]"
ImposeIsoMode                ::usage = "ImposeIsoMode[Wyckoff0, IsoMatrix, modeset, s]"
ImposeIsoStrainVariedDspMode ::usage = "ImposeIsoStrainVariedDspMode[Wyckoff0, IsoMatrix, modeset, LV]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)
(*Options[ImportIsodistortCIF]    = {Fractional->False, CorrectLabels->True, Tolerance->10^-6}*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)

ISODISTORT[R0_, pos0_, pos_, IsoMatrix_, label_] := Module[{imode, Amp, NN},
  Amp = Table[Chop[Flatten[R0.PbcDiff[#] & /@ (Transpose[pos][[1]] - Transpose[pos0][[1]])].Normalize@Normal[IsoMatrix[[;; , imode]]]], {imode, Length@label}];
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

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
