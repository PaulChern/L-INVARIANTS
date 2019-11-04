(* ::Package:: *)

BeginPackage["INVARIANT`Phonopy`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ReadQpoint      ::usage "ReadQpoint[filename]"
ReadBands       ::usage "ReadBands[filename]"
ModesSum        ::usage "ModesSum[list, Modes]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ReadQpoint[filename_] := Module[{QpointPhonon, EigenValuesVectors},
    QpointPhonon = Import[filename, "YAML"];
    EigenValuesVectors = QpointPhonon[[4]][[2]][[1]][[2]][[2]];
    Return[{"reciprocal_lattice"/.QpointPhonon,Table[{"frequency", "eigenvector", "group_velocity"} /. EigenValuesVectors[[i]], {i, Length@EigenValuesVectors}]}];
]

ReadBands[filename_] := Module[{PhonopyBands},
  PhonopyBands = Import[filename, "YAML"];
  Return[Table[{"distance" /. ph, "q-position" /. ph, #2 & @@@ Flatten["band" /. ph]}, {ph, PhonopyBands[[10]][[2]]}]];
]

ModesSum[list_, Modes_] := Module[{modes, i},
  Return[Sum[i[[2]] #\[Transpose][[1]] & /@ Modes[[i[[1]]]][[2]], {i, list}]];
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
