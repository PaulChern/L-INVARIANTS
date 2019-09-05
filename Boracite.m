BeginPackage["INVARIANT`Boracite`",{"INVARIANT`Structure`", "INVARIANT`ISODISTORT`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetPmode           ::usage = "GetPmode[DspMode]"
GetX5mode          ::usage = "GetX5mode[DspMode]"
GetStrain          ::usage = "GetStrain[parent, sub]"
GenTraningSets     ::usage = "GenTraningSets[AllModes, num, PhaseName, type]"
GenStrainSets      ::usage = "GenStrainSets[LatticeVector, type, num]"
GoalFunction       ::usage = "GoalFunctionEnergy[func, ts, coe, volume]"
GetElasticModuli   ::usage = "GetElasticModuli[file, volume]"
ReadInvariant      ::usage = "ReadInvariant[invariants_, \[CapitalGamma]4_]"
GetTensor          ::usage = "GetTensor[IvariantTerm]"
LandauInvariant    ::usage = "LandauInvariant[]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)
InvariantString = Association[
"3rd\[CapitalGamma]4" -> "3 n1 n2 n3",
"3rdX" -> "3   n4 n6 n9 - n5 n7 n8",
"3rd\[CapitalGamma]4X" -> "3   n1 n6 n7 + n2 n4 n5 + n3 n8 n9",
"4thX" -> 
"4   n4^4 + 2 n4^2 n5^2 + 2 n4^2 n6^2 + 2 n4^2 n7^2 + 2 n4^2 n8^2 + 2 n4^2 n9^2 + n5^4 + 2 n5^2 n6^2 + 2 n5^2 n7^2 + 2 n5^2 n8^2 + 2 n5^2 n9^2 + n6^4 + 2 n6^2 n7^2 + 2 n6^2 n8^2 + 2 n6^2 n9^2 + n7^4 + 2 n7^2 n8^2 + 2 n7^2 n9^2 + n8^4 + 2 n8^2 n9^2 + n9^4
4   n4^4 + n5^4 + n6^4 + n7^4 + n8^4 + n9^4
4   n4^2 n5^2 + n6^2 n7^2 + n8^2 n9^2
4   n4^2 n6^2 + n4^2 n9^2 + n5^2 n7^2 + n5^2 n8^2 + n6^2 n9^2 + n7^2 n8^2
4   n4^2 n7^2 + n5^2 n9^2 + n6^2 n8^2",
"4th\[CapitalGamma]4X" -> 
"4   n1^2 n4^2 + n1^2 n5^2 + n1^2 n6^2 + n1^2 n7^2 + n1^2 n8^2 + n1^2 n9^2 + n2^2 n4^2 + n2^2 n5^2 + n2^2 n6^2 + n2^2 n7^2 + n2^2 n8^2 + n2^2 n9^2 + n3^2 n4^2 + n3^2 n5^2 + n3^2 n6^2 + n3^2 n7^2 + n3^2 n8^2 + n3^2 n9^2
4   n1^2 n4^2 + n1^2 n8^2 + n2^2 n7^2 + n2^2 n9^2 + n3^2 n5^2 + n3^2 n6^2
4   n1^2 n5^2 + n1^2 n9^2 + n2^2 n6^2 + n2^2 n8^2 + n3^2 n4^2 + n3^2 n7^2
4   n1 n2 n8 n9 + n1 n3 n4 n5 + n2 n3 n6 n7
4   n1 n4 n7 n9 - n1 n5 n6 n8 - n2 n4 n7 n8 + n2 n5 n6 n9 + n3 n4 n6 n8 - n3 n5 n7 n9",
"6thX" -> 
"6   n4^6 + n5^6 + n6^6 + n7^6 + n8^6 + n9^6
6   n4^4 n5^2 + n4^2 n5^4 + n6^4 n7^2 + n6^2 n7^4 + n8^4 n9^2 + n8^2 n9^4
6   n4^4 n6^2 + n4^2 n9^4 + n5^4 n8^2 + n5^2 n7^4 + n6^4 n9^2 + n7^2 n8^4
6   n4^4 n7^2 + n4^2 n7^4 + n5^4 n9^2 + n5^2 n9^4 + n6^4 n8^2 + n6^2 n8^4
6   n4^4 n8^2 + n4^2 n8^4 + n5^4 n6^2 + n5^2 n6^4 + n7^4 n9^2 + n7^2 n9^4"];
x; y; z;
aa; AA; bb; BB; cc; CC;
n1; n2; n3; n4; n5; n6; n7; n8; n9;
Subscript[\[Epsilon], 1, 1]; Subscript[\[Epsilon], 2, 2]; Subscript[\[Epsilon], 3, 3]; 
Subscript[\[Epsilon], 2, 3]; Subscript[\[Epsilon], 1, 3]; Subscript[\[Epsilon], 1, 2];
Subscript[P, 1]; Subscript[P, 2]; Subscript[P, 3];
Subscript[\[Alpha], xx]; 
Subscript[\[Alpha], xxxx]; Subscript[\[Alpha], xxXX]; Subscript[\[Alpha], xxyy]; Subscript[\[Alpha], xxYY]; 
Subscript[\[Alpha], xxxxxx]; Subscript[\[Alpha], xxxxXX]; Subscript[\[Alpha], xxxxyy]; Subscript[\[Alpha], xxxxYY]; Subscript[\[Alpha], XXXXyy];
Subscript[\[Beta], 11];
Subscript[\[Beta], 1111]; Subscript[\[Beta], 1122]; 
Subscript[\[Beta], 111111]; Subscript[\[Beta], 111122]; Subscript[\[Beta], 112233];
Subscript[C, 1111]; Subscript[C, 1122]; Subscript[C, 1212];
Subscript[Q, 1111]; Subscript[Q, 1122]; Subscript[Q, 1212];
Subscript[G, 1111]; Subscript[G, 1122]; Subscript[G, 1212];
Subscript[\[Gamma], x5];
Subscript[\[Gamma], p];
Subscript[\[Gamma], c];
Subscript[\[Gamma], pc];
Subscript[\[Gamma], xc];
Subscript[\[Gamma], xp];
Subscript[\[Eta], xxpp];
(*--------------------------- Options ----------------------------*)
(*Options[ImportIsodistortCIF]    = {Fractional->False, CorrectLabels->True, Tolerance->10^-6}*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)

GetPmode[DspMode_] := Module[{GM4NoneZero, Pmode, Pindex, Px, Py, Pz},
  Pindex = Partition[{43, 44, 45, 55, 56, 57, 67, 68, 69, 313, 314, 315, 340, 341, 342, 343, 344, 345, 412, 413, 414, 415, 416, 417}, 3]\[Transpose];
  GM4NoneZero = If[StringCases[#[[2]], RegularExpression["GM4"]] != {}, #, ## &[]] & /@ DspMode;
  Px = Sign[DspMode[[340]][[4]]] Norm[(If[MemberQ[Pindex[[1]], #[[1]]], #, ## &[]] & /@ GM4NoneZero)\[Transpose][[4]]];
  Py = Sign[DspMode[[341]][[4]]] Norm[(If[MemberQ[Pindex[[2]], #[[1]]], #, ## &[]] & /@ GM4NoneZero)\[Transpose][[4]]];
  Pz = Sign[DspMode[[342]][[4]]] Norm[(If[MemberQ[Pindex[[3]], #[[1]]], #, ## &[]] & /@ GM4NoneZero)\[Transpose][[4]]];
  Return[{Px, Py, Pz}]
]

GetX5mode[DspMode_] := Module[{X5string, X5a, X5b, X5c, X5d, X5e, X5f},
  X5string = "X5(\\()[abcdef,]+(\\))(\[)[[:upper:]]([[:lower:]]?)\\d+:\\w:dsp(\])\\w+(\\()";
  X5a = Sign[DspMode[[379]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "a" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5b = Sign[DspMode[[380]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "b" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5c = Sign[DspMode[[381]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "c" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5d = Sign[DspMode[[382]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "d" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5e = Sign[DspMode[[383]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "e" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5f = Sign[DspMode[[384]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "f" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  Return[{X5a, X5b, X5c, X5d, X5e, X5f}]
]

GetStrain[parent_, sub_] := Module[{e, strain}, 
  e = N[sub[[6]] /. sub[[8]]].Inverse[N[parent[[6]] /. parent[[8]]]] - IdentityMatrix[3];
  strain = 1/2 (e + Transpose[e] + Transpose[e].e);
  Return[{strain[[1, 1]], strain[[2, 2]], strain[[3, 3]], strain[[1, 2]], strain[[2, 3]], strain[[1, 3]]}]
]

GenTraningSets[AllModes_, num_, PhaseName_, type_] := Module[{i, str, N, ModesNum, DirName, ModeTable, X5NoneZero, GM4NoneZero, Pmode, Pxmode, Pymode, Pzmode, Pindex, X5string, X5a1,X5a2, X5b1, X5b2, X5c1, X5c2},
  Pindex = Partition[{43, 44, 45, 55, 56, 57, 67, 68, 69, 313, 314, 315, 340, 341, 342, 343, 344, 345, 412, 413, 414, 415, 416, 417}, 3]\[Transpose];
  X5NoneZero = If[StringCases[#[[2]], RegularExpression["X5"]] != {} && Abs[#[[4]]] > 0.0001, #, ## &[]] & /@ AllModes;
  X5string = "X5(\\()[abcdef,]+(\\))(\[)[[:upper:]]([[:lower:]]?)\\d+:\\w:dsp(\])\\w+(\\()";
  {X5a1, X5a2, X5b1, X5b2, X5c1, X5c2} = Table[If[StringCases[#[[2]], RegularExpression[X5string <> str <> "(\\))"]] != {}, #, ## &[]] & /@ AllModes, {str, {"a", "b", "c", "d", "e", "f"}}];
  GM4NoneZero = If[StringCases[#[[2]], RegularExpression["GM4"]] != {} && Abs[#[[4]]] > 0.000, #, ## &[]] & /@ AllModes;
  Pmode = If[MemberQ[Flatten[Pindex], #[[1]]], #, ## &[]] & /@ GM4NoneZero;
  {Pxmode, Pymode, Pzmode} = Table[If[MemberQ[Pindex[[str]], #[[1]]], #, ## &[]] & /@ AllModes, {str, 3}];
  DirName = PhaseName <> "-" <> type;
  Print[DirName<>" Generated."];
  N = Quotient[num, 2];
  Which[
   type == "xp",
   ModesNum = Join[Pmode\[Transpose][[1]], X5NoneZero\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], (1 + i/100) #[[4]]}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "x",
   ModesNum = Join[X5NoneZero\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], (1 + i/100) #[[4]]}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p",
   ModesNum = Join[Pmode\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], (1 + i/100) #[[4]]}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p1p2",
   ModesNum = {Pxmode\[Transpose][[1]], Pymode\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]],Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] i/500}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p2",
   ModesNum = Join[Pymode\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], i/500}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p3",
   ModesNum = Join[Pzmode\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], i/500}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "a2b1",
   ModesNum = {X5a2\[Transpose][[1]], X5b1\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]], Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "a1b2",
   ModesNum = {X5a1\[Transpose][[1]], X5b2\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]], Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "c1c2",
   ModesNum = {X5c1\[Transpose][[1]], X5c2\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]], Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "a1b1c2",
   ModesNum = Join[X5a1\[Transpose][[1]], X5b1\[Transpose][[1]], X5c2\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}]
  ]; 
  Return[{ModeTable, DirName}]
]

GenStrainSets[LatticeVector_, PhaseName_, type_, num_] := Module[{i, N, R, DirName, StrainTable},
  DirName = type <> "-" <> PhaseName;
  N = Quotient[num, 2];
  R = Normal[Symmetrize[LatticeVector]];
  Which[
   type == "strain",
   StrainTable = Table[R + DiagonalMatrix[i/1000 (Diagonal@R)], {i, -N, N}],
   type == "Q1111",
   StrainTable = Table[DiagonalMatrix[{1.0, 1.0, 1.0 + i/1000}].R, {i, -N, N}],
   type == "Q1122",
   StrainTable = Table[DiagonalMatrix[{1.0 + i/1000, 1.0 + i/1000, 1.0}].R, {i, -N, N}],
   type == "Q1212",
   StrainTable = Table[{{1.0, i/1000, 0}, {i/1000, 1.0, 0}, {0, 0, 1.0}}.R, {i, -N, N}],
   type == "gammaxc",
   StrainTable = Table[{{1.0, i/1000, 0}, {i/1000, 1.0, 0}, {0, 0, 1.0}}.R, {i, -N, N}]
  ];
  Return[{StrainTable, DirName}]
]

GoalFunction[func_, ts_, coe_, volume_] := Module[{},
  Total[((func /. coe /. Thread[{aa, AA, bb, BB, CC, cc, Subscript[P, 1], Subscript[P, 2], Subscript[P, 3], Subscript[\[Epsilon], 1, 1], Subscript[\[Epsilon], 2, 2], Subscript[\[Epsilon], 3, 3], Subscript[\[Epsilon], 1, 2], Subscript[\[Epsilon], 2, 3], Subscript[\[Epsilon], 1, 3]} -> #[[3]]]) - #[[2]]/volume)^2 & /@ ts]
]

GetElasticModuli[file_, volume_] := Module[{KBar, KBar2meV, Cijkl},
  KBar2meV = 0.1*1000/160.21766208;
  KBar = Select[Import[file, {"Data"}], UnsameQ[#, {}] &];
  Cijkl = SparseArray[{{i_, i_, i_, i_} -> Subscript[C, 1111],
                       {i_, i_, j_, j_} /; i != j -> Subscript[C, 1122],
                       {i_, j_, i_, j_} /; i != j -> Subscript[C, 1212]}, {3, 3, 3, 3}];
  DeleteDuplicates@Flatten[Table[Which[i == j && k == l, Cijkl[[i, j, k, l]] -> KBar[[i, k]] KBar2meV, i == k && j == l && i != j, Cijkl[[i, j, k, l]] -> KBar[[i + 3, i + 3]] KBar2meV, True, ## &[]], {i, 3}, {j, 3}, {k, 3}, {l, 3}]]
]   

ReadInvariant[invariants_, \[CapitalGamma]4_] := Module[{R, X5, ops, IsoVarsX5, IsoVars\[CapitalGamma]4, IsoInvOrder, InvariantTerms},
  R = {x, y, z};
  X5 = {aa, AA, bb, BB, CC, cc};
  IsoVarsX5 = {n4, n5, n6, n7, n8, n9};
  IsoVars\[CapitalGamma]4 = {n1, n2, n3};
  Which[\[CapitalGamma]4 == "ShearStrains", ops = {Subscript[\[Epsilon], 2, 3], Subscript[\[Epsilon], 1, 3], Subscript[\[Epsilon], 1, 2]}, 
        \[CapitalGamma]4 == "Polarizations", ops = {Subscript[P, 1], Subscript[P, 2], Subscript[P, 3]}
       ];
  {IsoInvOrder, InvariantTerms} = ReadList[StringToStream[invariants], {Number, Expression}]\[Transpose];
  InvariantTerms = # /. Thread[IsoVarsX5 -> Through[X5 @@ R]] /. Thread[IsoVars\[CapitalGamma]4 -> Through[ops @@ R]] & /@ InvariantTerms;
  Return[InvariantTerms]
]

GetTensor[IvariantTerm_] := Module[{CoeffRules, X5, tensor, sa, XX},
  X5 = {aa, AA, bb, BB, CC, cc};
  tensor = Table[CoeffRules = CoefficientRules[ter, X5];
    sa = SparseArray[
      Table[Flatten[Join[ConstantArray[Position[Thread[ru[[1]] == 6], True], 3], 
                         ConstantArray[Position[Thread[ru[[1]] == 4], True], 2], 
                         ConstantArray[Position[Thread[ru[[1]] == 2], True], 1]]] -> ru[[2]], {ru, CoeffRules}], 6]; 
    1/2 (sa + Transpose[sa]), {ter, IvariantTerm /. {XX_[x, y, z] -> XX}}];
  Return[tensor]
]

LandauInvariant[] := Module[{Boracite, R, r, i, j, k, l, m, n, ndim, Pk, Pl, Pi2, Pj2, Pk2, Pij, Pkl, Pmn, \[Epsilon]ij, \[Epsilon]kl, \[Epsilon]mn, epsilon, X, X5, X5k, X5l, X5i2, X5j2, X5k2, X5ij, X5kl, X5mn, \[CapitalGamma]4P, \[CapitalGamma]4\[Epsilon], \[Alpha]i2, TensorX5, TensorX56th, \[Alpha]i2\[Alpha]j2, \[Alpha]i2\[Alpha]j2\[Alpha]k2, \[Beta]i2, \[Beta]iijj, \[Beta]iijjkk, Cijkl, Cijmn, Cklmn, Qmnkl, qijkl, Gijkl, EX5termsX5\[CapitalGamma]4P, fx5, fx56, fg, fp, fc, fpc, fxc, fxp},

(**********************************************-----Variables-----************************************************)
  R = {x, y, z};
  ndim = Length@R;
  Subscript[r, i_] := R[[i]];
  Pk = Pl = Array[Subscript[P, #] @@ R &, {3}];
  Pi2 = Pj2 = Pk2 = Array[(Subscript[P, #] @@ R)^2 &, {3}];
  Pij = Pkl = Pmn = Array[D[Subscript[P, #1] @@ R, Subscript[r, #2]] &, {3, ndim}];
  \[Epsilon]ij = \[Epsilon]kl = \[Epsilon]mn = SparseArray[{{i_, j_} /; i == j -> Subscript[\[Epsilon], i, j] @@ R, 
                                                            {i_, j_} /; i < j -> Subscript[\[Epsilon], i, j] @@ R, 
                                                            {i_, j_} /; i > j -> Subscript[\[Epsilon], j, i] @@ R}, {3, 3}];
  epsilon = Table[Subscript[\[Epsilon], i, j] @@ R -> Array[1/2 (D[Subscript[u, #1] @@ R, Subscript[r, #2]] + D[Subscript[u, #2] @@ R, Subscript[r, #1]]) &, {ndim, ndim}][[i, j]], {i, 1, ndim}, {j, 1, ndim}] // Flatten;
  X5 = {aa, AA, bb, BB, CC, cc};
  X5k = X5l = Array[X5[[#]] @@ R &, {6}];
  X5i2 = X5j2 = X5k2 = Array[(X5[[#]] @@ R)^2 &, {6}];
  X5ij = X5kl = X5mn = Array[D[X5[[#1]] @@ R, Subscript[r, #2]] &, {6, ndim}];
  X = Through[X5 @@ R];
  \[CapitalGamma]4P = {Subscript[P, 1], Subscript[P, 2], Subscript[P, 3]};
  \[CapitalGamma]4\[Epsilon] = {Subscript[\[Epsilon], 2, 3], Subscript[\[Epsilon], 1, 3], Subscript[\[Epsilon], 1, 2]};

(**********************************************-----Tensors-----************************************************)
  TensorX5 = GetTensor[ReadInvariant[InvariantString["4thX"], "Polarizations"]];
  TensorX56th = GetTensor[ReadInvariant[InvariantString["6thX"], "Polarizations"]];
  EX5termsX5\[CapitalGamma]4P = ReadInvariant[InvariantString["6thX"], "Polarizations"];
  \[Alpha]i2 = SparseArray[{{i_} -> Subscript[\[Alpha], xx]}, {6}];
  \[Alpha]i2\[Alpha]j2 = Total[Normal[#1*#2] & @@ {TensorX5, {2 Subscript[\[Alpha], xxxx], 0, 2 Subscript[\[Alpha], xxXX], 2 Subscript[\[Alpha], xxyy], 2 Subscript[\[Alpha], xxYY]}}];
  \[Alpha]i2\[Alpha]j2\[Alpha]k2 = Total[Normal[#1*#2] & @@ {TensorX56th, {Subscript[\[Alpha], xxxxxx], Subscript[\[Alpha], xxxxXX], Subscript[\[Alpha], xxxxyy], Subscript[\[Alpha], xxxxYY], Subscript[\[Alpha], XXXXyy]}}];
  \[Beta]i2 = SparseArray[{{i_} -> Subscript[\[Beta], 11]}, {3}];
  \[Beta]iijj = SparseArray[{{i_, i_} -> 2 Subscript[\[Beta], 1111], {i_, j_} /; i != j -> Subscript[\[Beta], 1122]}, {3, 3}];
  \[Beta]iijjkk = SparseArray[{{i_, i_, i_} -> 2 Subscript[\[Beta], 111111], 
                               {i_, i_, j_} /; i != j -> 2/3 Subscript[\[Beta], 111122], 
                               {i_, j_, i_} /; i != j -> 2/3 Subscript[\[Beta], 111122], 
                               {j_, i_, i_} /; i != j -> 2/3 Subscript[\[Beta], 111122], 
                               {i_, j_, k_} /; i != j && i != k && j != k -> 1/3 Subscript[\[Beta], 112233]}, {3, 3, 3}];
  Cijkl = Cijmn = Cklmn = SparseArray[{{i_, i_, i_, i_} -> Subscript[C, 1111], 
                                       {i_, i_, j_, j_} /; i != j -> Subscript[C, 1122], 
                                       {i_, j_, i_, j_} /; i != j -> Subscript[C, 1212]}, {3, 3, 3, 3}];
  Qmnkl = SparseArray[{{m_, m_, m_, m_} -> Subscript[Q, 1111], 
                       {m_, m_, n_, n_} /; m != n -> Subscript[Q, 1122], 
                       {m_, n_, m_, n_} /; m != n -> Subscript[Q, 1212]}, {3, 3, 3, 3}];
  qijkl = TensorContract[2 Cijmn\[TensorProduct]Qmnkl, {{3, 5}, {4, 6}}];
  Gijkl = SparseArray[{{i_, i_, i_, i_} -> Subscript[G, 1111], 
                       {i_, i_, j_, j_} /; i != j -> Subscript[G, 1122], 
                       {i_, j_, i_, j_} /; i != j -> Subscript[G, 1212]}, {6, 3,6, 3}];

(**********************************************-----Energy-----************************************************)
    fx5 = TensorContract[\[Alpha]i2\[TensorProduct]X5i2, {1, 2}] 
        + TensorContract[1/2 \[Alpha]i2\[Alpha]j2\[TensorProduct]X5i2\[TensorProduct]X5j2,{{1, 3}, {2, 4}}] 
        + 0 TensorContract[\[Alpha]i2\[Alpha]j2\[Alpha]k2\[TensorProduct]X5i2\[TensorProduct]X5j2\[TensorProduct]X5k2, {{1, 4}, {2, 5}, {3, 6}}] 
        + Subscript[\[Gamma], x5] (aa[x, y, z] bb[x, y, z] cc[x, y, z] - AA[x, y, z] BB[x, y, z] CC[x, y, z]); 
  fg = TensorContract[1/2 X5ij\[TensorProduct]Gijkl\[TensorProduct]X5kl, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}];
  fp = TensorContract[1/2 \[Beta]i2\[TensorProduct]Pi2, {{1, 2}}] 
     + TensorContract[1/2 \[Beta]iijj\[TensorProduct]Pi2\[TensorProduct]Pj2, {{1, 3}, {2, 4}}] 
     + 0 TensorContract[1/2 \[Beta]iijjkk\[TensorProduct]Pi2\[TensorProduct]Pj2\[TensorProduct]Pk2, {{1, 4}, {2, 5}, {3, 6}}] 
     + Subscript[\[Gamma], p] Subscript[P, 1][x, y, z] Subscript[P, 2][x, y, z] Subscript[P, 3][x, y, z];
  fx56 = Subscript[\[Alpha], xxxxXX] (aa[x, y, z]^2 AA[x, y, z]^2 (aa[x, y, z]^2 + AA[x, y, z]^2) + bb[x, y, z]^2 BB[x, y, z]^2 (bb[x, y, z]^2 + BB[x, y, z]^2) + CC[x, y, z]^2 cc[x, y, z]^2 (CC[x, y, z]^2 + cc[x, y, z]^2)) + Subscript[\[Alpha], xxxxxx] ((aa[x, y, z]^2 + AA[x, y, z]^2)^3 + (bb[x, y, z]^2 + BB[x, y, z]^2)^3 + (CC[x, y, z]^2 + cc[x, y, z]^2)^3);
  fc = TensorContract[1/2 \[Epsilon]ij\[TensorProduct]Cijkl\[TensorProduct]\[Epsilon]kl, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}] 
     + Subscript[\[Gamma], c] Subscript[\[Epsilon], 1, 2][x, y, z] Subscript[\[Epsilon], 2, 3][x, y, z] Subscript[\[Epsilon], 1, 3][x, y, z];
  fpc = TensorContract[-(1/2) \[Epsilon]ij\[TensorProduct]qijkl\[TensorProduct]Pk\[TensorProduct]Pl, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}] + Subscript[\[Gamma], pc] {Subscript[\[Epsilon], 2, 3][x, y, z], Subscript[\[Epsilon], 1, 3][x, y, z], Subscript[\[Epsilon], 1, 2][x, y, z]}.{Subscript[P, 1][x, y, z], Subscript[P, 2][x, y, z], Subscript[P, 3][x, y, z]};
  fxc = Subscript[\[Gamma], xc] (Subscript[\[Epsilon], 2, 3][x, y, z] bb BB + Subscript[\[Epsilon], 1, 3][x, y, z] aa AA + Subscript[\[Epsilon], 1, 2][x, y, z] cc CC) /. Thread[X5 -> X];
  fxp = (Subscript[\[Gamma], xp] (Subscript[P, 1][x, y, z] bb BB + Subscript[P, 2][x, y, z] aa AA + Subscript[P, 3][x, y, z] cc CC) /. Thread[X5 -> X]) 
      + EX5termsX5\[CapitalGamma]4P.{Subscript[\[Eta], xxpp], 0, 0, 0, 0};
  Boracite = fx5 + fp + fc + fpc + fxc + fxp;
  Return[Boracite]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
