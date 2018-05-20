(* ::Package:: *)

(*:Title: Semi-simple groups: roots and weights*)

(*:Author: Alexey A. Vladimirov, vladimirov.aleksey@googlemail.com *)

(*:Summary: This package provides
	1) The roots&weights systems for any semi-simple group (using the defentions of Bourbaki)
	2) The constructor of weight diagram for irredusable representations
	3) The destructor of representations with respect to given root system, including C-G decomposition of products.
	4) Various methods to plot multiplets
	5) Many others usefull tools.
*)

(*:Package Version: 0.96\[Beta] *)

(*:Mathetica version: 7.0 *)

(*:History: 
	15.01.2011: Package redesined and partially rewriten (from version 0.91), in more unified way and compact form.
	Additionally changed the algorithm of RepDecomposition and CGDecomposition onto more effective (much more ~4n times!!).
	6.02.2011: Some commands (BuildRep, CGDecomposition) complitelly rewritten, now they hardly use the GT-patterns (for work with them the `private` commands added : REPRESNTATION, GTs, CheckGTPattern
	Added command Weight, many commands obtain insertion of two new options InputMethod, OutputMethod. Added NormalizationMethod->"Canonical"
	Fixed mistake in SetGroup with D[3] definition.
	19.02.2011: The command CGCoefficints added, also the variables CGC, RulesForReps, RulesForStates, RulesForStates1, RulesForStates2 added and protected.
	The option OnlyRep added. 
	Begining of open \[Beta]-testing.
	6.05.2011: The mistake in CGCoefficients is corrected (was wrong chousing of equivalent representations)
*)


BeginPackage["SimpleGroups`"];


Unprotect[SetGroup,CartanMatrix,FundamentalWeight,RootSystem,Representation,RepProduct,DecomposeRep,PhysicalNormalizationMatrix,DynkinDiagram,CoxeterPlane,BuildRep,CGDecomposition,PhysicalNormalization,RepDimension,DynkinLabels,PlotRep2D,PlotRep3D,Weight];
Clear[SetGroup,CartanMatrix,FundamentalWeight,RootSystem,Representation,RepProduct,DecomposeRep,PhysicalNormalizationMatrix,DynkinDiagram,CoxeterPlane,BuildRep,CGDecomposition,PhysicalNormalization,RepDimension,DynkinLabels,PlotRep2D,PlotRep3D,Weight];
Unprotect[SimpleRoots,PositiveRoots,FundamentalWeights,MaximalRoots,Rank,CoxeterNumber,Group];
Clear[SimpleRoots,PositiveRoots,FundamentalWeights,MaximalRoots,Rank,CoxeterNumber,Group];
Unprotect[CGC,RulesForReps,RulesForStates,RulesForStates2,RulesForStates1];
Clear[CGC,RulesForReps,RulesForStates,RulesForStates2,RulesForStates1];


SetGroup::usage="Sets the global constants for the given group.";

CartanMatrix::usage="Returns the Cartan matrix of the root system (the argument is the system of simple roots).";
FundamentalWeight::usage="Returns fundamental weights for the given (simple)root system.";
RootSystem::usage="Builds the complete root system from the given set of simple roots ([[1]]  -- positive roots, [[2]] -- negative roots)";
Representation::usage="Returns the weight diagram for the group with given set of simple roots. The representation specified by Dynkin labels(default) of eldest weight(RepNotation->EldestWeight)";
RepProduct::usage="Returns the direct (Kronecker)product of two multiplets";
DecomposeRep::usage="Decompose given multiplet with respect to the given roots. Returns the list of Dynkin indices for eldest weights of sub-multiplets.";
PhysicalNormalizationMatrix::usage="Returns the weight-space rotation-rescaling matrix, which rotates the weight space to the physical normalized. It takes NormalizationMethod as an option";

DynkinDiagram::usage="Plots the Dynkin diagram for the set group. The numbers corresponds to the position of the roots in root-lists.";
CoxeterPlane::usage="Returns the vectors which forms the Coxeter plane for the simple groups of rank higher then 2";
BuildRep::usage="Returns the weight diagram of the representation specified by Dynkin labels(default) of eldest weight(RepNotation->EldestWeight)";
CGDecomposition::usage="Decompose the direct product of representations over irreducible representations. Returns the list of Dynkin indices for eldest weights of irreducible representations.";
PhysicalNormalization::usage="Transform the weight (list of weights or representation) from the working group weight-space to the physically normalized weight space. It takes NormalizationMethod as an option.";
RepDimension::usage="Returns the dimension of representation specified by Dynkin labels or eldest weight (for the last case use option RepNotation->EldestWeight).";
DynkinLabels::usage="Returns the Dynkin labels for the element of weight space.";
Weight::usage="Transforms the element of the package to the weight."
CGCoefficients::usage="Returns the Clebsh-Gordan coefficients for SU(N) groups."
PlotRep2D::usage="Plots the weigth diagram of representation, for group with rank 2 (2D).
It has the following options {PointStyle,LinkStyle,MultiplicitiesStyle,AxesStyle,ShowMultiplicities,ShowLinks,LinksList,ShowAxes}.
By defenition the Links to show are PositiveRoots of working group.";
PlotRep3D::usage="Plots the weigth diagram of representation, for group with rank 3 (3D).
It has the following options {PointStyle,LinkStyle,MultiplicitiesStyle,AxesStyle,ShowMultiplicities,ShowLinks,LinksList,ShowAxes}.
By defenition the Links to show are PositiveRoots of working group.";


SimpleRoots::usage="Simple roots of working group";
PositiveRoots::usage="Positive roots of working group";
MaximalRoots::usage="Maximal roots of the working group (for every group in direct product)";
FundamentalWeights::usage="Fundamental weights of working group";
WeylVector::usage="Weyl weight of working group";
Rank::usage="Rank of working group";
CoxeterNumber::usage="Coxeter number of the working group";

CGC::usage="The matrix of Clebsch-Gordan coefficients.";
RulesForReps::usage="The list of substitusions for the representation enumeration, in CGC";
RulesForStates::usage="The list of substitusions for the resulted states enumeration, in CGC";
RulesForStates2::usage="The list of substitusions for the first representation states enumeration, in CGC";
RulesForStates1::usage="The list of substitusions for the second representation states enumeration, in CGC";

Protect[A,B,F,G,SU,SP,SO];
Protect[RepNotation,InputMethod,OutputMethod,NormalizationMethod,PointStyle,LinkStyle,MultiplicitiesStyle,AxesStyle,ShowMultiplicities,ShowLinks,LinksList,ShowAxes,OnlyRep];


Begin["`Private`"]


Group::usage="The list of working group";
PackageElemType::usage="Returns the type of element; element of weight space, list of weights, multiplet or unknow.";
CheckGTPattern::usage="Checks the argument to be a GT-patern of the set group.";
REPRESENTATION::usage="Returns the multiplet (in dynkin indices) for the given group g, and dynI dynkin indices of the eldest weight. Use the special optimized algorithms for every non-exceptional group."
GTs::usage="Returns the GTpatterns for the representation specified by the dy for group g";


RootsWeightsA::usage="Canonical realization of roots for A-algebra.";
RootsWeightsB::usage="Canonical realization of roots for B-algebra.";
RootsWeightsC::usage="Canonical realization of roots for C-algebra.";
RootsWeightsD::usage="Canonical realization of roots for D-algebra.";
RootsWeightsE::usage="Canonical realization of roots for E-algebra.";
RootsWeightsF::usage="Canonical realization of roots for F-algebra.";
RootsWeightsG::usage="Canonical realization of roots for G-algebra.";


(* ::Section::Closed:: *)
(*Root and weights systems for standard type of groups*)


(* ::Subsubsection::Closed:: *)
(*A[l] roots & Weights*)


RootsWeightsA[rank_]:=Module[{\[Epsilon]=IdentityMatrix[rank+1]},
{Table[\[Epsilon][[i]]-\[Epsilon][[i+1]],{i,Length[\[Epsilon]]-1}],Flatten[Table[Table[\[Epsilon][[i]]-\[Epsilon][[j]],{i,j-1}],{j,Length[\[Epsilon]]}],1],Table[Sum[\[Epsilon][[j]],{j,i}]-i/(rank+1) Total[\[Epsilon]],{i,rank}],Join[{1},ConstantArray[0,rank-1],{-1}]}]


(* ::Subsubsection::Closed:: *)
(*B[l] roots & Weights*)


RootsWeightsB[rank_]:=Module[{\[Epsilon]=IdentityMatrix[rank]},
{Join[Table[\[Epsilon][[i]]-\[Epsilon][[i+1]],{i,Length[\[Epsilon]]-1}],{Last[\[Epsilon]]}],Join[Flatten[Table[Table[\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,Length[\[Epsilon]]}],1],Flatten[Table[Table[\[Epsilon][[i]]-\[Epsilon][[j]],{i,j-1}],{j,Length[\[Epsilon]]}],1],\[Epsilon]],Join[Table[Sum[\[Epsilon][[j]],{j,i}],{i,rank-1}],{1/2 Total[\[Epsilon]]}],Join[{1,1},ConstantArray[0,rank-2]]}]


(* ::Subsubsection::Closed:: *)
(*C[l] roots & weigths*)


RootsWeightsC[rank_]:=Module[{\[Epsilon]=IdentityMatrix[rank]},
{Join[Table[\[Epsilon][[i]]-\[Epsilon][[i+1]],{i,Length[\[Epsilon]]-1}],{2 Last[\[Epsilon]]}],Join[Flatten[Table[Table[\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,Length[\[Epsilon]]}],1],Flatten[Table[Table[\[Epsilon][[i]]-\[Epsilon][[j]],{i,j-1}],{j,Length[\[Epsilon]]}],1],2 \[Epsilon]],Table[Sum[\[Epsilon][[j]],{j,i}],{i,rank}],Join[{2},ConstantArray[0,rank-1]]}]


(* ::Subsubsection::Closed:: *)
(*D[l] roots & weigths*)


RootsWeightsD[rank_]:=Module[{\[Epsilon]=IdentityMatrix[rank]},
{Join[Table[\[Epsilon][[i]]-\[Epsilon][[i+1]],{i,Length[\[Epsilon]]-1}],{Last[\[Epsilon]]+\[Epsilon][[-2]]}],Join[Flatten[Table[Table[\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,Length[\[Epsilon]]}],1],Flatten[Table[Table[\[Epsilon][[i]]-\[Epsilon][[j]],{i,j-1}],{j,Length[\[Epsilon]]}],1]],Join[Table[Sum[\[Epsilon][[j]],{j,i}],{i,rank-2}],{1/2 Total[\[Epsilon]]-\[Epsilon][[rank]]},{1/2 Total[\[Epsilon]]}],Join[{1,1},ConstantArray[0,rank-2]]}]


(* ::Subsubsection::Closed:: *)
(*E[6, 7, 8] roots & weigths*)


RootsWeightsE[6]:=Module[{\[Epsilon]=IdentityMatrix[8],permut=Cases[Tuples[{0,1},5],x_/;EvenQ@Total[x]]},
{{1/2 (\[Epsilon][[1]]+\[Epsilon][[8]])-1/2 Sum[\[Epsilon][[i]],{i,2,7}],\[Epsilon][[1]]+\[Epsilon][[2]],\[Epsilon][[2]]-\[Epsilon][[1]],\[Epsilon][[3]]-\[Epsilon][[2]],\[Epsilon][[4]]-\[Epsilon][[3]],\[Epsilon][[5]]-\[Epsilon][[4]]},Join[Flatten[Table[Table[\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,5}],1],Flatten[Table[Table[-\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,5}],1],Table[1/2 (\[Epsilon][[8]]-\[Epsilon][[7]]-\[Epsilon][[6]]+Sum[(-1)^permut[[i,j]] \[Epsilon][[j]],{j,5}]),{i,Length[permut]}]],{2/3 (\[Epsilon][[8]]-\[Epsilon][[7]]-\[Epsilon][[6]]),1/2 Total[\[Epsilon]]-\[Epsilon][[6]]-\[Epsilon][[7]],5/6 (\[Epsilon][[8]]-\[Epsilon][[7]]-\[Epsilon][[6]])+1/2 (-\[Epsilon][[1]]+\[Epsilon][[2]]+\[Epsilon][[3]]+\[Epsilon][[4]]+\[Epsilon][[5]]),\[Epsilon][[3]]+\[Epsilon][[4]]+\[Epsilon][[5]]-\[Epsilon][[6]]-\[Epsilon][[7]]+\[Epsilon][[8]],2/3 (\[Epsilon][[8]]-\[Epsilon][[7]]-\[Epsilon][[6]])+\[Epsilon][[4]]+\[Epsilon][[5]],1/3 (\[Epsilon][[8]]-\[Epsilon][[7]]-\[Epsilon][[6]])+\[Epsilon][[5]]},1/2 {1,1,1,1,1,-1,-1,1}}]
RootsWeightsE[7]:=Module[{\[Epsilon]=IdentityMatrix[8],permut=Cases[Tuples[{0,1},6],x_/;OddQ@Total[x]]},
{{1/2 (\[Epsilon][[1]]+\[Epsilon][[8]])-1/2 Sum[\[Epsilon][[i]],{i,2,7}],\[Epsilon][[1]]+\[Epsilon][[2]],\[Epsilon][[2]]-\[Epsilon][[1]],\[Epsilon][[3]]-\[Epsilon][[2]],\[Epsilon][[4]]-\[Epsilon][[3]],\[Epsilon][[5]]-\[Epsilon][[4]],\[Epsilon][[6]]-\[Epsilon][[5]]},Join[Flatten[Table[Table[\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,6}],1],Flatten[Table[Table[-\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,6}],1],{\[Epsilon][[8]]-\[Epsilon][[7]]},Table[1/2 (\[Epsilon][[8]]-\[Epsilon][[7]]+Sum[(-1)^permut[[i,j]] \[Epsilon][[j]],{j,6}]),{i,Length[permut]}]],{(\[Epsilon][[8]]-\[Epsilon][[7]]),1/2 (\[Epsilon][[1]]+\[Epsilon][[2]]+\[Epsilon][[3]]+\[Epsilon][[4]]+\[Epsilon][[5]]+\[Epsilon][[6]]+2\[Epsilon][[8]]-2\[Epsilon][[7]]),1/2 (-\[Epsilon][[1]]+\[Epsilon][[2]]+\[Epsilon][[3]]+\[Epsilon][[4]]+\[Epsilon][[5]]+\[Epsilon][[6]]+3\[Epsilon][[8]]-3\[Epsilon][[7]]),(\[Epsilon][[3]]+\[Epsilon][[4]]+\[Epsilon][[5]]+\[Epsilon][[6]]+2\[Epsilon][[8]]-2\[Epsilon][[7]]),1/2 (2\[Epsilon][[4]]+2\[Epsilon][[5]]+2\[Epsilon][[6]]+3\[Epsilon][[8]]-3\[Epsilon][[7]]),\[Epsilon][[5]]+\[Epsilon][[6]]-\[Epsilon][[7]]+\[Epsilon][[8]],\[Epsilon][[6]]+1/2 (\[Epsilon][[8]]-\[Epsilon][[7]])},{0,0,0,0,0,0,-1,1}}]
RootsWeightsE[8]:=Module[{\[Epsilon]=IdentityMatrix[8],permut=Cases[Tuples[{0,1},7],x_/;EvenQ@Total[x]]},
{{1/2 (\[Epsilon][[1]]+\[Epsilon][[8]])-1/2 Sum[\[Epsilon][[i]],{i,2,7}],\[Epsilon][[1]]+\[Epsilon][[2]],\[Epsilon][[2]]-\[Epsilon][[1]],\[Epsilon][[3]]-\[Epsilon][[2]],\[Epsilon][[4]]-\[Epsilon][[3]],\[Epsilon][[5]]-\[Epsilon][[4]],\[Epsilon][[6]]-\[Epsilon][[5]],\[Epsilon][[7]]-\[Epsilon][[6]]},Join[Flatten[Table[Table[\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,8}],1],Flatten[Table[Table[-\[Epsilon][[i]]+\[Epsilon][[j]],{i,j-1}],{j,8}],1],Table[1/2 (\[Epsilon][[8]]+Sum[(-1)^permut[[i,j]] \[Epsilon][[j]],{j,7}]),{i,Length[permut]}]],{2\[Epsilon][[8]],1/2 Total[\[Epsilon]]+2\[Epsilon][[8]],1/2 Total[\[Epsilon]]-\[Epsilon][[1]]+3\[Epsilon][[8]],\[Epsilon][[3]]+\[Epsilon][[4]]+\[Epsilon][[5]]+\[Epsilon][[6]]+\[Epsilon][[7]]+5\[Epsilon][[8]],\[Epsilon][[4]]+\[Epsilon][[5]]+\[Epsilon][[6]]+\[Epsilon][[7]]+4\[Epsilon][[8]],\[Epsilon][[5]]+\[Epsilon][[6]]+\[Epsilon][[7]]+3\[Epsilon][[8]],\[Epsilon][[6]]+\[Epsilon][[7]]+2\[Epsilon][[8]],\[Epsilon][[7]]+\[Epsilon][[8]]},{0,0,0,0,0,0,1,1}}]


(* ::Subsubsection::Closed:: *)
(*F[4] roots & weigths*)


RootsWeightsF[4]:={{{0,1,-1,0},{0,0,1,-1},{0,0,0,1},{1/2,-(1/2),-(1/2),-(1/2)}},{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,1,0,0},{1,0,1,0},{0,1,1,0},{1,0,0,1},{0,1,0,1},{0,0,1,1},{1,-1,0,0},{1,0,-1,0},{0,1,-1,0},{1,0,0,-1},{0,1,0,-1},{0,0,1,-1},{1/2,1/2,1/2,1/2},{1/2,1/2,1/2,-(1/2)},{1/2,1/2,-(1/2),1/2},{1/2,1/2,-(1/2),-(1/2)},{1/2,-(1/2),1/2,1/2},{1/2,-(1/2),1/2,-(1/2)},{1/2,-(1/2),-(1/2),1/2},{1/2,-(1/2),-(1/2),-(1/2)}},{{1,1,0,0},{2,1,1,0},{3/2,1/2,1/2,1/2},{1,0,0,0}},{1,1,0,0}}


(* ::Subsubsection::Closed:: *)
(*G[2] roots & weigths*)


RootsWeightsG[2]:={{{1,-1,0},{-2,1,1}},{{1,-1,0},{-2,1,1},{-1,0,1},{0,-1,1},{1,-2,1},{-1,-1,2}},{{0,-1,1},{-1,-1,2}},{-1,-1,2}}


(* ::Section::Closed:: *)
(*SetGroup*)


SetAttributes[SetGroup,HoldAll];
SetGroup::nonPositiveIntRank="The rank `1` should be positive integer number.";
SetGroup::exceptionalGroupRank="The group `1` is an exceptional group. Its rank can be only `2`.";
SetGroup::u1Group="Rank of `1`-group should be greater then 1";
SetGroup::unknownName="The group `1` is unknow one. Please, use A, B, C, D, E, F, SU, SO or SP notations.";
SetGroup[group_]:=Module[{sGroup=Flatten@ReleaseHold[{Hold[group]}//.{Times->List,A->"A",B->"B",C->"C",D->"D",E->"E",F->"F",G->"G",SU->"SU",SO->"SO",SP->"SP"}],readG,rootWeight,tsDim,directProd},
(*....................................Addition functions.....................................................................*)
readG:=Module[{h=Head[#],r=First[#]},
If[Not[Positive[r]&&IntegerQ[r]],Message[SetGroup::nonPositiveIntRank,r];Abort[]];
Switch[Head[#],"A",{h,r},"B",If[r>=2,{h,r},{"A",r}],"C",If[r>=2,{h,r},{"A",r}],"D",If[r>3,{h,r},If[r==3,{"A",r},If[r==2,{{"A",1},{"A",1}},{"A",1}]]],"E",If[8>=r>=6,{h,r},Message[SetGroup::exceptionalGroupRank,h,"6,7,8"];Abort[]],"F",If[r==4,{h,r},Message[SetGroup::exceptionalGroupRank,h,"4"];Abort[]],"G",If[r==2,{h,r},Message[SetGroup::exceptionalGroupRank,h,"2"];Abort[]],"SU",If[r>=2,{"A",r-1},Message[SetGroup::u1Group,"SU"];Abort[]],"SP",If[r>=2,{"C",r},{"A",r}],"SO",Which[r==1,Message[SetGroup::u1Group,"SO"];Abort[],r==2||r==3,{"A",1},r==4,{{"A",1},{"A",1}},EvenQ[r],{"D",r/2},OddQ[r],{"B",(r-1)/2}],_,Message[SetGroup::unknownName,h];Abort[]]]&;
directProd:=Function[y,Flatten[Module[{i1=0,f},Module[{j=i1,l=Length[First[#]]},i1+=l;Function[x,ArrayPad[#,{j,tsDim-j-l}]&/@x][#]]&/@y],1]];
(*...................................Program...............................................................................*)
Unprotect[Group,Rank,SimpleRoots,PositiveRoots,FundamentalWeights,MaximalRoots,WeylVector,CoxeterNumber];
Clear[Group,Rank,SimpleRoots,PositiveRoots,FundamentalWeights,MaximalRoots,WeylVector,CoxeterNumber];
Group=Partition[Flatten[readG/@sGroup],2];
Rank=Total[Group[[All,2]]];
rootWeight=Switch[#[[1]],"A",RootsWeightsA[#[[2]]],"B",RootsWeightsB[#[[2]]],"C",RootsWeightsC[#[[2]]],"D",RootsWeightsD[#[[2]]],"E",RootsWeightsE[#[[2]]],"F",RootsWeightsF[#[[2]]],"G",RootsWeightsG[#[[2]]]]&/@Group;
tsDim=Total[Length/@rootWeight[[All,4]]];
SimpleRoots=directProd[rootWeight[[All,1]]];
PositiveRoots=directProd[rootWeight[[All,2]]];
FundamentalWeights=directProd[rootWeight[[All,3]]];
MaximalRoots=directProd[List/@rootWeight[[All,4]]];
WeylVector=1/2 Total[PositiveRoots];
CoxeterNumber=2 Length[PositiveRoots]/Rank;
Protect[Group,Rank,SimpleRoots,PositiveRoots,FundamentalWeights,MaximalRoots,WeylVector,CoxeterNumber];
Print["Set group is ",StringDrop[StringJoin@@(("\[Times] "<>ToString@#[[1]]<>"("<>ToString[#[[2]]]<>")")&/@Group),1]];
Print["Rank of the group = ",Rank];]


(* ::Section::Closed:: *)
(*`Private` commands*)


(* ::Subsection::Closed:: *)
(*PackageElemType*)


PackageElemType[elem_]:=Module[{l=Length@First[SimpleRoots]},Switch[elem,
{x__/;(And@@(NumericQ/@{x})&&Length[{x}]==l)},"weight",
{y__/;And@@(MatchQ[#,{x__/;(And@@(NumericQ/@{x})&&Length[{x}]==l)}]&/@{y})},"list of weights",
{y__/;And@@(MatchQ[#,{{x__/;(And@@(NumericQ/@{x})&&Length[{x}]==l)},_Integer}]&/@{y})},"multiplet",
_,"unknown"]]
PackageElemType[elem_,sRoots_]:=Module[{l=Length@First[sRoots]},Switch[elem,
{x__/;(And@@(NumericQ/@{x})&&Length[{x}]==l)},"weight",
{y__/;And@@(MatchQ[#,{x__/;(And@@(NumericQ/@{x})&&Length[{x}]==l)}]&/@{y})},"list of weights",
{y__/;And@@(MatchQ[#,{{x__/;(And@@(NumericQ/@{x})&&Length[{x}]==l)},_Integer}]&/@{y})},"multiplet",
_,"unknown"]]


(* ::Subsection::Closed:: *)
(*CheckGTPattern*)


CheckGTPattern::incorrectForm="Pattern has incorrect form";
CheckGTPattern::incorrectA="Elements of A-pattern have to be non-negative integers.";
CheckGTPattern::incorrectC="Elements of C-pattern have to be non-negative integers.";
CheckGTPattern::incorrectB1="The last element of B-pattern odd lines is 0 or 1";
CheckGTPattern::incorrectB2="Elements of B-pattern have to be non-negative integers or half-integers.";
CheckGTPattern::incorrectD="Elements of D-pattern have to be non-negative integers or half-integers (last one in can negative).";
CheckGTPattern::nonSnake="The pattern does not satisfy the snake-rule.";
CheckGTPattern[p_]:=Module[{n=Group[[1,2]],checkSnake,checkInt},
checkSnake=Function[x,Nand@@(GreaterEqual@@#&/@MapThread[Riffle[#1,#2]&,{Most[x],Rest[x]}])];
checkInt=Function[x,Nand@@((IntegerQ[#]&&NonNegative[#])&/@Flatten[x])];
Switch[Group[[1,1]],
"A",
If[Length/@p!=Range[n+1,1,-1],Message[CheckGTPattern::incorrectForm];Return[False]];
If[checkInt[p],Message[CheckGTPattern::incorrectA];Return[False]];
If[checkSnake[p],Message[CheckGTPattern::nonSnake];Return[False]];
Return[True];,
"B",
If[Length/@p!=Riffle[Range[n+1,2,-1],Range[n,1,-1]],Message[CheckGTPattern::incorrectForm];Return[False]];
If[Nand@@((#==1||#==0)&/@Flatten[Drop[p,{2,-1,2},{1,-2}]]),Message[CheckGTPattern::incorrectB1];Return[False]];
If[checkSnake[MapIndexed[If[OddQ@@#2,Most[#1],#1]&,p]],Message[CheckGTPattern::nonSnake];Return[False]];
If[And[checkInt[MapIndexed[If[OddQ@@#2,Most[#1],#1]&,p]],checkInt[MapIndexed[If[OddQ@@#2,Most[#1],#1]&,p]+1/2]],Message[CheckGTPattern::incorrectB2];Return[False]];
Return[True],
"C",
If[Length/@p!=Riffle[Range[n,1,-1],Range[n,1,-1]],Message[CheckGTPattern::incorrectForm];Return[False]];
If[checkInt[p],Message[CheckGTPattern::incorrectC];Return[False]];
If[checkSnake[p],Message[CheckGTPattern::nonSnake];Return[False]];
Return[True];,
"D",
If[Length/@p!=Riffle[Range[n,1,-1],Range[n-1,1,-1]],Message[CheckGTPattern::incorrectForm];Return[False]];
If[checkSnake[MapIndexed[If[OddQ@@#2,MapAt[Abs,#1,-1],#1]&,p]],Message[CheckGTPattern::nonSnake];Return[False]];
If[And[checkInt[MapIndexed[If[OddQ@@#2,MapAt[Abs,#1,-1],#1]&,p]],checkInt[MapIndexed[If[OddQ@@#2,MapAt[Abs,#1,-1],#1]&,p]+1/2]],Message[CheckGTPattern::incorrectD];Return[False]];
Return[True]]]


(* ::Subsection::Closed:: *)
(*REPRESENTATION*)


REPRESENTATION[g_,dynI_]:=Switch[g,
"A",(*..................Case of A algebra............................*)
Module[{upLine,calcLine,addLines,GTp,M,n=Length[dynI]+1},
(*................Addition functions......................*)
calcLine=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];
addLines=Function[x,Join[x,{#}]&/@calcLine[Last[x]]];
(*...............Program it-self..........................*)
upLine=Join[Reverse[Accumulate[Reverse[dynI]]],{0}];
M=Table[If[i==1,If[j==n-1,-1,0],If[i==n-j+1,2,If[i==n-j||i==n-j+2,-1,0]]],{j,n-1},{i,n}];
GTp=Nest[Flatten[addLines/@#,1]&,{{upLine}},Length[dynI]];
Tally[M.(Total/@#)&/@GTp]],
"B",(*..................Case of B algebra............................*)
Module[{calcLine1,calcLine2,addLines1,addLines2,GTp,M,n=Length[dynI],add\[Sigma],\[Sigma]},
(*................Addition functions......................*)
add\[Sigma]=Function[x,If[Last[x]!=0,{Join[x,{0}],Join[x,{1}]},{Join[x,{0}]}]];
If[EvenQ[Last[dynI]],
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];
calcLine2=Function[x,Flatten[add\[Sigma]/@calcLine1[x],1]];,
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x]~Join~{1/2},x}]]];
calcLine2=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];];
addLines1=Function[x,Join[x,{#}]&/@calcLine1[Last[x]]];
addLines2=Function[x,Join[x,{#}]&/@calcLine2[Last[x]]];
(*...............Program it-self..........................*)
M=Table[Which[j==n,Which[i==2n,4,i==2n-1,-2,True,0],i==2j-1,-1,i==2j,2,i==2j+2,-2,i==2j+3,1,True,0],{j,n},{i,2n}];
If[EvenQ[Last[dynI]],
GTp=Flatten[addLines1/@Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,Partition[add\[Sigma][Reverse[Accumulate[Reverse[MapAt[#/2&,dynI,-1]]]]],1],n-1],1];
\[Sigma]=Insert[(Most[#]-Rest[#]),2Last[#],-1]&[Flatten[Take[#,{1,-1,2},-1]]]&/@GTp;
Tally[MapThread[#1-#2&,{M.(Total/@(MapIndexed[If[OddQ@@#2,Most[#1],#1]&,#]))&/@GTp,\[Sigma]}]],
GTp=Flatten[addLines1/@Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,{{Reverse[Accumulate[Reverse[MapAt[#/2&,dynI,-1]]]]}},n-1],1];
\[Sigma]=Insert[(Most[#]-Rest[#]),2Last[#],-1]&/@Tuples[{0,1},n];
Tally[Flatten[Outer[#1-#2&,M.(Total/@#)&/@GTp,\[Sigma],1],1]]]],
"C",(*..................Case of C algebra............................*)
Module[{upLine,calcLine1,calcLine2,addLines1,addLines2,GTp,M,n=Length[dynI]},
(*................Addition functions......................*)
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x]~Join~{0},x}]]];
calcLine2=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];
addLines1=Function[x,Join[x,{#}]&/@calcLine1[Last[x]]];
addLines2=Function[x,Join[x,{#}]&/@calcLine2[Last[x]]];
(*...............Program it-self..........................*)
upLine=Reverse[Accumulate[Reverse[dynI]]];
M=Table[Which[i==2j-1,-1,i==2j,2,i==2j+2,-2,i==2j+3,1,True,0],{j,n},{i,2n}];
GTp=Flatten[addLines1/@Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,{{upLine}},Length[dynI]-1],1];
Tally[M.(Total/@#)&/@GTp]],
"D",(*..................Case of D algebra............................*)
Module[{upLine,calcLine1,calcLine2,addLines1,addLines2,GTp,M,tTotal,n=Length[dynI]},
(*................Addition functions......................*)
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{MapAt[Abs,Rest[x],-1],Most[x]}]]];
calcLine2=Function[x,Tuples[MapThread[Range[#1,#2]&,{Join[Rest[x],{-Last[x]}],x}]]];
addLines1=Function[x,Join[x,{#}]&/@calcLine1[Last[x]]];
addLines2=Function[x,Join[x,{#}]&/@calcLine2[Last[x]]];
tTotal=Function[y,Join[(Total/@y)+Function[x,MapIndexed[If[EvenQ@@#2,Min[x[[#2-1]],x[[#2+1]]],0]&,x]][y[[All,-1]]],{y[[-1,-1]]}]];
(*...............Program it-self..........................*)
upLine=MapAt[#-dynI[[-2]]/2&,Reverse[Accumulate[Reverse[MapAt[#/2&,dynI,{{-1},{-2}}]]]],-1];
M=Table[Which[j==n,Which[i==2n,2,i==2n-1,-2,i==2n-2,2,i==2n-3,-1,True,0],i==2j-1,-1,i==2j,2,i==2j+2,-2,i==2j+3,1,True,0],{j,n},{i,2n}];
GTp=Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,{{upLine}},Length[dynI]-1];
Tally[(M.tTotal[#])&/@GTp]],
_,(*..................Case of exceptional groups......................*)
Module[{\[Lambda],dim,f1,factor,getMultip,nextLayer,calcMultip,lvl,resultRep,curLayer,curDim,sRoot,fWeight,pRoot,wVector,toDynkI},
(*....................................Set root system......................................................................*)
Switch[g,"E",{sRoot,pRoot,fWeight}=Drop[RootsWeightsE[Length[dynI]],-1],
"F",{sRoot,pRoot,fWeight}=Drop[RootsWeightsF[Length[dynI]],-1],
"G",{sRoot,pRoot,fWeight}=Drop[RootsWeightsG[Length[dynI]],-1],
_,Print["Unknown group"];Abort[]];
wVector=1/2 Total[pRoot];
\[Lambda]=dynI.fWeight;
(*....................................Addition functions and constants.....................................................*)
dim=Times@@((\[Lambda].#/wVector.#+1)&/@pRoot);
f1=(\[Lambda]+wVector).(\[Lambda]+wVector);
factor:=(f1-(#+wVector).(#+wVector))&;
getMultip:=Module[{qq=Cases[resultRep,{#,_},1,1]},If[qq=={},0,qq[[1,2]]]]&;
nextLayer:=Module[{},lvl++;DeleteDuplicates[Flatten[Outer[(#1-#2)&,curLayer,sRoot,1],1]]];
calcMultip:=Function[\[Mu],Module[{f=factor[\[Mu]]},If[\[Mu].\[Mu]>\[Lambda].\[Lambda]||f==0,0,Total[2/f Flatten[Outer[(getMultip[\[Mu]+#2 #1](\[Mu]+#2 #1).#1)&,pRoot,Range[lvl],1],1]]]]];
(*....................................Program it-self.....................................................................*)
resultRep={{\[Lambda],1}};
curLayer={\[Lambda]};
lvl=0;
Monitor[NestWhile[Function[d,Module[{m=nextLayer,lCalc},lCalc=Cases[{#,calcMultip[#]}&/@m,{_,x_/;x!=0}];resultRep=Join[resultRep,lCalc];curLayer=lCalc[[All,1]];curDim=d+Total[lCalc[[All,2]]]]],1,#<dim&],ProgressIndicator[curDim,{1,dim}]];
(*.....................................Form output........................................................................*)
toDynkI=Function[x,2 sRoot[[#]].x/sRoot[[#]].sRoot[[#]]&/@Range[Length[dynI]]];
{toDynkI[#[[1]]],#[[2]]}&/@resultRep]]


(* ::Subsection::Closed:: *)
(*GTs*)


GTs[dynI_,g_]:=Switch[g,
"A",Module[{upLine,calcLine,addLines},
(*................Addition functions......................*)
calcLine=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];
addLines=Function[x,Join[x,{#}]&/@calcLine[Last[x]]];
(*...............Program it-self..........................*)
upLine=Join[Reverse[Accumulate[Reverse[dynI]]],{0}];
Nest[Flatten[addLines/@#,1]&,{{upLine}},Length[dynI]]],
"B",
Module[{upLine,calcLine1,calcLine2,addLines1,addLines2,add\[Sigma]},
(*................Addition functions......................*)
add\[Sigma]=Function[x,If[Last[x]!=0,{Join[x,{0}],Join[x,{1}]},{Join[x,{0}]}]];
If[EvenQ[Last[dynI]],
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];
calcLine2=Function[x,Flatten[add\[Sigma]/@calcLine1[x],1]];,
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x]~Join~{1/2},x}]]];
calcLine2=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];];
addLines1=Function[x,Join[x,{#}]&/@calcLine1[Last[x]]];
addLines2=Function[x,Join[x,{#}]&/@calcLine2[Last[x]]];
(*...............Program it-self..........................*)
If[EvenQ[Last[dynI]],upLine=Partition[add\[Sigma][Reverse[Accumulate[Reverse[MapAt[#/2&,dynI,-1]]]]],1];
Flatten[addLines1/@Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,upLine,Length[dynI]-1],1],upLine={{Reverse[Accumulate[Reverse[MapAt[#/2&,dynI,-1]]]]}};
Module[{gt=Flatten[addLines1/@Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,upLine,Length[dynI]-1],1],
t=Insert[Riffle[#,0],0,-1]&/@Tuples[{0,1},Length[dynI]],n=Range[2Length[dynI]]},
Flatten[Outer[Function[{x,y},MapThread[If[OddQ[#3],Insert[#1,#2,-1],#1]&,{x,y,n}]],gt,t,1],1]]]],
"C",Module[{upLine,calcLine1,calcLine2,addLines1,addLines2},
(*................Addition functions......................*)
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x]~Join~{0},x}]]];
calcLine2=Function[x,Tuples[MapThread[Range[#1,#2]&,{Rest[x],Most[x]}]]];
addLines1=Function[x,Join[x,{#}]&/@calcLine1[Last[x]]];
addLines2=Function[x,Join[x,{#}]&/@calcLine2[Last[x]]];
(*...............Program it-self..........................*)
upLine=Reverse[Accumulate[Reverse[dynI]]];
Flatten[addLines1/@Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,{{upLine}},Length[dynI]-1],1]],
"D",
Module[{upLine,calcLine1,calcLine2,addLines1,addLines2},
(*................Addition functions......................*)
calcLine1=Function[x,Tuples[MapThread[Range[#1,#2]&,{MapAt[Abs,Rest[x],-1],Most[x]}]]];
calcLine2=Function[x,Tuples[MapThread[Range[#1,#2]&,{Join[Rest[x],{-Last[x]}],x}]]];
addLines1=Function[x,Join[x,{#}]&/@calcLine1[Last[x]]];
addLines2=Function[x,Join[x,{#}]&/@calcLine2[Last[x]]];
(*...............Program it-self..........................*)
upLine=MapAt[#-dynI[[-2]]/2&,Reverse[Accumulate[Reverse[MapAt[#/2&,dynI,{{-1},{-2}}]]]],-1];
Nest[Flatten[addLines2/@(Flatten[addLines1/@#,1]),1]&,{{upLine}},Length[dynI]-1]],
_,Message[GTs::wrongArg];Abort[]];


(* ::Section:: *)
(*General commands*)


(* ::Subsection::Closed:: *)
(*CartanMatrix*)


CartanMatrix::badRoot="The roots are of the diferent size.";
CartanMatrix[SimpleRoots_]:=Check[Outer[2 #1.#2/ #2.#2&,SimpleRoots,SimpleRoots,1],Message[CartanMatrix::badRoot];Abort[],{Dot::dotsh}]


(* ::Input:: *)
(**)


(* ::Subsection::Closed:: *)
(*FundamentalWeight*)


FundamentalWeight::badRoot="The root system is degenerate.";
FundamentalWeight[SysRoots_]:=Check[Inverse[CartanMatrix[SysRoots]].SysRoots,Message[FundamentalWeight::badRoot];Abort[],{Inverse::sing}]


(* ::Subsection::Closed:: *)
(*RootSystem*)


RootSystem[sRoots_]:=Module[{L=FundamentalWeight[sRoots],iteration,curL=0,res,IsPositive},
iteration:=Function[x,Module[{i=DeleteDuplicates@Simplify[Flatten[Outer[#1-2 #1.#2/#2.#2 #2&,x,x,1],1]]},curL=Length[x];DeleteDuplicates[Join[i,-i]]]];
IsPositive:=Function[x,Catch[Module[{y=x.#},If[y>0,Throw[True],If[y<0,Throw[False]]]]&/@L;True]];
res=GatherBy[NestWhile[iteration,sRoots,Length[#]>curL&],IsPositive[#]&];
If[IsPositive[res[[1,1]]],res,{res[[2]],res[[1]]}]]


(* ::Subsection::Closed:: *)
(*Representation*)


Representation::nonIntDynkin="Dynkin labels `1` are not non-negative integers.";
Representation::notEldest="`1` is not an eldest weight for any representation of the current group.";
Representation::incorrectLength="`1` should have `2` number of components.";
Representation::unknownOpts="`1` is unknonw form of notation. Switch to DynkinLabels";
Representation[Rep_,sRoots_,opts___]:=Module[{rNot,\[Lambda],dim,f1,factor,getMultip,nextLayer,calcMultip,lvl,resultRep,curLayer,curDim,wVector,pRoots,rank,fWeight},
rank=Length[sRoots];
(*....................................Options desdecoding and input-check..........................................................*)
{rNot}=(({RepNotation}/.{opts})/.{RepNotation->"DynkinLabels"});
Switch[rNot,
"DynkinLabels",If[Length[Rep]==rank,If[And@@((NonNegative[#]&&IntegerQ[#])&/@Rep),Null,Message[Representation::nonIntDynkin,Rep];Abort[]],Message[Representation::incorrectLength,Rep,rank];Abort[]],
"EldestWeight",If[(Length[Rep]==Length[First[sRoots]])&&(And@@(NumericQ/@Rep)),If[And@@((NonNegative[#]&&IntegerQ[#])&/@(2 Rep.#/#.#&/@sRoots)),Null,Message[Representation::notEldest,Rep];Abort[]],Message[Representation::incorrectLength,Rep,Length[First[sRoots]]];Abort[]],
_,Message[Representation::unknownOpts,rNot];If[Length[Rep]==rank,If[And@@((NonNegative[#]&&IntegerQ[#])&/@Rep),rNot="DynkinLabels",Message[Representation::nonIntDynkin,Rep];Abort[]],Message[Representation::incorrectLength,Rep,rank];Abort[]]];
(*....................................Seting the root system..............................................................*)
PrintTemporary["Building root system."]; 
fWeight=Simplify[FundamentalWeight[sRoots]];
pRoots=Simplify[RootSystem[sRoots][[1]]];
wVector=Simplify[Total[pRoots]/2];
PrintTemporary["Done. Calculating weight diagram:"];
(*....................................Addition functions and constants.....................................................*)
\[Lambda]=Simplify[If[rNot=="DynkinLabels",Rep.fWeight,Rep]];
dim=Simplify[Times@@((\[Lambda].#/wVector.#+1)&/@pRoots)];
f1=Simplify[(\[Lambda]+wVector).(\[Lambda]+wVector)];
factor:=Simplify[(f1-(#+wVector).(#+wVector))]&;
getMultip:=Module[{qq=Cases[resultRep,{#,_},1,1]},If[qq=={},0,qq[[1,2]]]]&;
nextLayer:=Module[{},lvl++;DeleteDuplicates[Flatten[Outer[Simplify[(#1-#2)]&,curLayer,sRoots,1],1]]];
calcMultip:=Function[\[Mu],Module[{f=factor[\[Mu]]},If[\[Mu].\[Mu]>\[Lambda].\[Lambda]||f==0,0,Total[2/f Flatten[Outer[(getMultip[Simplify[\[Mu]+#2 #1]]Simplify[(\[Mu]+#2 #1)].#1)&,pRoots,Range[lvl],1],1]]]]];
(*....................................Program it-self.....................................................................*)
resultRep={{\[Lambda],1}};
curLayer={\[Lambda]};
lvl=0;
Monitor[NestWhile[Function[d,Module[{m=nextLayer,lCalc},lCalc=Cases[{#,calcMultip[#]}&/@m,{_,x_/;x!=0}];resultRep=Join[resultRep,lCalc];curLayer=lCalc[[All,1]];curDim=d+Total[lCalc[[All,2]]]]],1,#<dim&],ProgressIndicator[curDim,{1,dim}]];
resultRep]


(* ::Subsection::Closed:: *)
(*RepProduct*)


RepProduct::wrongArg="The arguments are not multiplets.";
RepProduct[R1_,R2_]:=Module[{},
If[Not[MatchQ[R1,{y__/;And@@(MatchQ[#,{{x__/;(And@@(NumericQ/@{x}))},_Integer}]&/@{y})}]&&MatchQ[R2,{y__/;And@@(MatchQ[#,{{x__/;(And@@(NumericQ/@{x}))},_Integer}]&/@{y})}]&&Equal@@((Length[#[[1]]]&/@R1)~Join~(Length[#[[1]]]&/@R2))],Message[RepProduct::wrongArg];Abort[];];
MapAt[Total,MapAt[First,Transpose[#],1],2]&/@Gather[Flatten[Outer[{Simplify[#1[[1]]+#2[[1]]],#1[[2]]#2[[2]]}&,R1,R2,1],1],#1[[1]]==#2[[1]]&]];


(* ::Input:: *)
(**)


(* ::Subsection::Closed:: *)
(*DecomposeRep*)


DecomposeRep::wrongArgRep="Wrong format of the weight diagram.";
DecomposeRep::wrongArgRoot="Wrong format of roots.";
DecomposeRep[Rep_,sysRoots_]:=Module[{fWeight,pRoots,wVector,dimension,elder,cutDomSector,getSubMultip,subtract,domwSub,i,tDim,cDim,startM},
(*....................................Options desdecoding and input-check..........................................................*)
If[PackageElemType[Rep,sysRoots]!="multiplet",Message[DecomposeRep::wrongArgRep];Abort[]];
If[PackageElemType[sysRoots,sysRoots]!="list of weights",Message[DecomposeRep::wrongArgRoot];Abort[]];
(*....................................Constructing the root system.................................................................*)
PrintTemporary["Initializing."]; 
fWeight=Simplify[FundamentalWeight[sysRoots]];
pRoots=Simplify[RootSystem[sysRoots][[1]]];
wVector=Simplify[Total[pRoots]/2];
PrintTemporary["Done. Calculating decomposition:"];
(*....................................Addition functions and constants.....................................................*)
elder:=Function[{w1,w2},Module[{W1=w1.#/#.#&/@sysRoots,W2=w2.#/#.#&/@sysRoots,c1,c2},c1=And@@(NonNegative/@W1);c2=And@@(NonNegative/@W2);If[c1==c2,w1.w1>w2.w2,c1]]];
cutDomSector=Select[#,Function[w,And@@(NonNegative/@(2 w[[1]].#/#.#&/@sysRoots))]]&;
subtract:=Function[{r1,r2},(Select[Replace[r1,({#[[1]],x_}:>Evaluate[{#[[1]],x-#[[2]]}])&/@r2,1],#[[2]]!=0&])];
dimension:=Function[x,Simplify[Times@@((x.#/wVector.#+1)&/@pRoots)]];
(*..............This long sub-program is actually cuted Representation, i.e. it calculates weights only in positive sector.........*)
getSubMultip=Function[\[Lambda],Module[{f1,factor,resultRep,getMultip,nextLayer,calcMultip,lvl},
f1=Simplify[(\[Lambda]+wVector).(\[Lambda]+wVector)];
factor:=Simplify[(f1-(#+wVector).(#+wVector))]&;
getMultip:=Module[{qq=Cases[resultRep,{#,_},1,1]},If[qq=={},0,qq[[1,2]]]]&;
nextLayer:=Function[layer,Module[{},lvl++;Select[DeleteDuplicates[Flatten[Outer[Simplify[(#1-#2)]&,layer,sysRoots,1],1]],Function[x,And@@(NonNegative/@(x.#&/@fWeight))]]]];
calcMultip:=Function[\[Mu],Module[{f=factor[\[Mu]]},If[\[Mu].\[Mu]>\[Lambda].\[Lambda]||f==0,0,Total[2/f Flatten[Outer[(getMultip[Simplify[\[Mu]+#2 #1]]Simplify[(\[Mu]+#2 #1).#1])&,pRoots,Range[lvl],1],1]]]]];
resultRep={{\[Lambda],1}};lvl=0;
NestWhile[Function[curL,Module[{m=nextLayer[curL],lCalc},lCalc=Cases[{#,calcMultip[#]}&/@m,{_,x_/;x!=0}];resultRep=Join[resultRep,lCalc];cDim+=Total[lCalc[[All,2]]];lCalc[[All,1]]]],{\[Lambda]},#!={}&];
resultRep]];
(*..............Program it-self........................................................................................*)
domwSub={};i=1;tDim=Simplify[Total[Select[Rep,Function[w,And@@(NonNegative/@(w[[1]].#&/@fWeight))]][[All,2]]]];cDim=0;
startM=Simplify[Sort[cutDomSector[Rep],elder[#1[[1]],#2[[1]]]&]];
Monitor[NestWhile[Function[restM,Module[{lll=restM[[1,1]]},AppendTo[domwSub,lll];Print["Sub-multiplet ",i++,":  ",Simplify[2 lll.#/#.#&/@sysRoots],"  dim=",dimension[lll]];subtract[restM,cutDomSector[getSubMultip[lll]]]]],startM,#!={}&],ProgressIndicator[cDim,{0,tDim}]];
Print["Total number of sub-mupltiplets: ",i-1];
Simplify[Outer[2 #1.#2/#2.#2&,domwSub,sysRoots,1]]]


(* ::Subsection::Closed:: *)
(*PhysicalNormalizationMatrix*)


PhysicalNormalizationMatrix::wrongMethod="The NormalizationMethod should be Coxeter, Unity or Half";
PhysicalNormalizationMatrix::wrongGroup="The label `1` is unknown label for the group. Please, use A,B,C,D,E,F,G -notation.";
PhysicalNormalizationMatrix::wrongRank="The rank `1` is not a positive integer";
PhysicalNormalizationMatrix[g_,n_,opts___?OptionQ]:=Module[{method,posInt},
posInt:=(IntegerQ[#]&&Positive[#])&;
method=(NormalizationMethod/.{opts})/.{NormalizationMethod->"Coxeter"};
If[method=="Canonical",Return[IdentityMatrix[n]]];
If[method!="Half"&&method!="Coxeter"&&method!="Unity",Message[PhysicalNormalizationMatrix::wrongMethod];Abort[]];
Switch[{g,n},
{"A",_?posInt},Switch[method,"Half",Sqrt[1/(n+1)],"Coxeter",1,"Unity",Sqrt[2/(n+1)]]Table[Which[i==j-1,-Sqrt[i/(2(i+1))],i==(n+1),Sqrt[1/(2(n+1))],i>=j,Sqrt[1/(2i (i+1))],True,0],{i,n+1},{j,n+1}],
{"B",_?posInt},Switch[method,"Half",Sqrt[1/(2n-1)],"Coxeter",Sqrt[(2n)/(2n-1)],"Unity",Sqrt[2/(2n-1)]]Table[Which[i==j-1,-Sqrt[i/(2(i+1))],i==n,Sqrt[1/(2n)],i>=j,Sqrt[1/(2i (i+1))],True,0],{i,n},{j,n}],
{"C",_?posInt},Switch[method,"Half",Sqrt[1/(2(n+1))],"Coxeter",Sqrt[n/(n+1)],"Unity",Sqrt[1/(n+1)]]Table[Which[i==j-1,-Sqrt[i/(2(i+1))],i==n,Sqrt[1/(2n)],i>=j,Sqrt[1/(2i (i+1))],True,0],{i,n},{j,n}],
{"D",_?posInt},Switch[method,"Half",Sqrt[1/(2(n-1))],"Coxeter",1,"Unity",Sqrt[1/(n-1)]]Table[Which[i==j-1,-Sqrt[i/(2(i+1))],i==n,Sqrt[1/(2n)],i>=j,Sqrt[1/(2i (i+1))],True,0],{i,n},{j,n}],
{"E",6},Switch[method,"Half",1/(2Sqrt[3]),"Coxeter",1,"Unity",Sqrt[1/6]]{{1/4,-(1/4),-(1/4),-(1/4),-(1/4),-(1/4),-(1/4),1/4},{-(Sqrt[3]/4),Sqrt[3]/4,-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),1/(4 Sqrt[3])},{-(Sqrt[(3/2)]/4),-(Sqrt[(3/2)]/4),5/(4 Sqrt[6]),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),1/(4 Sqrt[6])},{-(3/(4 Sqrt[10])),-(3/(4 Sqrt[10])),-(3/(4 Sqrt[10])),7/(4 Sqrt[10]),-(1/(4 Sqrt[10])),-(1/(4 Sqrt[10])),-(1/(4 Sqrt[10])),1/(4 Sqrt[10])},{-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),(3 Sqrt[3/5])/4,-(1/(4 Sqrt[15])),-(1/(4 Sqrt[15])),1/(4 Sqrt[15])},{1/4,1/4,1/4,1/4,1/4,-(1/4),-(1/4),1/4},{0,0,0,0,0,1/(2 Sqrt[3]),-(1/Sqrt[3]),-(1/(2 Sqrt[3]))},{0,0,0,0,0,-(1/2),0,-(1/2)}},
{"E",7},Switch[method,"Half",1/(3Sqrt[2]),"Coxeter",1,"Unity",1/3]{{1/4,-(1/4),-(1/4),-(1/4),-(1/4),-(1/4),-(1/4),1/4},{-(Sqrt[3]/4),Sqrt[3]/4,-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),1/(4 Sqrt[3])},{-(Sqrt[(3/2)]/4),-(Sqrt[(3/2)]/4),5/(4 Sqrt[6]),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),1/(4 Sqrt[6])},{-(3/(4 Sqrt[10])),-(3/(4 Sqrt[10])),-(3/(4 Sqrt[10])),7/(4 Sqrt[10]),-(1/(4 Sqrt[10])),-(1/(4 Sqrt[10])),-(1/(4 Sqrt[10])),1/(4 Sqrt[10])},{-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),(3 Sqrt[3/5])/4,-(1/(4 Sqrt[15])),-(1/(4 Sqrt[15])),1/(4 Sqrt[15])},{-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),11/(4 Sqrt[21]),-(1/(4 Sqrt[21])),1/(4 Sqrt[21])},{1/(2 Sqrt[7]),1/(2 Sqrt[7]),1/(2 Sqrt[7]),1/(2 Sqrt[7]),1/(2 Sqrt[7]),1/(2 Sqrt[7]),-(1/Sqrt[7]),1/Sqrt[7]},{0,0,0,0,0,0,1/2,1/2}},
{"E",8},Switch[method,"Half",1/Sqrt[30],"Coxeter",1,"Unity",Sqrt[1/15]]{{1/4,-(1/4),-(1/4),-(1/4),-(1/4),-(1/4),-(1/4),1/4},{-(Sqrt[3]/4),Sqrt[3]/4,-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),-(1/(4 Sqrt[3])),1/(4 Sqrt[3])},{-(Sqrt[(3/2)]/4),-(Sqrt[(3/2)]/4),5/(4 Sqrt[6]),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),-(1/(4 Sqrt[6])),1/(4 Sqrt[6])},{-(3/(4 Sqrt[10])),-(3/(4 Sqrt[10])),-(3/(4 Sqrt[10])),7/(4 Sqrt[10]),-(1/(4 Sqrt[10])),-(1/(4 Sqrt[10])),-(1/(4 Sqrt[10])),1/(4 Sqrt[10])},{-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),-(Sqrt[(3/5)]/4),(3 Sqrt[3/5])/4,-(1/(4 Sqrt[15])),-(1/(4 Sqrt[15])),1/(4 Sqrt[15])},{-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),-(Sqrt[(3/7)]/4),11/(4 Sqrt[21]),-(1/(4 Sqrt[21])),1/(4 Sqrt[21])},{-(3/(8 Sqrt[7])),-(3/(8 Sqrt[7])),-(3/(8 Sqrt[7])),-(3/(8 Sqrt[7])),-(3/(8 Sqrt[7])),-(3/(8 Sqrt[7])),13/(8 Sqrt[7]),1/(8 Sqrt[7])},{-(1/8),-(1/8),-(1/8),-(1/8),-(1/8),-(1/8),-(1/8),-(5/8)}},
{"F",4},Switch[method,"Half",1/3,"Coxeter",2/Sqrt[3],"Unity",Sqrt[2]/3]{{0,1/2,-(1/2),0},{0,1/(2 Sqrt[3]),1/(2 Sqrt[3]),-(1/Sqrt[3])},{0,1/Sqrt[6],1/Sqrt[6],1/Sqrt[6]},{-(1/Sqrt[2]),0,0,0}},
{"G",2},Switch[method,"Half",1/(2Sqrt[3]),"Coxeter",1/Sqrt[2],"Unity",1/Sqrt[6]]{{1/2,-(1/2),0},{-(1/(2 Sqrt[3])),-(1/(2 Sqrt[3])),1/Sqrt[3]},{-(1/Sqrt[6]),-(1/Sqrt[6]),-(1/Sqrt[6])}},
{_,_?IntegerQ},Message[PhysicalNormalizationMatrix::wrongGroup,g[n]],
{_,_},Message[PhysicalNormalizationMatrix::wrongRank,n]]]


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Commands for "working" group only*)


(* ::Subsection::Closed:: *)
(*RepDimension*)


RepDimension::notEldest="`1` is not an eldest weight for any representation of the current group.";
RepDimension::incorrectLength="`1` should have `2` number of components.";
RepDimension::unknownOpts="`1` is unknonw form of notation. Switch to DynkinLabels";
RepDimension[Rep_List,opts___]:=Module[{rNot,\[Lambda]},
{rNot}=(({RepNotation}/.{opts})/.{RepNotation->"DynkinLabels"});
Switch[rNot,
"DynkinLabels",If[Length[Rep]==Rank,\[Lambda]=Rep.FundamentalWeights,Message[RepDimension::incorrectLength,Rep,Rank];Abort[]],
"EldestWeight",If[Length[Rep]==Length[First[SimpleRoots]],If[And@@((NonNegative[#]&&IntegerQ[#])&/@(2 Rep.#/#.#&/@SimpleRoots)),\[Lambda]=Rep,Message[RepDimension::notEldest,Rep];Abort[]],Message[RepDimension::incorrectLength,Rep,Length[First[SimpleRoots]]];Abort[]],
_,Message[RepDimension::unknownOpts,rNot];If[Length[Rep]==Rank,\[Lambda]=Rep.FundamentalWeights,Message[RepDimension::notDynkin,Rep];Abort[]]];
Times@@((\[Lambda].#/WeylVector.#+1)&/@PositiveRoots)]


(* ::Subsection::Closed:: *)
(*DynkinDiagram*)


DynkinDiagram:=Module[{k=0,origin={0,0},g={},DyD,j=1},
DyD[{t_,n_,k_,j_}]:=
Module[{},Switch[{t,n},{"A",1},origin+={0,-0.85};Graphics[{Circle[{1,0}+origin,0.15],Text[1+k,{1,0}+origin]}],
{"A",_},origin+={0,-1.85};Graphics[{Table[{Circle[{i,0}+origin,0.15],Text[i+k,{i,0}+origin]},{i,n}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,n-1}],{Circle[{(n+1)/2,1}+origin,0.15],Text["max",{(n+1)/2,1.04}+origin],Text[j,{(n+1)/2,0.96}+origin],{Dashed,Line[{{(n+1)/2-0.15,1}+origin,{1,0.15}+origin}],Line[{{(n+1)/2+0.15,1}+origin,{n,0.15}+origin}]}}}],
{"B",2},origin+={0,-1.85};Graphics[{Table[{Circle[{i,0}+origin,0.15],Text[i+k,{i,0}+origin]},{i,2}],{Line[{{n-1+0.144,0.03}+origin,{n-0.144,0.03}+origin}],Line[{{n-1+0.144,-0.03}+origin,{n-0.144,-0.03}+origin}],Line[{{n-1+0.43,0.07}+origin,{n-1+0.57,0}+origin,{n-1+0.43,-0.07}+origin}]},{Circle[{3,0}+origin,0.15],Text["max",{3,0.04}+origin],Text[j,{3,-0.04}+origin],{Line[{{2+0.57,0.07}+origin,{2+0.43,0}+origin,{2+0.57,-0.07}+origin}]},Dashed,{Line[{{2+0.144,0.03}+origin,{3-0.144,0.03}+origin}],Line[{{2+0.144,-0.03}+origin,{3-0.144,-0.03}+origin}]}}}],
{"B",_},origin+={0,-1.85};Graphics[{Table[{Circle[{i,0}+origin,0.15],Text[i+k,{i,0}+origin]},{i,n}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,n-2}],{Line[{{n-1+0.144,0.03}+origin,{n-0.144,0.03}+origin}],Line[{{n-1+0.144,-0.03}+origin,{n-0.144,-0.03}+origin}],Line[{{n-1+0.43,0.07}+origin,{n-1+0.57,0}+origin,{n-1+0.43,-0.07}+origin}]},{Circle[{2,1}+origin,0.15],Text["max",{2,1.04}+origin],Text[j,{2,0.96}+origin],Dashed,Line[{{2,1-0.15}+origin,{2,0.15}+origin}]}}],
{"C",_},origin+={0,-0.85};Graphics[{Table[{Circle[{i,0}+origin,0.15],Text[i+k,{i,0}+origin]},{i,n}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,n-2}],{Line[{{n-1+0.144,0.03}+origin,{n-0.144,0.03}+origin}],Line[{{n-1+0.144,-0.03}+origin,{n-0.144,-0.03}+origin}],Line[{{n-1+0.57,0.07}+origin,{n-1+0.43,0}+origin,{n-1+0.57,-0.07}+origin}]},{Circle[{0,0}+origin,0.15],Text["max",{0,0.04}+origin],Text[j,{0,-.04}+origin],{Line[{{0.43,0.07}+origin,{0.57,0}+origin,{0.43,-0.07}+origin}]},Dashed,{Line[{{0.144,0.03}+origin,{1-0.144,0.03}+origin}],Line[{{0.144,-0.03}+origin,{1-0.144,-0.03}+origin}]}}}],
{"D",_},origin+={0,-1.35};Graphics[{Table[{Circle[{i,0}+origin,0.15],Text[i+k,{i,0}+origin]},{i,2,n-2}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,2,n-3}],{Circle[{1,-0.5}+origin,0.15],Circle[{n-1,-0.5}+origin,0.15],Circle[{n-1,0.5}+origin,0.15],Circle[{1,0.5}+origin,0.15],Text[1+k,{1,-0.5}+origin],Text[n+k,{n-1,-0.5}+origin],Text[n-1+k,{n-1,0.5}+origin],Text["max",{1,0.5+0.04}+origin],Text[j,{1,0.5-0.04}+origin],Line[{{1+0.15,-0.5}+origin,{2,-0.15}+origin}],{Dashed,Line[{{1+0.15,+0.5}+origin,{2,+0.15}+origin}]},Line[{{n-1-0.15,-0.5}+origin,{n-2,-0.15}+origin}],Line[{{n-1-0.15,+0.5}+origin,{n-2,+0.15}+origin}]}}],
{"E",6},origin+={0,-2.85};Graphics[{Circle[{1,0}+origin,0.15],Text[1+k,{1,0}+origin],Circle[{3,1}+origin,0.15],Text[2+k,{3,1}+origin],Table[{Circle[{i-1,0}+origin,0.15],Text[i+k,{i-1,0}+origin]},{i,3,6}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,4}],Line[{{3,1-0.15}+origin,{3,0.15}+origin}],{Circle[{3,2}+origin,0.15],Text["max",{3,2+0.04}+origin],Text[j,{3,2-0.04}+origin],{Dashed,Line[{{3,2-0.15}+origin,{3,1.15}+origin}]}}}],
{"E",_},origin+={0,-1.85};Graphics[{Circle[{1,0}+origin,0.15],Text[1+k,{1,0}+origin],Circle[{3,1}+origin,0.15],Text[2+k,{3,1}+origin],Table[{Circle[{i-1,0}+origin,0.15],Text[i+k,{i-1,0}+origin]},{i,3,n}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,n-2}],Line[{{3,1-0.15}+origin,{3,0.15}+origin}],If[n==7,{Circle[{0,0}+origin,0.15],Text["max",{0,0.04}+origin],Text[j,{0,-0.04}+origin],{Dashed,Line[{{0.15,0}+origin,{1-0.15,0}+origin}]}},{Circle[{8,0}+origin,0.15],Text["max",{8,+0.04}+origin],Text[j,{8,-0.04}+origin],{Dashed,Line[{{8-0.15,0}+origin,{7+0.15,0}+origin}]}}]}],
{"F",4},origin+={0,-0.85};Graphics[{Table[{Circle[{i,0}+origin,0.15],Text[i+k,{i,0}+origin]},{i,n}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,1}],{Line[{{2+0.144,0.03}+origin,{3-0.144,0.03}+origin}],Line[{{2+0.144,-0.03}+origin,{3-0.144,-0.03}+origin}],Line[{{3+0.15,0}+origin,{4-0.15,0}+origin}],Line[{{2+0.43,0.07}+origin,{2+0.57,0}+origin,{2+0.43,-0.07}+origin}]},{Circle[{0,0}+origin,0.15],Text["max",{0,0.04}+origin],Text[j,{0,-0.04}+origin],Dashed,Line[{{0.15,0}+origin,{1-0.15,0}+origin}]}}],
{"G",2},origin+={0,-0.85};Graphics[{Table[{Circle[{i,0}+origin,0.15],Text[i+k,{i,0}+origin]},{i,n}],Table[Line[{{i+0.15,0}+origin,{i+1-0.15,0}+origin}],{i,n-1}],{Line[{{n-1+0.144,0.03}+origin,{n-0.144,0.03}+origin}],Line[{{n-1+0.144,-0.03}+origin,{n-0.144,-0.03}+origin}],Line[{{n-1+0.57,0.07}+origin,{n-1+0.43,0}+origin,{n-1+0.57,-0.07}+origin}]},{Circle[{3,0}+origin,0.15],Text["max",{3,0.04}+origin],Text[j,{3,-0.04}+origin],Dashed,Line[{{2.15,0}+origin,{3-0.15,0}+origin}]}}]]];
(k+=#[[2]];g=Append[g,DyD[{#[[1]],#[[2]],k-#[[2]],j++}]];)&/@Group;
Show[g]];


(* ::Subsection::Closed:: *)
(*CoxeterPlane*)


CoxeterPlane::direcProduct="The Coxeter plane can be built only for simple groups (no direct products)";
CoxeterPlane::lowRank="The Coxeter plane has no sense for groups with rank one or two";
CoxeterPlane:=If[Length[Group]>1,Message[CoxeterPlane::direcProduct],
If[Group[[1,2]]<3,Message[CoxeterPlane::lowRank],
Module[{c1,JK},
JK:=Switch[Group[[1,1]],
"A",{Range[1,Group[[1,2]],2],Range[2,Group[[1,2]],2]},
"B",{Range[1,Group[[1,2]],2],Range[2,Group[[1,2]],2]},
"C",{Range[1,Group[[1,2]],2],Range[2,Group[[1,2]],2]},
"D",{Join[Range[1,Group[[1,2]]-1,2],{Group[[1,2]]}],Range[2,Group[[1,2]]-1,2]},
"E",{Join[{2},Range[3,Group[[1,2]],2]],Join[{1},Range[4,Group[[1,2]],2]]},
"F",{{1,3},{2,4}},
"G",{{1},{2}}];(*Subprogram of separation of roots onto ortogonal sub-spaces*)
c1=First@Select[Eigenvectors[Outer[N[#1.#2]&,SimpleRoots,SimpleRoots,1]],And@@(Positive/@#) || And@@(Negative/@#)&]; (*Fully positive (or negative) eigen vector of Coxeter matrix*)
{Total[c1[[#]]SimpleRoots[[#]]&/@JK[[1]]],Total[c1[[#]]SimpleRoots[[#]]&/@JK[[2]]]}]]]


(* ::Input:: *)
(**)


(* ::Subsection::Closed:: *)
(*BuildRep*)


BuildRep::nonIntDynkin="Dynkin labels `1` are not non-negative integers.";
BuildRep::notEldest="`1` is not an eldest weight for any representation of the current group.";
BuildRep::incorrectLength="`1` should have `2` number of components.";
BuildRep::unknownOpts="`1` is unknonw form of option `2`. `3` is used instead.";
BuildRep::notSimple="The GT-patterns have sence only for simple-classical groups.";
BuildRep[Rep_,opts___]:=Module[{rNot,outM,normM,\[Lambda],inForm,dProd,result,resList},
(*....................................Options desdecoding and input-check..........................................................*)
{rNot,outM,normM}=(({RepNotation,OutputMethod,NormalizationMethod}/.Flatten[{opts}])/.{RepNotation->"DynkinLabels",OutputMethod->"DynkinLabels",NormalizationMethod->"Canonical"});
Switch[rNot,
"DynkinLabels",If[Length[Rep]==Rank,If[And@@((NonNegative[#]&&IntegerQ[#])&/@Rep),\[Lambda]=Rep,Message[BuildRep::nonIntDynkin,Rep];Abort[]],Message[BuildRep::incorrectLength,Rep,Rank];Abort[]],
"EldestWeight",If[(Length[Rep]==Length[First[SimpleRoots]])&&(And@@(NumericQ/@Rep)),If[And@@((NonNegative[#]&&IntegerQ[#])&/@(2 Rep.#/#.#&/@SimpleRoots)),\[Lambda]=(2Rep.SimpleRoots[[#]])/SimpleRoots[[#]].SimpleRoots[[#]]&/@Range[Rank],Message[BuildRep::notEldest,Rep];Abort[]],Message[BuildRep::incorrectLength,Rep,Length[First[SimpleRoots]]];Abort[]],
_,Message[BuildRep::unknownOpts,rNot,"RepNotation","DynkinLabels"];If[Length[Rep]==Rank,If[And@@((NonNegative[#]&&IntegerQ[#])&/@Rep),\[Lambda]=Rep,Message[BuildRep::nonIntDynkin,Rep];Abort[]],Message[BuildRep::incorrectLength,Rep,Rank];Abort[]]];
(*....................................Additional functions.................................................................*)
dProd=Function[{x,y},Flatten[Outer[{Join[#1[[1]],#2[[1]]],#1[[2]]#2[[2]]}&,x,y,1],1]];
(*....................................Programm it-self.....................................................................*)
If[outM=="GTPatterns",If[Length[Group]==1&&MemberQ[{"A","B","C","D"},Group[[1,1]]],Return[GTs[\[Lambda],Group[[1,1]]]],Message[BuildRep::notSimple];Abort[]]];
inForm=Transpose[{Group[[All,1]],Function[x,Module[{t=Accumulate[Group[[All,2]]]},MapThread[Take[x,{#1+1,#2}]&,{Insert[Most[t],0,1],t}]]][\[Lambda]]}];
resList=REPRESENTATION@@#&/@inForm;
result=Fold[dProd,resList[[1]],Drop[resList,1]];(*....................................Formating output......................................................................*)
Switch[outM,"DynkinLabels",Return[result],
"Weights",
Which[normM=="Canonical",Return[{#[[1]].FundamentalWeights,#[[2]]}&/@result],
MemberQ[{"Unity","Coxeter","Half"},normM],
Module[{normMatrix,k={}},normMatrix=Module[{ph=(PhysicalNormalizationMatrix[#[[1]],#[[2]],NormalizationMethod->normM]&/@Group),m},ArrayFlatten[DiagonalMatrix[Array[m,Length[ph]]]/.{m[i_]:>ph[[i]]}]];
Module[{i=0},
Function[Switch[#1,
"A",AppendTo[k,{i+#2+1}];i+=(#2+1);,
"E",If[#2==6,AppendTo[k,i+{7}];AppendTo[k,i+{8}];i+=8;,If[#2==7,AppendTo[k,{i+8}];i+=8;,i+=8;]];,
"G",AppendTo[k,{i+3}];i+=3;,_,i+=#2;]]@@#&/@Group];
MapThread[{Simplify@(Delete[normMatrix.#1,k]),#2}&,Transpose[{#[[1]].FundamentalWeights,#[[2]]}&/@result]]],
True,Message[BuildRep::unknownOpts,normM,"NormalizationMethod","Canonical"];Return[{#[[1]].FundamentalWeights,#[[2]]}&/@result]],
_,Message[BuildRep::unknownOpts,outM,"OutputMethod","DynkinLabels"];Return[result]]]


(* ::Subsection::Closed:: *)
(*CGDecomposition*)


CGDecomposition::wrongArg="The argument at position `1` is not a correct multiplet.";
CGDecomposition::unknownOpts="`1` is unknonw form of option `2`. `3` is used instead.";
CGDecomposition[R1_,R2__List,opts___Rule]:=Module[{r1,r2,dimension,elder,cutDomSector,getSubMultip,subtract,domwSub,cDim,startM,addString,Rep,rNot},
(*....................................Options desdecoding and input-check..........................................................*)
{rNot}=({InputMethod}/.Flatten[{opts}])/.{InputMethod->"DynkinLabels"};
Switch[rNot,
"Weights",MapIndexed[If[PackageElemType[#1]!="multiplet",Message[CGDecomposition::wrongArg,#2];Abort[]]&,{R1,R2}];r1=DynkinLabels[R1];r2=DynkinLabels/@{R2},
"DynkinLabels",MapIndexed[If[Not[And@@({Length[#1[[1]]]==Rank,ArrayQ[#1[[1]],1,IntegerQ],NonNegative[#1[[2]]],IntegerQ[#1[[2]]]})],Message[CGDecomposition::wrongArg,#2];Abort[]]&,Flatten[{R1,R2},1]];r1=R1;r2={R2},
_,Message[CGDecomposition::unknownOpts,rNot,"InputMethod","DynkinLabels"];MapIndexed[If[Not[And@@({Length[#1[[1]]]==Rank,ArrayQ[#1[[1]],1,IntegerQ],NonNegative[#1[[2]]],IntegerQ[#1[[2]]]})],Message[CGDecomposition::wrongArg,#2];Abort[]]&,Flatten[{R1,R2},1]];r1=R1;r2={R2}];
(*....................................Addition functions and constants.....................................................*)
elder:=Function[{w1,w2},Norm[(w1+1).FundamentalWeights]>= Norm[(w2+1).FundamentalWeights]];
cutDomSector=Select[#,Function[w,And@@(NonNegative/@w[[1]])]]&;
subtract:=Function[{r1,r2},Select[Replace[r1,({#[[1]],x_}:>Evaluate[{#[[1]],x-#[[2]]}])&/@r2,1],#[[2]]!=0&]];
dimension:=Function[x,Times@@(((x.FundamentalWeights).#/WeylVector.#+1)&/@PositiveRoots)];
(*..............Program it-self........................................................................................*)
Rep=Fold[RepProduct,r1,r2];
startM=Sort[cutDomSector[Rep],elder[#1[[1]],#2[[1]]]&];
domwSub={};
Monitor[NestWhile[Function[restM,AppendTo[domwSub,restM[[1,1]]];addString[restM[[1,1]]];subtract[restM,cutDomSector[BuildRep[restM[[1,1]]]]]],startM,#!={}&],ProgressIndicator[Dynamic[Clock[Infinity]],Indeterminate]];
Print[ToString[CircleTimes@@(Total[#[[All,2]]]&/@Join[{r1},r2])]<>" = "<>ToString[CirclePlus@@(dimension/@domwSub)]];
domwSub]


(* ::Subsection::Closed:: *)
(*PhysicalNormalization*)


PhysicalNormalization::wrongArg="The argument should be weight, list of weights, or weight diagram.";
PhysicalNormalization[v_,opts___?OptionQ]:=Module[{normMatrix,k={},method},
method=(NormalizationMethod/.{opts})/.{NormalizationMethod->"Coxeter"};
If[method=="Canonical",Return[v]];
If[method!="Half"&&method!="Coxeter"&&method!="Unity",Message[PhysicalNormalizationMatrix::wrongMethod];Abort[]];
normMatrix=Module[{ph=(PhysicalNormalizationMatrix[#[[1]],#[[2]],NormalizationMethod->method]&/@Group),m},ArrayFlatten[DiagonalMatrix[Array[m,Length[ph]]]/.{m[i_]:>ph[[i]]}]];
Module[{i=0},Function[Switch[#1,"A",AppendTo[k,{i+#2+1}];i+=(#2+1);,"E",If[#2==6,AppendTo[k,i+{7}];AppendTo[k,i+{8}];i+=8;,If[#2==7,AppendTo[k,{i+8}];i+=8;,i+=8;]];,"G",AppendTo[k,{i+3}];i+=3;,_,i+=#2;]]@@#&/@Group];
Switch[PackageElemType[v],
"weight",Simplify@Delete[normMatrix.v,k],
"list of weights",Simplify@(Delete[normMatrix.#,k]&/@v),
"multiplet",MapThread[{Simplify@(Delete[normMatrix.#1,k]),#2}&,Transpose[v]],
_,Message[PhysicalNormalization::wrongArg]]]


(* ::Input:: *)
(**)


(* ::Subsection::Closed:: *)
(*Weight*)


Weight::wrongArg="The argument is neither weight, list of weights or weight diagram."
Weight::wrongArg2="The argument is GTpattern, or list of GTpatterns."
Weight::unknownOpts="`1` is unknown form of option `2`. `3` is used instead.";
Weight::wrongGT="Gelfand-Tsetlyn patterns are defined only for classical groups";
Weight[elem_,opts___]:=Module[{v,inM,pNorm,type},
{inM,pNorm}=({InputMethod,NormalizationMethod}/.Flatten[{opts}])/.{InputMethod->"DynkinLabels",NormalizationMethod->"Canonical"};
If[Not[MemberQ[{"DynkinLabels","GTPatterns"},inM]],Message[Weight::unknownOpts,inM,"InputMethod","DynkinLabels"];inM="DynkinLabels"];
If[Not[MemberQ[{"Coxeter","Half","Canonical","Unity"},pNorm]],Message[Weight::unknownOpts,pNorm,"NormalizationMethod","Canonical";pNorm="Canonical"]];
Switch[inM,"DynkinLabels",
Switch[elem,
{x__/;(And@@(IntegerQ/@{x})&&Length[{x}]==Rank)},v=elem.FundamentalWeights,
{y__/;And@@(MatchQ[#,{x__/;(And@@(IntegerQ/@{x})&&Length[{x}]==Rank)}]&/@{y})},v=(#.FundamentalWeights&/@elem),
{y__/;And@@(MatchQ[#,{{x__/;(And@@(IntegerQ/@{x})&&Length[{x}]==Rank)},_Integer}]&/@{y})},v=({#[[1]].FundamentalWeights,#[[2]]}&/@elem),
_,Message[Weight::wrongArg];Abort[]];,
"GTPatterns",If[Length[Group]!=1||Not[MemberQ[{"A","B","C","D"},Group[[1,1]]]],Message[Weight::wrongGT];Abort[]];
If[Nor[Depth[elem]==3,Depth[elem]==4],Message[Weight::wrongArg2];Abort[]];
If[Depth[elem]==3&&Not[CheckGTPattern[elem]],Message[Weight::wrongArg2];Abort[]];
If[Depth[elem]==4&&Nand@@(CheckGTPattern[#]&/@elem),Message[Weight::wrongArg2];Abort[]];
Switch[Group[[1,1]],
"A",Module[{n=Group[[1,2]]+1,M},M=Table[If[i==1,If[j==n-1,-1,0],If[i==n-j+1,2,If[i==n-j||i==n-j+2,-1,0]]],{j,n-1},{i,n}];If[Depth[elem]==3,v=M.(Total/@elem),v=(M.(Total/@#))&/@elem]],
"B",Module[{n=Group[[1,2]],M,\[Sigma]},\[Sigma]=Function[x,Insert[(Most[#]-Rest[#]),2Last[#],-1]&[Flatten[Take[x,{1,-1,2},-1]]]];M=Table[Which[j==n,Which[i==2n,4,i==2n-1,-2,True,0],i==2j-1,-1,i==2j,2,i==2j+2,-2,i==2j+3,1,True,0],{j,n},{i,2n}];If[Depth[elem]==3,v=M.(Total/@(MapIndexed[If[OddQ@@#2,Most[#1],#1]&,elem]))-\[Sigma][elem],v=(M.(Total/@(MapIndexed[If[OddQ@@#2,Most[#1],#1]&,#]))-\[Sigma][#])&/@elem]],
"C",Module[{n=Group[[1,2]],M},M=Table[Which[i==2j-1,-1,i==2j,2,i==2j+2,-2,i==2j+3,1,True,0],{j,n},{i,2n}];If[Depth[elem]==3,v=M.(Total/@elem),v=(M.(Total/@#))&/@elem]],
"D",Module[{n=Group[[1,2]],M,tTotal},tTotal=Function[y,Join[(Total/@y)+Function[x,MapIndexed[If[EvenQ@@#2,Min[x[[#2-1]],x[[#2+1]]],0]&,x]][y[[All,-1]]],{y[[-1,-1]]}]];M=Table[Which[j==n,Which[i==2n,2,i==2n-1,-2,i==2n-2,2,i==2n-3,-1,True,0],i==2j-1,-1,i==2j,2,i==2j+2,-2,i==2j+3,1,True,0],{j,n},{i,2n}];If[Depth[elem]==3,v=M.(tTotal/@elem),v=(M.(tTotal/@#))&/@elem]]];
If[Depth[elem]==3,v=v.FundamentalWeights,v=(#.FundamentalWeights&/@v)]];PhysicalNormalization[v,NormalizationMethod->pNorm]]


(* ::Subsection::Closed:: *)
(*CGCoeffecients*)


CGCoefficients::unknownOpts="`1` is unknown form of option `2`. `3` is used instead.";
CGCoefficients::noRep="There is no `1`-representation in the decomposition.";
CGCoefficients::wrongGroup="The CGCoefficents can deal only with SU-group.";
CGCoefficients[R1_,R2_,opts___]:=Module[{r1,r2,dimr1,dimr2,dimr,numberOfReps,numberList1,numberList2,numberList,listReps,r,eldestGT1,eldestGT2,eldestGTr,greater,CheckGTPatternA,w1,w2,w,cgc,Ep,Em,P1,P2,M1,M2,M,ccc,S,tNum,cNum,outM,tRep},
(*........................Input checks and options deshifration......................................*)
If[Length[Group]>1,Message[CGCoefficients::wrongGroup];Abort[]];
If[Group[[1,1]]=!="A",Message[CGCoefficients::wrongGroup];Abort[]];
{outM,tRep}=(({OutputMethod,OnlyRep}/.Flatten[{opts}])/.{OutputMethod->"Tables",OnlyRep->{}});
If[Not[MemberQ[{"List","Tables","None"},outM]],Message[CGCoefficients::unknownOpts,outM,"OutputMethod","Tables"];outM="Tables"];
PrintTemporary["Initialization..."];
(*..............................Intialization......................................................*)
r1=GTs[R1,"A"];(*The list of states for 1 rep*)
r2=GTs[R2,"A"];(*The list of states for 2 rep*)
dimr1=Length[r1];(*dimension of 1 rep*)
dimr2=Length[r2];(*dimension of 2 rep*)
eldestGT1=NestList[Most,r1[[1,1]],Rank];(*eldest pattern for 1 rep*)
eldestGT2=NestList[Most,r2[[1,1]],Rank];(*eldest pattern for 2 rep*)
Module[{lr=Sort[CGDecomposition[BuildRep[R1],BuildRep[R2]]]},
If[tRep==={},listReps=lr,listReps=Cases[lr,tRep];
If[listReps==={},Message[CGCoefficients::noRep,tRep];Return[]]]];(*list of decompositions reps*)
numberOfReps=Length[listReps];(*total number of reps*)
r=(GTs[#,"A"]&/@listReps);(*list of resulted reps*)
eldestGTr=(NestList[Most,#[[1,1]],Rank]&/@r);(*list of eldest patterns of resRep*)
dimr=(Length/@r);(*list of dimension of resRep*)
w1=DynkinLabels[Weight[r1,InputMethod->"GTPatterns"]];(*dynI of 1 rep*)
w2=DynkinLabels[Weight[r2,InputMethod->"GTPatterns"]];(*dynI of 2 rep*)
w=(DynkinLabels[Weight[#,InputMethod->"GTPatterns"]]&/@r);(*list of dynI of resRep*)
numberList1=Table[Rule[r1[[i]],i],{i,dimr1}];(*replacement list of rep1*)
numberList2=Table[Rule[r2[[i]],i],{i,dimr2}];(*replacement list of rep2*)
numberList=Table[Rule[r[[#,i]],i],{i,Length[r[[#]]]}]&/@Range[numberOfReps];(*replacement list of resRep*)
cgc=Table[Map[Table[If[w[[numr,#]]==w1[[i]]+w2[[j]],ccc[numr,#,i,j],0],{i,dimr1},{j,dimr2}]&,Range[dimr[[numr]]]],{numr,numberOfReps}];
(*the full list of Clebsch*)
tNum=Length[DeleteCases[Flatten[cgc],0]];
cNum=0;(*That is the total number of non-zeroCGCs, and number of calculated CGC, needed only for the visualization*)
(*..............................Additional functions...............................................*)
greater=Function[{x,y},Catch[If[#>0,Throw[True],If[#<0,Throw[False]]]&/@(w1[[x]]-w1[[y]]);Throw[False]]];(*Compare two weigths*)
CheckGTPatternA[p_]:=Module[{checkSnake,checkInt},
checkSnake=Function[x,Nand@@(GreaterEqual@@#&/@MapThread[Riffle[#1,#2]&,{Most[x],Rest[x]}])];
checkInt=Function[x,Nand@@((IntegerQ[#]&&NonNegative[#])&/@Flatten[x])];
If[Length/@p!=Range[Rank+1,1,-1],Return[False]];
If[checkInt[p],Return[False]];
If[checkSnake[p],Return[False]];
Return[True];];(*Check the paterns on self-consistent (the part of CheckGTPattern)*)
(*Action of the rising operator*)
Ep[Gt_,k_]:=Module[{gt=Reverse[Gt]},Sum[Module[{nnnn=Reverse[ReplacePart[gt,{k,i}->gt[[k,i]]+1]]},If[CheckGTPatternA[nnnn],S[nnnn]Sqrt[(-Product[(gt[[k+1,j]]-gt[[k,i]]+i-j),{j,1,k+1}]Product[(gt[[k-1,j]]-gt[[k,i]]+i-j-1),{j,1,k-1}])/(Product[(gt[[k,j]]-gt[[k,i]]+i-j)(gt[[k,j]]-gt[[k,i]]+i-j-1),{j,1,i-1}]Product[(gt[[k,j]]-gt[[k,i]]+i-j)(gt[[k,j]]-gt[[k,i]]+i-j-1),{j,i+1,k}])],0]],{i,1,k}]];
(*Action of the lowering operator*)
Em[Gt_,k_]:=Module[{gt=Reverse[Gt]},Sum[Module[{nnnn=Reverse[ReplacePart[gt,{k,i}->gt[[k,i]]-1]]},If[CheckGTPatternA[nnnn],S[nnnn]Sqrt[(-Product[(gt[[k+1,j]]-gt[[k,i]]+i-j+1),{j,1,k+1}]Product[(gt[[k-1,j]]-gt[[k,i]]+i-j),{j,1,k-1}])/(Product[(gt[[k,j]]-gt[[k,i]]+i-j+1)(gt[[k,j]]-gt[[k,i]]+i-j),{j,1,i-1}]Product[(gt[[k,j]]-gt[[k,i]]+i-j+1)(gt[[k,j]]-gt[[k,i]]+i-j),{j,i+1,k}])],0]],{i,1,k}]];
P1=Outer[Module[{res=Ep[#2,#1]},Coefficient[res,#]&/@(S/@r1)]&,Range[Rank],r1,1];(*Transformation matrix for rising rep1*)
P2=Outer[Module[{res=Ep[#2,#1]},Coefficient[res,#]&/@(S/@r2)]&,Range[Rank],r2,1];(*Transformation matrix for rising rep1*)
M1=Outer[Module[{res=Em[#2,#1]},Coefficient[res,#]&/@(S/@r1)]&,Range[Rank],r1,1];(*Transformation matrix for lowering rep1*)
M2=Outer[Module[{res=Em[#2,#1]},Coefficient[res,#]&/@(S/@r2)]&,Range[Rank],r2,1];(*Transformation matrix for lowering rep1*)
M=Function[x,Outer[Module[{res=Em[#2,#1]},Coefficient[res,#]&/@(S/@x)]&,Range[Rank],x,1]]/@r;(*Transformation matrix for lowering resRep*)
PrintTemporary["Done."];
Monitor[
(*..............................Eldest Clebchs calc...............................................*)
Function[numberW,Module[{nRep,LHS,\[Alpha]=Length[numberW],kkk},
nRep=First@First@Position[listReps,First[numberW],1,1];
kkk=eldestGTr[[nRep]]/.numberList[[nRep]];
LHS=DeleteCases[Flatten[Table[Sum[cgc[[nRep,kkk,i,j]](P1[[l,i,ii]]KroneckerDelta[j,jj]+KroneckerDelta[i,ii]P2[[l,j,jj]]),{i,dimr1},{j,dimr2}],{l,Rank},{ii,dimr1},{jj,dimr2}]],0];
If[\[Alpha]==1,If[LHS=={},(*if there are no equations to consider it the only cgc=1*)
Module[{i,j},{{i,j}}=Position[cgc[[nRep,kkk]],ccc[__],2,1];
cgc[[nRep,kkk,i,j]]=1;
cNum+=1;],
(*if there are several equation to consider, consider them*)
Module[{sol,vars,findConvention},
sol=Solve[((#==0)&/@LHS)~Join~{Sum[Power[cgc[[nRep,kkk,i,j]],2],{i,dimr1},{j,dimr2}]==1}];
vars=Sort[DeleteCases[Flatten[cgc[[nRep,kkk]]],0],greater[#1[[1]],#2[[1]]]&];
findConvention=Catch[FoldList[Function[{x,y},If[Length[#]==1,Throw[First[#]],#]&[Select[x,#[[y]]>=0&]]],vars/.sol,Range[Length[vars]]];Throw[First[vars/.sol]]];
cgc[[nRep,kkk]]=(cgc[[nRep,kkk]]/.(MapThread[Rule,{vars,findConvention}]));
cNum+=Length[sol[[1]]];]],
(*if \[Alpha]!=1 choose in the most stupid way*)
Module[{zeroEqn,complEqn,mainEqn,vars},
complEqn=(Sum[Power[cgc[[nRep,kkk,i,j]],2],{i,dimr1},{j,dimr2}]==1);
vars=Sort[DeleteCases[Flatten[cgc[[nRep,kkk]]],0],greater[#1[[1]],#2[[1]]]&];
zeroEqn=((#==0)&/@Take[vars,\[Alpha]]);
mainEqn=((#==0)&/@LHS);
Function[y,Module[{sol,findConvention},
sol=Solve[Join[mainEqn,ReplacePart[zeroEqn,{y}->complEqn]]];
findConvention=Catch[FoldList[Function[{x,y1},If[Length[#]==1,Throw[First[#]],#]&[Select[x,#[[y1]]>=0&]]],vars/.sol,Range[Length[vars]]];Throw[First[vars/.sol]]];
cgc[[nRep+y-1,kkk]]=(cgc[[nRep+y-1,kkk]]/.(MapThread[Rule,{ReplaceAll[vars,ccc[nRep,a_,b_,c_]:>ccc[nRep+y-1,a,b,c]],findConvention}]));(*Here was mistake*)
cNum+=Length[sol[[1]]];
]]/@Range[\[Alpha]];]]]]/@Split[listReps];
(*...........................The least coefficients..................................*)
Function[nRep,NestWhile[Function[layer,
Module[{RHS,LHS,sol,eqn},
RHS=Flatten[(Table[Sum[M[[nRep,l,#,i]]cgc[[nRep,i,k1,k2]],{i,dimr[[nRep]]}],{k1,dimr1},{k2,dimr2},{l,Rank}])&/@layer];
LHS=Flatten[(Table[Sum[cgc[[nRep,#,i,j]](M1[[l,i,k1]]KroneckerDelta[j,k2]+KroneckerDelta[i,k1]M2[[l,j,k2]]),{i,dimr1},{j,dimr2}],{k1,dimr1},{k2,dimr2},{l,Rank}])&/@layer];
eqn=DeleteDuplicates[DeleteCases[MapThread[#1==#2&,{RHS,LHS}],True]];
If[eqn!={},sol=Module[{lines=eqn[[All,1]],vars=Variables[Plus@@DeleteCases[RHS,0]]},MapThread[Rule,{vars,LinearSolve[Function[x,Coefficient[x,#]&/@vars]/@lines,eqn[[All,2]]]}]];
cgc[[nRep]]=(cgc[[nRep]]/.sol);
cNum+=Length[sol];];
DeleteDuplicates[DeleteCases[Flatten[MapIndexed[If[#1!=0,#2,0]&,(Total/@Transpose[M[[nRep,All,#]]])]&/@layer],0]]]],({eldestGTr[[nRep]]}/.numberList[[nRep]]),#!={}&];]/@Range[numberOfReps];
,TableForm[{ProgressIndicator[cNum,{1,tNum}],ToString[cNum]<>"/"<>ToString[tNum]}]];
(*........................................Formating output.............................................................*)
Unprotect[CGC,RulesForReps,RulesForStates,RulesForStates2,RulesForStates1];
Clear[CGC,RulesForReps,RulesForStates,RulesForStates2,RulesForStates1];
CGC=cgc;
RulesForReps=MapIndexed[Rule[First[#2],#1]&,listReps];
RulesForStates=MapIndexed[Rule[First[#2],#1]&,r[[#]]]&/@Range[numberOfReps];
RulesForStates1=MapIndexed[Rule[First[#2],#1]&,r1];
RulesForStates2=MapIndexed[Rule[First[#2],#1]&,r2];
Protect[CGC,RulesForReps,RulesForStates,RulesForStates2,RulesForStates1];
Switch[outM,
"None",Print["Calculation finished."],
"List",
Function[nRep,Print[TableForm[DeleteCases[Flatten[MapIndexed[If[#1!=0,{C[nRep,MatrixForm[r[[nRep,#2[[1]]]]],MatrixForm[r1[[#2[[2]]]]],MatrixForm[r2[[#2[[3]]]]]],"   ",#1}]&,cgc[[nRep]],{3}],2],Null],TableHeadings->{None,{"Representation","   ",listReps[[nRep]]}}]]]/@Range[numberOfReps];,
"Tables",
Print[Column[Table[Labeled[Framed[Row[Table[Framed[Grid[Replace[Select[Transpose[Select[Insert[MapIndexed[Insert[#1,MatrixForm[r1[[First[#2]]]],1]&,cgc[[nRep,nState]]],Insert[MatrixForm/@r2,MatrixForm[r[[nRep,nState]]],1],1],Function[x,Nand@@((#===0)&/@Rest[x])]]],Function[x,Nand@@((#===0)&/@Rest[x])]],{0->""},2],Dividers->{{2->Black},{2->Black}}],FrameMargins->20],{nState,dimr[[nRep]]}],"  "],FrameMargins->20],"Representation   "<>ToString[listReps[[nRep]]],Top,LabelStyle->{Bold,22}],{nRep,Length[listReps]}],Center,5]]
]]


(* ::Section:: *)
(*Mixed commands*)


(* ::Subsection::Closed:: *)
(*DynkinLabels*)


DynkinLabels::wrongArg="The argument is not a weight, list of weights of weight diagram.";
DynkinLabels[elem_]:=Switch[PackageElemType[elem],"weight",(2elem.#)/#.#&/@SimpleRoots,"list of weights",Outer[2 #1.#2/#2.#2&,elem,SimpleRoots,1],"multiplet",Function[x,{2 x[[1]].#/#.#&/@SimpleRoots,x[[2]]}]/@elem,"unknown",DynkinLabels::wrongArg]
DynkinLabels[elem_,sRoots_]:=Switch[PackageElemType[elem,sRoots],"weight",(2elem.#)/#.#&/@sRoots,"list of weights",Outer[2 #1.#2/#2.#2&,elem,sRoots,1],"multiplet",Function[x,{2 x[[1]].#/#.#&/@sRoots,x[[2]]}]/@elem,"unknown",DynkinLabels::wrongArg]


(* ::Subsection::Closed:: *)
(*PlotRep2D*)


PlotRep2D::nonMultiplet="The argument is not a multiplet";
PlotRep2D[Rep_,sRoots_List/;Length[sRoots]==2,opts___?OptionQ]:=
Module[{v=Orthogonalize[sRoots],points,multip,links,pStyle,lStyle,mStyle,aStyle,showM,showL,showA,l,vertx,axes},
{pStyle,lStyle,mStyle,aStyle,showM,showL,l,showA}=({PointStyle,LinkStyle,MultiplicitiesStyle,AxesStyle,ShowMultiplicities,ShowLinks,LinksList,ShowAxes}/.Flatten[{opts}]/.{PointStyle->{Medium,Black},LinkStyle->{Black},MultiplicitiesStyle->{Large,Black},AxesStyle->{Arrowheads[.02],Blue},ShowMultiplicities->True,ShowLinks->True,LinksList->PositiveRoots,ShowAxes->True});
vertx=Transpose[Rep][[1]];
points=Point[Outer[#1.#2&,vertx,v,1]];
If[showM===True,
multip=Text[Style@@Join[{If[#1[[2]]===1,"",#1[[2]]]},mStyle],{#1[[1]].v[[1]],#1[[1]].v[[2]]}+{0.1,0.1}]&/@Rep,multip=Sequence[]];
If[showL===True,
links=Module[{interV},
interV=Simplify[DeleteCases[Flatten[Outer[If[MemberQ[vertx,#1-#2],{{#1.v[[1]],#1.v[[2]]},{(#1-#2).v[[1]],(#1-#2).v[[2]]}}] &,vertx,l,1],1],Null]];
Line/@interV],links=Sequence[]];
If[showA===True,
Module[{directions=Outer[#1.#2&,sRoots,v,1],minimax},
minimax=1.2*Max[#]&/@Outer[#2.#1/#1.#1&,directions,points[[1]],1];
axes=Join[aStyle,{MapThread[{Thick,Arrow[{-#1 #2,#1 #2}]}&,{minimax,directions}],MapThread[Text[Style[#1,Large],#2 #3+{0.1,0.1}]&,{Range[2],minimax,directions}]}]],axes=Sequence[]];
Graphics[{axes,multip,Join[lStyle,links],{PointSize[pStyle[[1]]],pStyle[[2]],points}}]];


(* ::Subsection::Closed:: *)
(*PlotRep3D*)


PlotRep3D::nonMultiplet="The argument is not a multiplet";
PlotRep3D[Rep_,sRoots_List/;Length[sRoots]==3,opts___?OptionQ]:=
Module[{v=Orthogonalize[sRoots],points,multip,links,pStyle,lStyle,mStyle,aStyle,showM,showL,showA,l,vertx,axes},
{pStyle,lStyle,mStyle,aStyle,showM,showL,l,showA}=({PointStyle,LinkStyle,MultiplicitiesStyle,AxesStyle,ShowMultiplicities,ShowLinks,LinksList,ShowAxes}/.Flatten[{opts}]/.{PointStyle->{Medium,Red},LinkStyle->{Black},MultiplicitiesStyle->{Large},AxesStyle->{Arrowheads[.02],Blue},ShowMultiplicities->True,ShowLinks->True,LinksList->PositiveRoots,ShowAxes->False});
vertx=Transpose[Rep][[1]];
points=Point[Outer[#1.#2&,vertx,v,1]];
If[showM===True,
multip=Text[Style@@Join[{If[#1[[2]]===1," ",#1[[2]]]},mStyle],{#1[[1]].v[[1]],#1[[1]].v[[2]],#1[[1]].v[[3]]}+{0.1,0.1,0}]&/@Rep,multip=Sequence[]];
If[showL===True,
links=Module[{interV},
interV=DeleteCases[Flatten[Outer[If[MemberQ[vertx,#1-#2],{{#1.v[[1]],#1.v[[2]],#1.v[[3]]},{(#1-#2).v[[1]],(#1-#2).v[[2]],(#1-#2).v[[3]]}}] &,vertx,l,1],1],Null];
Line/@interV],links=Sequence[]];
If[showA===True,
Module[{directions=Outer[#1.#2&,sRoots,v,1],minimax},
minimax=1.2*Max[#]&/@Outer[#2.#1/#1.#1&,directions,points[[1]],1];
axes=Join[aStyle,{MapThread[{Thick,Arrow[{-#1 #2,#1 #2}]}&,{minimax,directions}],
MapThread[Text[Style[#1,Large],#2 #3+{0.1,0.1,0.1}]&,{Range[3],minimax,directions}]}]],
axes=Sequence[]];
Graphics3D[{multip,axes,Join[lStyle,links],PointSize[pStyle[[1]]],pStyle[[2]],points},Boxed->False]];


(* ::Input:: *)
(**)


(* ::Section:: *)
(*End Of Package*)


{SimpleRoots,PositiveRoots,FundamentalWeights,MaximalRoots,Rank,CoxeterNumber,Group}={{},{},{},{},0,0,{}};
{CGC,RulesForReps,RulesForStates,RulesForStates2,RulesForStates1}={{},{},{},{},{}};
Protect[SimpleRoots,PositiveRoots,FundamentalWeights,MaximalRoots,Rank,CoxeterNumber,Group];
Protect[CGC,RulesForReps,RulesForStates,RulesForStates2,RulesForStates1];
Protect[SetGroup,CartanMatrix,FundamentalWeight,RootSystem,Representation,RepProduct,DecomposeRep,PhysicalNormalizationMatrix,DynkinDiagram,CoxeterPlane,BuildRep,CGDecomposition,PhysicalNormalization,RepDimension,DynkinLabels,PlotRep2D,PlotRep3D,Weight];


Print["Package SimpleGroups (v.0.98\[Beta] (date:06.05.10)) is loaded"];
Print["Use SetGroup[] to set the group to work with."];


End[]
EndPackage[]
