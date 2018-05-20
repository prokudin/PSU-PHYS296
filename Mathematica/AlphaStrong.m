(* ::Package:: *)

Print["Package AlphaStrong for evaluation of as[\[Mu]]=\!\(\*FractionBox[SuperscriptBox[\(g\), \(2\)], SuperscriptBox[\((4\\\ \[Pi])\), \(2\)]]\), at 1-,2-,3- loops."];
Print["Contains the following functions: as[n,\[Mu]], asIt[n,\[Mu]], \[CapitalLambda]QCD[n,\[Mu]], ActiveN[\[Mu]], As[n,m,\[Mu],\[Mu]0], AsLog[n,m,\[Mu],\[Mu]0]"];
Print["Contains the following anomalous dimensions: \[Beta]0[Nf], \[Beta]1[Nf], \[Beta]2[Nf], \[Beta]3[Nf]"];
Print["Contains the following constants: mBOTTOM, mCHARM"];
Print["Copyright: Alexey Vladimirov, Regensburg University. Version 1.101 (1 Sep. 2016)"];
Print["___________________________________________________________________________"];

BeginPackage["AlphaStrong`"];
\[Beta]0::usage="\[Beta]0[Nf_] the 1-loop \[Beta]-function  RGE defintion: a'=-\[Beta]0 \!\(\*SuperscriptBox[\(a\), \(2\)]\)-\[Beta]1 \!\(\*SuperscriptBox[\(a\), \(3\)]\)-\[Beta]2 \!\(\*SuperscriptBox[\(a\), \(4\)]\)-\[Beta]3 \!\(\*SuperscriptBox[\(a\), \(5\)]\)+.. with \[Beta]0=\!\(\*FractionBox[\(11\), \(3\)]\)CA-\!\(\*FractionBox[\(2\), \(3\)]\)Nf";
\[Beta]1::usage="\[Beta]1[Nf_] the 1-loop \[Beta]-function  RGE defintion: a'=-\[Beta]0 \!\(\*SuperscriptBox[\(a\), \(2\)]\)-\[Beta]1 \!\(\*SuperscriptBox[\(a\), \(3\)]\)-\[Beta]2 \!\(\*SuperscriptBox[\(a\), \(4\)]\)-\[Beta]3 \!\(\*SuperscriptBox[\(a\), \(5\)]\)+.. with  \[Beta]0=\!\(\*FractionBox[\(11\), \(3\)]\)CA-\!\(\*FractionBox[\(2\), \(3\)]\)Nf";
\[Beta]2::usage="\[Beta]2[Nf_] the 1-loop \[Beta]-function  RGE defintion: a'=-\[Beta]0 \!\(\*SuperscriptBox[\(a\), \(2\)]\)-\[Beta]1 \!\(\*SuperscriptBox[\(a\), \(3\)]\)-\[Beta]2 \!\(\*SuperscriptBox[\(a\), \(4\)]\)-\[Beta]3 \!\(\*SuperscriptBox[\(a\), \(5\)]\)+.. with  \[Beta]0=\!\(\*FractionBox[\(11\), \(3\)]\)CA-\!\(\*FractionBox[\(2\), \(3\)]\)Nf";
\[Beta]3::usage="\[Beta]3[Nf_] the 1-loop \[Beta]-function  RGE defintion: a'=-\[Beta]0 \!\(\*SuperscriptBox[\(a\), \(2\)]\)-\[Beta]1 \!\(\*SuperscriptBox[\(a\), \(3\)]\)-\[Beta]2 \!\(\*SuperscriptBox[\(a\), \(4\)]\)-\[Beta]3 \!\(\*SuperscriptBox[\(a\), \(5\)]\)+.. with  \[Beta]0=\!\(\*FractionBox[\(11\), \(3\)]\)CA-\!\(\*FractionBox[\(2\), \(3\)]\)Nf";
\[CapitalLambda]QCDit::usage="\[CapitalLambda]QCDit[n_Integer,\[Mu]_Number]  returns \[CapitalLambda]QCD for n-loops, \[Mu]-scale obtained by double-iterative solution.
as[MZ]=asMSTW2008. + threashold matching.
";
\[CapitalLambda]QCD::usage="\[CapitalLambda]QCD[n_Integer,\[Mu]_Number]  returns \[CapitalLambda]QCD for n-loops, \[Mu]-scale obtained by HEP solution.
as[MZ]=asMSTW2008. + threashold matching.";
ActiveN::usage="ActiveN[\[Mu]_] the number of the active flavors, take into account threasholds mb and mb"
asIt::usage="asIt[n_Integer,\[Mu]_Number]  returns as=\!\(\*SuperscriptBox[\(g\), \(2\)]\)/(4 \[Pi])^2 for n-loops and \[Mu]-scale.
Here the double iteractive solution is used. The double iterative solution is very close to an exact solution, but it give non-common \[CapitalLambda]'s, since the standard fits uses the HEP solution.
asIt[MZ]=asMSTW2008. + threashold matching.";
as::usage="as[n_Integer,\[Mu]_Number]  returns as=\!\(\*SuperscriptBox[\(g\), \(2\)]\)/(4 \[Pi])^2 for n-loops and \[Mu]-scale.
Here the (standard) HEP solution is used. This solution is rather weakly approximate an exact solution, but it is commonly used. Since it is the standard, it is preferable, and better normalized.
as[MZ]=asMSTW2008. + threashold matching.";
mBOTTOM::usage="mass of bottom quark (from MSTW2008)";
mCHARM::usage="mass of charm quark (from MSTW2008)";
As::usage="As[n_,m_,  \[Mu]_, \[Mu]0_]  \!\(\*FormBox[\(is\\\ the\\\ function\\\ defined\\\ as\\\ \(\*OverscriptBox[SubscriptBox[\(\[Integral]\), \(\[Mu]0\)], \(\[Mu]\)]\*FractionBox[\(d\[Nu]\), \(\[Nu]\)] \*SuperscriptBox[\(as[m, \[Nu]]\), \(n\)]\)\),
TraditionalForm]\).
It uses default definition of as.
IMPORTANT: MEMORY-CONSUMING the function self-memorize the evaluated values. We have done it in order to speed up calculation, since these values are often-used functions. The procedure it sefl is not very memory-consuming, basically every point is 4-real numbers.
None the less we recommend to evaluate functions controlably, i.e. point-by-point, since uncontrolable evaluation (e.g. functions like Plot) can generate many useless data. It is also convinient to use Share[As] from time to time, it reduces the memory usage.
In the case of sirious memory-leak coused by function As use Clear[As], and restart  package.";
AsLog::usage="AsLog[n_,m_,  \[Mu]_, \[Mu]0_]  \!\(\*FormBox[\(is\\\ the\\\ function\\\ defined\\\ as\\\ \(\*OverscriptBox[SubscriptBox[\(\[Integral]\), \(\[Mu]0\)], \(\[Mu]\)]\*FractionBox[\(d\[Nu]\), \(\[Nu]\)] \*SuperscriptBox[\(as[m, \[Nu]]\), \(n\)]\)\),
TraditionalForm]\)Log[\!\(\*FractionBox[SuperscriptBox[\(\[Nu]\), \(2\)], SuperscriptBox[\(\[Mu]0\), \(2\)]]\)].
It uses default definition of as.
IMPORTANT: MEMORY-CONSUMING the function self-memorize the evaluated values. We have done it in order to speed up calculation, since these values are often-used functions. The procedure it sefl is not very memory-consuming, basically every point is 4-real numbers.
None the less we recommend to evaluate functions controlably, i.e. point-by-point, since uncontrolable evaluation (e.g. functions like Plot) can generate many useless data. It is also convinient to use Share[AsLog] from time to time, it reduces the memory usage.
In the case of sirious memory-leak coused by function As use Clear[AsLog], and restart  package.";


Begin["Private`"];


\[Beta]0[Nf_]:=\[Beta]0[Nf]=11-(2 Nf)/3
\[Beta]1[Nf_]:=\[Beta]1[Nf]=102-(38 Nf)/3
\[Beta]2[Nf_]:=\[Beta]2[Nf]=2857/2-(5033 Nf)/18+(325 Nf^2)/54
\[Beta]3[Nf_]:=\[Beta]3[Nf]=149753/6+(1093 Nf^3)/729+3564 Zeta[3]+Nf^2 (50065/162+(6472 Zeta[3])/81)-Nf (1078361/162+(6508 Zeta[3])/27)


mBOTTOM=4.75;
mCHARM=1.4;


ActiveN[\[Mu]_?NumberQ]:=Which[\[Mu]>=4.75`,5,\[Mu]>=1.4`,4,True,3]


\[CapitalLambda]QCDit[n_?IntegerQ,\[Mu]_?NumberQ]:=\[CapitalLambda]QCDitCompiled[n,\[Mu]]
\[CapitalLambda]QCDitCompiled=Compile[{{n,_Integer},{\[Mu],_Real}},
Which[n==1,
Which[\[Mu]>=4.75`,0.08311262598220348`,\[Mu]>=1.4`,0.11487595443310693`,True,0.13825030193630825`],
n==2,
Which[\[Mu]>=4.75`,0.23021053364096944`,\[Mu]>=1.4`,0.3240595819790062`,True,0.37909977715127163`],
n==3,
Which[\[Mu]>=4.75`,0.2737424902051772`,\[Mu]>=1.4`,0.4087221166574235`,True,0.4810477102019111`]]
];


\[CapitalLambda]QCD[n_?IntegerQ,\[Mu]_?NumberQ]:=\[CapitalLambda]QCDCompiled[n,\[Mu]]
\[CapitalLambda]QCDCompiled=Compile[{{n,_Integer},{\[Mu],_Real}},Which[n==1,
Which[\[Mu]>=4.75`,0.08311262598220348`,\[Mu]>=1.4`,0.11487595443310693`,True,0.13825030193630825`],
n==2,
Which[\[Mu]>=4.75`,0.214635561673419`,\[Mu]>=1.4`,0.3127952142217667`,True,0.3630391178804666`],
n==3,
Which[\[Mu]>=4.75`,0.19768896276830034`,\[Mu]>=1.4`,0.27783484947934445`,True,0.3200215629198203`]]
];


asIt[n_?IntegerQ,\[Mu]_?NumberQ]:=asItCompiled[n,\[Mu]]
asItCompiled=Compile[{{n,_Integer},{\[Mu],_Real}},Module[{\[CapitalLambda],B1,B2,ZZ,Nf},
Nf=Which[\[Mu]>=4.75`,5,\[Mu]>=1.4`,4,True,3];
B1=(306-38 Nf)/(33-2 Nf);
B2=(77139-15099 Nf+325 Nf^2)/(594-36 Nf);
\[CapitalLambda]=Which[n==1,
Which[\[Mu]>=4.75`,0.08311262598220348`,\[Mu]>=1.4`,0.11487595443310693`,True,0.13825030193630825`],
n==2,
Which[\[Mu]>=4.75`,0.23021053364096944`,\[Mu]>=1.4`,0.3240595819790062`,True,0.37909977715127163`],
n==3,
Which[\[Mu]>=4.75`,0.2737424902051772`,\[Mu]>=1.4`,0.4087221166574235`,True,0.4810477102019111`]];
ZZ=(22-(4 Nf)/3)Log[\[Mu]/\[CapitalLambda]];
If[\[Mu]<\[CapitalLambda],Return[0.]];
Which[n==1,1/ZZ,
n==2,1/(ZZ+B1 Log[1+ZZ/B1+Log[1+ZZ/B1]]),
n==3,With[{\[Rho]=ZZ+B1/2 Log[1+B1/B2 ZZ+ZZ^2/B2]-(B1^2-2B2)/Sqrt[4 B2-B1^2] (ArcTan[(Sqrt[-B1^2+4 B2] ZZ)/(2 B2+B1 ZZ)])},
(ZZ+B1/2 Log[1+B1/B2 \[Rho]+\[Rho]^2/B2]-(B1^2-2B2)/Sqrt[4 B2-B1^2] (ArcTan[(Sqrt[-B1^2+4 B2] \[Rho])/(2 B2+B1 \[Rho])]))^-1
]]]];


as[n_?IntegerQ,\[Mu]_?NumberQ]:=asCompiled[n,\[Mu]]
asCompiled=Compile[{{n,_Integer},{\[Mu],_Real}},Module[{\[CapitalLambda],B1,B2,ZZ,Nf,LL},
\[CapitalLambda]=Which[n==1,
Which[\[Mu]>=4.75`,0.08311262598220348`,\[Mu]>=1.4`,0.11487595443310693`,True,0.13825030193630825`],
n==2,
Which[\[Mu]>=4.75`,0.214635561673419`,\[Mu]>=1.4`,0.3127952142217667`,True,0.3630391178804666`],
n==3,
Which[\[Mu]>=4.75`,0.19768896276830034`,\[Mu]>=1.4`,0.27783484947934445`,True,0.3200215629198203`]];
If[\[Mu]<\[CapitalLambda],Return[0.]];
Nf=Which[\[Mu]>=4.75`,5,\[Mu]>=1.4`,4,True,3];
B1=(306-38 Nf)/(33-2 Nf);
B2=(77139-15099 Nf+325 Nf^2)/(594-36 Nf);
ZZ=(22-(4 Nf)/3)Log[\[Mu]/\[CapitalLambda]];
LL=Log[2Log[\[Mu]/\[CapitalLambda]]];
Which[n==1,1/ZZ,
n==2,1/ZZ (1-B1/ZZ LL),
n==3,1/ZZ (1-B1/ZZ LL)+1/ZZ^3 (B1^2 (LL^2-LL-1)+B2)
]]
];


AsNUM=Compile[{{n,_Integer},{m,_Integer},{\[Mu],_Real},{\[Mu]0,_Real}},NIntegrate[as[m,\[Nu]]^n/\[Nu],{\[Nu],\[Mu]0,\[Mu]}]];
AsLogNUM=Compile[{{n,_Integer},{m,_Integer},{\[Mu],_Real},{\[Mu]0,_Real}},NIntegrate[2 as[m,\[Nu]]^n/\[Nu] Log[\[Nu]/\[Mu]0],{\[Nu],\[Mu]0,\[Mu]}]];


As[n_?IntegerQ,m_?IntegerQ,\[Mu]_?NumberQ,\[Mu]0_?NumberQ]:=As[n,m,\[Mu],\[Mu]0]=AsNUM[n,m,\[Mu],\[Mu]0]
AsLog[n_?IntegerQ,m_?IntegerQ,\[Mu]_?NumberQ,\[Mu]0_?NumberQ]:=AsLog[n,m,\[Mu],\[Mu]0]=AsLogNUM[n,m,\[Mu],\[Mu]0]


End[];
EndPackage[];
