(* ::Package:: *)

Print["Package TMDevolution contains the set of functions related to the TMD-evolution kernel R"];
Print["Load package AlphaStrong"];
Print["Contains the following functions: Cf, bstar, Dboundary, \[CapitalGamma]evolutor, \[CapitalGamma]Logevolutor, \[Gamma]evolutor, RTMD, RTMD2"];
Print["Contains the following anomalous dimensions: \[CapitalGamma]cusp, \[CapitalGamma]c, \[Gamma]V, ddTMD, dd"];
Print["Copyright: Alexey Vladimirov, Regensburg University. Version 1.501 (1 Sep. 2016)"];
Print["___________________________________________________________________________"];

BeginPackage["TMDevolution`",{"AlphaStrong`"}];
C0::usage="C0 the constant often which appear in TMD business.
C0=2 Exp[-\!\(\*SubscriptBox[\(\[Gamma]\), \(E\)]\)]";
Cf::usage="Cf[f_] f is a flavour. It can be \"q\" for quark and anti-quark, or \"g\" for gluon.
Cf[q]=CF=4/3,  Cf[g]=CA=3";
\[CapitalGamma]cusp::usage="\[CapitalGamma]cusp[f_,n_,Nf_] the coefficeint of the cusp anomalous dimension at \!\(\*SuperscriptBox[\(as\), \(n\)]\).
Nf is the number of active quarks.
f is a flavour. It can be \"q\" for quark and anti-quark, or \"g\" for gluon";
\[CapitalGamma]c::usage="\[CapitalGamma]c[n_,Nf_] the flavour independent part of \!\(\*SuperscriptBox[\(as\), \(n\)]\) coefficient of the cusp anomalous dimension.
It is defined as \[CapitalGamma]cusp= 4 \!\(\*SubscriptBox[\(C\), \(f\)]\) \[CapitalGamma]c.
Nf is the number of active quarks.
The values are give up to 3 loops. 4-loop expression is the Pade.";
\[Gamma]V::usage="\[Gamma]V[f_,n_,Nf_]  the part coefficeint of the TMD anomalous dimension at \!\(\*SuperscriptBox[\(as\), \(n\)]\).
Anomalous dimension defined as \!\(\*FractionBox[\(d\\\ F\), \(d\\\ ln\\\ \*SuperscriptBox[\(\[Mu]\), \(2\)]\)]\)=\!\(\*FractionBox[\(\[Gamma]\), \(2\)]\)F,  where \[Gamma]=\[CapitalGamma]cusp ln(\!\(\*SuperscriptBox[\(mu\), \(2\)]\)/\[Zeta])-\[Gamma]V. This defines \[Gamma]V.
Nf is the number of active quarks.
f is a flavour. It can be \"q\" for quark and anti-quark, or \"g\" for gluon.";
ddTMD::usage="dTMD[f_,n_,k_,Nf_] the coefficeint for rapidity anomalous dimension D at \!\(\*SuperscriptBox[\(as\), \(n\)]\)\!\(\*SuperscriptBox[\(L\[Mu]\), \(k\)]\).
f is a flavour. It can be \"q\" for quark and anti-quark, or \"g\" for gluon.
It is defined as (see also dd) \!\(\*FractionBox[\(d\\\ F\), \(d\\\ ln\\\ \[Zeta]\)]\)=- D F, where D=\!\(\*SubscriptBox[\(C\), \(f\)]\) dd.";
dd::usage="dd[n_,k_,Nf_]  flavour independent the coefficeint for rapidity anomalous dimension D at \!\(\*SuperscriptBox[\(as\), \(n\)]\)\!\(\*SuperscriptBox[\(L\[Mu]\), \(k\)]\).
It is defined as \!\(\*FractionBox[\(d\\\ F\), \(d\\\ ln\\\ \[Zeta]\)]\)=- D F, where D=\!\(\*SubscriptBox[\(C\), \(f\)]\) dd.
Nf is the number of active quarks.";
\[CapitalGamma]evolutor::usage="\[CapitalGamma]evolutor[f_,n_,\[Mu]_,\[Mu]0_] the evolution integral defined as
\!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(\[Mu]0\)], \(\[Mu]\)]\)\!\(\*FractionBox[\(d\[Nu]\), \(\[Nu]\)]\)\[CapitalGamma][as[\[Nu]]] at n-loop order.
f is for flavor of D.
USES As";
Dboundary::usage="Dboundary[f_,n_,\[Mu]_,bt_] pertrubative value of the rapirity anomalous dimension D. 
Used as a boundary for anomalous dimension D at \[Mu]0.";
bstar::usage="bstar[bt_Real,bmax_Real] = \!\(\*FractionBox[\(bt\), \(Sqrt[1 + \*FractionBox[SuperscriptBox[\(bt\), \(2\)], SuperscriptBox[\(bmax\), \(2\)]]]\)]\)";
DTMD::usage="DTMD[f_,n_,\[Mu]_,bt_,\[Mu]0_] the rapirity anomalous dimension D. 
Defined as \!\(\*FractionBox[\(d\\\ F\), \(d\\\ ln\\\ \[Zeta]\)]\)=- D F. DTMD=\[CapitalGamma]evolutor+Dboundary.
n is order of evaluation, f is flavor, \[Mu] is high-energy parameter, bt is bt, \[Mu]0 is low-energy parameter (see documentatio)";
\[CapitalGamma]Logevolutor::usage="\[CapitalGamma]Logevolutor[f_,n_,\[Mu]_,\[Mu]0_] the evolution integral for double-logs.
\!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(\[Mu]0\)], \(\[Mu]\)]\)\!\(\*FractionBox[\(d\[Nu]\), \(\[Nu]\)]\)\[CapitalGamma][as[\[Nu]]] Log[\!\(\*FractionBox[SuperscriptBox[\(\[Nu]\), \(2\)], SuperscriptBox[\(\[Mu]0\), \(2\)]]\)] at n-loop order. 
f is for flavor.
USES AsLog";
\[Gamma]evolutor::usage="\[Gamma]evolutor[f_,n_,\[Mu]_,\[Mu]0_] the evolution integral for TMD anomalous dimension \[Gamma]V.
\!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(\[Mu]0\)], \(\[Mu]\)]\)\!\(\*FractionBox[\(d\[Nu]\), \(\[Nu]\)]\)\[Gamma]V[as[\[Nu]]] at n-loop order.
f is for flavor.
USES As.";
RTMD::usage="RTMD[f_,n_,bt_,\[Zeta]f_,\[Mu]f_,\[Zeta]i_,\[Mu]i_,\[Mu]0_] the TMD evolution kernel. 
Defined as Exp[\!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(\[Mu]i\)], \(\[Mu]f\)]\)\!\(\*FractionBox[\(d\[Mu]\), \(\[Mu]\)]\)\[Gamma][\[Mu],\[Zeta]f]]Exp[-D[\[Mu]i,bt]Log[\[Zeta]f/\[Zeta]i].
The parameter f(\"q\" of \"g\") is flavour.
The parameter \[Mu]0 is low-energy behaviour of function D.
The following verions of function exists:
RTMD[f,n,bt,\[Zeta]f,\[Mu]f,\[Mu]i,\[Mu]0] \[Zeta]i=\[Zeta]b=\!\(\*SuperscriptBox[\(C0\), \(2\)]\)/\!\(\*SuperscriptBox[\(bt\), \(2\)]\) 
RTMD[f,n,bt,\[Zeta]f,\[Mu]f,\[Mu]0] \[Zeta]i=\[Zeta]b=\!\(\*SuperscriptBox[\(C0\), \(2\)]\)/\!\(\*SuperscriptBox[\(bt\), \(2\)]\), \[Mu]i=\[Mu]0";
RTMD2::usage="RTMD2[f_,pn_,bt_,\[Zeta]f_,\[Mu]f_,\[Zeta]i_,\[Mu]i_,\[Mu]0_] the TMD evolution kernel over the alternative contour. 
Defined as Exp[\!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(\[Mu]i\)], \(\[Mu]f\)]\)\!\(\*FractionBox[\(d\[Mu]\), \(\[Mu]\)]\)\[Gamma][\[Mu],\[Zeta]i]]Exp[-D[\[Mu]f,bt]Log[\[Zeta]f/\[Zeta]i].
The parameter f(\"q\" of \"g\") is flavour.
The parameter \[Mu]0 is low-energy behaviour of function D.
The following verions of function exists:
RTMD2[f,n,bt,\[Zeta]f,\[Mu]f,\[Mu]i,\[Mu]0] \[Zeta]i=\[Zeta]b=\!\(\*SuperscriptBox[\(C0\), \(2\)]\)/\!\(\*SuperscriptBox[\(bt\), \(2\)]\) 
RTMD2[f,n,bt,\[Zeta]f,\[Mu]f,\[Mu]0] \[Zeta]i=\[Zeta]b=\!\(\*SuperscriptBox[\(C0\), \(2\)]\)/\!\(\*SuperscriptBox[\(bt\), \(2\)]\)";


Begin["`Private`"];


C0=2. Exp[-EulerGamma];


Cf::wrongf="The argument for flavour should be \"f\" or \"g\", but not `1`";
Cf["q"]=4/3;
Cf["g"]=3;
Cf[f_]:=Message[Cf::wrongf,f];


bstar=Compile[{{bt,_Real},{bmax,_Real}},bt/Sqrt[1+bt^2/bmax^2]];


\[CapitalGamma]cusp::wrongf="The argument for flavour should be \"f\" or \"g\", but not `1`";
\[CapitalGamma]cusp["q",n_,Nf_]:=16/3 \[CapitalGamma]c[n,Nf]
\[CapitalGamma]cusp["g",n_,Nf_]:=12 \[CapitalGamma]c[n,Nf]
\[CapitalGamma]cusp[f_,n_,Nf_]:=Message[\[CapitalGamma]cusp::wrongf,f]


\[CapitalGamma]c[1,Nf_]:=\[CapitalGamma]c[1,Nf]=1
\[CapitalGamma]c[2,Nf_]:=\[CapitalGamma]c[2,Nf]= 67/3-\[Pi]^2-10/9  Nf
\[CapitalGamma]c[3,Nf_]:=\[CapitalGamma]c[3,Nf]=735/2-(4 Nf^2)/27-(134 \[Pi]^2)/3+(11 \[Pi]^4)/5+1/9 Nf (-319+20 \[Pi]^2-156 Zeta[3])+66 Zeta[3]
\[CapitalGamma]c[4,Nf_]:=\[CapitalGamma]c[4,Nf]=3/16 Which[Nf== 3, 7849,Nf== 4,4313,Nf== 5,1553, True,0]


\[Gamma]V::wrongf="The argument for flavour should be \"f\" or \"g\", but not `1`";
\[Gamma]V["q",1,Nf_]:=\[Gamma]V["q",1,Nf]=-8
\[Gamma]V["q",2,Nf_]:=\[Gamma]V["q",2,Nf]=4/81 (2 Nf (65+9 \[Pi]^2)-3 (997+51 \[Pi]^2-828 Zeta[3]))
\[Gamma]V["q",3,Nf_]:=\[Gamma]V["q",3,Nf]=(8 Nf (45295+17115 \[Pi]^2+51 \[Pi]^4-84060 Zeta[3]))/3645-(8 Nf^2 (-2417+270 \[Pi]^2+216 Zeta[3]))/2187-(2 (983105+6198 \[Pi]^4-3693240 Zeta[3]+30 \[Pi]^2 (2531+2440 Zeta[3])+2069280 Zeta[5]))/1215
\[Gamma]V["g",1,Nf_]:=\[Gamma]V["g",1,Nf]=-22+(4 Nf)/3
\[Gamma]V["g",2,Nf_]:=\[Gamma]V["g",2,Nf]=-(1384/3)+11 \[Pi]^2+Nf (304/9-(2 \[Pi]^2)/3)+36 Zeta[3]
\[Gamma]V["g",3,Nf_]:=\[Gamma]V["g",3,Nf]=-(194372/27)-(319 \[Pi]^4)/5+1/243 Nf^2 (-1853+180 \[Pi]^2-3024 Zeta[3])+\[Pi]^2 (6109/9-120 Zeta[3])+2196 Zeta[3]+1/405 Nf (298175-19050 \[Pi]^2+1926 \[Pi]^4+41400 Zeta[3])-864 Zeta[5]
\[Gamma]V[f_,n_,Nf_]:=Message[\[Gamma]V::wrongf,f]


ddTMD[f_,n_,k_,Nf_]:=Cf[f]dd[n,k,Nf]


dd[1,1,Nf_]:=dd[1,1,Nf]=2\[CapitalGamma]c[1,Nf]
dd[1,0,Nf_]:=dd[1,0,Nf]=0
dd[2,2,Nf_]:=dd[2,2,Nf]=\[CapitalGamma]c[1,Nf]\[Beta]0[Nf]
dd[2,1,Nf_]:=dd[2,1,Nf]=2 \[CapitalGamma]c[2,Nf]
dd[2,0,Nf_]:=dd[2,0,Nf]=3(404/27-14 Zeta[3])-112/27 Nf/2
dd[3,3,Nf_]:=dd[3,3,Nf]=2 \[CapitalGamma]c[1,Nf]/3 \[Beta]0[Nf]^2
dd[3,2,Nf_]:=dd[3,2,Nf]=2\[CapitalGamma]c[2,Nf]\[Beta]0[Nf]+\[CapitalGamma]c[1,Nf]\[Beta]1[Nf]
dd[3,1,Nf_]:=dd[3,1,Nf]=2\[Beta]0[Nf]dd[2,0,Nf]+2 \[CapitalGamma]c[3,Nf]
dd[3,0,Nf_]:=dd[3,0,Nf]=297029/162-(77 \[Pi]^4)/30-(6164 Zeta[3])/3+\[Pi]^2 (-(1598/27)+44 Zeta[3])+16/729 Nf^2 (58+81 Zeta[3])+(Nf (-207895+3090 \[Pi]^2+9 \[Pi]^4+88380 Zeta[3]))/1215+864 Zeta[5]


\[CapitalGamma]evolutor[f_,n_,\[Mu]_,\[Mu]0_]:=(4Cf[f]Which[
(\[Mu]0<=mCHARM && \[Mu]<=mCHARM)||(mCHARM<\[Mu]0<=mBOTTOM && mCHARM<\[Mu]<=mBOTTOM)||(mBOTTOM<\[Mu]0 && mBOTTOM<\[Mu]),Sum[\[CapitalGamma]c[k,ActiveN[\[Mu]]]As[k,n,\[Mu],\[Mu]0],{k,1,n}],
\[Mu]0<=mCHARM && mCHARM<\[Mu]<=mBOTTOM,Sum[\[CapitalGamma]c[k,4]As[k,n,\[Mu],mCHARM]+\[CapitalGamma]c[k,3]As[k,n,mCHARM,\[Mu]0],{k,1,n}],
\[Mu]0<=mCHARM && mBOTTOM<\[Mu],Sum[\[CapitalGamma]c[k,5]As[k,n,\[Mu],mBOTTOM]+\[CapitalGamma]c[k,4]As[k,n,mBOTTOM,mCHARM]+\[CapitalGamma]c[k,3]As[k,n,mCHARM,\[Mu]0],{k,1,n}],
mCHARM<\[Mu]0<=mBOTTOM && mBOTTOM<=\[Mu],Sum[\[CapitalGamma]c[k,5]As[k,n,\[Mu],mBOTTOM]+\[CapitalGamma]c[k,4]As[k,n,mBOTTOM,\[Mu]0],{k,1,n}],
\[Mu]<=mCHARM && mCHARM<\[Mu]0<=mBOTTOM,Sum[-\[CapitalGamma]c[k,3]As[k,n,mCHARM,\[Mu]]-\[CapitalGamma]c[k,4]As[k,n,\[Mu]0,mCHARM],{k,1,n}],\[Mu]<=mCHARM && mBOTTOM<\[Mu]0,Sum[-\[CapitalGamma]c[k,3]As[k,n,mCHARM,\[Mu]]+\[CapitalGamma]c[k,4]As[k,n,mBOTTOM,mCHARM]-\[CapitalGamma]c[k,5]As[k,n,\[Mu]0,mBOTTOM],{k,1,n}],
mCHARM<\[Mu]<=mBOTTOM && mBOTTOM<=\[Mu]0,Sum[-\[CapitalGamma]c[k,4]As[k,n,mBOTTOM,\[Mu]]-\[CapitalGamma]c[k,5]As[k,n,\[Mu]0,mBOTTOM],{k,1,n}]])/;(n===1||n===2||n===3)


Dboundary[f_,n_,\[Mu]_,bt_]:=Cf[f]If[bt==0,-10^9,DboundaryNUM[n,as[n,\[Mu]],Log[(bt^2 \[Mu]^2)/C0^2]]]


DboundaryNUM[n_,a_,Lb_]:=Which[
n==1,a(dd[1,1,3]Lb),
n==2,a(dd[1,1,3]Lb)+a^2 (dd[2,2,3]Lb^2+dd[2,1,3]Lb+dd[2,0,3]),
n==3,a(dd[1,1,3]Lb)+a^2 (dd[2,2,3]Lb^2+dd[2,1,3]Lb+dd[2,0,3])+a^3 (dd[3,3,3]Lb^3+dd[3,2,3]Lb^2+dd[3,1,3]Lb+dd[3,0,3])]


DTMD[f_,n_,\[Mu]_,bt_,\[Mu]0_]:=Dboundary[f,n,\[Mu]0,bt]+\[CapitalGamma]evolutor[f,n,\[Mu],\[Mu]0]


\[CapitalGamma]Logevolutor[f_,n_,\[Mu]_,\[Mu]0_]:=(Which[
(\[Mu]0<=mCHARM && \[Mu]<=mCHARM)||(mCHARM<\[Mu]0<=mBOTTOM && mCHARM<\[Mu]<=mBOTTOM)||(mBOTTOM<\[Mu]0 && mBOTTOM<\[Mu]),4Cf[f]Sum[\[CapitalGamma]c[k,ActiveN[\[Mu]]]AsLog[k,n,\[Mu],\[Mu]0],{k,1,n}],
\[Mu]0<=mCHARM && mCHARM<\[Mu]<=mBOTTOM,4Cf[f]Sum[\[CapitalGamma]c[k,4]AsLog[k,n,\[Mu],mCHARM]+\[CapitalGamma]c[k,3]AsLog[k,n,mCHARM,\[Mu]0]+2\[CapitalGamma]c[k,4]As[k,n,\[Mu],mCHARM]Log[mCHARM/\[Mu]0],{k,1,n}],
\[Mu]0<=mCHARM && mBOTTOM<\[Mu],4Cf[f]Sum[\[CapitalGamma]c[k,5]AsLog[k,n,\[Mu],mBOTTOM]+\[CapitalGamma]c[k,4]AsLog[k,n,mBOTTOM,mCHARM]+\[CapitalGamma]c[k,3]AsLog[k,n,mCHARM,\[Mu]0]+2\[CapitalGamma]c[k,5]As[k,n,\[Mu],mBOTTOM]Log[mBOTTOM/\[Mu]0]+2\[CapitalGamma]c[k,4]As[k,n,mBOTTOM,mCHARM]Log[mCHARM/\[Mu]0],{k,1,n}],
mCHARM<\[Mu]0<=mBOTTOM && mBOTTOM<=\[Mu],4Cf[f]Sum[\[CapitalGamma]c[k,5]AsLog[k,n,\[Mu],mBOTTOM]+\[CapitalGamma]c[k,4]AsLog[k,n,mBOTTOM,\[Mu]0]+2\[CapitalGamma]c[k,5]As[k,n,\[Mu],mBOTTOM]Log[mBOTTOM/\[Mu]0],{k,1,n}],
\[Mu]==\[Mu]0,0,
\[Mu]<\[Mu]0,-\[CapitalGamma]Logevolutor[f,n,\[Mu]0,\[Mu]]+2\[CapitalGamma]evolutor[f,n,\[Mu]0,\[Mu]]Log[\[Mu]0/\[Mu]]
])/;(n===1||n===2||n===3)


\[Gamma]evolutor[f_,n_,\[Mu]_,\[Mu]0_]:=(Which[
(\[Mu]0<=mCHARM && \[Mu]<=mCHARM)||(mCHARM<\[Mu]0<=mBOTTOM && mCHARM<\[Mu]<=mBOTTOM)||(mBOTTOM<\[Mu]0 && mBOTTOM<\[Mu]),Sum[\[Gamma]V[f,k,ActiveN[\[Mu]]]As[k,n,\[Mu],\[Mu]0],{k,1,n}],
\[Mu]0<=mCHARM && mCHARM<\[Mu]<=mBOTTOM,Sum[\[Gamma]V[f,k,4]As[k,n,\[Mu],mCHARM]+\[Gamma]V[f,k,3]As[k,n,mCHARM,\[Mu]0],{k,1,n}],
\[Mu]0<=mCHARM && mBOTTOM<\[Mu],Sum[\[Gamma]V[f,k,5]As[k,n,\[Mu],mBOTTOM]+\[Gamma]V[f,k,4]As[k,n,mBOTTOM,mCHARM]+\[Gamma]V[f,k,3]As[k,n,mCHARM,\[Mu]0],{k,1,n}],
mCHARM<\[Mu]0<=mBOTTOM && mBOTTOM<=\[Mu],Sum[\[Gamma]V[f,k,5]As[k,n,\[Mu],mBOTTOM]+\[Gamma]V[f,k,4]As[k,n,mBOTTOM,\[Mu]0],{k,1,n}],
\[Mu]<=mCHARM && mCHARM<\[Mu]0<=mBOTTOM,Sum[-\[Gamma]V[f,k,3]As[k,n,mCHARM,\[Mu]]-\[Gamma]V[f,k,4]As[k,n,\[Mu]0,mCHARM],{k,1,n}],\[Mu]<=mCHARM && mBOTTOM<\[Mu]0,Sum[-\[Gamma]V[f,k,3]As[k,n,mCHARM,\[Mu]]+\[Gamma]V[f,k,4]As[k,n,mBOTTOM,mCHARM]-\[Gamma]V[f,k,5]As[k,n,\[Mu]0,mBOTTOM],{k,1,n}],
mCHARM<\[Mu]<=mBOTTOM && mBOTTOM<=\[Mu]0,Sum[-\[Gamma]V[f,k,4]As[k,n,mBOTTOM,\[Mu]]-\[Gamma]V[f,k,5]As[k,n,\[Mu]0,mBOTTOM],{k,1,n}]])/;(n===1||n===2||n===3)


RTMD[f_,n_?IntegerQ,bt_?NumberQ,\[Zeta]f_?NumberQ,\[Mu]f_?NumberQ,\[Zeta]i_?NumberQ,\[Mu]i_?NumberQ,\[Mu]0_?NumberQ]:=Exp[\[CapitalGamma]Logevolutor[f,n,\[Mu]f,\[Mu]i]-\[Gamma]evolutor[f,n,\[Mu]f,\[Mu]i]+\[CapitalGamma]evolutor[f,n,\[Mu]f,\[Mu]i]Log[\[Mu]i^2/\[Zeta]f]-(\[CapitalGamma]evolutor[f,n,\[Mu]i,\[Mu]0]+Dboundary[f,n,\[Mu]0,bt])Log[\[Zeta]f/\[Zeta]i]]
RTMD[f_,n_?IntegerQ,bt_?NumberQ,\[Zeta]f_?NumberQ,\[Mu]f_?NumberQ,\[Mu]i_?NumberQ,\[Mu]0_?NumberQ]:=RTMD[f,n,bt,\[Zeta]f,\[Mu]f,C0^2/bt^2,\[Mu]i,\[Mu]0]
RTMD[f_,n_?IntegerQ,bt_?NumberQ,\[Zeta]f_?NumberQ,\[Mu]f_?NumberQ,\[Mu]0_?NumberQ]:=Exp[\[CapitalGamma]Logevolutor[f,n,\[Mu]f,\[Mu]0]-\[Gamma]evolutor[f,n,\[Mu]f,\[Mu]0]+\[CapitalGamma]evolutor[f,n,\[Mu]f,\[Mu]0]Log[\[Mu]0^2/\[Zeta]f]-Dboundary[f,n,\[Mu]0,bt]Log[\[Zeta]f bt^2/C0^2]]


RTMD2[f_,n_?IntegerQ,bt_?NumberQ,\[Zeta]f_?NumberQ,\[Mu]f_?NumberQ,\[Zeta]i_?NumberQ,\[Mu]i_?NumberQ,\[Mu]0_?NumberQ]:=Exp[\[CapitalGamma]Logevolutor[f,n,\[Mu]f,\[Mu]i]-\[Gamma]evolutor[f,n,\[Mu]f,\[Mu]i]+\[CapitalGamma]evolutor[f,n,\[Mu]f,\[Mu]i]Log[\[Mu]i^2/\[Zeta]i]-(\[CapitalGamma]evolutor[f,n,\[Mu]f,\[Mu]0]+Dboundary[f,n,\[Mu]0,bt])Log[\[Zeta]f/\[Zeta]i]]
RTMD2[f_,n_?IntegerQ,bt_?NumberQ,\[Zeta]f_?NumberQ,\[Mu]f_?NumberQ,\[Mu]i_?NumberQ,\[Mu]0_?NumberQ]:=RTMD2[f,n,bt,\[Zeta]f,\[Mu]f,C0^2/bt^2,\[Mu]i,\[Mu]0]
RTMD2[f_,n_?IntegerQ,bt_?NumberQ,\[Zeta]f_?NumberQ,\[Mu]f_?NumberQ,\[Mu]0_?NumberQ]:=Exp[\[CapitalGamma]Logevolutor[f,n,\[Mu]f,\[Mu]0]-\[Gamma]evolutor[f,n,\[Mu]f,\[Mu]0]+\[CapitalGamma]evolutor[f,n,\[Mu]f,\[Mu]0]Log[\[Mu]0^2/\[Zeta]f]-Dboundary[f,n,\[Mu]0,bt]Log[\[Zeta]f bt^2/C0^2]]


End[];
EndPackage[];
