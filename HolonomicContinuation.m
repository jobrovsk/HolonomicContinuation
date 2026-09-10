(* ::Package:: *)

(* ::Text:: *)
(*Copyright (C) 2026  Abilio de Freitas, Carsten Schneider*)
(**)
(*This file is part of HolonomicContinuation.*)
(**)
(*HolonomicContinuation is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation; either version 3 of the License, or  (at your option) any later version.  See https://www.gnu.org/licenses/*)


(* ::Section:: *)
(*Start package*)


(* ::Input::Initialization:: *)
$HolonomicContinuationVersion="HolonomicContinuation Package by Abilio De Freitas, Jakob Obrovsky and Carsten Schneider; RISC Linz \[LongDash] V 1.2 (09/07/2026)";
If[TrueQ[$Notebooks],CellPrint[Cell[BoxData[$HolonomicContinuationVersion],"Print",FontColor->RGBColor[0,0,0],CellFrame->0.5,Background->RGBColor[0.796887,0.789075,0.871107]]],
Print[$HolonomicContinuationVersion]];


Get[FileNameJoin[{DirectoryName[$InputFileName],"RecToValuesFLINT_V3.3.m"}]];


(* ::Input::Initialization:: *)
BeginPackage["HolonomicContinuation`"]
ClearAll@@Names["HolonomicContinuation`*"];


(* ::Input::Initialization:: *)
FindIndicialShift::usage="FindIndicialShift[diffeqz, g, z, \[Alpha]] finds the shift in the power of z^(-r+shift-1+k+\[Alpha]) when we replace the function g[z] in the differential equation 'diffeqz' required to find the corresponding indicial equation. Here 'r' is the order of the differential equation, 'k' is a generic integer power (the summation index in the power expansion in 'z'), and '\[Alpha]' is an arbitrary power in 'z' multiplying the expansion, which ends up being the variable in which the indicial equation is expressed."

CheckIndicial::usage="CheckIndicial[diffeqz, g, z, \[Alpha]] computes the indicial equation associated with the differential equation 'diffeqz', which is satisfied by the function g[z]. The indicial equation is given in terms of the variable \[Alpha]. The output will be the indicial equation together with an printout indicating whether the indicial equation has only integer and/or half-integer solutions."


(* ::Input::Initialization:: *)
MakeIntegerDE::usage="MakeIntegerDe[diffeq,g[z]] clears integer denominators, i.e., produces a differential equation for g[z] with integer coefficients.";


(* ::Input::Initialization:: *)
GetDEQInf::usage="GetDEQInf[diffeqs, {s,spt}, g,z] obtains a differential equation at the point s=+/-infinity, starting from the differential equation 'diffeqs', which is the one that we obtain from the expansion at s=0. The resulting differential equation is written in terms of g[z] (that is, it's obyed by g[z]). The relation between the variable 'z' and the variable 's' depends on 'spt' :

If spt=-Infinity then s = -1/z
If spt=Infinity then s = 1/z

The resulting differential equation will be satisfied by a power/log expansion in z at z=+/-Infinity."


(* ::Input::Initialization:: *)
GetDEQspt::usage="GetDEQspt[diffeqs, {s, spt}, g, z, qlist] obtains a differential equation at the point s=spt, starting from the differential equation 'diffeqs', which is the one that we obtain from the expansion at s=0. The resulting differential equation is written in terms of g[z] (that is, it's obyed by g[z]). The relation between the variable 'z' and the variable 's' depends on 'spt' and the list of points in 'qlist' as follows:

If spt < 0 and spt \[NotElement] qlist  then s = z + spt
If spt > 0 and spt \[NotElement] qlist  then s = spt - z
If spt < 0 and spt \[Element] qlist  then s = \!\(\*SuperscriptBox[\(z\), \(2\)]\) + spt
If spt > 0 and spt \[Element] qlist  then s = spt - \!\(\*SuperscriptBox[\(z\), \(2\)]\) 

The resulting differential equation will be satisfied by a power/log expansion in z at z=0.

Options:
\"Print s-point\" can be set to \"yes\", in which case the value of 'spt' will be printed, or \"no\", in which case it won't. Default: \"no\".
\"Check indicial equation\" can be set either to \"yes\", in which case a printout appears indicating what type of solutions the indicial equation has (this is done using the function 'CheckIndicial') or \"no\", in which case it won't.  Default: \"no\"."



(* ::Input::Initialization:: *)
ParallelGetDEQspt::usage="ParallelGetDEQspt[diffeqs, {s, spt}, g, z, qlist, nsplit] computes the same as 'GetDEQspt', but it first splits the differential equation into several pieces, evaluates each one of then separately in parallel, and the adds up the results. This allows to obtain the result much faster than with 'GetDEQspt'. The splitting is done based on the following options.

If \"Method\"\[Rule]\"One by one\", then 'diffeq' is split in as many terms as it has, in other words, if the order of the differential equation is r, then 'diffeq' will be split in a list of r+2 terms (one for each derivative in g[s], another one for the g[s] term itself, and another one for the inhomogeneous part).

If \"Method\"\[Rule]\"Reduced split\", then 'diffeq' is split using the function 'ReducedSplit'.

If \"Method\"\[Rule]\"Optimal split\" (default), then 'diffeq' is split using the function 'OptimalSplit'

It also has the options of 'GetDEQspt'."


(* ::Input:: *)
(**)


(* ::Input::Initialization:: *)
GetAnsatzParameters::usage="GetAnsatzParameters[diffeq, g, z, a, teststart, testnlogs, nc] finds the appropriate lowest power in the expansion variable of an ansatz (power-log expansion) solution to a differential equation, as well as the appropriate power of the log, and the free coefficients (using the differential equation, all other coefficients can be written in terms of these). 

The input variables are 
'diffeq': the differential equation.
'g': the function that satisfies the differential equation. 
'z': the variable on which 'g' depends.
'a': a symbol to construct the coefficients of the ansatz. 
'teststart': a preliminary test lowest power of 'z' in the ansatz (which should be low enough to make sure nothing is missed). 
'testnlogs': a preliminary test highest power of the log (it should be chosen high enough). 
'nc': the highest power in 'z' in the ansatz. 

The output is given in the form {start,nlogs,freecoeffs}, where 'start' is the correct lowest power of 'z', 'nlogs' is the correct highest power of the log, and 'freecoeffs' is the list of free coefficients.\[IndentingNewLine]\[IndentingNewLine]If the option  \"Regular Point \" is set to True, the output can be given without any extra calculations."



(* ::Input::Initialization:: *)
GetCoeffSubs::usage="GetCoeffSubs[diffeq, g, z, a, teststart, testnlogs, nc] finds the list of substitutions of all coefficients in a power-log expansion ansatz solution to a differential equation in terms of unconstrained coefficients (to be determined by matching conditions).

The input variables are 
'diffeq': the differential equation.
'g': the function that satisfies the differential equation. 
'z': the variable on which 'g' depends.
'a': a symbol to construct the coefficients of the ansatz. 
'teststart': a preliminary test lowest power of 'z' in the ansatz (which should be low enough to make sure nothing is missed). 
'testnlogs': a preliminary test highest power of the log (it should be chosen high enough). 
'nc': the highest power in 'z' in the ansatz.

If in Mathematica several kernels are launched, the underlying solver will be executed in parallel mode.

Options: \"Number of coeffs to get parameters\" sets the number of coefficients used to determine the ansatz parameters. Default: 50.";


(* ::Input::Initialization:: *)
DEToRE::usage="DEToRE[diffeq,g,h,z,n] produces the underlying recurrences of the differential equation diffeq in g[z] and the inhomogenous part h[z] in n. The recurrence is needed for the command GetCoeffSubsFast.";

GetCoeffSubsFast::usage="GetCoeffSubsFast[initialL,diffeq,receq,g,h,z,n,nc,nlogs,start,precision,noKernels] computes the list of substitutiuons of all coefficients in a power-log expansion ansatz solution to a differential equation in terms of unconstrained coefficients (to be determined by matching conditions).

The input variables are 
'initialL: the number of initial values in a to prolong the expansion coefficients by the underlying recurrence receq.
'diffeq': the differential equation; if there are non-trivial log-contributions, the right hand side must be a generic function h[z].
'receq': the underlying recurrence of the differential equation where g[n] and h[n+v] with v an integers stands for the coefficients of the power series representations of g[z] and h[z]; it can be produced with the command DEToRE.
'g': the function that satisfies the differential equation. 
'h': the rhs function of the differential equation (arising also on the rhs of receq). 
'z': the variable on which 'g' and 'g' depends.
'a': a symbol to construct the coefficients of the ansatz. 
'start': the lowest power of 'z' of the solution expansion. 
'nlogs': the highest power of the log. 
'nc': the highest power in 'z' of the desired solution expansion.
'noKernels'n: If the option UseFlintByC is set to False, the Mathematica parallelization can be utilized. If a certain number of kernels are launched, they will be used. Alternatively, one define with 'noKernels' a positive integer. Then the kernels are launched as optimally needed but does not load more than noKernels. If UseFlintByC is set to True, then up to 'noKernels' many threads are used by the C-backend.

Remarks: 
'receq' can be derived by the command REtoDE.
'initialL', 'start' and 'nlogs' can be derived by the command GetCoeffSubs.
";


(* ::Input::Initialization:: *)
NumberOfInitialValues::usage="NumberOfInitialValues[receq,g[n]] determines the number of initial values needed to prolong the sequence (under the assumption that the sequence is two-sided and the negative entries are zero.";


(* ::Input::Initialization:: *)
UseFlintByC::usage="If this option is set to True, the package RecToValuesFLINT is used in addition. So it has to be loaded into the system beforehand and the underlying code has to be compliled accordingly.";


BackendC::usage="Which of the availbalbe C-backend should be used. Currently available are \"rec_to_val_V1\" and \"rec_to_val_V2\". The first,\"rec_to_val_V1\", needs only very little memory. 
The second, \"rec_to_val_V2\", is somewhat faster (up to 2 times) but needs more memory. This option is ignored if UseFlintByC->False.";


(* ::Input::Initialization:: *)
MatchExpansions::usage="MatchExpansions[func1, func2, s, freecoeffs, coeff, {point, delta}, workingprecision, mwp] matches two expansions, one of which is completely known, and one that contains free coefficients to be determined numerically. The matching is done by evaluating the expansions at a set of points and equating the results at each point. The number of points should be equal to the number of free coefficients plus an extra coefficient for which we already know the value. This extra coefficient is introduced artificially to be able to estimate the precision of the remaining free coefficients. In this way, we obtain a system of equations for the free coefficients that we can solve numerically.

The input variables are
'func1': the known expansion.
'func2': the expansion containing free coefficients to be determined.
's': the expansion variable.
'freecoeffs': the list of free coefficients.
'coeff': the extra coefficient with known value (often zero).
'point': first point where the matching is done.
'delta': separation between consecutive points.
'workingprecision': self-explanatory.
'wmp': used to reduce the working precision in case it's necessary ('NSolve' is used and it may demand it).

The output is a list. The first item in the list is the size of the extra coefficient calculated with the matching conditions, which we can then compare with the known value. The second item in the list is the solution to all free coefficients found using the matching conditions.

Remark: The functions 'func1' and 'func2' must be built by the function call 'BuildMatchExpansions'."






(* ::Input::Initialization:: *)
GetBestPointMatch::usage="GetBestPointMatch[s,freecoeffs,testcoeff,delta,workingprecision1,wmp,{lpoint,rpoint},nstartpts,prec,iterations] executes systematically MatchExpansions trough the interval [lpoint,rpoint] by the disection method using a certain number of iterations specified by 'iteration'. Here 'nstartpts' determines the number of check points within the specified interval in which one searches a good point. The inputs 's','freecoeffs','testcoeff','delta','workingprecision1','wmp' are the same as described for MatchExpansions.

The output is a list. The first item is the best point represented as a rational number and the second entry is the expected precision using 'testcoeff'. The value 'prec' determines how good the rational number value of the found point approximtes the floating point representation. The parameters 's','freecoeffs','testcoeff','delta','workingprecision1','wmp' are the same as described for MatchExpansions.

Remark: To make this command feasible, it is executed in parallel. Thus suffiently many kernels should be launched. In particlar, the definitions of the values 'freecoeffs','testcoeff','funcA','funcB' (for details see MatchExpansions) must be distributed to the subkernels, the value '$MaxExtraPrecision' must be set high enough in the subkernels and the package 'HolonomicContinuation.m' must be loaded into the subkernels. This can be carried out, e.g., with the calls

    DistributeDefinitions[freecoeffs,testcoeff,funcA,funcB]
    ParallelEvaluate[$MaxExtraPrecision =10200]
    ParallelEvaluate[Get[''HolonomicContinuation.m'']]"


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)


(* ::Input::Initialization:: *)
GetIniPoints::usage="GetIniPoints[lpoint, rpoint, nstartpts] produces a list of 'nstartpts' equally spaced values between 'lpoint' and 'rpoint'. The option \"Drop\" allows to drop some terms from the list. The default is \"Drop\"\[Rule]0 (nothing is dropped)."


(* ::Input::Initialization:: *)
TestCoeffInfo::usage="TestCoeffInfo[coeff] determines the power in the expansion variable and the power of the log associated to the coefficient 'coeff'. For example,

TestCoeffInfo[a2[-3]]

will give

{a2[-3],-3,2}"



(* ::Input::Initialization:: *)
BuildMatchExpansions::usage="BuildMatchExpansions[FF, i, z, s, {smA, startA, ncA, nlogsA, ruleA}, {smB, startB, ncB, nlogsB, ruleB}, {coeff, j, k}] constructs two expansions to be matched. The output is a list containing the two expansions. The inputs 'FF' and 'i' are the form factor and the specific case (an integer) under consideration, respectively. The expansions are built in the variable 'z' using the function 'BuildExpansion'. The first expansion is built using the input parameters {smA,startA,ncA,nlogsA,ruleA}, where startA, ncA and nlogsA will be the input for the function 'BuildExpansion'. The variable 'smA' indicates the point in 's' around which we are doing the expansion, and 'ruleA' is a rule for replacing 'z' in terms of 's'. The second expansion is built using the input parameters {smB,startB,ncB,nlogsB,ruleB} in a similar way, but an extra term is added using 'coeff'. This should be a coefficient we know beforehand to be equal to zero, which we add to the expansion in order to be able to determine the precision of our calculation. The values of 'j' and 'k' are such that the extra term added is given by coeff*z^j*Log[z]^k.

In order to use this function, the list of coefficient substitutions must be available for both expansions, as well as the numerical solution of the free coefficients of the first expansion, since 'BuildMatchExpansions' looks for this information in order to build the expansions, based on the values of 'smA' and 'smB'."




(* ::Input::Initialization:: *)
Begin["`Private`"]


$startTimeUsed=TimeUsed[];


(* ::Section:: *)
(*Transform to DEs at different points*)


(* ::Input::Initialization:: *)
RemoveRational[expr_]:=Module[{list,newlist},
list=Table[expr[[i]],{i,1,Length[expr]}];
newlist=DeleteCases[DeleteCases[list,_Rational],_Integer];
newlist=newlist /. List->Times;
If[expr===0,0,newlist] ]


(* ::Input::Initialization:: *)
Options[FindIndicialShift]={"Initial shift"->7,"Details"->"no"};

FindIndicialShift[diffeqz_,g_,z_,\[Alpha]_,OptionsPattern[]]:=Module[{r,j,Z,k,inishift,shift,i,n,test,indicial1,indicial2,indicial3,indicial4},
inishift=OptionValue["Initial shift"];
r=Exponent[diffeqz /. Derivative[j_][g][z]:>Z^j,Z]; (* Order of the differential equation *)
shift=inishift;
test=True;
indicial1=RemoveRational[Factor[FunctionExpand[Coefficient[Expand[diffeqz /. g[z]->z^(\[Alpha]+k) /. Derivative[n_][g][z]->D[z^(\[Alpha]+k),{z,n}]],z^(-r+shift-1+k+\[Alpha]) ,1]/. k->0]]];
indicial2=RemoveRational[Factor[FunctionExpand[Coefficient[Expand[diffeqz /. g[z]->z^(\[Alpha]+k) /. Derivative[n_][g][z]->D[z^(\[Alpha]+k),{z,n}]],z^(-r+shift-2+k+\[Alpha]) ,1]/. k->0]]];
indicial3=RemoveRational[Factor[FunctionExpand[Coefficient[Expand[diffeqz /. g[z]->z^(\[Alpha]+k) /. Derivative[n_][g][z]->D[z^(\[Alpha]+k),{z,n}]],z^(-r+shift-3+k+\[Alpha]) ,1]/. k->0]]];
indicial4=RemoveRational[Factor[FunctionExpand[Coefficient[Expand[diffeqz /. g[z]->z^(\[Alpha]+k) /. Derivative[n_][g][z]->D[z^(\[Alpha]+k),{z,n}]],z^(-r+shift-4+k+\[Alpha]) ,1]/. k->0]]];
test=Not[{indicial1,indicial2,indicial3,indicial4}==={0,0,0,0}];
While[test,
shift=shift-1;    If[OptionValue["Details"]=="yes",Print[shift]]; 
indicial1=indicial2;
indicial2=indicial3;
indicial3=indicial4;
indicial4=RemoveRational[Factor[FunctionExpand[Coefficient[Expand[diffeqz /. g[z]->z^(\[Alpha]+k) /. Derivative[n_][g][z]->D[z^(\[Alpha]+k),{z,n}]],z^(-r+shift-4+k+\[Alpha]) ,1]/. k->0]]];
 If[OptionValue["Details"]=="yes",Print[{indicial1,indicial2,indicial3,indicial4}]];
test=Not[{indicial1,indicial2,indicial3,indicial4}==={0,0,0,0}]];
shift ]


(* ::Input::Initialization:: *)
Options[CheckIndicial]={"Initial shift"->7,"Details"->"no"};

CheckIndicial[diffeqz_,g_,z_,\[Alpha]_,OptionsPattern[]]:=Module[{r,shift,j,Z,indicial,k,n,sol,solseq,integerq,halfintegerq},
r=Exponent[diffeqz /. Derivative[j_][g][z]:>Z^j,Z]; (* order of the differentia equation *) 
shift=FindIndicialShift[diffeqz,g,z,\[Alpha],"Initial shift"->OptionValue["Initial shift"],"Details"->OptionValue["Details"]];
indicial=RemoveRational[Factor[FunctionExpand[Coefficient[Expand[diffeqz /. g[z]->z^(\[Alpha]+k) /. Derivative[n_][g][z]->D[z^(\[Alpha]+k),{z,n}]],z^(-r+shift+k+\[Alpha]) ,1]/. k->0]]];
sol=Solve[indicial==0,\[Alpha]];
solseq=Table[sol[[i,1,2]],{i,1,Length[sol]}]; (* list of values of the solutions to the indicial equation *)
integerq=AllTrue[solseq,IntegerQ]; (* finding out if all solutions consist of only integers *) 
halfintegerq=AllTrue[2*solseq,IntegerQ] && AnyTrue[2*solseq,OddQ];(* finding out if all solutions consist of only integers or half-integers *) 
Which[
sol!={{}} && integerq,Print["The indicial equation has only integer solutions."],
sol!={{}} && halfintegerq,Print["The indicial equation has integer and half-integer solutions (only)."],
True,Print["The indicial equation has solutions outside of the set of integers and half-integers."]];
indicial ]


GetDEQInf[deIn_, {s_,spt_}, g_,z_]:=
Module[{de=deIn,derivsubs,ord,sign,inhompart,exp,j},
ord=Max[Cases[{de},Derivative[A_][_][_]->A,Infinity],0];
derivsubs={Derivative[1][g][s]->z^2*Derivative[1][g][z]};
sign=If[spt===-Infinity,-1,1];
Do[
	derivsubs=Append[derivsubs,D[g[s],{s,j+1}]->-sign*Expand[z^2*D[derivsubs[[j,2]],z]]];
,{j,ord}];
de=Collect[de /. derivsubs/. s->sign/z /.{g[_]->g[z]},{Derivative[_][g][z],g[z]},Expand];
inhompart=de /. Derivative[_][g][z]->0 /. g[z]->0;
If[inhompart=!=0,
de=de-inhompart;
de=Collect[D[-inhompart,z]*de-(-inhompart)*D[de,z],{g[z],Derivative[_][g][z]},Expand];
];
exp=Exponent[Expand[de /. z->1/z],z];
de=Collect[z^exp*de,{g[z], Derivative[_][g][z]},Expand];
MakeIntegerDE[de,g[z]]
]


(* ::Input::Initialization:: *)
Options[GetDEQspt]={"Print s-point"->"no","Check indicial equation"->"no","Initial shift"->7,"Details"->"no"};

GetDEQspt[diffeqs_,{s_,spt_},g_,z_,qlist_,opts:OptionsPattern[]]:=Module[{diffeqz,inhompart,hompart,homdiffeqz,res},
diffeqz=GetDEQsptInternal[diffeqs,{s,spt},g,z,qlist,opts];
inhompart=diffeqz /. Derivative[_][g][z]->0 /. g[z]->0; (* inhomogeneous part of the differential equation *) 

hompart=diffeqz-inhompart; (* homogeneous part of the differential equation *) 


res=If[inhompart===0,diffeqz,Collect[D[-inhompart,z]*hompart-(-inhompart)*D[hompart,z],{g[z],Derivative[_][g][z]},Expand]]; (* homogenization of the differential equation (only if needed) *)
If[OptionValue["Print s-point"]=="yes",Print[spt],Null];If[OptionValue["Check indicial equation"]=="no",Null,CheckIndicial[res,g,z,\[Alpha],"Initial shift"->OptionValue["Initial shift"],"Details"->OptionValue["Details"]]];
MakeIntegerDE[res,g[z]]
 ]


(* ::Input::Initialization:: *)
GetDEQsptInternal[diffeqs_,{s_,spt_},g_,z_,qlist_,OptionsPattern[]]:=Module[ {nderiv,diffeqz,n,derivsubs,Z,newderiv,j},
nderiv=Exponent[diffeqs /. Derivative[n_][g][s]->Z^n,Z]+3;
diffeqz=Which[
Not[MemberQ[qlist,spt]] && spt<0,   Collect[diffeqs /. s->z+spt /. Derivative[n_][g][z+spt]->Derivative[n][g][z] /. g[z+spt]->g[z] ,{Derivative[_][g][z],g[z]},Expand],
Not[MemberQ[qlist,spt]] && spt>0,   Collect[diffeqs /. s->spt-z /. Derivative[n_][g][spt-z]->(-1)^n Derivative[n][g][z] /. g[spt-z]->g[z] ,{Derivative[_][g][z],g[z]},Expand],
MemberQ[qlist,spt] && spt<0,     
derivsubs={Derivative[1][g][s]->1/2/z*Derivative[1][g][z]};
For[j=1,j<=nderiv,j++,
newderiv=D[g[s],{s,j+1}]->Expand[1/2/z*D[derivsubs[[j,2]],z]];
derivsubs=Append[derivsubs,newderiv];
];
Collect[diffeqs /. derivsubs/. s->z^2+spt /. g[z^2+spt]->g[z] ,{Derivative[_][g][z],g[z]},Expand],
MemberQ[qlist,spt] && spt>0,     derivsubs={Derivative[1][g][s]->-1/2/z*Derivative[1][g][z]};
For[j=1,j<=nderiv,j++,
newderiv=D[g[s],{s,j+1}]->Expand[-1/2/z*D[derivsubs[[j,2]],z]];
derivsubs=Append[derivsubs,newderiv];
];
Collect[diffeqs /. derivsubs/. s->spt-z^2 /. g[spt-z^2]->g[z] ,{Derivative[_][g][z],g[z]},Expand]  
]; 
diffeqz
]


(* ::Input::Initialization:: *)
FindReducedSplitValue[nterms_,nsplit_]:=Module[ {n=nsplit,q=1,r=2},
While[r>q,
{q,r}=QuotientRemainder[nterms,n];
n=n-1];
n+1]


(* ::Input::Initialization:: *)
ReducedSplit[list_,nsplit_]:=Module[{nterms,newnsplit,q,r,Q1,Q2},
nterms=Length[list];
newnsplit=FindReducedSplitValue[nterms,nsplit];
{q,r}=QuotientRemainder[nterms,newnsplit];
Q1=Table[Sum[list[[j]],{j,1+q*i,q*(i+1)}],{i,0,newnsplit-1}];
Q2={Sum[list[[j]],{j,1+newnsplit*q,nterms}]};
  Join[Q1,Q2] ]


(* ::Input::Initialization:: *)
OptimalSplit[list_,nsplit_]:=Module[{nterms,q,r,Q1,Q2},
nterms=Length[list];
{q,r}=QuotientRemainder[nterms,nsplit];
Q1=Table[Sum[list[[j]],{j,1+(q+1)*i,(q+1)*(i+1)}],{i,0,r-1}];
Q2=If[q==0,{},Table[Sum[list[[j]],{j,1+(q+1)*r+q*i,(q+1)*r+q*(i+1)}],{i,0,nsplit-r-1}]];
 Join[Q1,Q2] ]


(* ::Input::Initialization:: *)
Options[SplitDiffEq]={"Method"->"Optimal split"};

SplitDiffEq[diffeq_,g_,s_,nsplit_,OptionsPattern[]]:=Module[{terms={},nterms,order,Z,newnsplit,q,r,n,i,j,k,T0,T1,T,Q,Q0,Q1,Q2},
order=Exponent[diffeq /. Derivative[n_][g][s]->Z^n,Z]; (* order of the differential equation *)
T0={diffeq /. g[s]->0 /.  Derivative[n_][g][s]->0};
T1={Coefficient[diffeq,g[s],1]*g[s]};
T=Table[Coefficient[diffeq,Derivative[k][g][s],1]*Derivative[k][g][s],{k,1,order}];
T=Join[T0,T1,T];
Which[
OptionValue["Method"]=="One by one",          Table[T[[j]],{j,1,Length[T]}],
OptionValue["Method"]=="Reduced split",   ReducedSplit[T,nsplit] ,
OptionValue["Method"]=="Optimal split",   OptimalSplit[T,nsplit] ]  ]


(* ::Input::Initialization:: *)
Options[ParallelGetDEQspt]={"Method"->"Optimal split","Print s-point"->"no","Check indicial equation"->"no","Initial shift"->7,"Details"->"no"};

ParallelGetDEQspt[diffeq_,{s_,spt_},g_,z_,qlist_,nsplit_,OptionsPattern[]]:=Module[ {diffeqlist,l,i,T,diffeqz,inhompart,hompart,\[Alpha],res},
diffeqlist=SplitDiffEq[diffeq,g,s,nsplit,"Method"->OptionValue["Method"]];
l=Length[diffeqlist]; (* Print[l]; *)
DistributeDefinitions[diffeqlist,GetDEQsptInternal];

T=ParallelTable[GetDEQsptInternal[diffeqlist[[i]],{s,spt},g,z,qlist],{i,1,l}];
diffeqz=Plus@@T;
inhompart=diffeqz /. \!\(\*SuperscriptBox[\(g\), 
TagBox[
RowBox[{"(", "_", ")"}],
Derivative],
MultilineFunction->None]\)[z]->0 /. g[z]->0;
hompart=diffeqz-inhompart;
res=If[inhompart===0,diffeqz,Collect[D[-inhompart,z]*hompart-(-inhompart)*D[hompart,z],{g[z],\!\(\*SuperscriptBox[\(g\), 
TagBox[
RowBox[{"(", "_", ")"}],
Derivative],
MultilineFunction->None]\)[z]},Expand]];
If[OptionValue["Print s-point"]=="yes",Print[spt],Null];
If[OptionValue["Check indicial equation"]=="no",Null,CheckIndicial[res,g,z,\[Alpha],"Initial shift"->OptionValue["Initial shift"],"Details"->OptionValue["Details"]]];
MakeIntegerDE[res,g[z]]]


(* ::Section::Closed:: *)
(*Produce coefficient relations (find truncated series solutions)*)


(* ::Subsection::Closed:: *)
(*Get ansatz parameters*)


(* ::Input::Initialization:: *)
Options[GetAnsatzParameters]={"Regular Point"->False};


(* ::Input::Initialization:: *)
GetAnsatzParameters[diffeq_,g_,z_,a_,teststart_,testnlogs_,nc_,opts___Rule]:=GetAnsatzParametersInternal[diffeq,g,z,a,teststart,testnlogs,nc,opts]//ToStringForm[#,a]&


(* ::Input::Initialization:: *)
GetAnsatzParametersInternal[diffeq_,g_,z_,a_,teststart_,testnlogs_,nc_,opts___Rule]:=Module[{ncoeffs,ansatz,i,j,ord,n,r,deq,ansatzindeq,fv,lowestpower,coeffs,M,R,LogTerm,freepos,freecoeffs,coeffsubz,nlogs=testnlogs,start=teststart},

ord=Max[Cases[diffeq,Derivative[n_][g][z]->n,Infinity],0];

If[("Regular Point"/.{opts}/.Options[GetAnsatzParameters])===True,
Return[{0,0,Table[a[0,i],{i,0,ord-1}]}];
];


ncoeffs=nc+100;
ansatz =Collect[ Sum[ a[i,j]*z^j*Log[z]^i,{j,teststart,ncoeffs},{i,0,testnlogs}],{Log[z],a[_,_]},Expand];

deq=Collect[diffeq,{g[_],Derivative[_][g][z]},Expand];
deq=deq/.(Power[z,r_]/;r>ord+ncoeffs-teststart)->0;


ansatzindeq =deq /. g[z] -> ansatz /. {\!\(\*SuperscriptBox[\(g\), 
TagBox[
RowBox[{"(", "k_", ")"}],
Derivative],
MultilineFunction->None]\)[z]:>Expand[D[ansatz,{z,k}]]};

If[Head[ansatzindeq]=!=Plus,
fv=Expand[ansatzindeq];
fv=fv/.(Power[z,r_]/;r>ncoeffs)->0,
fv=0;
Do[
fv=fv+(Expand[ansatzindeq[[i]]]/.(Power[z,r_]/;r>ncoeffs)->0),
{i,Length[ansatzindeq]}
];
];
fv=Collect[fv,Log[z]];

lowestpower=-Exponent[fv /. z->1/z,z];
coeffs=Flatten[Table[Table[a[i,j],{i,0,testnlogs}],{j,nc,teststart,-1}]];

M=Table[
LogTerm=Coefficient[fv,Log[z],l];
Table[
Coefficient[LogTerm,z,i],
{i,lowestpower,nc-teststart+lowestpower}],{l,testnlogs,0,-1}] // Flatten;

M=Select[M,FreeQ[#,a[_,AA_]/;AA>nc]&];

M=Table[Coefficient[M,coeffs[[k]],1],{k,1,Length[coeffs]}]//Transpose; 
R=RowReduce[M];

freepos=Complement[Range[Length[coeffs]]-1,GetFixedCoeffs[R,nc,testnlogs,teststart]];
freecoeffs=Sort[Table[coeffs[[Length[coeffs]-freepos[[j]]]],{j,1,Length[freepos]}]];

If[Length[freecoeffs]==ord,Null,
Print["Warning: Amount of free coefficients DOES NOT match the order of the differential equation"]];
coeffsubz=Expand[Solve[MapThread[Equal,{DeleteCases[R . coeffs,0],Table[0,Length[DeleteCases[R . coeffs,0]]]}],Complement[Variables[DeleteCases[R . coeffs,0]],freecoeffs]]][[1]];
ansatz=Collect[ansatz  /. coeffsubz,Log[z],Expand] /. (z^r_/;r>nc)->0 ;
For[i=testnlogs,i>=0,i--,
nlogs=If[Coefficient[ansatz,Log[z],i]===0,nlogs-1,nlogs]];
start=-Exponent[ansatz /. z->1/z,z];
{start,nlogs,freecoeffs}]


(* ::Subsection::Closed:: *)
(*DE to RE and analysis of RE*)


(* ::Input::Initialization:: *)
NumberOfInitialValues[rec_,g_[n_]]:=
Module[{ord,lc,root},
ord=Cases[rec[[1]],g[A_]->A-n,Infinity]//Max;
lc=Coefficient[rec[[1]],g[ord+n]];
Max[(FindIntegerRoots[lc,n]//Max)-ord+1,0]
]


(* ::Input::Initialization:: *)
DEToRE[de_Equal,F_,H_,x_,N_]:=
Module[{rec,ord},
rec=de[[1]]-de[[2]];

rec=rec/.F->FF[1]/.H->FF[2];

ord=Max[Cases[{de},Derivative[A_][_][_]->A,Infinity],0];
rec=DEtoREPrep[rec,FF,x,N];


rec=DEtoREInternal[rec,FF,x,N,ord];
rec==0/.{FF[1]->F,FF[2]->H}
]


(* ::Input::Initialization:: *)
DEtoREPrep[fIn_,F_,x_,N_,start_:5]:=
Module[{f},
f=fIn/.Derivative[d_][F[B__]][x]:>MySum[F[B][N]Product[N-i,{i,0,d-1}]MyPower[x,N]/x^d,{N,d+start,Infinity}];
f=f/.F[B__][x]->MySum[F[B][N]MyPower[x,N],{N,start,Infinity}];
f
]


(* ::Input::Initialization:: *)
DEtoREInternal[equIn_,INT_,x_,N_,ordIn_]:=
Module[{equ,den,A,B,C,D,kk,ord},

equ=LinearCollectP[equIn,MySum[__],SpecFactor[#,0.2]&];
equ=If[Head[equ]===Plus,Apply[List,equ],{equ}];
equ=equ/.(MySum[A_,B_]/;Head[A]===Times):>Select[A,!FreeQ[#,x]&&FreeQ[#,B[[1]]]&]MySum[Select[A,!(!FreeQ[#,x]&&FreeQ[#,B[[1]]])&],B];


den=Apply[PolynomialLCM,Denominator[equ]];

ord=0;

equ=equ*den;
equ=Table[If[Head[equ[[i]]]===Times,Map[TimeConstrained[Factor[#],10,#]&,equ[[i]]],TimeConstrained[equ[[i]],10,equ[[i]]]],{i,Length[equ]}];


equ=Map[
(#/.MySum[A_,{B_,C_,D_}]:>Sum[(A/.N->N+kk)//MyNormalize,{kk,-Exponent[#,x],ord}])&,
equ
];
equ=equ/.MyPower[x,N]->1;
equ=Map[Collect[#,INT[__][_],ExpandTogether[Coefficient[#,x,0]]&]&,equ];
equ=Collect[Apply[Plus,equ],INT[__][_],ExpandTogether[#]&];

equ
]


(* ::Input::Initialization:: *)
LaurentPolynomialQ[f_]:=
Cases[{f},(Power[A_,B_Integer]/;(B<0&&Head[A]===Plus))->A,Infinity]==={}


(* ::Input::Initialization:: *)
ExpandTogether[f_]:=If[LaurentPolynomialQ[f],Expand[f],Together[f]]


(* ::Input::Initialization:: *)
LinearCollectP[fIn_, patt_, Fu_]:=
Module[{varL,res},
varL=Cases[{fIn},If[Head[patt]===List,Apply[Alternatives,patt],patt],Infinity]//Union;
res=LinearCollect[fIn, varL, Fu];
res
]


(* ::Input::Initialization:: *)
LinearCollect[fIn_, varL_, Fu_] :=
  Module[{subst, f, const,varLL},
If[Head[fIn]===List||Head[fIn]===Equal,Map[LinearCollect[#,varL,Fu]&,fIn],
varLL=Variables[varL];subst = Table[varLL[[i]] -> 0, {i, Length[varLL]}];const = fIn /. subst;f = fIn - const;

(*Fu[const]+Table[(Fu[f /. Append[Delete[subst, i], subst[[i, 1]] -> 1]])*subst[[i, 1]], {i, Length[subst]}]*)

f=Prepend[Table[(f /. Append[Delete[subst, i], subst[[i, 1]] -> 1]), {i, Length[subst]}],const];
subst=Prepend[Table[subst[[i, 1]],{i,Length[subst]}],1];
If[Fu=!=Identity,
If[System`Parallel`$SubKernel===True||Kernels[]==={}||ByteCount[f]/Length[f]<100000,
f=Map[Fu,f];,
f=ParallelMap[Fu,f]
];
];
f=f . subst;
f
]
]


(* ::Input::Initialization:: *)
SpecFactor[f__,time_:10]:=
Module[{res},
res=Together[f];
TimeConstrained[Factor[Numerator[res]],time,Numerator[res]]/Factor[Denominator[res]]
];


(* ::Input::Initialization:: *)
MyNormalize[f_]:=
Module[{var,varSubst,AA,BB},

var=Cases[{f},Power[AA_,BB_]/;!IntegerQ[BB]->MP[AA,BB],Infinity]//Union;

var=Join[var,Union[Cases[{f},MyPower[__],Infinity]]];

varSubst=var;
varSubst=varSubst/.{MyPower[A_,B_]:>MyPower[A,Expand[B]]};

varSubst=varSubst/.{MyPower[A_,B_Plus]:>Apply[Times,Map[MyPower[A,#]&,Apply[List,B]]]};

varSubst=varSubst/.{MyPower[AA_,(-1)BB_]->1/MyPower[AA,BB]};
varSubst=varSubst/.MyPower[AA_,BB_Integer]->Power[AA,BB];

f/.MapThread[Rule,{var,varSubst}]
];


(* ::Input::Initialization:: *)
FindIntegerRoots[f_,d_]:=
Module[{degrees},
degrees=Map[#[[1]]&,f//FactorList];
degrees=Select[degrees,(Exponent[#,d]===1)&];
degrees=d/.Map[Solve[(#==0),d]&,degrees]/. 
        d-> {}//Flatten;
degrees=Select[degrees, IntegerQ ];
degrees
]


(* ::Subsection::Closed:: *)
(*Compute general solutions (slow) for initial values*)


(* ::Input::Initialization:: *)
Options[GetCoeffSubs]={"Number of coeffs to get parameters"->50};

GetCoeffSubs[diffeqIn_,g_,z_,a_,teststart_,testnlogs_,nc_,OptionsPattern[]]:=Module[{diffeq,c,start,nlogs,ncoeffs,ansatz,i,j,ord,n,r,deq,ansatzindeq,fv,lowestpower,coeffs,shift,M,subst,vec,posL,posR,R,LogTerm,freepos,freecoeffs,coeffsubz,step,kernelL,ordD,const,pos,posZ,posExtra,MConst,varD,posVar,posVar2,varLS,blockSize,varCommon},
$startTimeUsed=TimeUsed[];

If[True||Global`PrintStep===True,
Print[{"Step1: Analyze the system",TimeUsed[]-$startTimeUsed,MaxMemoryUsed[]}];
];

diffeq=MakeIntegerDE[diffeqIn,g[z]];

{start,nlogs,freecoeffs}=GetAnsatzParametersInternal[diffeq,g,z,a,teststart,testnlogs,OptionValue["Number of coeffs to get parameters"]];

ord=Max[Cases[diffeq,Derivative[n_][g][z]->n,Infinity],0]; (* Order of the differential equation *)
ncoeffs=nc+If[nlogs>0||start<0,ord+Max[-start,0],0]+5; 

If[True||Global`PrintStep===True,
Print[{"Step2: Set up the matrix",TimeUsed[]-$startTimeUsed,MaxMemoryUsed[]}];
];
kernelL=Kernels[];

ansatz=Sum[ a[i,j]*z^j*Log[z]^i,{j,start,ncoeffs},{i,0,nlogs}];
ansatz=ansatz/.a[A_,B_]:>a[A]^(B-start+2);


deq=Collect[diffeq,{g[_],Derivative[_][g][z]},Expand];
deq=deq/.(Power[z,r_]/;r>ncoeffs-start)->0;

ansatz=Collect[ansatz,{Log[z],z},Expand];

ansatz={ansatz};
Do[
AppendTo[ansatz,Collect[D[ansatz[[-1]],z],{Log[z],z},Expand]],
{i,ord}
];

ordD=Min[Cases[ansatz,(Power[z,AA_]/;AA<0)->AA,Infinity]//Union,0];

ansatz=ansatz//Expand;


If[ordD<0,
ansatz=Map[If[Head[#]===Plus,Apply[Plus,Apply[List,#]z^(-ordD)],# z^(-ordD)]&,ansatz];
];

const=Map[If[Head[#]===Plus,Apply[List,#],#]&,ansatz]//Flatten//Denominator;
const=Apply[LCM,const];
If[const=!=1,
ansatz=Map[If[Head[#]===Plus,Apply[Plus,Apply[List,#]const],# const]&,ansatz];
];
ansatz=Table[Coefficient[ansatz,Log[z],i],{i,0,nlogs}]//Expand;

ansatzindeq =Table[deq /. g[z] -> ansatz [[i,1]]/. {\!\(\*SuperscriptBox[\(g\), 
TagBox[
RowBox[{"(", "k_", ")"}],
Derivative],
MultilineFunction->None]\)[z]:>(ansatz[[i,k+1]]/.(Power[z,r_]/;r>ncoeffs-start+k+10)->0)},
{i,Length[ansatz]}];



ansatzindeq=Map[Apply[List,#]&,ansatzindeq];


If[Global`PrintStep===True,Print[{T0,MaxMemoryUsed[],TimeUsed[]-$startTimeUsed}];
];

If[Global`MatrixGenNoKernels===True||System`Parallel`$SubKernel===True||Length[kernelL]===0||nlogs===0,

fv=0;
Do[

varLSubst=Cases[ansatzindeq[[i]],a[_],Infinity]//Union;
varLSubst=Table[varLSubst[[k]]->0,{k,Length[varLSubst]}];


(*do not normalize*)
const=1;
ExpandV=0;
Do[
If[Global`PrintStep===True,
Print[{i,k,TimeUsed[]-$startTimeUsed}];
];
exp=Table[0,{ncoeffs-ordD+10+1},{ncoeffs-start+20+1}];
Do[
If[Global`PrintStep===True,
Print[{i,k,j,TimeUsed[]-$startTimeUsed}];
];
exp=exp+CoefficientList[CoefficientList[ansatzindeq[[i,j]]/.Delete[varLSubst,k],varLSubst[[k,1]],ncoeffs-ordD+10+1],z,ncoeffs-start+20+1],
{j,Length[ansatzindeq[[i]]]}
];
exp=MapThread[Times,{exp,varLSubst[[k,1]]^Range[0,Length[exp]-1]}];
exp=Transpose[exp]*z^Range[0,Length[exp[[1]]]-1];
exp=Apply[Plus,exp//Flatten];
ExpandV=ExpandV+exp,
{k,Length[varLSubst]}
];
If[ordD<0,
ExpandV=Apply[Plus,Apply[List,ExpandV]*z^(ordD)/const];
];
fv=fv+Log[z]^(i-1)ExpandV,
{i,Length[ansatzindeq]}
],

DistributeDefinitions[ncoeffs,ordD,start];
ParallelEvaluate@SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->8];

ParallelEvaluate[$HistoryLength=0];
(*WithCleanup*)

fv=ParallelMap[
WithCleanup[(
varLSubst=Cases[#,a[_],Infinity]//Union;
varLSubst=Table[varLSubst[[k]]->0,{k,Length[varLSubst]}];


(*do not normalize*)
const=1;
ExpandV=0;
Do[
exp=Table[0,{ncoeffs-ordD+10+1},{ncoeffs-start+20+1}];
Do[
If[False&&Global`PrintStep===True,Print[{k,j,TimeUsed[]-$startTimeUsed}];
];
exp=exp+CoefficientList[CoefficientList[#[[j]]/.Delete[varLSubst,k],varLSubst[[k,1]],ncoeffs-ordD+10+1],z,ncoeffs-start+20+1],
{j,Length[#]}
];
exp=MapThread[Times,{exp,varLSubst[[k,1]]^Range[0,Length[exp]-1]}];
exp=Transpose[exp]*z^Range[0,Length[exp[[1]]]-1];
exp=Apply[Plus,exp//Flatten];
ExpandV=ExpandV+exp,
{k,Length[varLSubst]}
];
ExpandV),Clear[ExpandV];Clear[partL];ClearSystemCache[];]&,
ansatzindeq
];

If[Global`PrintStep===True,Print[{T0a,MaxMemoryUsed[],TimeUsed[]-$startTimeUsed}];
];

ClearSystemCache[];

(*do not normalize...*)
const=1;
If[ordD<0||const=!=1,
fv=Map[If[Head[#]===Plus,Apply[Plus,Apply[List,#]z^(ordD)/const],# z^(ordD)/const]&,fv];
];

If[Global`PrintStep===True,Print[{T0b,MaxMemoryUsed[],TimeUsed[]-$startTimeUsed}];
];


fv=Sum[fv[[i+1]]Log[z]^i,{i,0,Length[ansatzindeq]-1}];

If[Global`PrintStep===True,Print[{T0c,MaxMemoryUsed[],TimeUsed[]-$startTimeUsed}];
];


];

ansatzindeq=0;

fv=Collect[fv,Log[z]];

fv=fv/.Power[a[A_],B_]:>a[A,B+start-2];

lowestpower=-Exponent[fv /. z->1/z,z];

coeffs=Flatten[Table[Table[a[i,j],{i,0,nlogs}],{j,nc,start,-1}]];
shift=Which[ (* Possible shift in the powers of z for the contruction of the matrix M *)
(start==-1 || start==0) && nlogs==0,-ord,
(start<0 && nlogs==0) || (start==0 && nlogs==1),-1,
True,0];

If[nlogs===0&&start>=0,
shift=ord,
shift=0
];

M=Table[
LogTerm=Coefficient[fv,Log[z],l];
LogTerm=Apply[Plus,Apply[List,LogTerm]z^(-lowestpower)];
LogTerm=CoefficientList[LogTerm,z,nc+1-start];
LogTerm,
{l,nlogs,0,-1}] // Flatten;
M=Select[M,Complement[Variables[#],coeffs]==={}&];

M=Table[
GetMatrixRow[M[[k]],a,nlogs,start,nc],
{k,Length[M]}
];

M=PrepareRows/@M;

M=ReduceMatrix[M,ord];
If[True||Global`PrintStep===True,
Print[{"Step3: solve the system",TimeUsed[]-$startTimeUsed,MaxMemoryUsed[]}];
];

coeffsubz=MyNullSpaceQ[M];

If[True||Global`PrintStep===True,
Print[{"Step4: extract the relations ",TimeUsed[]-$startTimeUsed,MaxMemoryUsed[]}];
];

coeffsubz=FromNSToSubst[coeffsubz,coeffs,posZ];

If[True||Global`PrintStep===True,
Print[{"Step5: prepare output ",TimeUsed[]-$startTimeUsed,MaxMemoryUsed[]}];
];


{{nlogs,start},coeffsubz }//ToStringForm[#,a]&
]



(* ::Input::Initialization:: *)
FromNSToSubst[solIn_,coeffs_,posZ_]:=
Module[{R,pos,zeroV,mySubs},
R=Map[Reverse,Sort[Map[Reverse,solIn]]];
Do[
pos=Max[Position[R[[k]],(A_/;A=!=0),1,Heads->False]];
R[[k]]=R[[k]]/R[[k,pos]],
{k,Length[R]}];


zeroV=Table[0,{Length[R[[1]]]}];

Do[
If[Length[R]<k||R[[k,k]]=!=1||Complement[Drop[R[[k]],k],{0}]=!={},
R=Insert[R,zeroV,k]],
{k,Length[R[[1]]]}];


R=Transpose[IdentityMatrix[Length[R[[1]]]]-R];
R=Select[R,Union[#]=!={0}&];

If[ListQ[posZ],
R=R//Transpose;
Do[
R=Insert[R,Table[0,{Length[R[[1]]]}],posZ[[k]]],
{k,Length[posZ]}
];
R=R//Transpose;
(*coeffsubz=Join[coeffsubz,Table[UnitVector[Length[coeffsubz[[1]]],posZ[[k]]],{k,Length[posZ]}]]//Sort;*)
];

mySubs=Table[
pos=Position[R[[k]],1]//Min;
coeffs[[pos]]->-(R[[k]] . coeffs-coeffs[[pos]])//Expand,
{k,Length[R]}];
mySubs
]


(* ::Input::Initialization:: *)
MyParallelMap[IsParallel_,fu_,f_]:=
If[IsParallel===True,
ParallelMap[fu,f,Method->"FinestGrained"],
Map[fu,f]
]


Options[GenerateRelationFromRecC]={BackendC->"rec_to_val_V2"};
GenerateRelationFromRecC[{rec_,F_[n_],h_:dummy},{initialSubstIn_,a_},no_,IsInhom_,DigitPrec_,maxKernels_Integer,OptionsPattern[]]:=
Module[{result,check,maxNo,initialSubst,nu,extraValue,testvals,randomsubst,vars},

nu=NumberOfInitialValues[rec,F[n]];
maxNo=Map[#[[1,1]]&,initialSubstIn]//Max;
extraValue=maxNo-nu;

If[extraValue<5,
Print["Warning: I do not have suffently many initial values; I need at least :",nu," plus 5 extra values for safety/checking reasons."]'
];
extraValue=Min[extraValue,50];
initialSubst=Select[initialSubstIn,#[[1,1]]<maxNo-extraValue&];

result=RecToValuesFLINT`RecToValuesFLINT[{rec,F[n]},{initialSubst,a},no,DigitPrec,
	"Details"->True,
	"NumberOfThreads"->maxKernels,
	"Inhomogeneous"->IsInhom,
	"WriteOutputToFile"->False,
"Backend"->OptionValue[BackendC]
];
testvals=Table[a[k],{k,maxNo-extraValue,maxNo}];
vars=Union[Variables[Values[initialSubst]],Variables[Values[result]]];
randomsubst=Thread[vars->RandomInteger[{-100,100},Length[vars]]];
check=(testvals/.result/.randomsubst)-(testvals/.initialSubstIn/.randomsubst);
check=Max[Abs[check]];
If[!TrueQ[check<=1000 10^(-DigitPrec)],
	Print[check];
	Print[{(testvals/.result),(testvals/.initialSubstIn)}];
	Print["Error: the result does not agree with the initial values..."];
	Abort[];
];
Join[initialSubst,result]
]



(* ::Input::Initialization:: *)
GenerateRelationFromRec[{recIn_,F_[n_],h_:dummy},{initialSubst_,a_},no_,offsetNo_,DigitPrecIn_:Infinity]:=
Module[{IsParallel,SigmaFile,rec,ord,initialL,varL,DigitPrec,aNo,aNoMin,aNoMax,seqRes,k,initialP,seqP,testNo,test,check,zeroRecognition},

Clear[FP];
Clear[hP];


rec=recIn/.{F->FP,h->hP};


ord=Cases[rec[[1]],FP[A_]->A-n,Infinity];
ord=Max[ord]-Min[ord];
aNo=Cases[initialSubst,a[A_]->A,Infinity];
aNoMin=aNo//Min;
aNoMax=aNo//Max;
initialL=Table[a[k],{k,aNoMin,aNoMax}]/.initialSubst;
If[(IntegerQ[DigitPrecIn]&&DigitPrecIn>0),
DigitPrec=DigitPrecIn,
DigitPrec=Infinity
];
varL=Variables[initialL];
initialL=Map[Coefficient[initialL,#]&,varL];
initialL=MapThread[List,{Table[k,{k,1,Length[varL]}],initialL}];


If[Kernels[]=!={},
Print["I enter the parallel mode with ",Length[Kernels[]]," kernels."];
DistributeDefinitions[{JRecToListHorner,JEvalHorne}];
DistributeDefinitions[{rec,ord,aNoMin,DigitPrec,varL,h}];
IsParallel=True;
];


seqRes=MyParallelMap[IsParallel,
({k,initialP}=#;
Print["Consider: ",{k,varL[[k]]}];
Clear[FP];
Clear[hP];
hP[nn_Integer]:=Coefficient[h[nn],varL[[k]]];

seqP=JRecToListHorner[rec,initialP,aNoMin,no,FP[n]];
Clear[hP];
seqP=Join[initialP,seqP];
testNo=Min[no-aNoMin+1,Length[initialP]];
test=Take[seqP,testNo]-Take[initialP,testNo];
check=Union[test];
If[check=!={0},
Print["Error: the result does not agree with the initial values..."];
Abort[];
];
seqP)&,
initialL
];
seqRes=Transpose[seqRes] . varL;
seqRes=Table[a[k]->seqRes[[k-aNoMin+1]],{k,aNoMin,no}];
seqRes
]


(* ::Input::Initialization:: *)
GetFixedCoeffs[M_,ncoeffs_,nlogs_,start_]:=Module[{res={},i,j,k,n},
For[i=1,i<=Length[M],i++,
j=1;
n=0;
While[M[[i,j]]==0 && j<Length[M[[i]]],n++;j++];
res=Append[res,(nlogs+1)*(ncoeffs-start+1)-n-1]
];
res]


(* ::Input::Initialization:: *)
MakeIntegerDE[diffeqIn_,g_[z_]]:=
Module[{diffeq,c},
diffeq=Apply[List,diffeqIn];
c=diffeq/.{g[_]->1,Derivative[_][g][z]->1};
diffeq=diffeq/c;
c=Map[Apply[List,#]&,c];
c=c*Apply[PolynomialLCM,Denominator[Flatten[c]]];
c=Map[Apply[Plus,#]&,c];
diffeq=c . diffeq;
diffeq
];


(* ::Input::Initialization:: *)
GetMatrixRow[h_,a_,LogNo_,Start_,End_]:=
Module[{res,M,q,r},
If[h===0,Return[Table[0,{(End-Start+1)(LogNo+1)}]]];
res=h/.a[B_,A_]->a^((A-Start)*(LogNo+1)+(LogNo-B));
res=CoefficientList[res,a,(End-Start+1)(LogNo+1)];

res//Reverse
]


(* ::Input::Initialization:: *)
GetAVec[a_,A_,nlog_,start_]:=
Module[{qr},
qr=QuotientRemainder[A-2,nlog+1];
a[qr[[2]],qr[[1]]+start]
]



(* ::Input::Initialization:: *)
ReduceMatrix[MIn_,ord_]:=
Module[{M=MIn,MP,p,check,sol},
check=False;
p=NextPrime[Developer`$MaxMachineInteger/2];
MP=PolynomialMod[M,p];
While[check===False,
sol=NullSpace[MP,Modulus->p];
If[Length[sol]<ord,
M=Drop[M,-1];
MP=Drop[MP,-1],
check=True
];
];
M
];


(* ::Input::Initialization:: *)
ToStringForm[f_,a_]:=
Module[{var,varS},
var=Cases[f,a[_,_],Infinity]//Union;
varS=var/.a[A_,B_]:>ToExpression[(ToString[a]~~ToString[A])][B];
f/.MapThread[Rule,{var,varS}]
];


(* ::Input::Initialization:: *)
Clear[MyNullSpaceQ];


(* ::Input::Initialization:: *)
MyNullSpaceQ[MIn_]:=If[Global`AlwaysP===True||System`Parallel`$SubKernel===False&&Length[Kernels[]]>0,MyNullSpaceQParallel[MIn],MyNullSpaceQStandard[MIn]]


(* ::Input::Initialization:: *)
Clear[PrepareNSParallel];
PrepareNSParallel[M_,p_,no_,solNo_]:=
Module[{ML,pL,solL,solNoNew,pprod,sol},
pL=Table[NextPrime[p,k],{k,1,no}];

solL=Table[
Mod[M,pL[[k]]],
{k,Length[pL]}
];
solL=MapThread[List,{solL,pL}];
If[System`Parallel`$SubKernel===False&&Length[Kernels[]]>0,
solL=ParallelMap[{NullSpace[#[[1]],Modulus->#[[2]]],#[[2]],#[[2]]}&,solL],
solL=Map[{NullSpace[#[[1]],Modulus->#[[2]]],#[[2]],#[[2]]}&,solL];
];


If[solNo===Infinity,
Return[solL]
];
solNoNew=Map[Length[#[[1]]]&,solL]//Min;
solL=Select[solL,Length[#[[1]]]===solNoNew&];
If[solNoNew<solNo,
Print[{solNoNew,solNo}];
Print["All ealier solutions are wrong..."];
];
If[solL==={},
Print["Something is very weird: All primePs are wrong..."];
];
pL=Map[#[[2]]&,solL];
solL=Map[#[[1]]&,solL];

sol=Map[ChineseRemainder[##,pL]&,Transpose[solL,{3,1,2}],{2}];
{{sol,Apply[Times,pL],pL[[-1]]}}
];



(* ::Input::Initialization:: *)
Clear[MyNullSpaceQStandard];



(* ::Input::Initialization:: *)
MyNullSpaceQParallel[MIn_]:=
Module[{M,ML,solold,p,pLast,t=True,Mp,sol,solP,r,pprod,i,k,j,stepwise},
Print["Parallel version"];
M=PrepareRows/@MIn;
solold={};
stepwise=Max[Length[Kernels[]],1];
Which[
stepwise===1,stepwise=10,
stepwise===2,stepwise=12,
stepwise===3,stepwise=12,
stepwise===4,stepwise=12,
stepwise===5,stepwise=10,
stepwise===6,stepwise=12,
stepwise===7,stepwise=14,
stepwise===8,stepwise=16,
stepwise===9,stepwise=18
];
ML=PrepareNSParallel[M,Developer`$MaxMachineInteger^19,stepwise,Infinity];
k=1;
j=1;
{solP,p,pLast}=ML[[j]];
While[t,
If[solold=!={},sol={};
For[i=1,i<=Length[solold],i++,
r=ReconstructRationalNumber[#,pprod]&/@solold[[i]];
Mp=PolynomialMod[M,p];
If[Union[Mod[Mp . PolynomialMod[r,p],p]]==={0},sol=Append[sol,r];,sol={};
i=Length[solold]+1;
];
];
];
If[Length[sol]>0,t=False;,sol=solP;
If[sol==={},t=False;,If[Length[sol]<Length[solold]||solold==={},solold=sol;
pprod=p;,
If[Length[sol]>Length[solold],Print["Unlucky prime!"];,(*solold=Apply[ChineseRemainder[{##},{pprod,p}]&,Transpose[{solold,sol},{3,1,2}],{2}];*)
solold=Map[ChineseRemainder[##,{pprod,p}]&,Transpose[{solold,sol},{3,1,2}],{2}];
pprod=p*pprod;
];
];
];
];
Print["Number of primes: ",k," (total used time: ",TimeUsed[]-$startTimeUsed,", max used memory: ",MaxMemoryUsed[],", `prime` digits: ",Apply[Plus,DigitCount[p]],", CRA digits: ",Apply[Plus,DigitCount[pprod]],")"];
If[j<Length[ML],
j=j+1;
k=k+1,
Print["Compute the next homorphic images"];
ML=PrepareNSParallel[M,pLast,stepwise,Length[sol]];
k=k+stepwise;
j=1;
];
{solP,p,pLast}=ML[[j]];
];
sol
];




(* ::Input::Initialization:: *)
Clear[MyNullSpaceQStandard];


(* ::Input::Initialization:: *)
MyNullSpaceQStandard[MIn_]:=Module[{M,solold={},p,t=True,Mp,sol,r,pprod,i,k=0},
M=PrepareRows/@MIn;
p=NextPrime[Developer`$MaxMachineInteger^19];
While[t,
Mp=Mod[M,p];
If[solold=!={},sol={};
For[i=1,i<=Length[solold],i++,
r=ReconstructRationalNumber[#,pprod]&/@solold[[i]];
If[Union[Mod[Mp . PolynomialMod[r,p],p]]==={0},sol=Append[sol,r];,sol={};
i=Length[solold]+1;];];];
If[Length[sol]>0,t=False;,sol=NullSpace[Mp,Modulus->p];
If[sol==={},t=False;,If[Length[sol]<Length[solold]||solold==={},solold=sol;
pprod=p;,If[Length[sol]>Length[solold],Print["Unlucky prime!"];,(*solold=Apply[ChineseRemainder[{##},{pprod,p}]&,Transpose[{solold,sol},{3,1,2}],{2}];*)solold=Map[ChineseRemainder[##,{pprod,p}]&,Transpose[{solold,sol},{3,1,2}],{2}];
pprod=p*pprod;];];];];
Print["Number of primes: ",k++," (total used time: ",TimeUsed[]-$startTimeUsed,", prime digits: ",Apply[Plus,DigitCount[p]],", CRA digits: ",Apply[Plus,DigitCount[pprod]],")"];
p=NextPrime[p,-1];];
sol];



(* ::Input::Initialization:: *)
Clear[PrepareRows];


(* ::Input::Initialization:: *)
PrepareRows[r_]:=Module[{s},s=r*(LCM@@Denominator[r]);
s/Max[GCD@@s,1]]


(* ::Input::Initialization:: *)
Clear[ReconstructRationalNumber]


(* ::Input::Initialization:: *)
ReconstructRationalNumber[n_,p_]:=If[n===0,0,(((#[[2,2]]/#[[1,2,2]])&)[Internal`HGCD[p,Mod[n,p]]]*2)/2];


(* ::Subsection:: *)
(*Compute general solutions (fast) using the underlying recurrence and initial values*)


(* ::Input::Initialization:: *)
Options[GetCoeffSubsFast]={UseFlintByC->False,BackendC->"rec_to_val_V2"};


(* ::Input::Initialization:: *)
GetCoeffSubsFast[initialSIn_,de_,inputrec_,g_,h_,z_,a_,n_,noIn_,nlogs_,start_,precision_,maxKernels_Integer:0,opts:OptionsPattern[]]:=
Module[{initialS,no,ord,initialSN,varKnown,initialStep,recStep,resL,prec,res,aVar,UseCCode},
no=noIn;

initialS=Complement[Map[#[[2]]&,initialSIn]//Variables,Map[#[[1]]&,initialSIn]//Variables];
initialS=Join[MapThread[Rule,{initialS,initialS}],initialSIn];


UseCCode=OptionValue[UseFlintByC];
If[UseCCode===True&&Head[(?"RecToValuesFLINT`*")]===Missing,
Print["Warning: the package RecToValuesFLINT is not loaded; I use the standard Mathematica code..."];
UseCCode=False
];

If[Max[Map[#[[1,1]]&,initialS]]-30>no,
no=Max[Map[#[[1,1]]&,initialS]];
Print["Warning: the number of initial values are larger than the number of values to be calculated; this might cause problems. Thus I increase the number to ",no," to guarantee correctness."
];
];

If[Kernels[]=!={},
ord=Cases[de,Derivative[A_][_][_]->A,Infinity]//Max;
deL=Map[CoefficientList[#,z,no+250]&,Table[-Coefficient[de,Derivative[k][g][z]],{k,0,ord}]];
DistributeDefinitions[deL];
DistributeDefinitions[ExpandSeriesParaIn,ExpandDEParaIn];
];

resL={};
{initialStep,recStep}={initialS,inputrec};
prec=Infinity;
Do[
If[k===0,prec=precision];
If[k===nlogs,
res=TopLogCoeffSubs[initialStep,recStep,h,g,a,n,no,k,maxKernels,prec,UseFlintByC->UseCCode,BackendC->OptionValue[BackendC]],
res=LogCoeffSubs[initialStep,recStep,initialS,inputrec,de,h,g,z,a,n,no,k,start,maxKernels,prec,UseFlintByC->UseCCode,BackendC->OptionValue[BackendC]]
];
resL=Join[resL,res[[1]]];
{initialStep,recStep}={res[[2]],res[[3]]},
{k,nlogs,0,-1}
];

(*whyneeded???*)
initialSN=N[initialS,precision];
resL=N[resL,precision];

aVar=Table[ToExpression[ToString[a]~~ToString[k]],{k,0,nlogs}];
aVar=Table[aVar[[k]][A_]->aVar[[k]][Hold[Round[A]]],{k,1,nlogs+1}]/.Rule->RuleDelayed/.Hold->Identity;
{initialSN,resL}={initialSN,resL} /.aVar;

varKnown=Map[#[[1]]&,initialS];
resL=Select[resL,!MemberQ[varKnown,#[[1]]]&];
resL=Join[initialSN,resL];
resL
]


(* ::Code:: *)
(**)


Clear[JEvalHorner]
JEvalHorner[plist_?VectorQ,val_?IntegerQ]:=Fold[(val #1+#2)&,0,plist]

Clear[JRecToListHorner] 
JRecToListHorner[recIn_Equal,initial_?VectorQ,start_?IntegerQ,end_?IntegerQ,g_[n_Symbol]]:=
Module[{hh,rec,k,i,AAA,curNum,den,ord,values={},recK,oldDen,mul,vecGCD,lIndex,tIndex},
	{tIndex,lIndex}=MinMax[Cases[recIn,g[AAA_]->AAA-n,Infinity]];
	ord=lIndex-tIndex;
	rec=Table[Coefficient[recIn[[1]]-recIn[[2]],g[i+n]],{i,tIndex,lIndex}];
	rec=Reverse[CoefficientList[#,n]]&/@rec;
	hh=(recIn[[1]]-recIn[[2]]/.g[_]->0);
	den=LCM@@Denominator/@initial[[;;UpTo[ord]]];
	curNum=den initial[[;;UpTo[ord]]];
	Do[JPrint[k," : "];JPrint[" combinded: ",Timing[	
		recK=JEvalHorner[#,k-lIndex]&/@rec[[-Length[curNum]-1;;]];
		AppendTo[values,recK[[;;-2]] . curNum];
		values[[-1]]/=den;
		values[[-1]]+=(hh/.n->(k-lIndex));
		values[[-1]]/=-recK[[-1]];
		oldDen=den;
		den=LCM[oldDen,Denominator[values[[-1]]]];
		mul=den/Denominator[values[[-1]]];
		AppendTo[curNum,mul Numerator[values[[-1]]]];
		If[Length[curNum]>ord,curNum=Rest[curNum]];
		mul=den/oldDen;
		curNum[[;;-2]]*=mul;
		If[Mod[k,50]==0,
			vecGCD=GCD[GCD@@curNum,den];
			 (*vecGCD is here almost always 1*)
			curNum/=vecGCD;
			den/=vecGCD;
		]	
	][[1]]];,{k,start+Min[ord,Length[initial]],end}];
	
	Return[values];
]



(* ::Input::Initialization:: *)
Options[TopLogCoeffSubs]={UseFlintByC->False,BackendC->"rec_to_val_V2"}
TopLogCoeffSubs[initialS_,inputrec_,h_,g_,a_,n_,no_,LogDeg_,maxKernels_,precision_,OptionsPattern[]]:=Module[{substLogPart,inhomExpr,inhom,ordInhom,A,initialSubst,rec,h2,i,numberKernels,optimalKernels,newKernels,time,substLogN},

inhomExpr=0; 
inhom=Table[0,{no+1}]; 
ordInhom=0;
rec=inputrec /. h->h2; 
A=ToExpression[ToString[a]~~ToString[LogDeg]] ;
initialSubst=Select[initialS,!FreeQ[#[[1]],A]&]; 

h2[i_Integer]:=inhom[[i+1]] ;


If[OptionValue[UseFlintByC]===True,
{time,substLogPart}=GenerateRelationFromRecC[{rec,g[n],h2},{initialSubst,A},no,False,If[LogDeg===0,precision,Infinity],maxKernels,BackendC->OptionValue[BackendC]]//AbsoluteTiming,

numberKernels=Length[Kernels[]];

If[numberKernels<maxKernels,
optimalKernels=Variables[Map[#[[2]]&,initialSubst]]//Length;
If[optimalKernels>1&&optimalKernels>numberKernels,
newKernels=Min[optimalKernels-numberKernels,maxKernels-numberKernels];
Print["Launch ",newKernels," kernels"];
LaunchKernels[newKernels];
]
];
{time,substLogPart}=GenerateRelationFromRec[{rec,g[n],h2},{initialSubst,A},no,1,If[LogDeg===0,precision,Infinity]]//AbsoluteTiming;
];
Print["Time needed to get top contribution: ",time];

 substLogN=substLogPart; 
Do[
substLogN[[k,2]]=N[substLogN[[k,2]],precision]/.A[b_]:>A[Round[b]],
{k,Length[substLogN]}];

(* substRelationL=Join[substRelationL,substLogN]; *) 
{substLogN,substLogPart,inhomExpr} ]


(* ::Input::Initialization:: *)
ExpandSeriesIn[fIn_,gIn_,x_,len_]:=Module[{denF,gI,fI,denG,result},If[fIn===0||gIn===0,Return[0]];
fI=CoefficientList[fIn,x,len+1];
gI=CoefficientList[gIn,x,len+1];
denF=LCM@@Denominator/@fI;denG=LCM@@Denominator/@gI;
gI=denG*gI;fI=denF*fI;
result=ListConvolve[fI,gI,{1,1},0];
result/=denF denG;
result=result . Table[x^i,{i,0,Length[result]-1}];
result
]


(* ::Input::Initialization:: *)
ExpandSeries[f_,gIn_,x_,len_]:=
Module[{varL,g,res},
varL=Complement[Variables[gIn],{x}];
g=CoefficientRules[gIn,varL];
varL=Map[Apply[Times,MapThread[Power,{varL,#[[1]]}]]&,g];
res=Map[Apply[List,ExpandSeriesIn[f,#[[2]],x,len]]&,g];
res=Apply[Plus,res*varL//Flatten];
res
]


(* ::Input::Initialization:: *)
ExpandDE[deIn_,x_,no_]:=
Module[{de,f,g,part,res},
de=If[Head[deIn]===Plus,Apply[List,deIn],{deIn}];
res=Table[
part=de[[i]];
part=If[Head[part]===Times,Apply[List,part],{part}];
f=Apply[Times,Select[part,Complement[Variables[#],{x}]==={}&]];
g=Apply[Times,Select[part,Complement[Variables[#],{x}]=!={}&]];

ExpandSeries[f,g,x,no],
{i,Length[de]}
];
res=Apply[Plus,res];
res
]



(* ::Input::Initialization:: *)
MyCoefficientRules[f_List,varL_]:=
Module[{expL,res,coeffL},
coeffL=CoefficientRules[f,varL];
expL=Apply[Union,Map[#[[1]]&,coeffL,{2}]];
coeffL=Table[
res=Table[Select[coeffL[[i]],#[[1]]===expL[[k]]&],{i,Length[coeffL]}];
res=Map[If[#==={},0,#[[1,2]]]&,res];
expL[[k]]->res,
{k,Length[expL]}
];
coeffL
]



(* ::Input::Initialization:: *)
ExpandDEPara[DiffLIn_,x_,no_]:=
Module[{DiffL,varL,res},
varL=Complement[Variables[DiffLIn],{x}];
DiffL=MyCoefficientRules[DiffLIn,varL];
varL=Map[Apply[Times,MapThread[Power,{varL,#[[1]]}]]&,DiffL];
res=ParallelMap[Apply[List,ExpandDEParaIn[#[[2]],x,no]]&,DiffL];
res=Apply[Plus,res*varL//Flatten];
res
]



(* ::Input::Initialization:: *)
ExpandDEParaIn[LogContrD_,x_,no_]:=
Module[{inhomPart,res,part,f,g},
res=Table[
ExpandSeriesParaIn[deL[[i]],LogContrD[[i]],x,no],
{i,Length[deL]}
];
res=Apply[Plus,res];
res
]


(* ::Input::Initialization:: *)
ExpandSeriesParaIn[fIn_,gIn_,x_,len_]:=
Module[{denF,gI,fI,denG,result},If[fIn===0||gIn===0,Return[0]];
fI=fIn;
gI=CoefficientList[gIn,x,len+1];
denF=LCM@@Denominator/@fI;
denG=LCM@@Denominator/@gI;
gI=denG*gI;fI=denF*fI;
result=ListConvolve[fI,gI,{1,1},0];
result/=denF denG;
result=result . Table[x^i,{i,0,Length[result]-1}];
result
]


Options[LogCoeffSubs]={UseFlintByC->False,BackendC->"rec_to_val_V2"};
LogCoeffSubs[substLogPartA_,inhomExprA_,initialS_,inputrec_,de_,h_,g_,z_,a_,n_,no_,LogDeg_,start_,maxKernels_,precision_,OptionsPattern[]]:=Module[{Aabove,A,orderSubst,b,ordDE,c,LogContr,LogContrD,time,ord,inhomPart,inhomExpr,inhom,h2,rec,i,ordInhom,initialSubst,numberKernels,optimalKernels,newKernels,substLogPart,substLogN},
Aabove=ToExpression["a"~~ToString[LogDeg+1]];
A=ToExpression[ToString[a]~~ToString[LogDeg]];

orderSubst=start;
ordDE=Cases[de,Derivative[b_][g][z]->b,Infinity]; 
ordDE=Max[ordDE];

rec=inputrec /. h->h2; 

c=Apply[LCM,Denominator[Flatten[Map[If[Head[#]===Plus,Apply[List,#],#]&,Table[Aabove[n],{n,orderSubst,no}]/.substLogPartA]]]];
LogContr=Log[z]^(LogDeg+1)Collect[Sum[Expand[c Aabove[n]]z^(n),{n,orderSubst,no}]/.substLogPartA,z];
LogContrD={0};
Do[
{time,LogContr}=Collect[D[LogContr,z],{Log[z],z},Expand]//AbsoluteTiming;
AppendTo[LogContrD,Expand[(LogContr/.Power[Log[z],LogDeg+1]->0)]],
{i,1,ordDE}];



ord=Max[-Cases[LogContrD,(Power[z,b_]/;b<0)->b,Infinity]];
LogContrD=Map[Apply[Plus,Apply[List,#]z^ord]&,LogContrD];


If[Kernels[]=!={},
inhomPart=ExpandDEPara[LogContrD,z,ord+no+200],
inhomPart=-de/.{Derivative[A_][g][z]:>LogContrD[[A+1]],g[z]->LogContrD[[1]]};
inhomPart=ExpandDE[inhomPart,z,ord+no+200];
];

inhomPart=Apply[Plus,Apply[List,inhomPart]/z^ord/c];
inhomPart=inhomPart/.(Power[z,b_]/;(b>no+200))->0;

inhomExpr=inhomExprA+inhomPart; 
inhom=Coefficient[inhomExpr,Log[z],LogDeg]; 
ordInhom=Max[-Cases[inhom,(Power[z,b_]/;b<0)->b,Infinity],0] ;

inhom=CoefficientList[Apply[Plus,Apply[List,inhom]z^ordInhom],z]; 

h2[i_Integer]:=inhom[[i+1+ordInhom]];
initialSubst=Select[initialS,!FreeQ[#[[1]],A]&]; 
Print["ByteCount in initialSubst: ",initialSubst//ByteCount];

If[OptionValue[UseFlintByC]===True,
	{time,substLogPart}=GenerateRelationFromRecC[{rec,g[n],h2},{initialSubst,A},no,True,If[LogDeg===0,precision,Infinity],maxKernels,BackendC->OptionValue[BackendC]]//AbsoluteTiming
,
	numberKernels=Length[Kernels[]];
	If[numberKernels<maxKernels,
	optimalKernels=Variables[Map[#[[2]]&,initialSubst]]//Length;
	If[optimalKernels>1&&optimalKernels>numberKernels,
	newKernels=Min[optimalKernels-numberKernels,maxKernels-numberKernels];
	Print["Launch ",newKernels," kernels"];
	LaunchKernels[newKernels];
	]
	];	
	{time,substLogPart}=GenerateRelationFromRec[{rec,g[n],h2},{initialSubst,A},no,1,If[LogDeg===0,precision,Infinity]]//AbsoluteTiming;
];
Print["Time needed to get sub-log parts: ",time];

substLogN=substLogPart;
Do[
	substLogN[[k,2]]=N[substLogN[[k,2]],precision]/.A[b_]:>A[Round[b]],
{k,Length[substLogN]}];

(* substRelationL=Join[substRelationL,substLogN] *)
{substLogN,substLogPart,inhomExpr} ]


(* ::Section:: *)
(*Match points*)


(* ::Input::Initialization:: *)
FindBestInterval[list_,prec_]:=Module[{newlist,minval,p,p1,p2},
newlist=Abs[list];
minval=Min[Table[newlist[[i,2]],{i,1,Length[list]}]];
p=Position[newlist,minval][[1,1]];
{p1,p2}=Which[
p==1,{p,p+1},
p==Length[list],{p-1,p},
1<p<Length[list],{p-1,p+1}];
Map[Rationalize[#,prec]&,{list[[p1,1]],list[[p2,1]]}] ]


(* ::Input::Initialization:: *)
RefineMatchInterval[s_,freecoeffs_,testcoeff_,delta_,workingprecision1_,wmp_,{lpoint_,rpoint_},nstartpts_,prec_,iterations_]:=Module[{inipts,try,i,interval},
DistributeDefinitions[{freecoeffs,testcoeff,workingprecision1,wmp}];
inipts=GetIniPoints[lpoint,rpoint,nstartpts];
For[i=1,i<=iterations,i++,
try=ParallelMap[N[{#,MatchExpansions[Global`funcA,Global`funcB,s,freecoeffs,testcoeff,{#,delta},workingprecision1,wmp][[1]]},10]&,inipts];
interval=FindBestInterval[try,prec]; Print["Interval from interation ",i,": ",interval];
inipts=GetIniPoints[interval[[1]],interval[[2]],nstartpts] ];
interval ]


(* ::Input::Initialization:: *)
Options[GetBestPointMatch]={"Show interval"->"no"};

GetBestPointMatch[s_,freecoeffs_,testcoeff_,delta_,workingprecision1_,wmp_,{lpoint_,rpoint_},nstartpts_,prec_,iterations_,OptionsPattern[]]:=Module[{interval,inipts,list,newlist,minval,p,L},
interval=RefineMatchInterval[s,freecoeffs,testcoeff,delta,workingprecision1,wmp,{lpoint,rpoint},nstartpts,prec,iterations];
inipts=GetIniPoints[interval[[1]],interval[[2]],nstartpts];
list=ParallelMap[N[{#,MatchExpansions[Global`funcA,Global`funcB,s,freecoeffs,testcoeff,{#,delta},workingprecision1,wmp][[1]]},10]&,inipts]; Print[list];
newlist=Abs[list];
minval=Min[Table[newlist[[i,2]],{i,1,Length[list]}]];
p=Position[newlist,minval][[1,1]];
L={Rationalize[list[[p,1]],prec],list[[p,2]]};
If[OptionValue["Show interval"]=="yes",{interval,L},L] ]


(* ::Input::Initialization:: *)
BuildExpansion[z_,s_,a_,start_,nc_,nlogs_]:=Module[{expansion,i,j,k},
expansion=Sum[a[i,k]*z^k*Log[z]^i,{i,0,nlogs},{k,start,nc}];
ToStringForm[expansion,a]
]


(* ::Input::Initialization:: *)
PrepareHForm[fIn_,z_]:=
Module[{ord,A,f=fIn},
ord=Append[Union[Cases[{f},(Power[z,A_Integer]/;A<0)->A,Infinity]],0]//Min;
If[ord<0,
f=If[Head[f]===Plus,Apply[Plus,Apply[List,f]z^(-ord)],f z^(-ord)];
];
HornerFormH[CoefficientList[f,z]//Reverse,z]z^ord
];



(* ::Input::Initialization:: *)
BuildMatchExpansions[z_,s_,a_,{coeffsubsA_,solA_,startA_,ncA_,nlogsA_,ruleA_},{coeffsubsB_,startB_,ncB_,nlogsB_,ruleB_},{coeff_,j_,k_}]:=Module[{funcA,funcB,varA},
funcA=(BuildExpansion[z,s,a,startA,ncA,nlogsA]/. coeffsubsA /. solA);
funcA=Collect[funcA,Log[_],PrepareHForm[#,z]&];
funcA=funcA/. ruleA;

funcB=(BuildExpansion[z,s,a,startB,ncB,nlogsB] /. coeffsubsB)+coeff*z^j*Log[z]^k;
varA=Append[Table[ToExpression[ToString[a]~~ToString[kk]],{kk,0,nlogsB}],Log[_]];
funcB=Collect[funcB,varA,PrepareHForm[#,z]&];
funcB=funcB/. ruleB;
(* Output functions *)
{funcA,funcB}]


(* ::Input::Initialization:: *)
TestCoeffInfo[coeff_]:=Module[{i,j},
i=coeff[[1]];
j=ToExpression[StringSplit[ToString[Head[coeff]],"a"][[1]]];
{coeff,i,j} ]


(* ::Input::Initialization:: *)
EvHorner[plist_List,valIn_,prec_]:=
Module[{val=valIn},
If[!NumberQ[val],val=N[val,prec]];
Fold[(val #1+#2)&,0,plist]
];


MatchExpansions[func1_,func2_,s_,freecoeffs_,coeff_,{point_,delta_},workingprecision_,mwp_]:=Module[{slist,func1vals,lhs,eqsys,sol,i,j,k,l},
slist=Table[point+(j-1)*delta,{j,Length[freecoeffs]}];

func1vals=Map[(SetPrecision[func1/.s->#/.HornerFormH[A_,B_]:> EvHorner[A,B,workingprecision+200],workingprecision]//Expand)&,slist];
lhs=Map[(SetPrecision[func2 /. s->#/.HornerFormH[A_,B_]:> EvHorner[A,B,workingprecision+200],workingprecision+200]//Expand)&,slist] /.a_[r_]:>a[Round[r]];

eqsys=MapThread[Equal,{lhs,func1vals}];

eqsys=eqsys /. a_[r_]:>a[Round[r]];

sol=NSolve[eqsys,Join[freecoeffs],workingprecision-mwp];

If[sol==={},
{"Failed",{}},
sol=sol[[1]]/.A_[B_]:>A[Round[B]];
{coeff /. sol, sol}
]
]


(* ::Input::Initialization:: *)
Options[GetIniPoints]={"Drop"->0};

GetIniPoints[lpoint_,rpoint_,nstartpts_,OptionsPattern[]]:=Module[{T,k},
T=Table[lpoint+(rpoint-lpoint)/(nstartpts-1)*(k-1),{k,1,nstartpts}];
Drop[T,OptionValue["Drop"]]  ]


(* ::Section::Closed:: *)
(*End package*)


(* ::Input::Initialization:: *)
End[]
EndPackage[]
