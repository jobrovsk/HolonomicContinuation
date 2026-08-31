(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["RecToValuesFLINT`"]
ClearAll@@Names["RecToValuesFLINT`*"];


DeqToRecFlint::usage=
"DeqToRecFlint[{deq,g[z]},{a[n],p}]
gives the recurrence for a[n] of the solutions Sum[a[n](z-p)^n,{n,-\[Infinity],\[Infinity]}] for a rational point p. The differential equation is shifted ( /.z->z-p) and the recurrence is extracted.
The inhomogeneous part is not handled. Also, p must not be plus/minus infinity.
The options \"Details\", \"WriteOutputToFile\", \"PathToExecutible\" \"Overwrite\" work like for the function RecToValuesFLINT.
";


(* ::Input::Initialization:: *)
RecToValuesFLINT::usage="RecToValuesFLINT[{rec,g[n]},{smallestDerivative,ordD,a0},endComp,finalPrecision]
or
RecToValuesFLINT[{rec,g[n]},{initial,a0},endComp,finalPrecision]
rec: recurrence or a Path (String) (.m or .mx) containing a recurrence of the form 

	a_ordR g[n+l]+a_{ordR-1} g[n+l-1]+...+a_{0}g[n+t] +h ==0,

where a_i in Z[n] and t,l in Z with t<l and ordR=l-t. The h is ignored by default (if option \"Inhomogeneous\" not True).
g[n]: The function-symbol appearing in the recurrence

For giving the initial values, there are two options.
1) 0 is an ordinary point of the associated differential equation (here it is assumed that the recurrence was genereted by transforming this differential equation to the given recurrence rec)
	 a_tc \!\(\*SuperscriptBox[\(f\), TagBox[
RowBox[{\"(\", \"tc\", \")\"}],
Derivative],\nMultilineFunction->None]\)[z]+...+ a_ordD \!\(\*SuperscriptBox[\(f\), TagBox[
RowBox[{\"(\", \"ordD\", \")\"}],
Derivative],\nMultilineFunction->None]\)[z]
Then set {smallestDerivative,ordD}={tc,ordD}. Set a0 to the desired lhs of the output-substitutions
2) Else give \"initial\" as an explicit list of substitutions of a0[i]. The leading coefficient of the recurrence must not have any integer-zeros after the initial values. \"initial\" can be also given as a path to a file. Values smaller than the smallest given initial value are assumed to be 0.

endComp: index of the last value which should be computed (e.g. 15000)
finalPrecision: Relative precision in decimal digits. (has minimal influence on the running time)
	If finalPrecision is Infinity, then the output is exact.

Output: A list of substitutions for a0[startComp],...,a0[endComp], where the RHS are linear combinations of 
  1)a0[tc],...,a[ordD-1] and startComp=ordD .
  2)the variables in the rhs of \"initial\", and startComp is one larger than the given initial values.
 The coefficients have precision \"finalPrecision\". 

This function creates temporary files in Directory[]. Use SetDirectory[\"dir\"] to change this directory. The files are deleted if the function returns without error.

Options: 
\"Details\"-> (True/False), Whether Details should be printed. Default is no.
\"NumberOfThreads\"->n, Number of threads used by Flint. Default is 1.
\"Inhomogeneous\"->(True/False), Whether there is an inhomogeneous part in the equation (which is then evaluated via (h/.n->k) for k in the range implied by the initial values and endComp. Default is False i.e. inhomogeneous parts are ignored.
\"WriteOutputToFile\"->(False/PathToOutput), Whether output is returned (=default) or written to \"PathToOutput\" in the format implied by FileExtension[PathToOutput] (.m or .mx). The file must not exist.
\"Overwrite\" -> (True/False), Whether a file at PathToOutput should be overwritten. Default is no, then PathToOutput must not exist. Only relevant a path is specified via option \"WriteOutputToFile\" 
\"PathToExecutible\"->\"\", Give path to the directory, where the C-backend lives.
\"Backend\"->\"rec_to_val_V1\", Enter the C-backend which should be used. Curently rec_to_val_V1 and rec_to_val_V2 (faster, but needs more memory) are available.
";


(* ::Input::Initialization:: *)
Begin["`Private`"]


PathsToExecutibles=<||>;
BackendRecToValuesFLINT="rec_to_val_V2";


Clear[JPrint]
JPrint[debug_,stuff__]:=If[debug,Print[stuff]];


RecToValuesFLINT::execnotfound="Error: Executible `1` for `2` not found. Add path manually via option \"pathToExecutible\" .";
RecToValuesFLINT::filenotfound="Error: File `1` not found.";
RecToValuesFLINT::wrongextension="Error: File `1` must have extension \".m\" or \".mx\".";
RecToValuesFLINT::notempty="Error: `1` is not an empty directory.";
RecToValuesFLINT::cerror="Error: Program `1` did not terminate normally.";
RecToValuesFLINT::outputfileexists="Error: Output-file `1` already exists.";
RecToValuesFLINT::leadingcoeffvanish="Error: The leading coefficient of the recurrence is zero for some k in the range k>=`1`=startcomp.";
RecToValuesFLINT::inhomimpliesinitial="Error: If inhomogeneous recurrence is given, then initial values must be specified.";
RecToValuesFLINT::inhomtrivial="Warning: Option \"Inhomogeneous\"->True given, but recurrence seems to homogeneous. Computation continues with homogeneous recurrence.";
RecToValuesFLINT::initialtrivial="Warning: The given initial values are all zero.";
RecToValuesFLINT::recisdeq="Warning: The given recurrence contains Derivative. Are you sure this is a recurrence?";
RecToValuesFLINT::novars="Error: The right hand side of the given initial values is free of variables";
RecToValuesFLINT::startgeend="Warning: The given endComp=`1` is < startComp=`2`, i.e. the requested values are part of the initial values";
RecToValuesFLINT::cantwriteoutput="Error: Cannot write to Output-Path `1` ";
RecToValuesFLINT::optionmalformed="Error: Option[\"`1`\"] has malformed value `2`";
Options[RecToValuesFLINT]={
"PathToExecutible"->DirectoryName[$InputFileName],
"Overwrite"->False,
"Details"->False,
"NumberOfThreads"->1,
"Inhomogeneous"->False,
"WriteOutputToFile"->False,
"Backend"->"rec_to_val_V2"
};
Clear[RecToValuesFLINT]
RecToValuesFLINT[{recIn_,g_[n_Symbol]},{initialInfo__,a0_},endComp_Integer,finalPrecision_,OptionsPattern[]]:=Module[
	{inhomvalues,largesIntRoot,trivialInitial=False,trivialInhom=False,startComp,ordD,nVars,PathToOutput,k,debug=OptionValue["Details"],startInitial,hh,Rlcoeff,variables,precisionFlint,PathInhom,PathInitial,initialValues,
		trailingDcoeff,AAA,initialValGivenQ,outputPathGivenQ,inhomGivenQ,deq,rec,pathToExecutible,command,tempdir,PathOutputFlint,PathRec,PathTempfiles,tcoef,result,processOut},
initialValGivenQ=(Length[{initialInfo}]==1&&(Head[initialInfo]===List||StringQ[initialInfo]) );
outputPathGivenQ=(OptionValue["WriteOutputToFile"]=!=False);
inhomGivenQ=(OptionValue["Inhomogeneous"]=!=False);
If[inhomGivenQ&&!initialValGivenQ,Message[RecToValuesFLINT::inhomimpliesinitial];Abort[]];
If[OptionValue["PathToExecutible"]=!="",
	pathToExecutible=OptionValue["PathToExecutible"];
	If[DirectoryQ[pathToExecutible],pathToExecutible=FileNameJoin[{pathToExecutible,OptionValue["Backend"]}]];
,
	pathToExecutible=(FileNameJoin[{#,OptionValue["Backend"]}]&/@PathsToExecutibles)[SystemInformation["Kernel","MachineName"]];
];
If[Head[pathToExecutible]===Missing||!FileExistsQ[pathToExecutible],
	Message[RecToValuesFLINT::execnotfound,pathToExecutible , SystemInformation["Kernel","MachineName"]];
	Abort[];
];

If[outputPathGivenQ,
	If[!TrueQ[StringQ[OptionValue["WriteOutputToFile"]]],
		Message[RecToValuesFLINT::optionmalformed,"WriteOutputToFile",OptionValue["WriteOutputToFile"]];
		Abort[];
	];
	PathToOutput=OptionValue["WriteOutputToFile"];
	CheckFileDoesNotExist[PathToOutput,OptionValue["Overwrite"]];
	If[!MemberQ[{"m","mx"},FileExtension[PathToOutput]],Message[RecToValuesFLINT::wrongextension,PathToOutput];Abort[];];
	Export[PathToOutput,0];(*check if it is possible to write to PathToOutput*)
	If[!FileExistsQ[PathToOutput],Message[RecToValuesFLINT::cantwriteoutput,PathToOutput];Abort[];];
	DeleteFile[PathToOutput];
];
If[StringQ[recIn]&&!FileExistsQ[recIn],Message[RecToValuesFLINT::filenotfound,recIn];Abort[];];

tempdir=CreateDirectory["temp_"<>ToString[$KernelID]<>"_"<>ToString[RandomInteger[10^12-1]]];
JPrint[debug,"Temporary directory: "<>tempdir];

PathOutputFlint=FileNameJoin[{tempdir,"rec_to_val.data"}];
PathRec=FileNameJoin[{tempdir,"recM.data"}];
PathTempfiles=FileNameJoin[{tempdir,OptionValue["Backend"]}];
If[FileExistsQ[PathRec]||FileExistsQ[PathOutputFlint],Message[RecToValuesFLINT::notempty,tempdir];Abort[];];

JPrint[debug,"Imported and written Rec: ",AbsoluteTiming[
rec=If[StringQ[recIn],Import[recIn],recIn];
If[!FreeQ[rec,Derivative],Message[RecToValuesFLINT::recisdeq]];
rec=If[Head[rec]===Equal,rec[[1]]-rec[[2]],rec];
Rlcoeff=WriteRecToFileFlint[rec,g[n],PathRec];
][[1]]];
If[inhomGivenQ&&((rec/.g[_]->0)===0),inhomGivenQ=False;Message[RecToValuesFLINT::inhomtrivial];];
precisionFlint=If[finalPrecision===Infinity,0,Ceiling[Log2[10]*finalPrecision]];
command={pathToExecutible,PathOutputFlint,PathRec,"deqdummy",PathTempfiles,nVars,startComp,endComp,precisionFlint,OptionValue["NumberOfThreads"]};
If[initialValGivenQ,
	JPrint[debug,"Imported and written initial values:  ",AbsoluteTiming[
	PathInitial=FileNameJoin[{tempdir,"initialM.data"}];
	initialValues=If[StringQ[initialInfo],Import[initialInfo],initialInfo];	
	variables=(initialValues//Values//Variables);
	nVars=Length[variables];
	{startInitial,startComp}=MinMax[Cases[Join[Keys[initialValues],variables],a0[AAA_]->AAA,Infinity]];
	startComp+=1;
	If[startComp>endComp,Message[RecToValuesFLINT::startgeend,endComp,startComp];Return[{}]];
	If[inhomGivenQ,
		hh=(rec/.g[_]->0);
		inhomvalues=Table[hh/.n->(k-Rlcoeff),{k,startComp,endComp}];
		variables=Union[variables,Variables[inhomvalues]];
		nVars=Length[variables];
	];
	If[nVars==0,Message[RecToValuesFLINT::novars];Abort[];];
	trivialInitial=WriteRationalsFlint[PathInitial,Table[a0[k]/.initialValues,{k,startInitial,startComp-1}],variables,startInitial];
	][[1]]];
	If[trivialInitial,Message[RecToValuesFLINT::initialtrivial];];
	AppendTo[command,PathInitial];
	If[inhomGivenQ,
		JPrint[debug,"Written inhomogeneous part:  ",AbsoluteTiming[
		
		PathInhom=FileNameJoin[{tempdir,"inhomM.data"}];
		(*JPrint[debug,"inhom: ",hh];*);
		trivialInhom=WriteRationalsFlint[PathInhom,inhomvalues,variables,startComp];
		If[trivialInhom,
			inhomGivenQ=False;
			Message[RecToValuesFLINT::inhomtrivial];
			DeleteFile[PathInhom];
		,
			AppendTo[command,PathInhom];
		];
		][[1]]];
	];
,
	{trailingDcoeff,ordD}={initialInfo};
	nVars=ordD-trailingDcoeff;
	startComp=ordD;
	variables=Table[a0[i],{i,trailingDcoeff,ordD-1}];
];
largesIntRoot=Max[Append[SolveValues[Coefficient[rec,g[Rlcoeff+n]]==0,n,Integers],-Infinity]]+Rlcoeff;
If[largesIntRoot>=startComp,Message[RecToValuesFLINT::leadingcoeffvanish,startComp];Abort[];];
JPrint[debug,"Run: ",StringRiffle[command]];
JPrint[debug,"Running "<>FileBaseName[pathToExecutible]<>" ",AbsoluteTiming[
processOut=RunProcess[command];
][[1]]];
If[processOut["ExitCode"]=!=0,
        Message[RecToValuesFLINT::cerror,pathToExecutible];
        Print[If[StringLength[processOut["StandardOutput"]]>6000,
			StringTake[processOut["StandardOutput"] ,;;5000]<>"\n....\n"<>StringTake[processOut["StandardOutput"] ,-1000;;],processOut["StandardOutput"]]];
		Return[processOut];
];
(*Print[processOut["StandardOutput"]];*)
DeleteFile[PathRec];
JPrint[debug,"Readed in output: ",AbsoluteTiming[
	result=ReadOutputFlint[PathOutputFlint,finalPrecision];
][[1]]];
DeleteFile[PathOutputFlint];
If[initialValGivenQ,DeleteFile[PathInitial];];
If[inhomGivenQ,DeleteFile[PathInhom];];
DeleteDirectory[tempdir];
If[outputPathGivenQ,
	CheckFileDoesNotExist[PathToOutput,OptionValue["Overwrite"]];
	JPrint[debug,"Exported final output: ",AbsoluteTiming[
	Export[PathToOutput,Table[a0[j]->result[[j-startComp+1]] . variables,{j,startComp,endComp}]];
	][[1]]];
,
	Return[Table[a0[j]->result[[j-startComp+1]] . variables,{j,startComp,endComp}]];
];
Return[0];
]


DeqToRecFlint::execnotfound="Error: Executible `1` for `2` not found. Add path manually via option \"pathToExecutible\" .";
DeqToRecFlint::optionmalformed="Error: Option[\"`1`\"] has malformed value `2`";
DeqToRecFlint::wrongextension="Error: File `1` must have extension \".m\" or \".mx\".";
DeqToRecFlint::filenotfound="Error: File `1` not found.";
DeqToRecFlint::outputfileexists="Error: Output-file `1` already exists.";
DeqToRecFlint::cerror="Error: Program `1` did not terminate normally.";

DeqToRecFlint::notempty="Error: `1` is not an empty directory.";

DeqToRecFlint::cantwriteoutput="Error: Cannot write to Output-Path `1` ";

Options[DeqToRecFlint]={
"PathToExecutible"->DirectoryName[$InputFileName],
"Overwrite"->False,
"Details"->False,
"WriteOutputToFile"->False
};
BackendDeqToRecFlint="deq_to_rec";
Clear[DeqToRecFlint];
DeqToRecFlint[{deqIn_,g_[z_Symbol]},{r_[n_Symbol],point_?NumberQ},OptionsPattern[]]:=Module[
	{trailingDcoeff,ordR,shiftR,PathDeq,ordD,PathToOutput,k,debug=OptionValue["Details"],AAA,outputPathGivenQ,
		deq,pathToExecutible,command,tempdir,PathOutputFlint,result,processOut},
outputPathGivenQ=(OptionValue["WriteOutputToFile"]=!=False);
If[OptionValue["PathToExecutible"]=!="",
	pathToExecutible=OptionValue["PathToExecutible"];
	If[DirectoryQ[pathToExecutible],pathToExecutible=FileNameJoin[{pathToExecutible,OptionValue["Backend"]}]];
,
	pathToExecutible=(FileNameJoin[{#,BackendDeqToRecFlint}]&/@PathsToExecutibles)[SystemInformation["Kernel","MachineName"]];
];
If[Head[pathToExecutible]===Missing||!FileExistsQ[pathToExecutible],
	Message[DeqToRecFlint::execnotfound,pathToExecutible , SystemInformation["Kernel","MachineName"]];
	Abort[];
];

If[outputPathGivenQ,
	If[!TrueQ[StringQ[OptionValue["WriteOutputToFile"]]],
		Message[DeqToRecFlint::optionmalformed,"WriteOutputToFile",OptionValue["WriteOutputToFile"]];
		Abort[];
	];
	PathToOutput=OptionValue["WriteOutputToFile"];
	CheckFileDoesNotExist[PathToOutput,OptionValue["Overwrite"]];
	If[!MemberQ[{"m","mx"},FileExtension[PathToOutput]],Message[DeqToRecFlint::wrongextension,PathToOutput];Abort[];];
	Export[PathToOutput,0];(*check if it is possible to write to PathToOutput*)
	If[!FileExistsQ[PathToOutput],Message[DeqToRecFlint::cantwriteoutput,PathToOutput];Abort[];];
	DeleteFile[PathToOutput];
];
If[StringQ[deqIn]&&!FileExistsQ[deqIn],Message[DeqToRecFlint::filenotfound,deqIn];Abort[];];

tempdir=CreateDirectory["temp_"<>ToString[$KernelID]<>"_"<>ToString[RandomInteger[10^12-1]]];
JPrint[debug,"Temporary directory: "<>tempdir];

PathOutputFlint=FileNameJoin[{tempdir,"recFlint.data"}];
PathDeq=FileNameJoin[{tempdir,"deqM.data"}];
PathTempfiles=FileNameJoin[{tempdir,BackendDeqToRecFlint}];
If[FileExistsQ[PathDeq]||FileExistsQ[PathOutputFlint],Message[DeqToRecFlint::notempty,tempdir];Abort[];];

JPrint[debug,"Imported and written Deq: ",AbsoluteTiming[
deq=If[StringQ[deqIn],Import[deqIn],deqIn];
deq=If[Head[deq]===Equal,deq[[1]]-deq[[2]],deq];
{trailingDcoeff,ordD}=WriteDeqToFileFlint[deq,g[z],PathDeq];
][[1]]];

command={pathToExecutible,PathOutputFlint,PathDeq,Numerator[point],Denominator[point]};
JPrint[debug,"Run: ",StringRiffle[command]];
JPrint[debug,"Running "<>FileBaseName[pathToExecutible]<>" ",AbsoluteTiming[
processOut=RunProcess[command];
][[1]]];
If[processOut["ExitCode"]=!=0,
        Message[DeqToRecFlint::cerror,pathToExecutible];
        Print[If[StringLength[processOut["StandardOutput"]]>6000,
			StringTake[processOut["StandardOutput"] ,;;5000]<>"\n....\n"<>StringTake[processOut["StandardOutput"] ,-1000;;],processOut["StandardOutput"]]];
		Return[processOut];
];
(*Print[processOut["StandardOutput"]];*)
DeleteFile[PathDeq];
JPrint[debug,"Readed in output: ",AbsoluteTiming[
	{result,ordR,shiftR}=ReadOutputFlintPoly[PathOutputFlint,n];
][[1]]];
DeleteFile[PathOutputFlint];
DeleteDirectory[tempdir];
result=((result . Table[r[i+n],{i,-Length[result]+1+shiftR,shiftR}](*-h[n-ordD]*))==0);
If[outputPathGivenQ,
	CheckFileDoesNotExist[PathToOutput,OptionValue["Overwrite"]];
	JPrint[debug,"Exported final output: ",AbsoluteTiming[
	Export[PathToOutput,result];
	][[1]]];
	Return[0];
,
	Return[result];
];
]


Clear[WriteRationalsFlint]
WriteRationalsFlint::fileexists="Warning: The file `1` already exists, no overwrite performed.";
WriteRationalsFlint[PathOut_String,Rat_?VectorQ,variables_?VectorQ,start_Integer]:=Module[{trivial,numDen,File,i,j},
trivial=(Length[Intersection[variables,Variables[Rat]]]==0);
If[FileExistsQ[PathOut],Message[WriteRationalMatrixFlint::fileexists,PathOut];Return[1]];
File=OpenWrite[PathOut];
Write[File,start];
Do[
	Do[
		numDen=NumeratorDenominator[Coefficient[Rat[[i]],variables[[j]]]];
		WriteString[PathOut,ToString[numDen[[1]]]," ",ToString[numDen[[2]]],"\n"];
	,{j,Length[variables]}];
,{i,1,Length[Rat]}];
Close[File];Return[trivial];

]


Clear[CheckFileDoesNotExist]
CheckFileDoesNotExist[filename_,overwriteQ_]:=If[FileExistsQ[filename]&&!TrueQ[overwriteQ],
	Message[RecToValuesFLINT::outputfileexists,filename];
	Abort[];
];


Clear[WriteDeqToFileFlint];
WriteDeqToFileFlint[deq_,g_[z_Symbol],WritePathRec_String]:=Module[
{trailing,AAA,ordD,i},
{trailing,ordD}=MinMax[Join[Cases[deq, Derivative[AAA_][g][z] -> AAA, All], 
   Cases[deq, g[z] -> 0, All]]];
WritePolyListToFileFlint[Table[Coefficient[deq,Derivative[i][g][z]],{i,0,ordD}],trailing,WritePathRec,z];
Return[{trailing,ordD}];
]


Clear[WriteRecToFileFlint];
WriteRecToFileFlint[rec_,g_[n_Symbol],WritePathRec_String]:=Module[
{ordR,Rlcoeff,AAA,Rtcoeff,deqlist,ordD,tcoef},
{Rtcoeff,Rlcoeff}=MinMax[Cases[rec,g[AAA_]->AAA-n,Infinity]];
WritePolyListToFileFlint[Table[Coefficient[rec,g[i+n]],{i,Rtcoeff,Rlcoeff}],Rlcoeff,WritePathRec,n];
Return[Rlcoeff];
]


Clear[WritePolyListToFileFlint]
WritePolyListToFileFlint::fileexists="Warning: The file `1` already exists, no overwrite performed.";
WritePolyListToFileFlint::notinteger="Error: The coefficients must be polynomials in Z[`1`].";
WritePolyListToFileFlint[pols_?VectorQ,coeff_Integer,PathOut_String,z_Symbol]:=Module[
{File,ord=Length[pols]-1},
If[!FreeQ[pols,Rational],Message[WritePolyListToFileFlint::notinteger,z];Abort[];];
If[FileExistsQ[PathOut],Message[WritePolyListToFileFlint::fileexists,PathOut];Return[1]];
File=OpenWrite[PathOut];
Write[File,ord];
Write[File,coeff];
Do[
	WriteLine[File,ToString[
		Max[Exponent[pols[[i]],z]+1,0]]<>"  "<>
		StringTake[StringReplace[ToString[CoefficientList[pols[[i]],z]],","->""],2;;-2]]
,{i,1,ord+1}];
Close[File];
Return[0];
]


Clear[WriteInitialToFileFlint]
WriteInitialToFileFlint::fileexists="Warning: The file `1` already exists, no overwrite performed.";
WriteInitialToFileFlint[PathOut_String,ordD_Integer,tcoef_Integer]:=Module[{File,initial},
If[FileExistsQ[PathOut],Message[WriteInitialToFileFlint::fileexists,PathOut];Return[1]];
File=OpenWrite[PathOut];
Write[File,0];(*index of first initial value*)
(*Write[File,ordD];(*index of first value to compute*)
Write[File,ordD-tcoef];(*number of setes initial values*)*)
initial=IdentityMatrix[ordD][[;;,tcoef+1;;]];
Print[Dimensions[initial]];
Do[
	WriteLine[File,StringReplace[ToString[initial[[i]]],{"{"->"" , "}"->" 1",","->" 1"}]]
,{i,1,ordD}];
Close[File];
Return[0];
]




Clear[ReadOutputFlint]
ReadOutputFlint[PathOutput_String,finalPrecision_]:=Module[{file,result={},line,i},
file = OpenRead[PathOutput];
result = {};
While[True,
  line = ReadLine[file];
  If[line===EndOfFile,Break[]];
  line = ToExpression/@ StringSplit[line];
  AppendTo[result, 
  If[(finalPrecision===Infinity)||(finalPrecision<=0),
   line[[;;-2]]/line[[-1]]
   ,
   Table[N[line[[2 i - 1]],finalPrecision] 2^line[[2 i]],{i, Length[line]/2}]
   ]
  ];
 ];
Close[file];
Return[result];
]


Clear[ReadOutputFlintPoly];
ReadOutputFlintPoly[PathOutput_String,n_Symbol]:=Module[{file,result={},line,i,ordR,shiftR},
file = OpenRead[PathOutput];
ordR=Read[file,Number];
shiftR=Read[file,Number];
result = {};
Do[
  line = ReadLine[file];
  If[line===EndOfFile,Abort[]];
  line = ToExpression/@ StringSplit[line];
  AppendTo[result,FromDigits[Reverse[line[[2;;]]],n]];
 ,{ordR+1}];
Close[file];
Return[{result,ordR,shiftR}];
]


(* ::Input::Initialization:: *)
End[]
EndPackage[]
