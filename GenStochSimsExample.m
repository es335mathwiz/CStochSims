(* Wolfram Language Package *)

BeginPackage["GenStochSimsExample`", {"ProtectedSymbols`", "AMAModel`","AccelerateAMA`","Stack`","Stoch`","MmaModelToC`","JLink`"}]


cnstrctDModel::usage="cnstrctDModel[inFName_String,outFName_String]"
getCoeffs::usage="getCoeffs[inFile_String]"
getVarNames::usage="getVarNames[inFile_String]"




Begin["Private`"] (* Begin Private Context *) 
getCoeffs[inFile_String]:=
Module[{ds,conn,exp,result,sTemplate=StringTemplate[
"<modCoefficients>{
let $mod:=doc('`1`')
for $var in $mod//model/variable
return 
if ($var/mce_equation/coeff)
then (
for $coeffMCE in $var/mce_equation/coeff 
let $nmMCE := $coeffMCE/cf_name/text()
let $vlMCE := $coeffMCE/cf_value/number()
return  <coeffInfo theName=\"{$nmMCE}\" theValue=\"{$vlMCE}\"/>)
else(
for $coeff in $var/standard_equation/coeff 
let $nm := $coeff/cf_name/text()
let $vl := $coeff/cf_value/number()
return  <coeffInfo theName=\"{$nm}\" theValue=\"{$vl}\"/>)
}</modCoefficients>"]},
ds=JavaNew["com.saxonica.xqj.SaxonXQDataSource"];
conn = ds[getConnection[]];
exp = conn[prepareExpression[
sTemplate[inFile]]];
result = exp[executeQuery[]];
elemsToAssns[ImportString[
sequenceToListString[result][[1]],"XML"]
]]


elemsToSubs[xml_]:=
Cases[xml,XMLElement["coeffInfo", {"theName" -> tn_String, 
        "theValue" -> tv_String}, {}]->(ToExpression[tn]->ToExpression[tv]),Infinity]


elemsToAssns[xml_]:=
Cases[xml,XMLElement["coeffInfo", {"theName" -> tn_String, 
        "theValue" -> tv_String}, {}]:>{tn,(tn<>"="<>tv)},Infinity]



theNextString[aSeq_?JavaObjectQ]:=
With[{prime=aSeq[next[]]},If[prime,aSeq[getItemAsString[Null]],"$noMore"]]

sequenceToListString[aSeq_?JavaObjectQ]:=
Drop[Drop[FixedPointList[theNextString[aSeq]&,xx],1],{-2,-1}]


theNextNumber[aSeq_?JavaObjectQ]:=
With[{prime=aSeq[next[]]},If[prime,aSeq[getDouble[]],"$noMore"]]

sequenceToListNumber[aSeq_?JavaObjectQ]:=
Drop[Drop[FixedPointList[theNextNumber[aSeq]&,xx],1],{-2,-1}]


getVarNames[inFile_String]:=
Module[{ds,conn,exp,result,sTemplate=StringTemplate[
"let $mod:=doc('`1`')
for $var in $mod//model/variable/name/text() return $var"]},
ds=JavaNew["com.saxonica.xqj.SaxonXQDataSource"];
conn = ds[getConnection[]];
exp = conn[prepareExpression[
sTemplate[inFile]]];
result = exp[executeQuery[]];
sequenceToListString[result]
]

nameSubs[aStr_String]:=(ridParens[aStr]->ridUSSubs[ridParens[aStr]])

nameSubsWP[aStr_String]:=(aStr->ridUSSubs[ridParens[aStr]])


ridParens[str_String]:=StringReplace[str,{"("->"_",")"->""}]

ridUSSubs[str_String]:=StringReplace[
StringReplace[str,{"("->"_",")"->""}],"_"->"US"]





getEVEqns[inFile_String]:=
Module[{ds,conn,exp,result,sTemplate=StringTemplate[
"{let $mod:=doc('`1`')
for $eqn in $mod//model/variable/mce_equation/eviews_equation/text() return $eqn}"]},
ds=JavaNew["com.saxonica.xqj.SaxonXQDataSource"];
conn = ds[getConnection[]];
exp = conn[prepareExpression[
sTemplate[inFile]]];
result = exp[executeQuery[]];
sequenceToListString[result]
]


getDynareEqns[inFile_String]:=
Module[{ds,conn,exp,result,sTemplate=StringTemplate[
"let $mod:=doc('`1`')
for $var in $mod//model/variable 
return (if ($var/mce_equation) then $var/mce_equation/dynare_equation/text() else if ($var/standard_equation/dynare_equation) then $var/standard_equation/dynare_equation/text() else $var/name/text())
"]},
ds=JavaNew["com.saxonica.xqj.SaxonXQDataSource"];
conn = ds[getConnection[]];
exp = conn[prepareExpression[
sTemplate[inFile]]];
result = exp[executeQuery[]];
sequenceToListString[result]
]



$frbLocDesktop="g:/.m2/repository/";
$frbLoc="/msu/home/m1gsa00/.m2/repository/";
$garyMacLoc="/Users/garyanderson/.m2/repository/";
$reposLoc=
Which[
FileExistsQ[$frbLoc],$frbLoc,
FileExistsQ[$frbLocDesktop],$frbLocDesktop,
FileExistsQ[$garyMacLoc],$garyMacLoc
];


AddToClassPath[FileNameJoin[{$reposLoc,"gov/frb/ma/msu/dynareAntlr/1.0/dynareAntlr-1.0.jar"}]]
AddToClassPath[FileNameJoin[{$reposLoc,"org/antlr/antlr/3.5.2/antlr-3.5.2.jar"}]]
AddToClassPath[FileNameJoin[{$reposLoc,"org/antlr/antlr-runtime/3.5.2/antlr-runtime-3.5.2.jar"}]]
AddToClassPath[FileNameJoin[{$reposLoc,"SAXON-HE/saxon9-xqj/9-7-0-2/saxon9-xqj-9-7-0-2.jar"}]]
AddToClassPath[FileNameJoin[{$reposLoc,"SAXON-HE/saxon9he/9-7-0-2/saxon9he-9-7-0-2.jar"}]]
AddToClassPath[FileNameJoin[{$reposLoc,"xalan-j_2_7_1/xalan/2-7-1/xalan-2-7-1.jar"}]]




cnstrctDModel[inFName_String,outFName_String]:=
With[{eqns=getDynareEqns[inFName],
vars=getVarNames[inFName],
coeffs=getCoeffs[inFName]},
With[{addFctrs=#<>"_aerr"&/@vars},
With[{allVars=Join[addFctrs,vars]},
With[{vcSubs=
DeleteCases[Join[(nameSubs/@allVars),
(nameSubsWP/@(First/@coeffs)),
(nameSubs/@(First/@coeffs))],xx_->xx_]},
{eqns,vars,coeffs,doModStr[allVars,coeffs,eqns,vcSubs,inFName,outFName]}]]]]





doModStr[allVars_List,coeffs_List,eqns_List,vcSubs_List,inFName_String,outFName_String]:=
With[{modStr=StringReplace[
"var " <> (StringJoin @@ Riffle[allVars," "]) <> ";\n" <>
ridParens["parameters " <> (StringJoin @@ Riffle[First/@coeffs," "]) <> 
";\n"] <>
ridParens[(StringJoin @@ Riffle[Last/@coeffs,";\n"]) <> ";\n"] <>
"model;\n"<> (StringJoin @@ Riffle[eqns,";\n"]) <>";\n end;\n",vcSubs]},
WriteString[outFName,modStr]]




End[] (* End Private Context *)

EndPackage[]
