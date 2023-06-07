(* ::Package:: *)

(*  This function prints the version number for the installed verison of this package.  *)

PhylogeneticsVersion[]:=Print["Phylogenetics for Mathematica 6.8\n(c) P. David Polly, 7 June 2022\n"];


PhylogeneticsVersion[]


(* This functions reads one or more Newick format trees from a text file.  

   Usage:    ReadNewick[FileName] 

             where FilePath is the name (possibly including path) to the tree file.

	Created by P. David Polly, 21 March 2012
	Updated 14 December 2016 to set branch lengths to 1 if they are not already set and to remove PHYLIPs tree weights if they exist
*)

ReadNewick[TreeFile_]:=Block[{Tree,NumTrees,NumNodes,NumBranchLengths,NumTips},

Tree=Import[TreeFile,"Text"];
Tree=StringReplace[Tree,RegularExpression["\[[^\]]+]"]->""];
NumTrees=Length[StringCases[Tree,";"]];
If[NumTrees==0,Tree=Tree<>";"];
NumTrees=Length[StringCases[Tree,";"]];
If[NumTrees>1,(Tree=#<>";"&/@StringSplit[Tree,";"]);
Do[
Tree[[x]]=StringReplace[Tree[[x]],{RegularExpression["[^\\x21-\\x7E]"]->""}];
NumNodes=Length[StringCases[Tree[[x]],"("]];
NumBranchLengths=Length[StringCases[Tree[[x]],":"]];
If[NumBranchLengths==0,Tree[[x]]=StringReplace[Tree[[x]],{")"->":1)",","->":1,",";"->":1;"}];Print["All branch lengths set to 1"];NumBranchLengths=Length[StringCases[Tree[[x]],":"]]];
NumTips=(Length[StringCases[Tree[[x]],","]])+1;
If[NumNodes>NumBranchLengths,Print["Your tree does not have lengths for all branches"]],
{x,Length[Tree]}],

Tree=StringReplace[Tree,{RegularExpression["[^\\x21-\\x7E]"]->""}];
NumNodes=Length[StringCases[Tree,"("]];
NumBranchLengths=Length[StringCases[Tree,":"]];
If[NumBranchLengths==0,Tree=StringReplace[Tree,{")"->":1)",","->":1,",";"->":1;"}];Print["All branch lengths set to 1"];NumBranchLengths=Length[StringCases[Tree,":"]]];
NumTips=(Length[StringCases[Tree,","]])+1;
If[NumNodes>NumBranchLengths,Print["Your tree does not have lengths for all branches"]]
];

Return[Tree];
];



(* This function transforms a Newick format tree into a table where each row represents
   a branch in the tree, the first column giving the name of the branch end point, the
   second column giving the name of the branch ancestor, the third column giving the
   length of the branch, and the fourth column indicating whether the branch ends in a 
   tip (1) or a node (0).  The table is returned with a header line that must be discarded
   for functions that use the table as input. 

   Usage:    TreeToTable[tree] 

             where tree is a Newick format tree in a single string.  The tree must have
             branch lengths and must end in a semicolon. 
			 Example Tree:  "((A:1,B:1):1,C:1);"
 
	Created by P. David Polly, 21 March 2012
	Updated 14 December 2016 to fix new incompatibility with NumberForm formatting
	
*)

TreeToTable[Tree_]:=Block[{LastNode,SplitTree,Pos,TreeTable,Tip,Node,ParseTree},
(* This subfunction processes a tip, reading the name and length and inserting it into TreeTable *)

Tip[AncNum_]:=Block[{MyName, MyNum, MyAnc, MyLength},
MyAnc="Node "<>ToString[AncNum];
MyName="";
MyLength="";

While[StringMatchQ[SplitTree[[Pos]], RegularExpression["[^:]"]],
MyName=MyName<>SplitTree[[Pos]];
Pos=Pos+1;
];
Pos=Pos+1;
While[StringMatchQ[SplitTree[[Pos]], RegularExpression["[^,);]"]],
MyLength=MyLength<>SplitTree[[Pos]];
Pos=Pos+1;
];
MyLength=ToExpression[MyLength];
TreeTable=Append[TreeTable,{MyName,MyAnc,MyLength, 1}];
Pos=Pos-1;
Return[];
]; 

(* This subfunction processes a node, calling new node functions or tips as they are encountered within the node *)
Node[AncNum_]:=Block[{MyName, MyAnc,MyLength,MyNum},
MyLength="";
MyNum=LastNode;
LastNode=LastNode+1;
MyName="Node "<>ToString[MyNum];
MyAnc="Node "<>ToString[AncNum];

While[StringMatchQ[SplitTree[[Pos]],RegularExpression["[^);]+"]],
Pos=Pos+1;
If[StringMatchQ[SplitTree[[Pos]],"("],Node[MyNum],
If[StringMatchQ[SplitTree[[Pos]], RegularExpression["[A-Za-z0-9]"]],Tip[MyNum],
If[StringMatchQ[SplitTree[[Pos]],RegularExpression["[);]"]],Break[]]]];
];
Pos=Pos+1;
If[StringMatchQ[SplitTree[[Pos]],":"],Pos=Pos+1];
While[StringMatchQ[SplitTree[[Pos]], RegularExpression["[^,);]+"]],
MyLength=MyLength<>SplitTree[[Pos]];
Pos=Pos+1;
];
MyLength=ToExpression[MyLength];
TreeTable=Append[TreeTable,{MyName,MyAnc,MyLength,0}];
Return[];
];

(* main code *)
TreeTable={{"Descendant","Ancestor","Branch Length", "Tip?"}};
Pos=1;
LastNode=0;
SplitTree=StringSplit[Tree,""];
While[True, 
If[StringMatchQ[SplitTree[[Pos]],"("],Node[LastNode]];
If[StringMatchQ[SplitTree[[Pos]],";"],Break[]];
];
Return[TreeTable[[1;;-2]]];
];



(* This function calculates the three matrices used by Martins and Hansen (1997) for phylogenetic regressions.
   The matrices are:  (1)  Var[Y], which is the matrix of shared ancestry of tree tips;  (2) Var[AY], which is 
   the matrix of shared ancestry of tips; and (3) Var[A], which is the matrix of shared ancestry of tree nodes; and
   and nodes.  This function takes a tree table as its argument, in the form produced by the function TreeToTable[]
   (but without the row of column labels!).  The three matrices are returned with column and row labels.


   Usage:    {varY,varAY,varA}=PhylogeneticMatrices[tree];
 
             where varY, varAY, and varA are the three output matrices described above and tree is a string 
             containing a tree with branch lengths in Newick format or a tree table of  the sort produced by 
             TreeToTable[].
 
	Created by P. David Polly, updated 20 May 2012
	
*)

PhylogeneticMatrices[tree_]:=Block[{FindLineage,TipEntries,NodeEntries,TipLineages,NodeLineages,varY,varA,varAY,i,j,myTable},

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];

FindLineage[Point_]:=Block[{Lineage,NewBranch},
Lineage={Point};
While[!StringMatchQ[Lineage[[-1,2]],"Node 0"],
NewBranch=Flatten[myTable[[Flatten[Position[myTable[[1;;,1]],Lineage[[-1,2]]]]]]];
Lineage=Append[Lineage,NewBranch];
];
Return[Lineage];
];

TipEntries=myTable[[Flatten[Position[StringCases[myTable[[1;;,1]],"Node "~~__],{}]]]];
NodeEntries=myTable[[Flatten[Position[myTable[[1;;,1]],#]&/@Flatten[StringCases[myTable[[1;;,1]],"Node "~~__]]]]];
TipLineages=FindLineage[#]&/@TipEntries;
NodeLineages=FindLineage[#]&/@NodeEntries;
varY=Table[Table[Plus@@(Intersection[TipLineages[[i]],TipLineages[[j]]][[1;;,3]]),{i,Length[TipLineages]}],{j,Length[TipLineages]}];
varA=Table[Table[Plus@@(Intersection[NodeLineages[[i]],NodeLineages[[j]]][[1;;,3]]),{i,Length[NodeLineages]}],{j,Length[NodeLineages]}];
varA=Transpose[Prepend[Transpose[Prepend[varA,Table[1,{Length[varA[[1]]]}]]],Table[1,{Length[varA[[1]]]+1}]]];
varAY=Table[Table[Plus@@(Intersection[TipLineages[[i]],NodeLineages[[j]]][[1;;,3]]),{i,Length[TipLineages]}],{j,Length[NodeLineages]}];
varAY=Prepend[varAY,Table[0,{Length[varAY[[1]]]}]];
Return[{Prepend[Transpose[Prepend[varY,TipEntries[[1;;,1]]]],Prepend[TipEntries[[1;;,1]],"Var Y"]],Transpose[Prepend[Transpose[Prepend[varAY,TipEntries[[1;;,1]]]],Prepend[Prepend[NodeEntries[[1;;,1]],"Node 0"],"Var AY"]]],Prepend[Transpose[Prepend[varA,Prepend[NodeEntries[[1;;,1]],"Node 0"]]],Prepend[Prepend[NodeEntries[[1;;,1]],"Node 0"],"Var A"]]}];
];




(*  This function reconstructs ancestral quantitative values for traits on a phylogenetic tree following Martins & Hansen,
    1997, Am. Nat., 149:646-667.  Estimate of the standard error at the base of the tree is from Rohlf, 2001, Evolution, 
    55: 2143-2160.


    Usage:  {rate, nodes, se} = NodeReconstruct[tree, trait]
            where tree is a string containing a tree with branch lengths in Newick format or a tree table of the sort produced 
            by TreeToTable[].  The function returns the estimated per-step rate of the trait on the tree, the ancestral 
            values of Y at each node in the tree assuming a Brownian motion, and the uncertainty in the ancestral values due to 
            the Brownian motion evolutionary process expressed as 1 standard deviation of the variance due to the evolutionary 
            process (the "standard error" of some authors).  

    	Created by P. David Polly, updated 20 May 2012

*)

ReconstructNodes[tree_,trait_]:=Block[{J,GrandMean,Yprime,Nodes,Rate,x,contrasts, i,j,SDNode,X,SyGLS,varY,varAY,varA,TipLabels,NodeLabels,Y,VarY,VarAY,VarA,myTable,RateSquared,VarYNew},

If[StringQ[tree],myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];
];
];


{VarY,VarAY,VarA}=PhylogeneticMatrices[myTable];
varY=VarY[[2;;,2;;]];varAY=VarAY[[2;;,2;;]]; varA=VarA[[2;;,2;;]];
TipLabels=VarY[[2;;,1]];NodeLabels=VarA[[2;;,1]];

If[Length[trait[[1]]]!=2,Return["Trait must be a two-column matrix with tip labels in the first column and trait values in the second column"]];
If[Sort[TipLabels]!=Sort[trait[[1;;,1]]],Return["Traits labels do not match tip labels from Var Y."]];
Y=trait[[1;;,2]][[Flatten[Position[trait[[1;;,1]],#]&/@TipLabels]]];

J=Table[1,{Length[varY]}];
GrandMean=(J . Inverse[varY] . Y)/(J . Inverse[varY] . J);
Yprime=#-GrandMean&/@Y;
Nodes=#+GrandMean&/@(varAY . Inverse[varY] . Yprime);

contrasts=IndependentContrasts[tree,trait];
RateSquared=Mean[contrasts^2];
Rate=Sqrt[RateSquared];

SDNode={};

Do[
VarYNew=PhylogeneticMatrices[RerootTree[myTable,NodeLabels[[i]]]][[1,2;;,2;;]];
SyGLS=Sqrt[RateSquared/(Plus@@Plus@@Inverse[VarYNew])]//N;
SDNode=Append[SDNode,{NodeLabels[[i]],SyGLS}]
,{i,Length[NodeLabels]}];

Return[{{"Rate",Rate},Transpose[Prepend[{Nodes},#&/@NodeLabels]],SDNode}];

];





(* A function to estimate all the node values using simple matrix Algebra.  From Revell, 2009 Equation 3. *)
(*  This function reconstructs ancestral quantitative values for traits on a phylogenetic tree following Martins & Hansen,
    1997, Am. Nat..  It is the same as ReconstructNodes but does not return rate or standard errors. *)

(*
    Usage:  nodes = ReconstructNodesSimple[tree, trait]
            where tree is a string containing a tree with branch lengths in Newick format or a tree table of the sort produced 
            by TreeToTable[].  The function returns the ancestral 
            values of Y at each node in the tree assuming a Brownian motion.  

    	Created by P. David Polly, 5 March 2017

*)

ReconstructNodesSimple[tree_,trait_]:=Block[{J,GrandMean,Yprime,Nodes,Rate,x,contrasts, i,j,SDNode,X,SyGLS,varY,varAY,varA,TipLabels,NodeLabels,Y,VarY,VarAY,VarA,myTable,RateSquared,VarYNew},

If[StringQ[tree],myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];
];
];


{VarY,VarAY,VarA}=PhylogeneticMatrices[myTable];
varY=VarY[[2;;,2;;]];varAY=VarAY[[2;;,2;;]]; varA=VarA[[2;;,2;;]];
TipLabels=VarY[[2;;,1]];NodeLabels=VarA[[2;;,1]];

If[Length[trait[[1]]]!=2,Return["Trait must be a two-column matrix with tip labels in the first column and trait values in the second column"]];
If[Sort[TipLabels]!=Sort[trait[[1;;,1]]],Return["Traits labels do not match tip labels from Var Y."]];
Y=trait[[1;;,2]][[Flatten[Position[trait[[1;;,1]],#]&/@TipLabels]]];

J=Table[1,{Length[varY]}];
GrandMean=(J . Inverse[varY] . Y)/(J . Inverse[varY] . J);
Yprime=#-GrandMean&/@Y;
Nodes=#+GrandMean&/@(varAY . Inverse[varY] . Yprime);

Return[Transpose[{NodeLabels,Nodes}]];

];




(* This function draws a Newick format tree and labels the nodes and tips, making use of helpful suggestions given by Felsenstein (2004).


    Usage:  DrawNewickTree[tree (, aspectratio)]
            where tree tree is a string containing a tree with branch lengths in Newick format or a tree table of 
            the sort produced by TreeToTable[], and aspectratio is a number specifying the aspect ratio of the graphic output.

	The function requires that taxon names be unique.  It handles polytomies, but fails when the tree format is not specified correctly or when
	branch lengths are not given.  

    Created by P. David Polly, 5 April 2012
    22 December 2018:  Updated to reduce the size of the dot on Node 0.
    26 December 2018:  Updated by changing tip taxon text formatting for easier import into Illustrator.
    5 January 2018: Updated to make aspect ratio scaling responsive to number of taxa.
    29 May 2023: Updated to add user-controlled option for aspect ratio scaling.

*)

DrawNewickTree[tree_,aspect_:0]:=Block[{TreeTable,myTable,ExtNodeNames,ImmediateDescendants,FindLineage,TipEntries,NodeEntries,NodeNames,TipNames,TipY,TipLineages,NodeLineages,Descendants,NodeY,TreeX,TreeY,TreeXY,x},

If[!StringQ[tree]||!StringMatchQ[tree,"("~~___~~";"],Return["Input not in Newick format."]];
FindLineage[Point_]:=Block[{Lineage,NewBranch},
Lineage={Point};
While[!StringMatchQ[Lineage[[-1,2]],"Node 0"],
NewBranch=Flatten[TreeTable[[Flatten[Position[TreeTable[[1;;,1]],Lineage[[-1,2]]]]]]];
Lineage=Append[Lineage,NewBranch];
];
Return[Lineage];
];

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];
TreeTable=myTable;
Clear[myTable];


TipEntries=TreeTable[[Flatten[Position[StringCases[TreeTable[[1;;,1]],"Node "~~__],{}]]]];
NodeEntries=TreeTable[[Flatten[Position[TreeTable[[1;;,1]],#]&/@Flatten[StringCases[TreeTable[[1;;,1]],"Node "~~__]]]]];
NodeNames=NodeEntries[[1;;,1]];
TipNames=TipEntries[[1;;,1]];
TipY=Transpose[{TipNames,Table[x,{x,Length[TipNames],1,-1}]}];
TipLineages=FindLineage[#]&/@TipEntries;
NodeLineages=FindLineage[#]&/@NodeEntries;
Descendants=Table[Cases[If[MemberQ[#,{__,NodeNames[[x]],__}],#[[1,1]]]&/@TipLineages,Except[Null]],{x,Length[NodeNames]}];
NodeY=Transpose[{NodeNames,Table[Mean[TipY[[Flatten[Position[TipY[[1;;,1]],#]&/@Descendants[[x]]],2]]]//N,{x,Length[Descendants]}]}];
TreeX=Union[Sort[Partition[Flatten[Table[Transpose[{TipLineages[[x,1;;,1]],Reverse[Drop[FoldList[Plus,0,Reverse[TipLineages[[x,1;;,3]]]],1]]}],{x,Length[TipLineages]}]],2]]];
TreeY=Sort[Flatten[{TipY,NodeY},1]];
TreeXY=Transpose[{TreeX[[1;;,1]],TreeX[[1;;,2]],TreeY[[1;;,2]]}];
TreeXY=Sort[Append[TreeXY,{"Node 0",0,Mean[TreeY[[1;;,2]]]}]];

ExtNodeNames=Reverse[Sort[Append[NodeNames,"Node 0"]]];
ImmediateDescendants=TreeTable[[Flatten[Position[TreeTable[[1;;,2]],#]],1]]&/@ExtNodeNames;
Do[Do[
TreeXY[[Flatten[Position[TreeXY[[1;;,1]],ExtNodeNames[[x]]]],3]]=Mean[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#]&/@ImmediateDescendants[[x]]],3]]]//N,
{x,Length[ExtNodeNames]}],{2}];

If[aspect==0,
Return[
Graphics[
{
(* Horizontal lines *)
{Thick,Line[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{2,3}]]],Flatten[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2}]]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{3}]]]}]}]&/@TreeTable},
(* Vertical lines *)
{Thick,Line[{Flatten[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2}]]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{3}]]]}],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2,3}]]]}]&/@TreeTable},
(* tip labels *)
Table[{FontFamily->"Arial",FontSize->12,Text[StringReplace[TipNames[[x]],{"_"->" ","-"->" "}],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],TipNames[[x]]]],{2,3}]]],{-1.3,0}]},{x,Length[TipNames]}],
(* node labels *)
Table[{FontFamily->"Arial",FontSize->10,Text[NodeNames[[x]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],NodeNames[[x]]]],{2,3}]]],{-1.3,0}]},{x,Length[NodeNames]}],
(* Base node and label *)
{PointSize[0.01], Point[Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],"Node 0"]],{2,3}]]]]},
Text["Node 0",Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],"Node 0"]],{2,3}]]],{-1.3,0}]
}
,Axes->{True,False},AxesOrigin->{0,0*Length[TipNames]},AspectRatio->Log[Length[TipNames]]/GoldenRatio]

],

Return[
Graphics[
{
(* Horizontal lines *)
{Thick,Line[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{2,3}]]],Flatten[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2}]]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{3}]]]}]}]&/@TreeTable},
(* Vertical lines *)
{Thick,Line[{Flatten[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2}]]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{3}]]]}],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2,3}]]]}]&/@TreeTable},
(* tip labels *)
Table[{FontFamily->"Arial",FontSize->12,Text[StringReplace[TipNames[[x]],{"_"->" ","-"->" "}],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],TipNames[[x]]]],{2,3}]]],{-1.3,0}]},{x,Length[TipNames]}],
(* node labels *)
Table[{FontFamily->"Arial",FontSize->10,Text[NodeNames[[x]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],NodeNames[[x]]]],{2,3}]]],{-1.3,0}]},{x,Length[NodeNames]}],
(* Base node and label *)
{PointSize[0.01], Point[Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],"Node 0"]],{2,3}]]]]},
Text["Node 0",Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],"Node 0"]],{2,3}]]],{-1.3,0}]
}
,Axes->{True,False},AxesOrigin->{0,0*Length[TipNames]},AspectRatio->aspect]

]
]

];



(* This function takes a list of string labels and converts them to PHYLIP format by
   truncating labels greater than 10 characters or padding ones less than 10 characters.


    Usage:  phyliplabels = MakePHYLIPLabels[labels]
            where labels is a list of strings.

    Created by P. David Polly, 5 April 2012

*)

MakePHYLIPLabels[labels_]:=Block[{phyliplabels,x},
phyliplabels=labels;
Return[Table[If[StringLength[phyliplabels[[x]]]<= 10, phyliplabels[[1]]=phyliplabels[[x]]<>StringJoin[Table[" ",{10-StringLength[phyliplabels[[x]]]}]],StringDrop[phyliplabels[[x]], -1*(StringLength[phyliplabels[[x]]]-10)]],{x,Length[phyliplabels]}]];
]




(*  This function performs a brownian motion random walk for n generations with 
    a step rate of i.

	Usage:  walk = RandomWalk[100, 1];

	Created by P. David Polly, 13 April 2012

*)

RandomWalk[n_,i_]:=NestList[#+Random[NormalDistribution[0,i]]&,0,n];




(*  This function returns standardized Independent Contrasts (Felsenstein, 1985) for a
    continuous trait.   

	Usage:  contrasts = IndependentContrasts[tree, trait];

	where tree is a phylogentic tree with branch lengths in Newick format and trait is a 
	a two-column matrix with taxon labels in the first column (which must match the taxon
	names in tree) and trait values in the second column.

	Created by P. David Polly, Updated 19 May 2012

*)

IndependentContrasts[tree_,trait_]:=Block[{myTable, contrasts,myData0,x,stdcontrasts,FindDescendants,ProcessTip,TraverseNode},

FindDescendants[me1_]:=Block[{},
Return[Select[myTable,#[[2]]==me1&]]
];
ProcessTip[me2_]:=Block[{myLength,myID,myValue},
myLength=Select[myTable,#[[1]]==me2&][[1,3]];
myValue=Select[trait,#[[1]]==me2&][[1,2]];
Return[{myValue,myLength}];
];
TraverseNode[me3_]:=Block[{myDescendants,myLength,myValue,myData},
myLength=Select[myTable,#[[1]]==me3&][[1,3]];
myDescendants=FindDescendants[me3];
myData=If[#[[4]]==0,TraverseNode[#[[1]]],ProcessTip[#[[1]]]]&/@myDescendants;
Do[
contrasts=Append[contrasts,{myData[[x,1]]-myData[[x+1,1]],myData[[x,2]]+myData[[x+1,2]]}];
,{x,Length[myData]-1}];
myValue=(Plus@@((1/#[[2]])*#[[1]]&/@myData))/(Plus@@(1/#[[2]]&/@myData));
myLength=myLength+(Fold[Times, 1, myData[[1;;,2]]]/(Plus@@myData[[1;;,2]]));
(* Return[{myValue,myLength}]; *)
Return[{myValue,myLength}];
];

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Return["Tree variable must be a Newick tree string or a tree table."];];];
];

If[Length[Dimensions[trait]]!= 2, Return["Trait variable must have two columns"]];
If[Sort[trait[[1;;,1]]]!= Sort[Union[Select[myTable,#[[4]]==1&][[1;;,1]]]], Return["Trait labels must equal tree names"]];

contrasts={};
myData0=If[#[[4]]==0,TraverseNode[#[[1]]],ProcessTip[#[[1]]]]&/@FindDescendants["Node 0"];
Do[
contrasts=Append[contrasts,{myData0[[x,1]]-myData0[[x+1,1]],myData0[[x,2]]+myData0[[x+1,2]]}//N];
,{x,Length[myData0]-1}];

stdcontrasts=#[[1]]/Sqrt[#[[2]]]&/@contrasts;
Return[stdcontrasts];
];




(*  This function converts a tree in tabular form to a Newick format string.  .  

	Usage:  tree = TableToTree[treetable];

	where treetable is a table in the same format as produced by the 
    TreeToTable[] function.

	Created by P. David Polly, 3 May 2012

*)

TableToTree[myTable_]:=Block[{contrasts,myData0,x,stdcontrasts,newick,myValue0,myLength0,myEntry0,FindDescendants,ProcessTip,TraverseNode},

FindDescendants[me1_]:=Block[{},
Return[Select[myTable,#[[2]]==me1&]]
];

ProcessTip[me2_]:=Block[{myLength,myEntry},
myLength=Select[myTable,#[[1]]==me2&][[1,3]];
myEntry=me2<>":"<>ToString[myLength];
Return[myEntry];
];

TraverseNode[me3_]:=Block[{myDescendants,myLength,myData,myEntry},
myLength=Select[myTable,#[[1]]==me3&][[1,3]];
myDescendants=FindDescendants[me3];
myData=If[#[[4]]==0,TraverseNode[#[[1]]],ProcessTip[#[[1]]]]&/@myDescendants;
	   myEntry="("<>StringJoin[Riffle[myData,","]]<>"):"<>ToString[myLength];

Return[myEntry];
];



myData0=If[#[[4]]==0,TraverseNode[#[[1]]],ProcessTip[#[[1]]]]&/@FindDescendants["Node 0"];
   myEntry0="("<>StringJoin[Riffle[myData0,","]]<>");";

Return[myEntry0];
];



(*  This function reroots a Newick tree at the specified node. 

	Usage:  newtree = RerootTree[tree, node];

	where tree is a string containing a tree with branch lengths in Newick format 
    or a tree table of the sort produced by TreeToTable[], and node 
    is the name of the node where the new tree will be rooted.  The name of the
    node can be found by plotting the tree with DrawNewickTree[tree].  

	Created by P. David Polly, 10 May 2012

*)

RerootTree[tree_,root_]:=Block[{myData0,x,myLength0,myEntry0,TablePosition,myLine,myTable,FindLineage, FindDescendants,ProcessTip,TraverseNode},

FindLineage[Point_]:=Block[{Lineage,NewBranch},
Lineage={Point};
While[!StringMatchQ[Lineage[[-1,2]],"Node 0"],
TablePosition=Append[TablePosition,Flatten[Position[myTable[[1;;,1]],Lineage[[-1,2]]]]];
NewBranch=Flatten[myTable[[TablePosition[[-1]]]]];
Lineage=Append[Lineage,NewBranch];
];
Return[{Lineage,Flatten[TablePosition]}];
];

FindDescendants[me1_]:=Block[{},
Return[Select[myTable,#[[2]]==me1&]]
];

ProcessTip[me2_]:=Block[{myLength,myEntry},
myLength=Select[myTable,#[[1]]==me2&][[1,3]];
myEntry=me2<>":"<>ToString[myLength];
Return[myEntry];
];

TraverseNode[me3_]:=Block[{myDescendants,myLength,myData,myEntry},
myLength=Select[myTable,#[[1]]==me3&][[1,3]];
myDescendants=FindDescendants[me3];
myData=If[#[[4]]==0,TraverseNode[#[[1]]],ProcessTip[#[[1]]]]&/@myDescendants;
	   myEntry="("<>StringJoin[Riffle[myData,","]]<>"):"<>ToString[myLength];

Return[myEntry];
];

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];

TablePosition={};


myLine=FindLineage[myTable[[Flatten[Position[myTable[[1;;,2]],root]][[1]]    ]] ];

Do[
myTable[[TablePosition[[x]],{1,2}]]=myTable[[TablePosition[[x]],{2,1}]],
{x,Length[TablePosition]}];

myData0=If[#[[4]]==0,TraverseNode[#[[1]]],ProcessTip[#[[1]]]]&/@FindDescendants[root];
 myEntry0="("<>StringJoin[Riffle[myData0,","]]<>");";

Return[myEntry0];
];


(*  This function is deprecated, having been replaced by SimulateContinuousTraitsOnTree[].  It is
    retained for backward compatibility *)

(* This function simulates the evolution of a continuous trait on a tree using
   a Brownian motion model of evolution.  By default the trait starts with a 
   value of 0.0 at the root of the tree (this can be changed by specifying a 
   different value as the third input parameter of the function).  Evolution is
   simulated along each branch by Brownian motion using the rate specified by 
   the second input parameter.  The function returns a list of node and tip labels
   followed by the simulated trait value. 

   Usage:   SimulateContinuousTraitOnTree[tree, rate (,startvalue, includenodes)]

   where tree is a string containing a tree with branch lengths in Newick format or a tree table of 
   the sort produced by TreeToTable[], rate is the amount of trait change
   per unit of branch length (as specified in the tree), startvalue is the optional starting
   value of the trait at the root of the tree, and includenodes is the optional parameter
   that specifies whether simulated node values should be returned (1) or not returned (0, default).

   Created by P. David Polly, 20 April 2012

*)

SimulateContinuousTraitOnTree[tree_,Rate_,StartValue_:0,IncludeNodes_:0]:=Block[{myTable,x,myData0,TraitValues,FindDescendants,ProcessTip,TraverseNode},

FindDescendants[me_]:=Block[{},
Return[Select[myTable,#[[2]]==me&]]
];

ProcessTip[me1_,myAnc1_,myLength1_,myStart1_]:=Block[{myValue},
myValue=RandomReal[NormalDistribution[myStart1,Sqrt[myLength1]*Rate]];
TraitValues=Append[TraitValues,{me1,myValue}];
Return[];
];

TraverseNode[me2_,myAnc2_,myLength2_,myStart2_]:=Block[{myDescendants,myValue},
myValue=RandomReal[NormalDistribution[myStart2,Sqrt[myLength2]*Rate]];
TraitValues=Append[TraitValues,{me2,myValue}];
myDescendants=FindDescendants[me2];
If[#[[4]]==0,TraverseNode[#[[1]],#[[2]],#[[3]],myValue],ProcessTip[#[[1]],#[[2]],#[[3]],myValue]]&/@myDescendants;
Return[];
];

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];

TraitValues={{"Node 0",StartValue}};
If[#[[4]]==0,TraverseNode[#[[1]],#[[2]],#[[3]],StartValue],ProcessTip[#[[1]],#[[2]],#[[3]],StartValue]]&/@FindDescendants["Node 0"];

If[IncludeNodes==0,TraitValues=TraitValues[[Flatten[Position[TraitValues[[1;;,1]],#]&/@Complement[TraitValues[[1;;,1]],Union[Flatten[StringCases[TraitValues[[1;;,1]],"Node "~~DigitCharacter..]]]]]]]];

Return[TraitValues];
];





(* This function simulates the evolution of a continuous trait on a tree using
   a Brownian motion model of evolution.  By default the trait starts with a 
   value of 0.0 at the root of the tree (this can be changed by specifying a 
   different value as the third input parameter of the function).  Evolution is
   simulated along each branch by Brownian motion using the rate specified by 
   the second input parameter.  The function returns a list of node and tip labels
   followed by the simulated trait value. 

   Usage:   SimulateContinuousTraitsOnTree[tree, rate (,startvalue, includenodes)]

   where tree is a string containing a tree with branch lengths in Newick format or a tree table of 
   the sort produced by TreeToTable[], rate is the amount of trait change
   per unit of branch length (as specified in the tree), stepsperunit is an optional scaling factor for 
   the number of steps in the random walk process per unit of
   tree length (by default this is 1.0), startvalue is the optional starting
   value of the trait at the root of the tree, includenodes is the optional parameter
   that specifies whether simulated node values should be returned (1) or not returned (0, default).

   Created by P. David Polly, 20 April 2012
   Updated, 17 October 2014:  removed number of steps parameter.

*)

SimulateContinuousTraitsOnTree[tree_,Rate_,StartValue_:0,IncludeNodes_:0]:=Block[{myTable,x,myData0,TraitValues,FindDescendants,ProcessTip,TraverseNode,Rates,StartValues,i},

If[Length[Rate]==0,Rates=Flatten[{Rate}], Rates=Rate];
If[Length[StartValue]==0,StartValues=Table[StartValue,{Length[Rates]}],StartValues=StartValue];

FindDescendants[me_]:=Block[{},
Return[Select[myTable,#[[2]]==me&]]
];

ProcessTip[me1_,myAnc1_,myLength1_,myStart1_]:=Block[{myValue},
myValue=Table[Random[NormalDistribution[myStart1[[i]],Sqrt[myLength1]*Rates[[i]]]],{i,Length[Rates]}];
TraitValues=Append[TraitValues,{me1,myValue}];
Return[];
];

TraverseNode[me2_,myAnc2_,myLength2_,myStart2_]:=Block[{myDescendants,myValue},
myValue=Table[Random[NormalDistribution[myStart2[[i]],Sqrt[myLength2]*Rates[[i]]]],{i,Length[Rates]}];
TraitValues=Append[TraitValues,{me2,myValue}];
myDescendants=FindDescendants[me2];
If[#[[4]]==0,TraverseNode[#[[1]],#[[2]],#[[3]],myValue],ProcessTip[#[[1]],#[[2]],#[[3]],myValue]]&/@myDescendants;
Return[];
];

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];

TraitValues={{"Node 0",StartValues}};
If[#[[4]]==0,TraverseNode[#[[1]],#[[2]],#[[3]],StartValues],ProcessTip[#[[1]],#[[2]],#[[3]],StartValues]]&/@FindDescendants["Node 0"];

If[IncludeNodes==0,TraitValues=TraitValues[[Flatten[Position[TraitValues[[1;;,1]],#]&/@Complement[TraitValues[[1;;,1]],Union[Flatten[StringCases[TraitValues[[1;;,1]],"Node "~~DigitCharacter..]]]]]]]];

Return[Flatten[#]&/@TraitValues];
];





(* This function simulates the evolution of correlated continuous traits on a tree using
   a Brownian motion model of evolution.  By default the traits each start with a 
   value of 0.0 at the root of the tree (this can be changed by specifying a 
   vector of values as the third input parameter of the function).  Evolution is
   simulated along each branch by Brownian motion using the rate specified by 
   the second input parameter with covariance between the traits as specified in the matrix.  
 The number of traits is determined from the structure of the covariance matrix.  The rates
  of evolution are derived from the trace of the covariance matrix, with the rate at each 
  step of the simulation being equivalent to a random step whose variance is equal to the 
  variance of that trait in the covariance matrix.  
 The function returns a list of node and tip labels followed by the simulated trait values. 

   Usage:   SimulateCorrelatedTraitsOnTree[tree, covariancematrix (,startvector, includenodes)]

   where tree is a string containing a tree with branch lengths in Newick format or a tree table of 
   the sort produced by TreeToTable[], covariancematrix is the covariance matrix of the traits.
   per unit of branch length (as specified in the tree), startvector is the optional starting
   values of the traits at the root of the tree, and includenodes is the optional parameter
   that specifies whether simulated node values should be returned (1) or not returned (0, default).

   Created by P. David Polly, 23 July 2012
   Updated Januery 2014

*)

SimulateCorrelatedTraitsOnTree[tree_,P_,StartVector_:0,IncludeNodes_:0]:=Block[{myTable,x,myData0,TraitValues,FindDescendants,ProcessTip,TraverseNode,Chol},

FindDescendants[me_]:=Block[{},
Return[Select[myTable,#[[2]]==me&]]
];

ProcessTip[me1_,myAnc1_,myLength1_,myStart1_]:=Block[{myValue},
myValue=Table[Random[NormalDistribution[0,1]],{Length[Chol]}] . (Sqrt[myLength1]*Chol)+myStart1;
TraitValues=Append[TraitValues,{me1,myValue}];
Return[];
];

TraverseNode[me2_,myAnc2_,myLength2_,myStart2_]:=Block[{myDescendants,myValue},
myValue=Table[Random[NormalDistribution[0,1]],{Length[Chol]}] . (Sqrt[myLength2]*Chol)+myStart2;
TraitValues=Append[TraitValues,{me2,myValue}];
myDescendants=FindDescendants[me2];
If[#[[4]]==0,TraverseNode[#[[1]],#[[2]],#[[3]],myValue],ProcessTip[#[[1]],#[[2]],#[[3]],myValue]]&/@myDescendants;
Return[];
];

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];

Chol=CholeskyDecomposition[P];
TraitValues={{"Node 0",If[StartVector==0,Table[0,{Length[Chol]}],StartVector]}};
If[#[[4]]==0,TraverseNode[#[[1]],#[[2]],#[[3]],StartVector],ProcessTip[#[[1]],#[[2]],#[[3]],StartVector]]&/@FindDescendants["Node 0"];

If[IncludeNodes==0,TraitValues=TraitValues[[Flatten[Position[TraitValues[[1;;,1]],#]&/@Complement[TraitValues[[1;;,1]],Union[Flatten[StringCases[TraitValues[[1;;,1]],"Node "~~DigitCharacter..]]]]]]]];

Return[TraitValues];
];


(* 

	This function reads in a Phylip formatted cladistic character matrix and returns it as an array with the 
	taxon name in the first column and the characters separated into columns.
	
	Usage:    ReadPhylip[filename]

	where filename is the name of the phylip text file, with full path if necessary.

	Created by P. David Polly, 29 November 2012.

*)

ReadPhylip[fl_]:=Block[{data,numtaxa,numchars,Taxa,Chars,x},
data=OpenRead[fl];
numtaxa=ToExpression[Read[data,Word]];
numchars=ToExpression[Read[data,Word]];
Taxa={};
Chars={};
Do[
Taxa=Append[Taxa,Read[data,Word]]; 
Chars=Append[Chars,Read[data,Word]],
{x,numtaxa}];
Close[data];
Chars=Table[If[#=="?",#,ToExpression[#]]&/@StringSplit[Chars[[x]],""],{x,numtaxa}];
Return[Table[Prepend[Chars[[x]],Taxa[[x]]],{x,numtaxa}]];
]


(*  This function performs a Monte Carlo simulation of an evolving trait mean
    for n generations with a random component of r, a directional component 
    of m, and a starting value of s.  When OU constraints are modeled, the target 
    or optimal value is t and the strength of selection is a (0 to 1).  

	Usage:  mylineage = LineageEvolution[100, 1];

	Created by P. David Polly,  24 September 2014

*)

LineageEvolution[n_,OptionsPattern[{"r"->1,"d"->0,"s"->0,"a"->0,"t"->0}]]:=Block[{x},
If[OptionValue["a"]<0||OptionValue["a"]>1, Return["\"a\" must be between 0 and 1."]];
Return[Thread[{Table[x,{x,0,n}],NestList[#+Random[NormalDistribution[OptionValue["a"](OptionValue["t"]-#)+OptionValue["d"],OptionValue["r"]]]&,OptionValue["s"],n]}]];
];




(*  This function multiplies a tree matrix by a scaling factor Pagel's \[Lambda], which usually has a value between 0 and 1.  

	Usage:  NewC = LambdaMultiplication[C,lambda]
		where C is a phylogenetic tree covariance matrix and lambda is a value between 0 and 1 (can be higher up to the maximum value of Pagel's \[Lambda] for the tree)

	Created by P. David Polly, 26 November 2022
      
*)

LambdaMultiplication[cmat_,lambda_]:=Block[{TemplateMatrix},
TemplateMatrix=Table[1,{Length[cmat]},{Length[cmat]}]-IdentityMatrix[Length[cmat]];
Return[cmat*TemplateMatrix*lambda+IdentityMatrix[Length[cmat]]*cmat];
];



(*  This function estimates Pagel's Lambda using the Likelihood method described by Freckleton et 
	al. 2002 and partially following Revell's phylosig() function in R.

	Usage:  mylineage = PagelsLambda[tree, trait];

	where tree is a tree in Newick format or a tree table and trait is a vector trait values with tip labels in the first column.
	
	Created by P. David Polly,  24 September 2014
	Rewritten 26 November 2022
	Tweaked 8 December 2022
	Substituted PseudoInverse for Inverse to deal with singular adjusted phylogenetic matrices.

*)

PagelsLambda[mytree_,trait_]:=Block[{LambdaMultiplication,treetable,treelabels,treeorder,treetrait,y,Cmat,minbranch,treeheight,MaxLambda,One,n,LogLikelihoodLambda,loglike,lambda,LCmat,InvLCmat,aL,residsL,sigmaL,alllikes,maxlikes,warn},

LambdaMultiplication[cmat_,lambda_]:=Block[{TemplateMatrix},
TemplateMatrix=Table[1,{Length[cmat]},{Length[cmat]}]-IdentityMatrix[Length[cmat]];
Return[cmat*TemplateMatrix*lambda+IdentityMatrix[Length[cmat]]*cmat];
];

treetable=TreeToTable[mytree];
treelabels=PhylogeneticMatrices[mytree][[1,1,2;;]];
treeorder=Flatten[Position[trait[[1;;,1]],#]&/@treelabels];
treetrait=trait[[treeorder]];
y=treetrait[[1;;,2;;]];
Cmat=PhylogeneticMatrices[mytree][[1,2;;,2;;]];
minbranch=Min[Select[treetable,#[[-1]]==1&][[1;;,3]]];
treeheight=Max[Tr[Cmat,List]];
MaxLambda=treeheight/(treeheight-minbranch);
(* MaxLambda=1; *)
One=Table[1,{Length[Cmat]}];
n=Length[y];
warn="";

LogLikelihoodLambda[lam_,Cmat_,y_]:=Block[{LambdaCmat,a,sigma,InvLambdaCmat,resids,LogL},
LambdaCmat=LambdaMultiplication[Cmat,lam];
InvLambdaCmat=PseudoInverse[LambdaCmat];
a=((One . InvLambdaCmat . y)/(One . InvLambdaCmat . One))[[1]];
resids=#-a&/@y;
sigma=Flatten[((Transpose[resids] . InvLambdaCmat . resids)/n)][[1]];
LogL=Flatten[-1* (Transpose[resids] . (InvLambdaCmat/sigma) . resids)/2-n*Log[2 Pi]/2 -Det[sigma*LambdaCmat]/2][[1]];
Return[LogL]
];

 (*
{loglike,lambda}=Quiet[NMaximize[{LogLikelihoodLambda[lam,Cmat,y],0<=lam<MaxLambda},lam,WorkingPrecision->5]];
lambda=Chop[lambda[[1,2]]];
*)

alllikes=Table[{x,LogLikelihoodLambda[x,Cmat,y]},{x,0,MaxLambda,0.005}];
maxlikes=alllikes[[Flatten[Position[alllikes[[1;;,2]],Max[alllikes[[1;;,2]]]]]]];
If[Length[maxlikes]>1,{lambda,loglike}=Flatten[RandomSample[maxlikes,1]]; Print["Warning: Likelihood surface is flat over part of the parameter space, the value of \[Lambda] is ambiguous."],{lambda,loglike}=Flatten[maxlikes]];


LCmat=LambdaMultiplication[Cmat,lambda];
InvLCmat=PseudoInverse[LCmat];
aL=((One . InvLCmat . y)/(One . InvLCmat . One))[[1]];
residsL=#-aL&/@y;
sigmaL=Flatten[((Transpose[residsL] . InvLCmat . residsL)/n)][[1]];

Return[{"\[Lambda]"->Round[lambda,0.01],"LogLikelihood"->loglike,"Ancestral Value"->aL,"\[Sigma]2"->sigmaL}]
]




(*  This function fits three evolutionary models to data derived from an evolutionary lineage 
	following Hunt (2006). Using Hunt's terminology, it fits an unbiased random walk, a general
	random walk, and a stasis model to a set of trait data from a lineage sampled at different
	points along a lineage.  

	Usage:  ThreeModelTest[lineage];

	where lineage is a data matrix with time in the first column, trait values in the second column,
	sample error as a variance in the third column, and number of individuals in the sample in the
	fourth column.  If only the first two columns are given then the algorithm assumes there is 
	no sampling error.
	
	Created by P. David Polly,  21 January 2018
	Updated 7 July 2018 to remove minor inconsistencies with Hunt's function.

*)
BestEstimates[lineage_]:=Block[{MaxRand,MaxStasis,mu,sigma,theta,omega,i,j},
MaxRand=Maximize[{LogLikeRandWalk[i,j,lineage],0<j},{i,j}];
MaxStasis=Maximize[{LogLikeStasis[i,j,lineage],0<j},{i,j}];
mu=MaxRand[[2,1,2]];
sigma=MaxRand[[2,2,2]];
theta=MaxStasis[[2,1,2]];
omega=MaxStasis[[2,2,2]];
Return[{mu,sigma,theta,omega}];
]

LogLikeRandWalk[mu_,sigma_,lineage_]:=Block[{l, t,trait,SE,n,df,sv,svA,svD,svAD,Var,Mn},
t=Abs[Differences[lineage[[1;;,1]]]];
trait=Differences[lineage[[1;;,2]]];
df=Length[trait];
If[Length[lineage[[1]]]>2,SE=lineage[[1;;,3]],SE=Table[0,{df+1}]];
If[Length[lineage[[1]]]>3,n=lineage[[1;;,4]],n=Table[1,{df+1}]];
sv=SE/n;
svA=sv[[1;;df]];
svD=sv[[2;;]];
svAD=svA+svD;
l=Table[Mn=mu*t[[x]];Var=sigma*t[[x]]+svAD[[x]];-0.5Log[2Pi]-0.5*Log[Var]-((1/(2Var))*((trait[[x]]-Mn)^2)),{x,df}];
Return[Plus@@l];
]

LogLikeStasis[theta_,omega_,lineage_]:=Block[{trait,SE,n,l,df,sv,svD,anc,Mn,Var},
trait=Differences[lineage[[1;;,2]]];
df=Length[trait];
If[Length[lineage[[1]]]>2,SE=lineage[[1;;,3]],SE=Table[0,{df+1}]];
If[Length[lineage[[1]]]>3,n=lineage[[1;;,4]],n=Table[1,{df+1}]];
sv=SE/n;
svD=sv[[2;;]];
anc=lineage[[1;;df,2]];
l=Table[Mn=theta-anc[[x]];Var=omega+svD[[x]];-0.5Log[2Pi]-0.5*Log[Var]-((1/(2Var))*((trait[[x]]-Mn)^2)),{x,df}];
Return[Plus@@l];
]

AIC[l_,k_]:=2k-2l
AICc[l_,k_,n_]:=2k-2l+(2*k*(k+1))/(n-k-1) 

AICweights[AICs_]:=Block[{minAIC,deltaAIC,sumExp,w},
minAIC=Min[AICs];
deltaAIC=(#-minAIC)&/@AICs;
sumExp=Plus@@(Exp[-0.5*#]&/@deltaAIC);
w=(Exp[-0.5*#]/sumExp)&/@deltaAIC;
Return[w];
]
ThreeModelTest[lineage_]:=Block[{Xbar,Xsigma,Xtheta,Xomega,lGRW,lURW,lStasis,AICGRW,AICURW,AICStasis,w,myTable,bestAIC,bestModel},
{Xbar,Xsigma,Xtheta,Xomega}=BestEstimates[lineage];
lURW=LogLikeRandWalk[0,Xsigma,lineage];
lGRW=LogLikeRandWalk[Xbar,Xsigma,lineage];
lStasis=LogLikeStasis[Xtheta,Xomega,lineage];
AICURW=AICc[lURW,1,Length[lineage]-1];
AICGRW=AICc[lGRW,2,Length[lineage]-1];
AICStasis=AICc[lStasis,2,Length[lineage]-1];
w=AICweights[{AICURW,AICGRW,AICStasis}];
bestAIC=Min[{AICURW,AICGRW,AICStasis}];
bestModel=Which[AICGRW==bestAIC,"GRW",AICURW==bestAIC,"URW",AICStasis==bestAIC,"Stasis"];
myTable={
{"Best model",bestModel},
{"N",Length[lineage]},
{"\!\(\*SubscriptBox[OverscriptBox[\(\[Mu]\), \(^\)], \(step\)]\)",NumberForm[Xbar,{5,3}]},
{"\!\(\*SubscriptBox[SuperscriptBox[OverscriptBox[\(\[Sigma]\), \(^\)], \(2\)], \(step\)]\)",NumberForm[Xsigma,{5,3}]},
{"\!\(\*OverscriptBox[\(\[Theta]\), \(^\)]\)",NumberForm[Xtheta,{5,3}]},
{"\!\(\*OverscriptBox[\(\[Omega]\), \(^\)]\)",NumberForm[Xomega,{5,3}]},
{"\!\(\*SubscriptBox[\(AIC\), \(c\)]\) URW",NumberForm[Chop[AICURW],{5,3}]},
{"\!\(\*SubscriptBox[\(AIC\), \(c\)]\) GRW",NumberForm[Chop[AICGRW],{5,3}]},
{"\!\(\*SubscriptBox[\(AIC\), \(c\)]\) Stasis",NumberForm[Chop[AICStasis],{5,3}]},
{"wt URW",NumberForm[Chop[w[[1]]],{5,3}]},
{"wt GRW",NumberForm[Chop[w[[2]]],{5,3}]},
{"wt Stasis",NumberForm[Chop[w[[3]]],{5,3}]}
};
Return[TableForm[myTable[[1;;,2]],TableHeadings->{myTable[[1;;,1]],{"A","B"}}]];
]





(*	This function estimates the ancestor value for the base node in a tree
	  assuming Brownian motion based on the equations presented by O'Meara 
	  et al., Revell and Harmon, Adams, and others.  
	  
	  Usage:  OMearaAncestor[MyTree, MyTrait]
	  
	  where MyTree is a Newick-format tree in the format of ReadNewick[], MyTrait
	  is a matrix with a vector of taxon names in the first column and trait values in the second.
	  
	  Created by P. David Polly, 21 December 2018.
	  *)
OMearaAncestor[tree_,trait_]:=Block[{One,Cmatrix,a,resids,rate,treelabels,treedata,treeorder,data},
data=trait[[1;;,2;;]];
treelabels=PhylogeneticMatrices[tree][[1,1,2;;]];
treeorder=Flatten[Position[trait[[1;;,1]],#]&/@treelabels];
treedata=data[[treeorder]];
Cmatrix=PhylogeneticMatrices[tree][[1,2;;,2;;]];
One=Table[1,{Length[Cmatrix]}];
a=(One . Inverse[Cmatrix] . treedata)/(One . Inverse[Cmatrix] . One);
Return[a//N];
]






(*	This function estimates maps trait values onto a phylogenetic tree assuming a
	  Brownian motion model of evolution and using a color spectrum of the user's
	  choice.
	  
	  
	  Usage:  MapTraitsOntoTree[MyTree, MyTrait ("ColorPattern"-> MyColors)]
	  
	  where MyTree is a Newick-format tree in the format of ReadNewick[], MyTrait
	  is a vector of trait values, and DataLabels is a vector of taxon names for
	  MyTrait that match the taxon names in MyTree, and MyColors is an optional
	  list of three colors that will be used for minimum, base node, and maximum values of the 
	  trait.
	  
	  Created by P. David Polly, 21 December 2018
	  Updated 26 December 2018 to change text formating for easier import into Illustrator.
	  Updated 7 January 2019 to better differentiate 
	  *)
	  
	  
ColorLine[anc_,desc_,line_,rootvalue_,maxtrait_,mintrait_,colorpattern_]:=Block[{sortedline,newline,ancsign,descsign,ancresid,descresid,maxresid,minresid,line1,line2,prop},
ancresid=anc-rootvalue;
descresid=desc-rootvalue;
maxresid=maxtrait-rootvalue;
minresid=mintrait-rootvalue;
ancsign=Sign[ancresid];
descsign=Sign[descresid];
sortedline=Sort[line,Less];
Which[
(ancsign==1||ancsign==0)&&(descsign==1||descsign==0), newline=Line[sortedline,VertexColors->{Blend[colorpattern[[{2,3}]],ancresid/maxresid],Blend[colorpattern[[{2,3}]],descresid/maxresid]}],
(ancsign==-1||ancsign==0)&&(descsign==-1||descsign==0), newline=Line[sortedline,VertexColors->{Blend[colorpattern[[{2,1}]],ancresid/minresid],Blend[colorpattern[[{2,1}]],descresid/minresid]}],
(ancsign==-1)&&(descsign==1), 
prop=Abs[ancresid]/(Abs[ancresid]+Abs[descresid]);
line1=sortedline;
line1[[2,1]]=prop*(sortedline[[2,1]]-sortedline[[1,1]])+sortedline[[1,1]];
line2=sortedline;
line2[[1,1]]=prop*(sortedline[[2,1]]-sortedline[[1,1]])+sortedline[[1,1]];
newline={Line[line1,VertexColors->{Blend[colorpattern[[{2,1}]],ancresid/minresid],Blend[colorpattern[[{2,1}]],0]}],
Line[line2,VertexColors->{Blend[colorpattern[[{2,3}]],0],Blend[colorpattern[[{2,3}]],descresid/maxresid]}]},
(ancsign==1)&&(descsign==-1), 
prop=Abs[ancresid]/(Abs[ancresid]+Abs[descresid]);
line1=sortedline;
line1[[2,1]]=prop*(sortedline[[2,1]]-sortedline[[1,1]])+sortedline[[1,1]];
line2=sortedline;
line2[[1,1]]=prop*(sortedline[[2,1]]-sortedline[[1,1]])+sortedline[[1,1]];
newline={Line[line1,VertexColors->{Blend[colorpattern[[{2,3}]],ancresid/maxresid],Blend[colorpattern[[{2,3}]],0]}],
Line[line2,VertexColors->{Blend[colorpattern[[{2,1}]],0],Blend[colorpattern[[{2,1}]],descresid/minresid]}]}
];
Return[newline]
];
	  
MapTraitsOntoTree[tree_,traits_,OptionsPattern[{"ColorPattern"->{RGBColor[0.229801, 0.528862, 0.712718],RGBColor[0.992035, 0.876677, 0.310185],RGBColor[0.705653, 0.275837, 0.114107]},"TraitRange"->{"Default"},"NodeLabels"->False}]]:=Block[{LegendFactor,ancvalue,TraitValues,TreeTable,myTable,ExtNodeNames,ImmediateDescendants,FindLineage,TipEntries,NodeEntries,NodeNames,TipNames,TipY,TipLineages,NodeLineages,Descendants,NodeY,TreeX,TreeY,TreeXY,x,MaxTrait,MinTrait},
FindLineage[Point_]:=Block[{Lineage,NewBranch},
Lineage={Point};
While[!StringMatchQ[Lineage[[-1,2]],"Node 0"],
NewBranch=Flatten[TreeTable[[Flatten[Position[TreeTable[[1;;,1]],Lineage[[-1,2]]]]]]];
Lineage=Append[Lineage,NewBranch];
];
Return[Lineage];
];

If[StringQ[tree],
myTable=TreeToTable[tree][[2;;]],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
myTable=tree[[2;;]],
If[Length[tree[[1]]]==4,myTable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];
TreeTable=myTable;
Clear[myTable];

TipEntries=TreeTable[[Flatten[Position[StringCases[TreeTable[[1;;,1]],"Node "~~__],{}]]]];
NodeEntries=TreeTable[[Flatten[Position[TreeTable[[1;;,1]],#]&/@Flatten[StringCases[TreeTable[[1;;,1]],"Node "~~__]]]]];
NodeNames=NodeEntries[[1;;,1]];
TipNames=TipEntries[[1;;,1]];
TipY=Transpose[{TipNames,Table[x,{x,Length[TipNames],1,-1}]}];
TipLineages=FindLineage[#]&/@TipEntries;
NodeLineages=FindLineage[#]&/@NodeEntries;
Descendants=Table[Cases[If[MemberQ[#,{__,NodeNames[[x]],__}],#[[1,1]]]&/@TipLineages,Except[Null]],{x,Length[NodeNames]}];
NodeY=Transpose[{NodeNames,Table[Mean[TipY[[Flatten[Position[TipY[[1;;,1]],#]&/@Descendants[[x]]],2]]]//N,{x,Length[Descendants]}]}];
TreeX=Union[Sort[Partition[Flatten[Table[Transpose[{TipLineages[[x,1;;,1]],Reverse[Drop[FoldList[Plus,0,Reverse[TipLineages[[x,1;;,3]]]],1]]}],{x,Length[TipLineages]}]],2]]];
TreeY=Sort[Flatten[{TipY,NodeY},1]];
TreeXY=Transpose[{TreeX[[1;;,1]],TreeX[[1;;,2]],TreeY[[1;;,2]]}];
TreeXY=Sort[Append[TreeXY,{"Node 0",0,Mean[TreeY[[1;;,2]]]}]];

ExtNodeNames=Reverse[Sort[Append[NodeNames,"Node 0"]]];
ImmediateDescendants=TreeTable[[Flatten[Position[TreeTable[[1;;,2]],#]],1]]&/@ExtNodeNames;
Do[Do[
TreeXY[[Flatten[Position[TreeXY[[1;;,1]],ExtNodeNames[[x]]]],3]]=Mean[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#]&/@ImmediateDescendants[[x]]],3]]]//N,
{x,Length[ExtNodeNames]}],{2}];

TraitValues=Flatten[{#[[1]]->#[[2]]&/@traits,#[[1]]->#[[2]]&/@ReconstructNodesSimple[tree,traits]}];
ancvalue=OMearaAncestor[tree,traits][[1]];
If[OptionValue["TraitRange"][[1]]=="Default",
MaxTrait=Max[TraitValues[[1;;,2]]];
MinTrait=Min[TraitValues[[1;;,2]]],
MaxTrait=OptionValue["TraitRange"][[2]];
MinTrait=OptionValue["TraitRange"][[1]]
];
(* TraitColors[[1;;,2]]=#-MinTrait&/@TraitColors[[1;;,2]]/(MaxTrait-MinTrait);*)

Return[
Graphics[
{
(* Vertical lines *)
{AbsoluteThickness[3.5],CapForm["Square"],ColorLine[#[[2]]/.TraitValues,#[[2]]/.TraitValues,{Flatten[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2}]]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{3}]]]}],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2,3}]]]},ancvalue,MaxTrait,MinTrait,OptionValue["ColorPattern"]]&/@TreeTable},
(* Horizontal lines *)
{AbsoluteThickness[3.5],CapForm["Butt"],ColorLine[#[[1]]/.TraitValues,#[[2]]/.TraitValues,{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{2,3}]]],Flatten[{Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[2]]]],{2}]]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],#[[1]]]],{3}]]]}]},ancvalue,MaxTrait,MinTrait,OptionValue["ColorPattern"]]&/@TreeTable},
(* tip labels *)
Table[{Text[StringReplace[TipNames[[x]],{"_"->" ", "-"->" "}],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],TipNames[[x]]]],{2,3}]]],{-1.3,0}]},{x,Length[TipNames]}],

(* node labels *)
If[OptionValue["NodeLabels"]==True,
Table[{Gray,FontFamily->"Arial",FontSize->8,Text[NodeNames[[x]],Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],NodeNames[[x]]]],{2,3}]]],{-1.3,0}]},{x,Length[NodeNames]}],
(* Base node and label *)
{PointSize[0.01], Point[Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],"Node 0"]],{2,3}]]]]},
Text["Node 0",Flatten[TreeXY[[Flatten[Position[TreeXY[[1;;,1]],"Node 0"]],{2,3}]]],{-1.3,0}]],

(* Legend *)
If[Max[TreeXY[[1;;,3]]]<10,LegendFactor=1,LegendFactor=2];
{{CapForm["Butt"],AbsoluteThickness[6],ColorLine[MinTrait,MaxTrait,{{0,Max[TreeXY[[1;;,3]]]+LegendFactor},{Max[TreeXY[[1;;,2]]]/2,Max[TreeXY[[1;;,3]]]+LegendFactor}},ancvalue,MaxTrait,MinTrait,OptionValue["ColorPattern"]]},
Text[ToString[NumberForm[MinTrait,{3,2}]],{0,Max[TreeXY[[1;;,3]]]+LegendFactor},{-1,-2}],
Text[ToString[NumberForm[MaxTrait,{3,2}]],{Max[TreeXY[[1;;,2]]]/2,Max[TreeXY[[1;;,3]]]+LegendFactor},{1,-2}],
Text[ToString[NumberForm[ancvalue,{3,2}]],{((ancvalue-MinTrait)/(MaxTrait-MinTrait))*Max[TreeXY[[1;;,2]]]/2,Max[TreeXY[[1;;,3]]]+LegendFactor},{0,-2}]}
},
Axes->{True,False},AxesOrigin->{0,0},AspectRatio->Log[Length[TipNames]]/GoldenRatio]

]; 

];



(*  This function estimates Blomberg's Kappa for a trait given a tree.

	  Usage:  Kappa[MyTree, MyTrait]
	  
	  where MyTree is a Newick-format tree in the format of ReadNewick[]and MyTrait
	  is a vector of taxon names in the first column and trait values in the second column.
	  
	  Created by P. David Polly, 21 December 2018
	  14 February 2019: Fixed trait sorting bug

  *)

Kappa[tree_,trait_]:=Block[{y,varA,varY,varAY,C,A,EofY,num,denom,TipLabels,InvC},
{varY,varAY, varA}=PhylogeneticMatrices[tree];
C=varY[[2;;,2;;]];
InvC=Inverse[C];
TipLabels=varY[[2;;,1]];
y=trait[[Flatten[Position[trait[[1;;,1]],#]&/@TipLabels]]];
A=OMearaAncestor[tree,y][[1]];
EofY=Table[A,{Length[y]}];
num=((y[[1;;,2]]-EofY) . (y[[1;;,2]]-EofY))/((y[[1;;,2]]-EofY) . InvC . (y[[1;;,2]]-EofY));
denom=(Tr[C]-Length[trait]/(Table[1,{Length[y]}] . InvC . Table[1,{Length[y]}]))/(Length[y]-1);
Return[num/denom];
]




(*  This function estimates Adams's (2014)  multivariate Kappa from tree and univariate trait vector.

	  Usage:  KappaMultivariate[MyTree, MyTrait]
	  
	  where MyTree is a Newick-format tree in the format of ReadNewick[]and MyTrait
	  is a matrix containing a vector of taxon names in the first column and trait values in the second and subsequent columns.
	  
	  Created by P. David Polly, 21 December 2018

  *)

KappaMultivariate[tree_,traits_]:=Block[{varY,varAY,varA,Cmat,TipLabels,y,Imat,Ivec,a,EofY,Dya,eigvects,eigvals,Emat,Uy,PDuo,num, denom, kappamult},
{varY,varAY, varA}=PhylogeneticMatrices[tree];
Cmat=varY[[2;;,2;;]];
TipLabels=varY[[2;;,1]];
y=traits[[1;;,2;;]][[Flatten[Position[traits[[1;;,1]],#]&/@TipLabels]]];
Imat=IdentityMatrix[Length[y]];
Ivec=Table[1,{Length[y]}];
a=(Ivec . Inverse[Cmat] . y)/(Ivec . Inverse[Cmat] . Ivec);
EofY=Table[a,{Length[y]}];
Dya=Sqrt[Plus@@#]&/@((y-EofY)^2);
{eigvects,eigvals,v}=SingularValueDecomposition[Cmat];
Emat=Inverse[eigvects . Sqrt[eigvals] . Transpose[eigvects]];
Uy=Emat . (y-EofY);
PDuo=Sqrt[Plus@@(#^2)]&/@Uy;
num=(Dya . Dya)/(PDuo . PDuo);
denom=(Tr[Cmat]-(Length[y]/(Ivec . Inverse[Cmat] . Ivec)))/(Length[y]-1);
kappamult=num/denom;
Return[kappamult];
]




(* This function calculates a rate (or rate matrix) assuming Brownian motion based on the equations presented by 
	O'Meara et al., Revell and Harmon, Adams, and others.  
	
	Usage:  OMearaRate[MyTree, MyTrait]
	  
	  where MyTree is a Newick-format tree in the format of ReadNewick[]and MyTrait
	  is a matrix containing a vector of taxon names in the first column and trait values in the second and subsequent columns.
	  
	  Created by P. David Polly, 21 December 2018

	
	*)

OMearaRate[tree_,trait_]:=Block[{One,Cmatrix,a,resids,rate,treelabels,treedata,treeorder,data},
data=trait[[1;;,2;;]];
treelabels=PhylogeneticMatrices[tree][[1,1,2;;]];
treeorder=Flatten[Position[trait[[1;;,1]],#]&/@treelabels];
treedata=data[[treeorder]];
treedata=data[[treeorder]];
Cmatrix=PhylogeneticMatrices[tree][[1,2;;,2;;]];
One=Table[1,{Length[Cmatrix]}];
a=(One . Inverse[Cmatrix] . treedata)/(One . Inverse[Cmatrix] . One);
resids=#-a&/@treedata;
rate=(Transpose[resids] . Inverse[Cmatrix] . resids)/Length[treedata];
Return[rate//N];
]





(* This function traces the branches of a tip taxon to the base of the tree from a tree table.  
	It returns the branch entries of the table.  This function is used by other functions in the 
	package, but it may be useful on its own.
	
	Usage:  taxon is the name of a taxon as a string and treetable is a tree that has been converted to a four-column
	table using the TreeToTable[] function.
	
	Created by P. David Polly, 22 December 2018.
	Updated 27 December 2018 to adapt for tree tables where Node 0 is absent.

*)
GetFullLineage[taxon_,treetable_]:=Block[{workingbranches},
workingbranches=Select[treetable,#[[1]]==taxon&];
While[Length[Flatten[Select[treetable,#[[1]]==workingbranches[[-1,2]]&]]]>0,workingbranches=Append[workingbranches,Flatten[Select[treetable,#[[1]]==workingbranches[[-1,2]]&]]]];
Return[workingbranches];
]




(*   This function takes a set of taxon names and a tree and returns:  
	(1) the name of their last common ancestor; 
	(2) the number of taxa in the selected set as a proportion of the total in the clade, 
	(3) the time of their last common ancestor proportional to the total time between base 
	     of tree and latest occurring tip; and 
	(4) proportion of the lineage length shared by the taxon set and the total lineage length on the tree 
	
	Usage: taxa is a list of taxon names as strings (all taxa must be in the tree and none must 
	be entered more than once) and tree is a phylogenetic tree in Newick format.

	Created by P. David Polly, 22 December 2018
*)

SharedPhylogeneticHistory[taxa_,tree_]:=Block[{VarY,VarAY,VarA,mytable,totallength,totaldepth,workingbranches,mytaxa,totaltaxa,MyLineages,MySharedBranches,SortOrder,MySharedTime,MyCA,MyUniqueBranches,MyCladeLength,MissingTaxa,mytally},

If[StringQ[tree],mytable=TreeToTable[tree],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},mytable=tree,
If[Length[tree[[1]]]==4,mytable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];
];
];

mytaxa=taxa;
mytally=Tally[mytaxa];
If[Max[mytally[[1;;,2]]]>1,Return[Print["The following taxa were entered more than once: "<>(#[[1]]<>", "&/@Select[mytally,#[[2]]>1&])]]];
{VarY,VarAY,VarA}=PhylogeneticMatrices[tree];
totaltaxa=Length[Select[mytable,#[[4]]==1&]];
totallength=Plus@@mytable[[2;;,3]];
totaldepth=Max[Tr[VarY[[2;;,2;;]],List]];
workingbranches=Flatten[Table[Select[mytable,#[[1]]==mytaxa[[x]]&],{x,Length[mytaxa]}],1];
If[Length[workingbranches]!=Length[mytaxa],MissingTaxa=Complement[mytaxa,workingbranches[[1;;,1]]];Return[Print["The following taxa were not found: "<>(#<>", "&/@MissingTaxa)]]];
MyLineages=GetFullLineage[#,mytable]&/@mytaxa;
MySharedBranches=FoldList[Intersection,MyLineages[[1]],MyLineages[[2;;]]][[-1]];
SortOrder=Ordering[ToExpression[#]&/@StringSplit[Select[MySharedBranches, StringMatchQ[#[[1]],"Node*"]&][[1;;,1]]][[1;;,2]]];
MySharedBranches=MySharedBranches[[SortOrder]];
MySharedTime=Plus@@MySharedBranches[[1;;,3]];
If[Length[MySharedBranches]==0,MyCA="Node 0",MyCA=MySharedBranches[[-1,1]]];
MyUniqueBranches=Complement[Union[Flatten[MyLineages,1]],MySharedBranches];
MyCladeLength=Plus@@MyUniqueBranches[[1;;,3]];
Return[{MyCA,NumberForm[Length[mytaxa]/totaltaxa//N,{3,2}],NumberForm[MySharedTime/totaldepth//N,{3,2}],NumberForm[MyCladeLength/totallength//N,{3,2}]}];
];




(*  This function takes returns a Newick tree composed of a subset of taxa selected
	from a larger tree.  

	Usage: taxa is a list of taxon names as strings (all taxa must be in the tree and none must 
	be entered more than once) and tree is a phylogenetic tree in Newick format or as a
	tree table.

	Created by P. David Polly, 23 December 2018
	Updated 27 December 2018 to fix bug that caused function to fail with large trees.
*)

PruneTree[mytaxa_,tree_]:=Block[{MyTable,MyNodeNumbers, MyLineages,MySharedBranches,MyUniqueBranches,MyNewTree,MySingletonNodes,MySingletonNodePositions,MyAncestors,MyAncestorNodePositions,MyNodes,NewNodeNumbers,NodeChange,SortOrder,MissingTaxa},

If[StringQ[tree],MyTable=TreeToTable[tree],
If[Length[tree]<=1,Return["Tree variable must be a Newick tree string or a tree table."],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},MyTable=tree,
Return["Tree variable must be a Newick tree string or a tree table."];
];
];
];

MissingTaxa=Complement[mytaxa,MyTable[[2;;,1]]];If[Length[MissingTaxa]>0,Return["The following taxa were not found: "<>(#<>", "&/@MissingTaxa)]];
MyLineages=GetFullLineage[#,MyTable]&/@mytaxa;
MySharedBranches=FoldList[Intersection,MyLineages[[1]],MyLineages[[2;;]]][[-1]];
SortOrder=Ordering[ToExpression[#]&/@StringSplit[Select[MySharedBranches, StringMatchQ[#[[1]],"Node*"]&][[1;;,1]]][[1;;,2]]];
MySharedBranches=MySharedBranches[[SortOrder]];
MyUniqueBranches=Complement[Union[Flatten[MyLineages,1]],MySharedBranches];
MyNewTree=Prepend[MyUniqueBranches,MyTable[[1]]];

Do[
MySingletonNodes=Select[Tally[MyNewTree[[2;;,2]]],#[[2]]==1&][[1;;,1]];
If[Length[MySingletonNodes]==0,Break[]];
MySingletonNodePositions=Flatten[Position[MyNewTree[[1;;,2]],#]&/@MySingletonNodes];
MyAncestorNodePositions=Flatten[Position[MyNewTree[[1;;,1]],#]&/@MySingletonNodes];
MyNewTree[[MyAncestorNodePositions[[1]],1]]=MyNewTree[[MySingletonNodePositions[[1]],1]];
MyNewTree[[MyAncestorNodePositions[[1]],3]]=MyNewTree[[MyAncestorNodePositions[[1]],3]]+MyNewTree[[MySingletonNodePositions[[1]],3]];
MyNewTree=Drop[MyNewTree,{MySingletonNodePositions[[1]]}]
,{Infinity}];
MyNodes=Union[MyNewTree[[2;;,2]]];
MyNodeNumbers=ToExpression[#]&/@((StringSplit[#]&/@MyNodes)[[1;;,2]]);
NewNodeNumbers=MyNodeNumbers-Min[MyNodeNumbers];
NodeChange=Table[MyNodes[[x]]->"Node "<>ToString[NewNodeNumbers[[x]]],{x,Length[MyNodes]}];
MyNewTree=MyNewTree/.NodeChange;
Do[If[StringMatchQ[MyNewTree[[x,1]],"Node *"],MyNewTree[[x,4]]=0,MyNewTree[[x,4]]=1],{x,2,Length[MyNewTree],1}];
MyNewTree=TableToTree[MyNewTree]; 
Return[MyNewTree];
];




(* This function calculates mean tree distance for a subset of taxa in a larger clade.
	The metric is the mean pairwise distance (in branch length units) between the selected
	of taxa and is analogous to Webb's (2000) mean pairwise nodal distance.  This function
	assumes that the input tree is ultrametric because it uses depth of shared ancestry 
	of first taxon in the matrix as measure of depth of the tree as a whole.
	
	Input is a list of taxon names and the VarY matrix from the PhylogeneticMatrices[] function.    
	The function returns the size of the community and mean distance.
	
	Created by P. David Polly, 23 December 2018
  *)

MeanTreeDistance[community_,tree_]:=Block[{n,i,j,communitypairs,TreeDistance,taxa,mytable,VarY,mytally},

If[StringQ[tree],mytable=TreeToTable[tree],
If[tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},mytable=tree,
If[Length[tree[[1]]]==4,mytable=tree,Print["Tree variable must be a Newick tree string or a tree table."];];
];
];

mytally=Tally[community];
If[Max[mytally[[1;;,2]]]>1,Return[Print["The following taxa were entered more than once: "<>(#[[1]]<>", "&/@Select[mytally,#[[2]]>1&])]]];

VarY=PhylogeneticMatrices[mytable][[1]];

taxa=VarY[[1,2;;]];
TreeDistance=VarY[[2,2]]-VarY[[2;;,2;;]];
communitypairs=Flatten[Table[Table[Flatten[{Position[taxa, community[[i]]],Position[taxa,community[[j]]]}],{i,j+1,Length[community]}],{j,Length[community]}],1];
{Length[community],Mean[TreeDistance[[#[[1]],#[[2]]]]&/@communitypairs]//N}
]





(* This function performs a phylogenetic regression on a univariate trait using the algorithm of Martins and Hansen (1996).

	Input is a matrix of data with taxon names in the first column and trait values in the second and third columns
	and a tree in Newick format or as a tree table.
	
	Created by P. David Polly, 23 December 2018.
	Corrected error in MS caused by dividing by wrong degree of freedom, which 
	propogated the mistake into F-ration and P-value.   28 Sept 2019.
*)

PhylogeneticRegression[MyData_,MyTree_]:=Block[{X,Y,params,Cmatrix,treelabels,originalorder,treeorder,treedata,a,One,expected,residuals,
SSE,SST,SSR,SERegression,RSquare,DF,SS,MS,F,P,ANOVATable},
treelabels=PhylogeneticMatrices[MyTree][[1,1,2;;]];
treeorder=Flatten[Position[MyData[[1;;,1]],#]&/@treelabels];
originalorder=Flatten[Position[treelabels,#]&/@MyData[[1;;,1]]];
treedata=MyData[[treeorder,{2,3}]];
Cmatrix=PhylogeneticMatrices[MyTree][[1,2;;,2;;]];
X=Transpose[{Table[1,{Length[treedata]}],treedata[[1;;,1]]}];
Y=treedata[[1;;,2]];
params=Inverse[Transpose[X] . Inverse[Cmatrix] . X] . Transpose[X] . Inverse[Cmatrix] . Y;


expected=((params[[2]]*#)+params[[1]])&/@treedata[[1;;,1]];
residuals=Plus@@{-1*expected,treedata[[1;;,2]]};
SSE=Plus@@(residuals^2);
SST=Plus@@((#-Mean[Y]&/@treedata[[1;;,2]])^2);
SSR=SST-SSE;
SERegression=Sqrt[SSE/(Length[residuals]-2)];
RSquare=SSR/SST;
DF={Length[residuals]-1,Length[residuals]-2,1};
SS={SST,SSE,SSR};
MS=SS/DF;
F=MS[[3]]/MS[[2]];
P=1-CDF[FRatioDistribution[DF[[3]],DF[[2]]],F];
ANOVATable=TableForm[{{DF[[3]],SS[[3]],MS[[3]],NumberForm[F,{6,4}],NumberForm[P,{6,5}]},
{DF[[2]],SS[[2]],MS[[2]]},
{DF[[1]],SS[[1]]}},TableHeadings->{{"Model","Error","Total"},{"DF","SS","MS","F","P"}}];


Return[{#[[1]]->#[[2]]&/@Transpose[{{"Intercept","Slope"},params}],"R-Squared"->NumberForm[RSquare,{3,2}],"ANOVA Table"->ANOVATable,"SE Regression"->SERegression,"Predicted Values"->expected[[originalorder]]}];
]




(* This function transforms a phylogenetic covariance matrix (the C matrix of authors, VarY of Martins and Hansen) into a Newick format tree.

	The first input parameter contains a phylogenetic covariance matrix in the same format as produced in element 1 (VarY) by the PhylogeneticMatrices[] 
	function minus its labels.  In this format the diagonal contains the distances between root and repsecitve tip and the off diagonal elements contain 
	the shared branchlength of two tips.  The second input parameter is a list of taxa in the same order as the matrix.  
	
	Created by P. David Polly, 29 June 2019.
*)


CtoTree[BareC_,MyKin_]:=Block[{FindLineage,myTable,UsedNodes,StillToProcess,NodeNumCount,count,MyAnc,MyPositions,MyTransform,CladeC,MyDepth,MyClades,SubCladeC,MyLineageLength,MyNodeName,MyLine},
If[Length[MyKin]!=Length[BareC],Return["Taxon labels do not match phylogenetic covariance matrix."]];

FindLineage[Point_]:=Block[{Lineage,NewBranch},
Lineage={Point};
While[!StringMatchQ[Lineage[[-1,2]],"Node 0"],
NewBranch=Flatten[myTable[[Flatten[Position[myTable[[1;;,1]],Lineage[[-1,2]]]]]]];
Lineage=Append[Lineage,NewBranch];
];
Return[Lineage];
];

myTable=Transpose[{MyKin,Table["Node 0",{Length[MyKin]}],Tr[BareC,List],Table[1,{Length[MyKin]}]}];
UsedNodes={};
StillToProcess={"Node 0"};
NodeNumCount=0;

While[Length[StillToProcess]>=1,
count=1;
Clear[MyDepth];
MyAnc=StillToProcess[[1]];
MyPositions=Flatten[Position[MyKin,#]&/@Sort[Select[myTable,#[[2]]==MyAnc&&#[[4]]==1&][[1;;,1]]]];
MyTransform=#[[2]]->#[[1]]&/@Transpose[{MyPositions,Range[Length[MyPositions]]}];
CladeC=Chop[Table[BareC[[MyPositions[[i]],MyPositions[[j]]]],{i,Length[MyPositions]},{j,Length[MyPositions]}]];
MyDepth=Min[CladeC];
CladeC=Chop[Table[BareC[[MyPositions[[i]],MyPositions[[j]]]],{i,Length[MyPositions]},{j,Length[MyPositions]}]-MyDepth];
MyClades=Union[Table[Sort[Complement[Range[Length[CladeC]],Flatten[Position[CladeC[[n]],0]]]],{n,Length[CladeC]}]]/.MyTransform;
StillToProcess=Select[StillToProcess,#!=MyAnc&];
UsedNodes=Append[UsedNodes,MyAnc];

Do[
MyPositions=MyClades[[node]];
SubCladeC=Table[Chop[BareC[[MyPositions[[i]],MyPositions[[j]]]]],{i,Length[MyPositions]},{j,Length[MyPositions]}];
MyDepth=Min[SubCladeC];
MyLineageLength=Chop[If[MyAnc=="Node 0",0,Plus@@FindLineage[Flatten[Select[myTable,#[[1]]==MyAnc&]]][[1;;,3]]]];
SubCladeC=Chop[SubCladeC-MyDepth];

If[Length[MyClades[[node]]]==1,StillToProcess=Select[StillToProcess,#!=MyAnc&];UsedNodes=Append[UsedNodes,MyAnc];Continue[]];
NodeNumCount=NodeNumCount+1;
MyNodeName="Node "<>ToString[NodeNumCount];
StillToProcess=Append[StillToProcess,MyNodeName];
myTable=Append[myTable,{MyNodeName,MyAnc,MyDepth-MyLineageLength,0}];

Do[
MyLine=Flatten[myTable[[Flatten[Position[myTable[[1;;,1]],MyKin[[MyClades[[node,n]]]]]]]]];
MyLine[[2]]=MyNodeName;
MyLine[[3]]=MyLine[[3]]-MyDepth+MyLineageLength;
myTable[[Flatten[Position[myTable[[1;;,1]],MyKin[[MyClades[[node,n]]]]]]]]=MyLine,
{n,Length[MyClades[[node]]]}],
{node,Length[MyClades]}
];
count=count+1;
If[count>=Length[MyKin]^2,Break[]];
];
Return[TableToTree[myTable]];
]



(* This function performs a phylogenetic regression on a univariate trait using the algorithm of Martins and Hansen (1996), but allows the user to scale the tree branch lengths by a value of Pagel's lambda.  By looping the function into an exploration of all possible values of lambda, this function can be used to find the lambda value that minimizes the phylogenetic structure in the regression residuals.  

	Input is a matrix of data with taxon names in the first column and trait values in the second and third columns
	and a tree in Newick format or as a tree table and a single number for lambda, which normally ranges from 0 to 1 but can exceed one up to a maximum defined by the proportion of the length of the smallest branch over the tree height..
	
	Created by P. David Polly, 8 December 2022.
*)

PhylogeneticRegressionWithLambda[MyData_,MyTree_,lambda_]:=Block[{X,Y,params,Cmatrix,treelabels,originalorder,treeorder,treedata,a,One,expected,residuals,
SSE,SST,SSR,SERegression,RSquare,DF,SS,MS,F,P,ANOVATable},
treelabels=PhylogeneticMatrices[MyTree][[1,1,2;;]];
treeorder=Flatten[Position[MyData[[1;;,1]],#]&/@treelabels];
originalorder=Flatten[Position[treelabels,#]&/@MyData[[1;;,1]]];
treedata=MyData[[treeorder,{2,3}]];
Cmatrix=PhylogeneticMatrices[MyTree][[1,2;;,2;;]];
Cmatrix=LambdaMultiplication[Cmatrix,lambda];
X=Transpose[{Table[1,{Length[treedata]}],treedata[[1;;,1]]}];
Y=treedata[[1;;,2]];
params=Inverse[Transpose[X] . PseudoInverse[Cmatrix] . X] . Transpose[X] . PseudoInverse[Cmatrix] . Y;


expected=((params[[2]]*#)+params[[1]])&/@treedata[[1;;,1]];
residuals=Plus@@{-1*expected,treedata[[1;;,2]]};
SSE=Plus@@(residuals^2);
SST=Plus@@((#-Mean[Y]&/@treedata[[1;;,2]])^2);
SSR=SST-SSE;
SERegression=Sqrt[SSE/(Length[residuals]-2)];
RSquare=SSR/SST;
DF={Length[residuals]-1,Length[residuals]-2,1};
SS={SST,SSE,SSR};
MS=SS/DF;
F=MS[[3]]/MS[[2]];
P=1-CDF[FRatioDistribution[DF[[3]],DF[[2]]],F];
ANOVATable=TableForm[{{DF[[3]],SS[[3]],MS[[3]],NumberForm[F,{6,4}],NumberForm[P,{6,5}]},
{DF[[2]],SS[[2]],MS[[2]]},
{DF[[1]],SS[[1]]}},TableHeadings->{{"Model","Error","Total"},{"DF","SS","MS","F","P"}}];


Return[{#[[1]]->#[[2]]&/@Transpose[{{"Intercept","Slope"},params}],"R-Squared"->NumberForm[RSquare,{3,2}],"ANOVA Table"->ANOVATable,"SE Regression"->SERegression,"Predicted Values"->expected[[originalorder]],"Residuals"->residuals[[originalorder]]}];
]


(* This function takes as input a non-ultrametric tree and creates an ultrametric tree by lengthening all tip branches to be coeval with the latest occurring tip.  

	Input is a tree in Newick format or as a tree table.  Output is a Newick format ultrametric tree.
	
	Created by P. David Polly, 29 May 2023
*)

UltrametricizeTree[tree_]:=Block[{Cmat,tiplabels,myheight,myultrametrictree},
Cmat=PhylogeneticMatrices[tree][[1]];
tiplabels=Cmat[[1,2;;]];
Cmat=Cmat[[2;;,2;;]];
myheight=Max[Tr[Cmat,List]];
Do[Cmat[[x,x]]=myheight,{x,Length[Cmat]}];
myultrametrictree=CtoTree[Cmat,tiplabels];
Return[myultrametrictree];
]
