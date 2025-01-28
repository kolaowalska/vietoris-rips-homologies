(* ::Package:: *)

(* non-trivial examples for calculating homology groups *)

TORUS = {{1},{2},{3},{4},{5},{6},{7},{8},{9},
{1,2},{1,4},{1,3},{1,5},{1,7},{1,8},{2,6},{2,3},{2,9},{2,8},{2,4},{3,5},{3,6},{3,7},{3,9},{4,5},{4,6},{4,7},{4,9},{5,6},{5,8},{5,9},{6,7},{6,8},{7,8},{7,9},{8,9},
{1,2,4},{2,4,6},{2,3,6},{3,6,7},{1,3,7},{1,4,7},{4,5,6},{5,6,8},{6,7,8},{7,8,9},{4,7,9},{4,5,9},{1,5,8},{1,2,8},{2,8,9},{2,3,9},{3,5,9},{1,3,5}};

KLEINBOTTLE = {{1},{2},{3},{4},{5},{6},{7},{8},{9},
{1,2},{1,4},{1,3},{1,5},{1,6},{1,9},{2,6},{2,3},{2,7},{2,8},{2,5},{3,5},{3,8},{3,7},{3,9},{4,5},{4,6},{4,7},{4,9},{4,8},{5,7},{5,8},{6,7},{6,8},{6,9},{7,9},{8,9},
{1,2,6},{1,4,6},{2,3,7},{2,6,7},{1,3,5},{3,5,7},{4,5,8},{4,6,8},{6,7,9},{6,8,9},{4,7,9},{4,5,7},{1,5,2},{2,5,8},{3,8,9},{2,3,8},{1,3,9},{1,4,9}};

S2 ={{1},{2},{3},{4},
{1,3},{1,2},{1,4},{2,3},{2,4},{3,4},
{1,2,4},{2,3,4},{1,3,4},{1,3,2}};

RP = {{1},{2},{3},{4},{5},{6},
{1,2},{1,3},{1,4},{1,5},{1,6},{2,3},{2,4},{2,5},{2,6},{3,4},{3,5},{3,6},{4,5},{4,6},{5,6},
{1,2,5},{1,2,6},{1,3,4},{1,3,6},{1,4,5},{2,3,4},{2,3,5},{2,4,6},{3,5,6},{4,5,6}};

EXAMPLEDIM2 = {{17,15},{13,12},{19,10},{11,18},{16,8},{13,16},{11,16}};

EXAMPLE2DIM2 = {{13,8},{17,13},{14,13},{14,8},{17,15},{14,15},{16,13}};

EXAMPLEDIM3 = {{15,9,7},{14,8,5},{16,9,8},{15,8,8},{13,10,5},{14,11,8},{16,1,6},{9,4,5},{14,6,8},{9,10,5}};

EXAMPLEDIM5 = {{12,1,11,2,4},{11,1,12,4,3},{12,2,12,2,0},{12,1,17,4,0},{12,2,11,7,1},{11,1,14,3,3},{12,1,16,6,2},{12,2,15,7,3}};


(* a function used to generate k random points of dimension n *)

GenerateRandomSet[n_, k_, LowerLeft_, UpperRight_] := Module[
    {Points},
    If[Length[LowerLeft] != n || Length[UpperRight] != n,
        Return["Boundary points should be of given length n."]
    ];
    (* generating k points R^n *)
    Points = Table[
        (* using integers here for simplicity but this can be changed to RandomReal *)
        RandomInteger[{LowerLeft[[i]], UpperRight[[i]]}], {i, n}, {j, k}
    ];
    Transpose[Points]
];


(* a function which fills in Euclidean distances between point coordinates 
for later checking upon the Vietoris-Rips condition *)

GenerateDistanceMatrix[points_, r_] := Module[
  {n = Length[points], adjacency},
  adjacency = Table[
    EuclideanDistance[points[[i]], points[[j]]],
    {i, n}, {j, n}
  ];
  adjacency
];


(* a function for constructing the Vietoris-Rips simplicial complex with the following arguments:
- points_: a list of points which act as vertices for the complex
- distanceMatrix_: a matrix of precomputed Euclidean distances between coordinate pairs
- r_: the given radius parameter
- maxDim_: a rough estimate of the dimension up to which the calculations will be sensible
for the user input, it can obviously be increased but it is taken as an argument to avoid
redundant calculations for simpler examples *)

VietorisRipsComplex[points_, distanceMatrix_, r_, maxDim_ : 10] := Module[{
    n = Length[points],
    dims, k, combos
  },
  
  (* validating input *)
  If[Length[points] == 0, Return[Association[]]];
  
  (* using a dictionary to store k-simplices to speed up later calls *)
  dims = Association[];
  
  (* dim0 = all singletons {i} for i in 1..n *)
  dims["dim0"] = Table[{i}, {i, n}];
  
  (* for each dimension from 1 up to maxDim gathering all subsets of size 
  (k+1) whose pairwise distances are <= r using the distance matrix *)
  Do[
    combos = Subsets[ Range[n], {k + 1} ];
    combos = Select[
      combos,
      AllTrue[Subsets[#, {2}], distanceMatrix[[#[[1]], #[[2]]]] <= r &] &
    ];
    dims["dim" <> ToString[k]] = combos;
    ,
    {k, 1, maxDim}
  ];
  Return[dims];
];


(* a function analogous to the Vietoris-Rips function, but intead of utilizing point coordinates
it builds the complex based on a list of input simplices in curly-bracket representation. this
is used for examples in graph representation like the ones at the very top of this notebook *)

BuildSimplicialComplex[inputSimplices_, maxDim_ : 10] := Module[{
    dims, sortedSimplices, n
  },
  
  (* validating input *)
  If[Length[inputSimplices] == 0, Return[Association[]]];
  
  (* removing duplicates, just a precaution for long lists *)
  sortedSimplices = DeleteDuplicates[Sort /@ inputSimplices];
  
  (* finding the maximum vertex number *)
  n = Max[Flatten[sortedSimplices]];
  
  dims = Association[];
  
  dims["dim0"] = Select[sortedSimplices, Length[#] == 1 &];
  
  (* for each dimension from 1 up to maxDim *)
  Do[
    (* select simplices of current dimension *)
    dims["dim" <> ToString[k]] = Select[
      sortedSimplices, 
      Length[#] == k + 1 &
    ];
    ,
    {k, 1, maxDim}
  ];
  
  Return[dims];
];


(* a function to draw the neighborhoods around the given points to enable when toggling *)

DrawNeighborhoods[showNeighborhoods_, vertices_, r_, showIndices_] :=
  If[
    showNeighborhoods,
    Table[
      {Hue[0.9, 0.9, 0.4, 0.1], Disk[vertices[[i]], r/2]},
      {i, showIndices}
    ],
    {}
  ];
  
(* the functions below are utilized in Manipulate to visualize the complex,
drawing one- and two-simplices is the basis for the representation of 
higher-dimensional simplices for this 2D visualization. as arguments they
receive that same list of vertices, the distance matrix and the radius 
parameter as above, and additionaly a showIndices_ list which is the list of 
vertices that have been chosen to be shown with Manipulate by toggling *)

DrawOneSimplices[vertices_, distanceMatrix_, r_, showIndices_] := 
  Table[
    If[
    (* checking the Vietoris-Rips condition for distances between each pair of points *)
      i < j &&
      (**)
      MemberQ[showIndices, i] &&
      MemberQ[showIndices, j] &&
      distanceMatrix[[i, j]] <= r,
      {AbsoluteThickness[1.8], Black, Line[{vertices[[i]], vertices[[j]]}]},
      {}
    ],
    {i, Length[vertices]}, {j, i + 1, Length[vertices]}
  ] // Flatten;

DrawTwoSimplices[vertices_, distanceMatrix_, r_, showIndices_] := 
  Table[
    If[
    (* checking the Vietoris-Rips condition for distances between each triple of points *)
      distanceMatrix[[i, j]] <= r &&
      distanceMatrix[[i, k]] <= r &&
      distanceMatrix[[j, k]] <= r &&
      MemberQ[showIndices, i] &&
      MemberQ[showIndices, j] &&
      MemberQ[showIndices, k],
      {
        Hue[0.6, 0.9, 0.5, 0.1],
        EdgeForm[Black],
        Polygon[{vertices[[i]], vertices[[j]], vertices[[k]]}]
      },
      {}
    ],
    {i, Length[vertices]}, {j, i + 1, Length[vertices]}, {k, j + 1, Length[vertices]}
  ] // Flatten;


(* a Manipulate function enabling visualization utilizing the above functions *)

ManipulateComplex[vertices_, distanceMatrix_, R_, LowerLeftBound_, UpperRightBound_] := 
  Manipulate[
    Graphics[
      Flatten[
        {
          DrawNeighborhoods[showNeighborhoods, vertices, r, selectedPoints],
          DrawOneSimplices[vertices, distanceMatrix, r, selectedPoints],
          DrawTwoSimplices[vertices, distanceMatrix, r, selectedPoints],
          
          (* drawing vertices as points *)
          Table[{Black, Point[vertices[[i]]]}, {i, selectedPoints}]
          (* Table[{Black, PointSize[0.01], Point[vertices[[i]]]}, {i, selectedPoints}] *)
        }
      ],
      PlotRange -> {
        {LowerLeftBound[[1]] - R, UpperRightBound[[1]] + R},
        {LowerLeftBound[[2]] - R, UpperRightBound[[2]] + R}
      },
      PlotRangeClipping -> False,
      ImageSize -> 600
    ],
    
    {{showNeighborhoods, True, "Draw neighborhoods"}, {True, False}},
    {{r, R, "Neighborhood radius"}, 0, R + 5, Appearance -> "Labeled"},
    {{selectedPoints, Range[Length[vertices]], "Points to show"},
      Range[Length[vertices]],
      ControlType -> CheckboxBar
    },
    Delimiter,
    Style["Vietoris-Rips simplicial complex", 12, Bold]
  ];


(* counting simplices of each dimension using the previously implemented dictionary *)

CountSimplices[complex_] := Module[
  {dimKeys, counts},
  
  dimKeys = Keys[complex];
  counts = Table[
    {StringReplace[dim, "dim" -> ""], Length[complex[dim]]},
    {dim, dimKeys}
  ];

  Dataset[counts]
];


(* a function for drawing an indicative graph representation of the simplicial 
complex based on the vertices and edges; this was mainly used for double-checking 
whether the number of connected components (and therefore Betti0) is correct *)

GenerateGraph[complex_, n_] := Module[
  {edges, vertices},
  
  edges = complex["dim1"];
  vertices = Range[n];
  
  Graph[
    vertices,
    UndirectedEdge @@@ edges,
    VertexLabels -> Table[i -> Placed[i, Center], {i, vertices}],
    VertexSize -> 0.5
  ]
];


VisualizeComplex = GenerateGraph[SimplicialComplex, CardSn[[1, 2]]];


(* harvesting the connected components of the graph to later obtain 
the 0-th Betti number and avoid unnecessary matrix calculations *)

ShowConnectedComponents[complex_] := Module[
  {edges, vertices, graph},

  edges = complex["dim1"];
  vertices = Range[Length[complex["dim0"]]];
  
  graph = Graph[
    vertices,
    UndirectedEdge @@@ edges
  ];
  
  ConnectedComponents[graph]
];


ComplexConnectedComponents = ShowConnectedComponents[SimplicialComplex];
Betti0 = Length[ComplexConnectedComponents];


(* counting dimensions in which any simplices in the simplicial complex exist,
so that there are no redundant calculations for homology groups that are trivial *)

(* for example, if there are simplices only up to dimension 3, homology groups 
of higher dimensions needn't be calculated since they will always be 0 *)

CountPotentiallyNontrivialHomologies[complex_] := Module[
  {counts},
  
  counts = Count[Length[complex[#]] > 0 & /@ Keys[complex], True];
  counts
];


(* a helper function used to transform border matrices into column echelon form, 
since the built-in Mathematica function RowReduceAlternative[matrix, Modulus -> 2]
does not seem to work for non-square border matrices *)

ColumnEchelonForm[matrix_] := Module[
  {mat = matrix, rows, cols},
  {rows, cols} = Dimensions[mat];
  
  (* helper function to find the leading index in a column *)
  LeadingIndex[col_] := FirstPosition[col, 1, \[Infinity]][[1]];
  
  (* helper function to check whether the column is 0 *)
  IsZeroColumn[col_] := AllTrue[col, # == 0 &];
  
  Do[
    (* finding column with the minimum leading index among remaining columns *)
    Module[{
      remainingCols = Range[j, cols],
      leadingIndices,
      minLeadIndex,
      colWithMinLead
    },
      (* calculating leading indices for remaining columns *)
      leadingIndices = Map[LeadingIndex[mat[[All, #]]] &, remainingCols];
      
      (* finding the minimum leading index and its column *)
      minLeadIndex = Min[leadingIndices];
      
      (* if we find a non-zero column *)
      If[minLeadIndex < \[Infinity],
        (* locate the first column with minimum leading index *)
        colWithMinLead = remainingCols[[FirstPosition[leadingIndices, minLeadIndex][[1]]]];
        
        (* swap if not already in position *)
        If[colWithMinLead != j,
          mat[[All, {j, colWithMinLead}]] = mat[[All, {colWithMinLead, j}]];
        ];
        
        (* add column j to all other columns that have 1 in the same leading position *)
        Do[
          If[k != j && mat[[minLeadIndex, k]] == 1,
            mat[[All, k]] = Mod[mat[[All, k]] + mat[[All, j]], 2];
          ],
          {k, j + 1, cols}
        ];
      ];
    ];
    ,
    {j, 1, cols}
  ];
  
  (* moving all zero columns to the end *)
  Module[{
    nonZeroCols = Select[Range[cols], !IsZeroColumn[mat[[All, #]]] &],
    zeroCols = Select[Range[cols], IsZeroColumn[mat[[All, #]]] &]
  },
    mat = mat[[All, Join[nonZeroCols, zeroCols]]];
  ];
  
  mat
]


(* a function for constructing the simplicial complex border operator matrix for a given dimension *)

BorderMatrix[complex_, k_] := Module[
  {
    kSimplices, (* k-simplices *)  
    kMinusOneSimplices, (*(k-1)-simplices *)
    matrix, columnLabels, rowLabels, rank
  },
  
  (* indentifying the k-simplices and the (k-1)-simplices using the dictionary *)
  kSimplices = complex["dim" <> ToString[k]];
  kMinusOneSimplices = complex["dim" <> ToString[k - 1]];
  
  (* constructing the 0\[Dash]1 boundary matrix based on incidence *)
  matrix = Table[
    If[SubsetQ[kSimplices[[j]], kMinusOneSimplices[[i]]], 1, 0],
    {i, 1, Length[kMinusOneSimplices]},
    {j, 1, Length[kSimplices]}
  ];
  
  matrix = ColumnEchelonForm[matrix];
  
  (* building the string labels for row and column simplices *)
  columnLabels = Map[
    (StringJoin["\:3008", Riffle[ToString /@ #, ","], "\:3009"] &),
    kSimplices
  ];
  rowLabels = Map[
    (StringJoin["\:3008", Riffle[ToString /@ #, ","], "\:3009"] &),
    kMinusOneSimplices
  ];
  
  rank = MatrixRank[matrix];
  
  <|
    "RawMatrix" -> matrix,
    "LabeledMatrix" -> MatrixForm[
      matrix,
      TableHeadings -> {rowLabels, columnLabels}
    ],
    "Rank" -> rank,
    "k-Simplices" -> kMinusOneSimplices,
    "(k-1)-Simplices" -> kSimplices
  |>
];


(* an array which stores the ranks of the matrices for Betti number calculations *)
ranks = {}; 

(* a function to output border matrices which make sense 
and are non-trivial for the given simplicial complex *)

ShowBorderMatrices[SimplicialComplex_, PotentiallyNontrivialHomologies_] := Module[
  {res},
  ranks = Table[0, {PotentiallyNontrivialHomologies - 1}]; 
  
  Do[
    res = BorderMatrix[SimplicialComplex, k];
    Print[Indexed[\[Delta], k], ": "];
    Print[res["LabeledMatrix"]];
    Print["Matrix rank: ", res["Rank"]];
    Print[" "];
    Print[" "];
    ranks[[k]] = res["Rank"];
    ,
    {k, 1, PotentiallyNontrivialHomologies - 1}
  ]
]


(* function to calculate and display Betti numbers *)

DisplayBettiNumbers[SimplicialComplex_, PotentiallyNontrivialHomologies_, ranks_, CardSn_] := 
Module[{
    BettiNumbers, cardSn, rkCurr, rkNext, betti
  },
  
  (* initializing the list with Betti0 which we already have *)
  BettiNumbers = {Betti0};
  
  For[k = 1, k <= PotentiallyNontrivialHomologies - 1, k++,
  
    (* retrieving ardinality for the current dimension k *)
    cardSn = CardSn[[k + 1, 2]]; 
    
    (* retrieving the rank for (k-1) if k > 1; otherwise using 0 *)
    rkCurr = If[k <= Length[ranks], ranks[[k]], 0]; 
    
    (* retrieving the rank for (k+1) if within bounds; otherwise using 0 *)
    rkNext = If[k < Length[ranks], ranks[[k + 1]], 0];
    
    (* calculation of the Betti number using the twierdzenie 6.2.2 ze skryptu formula *)
    betti = cardSn - rkCurr - rkNext;
    
    AppendTo[BettiNumbers, Max[betti, 0]];
  ];
  
  (* displaying the numbers *)
  Column[
    Table[
      Row[{Subscript[\[ScriptCapitalB], i], " = ", BettiNumbers[[i + 1]]}],
      {i, 0, Length[BettiNumbers] - 1}
    ]
  ]
];


(* ::Title:: *)
(*Vietoris-Rips simplicial complex*)


n = 2; (* dimension *)
k = 7; (* number of points *)
r = 5; (* radius *)

(* picking the boundary points at random *)  

RandomBounds = RandomInteger[{0, 20}, {2, n}]; 
LowerLeftBound = Map[Min, Transpose[RandomBounds]];    
UpperRightBound = Map[Max, Transpose[RandomBounds]];
{LowerLeftBound, UpperRightBound}  


(* setting bounds

LowerLeftBound = {0, 0};
UpperRightBound = {10, 10};
*)


(* GENERATED VIETORIS RIPS *)
InputPoints = GenerateRandomSet[n, k, LowerLeftBound, UpperRightBound]
Edges = GenerateDistanceMatrix[InputPoints, r];
SimplicialComplex = VietorisRipsComplex[InputPoints, Edges, r, n];
SimplicialComplex



(*VIETORIS RIPS ENERATED FROM GIVEN SET
InputPoints = EXAMPLE2DIM2
Edges = GenerateDistanceMatrix[InputPoints, r];
SimplicialComplex = VietorisRipsComplex[InputPoints,Edges,r,2];
SimplicialComplex
*)

(*
COUNTING HOMOLOGIES TEST
SimplicialComplex = BuildSimplicialComplex[TORUS, 3]
SimplicialComplex = BuildSimplicialComplex[S2, 3]
SimplicialComplex = BuildSimplicialComplex[KLEINBOTTLE, 5]
SimplicialComplex = BuildSimplicialComplex[RP, 5]
*)


If[n == 2, ManipulateComplex[InputPoints, Edges, r, LowerLeftBound, UpperRightBound]]


Print["Number of simplices in each dimension:"]
NumberOfSimplices = CountSimplices[SimplicialComplex]
CardSn = Normal[NumberOfSimplices];


VisualizeComplex = GenerateGraph[SimplicialComplex, CardSn[[1, 2]]];
ComplexConnectedComponents = ShowConnectedComponents[SimplicialComplex];
Betti0 = Length[ComplexConnectedComponents];
PotentiallyNontrivialHomologies = CountPotentiallyNontrivialHomologies[SimplicialComplex];

Print["Connected components of the simplicial complex: ", ComplexConnectedComponents]
Print["0-th Betti number: ", Betti0]
Print["Potentially non-trivial homologies: ", PotentiallyNontrivialHomologies]


ShowBorderMatrices[SimplicialComplex, PotentiallyNontrivialHomologies]


Print["Betti numbers: "]
DisplayBettiNumbers[SimplicialComplex, PotentiallyNontrivialHomologies, ranks, CardSn]



