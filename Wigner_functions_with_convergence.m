(* ============================================================ *)
(* Absorbing continuous-time quantum walk                       *)
(* Wigner and conditionally normalized Wigner functions         *)
(*                                                              *)
(* This script generates the physical Wigner-function figure    *)
(* and the conditionally normalized Wigner-function figure used *)
(* in the accompanying manuscript.                              *)
(*                                                              *)
(* Wolfram Language / Mathematica                               *)
(* No external packages are required.                           *)
(* ============================================================ *)

ClearAll["Global`*"];

(*------------------------------------------------------------*)
(* 0. COMMON FIGURE STYLE                                     *)
(*------------------------------------------------------------*)
(* Formatting shared by both Wigner-function figures.         *)
plotWidth = 400;
plotHeight = 520;
barHeight = 360;

commonPlotStyle = {InterpolationOrder -> 2, MaxPlotPoints -> Infinity,
    PerformanceGoal -> "Quality", Frame -> True, FrameStyle -> Black, 
   FrameTicks -> {Automatic, Automatic}, 
   FrameTicksStyle -> Directive[Black, 34, FontFamily -> "Times"], 
   LabelStyle -> {FontSize -> 34, FontColor -> Black, 
     FontFamily -> "Times"}, 
   FrameLabel -> {Style["(s+s')/2", Italic, Black, 38, 
      FontFamily -> "Times"], 
     Style["k", Italic, Black, 38, FontFamily -> "Times"]}, 
   PlotRangePadding -> None, AspectRatio -> 1.15, 
   ImageSize -> {plotWidth, plotHeight}};

timeLabel[omegaT_] := 
  Style[Row[{"\[CapitalOmega] ", Style["t", Italic], " = ", 
     NumberForm[omegaT, {4, 1}]}], Black, 34, FontFamily -> "Times"];

legendStyle = {FontSize -> 26, FontColor -> Black, 
   FontFamily -> "Times"};

(*------------------------------------------------------------*)
(*1. Model and plotting parameters*)
(*------------------------------------------------------------*)
(* Parameters below reproduce the Wigner-function panels shown *)
(* in the manuscript. They may be modified to explore other     *)
(* regimes. eta = kappa/Omega is the dimensionless absorption   *)
(* strength used throughout the paper.                           *)
(*------------------------------------------------------------*)

Omega = 1.;
kappa = 1.5;
eta = kappa/Omega;

s0 = 3;
Nsites = 80;

times = {1.0, 1.5, 3.0, 5.0};

(* Maximum doubled-lattice coordinate m displayed at each time. *)
(* This controls only the displayed phase-space window; it is not *)
(* a truncation of an infinite Wigner sum.                        *)
mMaxList = {14, 16, 18, 22};

(* Number of momentum intervals used to sample k in [-Pi, Pi].   *)
(* This is a plotting-resolution parameter, not a physical cutoff. *)
nk = 900;

(*------------------------------------------------------------*)
(*2. Effective non-Hermitian Hamiltonian*)
(*------------------------------------------------------------*)
(*Construct the finite-lattice approximation to*)
(**)
(*H_eff=H_0-i kappa/2|1><1|,*)
(**)
(*where H_0 is the nearest-neighbor tight-binding Hamiltonian.*)
(*The imaginary term removes probability at the boundary*)
(*site s=1.*)
(*------------------------------------------------------------*)

heff[kappa_?NumericQ] := 
  Module[{H}, H = ConstantArray[0. + 0. I, {Nsites, Nsites}];
   (*On-site energy Omega.*)Do[H[[s, s]] = Omega, {s, 1, Nsites}];
   (*Nearest-neighbor hopping-Omega/2.*)Do[H[[s, s + 1]] = -Omega/2;
    H[[s + 1, s]] = -Omega/2, {s, 1, Nsites - 1}];
   (*Absorbing boundary at site s=1.*)H[[1, 1]] -= I kappa/2;
   SparseArray[H]];

Habs = heff[kappa];

(*------------------------------------------------------------*)
(*3. Surviving wavefunction*)
(*------------------------------------------------------------*)
(*Start from a walker localized at site s0 and propagate it*)
(*with the effective non-Hermitian Hamiltonian. Its norm is*)
(*the survival probability S(t).*)
(*------------------------------------------------------------*)

psiFullAt[t_?NumericQ] := Module[{psi0}, psi0 = UnitVector[Nsites, s0];
   MatrixExp[-I Habs t] . psi0];

(*------------------------------------------------------------*)
(*4. Analytic boundary-localized pole mode*)
(*------------------------------------------------------------*)
(*For eta>1 the Green's function contains a normalizable*)
(*boundary-localized mode. This routine evaluates its spatial*)
(*wavefunction directly from the analytic pole expression.*)
(*For eta<=1 no physical localized pole mode is included.*)
(*------------------------------------------------------------*)

qp = -I/eta;

zp = Omega - I Omega/2 (eta - 1/eta);

psiPoleAt[t_?NumericQ] := 
  Table[If[eta > 1, (1 - qp^2) qp^(s + s0 - 2) Exp[-I zp t], 0.], {s, 
    1, Nsites}];

(*------------------------------------------------------------*)
(*5. Doubled-lattice Wigner function*)
(*------------------------------------------------------------*)
(*Given any lattice wavefunction psi, evaluate*)
(**)
(*W_+(m,k)=1/(2 Pi) Sum_n*)
(*psi_n psi^*_(m-n) Exp[-i(2n-m)k].*)
(**)
(*The returned data have the form*)
(**)
(*{(s+s')/2,k,W}*)
(**)
(*and can be passed directly to ListDensityPlot.*)
(*------------------------------------------------------------*)

wignerFromPsi[psi_, mMax_Integer] := 
  Module[{psiLocal, amplitude, mValues, 
    kValues},(*Only amplitudes up to mMax-1 can enter the requested sum.*)
   psiLocal = Take[psi, Min[mMax - 1, Length[psi]]];
   (*Return the stored amplitude inside the lattice range and zero outside it.*)amplitude[s_Integer] := 
    If[1 <= s <= Length[psiLocal], psiLocal[[s]], 0];
   mValues = Range[2, mMax];
   (*Uniform momentum grid from-Pi to Pi.*)
   kValues = Subdivide[-Pi, Pi, nk];
   (*Table first produces one list for each m.Flatten[...,
   1] combines these into one list of {m/2,k,W} triplets for plotting.*)
   Flatten[Table[{m/2, k, 
      Re[1/(2 Pi) Sum[
         amplitude[n] Conjugate[
           amplitude[m - n]] Exp[-I (2 n - m) k], {n, 1, 
          m - 1}]]}, {m, mValues}, {k, kValues}], 1]];

(*------------------------------------------------------------*)
(* 6. PHYSICAL WIGNER SNAPSHOTS                              *)
(*------------------------------------------------------------*)
(*For each plotted time, compute*)
(**)
(*-the full surviving wavefunction,*)
(*-the analytic pole component,*)
(*-the Wigner function of each,*)
(*-the survival probability S(t),*)
(*-the maximum absolute Wigner values used for plotting.*)
(**)
(*Everything is stored in an Association for convenient access*)
(*when constructing the figures.*)
(*------------------------------------------------------------*)

physicalSnapshotAt[t_?NumericQ, mMax_Integer] := 
  Module[{psi, psiPole, fullWigner, poleWigner, survival}, 
   psi = N[psiFullAt[t]];
   psiPole = N[psiPoleAt[t]];
   fullWigner = wignerFromPsi[psi, mMax];
   poleWigner = wignerFromPsi[psiPole, mMax];
   survival = Total[Abs[psi]^2];
   <|"t" -> t, "OmegaT" -> Omega t, "survival" -> survival, 
    "fullWigner" -> fullWigner, "poleWigner" -> poleWigner, 
    "fullMax" -> Max[Abs[fullWigner[[All, 3]]]], 
    "poleMax" -> Max[Abs[poleWigner[[All, 3]]]]|>];

physicalSnapshots = MapThread[physicalSnapshotAt, {times, mMaxList}];


(*------------------------------------------------------------*)
(* 7. CONDITIONALLY NORMALIZED WIGNER SNAPSHOTS               *)
(*------------------------------------------------------------*)
(*Condition on the walker not having been absorbed:*)
(**)
(*        |psi_cond(t)> = |psi(t)> /Sqrt[S(t)].*)
(**)
(*Because the Wigner function is bilinear in the wavefunction,*)
(*this gives exactly*)
(**)
(*W_cond(m,k,t)=W_+(m,k,t)/S(t).*)
(**)
(* These snapshots remove the overall loss of probability and *)
(*show only the evolving phase-space structure of the*)
(*surviving walker.*)
(*------------------------------------------------------------*)

normalizedSnapshotAt[t_?NumericQ, mMax_Integer] := 
  Module[{psi, survival, psiConditional, conditionalWigner}, 
   psi = N[psiFullAt[t]];
   survival = Total[Abs[psi]^2];
   (*Normalize the surviving component to unit norm.*)
   psiConditional = psi/Sqrt[survival];
   conditionalWigner = wignerFromPsi[psiConditional, mMax];
   <|"t" -> t, "OmegaT" -> Omega t, "wigner" -> conditionalWigner, 
    "max" -> Max[Abs[conditionalWigner[[All, 3]]]]|>];

normalizedSnapshots = 
  MapThread[normalizedSnapshotAt, {times, mMaxList}];

(* ============================================================ *)
(* 8. PHYSICAL WIGNER FUNCTION                                  *)
(* Full Wigner function with localized-pole contours.            *)
(* ============================================================ *)

(*------------------------------------------------------------*)
(* 8.1 Global color scale                                     *)
(*------------------------------------------------------------*)
(* Use the same Wigner scale in every time panel so that the   *)
(* colors are quantitatively comparable across times.          *)
globalWMax = Max[physicalSnapshots[[All, "fullMax"]]];


(*------------------------------------------------------------*)
(*8.2 Wigner color map*)
(*------------------------------------------------------------*)
(*Negative values are blue, zero is white,and positive values*)
(*are red.The same map is used in every panel.*)

wignerColor = (Blend[{RGBColor[0.1, 0.2, 0.8], White, 
      RGBColor[0.85, 0.05, 0.05]}, 
     Rescale[Clip[#, {-globalWMax, globalWMax}], {-globalWMax, 
       globalWMax}]] &);


(*------------------------------------------------------------*)
(*8.3 Pole-contour visualization*)
(*------------------------------------------------------------*)
(*The localized pole contribution becomes much smaller than*)
(*the full Wigner function at later times. To keep its spatial*)
(*profile visible,only the contour overlay is rescaled.*)
(*The colored physical Wigner function is never rescaled.*)

poleContourScale = 0.35;


rescaledPoleData[j_Integer] := 
  Module[{poleData, poleMax, target}, 
   poleData = physicalSnapshots[[j, "poleWigner"]];
   poleMax = physicalSnapshots[[j, "poleMax"]];
   target = poleContourScale*globalWMax;
   If[poleMax == 0, poleData, 
    poleData /. {x_, k_, w_} :> {x, k, target Abs[w]/poleMax}]];


(*Four equally spaced contour strengths within the chosen visualization scale.*)
poleContourLevels = 
  poleContourScale*globalWMax*{0.20, 0.40, 0.60, 0.80};


(*------------------------------------------------------------*)(*8.4 Physical-Wigner-specific plot options*)(*------------------------------------------------------------*)
physicalWignerOptions = {ColorFunction -> wignerColor, 
   ColorFunctionScaling -> False, 
   PlotRange -> {-globalWMax, globalWMax}};


(*------------------------------------------------------------*)
(*8.5 Build each Wigner panel*)
(*------------------------------------------------------------*)

physicalWignerPanel[j_Integer] := 
  Module[{basePlot, poleContours}, 
   basePlot = 
    ListDensityPlot[physicalSnapshots[[j, "fullWigner"]], 
     Evaluate[Sequence @@ commonPlotStyle], 
     Evaluate[Sequence @@ physicalWignerOptions], 
     PlotLabel -> timeLabel[physicalSnapshots[[j, "OmegaT"]]]];
   poleContours = 
    ListContourPlot[rescaledPoleData[j], 
     Contours -> poleContourLevels, ContourShading -> False, 
     ContourStyle -> Directive[Black, Thick], InterpolationOrder -> 2,
      MaxPlotPoints -> Infinity, PlotRange -> All, Frame -> False, 
     Axes -> False, BoundaryStyle -> None];
   Show[basePlot, poleContours]];

wignerPlots = Table[physicalWignerPanel[j], {j, Length[physicalSnapshots]}];


(*------------------------------------------------------------*)
(*8.6 Shared color legend*)
(*------------------------------------------------------------*)

wignerBar = 
  BarLegend[{wignerColor, {-globalWMax, globalWMax}}, 
   LegendLabel -> None, LabelStyle -> legendStyle, 
   LegendMarkerSize -> {18, barHeight}];


(*------------------------------------------------------------*)
(*8.7 Assemble final figure*)
(*------------------------------------------------------------*)

physicalWignerFigure = 
  Grid[{{GraphicsRow[wignerPlots, Spacings -> 0.1, 
      ImageSize -> 4 plotWidth], wignerBar}}, 
   Alignment -> {Center, Center}, Spacings -> {0.4, 0}];

(* ============================================================ *)
(* 9. CONDITIONALLY NORMALIZED WIGNER FUNCTION                   *)
(* Wtilde_+(m,k,t) = W_+(m,k,t)/S(t)                             *)
(* ============================================================ *)
(* The state has already been normalized according to           *)
(*                                                               *)
(*   |psi_cond(t)> = |psi(t)>/Sqrt[S(t)],                        *)
(*                                                               *)
(* so, because the Wigner function is bilinear in psi,           *)
(*                                                               *)
(*   W_cond(m,k,t) = W_+(m,k,t)/S(t).                            *)
(*                                                               *)
(* These panels show the phase-space structure of the walker     *)
(* conditioned on nonabsorption.                                 *)

(*------------------------------------------------------------*)
(* 9.1 Global normalized-Wigner color scale                   *)
(*------------------------------------------------------------*)
(* Use one common scale for all four normalized snapshots so  *)
(* that the colors remain quantitatively comparable in time.   *)
globalWMaxNormalized = Max[normalizedSnapshots[[All, "max"]]];


(*------------------------------------------------------------*)
(*9.2 Normalized-Wigner color map*)
(*------------------------------------------------------------*)
(*Use the same blue-white-red convention as in the physical*)
(*Wigner figure:blue=negative,white=zero,red=positive.*)

normalizedWignerColor = (Blend[{RGBColor[0.1, 0.2, 0.8], White, 
      RGBColor[0.85, 0.05, 0.05]}, 
     Rescale[Clip[#, {-globalWMaxNormalized, 
        globalWMaxNormalized}], {-globalWMaxNormalized, 
       globalWMaxNormalized}]] &);


(*------------------------------------------------------------*)
(*9.3 Normalized-Wigner-specific plot options*)
(*------------------------------------------------------------*)
(*Only the color function and vertical color scale differ from*)
(*the physical Wigner figure. All typography and dimensions*)
(*are inherited from commonPlotStyle.*)

normalizedWignerOptions = {ColorFunction -> normalizedWignerColor, 
   ColorFunctionScaling -> False, 
   PlotRange -> {-globalWMaxNormalized, globalWMaxNormalized}};


(*------------------------------------------------------------*)
(*9.4 Build each normalized Wigner panel*)
(*------------------------------------------------------------*)

normalizedWignerPanel[j_Integer] := 
  ListDensityPlot[normalizedSnapshots[[j, "wigner"]], 
   Evaluate[Sequence @@ commonPlotStyle], 
   Evaluate[Sequence @@ normalizedWignerOptions], 
   PlotLabel -> timeLabel[normalizedSnapshots[[j, "OmegaT"]]]];


normalizedWignerPlots = 
  Table[normalizedWignerPanel[j], {j, Length[normalizedSnapshots]}];


(*------------------------------------------------------------*)
(*9.5 Shared normalized-Wigner color legend*)
(*------------------------------------------------------------*)
(*A single color bar is valid because all four panels use the*)
(*same global normalized-Wigner scale.*)

normalizedWignerBar = 
  BarLegend[{normalizedWignerColor, {-globalWMaxNormalized, 
     globalWMaxNormalized}}, LegendLabel -> None, 
   LabelStyle -> legendStyle, LegendMarkerSize -> {18, barHeight}];


(*------------------------------------------------------------*)
(*9.6 Assemble the final normalized-Wigner figure*)
(*------------------------------------------------------------*)

normalizedWignerFigure = 
  Grid[{{GraphicsRow[normalizedWignerPlots, Spacings -> 0.1, 
      ImageSize -> 4 plotWidth], normalizedWignerBar}}, 
   Alignment -> {Center, Center}, Spacings -> {0.4, 0}];

(* ------------------------------------------------------------ *)
(* 10. DISPLAY AND OPTIONAL EXPORT                              *)
(* ------------------------------------------------------------ *)
(* Display the two manuscript figures.                          *)

Print[physicalWignerFigure];
Print[normalizedWignerFigure];

(* Optional export to the current working directory:
Export["wigner_physical.pdf", physicalWignerFigure];
Export["wigner_conditionally_normalized.pdf", normalizedWignerFigure];
*)

(* ============================================================ *)
(* 11. NUMERICAL CONVERGENCE CHECKS                             *)
(* ============================================================ *)
(* The physical problem is defined on a semi-infinite lattice.  *)
(* Numerically, the lattice is truncated at Nsites = 80.        *)
(* The checks below verify that the displayed Wigner functions  *)
(* are insensitive to this artificial right boundary.           *)
(*                                                              *)
(* We compare the production cutoff N = 80 with N = 120 and    *)
(* N = 160, keeping all physical and plotting parameters fixed. *)
(* For each snapshot, we report the largest pointwise difference*)
(* over the complete displayed (m,k) domain.                    *)
(*                                                              *)
(* These checks do not test an infinite Wigner sum: for each m, *)
(* the Wigner definition contains exactly m-1 terms.            *)
(* ============================================================ *)

ClearAll[
  convergenceHamiltonian,
  convergencePsi,
  convergenceWignerValue,
  convergenceConditionalPsi,
  maxWignerCutoffDifference,
  maxConditionalWignerCutoffDifference
  ];

(*------------------------------------------------------------*)
(*11.1 Finite-chain Hamiltonian for arbitrary cutoff*)
(*------------------------------------------------------------*)

convergenceHamiltonian[nSites_Integer] :=
 Module[{H},
  H = ConstantArray[0. + 0. I, {nSites, nSites}];

  Do[H[[s, s]] = Omega, {s, 1, nSites}];

  Do[
   H[[s, s + 1]] = -Omega/2;
   H[[s + 1, s]] = -Omega/2,
   {s, 1, nSites - 1}
   ];

  H[[1, 1]] -= I kappa/2;

  SparseArray[H]
  ];

(*------------------------------------------------------------*)
(*11.2 Surviving wavefunction for arbitrary cutoff*)
(*------------------------------------------------------------*)

convergencePsi[t_?NumericQ, nSites_Integer] :=
 Module[{psi0},
  psi0 = UnitVector[nSites, s0];
  N[MatrixExp[-I convergenceHamiltonian[nSites] t] . psi0]
  ];

(*------------------------------------------------------------*)
(*11.3 Wigner function at one (m,k) point*)
(*------------------------------------------------------------*)

convergenceWignerValue[psi_List, m_Integer, k_?NumericQ] :=
 Re[
  1/(2 Pi)
   Sum[
    psi[[n]]
     Conjugate[psi[[m - n]]]
     Exp[-I (2 n - m) k],
    {n, 1, m - 1}
    ]
  ];

(*------------------------------------------------------------*)
(*11.4 Physical Wigner: lattice-cutoff convergence*)
(*------------------------------------------------------------*)

kConvergenceGrid = N@Subdivide[-Pi, Pi, nk];

maxWignerCutoffDifference[
   t_?NumericQ,
   mMax_Integer,
   nSites1_Integer,
   nSites2_Integer
   ] :=
 Module[{psi1, psi2},

  psi1 = convergencePsi[t, nSites1];
  psi2 = convergencePsi[t, nSites2];

  Max[
   Flatten[
    Table[
     Abs[
      convergenceWignerValue[psi1, m, k] -
       convergenceWignerValue[psi2, m, k]
      ],
     {m, 2, mMax},
     {k, kConvergenceGrid}
     ]
    ]
   ]
  ];

physicalWignerConvergence =
 Table[
  {
   times[[j]],
   mMaxList[[j]],
   maxWignerCutoffDifference[times[[j]], mMaxList[[j]], 80, 160],
   maxWignerCutoffDifference[times[[j]], mMaxList[[j]], 120, 160]
   },
  {j, Length[times]}
  ];

Print[""];
Print["Physical Wigner-function lattice-cutoff convergence:"];
Print[
 Grid[
  Prepend[
   Map[
    {#[[1]], #[[2]], ScientificForm[#[[3]], 6],
      ScientificForm[#[[4]], 6]} &,
    physicalWignerConvergence
    ],
   {"\[CapitalOmega] t", "mMax", "max |W80 - W160|",
    "max |W120 - W160|"}
   ],
  Frame -> All,
  Alignment -> Center,
  Spacings -> {2, 1}
  ]
 ];

(*------------------------------------------------------------*)
(*11.5 Conditionally normalized Wigner: cutoff convergence*)
(*------------------------------------------------------------*)

convergenceConditionalPsi[t_?NumericQ, nSites_Integer] :=
 Module[{psi, survival},
  psi = convergencePsi[t, nSites];
  survival = Total[Abs[psi]^2];
  psi/Sqrt[survival]
  ];

maxConditionalWignerCutoffDifference[
   t_?NumericQ,
   mMax_Integer,
   nSites1_Integer,
   nSites2_Integer
   ] :=
 Module[{psi1, psi2},

  psi1 = convergenceConditionalPsi[t, nSites1];
  psi2 = convergenceConditionalPsi[t, nSites2];

  Max[
   Flatten[
    Table[
     Abs[
      convergenceWignerValue[psi1, m, k] -
       convergenceWignerValue[psi2, m, k]
      ],
     {m, 2, mMax},
     {k, kConvergenceGrid}
     ]
    ]
   ]
  ];

normalizedWignerConvergence =
 Table[
  {
   times[[j]],
   mMaxList[[j]],
   maxConditionalWignerCutoffDifference[
    times[[j]], mMaxList[[j]], 80, 160],
   maxConditionalWignerCutoffDifference[
    times[[j]], mMaxList[[j]], 120, 160]
   },
  {j, Length[times]}
  ];

Print[""];
Print["Conditionally normalized Wigner-function lattice-cutoff convergence:"];
Print[
 Grid[
  Prepend[
   Map[
    {#[[1]], #[[2]], ScientificForm[#[[3]], 6],
      ScientificForm[#[[4]], 6]} &,
    normalizedWignerConvergence
    ],
   {"\[CapitalOmega] t", "mMax",
    "max |Wtilde80 - Wtilde160|",
    "max |Wtilde120 - Wtilde160|"}
   ],
  Frame -> All,
  Alignment -> Center,
  Spacings -> {2, 1}
  ]
 ];

(*------------------------------------------------------------*)
(*11.6 Population near the artificial right boundary*)
(*------------------------------------------------------------*)

rightBoundarySites = 10;

maxRightBoundaryPopulation[nSites_Integer] :=
 Max[
  Table[
   Module[{psi},
    psi = convergencePsi[t, nSites];
    Total[Abs[Take[psi, -rightBoundarySites]]^2]
    ],
   {t, times}
   ]
  ];

rightBoundaryPopulation80 = maxRightBoundaryPopulation[80];

Print[""];
Print[
 "Maximum population in the last ",
 rightBoundarySites,
 " sites of the Nsites = 80 chain over the four plotted times: ",
 ScientificForm[rightBoundaryPopulation80, 6]
 ];

(*------------------------------------------------------------*)
(* 11.7 Ballistic-front sanity check                           *)
(*------------------------------------------------------------*)
(* For the nearest-neighbor dispersion E(k)=Omega(1-cos k),    *)
(* the maximum group velocity is v_max = Omega. A conservative *)
(* estimate of the right-moving ballistic front is therefore   *)
(* s0 + Omega t. This is not used as a numerical approximation; *)
(* it is an independent sanity check on the finite-chain cutoff. *)

maxGroupVelocity = Omega;
maxPlottedTime = Max[times];
ballisticFrontEstimate = s0 + maxGroupVelocity maxPlottedTime;
rightBoundaryDistance = Nsites - s0;
ballisticSafetyMargin = Nsites - ballisticFrontEstimate;

Print[""];
Print["Ballistic-front sanity check:"];
Print[
 Grid[
  {
   {"maximum group velocity", maxGroupVelocity},
   {"maximum plotted time", maxPlottedTime},
   {"estimated right-moving front", ballisticFrontEstimate},
   {"production lattice cutoff Nsites", Nsites},
   {"distance from initial site to right boundary", rightBoundaryDistance},
   {"boundary minus estimated front", ballisticSafetyMargin}
  },
  Frame -> All,
  Alignment -> Left
 ]
];

(*------------------------------------------------------------*)
(* 11.8 Momentum-grid resolution check                        *)
(*------------------------------------------------------------*)
(* The Wigner sum is exact and finite for every (m,k). The k   *)
(* grid only samples this continuous function for plotting.     *)
(* We test nk=900 by comparing exact Wigner values at the       *)
(* midpoints of the coarse grid with linear interpolation from  *)
(* the neighboring coarse-grid values. The error is reported    *)
(* over the complete displayed m range for every snapshot.      *)

momentumResolutionError[t_?NumericQ, mMax_Integer, nIntervals_Integer] :=
 Module[{psi, dk, kLeft, kMid},
  psi = convergencePsi[t, Nsites];
  dk = 2 Pi/nIntervals;
  kLeft = N@Table[-Pi + j dk, {j, 0, nIntervals - 1}];
  kMid = kLeft + dk/2;
  Max[
   Flatten[
    Table[
     Abs[
      convergenceWignerValue[psi, m, kMid[[j + 1]]] -
       1/2 (
        convergenceWignerValue[psi, m, kLeft[[j + 1]]] +
        convergenceWignerValue[psi, m, kLeft[[j + 1]] + dk]
       )
     ],
     {m, 2, mMax},
     {j, 0, nIntervals - 1}
    ]
   ]
  ]
 ];

momentumGridConvergence =
 Table[
  {
   times[[j]],
   mMaxList[[j]],
   momentumResolutionError[times[[j]], mMaxList[[j]], 900],
   momentumResolutionError[times[[j]], mMaxList[[j]], 1800]
  },
  {j, Length[times]}
 ];

Print[""];
Print["Momentum-grid resolution check:"];
Print[
 Grid[
  Prepend[
   Map[
    {#[[1]], #[[2]], ScientificForm[#[[3]], 6],
      ScientificForm[#[[4]], 6]} &,
    momentumGridConvergence
   ],
   {"Omega t", "mMax", "midpoint interpolation error (nk=900)",
    "midpoint interpolation error (nk=1800)"}
  ],
  Frame -> All,
  Alignment -> Center,
  Spacings -> {2, 1}
 ]
];

(* ============================================================ *)
(* 12. CONVERGENCE SUMMARY                                      *)
(* ============================================================ *)

convergenceSummary = <|
  "PhysicalWigner" -> physicalWignerConvergence,
  "NormalizedWigner" -> normalizedWignerConvergence,
  "RightBoundaryPopulation" -> rightBoundaryPopulation80,
  "BallisticFrontEstimate" -> ballisticFrontEstimate,
  "BallisticSafetyMargin" -> ballisticSafetyMargin,
  "MomentumGrid" -> momentumGridConvergence
|>;

convergenceSummary
