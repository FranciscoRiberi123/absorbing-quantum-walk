(* ============================================================ *)
(* Absorbing continuous-time quantum walk                       *)
(* Figure 1: survival probability and first-passage density     *)
(*                                                              *)
(* This script generates the survival-probability and           *)
(* first-passage-density panels used in Fig. 1 of the           *)
(* accompanying manuscript.                                     *)
(*                                                              *)
(* Wolfram Language / Mathematica                               *)
(* No external packages are required.                           *)
(* ============================================================ *)

ClearAll["Global`*"];

(* ------------------------------------------------------------ *)
(* 1. Model and numerical parameters                            *)
(* ------------------------------------------------------------ *)
(* The parameters below reproduce Fig. 1.                       *)
(*                                                              *)
(* Omega  : coherent hopping rate                               *)
(* s0     : initial lattice site                                *)
(* eta    : dimensionless absorption strength, eta=kappa/Omega  *)
(* tMax   : largest plotted time                                *)
(* nT     : number of time intervals                            *)
(* Nsites : finite-chain cutoff used to approximate the         *)
(*          semi-infinite lattice                               *)
(* ------------------------------------------------------------ *)

Omega = 1.;
s0 = 8;

etas = {0.25, 1.0, 4.0};
etaLabels = {"\[Eta] = 0.25", "\[Eta] = 1", "\[Eta] = 4"};

tMax = 30.;
nT = 1800;

(* Uniform time grid containing nT+1 points. *)
times = Subdivide[0., tMax, nT];

(* The chain is chosen long enough that the wavepacket does not
   reach the artificial right boundary during the plotted time. *)
Nsites = 220;

(* ------------------------------------------------------------ *)
(* 2. Plot style                                                *)
(* ------------------------------------------------------------ *)
(* All three eta values use the same line thickness but         *)
(* different colors and dash patterns.                          *)
(* ------------------------------------------------------------ *)

curveStyles = {
   Directive[Red, AbsoluteThickness[3.2], Dashing[{0.035, 0.018}]],
   Directive[Blue, AbsoluteThickness[3.2]],
   Directive[Black, AbsoluteThickness[3.2], DotDashed]
   };

(* Common formatting for both panels. *)
commonPlotOptions = {
   Frame -> True,
   FrameStyle -> Black,
   LabelStyle -> {FontSize -> 24, FontColor -> Black, FontFamily -> "Times"},
   GridLines -> Automatic,
   ImageSize -> 500,
   ImagePadding -> {{82, 30}, {70, 38}},
   PlotRangePadding -> Scaled[0.03]
   };

(* ------------------------------------------------------------ *)
(* 3. Effective non-Hermitian Hamiltonian                       *)
(* ------------------------------------------------------------ *)
(* For a given eta, kappa = eta Omega and                       *)
(*                                                              *)
(*   H_eff = H_0 - i kappa/2 |1><1|,                            *)
(*                                                              *)
(* where H_0 is the nearest-neighbor tight-binding Hamiltonian. *)
(* The imaginary term removes probability from the absorbing    *)
(* boundary site s=1.                                           *)
(* ------------------------------------------------------------ *)

heff[eta_?NumericQ] :=
 Module[{H, kappa},
  kappa = eta*Omega;

  (* Start with an Nsites x Nsites complex matrix. *)
  H = ConstantArray[0. + 0. I, {Nsites, Nsites}];

  (* On-site energy Omega. *)
  Do[
   H[[s, s]] = Omega,
   {s, 1, Nsites}
   ];

  (* Nearest-neighbor hopping -Omega/2. *)
  Do[
   H[[s, s + 1]] = -Omega/2;
   H[[s + 1, s]] = -Omega/2,
   {s, 1, Nsites - 1}
   ];

  (* Absorbing boundary at site s=1. *)
  H[[1, 1]] -= I*kappa/2;

  (* Sparse storage is more efficient for this tridiagonal matrix. *)
  SparseArray[H]
  ];

(* ------------------------------------------------------------ *)
(* 4. Compute S(t) and F(t) for one absorption strength         *)
(* ------------------------------------------------------------ *)
(* The initial state is localized at site s0.                   *)
(*                                                              *)
(* For the surviving state                                      *)
(*                                                              *)
(*   |psi(t)> = exp[-i H_eff t] |s0>,                           *)
(*                                                              *)
(* the survival probability is                                  *)
(*                                                              *)
(*   S(t) = Sum_s |psi_s(t)|^2.                                 *)
(*                                                              *)
(* Probability leaves the lattice only through site s=1, so     *)
(* the first-passage density is                                 *)
(*                                                              *)
(*   F(t) = kappa |psi_1(t)|^2 = -dS(t)/dt.                     *)
(*                                                              *)
(* The function returns an Association containing plotting data *)
(* for both observables.                                        *)
(* ------------------------------------------------------------ *)

makeData[eta_?NumericQ] :=
 Module[{H, psi0, states, survivalValues, firstPassageValues, kappa},
  kappa = eta*Omega;
  H = heff[eta];

  (* Initial walker localized at site s0. *)
  psi0 = UnitVector[Nsites, s0];

  (* Surviving wavefunction at every time in the grid. *)
  states = (MatrixExp[-I*H*#] . psi0) & /@ times;

  (* Survival probability S(t). *)
  survivalValues = Total[Abs[#]^2] & /@ states;

  (* First-passage density F(t)=kappa |psi_1(t)|^2. *)
  firstPassageValues = kappa*Abs[#[[1]]]^2 & /@ states;

  <|
   "eta" -> eta,
   "Survival" -> Transpose[{times, survivalValues}],
   "FirstPassage" -> Transpose[{times, firstPassageValues}]
   |>
  ];

(* Evaluate the dynamics for all three absorption strengths. *)
data = makeData /@ etas;

(* ------------------------------------------------------------ *)
(* 5. Probability-conservation check                            *)
(* ------------------------------------------------------------ *)
(* The exact dynamics obey                                      *)
(*                                                              *)
(*   S(t) + Integral_0^t F(tau) d tau = 1.                      *)
(*                                                              *)
(* We verify this numerically at the final plotted time using   *)
(* the trapezoidal rule.                                        *)
(* ------------------------------------------------------------ *)

conservationCheck[j_Integer] :=
 Module[{survivalValues, firstPassageValues, absorbedProbability, finalSurvival},
  survivalValues = data[[j, "Survival"]][[All, 2]];
  firstPassageValues = data[[j, "FirstPassage"]][[All, 2]];

  (* Trapezoidal approximation to Integral F(t) dt. *)
  absorbedProbability =
   Total[
    Differences[times]*
     MovingAverage[firstPassageValues, 2]
    ];

  finalSurvival = Last[survivalValues];

  <|
   "eta" -> etas[[j]],
   "S(tMax)" -> finalSurvival,
   "IntegralF" -> absorbedProbability,
   "Total" -> finalSurvival + absorbedProbability
   |>
  ];

conservationResults =
 Table[
  conservationCheck[j],
  {j, Length[etas]}
  ];

conservationResults

(* ------------------------------------------------------------ *)
(* 6. Automatic plot ranges                                     *)
(* ------------------------------------------------------------ *)
(* Determine ranges from all three curves so no data are clipped *)
(* when parameters are modified.                                *)
(* ------------------------------------------------------------ *)

allSurvivalValues =
 Flatten[
  data[[All, "Survival"]][[All, All, 2]]
  ];

allFirstPassageValues =
 Flatten[
  data[[All, "FirstPassage"]][[All, All, 2]]
  ];

(* Small visual margins around the data. *)
survivalPadding =
 0.05*(Max[allSurvivalValues] - Min[allSurvivalValues]);

firstPassagePadding =
 0.08*(Max[allFirstPassageValues] - Min[allFirstPassageValues]);

survivalMin =
 Max[0, Min[allSurvivalValues] - survivalPadding];

survivalMax =
 Min[1.05, Max[allSurvivalValues] + survivalPadding];

firstPassageMin =
 Max[0, Min[allFirstPassageValues] - firstPassagePadding];

firstPassageMax =
 Max[allFirstPassageValues] + firstPassagePadding;

(* ------------------------------------------------------------ *)
(* 7. Legend                                                    *)
(* ------------------------------------------------------------ *)
(* The custom legend reproduces the line colors and dash styles *)
(* used in the manuscript figure.                               *)
(* ------------------------------------------------------------ *)

legendBox =
 Graphics[
  {
   EdgeForm[Black],
   FaceForm[White],
   Rectangle[{0, 0}, {1, 1}],

   curveStyles[[1]],
   Line[{{0.08, 0.78}, {0.22, 0.78}}],
   Black,
   Text[
    Style[etaLabels[[1]], 24, FontFamily -> "Times"],
    {0.62, 0.78}
    ],

   curveStyles[[2]],
   Line[{{0.08, 0.50}, {0.22, 0.50}}],
   Black,
   Text[
    Style[etaLabels[[2]], 24, FontFamily -> "Times"],
    {0.62, 0.50}
    ],

   curveStyles[[3]],
   Line[{{0.08, 0.22}, {0.22, 0.22}}],
   Black,
   Text[
    Style[etaLabels[[3]], 24, FontFamily -> "Times"],
    {0.62, 0.22}
    ]
   },
  PlotRange -> {{0, 1}, {0, 1}},
  ImageSize -> {300, 120},
  Background -> None,
  ImagePadding -> 0
  ];

(* ------------------------------------------------------------ *)
(* 8. Survival-probability panel                                *)
(* ------------------------------------------------------------ *)

survivalPlot =
 ListLinePlot[
  data[[All, "Survival"]],
  PlotStyle -> curveStyles,
  Joined -> True,
  Evaluate[Sequence @@ commonPlotOptions],
  PlotRange -> {{0, tMax}, {survivalMin, survivalMax}},
  FrameLabel -> {
    Style["\[CapitalOmega] t", Italic, Black, 34, FontFamily -> "Times"],
    Style["S(t)", Italic, Black, 34, FontFamily -> "Times"]
    },
  PlotLabel ->
   Style["Survival probability", Black, 22, FontFamily -> "Times"]
  ];

(* ------------------------------------------------------------ *)
(* 9. First-passage-density panel                               *)
(* ------------------------------------------------------------ *)

firstPassagePlot =
 ListLinePlot[
  data[[All, "FirstPassage"]],
  PlotStyle -> curveStyles,
  Joined -> True,
  Evaluate[Sequence @@ commonPlotOptions],
  PlotRange -> {{0, tMax}, {firstPassageMin, firstPassageMax}},
  (* Place the common legend inside the right-hand panel. *)
  Epilog ->
   Inset[
    legendBox,
    Scaled[{0.84, 0.70}]
    ],
  FrameLabel -> {
    Style["\[CapitalOmega] t", Italic, Black, 34, FontFamily -> "Times"],
    Style["F(t)", Italic, Black, 34, FontFamily -> "Times"]
    },
  PlotLabel ->
   Style["First-passage density", Black, 22, FontFamily -> "Times"]
  ];

(* ------------------------------------------------------------ *)
(* 10. Assemble Figure 1                                        *)
(* ------------------------------------------------------------ *)

figure1 =
 GraphicsRow[
  {survivalPlot, firstPassagePlot},
  Spacings -> 0.18,
  ImageSize -> 1100
  ];

figure1

(* ------------------------------------------------------------ *)
(* Optional export                                              *)
(* ------------------------------------------------------------ *)

(*
Export["Figure1.pdf", figure1];
*)
