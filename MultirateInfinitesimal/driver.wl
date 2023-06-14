(* ::Package:: *)

(* ::Title:: *)
(*MRI-GARK Analysis*)


(* ::Text:: *)
(*This script analyzes properties of an multirate infinitesimal generalized additive Runge\[Dash]Kutta (MRI-GARK) methods. It uses Integreat as a dependency.*)


PacletInstall["https://github.com/ComputationalScienceLaboratory/Integreat/releases/download/v0.4.0/Integreat.paclet"]


(* ::Text:: *)
(*Load the packages*)


SetDirectory[NotebookDirectory[]];
<<MRI`;
<<Integreat`RK`;


(* ::Text:: *)
(*View the usage notes for the MRI functions*)


?MRI*


mri = MRI["MRI-GARK-ERK45a"];
MRIDeltaC[mri]
MRIOrderConditions[mri, 4]
Norm[MRIOrderConditions[mri, {5}]]//N (* The principal error of the method *)
Limit[MRIStability[mri, zf, zs], zf->-\[Infinity]]


underlyingRK = RK[mri]
RKOrder[underlyingRK]
RKLinearStabilityPlot[underlyingRK]



