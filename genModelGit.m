(* ::Package:: *)

(* Simulates gene drive dynamics in multiple populations over time, similar to an SIR model for disease spread*)


(* ::Input::Initialization:: *)
Clear[t];
Clear[q1,q2,q3];
Clear[s];
Clear[c];
Clear[h];
Clear[b];
Clear[m];

dynamics[q10_, q20_, q30_, s_, c_, h_, m_,time_]:=
Module[{q1,q2,q3,list1,list2,list3,q1t,q2t,q3t,sn,sc},
sn = 1/2 (1-c)(1-h*s);
sc = c(1-s);
q1 = q10;
q2 =q20;
q3 = q30;
list1 = {q1};
list2 = {q2};
list3 = {q3};
Do[
q1t = (q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q1+m*q2)};
q2t = (q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-(2*m))q2+m*q1+m*q3)};
q3t = (q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q3+m*q2)};
q1 = q1t; q2 = q2t; q3 = q3t;
AppendTo[list1,Round[q1, 10^-4]];
AppendTo[list2,Round[q2, 10^-4]];
AppendTo[list3,Round[q3, 10^-4]];
, time];
{list1,list2,list3}
]

dynamics2[q10_, q20_, s_, c_, h_, m_,time_]:=
Module[{q1,q2,list1,list2,q1t,q2t,sn,sc},
sn = 1/2 (1-c)(1-h*s);
sc = c(1-s);
q1 = q10;
q2 =q20;
list1 = {q1};
list2 = {q2};
Do[
q1t = (q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q1+m*q2)};
q2t = (q^2 (1-s)+2q(1-q)(sn+sc))/(q^2 (1-s)+2q(1-q)(2sn +sc)+(1-q)^2)/.{q->((1-m)q2+m*q1)};
q1 = q1t; q2 = q2t;
AppendTo[list1,Round[q1, 10^-4]];
AppendTo[list2,Round[q2, 10^-4]];
, time];
{list1,list2}
]


m1 =Manipulate[Animate[ListPlot[dynamics[0,q2,0,s,1,0,m,t], PlotRange->{{0,102},{0,1.01}},PlotStyle->{Blue,Red,Green}, AxesLabel->{"Time","Allele Frequency"}, PlotLegends->{"Deme 1", "Deme 2", "Deme 3"}, PlotLabel->"Central Targeting"],{t, 1,100}], {q2, 0, 1, Appearance->"Labeled"}, {s,0.5,1, Appearance->"Labeled"}, {m, 0, 0.15, Appearance->"Labeled"}];
(*Export["Yonatan/centralAnimated.avi", m1];*)
m1

Manipulate[Animate[ListPlot[dynamics[q1,0,0,s,1,0,m,t], PlotRange->{{0,102},{0,1.01}},PlotStyle->{Blue,Red,Green}, AxesLabel->{"Time","Allele Frequency"},PlotLegends->{"Deme 1", "Deme 2", "Deme 3"},PlotLabel->"Peripheral Targeting"],{t, 1,100}], {q1, 0, 1, Appearance->"Labeled"}, {s,0.5,1, Appearance->"Labeled"}, {m, 0, 0.15, Appearance->"Labeled"}]
m2= Animate[ListPlot[dynamics[0.752,0,0,0.714,1,0,0.0636,t], PlotRange->{{0,102},{0,1.01}},PlotStyle->{Blue,Red,Green}, AxesLabel->{"Time","Allele Frequency"},PlotLegends->{"Deme 1", "Deme 2", "Deme 3"},PlotLabel->"Peripheral Targeting"],{t, 1,100}, AnimationDirection->Forward, AppearanceElements->None];
(*Export["Yonatan/periphSpillover.gif", m2];*)
m2


m1 =Manipulate[Animate[ListPlot[dynamics2[0,q2,s,1,0,m,t], PlotRange->{{0,102},{0,1.01}},PlotStyle->{Blue,Red}, AxesLabel->{"Time","Allele Frequency"}, PlotLegends->{"Deme 1", "Deme 2"}, PlotLabel->"Two Demes"],{t, 1,100}], {q1, 0, 1, Appearance->"Labeled"}, {s,0.5,1, Appearance->"Labeled"}, {m, 0, 0.15, Appearance->"Labeled"}];
(*Export["Yonatan/centralAnimated.avi", m1];*)
m1

m2= Animate[ListPlot[dynamics2[0.752,0,0.714,1,0,0.0636,t], PlotRange->{{0,102},{0,1.01}},PlotStyle->{Blue,Red}, AxesLabel->{"Time","Allele Frequency"},PlotLegends->{"Deme 1", "Deme 2"},PlotLabel->"Two Demes"],{t, 1,100}, AnimationDirection->Forward, AppearanceElements->None];
(*Export["Yonatan/periphSpillover.gif", m2];*)
m2


(* ::InheritFromParent:: *)
(**)
