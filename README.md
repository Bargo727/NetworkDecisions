# NetworkDecisions
The code provided here was used to produce the figures in the paper "Optimal Evidence Accumulation on Networks"

The name of the file corresponds to the figure in the paper that was produced by the file.  Running the code will produce the figures.  All simulations here were done in MATLAB

fig1&2.m
This file produces evidence accumulation trajectories for a two-agent unidirectional network.  The script was used to produce both figures 1 and 2 in the paper.

fig3a.m
This file produces plots showing the difference in the percent of trials for which both agents of a unidirectional two-agent network get the correct answers.  It is figure 3a.

fig3b.m
This file produces the decision time densities for agent 2 in a unidirectional two-agent network when it receives social information and when it does not.  It is figure 3b.

fig3c.m
This file produces plots that show the percent difference in decision times of an agent that receives social information versus an agent that does not receive social information.  It is figure 3c. 

fig5.m
This file produces evidence accumulation trajectories for a recurrent two-agent network.  Equilibration is necessary for proper transfer of decision information.  It is shown as figure 5 in the paper.

fig7.m
This file produces plots showing the decision time for all three agents of a unidirectional network to make a decision and a plot showing what percentage of trials all three agents choose the correct decision.  It is figure 7 in the paper.

fig8b.m
This file produces plots showing the stochastic trajectories of three agents' evidence that are connected all-to-all.  It is figure 8b in the paper.

fig9.m
This file produces plots showing the decision time for all three agents of a clique to make a decision and a plot showing what percentage of trials all three agents choose the correct decision.  It is figure 9 in the paper.

fig10b-d.m
This file produces plots showing various statistics for cliques of size N.  It is for figure 10 b-d  in the paper.

fig10e.m
This file produces a plot showing the disparity in the fraction of a clique of size N that is correct at time of decision versus at the end of first equilibration.  It is figure 10e in the paper.

cliqueCodeLowSNRLargeThreshold.m
This file produces a plot like the Fig. 8.1 in the text.  The SNR is lower and the thresholds are higher.  One can easily modify these parameters to investigate how they alter key statistics for large cliques.
