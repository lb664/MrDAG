
# News/Changelog 

Changes in **MrDAG** package version **0.1.3** (25/10/2026):
<ul>
    <li> Documentation updated: 'genetically predicted effects' replaces 'causal effects' in the text, PEIP (posterior edge inclusion probability) replaces PPEI, and the references are completed
</ul> 

Changes in **MrDAG** package version **0.1.2** (01/07/2025):
<ul>
    <li> Typos in the functions' descriptions have been corrected. Clarifications of some functions' outputs are also included
</ul> 

Changes in **MrDAG** package version **0.1.1** (01/05/2025):
<ul>
    <li> In MrDAG(), the option 'fastMCMC' has been removed. Now it is automatically executed
    <li> In MrDAG(), 'niter', the number of Markov chain Monte Carlo iterations, includes the burn-in phase
    <li> In MrDAG(), the names of some input parameters have changed ('w' -> 'pp', 'tempmax' -> 'tempMax', 'filename' -> 'fileName', and 'filepath' -> 'filePath')
    <li> In MrDAG(), the names of some output have changed ('logmarglik' -> 'logMargLik', 'validpropMrDAG' -> 'validPropMrDAG', 'acceptpropDAG' -> 'acceptPropDAG', 'hyperpar' -> 'hyperPar', and 'samplerpar' -> 'samplerPar')
    <li> In get_causaleffect(), a better description of the output is now provided
    <li> In get_causaleffect(), the names of some output have changed ('causaleffect' -> 'causalEffect', 'causaleffect_LL' -> 'causalEffect_LL', and 'causaleffect_UL -> 'causalEffect_UL)
    <li> In get_causaleffects(), a bug has been corrected when 'ord' is specified
    <li> In get_causaleffects(), the names of some output have changed ('causaleffects' -> 'causalEffects', 'causaleffects_LL' -> 'causalEffects_LL', and 'causaleffects_UL -> 'causalEffects_UL)
</ul>
