
# News/Changelog 

Changes in **MrDAG** package version **1.1.0**:
<ul>
    <li> In MrDAG(), the option 'fastMCMC' has been removed. Now it is automatically executed
    <li> In MrDAG(), 'niter', the number of Markov chain Monte Carlo iterations, includes the burn-in phase
    <li> In MrDAG(), the names of some input parameters have changed ('w' -> 'pp', 'tempmax' -> 'tempMax', 'filename' -> 'fileName', and 'filepath' -> 'filePath')
    <li> In MrDAG(), the names of some output have changed ('logmarglik' -> 'logMargLik', 'validpropMrDAG' -> 'validPropMrDAG', 'acceptpropDAG' -> 'acceptPropDAG', 'hyperpar' -> 'hyperPar', and 'samplerpar' -> 'samplerPar')
    <li> In get_causaleffect(), a better description of the output is now provided
    <li> In get_causaleffect(), the names of some output have changed ('causaleffect' -> 'causalEffect', 'causaleffect_LL' -> 'causalEffect_LL', and 'causaleffect_LL -> 'causalEffect_UL)
    <li> In get_causaleffects(), a bug has been corrected when 'ord' is specified
    <li> In get_causaleffects(), the names of some output have changed ('causaleffects' -> 'causalEffects', 'causaleffects_LL' -> 'causalEffects_LL', and 'causaleffects_LL -> 'causalEffects_UL)
</ul>
