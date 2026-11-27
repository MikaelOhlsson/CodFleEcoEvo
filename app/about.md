## About this app

This app explores the eco-evolutionary model of cod and flounder from the manuscript *Eco-evolutionary dynamics hinder recovery in harvested populations*. Cod mean body size evolves under selection from intrinsic growth, size-selective fishing, hypoxia, and competition with flounder. Fishing can push cod below the flounder's body size, after which competition from flounder may keep them small even once fishing stops: an eco-evolutionary trap.

The app has two tabs:

- **Model overview** shows where the equilibria are and how they depend on fishing effort, computed directly from the model equations.
- **Simulation** integrates the model over time, with fishing switched on during a chosen period, and animates the result.

All parameters are set in the sidebar; both tabs update as soon as a parameter changes. **Reset to default parameters** restores the defaults listed at the bottom of this page.


### The model

Cod (density *N*<sub>C</sub>, mean body size μ<sub>C</sub>) and flounder (density *N*<sub>F</sub>, mean body size μ<sub>F</sub>) have normally distributed body sizes with standard deviations σ<sub>C</sub> and σ<sub>F</sub>. Body sizes are scaled so that the flounder's mean is μ<sub>F</sub> = 0. Flounder do not evolve (their heritability is zero), so μ<sub>F</sub> stays at 0 throughout.

The per capita growth rate of cod is the sum of four terms (Eq. 3 of the manuscript):

- **Intrinsic growth**, highest when cod mean size equals the optimal size *z*\*: ϱ(θ<sup>2</sup> − (*z*\* − μ<sub>C</sub>)<sup>2</sup> − σ<sub>C</sub><sup>2</sup>)/θ<sup>2</sup>;
- **Fishing mortality**, a sigmoid function of body size with maximum η, inflection point φ and transition width τ;
- **Hypoxic mortality**, a sigmoid with maximum κ, inflection point ζ and transition width ν;
- **Competition**, −α<sub>CC</sub>*N*<sub>C</sub> − α<sub>CF</sub>*N*<sub>F</sub>. The competition coefficients decline with the difference in mean body size, with width ω (Eq. 5); α<sub>0</sub> sets their overall strength, and competition between species is further scaled by α<sub>I</sub>, the overlap in niche dimensions other than body size.

Cod mean body size changes as d*μ*<sub>C</sub>/d*t* = *h*<sup>2</sup><sub>C</sub> × (sum of the selection terms from these four processes) (Eq. 4).

**Two versions of the model** can be chosen at the top of the sidebar:

- **Fixed flounder density** (main text, Fig. 1): *N*<sub>F</sub> is a parameter.
- **Dynamic flounder density** (Supporting Information, Fig. S2): flounder density follows the same kind of equation as cod, with their own intrinsic growth, hypoxia, fishing (effort η<sub>F</sub>), competition among flounder, and competition from cod. The *N*<sub>F</sub> slider then only sets the initial flounder density of simulations.


### Model overview tab

**Rate of change of cod mean size.** d*μ*<sub>C</sub>/d*t* as a function of cod mean size, without fishing (grey) and at the current fishing effort η (red). Where a curve crosses zero, cod mean size is at an equilibrium: filled circles are stable equilibria, open circles unstable ones. Above zero, cod mean size increases; below zero, it decreases. The dotted line marks the flounder's mean size and the dashed line the optimal size *z*\*.

**Bifurcation diagram.** The equilibria of cod mean size as a function of fishing effort, from 0 up to the value set by *Largest fishing effort shown*. Solid lines are stable equilibria, dashed lines unstable ones, and grey dotted lines equilibria at which the cod population would be extinct (see the notes below). The red vertical line marks the current fishing effort η. Where they fall within the plotted range, two thresholds are marked:

- η<sub>crit</sub>, the critical fishing effort: above it, the large-cod state no longer exists, and cod tip to the small-bodied state.
- η<sub>rec</sub>, the recovery threshold: below it, the small-cod state no longer exists. If η<sub>rec</sub> < 0 (not shown), the small-cod state persists even without fishing, so that reducing fishing alone cannot reverse the shift (the irreversible trap). If 0 < η<sub>rec</sub> < η<sub>crit</sub>, cod return to large sizes once fishing is reduced below η<sub>rec</sub>.

**Terms of the model.** Each process at the current fishing effort, as a function of cod mean size, with the sign it enters the model with, and their sum (black):

- *Per capita growth rate of rare cod*: intrinsic growth minus mortality from fishing and hypoxia, minus competition with flounder. The total is the growth rate of cod when rare, i.e. without competition among cod themselves.
- *Selection on body size*: the selection on cod mean size from each process. The total equals d*μ*<sub>C</sub>/d*t* divided by *h*<sup>2</sup><sub>C</sub>, so it crosses zero exactly at the equilibria shown in the panel above. This panel shows which force dominates at each body size: for instance, whether flounder competition (pushing cod down when they are smaller than flounder) is strong enough to hold cod below the flounder's size against the pull of intrinsic growth towards *z*\*.

**How these are computed.** For a given cod mean size, the population densities settle to their ecological equilibrium: with fixed *N*<sub>F</sub>, the cod density follows from its growth rate; with dynamic *N*<sub>F</sub>, cod and flounder densities solve a two-species Lotka-Volterra system. Evaluated at these densities, d*μ*<sub>C</sub>/d*t* depends linearly on fishing effort, so the fishing effort at which any given body size is an equilibrium can be solved for directly. This gives the whole bifurcation diagram at once, without simulation. Stability is determined from the slope of this curve (fixed *N*<sub>F</sub>) or from the eigenvalues of the full three-variable system (dynamic *N*<sub>F</sub>).

In the dynamic model, the Model overview tab assumes that flounder are fished in proportion to cod, at the ratio η<sub>F</sub>/η set in the sidebar, so that what it shows at the current η corresponds to the fishing period of the simulation, and η = 0 means no fishing on either species.


### Simulation tab

The full model is integrated over time (from 0 to *Simulation length*), with fishing switched on during the *Fishing period*: cod are fished at effort η and, in the dynamic model, flounder at η<sub>F</sub>. Outside the fishing period there is no fishing.

If *Start from the equilibrium without fishing* is ticked (the default), the model is first run for 2000 time units without fishing, starting from the initial values in the sidebar, and the simulation starts from the state reached. Otherwise it starts from the initial values directly.

**Animation.** The body size distributions of cod (blue) and flounder (green) at the current time, each drawn as a normal distribution multiplied by population density, so that the area under each curve equals the population's density. The dashed line marks the optimal size *z*\*. The title turns red while fishing is on, and the line below it gives the current densities and cod mean size. Use **Play/Pause**, drag the time slider to jump to any time, and choose the playback speed on the right (at 1×, a full run takes about 20 seconds).

**Time series.** Mean body sizes (top) and population densities (bottom) over the whole simulation. The shaded region is the fishing period, and the vertical line follows the time shown in the animation.


### Parameters

| Parameter | Sidebar label | Default | Meaning |
|---|---|---|---|
| η | Fishing effort on cod | 3 | Maximum fishing mortality (for large cod) |
| η<sub>F</sub> | Fishing effort on flounder | 0 | Same for flounder (dynamic model only) |
| φ | Fishing inflection point | 1 | Body size at which fishing mortality is half its maximum |
| τ | Fishing transition width | 0.5 | Width of the transition from low to high fishing mortality |
| *N*<sub>F</sub> | Flounder density | 12 | Flounder density (fixed model) or its initial value (dynamic model) |
| ϱ | Maximum intrinsic growth | 5 | Intrinsic growth rate at the optimal body size |
| θ | Growth range width | 5 | Range of body sizes over which intrinsic growth is positive |
| *z*\* | Optimal body size | 1 | Body size at which cod intrinsic growth is highest |
| κ | Hypoxic effect intensity | 2 | Maximum mortality from hypoxia |
| ζ | Hypoxia inflection point | 1 | Body size at which hypoxic mortality is half its maximum |
| ν | Hypoxia transition width | 1.5 | Width of the transition from low to high hypoxic mortality |
| α<sub>0</sub> | Baseline competition | 0.707 | Overall strength of competition |
| α<sub>I</sub> | Non-size niche overlap | 0.5 | Scaling of competition between species, from overlap in traits other than body size |
| ω | Competition width | 0.33 | How quickly competition declines with difference in body size |
| σ<sub>C</sub>, σ<sub>F</sub> | Cod / flounder size SD | 0.5 | Standard deviation of body size within each species |
| *h*<sup>2</sup><sub>C</sub> | Cod heritability | 0.5 | Heritability of cod body size |

Simulation defaults: length 150, fishing during 0–50, initial cod mean size μ<sub>C</sub>(0) = 1.5 and density *N*<sub>C</sub>(0) = 5, starting from the equilibrium without fishing.

The defaults are the parameters of Fig. 1 of the manuscript. In the fixed model, α<sub>0</sub>, α<sub>I</sub> and *N*<sub>F</sub> only affect the dynamics through their product α<sub>0</sub>α<sub>I</sub>*N*<sub>F</sub>: for example, α<sub>I</sub> = 1 with *N*<sub>F</sub> = 6 gives exactly the same results. In the dynamic model, flounder settle at a density of about 12 regardless of α<sub>I</sub>, so α<sub>I</sub> = 0.5 gives about the same competitive pressure on cod as in the fixed model.


### Notes

- **Heritability** *h*<sup>2</sup><sub>C</sub> sets how fast cod body size evolves, but not where the equilibria and thresholds are: it does not affect the Model overview tab apart from scaling d*μ*<sub>C</sub>/d*t*.
- **Extinct cod.** In the fixed model, cod mean size has the same equilibria whether or not the cod population can persist there; the bifurcation diagram shows such equilibria as "cod extinct". In the dynamic model, the "cod extinct" branches are equilibria with flounder alone, which cod cannot invade.
