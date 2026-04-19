## Equilibrium
Given sequences:
    ROW discount shocks {Bt}
    monetary shocks {et}
    initial wealth distribution dDi0(s,BH,BF,e)
Eq is a path of:
    Policies {cHit(ap,e), cFit(ap,e), cit(ap, e), ait+1(ap, e)}
    distributions Dit(ap,e)
    prices {Et, Qt, Pt, PHt, PFt, Wt, rt, it}
    aggregate quantities {Ct, CHt, CFt,Yt, At, pt, dt, nfat}
such that:
    HH opotimize
    distributions evolve consistently w policies and initial condition
    firms optimize
    domestic goods mkt clears

## Steady state
zero NFA
zero inflation
normalize prices, output, consumption to 1

## Parameters
beta = discount factor
sigma = inverse EIS
phi_labor = inverse Frisch (labor supply)
eta = h/f goods elas
alpha = openness/home bias
psi = scaling param in labor supply
kappaw = phillips slope, depends on beta, :
    thetaw = Calvo prob
muw = wage markdown

**not calibrated directly
chi = eta*(1-alpha) + gamma
    with gamma = eta, chi = eta*(2-alpha)

phi_pi = Taylor rule coefficient
rhom = MP shock persistence

clowerbar = subsistence consumption
zeta = weird labor thing

mu = markup
thetah = goods Calvo prob
kappah = home goods Phillips slope
theta, kappa for foreign goods, home exports

theta_sub = consumption reallocation parameter

eps = domestic MP shock

## Equations

### Base equations (1)-(23)

(1) HH value function -- note three assets but one borrowing constraint -- unused
    not needed (in HH)
(2) utility function
    not needed (in HH/implicit)
(3) CES consumption aggregator cF, cH -> c -- replaced by (47)
    not needed if we only solve for H/F aggregate
(4) price index PH_t, PF_t -> P_t
    simple: price_index
(5) stock, bond, foreign bond return equalization (nominal UIP)
    not needed (real instead)
(6) real exchange rate E_t, Pstar_t, P_t -> Q_t
    simple: exchange_rate_identity (target residual)
(unnumbered) Fisher eqs -- altered to match constant-r MP rule
    simple: fisher_equation (written as target)
(7) real UIP + equalization with stock return
    simple: uip_real (target residual)
(8) real asset positions apre_t, apost_t
    not needed (in HH/implicit)
(9) HH value function -- one asset state variable -- usable version
    already in SSJ package sj.hetblocks.hh_sim.hh
(10)-(11) c policy function -> cH, cF policy function mapping
    not needed if we only solve for H/F aggregate
(unnumbered) C_t, A_t aggregates; distribution D_t; law of motion
    not needed if we work with these directly
(12) ROW H export demand CstarH_t, Cstar_t, Pstar_t, PstarH_t
    not needed if we only solve for H/F aggregate
(13) production function Y_t = N_t
    simple: production
(14) pricing constant markup mu, W_t -> PH_t -- relaxed?
    not needed; only Phillips
(15) dividends (profits) d_t
    can include but only needed for IRFs
(16) foreign MP rule (constant Pstar_t, PstarF_t, Cstar) istar_t, rstar_t --(varies with beta shock) -- is PstarH not constant then?
    simple: foreign rate block
(unnumbered) PF = E
    simple: import_price_pass_through
(17) LOOP for home good PstarH_t, PH_t, E_t -- (could relax with DCP)
    not needed if we only solve for H/F aggregate
(18) (wage) Phillips curve piw_t, piw_t+1 -- simpler than IKC version
    simple: union_wage_nkpc (target residual)
(19)-(20) domestic MP rule pi_t+1, r_ss -> i_t -- use 20. shock here!
    simple: monetary_rule_taylor (write as target)
(21) domestic goods mkt clearing CH_t CHstar_t Y_t
    simple: goods_market (write as target)
(22) net foreing asset position def A_t, p_t -> nfa_t
    only used in IRFs
(23) currrent account (nfa law of motion) -- note (PH/P)Y-C is net exports
    only used in IRFs
extra:
    simple: inflation_defs


### Quantitative equations (47)-(53)
(47) Stone-Geary non-homothetic CES aggregator -- replaces (3)
    not needed if we only solve for H/F aggregate
(48) crazy labor shock formulation
    used as hetinputs to HH block
(49) home goods Phillips curve piH_t
    simple: piH_nkpc (written as target)
(50) foreign goods Phillips curve piF_t
    not used in paper calibration (thetaF=0, no inflation)
(51) home goods foreign Phillips curve
    simple: piHstar_nkpc (written as target)
(52)-(53)-(54)-(55) delayed substitution (sticky consumption)
    simple: xH_target_law, xHstar_target_law, delayed_substitution
    see appendix D3

## 4) notes

portfolio choice is indeterminate -- this is one-asset HANK (so we are missing adjustment friction ...)
ROW is RA with exognenous patience shock (can think of as interest rate shock)
PCP, no foreign price rigidity
no government spending, no taxation, no debt
only solve for aggregate (not H/F)
quantitative modifications:
    non-homothetic consumption (poor consume more imports)
    heterogeneous labor income incidence
    price stickiness -- how do we deal with both
    delayed substitution (sticky consumption shares)


is delayed substitution alright? is HH block hetinputs alright
how did we get grids and SS anchors?
why do we have so many targets?
why do we not use so many equations involving cF, cH, pH, etc (eg the consumption aggregator)?
in general, how does the model as currently written differ from the DAG plan you made?
delayed substitution**
does eq 50 matter?


1. In implementation, equation (9) is the key HH dynamic block; (1) and most of portfolio detail become accounting identities unless you explicitly model endogenous portfolio composition.
2. Treat forward-looking equations as residual blocks and solve with model-level unknowns/targets (as in `main_sec67`).
3. Use equation switches:
   - baseline demand: (3), baseline pricing: (14), policy: (19)
   - quantitative extensions: replace with (47), (49)-(53), and optionally (20)
