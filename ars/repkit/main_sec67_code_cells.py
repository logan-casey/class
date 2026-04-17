
# ==== CODE CELL 0 ====

import numpy as np
from scipy import optimize
import json

# ==== CODE CELL 1 ====

import sequence_jacobian as sj  
from base import calibration, fiscal, winding_number, capital_sticky_prices
from models import models_analytical, models_heterogeneous
from plotting import sec67_plots as plots

opts = {'texfig': True, 'savefig': True}

# ==== CODE CELL 2 ====

@sj.simple
def fiscalrule(B, r, G, w, N):
    T = (1 + r(-1)) * B(-1) + G - B # taxes
    Z = w*N - T # after-tax income
    return T, Z

@sj.simple
def wage_inflation(pi, w):
    piw = (1 + pi) * w / w(-1) - 1
    return piw

@sj.simple
def union(piw, N, UCE, Z, kappaw, vphi, frisch, muw, beta, theta):
    wnkpc = (kappaw * (muw * vphi * N ** (1 / frisch)/(UCE * (1-theta)*Z/N) - 1) +
             beta * piw(+1)*(1 + piw(+1)) - piw*(1 + piw))
    return wnkpc

@sj.simple
def taylor(r, pi, phi_pi):
    i = r.ss + phi_pi * pi
    return i

@sj.simple
def share_bond(B, J_end_period):
    bondshare = B / ( B + J_end_period)
    return bondshare

@sj.simple
def finance(r, J, J_end_period, bondshare, i, pi):
    rpost = bondshare(-1) * (1+r(-1)) + (1 - bondshare(-1)) * J/J_end_period(-1) - 1
    fisher = 1 + i - (1 + r) * (1 + pi(+1))
    return rpost, fisher

@sj.simple
def mkt_clearing(A, B, J_end_period, C, I, G, p_adjust, k_adjust, Y):
    asset_mkt = A - B - J_end_period
    goods_mkt = C + I + G + p_adjust + k_adjust - Y
    return asset_mkt, goods_mkt

# ==== CODE CELL 3 ====

other_blocks = [fiscalrule, wage_inflation, union, taylor, 
                capital_sticky_prices.production_pricesetting, share_bond, finance, mkt_clearing]

# ==== CODE CELL 4 ====

with open('intermediates/solved_params.json', 'r') as f:
    params = json.load(f)

# ==== CODE CELL 5 ====

hh_analytical, ss_analytical = models_analytical.get_all_quant(params)
hh_het, ss_het = models_heterogeneous.get_all_quant(params)
hh_all = {**hh_analytical, **hh_het}
ss_all = {**ss_analytical, **ss_het}
hh_names = {m: hh_all[m].name for m in hh_all}

# ==== CODE CELL 6 ====

calib = calibration.get_calibration_quant()
muw, N, frisch, theta, Z = [calib[k] for k in ('muw', 'N', 'frisch', 'theta', 'Z')]
vphi_base = (1-theta)*Z/N**(1+1/frisch)/muw
for m in ss_all:
    ss_all[m]['vphi'] = vphi_base*ss_all[m]['UCE']

# ==== CODE CELL 7 ====

models_all = {m: sj.combine([*other_blocks, hh_all[m]]) for m in hh_all}
for m in models_all:
    if m in ('RA', 'TA'):
        ss_all[m] = models_all[m].steady_state({**ss_all[m], **calib}, dissolve=[f'hh_{m.lower()}'])
    else:
        ss_all[m] = models_all[m].steady_state({**ss_all[m], **calib})

# ==== CODE CELL 8 ====

for m in models_all:
    for k in ('wnkpc', 'pi', 'piw', 'fisher', 'asset_mkt'):
        assert np.isclose(ss_all[m][k], 0)
    assert np.isclose(ss_all[m]['goods_mkt'], 0, atol=1E-5) # less exact for HA-two

# ==== CODE CELL 9 ====

plots.table3({**calib, 'rhoG': calibration.rhoG, 'rhoB': calibration.rhoB})

# ==== CODE CELL 10 ====

T = 500
Js_cache = {m: hh_all[m].partial_jacobians(ss_all[m], inputs=['Z', 'rpost'], outputs=['C', 'A', 'UCE'], T=T) for m in models_all}
Js = {m: hh_all[m].jacobian(ss_all[m], inputs=['Z', 'rpost'], outputs=['C', 'A'], T=T, Js=Js_cache[m]) for m in models_all}

# ==== CODE CELL 11 ====

Ms_analytical, As_analytical = models_analytical.MA_all(params, calib['r'], T)
for m in ('BU', 'TABU'):
    assert np.allclose(Js[m]['C', 'Z'][:-50, :-50], Ms_analytical[m][:-50, :-50])
    assert np.allclose(Js[m]['A', 'Z'][:-50, :-50], As_analytical[m][:-50, :-50])

# ==== CODE CELL 12 ====

mcap = {m: Js[m]['C', 'rpost'][:, 0] / calib['A'] for m in Js}
mcap['RA'] = mcap['TA'] = Ms_analytical['RA'][:, 0] # numerical Jacobian misses unit root

mlab = {m: np.abs(Js[m]['C', 'Z'][:, 0]) for m in Js}
mlab['RA'] = Ms_analytical['RA'][:, 0]
mlab['TA'] = Ms_analytical['TA'][:, 0]

# ==== CODE CELL 13 ====

plots.figure7(mcap, mlab, **opts)

# ==== CODE CELL 14 ====

# specify unknowns to solve for and GE targets on DAG to hit
unknowns = ['Y', 'w', 'r']
targets = ['asset_mkt', 'wnkpc', 'fisher']
inputs = ['B', 'G']
outputs = ['C', 'I', 'N', 'pi', 'Y', 'r', 'Z', 'rpost']

Gs = {m: models_all[m].solve_jacobian(ss_all[m], unknowns, targets, inputs, outputs, T=T, Js=Js_cache[m]) for m in models_all}

# ==== CODE CELL 15 ====

rhoG = 0.76
rhoB = 0.93

dG = rhoG**np.arange(T)
dB = fiscal.Bplan(dG, rhoB)
shock = {'G': dG, 'B': dB}

impulses = {m: Gs[m] @ shock for m in ('HA-two', 'RA', 'TA', 'TABU')}

# ==== CODE CELL 16 ====

plots.figure8(impulses, **opts)

# ==== CODE CELL 17 ====

decomp = {}
for m in ('HA-two', 'TABU'):
    # load impulse responses that are outputs or inputs to the household
    dC, dZ = impulses[m]['C'], impulses[m]['Z']
    dcap = impulses[m]['rpost']*(np.arange(T) == 0) # change from cap gain: only rpost_0
    dr = impulses[m]['rpost']*(np.arange(T) != 0)   # change from r: everything except rpost_0

    # now store the contributions from each component, check that they sum to dC
    decomp[m] = {}
    for dx, name, i in zip((dZ, dcap, dr), ('Income', 'Capital gain', 'Rate'), ('Z', 'rpost', 'rpost')):
        decomp[m][name] = Js[m]['C', i] @ dx

    assert np.allclose(sum(decomp[m][x] for x in decomp[m]), dC)
    decomp[m] = {'Overall': dC, **decomp[m]}
    

# ==== CODE CELL 18 ====

plots.figure10(decomp, **opts)

# ==== CODE CELL 19 ====

rho_Bs_short = [0, 0.5, 0.76, 0.93]
rho_Bs_labels = [rf'$\rho_B = {rho_B}$' for rho_B in rho_Bs_short]

# ==== CODE CELL 20 ====

impulses = {}
for rhoB, label in zip(rho_Bs_short, rho_Bs_labels):
    dB = fiscal.Bplan(dG, rhoB)
    shock = {'G': dG, 'B': dB}
    impulses[label] = Gs['HA-two'] @ shock
plots.figureG(impulses, 'green', '1', **opts)

# ==== CODE CELL 21 ====

impulses = {}
for rhoB, label in zip(rho_Bs_short, rho_Bs_labels):
    dB = fiscal.Bplan(dG, rhoB)
    shock = {'G': dG, 'B': dB}
    impulses[label] = Gs['TABU'] @ shock
plots.figureG(impulses, 'purple', '2', **opts)

# ==== CODE CELL 22 ====

NrhoB = 20
rho_Bs = np.linspace(0, 0.93, NrhoB)
mult_impact, mult_cumul = {}, {}
for m in ('TABU', 'HA-two', 'RA', 'TA'):
    mult_impact[m], mult_cumul[m] = np.empty(NrhoB), np.empty(NrhoB)
    G_Y = Gs[m][['Y'], :] # only interested in multiplier on output here
    for i, rhoB in enumerate(rho_Bs):
        dB = fiscal.Bplan(dG, rhoB)
        shock = {'G': dG, 'B': dB}
        dY = (G_Y @ shock)['Y']
        mult_impact[m][i], _, mult_cumul[m][i] = fiscal.compute_multipliers(dY, dG, calib['r'])

# ==== CODE CELL 23 ====

plots.figure9(rho_Bs, mult_impact, mult_cumul, **opts)

# ==== CODE CELL 24 ====

plots.table4(mult_impact, mult_cumul, ikc=False).round(1)

# ==== CODE CELL 25 ====

dB = fiscal.Bplan(dG, 0.93)
shock = {'G': dG, 'B': dB}
phi_pis = [1.1, 1.5, 1.7, 2]
phi_pi_labels = [rf'$\phi_\pi = {phi_pi}$' for phi_pi in phi_pis]
outputs=['C', 'I', 'N', 'pi', 'Y', 'r']

impulses = {}
for phi_pi, label in zip(phi_pis, phi_pi_labels):
    ss_alt = ss_all['HA-two'].copy()
    ss_alt['phi_pi'] = phi_pi
    impulses[label] = models_all['HA-two'].solve_jacobian(ss_alt, unknowns, targets, inputs, outputs, T=T, Js=Js_cache['HA-two']) @ shock

plots.figureG(impulses, 'green', '3', **opts)

# ==== CODE CELL 26 ====

eps_Is = [1, 4, 10, 20]
eps_I_labels = [rf'$\epsilon_I = {eps_I}$' for eps_I in eps_Is]

impulses = {}
for eps_I, label in zip(eps_Is, eps_I_labels):
    ss_alt = ss_all['HA-two'].copy()
    ss_alt['epsI'] = eps_I
    impulses[label] = models_all['HA-two'].solve_jacobian(ss_alt, unknowns, targets, inputs, outputs, T=T, Js=Js_cache['HA-two']) @ shock

plots.figureG(impulses, 'green', '4', **opts)

# ==== CODE CELL 27 ====

kappa_ps = [calib['kappap']*6 / (1 + Gammap) for Gammap in (9, 5, 1, 0)]
kappa_p_labels = [rf'$\kappa_p = {kappa_p:.3f}$' for kappa_p in kappa_ps]

impulses = {}
for kappa_p, label in zip(kappa_ps, kappa_p_labels):
    ss_alt = ss_all['HA-two'].copy()
    ss_alt['kappap'] = kappa_p
    impulses[label] = models_all['HA-two'].solve_jacobian(ss_alt, unknowns, targets, inputs, outputs, T=T, Js=Js_cache['HA-two']) @ shock

plots.figureG(impulses, 'green', '5', **opts)

# ==== CODE CELL 28 ====

kappa_ws = [calib['kappaw']*6 / (1 + Gammap) for Gammap in (9, 5, 1, 0)]
kappa_w_labels = [rf'$\kappa_w = {kappa_w:.3f}$' for kappa_w in kappa_ws]

impulses = {}
for kappa_w, label in zip(kappa_ws, kappa_w_labels):
    ss_alt = ss_all['HA-two'].copy()
    ss_alt['kappaw'] = kappa_w
    impulses[label] = models_all['HA-two'].solve_jacobian(ss_alt, unknowns, targets, inputs, outputs, T=T, Js=Js_cache['HA-two']) @ shock

plots.figureG(impulses, 'green', '6', **opts)

# ==== CODE CELL 29 ====

J_taylor = taylor.jacobian(ss_all['HA-two'], inputs=['pi'], outputs=['i'], T=T)['i', 'pi'].matrix(T)
J_taylor[:5, :5]

# ==== CODE CELL 30 ====

J_constant_r = np.diag(np.full(T-1, 1+calib['r']), 1)
J_constant_r[3:] = J_taylor[3:]

# ==== CODE CELL 31 ====

J_zlb = J_taylor.copy()
J_zlb[:3] = 0

# ==== CODE CELL 32 ====

names = ['Taylor rule', 'Constant-r (3 years)', 'ZLB (3 years)']
models_alt = {}
other_blocks_alt = other_blocks.copy()
for name, Ji in zip(names, (J_taylor, J_constant_r, J_zlb)):
    other_blocks_alt[3] = sj.JacobianDict({'i': {'pi': Ji}}, name='monetary')
    models_alt[name] = sj.combine([*other_blocks_alt, hh_all['HA-two']])

# ==== CODE CELL 33 ====

Gs_alt = {m: models_alt[m].solve_jacobian(ss_all['HA-two'], unknowns, targets, inputs, outputs=['Y'], T=T, Js=Js_cache['HA-two']) for m in models_alt}

# ==== CODE CELL 34 ====

shock = {'G': dG, 'B': fiscal.Bplan(dG, 0.76)}
dYs_rho_76 = {m: (Gs_alt[m] @ shock)['Y'] for m in models_alt}

shock = {'G': dG, 'B': fiscal.Bplan(dG, 0.93)}
dYs_rho_93 = {m: (Gs_alt[m] @ shock)['Y'] for m in models_alt}

# ==== CODE CELL 35 ====

plots.figure11(dYs_rho_76, dYs_rho_93, **opts)

# ==== CODE CELL 36 ====

ss = ss_all['TABU'].copy()
ss['phi_pi'] = 1.2
H_U = models_all['TABU'].jacobian(ss, unknowns, targets, T=T, Js=Js_cache['TABU'])
winding_number.from_H_U(H_U, unknowns, targets)

# ==== CODE CELL 37 ====

ss['phi_pi'] = 1.0
H_U = models_all['TABU'].jacobian(ss, unknowns, targets, T=T, Js=Js_cache['TABU'])
winding_number.from_H_U(H_U, unknowns, targets)

# ==== CODE CELL 38 ====

for m in ('BU', 'TABU', 'HA-hi-liq', 'HA-two'):
    H_U = models_all[m].jacobian(ss_all[m], unknowns, targets + ['pi'], T=T, Js=Js_cache[m])
    A = winding_number.pack_jacdict_center(H_U, unknowns, targets+['pi'])

    def winding(phi_pi):
        A0 = A[:, :-1, :].copy()                                     # copy everything except pi output
        A0[:, -1, :] += (phi_pi - ss_all[m]['phi_pi']) * A[:, -1, :] # manually correct 'fisher' for phi_pi change
        return winding_number.winding_number(A0)
    
    phi_pi = optimize.bisect(lambda x: winding(x)+0.5, 0.9, 1.2, xtol=1E-5) # where does winding go from -1 to 0?
    print(f'Threshold for {m}: {phi_pi:.2f}')


# ==== CODE CELL 39 ====

ss_two_alt = ss_all['HA-two'].copy()
ss_two_alt['liquidshare'] = ss_two_alt['Aliq'] / ss_two_alt['A']
hh_two_alt = models_heterogeneous.ha_two_alt
model_two_alt = sj.combine([*other_blocks, hh_two_alt])

# ==== CODE CELL 40 ====

mcap['HA-two-alt'] = hh_two_alt.jacobian(ss_two_alt, inputs=['rpost'],
                            outputs=['C'], T=25)['C','rpost'][:, 0] / ss_two_alt['A']

# ==== CODE CELL 41 ====

plots.figureF1a(mcap, **opts)

# ==== CODE CELL 42 ====

mcap['HA-two-alt'][:10]

# ==== CODE CELL 43 ====

G_C_alt = model_two_alt.solve_jacobian(ss_two_alt, unknowns, targets, inputs, ['C'], T=T, Js=Js_cache['HA-two'])
G_C = Gs['HA-two']

# ==== CODE CELL 44 ====

shock_bb = {'G': dG, 'B': 0*dG}
shock_df = {'G': dG, 'B': fiscal.Bplan(dG, 0.93)}

# ==== CODE CELL 45 ====

dC_bb = (G_C @ shock_bb)['C']
dC_bb_alt = (G_C_alt @ shock_bb)['C']
dC_df = (G_C @ shock_df)['C']
dC_df_alt = (G_C_alt @ shock_df)['C']

# ==== CODE CELL 46 ====

plots.figureF1b(dC_bb, dC_bb_alt, dC_df, dC_df_alt, **opts)