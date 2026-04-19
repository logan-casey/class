import sequence_jacobian as sj
hh = sj.hetblocks.hh_sim.hh

print(hh)
print(f'Inputs: {hh.inputs}')
print(f'Macro outputs: {hh.outputs}')
print(f'Micro outputs: {hh.internals}')