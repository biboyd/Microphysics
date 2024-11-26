import pynucastro as pyna
from pynucastro.screening import chugunov_2009
import pandas as pd


def get_states(network):
    states = []

    #set composotion
    comp = pyna.Composition(network.unique_nuclei)
    with open("react_large_comp_final.txt", "r") as f:
        for line in f:
            spline = line.split(' ')
            compX = float(spline[2])
            compName = spline[-1].removeprefix("#").removesuffix('\n')

            comp.set_nuc(compName, compX)

    comp.normalize()

    rho_arr = [4.5e9, 4.023242e+09, 3.350290e+09, 1.629665e+09, 1.304853e+09]
    temp_arr = [5.5e8, 5.188529e+08, 4.716727e+08, 3.240484e+08, 2.886596e+08]

    for r, T in zip(rho_arr, temp_arr):
        states.append((r, T, comp))

    return states


nuclei = ["p", 'he4', 'n',
          'c12', 'c13', #'c14',
          'n13', 'n14', #'n15',
          'o16', 'o17', 'o18', 'o19', #'o20',
          'f18', 'f19', #'f20', 
          'f21', #'f22',
          'ne20', 'ne21', 'ne22', 'ne23',# 'ne24',
          #'na22', 
          'na23', #'na24', 
          'na25',
          #'mg23', 
          'mg24', 'mg25', #'mg26', 
          'mg27', #'mg28',
          #'al25', 
          #'al26', 
           'al27', #'al28',
          #'si27', 'si28', 'si29', 'si30', 'si31', 'si32',
          #'p30', 'p31', 'p32'
         ]

rl = pyna.ReacLibLibrary()
rl_rates = rl.linking_nuclei(nuclei, with_reverse=False)

tl = pyna.TabularLibrary()
tl_rates = tl.linking_nuclei(nuclei)

complete_lib = rl_rates + tl_rates # don't edit this one

all_lib = rl_rates + tl_rates # will reduce this one
dupes = all_lib.find_duplicate_links()

rates_to_remove = []
for d in dupes:
    for r in d:
        if isinstance(r, pyna.rates.ReacLibRate):
            rates_to_remove.append(r)

for r in rates_to_remove:
    all_lib.remove_rate(r)

cutoff = 1e-11
urca_big_net = pyna.AmrexAstroCxxNetwork(libraries=[all_lib])
states = get_states(urca_big_net)
unimportant_rates = urca_big_net.find_unimportant_rates(states, cutoff, chugunov_2009)

all_lib.remove_rate(all_lib.get_rate_by_name("o18(a,n)ne21"))
all_lib.remove_rate(all_lib.get_rate_by_name("o17(a,n)ne20"))
all_lib.remove_rate(all_lib.get_rate_by_name("f19(e,)o19"))
all_lib.remove_rate(all_lib.get_rate_by_name("o19(,e)f19"))
all_lib.remove_rate(all_lib.get_rate_by_name("ne22(a,n)mg25"))

for r in unimportant_rates:
    all_lib.remove_rate(r)
print(all_lib)

urca_large_net = pyna.AmrexAstroCxxNetwork(libraries=[all_lib])
print(f"Total Nuclei: {len(urca_large_net.unique_nuclei)}")
print(f"Total Rates Included: {len(urca_large_net.rates)}")
urca_large_net.write_network()
urca_large_net.plot(outfile="urca_large.png", rotated=True)

#print(urca_big_net)
