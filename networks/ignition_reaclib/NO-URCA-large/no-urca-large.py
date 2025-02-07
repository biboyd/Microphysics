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
          'c12', 'c13', 'c14',
          'n13', 'n14', #'n15',
          'o16', 'o17', 
          'o18', #'o19', #'o20',
          'f18', #'f19', #'f20', 
          'f21', #'f22',
          'ne20', 'ne21', 'ne22', 'ne23',# 'ne24',
          #'na22', 
          'na23', #'na24', 
          'na25',
          #'mg23', 
          'mg24', 'mg25', #'mg26', 'mg27', #'mg28',
          #'al25', 
          #'al26', 'al27', #'al28',
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
#all_lib.remove_rate(all_lib.get_rate_by_name("o17(a,n)ne20"))
#all_lib.remove_rate(all_lib.get_rate_by_name("f19(e,)o19"))
#all_lib.remove_rate(all_lib.get_rate_by_name("o19(,e)f19"))
all_lib.remove_rate(all_lib.get_rate_by_name("ne22(a,n)mg25"))

# testing
all_lib.remove_rate(all_lib.get_rate_by_name("o16(c12,a)mg24"))
#all_lib.remove_rate(all_lib.get_rate_by_name("o16(c12,p)al27")) already removed
#all_lib.remove_rate(all_lib.get_rate_by_name("al27(e,)mg27"))
#all_lib.remove_rate(all_lib.get_rate_by_name("mg27(p,n)al27"))
all_lib.remove_rate(all_lib.get_rate_by_name("na25(p,a)ne22"))
all_lib.remove_rate(all_lib.get_rate_by_name("na25(p,n)mg25"))

#quick test
#all_lib.remove_rate(all_lib.get_rate_by_name("n13(e,)c13"))
#all_lib.remove_rate(all_lib.get_rate_by_name("na23(e,)ne23"))
#all_lib.remove_rate(all_lib.get_rate_by_name("ne23(,e)na23"))

for r in unimportant_rates:
    all_lib.remove_rate(r)

#all_lib.add_rate(tl.get_rate_by_name("mg27(,)al27"))
#all_lib.add_rate(rl.get_rate_by_name("ne23(p,)na24"))
#all_lib.add_rate(rl.get_rate_by_name("ne23(a,)mg27"))
#all_lib.add_rate(rl.get_rate_by_name("ne23(n,)ne24"))
#all_lib.add_rate(rl.get_rate_by_name("c12(c12,n)mg23"))
#all_lib.add_rate(rl.get_rate_by_name("mg23(n,)mg24"))

# remove select Urca rates                                                                                                                                                                      
all_lib.remove_rate(all_lib.get_rate_by_name("na25(,)mg25"))                                                                                                                                    
#all_lib.remove_rate(all_lib.get_rate_by_name("mg25(,)na25"))                                                                                                                                   
all_lib.remove_rate(all_lib.get_rate_by_name("ne23(,)na23"))                                                                                                                                    
#all_lib.remove_rate(all_lib.get_rate_by_name("na23(,)ne23"))                                                                                                                                   
all_lib.remove_rate(all_lib.get_rate_by_name("f21(,)ne21"))                                                                                                                                     
#all_lib.remove_rate(all_lib.get_rate_by_name("ne21(,)f21"))                                                                                                                                    

urca_large_net = pyna.AmrexAstroCxxNetwork(libraries=[all_lib])
print(all_lib)
print(f"Total Nuclei: {len(urca_large_net.unique_nuclei)}")
print(f"Total Rates Included: {len(urca_large_net.rates)}")
urca_large_net.write_network()

for i, (rho, T, comp) in enumerate(states):
    urca_large_net.plot(outfile=f"urca_large_{i}.png", rho=rho, T=T, comp=comp,screen=chugunov_2009, 
                    rotated=True, curved_edges=True, ydot_cutoff_value=1e-20)#, always_show_p=True, always_show_alpha=True)

urca_large_net.plot(outfile="urca_large.png", #rho=rho, T=T, comp=comp,screen=chugunov_2009, 
                    rotated=True, curved_edges=True)#, ydot_cutoff_value=1e-20, always_show_p=True, always_show_alpha=True)

#print(urca_big_net)
