import pynucastro as pyna

# C-burning with A=23 URCA rate module generator
rl = pyna.ReacLibLibrary()
rl_lib = rl.linking_nuclei(["n", "p", "he4", "c12", "n13", "c13", "o16",
                              "ne20", "na23", "mg24", "mg25", "na25"])

# let's remove some unimportant reverse rates (we are not hot
# enough to drive these
rates_to_remove = ["mg24(a,c12)o16",
                   "mg24(g,a)ne20",
                   "mg24(,p)na23",
                   "c12(g,aa)he4",
                   "n13(g,p)c12",
                   "o16(p,a)n13",
                   "o16(g,a)c12",
                   "ne20(g,a)o16",
                   "na23(p,c12)c12",
                   "ne20(a,c12)c12",
                   "c13(p,n)n13",
                   "c13(,n)c12",
                   "o16(n,a)c13",
                   "ne20(a,p)na23",
                   "mg25(n,p)na25",
                   "mg25(,n)mg24",
                   "n(,)p",
                   "aa(a,)c12",
                   "na23(p,a)ne20",
                   "na25(p,n)mg25",
                   "c12(o16,a)mg24",
                   "n13(a,p)o16",
                   "o16(a,)ne20",
                   "ne20(a,)mg24",
                   #
                   # actually maybe unnecessary "c12(a,)o16",
                   
                   #"n13(n,p)c13",
                   #"c13(a,n)o16",
                   #"c12(p,)n13",
                   #"c12(n,)c13",
                   #"na23(p,)mg24",
                   #"mg24(n,)mg25",
                   #"n13(,)c13",
                  ]


for r in rates_to_remove:
    print(f"removing {r:15s} : {rl.get_rate_by_name(r).Q}")
    rl_lib.remove_rate(r)

tl = pyna.TabularLibrary()
tl_rates = tl.get_rate_by_name(["na23(,)ne23",
                                "ne23(,)na23",
                                "na25(,)mg25",
                                "mg25(,)na25"])
print("tl_rates: ")
print(tl_rates)

tl_lib = pyna.Library(rates=tl_rates)

#print(tl_lib)

all_lib = rl_lib + tl_lib
# we will have duplicate rates -- we want to remove any ReacLib rates
# that we have tabular rates for

dupes = all_lib.find_duplicate_links()

rates_to_remove = []
for d in dupes:
    for r in d:
        if isinstance(r, pyna.rates.ReacLibRate):
            rates_to_remove.append(r)

for r in rates_to_remove:
    all_lib.remove_rate(r)

print(all_lib)

urca23_25_reduced_net = pyna.AmrexAstroCxxNetwork(libraries=[all_lib])
urca23_25_reduced_net .write_network()
