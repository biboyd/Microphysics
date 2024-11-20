import pynucastro as pyna

rate_ids = [
            "N13 --> p + C12 <ls09_reaclib__reverse>",
            "C12 + n --> C13 <ks03_reaclib__>",
            "C12 + p --> N13 <ls09_reaclib__>",
            "C12 + He4 --> O16 <nac2_reaclib__>",
            "O16 + n --> O17 <ks03_reaclib__>",
            "O16 + He4 --> Ne20 <co10_reaclib__>",
            "Ne22 + n --> Ne23 <ks03_reaclib__>",
            "Ne22 + p --> Na23 <ke17_reaclib__>",
            "Ne23 + He4 --> Mg27 <rath_reaclib__>",
            "Na23 + He4 --> Al27 <ths8_reaclib__>",
            "C12 + C12 --> p + Na23 <cf88_reaclib__>",
            "C12 + C12 --> He4 + Ne20 <cf88_reaclib__>",
            "C13 + He4 --> n + O16 <gl12_reaclib__>",
            "N13 + n --> p + C13 <nacr_reaclib__reverse>",
            "N13 + He4 --> p + O16 <cf88_reaclib__>",
            "O16 + C12 --> p + Al27 <cf88_reaclib__>",
            "O17 + He4 --> n + Ne20 <nacr_reaclib__>",
            "Ne20 + n --> He4 + O17 <nacr_reaclib__reverse>",
            "Ne23 + p --> n + Na23 <rath_reaclib__>",
            "Na23 + p --> He4 + Ne20 <il10_reaclib__>",
            "Mg27 + p --> n + Al27 <rath_reaclib__>",
            "N13 --> C13 <tabular_tabular>",
            "p --> n <tabular_tabular>",
            "Al27 --> Mg27 <tabular_tabular>",
            "Na23 --> Ne23 <tabular_tabular>",
            "Ne23 --> Na23 <tabular_tabular>"]

tl = pyna.TabularLibrary()
full_lib = pyna.ReacLibLibrary() + tl
final_lib = pyna.Library()

for rid in rate_ids:
    final_lib.add_rate(full_lib.get_rate(rid))

#explicitly ad A=21,25,27 pairs
final_lib.add_rate(tl.get_rate_by_name("ne21(,)f21"))
final_lib.add_rate(tl.get_rate_by_name("f21(,)ne21"))
final_lib.add_rate(tl.get_rate_by_name("na25(,)mg25"))
final_lib.add_rate(tl.get_rate_by_name("mg25(,)na25"))
final_lib.add_rate(tl.get_rate_by_name("mg27(,)al27"))

urca_net = pyna.AmrexAstroCxxNetwork(libraries=[final_lib])
urca_net.write_network()
