# C-burning with A=23 URCA rate module generator
import pynucastro as pyna

tl = pyna.TabularLibrary()
rl = pyna.ReacLibLibrary()

# get Tabular rates
tl_names = ["ne21(e,_)f21",
            "f21(_,e)ne21",
            "na23(e,_)ne23",
            "ne23(_,e)na23",
            "mg25(e,_)na25",
            "na25(_,e)mg25"]
tl_rates = tl.get_rate_by_name(tl_names)

# get ReacLib rates
rl_names = [
            "c12(c12,a)ne20",
            "c12(c12,n)mg23",
            "c12(c12,p)na23",
            "c12(a,g)o16",
            "n(_,e)p",]
rl_rates = rl.get_rate_by_name(rl_names)

rates = tl_rates + rl_rates
urca_net =  pyna.AmrexAstroCxxNetwork(rates=rates)
urca_net.write_network()
