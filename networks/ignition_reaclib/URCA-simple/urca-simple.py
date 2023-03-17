# C-burning with A=23 URCA rate module generator
import pynucastro as pyna

tl = pyna.TabularLibrary()
rl = pyna.ReacLibLibrary()

files = ["c12-c12a-ne20-cf88",
         "c12-c12n-mg23-cf88",
         "c12-c12p-na23-cf88",
         "c12-ag-o16-nac2",
         "na23--ne23-toki",
         "ne23--na23-toki",
         "n--p-wc12"]
# get Tabular rates
tl_names = ["na23(e,)ne23", "ne23(,e)na23"]
tl_rates = tl.get_rate_by_name(tl_names)

# get ReacLib rates
rl_names = [
            "c12(c12,a)ne20",
            "c12(c12,n)mg23",
            "c12(c12,p)na23",
            "c12(a,g)o16",
            "n(,e)p",]
rl_rates = rl.get_rate_by_name(rl_names)

rates = tl_rates + rl_rates
urca_net =  pyna.AmrexAstroCxxNetwork(rates=rates)
urca_net.write_network()
