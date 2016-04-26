"""
Script to finally get my head around putting things in terms of dBm!
"""

def main():

    # "radio frequency work dBm is typically referenced relative to a 50 ohm impedance"
    # from https://en.wikipedia.org/wiki/DBm
    # this also agrees with what Abby told me about the scope (50 ohm impedance)

    # x [dBm] = 10*math.log10( P[mW] / 1mW )

    # So if I have a power of 1mW, x is 0.
    # i.e. 0 dBm imples 1mW of power 

    # Now, how to define power in a wave?
    # Ryan has http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf

    # If I define power in wave by
    # time_integrated_power = deltaT * V^2/R
    # t_i_p has the desired properties for characterizing impulses.

    # Let's calculate the power for a 'top hat' function over a 100 ns period
    total_time = 100 #ns
    dt_low = 10 #ns
    dt_med = 1 #ns
    dt_high = 0.1 #ns

    times_low = [dt_low * i for i in range(total_time/dt_low)]
    times_med = [dt_med * i for i in range(total_time/dt_med)]
    times_high = [dt_high * i for i in range(int(total_time/dt_high))]

    # Now suppose I'm sending in Z mW of power...
    Z_mW = 1 # to start with
    Z_W = Z_mW*1e-3
    # in a top hat pulse of 10 ns, which makes the sums easy
    # If you deliver a pulse of energy, power over a short time
    # What does it mean to talk about power if you can aritrarily increase the time?
    # (The zero padding problem?)

    # Talk about energy - which is Ryan's normalization
    # mW/Hz = mWs = mJ
    
    # So count the mJ of power in the wave..
    # sum(V**2/R)
    
    return 0


if __name__ == '__main__':
    main()
