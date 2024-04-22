def poly_circuit(volts, ph, x1, x2, xM, r1, r2, rfe, f_w, sync):
    '''
    This function computes the performance of an induction motor given its equivalent
    circuit parameters aka motor constants. 
    Args:
        volts: Phase-voltage value [volts]
        ph:    No. of phases
        x1:    Primary leakage reactance [ohms]
        x2:    Secondary leakage reactance [ohms]
        xM:    Mutual reactance [ohms]
        r1:    Primary resistance [ohms] 
        r2:    Secondary resistance [ohms] 
        rfe:   Resistance for iron losses [ohms]
        f_w:   Friction and windage loss [watts]
        sync:  Synchronous speed of machine [rpm]
    Returns:
        perf:  Dictionary containing performance over range of operating points
    '''
    import numpy as np

    # series core-loss resistance
    rM = xM**2 / (rfe * (1 + (xM/rfe)**2))

    # F-CONSTANTS 
    F1 = x1 + x2 * xM / (xM + x2) - r1 * rM / (xM + x2)
    F2 = r2 * (xM + x1) / (xM + x2)
    F3 = r2 * (r1 + rM) / (xM + x2)
    F4 = rM * (x1 + x2) / (xM + x2) + r1
    F5 = volts * rM / (xM + x2)
    F6 = volts * r2 / (xM + x2)
    F7 = volts * np.sqrt(rM**2 + xM**2) / (xM + x2)
    A10 = volts * x2 / (xM + x2)

    # points calculated: no-load plus operation
    slip_points = [0.0001, 0.0005, 0.001, 0.005] + [abs(round(i, 2)) for i in np.arange(0.01, 1.01, 0.01)] 
    
    # dictionary to hold expected performance data
    perf = {'slip': slip_points,
            'current': [],
            'pri_loss': [],
            'sec_loss': [],
            'core_loss': [],
            'input': [],
            'torque': [],
            'eff': [],
            'pf': []}
    
    for s in slip_points: 
        
        part1 = np.sqrt((F5 + F6 / s)**2 + volts**2)
        part2 = np.sqrt((F3 / s - F1)**2 + (F2 / s + F4)**2)
        
        I = part1 / part2 # primary current
        I2 = F7 / part2 # secondary current
        IM = np.sqrt(((F6/s)**2 + A10**2) / part2**2) # magnetizing current

        pri_loss = I**2 * r1 * ph 
        sec_loss = I2**2 * r2 * ph
        core_loss = IM**2 * rM * ph
        
        if s != 1: 
            sec_output = sec_loss * (1 - s) / s  
            if sec_output > f_w:
                output = sec_output - f_w
            else: output = 0
        else: 
            sec_output = 0
            output = 0
        
        input_ = pri_loss + sec_loss + core_loss + sec_output 
        
        torque = 7.04 * (sec_loss - f_w) / (s * sync) 
            
        eff = output / input_ * 100
        
        pf = input_ / (volts * I * ph) * 100
        
        # add perf for this speed in results dict
        perf['current'].append(I)
        perf['pri_loss'].append(pri_loss)
        perf['sec_loss'].append(sec_loss)
        perf['core_loss'].append(core_loss)
        perf['input'].append(input_)
        perf['torque'].append(torque)
        perf['eff'].append(eff)
        perf['pf'].append(pf)

    return perf