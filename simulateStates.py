import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numba import jit
import json


@jit(nopython=True, parallel=True)
def gauss(x, amp, mu, sig):
    '''Numba gaussian function'''
    return amp*np.exp(-(x-mu)**2/(2.*sig**2))


@jit(nopython=True, parallel=True)
def seno(x, a1, a2, a3, a4):
    '''Numba sine function'''
    return a1*np.sin(a2*x+a3)+a4


def gausshist(n, bins):
    '''Fits a gaussian to a histogram'''
    # Get bins centers
    x = (np.abs(bins[1:]-bins[:-1])/2)+bins[:-1]
    # Perform gaussian fit over histogram
    arguments, variance = curve_fit(
        gauss, x, n, p0=[0.5, np.mean(x), np.std(x)])
    return arguments


def gen(mus, sig, density=100, v1=-12, v2=12):
    '''Generates sinograms using mu and sigma'''
    m = np.zeros((np.size(mus), density))
    y = np.linspace(v1, v2, density)
    for i in range(np.size(mus)):
        r = np.random.normal(mus[i], sig[i], 10000)
        n, bins = np.histogram(r, bins=density, density=True)
        arguments = gausshist(n, bins)
        m[i, :] = gauss(y, arguments[0], arguments[1], arguments[2])
    return m, y


def isAngle(s):
    '''Checks if s is a valid angle '''
    try:
        number = float(s)
        if number >= 0 and number < 360:
            return True
        return False
    except ValueError:
        return False


def isVoltage(s):
    '''Checks if s is a valid voltage '''
    try:
        number = float(s)
        if number > -120 and number < 120:
            return True
        return False
    except ValueError:
        return False


def isAlpha(s):
    '''Checks if s is a valid alpha '''
    try:
        number = float(s)
        if number >= 0 and number <= 10:
            return True
        return False
    except ValueError:
        return False


def saveData(data, filename):
    '''Saves data in a json file'''
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)


def musForAngle(angle, phi, alpha):
    '''Returns mus for a given angle'''
    angle = np.deg2rad(angle)
    a = np.cos(angle)
    b = np.sin(angle)
    value = a*np.cos(phi) - b*np.sin(phi)
    return value*alpha


def coherentWavefunction(x, x0, p0):
    '''Wavefunction of the coherent state D(x0,p0)|0> in the quadrature basis (vacuum variance 1/2)'''
    return np.pi**(-0.25)*np.exp(-0.5*(x-x0)**2 + 1j*p0*x - 0.5j*x0*p0)


def catMarginals(phi, x, alpha, angle, parity):
    '''Exact quadrature distributions pr(x,phi) of the cat state N(|beta> +- |-beta>), beta = alpha*exp(i*angle).
    alpha is the coherent amplitude, so the two components are displaced by sqrt(2)*alpha.'''
    s = 1 if parity == "even" else -1
    norm = 2*(1 + s*np.exp(-2*alpha**2))
    if norm < 1e-12:
        raise ValueError("The odd cat state needs alpha > 0")
    beta = alpha*np.exp(1j*np.deg2rad(angle))
    pr = np.empty((np.size(phi), np.size(x)))
    for i in range(np.size(phi)):
        # Rotating the quadrature by phi is equivalent to rotating beta by -phi
        b = beta*np.exp(-1j*phi[i])
        x0, p0 = np.sqrt(2)*b.real, np.sqrt(2)*b.imag
        psi = coherentWavefunction(x, x0, p0) + s*coherentWavefunction(x, -x0, -p0)
        pr[i, :] = np.abs(psi)**2/norm
    return pr


def sampleMarginals(pr, x, samples):
    '''Replaces each exact distribution by a histogram of a finite number of measurements (shot noise)'''
    rng = np.random.default_rng()
    dx = np.abs(x[1]-x[0])
    out = np.empty_like(pr)
    for i in range(np.shape(pr)[0]):
        p = np.clip(pr[i, :], 0, None)
        counts = rng.multinomial(samples, p/np.sum(p))
        out[i, :] = counts/(samples*dx)
    return out


def generateData(mode, v1, v2, alpha, angle, eta, parity="even", samples=0):
    data = {}
    name=""

    # Endpoint excluded so that phi=0 and phi=2pi are not counted twice
    phi = np.linspace(0, 2*np.pi, 180, endpoint=False)
    if mode == "1":
        data["state"] = "vacuum"
        name = f'sim_{data["state"]}'
        mus = np.zeros_like(phi)
        sig = np.zeros_like(phi)+(1/np.sqrt(2))
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    elif mode == "2":
        data["state"] = "coherent"
        data["alpha"] = alpha
        data["angle"] = angle
        name = f'sim_{data["state"]}_{alpha}_{angle}'
        sig = np.zeros_like(phi)+(1/np.sqrt(2))
        mus = musForAngle(angle, phi, alpha)
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    elif mode == "3":
        data["state"] = "thermal"
        name = f'sim_{data["state"]}'
        mus = np.zeros_like(phi)
        sig = np.zeros_like(phi)+3
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    elif mode == "4":
        data["state"] = "squeezed"
        data["alpha"] = alpha
        data["angle"] = angle
        data["eta"] = eta
        name = f'sim_{data["state"]}_{alpha}_{angle}_{eta}'
        sig = ((np.sin(2*(phi+(np.deg2rad(eta))))+1)*0.5 *
               (1/np.sqrt(2)))+(1/np.sqrt(2)) - 0.25
        mus = musForAngle(angle, phi, alpha)
        pr, x = gen(mus, sig, v1=v1, v2=v2)

    elif mode == "5":
        data["state"] = "cat"
        data["alpha"] = alpha
        data["angle"] = angle
        data["parity"] = parity
        data["samples"] = samples
        name = f'sim_{data["state"]}_{parity}_{alpha}_{angle}'
        x = np.linspace(v1, v2, 256)
        pr = catMarginals(phi, x, alpha, angle, parity)
        if samples > 0:
            pr = sampleMarginals(pr, x, samples)

    name = name.replace('.', '-')+".json"
    data.update({
        "voltage_range": [v1,v2],
        "phi": phi.tolist(),
        "x": x.tolist(),
        "pr": pr.tolist()
    })
    
    return data, name, phi, pr, x 


def main():
    print("\nWelcome to the quadrature distribution data generator!\n")
    v1, v2 = -8, 8
    print("To begin, please enter the voltage range of the simulated scope")
    while (voltage := input("Enter the voltage range (default -8,8): ")):
        if voltage == "":
            break
        elif voltage.count(',') == 1:
            v1I, v2I = voltage.split(',')
            print(v1I, v2I)
            if isVoltage(v1I) and isVoltage(v2I) and (v1I := float(v1I)) < (v2I := float(v2I)):
                v1 = v1I
                v2 = v2I
                break
            else:
                print("Invalid voltage range")
        else:
            print("Invalid voltage range")

    print(f'Voltage range set to {v1} to {v2} V\n')

    print("Now, please select the state to simulate")
    while (mode := input("Choose state (1: Vacuum, 2: Coherent, 3: Thermal, 4: Squeezed, 5: Cat): ")) not in ["1", "2", "3", "4", "5"]:
        print("Invalid mode")
    print(f'State {mode} selected\n')

    alpha = 3
    angle = 0
    if mode == "2" or mode == "4" or mode == "5":
        print("Now, please enter the alpha value")
        while (alphaI := input("Choose alpha (0-10) (default 3): ")):
            if alphaI == "":
                break
            elif isAlpha(alphaI):
                alpha = float(alphaI)
                break
            else:
                print("Invalid alpha")
        print(f'Alpha set to {alpha}\n')

        print("Now, please enter the angle value")

        while (angleI := input("Choose angle in degrees (default 0°): ")):
            if angleI == "":
                break
            elif isAngle(angleI):
                angle = float(angleI)
                break
            else:
                print("Invalid angle option")

    eta = 0
    if mode == "4":
        print("Now, please enter the eta value (squeezing direction)")
        while (etaI := input("Choose eta in degrees (default 0°): ")):
            if etaI == "":
                break
            elif isAngle(etaI):
                eta = float(etaI)
                break
            else:
                print("Invalid eta option")

    parity = "even"
    samples = 0
    if mode == "5":
        while (parityI := input("Choose parity (even/odd) (default even): ")):
            if parityI in ["even", "odd"]:
                parity = parityI
                break
            print("Invalid parity")
        while (samplesI := input("Measurements per angle (0 = exact distribution) (default 0): ")):
            if samplesI.isdigit():
                samples = int(samplesI)
                break
            print("Invalid number of measurements")

    data, name, phi, pr, x = generateData(mode, v1, v2, alpha, angle, eta, parity, samples)
    
    print("Data generated successfully!")
    print("Please select a filename to save the data (if the file already exists, it will be overwritten)")
    
    while (filename := input(f"Enter filename (default {name}): ")):
        if filename == "":
            break
        elif filename.count('.') == 1:
            name, ext = filename.split('.')
            if ext == "json":
                if name == "":
                    print("Invalid filename")
                else:
                    name = filename
                    break
            else:
                print("Invalid file extension (use .json)")
        else:
            print("Invalid file name (use .json and no additional dots)")
    
    saveData(data, name)

    print("\nThe data has been generated and saved to the current directory\n")

    view = input("If you want to see the generated data, enter 'y': ")
    if view == "y":
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i in range(np.size(phi)-1, 0, -1):
            ax.plot(x, np.zeros_like(x)+phi[i], pr[i, :], c='mediumvioletred')
        ax.set_ylabel(r"$\phi$")
        ax.set_xlabel(r"$x(v)$")
        ax.set_zlabel(r"$pr(x,\phi)$")
        ax.view_init(elev=0, azim=0)
        plt.show()
        plt.close()


if __name__ == "__main__":
    main()
