# Quantum Tomography
This program allows to perform quantum tomography on quadrature marginal distributions to obtain the Wigner function. It also allows to calculate the density matrix in quadrature and the number representation from the generated Wigner function. In addition, it offers the option to simulate quadrature measurements to use as input.
These functionalities are presented in a user-friendly graphical interface and take advantage of multiprocessing techniques to provide fast results.

[//]: <img width="80%" alt="Screenshot 2022-10-30 at 7 52 04 PM" src="https://user-images.githubusercontent.com/44838230/198908531-4659e160-90d2-48f7-aea1-8a06d42f3ec7.png">
<img width="80%" alt="Main window of the tomography displaying the Wigner function of a coherent state"  src="https://github.com/amartinez1224/quantum-tomography/assets/44838230/af51d5f4-d2c4-4e04-92bb-a74979101869">

### Density matrix
<img width="40%" alt="Window displaying the Density matrix in Fock representation of a coherent state"  src="https://github.com/amartinez1224/quantum-tomography/assets/44838230/0c67ef93-97be-43c2-b2b6-678db7113459">
<img width="59%" alt="Window displaying the Density matrix in quadrature representation of a coherent state"  src="https://github.com/amartinez1224/quantum-tomography/assets/44838230/3df5bb2d-bb73-4964-bca7-247305fc8cb2">

### Data simulation (marginal distributions)
<img width="40%" alt="Data simulation window" src="https://github.com/amartinez1224/quantum-tomography/assets/44838230/d6a31423-c3d0-484b-b2a8-a203d882ef27">
<img width="49%" alt="Squeezed state simulated data" src="https://github.com/amartinez1224/quantum-tomography/assets/44838230/35233e81-7b70-4018-9f57-48b58ac6afc2">






# Installation

## Install using binary files
You can use the installers for your system to easily install the software. 
- **[MacOS (Ventura 13.0 - Intel)](https://github.com/amartinez1224/quantum-tomography/releases/download/v3.0/QuantumTomography_MacOS_Intel.app.zip)**
- **[MacOS (Ventura 13.0 - Silicon)](https://github.com/amartinez1224/quantum-tomography/releases/download/v3.0/QuantumTomography_MacOS_Silicon.app.zip)**
- **[Windows](https://github.com/amartinez1224/quantum-tomography/releases/download/v3.0/QuantumTomography_Windows.zip)**

If you get a Launch error or cannot install the application using the binaries, try to clone the repository following the instructions below.

## Clone repo
As an alternative you can directly clone the repo and run the `gui.py` file (works in Windows, Linux and MacOS).
```
$ git clone https://github.com/amartinez1224/quantum-tomography.git
$ cd quantum-tomography
$ pip install -r requirements.txt
$ python gui.py
```

### Use a virtual environment (optional)
It's common practice to use a virtual environment to run Python applications. Virtual environments allow you to use an specific Python version (or the same one as the OS) and install specific versions of libraries that match those required by the application without creating version conflicts with the global libraries in your OS. You can find more information on virtual environments [here](https://docs.python.org/3/library/venv.html).

First clone the repo and navigate to the main folder
```
$ git clone https://github.com/amartinez1224/quantum-tomography.git
$ cd quantum-tomography
```
Once there create a virtual environment using *venv* run
```
$ python -m venv DIR_OF_NEW_ENV
```
Then activate the environment, install the requirements and run the `gui.py` file on the *venv*: (Mac)
```
$ source DIR_OF_NEW_ENV/bin/activate
$ pip install -r requirements.txt
$ python gui.py
```
Then activate the environment, install the requirements and run the `gui.py` file on the *venv*: (Windows)
```
$ .\DIR_OF_NEW_ENV\Scripts\activate
$ python -m pip install -r requirements.txt
$ python gui.py
```

To exit the *venv* run
```
deactivate
```

# Usage
To perform a tomography, select an input file using the button in the lower right corner. Then click on the *Tomography* button. 

To tune the tomography, use the values on the right and then regenerate the tomography by clicking on the *Tomography* button:
- X min, X max: set the range of the X axis.
- Y min, Y max: set the range of the Y axis.
- Density: changes the amount of points generated for the Wigner function. The size of the generated Wigner function function is *density*\**density*.
- Kc: the value for the cutoff frequency.
- Filter: the reconstruction filter (see [Filters](#filters) below).

You can also change the color of the graph and the angle from which it is viewed. For this you don't have to regenerate the tomography. In addition, the graph is fully interactive, so you can use the mouse to scroll.

At the bottom of the image you will find three options:
- Save data: saves the generated Wigner function in a JSON file.
- Density matrix quadrature: shows the density matrix in the quadrature space.
- Density matrix Fock: shows the density matrix in the Fock space. You can select the number of photons *n* and *m*.

## Filters
The Wigner function is obtained from the marginal distributions $pr(x,\theta)$ through the inverse Radon transform

$$W(q,p)=\frac{1}{2\pi^2}\int_0^{\pi}d\theta\int dx\ pr(x,\theta)\ K(q\cos\theta+p\sin\theta-x),\qquad K(x)=\frac{1}{2}\int_{-k_c}^{k_c}|k|\ e^{ikx}\ dk .$$

The kernel $K$ is a ramp filter $|k|$ cut at the frequency $k_c$ (the *Kc* value in the GUI). `tomo.py` has two ways of computing this, and you choose one with the *Filter* menu.

### Filtered back-projection (`tomo.tomo_fbp`, default)
Each marginal distribution is filtered **once** in Fourier space and then back-projected onto the $(q,p)$ grid:
1. The voltage axis is reversed if it is decreasing, and each projection is zero-padded to a power of two $\geq 4n$, where $n$ is the number of voltage points. The padding prevents circular wrap-around in the FFT convolution.
2. The FFT of each projection is multiplied by the filter $H(k)=|k|\ w(|k|/k_c)$ for $|k|\le k_c$ and $H(k)=0$ above $k_c$, and then transformed back.
3. For every angle, the filtered projection is linearly interpolated at $x=q\cos\theta+p\sin\theta$ and added to $W$. Points outside the measured voltage range contribute 0.
4. The result is normalized by $1/(2N)$, where $N$ is the number of angles.

The available windows $w(u)$, with $u=|k|/k_c$, are listed in `tomo.FILTERS`:

| Filter | Window $w(u)$ | Notes |
|---|---|---|
| `shepp-logan` (default) | $\mathrm{sinc}(u/2)=\dfrac{\sin(\pi u/2)}{\pi u/2}$ | Mild smoothing; keeps most of the resolution. |
| `ramlak` | $1$ | Sharp cutoff. Same as the original kernel; gives the most Gibbs ringing. |
| `hamming` | $0.54+0.46\cos(\pi u)$ | Strong smoothing. |
| `hann` | $0.5\,(1+\cos(\pi u))$ | Strongest smoothing; reaches 0 at $k_c$. |

Smoother windows reduce ringing and noise but blur fine features. A larger *Kc* gives more resolution and more noise. The default is $k_c=12$.

**Requirement:** the angles must be uniformly sampled over $[0,\pi)$ or $[0,2\pi)$, without repeating the endpoint.

### Kernel (legacy) (`tomo.tomo`)
This is the original method. For every $(q,p)$ point it evaluates the closed-form real-space kernel and sums over all angles and voltages, running in parallel with `multiprocessing`:

$$K(x)=\frac{\cos(k_c x)+k_c x\sin(k_c x)-1}{x^2}.$$

When $|k_c x|<0.01$, the code uses a 5-term Taylor expansion instead, $K\approx\frac{k_c^2}{2}\left(1-\frac{(k_cx)^2}{4}+\frac{(k_cx)^4}{72}-\frac{(k_cx)^6}{2880}+\frac{(k_cx)^8}{201600}\right)$, to avoid the $0/0$ cancellation (functions `K`, `Kor` and `Kcomp`). This method is equivalent to the `ramlak` filter but costs $O(D^2 N n)$ operations, compared with about $O(N n\log n + N D^2)$ for filtered back-projection ($D$ = *Density*). It is much slower.

The selected filter and *Kc* are stored in the JSON written by *Save data*. The file `docs/filters.tex` has the full derivation.

## Input file
The input file has to be in the following JSON format:
```
{
  "phi": [],
  "x": [],
  "pr": [
    []
  ]
}
```
where *phi* is a one dimensional array that contains the angle of each marginal distribution, *x* another one dimensional array that contains the voltage range of the distributions and *pr* is a two dimensional array of size *phi*\**x* that contains the distribution.

## Samples
The folder *SampleStates* contains samples of simulated marginal distributions for the vacuum, coherent, squeezed and thermal states.

**sim_vacuum**
- Voltage range: -8V, 8V

**sim_thermal.json**:
- Voltage range: -8V, 8V

**sim_coherent.json**
- Voltage range: -8V, 8V
- Alpha: 3
- Angle: 180°

**sim_squeezed.json**
- Voltage range: -8V, 8V
- Alpha: 0
- Angle: 0°
- Eta: 0°

**sim_squeezed_vacuum.json**
- Voltage range: -8V, 8V
- Alpha: 3
- Angle: 90°
- Eta: 45°


## Simulate data
You can simulate your own data using the *simulate data* button or by running the `simulateStates.py` script. The script has an inline interface for guidance. You can also preview the simulated data in a graphic.

## Read scope data
Data taken from a scope can be processed with the help of the `readScopeData.py` script to fit the input file format.

## Wigner functions online gallery
Coming soon [here](https://amartinez1224.github.io/wigner-functions-gallery/).

# Documents
## Monograph
You can find the monograph that lead to this development [here](https://github.com/amartinez1224/quantum-tomography/blob/567bfe7e2c01903786542dab36ff7c925361e157/docs/monograph.pdf).
## Lab guides
- [Lab guide 1](https://github.com/amartinez1224/quantum-tomography/blob/567bfe7e2c01903786542dab36ff7c925361e157/docs/Lab_guide_1.pdf)
- [Lab guide 2](https://github.com/amartinez1224/quantum-tomography/blob/567bfe7e2c01903786542dab36ff7c925361e157/docs/Lab_guide_2.pdf)

# Contact
If you find any bugs or are further interested in the project, please do not hesitate to contact the team or the [Quantum Optics Group ](https://opticacuantica.uniandes.edu.co/en/) at [Uniandes](https://uniandes.edu.co/en).
- Prof. Alejandra Valencia, Ph.D ([ac.valencia@uniandes.edu.co](mailto:ac.valencia@uniandes.edu.co))
- Andrés Martínez S. ([a.martinez@uniandes.edu.co](mailto:a.martinez@uniandes.edu.co))
- Juan-Rafael Álvarez, Ph.D ([jr.alvarez2101@uniandes.edu.co](mailto:jr.alvarez2101@uniandes.edu.co))
- Jose Ricardo Mejía, Lab Analyst ([jr.mejia1228@uniandes.edu.co](mailto:jr.mejia12281@uniandes.edu.co))
