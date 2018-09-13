### Dependencies:
  
numpy 
matplotlib.pyplot
seaborn
pyfftw

# Sample Normalized Output
## Initial Conditions
```python
y = (i * dx) - (L/3.8) # or 3.8
x = (j*dx) - (L/3.8) # or 3.8
kx = ((0.7 * math.pi) / L)
ky = ((1.0 * math.pi) / L)
```
N = 128, L = 40.0, U = 80.0, dt = 0.01 tf = 50.0 undersampled = 2 twice
<p align="center">
<img src="https://github.com/mauckc/2D-Quantum-Free-Particle/blob/master/media/sample.gif"/>
</p>

## Usage:
  
  python dev-2Dparticle.py

### Sample output:
  
```
printing field number 128 output number64
at time: 1.27
(256L, 256L)
Phi sum:1019.97323996
Phi sum1:1019.97323996
Phi max:1.99920403047
Phi avg:0.015563556518
(256L, 256L)
Phi sum:22330354042.8
Phi sum1:22330354042.8
Phi max:43704420.4754
Phi avg:340734.162029

printing field number 130 output number65
at time: 1.29
(256L, 256L)
Phi sum:21871780457.9
Phi sum1:21871780457.9
Phi max:42917529.9614
Phi avg:333736.884429
(256L, 256L)
Phi sum:1011.62579958
Phi sum1:1011.62579958
Phi max:2.0
Phi avg:0.015436184686

printing field number 132 output number66
at time: 1.31
(256L, 256L)
Phi sum:21526789148.8
Phi sum1:21526789148.8
Phi max:42346999.5531
Phi avg:328472.734815
(256L, 256L)
```
