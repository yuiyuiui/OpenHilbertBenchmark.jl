Test (vectors input) Open Hilbert transform algrorithm in OpenHilbert.jl

# 1. Accuracy test
Points gap $h=2^{-5}\sim 2^{-4}$ 


|Number of Points $\sim$ accuracy (max error)|FFT|FIR|ASY|AAA|
|---|---|---|---|---|
|Type1: Schwartz Function|$10^6\sim10^{-5.5}$ | $10^2\sim10^{-16}$| $10^{2.5}\sim10^{-16}$|$10^{5.7}\sim10^{-12}$ |
|Type2: Rantional Function (order number = 1)|$10^6\sim10^{-4}$ | $10^6\sim10^{-4}$|$10^4\sim 10^{-15}$ |$10^6\sim10^{-5}$ |
|Type3: Rational Function (order number > 1)|$10^6\sim10^{-4.5}$ | $10^6\sim10^{-4.5}$ |$10^5\sim10^{-15.5}$ |$10^6\sim10^{-6.5}$ |
|Type4: Rational (order number > 1) + Schwartz|$10^6\sim10^{-4.5}$ | $10^6\sim10^{-4.5}$ |$10^5\sim10^{-15.5}$ |$10^{5.5}\sim10^{-11}$ |
|Type5: Log-Rational |$10^6\sim10^{-5.5}$ | $10^6\sim10^{-5.5}$ |$10^6\sim10^{-10}$ |$10^{5.7}\sim10^{-6}$|
|Type6: Log-Rational + Rational (order number > 1) + Schwartz|$10^6\sim10^{-4.5}$ | $10^6\sim10^{-4.5}$ |$10^6\sim10^{-10}$ |$10^{5.7}\sim10^{-6}$ |
|Type7: d-Rational|$10^6\sim10^{-2}$ | $10^6\sim10^{-2}$ |$10^2\sim10^{-15.5}$ |$10^{5.5}\sim10^{-2}$ |
|Type8: d-Rational + Rational (order number = 1) + Schwartz|$10^6\sim10^{-2}$ |$10^6\sim10^{-2}$ | $10^{5.7}\sim10^{-11}$ |$10^{5.5}\sim10^{-2}$ |
|Type9: d-Rational + Log-Rational + Rational (order number > 1) + Schwartz|$10^6\sim10^{-2}$ | $10^6\sim10^{-2}$ |$10^{6.3}\sim10^{-10}$ |$10^{5.7}\sim10^{-2}$ |

# 2. Speed test
|Objects/Methods|FFT|FIR|ASY|AAA|
|---|---|---|---|---|
|Type1: Schwartz Function|
|Type2: Rantional Function (order number = 1)|
|Type3: Rational Function (order number > 1)|
|Type4: Rational (order number > 1) + Schwartz|
|Type5: Log-Rational |
|Type6: Log-Rational + Rational (order number > 1) + Schwartz|
|Type7: d-Rational|
|Type8: d-Rational + Rational (order number = 1) + Schwartz|
|Type9: d-Rational + Log-Rational + Rational (order number > 1) + Schwartz|

# 3. To Do

### 1. Add correctness test of hilbert transform by grid integral
### 2. Compaare with current algorithms
### 3. Add a analytic calculated hilbert transform or a $O(1/x\ln x )$ function
