---
title: GSoC 2021 - Final Report
tags: GSoC
excerpt: false
layout: article
---

This is the final report for my GSoC project with SymPy.

<h2>Project Synopsis</h2>

My project consisted of 2 parts -

* Rational Riccati Solver - I implemented an ODE solver for Riccati equations with rational solutions. This solver is quite important as it can be used as a subroutine for other solvers like the Kovacic solver (for second order linear homogeneous ODEs), Lie group solver (for 1st order nonlinear ODEs), etc.

* ODE test suite - I created a test suite of all ODEs present in _Differential Gleichungen_ by Kamke in SymPy. This book consists of 1940 ODEs of various orders and linearities including systems of ODEs as well. It will be useful to check how many of these SymPy can solve and improve/add solvers based on the errors we encounter in these set of ODEs.

Additionally, I worked for sometime on `dsubs`, a function to facilitate substitutions in differential equations.

<h2>Work Accomplished</h2>

[#21459](https://github.com/sympy/sympy/pull/21459) - Rational Riccati Solver

[Test Suite Repository](https://github.com/naveensaigit/kamke_test_suite/) - A repository which contains all examples from Kamke. The repo has a way to run all the examples and generate a report in both human/machine readable format.

<h2>Work in Progress</h2>

[#20979](https://github.com/sympy/sympy/pull/20979) - PR for `dsubs`, a function for substitutions in differential equations.

<h2>Examples</h2>

<h3>Riccati ODEs</h3>

Any Riccati ODE with atleast 1 rational particular solution can now be solved in SymPy. Some examples are -
``` python
>>> from sympy import Symbol, Function, dsolve, checkodesol
>>> f = Function('f')
>>> x = Symbol('x')

>>> eq = -x**4*f(x)**2 + x**3*f(x).diff(x) + x**2*f(x) + 20
>>> sol = dsolve(eq, hint="1st_rational_riccati")
>>> sol
Eq(f(x), (4*C1 - 5*x**9 - 4)/(x**2*(C1 + x**9 - 1)))
>>> checkodesol(eq, sol)
(True, 0)

>>> eq = f(x).diff(x) + f(x)**2 - (2*x + 1/x)*f(x) + x**2
>>> sol = dsolve(eq, hint="1st_rational_riccati")
>>> sol
Eq(f(x), x*(C1 + x**2 + 1)/(C1 + x**2 - 1))
>>> checkodesol(eq, sol)
(True, 0)
```

<h3>Kamke Test Suite</h3>

To view the results of the test suite run against the latest version of SymPy, go to the [home page](https://naveensaigit.github.io/kamke_test_suite). To use the Kamke Test Suite, please clone the [repository](https://github.com/naveensaigit/kamke_test_suite). To run the entire test suite, run -
```
python test_kamke.py
```

This will generate a folder with the results of the run in JSON format. To generate HTML pages from these files, run -

```
python test_kamke.py --html
```

There are different ways to run the test suite using command line arguments. To view all the different ways to run the test suite, please refer to the [README](https://www.github.com/naveensaigit/kamke_test_suite/tree/gh-pages/README.md) of the repository.

<h3>Substitutions in ODEs</h3>

A basic transformation of `x -> t^2` can be achieved by -

``` python
>>> from sympy import symbols, dsubs, Function, cos
>>> f, g, h = symbols('f g h', cls=Function)
>>> x, a, t = symbols('x a t')
>>> eq = f(x).diff(x) + f(x)
>>> dsubs(eq, {x: t**2})
f(t**2) + Derivative(f(t**2), t)/(2*t)
```

Functions in the equation can be transformed too

``` python
>>> eq = f(x).diff(x, 2) + cos(x)
>>> dsubs(eq, {x: t**2, f(x): g(t)})
cos(t**2) + (Derivative(g(t), (t, 2))/(2*t)
- Derivative(g(t), t)/(2*t**2))/(2*t)

>>> dsubs(eq, {x: t, f(x): g(t)**2})
2*g(t)*Derivative(g(t), (t, 2)) + cos(t)
+ 2*Derivative(g(t), t)**2
```

To apply a transformation with number of new symbols greater than number of old symbols, a list of all new variables must be given. For example, a transformation `x = a t` where `t` is the new variable and `a` is a symbolic constant, would be done as -

``` python
>>> dsubs(eq, {x: a*t, f(x): g(t)}, [t])
cos(a*t) + Derivative(g(t), (t, 2))/a**2
```

We can do it similarly for functions as well.

``` python
>>> dsubs(eq, {x: t, f(x): h(a)*g(t)**2}, [t, g(t)])
2*g(t)*h(a)*Derivative(g(t), (t, 2))
+ 2*h(a)*Derivative(g(t), t)**2 + cos(t)
```

Let us see how `dsubs` works for some typical cases using an example of a homogeneous ODE which can be solved using the substitution `v = y(x)/x`.

``` python
>>> from sympy import symbols, Function, dsolve
>>> x, t = symbols('x t')
>>> y, v = symbols('y v', cls=Function)
>>> eq = y(x).diff(x) - (x**2 + y(x)**2)/(x*y(x))
```

Transform the homogeneous equation to a linear equation

``` python
>>> homeq = dsubs(eq, {x:t, y(x): t*v(t)}, [t, v(t)]).simplify()
>>> homeq
t*Derivative(v(t), t) - 1/v(t)
```

Find the solutions for the transformed equation

``` python
>>> solt1, solt2 = dsolve(homeq)
>>> solt1
Eq(v(t), -sqrt(C1 + 2*log(t)))
```

Re-substitute `v = y(x)/x` to get the solutions to the original equation

``` python
>>> sol = dsubs(solt1, {v(t): y(x)/x, t: x})
>>> sol
Eq(y(x)/x, -sqrt(C1 + 2*log(x)))
```

Let us look at another example - The 2nd order Cauchy Euler equation.

``` python
>>> x, a, b, t = symbols('x a b t')
>>> y, phi = symbols('y phi', cls=Function)
>>> eq = x**2*y(x).diff(x, 2) + a*x*y(x).diff(x) + b*y(x)
```

Transform the equation using the substitutions `x = e^t` and `y(x) = phi(t)`

``` python
>>> from sympy import exp, log, logcombine
>>> eqtrans = dsubs(eq, {x: exp(t), y(x): phi(t)}, [t, phi(t)]).simplify()
>>> eqtrans
a*Derivative(phi(t), t) + b*phi(t) - Derivative(phi(t), t) + Derivative(phi(t), (t, 2))
```

The equation is now transformed to a constant-coefficient equation.

``` python
>>> solt = dsolve(eqtrans)
>>> solt
Eq(phi(t), C1*exp(t*(-a - sqrt(a**2 - 2*a - 4*b + 1) + 1)/2) + C2*exp(t*(-a + sqrt(a**2 - 2*a - 4*b + 1) + 1)/2))
```

Re-substitute for `phi(t)` and `t` to get the solution for the original equation.

``` python
>>> sol = logcombine(dsubs(solt, {phi(t): y(x), t: log(x)}), force=True)
>>> sol
Eq(y(x), C1*x**(-a/2 - sqrt(a**2 - 2*a - 4*b + 1)/2 + 1/2) + C2*x**(-a/2 + sqrt(a**2 - 2*a - 4*b + 1)/2 + 1/2))
```

<h2>Future Work</h2>

[#21770](https://github.com/sympy/sympy/pull/21770/) - Solver for Part-1 of Kovacic Algorithm. I stopped working on this because it was not clear how to implement parts 2 and 3 of the Kovacic algorithm with the Riccati solver I added.

[#21704](https://github.com/sympy/sympy/pull/21704/) - I initally wrote a parser for Maple to parse the examples from Maple to SymPy. However, most examples could be parsed using `sympify`, so I didn't work on this during the project. However, it will be useful to have a parser for Maple just like how we have one for Mathematica.

<b>Test Suite Repository</b> - Currently, the test suite only contains ODEs from Kamke. There are about 2300 more ODEs from the book _Ordinary Differential Equations and Their Solutions_ by Murphy that can also be added to the test suite. They are available in Mathematica syntax [here](https://drive.google.com/file/d/1kwDlhBQwUh2H8KeBmovUzkMY3SCFCLPa/view?usp=sharing). We would have to use the `MmaTranslator` function in Maple to convert these ODEs into Maple syntax so that they can be easily parsed with `sympify`.

<h2>Conclusion</h2>

I have definitely learnt a lot during this summer! The journey from being new to open source to completing an entire project gives me a lot of satisfaction! I would like to sincerely thank my mentors [Oscar Benjamin](https://github.com/oscarbenjamin) and [Aaron Meurer](https://github.com/asmeurer) for helping me throughout the project by continuously reviewing my code and suggesting changes to make it more efficient! I would also like to thank [Nijso](https://github.com/bigfooted) for suggesting the idea for this project and his feedback during the initial days. Lastly, thanks to [Nasser M. Abbasi](https://github.com/nasser1) for his [website](https://www.12000.org) where he published a list of all ODEs from Kamke in Maple and Mathematica syntax. The test suite could not have been created without his work.