# Numerical Methods

This project was designed for the Computer Engineering Numerical Methods Class 2019.2

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

First clone this repository

```
git clone https://github.com/zsmn/numerical_methods.git
```

You also will need:

```
python >= 3
virtualenv
```

### Installing

Now after clone the repository you will need to run the RUNME file

```
source RUNME
```

And now after the end of installation of dependences, you can run the project

```
python3 methods.py
```

Now you can open the input.txt file and insert some methods to be calculated. You can see how you can do this on the following topic.

## How to insert a method in input.txt

Actually this project can do the following methods:

* Euler
* Aprimorated Euler
* Reverse Euler
* Runge-Kutta (order = 4)
* Adam Multon (order = 1 to order = 8)
* Adam Bashforth (order = 1 to order = 8)
* Inverse Formule (order = 1 to order = 6)

### Euler

```
euler y0 t0 steps function
```

### Aprimorated Euler

```
euler_aprimorado y0 t0 steps function
```

### Reverse Euler

```
euler_inverso y0 t0 steps function
```

### Runge Kutta

```
runge_kutta y0 t0 steps function
```

### Adam Multon

```
adam_multon y0 y1 ... yn-1 t0 steps function order
```

You can also use this method predicting the points using Euler, Aprimorated Euler, Reverse Euler or Runge Kutta methods

```
adam_multon_by_euler y0 t0 steps function order
```
```
adam_multon_by_euler_aprimorado y0 t0 steps function order
```
```
adam_multon_by_euler_inverso y0 t0 steps function order
```
```
adam_multon_by_runge_kutta y0 t0 steps function order
```

### Adam Bashforth

```
adam_bashforth y0 y1 ... yn t0 steps function order
```

You can also use this method predicting the points using Euler, Aprimorated Euler, Reverse Euler or Runge Kutta methods

```
adam_bashforth_by_euler y0 t0 steps function order
```
```
adam_bashforth_by_euler_aprimorado y0 t0 steps function order
```
```
adam_bashforth_by_euler_inverso y0 t0 steps function order
```
```
adam_bashforth_by_runge_kutta y0 t0 steps function order
```

### Inverse Formule

```
formula_inversa y0 y1 ... yn t0 steps function order
```

You can also use this method predicting the points using Euler, Aprimorated Euler, Reverse Euler or Runge Kutta methods

```
formula_inversa_by_euler y0 t0 steps function order
```
```
formula_inversa_by_euler_aprimorado y0 t0 steps function order
```
```
formula_inversa_by_euler_inverso y0 t0 steps function order
```
```
formula_inversa_by_runge_kutta y0 t0 steps function order
```

## Example of input and output

Input:

```
test
```

Output:

```
output
```

## Authors

* **Zilde Neto** - *Development* - [zsmn](https://github.com/zsmn)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
