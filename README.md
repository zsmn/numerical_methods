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
euler 1.0 0 0.01 50  1-t+4*y
runge_kutta 1.0 0 0.025 40 -100*y+100*t+1
adam_bashforth 0.5 0.8292986 1.2140877 1.6489406 0 0.2 10 y-t*t+1 4
```

Output:

```
Euler
y(0) = 1.0
it = 50
h = 0.01
0   1.0
1   1.05000000000000
2   1.10190000000000
3   1.15577600000000
4   1.21170704000000
5   1.26977532160000
6   1.33006633446400
7   1.39266898784256
8   1.45767574735626
9   1.52518277725051
10   1.59529008834053
11   1.66810169187415
12   1.74372575954912
13   1.82227478993109
14   1.90386578152833
15   1.98862041278946
16   2.07666522930104
17   2.16813183847308
18   2.26315711201201
19   2.36188339649249
20   2.46445873235219
21   2.57103708164627
22   2.68177856491212
23   2.79684970750861
24   2.91642369580895
25   3.04068064364131
26   3.16980786938696
27   3.30400018416244
28   3.44346019152894
29   3.58839859919010
30   3.73903454315770
31   3.89559592488401
32   4.05831976187937
33   4.22745255235454
34   4.40325065444873
35   4.58598068062667
36   4.77591990785174
37   4.97335670416581
38   5.17859097233244
39   5.39193461122574
40   5.61371199567477
41   5.84426047550176
42   6.08393089452183
43   6.33308813030271
44   6.59211165551481
45   6.86139612173541
46   7.14135196660482
47   7.43240604526902
48   7.73500228707978
49   8.04960237856297
50   8.37668647370549
Runge Kutta
y(0) = 1.0
it = 40
h = 0.025
0   1.0
1   0.673437500000000
2   0.470471191406250
3   0.347649288177491
4   0.276796022802592
5   0.239641171036056
6   0.224337634343692
7   0.223203309769738
8   0.231256833678814
9   0.245268103088606
10   0.263142598096518
11   0.283522153453211
12   0.305526083879816
13   0.328583320015818
14   0.352323559072757
15   0.376506682836241
16   0.400976989651625
17   0.425633516727226
18   0.450410796002811
19   0.475266375533073
20   0.500172727884727
21   0.525112003237752
22   0.550072627099480
23   0.575047094134819
24   0.600030537603047
25   0.625019801726976
26   0.650012840182336
27   0.675008326055733
28   0.700005398926765
29   0.725003500866574
30   0.750002270093169
31   0.775001472013540
32   0.800000954508780
33   0.825000618939287
34   0.850000401343444
35   0.875000260246140
36   0.900000168753356
37   0.925000109426005
38   0.950000070955925
39   0.975000046010483
40   1.00000002983492
Adam Bashforth 
it = 10
h = 0.2
0   0.5
1   0.8292986
2   1.2140877
3   1.6489406
4   2.06331232416667
5   2.46774765690972
6   2.86303688111140
7   3.22380360859157
8   3.52119423255093
9   3.72506011232037
10   3.79773717936204
```

## Authors

* **Zilde Neto** - *Development* - [zsmn](https://github.com/zsmn)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
