# Notes on runs with k=0.1

We use n=20 for the frequency perturbation, and $\alpha = 0.01$, with $\theta = 0.001$.
This corresponds to $\hat{k} = 0.1$, so the frequency should be equal to $1.015$
($\omega = 1 + \frac{3}{2}\hat{k}^2 = 1 + 1.5*0.01$). The imaginary part should still be negligible.


## run: 
```
run_name = "norel_norm_undamped_rel_t400_L39.738_N85_n20_J3522_M16001_o5_Th3_alp2"
(N, L, J, dx, order) = par_grid = (800000, 39.738, 3522, 0.011282793867120954, 5)
(t_i, t_f, M, M_g, dt) = par_evolv = (0.0, 400.0, 16001, 1601, 0.025)
(θ, nm, k) = par_f = (0.001, 20, 0.1581152878146758)
```
![some values](Images/norel_norm_undamped_rel_t400_L39.738_N85_n20_J3522_M16001_o5_Th3_alp2_total_run.png)

The current is growing linearly, but very slowly. Is it just the error?
The energy is not growing, so the system seems stable.

![Energy conservation](Images/norel_norm_undamped_rel_t400_L39.738_N85_n20_J3522_M16001_o5_Th3_alp2energy_conservation.png)


![Energy fit](Images/norel_norm_undamped_rel_t400_L39.738_N85_n20_J3522_M16001_o5_Th3_alp2_energy_fit.png)

The fit seems to be ok, the parameters are: 
```
@. model_e1(x,p) = p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x) + p[6]*cos(p[7]*x + p[8])*exp(-p[9]*x)
@. model_e2(x,p) = p[1] + p[2]*(cos(p[3]*x + p[4])^2-p[6])*exp(-p[5]*x)

```
```
p_e_1 = [0.0020502155276045917; 5.039024594434988e-5; 2.027635371154509; 6.345690430337964; -0.0003358830579111077; 0.0019667201803559324; 1.99859474123304; 0.017465375859802772; 0.00020549260871791195]
```

So, there are two frequencies, $\omega_2 = 2.02763$ and $\omega_1 = 1.99859$, the one close to 2 seems to be the dominant. I don't understand this. 

## Temperature fit.

```
@. model_tl001(x,p) = p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x) + p[6]*cos(p[7]*x + p[8])*exp(-p[9]*x)
```

![Temperature fit](Images/norel_norm_undamped_rel_t400_L39.738_N85_n20_J3522_M16001_o5_Th3_alp2temperature_fit.png)

p_t_1 = [0.0011075629178367128; -2.539600079603148e-6; 2.027635613940527; 0.06243459044499647; -0.0003363494104947458; -9.914225002085228e-5; 1.9985947350845683; 0.017467277743187723; 0.0002055007935678624]

The fit seems to be very good also, and the frequencies are the same as before, the relative contributions are different, but in both cases the lower frequency is the dominant.

## More sampling points.

We made another runs, in particular one with more sampling points.

### run: norel_norm_undamped_rel_t400_L4_N85_n2_J80_M16001_o5_Th3_alp2_ave.jld2 

```
run_name = "norel_norm_undamped_rel_t400_L4_N85_n2_J80_M16001_o5_Th3_alp2"
(N, L, J, dx, order) = par_grid = (800000, 4, 80, 0.05, 5)
(t_i, t_f, M, M_g, dt) = par_evolv = (0.0, 400.0, 16001, 8001, 0.025)
(θ, nm, k) = par_f = (0.001, 2, 1.5707963267948966)
```

**Notice that here L=4, so the wavelength is a bit different.**

![averages](Images/norel_norm_undamped_rel_t400_L4_N85_n2_J80_M16001_o5_Th3_alp2_total_run.png)

![Energy Conservation](Images/norel_norm_undamped_rel_t400_L4_N85_n2_J80_M16001_o5_Th3_alp2_energy_conservation.png
)

![Energy fit](Images/norel_norm_undamped_rel_t400_L4_N85_n2_J80_M16001_o5_Th3_alp2_energy_fit.png)

```
@. model_e1(x,p) =  p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x) + p[6]*cos(p[7]*x + p[8])*exp(-p[9]*x)
```

```
p_e = [6.482617995271001e-6; 3.840406038532317e-6; 2.012904971545909; 0.18051577914897698; -0.0006297953349135141; -2.7304791719562594e-6;  2.0033087748744984; 2.912384672766882;  9.485628705397978e-5]
```

The temperature fit is also good:

![Temperature Fit](Images/norel_norm_undamped_rel_t400_L4_N85_n2_J80_M16001_o5_Th3_alp2temperature_fit.png)

```
#model_tl001(x,p) = p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x) + p[6]*cos(p[7]*x + p[8])*exp(-p[9]*x)
fit_T_1.param = [0.0010060110877331346; -1.9068792348543928e-6;  2.012907059060496;  0.17919945844298754; -0.0006247012685843706; -1.3547776221701724e-6;  2.003298278286806; -0.22784260925642166;  9.178632138800702e-5]
```

The temperature spectrum shows several exited modes:

![Temperature spectrum](Images/norel_norm_undamped_rel_t400_L4_N85_n2_J80_M16001_o5_Th3_alp2_temperature_spectrum.png)

For this case we have:

$\omega_{energy} = (1.00645, 1.00165)$ (the dominant first)

$\omega_{temp} = (1.00645,1.00165)$




## run with 8e6 points, same as before but the time sampling is as the first one.

### run norel_norm_undamped_rel_t400_L4_N86_n2_J80_M16001_o5_Th3_alp2_ave.jld2

```
run_name = "norel_norm_undamped_rel_t400_L4_N86_n2_J80_M16001_o5_Th3_alp2"
(N, L, J, dx, order) = par_grid = (8000000, 4, 80, 0.05, 5)
(t_i, t_f, M, M_g, dt) = par_evolv = (0.0, 400.0, 16001, 1601, 0.025)
(θ, nm, k) = par_f = (0.001, 2, 1.5707963267948966)
```

![averages](Images/norel_norm_undamped_rel_t400_L4_N86_n2_J80_M16001_o5_Th3_alp2_total_run.png)

![Energy Conservation](Images/norel_norm_undamped_rel_t400_L4_N86_n2_J80_M16001_o5_Th3_alp2_energy_conservation.png)

![Energy Fit](Images/norel_norm_undamped_rel_t400_L4_N86_n2_J80_M16001_o5_Th3_alp2_energy_fit.png)

```
@. model_e1(x,p) =  p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x) + p[6]*cos(p[7]*x + p[8])*exp(-p[9]*x)
```
```
p_e = [5.711450907054283e-6; 1.3608018340426746e-5; 2.012575349105527; 0.1350197841144081; -0.0002551407684421814; -7.968990739872489e-6; 2.012123942134308; 0.2307122079967403; -0.00044981767181786797]
```

**Notice that here there is not the large oscillation, but rather a small one. We should repeat it with more sampling points...**


The temperature fit is also good:

![Temperature Fit](Images/norel_norm_undamped_rel_t400_L4_N86_n2_J80_M16001_o5_Th3_alp2temperature_fit.png)

```
#model_tl001(x,p) = p[1] + p[2]*cos(p[3]*x + p[4])*exp(-p[5]*x) + p[6]*cos(p[7]*x + p[8])*exp(-p[9]*x)
fit_T_1.param = [0.0010034058638794707; -2.842954241819258e-6;  2.013166365765038;  0.01628772065320813; -1.3037764963152664e-5;  4.122742038753967e-8;  2.0063015016254666;  1.7863027971973944; -0.0017475811770664323]
```

The temperature spectrum shows several exited modes:

![Temperature spectrum](Images/norel_norm_undamped_rel_t400_L4_N86_n2_J80_M16001_o5_Th3_alp2_temperature_spectrum.png)


For this case we have:

$\omega_{energy} = (1.00628, 1.0060)$ (the dominant first)

$\omega_{temp} = (1.00658,1.00315)$









