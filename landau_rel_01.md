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










