Notes for the current system


Clamping [0, 1] in the forward function makes loss always 0
Clipping gradient in the train didn’t work

Forward types:

relu very unstable, after iter 4 it would go to crazy numbers (i.e -7458534.000)
sigmoid after iter 3 loss would go to 0

Changed to softplus with residual and it stabilised to low values but still oscillating

Could be changed:
lr, gradient clippping, dt
collocation point generation not be uniform
weight initialization
forward function can be changed


LR

with lr = 1e-4, we have constant oscillations between 0.001.. and 0.002.. and an overal negative residual for the PDE of -0.103 till -0.106
with lr = 1e-5, it starts in the 40.00 and then decreases fast to 0.00 we have oscillations between