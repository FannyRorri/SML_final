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