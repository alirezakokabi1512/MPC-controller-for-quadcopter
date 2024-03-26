function [Aa,Ba,Ca]= augment_mimo(Ad,Bd,Cd,n_states,n_inputs,n_outputs)

Aa=[Ad zeros(n_states,n_inputs); Cd*Ad eye(n_outputs)];
Ba=[Bd;Cd*Bd];
Ca=[zeros(n_inputs,n_states) eye(n_inputs)];




