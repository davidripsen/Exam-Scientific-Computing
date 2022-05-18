function f = CSTR_3D( [beta; k0; Ea])
% Init
% - se en masse konstanter i paperet
% gang k0 med 60 for at få i minutter
% Lav en F(tid) funktion hardcodet step-wise der passer med grafen. Husk at
% omregn fra mL -> L.
        % Mål: Temperatur grafens dynamik.
f = zeros(3,1);

kT = k0*exp(-Ea / (R*T));
r = kt * Ca * Cb

Ra = -r;
Rb = -2*r;
RT = beta*r;

Cdota = F/V * (Ca_in - Ca) + Ra;
Cdotb = F/V * (Cb_in - Cb) + Rb;
Tdot =  F/V * (T_in - T) + RT;

f(1) = Cdota;
f(2) = Cdotb;
f(3) = Tdot;
end