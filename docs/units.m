kpc := 3.08568025e19 * m;
c := 3e8 * m/s;

eta := 1.22e-2;
M := (2.5e-22) * ev / c^2;
H0 := 1./(13.7e9 * 3600*24*365.25) / s;
hbar := 6.58e-16*ev*s;

solve(eta = M*(q^2)*H0/hbar,q);
Delta := 0.3534731088e20 * m;

LJ := Delta / ((6*a)^(0.25) * sqrt(eta)) / kpc;

