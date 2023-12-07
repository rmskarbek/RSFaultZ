function x = bclog(z, theta)
%%% Zp locates the branch cut. Zp2 shifts the displacements.


% arg[z_, σ_: - Pi] := Arg[z Exp[-I (σ + Pi)]] + σ + Pi;
% log[z_, σ_: - Pi] := Log[Abs[z]] + I arg[z, σ]

    x = log(abs(z)) + 1i*angle(z.*exp(-1i*(theta + pi))) + 1i*pi;

%    x = log(abs(z)) + 1i*angle(z.*exp(-1i*(theta + pi))) + theta + pi;% + 1i*pi;
%    x = log(abs(z)) + 1i*angle(z.*exp(1i*(theta+pi))) - 1i*theta + 1i*pi;
%    x = log(abs(z)) + 1i*(angle(z) + pi);
end