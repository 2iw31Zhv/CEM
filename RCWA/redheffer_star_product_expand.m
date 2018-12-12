function S = redheffer_star_product_expand(SA11, SA12, SA21, SA22,...
    SB11, SB12, SB21, SB22)

[a, ~] = size(SA11);
[b, ~] = size(SA22);
% [b, ~] = size(SB11);
[c, ~] = size(SB22);

S = zeros(a+c, a+c);
I = eye(b);

S(1:a, 1:a) = SA11 + SA12 * (I - SB11 * SA22)^(-1) * SB11 * SA21;
S(1:a, (a+1):(a+c)) = SA12 * (I - SB11 * SA22)^(-1) * SB12;
S((a+1):(a+c), 1:a) = SB21 * (I - SA22 * SB11)^(-1) * SA21;
S((a+1):(a+c), (a+1):(a+c)) = SB22 + SB21 * (I - SA22 * SB11)^(-1) * SA22 * SB12;

end

