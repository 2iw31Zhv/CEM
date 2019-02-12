function S = redheffer_star_product(SA, SB)

[b, ~] = size(SA.S22);

S = struct;
I = eye(b);

S.S11 = SA.S11 + SA.S12 / (I - SB.S11 * SA.S22) * SB.S11 * SA.S21;
S.S12 = SA.S12 / (I - SB.S11 * SA.S22) * SB.S12;
S.S21 = SB.S21 / (I - SA.S22 * SB.S11) * SA.S21;
S.S22 = SB.S22 + SB.S21 / (I - SA.S22 * SB.S11) * SA.S22 * SB.S12;

end
