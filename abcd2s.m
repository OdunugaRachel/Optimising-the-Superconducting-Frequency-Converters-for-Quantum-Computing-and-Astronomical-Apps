function s = abcd2s( abcd, Z0) 

a = abcd( 1, 1);
b = abcd( 1, 2);
c = abcd( 2, 1);
d = abcd( 2, 2);

denom = (a + b/Z0 + c*Z0 + d);

s11 = (a + b/Z0 - c*Z0 - d) / denom;
s12 = 2 * (a*d - b*c) / denom;
s21 = 2 / denom;
s22 = (-a + b/Z0 - c*Z0 + d) / denom;

s = [s11 s12; s21 s22];
