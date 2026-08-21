restart

-- Setting up K and L

K = QQ; L = K[a]/ideal(a^24-2);

-- Defining a function to compute 
-- the standard coordinates of a polynomial in L

polyToList = p -> ( apply(24, i -> coefficient(a^i,p)) );

-- Setting up the subspace V of L

gensV = {1_L, a - a^15, a^2 - a^12 - a^16, a^4 - a^14, 
    a^6 - a^13 + a^20, a^9 + a^16 + a^19};

gensV = {a, a^2 - a^16, a^3 - a^13 - a^17, a^5 - a^15, 
    a^7 - a^14 + a^21, a^19 + a^17 + a^20};

gensV = {1_L, a - a^3, a^2 - a^6, a^4 - a^9, a^5 - a^12, a^7 - a^15};

V = image transpose matrix apply( gensV, b -> polyToList(b) );
rank V -- rank is 6
 
-- Setting up alpha times V and alpha^2 times V
 
gensaV  = apply( gensV, b-> a*b );
aV = image transpose matrix apply( gensaV, b -> polyToList(b) );
rank aV --  ourput: 6

gensa2V = apply( gensV, b-> a^2*b );
a2V = image transpose matrix apply( gensa2V, b -> polyToList(b) );
rank a2V -- ourput: 6

-- Verifying condition 3.1(1)

rank(V+aV+a2V) == rank V + rank aV + rank a2V -- ourput: true

-- Verifying condition 3.1(2)

gensV2 = flatten apply( 6,i -> apply( i+1, j -> (gensV#i)*(gensV#j) ));
V2 = image transpose matrix apply( gensV2, b -> polyToList(b) );
rank V2 -- ourput: 21

basisL20 = apply( 21, i -> a^i );
L20 = image transpose matrix apply( basisL20, b -> polyToList(b) );
V2 == L20 -- ourput: true

-- Verifying condition 3.1(3)

gensVL21 = flatten apply( gensV, b -> apply( 22, i -> a^i*b ));
VL21 = image transpose matrix apply( gensVL21, b -> polyToList(b) );
rank VL21 -- ourput: 24, so it is all of L

-- Realizing R as quotient of polynomial algebra

-- Setting up the polynomial algebra

Q = QQ[x_1..x_6,y_1..y_6];

-- Setting up the basis for Wt = (V + aV)t 

Lt = L[t];

gensWt = apply( gensV|gensaV, b -> b*t );

-- Define the canonical map phi: Q -> L[t] that maps the 
-- 12 indeterminates to the basis for Wt

phi = map(Lt,Q,gensWt);

I = ker phi;

-- Set up the 78 - 23 = 55 relations in Q to define R

rels = { x_2^2 - y_1^2 - 2*y_1*y_4,
  
  x_2*x_6 + 2*x_1^2 - 36*x_1*x_2 + 72*x_1*x_3 + 42*x_1*x_4 - 2*x_1*x_5
  + 36*x_1*y_1 - 74*x_1*y_2 + 4*x_1*y_3 + 72*x_1*y_4 + 12*x_1*y_5 +
  x_1*y_6 + 2*y_1^2 + 2*y_1*y_2 - 48*y_1*y_3 + 2*y_1*y_4 - 12*y_1*y_6
  - 2*x_2*x_3 - 6*x_2*x_5 - 36*x_2*y_5 + 18*y_2^2,
  
  x_2*y_1 - x_1*y_2, 
  
  x_2*y_2 + 24*x_1*x_2 - 48*x_1*x_3 - 28*x_1*x_4 - 24*x_1*y_1 +
  48*x_1*y_2 - x_1*y_3 - 48*x_1*y_4 - 8*x_1*y_5 - 3*y_1*y_2 +
  32*y_1*y_3 + 8*y_1*y_6 + x_2*x_3 + 4*x_2*x_5 + 24*x_2*y_5 -
  12*y_2^2,
  
  x_2*y_3 + 2*x_1*x_2 - 4*x_1*x_3 - 3*x_1*x_4 - 2*x_1*y_1 + 4*x_1*y_2
  - 4*x_1*y_4 + 2*y_1*y_3 + 2*x_2*y_5 - 2*y_2^2,
  
  x_2*y_4 + x_1*x_5 + 3*x_1*y_2 - x_1*y_3 - 3*y_1^2 + y_1*y_2 -
  4*y_1*y_4,
  
  x_2*y_6 - 6*x_1*x_2 + 12*x_1*x_3 + 10*x_1*x_4 + 8*x_1*y_1 -
  12*x_1*y_2 + 12*x_1*y_4 + 2*x_1*y_5 - 10*y_1*y_3 - 3*y_1*y_6 -
  2*x_2*x_5 - 6*x_2*y_5 + 4*y_2^2,
  
  x_3^2 - 2*x_1^2 + 4*x_1*x_2 - 8*x_1*x_3 - 6*x_1*x_4 - 4*x_1*y_1 +
  8*x_1*y_2 - 8*x_1*y_4 + 4*y_1*y_3 + 4*x_2*y_5 - 3*y_2^2,

  x_3*x_4 + x_1*x_5 + 2*x_1*y_2 - x_1*y_3 - 4*y_1^2 + y_1*y_2 -
  4*y_1*y_4,
  
  x_3*x_5 + 2*x_1*x_2 - 4*x_1*x_3 + x_1*x_4 - 4*x_1*y_1 + 4*x_1*y_2 -
  4*x_1*y_4 - y_1*y_3 + x_2*y_5,
  
  x_3*x_6 - 18*x_1*x_2 + 36*x_1*x_3 + 22*x_1*x_4 + 20*x_1*y_1 -
  36*x_1*y_2 + 36*x_1*y_4 + 5*x_1*y_5 - 24*y_1*y_3 - 5*y_1*y_6 -
  3*x_2*x_5 - 18*x_2*y_5 + 10*y_2^2,
  
  x_3*y_1 - x_1*y_3,
  
  x_3*y_2 + 2*x_1*x_2 - 4*x_1*x_3 - 3*x_1*x_4 - 2*x_1*y_1 + 4*x_1*y_2
  - 4*x_1*y_4 + 2*y_1*y_3 + 2*x_2*y_5 - 2*y_2^2,
  
  x_3*y_3 - 5*x_1*x_2 - 2*x_1*x_6 + 3*x_1*y_1 - 2*x_1*y_2 + 7*x_1*y_4
  + 2*y_1^2 - 4*x_2*x_4,
  
  x_3*y_4 + 22*x_1*x_2 - 44*x_1*x_3 - 27*x_1*x_4 - 22*x_1*y_1 +
  44*x_1*y_2 - 44*x_1*y_4 - 7*x_1*y_5 - 2*y_1*y_2 + 30*y_1*y_3 +
  8*y_1*y_6 + 4*x_2*x_5 + 22*x_2*y_5 - 11*y_2^2,
  
  x_3*y_5 + x_1*x_6 + 2*x_1*y_2 - 2*x_1*y_3 - 4*y_1^2 + 2*y_1*y_2 -
  2*y_1*y_4,
  
  x_3*y_6 + 6*x_1*x_2 - 5*x_1*x_3 + 2*x_1*x_6 - 6*x_1*y_1 + 7*x_1*y_2
  - 11*x_1*y_4 + 3*x_2*x_4 + 2*x_2*y_5,
  
  x_4^2 + 3*x_1*x_2 - 6*x_1*x_3 - 2*x_1*x_4 - 3*x_1*y_1 + 6*x_1*y_2 -
  6*x_1*y_4 + 2*y_1*y_3 + 3*x_2*y_5 - 2*y_2^2,
  
  x_4*x_5 - 2*x_1^2 - 24*x_1*x_2 + 48*x_1*x_3 + 28*x_1*x_4 +
  24*x_1*y_1 - 48*x_1*y_2 + 2*x_1*y_3 + 48*x_1*y_4 + 8*x_1*y_5 +
  x_1*y_6 + 2*y_1*y_2 - 32*y_1*y_3 - 8*y_1*y_6 - 2*x_2*x_3 - 4*x_2*x_5
  - 24*x_2*y_5 + 12*y_2^2,
  
  x_4*x_6 + 4*x_1*x_2 - x_1*x_5 + 2*x_1*x_6 - 4*x_1*y_1 - x_1*y_2 +
  2*x_1*y_3 - 6*x_1*y_4 + y_1^2 - 2*y_1*y_2 + 3*y_1*y_4 + 2*x_2*x_4,
  
  x_4*y_1 - x_1*y_4,
  
  x_4*y_2 + x_1*x_5 + 3*x_1*y_2 - x_1*y_3 - 3*y_1^2 + y_1*y_2 -
  4*y_1*y_4,
  
  x_4*y_3 + 22*x_1*x_2 - 44*x_1*x_3 - 27*x_1*x_4 - 22*x_1*y_1 +
  44*x_1*y_2 - 44*x_1*y_4 - 7*x_1*y_5 - 2*y_1*y_2 + 30*y_1*y_3 +
  8*y_1*y_6 + 4*x_2*x_5 + 22*x_2*y_5 - 11*y_2^2,
  
  x_4*y_4 - 4*x_1*x_2 - x_1*x_6 + 4*x_1*y_1 - x_1*y_2 + 7*x_1*y_4 +
  y_1^2 - 3*x_2*x_4,
  
  x_4*y_5 + 4*x_1*x_2 - 8*x_1*x_3 - 4*x_1*x_4 - 6*x_1*y_1 + 8*x_1*y_2
  - 8*x_1*y_4 + 4*y_1*y_3 + y_1*y_6 + 4*x_2*y_5 - 2*y_2^2,
  
  x_4*y_6 - 38*x_1*x_2 + 76*x_1*x_3 + 47*x_1*x_4 - 2*x_1*x_5 +
  38*x_1*y_1 - 78*x_1*y_2 + 4*x_1*y_3 + 76*x_1*y_4 + 13*x_1*y_5 +
  2*x_1*y_6 + 2*y_1^2 + 2*y_1*y_2 - 52*y_1*y_3 + 2*y_1*y_4 -
  14*y_1*y_6 - 2*x_2*x_3 - 7*x_2*x_5 - 38*x_2*y_5 + 19*y_2^2,
  
  x_5^2 + 4*x_1*x_2 + x_1*x_3 + 4*x_1*x_6 - 4*x_1*y_1 + 5*x_1*y_2 -
  6*x_1*y_4 - 12*y_1^2 + 2*x_2*x_4,
  
  x_5*x_6 + 5*x_1*x_2 - 2*x_1*x_3 - 7*x_1*y_1 + 2*x_1*y_2 - 4*x_1*y_4
  + 2*x_2*y_5,
  
  x_5*y_1 - x_1*y_5,
  
  x_5*y_2 - x_2*y_5, 
  
  x_5*y_3 + x_1*x_6 + 2*x_1*y_2 - 2*x_1*y_3 - 4*y_1^2 + 2*y_1*y_2 -
  2*y_1*y_4,

  x_5*y_4 + 4*x_1*x_2 - 8*x_1*x_3 - 4*x_1*x_4 - 6*x_1*y_1 + 8*x_1*y_2
  - 8*x_1*y_4 + 4*y_1*y_3 + y_1*y_6 + 4*x_2*y_5 - 2*y_2^2,
  
  x_5*y_5 - 144*x_1*x_2 + 288*x_1*x_3 + 168*x_1*x_4 - 2*x_1*x_5 +
  144*x_1*y_1 - 290*x_1*y_2 + 15*x_1*y_3 + 288*x_1*y_4 + 48*x_1*y_5 +
  4*x_1*y_6 + 2*y_1^2 + 15*y_1*y_2 - 192*y_1*y_3 + 2*y_1*y_4 -
  48*y_1*y_6 - 12*x_2*x_3 - 24*x_2*x_5 - 144*x_2*y_5 + 72*y_2^2,
  
  x_5*y_6 + 4*x_1*x_2 + 2*x_1*x_6 - 4*x_1*y_1 + 5*x_1*y_2 + 2*x_1*y_3
  - 6*x_1*y_4 - 7*y_1^2 - 2*y_1*y_2 + 2*x_2*x_4,
  
  x_6^2 + 26*x_1*x_2 - 52*x_1*x_3 - 31*x_1*x_4 - 30*x_1*y_1 +
  52*x_1*y_2 - 52*x_1*y_4 - 4*x_1*y_5 + 33*y_1*y_3 + 4*y_1*y_6 +
  4*x_2*x_5 + 26*x_2*y_5 - 14*y_2^2,
  
  x_6*y_1 - x_1*y_6, 
  
  x_6*y_2 - 6*x_1*x_2 + 12*x_1*x_3 + 10*x_1*x_4 + 8*x_1*y_1 -
  12*x_1*y_2 + 12*x_1*y_4 + 2*x_1*y_5 - 10*y_1*y_3 - 3*y_1*y_6 -
  2*x_2*x_5 - 6*x_2*y_5 + 4*y_2^2,
  
  x_6*y_3 + 6*x_1*x_2 - 5*x_1*x_3 + 2*x_1*x_6 - 6*x_1*y_1 + 7*x_1*y_2
  - 11*x_1*y_4 + 3*x_2*x_4 + 2*x_2*y_5,
  
  x_6*y_4 - 38*x_1*x_2 + 76*x_1*x_3 + 47*x_1*x_4 - 2*x_1*x_5 +
  38*x_1*y_1 - 78*x_1*y_2 + 4*x_1*y_3 + 76*x_1*y_4 + 13*x_1*y_5 +
  2*x_1*y_6 + 2*y_1^2 + 2*y_1*y_2 - 52*y_1*y_3 + 2*y_1*y_4 -
  14*y_1*y_6 - 2*x_2*x_3 - 7*x_2*x_5 - 38*x_2*y_5 + 19*y_2^2,
  
  x_6*y_5 + 4*x_1*x_2 + 2*x_1*x_6 - 4*x_1*y_1 + 5*x_1*y_2 + 2*x_1*y_3
  - 6*x_1*y_4 - 7*y_1^2 - 2*y_1*y_2 + 2*x_2*x_4,
  
  x_6*y_6 + 4*x_1*x_2 + 4*x_1*x_3 - 2*x_1*x_6 - 4*x_1*y_1 - 6*x_1*y_2
  - x_1*y_4 - 2*y_1^2 - x_2*x_4,
  
  y_2*y_3 - 6*x_1*x_2 - 2*x_1*x_6 + 6*x_1*y_1 - 2*x_1*y_2 + 9*x_1*y_4
  + 2*y_1^2 - 4*x_2*x_4,

  y_2*y_4 + 34*x_1*x_2 - 68*x_1*x_3 - 41*x_1*x_4 - 34*x_1*y_1 +
  68*x_1*y_2 - x_1*y_3 - 68*x_1*y_4 - 11*x_1*y_5 - 2*y_1*y_2 +
  46*y_1*y_3 + 12*y_1*y_6 + x_2*x_3 + 6*x_2*x_5 + 34*x_2*y_5 -
  17*y_2^2,
  
  y_2*y_5 - 2*x_1*x_2 - x_1*x_6 + 2*x_1*y_1 - 2*x_1*y_3 + 3*x_1*y_4 +
  2*y_1*y_2 - 2*y_1*y_4 - x_2*x_4,
  
  y_2*y_6 + 8*x_1*x_2 - 3*x_1*x_3 + 2*x_1*x_6 - 8*x_1*y_1 + 5*x_1*y_2
  - 11*x_1*y_4 - y_1*y_5 + 3*x_2*x_4 + x_2*y_5,
  
  y_3^2 + 24*x_1*x_2 - 48*x_1*x_3 - 28*x_1*x_4 + 4*x_1*x_5 -
  24*x_1*y_1 + 55*x_1*y_2 - 6*x_1*y_3 - 48*x_1*y_4 - 8*x_1*y_5 -
  2*x_1*y_6 - 9*y_1^2 + 32*y_1*y_3 - 9*y_1*y_4 + 8*y_1*y_6 + 2*x_2*x_3
  + 4*x_2*x_5 + 24*x_2*y_5 - 12*y_2^2,
  
  y_3*y_4 - x_1*x_2 - 2*x_1*x_4 + x_1*y_1 + 2*y_1*y_3 + y_1*y_5 -
  2*y_2^2,
  
  y_3*y_5 - 16*x_1*x_2 + 32*x_1*x_3 + 16*x_1*x_4 + 16*x_1*y_1 -
  32*x_1*y_2 + 2*x_1*y_3 + 32*x_1*y_4 + 4*x_1*y_5 + x_1*y_6 +
  2*y_1*y_2 - 20*y_1*y_3 - 4*y_1*y_6 - 2*x_2*x_3 - 2*x_2*x_5 -
  16*x_2*y_5 + 8*y_2^2,
  
  y_3*y_6 + 4*x_1*x_2 - 3*x_1*x_5 + 2*x_1*x_6 - 4*x_1*y_1 - 3*x_1*y_2
  + 2*x_1*y_3 - 6*x_1*y_4 + 2*x_1*y_6 + 3*y_1^2 + 5*y_1*y_4 +
  2*x_2*x_4,
  
  y_4^2 + 12*x_1*x_2 - 24*x_1*x_3 - 14*x_1*x_4 + 3*x_1*x_5 -
  12*x_1*y_1 + 29*x_1*y_2 - 4*x_1*y_3 - 24*x_1*y_4 - 4*x_1*y_5 -
  x_1*y_6 - 5*y_1^2 + y_1*y_2 + 16*y_1*y_3 - 5*y_1*y_4 + 4*y_1*y_6 +
  x_2*x_3 + 2*x_2*x_5 + 12*x_2*y_5 - 6*y_2^2,
  
  y_4*y_5 - 2*x_1*x_2 + x_1*x_3 + 2*x_1*y_1 - x_1*y_2 + 3*x_1*y_4 -
  2*y_1^2 + y_1*y_5 - x_2*x_4 - x_2*y_5,
  
  y_4*y_6 - 3*x_1*x_2 + 10*x_1*x_3 + 10*x_1*x_4 + 3*x_1*y_1 -
  10*x_1*y_2 + 10*x_1*y_4 + 2*x_1*y_5 - 10*y_1*y_3 - y_1*y_5 -
  2*y_1*y_6 - 2*x_2*x_5 - 5*x_2*y_5 + 4*y_2^2,

  y_5^2 + 6*x_1*x_2 - 12*x_1*x_3 - 7*x_1*x_4 - 6*x_1*y_1 + 12*x_1*y_2
  - 12*x_1*y_4 + 2*x_1*y_5 + 8*y_1*y_3 - 2*x_2*x_5 + 6*x_2*y_5 -
  3*y_2^2,
  
  y_5*y_6 - 80*x_1*x_2 + 160*x_1*x_3 + 96*x_1*x_4 - 2*x_1*x_5 +
  80*x_1*y_1 - 162*x_1*y_2 + 9*x_1*y_3 + 160*x_1*y_4 + 28*x_1*y_5 +
  2*x_1*y_6 + 2*y_1^2 + 10*y_1*y_2 - 108*y_1*y_3 + 2*y_1*y_4 -
  28*y_1*y_6 - 7*x_2*x_3 - 14*x_2*x_5 - 80*x_2*y_5 + 40*y_2^2,
  
  y_6^2 - 24*x_1*x_2 + 48*x_1*x_3 + 28*x_1*x_4 + x_1*x_5 + 24*x_1*y_1
  - 41*x_1*y_2 + 5*x_1*y_3 + 48*x_1*y_4 + 8*x_1*y_5 - 2*x_1*y_6 -
  7*y_1^2 - 3*y_1*y_2 - 32*y_1*y_3 - 5*y_1*y_4 - 8*y_1*y_6 - 2*x_2*x_3
  - 4*x_2*x_5 - 24*x_2*y_5 + 12*y_2^2 };

-- Verify the number of relations

#rels -- ourput: 55

-- Verify that the relations generate the kernel of phi

ideal rels == ker phi -- output: true

R = Q/(ideal rels)

-- Verify that the Hilbert function of R is as expected

apply(6,d -> hilbertFunction(d,R)) -- output: {1, 23, 23, 24, 24, 24}

-- Set up the ideals in R corresponding to t and at
 
X = ideal(x_1); 
Y = ideal(y_1);

-- Verify the condition on these ideals from Huneke-Iyengar-Wiegand

(X:Y)*(Y:X)  == intersect(X:Y,Y:X) -- output: true

-- Notice the relation between the colon ideals observed in proof of 3.2

mingens (X:Y) -- output: | x_6 x_5 x_4 x_3 x_2 x_1 |
mingens (Y:X) -- output: | y_6 y_5 y_4 y_3 y_2 y_1 |

