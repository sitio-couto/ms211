fn x = x - (x^5 - (10/9)*x^3 + (5/21)*x)/(5*(x^4) - 10*(x^2)/3 + 5/21)
f x = (x-(sqrt $ 64/(900/x^2-1)))^2 +64-(400/x^2)*((x-(sqrt $ 64/(900/x^2-1)))^2)
ffp a b = (a*(f b) - b*(f a))/((f b) - (f a))
fix x = (x^3)/3
fs k k1 = k1 - (f k1)*(k1 - k)/((f k1) - (f k))

newton e x0 = newton' e x0 $ (abs $ x0)+e+1
newton' e xk1 xk
  | ((abs $ xk - xk1) < e) = []
  | otherwise = xk2:newton' e xk2 xk1
  where xk2 = fn xk1

bissect e [a,b]
  | ((abs $ a-b) < e) = []
  | otherwise = cut:bissect e i
  where i = if (f cut)*(f a)<0 then [a,cut] else [cut,b]
        cut = (a+b)/2

falsePosition e [a,b]
  | ((abs $ a-b) < e) = []
  | otherwise = cut:bissect e i
  where i = if (f cut)*(f a)<0 then [a,cut] else [cut,b]
        cut = ffp a b

fixatedPoint e x0 = fixatedPoint' e x0 $ (abs $ x0)+e+1
fixatedPoint' e xk1 xk
  | ((abs $ xk - xk1) < e) = []
  | otherwise = xk2:fixatedPoint' e xk2 xk1
  where xk2 = fix xk1

sec e k k1
  | ((abs $ k - k1) < e) = []
  | otherwise = k1:sec e k1 k2
  where k2 = fs k k1
