xi = 1.5          -- X inicial
err = 0.00001     -- Erro predefinido
f x = x^2 + x - 6 -- Função F(x)=0

-- Função iterativa do método de newton
fn x = x - (x^2 + x - 6)/(2*x + 1)
-- Função falsa posição
ffp a b = (a*(f b) - b*(f a))/((f b) - (f a))
-- Função iterativa do ponto fixo
fix x = sqrt(6 - x)
-- Funcao para metodo da secante
fs k k1 = k1 - (f k1)*(k1 - k)/((f k1) - (f k))

-- Critérios de parada
relative e xk xk1 = (abs xk1)*e > (abs $ xk1 - xk)
absolute e xk xk1 = e > (abs $ xk1 - xk)
residual e f xk = e > abs (f xk)

-- Chamadas para aplicação dos métodos:
newton e x0 = newton' e x0 (fn x0)
newton' e xk xk1
  | (residual e (f) xk) = []
  | otherwise = xk:newton' e xk1 xk2
  where xk2 = fn xk1

bissect e [a,b]
  | (absolute e a b) = []
  | otherwise = cut:bissect e i
  where i = if (f cut)*(f a)<0 then [a,cut] else [cut,b]
        cut = (a+b)/2

falsePosition e [a,b]
  | (absolute e a b) = []
  | otherwise = cut:bissect e i
  where i = if (f cut)*(f a)<0 then [a,cut] else [cut,b]
        cut = ffp a b

fixedPoint e x0 = fixedPoint' e x0 $ fix x0
fixedPoint' e xk xk1
  | (absolute e xk xk1) = []
  | otherwise = xk2:fixedPoint' e xk1 xk2
  where xk2 = fix xk1

sec e xk xk1
  | (absolute e xk xk1) = []
  | otherwise = xk1:sec e xk1 xk2
  where xk2 = fs xk xk1
