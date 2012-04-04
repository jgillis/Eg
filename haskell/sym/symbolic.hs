data Expr = Val Double
          | Sym String
          | Mul Expr Expr
          | Add Expr Expr
          | Sub Expr Expr
          | Sin Expr
          | Cos Expr
          deriving Eq

paren x = "( " ++ x ++ " )"

instance Show Expr where
  show (Val x) = show x
  show (Sym n) = n
  show (Mul x y) = paren $ show x ++ "*" ++ show y
  show (Add x y) = paren $ show x ++ "+" ++ show y
  show (Sub x y) = paren $ show x ++ "-" ++ show y
  show (Sin x) = "sin( " ++ show x ++ " )"
  show (Cos x) = "cos( " ++ show x ++ " )"

instance Num Expr where
  (+) = Add
  (-) = Sub
  (*) = Mul
  fromInteger = Val . fromInteger

instance Fractional Expr where
  fromRational = Val . fromRational
  
instance Floating Expr where
  sin = Sin
  cos = Cos

dot :: Expr -> Expr
dot (Val x) = 0
dot (Sym name) = Sym $ "d"++name++"/dt"
dot (Mul x y) = x*(dot y) + (dot x)*y
dot (Add x y) = dot x + dot y
dot (Sub x y) = dot x - dot y
--dot (Div x y) = dot 
dot (Sin x) = cos x * dot x
dot (Cos x) = -sin x * dot x


main = do
  let x = Sym "x"
      y = Sym "y"
      
      f = x*y*cos(x+2*y)
  print f
  print $ dot f
