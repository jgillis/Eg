module Main where

data Dual a = Dual (a,a) deriving (Show , Eq)

instance Num a => Num (Dual a) where
  Dual (x,x') + Dual (y,y') = Dual (x + y, x' + y')
  Dual (x,x') - Dual (y,y') = Dual (x - y, x' - y')
  Dual (x,x') * Dual (y,y') = Dual (x*y, x*y' + x'*y)
  fromInteger k = Dual (fromInteger k, fromInteger 0)

fad f x = map (\(Dual (_,y)) -> y) $ f (Dual (x,1))
  
f x = [x*x, x * x -2*x]

main = do
  print $ f 7
  print $ (fad f 2)
