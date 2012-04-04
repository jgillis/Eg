cartprod :: (Num a) => [[a]] -> [a]

cartprod [x] = x

cartprod [x:xs] = [i*j | i<-[x] , j<-cartprod [xs]]

base = [3,5]
