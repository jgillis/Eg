module Main where
c = putStrLn "C!"
 
combine before after =
  do before
     putStrLn "In the middle"
     after
 
main = do combine c c
          let b = combine (putStrLn "Hello!") (putStrLn "Bye!")
          let d = combine (b) (combine c c)
	  f <- d
          putStrLn "So long!"
