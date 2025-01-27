import Data.Complex
import Data.List (sortBy)
import Text.Printf
import Data.Time.Clock
import Debug.Trace (trace)
-- Calculates derivative coefficients of a polynomial
-- Input: List of coefficients [a_n, a_{n-1}, ..., a_1, a_0]
-- Output: List of derivative coefficients [na_n, (n-1)a_{n-1}, ..., a_1]
derivative :: Num a => [a] -> [a]
derivative [] = []
derivative [_] = []
derivative (p:ps) = (fromIntegral (length ps) * p) : derivative ps

-- Evaluates polynomial at complex point x
-- Uses Horner's method for numerical stability
-- iterate generates: [1, x, x², x³, ...]
-- reverpse p [1,2,3] -> reverse p [3,2,1]
-- zipWith applay function on two arrays
-- map (:+0) p convert to complex numbers
polyval :: [Double] -> Complex Double -> Complex Double
polyval p x = sum $ zipWith (*) (map (:+ 0) (reverse p)) (iterate (*x) (1 :+ 0))

-- Divides polynomial p by polynomial q at point x
-- Uses different methods based on magnitude of x to avoid overflow
divide :: [Double] -> [Double] -> Complex Double -> Complex Double
divide p q x
 | magnitude x <= 1 = px / qx       -- Standard division for small x
 | otherwise = x * (pr / qr)        -- Modified for large x to avoid overflow
 where
   px = polyval p x
   qx = polyval q x
   pr = polyval (reverse p) (1/x)
   qr = polyval (reverse q) (1/x)


eulerEquation :: Double -> Double -> Complex Double
eulerEquation r theta =
  let realPart = r * cos theta
      imagPart = r * sin theta
  in realPart :+ imagPart

initRoots :: [Double] -> [Complex Double]
initRoots p = [ eulerEquation r  (2 * pi * fromIntegral k / fromIntegral n) | k <- [0..n-2] ]
  where
    r = (1 + maximum (map abs (drop 1 p)) / (abs (head p))) / 2
    n = length p

-- Calculates correction terms for Aberth-Ehrlich iteration
getOffsets :: [Double] -> [Double] -> [Complex Double] -> [Complex Double]
getOffsets p pTag roots = zipWith (/) nums denoms
 where
   nums = map (divide p pTag) roots                                                            -- Newton correction terms
   sums = [sum [1/(z1 - z2) | (z2, j) <- zip roots [1..], i /= j] | (z1, i) <- zip roots [1..]]-- Deflation terms
   denoms = zipWith (\n s -> 1 - n * s) nums sums


-- Main iteration function for Aberth-Ehrlich method
abertErlich :: [Double] -> [Double] -> [Complex Double] -> Double -> Int -> [Complex Double]
abertErlich p pTag roots epsilon maxTries = go roots 0
 where
   go currentRoots tries
     | maxW < epsilon || tries >= maxTries = newRoots           -- Convergence check
     | otherwise = go newRoots (tries + 1)                      -- Continue iteration
     where
       w = getOffsets p pTag currentRoots
       newRoots = zipWith (-) currentRoots w                    -- Update approximations
       maxW = maximum $ map magnitude w
       
       
readFileToList :: FilePath -> IO [Double]
readFileToList filePath = do
    contents <- readFile filePath
    let numbers = map read (lines contents) :: [Double]
    return numbers

sortByRealThenImag :: [Complex Double] -> [Complex Double]
sortByRealThenImag = sortBy compareComplex
  where
    compareComplex a b = case compare (realPart a) (realPart b) of
        EQ -> compare (imagPart a) (imagPart b)
        other -> other


main :: IO ()
main = do
    p <- readFileToList "poly_coeff_alberth.txt"
    start <- getCurrentTime
    let dp = derivative p
    let roots = initRoots p
    let abertErlichRoots = abertErlich p dp roots 1e-6 600
    mapM_ print (sortByRealThenImag abertErlichRoots)
    putStrLn "Found roots:"
    end <- getCurrentTime
    let duration = diffUTCTime end start
    putStrLn( "The operation took: " ++ show duration ++ " seconds" )



 
