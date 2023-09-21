import Numeric.LinearAlgebra
import Chart

-- Parameters

-- Δu(x, y) = f(x, y)
f x y = -sin(x) - sin(y)

-- u(x,0)
p x = sin x
-- u(0,y)
q y = sin y
-- u(x,2π)
r x = sin x
-- u(2π,y)
s y = sin y

xDomain = (0.0, 2*pi)
yDomain = (0.0, 2*pi)

-- Define finite difference matrices Au = v

a :: Int -> Matrix Double
a n = fromLists [[ if (i == j) then -4 else if (i == j-1 || j == i-1 || i == j+n || j == i+n) then 1 else 0 | i <- [1 .. n^2]] | j <- [1 .. n^2]]

fdm :: (Double -> Double -> Double) -- f : Δu = f
    -> (Double, Double)             -- (x0, xt) : x domain
    -> (Double, Double)             -- (y0, yt) : y domain
    -> (Double -> Double)           -- p : p(x) = u(x,0)
    -> (Double -> Double)           -- q : q(y) = u(0,y)
    -> (Double -> Double)           -- r : r(x) = u(x,yt)
    -> (Double -> Double)           -- s : s(y) = u(xt,y)
    -> Int                          -- n : number of steps
    -> Matrix Double
fdm (f) (x0, xt) (y0, yt) p q r s numSteps =
    let h = (xt - x0) / fromIntegral numSteps -- for now, calculate h based on x range
        c = 1 / h^2

        v :: (Double, Double) -> Double -> Double
        v (i, j) n | (i, j) == (1, 1)   = f i j - c * (p i + q j) 
                   | (i, j) == (1, j)   = f i j - c * q j
                   | (i, j) == (1, n)   = f i j - c * (q j + r i) 
                   | (i, j) == (i, 1)   = f i j - c * p i
                   | (i, j) == (i, j)   = f i j
                   | (i, j) == (i, n)   = f i j - c * r i
                   | (i, j) == (n, 1)   = f i j - c * (p i + s j)
                   | (i, j) == (n, j)   = f i j - c * s j
                   | (i, j) == (n, n)   = f i j - c * (r i + s j)

        -- define mesh
        xs = [x0, x0 + h .. xt]; ys = [y0, y0 + h .. yt]
        mesh = [[(x, y) | x <- xs] | y <- ys]

    in col [ v (x, y) (fromIntegral numSteps) | row <- mesh, (x, y) <- row ]

u numSteps = luSolve (luPacked $ a (numSteps+1)) (fdm (f) xDomain yDomain p q r s numSteps)

save numSteps = do 
    let us = ((numSteps+1)><(numSteps+1)) (concat $ toLists $ u numSteps)
        plotU = contour $ toLists $ us
    writeFile "data/eulerData.txt" (show us) 
    plotU
