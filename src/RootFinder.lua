-- Root finding algorithm to find the roots of any degree polynomial.
-- Put this code in a module script and require it. The module will return a 
-- function that takes an array of coefficients of a polynomial as an argument 
-- (ordered with constant coefficient first and leading coefficient last, e.g.
-- for a polynomial ax^3+bx+2+cx+d, you would pass:
-- {d, c, b, a}
-- as the array of coefficients).

-- The variable below controls the number of iterations for the bisection
-- method. More iterations means more precision, but slower performance
local ITERATIONS = 100

-- The minimum and maximum bounds in which solutions are allowed to exist.
local SOLUTION_BOUNDS = 10

-- Error margin within which a root is considered to be found
local EPSILON = 1e-3

-- Get value of polynomial
local function getPolynomialValue(coefficients, x)
	local y = 0
	for i = 1, #coefficients do
		y += coefficients[i] * x^(i-1)
	end
	return y
end

-- Compute and return the derivative of polynomial
local function getPolynomialDerivative(coefficients)
	local derivative = {}
	for i = 1, #coefficients - 1 do
		derivative[i] = coefficients[i+1] * i
	end
	return derivative
end

-- Compute polynomial roots from a given array of coefficients
local function FindPolynomialRoots(coefficients)
	if #coefficients == 2 then
		-- Polynomial is a line. Solve it
		local a = coefficients[2]
		local b = coefficients[1]
		local x = -b / a

		return {x}
	end
	
	-- Compute the derivative, and find its roots
	local derivative = getPolynomialDerivative(coefficients)
	local droots = FindPolynomialRoots(derivative)
	
	local roots = {}
	
	for i = 1, #droots + 1 do
		local xMin = -SOLUTION_BOUNDS
		local xMax = SOLUTION_BOUNDS
		if i > 1 then
			xMin = droots[i - 1]
		end
		if i <= #droots then
			xMax = droots[i]
		end
		
		local x = 0
		local y = 0
		local sign = (#coefficients - i) % 2 == 0 and -1 or 1
		
		-- Bisection method
		for n = 1, ITERATIONS do
			x = (xMin + xMax) * 0.5
			y = getPolynomialValue(coefficients, x) * sign
			if y < 0 then
				xMin = x
			else
				xMax = x
			end
		end
		
		if math.abs(y) < EPSILON then
			roots[#roots + 1] = x
		end
	end
	
	return roots
end

return FindPolynomialRoots
