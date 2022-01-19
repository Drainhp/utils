local v3 = Vector3
local cf = CFrame
local unpack = unpack
local typeof = typeof

local sqrt = math.sqrt
local cos = math.cos
local sin = math.sin
local tan2 = math.atan2

-- class

local matrix3 = {}
local mt = {}
local ref = {}

local function sumType(t, x)
	local sum = 0
	for i, v in pairs(t) do
		sum += typeof(v) == x and 1 or 0
	end
	return sum
end

-- metamethods

function mt.__index(cf, k)
	return matrix3[k]
end

-- public constructors

function matrix3.new(...)
	local self = {}
	local components = {...}

	if #components < 0 or #components > 9 then
		error("Invalid number of arguments: " .. #components)
	elseif #components == 0 then
		components = {1, 0, 0, 0, 1, 0, 0, 0, 1}
	elseif #components == 4 and sumType(components, "number") == 4 then
		-- from quaternion
		local x, y, z, w = unpack(components)
		
		local m = 1 / sqrt(x*x + y*y + z*z + w*w)
		x, y, z, w = x*m, y*m, z*m, w*m
		components = {
			1-2*(y*y+z*z), 2*(y*x-w*z), 2*(w*y+z*x),
			2*(w*z+x*y), 1-2*(x*x+z*z), 2*(z*y-w*x),
			2*(x*z-w*y), 2*(w*x+z*y), 1-2*(y*y+z*z)
		}
	elseif #components == 9 and sumType(components, "number") == 9 then
		-- from components
	else
		error("Invalid arguments")
	end

	ref[self] = components
	return setmetatable(self, mt)
end

function matrix3.lookAt(eye, target, upDir)
	upDir = upDir.Unit or v3.yAxis
	local forward = (target - eye).Unit
	local right = forward:Cross(upDir).Unit
	local up = right:Cross(forward).Unit

	return matrix3.new(
		right.X, up.X, forward.X,
		right.Y, up.Y, forward.Y,
		right.Z, up.Z, forward.Z
	)
end

function matrix3.fromEulerAnglesXYZ(rx, ry, rz)
	local cx, sx = cos(rx), sin(rx)
	local cy, sy = cos(ry), sin(ry)
	local cz, sz = cos(rz), sin(rz)
	
	return matrix3.new(
		cy*cz, -cy*sz, sy,
		sx*sy*cz+cx*sz, -sx*sy*sz+cx*cz, -sx*cy,
		-cx*sy*cz+sx*sz, cx*sy*sz+sx*cz, cx*cy
	)
end

function matrix3.fromEulerAnglesYXZ(rx, ry, rz)
	local cx, sx = cos(rx), sin(rx)
	local cy, sy = cos(ry), sin(ry)
	local cz, sz = cos(rz), sin(rz)
	
	return matrix3.new(
		cy*cz+sy*sx*sz, -cy*sz+sy*sx*cz, sy*cx,
		cx*sz, cx*cz, -sx,
		-sy*cz+cy*sx*sz, sy*sz+cy*sx*cz, cy*cx
	)
end

function matrix3.fromAxisAngle(v, r)
	v = v.Unit
	local vx, vy, vz = v.X, v.Y, v.Z
	local c, s = cos(r), sin(r)

	return matrix3.new(
		c + (1-c) * vx*vx, (1-c)*vx*vy - s*vz, (1-c)*vx*vz + s*vy,
		(1-c)*vy*vx + s*vz, c + (1-c)*vy*vy, (1-c)*vy*vz - s*vx,
		(1-c)*vz*vx - s*vy, (1-c)*vz*vy + s*vx, c + (1-c)*vz*vz
	)
end

function matrix3.fromMatrix(vx, vy, vz)
	if not vz then
		vz = vx:Cross(vy).Unit
	end
	return matrix3.new(vx.X, vx.Y, vx.Z, vy.X, vy.Y, vy.Z, vz.X, vz.Y, vz.Z)
end

-- public methods

function matrix3:GetComponents()
	return unpack(ref[self])
end

function matrix3:Determinant()
	local m = ref[self]

	return m[1] * (m[5]*m[9] - m[8]*m[6]) - m[2] * (m[4]*m[9] - m[7]*m[6]) + m[3] * (m[4]*m[8] - m[7]*m[5])
end

matrix3.components = matrix3.GetComponents

return matrix3