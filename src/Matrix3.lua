local v3 = Vector3
local cf = CFrame
local unpack = unpack
local typeof = typeof

local sqrt = math.sqrt
local cos = math.cos
local sin = math.sin
local acos = math.acos
local asin = math.asin
local atan2 = math.atan2

local format = string.format

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

local function matrix3ToQuaternion(m)
	local m11, m12, m13, m21, m22, m23, m31, m32, m33 = m:GetComponents()
	local qw, qx, qy, qz
	
	if (m11 + m22 + m33 > 0) then
		qw = sqrt(1 + m11 + m22 + m33) * 0.5
		qx = (m32-m23) / (4*qw)
		qy = (m13-m31) / (4*qw)
		qz = (m21-m12) / (4*qw)
	elseif (m11 > m22 and m11 > m33) then
		qx = sqrt(1 + m11 - m22 - m33) * 0.5
		qy = (m21+m12) / (4*qx)
		qz = (m31+m13) / (4*qx)
		qw = (m32-m23) / (4*qx)
	elseif (m22 > m33) then
		qy = sqrt(1 + m22 - m11 - m33) * 0.5
		qx = (m21+m12) / (4*qy)
		qz = (m32+m23) / (4*qy)
		qw = (m13-m31) / (4*qy)
	else
		qz = sqrt(1 + m33 - m11 - m22) * 0.5
		qx = (m31+m13) / (4*qz)
		qy = (m32+m23) / (4*qz)
		qw = (m21-m12) / (4*qz)
	end
	
	return qw, qx, qy, qz
end

-- metamethods

function mt.__index(matrix, property)
	local m = ref[matrix]

	if property == "XVector" then
		return v3.new(m[1], m[4], m[7])
	elseif property == "YVector" then
		return v3.new(m[2], m[5], m[8])
	elseif property == "ZVector" then
		return v3.new(m[3], m[6], m[9])
	elseif matrix3[property] then
		return matrix3[property]
	else
		error(property .. " is not a valid member of Matrix3")
	end
end

function mt.__newindex(matrix, property, value)
	error(property .. " cannot be assigned to")
end

function mt.__eq(a, b)
	if ref[a] and ref[b] then
		a, b = ref[a], ref[b]
		for i = 1, 9 do
			if a[i] ~= b[i] then return false end
		end
		return true
	else
		return false
	end
end

function mt.__mul(a, b)
	if ref[a] and ref[b] then
		-- Matrix3 * Matrix3
		a, b = ref[a], ref[b]

		return matrix3.new(
			a[1]*b[1] + a[2]*b[4] + a[3]*b[7], a[1]*b[2] + a[2]*b[5] + a[3]*b[8], a[1]*b[3] + a[2]*b[6] + a[3]*b[9], 
			a[4]*b[1] + a[5]*b[4] + a[6]*b[7], a[4]*b[2] + a[5]*b[5] + a[6]*b[8], a[4]*b[3] + a[5]*b[6] + a[6]*b[9], 
			a[7]*b[1] + a[8]*b[4] + a[9]*b[7], a[7]*b[2] + a[8]*b[5] + a[9]*b[8], a[7]*b[3] + a[8]*b[6] + a[9]*b[9]
		)
	elseif ref[a] and typeof(b) == "CFrame" then
		-- Matrix3 * CFrame
		a = ref[a]
		local x, y, z, c11, c12, c13, c21, c22, c23, c31, c32, c33 = b:GetComponents()

		return matrix3.new(
			a[1]*c11 + a[2]*c21 + a[3]*c31, a[1]*c12 + a[2]*c22 + a[3]*c32, a[1]*c13 + a[2]*c23 + a[3]*c33, 
			a[4]*c11 + a[5]*c21 + a[6]*c31, a[4]*c12 + a[5]*c22 + a[6]*c32, a[4]*c13 + a[5]*c23 + a[6]*c33, 
			a[7]*c11 + a[8]*c21 + a[9]*c31, a[7]*c12 + a[8]*c22 + a[9]*c32, a[7]*c13 + a[8]*c23 + a[9]*c33
		)
	elseif ref[a] and typeof(b) == "Vector3" then
		-- Matrix3 * Vector3
		a = ref[a]
		local x, y, z = b.X, b.Y, b.Z

		return v3.new(
			a[1]*x + a[2]*y + a[3]*z,
			a[4]*x + a[5]*y + a[6]*z,
			a[7]*x + a[8]*y + a[9]*z
		)
	elseif ref[a] and typeof(b) == "number" then
		-- Matrix3 * number
		a = ref[a]

		return matrix3.new(
			a[1]*b + a[2]*b + a[3]*b,
			a[4]*b + a[5]*b + a[6]*b,
			a[7]*b + a[8]*b + a[9]*b
		)
	elseif typeof(a) == "number" and ref[b] then
		-- number * Matrix3
		b = ref[b]

		return matrix3.new(
			a*b[1] + a*b[2] + a*b[3],
			a*b[4] + a*b[5] + a*b[6],
			a*b[7] + a*b[8] + a*b[9]
		)
	elseif ref[a] then
		error("bad argument #2 to '?' (Vector3 expected, got " .. typeof(b) .. ")")
	else
		error("bad argument #1 to '?' (Matrix3 expected, got " .. typeof(a) .. ")")
	end
end

function mt.__tostring(matrix)
	return format("%s, %s, %s, %s, %s, %s, %s, %s, %s", unpack(ref[matrix]))
end

mt.__metatable = false

-- public constructors

function matrix3.new(...)
	local self = {}
	local components = {...}

	if #components < 0 or #components > 9 then
		error("Invalid number of arguments: " .. #components)
	elseif #components == 0 then
		components = {1, 0, 0, 0, 1, 0, 0, 0, 1}
	elseif #components == 1 and typeof(components[1] == "CFrame") then
		-- from cframe
		x, y, z, m11, m12, m13, m21, m22, m23, m31, m32, m33 = components[1]:GetComponents()
		components = {m11, m12, m13, m21, m22, m23, m31, m32, m33}
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

matrix3.Angles = matrix3.fromEulerAnglesXYZ;
matrix3.fromOrientation = matrix3.fromEulerAnglesYXZ;

-- public methods

function matrix3:GetComponents()
	return unpack(ref[self])
end

function matrix3:Determinant()
	local m = ref[self]

	return m[1] * (m[5]*m[9] - m[8]*m[6]) - m[2] * (m[4]*m[9] - m[7]*m[6]) + m[3] * (m[4]*m[8] - m[7]*m[5])
end

function matrix3:Transpose()
	local m = ref[self]

	return matrix3.new(
		m[1], m[4], m[7],
		m[2], m[5], m[8],
		m[3], m[6], m[9]
	)
end

function matrix3:Inverse()
	local m = ref[self]

	local invDet = 1 / self:Determinant()

	return matrix3.new(
		(m[5]*m[9] - m[7]*m[6]) * invDet,
		(m[3]*m[7] - m[2]*m[9]) * invDet,
		(m[2]*m[6] - m[3]*m[5]) * invDet,
		(m[6]*m[7] - m[4]*m[9]) * invDet,
		(m[1]*m[9] - m[3]*m[7]) * invDet,
		(m[4]*m[3] - m[1]*m[6]) * invDet,
		(m[4]*m[7] - m[7]*m[5]) * invDet,
		(m[7]*m[2] - m[1]*m[7]) * invDet,
		(m[1]*m[5] - m[4]*m[2]) * invDet
	)
end

function matrix3:Transform(m2)
	return self * m2
end

function matrix3:InverseTransform(m2)
	return self:Inverse() * m2
end

function matrix3:VectorTransform(v)
	return self * v
end

function matrix3:VectorInverseTransform(v)
	return self:Inverse() * v
end

function matrix3:ToEulerAnglesXYZ()
	local m = ref[self]
	local rx = atan2(-m[9], m[12])
	local ry = asin(m[6])
	local rz = atan2(-m[5], m[4])
	return rx, ry, rz
end

function matrix3:ToEulerAnglesYXZ()
	local m = ref[self]
	local rx = asin(-m[9])
	local ry = atan2(m[6], m[12])
	local rz = atan2(m[7], m[8])
	return rx, ry, rz
end

function matrix3:ToAxisAngle()
	local qw, qx, qy, qz = matrix3ToQuaternion(self);
	
	-- pick the twin closest to identity quaternion
	if (qw <= 0) then
		qw, qx, qy, qz = -qw, -qx, -qy, -qz
	end	
	
	local theta = acos(qw) * 2
	local axis = v3(qx, qy, qz) / sin(theta*0.5)
	
	if (axis:Dot(axis) > 0) then
		return axis.unit, theta
	else
		return v3(1, 0, 0), theta
	end
end

function matrix3:Lerp(m2, t)
	local diff = self:inverse() * m2
	local axis, theta = diff:toAxisAngle()
	local m = ref[self * matrix3.fromAxisAngle(axis, theta * t)]
	return cframe.new(m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9])
end

function matrix3:ToCFrame()
	local m = ref[self]

	return cf.new(0, 0, 0,
		m[1], m[2], m[3],
		m[4], m[5], m[6],
		m[7], m[8], m[9]
	)
end

matrix3.components = matrix3.GetComponents
matrix3.inverse = matrix3.Inverse
matrix3.toOrientation = matrix3.ToEulerAnglesYXZ
matrix3.toEulerAnglesYXZ = matrix3.ToEulerAnglesYXZ
matrix3.toEulerAnglesXYZ = matrix3.ToEulerAnglesXYZ
matrix3.toAxisAngle = matrix3.ToAxisAngle
matrix3.lerp = matrix3.Lerp
matrix3.toCFrame = matrix3.ToCFrame

return matrix3