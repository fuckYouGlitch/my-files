--Taro template \o/


print('LOADED TEMPLATE')

TEMPLATE = {}

-- Taro's janky "TEMPLATE 1" implementation for FNF
-- shoutouts to Kade for letting me do this

-- PSYCH IMPROVED VERSION [README]
-- use psych engine version >= 0.6.3
-- there are two callbacks:init,update
-- write mods in init
-- some stuffs were stolen from troll engine
-- EASING EQUATIONS

---------------------------------------------------------------------------------------
----------------------DON'T TOUCH IT KIDDO---------------------------------------------
---------------------------------------------------------------------------------------
			
-- Adapted from
-- Tweener's easing functions (Penner's Easing Equations)
-- and http://code.google.com/p/tweener/ (jstweener javascript version)
--

--[[
Disclaimer for Robert Penner's Easing Equations license:

TERMS OF USE - EASING EQUATIONS

Open source under the BSD License.

Copyright © 2001 Robert Penner
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
    * Neither the name of the author nor the names of contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
]]

-- For all easing functions:
-- t = elapsed time
-- b = begin
-- c = change == ending - beginning
-- d = duration (total time)


local pow = math.pow
local sin = math.sin
local cos = math.cos
local pi = math.pi
local sqrt = math.sqrt
local abs = math.abs
local asin  = math.asin

function linear(t, b, c, d)
  return c * t / d + b
end

function inQuad(t, b, c, d)
  t = t / d
  return c * pow(t, 2) + b
end

function outQuad(t, b, c, d)
  t = t / d
  return -c * t * (t - 2) + b
end

function inOutQuad(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(t, 2) + b
  else
    return -c / 2 * ((t - 1) * (t - 3) - 1) + b
  end
end

function outInQuad(t, b, c, d)
  if t < d / 2 then
    return outQuad (t * 2, b, c / 2, d)
  else
    return inQuad((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inCubic (t, b, c, d)
  t = t / d
  return c * pow(t, 3) + b
end

function outCubic(t, b, c, d)
  t = t / d - 1
  return c * (pow(t, 3) + 1) + b
end

function inOutCubic(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * t * t * t + b
  else
    t = t - 2
    return c / 2 * (t * t * t + 2) + b
  end
end

function outInCubic(t, b, c, d)
  if t < d / 2 then
    return outCubic(t * 2, b, c / 2, d)
  else
    return inCubic((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inQuart(t, b, c, d)
  t = t / d
  return c * pow(t, 4) + b
end

function outQuart(t, b, c, d)
  t = t / d - 1
  return -c * (pow(t, 4) - 1) + b
end

function inOutQuart(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(t, 4) + b
  else
    t = t - 2
    return -c / 2 * (pow(t, 4) - 2) + b
  end
end

function outInQuart(t, b, c, d)
  if t < d / 2 then
    return outQuart(t * 2, b, c / 2, d)
  else
    return inQuart((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inQuint(t, b, c, d)
  t = t / d
  return c * pow(t, 5) + b
end

function outQuint(t, b, c, d)
  t = t / d - 1
  return c * (pow(t, 5) + 1) + b
end

function inOutQuint(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(t, 5) + b
  else
    t = t - 2
    return c / 2 * (pow(t, 5) + 2) + b
  end
end

function outInQuint(t, b, c, d)
  if t < d / 2 then
    return outQuint(t * 2, b, c / 2, d)
  else
    return inQuint((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inSine(t, b, c, d)
  return -c * cos(t / d * (pi / 2)) + c + b
end

function outSine(t, b, c, d)
  return c * sin(t / d * (pi / 2)) + b
end

function inOutSine(t, b, c, d)
  return -c / 2 * (cos(pi * t / d) - 1) + b
end

function outInSine(t, b, c, d)
  if t < d / 2 then
    return outSine(t * 2, b, c / 2, d)
  else
    return inSine((t * 2) -d, b + c / 2, c / 2, d)
  end
end

function inExpo(t, b, c, d)
  if t == 0 then
    return b
  else
    return c * pow(2, 10 * (t / d - 1)) + b - c * 0.001
  end
end

function outExpo(t, b, c, d)
  if t == d then
    return b + c
  else
    return c * 1.001 * (-pow(2, -10 * t / d) + 1) + b
  end
end

function inOutExpo(t, b, c, d)
  if t == 0 then return b end
  if t == d then return b + c end
  t = t / d * 2
  if t < 1 then
    return c / 2 * pow(2, 10 * (t - 1)) + b - c * 0.0005
  else
    t = t - 1
    return c / 2 * 1.0005 * (-pow(2, -10 * t) + 2) + b
  end
end

function outInExpo(t, b, c, d)
  if t < d / 2 then
    return outExpo(t * 2, b, c / 2, d)
  else
    return inExpo((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inCirc(t, b, c, d)
  t = t / d
  return(-c * (sqrt(1 - pow(t, 2)) - 1) + b)
end

function outCirc(t, b, c, d)
  t = t / d - 1
  return(c * sqrt(1 - pow(t, 2)) + b)
end

function inOutCirc(t, b, c, d)
  t = t / d * 2
  if t < 1 then
    return -c / 2 * (sqrt(1 - t * t) - 1) + b
  else
    t = t - 2
    return c / 2 * (sqrt(1 - t * t) + 1) + b
  end
end

function outInCirc(t, b, c, d)
  if t < d / 2 then
    return outCirc(t * 2, b, c / 2, d)
  else
    return inCirc((t * 2) - d, b + c / 2, c / 2, d)
  end
end

function inElastic(t, b, c, d, a, p)
  if t == 0 then return b end

  t = t / d

  if t == 1  then return b + c end

  if not p then p = d * 0.3 end

  local s

  if not a or a < abs(c) then
    a = c
    s = p / 4
  else
    s = p / (2 * pi) * asin(c/a)
  end

  t = t - 1

  return -(a * pow(2, 10 * t) * sin((t * d - s) * (2 * pi) / p)) + b
end

-- a: amplitud
-- p: period
function outElastic(t, b, c, d, a, p)
  if t == 0 then return b end

  t = t / d

  if t == 1 then return b + c end

  if not p then p = d * 0.3 end

  local s

  if not a or a < abs(c) then
    a = c
    s = p / 4
  else
    s = p / (2 * pi) * asin(c/a)
  end

  return a * pow(2, -10 * t) * sin((t * d - s) * (2 * pi) / p) + c + b
end

-- p = period
-- a = amplitud
function inOutElastic(t, b, c, d, a, p)
  if t == 0 then return b end

  t = t / d * 2

  if t == 2 then return b + c end

  if not p then p = d * (0.3 * 1.5) end
  if not a then a = 0 end

  local s

  if not a or a < abs(c) then
    a = c
    s = p / 4
  else
    s = p / (2 * pi) * asin(c / a)
  end

  if t < 1 then
    t = t - 1
    return -0.5 * (a * pow(2, 10 * t) * sin((t * d - s) * (2 * pi) / p)) + b
  else
    t = t - 1
    return a * pow(2, -10 * t) * sin((t * d - s) * (2 * pi) / p ) * 0.5 + c + b
  end
end

-- a: amplitud
-- p: period
function outInElastic(t, b, c, d, a, p)
  if t < d / 2 then
    return outElastic(t * 2, b, c / 2, d, a, p)
  else
    return inElastic((t * 2) - d, b + c / 2, c / 2, d, a, p)
  end
end

function inBack(t, b, c, d, s)
  if not s then s = 1.70158 end
  t = t / d
  return c * t * t * ((s + 1) * t - s) + b
end

function outBack(t, b, c, d, s)
  if not s then s = 1.70158 end
  t = t / d - 1
  return c * (t * t * ((s + 1) * t + s) + 1) + b
end

function inOutBack(t, b, c, d, s)
  if not s then s = 1.70158 end
  s = s * 1.525
  t = t / d * 2
  if t < 1 then
    return c / 2 * (t * t * ((s + 1) * t - s)) + b
  else
    t = t - 2
    return c / 2 * (t * t * ((s + 1) * t + s) + 2) + b
  end
end

function outInBack(t, b, c, d, s)
  if t < d / 2 then
    return outBack(t * 2, b, c / 2, d, s)
  else
    return inBack((t * 2) - d, b + c / 2, c / 2, d, s)
  end
end

function outBounce(t, b, c, d)
  t = t / d
  if t < 1 / 2.75 then
    return c * (7.5625 * t * t) + b
  elseif t < 2 / 2.75 then
    t = t - (1.5 / 2.75)
    return c * (7.5625 * t * t + 0.75) + b
  elseif t < 2.5 / 2.75 then
    t = t - (2.25 / 2.75)
    return c * (7.5625 * t * t + 0.9375) + b
  else
    t = t - (2.625 / 2.75)
    return c * (7.5625 * t * t + 0.984375) + b
  end
end

function inBounce(t, b, c, d)
  return c - outBounce(d - t, 0, c, d) + b
end

function inOutBounce(t, b, c, d)
  if t < d / 2 then
    return inBounce(t * 2, 0, c, d) * 0.5 + b
  else
    return outBounce(t * 2 - d, 0, c, d) * 0.5 + c * .5 + b
  end
end

function outInBounce(t, b, c, d)
  if t < d / 2 then
    return outBounce(t * 2, b, c / 2, d)
  else
    return inBounce((t * 2) - d, b + c / 2, c / 2, d)
  end
end
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pi = 3.14159265358979323846264338
rad = pi/180
function instant()
	return 1
end

function scale(x, l1, h1, l2, h2)
	return (((x) - (l1)) * ((h2) - (l2)) / ((h1) - (l1)) + (l2))
end

function math.clamp(val,min,max)
	if val < min then return min end
	if val > max then return max end
	return val
end

function lerp(a, b, t)
    return a + t * (b-a)
end
function square(angle)
	local fAngle = angle % (math.pi * 2)
		--Hack: This ensures the hold notes don't flicker right before they're hit.
		if fAngle < 0.01 then
		    fAngle = fAngle + math.pi * 2;
		end
	return fAngle >= math.pi and -1.0 or 1.0;
end

function triangle( angle )
	local fAngle= angle % math.pi * 2.0;
	if fAngle < 0.0 then
		fAngle= fAngle+math.pi * 2.0;
	end
	local result= fAngle * (1 / math.pi);
	if result < .5 then
		return result * 2.0;
	elseif result < 1.5 then
		return 1.0 - ((result - .5) * 2.0);
	else
		return -4.0 + (result * 2.0);
	end
	
end

function quantize(f,interval)
return math.floor((f+interval/2)/interval)*interval;
end
function round(v)
return math.floor(v+0.5)
end
function scale(x, l1, h1, l2, h2)
return ((x - l1) * (h2 - l2) / (h1 - l1) + l2)
end
function BeatToNoteRow(beat)
return round(beat * 48)
end
function fastTan(x)
return math.sin(x) / math.cos(x)
end
function fastCsc(x)
return 1/math.sin(x)
end
function selectTanType(angle,is_csc)
if is_csc ~= 0 then
return fastCsc(angle)
else
return fastTan(angle)
end
end
function cosecant(a)
return fastCsc(a)
end
function CalculateBumpyAngle(y_offset, offset, period)
  return (y_offset + (100.0 * offset)) / ((period * 16.0) + 16.0)
end
---------------vector?
function vec3(xp,yp,zp)--Vector3D
  v = {x=xp, y=yp, z=zp, w=0}
 return v
end
function vec4(xp,yp,zp,wp)--Vector4D
  v = {x=xp, y=yp, z=zp, w=wp}
 return v
end
function fromEuler(roll, pitch, yaw)
 cr = math.cos(roll * rad);
 sr = math.sin(roll * rad);
 cp = math.cos(pitch * rad);
 sp = math.sin(pitch * rad);
 cy = math.cos(yaw * rad);
 sy = math.sin(yaw * rad);
 q= {x=0, y=0, z=0, w=0 }
 
   q.w = cr * cp * cy + sr * sp * sy;
   q.x = sr * cp * cy - cr * sp * sy;
   q.y = cr * sp * cy + sr * cp * sy;
   q.z = cr * cp * sy - sr * sp * cy;
 return q;--fuck
end
function normalize(q)
   length = math.sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
   
    q.w = q.w / length;
    q.x = q.x / length;
    q.y = q.y / length;
    q.z = q.z / length;
   return q;
end
function CalculateDigitalAngle(y_offset, offset, period)
ARROW_SIZE = 112
return math.pi * (y_offset + (1.0 * offset)) / (ARROW_SIZE + (period * ARROW_SIZE))
end
function getCartesianCoords3D(theta, phi, radius)
local pos = {x=0,y=0,z=0}
rad = math.pi/180
pos.x = math.cos(theta*rad)*math.sin(phi*rad);
pos.y = math.cos(phi*rad);
pos.z = math.sin(theta*rad)*math.sin(phi*rad);
pos.x =pos.x* radius;
pos.y =pos.y* radius;
pos.z =pos.z* radius;
return pos;
end
function rotateXYZ( rX, rY, rZ )
    local PI=pi
	rX = rX*(PI/180)
	rY = rY*(PI/180)
	rZ = rZ*(PI/180)

	local cX = math.cos(rX)
	local sX = math.sin(rX)
	local cY = math.cos(rY)
	local sY = math.sin(rY)
	local cZ = math.cos(rZ)
	local sZ = math.sin(rZ)

	 return {
	 	m00=cZ*cY,m01=cZ*sY*sX+sZ*cX, m02=cZ*sY*cX+sZ*(-sX), m03=0,
	 	m10=(-sZ)*cY, m11=(-sZ)*sY*sX+cZ*cX, m12=(-sZ)*sY*cX+cZ*(-sX), m13=0,
	 	m20=-sY, m21=cY*sX, m22=cY*cX, m23=0,
	 	m30=0, m31=0, m32=0, m33=1
	  }--fuck mat4
end
function rotationXYZ( rX, rY, rZ )rotateXYZ( rX, rY, rZ )end
function rotation3d(axis, angle) 
	axis = normalize(axis);
	 s = math.sin(angle);
	 c = math.cos(angle);
	 oc = 1.0 - c;

	return {
		oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
		oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
		oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
		0.0,                                0.0,                                0.0,                                1.0
	}---又是mat 4 lol
end
function rotate(xx,yy,angle)
return {x=xx*math.cos(angle)-yy*math.sin(angle),y=xx*math.sin(angle)+yy*math.cos(angle)}
end
function rotateV3(vec, xA, yA, zA)
local rotateZ = rotate(vec.x, vec.y, zA);
local rotateY = rotate(rotateZ.x, vec.z, yA);
local rotateX = rotate(rotateY.y, rotateZ.y, xA);
return {x=rotateY.x, y=rotateX.y, z=rotateX.x}
end
function arotateV3(vec, xA, yA, zA)
local rotateZ = rotate(vec.x, vec.y, zA);
        local offZ = {x=rotateZ.x, y=rotateZ.y, z=vec.z}

        local rotateX = rotate(offZ.z, offZ.y, xA);
        local offX = {x=offZ.x, y=rotateX.y, z=rotateX.x}

        local rotateY = rotate(offX.x, offX.z, yA);
        local offY = {x=rotateY.x, y=offX.y, z=rotateY.y}
return offY
end
---------------------------------------------------------------------------------------
----------------------END DON'T TOUCH IT KIDDO-----------------------------------------
---------------------------------------------------------------------------------------
strumDefaultX = {92,204,316,428,732,844,956,1068}
defaultY = 50
CAMERA_COUNT = 1
beat = 0
ARROW_SIZE = 112
rainbow = false
--all of our mods, with default values
modList = {
	beat = 0,
	flip = 0,
	invert = 0,
	drunk = 0,
	drunky = 0,
	drunkspeed = 1,
	drunkoffset = 1,
	drunkperiod = 1,
	tipsy = 0,
	drunkyspeed = 1,
	drunkyoffset = 1,
	drunkyperiod = 1,
	tipsyspeed = 1,
	tipsyoffset = 1,
	tandrunk = 0,
	tandrunkspeed = 1,
	tandrunkoffset = 1,
	tandrunkperiod = 1,
	tandrunky = 0,
	tandrunkyspeed = 1,
	tandrunkyoffset = 1,
	tandrunkyperiod = 1,
	tipsyx = 0,
	tipsyxspeed = 1,
	tipsyxoffset = 1,
	tantipsy = 0,
	tantipsyspeed = 1,
	tantipsyoffset = 1,
	tantipsyx = 0,
	tantipsyxspeed = 1,
	tantipsyxoffset = 1,
	adrunk = 0, --non conflict accent mod
	atipsy = 0, --non conflict accent mod
	bumpyx = 0,
	bumpyxoffset = 1,
	bumpyxperiod = 1,
	beaty = 0,
	sawtooth = 0,
	sawtoothperiod = 0,
	sawtoothangle = 0,
	sawtoothangleperiod = 0,
	digital = 0,
	digitalsteps = 0,
	digitaloffset = 0,
	digitalperiod = 0,
	square = 0,
	squareoffset = 1,
	squareperiod = 1,
	squarey = 0,
	squareyoffset = 1,
	squareyperiod = 1,
	squareangle = 0,
	squareangleoffset = 1,
	squareangleperiod = 1,
	bounce = 0,
	bounceY = 0,
	bounceoffset = 0,
	bounceperiod = 0,
	xmode = 0,
	tiny = 0,
	zigzag = 0,
	zigzagoffset = 0,
	zigzagperiod = 0,
	swap = 0,
	parabolax = 0,
	parabolay = 0,
	scalex = 0,
	scaley = 0,
	squish = 0,
	stretch = 0,
	movex = 0,
	movey = 0,
	amovex = 0,
	amovey = 0,
	tornado = 0,
	reverse = 0,
	split = 0,
	cross = 0,
	centered = 0,
	dark = 0,
	stealth = 0,
	alpha = 1,
	randomvanish = 0,
	blink = 0,
	confusion = 0,
	dizzy = 0,
	wave = 0,
	brake = 0,
	boost = 0,
	boomerang = 0,
	hidden = 0,
	hiddenoffset = 0,
	sudden = 0,
	suddenoffset = 0,
	alternate = 0,
	camx = 0,
	scale = 0,
	zoom = 0,
	camy = 0,
	rotationz = 0,
	rotationx = 0,
	rotationy = 0,
	camwag = 0,
	xmod = 1, 
	movez = 0,
	amovez = 0,
	inside = 0,
	shrink = 0,
	curveX =0,
	pulseinner = 0,
	pulseouter = 0,
	pulseoffset = 0,
	pulseperiod = 0,
	shrinkmult = 0,
	shrinklinear = 0,
	attenuatey = 0,
	attenuatex = 0,
	wiggle = 0,
	cosecant = 1,
	expand = 0,
	expandperiod = 0,
	tanexpand = 0,
	tanexpandperiod = 0,
	tornadoy = 0,
	randomspeed =0,--用不了
	tanwavex = 0,
	tanwavespeed = 1,
	tanwavey = 0,
	twirl =0,
	roll =0,
	twirlspeed=1,
	rollspeed=1,
	shakynotes=0,
	shakynotesspeed=1,
	shakenotes=0,
	wavex=0,
	waveperiod =0,
	cosecantx=0,
	cosecantsize=1,
	cosecantspeed=1,
	cosecantspacing=1,
	cosecantoffset=1,
	cosecantperiod=1,
	cosecanty=0,
	jump = 0,
	bumpyy = 0,
	bumpyyoffset = 1,
	bumpyyperiod = 1,
	tanbumpyy = 0,
	tanbumpyyoffset =.1,
	tanbumpyyperiod =.1,
	tanbumpyx = 0,
	tanbumpyxoffset =.1,
	tanbumpyxperiod =.1,
    big = 0,
	tandigital = 0,
	tandigitalsteps = 0,
	tandigitaloffset = 0,
	tandigitalperiod = 0,
	bumpy = 0,
	bumpyoffset = 0,
	bumpyperiod = 0,
	tanbumpy = 0,
	tanbumpyoffset = 0,
	tanbumpyperiod = 0,
	attenuatez = 0,
	tornadoz = 0,
	parabolaz = 0,
	sawtoothz = 0,
	sawtoothzperiod = 0,
	sawtoothscalex = 0,
	sawtoothscalexperiod = 0,
	sawtoothscaley = 0,
	sawtoothscaleyperiod = 0,
	squarez = 0,
	squarezoffset = 1,
	squarezperiod = 1,
	bouncez = 0,
	bouncezoffset = 0,
	bouncezperiod = 0,
	digitalz = 0,
	digitalzsteps = 0,
	digitalzoffset = 0,
	digitalzperiod = 0,
	tandigitalz = 0,
	tandigitalzsteps = 10,
	tandigitalzoffset = 0,
	tandigitalzperiod = 10,
	zigzagz = 0,
	zigzagzoffset = 0,
	zigzagzperiod = 0,
	zigzagscalex = 0,
	zigzagscalexoffset = 0,
	zigzagscalexperiod = 0,
	zigzagscaley = 0,
	zigzagscaleyoffset = 0,
	zigzagscaleyperiod = 0,
	sawtoothy= 0,
	sawtoothyperiod = 0,
	beatz = 0,
	drunkz = 0,
	drunkzspeed = 0,
	drunkzoffset = 0,
	drunkzperiod = 0,
	tandrunkz = 0,
	tandrunkzspeed = 0,
	tandrunkzoffset = 0,
	tandrunkzperiod = 0,
	tipsyz = 0,
	tipsyzspeed = 0,
	tipsyzoffset = 0,
	tantipsyz = 0,
	tantipsyzspeed = 0,
	tantipsyzoffset = 0,
	zangle=0,
	------这些用不了--------
	confusionx = 0,
	confusionxoffset = 1,
	confusiony = 0,
	confusionyoffset = 1,
	-----------------------
	confusionoffset = 1,
	rotatex = 0,
	rotatey = 0,
	rotatez = 0,
	rotatePointY= (720/2)-(ARROW_SIZE/2),
	rotatePointX= (1280/2)-(ARROW_SIZE/2),
	incomingAngleX = 0,
	incomingAngleY= 0,
	sine = 0,
	distortwiggle = 0,
	distortwigglescratch = 0,
	distortwiggleperiod = 0,
	receptorscroll = 0,-----用不了
	vibrate = 0,
	outside = 0,
	zigzagy = 0,
	zigzagyoffset = 0,
	zigzagyperiod = 0,
	zigzagangle = 0,
	zigzagangleoffset = 0,
	zigzagangleperiod = 0,
	arotatex = 0,
	arotatey = 0,-----有漏洞
	arotatez = 0,
	fieldroll = 0,
	fieldpitch = 0,
	fieldyaw = 0,
	fieldx = 0,
	fieldy = 0,
	fieldz = 0,
	centerrotatex = 0,
	centerrotatey = 0,
	centerrotatez = 0,
	InvertIncomingAngle=0,
	IncomingAngleSmooth=0,
	IncomingAngleCurve=0,
	MoveYWaveShit=0,
	beatscalex=0,
	beatscaley=0,
	beatangle=0,    
	drunkscalexoffset=0,
    drunkscalex=0,
	drunkscalexperiod=0,
	drunkscalexspeed=0,
	drunkscaleyoffset=0,
    drunkscaley=0,
	drunkscaleyperiod=0,
	drunkscaleyspeed=0,
	drunkangleoffset=0,
    drunkangle=0,
	drunkangleperiod=0,
	drunkanglespeed=0,
	Paralysis=0,
	Paralysisamplitude=1,
	tipsyscalex =0,
	tipsyscaley =0,
	tipsyangle=0,
    tipsyanglespeed = 1,
	tipsyangleoffset = 1,
	tipsyscaleyspeed = 1,
	tipsyscaleyoffset = 1,
	tipsyscalexspeed = 1,
	tipsyscalexoffset = 1,
	camalpha = 1,
	camzoom = 0,
	localrotatex = 0,
	localrotatey = 0,
	localrotatez = 0,
	danceshit = 0,
	xzfollowstrums=0,
	xzfollowstrumsspeed=1,
}

--column specific mods
for i=0,3 do
	modList['movex'..i] = 0
	modList['movey'..i] = 0
	modList['amovex'..i] = 0
	modList['amovey'..i] = 0
	modList['dark'..i] = 0
	modList['stealth'..i] = 0
	modList['confusion'..i] = 0
	modList['reverse'..i] = 0
	modList['scalex'..i] = 0
	modList['scaley'..i] = 0
	modList['squish'..i] = 0
	modList['rotatex'..i] = 0
	modList['rotatey'..i] = 0
	modList['rotatez'..i] = 0
	modList['arotatex'..i] = 0
	modList['arotatey'..i] = 0
	modList['arotatez'..i] = 0
	modList['stretch'..i] = 0
	modList['xmod'..i] = 1 
	modList['tiny'..i] = 0
	modList['scale'..i] = 0
	modList['zoom'..i] = 0
	modList['movez'..i] = 0
	modList['bumpy'..i] = 0
	modList['amovez'..i] = 0
	modList['confusionx'..i] = 0
	modList['confusionxoffset'..i] = 0
	modList['confusiony'..i] = 0
	modList['confusionyoffset'..i] = 0
	modList['confusionoffset'..i] = 0
	modList['centerrotatex'..i] = 0
	modList['centerrotatey'..i] = 0
	modList['centerrotatez'..i] = 0
	modList['localrotatex'..i] = 0
	modList['localrotatey'..i] = 0
	modList['localrotatez'..i] = 0
end

activeMods = {{},{}}

for pn=1,2 do
	for k,v in pairs(modList) do
		activeMods[pn][k] = v
	end
end

storedMods = {{},{}}
targetMods = {{},{}}
isTweening = {{},{}}
tweenStart = {{},{}}
tweenLen = {{},{}}
tweenCurve = {{},{}}
tweenEx1 = {{},{}}
tweenEx2 = {{},{}}
modnames = {}

function definemod(t)
	local k,v = t[1],t[2]
	if not v then v = 0 end
	for pn=1,2 do
		storedMods[pn][k] = v
		targetMods[pn][k] = v
		isTweening[pn][k] = false
		tweenStart[pn][k] = 0
		tweenLen[pn][k] = 0
		tweenCurve[pn][k] = linear
		tweenEx1[pn][k] = nil
		tweenEx2[pn][k] = nil
		if pn == 1 then
			--print('registered modifier: '..k)
			table.insert(modnames,k)
		end
	end
end

function TEMPLATE.InitMods()
	for pn=1,2 do
		for k,v in pairs(activeMods[pn]) do
			definemod{k,v}
		end
	end
end

function TEMPLATE.setup()
	--sort tables, optimization step
	function modtable_compare(a,b)
		return a[1] < b[1]
	end
	
	if table.getn(event) > 1 then
		table.sort(event, modtable_compare)
	end
	if table.getn(mods) > 1 then
		table.sort(mods, modtable_compare)
	end
end
function receptorAlpha(iCol,pn)
	local alp = 1
	
	local m = activeMods[pn]
	
	if m.alpha ~= 1 then
		alp = alp*m.alpha
	end
	if m.dark ~= 0 or m['dark'..iCol] ~= 0 then
		alp = alp*(1-m.dark)*(1-m['dark'..iCol])
	end
	
	return alp
end

function arrowAlpha(fYOffset, iCol,pn)
	local alp = 1
	
	local m = activeMods[pn]
	
	if m.alpha ~= 1 then
		alp = alp*m.alpha
	end
	if m.stealth ~= 0 or m['stealth'..iCol] ~= 0 then
		alp = alp*(1-m.stealth)*(1-m['stealth'..iCol])
	end
	if m.hidden ~= 0 then
		if fYOffset < m.hiddenoffset and fYOffset >= m.hiddenoffset-200 then
			local hmult = -(fYOffset-m.hiddenoffset)/200
			alp = alp*(1-hmult)*m.hidden
		elseif fYOffset < m.hiddenoffset-100 then
			alp = alp*(1-m.hidden)
		end
	end
	if m.sudden ~= 0 then
		if fYOffset > m.suddenoffset and fYOffset <= m.suddenoffset+200 then
			local hmult = -(fYOffset-m.suddenoffset)/200
			alp = alp*(1-hmult)*m.sudden
		elseif fYOffset > m.suddenoffset+100 then
			alp = alp*(1-m.sudden)
		end
	end
	if m.blink ~= 0 then
		local time = getSongPosition()/1000
		local f = math.sin(time*10)
		f=quantize(f,0.3333)
		alp = alp + scale( f, 0, 1, -1, 0 );
	end
	if m.randomvanish ~= 0 then
		local fRealFadeDist = 80;
		alp = alp + scale( math.abs(fYOffset-360), fRealFadeDist, 2*fRealFadeDist, -1, 0 ) * m.randomvanish;
	end
	return alp
end
function frameCounter(start,ends,interval,func) -- timing must be "end"
local index = start;
while index < ends do
func(index)
index = index + interval
end
end
function getReverseForCol(iCol,pn)
	local m = activeMods[pn]
	local val = 0
	
	val = val+m.reverse+m['reverse'..iCol]
	
	if m.split ~= 0 and iCol > 1 then val = val+m.split end
	if m.cross ~= 0 and iCol == 1 or iCol == 2 then val = val+m.cross end
	if m.alternate ~= 0 and iCol % 2 == 1 then val = val+m.alternate end
	if m.centered ~= 0 then val = scale( m.centered, 0, 1, val, 0 ) end

	return val
end

function getYAdjust(fYOffset, iCol, pn)
	s =getSongPosition()
	local m = activeMods[pn]
	local fScrollSpeed = 1
	local yadj = 0

	
	if m.brake ~= 0 then

		local fEffectHeight = 500;
		local fScale = scale( fYOffset, 0, fEffectHeight, 0, 1 )
		local fNewYOffset = fYOffset * fScale; 
		local fBrakeYAdjust = m.brake * (fNewYOffset - fYOffset)
		
		fBrakeYAdjust = math.clamp( fBrakeYAdjust, -400, 400 )
		yadj = yadj+fBrakeYAdjust;
	
	end

	if m.boost ~= 0 then

		local fEffectHeight = 500;
		local fNewYOffset = fYOffset * 1.5 / ((fYOffset+fEffectHeight/1.2)/fEffectHeight); 
		local fAccelYAdjust = m.boost * (fNewYOffset - fYOffset)
		
		fAccelYAdjust = math.clamp( fAccelYAdjust, -400, 400 )
		yadj = yadj+fAccelYAdjust;
	
	end
	
    fYOffset = fYOffset+yadj
    
	if m.boomerang ~= 0 then
		fYOffset = ((-1*fYOffset*fYOffset/500) + 1.5*fYOffset)*m.boomerang
	end
	if m.expand ~= 0 then
	local last = 0
	local time = getSongPosition() / 1000
	local expandSeconds = 0
    expandSeconds = expandSeconds + (time - last);
    expandSeconds = expandSeconds % ((math.pi * 2) / (m.expandperiod + 1));
    last = time
	    local fExpandMultiplier = scale(math.cos(expandSeconds * 3 * (m.expandperiod + 1)), -1, 1, 0.75, 1.75);
      fScrollSpeed = fScrollSpeed * scale(m.expand, 0, 1, 1, fExpandMultiplier);
	end
	if m.tanexpand ~= 0 then
	local last = 0
	local time = getSongPosition() / 1000
	local tanExpandSeconds = 0
    tanExpandSeconds = tanExpandSeconds + (time - last);
    tanExpandSeconds = tanExpandSeconds % ((math.pi * 2) / (m.tanexpandperiod + 1));
    last = time
	    local fTanExpandMultiplier = scale(selectTanType(tanExpandSeconds * 3 * (m.tanexpandperiod + 1),m.cosecant), -1, 1, 0.75, 1.75);
      fScrollSpeed = fScrollSpeed * scale(m.tanexpand, 0, 1, 1, fTanExpandMultiplier);
	end
	if m.randomspeed ~= 0 then
fr = runHaxeCode([[
var seed:Int = (btr(]]..tostring(curBeat)..[[) << 8) + (]]..tostring(iCol)..[[ * 100);
for (i in 0...3){
seed = ((seed * 1664525) + 1013904223) & 0xFFFFFFFF;
}
var fRandom:Float = seed / 4294967296.0;
return fRandom;]])
      fScrollSpeed = fScrollSpeed*scale(fr, 0.0, 1.0, 1.0, m.randomspeed+ 1.0);
     -- debugPrint(fr)
	end
	if m.Paralysis ~=0 then
	 beat = (s/crochet/2);
       fixedperiod = (math.floor(beat)*crochet*2);
        strumTime = (s - (fYOffset / scrollSpeed));
      fYOffset=((fixedperiod - strumTime)*scrollSpeed/4)*m.Paralysisamplitude;
        
    end
	
	
	fYOffset = fYOffset * fScrollSpeed
	return fYOffset
end

function getZoom(fYOffset, iCol, pn)
	local m = activeMods[pn]
    local fZoom = 1
    
	local fPulseInner = 1.0
	if m.pulseinner ~= 0 or m.pulseouter ~= 0 then
		fPulseInner = ((m.pulseinner*0.5)+1)
		if fPulseInner == 0 then
			fPulseInner = 0.01
		end
	end
    if m.pulseinner ~= 0 or m.pulseouter ~= 0 then
    local sine = math.sin(((fYOffset+(100.0*m.pulseoffset))/(0.4*(ARROW_SIZE+(m.pulseperiod*ARROW_SIZE)))))

		fZoom = fZoom*((sine*(m.pulseouter*0.5))+fPulseInner)
	end
	if m.shrinkmult ~= 0 and fYOffset >= 0 then
	fZoom = fZoom * (1/(1+(fYOffset*(m.shrinkmult/100.0))))
    end
    if m.shrinklinear ~= 0 and fYOffset >= 0 then
    fZoom = fZoom + (fYOffset*(0.5*m.shrinklinear/ARROW_SIZE))
    end
	local fTinyPercent = m.tiny
	if fTinyPercent ~= 0 then
		fTinyPercent = math.pow( 0.5, fTinyPercent )
		fZoom = fZoom * fTinyPercent
	end
	if m['tiny'..iCol] ~= 0 then
		fTinyPercent = math.pow( 0.5,  m['tiny'..iCol] )
		fZoom = fZoom * fTinyPercent
	end
	fZoom = fZoom + m.zoom + m['zoom'..iCol]
	return fZoom
end
function getScale(fYOffset, iCol, pn, sx, sy)
    local x = sx
    local y = sy
    local m = activeMods[pn]
x = x + m.scalex + m['scalex'..iCol] + m.scale + m['scale'..iCol] + m.big
y = y + m.scaley + m['scaley'..iCol] + m.scale + m['scale'..iCol] + m.big
    local angle = 0;
    local stretch = m.stretch + m['stretch'..iCol]
    local squish = m.squish + m['squish'..iCol]
    local stretchX = lerp(1, 0.5, stretch);
    local stretchY = lerp(1, 2, stretch);
    local squishX = lerp(1, 2, squish);
    local squishY = lerp(1, 0.5, squish);
x = x * ((math.sin(angle * math.pi / 180) * squishY) + (math.cos(angle * math.pi / 180) * squishX));
x = x * ((math.sin(angle * math.pi / 180) * stretchY) + (math.cos(angle * math.pi / 180) * stretchX));
y = y * ((math.cos(angle * math.pi / 180) * stretchY) + (math.sin(angle * math.pi / 180) * stretchX));
y = y * ((math.cos(angle * math.pi / 180) * squishY) + (math.sin(angle * math.pi / 180) * squishX));
if m.twirl ~=0 then
x =x*(0+(m.twirl*math.cos(((fYOffset*0.001)*(5*m.twirlspeed)))));
end
if m.roll ~=0 then
y =y*(0+(m.roll*math.cos(((fYOffset*0.001)*(5*m.rollspeed)))));
end
if m.beatscalex ~= 0 then

		local fBeatStrength = m.beatscalex;

		local fAccelTime = 0.3;
		local fTotalTime = 0.7;

		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = beat + fAccelTime;

		local bEvenBeat = false;
		if math.floor(fBeat) % 2 ~= 0 then
			bEvenBeat = true;
		end

		fBeat = fBeat-math.floor( fBeat );
		fBeat = fBeat+1;
		fBeat = fBeat-math.floor( fBeat );

		if fBeat<fTotalTime then

			local fAmount = 0;
			if fBeat < fAccelTime then
				fAmount = scale( fBeat, 0.0, fAccelTime, 0.0, 1.0);
				fAmount = fAmount*fAmount;
			else 
				--fBeat < fTotalTime
				fAmount = scale( fBeat, fAccelTime, fTotalTime, 1.0, 0.0);
				fAmount = 1 - (1-fAmount) * (1-fAmount);
			end

			if bEvenBeat then
				fAmount = fAmount*-1;
			end

			local fShift = 40.0*fAmount*math.sin( ((fYOffset/30.0)) + (math.pi/2) );

			x = x + fBeatStrength * fShift/100

		end

	end
	if m.beatscaley ~= 0 then

		local fBeatStrength = m.beatscaley;

		local fAccelTime = 0.3;
		local fTotalTime = 0.7;

		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = beat + fAccelTime;

		local bEvenBeat = false;
		if math.floor(fBeat) % 2 ~= 0 then
			bEvenBeat = true;
		end

		fBeat = fBeat-math.floor( fBeat );
		fBeat = fBeat+1;
		fBeat = fBeat-math.floor( fBeat );

		if fBeat<fTotalTime then

			local fAmount = 0;
			if fBeat < fAccelTime then
				fAmount = scale( fBeat, 0.0, fAccelTime, 0.0, 1.0);
				fAmount = fAmount*fAmount;
			else 
				--fBeat < fTotalTime
				fAmount = scale( fBeat, fAccelTime, fTotalTime, 1.0, 0.0);
				fAmount = 1 - (1-fAmount) * (1-fAmount);
			end

			if bEvenBeat then
				fAmount = fAmount*-1;
			end

			local fShift = 40.0*fAmount*math.sin( ((fYOffset/30.0)) + (math.pi/2) );

		y = y + fBeatStrength * fShift/100

		end

	end
	if m.drunkscalex ~= 0 then
        x = x+ m.drunkscalex * math.cos(getSongPosition()*0.001 * (1 + m.drunkscalexspeed) + iCol * ((m.drunkscalexoffset * 0.2) + 0.2) + fYOffset * ((m.drunkscalexperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5/100
    end
    if m.drunkscaley ~= 0 then
        y = y+ m.drunkscaley * math.cos(getSongPosition()*0.001 * (1 + m.drunkscaleyspeed) + iCol * ((m.drunkscaleyoffset * 0.2) + 0.2) + fYOffset * ((m.drunkscaleyperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5/100
    end
    if m.tipsyscalex ~= 0 then
        x = x + m.tipsyscalex * math.cos( getSongPosition() * 0.001 * ((m.tipsyscalexspeed * 1.2) + 1.2) + (iCol * ((m.tipsyscalexoffset * 1.8) + 1.8))) * ARROW_SIZE * 0.4/100
    end
    if m.tipsyscaley ~= 0 then
        y = y + m.tipsyscaley * math.cos( getSongPosition() * 0.001 * ((m.tipsyscaleyspeed * 1.2) + 1.2) + (iCol * ((m.tipsyscaleyoffset * 1.8) + 1.8))) * ARROW_SIZE * 0.4/100
    end
    if m.zigzagscalex ~= 0 then
		local fResult = triangle( (math.pi * (1/(m.zigzagscalexperiod+1)) * ((fYOffset+(100.0*(m.zigzagscalexoffset)))/ARROW_SIZE) ) );
		x = x + (m.zigzagscalex*ARROW_SIZE/2) * fResult/100;
	end
	if m.zigzagscaley ~= 0 then
		local fResult = triangle( (math.pi * (1/(m.zigzagscaleyperiod+1)) * ((fYOffset+(100.0*(m.zigzagscaleyoffset)))/ARROW_SIZE) ) );
		y = y + (m.zigzagscaley*ARROW_SIZE/2) * fResult/100;
	end
	if m.sawtoothscaley ~= 0 then
		y = y + (m.sawtoothscaley*ARROW_SIZE) * ((0.5 / (m.sawtoothscaleyperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothscaleyperiod+1) * fYOffset) / ARROW_SIZE) )/100;
	end
	if m.sawtoothscalex ~= 0 then
		x = x + (m.sawtoothscalex*ARROW_SIZE) * ((0.5 / (m.sawtoothscalexperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothscalexperiod+1) * fYOffset) / ARROW_SIZE) )/100;
	end
return x,y
end
function receptorRotation(fYOffset,iCol,pn)
    local fRotationX, fRotationY, fRotationZ = 0, 0, 0
    local m = activeMods[pn]
    if m['confusionx'..iCol]~= 0 then
		fRotationX = fRotationX +m['confusionx'..iCol] * 180.0/math.pi;
	end
	if m.confusionxoffset ~= 0 then
		fRotationX = fRotationX +m.confusionxoffset * 180.0/math.pi;
	end
	if m['confusionxoffset'..iCol] ~= 0 then
		fRotationX = fRotationX +m['confusionxoffset'..iCol] * 180.0/math.pi;
	end
	if m.confusionx ~= 0 then
		local fConfRotation = beat
		local PI = math.pi
		fConfRotation = fConfRotation * m.confusionx
		fConfRotation = fConfRotation % ( 2*PI )
		fConfRotation = fConfRotation*(-180/PI)
		fRotationX = fRotationX + fConfRotation;
	end
	if m['confusiony'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusiony'..iCol] * 180.0/math.pi)
	end
	if m['confusionyoffset'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusionyoffset'..iCol] * 180.0/math.pi)
	end
	if m.confusionyoffset ~= 0 then
		fRotationY = fRotationY + m.confusionyoffset * 180.0/math.pi
	end
	if m.confusiony ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusiony
		fConfRotation = fConfRotation%(2*PI )
		fConfRotation = fConfRotation *(-180/PI)
		fRotationY = fRotationY + fConfRotation
	end
	if m['confusion'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusion'..iCol] * 180.0/math.pi
	end
	if m.confusionoffset ~= 0 then
		fRotationZ = fRotationZ+m.confusionoffset * 180.0/math.pi
	end
	if m['confusionoffset'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusionoffset'..iCol] * 180.0/math.pi
	end
	if m.confusion ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusion
		fConfRotation = fConfRotation%( 2*PI );
		fConfRotation = fConfRotation*(-180/PI)
		fRotationZ = fRotationZ*fConfRotation;
	end
    return fRotationX, fRotationY, fRotationZ
end
function arrowRotation(fYOffset,iCol,pn,noteBeat)
    local fRotationX, fRotationY, fRotationZ = 0, 0, 0
    local m = activeMods[pn]
    if m['confusionx'..iCol]~= 0 then
		fRotationX = fRotationX +m['confusionx'..iCol] * 180.0/math.pi;
	end
	if m.confusionxoffset ~= 0 then
		fRotationX = fRotationX +confusionxoffset * 180.0/math.pi;
	end
	if m['confusionxoffset'..iCol] ~= 0 then
		fRotationX = fRotationX +m['confusionxoffset'..iCol] * 180.0/math.pi;
	end
	if m.confusionx ~= 0 then
		local fConfRotation = beat
		fConfRotation = fConfRotation * m.confusionx
		fConfRotation = fConfRotation % ( 2*PI )
		fConfRotation = fConfRotation*(-180/PI)
		fRotationX = fRotationX + fConfRotation;
	end
	if m['confusiony'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusiony'..iCol] * 180.0/math.pi)
	end
	if m['confusionyoffset'..iCol] ~= 0 then
		fRotationY = fRotationY + (m['confusionyoffset'..iCol] * 180.0/math.pi)
	end
	if m.confusionyoffset ~= 0 then
		fRotationY = fRotationY + m.confusionyoffset * 180.0/math.pi
	end
	if m.confusiony ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusiony
		fConfRotation = fConfRotation%(2*PI )
		fConfRotation = fConfRotation *(-180/PI)
		fRotationY = fRotationY + fConfRotation
	end
	if m.dizzy ~= 0 then
	local fSongBeat = beat
	local PI = pi
	local fDizzyRotation = noteBeat - fSongBeat
		fDizzyRotation = fDizzyRotation * m.dizzy
		fDizzyRotation = fDizzyRotation % (2*PI)
		fDizzyRotation = fDizzyRotation*(180/PI)
		fRotationZ = fRotationZ+fDizzyRotation
	end
		if m['confusion'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusion'..iCol] * 180.0/math.pi
	end
	if m.confusionoffset ~= 0 then
		fRotationZ = fRotationZ+m.confusionoffset * 180.0/math.pi
	end
	if m['confusionoffset'..iCol] ~= 0 then
		fRotationZ = fRotationZ+m['confusionoffset'..iCol] * 180.0/math.pi
	end
	if m.confusion ~= 0 then
		local fConfRotation = beat
		local PI=math.pi
		fConfRotation = fConfRotation * m.confusion
		fConfRotation = fConfRotation%( 2*PI );
		fConfRotation = fConfRotation*(-180/PI)
		fRotationZ = fRotationZ*fConfRotation;
	end
    return fRotationX, fRotationY, fRotationZ
end
function addModchartingToolsmod(curPos, lane, pf)
local m = activeMods[pf]
local noteData = {x=0,y=0,z=0,angle=0}
s = getSongPosition()
-------------------
if m.danceshit~= 0 then
local shift = m.danceshit*30;
if (lane % 2 == 0) then
shift =shift* -1;
end
noteData.angle =noteData.angle+ m.danceshit*shift*0.35;
noteData.y =noteData.y+ m.danceshit*shift;    
end
-------------------
return noteData
end
function arrowEffects(fYOffset, iCol, pn)
    local m = activeMods[pn]
	local xpos, ypos, rotz, zpos = 0, 0, 0, 0
	s = getSongPosition()
	if m.rotatex ~= 0 or m.rotatey~= 0 or m.rotatez ~= 0 then	
    local laneShit = iCol%4;
    local offsetThing = 0.5
        if (iCol < 2) then
            offsetThing = -0.5;
            laneShit = iCol+1;
        end
    local distFromCenter = ((laneShit)-2)+offsetThing;
    xpos = xpos-distFromCenter*ARROW_SIZE;
    local q = fromEuler(90+m.rotatez, m.rotatex, (downscroll and m.rotatey or -m.rotatey))
    xpos= xpos+q.x * distFromCenter*ARROW_SIZE;
    ypos= ypos+q.y * distFromCenter*ARROW_SIZE;
    zpos= zpos+q.z * distFromCenter*ARROW_SIZE;
    end
    if m['confusion'..iCol] ~= 0 or m.confusion ~= 0 ~= 0 then
		rotz = rotz + m['confusion'..iCol] + m.confusion
	end
	if m.dizzy ~= 0 then
		rotz = rotz + m.dizzy*fYOffset
	end
	if m.beatangle ~= 0 then

		local fBeatStrength = m.beatangle;

		local fAccelTime = 0.3;
		local fTotalTime = 0.7;

		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = beat + fAccelTime;

		local bEvenBeat = false;
		if math.floor(fBeat) % 2 ~= 0 then
			bEvenBeat = true;
		end

		fBeat = fBeat-math.floor( fBeat );
		fBeat = fBeat+1;
		fBeat = fBeat-math.floor( fBeat );

		if fBeat<fTotalTime then

			local fAmount = 0;
			if fBeat < fAccelTime then
				fAmount = scale( fBeat, 0.0, fAccelTime, 0.0, 1.0);
				fAmount = fAmount*fAmount;
			else 
				--fBeat < fTotalTime
				fAmount = scale( fBeat, fAccelTime, fTotalTime, 1.0, 0.0);
				fAmount = 1 - (1-fAmount) * (1-fAmount);
			end

			if bEvenBeat then
				fAmount = fAmount*-1;
			end

			local fShift = 40.0*fAmount*math.sin( ((fYOffset/30.0)) + (math.pi/2) );

			rotz = rotz + fBeatStrength * fShift

		end

	end
	if m.tipsyangle ~= 0 then
        rotz = rotz + m.tipsyangle * math.cos( getSongPosition() * 0.001 * ((m.tipsyanglespeed * 1.2) + 1.2) + (iCol * ((m.tipsyangleoffset * 1.8) + 1.8))) * ARROW_SIZE * 0.4
    end
	    if m.tornado ~= 0 then
		local iTornadoWidth = 2

		local iStartCol = iCol - iTornadoWidth;
		local iEndCol = iCol + iTornadoWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMinX = 3.402823466*(10^38)
		local fMaxX = 1.175494351*(10^-38)
		
		-- TODO: Don't index by PlayerNumber.

		for i=iStartCol,iEndCol do
		
			fMinX = math.min( fMinX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
			fMaxX = math.max( fMaxX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
	end

		local fRealPixelOffset = getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + (fYOffset * 6 / screenHeight)
		
		local fAdjustedPixelOffset = scale( math.cos(fRads), -1, 1, fMinX, fMaxX );

		xpos = xpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornado
    end
    if m.tornadoy ~= 0 then
		local iTornadoyWidth = 2

		local iStartCol = iCol - iTornadoyWidth;
		local iEndCol = iCol + iTornadoyWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMiny = 3.402823466*(10^38)
		local fMaxy = 1.175494351*(10^-38)
		
		-- TODO: Don't index by PlayerNumber.

		for i=iStartCol,iEndCol do
		
			fMiny = math.min( fMiny, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
			fMaxy = math.max( fMaxy, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
	end

		local fRealPixelOffset = getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
		local fPositionBetween = scale( fRealPixelOffset, fMiny, fMaxy, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + (fYOffset * 6 / screenHeight)
		
		local fAdjustedPixelOffset = scale( math.cos(fRads), -1, 1, fMiny, fMaxy );

		ypos = ypos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornadoy
    end
    if m['rotatex'..iCol] ~= 0 or m['rotatey'..iCol] ~= 0 or m['rotatez'..iCol] ~= 0 then	
    local laneShit = iCol%4;
    local offsetThing = 0.5
        if (iCol < 2) then
            offsetThing = -0.5;
            laneShit = iCol+1;
        end
    local distFromCenter = ((laneShit)-2)+offsetThing;
    xpos = xpos-distFromCenter*ARROW_SIZE;
    local q = fromEuler(90+m['rotatez'..iCol], m['rotatex'..iCol], (downscroll and m['rotatey'..iCol] or -m['rotatey'..iCol]))
    xpos= xpos+q.x * distFromCenter*ARROW_SIZE;
    ypos= ypos+q.y * distFromCenter*ARROW_SIZE;
    zpos= zpos+q.z * distFromCenter*ARROW_SIZE;
    end
    if m.drunk ~= 0 then
        xpos = xpos + m.drunk * math.cos(getSongPosition()*0.001 * (1 + m.drunkspeed) + iCol * ((m.drunkoffset * 0.2) + 0.2) + fYOffset * ((m.drunkperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.drunky ~= 0 then
        ypos = ypos + m.drunky * math.cos(getSongPosition()*0.001 * (1 + m.drunkyspeed) + iCol * ((m.drunkyoffset * 0.2) + 0.2) + fYOffset * ((m.drunkyperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.tandrunk ~= 0 then
        xpos = xpos + m.tandrunk * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkspeed) + iCol * ((m.tandrunkoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkperiod * 10) + 10) / screenHeight,m.cosecant) * ARROW_SIZE * 0.5;
    end
    if m.tandrunky ~= 0 then
        ypos = ypos + m.tandrunky * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkyspeed) + iCol * ((m.tandrunkyoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkyperiod * 10) + 10) / screenHeight,m.cosecant) * ARROW_SIZE * 0.5;
    end
    if m.drunkangle ~= 0 then
        rotz = rotz+ m.drunkangle * math.cos(getSongPosition()*0.001 * (1 + m.drunkanglespeed) + iCol * ((m.drunkangleoffset * 0.2) + 0.2) + fYOffset * ((m.drunkangleperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5
    end
    if m.tipsy ~= 0 then
        ypos = ypos + m.tipsy * math.cos( getSongPosition() * 0.001 * ((m.tipsyspeed * 1.2) + 1.2) + (iCol * ((m.tipsyoffset * 1.8) + 1.8)))* ARROW_SIZE * 0.4
    end
    if m.tantipsy ~= 0 then
        ypos = ypos + m.tantipsy * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyoffset * 1.8) + 1.8)),m.cosecant)* ARROW_SIZE * 0.4
    end
    if m.tipsyx ~= 0 then
        xpos = xpos + m.tipsyx * math.cos( getSongPosition() * 0.001 * ((m.tipsyxspeed * 1.2) + 1.2) + (iCol * ((m.tipsyxoffset * 1.8) + 1.8))) * ARROW_SIZE * 0.4
    end
     if m.tantipsyx ~= 0 then
        xpos = xpos + m.tantipsyx * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyxspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyxoffset * 1.8) + 1.8)),m.cosecant)* ARROW_SIZE * 0.4
    end
    if m.adrunk ~= 0 then
        xpos = xpos + m.adrunk * ( math.cos( getSongPosition()*0.001 + iCol*(0.2) + 1*(0.2) + fYOffset*(10)/720) * ARROW_SIZE*0.5 )
    end
    if m.atipsy ~= 0 then
        ypos = ypos + m.atipsy * ( math.cos( getSongPosition()*0.001 *(1.2) + iCol*(2.0) + 1*(0.2) ) * ARROW_SIZE*0.4 )
    end
    if m.drunkz ~= 0 then
        zpos = zpos + m.drunkz * math.cos(getSongPosition()*0.001 * (1 + m.drunkzspeed) + iCol * ((m.drunkzoffset * 0.2) + 0.2) + fYOffset * ((m.drunkzperiod * 10) + 10) / screenHeight) * ARROW_SIZE * 0.5;
    end
    if m.tandrunkz ~= 0 then
        zpos = zpos + m.tandrunkz * selectTanType(getSongPosition()*0.001 * (1 + m.tandrunkzspeed) + iCol * ((m.tandrunkzoffset * 0.2) + 0.2) + fYOffset * ((m.tandrunkzperiod * 10) + 10) / screenHeight,m.cosecant) * ARROW_SIZE * 0.5;
    end
    if m.tipsyz ~= 0 then
        zpos = zpos + m.tipsyz * math.cos( getSongPosition() * 0.001 * ((m.tipsyzspeed * 1.2) + 1.2) + (iCol * ((m.tipsyzoffset * 1.8) + 1.8)))* ARROW_SIZE * 0.4
    end
    if m.tantipsyz ~= 0 then
        zpos = zpos + m.tantipsyz * selectTanType( getSongPosition() * 0.001 * ((m.tantipsyzspeed * 1.2) + 1.2) + (iCol * ((m.tantipsyzoffset * 1.8) + 1.8)),m.cosecant)* ARROW_SIZE * 0.4
    end
    if m.tandigital ~= 0 then
		xpos = xpos + (m.tandigital * ARROW_SIZE * 0.5) * round((m.tandigitalsteps+1) * selectTanType(math.pi * (fYOffset + (1.0 * m.tandigitaloffset ) ) / (ARROW_SIZE + (m.tandigitalperiod * ARROW_SIZE) ),m.cosecant) )/(m.tandigitalsteps+1);
	end
    if m.wave ~= 0 then
		ypos =ypos + m.wave * 20 *math.sin( fYOffset/((m.waveperiod*38)+38) );
	end
	if m.wavex ~= 0 then
		xpos = xpos + m.wavex * 20*math.sin( (fYOffset+250)/76 )
	end
	if m['movex'..iCol] ~= 0 or m.movex ~= 0 then
		xpos = xpos + m['movex'..iCol] + m.movex
	end
	if m['amovex'..iCol] ~= 0 or m.amovex ~= 0 then
		xpos = xpos + m['amovex'..iCol] + m.amovex
	end
	if m['movey'..iCol] ~= 0 or m.movey ~= 0 then
		ypos = ypos + m['movey'..iCol] + m.movey
	end
	if m['amovey'..iCol] ~= 0 or m.amovey ~= 0 then
		ypos = ypos + m['amovey'..iCol] + m.amovey
	end
	if m['movez'..iCol] ~= 0 or m.movez ~= 0 then
	    zpos = zpos + m['movez'..iCol] + m.movez
	end
	if m['amovez'..iCol] ~= 0 or m.amovez ~= 0 then
		zpos = zpos + m['amovez'..iCol] + m.amovez
	end
	if m['reverse'..iCol] ~= 0 or m.reverse ~= 0 or m.split ~= 0 or m.cross ~= 0 or m.alternate ~= 0 or m.centered ~= 0 then
		ypos = ypos + getReverseForCol(iCol,pn) * 470
	end
	if m.bumpy ~= 0 then
		zpos = zpos + m.bumpy * 40*math.sin((fYOffset+(100.0*m.bumpyoffset))/((m.bumpyperiod*16.0)+16.0))
    end
	if m['bumpy'..iCol] ~= 0 then
	    zpos = zpos + m['bumpy'..iCol] * 40*math.sin((fYOffset+(100.0*m.bumpyoffset))/((m.bumpyperiod*16.0)+16.0))
	end
	if m.tanbumpy ~= 0 then
		zpos = zpos + m.tanbumpy * 40*selectTanType((fYOffset+(100.0*m.tanbumpyoffset))/((m.tanbumpyperiod*16.0)+16.0),m.cosecant)
    end
    if m.attenuatez ~= 0 then
    local fXOffset =getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
    zpos = zpos +m.attenuatez * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (fXOffset / ARROW_SIZE);
    end
    if m.zigzagangle ~= 0 then
		local fResult = triangle( (math.pi * (1/(m.zigzagangleperiod+1)) * ((fYOffset+(100.0*(m.zigzagangleoffset)))/ARROW_SIZE) ) );
		rotz = rotz + (m.zigzagangle*ARROW_SIZE/2) * fResult/3.6;
	end
    if m.tornadoz ~= 0 then
		local iTornadoWidth = 2

		local iStartCol = iCol - iTornadoWidth;
		local iEndCol = iCol + iTornadoWidth;
		iStartCol = math.clamp( iStartCol, 0, 4-1 );
		iEndCol = math.clamp( iEndCol, 0, 4-1 );

		local fMinX = 3.402823466*(10^38)
		local fMaxX = 1.175494351*(10^-38)

		-- TODO: Don't index by PlayerNumber.

		for i=iStartCol,iEndCol do

			fMinX = math.min( fMinX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
			fMaxX = math.max( fMaxX, getPropertyFromGroup('strumLineNotes',(pn==2 and i+4 or i),"x") );
	end

		local fRealPixelOffset = getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x")
		local fPositionBetween = scale( fRealPixelOffset, fMinX, fMaxX, -1, 1 );
		local fRads = math.acos( fPositionBetween );
		fRads = fRads + (fYOffset * 6 / screenHeight)

		local fAdjustedPixelOffset = scale( math.cos(fRads), -1, 1, fMinX, fMaxX );

		zpos = zpos + (fAdjustedPixelOffset - fRealPixelOffset) * m.tornadoz
    end
    if m.parabolaz ~= 0 then
        zpos = zpos + m.parabolaz * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE)
    end
	if m.sawtoothz ~= 0 then
		zpos = zpos + (m.sawtoothz*ARROW_SIZE) * ((0.5 / (m.sawtoothzperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothzperiod+1) * fYOffset) / ARROW_SIZE) );
	end
	if m.sawtoothy ~= 0 then
		ypos = ypos + (m.sawtoothy*ARROW_SIZE) * ((0.5 / (m.sawtoothyperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothyperiod+1) * fYOffset) / ARROW_SIZE) );
	end
	if m.digitalz ~= 0 then
		zpos = zpos + (m.digitalz * ARROW_SIZE * 0.5) * round((m.digitalzsteps+1) * math.sin(math.pi * (fYOffset + (1.0 * m.digitalzoffset ) ) / (ARROW_SIZE + (m.digitalzperiod * ARROW_SIZE) )) )/(m.digitalzsteps+1);
	end
	if m.tandigitalz ~= 0 then
		zpos = zpos + (m.tandigitalz * ARROW_SIZE * 0.5) * round((m.tandigitalzsteps+1) * selectTanType(math.pi * (fYOffset + (1.0 * m.tandigitalzoffset ) ) / (ARROW_SIZE + (m.tandigitalzperiod * ARROW_SIZE) ),m.cosecant) )/(m.tandigitalzsteps+1);
	end
	if m.squarez ~= 0 then
		local fResult = square( (math.pi * (fYOffset+(1.0*(m.squarezoffset))) / (ARROW_SIZE+(m.squarezperiod*ARROW_SIZE))) );
		zpos = zpos + (m.squarez * ARROW_SIZE * 0.5) * fResult;
	end
	if m.squareangle ~= 0 then
		local fResult = square( (pi * (fYOffset+(1.0*(m.squareangleoffset))) / (ARROW_SIZE+(m.squareangleperiod*ARROW_SIZE))) );
		rotz = rotz + (m.squareangle * ARROW_SIZE * 0.5) * fResult;
	end
	if m.sawtoothangle ~= 0 then
		rotz = rotz + (m.sawtoothangle*ARROW_SIZE) * ((0.5 / (m.sawtoothangleperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothangleperiod+1) * fYOffset) / ARROW_SIZE) )/3.6;
	end
	if m.bouncez ~= 0 then
		local fBounceAmt = math.abs( math.sin( ( (fYOffset + (1.0 * (m.bouncezoffset) ) ) / ( 60 + m.bouncezperiod*60) ) ) );
		zpos = zpos + m.bouncez * ARROW_SIZE * 0.5 * fBounceAmt;
	end
	if m.zigzagz ~= 0 then
		local fResult = triangle( (math.pi * (1/(m.zigzagzperiod+1)) * ((fYOffset+(100.0*(m.zigzagzoffset)))/ARROW_SIZE) ) );

		zpos = zpos + (m.zigzagz*ARROW_SIZE/2) * fResult;
	end
	if m.zigzagy ~= 0 then
		local fResult = triangle( (math.pi * (1/(m.zigzagyperiod+1)) * ((fYOffset+(100.0*(m.zigzagyoffset)))/ARROW_SIZE) ) );

		ypos = ypos + (m.zigzagy*ARROW_SIZE/2) * fResult;
	end
	if m.rotationx ~= 0 or m.rotationy ~= 0 or rotatePointX ~=(1280/2)-(ARROW_SIZE/2) or rotatePointY~=(720/2)-(ARROW_SIZE/2) then
       x = strumDefaultX[(pn==2 and iCol+5 or iCol+1)]
       y = defaultY
       rotX = getCartesianCoords3D(m.rotationx, 90, x-m.rotatePointX);
       xpos =xpos+ rotX.x+m.rotatePointX-x;
       rotY = getCartesianCoords3D(90,m.rotationy, y-m.rotatePointY);
       ypos =ypos+ rotY.y+m.rotatePointY-y;
       zpos=zpos+ rotX.z + rotY.z;
	end
	if m.beatz ~= 0 then

		local fBeatStrength = m.beatz;

		local fAccelTime = 0.3;
		local fTotalTime = 0.7;

		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = beat + fAccelTime;

		local bEvenBeat = false;
		if math.floor(fBeat) % 2 ~= 0 then
			bEvenBeat = true;
		end

		fBeat = fBeat-math.floor( fBeat );
		fBeat = fBeat+1;
		fBeat = fBeat-math.floor( fBeat );

		if fBeat<fTotalTime then

			local fAmount = 0;
			if fBeat < fAccelTime then
				fAmount = scale( fBeat, 0.0, fAccelTime, 0.0, 1.0);
				fAmount = fAmount*fAmount;
			else 
				--fBeat < fTotalTime
				fAmount = scale( fBeat, fAccelTime, fTotalTime, 1.0, 0.0);
				fAmount = 1 - (1-fAmount) * (1-fAmount);
			end

			if bEvenBeat then
				fAmount = fAmount*-1;
			end

			local fShift = 40.0*fAmount*math.sin( ((fYOffset/30.0)) + (math.pi/2) );

			zpos = zpos + fBeatStrength * fShift

		end

	end
	if m.flip ~= 0 then
		local fDistance = ARROW_SIZE * 2 * (1.5 - iCol);
		xpos = xpos + fDistance * m.flip;
	end

	if m.invert ~= 0 then
		local fDistance = ARROW_SIZE * (iCol%2 == 0 and 1 or -1);
		xpos = xpos + fDistance * m.invert;
	end
	
	if m.beat ~= 0 then
			
		local fBeatStrength = m.beat;
		
		local fAccelTime = 0.3;
		local fTotalTime = 0.7;
		
		-- If the song is really fast, slow down the rate, but speed up the
		-- acceleration to compensate or it'll look weird.
		fBeat = beat + fAccelTime;
		
		local bEvenBeat = false;
		if math.floor(fBeat) % 2 ~= 0 then
			bEvenBeat = true;
		end
		
		fBeat = fBeat-math.floor( fBeat );
		fBeat = fBeat+1;
		fBeat = fBeat-math.floor( fBeat );
		
		if fBeat<fTotalTime then
		
			local fAmount = 0;
			if fBeat < fAccelTime then
				fAmount = scale( fBeat, 0.0, fAccelTime, 0.0, 1.0);
				fAmount = fAmount*fAmount;
			else 
				--fBeat < fTotalTime
				fAmount = scale( fBeat, fAccelTime, fTotalTime, 1.0, 0.0);
				fAmount = 1 - (1-fAmount) * (1-fAmount);
			end

			if bEvenBeat then
				fAmount = fAmount*-1;
			end

			local fShift = 40.0*fAmount*math.sin( ((fYOffset/30.0)) + (pi/2) );
			
			xpos = xpos + fBeatStrength * fShift
			
		end
	
	end

	if m.beaty ~= 0 then
			
		local fBeatStrength = m.beaty;
		
		local fAccelTime = 0.3;
		local fTotalTime = 0.7;
		fBeat = beat + fAccelTime;
		
		local bEvenBeat = false;
		if math.floor(fBeat) % 2 ~= 0 then
			bEvenBeat = true;
		end
		
		fBeat = fBeat-math.floor( fBeat );
		fBeat = fBeat+1;
		fBeat = fBeat-math.floor( fBeat );
		
		if fBeat<fTotalTime then
		
			local fAmount = 0;
			if fBeat < fAccelTime then
				fAmount = scale( fBeat, 0.0, fAccelTime, 0.0, 1.0);
				fAmount = fAmount*fAmount;
			else 
				--fBeat < fTotalTime
				fAmount = scale( fBeat, fAccelTime, fTotalTime, 1.0, 0.0);
				fAmount = 1 - (1-fAmount) * (1-fAmount);
			end

			if bEvenBeat then
				fAmount = fAmount*-1;
			end

			local fShift = 40.0*fAmount*math.sin( ((fYOffset/30.0)) + (pi/2) );
			
			ypos = ypos + fBeatStrength * fShift
			
		end
	
	end

	if m.sawtooth ~= 0 then
		xpos = xpos + (m.sawtooth*ARROW_SIZE) * ((0.5 / (m.sawtoothperiod+1) * fYOffset) / ARROW_SIZE - math.floor((0.5 / (m.sawtoothperiod+1) * fYOffset) / ARROW_SIZE) );
	end

	if m.digital ~= 0 then
		xpos = xpos + (m.digital * ARROW_SIZE * 0.5) * round((m.digitalsteps+1) * math.sin(pi * (fYOffset + (1.0 * m.digitaloffset ) ) / (ARROW_SIZE + (m.digitalperiod * ARROW_SIZE) )) )/(m.digitalsteps+1);
	end

	if m.bumpyx ~= 0 then
		xpos = xpos + m.bumpyx * 40*math.sin((fYOffset+(100.0*m.bumpyxoffset))/((m.bumpyxperiod*16.0)+16.0));
	end
    if m.bumpyy ~= 0 then
		ypos = ypos + m.bumpyy * 40*math.sin((fYOffset+(100.0*m.bumpyyoffset))/((m.bumpyyperiod*16.0)+16.0));
	end
	if m.square ~= 0 then
		local fResult = square( (pi * (fYOffset+(1.0*(m.squareoffset))) / (ARROW_SIZE+(m.squareperiod*ARROW_SIZE))) );
		xpos = xpos + (m.square * ARROW_SIZE * 0.5) * fResult;
	end
	if m.squarey ~= 0 then
		local fResult = square( (pi * (fYOffset+(1.0*(m.squareyoffset))) / (ARROW_SIZE+(m.squareyperiod*ARROW_SIZE))) );
		ypos = ypos + (m.squarey * ARROW_SIZE * 0.5) * fResult;
	end
	if m.bounce ~= 0 then
		local fBounceAmt = math.abs( math.sin( ( (fYOffset + (1.0 * (m.bounceoffset) ) ) / ( 60 + m.bounceperiod*60) ) ) );
		xpos = xpos + m.bounce * ARROW_SIZE * 0.5 * fBounceAmt;
	end

    if m.bounceY ~= 0 then
		local fBounceAmt = math.abs( math.sin( ( (fYOffset + (1.0 * (m.bounceoffset) ) ) / ( 60 + m.bounceperiod*60) ) ) );
		ypos = ypos + m.bounceY * ARROW_SIZE * 0.5 * fBounceAmt;
	end
	if m.xmode ~= 0 then
		xpos = xpos + m.xmode * (pn == 2 and -fYOffset or fYOffset)
	end
    if m.tanbumpyx ~= 0 then
    xpos=xpos+ m.tanbumpyx * 10 * selectTanType(CalculateBumpyAngle(fYOffset/10, m.tanbumpyxoffset,m.tanbumpyxperiod), m.cosecant);
    end
    if m.tanbumpyy ~= 0 then
    ypos=ypos+ m.tanbumpyy * 10 * selectTanType(CalculateBumpyAngle(fYOffset/10, m.tanbumpyyoffset,m.tanbumpyyperiod), m.cosecant);
    end
	if m.tiny ~= 0 then
		local fTinyPercent = m.tiny
		fTinyPercent = math.min( math.pow(0.5, fTinyPercent), 1 );
		xpos = xpos * fTinyPercent
	end

	if m.zigzag ~= 0 then
		local fResult = triangle( (pi * (1/(m.zigzagperiod+1)) * 
		((fYOffset+(100.0*(m.zigzagoffset)))/ARROW_SIZE) ) );
	    
		xpos = xpos + (m.zigzag*ARROW_SIZE/2) * fResult;
	end

    if m.swap ~= 0 then
        xpos = xpos + screenWidth / 2 * m.swap * (pn == 2 and -1 or 1)
    end
    
    if m.parabolax ~= 0 then
        xpos = xpos + m.parabolax * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE)
    end
    
    if m.parabolay ~= 0 then
        ypos = ypos + m.parabolay * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE)
    end
    if m.inside ~= 0 then
        xpos = xpos +math.sin(0 + (fYOffset*0.004))*(ARROW_SIZE* (iCol % 2 == 0 and 1 or -1) * m.inside*0.5);
    end
    if m.curveX ~= 0 then
        xpos = xpos +fYOffset*m.curveX*0.01
    end
    if m.attenuatey ~= 0 then
    ypos = ypos +m.attenuatey * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x") / ARROW_SIZE);
    end
    if m.attenuatex ~= 0 then
    xpos = xpos +m.attenuatex * (fYOffset / ARROW_SIZE) * (fYOffset / ARROW_SIZE) * (getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"x") / ARROW_SIZE);
    end
    if m.outside ~= 0 then
    multIn = iCol%4 <= 1 and -1 or 1
    xpos =xpos+math.pow(math.max((fYOffset*0.01)-1, 0)*m.outside, 2) * multIn
    end
    if m.vibrate ~= 0 then
		xpos = xpos + (math.random() - 0.5) * m.vibrate * 20;
		ypos = ypos + (math.random() - 0.5) * m.vibrate * 20;
    end
    if m.tanwavex ~= 0 then
    xpos = xpos+ 260*m.tanwavex*math.tan((s* (m.tanwavespeed)*0.0008)+(iCol/4))*0.2;
    end
    if m.tanwavey ~= 0 then
    ypos = ypos+ 260*m.tanwavey*math.tan((s* (m.tanwavespeed)*0.0008)+(iCol/4))*0.2;
    end
    if m.shakynotes ~= 0 then
    xpos = xpos+math.sin(500)+m.shakynotes * (math.cos(s * 4*0.2) + ((iCol%4)*0.2) - 0.002)* (math.sin(100 - (120 * m.shakynotesspeed * 0.4)))
    ypos = ypos+math.sin(500)+m.shakynotes * (math.cos(s * 8*0.2) + ((iCol%4)*0.2) - 0.002)* (math.sin(100 - (120 * m.shakynotesspeed * 0.4)))
    end
    if m.shakenotes ~=0 then
    xpos = xpos+math.sin(0.1)*(m.shakenotes * getRandomInt(1, 20));
    ypos = ypos+math.sin(0.1)*(m.shakenotes * getRandomInt(1, 20));
    end
if m.cosecantx~=0 then
xpos = xpos+m.cosecantx * (cosecant( ((s*(0.001*m.cosecantperiod)) + ((iCol%4)*0.2) + (fYOffset*(0.225*m.cosecantoffset))*((m.cosecantspacing*10)/720)) * (m.cosecantoffset*0.2)) * 108*(0.5*m.cosecantsize));
end
if m.cosecanty~=0 then
ypos = ypos+m.cosecanty * (cosecant( ((s*(0.001*m.cosecantperiod)) + ((iCol%4)*0.2) + (fYOffset*(0.225*m.cosecantoffset))*((m.cosecantspacing*10)/720)) * (m.cosecantoffset*0.2)) * 108*(0.5*m.cosecantsize));
end
if m.zangle~=0 then
zpos = fYOffset*4;
end
if m.jump ~= 0 then
        local beatVal = beat - math.floor(beat);
        ypos=ypos+ (beatVal*(stepCrochet*m.jump))*scrollSpeed;
    end
   if m.sine ~= 0 then
	local len=4
		ypos = ypos+len * math.sin(len / 100.0 + beat * 2.0 * math.pi) * m.sine
	end
	if m.distortwiggle ~= 0 then
	local y=getPropertyFromGroup('strumLineNotes',(pn==2 and iCol+4 or iCol),"y")+fYOffset
		xpos = xpos + (math.sin(y / (200+m.distortwiggleperiod) + beat + m.distortwigglescratch) * m.distortwiggle * 20)
	end
	if m.wiggle ~= 0 then
		xpos = xpos + math.sin(beat) * m.wiggle * 20;
		ypos = ypos+ math.sin(beat + 1) *m.wiggle * 20;
		if m.wiggle > 0 then
		m.rotatez = math.sin(beat)*0.2*m.wiggle
		end
    end

   xpos =xpos+ addModchartingToolsmod(fYOffset, iCol, pn).x
   ypos =ypos+ addModchartingToolsmod(fYOffset, iCol, pn).y
   zpos =zpos+ addModchartingToolsmod(fYOffset, iCol, pn).z
   rotz = rotz+ addModchartingToolsmod(fYOffset, iCol, pn).angle
    return xpos, ypos, rotz, zpos
end
--



defaultPositions = {{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0}}
defaultscale = {{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0},{x=0,y=0}}


-- events

mods,curmod = {},1
perframe = {}
event,curevent = {},1
songStarted = false


function set(t)
	table.insert(mods,{t[1],0,linear,t[2],t[3],pn=t.pn})
end
function me(t)
	table.insert(mods,t)
end
function mpf(t)
	table.insert(perframe,t)
end
function m2(t)
	table.insert(event,t)
end
flx=0
function TEMPLATE.songStart()
    
    downscroll = false

	for i=0,7 do
        defaultPositions[i+1].x = getPropertyFromGroup("strumLineNotes",i,"x")
        defaultPositions[i+1].y = getPropertyFromGroup("strumLineNotes",i,"y")
        defaultscale[i+1].x = getPropertyFromGroup("strumLineNotes",i,"scale.x")
        defaultscale[i+1].y = getPropertyFromGroup("strumLineNotes",i,"scale.y")
        --print(i .. ": " .. defaultPositions[i+1].x .. " " .. defaultPositions[i+1].y)
    end
	
	--fuck it, it's mods. You don't get a say here.
	--(this is done to prevent a lot of bugs and weird cases)
	storedScrollSpeed = 1.8
	--storedScrollSpeed = scrollSpeed
	
	for pn=1,2 do
		activeMods[pn].xmod = storedScrollSpeed
	end
	
	songStarted = true
	
end
local s = 0
function TEMPLATE.update(elapsed)
s = getSongPosition()
    beat = ( s/ 1000) * (curBpm/60)
	
	while curmod <= #mods and beat > mods[curmod][1] do
		local v = mods[curmod]
		
		local mn = v[5]
		local dur = v[2]
		if v.timing and v.timing == 'end' then
			dur = v[2]-v[1]
		end
		
		
		if v.plr and not v.pn then v.pn = v.plr end
		
		for pn=1,2 do
			if not v.pn or pn == v.pn then
				tweenStart[pn][mn] = v[1]
				tweenLen[pn][mn] = dur
				tweenCurve[pn][mn] = v[3]
				if v.startVal then
					storedMods[pn][mn] = v.startVal
				else
					storedMods[pn][mn] = activeMods[pn][mn]
				end
				targetMods[pn][mn] = v[4]
				tweenEx1[pn][mn] = v.ex1
				tweenEx2[pn][mn] = v.ex2
				isTweening[pn][mn] = true
			end
		end
		curmod = curmod+1
	end
	
	for pn=1,2 do
		for _,v in pairs(modnames) do
			
			if isTweening[pn][v] then
				local curtime = beat - tweenStart[pn][v]
				local duration = tweenLen[pn][v]
				local startstrength = storedMods[pn][v]
				local diff = targetMods[pn][v] - startstrength
				local curve = tweenCurve[pn][v]
				local strength = curve(curtime, startstrength, diff, duration, tweenEx1[pn][v], tweenEx2[pn][v])
				activeMods[pn][v] = strength
				if beat > tweenStart[pn][v]+duration then
					isTweening[pn][v] = false
				end
			else
				activeMods[pn][v] = targetMods[pn][v]
			end
			
		end
	end
	
	----------------------------------------
	-- do this stuff every frame --
	----------------------------------------
	if #perframe>0 then
		for i=1,#perframe do
			local a = perframe[i]
			if beat > a[1] and beat < a[2] then
				a[3](beat)
			end
		end
	end
	
	-----------------------------------------
	-- event queue --event,curevent = {},1 --
	-----------------------------------------
	while curevent <= #event and beat>=event[curevent][1] do
		if event[curevent][3] or beat < event[curevent][1]+2 then
			event[curevent][2]()
		end
		curevent = curevent+1;
	end
	
	---------------------------------------
	-- ACTUALLY APPLY THE RESULTS OF THE ABOVE CALCULATIONS TO THE NOTES
	---------------------------------------
runHaxeCode([[
var c0 = getVar("c0");
var c1 = getVar("c1");
var cs = getVar("cs");
var cs1 = getVar("cs1");
for (p in game.notes){
if (p.mustPress){
p.cameras = cs1;
}else{//fuck修复了
p.cameras = cs;}}
game.strumLineNotes.members[0].cameras = cs;
game.strumLineNotes.members[1].cameras = cs;
game.strumLineNotes.members[2].cameras = cs;
game.strumLineNotes.members[3].cameras = cs;
game.strumLineNotes.members[4].cameras = cs1;
game.strumLineNotes.members[5].cameras = cs1;
game.strumLineNotes.members[6].cameras = cs1;
game.strumLineNotes.members[7].cameras = cs1;
for (i in 0...]]..tostring(CAMERA_COUNT)..[[){
var v = ]]..tostring(CAMERA_COUNT)..[[;
var ag = ]]..tostring(activeMods[1].rotationz+activeMods[1].camwag * math.sin(beat*pi))..[[;
var ag2 = ]]..tostring(activeMods[2].rotationz+activeMods[2].camwag *math.sin(beat*pi))..[[;
var cx = ]]..tostring(activeMods[1].camx)..[[;
var cy = ]]..tostring(activeMods[1].camy)..[[;
var cx2 = ]]..tostring(activeMods[2].camx)..[[;
var cy2 = ]]..tostring(activeMods[2].camy)..[[;
c0.x = cx;c0.y = cy;c0.angle = ag;
c1.x = cx2;c1.y = cy2;c1.angle = ag2;
c0.alpha = ]]..tostring(activeMods[1].camalpha)..[[;
c1.alpha = ]]..tostring(activeMods[2].camalpha)..[[;
c0.zoom+=]]..tostring(activeMods[1].camzoom)..[[;
c1.zoom+=]]..tostring(activeMods[2].camzoom)..[[;
}]])--fuck

	if songStarted then
		for pn=1,2 do
		local m = activeMods[pn]

			local xmod = m.xmod
			for col=0,3 do
			
				local c = (pn-1)*4 + col
				local xp, yp, rz, zp = arrowEffects(0, col, pn)
				local alp = receptorAlpha(col,pn)
				local defaultx, defaulty = defaultPositions[c+1].x,
				defaultPositions[c+1].y
                if m.MoveYWaveShit ~= 0 then
                yp =yp+ 260*math.sin((s*0.0008)+(col/4))*m.MoveYWaveShit
                end
                if m.xzfollowstrums ~= 0 then
                v1 = m.xzfollowstrums;
                moveSpeed = m.xzfollowstrumsspeed
                xp =xp+ math.sin(s*0.001*3*moveSpeed)*320*v1;
                zp =zp+ math.cos(s*0.001*3*moveSpeed)*320*v1;
                end
                if m.receptorscroll~= 0 then
                --sc = stepCrochet;
                --beat = m.receptorscroll
                --a = sc*beat;
                --if (math.floor(s/sc/beat) % 2 == 0) then
                --yp =45 + ((s%a)/(a/550));
                --else
                --yp =620 - ((s%a)/(a/550));
                --end
                end
				setPropertyFromGroup("strumLineNotes",c,"x",defaultx + xp)
                setPropertyFromGroup("strumLineNotes",c,"y",defaulty + yp)
                setPropertyFromGroup("strumLineNotes",c,"angle",rz)
                setPropertyFromGroup("strumLineNotes",c,"alpha",alp)
if m.arotatex ~= 0 or m.arotatey ~= 0 or m.arotatez ~= 0 or m["arotatex"..col] ~= 0 or m["arotatey"..col] ~= 0 or m["arotatez"..col] ~= 0 then
local origin = {defaultx, screenHeight/2}
local diff = {x=getPropertyFromGroup("strumLineNotes",c,"x")-origin[1],y=getPropertyFromGroup("strumLineNotes",c,"y")-origin[2],
    z=zp}
    local out = rotateV3(diff,math.rad(m.arotatex+m["arotatex"..col]),math.rad(m.arotatey+m["arotatey"..col]),math.rad(m.arotatez+m["arotatez"..col]))
    setPropertyFromGroup('strumLineNotes', c, 'x', origin[1]+out.x)
			setPropertyFromGroup('strumLineNotes', c, 'y', origin[2]+out.y)
zp=zp+out.z
end
if m.centerrotatex ~= 0 or m.centerrotatey ~= 0 or m.centerrotatez ~= 0 or m["centerrotatex"..col] ~= 0 or m["centerrotatey"..col] ~= 0 or m["centerrotatez"..col] ~= 0 then
	local origin = {screenWidth*0.5, screenHeight* 0.5}
    local diff = {x=getPropertyFromGroup("strumLineNotes",c,"x")-origin[1],y=getPropertyFromGroup("strumLineNotes",c,"y")-origin[2],
    z=zp}
    local out = rotateV3(diff,math.rad(m.centerrotatex+m["centerrotatex"..col]),math.rad(m.centerrotatey+m["centerrotatey"..col]),math.rad(m.centerrotatez+m["centerrotatez"..col]))
    setPropertyFromGroup('strumLineNotes', c, 'x', origin[1]+out.x)
			setPropertyFromGroup('strumLineNotes', c, 'y', origin[2]+out.y)
zp=zp+out.z
end
if m.localrotatex ~= 0 or m.localrotatey ~= 0 or m.localrotatez ~= 0 or m["localrotatex"..col] ~= 0 or m["localrotatey"..col] ~= 0 or m["localrotatez"..col] ~= 0 then
    local x = (screenWidth* 0.5) - ARROW_SIZE - 54 + ARROW_SIZE * 1.5;
    if pn == 2 then
    x = x + (screenWidth* 0.5 - ARROW_SIZE * 2 - 100)
    elseif pn == 1 then
    x = x - (screenWidth* 0.5 - ARROW_SIZE * 2 - 100)
    end
		x = x-56;

		local origin = {x, screenHeight* 0.5}
		local diff = {x=getPropertyFromGroup("strumLineNotes",c,"x")-origin[1],y=getPropertyFromGroup("strumLineNotes",c,"y")-origin[2],
    z=zp}
        diff.z = diff.z * screenHeight;
		local out = arotateV3(diff, math.rad(m.localrotatex + m['localrotatex'..col]),math.rad(m.localrotatey+m['localrotatey'..col]),math.rad(m.localrotatez + m['localrotatez'..col]));
        out.z = out.z / screenHeight;
    setPropertyFromGroup('strumLineNotes', c, 'x', origin[1]+out.x)
			setPropertyFromGroup('strumLineNotes', c, 'y', origin[2]+out.y)
zp=zp+out.z
end

--local rotx,roty,rotz = receptorRotation(0,col,pn)
--local rotation = rotateXYZ(rotx,roty,rotz)
--local anglepos={x=rotation.m00+rotation.m01+rotation.m02+rotation.m03,
--setPropertyFromGroup('strumLineNotes', c, 'x',getPropertyFromGroup('strumLineNotes', c, 'x')+anglepos.x)
--	setPropertyFromGroup('strumLineNotes', c, 'y', getPropertyFromGroup("strumLineNotes",c,"y")+anglepos.y)
--	zp=zp+anglepos.z
local scalex, scaley = getScale(0, col, pn, defaultscale[c+1].x, defaultscale[c+1].y)
local zNear,zFar = 0,100
	local zRange = zNear - zFar
	local fov = 90
	local tanHalfFOV = math.tan(math.rad(fov/2))
	local origin={x=getPropertyFromGroup("strumLineNotes",c,"x") - (screenWidth/2),y=getPropertyFromGroup("strumLineNotes",c,"y") - (screenHeight/2),z=zp}

    local r=rotateV3(origin,math.rad(m.fieldpitch),math.rad(m.fieldyaw),math.rad(m.fieldroll))
			--local pos={x=_G[strum..'X'..c%4]+xp - (screenWidth/2),y=_G[strum..'Y'..c%4]+yp - (screenHeight/2),z=zp/1000-1}
			local pos={x=r.x+m.fieldx,y=r.y+m.fieldy,z=(r.z-(1000+m.fieldz))/1000}
			local X = pos.x*(1/tanHalfFOV)/-pos.z+(screenWidth/2)
			local Y = pos.y/(1/tanHalfFOV)/-pos.z+(screenHeight/2)
			setPropertyFromGroup('strumLineNotes', c, 'x', X)
			setPropertyFromGroup('strumLineNotes', c, 'y', Y)
			local scale = -pos.z
			scale = 1 / scale
setPropertyFromGroup('strumLineNotes', c, 'scale.x', scalex * scale)--z lol
setPropertyFromGroup('strumLineNotes', c, 'scale.y', scaley * scale)			
 local zoom = getZoom(0,col,pn)
setPropertyFromGroup("strumLineNotes",c,"scale.x",getPropertyFromGroup("strumLineNotes",c,"scale.x")*zoom)
setPropertyFromGroup("strumLineNotes",c,"scale.y",getPropertyFromGroup("strumLineNotes",c,"scale.y")*zoom)
			end
		end
		end
		for v = 0, getProperty("notes.length")-1 do
			if getPropertyFromGroup('notes',v,"alive") then
				local pn = 1
				if getPropertyFromGroup('notes',v,"mustPress") then pn = 2 end
				local m =activeMods[pn]				
				local xmod = activeMods[pn].xmod				
				local isSus = getPropertyFromGroup('notes',v,"isSustainNote")
				local col = getPropertyFromGroup('notes',v,"noteData")
				local c = (pn-1)*4 + col				
				local targTime = getPropertyFromGroup('notes',v,"strumTime")
				local defaultx, defaulty = defaultPositions[c+1].x, defaultPositions[c+1].y
			   
				local scrollSpeeds = xmod * activeMods[pn]['xmod'..col] * (1 - 2*getReverseForCol(col,pn))*scrollSpeed

				local off = (1 - 2*getReverseForCol(col,pn))

				local ypos = getYAdjust(defaulty - (getSongPosition() - targTime),col,pn) * scrollSpeeds * 0.45 - off + ARROW_SIZE/2
				local zoom = getZoom(ypos-defaulty,col,pn)
				local xa, ya, rz, za = arrowEffects(ypos-defaulty, col, pn)
				local scalex, scaley = getScale(ypos-defaulty, col, pn, defaultscale[c+1].x, defaultscale[c+1].y)
				local alp = arrowAlpha(ypos-defaulty, col, pn)
				local fYOffset = ypos-defaulty
				if m.MoveYWaveShit ~= 0 then
				ya =ya+ 260*math.sin(((targTime+fYOffset)*0.0008)+(col/4))*m.MoveYWaveShit
				end
				if m.xzfollowstrums ~= 0 then
                v1 = m.xzfollowstrums;
                moveSpeed = m.xzfollowstrumsspeed
                time = targTime;
                xa =xa+ math.sin(time*0.001*3*moveSpeed)*320*v1;
                za =za+ math.cos(time*0.001*3*moveSpeed)*320*v1;
                end
                setPropertyFromGroup("notes",v,"x",defaultx + xa + (getPropertyFromGroup('notes',v,"isSustainNote") and 35.9 or 0))
                setPropertyFromGroup("notes",v,"y",ypos + ya-35.9*(1 - 2*getReverseForCol(col,pn)))
                if isSus then
                setPropertyFromGroup("notes",v,"alpha",alp*0.6)
                else
                setPropertyFromGroup("notes",v,"alpha",alp)
                end

if m.shrink ~= 0 then
scaleMult = 1 + ((ypos-defaulty)*0.001*m.shrink);
scalex =scalex* scaleMult;
scaley =scaley* scaleMult;
end
if rainbow == true
and version == '0.6.3' then
setPropertyFromGroup("notes",v,"colorSwap.hue",getPropertyFromGroup('notes',v,"colorSwap.hue")+0.01)
end
local susend = string.find(string.lower(getPropertyFromGroup('notes', v, 'animation.curAnim.name')), 'end') or string.find(string.lower(getPropertyFromGroup('notes', v, 'animation.curAnim.name')), 'tail')
				if getPropertyFromGroup('notes',v,"isSustainNote") then
				local ypos2 = getYAdjust(defaulty - ((getSongPosition()+.1) - targTime),col,pn) * scrollSpeeds * 0.45 - off --+ ARROW_SIZE/2
				--debugPrint(ypos2)
					local xa2, ya2 = arrowEffects(ypos2-defaulty, col, pn)
					setPropertyFromGroup("notes",v,"angle",math.deg(math.atan2(((ypos2 + ya2)-(ypos + ya))*100,(xa2-xa)*100) + math.pi/2))
				else
					setPropertyFromGroup("notes",v,"angle",rz)
                end
                if getPropertyFromGroup('notes',v,"scale.y") <=0 or getPropertyFromGroup('notes',v,"scale.x") <=0 then
                    setPropertyFromGroup("notes",v,"alpha",0)
                end

if m.arotatex ~= 0 or m.arotatey ~= 0 or m.arotatez ~= 0 or m["arotatex"..col] ~= 0 or m["arotatey"..col] ~= 0 or m["arotatez"..col] ~= 0 then
	local origin = {defaultx, screenHeight/2}
    local diff = {x=getPropertyFromGroup("notes",v,"x")-origin[1],y=getPropertyFromGroup("notes",v,"y")-origin[2],
    z=za}
    local out = rotateV3(diff,math.rad(m.arotatex+m["arotatex"..col]),math.rad(m.arotatey+m["arotatey"..col]),math.rad(m.arotatez+m["arotatez"..col]))
    setPropertyFromGroup('notes', v, 'x', origin[1]+out.x)
	setPropertyFromGroup('notes', v, 'y', origin[2]+out.y)
    za=za+out.z
end
if m.centerrotatex ~= 0 or m.centerrotatey ~= 0 or m.centerrotatez ~= 0 or m["centerrotatex"..col] ~= 0 or m["centerrotatey"..col] ~= 0 or m["centerrotatez"..col] ~= 0 then
	local origin = {screenWidth, screenHeight}
    local diff = {x=getPropertyFromGroup("notes",v,"x")-origin[1],y=getPropertyFromGroup("notes",v,"y")-origin[2],
    z=za}
    local out = rotateV3(diff,math.rad(m.centerrotatex+m["centerrotatex"..col]),math.rad(m.centerrotatey+m["centerrotatey"..col]),math.rad(m.centerrotatez+m["centerrotatez"..col]))
    setPropertyFromGroup('notes', v, 'x', origin[1]+out.x)
	setPropertyFromGroup('notes', v, 'y', origin[2]+out.y)
    za=za+out.z
end
if m.incomingAngleY ~= 0 then --incoming shit lol
    local diff = {x=0,y=getPropertyFromGroup("notes",v,"y")-getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y"),z=za}
    local out = rotateV3(diff,math.rad(m.incomingAngleY),0,0)
	setPropertyFromGroup('notes', v, 'y', getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y")+out.y)
    za=za+out.z
end
if m.incomingAngleX~= 0 then
    local diff = {x=getPropertyFromGroup("notes",v,"x")-getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"x"),y=getPropertyFromGroup("notes",v,"y")-getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y"),z=0}
    local out = rotateV3(diff,0,0,math.rad(m.incomingAngleX))
    setPropertyFromGroup('notes', v, 'x', getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"x")+out.x)
	setPropertyFromGroup('notes', v, 'y', getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y")+out.y)
end
if m.InvertIncomingAngle ~= 0 then
local diff = {x=0,y=getPropertyFromGroup("notes",v,"y")-getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y"),z=za}

if col %2==0 then
local out = rotateV3(diff,math.rad(m.incomingAngleY+m.InvertIncomingAngle+(fYOffset*0.015)),0,0)
	setPropertyFromGroup('notes', v, 'y', getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y")+out.y)
    za=za+out.z
else
local out = rotateV3(diff,math.rad(m.incomingAngleY-m.InvertIncomingAngle-(fYOffset*0.015)),0,0)
	setPropertyFromGroup('notes', v, 'y', getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y")+out.y)
    za=za+out.z
end
end
if m.IncomingAngleSmooth ~=0 then
local diff = {x=0,y=getPropertyFromGroup("notes",v,"y")-getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y"),z=za}
local out = rotateV3(diff,math.rad(m.incomingAngleY+m.IncomingAngleSmooth+(fYOffset*0.015)),0,0)
	setPropertyFromGroup('notes', v, 'y', getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y")+out.y)
    za=za+out.z
end
if m.IncomingAngleCurve ~=0 then
local diff = {x=0,y=getPropertyFromGroup("notes",v,"y")-getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y"),z=za}
local out = rotateV3(diff,math.rad(m.incomingAngleY+m.IncomingAngleCurve*(fYOffset*0.015)),0,0)
	setPropertyFromGroup('notes', v, 'y', getPropertyFromGroup('strumLineNotes',(pn==2 and col+4 or col),"y")+out.y)
    za=za+out.z
end
if m.localrotatex ~= 0 or m.localrotatey ~= 0 or m.localrotatez ~= 0 or m["localrotatex"..col] ~= 0 or m["localrotatey"..col] ~= 0 or m["localrotatez"..col] ~= 0 then
    local x = (screenWidth* 0.5) - ARROW_SIZE - 54 + ARROW_SIZE * 1.5;
    if pn == 2 then
    x = x + (screenWidth* 0.5 - ARROW_SIZE * 2 - 100)
    elseif pn == 1 then
    x = x - (screenWidth* 0.5 - ARROW_SIZE * 2 - 100)
    end
    x = x-56;
	local origin = {x, screenHeight* 0.5}
	local diff = {x=getPropertyFromGroup("notes",v,"x")-origin[1],y=getPropertyFromGroup("notes",v,"y")-origin[2],z=za}
    diff.z = diff.z * screenHeight
	local out = arotateV3(diff, math.rad(m.localrotatex + m['localrotatex'..col]),math.rad(m.localrotatey+m['localrotatey'..col]),math.rad(m.localrotatez + m['localrotatez'..col]));
    out.z = out.z / screenHeight;
    setPropertyFromGroup('notes', v, 'x', origin[1]+out.x)
	setPropertyFromGroup('notes', v, 'y', origin[2]+out.y)
za=za+out.z
end

            local xoff=(getPropertyFromGroup('notes',v,"isSustainNote") and 35 or 0)
            local zNear,zFar = 0,100
         	local zRange = zNear - zFar
        	local fov = 90
        	local tanHalfFOV = math.tan(math.rad(fov/2))
            local origin={x=getPropertyFromGroup('notes',v,"x") - (screenWidth/2),y= getPropertyFromGroup('notes',v,"y") - (screenHeight/2),z=za}
            local r=rotateV3(origin,math.rad(m.fieldpitch),math.rad(m.fieldyaw),math.rad(m.fieldroll))
			local pos={x=r.x+m.fieldx,y=r.y+m.fieldy,z=(r.z-(1000+m.fieldz))/1000}
			local X = pos.x*(1/tanHalfFOV)/-pos.z+(screenWidth/2)
			local Y = pos.y/(1/tanHalfFOV)/-pos.z+(screenHeight/2)
			setPropertyFromGroup('notes', v, 'x', X)
			setPropertyFromGroup('notes', v, 'y', Y)
			local scale = -pos.z
			scale = 1 / scale
			local scalenewy=getPropertyFromGroup('notes',v,"isSustainNote") and 1 or scaley
            setPropertyFromGroup('notes', v, 'scale.x', scalex * scale)
            setPropertyFromGroup('notes', v, 'scale.y', scalenewy * scale)
            setPropertyFromGroup("notes",v,"scale.x",(getPropertyFromGroup('notes',v,"scale.x")*zoom))
            if susend then
            if m.reverse ~=0 then
            setPropertyFromGroup("notes",v,"flipY",true)
            else setPropertyFromGroup("notes",v,"flipY",false)
            end
            end
            if isSus and not susend then
                setPropertyFromGroup("notes",v,"scale.y",scrollSpeed*m.xmod*(stepCrochet / 100 * 1.05))
                updateHitboxFromGroup("notes", i)
            else
                setPropertyFromGroup("notes",v,"scale.y",(getPropertyFromGroup('notes',v,"scale.y")* zoom))
            end
        end
   	end
end

function string:split( inSplitPattern, outResults ) -- from here, code isnt mine, https://stackoverflow.com/questions/19262761/lua-need-to-split-at-comma
    if not outResults then
      outResults = { }
    end
    local theStart = 1
    local theSplitStart, theSplitEnd = string.find( self, inSplitPattern, theStart )
    while theSplitStart do
      table.insert( outResults, string.sub( self, theStart, theSplitStart-1 ) )
      theStart = theSplitEnd + 1
      theSplitStart, theSplitEnd = string.find( self, inSplitPattern, theStart )
    end
    table.insert( outResults, string.sub( self, theStart ) )
    return outResults
end
function onCreatePost()
function print(...)end
addHaxeLibrary('Math')
luaDebugMode = true
runHaxeCode([[
var cams = [];
var cams1 = [];
//list千万不要加进for i里！！
for (i in 0...]]..tostring(CAMERA_COUNT)..[[){
var c0 = new FlxCamera();
c0.bgColor = 0x00;
FlxG.cameras.add(c0, false);
cams.push(c0);

var c1 = new FlxCamera();
c1.bgColor = 0x00;
FlxG.cameras.add(c1, false);
cams1.push(c1);

game.strumLineNotes.members[7].cameras = cams1;
c0.zoom = i*1/]]..tostring(CAMERA_COUNT)..[[;
c1.zoom = i*1/]]..tostring(CAMERA_COUNT)..[[;
//lol不要用onCreatePost改notes照相机
setVar("c0",c0);
setVar("c1",c1);
setVar("cs",cams);
setVar("cs1",cams1);
}
function btr(beat){return Math.round(beat * 48);}
]])
    function hide(t)
		local bt,tpn = t[1],t.pn
		for i=0,3 do
			me{bt+i*.125-1,.5,outExpo,-70,'movey'..i,pn=tpn}
			me{bt+i*.125-.5,1.25,inExpo,650,'movey'..i,pn=tpn}
			set{bt+i*.125+1.75,1,'stealth',pn=tpn}
			set{bt+i*.125+1.75,1,'dark',pn=tpn}
		end
	end
	function unhide(t)
		local bt,tpn = t[1],t.pn
		for i=0,3 do
			set{bt+i*.125-2,0,'stealth',pn=tpn}
			set{bt+i*.125-2,0,'dark',pn=tpn}
			me{bt+i*.125-2,1,outExpo,-70,'movey'..i,pn=tpn}
			me{bt+i*.125-1,1,inExpo,50,'movey'..i,pn=tpn}
			me{bt+i*.125-0,1.25,outElastic,0,'movey'..i,pn=tpn}
		end
	end
	--wiggle(beat,num,div,ease,amt,mod)
	function wig(t)
		local b,num,div,ea,am,mo = t[1],t[2],t[3],t[4],t[5],t[6]
		local f = 1
		for i=0,num do
			local smul = i==0 and 1 or 0
			local emul = i==num and 0 or 1
			
			me{b+i*(1/div),1/div,ea,startVal = am*smul*f, am*emul*-f,mo,pn=t.pn}
			
			f = f*-1
		end
	end
	--simple mod 2
    function sm2(tab)
		local b,len,eas,amt,mods,intime = tab[1],tab[2],tab[3],tab[4],tab[5],tab.intime
		if not intime then intime = .1 end
		if intime <= 0 then intime = .001 end
		me{b-intime,intime,linear,amt,mods,pn=tab.pn}
		me{b,len-intime,eas,0,mods,pn=tab.pn}
	end
	function halo(t)
	mpf{t[1],t[2],function(beat)
			for pn=1,2 do
			local h = t[3]
				for col=0,3 do
					local ang = 2*pi*((pn-1)*4 + col)/8 + ((beat-476)*0.5*pi)
					activeMods[pn]['movex'..col] = 280 * h * math.sin(ang)
					activeMods[pn]['movey'..col] = 70 * (1-(activeMods[pn].amovey/300)) * h * math.cos(ang)
				end
			end
		end}
	end
	function ease(beat,time,ease,str)	
	if string.find(str,"\n") then
	lol = stringSplit(str, '\n')
	end
    for i =0,#lol-1 do
    if i ~=0 then
    t = lol[i]:split(",")
	local val = t[1]
	local mod = t[2]
	me{beat,time,_G[ease],val,mod}
	end
    end
	end
	TEMPLATE.InitMods()	
	definemod{'colspacing',112}
	definemod{'spacing',620}
	definemod{'addx',0}
	definemod{'waveamp',0}
	init()--me{step,秒数,ease,数值,mod名字,pn=多少}
	TEMPLATE.setup()
end
function onSongStart()
TEMPLATE.songStart()
end
function onUpdatePost(elapsed)
TEMPLATE.update(elapsed)
end

function init()
me{1,1,linear,360,'arotatez'}

end
