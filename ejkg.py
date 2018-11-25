#
# VPython shell program for simulating Rutherford scattering of an alpha particle
# off of a gold nucleus.
#
from visual import *
from visual.graph import *
from random import *
#
G = 6.67e-11
earth_mass = 5.972e24
earth_radius = 6.371e6
rocket_maxthrust = 845e3 * 3 #3 active thrusters
specific_impulse = 290 
kerosene_max_mass = 123e3
LOX_max_mass = 287e3
empty_rocket_mass = 22.2e3
ballast_mass = 10e3
atm_height = 1200       #when "reentry" begins
initial_velocity = 0      #approx. Mach 6 on reentry
    
# Set up the display window.
#
scene = display(title = 'Self-Landing Rocket', width = 900, height = 1000, 
        range = 50, center = (0, atm_height, 0))
#
scene.autoscale = 0                             # Turn off auto scaling.
scene.autocenter = 0

#
# Set up rocket parameters
rocket = cylinder(pos=(0,0,0),axis=(0,42.6,0),radius=(20))
rocket.length = 42.6
rocket.radius = 3.7 / 2
rocket.fuel_pct = 0.1  #10% fuel remaining
rocket.mass = empty_rocket_mass + rocket.fuel_pct*(kerosene_max_mass + LOX_max_mass) + ballast_mass

rocket.theta = random() * pi/2     #angle from the vertical (y axis)
rocket.phi = random() * 2*pi       #rotation angle in xz plane

x = rocket.length * cos(rocket.phi) * sin(rocket.theta)
y = rocket.length * cos(rocket.theta)
z = rocket.length * sin(rocket.phi) * sin(rocket.theta)

rocket.axis = vector(x,y,z)

rocket.com = vector(0,atm_height, 0)

#calculating the distance from the bottom of the rocket to the center of mass
com_dist_from_bot =  (empty_rocket_mass * rocket.length / 2) + ((kerosene_max_mass + LOX_max_mass) * rocket.length * rocket.fuel_pct**2 / 2)
com_dist_from_bot /= rocket.mass

rocket.momentum = vector(0,-initial_velocity,0) * rocket.mass
rocket.ang_momentum = vector(0,0,0)

rocket.pos = rocket.com - rocket.axis * (com_dist_from_bot/rocket.length)
trail = curve(color = (0,0,1))
trail2 = curve(color = (0,0,0.5))

#rocket = cylinder(pos = rocket.pos, axis = rocket.axis, radius = rocket.radius, color = (1,0,0))

height_label = label(pos = rocket.com, height = 16, text = "height")

#gravity_arrow = arrow(pos = rocket.com, axis = (0,-30,0), shaftwidth = 1)
thrust_cone = cone(pos = rocket.pos, axis = -rocket.axis * 0.5, radius=rocket.radius, color=(1,0,0))
momentum_arrow = arrow(pos = rocket.pos + rocket.axis, axis = (0,0,0), shaftwidth = 1)

#landing_pad = cylinder(pos=(0,0,0), axis = (0,-1,0), radius = 25)

dt = 0.01

running = True

zoom_out_triggered = False
thrusting = False

count = 0

while(running):
    rate(100)   #simulation runs in real time

    fuel_mass = rocket.fuel_pct * (kerosene_max_mass + LOX_max_mass)
    com_fuel = rocket.pos + (rocket.fuel_pct * rocket.length * norm(rocket.axis) / 2)
    
    rocket.mass = fuel_mass + empty_rocket_mass + ballast_mass
    
    #calculating moment of inertia
    I_rocket = empty_rocket_mass * rocket.length**2 /12
    I_fuel = fuel_mass * (rocket.fuel_pct * rocket.length)**2 /12
    I_fuel += fuel_mass * mag2(rocket.com - com_fuel) #parallel axis theorem
    I_ballast = ballast_mass * (0.95 * rocket.length / 2)**2
    rocket.moment_I = I_rocket + I_fuel + I_ballast

    # Calculating effects of gravity
    force_gravity = G*(rocket.mass)*(earth_mass)/((rocket.com.y + earth_radius)**2) * vector(0, -1, 0)
    
    # Calculating effects of air resistance
    angle = acos(abs(dot(norm(rocket.momentum), vector(0,1,0)))) # angle between momentum and vertical
    frontal_area = (pi * rocket.radius**2)*abs(cos(angle)) + (rocket.length * 2*rocket.radius)*sin(angle)
    force_drag = .005 * frontal_area * 0.5 * 1.225 * mag2(rocket.momentum/rocket.mass) * vector(0, 1, 0)
    torque_drag = cross(rocket.com - rocket.pos, force_drag)

    # Adding minor random turbulence to make things interesting
    force_turbulence = vector(random()-0.5, random()-0.5, random()-0.5) * 5e4
    torque_turbulence = cross(rocket.com - rocket.pos, force_turbulence)

    #
    #   Figuring out when to activate engines
    #
    r = sqrt((rocket.com.x - rocket.x)**2 + (rocket.com.z - rocket.z)**2 + (rocket.com.y - rocket.y)**2)
    rocket.theta = acos((rocket.com.y-rocket.y) / r)
    thrust_theta = rocket.theta * 0.5
    if (abs(thrust_theta - rocket.theta) > pi/3): thrust_theta = rocket.theta - pi/3

    thrust_x = sin(thrust_theta) * cos(rocket.phi)
    thrust_y = cos(thrust_theta)
    thrust_z = sin(thrust_theta) * sin(rocket.phi)
    force_righting = vector(thrust_x, thrust_y, thrust_z)
    
    rocket.phi = atan2((rocket.com.z - rocket.z), (rocket.com.x - rocket.x))
    righting_factor = rocket.theta

    #restoring force to get rocket back on target
    r = sqrt(rocket.x**2 + rocket.y**2 + rocket.z**2)
    theta = acos(rocket.y / r)
    length = sqrt(rocket.x**2 + rocket.y**2) * cos(theta)
    length = vector(0, length, 0)
    force_targeting = norm(length - rocket.pos)
    targeting_factor = theta / (2*pi)

    #stabilizing force to act against angular momentum
    test = cylinder(pos = rocket.com, axis = norm(rocket.ang_momentum), radius = 0)
    test.rotate(angle = -pi/2, axis = rocket.axis, origin = rocket.com)
    force_stab = norm(-test.axis)
    stab_factor = (1/righting_factor)*(mag(rocket.ang_momentum)/rocket.moment_I / pi)

    #if(dot(force_stab, rocket.pos) < 0): force_stab = vector(0,0,0)

    #making sure rocket is not falling too quickly
    falling_factor = abs(rocket.momentum.y / rocket.mass)**2 / (rocket.y)
    force_falling = vector(0,1,0)

    norm_thrust = norm(righting_factor * force_righting + targeting_factor*norm(force_targeting) + stab_factor*force_stab + falling_factor*force_falling)
    
    force_thrust = rocket_maxthrust * norm_thrust

    if(thrusting == False): force_thrust *= 0

    thrust_cone.axis = -norm(force_thrust) * 20


    torque_thrust = cross(rocket.com - rocket.pos, force_thrust)

    #calculates whether rocket needs to start thrusting or not
    #if(mag(rocket.ang_momentum)/rocket.moment_I > 0.1 or rocket.theta > 1): thrusting = True

    y_vel = rocket.momentum.y / rocket.mass
    accel = (-9.81 + (rocket_maxthrust / rocket.mass) * abs(norm_thrust.y))
    delta_y = -(rocket.y)
    x_vel = rocket.momentum.x / rocket.mass
    if(y_vel**2 - 4*(accel/2)*(-delta_y) > 0):
        time = (-(y_vel) + sqrt(y_vel**2 - 4*(accel/2)*(-delta_y)))/(accel)
    else:
        time = 0

    if(rocket.momentum.y < 0):
        arrival_speed_squared = y_vel**2 + 2 * accel * delta_y
        if(arrival_speed_squared > -400):  thrusting = True
        else:
            #if(rocket.theta > 0.1): thrusting = True
            #else: thrusting = False
            thrusting = False
    if(sqrt(rocket.x**2 + rocket.z**2) > 100):
        if(abs(x_vel + rocket.x) > abs(rocket.x)):
            thrusting = True
        elif((x_vel < 0 and rocket.x > 0) and (x_vel * time)+(rocket.x + 100) > 0):
            thrusting = True
        elif((x_vel > 0 and rocket.x < 0)and (x_vel * time)+(rocket.x - 100) < 0):
            thrusting = True
        else:
            thrusting = False
    elif(rocket.momentum.y > 0):
        thrusting = False

    # Changing linear momentum and rotation
    rocket.momentum += (force_gravity + force_drag + force_turbulence + force_thrust) * dt
    rocket.com += (rocket.momentum / rocket.mass) * dt

    #gravity_arrow.pos = rocket.com

    # Changing angular momentum
    rocket.ang_momentum += (torque_drag + torque_turbulence + torque_thrust) * dt

    delta_theta = mag(rocket.ang_momentum) / rocket.moment_I * dt
    rocket.rotate(angle = delta_theta, axis = norm(rocket.ang_momentum), origin = rocket.com)

    #calculating new weight distribution of rocket due to loss of fuel
    com_dist_from_bot = (empty_rocket_mass * rocket.length / 2) + ((kerosene_max_mass + LOX_max_mass) * rocket.length * rocket.fuel_pct**2 / 2)
    com_dist_from_bot /= rocket.mass

    rocket.pos = rocket.com - rocket.axis * (com_dist_from_bot/rocket.length)
    trail.append(rocket.pos)
    trail2.append(rocket.pos + rocket.axis)
    thrust_cone.pos = rocket.pos
    #momentum_arrow.axis = rocket.momentum / rocket.mass / 2
    #momentum_arrow.pos = rocket.pos + rocket.axis

    #zooms out when rocket is about to land so that the landing can be viewed better
    if(rocket.pos.y <= 1000 and not zoom_out_triggered):
        zoom_out_triggered = True
        scene.center = (0,500,0)
        scene.range = 600
        scene.forward = vector(0,0,-1)
        landing_pad = cylinder(pos = (0,0,0), axis = (0,-2,0), radius = 100, color = (0.5,0.5,0.5))
    elif((rocket.pos.y >= 1000 or abs(rocket.pos.x) > 500) and zoom_out_triggered):
        zoom_out_triggered = False

    #stops once rocket has landed (or crashed)
    if(rocket.pos.y <= 1):
        running = False
        
    if(not zoom_out_triggered): scene.center = rocket.com
    height_label.pos = rocket.com + vector(scene.width * 0.4, scene.height * 0.4, 0)
    height_label.text = "speed: %s m/s height %s" %(int(mag(rocket.momentum/rocket.mass)), int(rocket.com.y))

    count += 1

    if(count % 100 == 0):
        cyl = cylinder(pos = rocket.pos, axis = rocket.axis, radius = rocket.radius, color = (.5,.5,.5))

angle = sqrt(rocket.axis.x **2 + rocket.axis.z ** 2)

x = 1
while(True):
    x += 1

