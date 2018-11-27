#
# Simulation of a self landing rocket
# Jason Chadwick, Thomas Garrison
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
atm_height = 5000       #when "reentry" begins
initial_velocity = 0      #approx. Mach 6 on reentry
limiting_angle = pi/100
    
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

rocket.theta = random() * pi/4     #angle from the vertical (y axis)
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
info_label = label(pos = rocket.com, height = 16, text = "placeholder")

#gravity_arrow = arrow(pos = rocket.com, axis = (0,-30,0), shaftwidth = 1)
thrust_cone = cone(pos = rocket.pos, axis = -rocket.axis * 0.5, radius=rocket.radius, color=(1,0,0))
momentum_arrow = arrow(pos = rocket.pos + rocket.axis, axis = (0,0,0), shaftwidth = 1)

#landing_pad = cylinder(pos=(0,0,0), axis = (0,-1,0), radius = 25)

dt = 0.01

running = True

zoom_out_triggered = False
thrusting = False

count = 0

righting_factor = 0
stab_factor = 0
targeting_factor = 0
falling_factor = 0

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
    angle = acos(abs(dot(norm(rocket.momentum), norm(rocket.axis)))) # angle between momentum and vertical
    frontal_area = (pi * rocket.radius**2)*abs(cos(angle)) + (rocket.length * 2*rocket.radius)*sin(angle)
    force_drag = .5 * frontal_area * 0.5 * 1.225 * mag2(rocket.momentum/rocket.mass) * -norm(rocket.momentum)
    force_drag *= 0
    torque_drag = cross((rocket.pos + 0.5*rocket.axis) - rocket.com, force_drag)

    # Adding minor random turbulence to make things interesting
    force_turbulence = vector(random()-0.5, random()-0.5, random()-0.5) * 0
    torque_turbulence = cross(rocket.com - rocket.pos, force_turbulence)

    #
    #   Calculating amount of thrust in each direction
    #

    #thrust to correct rocket angle
    r = sqrt((rocket.com.x - rocket.x)**2 + (rocket.com.z - rocket.z)**2 + (rocket.com.y - rocket.y)**2)
    rocket.theta = acos((rocket.com.y-rocket.y) / r)
    rocket.phi = atan2((rocket.com.z - rocket.z), (rocket.com.x - rocket.x))

    thrust_theta = rocket.theta * 0.1

    thrust_x = sin(thrust_theta) * cos(rocket.phi)
    thrust_y = cos(thrust_theta)
    thrust_z = sin(thrust_theta) * sin(rocket.phi)
    force_righting = vector(thrust_x, thrust_y, thrust_z)
    
    
    #restoring thrust to get rocket back towards target
    r = sqrt(rocket.x**2 + rocket.y**2 + rocket.z**2)
    theta = abs(acos(rocket.y / r))
    length = mag(rocket.pos) / cos(theta)
    length = vector(0, length, 0)
    force_targeting = norm(length - rocket.pos)
    
    
    #stabilizing thrust to slow rotation (acting against angular momentum)
    test = cylinder(pos = rocket.com, axis = norm(rocket.ang_momentum), radius = 0)
    test.rotate(angle = -pi/2, axis = rocket.axis, origin = rocket.com)
    force_stab = norm(-test.axis) * 0
    

    #thrust to keep rocket from crashing into the ground
    force_falling = vector(0,1,0)


    #calculates whether rocket needs to start thrusting or not
    y_vel = rocket.momentum.y / rocket.mass
    accel = (-9.81 + (rocket_maxthrust / rocket.mass) * norm(rocket.axis).y)
    delta_y = -(rocket.com.y - com_dist_from_bot) + 1
    horizontal_pos = vector(rocket.x, 0, rocket.z)
    horizontal_momentum = vector(rocket.momentum.x, 0, rocket.momentum.z)
    horizontal_dist = mag(horizontal_pos)
    horizontal_vel = mag(horizontal_momentum) / rocket.mass

    #factors which decide how much each component of the thrust is weighted
    #generally between 0 and 1, but sometimes greater if it is very important
    righting_factor = rocket.theta
    targeting_factor = theta / (pi)
    #if(dot(horizontal_momentum, horizontal_pos) > 0):
    #    targeting_factor *= mag(horizontal_momentum)
    #stab_factor = (1/(righting_factor+0.1))*(mag(rocket.ang_momentum)/rocket.moment_I / pi)
    stab_factor = (1/righting_factor)*(mag(rocket.ang_momentum)/rocket.moment_I / pi)
    falling_factor = abs(rocket.momentum.y / rocket.mass)**2 / (rocket.y)

    
    #calculating time to arrival
    if(y_vel**2 - 4*(accel/2)*(-delta_y) > 0):
        time = (-(y_vel) - sqrt(y_vel**2 - 4*(accel/2)*(-delta_y)))/(accel)
    else:
        time = 1000
    #if(count % 1 < 0.01): print time

    thrusting = False
    info_label.text = ""

    if(rocket.y < 1):
        thrusting = False
    elif(rocket.momentum.y < 0):
        arrival_speed_squared = y_vel**2 + 2 * accel * delta_y
        if(arrival_speed_squared > -(rocket.y/2)):
            falling_factor = 1000
            thrusting = True
            info_label.text = "Correcting for falling"
        else:
            if(rocket.theta > 0.1):
                righting_factor = 50
                thrusting = True
                info_label.text = "Correcting for angle"
            elif(horizontal_dist > 10):
                if(dot(horizontal_momentum, horizontal_pos) > 0) or (abs(horizontal_vel * time) - horizontal_dist > 10):   #if momentum is away from (0,0,0)
                    #targeting_factor *= 10
                    thrusting = True
                    info_label.text = "Correcting for displacement"
            else: thrusting = False
    
    
    #targeting_factor *= 0

    ang_velocity = mag(rocket.ang_momentum) / rocket.moment_I
    thrust_cone.color = (1,0,0)
    '''
    #if rocket is currently righting itself, slows it down if it will overshoot
    if thrusting and (dot(rocket.ang_momentum, cross(rocket.axis, vector(0,1,0))) > 0) and ang_velocity > 0.05:
        thrust_cone.color = (0,1,0)
        #print("ahhh")
        max_torque = -rocket_maxthrust * sin(limiting_angle) * com_dist_from_bot
        if(ang_velocity**2 - 4*(max_torque/2)*(rocket.theta) > 0):
            time = (-(ang_velocity) + sqrt((ang_velocity)**2 - 4*(max_torque/2)*(rocket.theta)))/(max_torque)
        else:
            time = 1000
        #print(time)
        if(-max_torque * abs(time) < mag(rocket.ang_momentum)):
            righting_phi = atan2(force_righting.z, force_righting.x)
            righting_phi += pi
            force_righting.x = sin(limiting_angle) * cos(righting_phi)
            force_righting.y = cos(limiting_angle)
            force_righting.z = sin(limiting_angle) * sin(righting_phi)
            
            force_righting = -(rocket.axis).rotate(angle = limiting_angle, axis = -rocket.ang_momentum)
            force_righting = norm(force_righting)
            thrust_cone.color = (0,0,1)
            #righting_factor = 50
            #temp = force_righting.x
            #force_righting.x = -force_righting.z
            #force_righting.z = -temp
    '''
    
    norm_thrust = norm(righting_factor * force_righting + targeting_factor*force_targeting + stab_factor*force_stab + falling_factor*force_falling)
    #norm_thrust = vector(1,0,0)
    
    if(dot(norm_thrust, rocket.axis) < cos(limiting_angle)):
        thrust_phi = atan2(norm_thrust.z - rocket.axis.z, norm_thrust.x - rocket.axis.x)

        axis = cross(norm_thrust, rocket.axis)
        

        theta = rocket.theta 
        
        thrust_x = sin(limiting_angle) * cos(thrust_phi)
        thrust_y = cos(limiting_angle)
        thrust_z = sin(limiting_angle) * sin(thrust_phi)

        norm_thrust = vector(thrust_x, thrust_y, thrust_z)

        norm_thrust = (-rocket.axis).rotate(angle = limiting_angle, axis = axis)
    
    
    force_thrust = rocket_maxthrust * norm_thrust
    if(thrusting == False): force_thrust *= 0

    thrust_cone.axis = -norm(force_thrust) * 20 #visible cone showing thrust

    torque_thrust = cross(rocket.com - rocket.pos, force_thrust)
    

    # Changing linear momentum and rotation
    rocket.momentum += (force_gravity + force_drag + force_turbulence + force_thrust) * dt
    rocket.com += (rocket.momentum / rocket.mass) * dt

    #gravity_arrow.pos = rocket.com

    '''print("new")
    print(torque_drag)
    print(torque_turbulence)
    print(torque_thrust)'''
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

    #zooms out when rocket is about to land so that the landing can be viewed better
    if(rocket.pos.y <= 1000 and not zoom_out_triggered):
        zoom_out_triggered = True
        scene.center = (0,500,0)
        scene.range = 600
        scene.forward = vector(0,0,-1)
        landing_pad = cylinder(pos = (0,0,0), axis = (0,-2,0), radius = 100, color = (0.5,0.5,0.5))
    elif((rocket.pos.y >= 1000 or abs(rocket.pos.x) > 500) and zoom_out_triggered):
        zoom_out_triggered = False



    #centers screen on rocket
    if(not zoom_out_triggered): scene.center = rocket.com

    #moves height label to be right near rocket
    if(zoom_out_triggered):
        height_label.pos = rocket.com + vector(300,0,0)
        info_label.pos = rocket.com + vector(300, -50, 0)
    else:
        height_label.pos = rocket.com + vector(30,0,0)
        info_label.pos = rocket.com + vector(30, -5, 0)
    height_label.text = "speed: %s m/s height %s" %(int(mag(rocket.momentum/rocket.mass)), int(rocket.y))

    #creates intermediate cylinders to view progress
    count += dt
    #if(count % 100 == 0):
        #cyl = cylinder(pos = rocket.pos, axis = rocket.axis, radius = rocket.radius, color = (.5,.5,.5))

    #stops once rocket has landed (or crashed)
    if(rocket.y < 1):
        running = False

angle = sqrt(rocket.axis.x **2 + rocket.axis.z ** 2)
speed = mag(rocket.momentum) / rocket.mass
print("angle: %s") %(angle)
print("speed: %s") %(speed)
print(rocket.y)
