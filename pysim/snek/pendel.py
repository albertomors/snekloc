### PENDULUM EQUATION ###
#
# I = ml^2
# theta_ddot = -b/I*theta_dot -g/l*sin(theta) + u/I
# theta = angle from vertical > 0 ccw
#
#########################

import pygame
import sys
import numpy as np
import math
from buttons import StandardSlider, SquareButton, StageButton
from buttons import slider_size, btn_size
from plotter import PlotFigure

#===============================================================
# SCREEN
#===============================================================

SCREENW, SCREENH = 800, 600
pygame.init()
screen = pygame.display.set_mode((SCREENW, SCREENH))
fg = pygame.Surface((SCREENW, SCREENH), pygame.SRCALPHA)
db_screen = pygame.Surface((SCREENW, SCREENH), pygame.SRCALPHA)
pygame.display.set_caption("pendel")
clock = pygame.time.Clock()

small_font = pygame.font.SysFont('bookantiqua', 8)
font = pygame.font.SysFont('bookantiqua', 18)

SCREEN_SCALE = 300 #px = 1 meter
screen_scale = SCREEN_SCALE

#===============================================================
# COLORS
#===============================================================

DARK = np.array([30, 30, 30])
WHITE = np.array([255, 255, 255])
RED = np.array([255, 0, 0])
GREEN = np.array([0, 255, 0])

#===============================================================
# DRAWING METHODS
#===============================================================

"""convert xy axis coordinates to screen coordinates tuple (,2)"""
def xy2screen(x,y): 
    sx = int(x * screen_scale)
    sy = int(SCREENH - y * screen_scale)
    return sx, sy

"""convert screen coordinates to xy axis coordinates np.array[(2,1)"""
def screen2xy(sx, sy): 
    x = sx / screen_scale
    y = (SCREENH - sy) / screen_scale
    return np.array([x,y]).reshape(2,1)

"""flip sign of angle from screen to xy coordinates: ccw = >0 and viceversa"""
def a2a(angle):
    return -angle

#===============================================================
# SNAKE
#===============================================================

# pendel params

m = 1         #mass
l = 0.8         #length
g = 9.81        #gravity
I = m*l**2      #moment of inertia
b = 0.3         #viscous friction
r = 0.05        #radius of a disc

# sim params
dt = 0.01

# derived params
nx = 2
nu = 1

# initial conditions
# theta_0 = np.random.uniform(-np.pi/4, np.pi/4, 1)
theta_0 = 0
theta_dot_0 = 0
x0 = np.array([theta_0, theta_dot_0]).reshape(-1,1)
x = x0.copy()

#===============================================================
# REFERENCE
#===============================================================

theta_ref = np.pi/4

""" smallest signed angle function"""
def ssa(angle):
    return np.arctan2(np.sin(angle), np.cos(angle))

#===============================================================
# DYNAMICS
#===============================================================

# dynamics using (2.33) model
def dyn(x,u):

    theta = x[0]
    theta_dot = x[1]
    theta_ddot = -b/I*theta_dot -g/l*math.sin(theta) + u/I

    dx = np.array([theta_dot, theta_ddot]).reshape(-1,1)

    return dx

#===============================================================
# INTEGRATION METHODS
#===============================================================

def RK4(x,u,law,dt):
    k1 = law(x, u)
    k2 = law(x+dt/2*k1, u)
    k3 = law(x+dt/2*k2, u)
    k4 = law(x+dt*k3, u)

    x_next = x + dt/6*(k1+2*k2+2*k3+k4)
    return x_next

def euler_step(x,u,law,dt):
    return x + dt*law(x,u)

#===============================================================
# REFERENCE FILTERS
#===============================================================

dtheta_ref_filt, theta_ref_filt, theta_ref_dot, theta_ref_ddot = 0,0,0,0

xi = 1
deltaT = 0.1

def third_order_filter(x, r, xi, deltaT):
    w = 2*np.pi/deltaT
    F = np.array([[0, 1, 0], [0, 0, 1], [-w**3, -(2*xi+1)*(w**2), -(2*xi+1)*w]])
    tf = lambda x,r: F@x + np.array([0, 0, w**3]).reshape(3,1)*r
    z = euler_step(x, r, tf, dt)
    x=z[0]; xdot=z[1]; xddot=z[2]
    
    return x, xdot, xddot

#===============================================================
# BODY SHAPE CONTROLLER
#===============================================================

MAX_TORQUE = 2 #Nm

"""control saturation to [-treshold, treshold]"""
def sat(u, treshold):
    return np.clip(u, -treshold, treshold)

"""simple pd controller"""
class pd():
    def __init__(self, kp=16, kd=3.2):
        self.kp = kp
        self.kd = kd

    def u(self, x, theta_ref, theta_ref_dot, theta_ref_ddot):
        
        theta = x[0]
        theta_dot = x[1]

        u_bar = theta_ref_ddot - self.kd * (theta_dot - theta_ref_dot) - self.kp * ssa(theta - theta_ref)
        u = u_bar

        return sat(u, MAX_TORQUE)
        # return u

"""fblin controller"""
class fblin():
    def __init__(self, kp=16, kd=3.2):
        self.kp = kp
        self.kd = kd

    def u(self, x, theta_ref, theta_ref_dot, theta_ref_ddot):
        
        theta = x[0]
        theta_dot = x[1]

        u_bar = theta_ref_ddot - self.kd * (theta_dot - theta_ref_dot) - self.kp * ssa(theta - theta_ref)
        # theta_ddot = -b/I*theta_dot -g/l*sin(theta) + u/I
        u = I*(u_bar + g/l*math.sin(theta) + b/I*theta_dot)

        return sat(u, MAX_TORQUE)
        # return u

"""swingup controller"""
class swingup():
    def __init__(self, kp=16):
        self.kp = kp

    def u(self, x, theta_ref, theta_ref_dot, theta_ref_ddot):
        
        theta = x[0]
        theta_dot = x[1]

        #E = 0.5*I*theta_dot**2 - m*g*l*(1 - math.cos(theta))
        E_tilde = 0.5*I*theta_dot**2 - m*g*l*(1 + math.cos(theta)) #valid only for upright equilibrium
        u = -self.kp*E_tilde*theta_dot + b*theta_dot

        return sat(u, MAX_TORQUE)
        # return u

class hybrid():
    def __init__(self, kp=16, kd=3.2, threshold=np.deg2rad(20)):
        self.kp = kp
        self.kd = kd
        self.fblin_ctrl = fblin()
        self.swingup_ctrl = swingup()
        self.threshold = threshold
        self.update_gains()
    
    def update_gains(self):
        self.fblin_ctrl.kp = self.kp
        self.fblin_ctrl.kd = self.kd
        self.swingup_ctrl.kp = self.kp

    def u(self, x, theta_ref, theta_ref_dot, theta_ref_ddot):
        
        theta = x[0]
        theta_dot = x[1]

        self.update_gains()

        if abs(ssa(theta_ref-np.pi)) < self.threshold and abs(ssa(theta-theta_ref)) > self.threshold:
            u = self.swingup_ctrl.u(x, theta_ref, theta_ref_dot, theta_ref_ddot)
        else:
            u = self.fblin_ctrl.u(x, theta_ref, theta_ref_dot, theta_ref_ddot)

        return sat(u, MAX_TORQUE)
        # return u

#===============================================================
# DRAWING
#===============================================================

def draw_pendel(x):

    theta = x[0]
    theta_dot = x[1]

    jx, jy = screen2xy(SCREENW//2, SCREENH//2)  #find center of the screen
    jy += l/2                                   #move it up by half pendel length
    bx = jx + l*math.sin(theta)
    by = jy - l*math.cos(theta)

    pygame.draw.line(fg, WHITE, xy2screen(jx, jy), xy2screen(bx, by), 1)
    pygame.draw.circle(fg, WHITE, xy2screen(bx, by), r*screen_scale)

#===============================================================
# MAIN
#===============================================================

sliders = {
    "kp": StandardSlider("kp", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*3,    0, 16, 50),
    "kd": StandardSlider("kd", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*2,    0, 3.2, 50),
    "scale": StandardSlider("scale", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*1, 1, SCREEN_SCALE, 1000),
}
buttons = {
    "free" : SquareButton("free", SCREENW-btn_size-10, SCREENH-(slider_size[1]+1)*5-btn_size-10),
    "debug" : SquareButton("debug", SCREENW-btn_size-10, SCREENH-(slider_size[1]+1)*4-btn_size-10),
}

for btn in buttons.values():
    btn.draw_once(font, screen)
for slider in sliders.values():
    slider.draw_once(small_font, screen)

ctrl = hybrid()

debug = False
running = True
free = False
sim_time = 0
t_offset = pygame.time.get_ticks()

theta_hist = []
theta_ref_hist = []
u_hist = []
timestamp = []
plotf_theta = PlotFigure(100, 100, font=font, axiscolor=WHITE, margin=5)
plotf_u = PlotFigure(100, 100, font=font, axiscolor=WHITE, margin=5)

while running:

    #===============================================================
    # EVENTS HANDLING
    #===============================================================

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

        for slider in sliders.values():
            slider.handle_event(event)
        for button in buttons.values():
            button.handle_event(event)
            button.update()

    real_time = (pygame.time.get_ticks() - t_offset) / 1000
    keys = pygame.key.get_pressed()
    mouse_pos = pygame.mouse.get_pos()

    # get values from sliders
    ctrl.kp = sliders["kp"].get_value()
    ctrl.kd = sliders["kd"].get_value()
    debug = buttons["debug"].get_value()
    free = buttons["free"].get_value()
    screen_scale = round(sliders["scale"].get_value())

    #===============================================================
    # UPDATE
    #===============================================================

    jx, jy = screen2xy(SCREENW//2, SCREENH//2)
    sjx, sjy = xy2screen(jx, jy+l/2)
    angle = math.atan2(mouse_pos[1] - sjy, mouse_pos[0] - sjx)

    db_screen.fill((0,0,0,0))
    pygame.draw.line(db_screen, GREEN, (sjx, sjy), mouse_pos, 1)

    #fix screen flipped coordinates and ssa for continuity rotating aroung 360° [instead of making a step from -180° to +180°]
    theta_ref = a2a(angle)+np.pi/2
    theta_ref = theta_ref_filt + ssa(theta_ref-theta_ref_filt)
    theta_ref_filt, theta_ref_dot, theta_ref_ddot = third_order_filter(np.array([theta_ref_filt, theta_ref_dot, theta_ref_ddot]).reshape(-1,1), theta_ref, xi, deltaT)

    u = ctrl.u(x, theta_ref, 0, 0)
    if free:
        u = 0
    x = RK4(x, u, dyn, dt)

    #===============================================================
    # PLOTS
    #===============================================================

    if debug:
        theta_hist.append(x[0]%(2*np.pi))
        if len(theta_hist) > 100:
            theta_hist.pop(0)
        theta_ref_hist.append(theta_ref_filt%(2*np.pi))
        if len(theta_ref_hist) > 100:
            theta_ref_hist.pop(0)
        u_hist.append(u)
        if len(u_hist) > 100:
            u_hist.pop(0)

        timestamp.append(sim_time)
        if len(timestamp) > 100:
            timestamp.pop(0)

        plotf_theta.plot(timestamp, theta_ref_hist, xlim=(timestamp[0], timestamp[-1]), ylim=(0,2*np.pi), line_width=1, color=(128,128,128))
        plotf_theta.hold_on()
        plotf_theta.plot(timestamp, theta_hist, xlim=(timestamp[0], timestamp[-1]), ylim=(0,2*np.pi), line_width=1, color=(255,255,255))
        plotf_theta.hold_off()

        plotf_u.plot(timestamp, u_hist, xlim=(timestamp[0], timestamp[-1]), ylim=(-MAX_TORQUE-1, MAX_TORQUE+1), line_width=1, color=(255,255,255))

        db_screen.blit(plotf_theta.surf, (10, SCREENH-150))
        db_screen.blit(plotf_u.surf, (10+100+10, SCREENH-150))

    #===============================================================
    # DRAWING
    #===============================================================

    screen.fill(DARK)
    fg.fill((0,0,0,0))

    draw_pendel(x)
    screen.blit(fg, (0,0))
    screen.blit(db_screen, (0,0)) if debug else None

    pygame.draw.line(screen, WHITE, (10,SCREENH-20), (10, SCREENH-10))
    pygame.draw.line(screen, WHITE, (10,SCREENH-10), (round(10+screen_scale), SCREENH-10))
    pygame.draw.line(screen, WHITE, (round(10+screen_scale), SCREENH-10), (round(10+screen_scale), SCREENH-20))
    screen.blit(font.render(f"1 m", True, WHITE), (round(10+screen_scale//2-15), SCREENH-30))

    screen.blit(font.render(f"sim_t: {sim_time:.2f} s / t: {real_time:.2f} s", True, WHITE), (10, 10))
    screen.blit(font.render(f"dt: {dt:.3f} s", True, WHITE), (10, 30))

    for slider in sliders.values():
        slider.draw(font,screen)
    for button in buttons.values():
        button.draw(font, screen)

    pygame.display.flip()
    sim_time += dt
    #sinchronize with real time if sim is running faster
    while sim_time > real_time:
        real_time = (pygame.time.get_ticks() - t_offset) / 1000

pygame.quit()
sys.exit()