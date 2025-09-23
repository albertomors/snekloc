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
bg = pygame.Surface((SCREENW, SCREENH), pygame.SRCALPHA)
lux_surf = pygame.Surface((SCREENW, SCREENH), pygame.SRCALPHA)
fg = pygame.Surface((SCREENW, SCREENH), pygame.SRCALPHA)
db_screen = pygame.Surface((SCREENW, SCREENH), pygame.SRCALPHA)
pygame.display.set_caption("snek simulator")
clock = pygame.time.Clock()

small_font = pygame.font.SysFont('bookantiqua', 8)
font = pygame.font.SysFont('bookantiqua', 18)

SCREEN_SCALE = 100 #px = 1 meter

#===============================================================
# COLORS
#===============================================================

DEEP_SEA = np.array([0, 0, 50])
WHITE = np.array([255, 255, 255])
RED = np.array([255, 0, 0])
GREEN = np.array([0, 255, 0])
MAGENTA = np.array([255, 0, 255])

def sin_color(color, t, scale):
    dc = (np.ones(len(color)) * np.sin(t) * scale)
    color = DEEP_SEA + dc.astype(int)
    return np.clip(color, 0, 255)

#===============================================================
# DRAWING METHODS
#===============================================================

"""convert xy axis coordinates to screen coordinates tuple (,2)"""
def xy2screen(x,y): 
    sx = int(x * SCREEN_SCALE)
    sy = int(SCREENH - y * SCREEN_SCALE)
    return sx, sy

"""convert screen coordinates to xy axis coordinates np.array[(2,1)"""
def screen2xy(sx, sy): 
    x = sx / SCREEN_SCALE
    y = (SCREENH - sy) / SCREEN_SCALE
    return np.array([x,y]).reshape(2,1)

"""flip sign of angle from screen to xy coordinates: ccw = >0 and viceversa"""
def a2a(angle):
    return -angle

#===============================================================
# SNAKE
#===============================================================

# snek params
Nl = 10
l = 0.09
m = 1.56
J = (1/3)*m*l**2

dt = 0.01

# env params
ct = 0.5
cn = 3

# derived params
nx = 2*Nl+4
nu = Nl-1
e = np.ones((Nl,1))
e_bar = np.ones((Nl-1,1))
A = np.diag(np.ones(Nl)) + np.diag(np.ones(Nl-1), 1)
A = A[:Nl-1, :]
D = np.diag(np.ones(Nl)) + np.diag(-np.ones(Nl-1), 1)
D = D[:Nl-1, :]
D_bar = D.T @ np.linalg.inv(D @ D.T)
K = A.T @ np.linalg.inv(D @ D.T) @ D
V = A.T @ np.linalg.inv(D @ D.T) @ A
H = np.triu(np.ones((Nl, Nl)))

# initial conditions
# theta_0 = np.random.uniform(-np.pi/4, np.pi/4, (Nl, 1))
theta_0 = np.zeros((Nl, 1))
theta_dot_0 = np.zeros((Nl, 1))
sp0 = screen2xy(SCREENW//2, SCREENH//2)
p0 = np.array([sp0[0], sp0[1]]).reshape(2, 1)
v0 = np.zeros((2, 1))
x0 = np.vstack((theta_0, p0, theta_dot_0, v0))
x = x0.copy()

#===============================================================
# REFERENCE GAIT
#===============================================================

alpha = np.deg2rad(13)
omega = np.deg2rad(120)
delta = np.deg2rad(50)

def lat_ond(t):

    N = np.arange(0, Nl-1).reshape(-1,1)
    phi_ref = alpha * np.sin(omega * t + N * delta)
    phi_ref_dot = omega * alpha * np.cos(omega * t + N * delta)
    phi_ref_ddot = -(omega**2) * alpha * np.sin(omega * t + N * delta)

    return phi_ref, phi_ref_dot, phi_ref_ddot

def eel_like(t):

    N = np.arange(0, Nl-1).reshape(-1,1)
    g = (Nl-N+1)/(Nl+1)
    phi_ref = g * (alpha * np.sin(omega * t + N * delta))
    phi_ref_dot = g * (omega * alpha * np.cos(omega * t + N * delta))
    phi_ref_ddot = g * (-(omega**2) * alpha * np.sin(omega * t + N * delta))

    return phi_ref, phi_ref_dot, phi_ref_ddot

#===============================================================
# DYNAMICS
#===============================================================

# dynamics using (2.33) model
def dyn233(x,u):

    theta = x[0:Nl]
    # px = x[Nl]
    # py = x[Nl+1]
    theta_dot = x[Nl+2:-2]
    px_dot = x[-2]
    py_dot = x[-1]

    #(2.12)
    St = np.diag(np.sin(theta).flatten())
    Ct = np.diag(np.cos(theta).flatten())
    Xd = l*K.T@St@theta_dot + e*px_dot
    Yd = -l*K.T@Ct@theta_dot + e*py_dot

    #(2.34)
    Mt = J*np.eye(Nl) + m*(l**2)*St@V@St + m*(l**2)*Ct@V@Ct
    Mtinv = np.linalg.inv(Mt)
    Wt = m*(l**2)*St@V@Ct - m*(l**2)*Ct@V@St

    #(2.25)
    Fx = -((ct*(Ct**2) + cn*(St**2))@Xd + (ct-cn)*(St@Ct)@Yd)
    Fy = -((ct-cn)*St@Ct@Xd + (ct*(St**2) + cn*(Ct**2))@Yd)

    #(2.33)
    dx = np.zeros((nx, 1))
    dx[0:Nl] = theta_dot
    dx[Nl:Nl+2] = [px_dot, py_dot]
    dx[Nl+2:-2] = Mtinv@(-Wt@(theta_dot**2) + l*St@K@Fx - l*Ct@K@Fy + D.T@u)
    dx[-2] = e.T @ Fx / (Nl*m)
    dx[-1] = e.T @ Fy / (Nl*m)

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
phi0_filt, phi0_filt_dot, phi0_filt_ddot = 0,0,0

xi = 1
deltaT = 4

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

MAX_TORQUE = 100

"""control saturation to [-treshold, treshold]"""
def sat(u, treshold):
    return np.clip(u, -treshold, treshold)

"""partial feedback linearization body shape controller"""
class fblin():
    def __init__(self, kp=16, kd=3.2):
        self.kp = kp
        self.kd = kd

    def u(self, x, phi_ref, phi_ref_dot, phi_ref_ddot):
        theta = x[0:Nl]
        # px = x[Nl]
        # py = x[Nl+1]
        theta_dot = x[Nl+2:-2]
        px_dot = x[-2]
        py_dot = x[-1]

        # theta->phi
        phi = D@theta
        phi_dot = D@theta_dot

        # state-dependent stuff
        St = np.diag(np.sin(theta).flatten())
        Ct = np.diag(np.cos(theta).flatten())
        Xd = l*K.T@St@theta_dot + e*px_dot
        Yd = -l*K.T@Ct@theta_dot + e*py_dot
        Mt = J*np.eye(Nl) + m*(l**2)*St@V@St + m*(l**2)*Ct@V@Ct
        # Mtinv = np.linalg.inv(Mt)
        Wt = m*(l**2)*St@V@Ct - m*(l**2)*Ct@V@St
        Fx = -((ct*(Ct**2) + cn*(St**2))@Xd + (ct-cn)*(St@Ct)@Yd)
        Fy = -((ct-cn)*St@Ct@Xd + (ct*(St**2) + cn*(Ct**2))@Yd)
        Fr = np.vstack([Fx,Fy])

        #assemble the big matrices.. (2.40)
        M_bar = np.block([[H.T@Mt@H, np.zeros((Nl,2))],
                        [np.zeros((2,Nl)), Nl*m*np.eye(2)]])
        W_bar = np.vstack([H.T@Wt@(theta_dot**2), np.zeros((2,1))])
        G_bar = np.block([[-l*H.T@St@K, l*H.T@Ct@K],
                        [-e.T, np.zeros((1,Nl))],
                        [np.zeros((1,Nl)), -e.T]])
        
        #..and split them (2.41)
        M11 = M_bar[0:Nl-1, 0:Nl-1]
        M12 = M_bar[0:Nl-1, Nl-1:]
        M21 = M_bar[Nl-1:, 0:Nl-1]
        M22 = M_bar[Nl-1:, Nl-1:]
        W1 = W_bar[0:Nl-1]
        W2 = W_bar[Nl-1:]
        G1 = G_bar[0:Nl-1,:]
        G2 = G_bar[Nl-1:,:]
        M22inv = np.linalg.inv(M22)

        # controller matrix terms
        linearizing_term = M11 - M12@M22inv@M21
        friction_term = G1 - M12@M22inv@G2
        velocity_term = W1 - M12@M22inv@W2

        # fb linearizing control law
        u_bar = phi_ref_ddot - self.kd * (phi_dot - phi_ref_dot) - self.kp * (phi - phi_ref)
        u = linearizing_term@u_bar + velocity_term + friction_term@Fr

        return sat(u, MAX_TORQUE)
        # return u

#===============================================================
# HEADING CONTROLLER
#===============================================================

MAX_PHI0 = np.pi/4
HALF_LENGTH_ON_SCREEN = Nl*2*l*SCREEN_SCALE/2

""" smallest signed angle function"""
def ssa(angle):
    return np.arctan2(np.sin(angle), np.cos(angle))

"""heading controller"""
class pd_heading():
    def __init__(self, kp=0.3, kd=0.1):
        self.kp = kp
        self.kd = kd

    def u(self, x, theta_ref, theta_ref_dot):

        heading = np.mean(x[0:Nl])
        heading_dot = np.mean(x[Nl+2:-2])
        # heading = x[Nl-1]
        # heading_dot = x[-3]

        px = x[Nl]
        py = x[Nl+1]
        spx, spy = xy2screen(px, py)
        sheading = a2a(heading)

        pygame.draw.line(db_screen, MAGENTA, (spx, spy), (spx + HALF_LENGTH_ON_SCREEN * np.cos(sheading), spy + HALF_LENGTH_ON_SCREEN * np.sin(sheading)), 1)
        theta_ref_filt_screen = a2a(theta_ref_filt)
        pygame.draw.line(db_screen, GREEN, (spx, spy), (spx + HALF_LENGTH_ON_SCREEN * np.cos(theta_ref_filt_screen).item(), spy + HALF_LENGTH_ON_SCREEN * np.sin(theta_ref_filt_screen).item()), 1)

        dtheta = heading - theta_ref
        dtheta = ssa(dtheta)
        phi0 = + self.kp*(dtheta) + self.kd*(heading_dot - theta_ref_dot)

        return sat(phi0, MAX_PHI0)

#===============================================================
# SNEK DRAWING
#===============================================================

LINKW = 0.01 #meters

def draw_snek(x):

    theta = x[0:Nl]
    px = x[Nl]
    py = x[Nl+1]
    # theta_dot = x[Nl+2:-2]
    # px_dot = x[-2]
    # py_dot = x[-1]

    #links CM
    X = -l*K.T@np.cos(theta) + e*px
    Y = -l*K.T@np.sin(theta) + e*py

    # compute start and end of every link
    pis = np.zeros((Nl+1, 2))
    pis[0, :] = np.array([X[0], Y[0]]).reshape(1,-1) - l*np.array([np.cos(theta[0]), np.sin(theta[0])]).reshape(1,-1)
    for i in range(1, Nl+1):
        pis[i, :] = pis[i-1, :] + 2*l*np.array([np.cos(theta[i-1]), np.sin(theta[i-1])]).reshape(1,-1)

    # draw px,py
    spx, spy = xy2screen(px, py)
    pygame.draw.rect(fg, WHITE, (spx-5, spy-5, 10, 10), 1)

    for i in range(Nl):
        # draw links
        x1, y1 = xy2screen(pis[i, 0], pis[i, 1])
        x2, y2 = xy2screen(pis[i+1, 0], pis[i+1, 1])
        pygame.draw.line(fg, WHITE, (x1, y1), (x2, y2), int(LINKW*SCREEN_SCALE))

        # draw CM
        pix, piy = xy2screen(X[i], Y[i])
        pygame.draw.circle(fg, RED, (pix, piy), 2)

        # draw lux on an addictional surf and overlaps it with lux_surf
        # addict_surf = pygame.Surface(lux_surf.get_size(), pygame.SRCALPHA)
        # pygame.draw.circle(addict_surf, (255,255,255,64), (pix, piy), 0.9*SCREEN_SCALE)
        # pygame.draw.circle(addict_surf, (255,255,255,128), (pix, piy), 0.3*SCREEN_SCALE)
        # pygame.draw.circle(addict_surf, (255,255,255,255), (pix, piy), 0.1*SCREEN_SCALE)
        # lux_surf.blit(addict_surf, (0, 0), special_flags=pygame.BLEND_RGBA_MAX)

    # draw frontal light
    a1 = theta[-1] - np.deg2rad(10)
    a2 = theta[-1] + np.deg2rad(10)
    sp1 = xy2screen(X[-1], Y[-1])
    sp2 = sp1 + max(SCREENW, SCREENH) * np.hstack([np.cos(a1), -np.sin(a1)])
    sp3 = sp1 + max(SCREENW, SCREENH) * np.hstack([np.cos(a2), -np.sin(a2)])
    pygame.draw.polygon(lux_surf, (255, 255, 255, 255), [sp1, sp2, sp3])

#===============================================================
# MAIN
#===============================================================

sliders = {
    "kp": StandardSlider("kp", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*5,    0, 16, 50),
    "kd": StandardSlider("kd", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*4,    0, 3.2, 50),
    "alpha": StandardSlider("α", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*3,  0, np.rad2deg(alpha), 90),
    "omega": StandardSlider("ω", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*2,  0, np.rad2deg(omega),  3*120),
    "delta": StandardSlider("δ", SCREENW-slider_size[0]-50, SCREENH-(slider_size[1]+1)*1,  0, np.rad2deg(delta), 90),
}
buttons = {
    "<<": StageButton("<<", 25, SCREENH//2, 50, 200),
    ">>": StageButton(">>", SCREENW-25, SCREENH//2, 50, 200),
    "debug" : SquareButton("debug", SCREENW-btn_size-10, SCREENH-(slider_size[1]+1)*7-btn_size-10),
    "gait" : SquareButton("gait", SCREENW-btn_size-10, SCREENH-(slider_size[1]+1)*6-btn_size, val1="lat", val2="eel")
}

for btn in buttons.values():
    btn.draw_once(font, screen)
for slider in sliders.values():
    slider.draw_once(small_font, screen)

body_ctrl = fblin()
heading_ctrl = pd_heading()
gait = lat_ond

debug = False
running = True
sim_time = 0
t_offset = pygame.time.get_ticks()

stage = 1

phi5_hist = []
phi5_ref_hist = []
timestamp = []
plotf = PlotFigure(100, 100, font=font, axiscolor=WHITE, margin=5)

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
            print(button.title, button.get_value())

    real_time = (pygame.time.get_ticks() - t_offset) / 1000
    keys = pygame.key.get_pressed()
    mouse_pos = pygame.mouse.get_pos()

    # get values from sliders
    body_ctrl.kp = sliders["kp"].get_value()
    body_ctrl.kd = sliders["kd"].get_value()
    alpha = np.deg2rad(sliders["alpha"].get_value())
    omega = np.deg2rad(sliders["omega"].get_value())
    delta = np.deg2rad(sliders["delta"].get_value())
    debug = buttons["debug"].get_value()
    gait = lat_ond if buttons["gait"].get_value() == "lat" else eel_like

    if buttons[">>"].get_value() == 1:
        stage += 1
        buttons[">>"].off()
    elif buttons["<<"].get_value() == 1:
        stage -= 1
        buttons["<<"].off()
    stage = max(1, min(3, stage))

    #===============================================================
    # UPDATE
    #===============================================================

    px = x[Nl]
    py = x[Nl+1]
    spx, spy = xy2screen(px, py)
    angle = math.atan2(mouse_pos[1] - spy, mouse_pos[0] - spx)

    if stage == 1:
        angle += np.mean(x[0:Nl-1])
    elif stage == 2:
        angle = 0

    db_screen.fill((0,0,0,0))
    pygame.draw.line(db_screen, GREEN, (spx, spy), mouse_pos, 1)

    #fix screen flipped coordinates and ssa for continuity rotating aroung 360° [instead of making a step from -180° to +180°]
    theta_ref = a2a(angle)
    theta_ref = theta_ref_filt + ssa(theta_ref-theta_ref_filt)
    theta_ref_filt, theta_ref_dot, theta_ref_ddot = third_order_filter(np.array([theta_ref_filt, theta_ref_dot, theta_ref_ddot]).reshape(-1,1), theta_ref, xi, deltaT)

    phi_ref, phi_ref_dot, phi_ref_ddot = gait(sim_time)
    if stage == 1:
        phi_ref = np.vstack([np.zeros((Nl-2,1)),-theta_ref])
        phi_ref_dot = np.zeros_like(phi_ref)
        phi_ref_ddot = np.zeros_like(phi_ref)

    if stage == 1:
        phi0 = 0
    else:
        phi0 = heading_ctrl.u(x, theta_ref_filt, theta_ref_dot)
    
    # if stage == 3:
    #     if x[Nl] > screen2xy(SCREENW/2, SCREENH/2)[0]:
    #         ct = 4
    #         cn = 4.5
    #     else:
    #         ct = 0.5
    #         cn = 3

    phi0_filt, phi0_filt_dot, phi0_filt_ddot = third_order_filter(np.array([phi0_filt, phi0_filt_dot, phi0_filt_ddot]).reshape(-1,1), phi0, xi, deltaT)
    
    u = body_ctrl.u(x, phi_ref + phi0_filt, phi_ref_dot + phi0_filt_dot, phi_ref_ddot + phi0_filt_ddot)
    x = RK4(x, u, dyn233, dt)

    if stage in (1, 2):
        x[Nl:Nl+2] = screen2xy(SCREENW//2, SCREENH//2) #lock it

    bg_color = sin_color(DEEP_SEA, real_time, 10)

    #===============================================================
    # PLOTS
    #===============================================================

    if debug:
        phis = D @ x[0:Nl]
        phi5_hist.append(phis[Nl-2])
        if len(phi5_hist) > 100:
            phi5_hist.pop(0)
        phi5_ref_hist.append(phi_ref[Nl-2])
        if len(phi5_ref_hist) > 100:
            phi5_ref_hist.pop(0)
        timestamp.append(sim_time)
        if len(timestamp) > 100:
            timestamp.pop(0)

        plotf.plot(timestamp, phi5_ref_hist, xlim=(timestamp[0], timestamp[-1]), ylim=(-np.pi/2,np.pi/2), line_width=1, color=(128,128,128))
        plotf.hold_on()
        plotf.plot(timestamp, phi5_hist, xlim=(timestamp[0], timestamp[-1]), ylim=(-np.pi/2,np.pi/2), line_width=1, color=(255,255,255))
        plotf.hold_off()
        # plotf._xticks(10)
        # plotf._yticks(13)
        # plotf._draw_ticks()
        # plotf.label('Time (s)', 'Value')
        # plotf._draw_labels()

        db_screen.blit(plotf.surf, (10, SCREENH-150))

    #===============================================================
    # DRAWING
    #===============================================================

    bg.fill(bg_color)
    pygame.draw.line(bg, (128,128,128), (SCREENW // 2, 0), (SCREENW // 2, SCREENH), 1)
    lux_surf.fill((0,0,0,0))
    fg.fill((0,0,0,0))
    screen.fill((0,0,0))

    draw_snek(x)
    bg.blit(lux_surf, (0, 0), special_flags=pygame.BLEND_RGBA_MULT)
    screen.blit(bg, (0,0))
    screen.blit(fg, (0,0))
    screen.blit(db_screen, (0,0)) if debug else None

    if stage == 1:
        text = "STAGE 1: tune kp, kd for body shape control"
    elif stage == 2:
        text = "STAGE 2: tune α, ω, δ for gait pattern"
    elif stage == 3:
        text = "STAGE 3: play!"
    text = font.render(text, True, WHITE)
    screen.blit(text, ((SCREENW- text.get_width())//2, 100))

    pygame.draw.line(screen, WHITE, (10,SCREENH-20), (10, SCREENH-10))
    pygame.draw.line(screen, WHITE, (10,SCREENH-10), (10+SCREEN_SCALE, SCREENH-10))
    pygame.draw.line(screen, WHITE, (10+SCREEN_SCALE, SCREENH-10), (10+SCREEN_SCALE, SCREENH-20))
    screen.blit(font.render(f"1 m", True, WHITE), (10+SCREEN_SCALE//2-15, SCREENH-30))

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