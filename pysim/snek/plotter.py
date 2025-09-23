import pygame
import sys
import numpy as np

WHITE = (255, 255, 255)

class PlotFigure:
    def __init__(self, w, h, font=None, axiscolor=WHITE, margin=5):
        # assert isinstance(w, int) and w > 0, "width must be a positive integer"
        # assert isinstance(h, int) and h > 0, "height must be a positive integer"
        # assert font is None or isinstance(font, pygame.font.Font), "font must be a pygame font or None"
        # assert isinstance(axiscolor, (tuple, list, np.ndarray)) and len(axiscolor) == 3, "axiscolor must be an RGB tuple"
        # assert all(0 <= c <= 255 for c in axiscolor), "axiscolor value is invalid, must be in range 0-255"
        # assert isinstance(margin, int) and margin >= 0, "margin must be a non-negative integer"

        self.w = w
        self.h = h
        self.surf = pygame.Surface((self.w, self.h))

        self.xlim = None
        self.ylim = None

        self.xlabel = ''
        self.ylabel = ''
        self.font = font
        self.axiscolor = axiscolor

        self.xticks = None
        self.yticks = None

        self.margin = margin
        self.hold = False
    
    """ Set axis limits for both x and y axis. If 'auto', it will be set based on data range. """
    def _xlim(self, lim, X):
        if lim == 'auto':
            self.xlim = (min(X), max(X))
        else:
            # assert isinstance(lim, (list, tuple)) and len(lim) == 2 and lim[0] < lim[1], "xlim must be an increasing tuple"
            self.xlim = lim

    def _ylim(self, lim, Y):
        if lim == 'auto':
            self.ylim = (min(Y), max(Y))
        else:
            # assert isinstance(lim, (list, tuple)) and len(lim) == 2 and lim[0] < lim[1], "ylim must be an increasing tuple"
            self.ylim = lim
    
    def hold_on(self):
        self.hold = True
    def hold_off(self):
        self.hold = False

    """ Plot XY data on the surface. """
    def plot(self, X,Y, xlim='auto', ylim='auto', line_width=1, color=WHITE):

        #===============================================================
        # checks on XY data
        #===============================================================

        # assert isinstance(X, (list, tuple, np.ndarray)), "X must be a list or tuple"
        # assert isinstance(Y, (list, tuple, np.ndarray)), "Y must be a list or tuple"
        # assert len(X) == len(Y), "X and Y must be the same length"
        if len(X) < 2: #nothing to plot
            return

        #===============================================================
        # set axis limits
        #===============================================================

        self._xlim(xlim, X)
        self._ylim(ylim, Y)

        if self.xlim[0] == self.xlim[1] or self.ylim[0] == self.ylim[1]:
            return #nothing to plot
        
        #===============================================================
        # checks on line parameters
        #===============================================================
        
        # assert isinstance(line_width, int) and line_width > 0, "line_width must be a positive integer"
        # assert isinstance(color, tuple) and len(color) == 3, "color must be an RGB tuple"
        # assert all(0 <= c <= 255 for c in color), "color value is invalid, must be in range 0-255"

        # scale XY data based on surface size and axis limits

        sx = self.w/(self.xlim[1]-self.xlim[0])
        sy = self.h/(self.ylim[1]-self.ylim[0])

        px1 = int((X[0]-self.xlim[0])*sx)
        py1 = self.h - int((Y[0]-self.ylim[0])*sy)
    
        if not self.hold:
            self.surf.fill((0,0,0))
        
        for i in range(1,len(X)):
            px2 = int((X[i]-self.xlim[0])*sx)
            py2 = self.h - int((Y[i]-self.ylim[0])*sy)

            pygame.draw.line(self.surf, color, (px1, py1), (px2, py2), line_width)
            # and then go to next point
            px1, py1 = px2, py2

        pygame.draw.rect(self.surf, (255, 255, 255), (0, 0, self.w, self.h), 1)

    def label(self, xlabel=None, ylabel=None):
        # assert isinstance(xlabel, str), "xlabel must be a string"
        # assert isinstance(ylabel, str), "ylabel must be a string"
        self.xlabel = xlabel
        self.ylabel = ylabel

    def _xticks(self, nums):
        # assert isinstance(nums, int) and nums > 1, "nums must be a positive integer greater than 1"
        self.xticks = np.linspace(self.xlim[0], self.xlim[1], nums)

    def _yticks(self, nums):
        # assert isinstance(nums, int) and nums > 1, "nums must be a positive integer greater than 1"
        self.yticks = np.linspace(self.ylim[0], self.ylim[1], nums)
    
    def _draw_ticks(self):
        if self.xticks is None and self.yticks is None:
            return
        
        # assert self.font is not None, "font must be set to draw ticks"

        sx = self.w/(self.xlim[1]-self.xlim[0])
        sy = self.h/(self.ylim[1]-self.ylim[0])

        dummy_tick = self.font.render(f"{-88.88:.2f}", True, self.axiscolor)
        w = dummy_tick.get_size()[0]
        h = dummy_tick.get_size()[1]

        if self.xticks is not None:
            for x in self.xticks:
                px = int((x - self.xlim[0]) * sx)
                pygame.draw.line(self.surf, self.axiscolor, (px, self.h), (px, self.h - self.margin), 1)

        if self.yticks is not None:
            for y in self.yticks:
                py = self.h - int((y - self.ylim[0]) * sy)
                pygame.draw.line(self.surf, self.axiscolor, (0, py), (self.margin, py), 1)
        
        bigger_surf = pygame.Surface((self.w+w+self.margin, self.h+h+self.margin))
        bigger_surf.blit(self.surf, (w+self.margin, 0))

        if self.xticks is not None:
            for x in self.xticks:
                px = int((x - self.xlim[0]) * sx)
                text = self.font.render(f"{x:.2g}", True, self.axiscolor)
                text_size = text.get_size()
                bigger_surf.blit(text, (w+self.margin + px - text_size[0]//2, self.h+self.margin))
        if self.yticks is not None:
            for y in self.yticks:
                py = self.h - int((y - self.ylim[0]) * sy)
                text = self.font.render(f"{y:.2g}", True, self.axiscolor)
                text_size = text.get_size()
                bigger_surf.blit(text, (w - text_size[0], py - text_size[1]//2))
        
        pygame.draw.rect(bigger_surf, (255,0,0), (w, 0, self.w+w, self.h+self.margin), 1) #box around plot
        pygame.draw.rect(bigger_surf, (0,255,0), (w+self.margin, 0, self.w, self.h+1), 1) #box around plot
        self.surf = bigger_surf

    def _draw_labels(self):
        if not self.xlabel and not self.ylabel:
            return
        
        # assert self.font is not None, "font must be set to draw labels"
        
        sizex = (0,0)
        sizey = (0,0)
        if self.xlabel:
            textx = self.font.render(self.xlabel, True, self.axiscolor)
            sizex = textx.get_size()
        if self.ylabel:
            texty = self.font.render(self.ylabel, True, self.axiscolor)
            sizey = texty.get_size()

        bigger_surf = pygame.Surface((self.surf.get_width()+sizey[0]+self.margin, self.surf.get_height()+sizex[1]+self.margin))
        bigger_surf.blit(self.surf, (sizey[0]+self.margin, 0)) #plot it topright so theres space for labelx on the left and labely on the bottom
        bigger_surf.blit(textx, (sizey[0]+self.margin + (self.surf.get_width()-sizex[0])//2, self.surf.get_height()+self.margin))
        bigger_surf.blit(texty, (0, (self.surf.get_height()-sizey[1])//2))
        
        pygame.draw.rect(bigger_surf, (255,0,0), (sizey[0], 0, self.surf.get_width()+sizey[0], self.surf.get_height()+self.margin), 1) #box around plot
        pygame.draw.rect(bigger_surf, self.axiscolor, (sizey[0]+self.margin, 0, self.surf.get_width(), self.surf.get_height()), 1) #box around plot

        self.surf = bigger_surf

#===============================================================
# USAGE
#===============================================================

# pygame.init()
# screen = pygame.display.set_mode((800, 600))
# font = pygame.font.Font(None, 24)
# clock = pygame.time.Clock()

# plotf = PlotFigure(600, 400, font=font, axiscolor=WHITE, margin=5)

# while True:
#     for event in pygame.event.get():
#         if event.type == pygame.QUIT:
#             pygame.quit()
#             sys.exit()
        
#     time = pygame.time.get_ticks() / 1000
#     X = np.arange(time-6, time, 0.01)

#     plotf.plot(X, np.sin(X), xlim=(time-6, time), ylim=(-1.5, 1.5), color=(255,0,0))
#     plotf.hold_on()
#     plotf.plot(X, np.cos(X), xlim=(time-6, time), ylim=(-1.5, 1.5), color=(255,0,0))
#     plotf.plot(X, np.zeros_like(X), xlim=(time-6, time), ylim=(-1.5, 1.5), color=(255,0,0))
#     plotf.hold_off()
#     plotf._xticks(10)
#     plotf._yticks(13)
#     plotf._draw_ticks()
#     plotf.label('Time (s)', 'Value')
#     plotf._draw_labels()
    

#     screen.fill((60,60,60))
#     screen.blit(plotf.surf, (0, 0))
#     pygame.display.flip()