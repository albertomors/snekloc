import pygame
import sys

GRAY = (150,150,150)
WHITE = (255, 255, 255)

#===============================================================
# BUTTONS
#===============================================================

class Button():
    def __init__(self, title, cx, cy, w, h, val1, val2):
        self.title = title
        self.cx = cx
        self.cy = cy
        self.w = w
        self.h = h
        self.val1 = val1
        self.val2 = val2

        #calculate sizes
        self.rect = pygame.Rect(self.cx - self.w//2, self.cy - self.h // 2, self.w, self.h)
        self.surf = pygame.Surface((self.w, self.h), pygame.SRCALPHA)

        self.clicked = False
        self.hovered = False
        self.off()
    
    """ draw at initialization static stuff """
    def draw_once(self, font, screen):
        pygame.draw.rect(self.surf, GRAY, (0, 0, self.w, self.h))

    """ draw dynamic stuff """
    def draw(self, font, screen):
        screen.blit(self.surf, (self.rect.x, self.rect.y))

        text = font.render(f"{self.title}", True, WHITE)
        screen.blit(text, (self.rect.x - text.get_width() - 5, self.rect.y + (self.h - text.get_height())//2))

        text = font.render(f"{self.val}", True, WHITE)
        screen.blit(text, (self.rect.x + (self.w-text.get_width())//2, self.rect.y + (self.h-text.get_height())//2))

    """ update """
    def handle_event(self, event):
        if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1 and self.rect.collidepoint(event.pos):
            self.clicked = True
        elif event.type == pygame.MOUSEBUTTONUP and event.button == 1 and self.rect.collidepoint(event.pos):
            self.clicked = False
        
        if event.type == pygame.MOUSEMOTION:
            if self.rect.collidepoint(event.pos):
                self.hovered = True
            else:
                self.hovered = False

    def update(self):
        if self.clicked:
            self.toggle()

    """ clicks """
    def toggle(self):
        self.val = self.val2 if self.val == self.val1 else self.val1
    def on(self):
        self.val = self.val2
    def off(self):
        self.val = self.val1

    def get_value(self):
        return self.val

#===============================================================

class SquareButton(Button):
    def __init__(self, title, cx, cy, size=20, val1=0, val2=1):
        super().__init__(title, cx, cy, size, size, val1=val1, val2=val2)
        
#===============================================================

class StageButton(Button):
    def __init__(self, title, cx, cy, w, h, val1=0, val2=1):
        super().__init__(title, cx, cy, w, h, val1=val1, val2=val2)
    
    def update(self):
        if self.clicked:
            self.surf.set_alpha(255)
            self.on() # saturate to on, turned off from outside to avoid multiple toggles
        elif self.hovered and not self.clicked:
            self.surf.set_alpha(128)
            self.off()
        else:
            self.surf.set_alpha(32)
    
    def draw_once(self, font, screen):
        super().draw_once(font, screen)
        text = font.render(f"{self.title}", True, WHITE)
        self.surf.blit(text, (max(0, (self.w-text.get_width())//2), max(0, (self.h-text.get_height())//2)))

    def draw(self, font, screen):
        screen.blit(self.surf, (self.rect.x, self.rect.y))

#===============================================================

class SliderButton(StageButton):
    def __init__(self, title, cx, cy, size, val1=0, val2=1):
        super().__init__(title, cx, cy, size, size, val1=val1, val2=val2)
    
    def update(self):
        if self.clicked:
            self.on()
        else:
            self.off()

    def draw_once(self, font, screen):
        super().draw_once(font, screen)
        text = font.render(f"{self.title}", True, WHITE)
        self.surf.blit(text, (max(0, (self.w-text.get_width())//2), max(0, (self.h-text.get_height())//2)))

    def draw(self, font, screen):
        screen.blit(self.surf, (self.rect.x, self.rect.y))

#===============================================================
# SLIDER
#===============================================================

class Slider:
    def __init__(self, title, x, y, w, h, margin, min=0, init_val=50, max=100):
        self.title = title
        self.x = x
        self.y = y
        self.w = w
        self.h = h
        self.margin = margin

        self.min = min
        self.max = max
        self.val = init_val
        self.step = (self.max - self.min) / 100  # default step is 1% of range

        self.dragging = False

        #calculate sizes
        self.l_btn = SliderButton("<", self.x - self.margin - self.h // 2, self.y + self.h // 2, size=self.h)
        self.r_btn = SliderButton(">", self.x + self.w + self.margin + self.h//2, self.y + self.h // 2, size=self.h)
        self.s_pos = self.x + int((self.val - self.min) / (self.max - self.min) * self.w)
        self.s_rect = pygame.Rect(self.s_pos - self.margin//2, self.y - self.margin//2, self.margin, self.h + self.margin)
    
    def draw_once(self, small_font, screen):
        self.l_btn.draw_once(small_font, screen)
        self.r_btn.draw_once(small_font, screen)

    def draw(self, font, screen):
        # slider track
        pygame.draw.rect(screen, (100, 100, 100), (self.x, self.y, self.w, self.h))
        
        # slider handle
        self.s_pos = self.x + int((self.val - self.min) / (self.max - self.min) * self.w)
        self.s_rect.x = self.s_pos - self.margin // 2
        pygame.draw.rect(screen, (200, 200, 200), self.s_rect)
        
        # buttons
        # create a slightly smaller font based on the passed-in font and use it for the buttons
        self.l_btn.draw(font, screen)
        self.r_btn.draw(font, screen)

        # Draw val text
        text = font.render(f"{self.title}: {self.val:.2f}", True, (255, 255, 255))
        screen.blit(text, (self.x + self.w  + 2*self.margin + self.h, self.y + self.h // 2 - text.get_height() // 2))

    def handle_event(self, event):
        self.l_btn.handle_event(event)
        self.r_btn.handle_event(event)
        self.l_btn.update()
        self.r_btn.update()

        if self.l_btn.get_value() == self.l_btn.val2:  # if button is "on"
            self.decrement()
        elif self.r_btn.get_value() == self.r_btn.val2:
            self.increment()

        if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:  # Left mouse button
        # Check if slider is clicked
            if self.x <= event.pos[0] <= self.x + self.w and self.y - self.margin//2 <= event.pos[1] <= self.y + self.h + self.margin//2:
                self.dragging = True
                # Update slider position directly
                self.set_value_from_mouse_pos(event.pos[0])
        
        elif event.type == pygame.MOUSEBUTTONUP:
            if event.button == 1:  # Left mouse button
                self.dragging = False
        
        elif event.type == pygame.MOUSEMOTION:
            if self.dragging:
                # Update slider position
                self.set_value_from_mouse_pos(event.pos[0])

    def set_value_from_mouse_pos(self, x_pos):
        # Convert mouse position to val
        relative_pos = max(0, min(x_pos - self.x, self.w))
        self.val = self.min + (relative_pos / self.w) * (self.max - self.min)
        # Round to nearest step
        # self.val = round(self.val / self.step) * self.step
        # Ensure val is within range
        self.val = max(self.min, min(self.max, self.val))

    def increment(self):
        self.val = min(self.max, self.val + self.step)

    def decrement(self):
        self.val = max(self.min, self.val - self.step)

    def get_value(self):
        return self.val

class StandardSlider(Slider):
    def __init__(self, title, cx, cy, min=0, init_val=0.5, max=1):
        w, h = 50, 10
        margin = h//2
        x = cx - w//2
        y = cy - h//2
        super().__init__(title, x, y, w, h, margin, min, init_val, max)

#===============================================================

def get_size(Slider):
    tot_w = Slider.w + 2 * Slider.margin + 2 * Slider.h
    tot_h = Slider.h + Slider.margin
    return tot_w, tot_h

slider_size = get_size(StandardSlider("test", 0, 0))
btn_size = SquareButton("test", 0, 0).w