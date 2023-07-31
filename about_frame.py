#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 15:22:49 2023

@author: atilio
"""

from tkinter import ttk
import tkinter as tk
import webbrowser

def create_about():
    
    root = tk.Tk()
    root.title('About')
    
    window_width = 480
    window_height = 350
    
    # get the screen dimension
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    
    # find the center point
    center_x = int(screen_width/2 - window_width / 2)
    center_y = int(screen_height/2 - window_height / 2)
    
    # set the position of the window to the center of the screen
    root.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')
    
    #root.geometry("%dx%d" % (screen_width-2, screen_height-10))
    
    #root.resizable(False, False)
    
    title = ttk.Label(root, text='PolyWeight',font=("Arial", 32, "bold"))
    title.pack(pady=(28,0))
    
    version = ttk.Label(root, text='Version 0.5.0',font=("Arial", 16))
    version.pack(pady=(18,0))
    
    dev = ttk.Label(root, text='Developed by Atilio Minotto',font=("Arial", 16))
    dev.pack(pady=((32,0)))
    
    
    def callback(url):
        webbrowser.open_new(url)
        
    email_link = ttk.Label(root, text="minotto93@gmail.com | amneto2@ucs.br",font=("Arial", 16))
    email_link.pack(pady=(18,18))


    link1 = ttk.Label(root, text="GitHub",cursor="hand2",font=("Arial", 16))
    link1.pack(pady=(18,0))
    link1.bind("<Button-1>", lambda e: callback("https://github.com/a-minotto/PolyWeight"))
    link1.config(underline=True,foreground='blue')

    link2 = ttk.Label(root, text="License",cursor="hand2",font=("Arial", 16))
    link2.pack(pady=(18,0))
    link2.bind("<Button-1>", lambda e: callback("https://www.gnu.org/licenses/gpl-3.0.en.html"))
    link2.config(underline=True,foreground='blue')

    
    root.mainloop()