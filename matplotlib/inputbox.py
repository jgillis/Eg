#! /usr/bin/env python
# from the web
from Tkinter import *

class Window(Frame):
	def __init__(self, parent=None):
		Frame.__init__(self, parent)
		Label(self, text="Enter the path to a gif").pack(padx=5, pady=5)
		self.label = Label(self, text="")
		self.label.pack(padx=5, pady=5)
		self.entry = Entry(self, text="")
		self.entry.pack(padx=5, pady=5)

		self.pack()
		self.update()

	def update(self):
		try:
			self.image = PhotoImage(file=self.entry.get())
			self.label.config(image=self.image)
		except TclError:
			self.label.config(text=self.entry.get(), image="")
		self.after(20, self.update)

if __name__ == '__main__':
	root = Tk()
	Window().mainloop()