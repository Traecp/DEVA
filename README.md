DEVA
==========================
D2AM EDF image Visualisation and Analysis - A software for 2D data analysis

INSTALLATION on LINUX
==========================
+ Requirements:
	- Numpy >= 1.8
	- Scipy
	- Matplotlib >= 1.4
	- fabio
	- Cython
+ Installation:
	sudo python setup.py install
	
INSTALLATION on WINDOWS
==========================
+ Requirements:
	- pythonxy (https://code.google.com/p/pythonxy/)
	
	and 
	
	- fabio (https://pypi.python.org/pypi/fabio), Cython (https://pypi.python.org/pypi/Cython/), lmfit (https://pypi.python.org/pypi/lmfit/), matplotlib (you have to reinstall matplotlib unless an error of backend-gdk will rise), PyGTK all in one win32 py27 (http://ftp.gnome.org/pub/GNOME/binaries/win32/pygtk/2.24/pygtk-all-in-one-2.24.2.win32-py2.7.msi)
	
	Optional:
	
	- pyFAI, xrayutilities
	
	
	If OpenMP not found: Try to install these Microsoft packages:
	- Windows SDK for Windows Server 2008 and .NET Framework 3.5: http://www.microsoft.com/en-us/download/details.aspx?id=11310
	- Microsoft Visual C++ 2008 Redistributable Package (x86): http://www.microsoft.com/en-us/download/details.aspx?id=29
	
	Then retry the installation as below.
	
+ Installation:
	
	A binary package will be released soon.

	For the moment let's try to install DEVA with command line on windows: python setup.py install -c mingw32
	
	RUNNING: to run DEVA, just open a terminal (or command prompt) and type: Deva.py
	
