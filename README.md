DEVA
==========================
D2AM EDF image Visualisation and Analysis - A software for 2D data analysis. This software deals with special geometry of XPAD area detectors.
If you find any bug, please inform me at: thanhtra0104 @ gmail.com

Thank you.

INSTALLATION on LINUX and MAC OSX
==========================
+ Requirements:
	- Numpy >= 1.8
	- Scipy
	- Matplotlib >= 1.4
	- fabio
	- Cython
	- gcc (C compiler)
	- pyFAI
	- pyGTK -- for MAC, you can find it here: http://sourceforge.net/projects/macpkg/files/PyGTK/
	- lmfit
	- xrayutilities
+ Installation:
	sudo python setup.py install
	
INSTALLATION on WINDOWS
==========================
+ Requirements (all of these packages must be in 32 bits for python 2.7):

	- pythonxy (https://code.google.com/p/pythonxy/)
	
	- fabio (https://pypi.python.org/pypi/fabio), Cython (https://pypi.python.org/pypi/Cython/), lmfit (https://pypi.python.org/pypi/lmfit/), matplotlib (you have to reinstall matplotlib unless an error of backend-gdk will rise) http://matplotlib.org/downloads.html, PyGTK all in one win32 py27 (http://ftp.gnome.org/pub/GNOME/binaries/win32/pygtk/2.24/pygtk-all-in-one-2.24.2.win32-py2.7.msi)
	
	*Optional, but strongly recommended:
	
	- pyFAI (https://pypi.python.org/pypi/pyFAI)
	
	- xrayutilities (https://pypi.python.org/pypi/xrayutilities) -- See help on this website for installation
	
	
	If OpenMP not found or any error of VC++: Try to install these Microsoft packages:
	- Windows SDK for Windows Server 2008 and .NET Framework 3.5: http://www.microsoft.com/en-us/download/details.aspx?id=11310
	- Microsoft Visual C++ 2008 Redistributable Package (x86): http://www.microsoft.com/en-us/download/details.aspx?id=29
	
	Then retry the installation as below.
	
+ Installation:
	
	A binary package will be released soon.

	For the moment let's try to install DEVA with command line on windows: python setup.py install
	
	
RUNNING: 
==========================
to run DEVA, just open a terminal (or command prompt on Windows) and type: Deva.py
	
