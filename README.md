I have added a code. In that basically I have used Lax Friedrich scheme to sole shock tube problem. only visualization is an issue. GPT and CLaude gave an HTML based method to visualize. It works very well. But Sir recommended GNUPlot. I want someone to try to get the visualization by GNUPlot.
-Aditya Kini

New changes are to be made in the shock tube with reflection. @Aditya Kumar do the visualization using GNU Plot for the shock tube with reflection. We need graphs density, velocity, pressure and temperature graphs and contours(for 2D case later on). @Jaga and @Archit do the validation using the shock tube with reflection case. Validate the schemes by comparing with existing numerical or analytical solutions for shock tube problems after Aditya has put the GNU Plot code. 

After this we extend the code to 2D and implement immersed boundary method.
-Aditya Kini


I have added the gnu-plot changes for only shock-tube everyone have a look once i think there is a issue in the density plot may be issue with the dat file i will try to rectify that everyone have a look once - Aditya Kumar




I have added 2D code. Archit, Jaga and Aditya Kumar, I want you guys to validate the 2D code for Wedge angles of 10 15 20 degrees for mesh sizes of 300x300,500x500,1000x1000. Be sure to change the sampling points in geometry.cpp appropriately. I would suggest doing it in CADlab as it took a while to run on my computer. In VS studio go to main.cpp and press f5 to compile all together. No need to change inlet mach number.

For results I want the following: Mach number contours, and pressure contours. Also I want you guys to measure the shock angle and cross verify shock angle,pressure and mach numbers by theta- beta -M tables from last sem gas tables. In the contours be sure to demarcate the wedge properly in the final image.
-Aditya Kini
