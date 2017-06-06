The Sample Optimization 3DCoils 

Unpack all the files in the archive to a directory in the MATLAB path.
Start MATLAB and start the bnb user interface by typing 'bnbgui' at the MATLAB prompt.
Load '3DCoils.mat' in bnbgui ('file'->'load').
Update the function ('function'->'fun') to make sure it points to the file 'F3DCoils.m'.
Update the nonlinear constraints ('function'->'nonlincon') to point to the file 'C3DCoils.m'.
Start the optimization by pressing the 'optimize' button.

3dspoelen concerns the optimization of a device for measuring the magnetic field in 3 dimensions. The device is made up of 3 coils (for the 3 dimensions) and an op-amp for each coil (for amplification of the signal). See the picture '3DCoils.bmp'. Goal was to maximize the signal-to-noise-ratio. The same wire and the same op-amp had to be used for the three coils.

The number of possible op-amps and possible copper wires was limited and I had to use integer variables to force the optimizer to make the choice. Whenever there was a choice between different options (6 op-amps and 10 wires in this case) I used the word 'set' to denote the set of different options.
 


