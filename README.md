# 3D-points-renderer
Here's a 3D points renderer I made from scratch which runs single threaded on the CPU. I tried to make up my own techniques so it's using weird perspective math (not sure what/if there's a name for it). It only renders points now, but I want to add full 3d models, and light and physics simulation when I have the time. The camera control is weird in that it is click and drag, there are other (more conventional) options commented out in the handle_camera_inputs function. 

Try messing with the camera fov constant at the top of the file and how each cube is created/which cubes are rendered. 
WASD movement, R and F to go directly up and down. Have fun! :)
![Screenshot](two_cubes.png)
![Screenshot](360_cube.png)
