# 3D-points-renderer
Here's a 3D points renderer I made from scratch in Rust which runs single threaded on the CPU. I tried to make up my own techniques so it's using weird perspective math (not sure what it's called/if there's a name for it). It only renders points now, but I want to add full 3d models, and light and physics simulation when I have the time. The camera control is click and drag, there are other (more conventional) options commented out in the handle_camera_inputs function. 

Try messing with the camera fov constant at the top of the file and how each cube is created/which cubes are rendered in main. 
WASD movement, R and F to go directly up and down. Have fun! :)

Here's a filled-in and non-filled-in cube, the patterns on the filled-in cube are cool to see in motion, especially as you move far away.
![two_cubes](https://user-images.githubusercontent.com/64884461/208126102-2be0d7c8-9827-46f0-a169-5a28868ed2f1.PNG)

Here's what it's like inside a non-filled-in cube with a 360 degree camera, looks really cool in motion.
![360_cube](https://user-images.githubusercontent.com/64884461/208126077-1e7b9f3f-25bd-43ad-92ae-92e386a50231.PNG)

Randomized cube and color support
![red_white_random_cubes](https://user-images.githubusercontent.com/64884461/208986667-cb2f7823-6c49-43e4-939f-ad09f774be7f.PNG)
