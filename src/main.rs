// Wells Crosby 2022

use cgmath::{prelude::*, Euler, Point3, Vector3};
use minifb::{Key, Window, WindowOptions};
use std::f32::consts::PI;
use std::time::Instant;

// window dimensions
const WIDTH: usize = 1400;
const WIDTH_F32: f32 = WIDTH as f32;
const WIDTH_DIV_2_F32: f32 = WIDTH_F32 / 2.0;

const HEIGHT: usize = 800;
const HEIGHT_F32: f32 = HEIGHT as f32;
const HEIGHT_DIV_2_F32: f32 = HEIGHT_F32 / 2.0;

// camera properties
const CAMERA_KEYS_SENSITIVITY: f32 = 0.03;
const CAMERA_DRAG_SENSITIVITY: f32 = 2.0;

// MESS WITH THIS ONE!
// horizontal fov in radians, try setting this to a higher value, even beyond 2 PI
const HORIZONTAL_FOV_RAD: f32 = 0.25 * (PI * 2.0);

// in units per second
const MOVEMENT_SPEED: f32 = 50.0;

// rotation to be added later
struct PointCloud {
    points: Vec<Point3<f32>>,
    center: Point3<f32>,
    orientation: Euler<f32>,
}

// horizontal_fov there for later use with multiple cameras
struct Camera {
    position: Point3<f32>,
    orientation: Euler<f32>,
    horizontal_fov: f32,
    pixels_per_radian: f32,
}

fn main() {
    // our main point cloud, try messing with the create_cube parameters
    let point_cloud_cube: PointCloud = PointCloud {
        points: create_cube(50, true),
        center: Point3::new(0.0, 0.0, 0.0),
        orientation: Euler::new(0.0, 0.0, 0.0),
    };

    // our second point cloud, try messing with the create_cube parameters
    let point_cloud_cube_2: PointCloud = PointCloud {
        points: create_cube(50, false),
        center: Point3::new(0.0, 100.0, 0.0),
        orientation: Euler::new(0.0, 0.0, 0.0),
    };

    // for testing
    let point_cloud_origin: PointCloud = PointCloud {
        points: vec![Point3::new(0.0, 0.0, 0.0)],
        center: Point3::new(0.0, 0.0, 0.0),
        orientation: Euler::new(0.0, 0.0, 0.0),
    };

    // choose which objects to render
    let rendered_objects = vec![point_cloud_cube, point_cloud_cube_2];

    // create our camera
    let mut camera = Camera {
        position: Point3::new(-20.0, 0.0, 30.0),
        orientation: Euler::new(-1.0, 0.0, 0.0),
        horizontal_fov: HORIZONTAL_FOV_RAD,
        pixels_per_radian: WIDTH as f32 / HORIZONTAL_FOV_RAD,
    };

    // pixel buffer and window
    let mut buffer = vec![0; WIDTH * HEIGHT];
    let mut window = Window::new(
        "Test - ESC to exit",
        WIDTH,
        HEIGHT,
        WindowOptions::default(),
    )
    .unwrap_or_else(|e| {
        panic!("{}", e);
    });

    // Colors
    let black = color_val_from_rgb(0, 0, 0);
    let white = color_val_from_rgb(255, 255, 255);
    let red = color_val_from_rgb(255, 0, 0);

    // Limit to max ~60 fps update rate
    window.limit_update_rate(Some(std::time::Duration::from_micros(16600)));

    // timer variable
    let mut frame_time = 0.0;
    let mut prev_mouse_position = Some((0.0, 0.0));
    while window.is_open() && !window.is_key_down(Key::Escape) {
        // start frame timer
        let start = Instant::now();

        // reset screen to all black, this is more efficient than copying an existing one for some reason.
        for i in buffer.iter_mut() {
            *i = black;
        }

        let mouse_position = window.get_mouse_pos(minifb::MouseMode::Discard);
        // handle inputs and get the unit vector of the camera orientation
        let (camera_dir_vec, camera_plane_up) = handle_camera_inputs(
            &window,
            &mut camera,
            frame_time,
            prev_mouse_position,
            mouse_position,
        );
        prev_mouse_position = mouse_position;

        // render the desired point clouds
        for point_cloud in &rendered_objects {
            // figure out where each point should be on the screen
            for point in &point_cloud.points {
                // previous rendering technique, has strong distortion when looking upwards or downwards
                /*
                let angles = point_to_angles_old(&camera, point, camera_dir_vec);
                //println!("{:?}", point_to_angles_2(&camera, point, camera_dir_vec));
                //println!("{:?}", angles);
                if angles.0.abs() < (camera.horizontal_fov / 2.0)
                    && angles.1.abs() < (((camera.horizontal_fov / WIDTH_F32) * HEIGHT_F32) / 2.0)
                {
                    let x = (angles.0 * camera.pixels_per_radian) as isize + WIDTH as isize / 2;
                    let y = -(angles.1 * camera.pixels_per_radian) as isize + HEIGHT as isize / 2;
                    //println!("{}, {}", x, y);
                    buffer[x_y_to_index((x as usize, y as usize))] = white;
                }
                */

                // get the position of the point on the screen
                let (distance, angle) = point_to_position(
                    &camera,
                    add_points(point, &point_cloud.center),
                    camera_dir_vec,
                    camera_plane_up,
                );

                // convert distance and angle to x y coordinates on the window
                let x = (angle.sin() * distance) + WIDTH_DIV_2_F32;
                let y = (-angle.cos() * distance) + HEIGHT_DIV_2_F32;

                // check that the point is in the window
                if x > 0.0 && x < WIDTH_F32 && y > 0.0 && y < HEIGHT_F32 {
                    buffer[x_y_to_index((x as usize, y as usize))] = white;
                }
            }
        }

        // put red dot crosshair in center of screen
        buffer[x_y_to_index((WIDTH / 2, HEIGHT / 2))] = red;
        buffer[x_y_to_index(((WIDTH / 2) + 1, (HEIGHT / 2)))] = red;
        buffer[x_y_to_index(((WIDTH / 2), (HEIGHT / 2) + 1))] = red;
        buffer[x_y_to_index(((WIDTH / 2) - 1, (HEIGHT / 2)))] = red;
        buffer[x_y_to_index(((WIDTH / 2), (HEIGHT / 2) - 1))] = red;

        // We unwrap here as we want this code to exit if it fails
        window.update_with_buffer(&buffer, WIDTH, HEIGHT).unwrap();

        //frames per second
        frame_time = start.elapsed().as_secs_f32();
        println!("{:?}", 1.0 / frame_time);
    }
}

// point addition
fn add_points(first: &Point3<f32>, second: &Point3<f32>) -> Point3<f32> {
    Point3::new(first.x + second.x, first.y + second.y, first.z + second.z)
}

// creates a vec of points in the shape of a cube, can be either filled or hollow
// this is pretty inefficient for the non-filled in version
fn create_cube(size: usize, fill_in: bool) -> Vec<Point3<f32>> {
    let mut my_vec: Vec<Point3<f32>> = Vec::new();
    for i in 0..size {
        for j in 0..size {
            for k in 0..size {
                if fill_in {
                    my_vec.push(Point3::new(i as f32, j as f32, k as f32));
                } else if i == 0
                    || j == 0
                    || k == 0
                    || i == size - 1
                    || j == size - 1
                    || k == size - 1
                {
                    my_vec.push(Point3::new(i as f32, j as f32, k as f32));
                }
            }
        }
    }

    my_vec
}

// get the inputs from the user and manipulate the camera accordingly, return the unit vector
// in the direction of the camera's orientation and the vector pointing 90 deg up from the camera vector
fn handle_camera_inputs(
    window: &Window,
    camera: &mut Camera,
    frame_time: f32,
    prev_mouse_position: Option<(f32, f32)>,
    mouse_position: Option<(f32, f32)>,
) -> (Vector3<f32>, Vector3<f32>) {
    // control camera orientation with click and drag method
    // check that mouse cursor in window
    if mouse_position.is_some() {
        // check if left mouse button is down
        if window.get_mouse_down(minifb::MouseButton::Left) {
            // check that previous frame's mouse
            if prev_mouse_position.is_some() {
                // find out how many x and y pixels the mouse has moved since the last frame scaled by CAMERA_DRAG_SENSITIVITY
                let pixel_movement = (
                    (mouse_position.unwrap().0 - prev_mouse_position.unwrap().0)
                        * CAMERA_DRAG_SENSITIVITY,
                    -(mouse_position.unwrap().1 - prev_mouse_position.unwrap().1)
                        * CAMERA_DRAG_SENSITIVITY,
                );
                //println!("{:?}", (pixel_movement.0 / camera.pixels_per_radian));

                // convert the pixel movement to a new x and y camera orientation
                // find angles by pixels moved converted to angle moved using camera_pixels_per_radian plus the old orientation
                let moved = (
                    (pixel_movement.0 / camera.pixels_per_radian) + camera.orientation.x,
                    (pixel_movement.1 / camera.pixels_per_radian) + camera.orientation.y,
                );
                //println!("{:?}", movement);

                // move y component of camera orientation, keeping it in range
                if moved.1 > PI / 2.0 {
                    camera.orientation.y = (PI / 2.0) - 0.0001;
                } else if moved.1 < -PI / 2.0 {
                    camera.orientation.y = (-PI / 2.0) + 0.0001;
                } else {
                    camera.orientation.y = moved.1;
                }

                // move x component of camera orientation,
                if moved.0 > PI {
                    camera.orientation.x = moved.0 - (2.0 * PI);
                } else if moved.0 < -PI {
                    camera.orientation.x = moved.0 + (2.0 * PI);
                } else {
                    camera.orientation.x = moved.0;
                }
            }
        }
    }
    //println!("{:?}", mouse_position);

    // control camera orientation with the position of the mouse on the window
    /*
    if mouse_position.is_some() {
        camera.orientation.x = (mouse_position.unwrap().0 * ((PI * 2.0) / WIDTH as f32)) - PI;
        camera.orientation.y = (PI / 2.0) - (mouse_position.unwrap().1 * (PI / HEIGHT as f32));
    }
    */

    // control camera orientation with arrow keys
    /*
    if window.is_key_down(Key::Up) {
        if camera.orientation.y < (PI / 2.0) - CAMERA_KEYS_SENSITIVITY {
            camera.orientation.y += CAMERA_KEYS_SENSITIVITY;
        } else {
            camera.orientation.y = PI / 2.0 - 0.0001;
        }
    }
    if window.is_key_down(Key::Down) {
        if camera.orientation.y > CAMERA_KEYS_SENSITIVITY - (PI / 2.0) {
            camera.orientation.y -= CAMERA_KEYS_SENSITIVITY;
        } else {
            camera.orientation.y = -PI / 2.0 + 0.0001;
        }
    }
    if window.is_key_down(Key::Left) {
        if camera.orientation.x - CAMERA_KEYS_SENSITIVITY < -PI {
            //println!("left rollover");
            camera.orientation.x = camera.orientation.x - CAMERA_KEYS_SENSITIVITY + (PI * 2.0);
        } else {
            camera.orientation.x -= CAMERA_KEYS_SENSITIVITY;
        }
    }
    if window.is_key_down(Key::Right) {
        if camera.orientation.x + CAMERA_KEYS_SENSITIVITY > PI {
            camera.orientation.x = camera.orientation.x + CAMERA_KEYS_SENSITIVITY - (PI * 2.0);
        } else {
            camera.orientation.x += CAMERA_KEYS_SENSITIVITY;
        }
    }
    */

    // unit vector of camera direction
    let camera_dir_vec = euler_to_vector(&camera.orientation);

    // the vector pointing 90 deg up from the camera vector
    let camera_plane_up = Vector3::new(
        -camera_dir_vec.x * camera_dir_vec.z,
        -camera_dir_vec.y * camera_dir_vec.z,
        camera_dir_vec.x.powi(2) + camera_dir_vec.y.powi(2),
    );

    // move camera forwards and backwards
    if window.is_key_down(Key::W) {
        camera.position += camera_dir_vec * MOVEMENT_SPEED * frame_time;
    }
    if window.is_key_down(Key::S) {
        camera.position -= camera_dir_vec * MOVEMENT_SPEED * frame_time;
    }

    // move camera left and right (this is not the most efficient way possible)
    // (also doesn't take into account camera rotation)
    let unit_vec_left_of_camera = camera_plane_up.cross(camera_dir_vec)
        + camera_plane_up * (camera_plane_up.dot(camera_dir_vec));

    if window.is_key_down(Key::A) {
        camera.position += unit_vec_left_of_camera * MOVEMENT_SPEED * frame_time;
    }
    if window.is_key_down(Key::D) {
        camera.position -= unit_vec_left_of_camera * MOVEMENT_SPEED * frame_time;
    }

    // move camera up or down
    if window.is_key_down(Key::R) {
        camera.position.z += MOVEMENT_SPEED * frame_time;
    }
    if window.is_key_down(Key::F) {
        camera.position.z -= MOVEMENT_SPEED * frame_time;
    }

    // rotate camera
    if window.is_key_down(Key::Q) {
        if camera.orientation.z < PI - CAMERA_KEYS_SENSITIVITY {
            camera.orientation.z += CAMERA_KEYS_SENSITIVITY;
        } else {
            camera.orientation.z = PI - 0.0001;
        }
    }
    if window.is_key_down(Key::E) {
        if camera.orientation.z > CAMERA_KEYS_SENSITIVITY - PI {
            camera.orientation.z -= CAMERA_KEYS_SENSITIVITY;
        } else {
            camera.orientation.z = PI - 0.0001;
        }
    }

    // these must be calculated here for most accuracy
    (camera_dir_vec, camera_plane_up)
}

// x and y position to buffer index
fn x_y_to_index(x_y: (usize, usize)) -> usize {
    x_y.0 + (WIDTH * x_y.1)
}

// rgb values to color value, taken from example code
fn color_val_from_rgb(r: u8, g: u8, b: u8) -> u32 {
    let (r, g, b) = (r as u32, g as u32, b as u32);
    (r << 16) | (g << 8) | b
}

// convert euler coordinates to a unit vector
fn euler_to_vector(euler: &Euler<f32>) -> Vector3<f32> {
    Vector3::new(
        euler.x.cos() * euler.y.cos(),
        -euler.x.sin() * euler.y.cos(),
        euler.y.sin(),
    )
}

// finds the position of a point in the window by finding how far it is from
// the center and at what angle it is from facing upwards
fn point_to_position(
    camera: &Camera,
    point: Point3<f32>,
    camera_dir_vec: Vector3<f32>,
    camera_plane_up: Vector3<f32>,
) -> (f32, f32) {
    // vector from camera position to point
    let camera_to_point = point - camera.position;

    // cosine of the angle between the camera unit vec and camera to point vec
    let cos_angle = camera_to_point.dot(camera_dir_vec) / camera_to_point.magnitude();

    // smallest vector from line extending from camera_dir_vec to the point
    let point_on_camera_plane = camera_to_point - cos_angle * camera_dir_vec;

    // find the signed angle between camera_plane_up and point_on_camera_plane
    let plane_angle = (camera_plane_up.cross(point_on_camera_plane))
        .dot(camera_dir_vec)
        .atan2(camera_plane_up.dot(point_on_camera_plane));

    // (distance from center of window in pixels, angle from center pointing up)
    (
        cos_angle.acos() * camera.pixels_per_radian,
        plane_angle + camera.orientation.z,
    )
}

// old function for finding position of point on window
fn point_to_angles_old(
    camera: &Camera,
    point: &Point3<f32>,
    camera_dir_vec: Vector3<f32>,
) -> (f32, f32) {
    let camera_to_point = point - camera.position;

    //println!("{:?}", camera_to_point);
    let vertical_angle =
        (camera_to_point.z / (camera_to_point.x.powi(2) + camera_to_point.y.powi(2)).sqrt()).atan()
            - camera.orientation.y;

    let camera_to_point_angle = camera_to_point.angle(camera_dir_vec).0;

    let mut horizontal_angle = (camera_to_point_angle.powi(2) - vertical_angle.powi(2)).sqrt();

    if camera_dir_vec.x > 0.0 {
        if point.y
            > (camera_dir_vec.y / camera_dir_vec.x) * (point.x - (camera.position.x))
                + camera.position.y
        {
            horizontal_angle = -horizontal_angle;
        }
    } else {
        if point.y
            < (camera_dir_vec.y / camera_dir_vec.x) * (point.x - (camera.position.x))
                + camera.position.y
        {
            horizontal_angle = -horizontal_angle;
        }
    }

    //println!("{:?}", vertical_angle);
    (horizontal_angle, vertical_angle)
}
