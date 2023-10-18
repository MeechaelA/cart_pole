use std::fs::OpenOptions;
use std::io::Write;
use nalgebra::{Matrix2x1, Matrix2};

fn integrate_local(x_k: f64, x_k_m1: f64, time_delta: f64) -> f64{
    let y = ((x_k + x_k_m1) * 0.5) * time_delta; 
    return y;
}

fn main() {

    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .open("post/output.csv")
        .unwrap();


    let mut cart_pos: f64 = 0.0;
    let mut prev_cart_pos: f64 = 0.0;

    let mut cart_vel: f64 = 0.0;
    let mut prev_cart_vel: f64 = 0.0;
    
    let mut cart_accel: f64 = 0.0;
    let mut prev_cart_accel: f64 = 0.0;

    let cart_mass: f64 = 50.0; // Grams 

    let mut pole_angle_pos: f64 = 0.0; // Radians
    let mut prev_pole_angle_pos: f64 = pole_angle_pos;
    
    let mut pole_angle_vel: f64 = 0.0;
    let mut prev_pole_angle_vel: f64 = 0.0;
    
    let mut pole_angle_accel: f64 = 0.0;
    let mut prev_pole_angle_accel: f64 = 0.0;

    let pole_length: f64 = 0.5; // Meters
    let pole_mass: f64 = 25.0; // Grams


    let gravity = 9.81;
    let mut force = 0.0;

    let mut time = 0.0;
    let mut time_prev = time;
    let time_delta = 0.00001;
    let mut times: Vec<f64> = vec![];
    times.push(time);

    let mut sim_index = 0;

    let mut solved = false;
    
    // let output = format!("{:?},{:?},{:?},{:?},{:?},{:?},{:?}", time, cart_accel, cart_vel, cart_pos, pole_angle_accel, pole_angle_vel, pole_angle_pos);   
    let output = format!("{:?},{:?},{:?},{:?},{:?},{:?},{:?}", "TIME (s)", "Cart Acceleration (m/s/s)", "Cart Velocity (m/s)", "Cart Position (m)", "Pole Angular Acceleration (rad/s/s)", "Pole Angular Velocity (rad/s)", "Pole Angular Position (rad)");
    writeln!(file, "{}", output);

    
    while (solved != true){

        
        let a_matrix: nalgebra::Matrix<f64, nalgebra::Const<2>, nalgebra::Const<2>, nalgebra::ArrayStorage<f64, 2, 2>> = Matrix2::new(pole_angle_pos.cos(),   pole_length, cart_mass + pole_mass, pole_mass * pole_length * pole_angle_pos.cos());

        let b_matrix: nalgebra::Matrix<f64, nalgebra::Const<2>, nalgebra::Const<1>, nalgebra::ArrayStorage<f64, 2, 1>> = Matrix2x1::new(-1.0 * gravity * pole_angle_pos.sin(), force + pole_mass * pole_length * cart_vel.powf(2.0) * pole_angle_pos.sin());

        let a_matrix_inverse = a_matrix.try_inverse();

        match a_matrix_inverse{
            Some(a_matrix_solution)=>{

                let x = a_matrix_solution * b_matrix;

                cart_accel = x[0];
                pole_angle_accel = x[1];
                if sim_index > 0{
                    cart_vel += (cart_accel + prev_cart_accel) / 2.0 * (time - time_prev);
                    cart_pos += (cart_vel + prev_cart_vel) / 2.0 * (time - time_prev);

                    pole_angle_vel += (pole_angle_accel + prev_pole_angle_accel) / 2.0 * (time - time_prev);
                    pole_angle_pos += (pole_angle_vel + prev_pole_angle_vel) / 2.0 * (time - time_prev);

                    force = 0.0; 

                    prev_cart_accel = cart_accel;
                    prev_cart_vel = cart_vel;
                    prev_cart_pos = cart_pos;
            
                    prev_pole_angle_accel = pole_angle_accel;
                    prev_pole_angle_vel = pole_angle_vel;
                    prev_pole_angle_pos = pole_angle_pos;      
    
                }
            }
            None=>{
                println!("Could not invert A Matrix");
            }
        }
        


        if time > 60.0 {
            solved = true;
        }

        time_prev = time;

        let output = format!("{:?},{:?},{:?},{:?},{:?},{:?},{:?}", time, cart_accel, cart_vel, cart_pos, pole_angle_accel, pole_angle_vel, pole_angle_pos);
        writeln!(file, "{}", output);

        time += time_delta;
        times.push(time);
        sim_index += 1;
    }

    println!("Solution Found!");
}